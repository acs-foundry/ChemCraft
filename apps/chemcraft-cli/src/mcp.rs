use std::io::{BufRead, BufReader, Write};
use std::net::{SocketAddr, TcpListener, TcpStream};
use std::path::Path;
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Duration;

use mol_admet_bridge::{
    AdmetCraftProvider, AdmetProvider, AdmetReportFormat, compute_admet_delta, predict_nmr,
};
use mol_engine::{
    BuiltinLibraryPreset, LibraryFilterConfig, PathSandbox, ScreeningFunnel,
    filter_and_prepare_library, load_builtin_library, stream_load_library_from_file,
};
use mol_graph::{
    MolGraph, SmartsQuery, analyze_molecule, attach_superatom_group, compute_stoichiometry,
    generate_3d_coords, mutate_molecule, parse_cdxml, parse_reaction_smiles, parse_smiles,
    resolve_chemical_name, write_cdxml, write_pdb, write_rxn_v2000, write_sdf_v2000,
    write_sdf_v2000_3d,
};
use mol_sar::{MorganFp, SarFilter, SarTable};
use mol_vector::{
    SceneStyle, build_mol_scene, export_sar_grid_eps, export_sar_grid_pdf_bytes,
    export_sar_grid_png_bytes, export_sar_grid_svg, generate_2d_coords, scene_to_eps,
    scene_to_pdf_bytes, scene_to_png_bytes, scene_to_svg,
};
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};

pub const MAX_MCP_MESSAGE_BYTES: usize = 1024 * 1024;
pub const MAX_MCP_INPUT_FILE_BYTES: usize = 2 * 1024 * 1024;
pub const MAX_MCP_MOLECULES: usize = 1_000;
pub const MAX_TCP_CONNECTIONS: usize = 32;
pub const TCP_SOCKET_TIMEOUT: Duration = Duration::from_secs(30);

#[derive(Debug, Deserialize)]
struct JsonRpcRequest {
    jsonrpc: Option<String>,
    id: Option<Value>,
    method: String,
    #[serde(default)]
    params: Value,
}

#[derive(Debug, Serialize)]
struct JsonRpcResponse {
    jsonrpc: &'static str,
    id: Value,
    #[serde(skip_serializing_if = "Option::is_none")]
    result: Option<Value>,
    #[serde(skip_serializing_if = "Option::is_none")]
    error: Option<JsonRpcError>,
}

#[derive(Debug, Serialize)]
struct JsonRpcError {
    code: i64,
    message: String,
}

struct SessionAuth {
    initialized: bool,
    required_token: Option<String>,
}

struct ConnGuard {
    active: Arc<AtomicUsize>,
}

impl Drop for ConnGuard {
    fn drop(&mut self) {
        self.active.fetch_sub(1, Ordering::AcqRel);
    }
}

pub fn run_mcp_stdio(sandbox: &PathSandbox) -> Result<(), Box<dyn std::error::Error>> {
    let stdin = std::io::stdin();
    let stdout = std::io::stdout();
    let reader = BufReader::new(stdin.lock());
    let mut writer = stdout.lock();
    serve_stream(reader, &mut writer, sandbox, None)
}

pub fn run_mcp_tcp(
    addr_str: &str,
    sandbox: &PathSandbox,
    allow_remote: bool,
    auth_token: Option<String>,
) -> Result<(), Box<dyn std::error::Error>> {
    let addr: SocketAddr = addr_str.parse()?;
    if !addr.ip().is_loopback() && !allow_remote {
        return Err(format!(
            "refusing to bind non-loopback TCP address {addr} without --allow-remote"
        )
        .into());
    }
    let token = auth_token
        .map(|t| t.trim().to_string())
        .filter(|t| !t.is_empty())
        .ok_or("TCP MCP mode requires non-empty --mcp-token or CHEMCRAFT_MCP_TOKEN")?;

    let listener = TcpListener::bind(addr)?;
    eprintln!("ChemCraft MCP TCP listening on {addr} (auth required)");
    let active_conns = Arc::new(AtomicUsize::new(0));
    let shared_sandbox = Arc::new(sandbox.clone());

    for stream_res in listener.incoming() {
        let Ok(stream) = stream_res else {
            continue;
        };
        handle_accepted_tcp(
            stream,
            Arc::clone(&shared_sandbox),
            token.clone(),
            Arc::clone(&active_conns),
            TCP_SOCKET_TIMEOUT,
        );
    }
    Ok(())
}

fn handle_accepted_tcp(
    mut stream: TcpStream,
    sandbox: Arc<PathSandbox>,
    token: String,
    active_conns: Arc<AtomicUsize>,
    timeout: Duration,
) {
    let prev = active_conns.fetch_add(1, Ordering::AcqRel);
    let guard = ConnGuard {
        active: Arc::clone(&active_conns),
    };
    if prev >= MAX_TCP_CONNECTIONS {
        let err_resp = JsonRpcResponse {
            jsonrpc: "2.0",
            id: Value::Null,
            result: None,
            error: Some(JsonRpcError {
                code: -32001,
                message: "server busy: too many concurrent TCP connections".to_string(),
            }),
        };
        let _ = write_json_line(&mut stream, &err_resp);
        drop(guard);
        return;
    }
    let _ = stream.set_read_timeout(Some(timeout));
    let _ = stream.set_write_timeout(Some(timeout));
    let Ok(reader_stream) = stream.try_clone() else {
        drop(guard);
        return;
    };
    std::thread::spawn(move || {
        let _guard = guard;
        let reader = BufReader::new(reader_stream);
        let mut writer = stream;
        let _ = serve_stream(reader, &mut writer, &sandbox, Some(token));
    });
}

#[cfg(test)]
pub fn serve_stream_for_test<R: BufRead, W: Write>(
    reader: R,
    writer: &mut W,
    sandbox: &PathSandbox,
    required_token: Option<String>,
) -> Result<(), Box<dyn std::error::Error>> {
    serve_stream(reader, writer, sandbox, required_token)
}

enum LineRead {
    Eof,
    Line(String),
    TooLarge,
}

fn read_limited_line<R: BufRead>(reader: &mut R, limit: usize) -> std::io::Result<LineRead> {
    let mut buf = Vec::new();
    let mut exceeded = false;
    loop {
        let avail = reader.fill_buf()?;
        if avail.is_empty() {
            if buf.is_empty() && !exceeded {
                return Ok(LineRead::Eof);
            }
            break;
        }
        if let Some(pos) = avail.iter().position(|&b| b == b'\n') {
            if !exceeded && buf.len() + pos <= limit {
                buf.extend_from_slice(&avail[..pos]);
            } else {
                exceeded = true;
            }
            reader.consume(pos + 1);
            break;
        }
        if !exceeded && buf.len() + avail.len() <= limit {
            buf.extend_from_slice(avail);
        } else {
            exceeded = true;
        }
        let consumed = avail.len();
        reader.consume(consumed);
    }
    if exceeded {
        return Ok(LineRead::TooLarge);
    }
    let s = String::from_utf8(buf)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
    Ok(LineRead::Line(s))
}

fn serve_stream<R: BufRead, W: Write>(
    mut reader: R,
    writer: &mut W,
    sandbox: &PathSandbox,
    required_token: Option<String>,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut auth = SessionAuth {
        initialized: false,
        required_token,
    };
    loop {
        match read_limited_line(&mut reader, MAX_MCP_MESSAGE_BYTES)? {
            LineRead::Eof => break,
            LineRead::TooLarge => {
                let err_resp = JsonRpcResponse {
                    jsonrpc: "2.0",
                    id: Value::Null,
                    result: None,
                    error: Some(JsonRpcError {
                        code: -32600,
                        message: format!(
                            "request exceeds maximum message size ({MAX_MCP_MESSAGE_BYTES} bytes)"
                        ),
                    }),
                };
                write_json_line(writer, &err_resp)?;
            }
            LineRead::Line(line) => {
                let trimmed = line.trim();
                if trimmed.is_empty() {
                    continue;
                }
                let raw_val: Value = match serde_json::from_str(trimmed) {
                    Ok(v) => v,
                    Err(e) => {
                        let resp = JsonRpcResponse {
                            jsonrpc: "2.0",
                            id: Value::Null,
                            result: None,
                            error: Some(JsonRpcError {
                                code: -32700,
                                message: format!("parse error: {e}"),
                            }),
                        };
                        write_json_line(writer, &resp)?;
                        continue;
                    }
                };
                let has_id = raw_val.get("id").is_some();
                let req: JsonRpcRequest = match serde_json::from_value(raw_val.clone()) {
                    Ok(r) => r,
                    Err(e) => {
                        let resp = JsonRpcResponse {
                            jsonrpc: "2.0",
                            id: raw_val.get("id").cloned().unwrap_or(Value::Null),
                            result: None,
                            error: Some(JsonRpcError {
                                code: -32600,
                                message: format!("invalid JSON-RPC request: {e}"),
                            }),
                        };
                        write_json_line(writer, &resp)?;
                        continue;
                    }
                };

                if req.jsonrpc.as_deref() != Some("2.0") {
                    if has_id {
                        let resp = JsonRpcResponse {
                            jsonrpc: "2.0",
                            id: req.id.unwrap_or(Value::Null),
                            result: None,
                            error: Some(JsonRpcError {
                                code: -32600,
                                message: "invalid request: 'jsonrpc' must equal \"2.0\""
                                    .to_string(),
                            }),
                        };
                        write_json_line(writer, &resp)?;
                    }
                    continue;
                }
                if !has_id {
                    // JSON-RPC 2.0 notification -> do not respond
                    continue;
                }
                let id = req.id.unwrap_or(Value::Null);
                let resp = match std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                    handle_request(id.clone(), &req.method, &req.params, sandbox, &mut auth)
                })) {
                    Ok(r) => r,
                    Err(_) => JsonRpcResponse {
                        jsonrpc: "2.0",
                        id,
                        result: None,
                        error: Some(JsonRpcError {
                            code: -32603,
                            message: "internal error during MCP request execution".to_string(),
                        }),
                    },
                };
                write_json_line(writer, &resp)?;
            }
        }
    }
    Ok(())
}

fn write_json_line<W: Write, T: Serialize>(
    writer: &mut W,
    val: &T,
) -> Result<(), Box<dyn std::error::Error>> {
    let s = serde_json::to_string(val)?;
    writer.write_all(s.as_bytes())?;
    writer.write_all(b"\n")?;
    writer.flush()?;
    Ok(())
}

fn constant_time_eq(a: &[u8], b: &[u8]) -> bool {
    let max_len = a.len().max(b.len());
    let mut diff = a.len() ^ b.len();
    for i in 0..max_len {
        let ba = a.get(i).copied().unwrap_or(0);
        let bb = b.get(i).copied().unwrap_or(0);
        diff |= usize::from(ba ^ bb);
    }
    diff == 0
}

fn handle_request(
    id: Value,
    method: &str,
    params: &Value,
    sandbox: &PathSandbox,
    auth: &mut SessionAuth,
) -> JsonRpcResponse {
    if method == "initialize" {
        if let Some(ref req_tok) = auth.required_token {
            let provided = params.get("token").and_then(Value::as_str).unwrap_or("");
            if !constant_time_eq(provided.as_bytes(), req_tok.as_bytes()) {
                return JsonRpcResponse {
                    jsonrpc: "2.0",
                    id,
                    result: None,
                    error: Some(JsonRpcError {
                        code: -32001,
                        message: "unauthorized: invalid or missing token".to_string(),
                    }),
                };
            }
        }
        auth.initialized = true;
        return JsonRpcResponse {
            jsonrpc: "2.0",
            id,
            result: Some(json!({
                "protocolVersion": "2024-11-05",
                "serverInfo": { "name": "chemcraft-mcp", "version": env!("CARGO_PKG_VERSION") },
                "capabilities": { "tools": {} }
            })),
            error: None,
        };
    }

    if !auth.initialized {
        return JsonRpcResponse {
            jsonrpc: "2.0",
            id,
            result: None,
            error: Some(JsonRpcError {
                code: -32002,
                message: "server not initialized: call 'initialize' first".to_string(),
            }),
        };
    }

    match method {
        "ping" => JsonRpcResponse {
            jsonrpc: "2.0",
            id,
            result: Some(json!({})),
            error: None,
        },
        "tools/list" => JsonRpcResponse {
            jsonrpc: "2.0",
            id,
            result: Some(json!({ "tools": mcp_tool_catalog() })),
            error: None,
        },
        "tools/call" => {
            let Some(name) = params.get("name").and_then(Value::as_str) else {
                return JsonRpcResponse {
                    jsonrpc: "2.0",
                    id,
                    result: None,
                    error: Some(JsonRpcError {
                        code: -32602,
                        message: "invalid params: missing 'name'".to_string(),
                    }),
                };
            };
            let args = params
                .get("arguments")
                .cloned()
                .unwrap_or_else(|| json!({}));
            if !args.is_object() {
                return JsonRpcResponse {
                    jsonrpc: "2.0",
                    id,
                    result: None,
                    error: Some(JsonRpcError {
                        code: -32602,
                        message: "invalid params: 'arguments' must be a JSON object".to_string(),
                    }),
                };
            }
            match dispatch_tool(name, &args, sandbox) {
                Ok(res) => JsonRpcResponse {
                    jsonrpc: "2.0",
                    id,
                    result: Some(res),
                    error: None,
                },
                Err(invalid_param_msg) => JsonRpcResponse {
                    jsonrpc: "2.0",
                    id,
                    result: None,
                    error: Some(JsonRpcError {
                        code: -32602,
                        message: invalid_param_msg,
                    }),
                },
            }
        }
        other => JsonRpcResponse {
            jsonrpc: "2.0",
            id,
            result: None,
            error: Some(JsonRpcError {
                code: -32601,
                message: format!("method not found: {other}"),
            }),
        },
    }
}

fn mcp_tool_catalog() -> Value {
    json!([
        {
            "name": "draw_from_smiles",
            "description": "Generate clean 2D coordinates, valence validation, structural key, and vector SVG preview from SMILES.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "smiles": { "type": "string" },
                    "name": { "type": "string" }
                },
                "required": ["smiles"]
            },
            "outputSchema": {
                "type": "object",
                "properties": {
                    "canonical_smiles": { "type": "string" },
                    "formula": { "type": "string" },
                    "warnings": { "type": "array" },
                    "svg": { "type": "string" }
                }
            }
        },
        {
            "name": "mutate_scaffold",
            "description": "Modify a scaffold/functional group (e.g. 'tetrazole' bioisostere replacing -COOH, or 'fluoro', 'chloro', 'methyl', 'trifluoromethyl', 'cyano') and compute the AdmetCraft delta.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "smiles": { "type": "string" },
                    "group": { "type": "string" },
                    "atom_index": { "type": "integer", "minimum": 0 },
                    "name": { "type": "string" }
                },
                "required": ["smiles", "group"]
            },
            "outputSchema": {
                "type": "object",
                "properties": {
                    "old_smiles": { "type": "string" },
                    "new_smiles": { "type": "string" },
                    "delta_logp": { "type": "number" },
                    "delta_tpsa": { "type": "number" },
                    "new_gi": { "type": "string" }
                }
            }
        },
        {
            "name": "search_substructure",
            "description": "Filter an SDF/CSV/DWAR file inside --read-root using a SMARTS substructure pattern and optional Tanimoto threshold.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "file": { "type": "string" },
                    "smarts": { "type": "string" },
                    "reference_smiles": { "type": "string" },
                    "min_tanimoto": { "type": "number", "minimum": 0.0, "maximum": 1.0 },
                    "output_sdf": { "type": "string" }
                },
                "required": ["file", "smarts"]
            },
            "outputSchema": {
                "type": "object",
                "properties": {
                    "total_loaded": { "type": "integer" },
                    "matched_count": { "type": "integer" },
                    "matches": { "type": "array" }
                }
            }
        },
        {
            "name": "export_vector_figure",
            "description": "Export a single molecule or multi-molecule SAR grid to SVG, EPS, PDF, or PNG inside --write-root.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "smiles": { "type": "array", "items": { "type": "string" } },
                    "output_path": { "type": "string" },
                    "format": { "type": "string", "enum": ["svg", "eps", "pdf", "png"] },
                    "columns": { "type": "integer", "minimum": 1, "maximum": 8 }
                },
                "required": ["smiles", "output_path", "format"]
            },
            "outputSchema": {
                "type": "object",
                "properties": {
                    "written_path": { "type": "string" },
                    "count": { "type": "integer" },
                    "overlap_count": { "type": "integer" }
                }
            }
        },
        {
            "name": "mol.analyze_admet",
            "description": "Invoke the in-process AdmetCraft engine on molecules (SMILES or file inside --read-root) and optionally write a CSV/JSON/PDF report inside --write-root.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "smiles": { "type": "array", "items": { "type": "string" } },
                    "output_report": { "type": "string" },
                    "format": { "type": "string", "enum": ["csv", "json", "pdf"] }
                },
                "required": ["smiles"]
            },
            "outputSchema": {
                "type": "object",
                "properties": {
                    "count": { "type": "integer" },
                    "evaluations": { "type": "array" }
                }
            }
        },
        {
            "name": "export_3d",
            "description": "Generate 3D coordinates from SMILES and write PDB or 3D SDF V2000 inside --write-root.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "smiles": { "type": "string" },
                    "output_path": { "type": "string" },
                    "format": { "type": "string", "enum": ["pdb", "sdf3d"] }
                },
                "required": ["output_path", "format"]
            }
        },
        {
            "name": "export_chimerax",
            "description": "Export a self-contained UCSF ChimeraX 3D bundle (.pdb, .sdf, and .cxc script) inside --write-root without leaving the sandbox.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "smiles": { "type": "string" },
                    "stem": { "type": "string" }
                },
                "required": ["stem"]
            }
        },
        {
            "name": "screen_library",
            "description": "Stream/filter a built-in ('fda'|'fragments'|'kinase') or sandboxed --read-root chemical library (.sdf/.sdf.gz/.smi/.smi.gz/.csv/.dwar) using an ADMET/SMARTS/MaxMin funnel, ionize at pH 7.4, generate 3D conformers, and write a DockCraft-ready 3D SDF inside --write-root.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "preset": { "type": "string", "enum": ["fda", "fragments", "kinase"] },
                    "file": { "type": "string" },
                    "funnel": { "type": "string", "enum": ["all", "ro5", "lead", "ro3", "bbb"] },
                    "smarts": { "type": "string" },
                    "diverse_n": { "type": "integer", "minimum": 1, "maximum": 500 },
                    "protonate_ph74": { "type": "boolean" },
                    "output_sdf": { "type": "string" }
                }
            }
        },
        {
            "name": "analyze_molecule",
            "description": "ChemDraw Analysis Window: exact monoisotopic mass, average MW, charge-annotated Hill formula, elemental analysis (%C/%H/...), DBE, and ESI+/ESI- m/z adducts.",
            "inputSchema": {
                "type": "object",
                "properties": { "smiles": { "type": "string" } },
                "required": ["smiles"]
            }
        },
        {
            "name": "predict_nmr",
            "description": "Predict 1H and 13C NMR chemical shifts (ppm), multiplicities, and integrations for a SMILES structure.",
            "inputSchema": {
                "type": "object",
                "properties": { "smiles": { "type": "string" } },
                "required": ["smiles"]
            }
        },
        {
            "name": "name_to_structure",
            "description": "Convert a common drug/solvent/reagent name or systematic IUPAC name into canonical SMILES and exact mass.",
            "inputSchema": {
                "type": "object",
                "properties": { "name": { "type": "string" } },
                "required": ["name"]
            }
        },
        {
            "name": "cdxml_io",
            "description": "Export SMILES to ChemDraw `.cdxml` inside --write-root, or parse an input `.cdxml` file inside --read-root.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "smiles": { "type": "string" },
                    "output_cdxml": { "type": "string" },
                    "input_cdxml": { "type": "string" }
                }
            }
        },
        {
            "name": "rxn_stoichiometry",
            "description": "Parse Reaction SMILES, compute atom economy & theoretical/isolated yield, and optionally export MDL `.rxn` V2000 inside --write-root.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "reaction_smiles": { "type": "string" },
                    "mmol": { "type": "number" },
                    "actual_mg": { "type": "number" },
                    "output_rxn": { "type": "string" }
                },
                "required": ["reaction_smiles"]
            }
        },
        {
            "name": "attach_superatom",
            "description": "Attach a ChemDraw superatom/protecting group (Boc, Fmoc, Cbz, Ph, Bn, CF3, CO2Me, SO2Me, ...) onto a SMILES.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "smiles": { "type": "string" },
                    "group": { "type": "string" },
                    "atom_index": { "type": "integer" }
                },
                "required": ["smiles", "group"]
            }
        },
        {
            "name": "plan_synthesis",
            "description": "Run Buy-Before-Make building-block lookup + Retro* multi-step retrosynthesis planning with reagents, stoichiometry, chemoselectivity warnings, and optional ELN Markdown export.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "smiles": { "type": "string" },
                    "name": { "type": "string" },
                    "max_steps": { "type": "integer", "minimum": 1, "maximum": 6 },
                    "top": { "type": "integer", "minimum": 1, "maximum": 10 },
                    "scale_mg": { "type": "number", "minimum": 1.0 },
                    "out_md_path": { "type": "string" }
                },
                "required": ["smiles"]
            }
        },
    ])
}

fn dispatch_tool(name: &str, args: &Value, sandbox: &PathSandbox) -> Result<Value, String> {
    match name {
        "draw_from_smiles" => tool_draw_from_smiles(args),
        "mutate_scaffold" => tool_mutate_scaffold(args),
        "search_substructure" => tool_search_substructure(args, sandbox),
        "export_vector_figure" => tool_export_vector_figure(args, sandbox),
        "mol.analyze_admet" => tool_analyze_admet(args, sandbox),
        "export_3d" => tool_export_3d(args, sandbox),
        "export_chimerax" => tool_export_chimerax(args, sandbox),
        "screen_library" => tool_screen_library(args, sandbox),
        "analyze_molecule" => tool_analyze_molecule(args),
        "predict_nmr" => tool_predict_nmr(args),
        "name_to_structure" => tool_name_to_structure(args),
        "cdxml_io" => tool_cdxml_io(args, sandbox),
        "rxn_stoichiometry" => tool_rxn_stoichiometry(args, sandbox),
        "attach_superatom" => tool_attach_superatom(args),
        "plan_synthesis" => tool_plan_synthesis(args, sandbox),
        other => Err(format!("unknown tool '{other}'")),
    }
}

fn tool_ok(summary: String, content: Value) -> Value {
    json!({
        "isError": false,
        "content": [{ "type": "text", "text": summary }],
        "structuredContent": content
    })
}

fn tool_err(msg: impl Into<String>) -> Value {
    let m = msg.into();
    json!({
        "isError": true,
        "content": [{ "type": "text", "text": m }]
    })
}

fn tool_draw_from_smiles(args: &Value) -> Result<Value, String> {
    let smi = args
        .get("smiles")
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| "missing required non-empty 'smiles'".to_string())?;
    let name = args.get("name").and_then(Value::as_str);

    let mut mol = match parse_smiles(smi) {
        Ok(m) => m,
        Err(e) => return Ok(tool_err(format!("SMILES parse error: {e}"))),
    };
    if let Some(n) = name {
        mol.name = Some(n.to_string());
    }
    generate_2d_coords(&mut mol);
    let scene = build_mol_scene(&mol, &SceneStyle::default(), &[]);
    let svg = scene_to_svg(&scene);
    let can = mol.canonical_smiles();
    let key = mol.structural_key();
    let formula = mol.formula();
    let warnings = mol.warnings.clone();
    let coords: Vec<[f64; 2]> = mol
        .atoms
        .iter()
        .map(|a| a.coords_2d.unwrap_or([0.0, 0.0]))
        .collect();

    Ok(tool_ok(
        format!("Drew {can} ({formula})"),
        json!({
            "canonical_smiles": can,
            "formula": formula,
            "structural_key": key,
            "warnings": warnings,
            "coords_2d": coords,
            "label_overlaps": scene.label_overlap_count(),
            "svg": svg,
        }),
    ))
}

fn tool_mutate_scaffold(args: &Value) -> Result<Value, String> {
    let smi = args
        .get("smiles")
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| "missing required 'smiles'".to_string())?;
    let group = args
        .get("group")
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| "missing required 'group'".to_string())?;
    let atom_index = match args.get("atom_index") {
        None | Some(Value::Null) => None,
        Some(v) => Some(
            usize::try_from(
                v.as_u64()
                    .ok_or_else(|| "'atom_index' must be a non-negative integer".to_string())?,
            )
            .map_err(|_| "'atom_index' overflow".to_string())?,
        ),
    };
    let name = args.get("name").and_then(Value::as_str);

    let mol = match parse_smiles(smi) {
        Ok(m) => m,
        Err(e) => return Ok(tool_err(format!("SMILES parse error: {e}"))),
    };
    let mut mutated = match mutate_molecule(&mol, atom_index, group) {
        Ok(m) => m,
        Err(e) => return Ok(tool_err(format!("mutation error: {e}"))),
    };
    generate_2d_coords(&mut mutated);
    let new_smi = mutated.canonical_smiles();
    let old_smi = mol.canonical_smiles();
    let (delta, eval) =
        match compute_admet_delta(&AdmetCraftProvider::default(), &old_smi, &new_smi, name) {
            Ok(res) => res,
            Err(e) => return Ok(tool_err(format!("ADMET evaluation error: {e}"))),
        };
    let scene = build_mol_scene(&mutated, &SceneStyle::default(), &[]);
    let svg = scene_to_svg(&scene);

    Ok(tool_ok(
        format!(
            "Mutated {old_smi} -> {new_smi} (ΔlogP={:+.2}, ΔTPSA={:+.2}, GI={})",
            delta.delta_logp, delta.delta_tpsa, delta.new_gi
        ),
        json!({
            "old_smiles": delta.old_smiles,
            "new_smiles": delta.new_smiles,
            "delta_mw": delta.delta_mw,
            "delta_logp": delta.delta_logp,
            "delta_tpsa": delta.delta_tpsa,
            "delta_esol_log_s": delta.delta_esol_log_s,
            "old_gi": delta.old_gi,
            "new_gi": delta.new_gi,
            "old_bbb": delta.old_bbb,
            "new_bbb": delta.new_bbb,
            "new_admet": eval,
            "svg": svg,
        }),
    ))
}

fn tool_search_substructure(args: &Value, sandbox: &PathSandbox) -> Result<Value, String> {
    let file_str = args
        .get("file")
        .and_then(Value::as_str)
        .ok_or_else(|| "missing required 'file'".to_string())?;
    let smarts = args
        .get("smarts")
        .and_then(Value::as_str)
        .ok_or_else(|| "missing required 'smarts'".to_string())?;
    if SmartsQuery::parse(smarts).is_err() {
        return Err(format!("invalid SMARTS pattern '{smarts}'"));
    }

    let ref_smi = args.get("reference_smiles").and_then(Value::as_str);
    let min_tanimoto = args
        .get("min_tanimoto")
        .and_then(Value::as_f64)
        .unwrap_or(0.0);
    if !(0.0..=1.0).contains(&min_tanimoto) {
        return Err("'min_tanimoto' must be inside [0.0, 1.0]".to_string());
    }

    let (resolved, text) =
        match sandbox.read_to_string_limited(Path::new(file_str), MAX_MCP_INPUT_FILE_BYTES) {
            Ok(res) => res,
            Err(e) => return Ok(tool_err(format!("read sandbox error: {e}"))),
        };

    let mut table = SarTable::new();
    let ext = resolved
        .extension()
        .and_then(|s| s.to_str())
        .unwrap_or("")
        .to_ascii_lowercase();
    let load_res = if ext == "dwar" {
        table.load_dwar_text_limited(&text, MAX_MCP_MOLECULES)
    } else if ext == "csv" {
        table.load_csv_text_limited(&text, MAX_MCP_MOLECULES)
    } else {
        table.load_sdf_text_limited(&text, MAX_MCP_MOLECULES)
    };
    if let Err(e) = load_res {
        return Ok(tool_err(format!("parse error: {e}")));
    }

    let filter = SarFilter {
        smarts_substructure: Some(smarts.to_string()),
        reference_smiles: ref_smi.map(str::to_string),
        min_tanimoto,
        ..SarFilter::default()
    };
    let hit_indices = table
        .matching_indices(&filter)
        .map_err(|e| format!("filter error: {e}"))?;

    let ref_fp = ref_smi
        .and_then(|s| parse_smiles(s).ok())
        .map(|m| MorganFp::ecfp4(&m));

    let hits_json: Vec<Value> = hit_indices
        .iter()
        .map(|&idx| {
            let r = &table.rows[idx];
            let sim = ref_fp.as_ref().map(|rfp| r.fp.tanimoto(rfp));
            json!({
                "name": r.name,
                "canonical_smiles": r.canonical_smiles,
                "formula": r.formula,
                "mw": r.mw,
                "tanimoto": sim,
            })
        })
        .collect();

    if let Some(out_sdf) = args.get("output_sdf").and_then(Value::as_str) {
        let refs: Vec<&MolGraph> = hit_indices
            .iter()
            .map(|&idx| &table.rows[idx].mol)
            .collect();
        let sdf_out = write_sdf_v2000(&refs);
        match sandbox.open_write(Path::new(out_sdf)) {
            Ok((_, mut f)) => {
                if let Err(e) = f.write_all(sdf_out.as_bytes()) {
                    return Ok(tool_err(format!("write failed: {e}")));
                }
            }
            Err(e) => return Ok(tool_err(format!("write sandbox error: {e}"))),
        }
    }

    Ok(tool_ok(
        format!(
            "Matched {}/{} compounds against '{}'",
            hit_indices.len(),
            table.rows.len(),
            smarts
        ),
        json!({
            "total_loaded": table.rows.len(),
            "matched_count": hit_indices.len(),
            "matches": hits_json,
        }),
    ))
}

fn tool_export_vector_figure(args: &Value, sandbox: &PathSandbox) -> Result<Value, String> {
    let smi_arr = args
        .get("smiles")
        .and_then(Value::as_array)
        .filter(|a| !a.is_empty())
        .ok_or_else(|| "missing non-empty 'smiles' array".to_string())?;
    if smi_arr.len() > MAX_MCP_MOLECULES {
        return Err(format!("too many molecules (max {MAX_MCP_MOLECULES})"));
    }
    let out_str = args
        .get("output_path")
        .and_then(Value::as_str)
        .ok_or_else(|| "missing required 'output_path'".to_string())?;
    let fmt = args
        .get("format")
        .and_then(Value::as_str)
        .ok_or_else(|| "missing required 'format' (svg|eps|pdf|png)".to_string())?
        .to_ascii_lowercase();
    if !matches!(fmt.as_str(), "svg" | "eps" | "pdf" | "png") {
        return Err(format!("unsupported format '{fmt}'"));
    }
    let cols = args
        .get("columns")
        .and_then(Value::as_u64)
        .map(|v| (v as usize).clamp(1, 8))
        .unwrap_or(3);

    let mut mols = Vec::with_capacity(smi_arr.len());
    for (idx, v) in smi_arr.iter().enumerate() {
        let raw = v
            .as_str()
            .ok_or_else(|| format!("smiles[{idx}] must be a string"))?;
        let mut parts = raw.trim().splitn(2, char::is_whitespace);
        let smi = parts.next().unwrap_or("");
        let label = parts
            .next()
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .map(str::to_string)
            .unwrap_or_else(|| format!("Mol_{}", idx + 1));
        let mut m = match parse_smiles(smi) {
            Ok(m) => m,
            Err(e) => return Ok(tool_err(format!("invalid SMILES at index {idx}: {e}"))),
        };
        m.name = Some(label);
        generate_2d_coords(&mut m);
        mols.push(m);
    }

    let (resolved, mut file) = match sandbox.open_write(Path::new(out_str)) {
        Ok(pair) => pair,
        Err(e) => return Ok(tool_err(format!("write sandbox error: {e}"))),
    };

    let mut total_overlaps = 0usize;
    let bytes: Vec<u8> = if mols.len() == 1 {
        let scene = build_mol_scene(&mols[0], &SceneStyle::default(), &[]);
        total_overlaps = scene.label_overlap_count();
        match fmt.as_str() {
            "svg" => scene_to_svg(&scene).into_bytes(),
            "eps" => scene_to_eps(&scene).into_bytes(),
            "pdf" => scene_to_pdf_bytes(&scene).map_err(|e| e.to_string())?,
            "png" => scene_to_png_bytes(&scene, 2.0).map_err(|e| e.to_string())?,
            _ => unreachable!(),
        }
    } else {
        let refs: Vec<(&str, &MolGraph)> = mols
            .iter()
            .map(|m| (m.name.as_deref().unwrap_or("Mol"), m))
            .collect();
        match fmt.as_str() {
            "svg" => export_sar_grid_svg(&refs, cols, 250.0, 210.0, false).into_bytes(),
            "pdf" => {
                export_sar_grid_pdf_bytes(&refs, cols, 250.0, 210.0).map_err(|e| e.to_string())?
            }
            "eps" => export_sar_grid_eps(&refs, cols, 250.0, 210.0, false).into_bytes(),
            "png" => export_sar_grid_png_bytes(&refs, cols, 250.0, 210.0, false, 2.0)
                .map_err(|e| e.to_string())?,
            _ => unreachable!(),
        }
    };

    if let Err(e) = file.write_all(&bytes) {
        return Ok(tool_err(format!("export write failed: {e}")));
    }

    Ok(tool_ok(
        format!(
            "Exported {} structure(s) as {fmt} to {}",
            mols.len(),
            resolved.display()
        ),
        json!({
            "written_path": resolved.display().to_string(),
            "format": fmt,
            "count": mols.len(),
            "overlap_count": total_overlaps,
            "bytes": bytes.len(),
        }),
    ))
}

fn tool_analyze_admet(args: &Value, sandbox: &PathSandbox) -> Result<Value, String> {
    let smi_arr = args
        .get("smiles")
        .and_then(Value::as_array)
        .filter(|a| !a.is_empty())
        .ok_or_else(|| "missing non-empty 'smiles' array".to_string())?;
    if smi_arr.len() > MAX_MCP_MOLECULES {
        return Err(format!("too many molecules (max {MAX_MCP_MOLECULES})"));
    }

    let mut items = Vec::with_capacity(smi_arr.len());
    for (idx, v) in smi_arr.iter().enumerate() {
        let raw = v
            .as_str()
            .ok_or_else(|| format!("smiles[{idx}] must be a string"))?;
        let mut parts = raw.trim().splitn(2, char::is_whitespace);
        let smi = parts.next().unwrap_or("").to_string();
        let name = parts
            .next()
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .map(str::to_string)
            .or_else(|| Some(format!("Mol_{}", idx + 1)));
        items.push((smi, name));
    }

    let provider = AdmetCraftProvider::default();
    let evals: Vec<_> = provider
        .evaluate_batch(&items)
        .into_iter()
        .filter_map(Result::ok)
        .collect();
    if evals.is_empty() {
        return Ok(tool_err(
            "no valid molecules could be evaluated by AdmetCraft",
        ));
    }

    let mut report_path_written = None;
    if let Some(out_str) = args.get("output_report").and_then(Value::as_str) {
        let fmt_str = args
            .get("format")
            .and_then(Value::as_str)
            .ok_or_else(|| "'format' is required when 'output_report' is set".to_string())?
            .to_ascii_lowercase();
        let fmt = match fmt_str.as_str() {
            "csv" => AdmetReportFormat::Csv,
            "json" => AdmetReportFormat::Json,
            "pdf" => AdmetReportFormat::Pdf,
            _ => return Err(format!("unsupported ADMET report format '{fmt_str}'")),
        };
        let (resolved, mut file) = match sandbox.open_write(Path::new(out_str)) {
            Ok(p) => p,
            Err(e) => return Ok(tool_err(format!("write sandbox error: {e}"))),
        };
        if let Err(e) = provider.export_report_to_writer(&evals, &mut file, fmt) {
            return Ok(tool_err(format!("report error: {e}")));
        }
        report_path_written = Some(resolved.display().to_string());
    }

    Ok(tool_ok(
        format!(
            "Evaluated {} molecule(s) via AdmetCraft engine",
            evals.len()
        ),
        json!({
            "count": evals.len(),
            "report_written": report_path_written,
            "evaluations": evals,
        }),
    ))
}

fn resolve_mol_from_args(args: &Value) -> Result<MolGraph, String> {
    let smi = args
        .get("smiles")
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| "provide non-empty 'smiles'".to_string())?;
    let mut mol = parse_smiles(smi).map_err(|e| format!("SMILES parse error: {e}"))?;
    generate_3d_coords(&mut mol);
    Ok(mol)
}

fn tool_export_3d(args: &Value, sandbox: &PathSandbox) -> Result<Value, String> {
    let mol = match resolve_mol_from_args(args) {
        Ok(m) => m,
        Err(msg) => return Ok(tool_err(msg)),
    };
    let out_str = args
        .get("output_path")
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| "missing required 'output_path'".to_string())?;
    let fmt = args
        .get("format")
        .and_then(Value::as_str)
        .unwrap_or("pdb")
        .to_ascii_lowercase();
    let payload = match fmt.as_str() {
        "pdb" => write_pdb(&mol),
        "sdf3d" => write_sdf_v2000_3d(&[&mol]),
        _ => return Err(format!("unsupported 3D format '{fmt}'")),
    };
    let (resolved, mut f) = match sandbox.open_write(Path::new(out_str)) {
        Ok(p) => p,
        Err(e) => return Ok(tool_err(format!("write sandbox error: {e}"))),
    };
    f.write_all(payload.as_bytes()).map_err(|e| e.to_string())?;
    Ok(tool_ok(
        format!("Exported 3D {fmt} -> {}", resolved.display()),
        json!({
            "written_path": resolved.display().to_string(),
            "format": fmt,
            "atom_count": mol.atoms.len(),
            "bytes": payload.len(),
        }),
    ))
}

fn tool_export_chimerax(args: &Value, sandbox: &PathSandbox) -> Result<Value, String> {
    let mol = match resolve_mol_from_args(args) {
        Ok(m) => m,
        Err(msg) => return Ok(tool_err(msg)),
    };
    let raw_stem = args
        .get("stem")
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| "missing required non-empty 'stem'".to_string())?;
    let safe_stem: String = raw_stem
        .chars()
        .filter(|c| c.is_ascii_alphanumeric() || *c == '_' || *c == '-')
        .take(40)
        .collect();
    if safe_stem.is_empty() {
        return Err("stem contains no safe ASCII characters".to_string());
    }

    let is_peptide = mol
        .atoms
        .iter()
        .any(|a| a.pdb_meta.as_ref().is_some_and(|m| !m.is_hetatm));
    let pdb_rel = format!("{safe_stem}.pdb");
    let sdf_rel = format!("{safe_stem}.sdf");
    let cxc_rel = format!("{safe_stem}.cxc");

    let (pdb_resolved, mut pdb_f) = match sandbox.open_write(Path::new(&pdb_rel)) {
        Ok(p) => p,
        Err(e) => return Ok(tool_err(format!("write sandbox error: {e}"))),
    };
    pdb_f
        .write_all(write_pdb(&mol).as_bytes())
        .map_err(|e| e.to_string())?;

    let (sdf_resolved, mut sdf_f) = match sandbox.open_write(Path::new(&sdf_rel)) {
        Ok(p) => p,
        Err(e) => return Ok(tool_err(format!("write sandbox error: {e}"))),
    };
    sdf_f
        .write_all(write_sdf_v2000_3d(&[&mol]).as_bytes())
        .map_err(|e| e.to_string())?;

    let cxc = if is_peptide {
        format!(
            "# ChemCraft -> UCSF ChimeraX\nopen {pdb_rel}\npreset \"overall look\" publication\ncartoon\nshow atoms\nstyle stick\ncolor byhetero\nhbonds reveal true\nlighting soft\nset bgColor white\nview\n"
        )
    } else {
        format!(
            "# ChemCraft -> UCSF ChimeraX\nopen {sdf_rel}\npreset \"overall look\" publication\nstyle ball\ncolor byhetero\nlighting soft\nset bgColor white\nview\n"
        )
    };
    let (cxc_resolved, mut cxc_f) = match sandbox.open_write(Path::new(&cxc_rel)) {
        Ok(p) => p,
        Err(e) => return Ok(tool_err(format!("write sandbox error: {e}"))),
    };
    cxc_f.write_all(cxc.as_bytes()).map_err(|e| e.to_string())?;

    Ok(tool_ok(
        format!("Exported ChimeraX bundle -> {}", cxc_resolved.display()),
        json!({
            "script_path": cxc_resolved.display().to_string(),
            "pdb_path": pdb_resolved.display().to_string(),
            "sdf_path": sdf_resolved.display().to_string(),
            "is_peptide": is_peptide,
        }),
    ))
}

fn tool_screen_library(args: &Value, sandbox: &PathSandbox) -> Result<Value, String> {
    let (label, mols) = if let Some(rel_file) = args
        .get("file")
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
    {
        let (resolved, file) = match sandbox.open_read(Path::new(rel_file)) {
            Ok(p) => p,
            Err(e) => return Ok(tool_err(format!("read sandbox error: {e}"))),
        };
        let hint = resolved
            .file_name()
            .map(|s| s.to_string_lossy().to_string())
            .unwrap_or_else(|| rel_file.to_string());
        let loaded = match stream_load_library_from_file(file, &hint, 5_000) {
            Ok(v) => v,
            Err(e) => return Ok(tool_err(format!("library parse error: {e}"))),
        };
        (rel_file.to_string(), loaded)
    } else {
        let pre_str = args.get("preset").and_then(Value::as_str).unwrap_or("fda");
        let preset = BuiltinLibraryPreset::from_str_lossy(pre_str)
            .ok_or_else(|| format!("unknown library preset '{pre_str}'"))?;
        (preset.label().to_string(), load_builtin_library(preset))
    };

    let funnel_str = args.get("funnel").and_then(Value::as_str).unwrap_or("ro5");
    let funnel = ScreeningFunnel::from_str_lossy(funnel_str)
        .ok_or_else(|| format!("invalid funnel '{funnel_str}'"))?;
    let cfg = LibraryFilterConfig {
        funnel,
        smarts_substructure: args
            .get("smarts")
            .and_then(Value::as_str)
            .map(str::to_string),
        diversity_pick_n: args
            .get("diverse_n")
            .and_then(Value::as_u64)
            .map(|v| (v as usize).clamp(1, 500)),
        max_hits: 500,
        protonate_ph74: args
            .get("protonate_ph74")
            .and_then(Value::as_bool)
            .unwrap_or(true),
        generate_3d: true,
        ..LibraryFilterConfig::default()
    };

    let res = match filter_and_prepare_library(&label, mols, &cfg) {
        Ok(r) => r,
        Err(e) => return Ok(tool_err(format!("library screening failed: {e}"))),
    };

    let mut sdf_written = None;
    if let Some(out_rel) = args.get("output_sdf").and_then(Value::as_str) {
        let refs: Vec<&MolGraph> = res.hits.iter().map(|h| &h.mol).collect();
        let sdf = write_sdf_v2000_3d(&refs);
        let (resolved, mut f) = match sandbox.open_write(Path::new(out_rel)) {
            Ok(p) => p,
            Err(e) => return Ok(tool_err(format!("write sandbox error: {e}"))),
        };
        f.write_all(sdf.as_bytes()).map_err(|e| e.to_string())?;
        sdf_written = Some(resolved.display().to_string());
    }

    let hit_summaries: Vec<Value> = res
        .hits
        .iter()
        .take(50)
        .map(|h| {
            json!({
                "name": h.name,
                "canonical_smiles": h.canonical_smiles,
                "formula": h.eval.physicochemical.formula,
                "mw": h.eval.physicochemical.molecular_weight,
                "clogp": h.eval.lipophilicity.consensus_log_p,
                "tpsa": h.eval.physicochemical.tpsa,
                "lipinski_violations": h.eval.druglikeness.lipinski.violations,
            })
        })
        .collect();

    Ok(tool_ok(
        format!(
            "Screened '{}': {}/{} passed funnel -> {} selected in {} ms",
            res.source_label,
            res.passed_funnel,
            res.total_loaded,
            res.hits.len(),
            res.elapsed_ms
        ),
        json!({
            "source": res.source_label,
            "total_loaded": res.total_loaded,
            "passed_funnel": res.passed_funnel,
            "selected_count": res.hits.len(),
            "elapsed_ms": res.elapsed_ms,
            "output_sdf": sdf_written,
            "hits": hit_summaries,
        }),
    ))
}

fn tool_analyze_molecule(args: &Value) -> Result<Value, String> {
    let smiles = args
        .get("smiles")
        .and_then(Value::as_str)
        .ok_or_else(|| "missing 'smiles'".to_string())?;
    let mol = parse_smiles(smiles).map_err(|e| e.to_string())?;
    let analysis = analyze_molecule(&mol);
    Ok(tool_ok(
        format!(
            "{} | Exact Mass {:.4} | MW {:.2} | DBE {:.1} | EA: {}",
            analysis.formula_with_charge,
            analysis.exact_mass,
            analysis.molecular_weight,
            analysis.dbe,
            analysis.elemental_analysis_string
        ),
        serde_json::to_value(&analysis).map_err(|e| e.to_string())?,
    ))
}

fn tool_predict_nmr(args: &Value) -> Result<Value, String> {
    let smiles = args
        .get("smiles")
        .and_then(Value::as_str)
        .ok_or_else(|| "missing 'smiles'".to_string())?;
    let mol = parse_smiles(smiles).map_err(|e| e.to_string())?;
    let spec = predict_nmr(&mol);
    Ok(tool_ok(
        format!(
            "Predicted {} 1H NMR groups & {} 13C NMR groups",
            spec.h1_peaks.len(),
            spec.c13_peaks.len()
        ),
        serde_json::to_value(&spec).map_err(|e| e.to_string())?,
    ))
}

fn tool_name_to_structure(args: &Value) -> Result<Value, String> {
    let name = args
        .get("name")
        .and_then(Value::as_str)
        .ok_or_else(|| "missing 'name'".to_string())?;
    let mol = resolve_chemical_name(name).map_err(|e| e.to_string())?;
    let analysis = analyze_molecule(&mol);
    let smi = mol.canonical_smiles();
    Ok(tool_ok(
        format!(
            "Resolved '{name}' -> {smi} ({})",
            analysis.formula_with_charge
        ),
        json!({
            "name": name,
            "canonical_smiles": smi,
            "analysis": analysis,
        }),
    ))
}

fn tool_cdxml_io(args: &Value, sandbox: &PathSandbox) -> Result<Value, String> {
    if let Some(smi) = args.get("smiles").and_then(Value::as_str) {
        let mut mol = parse_smiles(smi).map_err(|e| e.to_string())?;
        generate_2d_coords(&mut mol);
        let xml = write_cdxml(&mol);
        let mut written = None;
        if let Some(out_rel) = args.get("output_cdxml").and_then(Value::as_str) {
            let (resolved, mut f) = match sandbox.open_write(Path::new(out_rel)) {
                Ok(p) => p,
                Err(e) => return Ok(tool_err(format!("write sandbox error: {e}"))),
            };
            f.write_all(xml.as_bytes()).map_err(|e| e.to_string())?;
            written = Some(resolved.display().to_string());
        }
        return Ok(tool_ok(
            format!("Generated CDXML ({} bytes)", xml.len()),
            json!({
                "output_cdxml": written,
                "bytes": xml.len(),
            }),
        ));
    }
    if let Some(in_rel) = args.get("input_cdxml").and_then(Value::as_str) {
        let xml = match sandbox.read_to_string_limited(Path::new(in_rel), MAX_MCP_INPUT_FILE_BYTES)
        {
            Ok((_, s)) => s,
            Err(e) => return Ok(tool_err(format!("read sandbox error: {e}"))),
        };
        let mol = parse_cdxml(&xml).map_err(|e| e.to_string())?;
        let analysis = analyze_molecule(&mol);
        let smi = mol.canonical_smiles();
        return Ok(tool_ok(
            format!("Parsed CDXML -> {smi} ({})", analysis.formula_with_charge),
            json!({
                "canonical_smiles": smi,
                "analysis": analysis,
                "warnings": mol.warnings,
            }),
        ));
    }
    Ok(tool_err(
        "provide either 'smiles' (export) or 'input_cdxml' (import)".to_string(),
    ))
}

fn tool_rxn_stoichiometry(args: &Value, sandbox: &PathSandbox) -> Result<Value, String> {
    let rsmiles = args
        .get("reaction_smiles")
        .and_then(Value::as_str)
        .ok_or_else(|| "missing 'reaction_smiles'".to_string())?;
    let mmol = args.get("mmol").and_then(Value::as_f64).unwrap_or(1.0);
    let actual_mg = args.get("actual_mg").and_then(Value::as_f64);
    let rxn = parse_reaction_smiles(rsmiles).map_err(|e| e.to_string())?;
    let stoich = compute_stoichiometry(&rxn, mmol, actual_mg);
    let mut written = None;
    if let Some(out_rel) = args.get("output_rxn").and_then(Value::as_str) {
        let rxn_text = write_rxn_v2000(&rxn);
        let (resolved, mut f) = match sandbox.open_write(Path::new(out_rel)) {
            Ok(p) => p,
            Err(e) => return Ok(tool_err(format!("write sandbox error: {e}"))),
        };
        f.write_all(rxn_text.as_bytes())
            .map_err(|e| e.to_string())?;
        written = Some(resolved.display().to_string());
    }
    Ok(tool_ok(
        format!(
            "Reaction atom economy {:.1}%, theoretical yield {:.2} mg",
            stoich.atom_economy_pct, stoich.theoretical_yield_mg
        ),
        json!({
            "reaction_smiles": rxn.reaction_smiles(),
            "output_rxn": written,
            "stoichiometry": stoich,
        }),
    ))
}

fn tool_attach_superatom(args: &Value) -> Result<Value, String> {
    let smiles = args
        .get("smiles")
        .and_then(Value::as_str)
        .ok_or_else(|| "missing 'smiles'".to_string())?;
    let group = args
        .get("group")
        .and_then(Value::as_str)
        .ok_or_else(|| "missing 'group'".to_string())?;
    let atom_idx = args
        .get("atom_index")
        .and_then(Value::as_u64)
        .map(|v| v as usize);
    let mol = parse_smiles(smiles).map_err(|e| e.to_string())?;
    let updated = attach_superatom_group(&mol, atom_idx, group).map_err(|e| e.to_string())?;
    let analysis = analyze_molecule(&updated);
    let smi = updated.canonical_smiles();
    Ok(tool_ok(
        format!(
            "Attached {group} -> {smi} ({})",
            analysis.formula_with_charge
        ),
        json!({
            "canonical_smiles": smi,
            "analysis": analysis,
        }),
    ))
}

fn tool_plan_synthesis(args: &Value, sandbox: &PathSandbox) -> Result<Value, String> {
    let smiles = args
        .get("smiles")
        .and_then(Value::as_str)
        .ok_or_else(|| "missing 'smiles'".to_string())?;
    let name = args.get("name").and_then(Value::as_str);
    let max_steps = args
        .get("max_steps")
        .and_then(Value::as_u64)
        .unwrap_or(4)
        .clamp(1, 6) as usize;
    let top = args
        .get("top")
        .and_then(Value::as_u64)
        .unwrap_or(5)
        .clamp(1, 10) as usize;
    let scale_mg = args
        .get("scale_mg")
        .and_then(Value::as_f64)
        .unwrap_or(100.0)
        .max(1.0);
    let cfg = mol_synth::RetroConfig {
        max_steps,
        max_routes: top,
        target_scale_mg: scale_mg,
        ..Default::default()
    };
    let plan = mol_synth::plan_retrosynthesis_smiles(smiles, name, &cfg)
        .map_err(|e| format!("invalid SMILES: {e}"))?;
    let mut out_md_written = None;
    if let Some(out_rel) = args.get("out_md_path").and_then(Value::as_str) {
        let md = plan.render_markdown();
        let (resolved, mut f) = match sandbox.open_write(Path::new(out_rel)) {
            Ok(p) => p,
            Err(e) => return Ok(tool_err(format!("write sandbox error: {e}"))),
        };
        f.write_all(md.as_bytes()).map_err(|e| e.to_string())?;
        out_md_written = Some(resolved.display().to_string());
    }
    let summary = format!(
        "Planned {} synthetic route(s) for {} ({} nodes, {} ms)",
        plan.routes.len(),
        plan.target_name,
        plan.nodes_expanded,
        plan.elapsed_ms
    );
    let mut payload = serde_json::to_value(&plan).unwrap_or_else(|_| json!({}));
    if let Some(obj) = payload.as_object_mut()
        && let Some(p) = out_md_written
    {
        obj.insert("written_md_path".to_string(), json!(p));
    }
    Ok(tool_ok(summary, payload))
}
