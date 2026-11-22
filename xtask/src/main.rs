#![forbid(unsafe_code)]

//! `xtask` — ChemCraft Architectural Guardrails (`storytold/craft-anything` discipline).
//!
//! Run `cargo xtask layers` to verify strict acyclic `L0 -> L4` crate dependency layering,
//! `#![forbid(unsafe_code)]` across every crate entry point, and workspace panic-denial lints.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

#[derive(Debug, Clone, Copy)]
struct LayerSpec {
    name: &'static str,
    rel_dir: &'static str,
    layer: u8,
    role: &'static str,
}

const WORKSPACE_LAYERS: &[LayerSpec] = &[
    LayerSpec {
        name: "mol-graph",
        rel_dir: "crates/mol-graph",
        layer: 0,
        role: "L0 Foundation (SMILES/SDF/CDXML/RXN/PDB Graph IR & 3D EtKDG Coords)",
    },
    LayerSpec {
        name: "mol-vector",
        rel_dir: "crates/mol-vector",
        layer: 1,
        role: "L1 2D ACS-1996 Vector Depiction (SVG/EPS/PDF/PNG)",
    },
    LayerSpec {
        name: "mol-admet-bridge",
        rel_dir: "crates/mol-admet-bridge",
        layer: 1,
        role: "L1 Live AdmetCraft IPC & 1H/13C NMR Shift Predictor",
    },
    LayerSpec {
        name: "mol-sar",
        rel_dir: "crates/mol-sar",
        layer: 2,
        role: "L2 DataWarrior-style SAR Table, MMPs & Activity Cliffs",
    },
    LayerSpec {
        name: "mol-peptide",
        rel_dir: "crates/mol-peptide",
        layer: 2,
        role: "L2 Peptide, Macrocycle, HELM & NRPS Biosynthesis Studio",
    },
    LayerSpec {
        name: "mol-synth",
        rel_dir: "crates/mol-synth",
        layer: 3,
        role: "L3 Buy-Before-Make Retro* Retrosynthesis & MedChem Stock Catalog",
    },
    LayerSpec {
        name: "mol-engine",
        rel_dir: "crates/mol-engine",
        layer: 3,
        role: "L3 CoW CommandRegistry, .chemcraft ZIP, Screening Funnel & PathSandbox",
    },
    LayerSpec {
        name: "chemcraft-cli",
        rel_dir: "apps/chemcraft-cli",
        layer: 4,
        role: "L4 Headless CLI, Hardened MCP Server, Craft Link & ArtCraft Export",
    },
    LayerSpec {
        name: "chemcraft",
        rel_dir: "apps/chemcraft",
        layer: 4,
        role: "L4 Native Desktop 2D Chemistry & Peptide/NRP GUI",
    },
];

fn workspace_root() -> PathBuf {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    manifest_dir
        .parent()
        .map_or(manifest_dir.clone(), Path::to_path_buf)
}

fn parse_dependencies_section(toml_text: &str) -> (Vec<String>, bool) {
    let mut deps = Vec::new();
    let mut in_deps = false;
    let mut in_lints = false;
    let mut has_workspace_lints = false;

    for raw in toml_text.lines() {
        let line = raw.split('#').next().unwrap_or("").trim();
        if line.is_empty() {
            continue;
        }
        if line.starts_with('[') && line.ends_with(']') {
            let hdr = line.trim_start_matches('[').trim_end_matches(']').trim();
            in_deps = hdr == "dependencies" || hdr == "build-dependencies";
            in_lints = hdr == "lints";
            continue;
        }
        if in_lints && line.replace(' ', "") == "workspace=true" {
            has_workspace_lints = true;
        }
        if in_deps && let Some((lhs, _)) = line.split_once('=') {
            let key = lhs.trim();
            let dep_name = key.split('.').next().unwrap_or(key).trim();
            if !dep_name.is_empty() {
                deps.push(dep_name.to_string());
            }
        }
    }
    (deps, has_workspace_lints)
}

fn check_layers() -> Result<(), String> {
    let root = workspace_root();
    let root_toml = std::fs::read_to_string(root.join("Cargo.toml"))
        .map_err(|e| format!("failed to read root Cargo.toml: {e}"))?;

    let mut violations = Vec::new();
    for required in [
        "unsafe_code = \"forbid\"",
        "unwrap_used = \"deny\"",
        "expect_used = \"deny\"",
        "panic = \"deny\"",
    ] {
        if !root_toml.contains(required) {
            violations.push(format!("root Cargo.toml missing lint policy: `{required}`"));
        }
    }

    let layer_map: BTreeMap<&str, u8> =
        WORKSPACE_LAYERS.iter().map(|s| (s.name, s.layer)).collect();

    println!("# ChemCraft Architecture Layer Verification (`cargo xtask layers`)\n");
    println!(
        "| Layer | Crate | Role | Internal Deps | `#![forbid(unsafe_code)]` | Workspace Lints | Status |"
    );
    println!("|---|---|---|---|---:|---:|---|");

    for spec in WORKSPACE_LAYERS {
        let crate_dir = root.join(spec.rel_dir);
        let cargo_toml_path = crate_dir.join("Cargo.toml");
        let toml_text = std::fs::read_to_string(&cargo_toml_path)
            .map_err(|e| format!("failed to read {}: {e}", cargo_toml_path.display()))?;
        let (deps, has_ws_lints) = parse_dependencies_section(&toml_text);

        let mut internal_deps = Vec::new();
        let mut crate_ok = true;
        for dep in &deps {
            if let Some(&dep_layer) = layer_map.get(dep.as_str()) {
                internal_deps.push(format!("{dep}(L{dep_layer})"));
                if dep_layer >= spec.layer {
                    crate_ok = false;
                    violations.push(format!(
                        "LAYER VIOLATION: '{}' (L{}) depends on '{}' (L{})",
                        spec.name, spec.layer, dep, dep_layer
                    ));
                }
            }
        }

        if !has_ws_lints {
            crate_ok = false;
            violations.push(format!(
                "LINT VIOLATION: '{}' missing `[lints] workspace = true`",
                spec.name
            ));
        }

        let entry_path = if crate_dir.join("src/lib.rs").is_file() {
            crate_dir.join("src/lib.rs")
        } else {
            crate_dir.join("src/main.rs")
        };
        let entry_src = std::fs::read_to_string(&entry_path)
            .map_err(|e| format!("failed to read {}: {e}", entry_path.display()))?;
        let has_forbid_unsafe = entry_src
            .lines()
            .take(10)
            .any(|l| l.trim() == "#![forbid(unsafe_code)]");
        if !has_forbid_unsafe {
            crate_ok = false;
            violations.push(format!(
                "SAFETY VIOLATION: '{}' missing `#![forbid(unsafe_code)]`",
                spec.name
            ));
        }

        let deps_str = if internal_deps.is_empty() {
            "—".to_string()
        } else {
            internal_deps.join(", ")
        };
        println!(
            "| L{} | `{}` | {} | {} | {} | {} | {} |",
            spec.layer,
            spec.name,
            spec.role,
            deps_str,
            if has_forbid_unsafe { "PASS" } else { "FAIL" },
            if has_ws_lints { "PASS" } else { "FAIL" },
            if crate_ok { "PASS" } else { "FAIL" }
        );
    }

    println!();
    if violations.is_empty() {
        println!(
            "PASS: {} ChemCraft crates verified across L0..L4 (0 layer violations, 0 unsafe_code, 100% workspace lints).",
            WORKSPACE_LAYERS.len()
        );
        Ok(())
    } else {
        for v in &violations {
            eprintln!("ERROR: {v}");
        }
        Err(format!(
            "{} architectural guardrail violation(s) detected",
            violations.len()
        ))
    }
}

fn main() -> ExitCode {
    let arg = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "layers".to_string());
    let res = match arg.as_str() {
        "layers" | "check" => check_layers(),
        _ => Err(format!("unknown xtask subcommand '{arg}' (use `layers`)")),
    };
    match res {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("xtask failed: {e}");
            ExitCode::FAILURE
        }
    }
}
