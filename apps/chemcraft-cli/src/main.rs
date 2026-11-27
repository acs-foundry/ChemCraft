#![forbid(unsafe_code)]

mod mcp;

use std::path::{Path, PathBuf};

use clap::{Parser, Subcommand, ValueEnum};
use mol_admet_bridge::{AdmetCraftProvider, AdmetProvider, AdmetReportFormat, predict_nmr};
use mol_engine::{
    BuiltinLibraryPreset, DEFAULT_MAX_LIBRARY_RECORDS, LibraryFilterConfig, PathSandbox,
    ScreeningFunnel, decode_chemcraft_zip, default_chimerax_export_dir,
    default_dockcraft_handoff_dir, export_library_hits_sdf3d, filter_and_prepare_library,
    launch_chimerax_for_mol, launch_dockcraft_for_mols, load_builtin_library,
    merge_dockcraft_results, read_dockcraft_results, save_chemcraft_atomic, stream_load_library,
    write_chimerax_bundle, write_dockcraft_bundle,
};
use mol_graph::{
    MolGraph, analyze_molecule, attach_superatom_group, compute_stoichiometry, ensure_3d_coords,
    generate_3d_coords, parse_cdxml, parse_reaction_smiles, parse_rxn_v2000, parse_smiles,
    resolve_chemical_name, write_cdxml, write_pdb, write_rxn_v2000, write_sdf_v2000,
    write_sdf_v2000_3d,
};
use mol_sar::{
    SarFilter, SarTable, find_activity_cliffs, find_duplicate_compounds,
    find_matched_molecular_pairs, ingest_experimental_assay_text, render_dbtl_sar_csv,
    render_next_round_seeds_smi, render_qsar_recalibration_csv,
};
use mol_synth::{
    E3RecruiterPreset, RetroConfig, StockCatalog, default_stock_catalog, default_stock_tsv_path,
    enumerate_protac_library, enumerate_rgroup_library, plan_retrosynthesis,
};
use mol_vector::{
    DepictionPreset, SceneStyle, build_mol_scene, generate_2d_coords, scene_to_eps,
    scene_to_pdf_bytes, scene_to_png_bytes, scene_to_svg,
};

#[derive(Debug, Parser)]
#[command(
    name = "chemcraft-cli",
    version,
    about = "ChemCraft CLI & MCP server — Clean-room 2D vector chemistry, DataWarrior-style SAR, and live AdmetCraft integration"
)]
struct Cli {
    #[command(subcommand)]
    cmd: SubCmd,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum RenderFormat {
    Svg,
    Eps,
    Pdf,
    Png,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum EnrichFormat {
    Csv,
    Json,
    Pdf,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum Export3dFormat {
    Pdb,
    Sdf3d,
}

#[derive(Debug, Subcommand)]
enum SubCmd {
    /// Render a SMILES string to vector SVG, EPS, PDF, or raster PNG (with optional ACS Document 1996 preset).
    Render {
        smiles: String,
        #[arg(long)]
        out: PathBuf,
        #[arg(long, value_enum, default_value_t = RenderFormat::Svg)]
        format: RenderFormat,
        #[arg(long, default_value_t = 520.0)]
        width: f64,
        #[arg(long, default_value_t = 380.0)]
        height: f64,
        #[arg(long)]
        dark: bool,
        #[arg(long)]
        acs1996: bool,
        #[arg(long)]
        monochrome: bool,
    },
    /// Filter a chemical library (.sdf, .dwar, .csv) by SMARTS substructure, Tanimoto, and Lipinski/PAINS.
    Filter {
        input: PathBuf,
        #[arg(long)]
        substructure: Option<String>,
        #[arg(long)]
        reference: Option<String>,
        #[arg(long, default_value_t = 0.0)]
        min_tanimoto: f64,
        #[arg(long)]
        lipinski_only: bool,
        #[arg(long)]
        no_pains: bool,
        #[arg(long)]
        out: PathBuf,
    },
    /// Enrich an SDF/DWAR/CSV library with live AdmetCraft predictions and export CSV/JSON/PDF.
    Enrich {
        input: PathBuf,
        #[arg(long)]
        with_admet: bool,
        #[arg(long)]
        out: PathBuf,
        #[arg(long, value_enum, default_value_t = EnrichFormat::Csv)]
        format: EnrichFormat,
    },
    /// Compute Matched Molecular Pairs (MMP) and Activity Cliffs (SALI).
    Mmp {
        input: PathBuf,
        #[arg(long, default_value_t = 0.40)]
        min_tanimoto: f64,
        #[arg(long)]
        out: Option<PathBuf>,
    },
    /// Generate a 3D conformer from SMILES and export PDB or 3D SDF V2000.
    Export3d {
        smiles: String,
        #[arg(long)]
        out: PathBuf,
        #[arg(long, value_enum)]
        format: Option<Export3dFormat>,
    },
    /// Export a self-contained ChimeraX 3D bundle (.pdb/.sdf + .cxc) and optionally launch ChimeraX.
    Chimerax {
        #[arg(long)]
        smiles: Option<String>,
        #[arg(long)]
        out_dir: Option<PathBuf>,
        #[arg(long)]
        launch: bool,
    },
    /// Export a 3D ligand or chemical library (.sdf/.smi/.pdb) for DockCraft virtual screening and optionally launch DockCraft.
    Dockcraft {
        #[arg(long)]
        smiles: Option<String>,
        #[arg(long)]
        library: Option<PathBuf>,
        #[arg(long)]
        out_dir: Option<PathBuf>,
        #[arg(long)]
        launch: bool,
    },
    /// Fold a DockCraft `*.dockcraft.json` results sidecar back into a .chemcraft project's SAR
    /// table, so docking scores sit beside the measured ones and MMP / activity cliffs can be
    /// computed on either.
    DockResults {
        /// Project to merge into
        #[arg(long)]
        project: PathBuf,
        /// DockCraft results sidecar (*.dockcraft.json)
        #[arg(long)]
        results: PathBuf,
        /// Where to write the merged project. Defaults to overwriting --project.
        #[arg(long)]
        out: Option<PathBuf>,
    },
    /// Stream, pre-filter (Ro5/Lead-Like/Ro3/BBB/SMARTS/MaxMin), protonate (pH 7.4), and prepare a 3D screening library for DockCraft.
    Library {
        /// Built-in preset ('fda', 'fragments', 'kinase')
        #[arg(long)]
        preset: Option<String>,
        /// Local library file (.sdf, .sdf.gz, .smi, .smi.gz, .csv, .dwar)
        #[arg(long)]
        input: Option<PathBuf>,
        /// Funnel ('all', 'ro5', 'lead', 'ro3', 'bbb')
        #[arg(long, default_value = "ro5")]
        funnel: String,
        #[arg(long)]
        smarts: Option<String>,
        #[arg(long)]
        reference: Option<String>,
        #[arg(long, default_value_t = 0.0)]
        min_tanimoto: f64,
        /// Optional MaxMin ECFP4 diversity subset size
        #[arg(long)]
        diverse_n: Option<usize>,
        #[arg(long, default_value_t = 500)]
        max_hits: usize,
        #[arg(long, default_value_t = true)]
        ph74: bool,
        #[arg(long)]
        out: PathBuf,
        #[arg(long)]
        send_dockcraft: bool,
    },
    /// Inspect SMILES structure, key, coordinates, and AdmetCraft profile as JSON.
    Inspect { smiles: String },
    /// ChemDraw Analysis Window: exact monoisotopic mass, average MW, elemental analysis (%C/%H/...), DBE, and ESI+/ESI- adducts.
    Analyze { smiles: String },
    /// Predict 1H and 13C NMR chemical shifts, multiplicities, and integrations as JSON.
    Nmr { smiles: String },
    /// Convert a common drug/solvent/reagent name or systematic IUPAC name into SMILES & analysis JSON.
    Name { name: String },
    /// Attach a ChemDraw superatom/protecting group (Boc, Fmoc, Cbz, Ph, Bn, CF3, CO2Me, SO2Me, ...) to a SMILES.
    Group {
        smiles: String,
        #[arg(long)]
        attach: String,
        #[arg(long)]
        atom: Option<usize>,
    },
    /// Export SMILES to ChemDraw `.cdxml` or import `.cdxml` to canonical SMILES & analysis JSON.
    Cdxml {
        #[arg(long)]
        smiles: Option<String>,
        #[arg(long)]
        input: Option<PathBuf>,
        #[arg(long)]
        out: Option<PathBuf>,
    },
    /// Parse/export MDL `.rxn` V2000 reactions and calculate stoichiometry, atom economy, and yield.
    Rxn {
        #[arg(long)]
        rsmiles: Option<String>,
        #[arg(long)]
        input: Option<PathBuf>,
        #[arg(long)]
        out: Option<PathBuf>,
        #[arg(long, default_value_t = 1.0)]
        mmol: f64,
        #[arg(long)]
        actual_mg: Option<f64>,
    },
    /// Multi-step Buy-Before-Make + Retro* retrosynthesis planning, stoichiometry, ELN Markdown & MDL `.rxn` export.
    Retro {
        /// Target SMILES string (optional if `--init-stock` is used).
        smiles: Option<String>,
        /// Target display name.
        #[arg(long)]
        name: Option<String>,
        /// Maximum synthetic steps in longest linear sequence (`1..=6`).
        #[arg(long, default_value_t = 4)]
        max_steps: usize,
        /// Top ranked synthetic routes to return (`1..=10`).
        #[arg(long, default_value_t = 5)]
        top: usize,
        /// Target synthesis scale in mg for stoichiometry calculations.
        #[arg(long, default_value_t = 100.0)]
        scale_mg: f64,
        /// Optional extra building-block catalog file (`.tsv` or `.smi`) merged on top of default stock.
        #[arg(long)]
        stock: Option<PathBuf>,
        /// Write complete retrosynthesis plan as JSON.
        #[arg(long)]
        out_json: Option<PathBuf>,
        /// Write ELN-ready Markdown lab protocol.
        #[arg(long)]
        out_md: Option<PathBuf>,
        /// Export every step of Route #1 as MDL `.rxn` V2000 files inside this directory.
        #[arg(long)]
        out_rxn_dir: Option<PathBuf>,
        /// Initialize/export the default MedChem building-block TSV catalog in CRAFT_DATA_DIR/building_blocks.
        #[arg(long)]
        init_stock: bool,
    },
    /// Push a 3D PDB (or build a 3D peptide/NRP/SMILES) directly into a running BiomolCraft window over Craft Link v1 (`craft-ipc`).
    LinkBiomol {
        #[arg(long)]
        pdb: Option<PathBuf>,
        #[arg(long)]
        smiles: Option<String>,
        #[arg(long)]
        module: Option<String>,
    },
    /// Export a 2D ACS-1996 chemical/NRP publication bundle (PNG + WordCraft .docx + DeckCraft .pptx + GridCraft .xlsx) and optionally live-push over Craft Link v1.
    ExportArtcraft {
        #[arg(long)]
        smiles: Option<String>,
        #[arg(long)]
        name: Option<String>,
        #[arg(long, default_value = "wordcraft")]
        target: String,
        #[arg(long)]
        out_dir: PathBuf,
        #[arg(long)]
        push_live: bool,
        #[arg(long)]
        json: bool,
    },
    /// Enumerate a combinatorial R-group/Markush library or assemble a bifunctional PROTAC library (CRBN/VHL/IAP + alkyl/PEG/rigid linkers) into 3D SDF / SMILES / JSON.
    Enumerate {
        /// Mode: `rgroup` (Markush `[*]`/`[R1]`/`[R2]` scaffold) or `protac` (warhead + linkers + E3 recruiter).
        #[arg(long, default_value = "rgroup")]
        mode: String,
        /// Scaffold SMILES with `[*]` / `[R1]` / `[R2]` (for `--mode rgroup`).
        #[arg(long)]
        scaffold: Option<String>,
        /// POI warhead SMILES with optional `[*]` exit vector (for `--mode protac`).
        #[arg(long)]
        warhead: Option<String>,
        /// E3 ligase recruiter preset for PROTAC mode (`crbn` | `vhl` | `iap` | `all`).
        #[arg(long, default_value = "crbn")]
        e3: String,
        /// Optional custom R1 or linker file (`NAME\tSMILES` or one SMILES per line).
        #[arg(long)]
        r1_or_linkers: Option<PathBuf>,
        /// Optional custom R2 file (`NAME\tSMILES` or one SMILES per line).
        #[arg(long)]
        r2: Option<PathBuf>,
        /// Maximum enumerated compounds to keep.
        #[arg(long, default_value_t = 250)]
        max_compounds: usize,
        /// Output path (`.sdf` exports 3D conformers for DockCraft; `.smi` exports SMILES; `.json` exports full metadata).
        #[arg(long)]
        out: PathBuf,
        /// Optional JSON summary path.
        #[arg(long)]
        out_json: Option<PathBuf>,
    },
    /// DBTL Experimental Assay Ingestion (`Test -> Learn`): ingest lab IC50/EC50/Kd CSV/TSV, compute ADMET + LE/LLE + MMPs + Activity Cliffs, and emit Round-2 seeds & QSAR calibration CSV.
    AssayIngest {
        /// Experimental assay CSV or TSV (`compound_id,smiles,ic50_nm|pic50[,series]`).
        #[arg(long)]
        assay: PathBuf,
        /// Optional existing SAR CSV/SDF/DWAR table to merge with.
        #[arg(long)]
        existing_sar: Option<PathBuf>,
        /// Minimum ECFP4 Tanimoto for SALI activity cliff detection.
        #[arg(long, default_value_t = 0.40)]
        min_tanimoto: f64,
        /// Number of top high-potency/LLE seeds to export for Round 2.
        #[arg(long, default_value_t = 10)]
        top_seeds: usize,
        /// Output directory receiving `dbtl_report.json`, `sar_updated.csv`, `qsar_train_ready.csv`, and `next_round_seeds.smi`.
        #[arg(long)]
        out_dir: PathBuf,
    },
    /// Start the hardened MCP server (stdio or TCP with token + sandbox).
    Mcp {
        #[arg(long)]
        read_root: Option<PathBuf>,
        #[arg(long)]
        write_root: Option<PathBuf>,
        #[arg(long)]
        tcp: Option<String>,
        #[arg(long)]
        allow_remote: bool,
        #[arg(long)]
        mcp_token: Option<String>,
    },
}

fn load_table_from_path(path: &Path) -> Result<SarTable, Box<dyn std::error::Error>> {
    let text = std::fs::read_to_string(path)?;
    let ext = path
        .extension()
        .and_then(|s| s.to_str())
        .unwrap_or("")
        .to_ascii_lowercase();
    let mut table = SarTable::new();
    if ext == "dwar" {
        table.load_dwar_text(&text)?;
    } else if ext == "csv" {
        table.load_csv_text(&text)?;
    } else {
        table.load_sdf_text(&text)?;
    }
    Ok(table)
}

fn run() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();
    match cli.cmd {
        SubCmd::Render {
            smiles,
            out,
            format,
            width,
            height,
            dark,
            acs1996,
            monochrome,
        } => {
            let mut mol = parse_smiles(&smiles)?;
            generate_2d_coords(&mut mol);
            let mut style = if acs1996 {
                SceneStyle::from_preset(
                    if monochrome {
                        DepictionPreset::Acs1996Bw
                    } else {
                        DepictionPreset::Acs1996Color
                    },
                    width,
                    height,
                    dark,
                )
            } else {
                SceneStyle {
                    width,
                    height,
                    dark_mode: dark,
                    monochrome,
                    ..SceneStyle::default()
                }
            };
            style.monochrome = monochrome || (acs1996 && monochrome);
            let scene = build_mol_scene(&mol, &style, &[]);
            let bytes = match format {
                RenderFormat::Svg => scene_to_svg(&scene).into_bytes(),
                RenderFormat::Eps => scene_to_eps(&scene).into_bytes(),
                RenderFormat::Pdf => scene_to_pdf_bytes(&scene)?,
                RenderFormat::Png => scene_to_png_bytes(&scene, 2.0)?,
            };
            std::fs::write(&out, bytes)?;
            eprintln!("Rendered {} -> {}", mol.canonical_smiles(), out.display());
        }
        SubCmd::Filter {
            input,
            substructure,
            reference,
            min_tanimoto,
            lipinski_only,
            no_pains,
            out,
        } => {
            let mut table = load_table_from_path(&input)?;
            if lipinski_only || no_pains {
                table.enrich_with_admet(&AdmetCraftProvider::default());
            }
            let filter = SarFilter {
                smarts_substructure: substructure,
                reference_smiles: reference,
                min_tanimoto,
                lipinski_pass_only: lipinski_only,
                no_pains_only: no_pains,
                ..SarFilter::default()
            };
            let hit_indices = table.matching_indices(&filter)?;
            let refs: Vec<&mol_graph::MolGraph> = hit_indices
                .iter()
                .map(|&idx| &table.rows[idx].mol)
                .collect();
            let sdf_out = write_sdf_v2000(&refs);
            std::fs::write(&out, sdf_out)?;
            eprintln!(
                "Filtered {}/{} compounds -> {}",
                hit_indices.len(),
                table.rows.len(),
                out.display()
            );
        }
        SubCmd::Enrich {
            input,
            with_admet,
            out,
            format,
        } => {
            let mut table = load_table_from_path(&input)?;
            run_enrich_command(&mut table, with_admet, &out, format)?;
        }
        SubCmd::Mmp {
            input,
            min_tanimoto,
            out,
        } => {
            let mut table = load_table_from_path(&input)?;
            table.enrich_with_admet(&AdmetCraftProvider::default());
            let duplicates = find_duplicate_compounds(&table.rows);
            let cliffs = find_activity_cliffs(&table.rows, min_tanimoto);
            let mmps = find_matched_molecular_pairs(&table.rows);
            let payload = serde_json::to_string_pretty(&serde_json::json!({
                "compounds": table.rows.len(),
                "duplicates": duplicates,
                "activity_cliffs": cliffs,
                "matched_molecular_pairs": mmps,
            }))?;
            if let Some(p) = out {
                std::fs::write(p, payload)?;
            } else {
                println!("{payload}");
            }
        }
        SubCmd::Inspect { smiles } => {
            let mut mol = parse_smiles(&smiles)?;
            generate_2d_coords(&mut mol);
            let eval = AdmetCraftProvider::default().evaluate_smiles(&smiles, None)?;
            let analysis = analyze_molecule(&mol);
            println!(
                "{}",
                serde_json::to_string_pretty(&serde_json::json!({
                    "canonical_smiles": mol.canonical_smiles(),
                    "formula": mol.formula(),
                    "formula_with_charge": analysis.formula_with_charge,
                    "exact_mass": analysis.exact_mass,
                    "structural_key": mol.structural_key(),
                    "warnings": mol.warnings,
                    "analysis": analysis,
                    "admet": eval,
                }))?
            );
        }
        SubCmd::Analyze { smiles } => {
            let mol = parse_smiles(&smiles)?;
            let analysis = analyze_molecule(&mol);
            println!("{}", serde_json::to_string_pretty(&analysis)?);
        }
        SubCmd::Nmr { smiles } => {
            let mol = parse_smiles(&smiles)?;
            let spectrum = predict_nmr(&mol);
            println!("{}", serde_json::to_string_pretty(&spectrum)?);
        }
        SubCmd::Name { name } => {
            let mol = resolve_chemical_name(&name)?;
            let analysis = analyze_molecule(&mol);
            println!(
                "{}",
                serde_json::to_string_pretty(&serde_json::json!({
                    "name": mol.name,
                    "canonical_smiles": mol.canonical_smiles(),
                    "formula": analysis.formula_with_charge,
                    "exact_mass": analysis.exact_mass,
                    "molecular_weight": analysis.molecular_weight,
                }))?
            );
        }
        SubCmd::Group {
            smiles,
            attach,
            atom,
        } => {
            let mol = parse_smiles(&smiles)?;
            let updated = attach_superatom_group(&mol, atom, &attach)?;
            let analysis = analyze_molecule(&updated);
            println!(
                "{}",
                serde_json::to_string_pretty(&serde_json::json!({
                    "attached_group": attach,
                    "canonical_smiles": updated.canonical_smiles(),
                    "formula": analysis.formula_with_charge,
                    "exact_mass": analysis.exact_mass,
                    "molecular_weight": analysis.molecular_weight,
                }))?
            );
        }
        SubCmd::Cdxml { smiles, input, out } => {
            if let Some(smi) = smiles {
                let mut mol = parse_smiles(&smi)?;
                generate_2d_coords(&mut mol);
                let xml = write_cdxml(&mol);
                if let Some(out_path) = out {
                    std::fs::write(&out_path, xml.as_bytes())?;
                    eprintln!("Exported CDXML -> {}", out_path.display());
                } else {
                    print!("{xml}");
                }
            } else if let Some(in_path) = input {
                let xml = std::fs::read_to_string(&in_path)?;
                let mol = parse_cdxml(&xml)?;
                let analysis = analyze_molecule(&mol);
                println!(
                    "{}",
                    serde_json::to_string_pretty(&serde_json::json!({
                        "canonical_smiles": mol.canonical_smiles(),
                        "formula": analysis.formula_with_charge,
                        "exact_mass": analysis.exact_mass,
                        "atoms": mol.atoms.len(),
                        "bonds": mol.bonds.len(),
                        "warnings": mol.warnings,
                    }))?
                );
            } else {
                return Err(
                    "provide --smiles (to export CDXML) or --input <file.cdxml> (to parse CDXML)"
                        .into(),
                );
            }
        }
        SubCmd::Rxn {
            rsmiles,
            input,
            out,
            mmol,
            actual_mg,
        } => {
            let rxn = if let Some(rs) = rsmiles {
                parse_reaction_smiles(&rs)?
            } else if let Some(in_path) = input {
                let text = std::fs::read_to_string(&in_path)?;
                parse_rxn_v2000(&text)?
            } else {
                return Err("provide --rsmiles or --input <file.rxn>".into());
            };
            if let Some(out_path) = out {
                let rxn_text = write_rxn_v2000(&rxn);
                std::fs::write(&out_path, rxn_text.as_bytes())?;
                eprintln!("Exported RXN V2000 -> {}", out_path.display());
            }
            let stoich = compute_stoichiometry(&rxn, mmol, actual_mg);
            println!(
                "{}",
                serde_json::to_string_pretty(&serde_json::json!({
                    "name": rxn.name,
                    "reaction_smiles": rxn.reaction_smiles(),
                    "reactants": rxn.reactants.len(),
                    "agents": rxn.agents.len(),
                    "products": rxn.products.len(),
                    "stoichiometry": stoich,
                }))?
            );
        }
        SubCmd::Enumerate {
            mode,
            scaffold,
            warhead,
            e3,
            r1_or_linkers,
            r2,
            max_compounds,
            out,
            out_json,
        } => {
            let load_frags =
                |p: &Path| -> Result<Vec<(String, String)>, Box<dyn std::error::Error>> {
                    let txt = std::fs::read_to_string(p)?;
                    let mut v = Vec::new();
                    for (idx, line) in txt.lines().enumerate() {
                        let s = line.trim();
                        if s.is_empty() || s.starts_with('#') {
                            continue;
                        }
                        if let Some((a, b)) = s.split_once('\t') {
                            v.push((a.trim().to_string(), b.trim().to_string()));
                        } else if let Some((a, b)) = s.split_once(',') {
                            v.push((a.trim().to_string(), b.trim().to_string()));
                        } else {
                            v.push((format!("F{}", idx + 1), s.to_string()));
                        }
                    }
                    Ok(v)
                };
            let custom1 = match r1_or_linkers {
                Some(ref p) => Some(load_frags(p)?),
                None => None,
            };
            let custom2 = match r2 {
                Some(ref p) => Some(load_frags(p)?),
                None => None,
            };
            let mut lib = if mode.eq_ignore_ascii_case("protac") {
                let wh = warhead
                    .as_deref()
                    .ok_or("PROTAC mode requires --warhead <SMILES>")?;
                let preset = match e3.to_ascii_lowercase().as_str() {
                    "crbn" => E3RecruiterPreset::Crbn,
                    "vhl" => E3RecruiterPreset::Vhl,
                    "iap" => E3RecruiterPreset::Iap,
                    "all" => E3RecruiterPreset::All,
                    other => {
                        return Err(format!(
                            "unknown --e3 preset '{other}' (expected crbn|vhl|iap|all)"
                        )
                        .into());
                    }
                };
                enumerate_protac_library(wh, preset, custom1.as_deref(), max_compounds)
            } else {
                let scaf = scaffold
                    .as_deref()
                    .or(warhead.as_deref())
                    .ok_or("R-group mode requires --scaffold <SMILES>")?;
                enumerate_rgroup_library(
                    scaf,
                    custom1.as_deref(),
                    custom2.as_deref(),
                    max_compounds,
                )
            };

            let ext = out
                .extension()
                .and_then(|s| s.to_str())
                .unwrap_or("")
                .to_ascii_lowercase();
            if let Some(parent) = out.parent()
                && !parent.as_os_str().is_empty()
            {
                std::fs::create_dir_all(parent)?;
            }
            if ext == "smi" {
                let mut smi_txt = String::new();
                for c in &lib.compounds {
                    smi_txt.push_str(&format!("{}\t{}\n", c.canonical_smiles, c.id));
                }
                std::fs::write(&out, smi_txt)?;
            } else if ext == "json" {
                std::fs::write(&out, serde_json::to_string_pretty(&lib)?)?;
            } else {
                let mut mols = Vec::with_capacity(lib.compounds.len());
                for c in &mut lib.compounds {
                    if let Some(mut m) = c.mol.take() {
                        m.properties
                            .insert("MW".to_string(), format!("{:.2}", c.mw));
                        m.properties
                            .insert("LinkerOrR1".to_string(), c.linker_or_r1.clone());
                        m.properties
                            .insert("RecruiterOrR2".to_string(), c.recruiter_or_r2.clone());
                        m.properties
                            .insert("bRo5_OK".to_string(), c.bro5_ok.to_string());
                        ensure_3d_coords(&mut m);
                        mols.push(m);
                    }
                }
                let refs: Vec<&MolGraph> = mols.iter().collect();
                let sdf = write_sdf_v2000_3d(&refs);
                std::fs::write(&out, sdf)?;
            }
            if let Some(jpath) = out_json {
                if let Some(parent) = jpath.parent()
                    && !parent.as_os_str().is_empty()
                {
                    std::fs::create_dir_all(parent)?;
                }
                std::fs::write(&jpath, serde_json::to_string_pretty(&lib)?)?;
            }
            println!("{}", serde_json::to_string_pretty(&lib)?);
        }
        SubCmd::AssayIngest {
            assay,
            existing_sar,
            min_tanimoto,
            top_seeds,
            out_dir,
        } => {
            let assay_txt = std::fs::read_to_string(&assay)?;
            let base = match existing_sar {
                Some(ref p) => Some(load_table_from_path(p)?),
                None => None,
            };
            let (table, report) = ingest_experimental_assay_text(&assay_txt, base, min_tanimoto)?;
            std::fs::create_dir_all(&out_dir)?;
            let report_json = serde_json::to_string_pretty(&report)?;
            std::fs::write(out_dir.join("dbtl_report.json"), &report_json)?;
            std::fs::write(out_dir.join("sar_updated.csv"), render_dbtl_sar_csv(&table))?;
            std::fs::write(
                out_dir.join("qsar_train_ready.csv"),
                render_qsar_recalibration_csv(&table),
            )?;
            std::fs::write(
                out_dir.join("next_round_seeds.smi"),
                render_next_round_seeds_smi(&report, top_seeds),
            )?;
            println!("{report_json}");
        }
        SubCmd::Mcp {
            read_root,
            write_root,
            tcp,
            allow_remote,
            mcp_token,
        } => {
            let sandbox = PathSandbox::new(read_root.as_deref(), write_root.as_deref())?;
            if let Some(addr) = tcp {
                let tok = mcp_token.or_else(|| std::env::var("CHEMCRAFT_MCP_TOKEN").ok());
                mcp::run_mcp_tcp(&addr, &sandbox, allow_remote, tok)?;
            } else {
                mcp::run_mcp_stdio(&sandbox)?;
            }
        }
        SubCmd::Export3d {
            smiles,
            out,
            format,
        } => {
            let mut mol = parse_smiles(&smiles)?;
            generate_3d_coords(&mut mol);
            let resolved_fmt = format.unwrap_or_else(|| {
                let ext = out
                    .extension()
                    .and_then(|e| e.to_str())
                    .map(str::to_ascii_lowercase);
                match ext.as_deref() {
                    Some("sdf" | "mol") => Export3dFormat::Sdf3d,
                    _ => Export3dFormat::Pdb,
                }
            });
            let txt = match resolved_fmt {
                Export3dFormat::Pdb => write_pdb(&mol),
                Export3dFormat::Sdf3d => write_sdf_v2000_3d(&[&mol]),
            };
            std::fs::write(&out, txt.as_bytes())?;
            eprintln!("Exported 3D -> {}", out.display());
        }
        SubCmd::Chimerax {
            smiles,
            out_dir,
            launch,
        } => {
            let smi = smiles.ok_or_else(|| "provide --smiles".to_string())?;
            let mut m = parse_smiles(&smi)?;
            generate_3d_coords(&mut m);
            let mol = m;
            let target_dir = out_dir.unwrap_or_else(default_chimerax_export_dir);
            let bundle = if launch {
                launch_chimerax_for_mol(&mol, Some(&target_dir))?
            } else {
                write_chimerax_bundle(&mol, &target_dir, None)?
            };
            let info = serde_json::json!({
                "dir": bundle.dir.display().to_string(),
                "script": bundle.script_path.display().to_string(),
                "structure": bundle.structure_path.display().to_string(),
                "is_peptide": bundle.is_peptide,
                "launched": launch,
            });
            println!("{}", serde_json::to_string_pretty(&info)?);
        }
        SubCmd::DockResults {
            project,
            results,
            out,
        } => {
            let bytes = std::fs::read(&project)?;
            let mut doc = decode_chemcraft_zip(&bytes)?;
            let sidecar = read_dockcraft_results(&results)?;
            let pose_dir = results.parent().unwrap_or_else(|| Path::new("."));
            let before = doc.sar_table.rows.len();
            let outcome = merge_dockcraft_results(&mut doc.sar_table, &sidecar, pose_dir);
            let target = out.unwrap_or(project);
            save_chemcraft_atomic(&target, &doc)?;
            let info = serde_json::json!({
                "project": target.display().to_string(),
                "sar_rows": before,
                "merged": outcome.merged,
                "unmatched": outcome.unmatched,
                "receptor": sidecar.run.receptor_name,
                "scoring": sidecar.run.scoring,
            });
            println!("{}", serde_json::to_string_pretty(&info)?);
        }
        SubCmd::Dockcraft {
            smiles,
            library,
            out_dir,
            launch,
        } => {
            let mut mols: Vec<MolGraph> = if let Some(lib_path) = library {
                let table = load_table_from_path(&lib_path)?;
                table
                    .rows
                    .into_iter()
                    .map(|r| {
                        let mut m = r.mol;
                        ensure_3d_coords(&mut m);
                        m
                    })
                    .collect()
            } else if let Some(smi) = smiles {
                let mut m = parse_smiles(&smi)?;
                generate_3d_coords(&mut m);
                vec![m]
            } else {
                return Err("provide --smiles or --library".into());
            };
            for m in &mut mols {
                ensure_3d_coords(m);
            }
            let refs: Vec<&MolGraph> = mols.iter().collect();
            let target_dir = out_dir.unwrap_or_else(default_dockcraft_handoff_dir);
            let bundle = if launch {
                launch_dockcraft_for_mols(&refs, Some(&target_dir), None)?
            } else {
                write_dockcraft_bundle(&refs, &target_dir, None)?
            };
            let info = serde_json::json!({
                "dir": bundle.dir.display().to_string(),
                "sdf_3d": bundle.sdf_path.display().to_string(),
                "smi": bundle.smi_path.display().to_string(),
                "pdb": bundle.pdb_path.as_ref().map(|p| p.display().to_string()),
                "molecule_count": bundle.molecule_count,
                "launched": launch,
            });
            println!("{}", serde_json::to_string_pretty(&info)?);
        }
        SubCmd::Library {
            preset,
            input,
            funnel,
            smarts,
            reference,
            min_tanimoto,
            diverse_n,
            max_hits,
            ph74,
            out,
            send_dockcraft,
        } => {
            let (source_label, mols) = if let Some(p) = input {
                (
                    p.display().to_string(),
                    stream_load_library(&p, DEFAULT_MAX_LIBRARY_RECORDS)?,
                )
            } else if let Some(ref pre) = preset {
                let p = BuiltinLibraryPreset::from_str_lossy(pre).ok_or_else(|| {
                    format!("unknown library preset '{pre}' (expected fda|fragments|kinase)")
                })?;
                (p.label().to_string(), load_builtin_library(p))
            } else {
                return Err("provide --preset or --input".into());
            };
            let funnel_enum = ScreeningFunnel::from_str_lossy(&funnel).ok_or_else(|| {
                format!("unknown funnel '{funnel}' (expected all|ro5|lead|ro3|bbb)")
            })?;
            let cfg = LibraryFilterConfig {
                funnel: funnel_enum,
                smarts_substructure: smarts,
                reference_smiles: reference,
                min_tanimoto,
                diversity_pick_n: diverse_n,
                max_hits,
                protonate_ph74: ph74,
                generate_3d: true,
            };
            let res = filter_and_prepare_library(&source_label, mols, &cfg)?;
            let count = export_library_hits_sdf3d(&res, &out)?;
            if send_dockcraft {
                let refs: Vec<&MolGraph> = res.hits.iter().map(|h| &h.mol).collect();
                let _ = write_dockcraft_bundle(
                    &refs,
                    &default_dockcraft_handoff_dir(),
                    Some("chemcraft_screen_lib"),
                )?;
            }
            let summary = serde_json::json!({
                "source": res.source_label,
                "total_loaded": res.total_loaded,
                "passed_funnel": res.passed_funnel,
                "selected_count": count,
                "elapsed_ms": res.elapsed_ms,
                "out_sdf_3d": out.display().to_string(),
            });
            println!("{}", serde_json::to_string_pretty(&summary)?);
        }
        SubCmd::Retro {
            smiles,
            name,
            max_steps,
            top,
            scale_mg,
            stock,
            out_json,
            out_md,
            out_rxn_dir,
            init_stock,
        } => {
            if init_stock {
                let path = default_stock_tsv_path();
                let cat = StockCatalog::builtin_medchem();
                cat.write_tsv(&path)?;
                eprintln!(
                    "Initialized {} commercial MedChem building blocks -> {}",
                    cat.len(),
                    path.display()
                );
                if smiles.is_none() {
                    return Ok(());
                }
            }
            let smi = smiles
                .as_deref()
                .ok_or("provide target SMILES (or pass --init-stock)")?;
            let mut mol = parse_smiles(smi)?;
            if let Some(n) = name {
                mol.name = Some(n);
            }
            let cfg = RetroConfig {
                max_steps: max_steps.clamp(1, 6),
                max_routes: top.clamp(1, 10),
                target_scale_mg: scale_mg.max(1.0),
                ..Default::default()
            };
            let plan = if let Some(extra_stock) = stock {
                let mut cat = StockCatalog::load_default();
                let added = cat.merge_file(&extra_stock)?;
                eprintln!(
                    "Merged {added} custom building block(s) from {}",
                    extra_stock.display()
                );
                plan_retrosynthesis(&mol, &cat, &cfg)
            } else {
                plan_retrosynthesis(&mol, default_stock_catalog(), &cfg)
            };

            if let Some(jp) = out_json {
                std::fs::write(&jp, serde_json::to_string_pretty(&plan)?)?;
            }
            if let Some(mp) = out_md {
                std::fs::write(&mp, plan.render_markdown())?;
            }
            if let Some(rdir) = out_rxn_dir {
                std::fs::create_dir_all(&rdir)?;
                if let Some(best) = plan.routes.first() {
                    for st in &best.steps {
                        if let Some(rxn_txt) = st.to_rxn_v2000() {
                            let fname = format!(
                                "step_{:02}_{}.rxn",
                                st.step_number,
                                st.template_id.to_ascii_lowercase()
                            );
                            std::fs::write(rdir.join(fname), rxn_txt)?;
                        }
                    }
                }
            }
            println!("{}", plan.render_markdown());
        }
        SubCmd::LinkBiomol {
            pdb,
            smiles,
            module,
        } => {
            let target_path = if let Some(p) = pdb {
                p
            } else if let Some(smi) = smiles {
                let mut mol = parse_smiles(&smi)?;
                ensure_3d_coords(&mut mol);
                let out = std::env::temp_dir().join("chemcraft_link_mol.pdb");
                std::fs::write(&out, write_pdb(&mol))?;
                out
            } else {
                return Err("provide --pdb or --smiles".into());
            };
            let reply = craftbio_scene::push_pdb_to_live_biomolcraft(
                "chemcraft",
                &target_path.display().to_string(),
                module.as_deref(),
            )?;
            println!("Craft Link v1 → BiomolCraft: {reply}");
        }
        SubCmd::ExportArtcraft {
            smiles,
            name,
            target,
            out_dir,
            push_live,
            json,
        } => {
            let smi = smiles.ok_or_else(|| "provide --smiles <SMILES>".to_string())?;
            let mut m = parse_smiles(&smi)?;
            let t = name.unwrap_or_else(|| {
                m.name
                    .clone()
                    .unwrap_or_else(|| "ChemCraft Compound".to_string())
            });
            m.name = Some(t.clone());
            let canon = m.canonical_smiles();
            let an = analyze_molecule(&m);
            let md = format!(
                "## ChemCraft Structural Profile — {t}\n\n- **Canonical SMILES:** `{canon}`\n- **Formula:** `{}`\n- **MW:** `{:.2} Da` | **Exact Mass:** `{:.4} Da` | **DBE:** `{:.1}`\n- **Heavy Atoms:** `{}` | **Rings:** `{}` | **Chiral Centers:** `{}`\n- **Elemental Analysis:** `{}`\n",
                an.formula,
                an.molecular_weight,
                an.exact_mass,
                an.dbe,
                an.heavy_atoms,
                an.ring_count,
                an.chiral_centers,
                an.elemental_analysis_string
            );
            let r = vec![
                ("Compound", t.clone()),
                ("Canonical SMILES", canon),
                ("Formula", an.formula),
                ("MW (Da)", format!("{:.2}", an.molecular_weight)),
                ("Exact Mass (Da)", format!("{:.4}", an.exact_mass)),
                ("DBE", format!("{:.1}", an.dbe)),
                (
                    "Heavy Atoms / Rings",
                    format!("{} / {}", an.heavy_atoms, an.ring_count),
                ),
                ("Chiral Centers", an.chiral_centers.to_string()),
            ];
            let mut mol = m;
            generate_2d_coords(&mut mol);
            let style =
                SceneStyle::from_preset(DepictionPreset::Acs1996Color, 1200.0, 900.0, false);
            let scene = build_mol_scene(&mol, &style, &[]);
            let png_bytes = scene_to_png_bytes(&scene, 2.0)?;
            let mut rows = r;
            rows.push((
                "Generator",
                format!("ChemCraft v{}", env!("CARGO_PKG_VERSION")),
            ));
            let sum = craftbio_report::export_artcraft_bundle(&t, &png_bytes, &rows, &md, &out_dir)
                .map_err(|e| format!("ArtCraft export failed: {e}"))?;

            let mut push_res = None;
            if push_live {
                let tgt = craftbio_scene::ArtCraftTarget::parse(&target).ok_or_else(|| {
                    format!(
                        "unknown ArtCraft target '{target}' (use wordcraft|deckcraft|gridcraft)"
                    )
                })?;
                let md_path = match tgt {
                    craftbio_scene::ArtCraftTarget::DeckCraft => sum.deckcraft_slides_md.clone(),
                    _ => sum.wordcraft_md.clone(),
                };
                let native_doc = match tgt {
                    craftbio_scene::ArtCraftTarget::WordCraft => sum.wordcraft_docx.clone(),
                    craftbio_scene::ArtCraftTarget::DeckCraft => sum.deckcraft_pptx.clone(),
                    craftbio_scene::ArtCraftTarget::GridCraft => sum.gridcraft_xlsx.clone(),
                };
                let payload = craftbio_scene::ArtCraftPushPayload {
                    target: tgt,
                    title: t.clone(),
                    source_app: "chemcraft".to_string(),
                    bundle_dir: sum.out_dir.clone(),
                    figure_png_path: sum.figure_png.clone(),
                    markdown_path: md_path,
                    table_tsv_path: sum.gridcraft_tsv.clone(),
                    native_doc_path: native_doc,
                    summary_markdown: md.clone(),
                };
                push_res = Some(craftbio_scene::push_to_artcraft_via_craft_link(
                    None, &payload,
                )?);
            }

            if json {
                println!(
                    "{}",
                    serde_json::to_string_pretty(&serde_json::json!({
                        "bundle": sum,
                        "push": push_res,
                    }))?
                );
            } else {
                println!(
                    "Exported ChemCraft ArtCraft bundle -> {} (DOCX={:?}, PPTX={:?}, XLSX={:?})",
                    sum.out_dir, sum.wordcraft_docx, sum.deckcraft_pptx, sum.gridcraft_xlsx
                );
                if let Some(pr) = push_res {
                    println!("{}", pr.message);
                }
            }
        }
    }
    Ok(())
}

fn escape_csv_cell(raw: &str) -> String {
    let trimmed = raw.trim_start();
    let safe = if trimmed.starts_with(['=', '+', '-', '@']) {
        format!("'{raw}")
    } else {
        raw.to_string()
    };
    if safe.contains([',', '"', '\n', '\r']) {
        format!("\"{}\"", safe.replace('"', "\"\""))
    } else {
        safe
    }
}

fn run_enrich_command(
    table: &mut SarTable,
    with_admet: bool,
    out: &Path,
    format: EnrichFormat,
) -> Result<(), Box<dyn std::error::Error>> {
    if with_admet {
        let provider = AdmetCraftProvider::default();
        table.enrich_with_admet(&provider);
        let evals: Vec<_> = table
            .rows
            .iter()
            .filter_map(|r| r.admet_full.clone())
            .collect();
        let fmt = match format {
            EnrichFormat::Csv => AdmetReportFormat::Csv,
            EnrichFormat::Json => AdmetReportFormat::Json,
            EnrichFormat::Pdf => AdmetReportFormat::Pdf,
        };
        let mut f = std::fs::File::create(out)?;
        provider.export_report_to_writer(&evals, &mut f, fmt)?;
        eprintln!(
            "Wrote ADMET report ({} compounds) -> {}",
            evals.len(),
            out.display()
        );
        return Ok(());
    }

    match format {
        EnrichFormat::Csv => {
            use std::fmt::Write as _;
            let mut csv_out = String::from(
                "Molecule,Canonical SMILES,Structural Key,Formula,MW,HeavyAtoms,pIC50,Series\n",
            );
            for r in &table.rows {
                let pic50 = r.pic50.map_or_else(String::new, |v| format!("{v:.3}"));
                let _ = writeln!(
                    csv_out,
                    "{},{},{},{},{:.2},{},{},{}",
                    escape_csv_cell(&r.name),
                    escape_csv_cell(&r.canonical_smiles),
                    escape_csv_cell(&r.structural_key),
                    escape_csv_cell(&r.formula),
                    r.mw,
                    r.heavy_atoms,
                    pic50,
                    escape_csv_cell(&r.series)
                );
            }
            std::fs::write(out, csv_out)?;
        }
        EnrichFormat::Json => {
            let rows: Vec<_> = table
                .rows
                .iter()
                .map(|r| {
                    serde_json::json!({
                        "name": r.name,
                        "canonical_smiles": r.canonical_smiles,
                        "structural_key": r.structural_key,
                        "formula": r.formula,
                        "mw": r.mw,
                        "heavy_atoms": r.heavy_atoms,
                        "pic50": r.pic50,
                        "series": r.series,
                    })
                })
                .collect();
            std::fs::write(out, serde_json::to_string_pretty(&rows)?)?;
        }
        EnrichFormat::Pdf => {
            return Err("PDF ADMET report requires --with-admet".into());
        }
    }
    eprintln!(
        "Wrote table metadata ({} compounds) -> {}",
        table.rows.len(),
        out.display()
    );
    Ok(())
}

fn main() {
    if let Err(e) = run() {
        eprintln!("chemcraft-cli error: {e}");
        std::process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::BufReader;

    #[test]
    fn e2e_50_compounds_sdf_substructure_filter_tetrazole_bioisostere_svg_and_admet()
    -> Result<(), Box<dyn std::error::Error>> {
        let base = std::env::temp_dir().join("chemcraft_e2e_integration_test");
        let _ = std::fs::remove_dir_all(&base);
        let rdir = base.join("read");
        let wdir = base.join("write");
        std::fs::create_dir_all(&rdir)?;
        std::fs::create_dir_all(&wdir)?;

        // 1. Export realistic 50-compound library as SDF into read sandbox
        let demo = SarTable::demo_50_compounds();
        assert_eq!(demo.rows.len(), 50);
        let refs: Vec<&mol_graph::MolGraph> = demo.rows.iter().map(|r| &r.mol).collect();
        std::fs::write(rdir.join("library_50.sdf"), write_sdf_v2000(&refs))?;

        // 2. Drive full workflow via hardened MCP server: initialize -> search_substructure -> mutate_scaffold (tetrazole) -> export_vector_figure -> mol.analyze_admet
        let sb = PathSandbox::new(Some(&rdir), Some(&wdir))?;
        let requests = [
            serde_json::json!({"jsonrpc":"2.0","id":1,"method":"initialize","params":{}}).to_string(),
            serde_json::json!({
                "jsonrpc":"2.0","id":2,"method":"tools/call",
                "params":{
                    "name":"search_substructure",
                    "arguments":{"file":"library_50.sdf","smarts":"c1ccccc1","output_sdf":"aromatic_hits.sdf"}
                }
            }).to_string(),
            serde_json::json!({
                "jsonrpc":"2.0","id":3,"method":"tools/call",
                "params":{
                    "name":"mutate_scaffold",
                    "arguments":{"smiles":"CC(C)Cc1ccc(C(C)C(=O)O)cc1","group":"tetrazole","name":"Ibuprofen-Tetrazole"}
                }
            }).to_string(),
            serde_json::json!({
                "jsonrpc":"2.0","id":4,"method":"tools/call",
                "params":{
                    "name":"export_vector_figure",
                    "arguments":{"smiles":["CC(C)Cc1ccc(C(C)c2nnn[nH]2)cc1 Ibuprofen-Tetrazole"],"output_path":"bioisostere.svg","format":"svg"}
                }
            }).to_string(),
            serde_json::json!({
                "jsonrpc":"2.0","id":5,"method":"tools/call",
                "params":{
                    "name":"mol.analyze_admet",
                    "arguments":{"smiles":["CC(C)Cc1ccc(C(C)C(=O)O)cc1 Ibuprofen","CC(C)Cc1ccc(C(C)c2nnn[nH]2)cc1 Ibuprofen-Tetrazole"],"output_report":"admet_delta.csv","format":"csv"}
                }
            }).to_string(),
            // Strict JSON-RPC 2.0 check: missing jsonrpc or invalid params -> -32602 / -32600
            serde_json::json!({"jsonrpc":"2.0","id":6,"method":"tools/call","params":{"name":"draw_from_smiles","arguments":{}}}).to_string(),
        ].join("\n") + "\n";

        let mut out_buf = Vec::new();
        mcp::serve_stream_for_test(BufReader::new(requests.as_bytes()), &mut out_buf, &sb, None)?;
        let out_str = String::from_utf8(out_buf)?;
        let lines: Vec<serde_json::Value> = out_str
            .lines()
            .filter(|l| !l.trim().is_empty())
            .map(serde_json::from_str)
            .collect::<Result<_, _>>()?;
        assert_eq!(lines.len(), 6);

        // Aromatic hits from 50-compound SDF
        let total_loaded = lines[1]["result"]["structuredContent"]["total_loaded"]
            .as_u64()
            .ok_or("missing total_loaded")?;
        let matched = lines[1]["result"]["structuredContent"]["matched_count"]
            .as_u64()
            .ok_or("missing matched_count")?;
        assert!(
            matched >= 45,
            "expected >=45 aromatic hits, got {matched}/{total_loaded}; full resp: {}",
            lines[1]
        );
        assert!(wdir.join("aromatic_hits.sdf").exists());

        // Tetrazole bioisostere mutation changes logP, TPSA (+17.16 A^2), and preserves GI absorption
        let delta_tpsa = lines[2]["result"]["structuredContent"]["delta_tpsa"]
            .as_f64()
            .ok_or("missing delta_tpsa")?;
        assert!(
            delta_tpsa > 15.0,
            "expected tetrazole to raise TPSA, got {delta_tpsa}"
        );
        let new_gi = lines[2]["result"]["structuredContent"]["new_gi"]
            .as_str()
            .ok_or("missing new_gi")?;
        assert_eq!(new_gi, "High");

        // SVG valid + zero overlaps + AdmetCraft CSV generated
        let svg_text = std::fs::read_to_string(wdir.join("bioisostere.svg"))?;
        assert!(svg_text.contains("<svg") && svg_text.contains("</svg>"));
        assert_eq!(lines[3]["result"]["structuredContent"]["overlap_count"], 0);
        assert!(wdir.join("admet_delta.csv").exists());
        assert_eq!(lines[5]["error"]["code"], -32602);

        let _ = std::fs::remove_dir_all(&base);
        Ok(())
    }

    #[test]
    fn enrich_with_admet_differs_from_without_admet() -> Result<(), Box<dyn std::error::Error>> {
        let base = std::env::temp_dir().join("chemcraft_enrich_flag_test");
        let _ = std::fs::remove_dir_all(&base);
        std::fs::create_dir_all(&base)?;

        let mut table_plain = SarTable::new();
        table_plain.load_csv_text(
            "name,smiles,pIC50\nAspirin,CC(=O)Oc1ccccc1C(=O)O,5.2\nCaffeine,Cn1cnc2c1c(=O)n(C)c(=O)n2C,4.8\n",
        )?;
        let mut table_admet = table_plain.clone();

        let plain_csv = base.join("plain.csv");
        let admet_csv = base.join("admet.csv");
        run_enrich_command(&mut table_plain, false, &plain_csv, EnrichFormat::Csv)?;
        run_enrich_command(&mut table_admet, true, &admet_csv, EnrichFormat::Csv)?;

        let plain_text = std::fs::read_to_string(&plain_csv)?;
        let admet_text = std::fs::read_to_string(&admet_csv)?;
        assert_ne!(plain_text, admet_text);
        assert!(admet_text.contains("Consensus LogP"));
        assert!(!plain_text.contains("Consensus LogP"));

        let _ = std::fs::remove_dir_all(&base);
        Ok(())
    }

    #[test]
    fn mcp_export_sar_grid_eps_and_png_include_multiple_compounds()
    -> Result<(), Box<dyn std::error::Error>> {
        let base = std::env::temp_dir().join("chemcraft_mcp_grid_eps_png_test");
        let _ = std::fs::remove_dir_all(&base);
        let wdir = base.join("write");
        std::fs::create_dir_all(&wdir)?;

        let sb = PathSandbox::new(None, Some(&wdir))?;
        let smiles = [
            "CC(=O)Oc1ccccc1C(=O)O Aspirin",
            "Cn1cnc2c1c(=O)n(C)c(=O)n2C Caffeine",
            "CC(C)Cc1ccc(C(C)C(=O)O)cc1 Ibuprofen",
            "CC(=O)Nc1ccc(O)cc1 Paracetamol",
            "c1ccccc1 Benzene",
            "CCO Ethanol",
        ];
        let requests = [
            serde_json::json!({"jsonrpc":"2.0","id":0,"method":"initialize","params":{}})
                .to_string(),
            serde_json::json!({
                "jsonrpc":"2.0","id":1,"method":"tools/call",
                "params":{
                    "name":"export_vector_figure",
                    "arguments":{"smiles":smiles,"output_path":"grid.eps","format":"eps"}
                }
            })
            .to_string(),
            serde_json::json!({
                "jsonrpc":"2.0","id":2,"method":"tools/call",
                "params":{
                    "name":"export_vector_figure",
                    "arguments":{"smiles":smiles,"output_path":"grid.png","format":"png"}
                }
            })
            .to_string(),
        ]
        .join("\n")
            + "\n";

        let mut out_buf = Vec::new();
        mcp::serve_stream_for_test(BufReader::new(requests.as_bytes()), &mut out_buf, &sb, None)?;
        let eps = std::fs::read_to_string(wdir.join("grid.eps"))?;
        for label in [
            "Aspirin",
            "Caffeine",
            "Ibuprofen",
            "Paracetamol",
            "Benzene",
            "Ethanol",
        ] {
            assert!(
                eps.contains(label),
                "EPS grid missing compound label {label}"
            );
        }
        let png = std::fs::read(wdir.join("grid.png"))?;
        assert!(png.starts_with(b"\x89PNG\r\n\x1a\n"));
        assert!(png.len() > 2000);

        let _ = std::fs::remove_dir_all(&base);
        Ok(())
    }

    #[test]
    fn mcp_search_substructure_rejects_oversized_csv() -> Result<(), Box<dyn std::error::Error>> {
        let base = std::env::temp_dir().join("chemcraft_mcp_oversized_csv_test");
        let _ = std::fs::remove_dir_all(&base);
        let rdir = base.join("read");
        std::fs::create_dir_all(&rdir)?;

        let mut big = String::from("name,smiles\n");
        for i in 0..80_000 {
            big.push_str(&format!("Mol_{i},CC(=O)Oc1ccccc1C(=O)O\n"));
        }
        assert!(big.len() > 2 * 1024 * 1024);
        std::fs::write(rdir.join("huge.csv"), big)?;

        let sb = PathSandbox::new(Some(&rdir), None)?;
        let req = [
            serde_json::json!({"jsonrpc":"2.0","id":0,"method":"initialize","params":{}})
                .to_string(),
            serde_json::json!({
                "jsonrpc":"2.0","id":1,"method":"tools/call",
                "params":{"name":"search_substructure","arguments":{"file":"huge.csv","smarts":"c1ccccc1"}}
            })
            .to_string(),
        ]
        .join("\n")
            + "\n";
        let mut out_buf = Vec::new();
        mcp::serve_stream_for_test(BufReader::new(req.as_bytes()), &mut out_buf, &sb, None)?;
        let lines: Vec<serde_json::Value> = String::from_utf8(out_buf)?
            .lines()
            .filter(|l| !l.trim().is_empty())
            .map(serde_json::from_str)
            .collect::<Result<_, _>>()?;
        assert_eq!(lines.len(), 2);
        assert_eq!(lines[1]["result"]["isError"], true);

        let _ = std::fs::remove_dir_all(&base);
        Ok(())
    }

    #[test]
    fn mcp_draw_from_smiles_1000_atoms_avoids_quadratic_exponential_dos()
    -> Result<(), Box<dyn std::error::Error>> {
        let base = std::env::temp_dir().join("chemcraft_mcp_1000_atom_dos_test");
        let _ = std::fs::remove_dir_all(&base);
        std::fs::create_dir_all(&base)?;

        // 1,000-atom poly-chiral/diene branched cyclic SMILES (must finish in < 2s)
        let mut big_smi = String::from("N[C@@H](C(=O)O)");
        for _ in 0..80 {
            big_smi.push_str("/C=C/C[C@@H](c1ccccc1)");
        }
        let sb = PathSandbox::new(None, None)?;
        let req = [
            serde_json::json!({"jsonrpc":"2.0","id":0,"method":"initialize","params":{}})
                .to_string(),
            serde_json::json!({
                "jsonrpc":"2.0","id":1,"method":"tools/call",
                "params":{"name":"draw_from_smiles","arguments":{"smiles":big_smi}}
            })
            .to_string(),
        ]
        .join("\n")
            + "\n";

        let t0 = std::time::Instant::now();
        let mut out_buf = Vec::new();
        mcp::serve_stream_for_test(BufReader::new(req.as_bytes()), &mut out_buf, &sb, None)?;
        let elapsed = t0.elapsed();
        let lines: Vec<serde_json::Value> = String::from_utf8(out_buf)?
            .lines()
            .filter(|l| !l.trim().is_empty())
            .map(serde_json::from_str)
            .collect::<Result<_, _>>()?;
        assert_eq!(lines.len(), 2);
        assert_eq!(lines[1]["result"]["isError"], false);
        assert!(
            elapsed.as_secs_f64() < 2.5,
            "1000-atom MCP draw_from_smiles took too long: {elapsed:?}"
        );

        let _ = std::fs::remove_dir_all(&base);
        Ok(())
    }
}
