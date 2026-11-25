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
    MolGraph, SecondaryStructure, analyze_molecule, attach_superatom_group, build_peptide,
    compute_stoichiometry, ensure_3d_coords, generate_3d_coords, parse_cdxml,
    parse_reaction_smiles, parse_rxn_v2000, parse_smiles, resolve_chemical_name, write_cdxml,
    write_pdb, write_rxn_v2000, write_sdf_v2000, write_sdf_v2000_3d,
};
use mol_peptide::{
    NrpsBgc, PeptideConformation, ScanMode, analyze_peptide_developability, build_peptide_or_nrp,
    compile_nrps_bgc, default_monomer_catalog, design_nrps_from_peptide, landmark_nrp_presets,
    parse_antismash_gbk, parse_peptide_spec, plan_spps_vs_nrps, predict_a_domain_specificity,
    preset_nrps_bgc, propose_nrps_swap, render_developability_markdown, render_nrps_markdown,
    render_scan_csv, render_spps_markdown, run_peptide_scan,
};
use mol_sar::{
    SarFilter, SarTable, find_activity_cliffs, find_duplicate_compounds,
    find_matched_molecular_pairs,
};
use mol_synth::{
    RetroConfig, StockCatalog, default_stock_catalog, default_stock_tsv_path, plan_retrosynthesis,
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
enum PeptideSsCli {
    Beta,
    Alpha,
    Extended,
}

impl From<PeptideSsCli> for SecondaryStructure {
    fn from(value: PeptideSsCli) -> Self {
        match value {
            PeptideSsCli::Beta => SecondaryStructure::BetaStrand,
            PeptideSsCli::Alpha => SecondaryStructure::AlphaHelix,
            PeptideSsCli::Extended => SecondaryStructure::Extended,
        }
    }
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum Export3dFormat {
    Pdb,
    Sdf3d,
}

#[derive(Debug, Subcommand)]
enum NrpCmd {
    /// Build any linear, cyclic, depsipeptide, lariat, disulfide, or stapled Peptide / NRP (notation, preset ID, or HELM) with bRo5 developability, proteolysis, ESI-MS, and SPPS plan.
    Build {
        /// Extended sequence (`cyclo(Arg-Gly-Asp-D-Phe-NMe-Val)`), preset ID (`cyclosporin_a`, `daptomycin`, `polymyxin_b1`, `gramicidin_s`, `cilengitide`, `octreotide`, `acv_tripeptide`, `semaglutide_core`), or HELM string.
        spec: String,
        /// 3D conformation (`ring` | `alpha` | `beta` | `extended` | `ppii` | `hairpin`).
        #[arg(long, default_value = "beta")]
        ss: String,
        /// Optional SVG output path (2D depiction).
        #[arg(long)]
        svg: Option<PathBuf>,
        /// Optional PDB output path (with SSBOND + CONECT records).
        #[arg(long)]
        pdb: Option<PathBuf>,
        /// Optional 3D SDF V2000/V3000 output path.
        #[arg(long)]
        sdf3d: Option<PathBuf>,
        /// Optional Pistoia HELM output path.
        #[arg(long)]
        helm_out: Option<PathBuf>,
        /// Optional Markdown developability + SPPS dossier output path.
        #[arg(long)]
        report_md: Option<PathBuf>,
        /// Optional JSON report output path.
        #[arg(long)]
        out_json: Option<PathBuf>,
        /// Export 3D DockCraft handoff bundle.
        #[arg(long)]
        send_dockcraft: bool,
        /// Export ChimeraX `.pdb`/`.sdf` + `.cxc` bundle.
        #[arg(long)]
        send_chimerax: bool,
        /// Print raw JSON summary to stdout instead of Markdown.
        #[arg(long)]
        json: bool,
    },
    /// Analyze, compile, reverse-design, or engineer (domain-swap) an NRPS Biosynthetic Gene Cluster (preset, JSON, antiSMASH `.gbk`, or Stachelhaus 10-aa code).
    Bgc {
        /// Curated BGC preset (`gramicidin_s`, `acv`, `polymyxin_b1`, `daptomycin`).
        #[arg(long)]
        preset: Option<String>,
        /// Input NRPS BGC JSON file.
        #[arg(long)]
        input_json: Option<PathBuf>,
        /// Input antiSMASH GenBank `.gbk` file.
        #[arg(long)]
        antismash_gbk: Option<PathBuf>,
        /// Reverse-design a minimal NRPS BGC blueprint from a peptide/NRP notation or preset.
        #[arg(long)]
        from_seq: Option<String>,
        /// Predict probabilistic substrate candidates from a 10-aa Stachelhaus A-domain code.
        #[arg(long)]
        stachelhaus: Option<String>,
        /// Evaluate a module A-domain swap (`MOD:MONOMER`, e.g. `3:Ile`).
        #[arg(long)]
        swap: Option<String>,
        /// Write NRPS compilation/audit JSON.
        #[arg(long)]
        out_json: Option<PathBuf>,
        /// Write NRPS report Markdown.
        #[arg(long)]
        out_md: Option<PathBuf>,
        /// Write compiled 3D product PDB.
        #[arg(long)]
        pdb: Option<PathBuf>,
        /// Write compiled 3D product SDF.
        #[arg(long)]
        sdf3d: Option<PathBuf>,
        /// Write compiled 2D SVG depiction.
        #[arg(long)]
        svg: Option<PathBuf>,
        /// Print raw JSON on stdout instead of Markdown.
        #[arg(long)]
        json: bool,
    },
    /// Run a combinatorial peptide optimization scan (`ala` | `d-scan` | `nme` | `staple`).
    Scan {
        /// Parent peptide/NRP notation, preset ID, or HELM string.
        spec: String,
        /// Scan mode (`ala` | `d-scan` | `nme` | `staple`).
        #[arg(long, default_value = "nme")]
        mode: String,
        /// Write scan CSV summary table.
        #[arg(long)]
        out_csv: Option<PathBuf>,
        /// Write 3D multi-molecule SDF of all scanned analogs.
        #[arg(long)]
        out_sdf: Option<PathBuf>,
        /// Write JSON summary.
        #[arg(long)]
        out_json: Option<PathBuf>,
        /// Export DockCraft 3D handoff bundle for all variants.
        #[arg(long)]
        send_dockcraft: bool,
        /// Print raw JSON on stdout instead of compact table.
        #[arg(long)]
        json: bool,
    },
    /// List all curated landmark NRP presets and monomer catalog entries.
    Presets {
        /// Also list all built-in monomer codes and N-terminal caps.
        #[arg(long)]
        monomers: bool,
        /// Output as JSON.
        #[arg(long)]
        json: bool,
    },
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
    /// Build a peptide in 2D + 3D from 1-letter or 3-letter sequence (`Ac-Ala-Gly-Ser-NH2`, `AGYS`) and export SVG/PDB/3D-SDF.
    Peptide {
        sequence: String,
        #[arg(long, value_enum, default_value_t = PeptideSsCli::Beta)]
        ss: PeptideSsCli,
        #[arg(long)]
        svg: Option<PathBuf>,
        #[arg(long)]
        pdb: Option<PathBuf>,
        #[arg(long)]
        sdf3d: Option<PathBuf>,
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
        peptide: Option<String>,
        #[arg(long, value_enum, default_value_t = PeptideSsCli::Alpha)]
        ss: PeptideSsCli,
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
        peptide: Option<String>,
        #[arg(long)]
        library: Option<PathBuf>,
        #[arg(long, value_enum, default_value_t = PeptideSsCli::Beta)]
        ss: PeptideSsCli,
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
    /// Peptide & Non-Ribosomal Peptide (NRP / NRPS) Design Studio CLI (`build`, `bgc`, `scan`, `presets`).
    Nrp {
        #[command(subcommand)]
        action: NrpCmd,
    },
    /// Push a 3D PDB (or build a 3D peptide/NRP/SMILES) directly into a running BiomolCraft window over Craft Link v1 (`craft-ipc`).
    LinkBiomol {
        #[arg(long)]
        pdb: Option<PathBuf>,
        #[arg(long)]
        smiles: Option<String>,
        #[arg(long)]
        nrp: Option<String>,
        #[arg(long)]
        module: Option<String>,
    },
    /// Export a 2D ACS-1996 chemical/NRP publication bundle (PNG + WordCraft .docx + DeckCraft .pptx + GridCraft .xlsx) and optionally live-push over Craft Link v1.
    ExportArtcraft {
        #[arg(long)]
        smiles: Option<String>,
        #[arg(long)]
        nrp: Option<String>,
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
        SubCmd::Peptide {
            sequence,
            ss,
            svg,
            pdb,
            sdf3d,
        } => {
            let mol = build_peptide(&sequence, ss.into())?;
            if let Some(svg_out) = svg {
                let scene = build_mol_scene(&mol, &SceneStyle::default(), &[]);
                std::fs::write(&svg_out, scene_to_svg(&scene).as_bytes())?;
            }
            if let Some(pdb_out) = pdb {
                std::fs::write(&pdb_out, write_pdb(&mol).as_bytes())?;
            }
            if let Some(sdf_out) = sdf3d {
                std::fs::write(&sdf_out, write_sdf_v2000_3d(&[&mol]).as_bytes())?;
            }
            let summary = serde_json::json!({
                "name": mol.name,
                "sequence": mol.properties.get("PEPTIDE_SEQUENCE"),
                "one_letter": mol.properties.get("PEPTIDE_1LETTER"),
                "secondary_structure": mol.properties.get("SECONDARY_STRUCTURE"),
                "residues": mol.properties.get("PEPTIDE_RESIDUES"),
                "est_charge_ph74": mol.properties.get("EST_CHARGE_PH74"),
                "formula": mol.formula(),
                "mw": mol.molecular_weight(),
                "canonical_smiles": mol.canonical_smiles(),
            });
            println!("{}", serde_json::to_string_pretty(&summary)?);
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
            peptide,
            ss,
            out_dir,
            launch,
        } => {
            let mol = if let Some(seq) = peptide {
                build_peptide(&seq, ss.into())?
            } else if let Some(smi) = smiles {
                let mut m = parse_smiles(&smi)?;
                generate_3d_coords(&mut m);
                m
            } else {
                return Err("provide --smiles or --peptide".into());
            };
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
            peptide,
            library,
            ss,
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
            } else if let Some(seq) = peptide {
                vec![build_peptide(&seq, ss.into())?]
            } else if let Some(smi) = smiles {
                let mut m = parse_smiles(&smi)?;
                generate_3d_coords(&mut m);
                vec![m]
            } else {
                return Err("provide --smiles, --peptide, or --library".into());
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
        SubCmd::Nrp { action } => {
            run_nrp_command(action)?;
        }
        SubCmd::LinkBiomol {
            pdb,
            smiles,
            nrp,
            module,
        } => {
            let target_path = if let Some(p) = pdb {
                p
            } else if let Some(spec) = nrp {
                let cat = default_monomer_catalog();
                let conf = PeptideConformation::from_str_lossy("ring")
                    .unwrap_or(PeptideConformation::BetaStrand);
                let built = build_peptide_or_nrp(&spec, conf, &cat)?;
                let out = std::env::temp_dir().join("chemcraft_link_nrp.pdb");
                std::fs::write(&out, write_pdb(&built.mol))?;
                out
            } else if let Some(smi) = smiles {
                let mut mol = parse_smiles(&smi)?;
                ensure_3d_coords(&mut mol);
                let out = std::env::temp_dir().join("chemcraft_link_mol.pdb");
                std::fs::write(&out, write_pdb(&mol))?;
                out
            } else {
                return Err("provide --pdb, --nrp, or --smiles".into());
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
            nrp,
            name,
            target,
            out_dir,
            push_live,
            json,
        } => {
            let (mut mol, title, summary_md, mut rows) = if let Some(spec) = nrp {
                let cat = default_monomer_catalog();
                let conf = PeptideConformation::from_str_lossy("ring")
                    .unwrap_or(PeptideConformation::BetaStrand);
                let built = build_peptide_or_nrp(&spec, conf, &cat)?;
                let dev = analyze_peptide_developability(&built, &cat);
                let spps = plan_spps_vs_nrps(&built.spec, &cat);
                let mut md = render_developability_markdown(&dev);
                md.push_str("\n---\n\n");
                md.push_str(&render_spps_markdown(&spps));
                let t = name.unwrap_or_else(|| dev.name.clone());
                let r = vec![
                    ("Compound", t.clone()),
                    ("Formula", dev.formula.clone()),
                    (
                        "Monoisotopic Mass (Da)",
                        format!("{:.4}", dev.monoisotopic_mass),
                    ),
                    ("Average MW (Da)", format!("{:.2}", dev.average_mw)),
                    (
                        "WLOGP / TPSA",
                        format!("{:.2} / {:.1} Å²", dev.wlogp, dev.tpsa),
                    ),
                    ("bRo5 Oral Score", format!("{:.3}", dev.bro5_oral_score)),
                    ("Oral Verdict", dev.bro5_verdict.clone()),
                    ("HELM", built.helm.clone()),
                ];
                (built.mol, t, md, r)
            } else if let Some(smi) = smiles {
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
                (m, t, md, r)
            } else {
                return Err("provide --smiles <SMILES> or --nrp <SPEC|PRESET>".into());
            };
            generate_2d_coords(&mut mol);
            let style =
                SceneStyle::from_preset(DepictionPreset::Acs1996Color, 1200.0, 900.0, false);
            let scene = build_mol_scene(&mol, &style, &[]);
            let png_bytes = scene_to_png_bytes(&scene, 2.0)?;
            rows.push((
                "Generator",
                format!("ChemCraft v{}", env!("CARGO_PKG_VERSION")),
            ));
            let sum = craftbio_report::export_artcraft_bundle(
                &title,
                &png_bytes,
                &rows,
                &summary_md,
                &out_dir,
            )
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
                    title: title.clone(),
                    source_app: "chemcraft".to_string(),
                    bundle_dir: sum.out_dir.clone(),
                    figure_png_path: sum.figure_png.clone(),
                    markdown_path: md_path,
                    table_tsv_path: sum.gridcraft_tsv.clone(),
                    native_doc_path: native_doc,
                    summary_markdown: summary_md,
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

fn run_nrp_command(action: NrpCmd) -> Result<(), Box<dyn std::error::Error>> {
    let cat = default_monomer_catalog();
    match action {
        NrpCmd::Presets { monomers, json } => {
            let presets = landmark_nrp_presets();
            if json {
                let mut payload = serde_json::json!({
                    "presets": presets,
                    "base_monomer_count": cat.monomers.len(),
                    "stereomonomer_count": cat.total_stereomonomer_count(),
                    "n_cap_count": cat.n_caps.len(),
                });
                if monomers && let Some(obj) = payload.as_object_mut() {
                    obj.insert("monomers".to_string(), serde_json::to_value(&cat.monomers)?);
                    obj.insert("n_caps".to_string(), serde_json::to_value(&cat.n_caps)?);
                }
                println!("{}", serde_json::to_string_pretty(&payload)?);
            } else {
                println!(
                    "ChemCraft v{} — Peptide & NRP Presets ({} base monomers / {} L+D stereomonomers, {} N-caps):",
                    env!("CARGO_PKG_VERSION"),
                    cat.monomers.len(),
                    cat.total_stereomonomer_count(),
                    cat.n_caps.len()
                );
                for p in &presets {
                    println!(
                        "  - {:<18} | {:<15} | {:>9.4} Da | {}",
                        p.id, p.expected_formula, p.expected_monoisotopic_mass, p.name
                    );
                }
                if monomers {
                    println!("\nBuilt-in Monomers:");
                    for m in &cat.monomers {
                        println!(
                            "  - {:<10} [{}] {} (chiral={}, anchor={:?})",
                            m.code,
                            m.class.label(),
                            m.name,
                            m.has_alpha_chiral,
                            m.anchor_kind
                        );
                    }
                    println!("\nN-Terminal Caps & Fatty-Acyl Starters:");
                    for c in &cat.n_caps {
                        println!("  - {:<16} {}", c.code, c.name);
                    }
                }
            }
        }
        NrpCmd::Build {
            spec,
            ss,
            svg,
            pdb,
            sdf3d,
            helm_out,
            report_md,
            out_json,
            send_dockcraft,
            send_chimerax,
            json,
        } => {
            let conf = PeptideConformation::from_str_lossy(&ss)
                .ok_or_else(|| format!("unknown conformation '{ss}'"))?;
            let built = build_peptide_or_nrp(&spec, conf, &cat)?;
            let dev = analyze_peptide_developability(&built, &cat);
            let spps = plan_spps_vs_nrps(&built.spec, &cat);

            if let Some(p) = svg {
                let scene = build_mol_scene(&built.mol, &SceneStyle::default(), &[]);
                std::fs::write(&p, scene_to_svg(&scene))?;
            }
            if let Some(p) = pdb {
                std::fs::write(&p, write_pdb(&built.mol))?;
            }
            if let Some(p) = sdf3d {
                std::fs::write(&p, write_sdf_v2000_3d(&[&built.mol]))?;
            }
            if let Some(p) = helm_out {
                std::fs::write(&p, format!("{}\n", built.helm))?;
            }
            let mut full_md = render_developability_markdown(&dev);
            full_md.push_str("\n---\n\n");
            full_md.push_str(&render_spps_markdown(&spps));
            if let Some(p) = report_md {
                std::fs::write(&p, &full_md)?;
            }
            if send_dockcraft {
                let dir = default_dockcraft_handoff_dir();
                let bundle = write_dockcraft_bundle(&[&built.mol], &dir, Some("nrp_target"))?;
                eprintln!("DockCraft handoff bundle -> {}", bundle.sdf_path.display());
            }
            if send_chimerax {
                let dir = default_chimerax_export_dir();
                let bundle = write_chimerax_bundle(&built.mol, &dir, Some("nrp_target"))?;
                eprintln!("ChimeraX 3D bundle -> {}", bundle.script_path.display());
            }

            let payload = serde_json::json!({
                "developability": dev,
                "spps_plan": spps,
            });
            if let Some(p) = out_json {
                std::fs::write(&p, serde_json::to_string_pretty(&payload)?)?;
            }
            if json {
                println!("{}", serde_json::to_string_pretty(&payload)?);
            } else {
                println!("{full_md}");
            }
        }
        NrpCmd::Bgc {
            preset,
            input_json,
            antismash_gbk,
            from_seq,
            stachelhaus,
            swap,
            out_json,
            out_md,
            pdb,
            sdf3d,
            svg,
            json,
        } => {
            if let Some(code) = stachelhaus {
                let pred = predict_a_domain_specificity(&code);
                if let Some(p) = out_json {
                    std::fs::write(&p, serde_json::to_string_pretty(&pred)?)?;
                }
                println!("{}", serde_json::to_string_pretty(&pred)?);
                return Ok(());
            }

            let bgc: NrpsBgc = if let Some(pid) = preset {
                preset_nrps_bgc(&pid).ok_or_else(|| {
                    format!("unknown NRPS BGC preset '{pid}' (expected gramicidin_s, acv, polymyxin_b1, daptomycin)")
                })?
            } else if let Some(jp) = input_json {
                let txt = std::fs::read_to_string(&jp)?;
                serde_json::from_str(&txt)?
            } else if let Some(gp) = antismash_gbk {
                let txt = std::fs::read_to_string(&gp)?;
                parse_antismash_gbk(&txt)?
            } else if let Some(seq) = from_seq {
                let spec = parse_peptide_spec(&seq, &cat, PeptideConformation::BetaStrand)?;
                design_nrps_from_peptide(&spec)
            } else {
                return Err(
                    "provide one of --preset, --input-json, --antismash-gbk, --from-seq, or --stachelhaus"
                        .into(),
                );
            };

            if let Some(sw) = swap {
                let (mod_s, mon_s) = sw
                    .split_once(':')
                    .ok_or("expected --swap MOD:MONOMER (e.g. 3:Ile)")?;
                let mod_idx: usize = mod_s.trim().parse()?;
                let (ev, built) = propose_nrps_swap(&bgc, mod_idx, mon_s.trim(), &cat)?;
                if let Some(p) = pdb {
                    std::fs::write(&p, write_pdb(&built.mol))?;
                }
                if let Some(p) = sdf3d {
                    std::fs::write(&p, write_sdf_v2000_3d(&[&built.mol]))?;
                }
                if let Some(p) = svg {
                    let scene = build_mol_scene(&built.mol, &SceneStyle::default(), &[]);
                    std::fs::write(&p, scene_to_svg(&scene))?;
                }
                if let Some(p) = out_json {
                    std::fs::write(&p, serde_json::to_string_pretty(&ev)?)?;
                }
                println!("{}", serde_json::to_string_pretty(&ev)?);
                return Ok(());
            }

            let (rep, built) = compile_nrps_bgc(&bgc, &cat)?;
            let md = render_nrps_markdown(&rep);
            if let Some(p) = out_json {
                std::fs::write(&p, serde_json::to_string_pretty(&rep)?)?;
            }
            if let Some(p) = out_md {
                std::fs::write(&p, &md)?;
            }
            if let Some(p) = pdb {
                std::fs::write(&p, write_pdb(&built.mol))?;
            }
            if let Some(p) = sdf3d {
                std::fs::write(&p, write_sdf_v2000_3d(&[&built.mol]))?;
            }
            if let Some(p) = svg {
                let scene = build_mol_scene(&built.mol, &SceneStyle::default(), &[]);
                std::fs::write(&p, scene_to_svg(&scene))?;
            }
            if json {
                println!("{}", serde_json::to_string_pretty(&rep)?);
            } else {
                println!("{md}");
            }
        }
        NrpCmd::Scan {
            spec,
            mode,
            out_csv,
            out_sdf,
            out_json,
            send_dockcraft,
            json,
        } => {
            let scan_mode = ScanMode::from_str_lossy(&mode).ok_or_else(|| {
                format!("unknown scan mode '{mode}' (expected ala|d-scan|nme|staple)")
            })?;
            let parsed = parse_peptide_spec(&spec, &cat, PeptideConformation::BetaStrand)?;
            let bundle = run_peptide_scan(&parsed, scan_mode, &cat)?;
            if let Some(p) = out_csv {
                std::fs::write(&p, render_scan_csv(&bundle))?;
            }
            if let Some(p) = out_sdf {
                let refs: Vec<&MolGraph> = bundle.built_mols.iter().map(|b| &b.mol).collect();
                std::fs::write(&p, write_sdf_v2000_3d(&refs))?;
            }
            if send_dockcraft {
                let dir = default_dockcraft_handoff_dir();
                let mols: Vec<&MolGraph> = bundle.built_mols.iter().map(|b| &b.mol).collect();
                let hb = write_dockcraft_bundle(&mols, &dir, Some("nrp_scan"))?;
                eprintln!("DockCraft scan bundle -> {}", hb.sdf_path.display());
            }
            let payload = serde_json::json!({
                "mode": bundle.mode.label(),
                "parent": bundle.parent,
                "variant_count": bundle.variants.len(),
                "variants": bundle.variants,
            });
            if let Some(p) = out_json {
                std::fs::write(&p, serde_json::to_string_pretty(&payload)?)?;
            }
            if json {
                println!("{}", serde_json::to_string_pretty(&payload)?);
            } else {
                println!("{}", render_scan_csv(&bundle));
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

    #[test]
    fn mcp_build_peptide_export_3d_and_chimerax_bundle_work_inside_sandbox()
    -> Result<(), Box<dyn std::error::Error>> {
        let base = std::env::temp_dir().join("chemcraft_mcp_peptide_chimerax_test");
        let _ = std::fs::remove_dir_all(&base);
        let wdir = base.join("write");
        std::fs::create_dir_all(&wdir)?;

        let sb = PathSandbox::new(None, Some(&wdir))?;
        let req = [
            serde_json::json!({"jsonrpc":"2.0","id":0,"method":"initialize","params":{}})
                .to_string(),
            serde_json::json!({
                "jsonrpc":"2.0","id":1,"method":"tools/call",
                "params":{
                    "name":"build_peptide",
                    "arguments":{
                        "sequence":"Ac-Ala-Gly-Ser-NH2",
                        "secondary_structure":"alpha-helix",
                        "output_pdb":"pep.pdb",
                        "output_svg":"pep.svg"
                    }
                }
            })
            .to_string(),
            serde_json::json!({
                "jsonrpc":"2.0","id":2,"method":"tools/call",
                "params":{
                    "name":"export_3d",
                    "arguments":{
                        "smiles":"CC(=O)Oc1ccccc1C(=O)O",
                        "output_path":"aspirin_3d.sdf",
                        "format":"sdf3d"
                    }
                }
            })
            .to_string(),
            serde_json::json!({
                "jsonrpc":"2.0","id":3,"method":"tools/call",
                "params":{
                    "name":"export_chimerax",
                    "arguments":{
                        "peptide":"Tyr-Gly-Gly-Phe-Met",
                        "secondary_structure":"beta-strand",
                        "stem":"enkephalin"
                    }
                }
            })
            .to_string(),
            serde_json::json!({
                "jsonrpc":"2.0","id":4,"method":"tools/call",
                "params":{
                    "name":"export_3d",
                    "arguments":{
                        "smiles":"CCO",
                        "output_path":"../outside.pdb",
                        "format":"pdb"
                    }
                }
            })
            .to_string(),
            serde_json::json!({
                "jsonrpc":"2.0","id":5,"method":"tools/call",
                "params":{
                    "name":"screen_library",
                    "arguments":{
                        "preset":"fragments",
                        "funnel":"ro3",
                        "diverse_n":8,
                        "output_sdf":"frags_3d.sdf"
                    }
                }
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
        assert_eq!(lines.len(), 6);
        assert_eq!(lines[1]["result"]["isError"], false);
        assert_eq!(lines[2]["result"]["isError"], false);
        assert_eq!(lines[3]["result"]["isError"], false);
        // Path traversal rejected by sandbox
        assert_eq!(lines[4]["result"]["isError"], true);
        assert_eq!(lines[5]["result"]["isError"], false);
        assert_eq!(lines[5]["result"]["structuredContent"]["selected_count"], 8);

        assert!(wdir.join("pep.pdb").exists());
        assert!(wdir.join("pep.svg").exists());
        assert!(wdir.join("aspirin_3d.sdf").exists());
        assert!(wdir.join("frags_3d.sdf").exists());
        assert!(wdir.join("enkephalin.cxc").exists());
        assert!(wdir.join("enkephalin.pdb").exists());
        let cxc = std::fs::read_to_string(wdir.join("enkephalin.cxc"))?;
        assert!(cxc.contains("open enkephalin.pdb"));
        assert!(cxc.contains("cartoon"));

        let _ = std::fs::remove_dir_all(&base);
        Ok(())
    }

    #[test]
    fn mcp_plan_synthesis_exports_eln_markdown_and_enforces_sandbox()
    -> Result<(), Box<dyn std::error::Error>> {
        let base = std::env::temp_dir().join("chemcraft_mcp_plan_synthesis_test");
        let _ = std::fs::remove_dir_all(&base);
        let wdir = base.join("write");
        std::fs::create_dir_all(&wdir)?;

        let sb = PathSandbox::new(None, Some(&wdir))?;
        let req = [
            serde_json::json!({"jsonrpc":"2.0","id":0,"method":"initialize","params":{}})
                .to_string(),
            serde_json::json!({
                "jsonrpc":"2.0","id":1,"method":"tools/call",
                "params":{
                    "name":"plan_synthesis",
                    "arguments":{
                        "smiles":"CC(=O)Nc1ccc(O)cc1",
                        "name":"Paracetamol",
                        "max_steps":3,
                        "scale_mg":250.0,
                        "out_md_path":"paracetamol_eln.md"
                    }
                }
            })
            .to_string(),
            serde_json::json!({
                "jsonrpc":"2.0","id":2,"method":"tools/call",
                "params":{
                    "name":"plan_synthesis",
                    "arguments":{
                        "smiles":"CC(=O)Nc1ccc(O)cc1",
                        "out_md_path":"../outside_eln.md"
                    }
                }
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
        assert_eq!(lines.len(), 3);
        assert_eq!(lines[1]["result"]["isError"], false);
        let routes = lines[1]["result"]["structuredContent"]["routes"]
            .as_array()
            .ok_or("missing routes array")?;
        assert!(!routes.is_empty());
        assert_eq!(routes[0]["fully_resolved_in_stock"], true);
        assert_eq!(lines[2]["result"]["isError"], true);

        let md = std::fs::read_to_string(wdir.join("paracetamol_eln.md"))?;
        assert!(md.contains("Paracetamol"));
        assert!(md.contains("Shopping List"));

        let _ = std::fs::remove_dir_all(&base);
        Ok(())
    }

    #[test]
    fn mcp_nrp_peptide_bgc_and_scan_tools_enforce_sandbox() -> Result<(), Box<dyn std::error::Error>>
    {
        let base = std::env::temp_dir().join("chemcraft_mcp_nrp_test");
        let _ = std::fs::remove_dir_all(&base);
        let wdir = base.join("write");
        std::fs::create_dir_all(&wdir)?;

        let sb = PathSandbox::new(None, Some(&wdir))?;
        let req = [
            serde_json::json!({"jsonrpc":"2.0","id":0,"method":"initialize","params":{}})
                .to_string(),
            serde_json::json!({
                "jsonrpc":"2.0","id":1,"method":"tools/call",
                "params":{
                    "name":"design_peptide_nrp",
                    "arguments":{
                        "spec":"cyclosporin_a",
                        "conformation":"ring",
                        "output_pdb":"csa.pdb",
                        "output_md":"csa.md"
                    }
                }
            })
            .to_string(),
            serde_json::json!({
                "jsonrpc":"2.0","id":2,"method":"tools/call",
                "params":{
                    "name":"analyze_nrps_bgc",
                    "arguments":{
                        "preset":"gramicidin_s",
                        "swap_module":3,
                        "swap_monomer":"Ile",
                        "output_pdb":"grs_ile.pdb"
                    }
                }
            })
            .to_string(),
            serde_json::json!({
                "jsonrpc":"2.0","id":3,"method":"tools/call",
                "params":{
                    "name":"generate_peptide_scan",
                    "arguments":{
                        "spec":"Ac-Arg-Gly-Asp-Phe-Val-NH2",
                        "mode":"nme",
                        "output_csv":"nme_scan.csv",
                        "output_sdf3d":"nme_scan.sdf"
                    }
                }
            })
            .to_string(),
            serde_json::json!({
                "jsonrpc":"2.0","id":4,"method":"tools/call",
                "params":{
                    "name":"design_peptide_nrp",
                    "arguments":{
                        "spec":"cilengitide",
                        "output_pdb":"../escape.pdb"
                    }
                }
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
        assert_eq!(lines.len(), 5);
        assert_eq!(lines[1]["result"]["isError"], false);
        assert_eq!(
            lines[1]["result"]["structuredContent"]["developability"]["formula"],
            "C62H111N11O12"
        );
        assert_eq!(lines[2]["result"]["isError"], false);
        assert_eq!(lines[3]["result"]["isError"], false);
        assert_eq!(lines[3]["result"]["structuredContent"]["variant_count"], 5);
        // Sandbox escape blocked
        assert_eq!(lines[4]["result"]["isError"], true);
        assert!(wdir.join("csa.pdb").exists());
        assert!(wdir.join("csa.md").exists());
        assert!(wdir.join("grs_ile.pdb").exists());
        assert!(wdir.join("nme_scan.csv").exists());
        assert!(wdir.join("nme_scan.sdf").exists());

        let _ = std::fs::remove_dir_all(&base);
        Ok(())
    }
}
