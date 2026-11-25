# ChemCraft v2.0.0

**Local-First, Clean-Room 2D Vector + 3D Molecular Design, Retrosynthesis & SAR Workbench in Pure Rust**  
*(Works side-by-side with [`AdmetCraft`](../AdmetCraft) for deep ADMET/PK, [`DockCraft`](../DockCraft) for 3D SBDD/MD/FEP & [`ChimeraX`](https://www.cgl.ucsf.edu/chimerax/) for macromolecular 3D)*

`ChemCraft v2.0.0` combines ChemDraw/VectorCraft-style 2D chemical drawing, interactive pure-Rust 3D molecular visualization, DataWarrior-style SAR exploration (ECFP4 Tanimoto, SMARTS, MMP & Activity Cliffs), Retro\* synthesis planning, one-click **UCSF ChimeraX / DockCraft / AdmetCraft** handoff — all written in **100% safe, offline pure Rust** (`#![forbid(unsafe_code)]`, zero `.unwrap()` / `.expect()` / `panic!()`).

> **v2.0.0 breaking change:** peptide / NRP / NRPS design has been extracted from ChemCraft. It now lives in the private [`GeneCraft`](https://github.com/acs-foundry/GeneCraft) repository, together with the downstream synthesis-module and gene-design layers. ChemCraft deliberately has no private dependency: interoperability uses public handoff formats (SMILES, HELM, JSON bundles and Craft Link).

---

## ChemCraft vs AdmetCraft vs ChimeraX

| Tool | Primary Focus | Interoperability |
|---|---|---|
| **ChemCraft** | 2D vector drawing, ring snapping, 3D interactive preview, MMP/SAR exploration, SVG/EPS/PDF/PNG/PDB export | Sends SMILES/SAR libraries to **AdmetCraft** & exports 3D `.pdb`/`.sdf` + `.cxc` bundles to **ChimeraX** |
| **AdmetCraft** | Dedicated ADMET, pharmacokinetics (GI/BBB/P-gp/CYP), ESOL/Ali/Silicos-IT solubility, Lipinski/Ghose/Veber/Egan/Muegge, PAINS/Brenk/SA, Bioavailability Radar & BOILED-Egg SVGs | Invoked optionally from ChemCraft via `admetcraft-cli eval/batch` or launched standalone (`admetcraft`) |
| **ChimeraX** | Full macromolecular 3D rendering, publication cartoons, molecular surfaces, electrostatic potentials & H-bond networks | Opened directly from ChemCraft (`Visualizar en ChimeraX` / `chemcraft-cli chimerax`) using sandbox-safe relative `.cxc` bundles |

---

## Architecture (8-Crate Cargo Workspace)

| Layer | Crate / App | Role |
|---|---|---|
| L0 | `crates/mol-graph` | Clean-room editable `MolGraph`, unbounded GF(2) SSSR macrocycle ring basis, iterative Hückel aromaticity, valence diagnostics, Cahn-Ingold-Prelog `R/S` & `E/Z` stereo, canonical SMILES, SMARTS subgraph matching, SDF **V2000 & V3000**, **ChemDraw `.cdxml`** bidirectional reader/writer, **MDL `.rxn` V2000** reaction & stoichiometry/yield engine, **Analysis Window** (exact monoisotopic mass, elemental analysis `%C/%H/...`, DBE, ESI+/ESI- MS adducts), **30 Superatom/Protecting Groups** (`Boc`, `Fmoc`, `Cbz`, `Ph`, `Bn`, `CF3`, `CO2Me`, `SO2Me`, ...), **Name-to-Structure** (110+ drugs/solvents + systematic IUPAC parser), DataWarrior `.dwar`, CSV, bioisostere mutation, pure-Rust **3D conformer generator** (chiral-volume & Jacobi MDS cage embedding, PDB with `SSBOND` & `CONECT`, 3D-SDF V2000+V3000) |
| L1 | `crates/mol-vector` | 2D layout generator (fused/spiro/bridged rings + `Clean Structure 2D`), **ACS Document 1996** (`B&W` / `Color`) & Modern depiction presets with optional terminal `CH3`, ring-snapping templates, ring-aligned double/triple/aromatic bonds, stereo wedges/hashes/wavy/dative arrows, label collision resolution, and **SVG / EPS / FlateDecode PDF / PNG** export (single or multi-molecule grid) |
| L2 | `crates/mol-admet-bridge` | Decoupled `AdmetProvider` trait + standalone pure-Rust ADMET preview engine (**exact Ertl 2000 TPSA**, **Wildman-Crippen 1999 WLOGP/MR**, strict rotatable bonds), **`1H` & `13C` NMR chemical shift/multiplicity predictor**, external `admetcraft-cli` delegation (`eval`/`batch` + BOILED-Egg SVG handoff), and `admetcraft` GUI launcher — zero sibling path dependencies |
| L2 | `crates/mol-sar` | In-memory SAR compound table, 2048-bit ECFP4 (`MorganFp`) fingerprints, `rayon`-parallel Tanimoto/SMARTS/numeric filtering, deterministic Matched Molecular Pairs (MMP), duplicate detection, Activity Cliffs (`SALI`), and Ligand Efficiency (`LE`, `LLE`) |
| L2 | `crates/mol-synth` | **Buy-Before-Make & Retro\* Multi-Step Retrosynthesis Engine**: built-in starter MedChem building-block index (172 starter BBs for offline tests & demos, extensible via user-supplied real supplier `.tsv`/`.smi` catalogs such as Enamine/Mcule/eMolecules in `CRAFT_DATA_DIR/building_blocks`, default `$XDG_DATA_HOME/CraftBio/building_blocks/`), ECFP4 Tanimoto nearest-stock analog finder, 25+ graph-level MedChem disconnections, chemoselectivity & protecting-group alerts, memoized Retro\* AND-OR search, backward-propagated forward stoichiometry (`mg`/`mmol`), MDL `.rxn` V2000 & ELN Markdown protocol export |
| L3 | `crates/mol-engine` | Session state, Copy-on-Write structural Undo/Redo (`Arc` snapshots), `CommandRegistry`, bounded `.chemcraft` ZIP+JSON archives (`O_EXCL` atomic save), **ChimeraX `.pdb`/`.sdf`+`.cxc` bundle launcher**, **DockCraft 3D SDF/SMILES handoff & pH 7.4-aware sidecar merger**, and kernel-confined `PathSandbox` (`openat2` with `RESOLVE_BENEATH \| RESOLVE_NO_SYMLINKS`) |
| L4 | `apps/chemcraft-cli` | Headless CLI (`render --acs1996`, `analyze`, `nmr`, `name`, `group`, `cdxml`, `rxn`, `retro`, `export3d`, `chimerax`, `dockcraft`, `dock-results`, `library`, `filter`, `enrich`, `mmp`, `inspect`) + strict JSON-RPC 2.0 **MCP server** (15 tools over `stdio` & `--tcp` with constant-time bearer auth) |
| L5 | `apps/chemcraft` | Native `egui` / `eframe` ChemDraw-grade desktop studio with `2D / 3D / Split 2D+3D` viewports, optional **SAR Explorer**, on-demand **`Synth Module`** window (hidden by default), and **ChimeraX / DockCraft / AdmetCraft** launchers |

---

## Building & Running

```bash
cargo test --workspace --offline
cargo build --workspace --release --offline
./packaging/install.sh
```

### Desktop Studio
```bash
chemcraft
```

### CLI & MCP Examples
```bash
# Render vector SVG, PDF, EPS, or PNG (including ACS Document 1996 B&W preset)
chemcraft-cli render "CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C" --out testosterone.svg --format svg --acs1996 --monochrome
chemcraft-cli render "CCCc1nn(C)c2c(=O)[nH]c(-c3cc(S(=O)(=O)N4CCN(C)CC4)ccc3OCC)nc12" --out sildenafil.pdf --format pdf

# ChemDraw Analysis Window (exact mass, elemental analysis, DBE, ESI+/ESI- adducts) & 1H/13C NMR
chemcraft-cli analyze "CC(=O)Oc1ccccc1C(=O)O"
chemcraft-cli nmr "CC(=O)Oc1ccccc1C(=O)O"

# Name-to-Structure & Superatom / Protecting Group Attachment
chemcraft-cli name "2-methylpropan-1-ol"
chemcraft-cli group "Nc1ccccc1" --attach Boc --atom 0

# ChemDraw .cdxml bidirectional I/O & MDL .rxn V2000 + stoichiometry/yield
chemcraft-cli cdxml --smiles "CC(C)Cc1ccc([C@@H](C)C(=O)O)cc1" --out ibuprofen.cdxml
chemcraft-cli cdxml --input ibuprofen.cdxml
chemcraft-cli rxn --rsmiles "CC(=O)O.CCO>>CC(=O)OCC.O" --mmol 2.0 --actual-mg 140 --out ester.rxn

# Substructure + Lipinski filter
chemcraft-cli filter library.sdf --substructure "c1ccccc1" --lipinski-only --out hits.sdf

# Enrich an SDF/DWAR/CSV dataset with AdmetCraft and export PDF/CSV
chemcraft-cli enrich library.sdf --with-admet --out admet_report.pdf --format pdf

# Export a 3D conformer from SMILES to PDB or 3D SDF V2000
chemcraft-cli export3d "CC(=O)Oc1ccccc1C(=O)O" --out aspirin_3d.pdb --format pdb

# Export a ChimeraX bundle (.pdb/.sdf + .cxc) and launch ChimeraX
chemcraft-cli chimerax --smiles "CC(C)Cc1ccc(C(C)C(=O)O)cc1" --launch

# Export 3D ligand or SDF library to DockCraft and optionally launch DockCraft
chemcraft-cli dockcraft --smiles "CC(C)Cc1ccc(C(C)C(=O)O)cc1" --launch
chemcraft-cli dockcraft --library hits.sdf --out-dir ./dock_handoff

# Stream, pre-filter (Ro5/Lead-Like/Ro3/BBB/MaxMin), protonate pH 7.4, and prepare 3D SDF libraries
chemcraft-cli library --preset fragments --funnel ro3 --diverse-n 15 --out frags_3d.sdf --send-dockcraft
chemcraft-cli library --input ~/.local/share/chemcraft/libraries/zinc_subset.sdf.gz --funnel lead --out lead_hits_3d.sdf

# Retrosynthesis & Buy-Before-Make (Synth Module CLI)
chemcraft-cli retro --init-stock
chemcraft-cli retro "CC(=O)Nc1ccc(O)cc1" --name Paracetamol --scale-mg 250 --out-md paracetamol_eln.md --out-rxn-dir ./rxn_steps

# Hardened MCP Server (stdio or loopback TCP with token)
chemcraft-cli mcp --read-root ./data --write-root ./exports
chemcraft-cli mcp --read-root ./data --write-root ./exports --tcp 127.0.0.1:7878 --mcp-token "$TOKEN"
```

### MCP Tools Exposed
- `draw_from_smiles`: Computes clean 2D coordinates, valence warnings, deterministic open connectivity key (`CCKEY-...`), and SVG output.
- `mutate_scaffold`: Replaces carboxylic acids (`-COOH` -> `tetrazole`) or attaches substituents and returns the **AdmetCraft** delta (`ΔMW`, `ΔLogP`, `ΔTPSA`, `ΔLogS`, `GI`, `BBB`).
- `search_substructure`: Filters SDF/DWAR/CSV libraries inside `--read-root` by SMARTS and ECFP4 Tanimoto.
- `export_vector_figure`: Exports single molecules or multi-molecule SAR grids (`svg`, `eps`, `pdf`, `png`) inside `--write-root`.
- `mol.analyze_admet`: Runs the built-in ADMET engine (or optional `admetcraft-cli`) and optionally writes CSV/JSON/PDF ADMET reports.
- `export_3d`: Generates 3D conformers from SMILES and writes PDB or 3D SDF inside `--write-root`.
- `export_chimerax`: Exports self-contained ChimeraX `.pdb`/`.sdf` + `.cxc` scripts inside `--write-root`.
- `screen_library`: Streams built-in presets or sandboxed `.sdf`/`.sdf.gz`/`.smi`/`.smi.gz`/`.csv`/`.dwar` chemical libraries through ADMET/SMARTS/MaxMin funnels, ionizes at pH 7.4, generates 3D coordinates, and writes DockCraft-ready 3D SDF files inside `--write-root`.
- `plan_synthesis`: Runs Buy-Before-Make starter/custom catalog lookup + Retro\* multi-step retrosynthesis planning with reagents, backward-propagated stoichiometry, chemoselectivity/PG alerts, and optional ELN Markdown export inside `--write-root`.

---

## Clean-Room & Scientific Scope Notes
- Zero third-party proprietary code (`ChemDraw`, `MarvinSketch`, `DataWarrior`).
- `structural_key()` produces a deterministic topological `CCKEY-...` digest from canonical SMILES (does not invoke external C InChI binaries).
- Live ADMET descriptors run cleanly standalone inside `mol-admet-bridge` (with optional `admetcraft-cli` IPC discovery when enabled).
- Peptide / NRPS design and gene synthesis live in the private `GeneCraft` repository; ChemCraft only exchanges open formats with it.
