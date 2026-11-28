#![forbid(unsafe_code)]

mod chem_panels;
mod library_ui;
mod miniapps;
mod synth_ui;
mod tutorials;
mod viewer3d;
mod widgets;

use std::path::{Path, PathBuf};
use std::time::Instant;

use chem_panels::{
    InspectorTab, ReactionPanelState, render_analysis_tab, render_nmr_tab, render_reaction_tab,
};
use craftbio_ui::craft_dialogs::{Answer, Dialogs, FileDialog};
use craftbio_ui::craft_settings::{Location, Settings, SettingsStore};
use craftbio_ui::{
    AppAccent, CommonUiPrefs, CraftBioTheme, ThemeMode, WindowConfig, build_viewport,
    choose_folder, dialogs_for_egui, open_settings, save_file, status_bar, theme_toggle_button,
    titlebar_pref,
};
use eframe::egui;
use egui::{Color32, Pos2, RichText, Sense, Stroke, Vec2};
use library_ui::{
    LibraryAction, LibraryWindowState, poll_library_docking_job, render_library_contents,
};
use miniapps::{
    ALL_MINI_APPS, CompoundSubTab, EnumeratorState, MiniAppId, WorkspaceWindows, mini_app_meta,
    render_mini_app_header,
};
use mol_admet_bridge::{
    AdmetCraftProvider, AdmetProvider, AdmetReportFormat, launch_admetcraft_gui,
};
use mol_engine::{
    BuiltinLibraryPreset, ChemCraftDocument, ChemSession, Command, CommandRegistry, MAX_READ_BYTES,
    PathSandbox, ScreeningFunnel, decode_chemcraft_zip, is_chimerax_available,
    is_dockcraft_available, launch_chimerax_for_mol, launch_dockcraft_for_mols,
    save_chemcraft_atomic,
};
use mol_graph::{
    BondOrder, BondStereo, SUPERATOM_GROUPS, analyze_molecule, compute_3d_coords,
    generate_3d_coords, parse_csv_records, parse_dwar, parse_rxn_v2000, parse_sdf, write_cdxml,
    write_pdb, write_sdf_v2000_3d,
};
use mol_synth::{E3RecruiterPreset, EnumerateMode};
use mol_vector::{
    DepictionPreset, ReactionArrow, RingTemplate, STANDARD_BOND_LEN, SceneStyle, VectorScene,
    build_mol_scene, scene_to_eps, scene_to_pdf_bytes, scene_to_png_bytes, scene_to_svg,
};
use serde::{Deserialize, Serialize};
use synth_ui::{SynthAction, SynthModuleState, render_synth_contents};
use tutorials::{TutorialAction, TutorialState, TutorialViewMode, built_in_tutorials};
use viewer3d::{Viewer3dState, render_3d_viewport};
use widgets::{
    ScatterXAxis, draw_admet_radar, draw_sar_scatter, paint_mol_thumbnail, paint_vector_scene,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ViewMode {
    View2D,
    View3D,
    Split2D3D,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DrawTool {
    SelectMove,
    BondOrder(BondOrder),
    BondStereo(BondStereo),
    Ring(RingTemplate),
    ChargePlus,
    ChargeMinus,
    Superatom(&'static str),
    ReactionArrow,
    Eraser,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SarTab {
    Table,
    Scatter,
    MmpCliffs,
    Filters,
}

#[derive(Debug, Clone, PartialEq, Default, Serialize, Deserialize)]
#[serde(default)]
struct ChemPrefs {
    ui: CommonUiPrefs,
    export_dir: String,
    show_terminal_methyls: bool,
}

impl Settings for ChemPrefs {
    const APP_ID: &'static str = "chemcraft";
    const SCHEMA_VERSION: u32 = 1;
}

enum ChemDialogAction {
    OpenFile,
    SaveProjectAs,
    ChooseExportDir,
    SaveBytes(Vec<u8>),
}

pub(crate) struct ChemCraftApp {
    pub(crate) session: ChemSession,
    settings: SettingsStore<ChemPrefs>,
    dialogs: Dialogs<ChemDialogAction>,
    tool: DrawTool,
    view_mode: ViewMode,
    depiction_preset: DepictionPreset,
    show_terminal_methyls: bool,
    inspector_tab: InspectorTab,
    reaction_state: ReactionPanelState,
    viewer_3d: Viewer3dState,
    pub(crate) library_win: LibraryWindowState,
    pub(crate) synth_module: SynthModuleState,
    pub(crate) workspace: WorkspaceWindows,
    pub(crate) tutorial_state: TutorialState,
    pub(crate) enumerator: EnumeratorState,
    sar_panel_open: bool,
    chimerax_available: bool,
    dockcraft_available: bool,
    admet_cli_available: bool,
    admet_gui_available: bool,
    active_element: &'static str,
    smiles_bar: String,
    name_bar: String,
    file_path_bar: String,
    export_dir: String,
    current_project_path: Option<PathBuf>,
    dirty: bool,
    status_msg: String,
    dark_mode: bool,
    use_system_titlebar: bool,
    applied_theme: Option<ThemeMode>,
    sar_tab: SarTab,
    scatter_x: ScatterXAxis,
    palette_open: bool,
    palette_query: String,
    hover_atom: Option<usize>,
    drag_start_atom: Option<usize>,
    bond_drag_from: Option<usize>,
    last_canvas_scene: Option<VectorScene>,
    arrow_start: Option<[f64; 2]>,
}

impl ChemCraftApp {
    #[allow(dead_code)]
    pub(crate) fn new() -> Self {
        let ctx = egui::Context::default();
        let (settings, _report) = SettingsStore::<ChemPrefs>::open(Location::Memory);
        Self::with_context(&ctx, settings)
    }

    fn with_context(ctx: &egui::Context, settings: SettingsStore<ChemPrefs>) -> Self {
        let saved = settings.get().clone();
        let initial_mode = ThemeMode::from_dark(saved.ui.dark_mode);
        CraftBioTheme::new(initial_mode, AppAccent::CHEM_AMBER).apply(ctx);

        let mut session = ChemSession::new();
        session.auto_refresh_admet = false;
        let _ = CommandRegistry::execute(
            &mut session,
            Command::LoadCanvasSmiles {
                smiles: "CC(C)Cc1ccc(C(C)C(=O)O)cc1".to_string(),
                name: Some("Ibuprofen".to_string()),
            },
        );
        let export_dir = {
            let saved_dir = saved.export_dir.trim();
            if !saved_dir.is_empty() && Path::new(saved_dir).is_dir() {
                saved_dir.to_string()
            } else {
                std::env::current_dir()
                    .unwrap_or_else(|_| PathBuf::from("."))
                    .display()
                    .to_string()
            }
        };
        let smiles_bar = session.canvas_mol.canonical_smiles();
        let chimerax_available = is_chimerax_available();
        let dockcraft_available = is_dockcraft_available();
        let admet_cli_available = AdmetCraftProvider::is_external_cli_available();
        let admet_gui_available = AdmetCraftProvider::is_external_gui_available();
        session.prefer_external_admet = admet_cli_available;
        Self {
            session,
            settings,
            dialogs: dialogs_for_egui(ctx),
            tool: DrawTool::BondOrder(BondOrder::Single),
            view_mode: ViewMode::Split2D3D,
            depiction_preset: DepictionPreset::ModernColor,
            show_terminal_methyls: saved.show_terminal_methyls,
            inspector_tab: InspectorTab::Analysis,
            reaction_state: ReactionPanelState::default(),
            viewer_3d: Viewer3dState::default(),
            library_win: LibraryWindowState::default(),
            synth_module: SynthModuleState::default(),
            workspace: WorkspaceWindows::default(),
            tutorial_state: TutorialState::default(),
            enumerator: EnumeratorState::default(),
            sar_panel_open: false,
            chimerax_available,
            dockcraft_available,
            admet_cli_available,
            admet_gui_available,
            active_element: "C",
            smiles_bar,
            name_bar: String::new(),
            file_path_bar: String::new(),
            export_dir,
            current_project_path: None,
            dirty: false,
            status_msg: format!(
                "Listo — ChemCraft v{} (Editor 2D/3D · Mini-Apps Móviles · Tutoriales · SAR/PROTAC · Retro*)",
                env!("CARGO_PKG_VERSION")
            ),
            dark_mode: saved.ui.dark_mode,
            use_system_titlebar: saved.ui.use_system_titlebar,
            applied_theme: Some(initial_mode),
            sar_tab: SarTab::Table,
            scatter_x: ScatterXAxis::ConsensusLogP,
            palette_open: false,
            palette_query: String::new(),
            hover_atom: None,
            drag_start_atom: None,
            bond_drag_from: None,
            last_canvas_scene: None,
            arrow_start: None,
        }
    }

    fn process_dialog_answers(&mut self) {
        for (action, answer) in self.dialogs.take_answers() {
            match (action, answer) {
                (ChemDialogAction::OpenFile, Answer::Paths(paths)) => {
                    if let Some(p) = paths.into_iter().next() {
                        let s = p.display().to_string();
                        self.file_path_bar.clone_from(&s);
                        self.import_file(&s);
                    }
                }
                (ChemDialogAction::SaveProjectAs, Answer::Paths(paths)) => {
                    if let Some(p) = paths.into_iter().next() {
                        self.save_project_to_path(&p);
                    }
                }
                (ChemDialogAction::ChooseExportDir, Answer::Paths(paths)) => {
                    if let Some(dir) = paths.into_iter().next() {
                        self.export_dir = dir.display().to_string();
                        self.status_msg = format!("Default folder -> {}", self.export_dir);
                    }
                }
                (ChemDialogAction::SaveBytes(bytes), Answer::Paths(paths)) => {
                    if let Some(target_path) = paths.into_iter().next() {
                        self.write_bytes_to_path(&target_path, &bytes);
                    }
                }
                (_, Answer::Cancelled) => {
                    self.status_msg = "Dialog cancelled".to_string();
                }
                (_, Answer::Unavailable(why)) => {
                    self.status_msg = format!("Native dialog unavailable: {why}");
                }
                (_, Answer::Message(_)) => {}
            }
        }
    }

    fn sync_theme_and_prefs(&mut self, ctx: &egui::Context) {
        let mode = ThemeMode::from_dark(self.dark_mode);
        if self.applied_theme != Some(mode) {
            CraftBioTheme::new(mode, AppAccent::CHEM_AMBER).apply(ctx);
            self.applied_theme = Some(mode);
        }
        self.settings.update(|prefs| {
            prefs.ui.dark_mode = self.dark_mode;
            prefs.ui.use_system_titlebar = self.use_system_titlebar;
            prefs.export_dir.clone_from(&self.export_dir);
            prefs.show_terminal_methyls = self.show_terminal_methyls;
        });
        let now = Instant::now();
        let tick = self.settings.tick(now);
        if let Some(due) = tick.next_deadline {
            ctx.request_repaint_after(due.saturating_duration_since(now));
        }
    }

    fn default_stem(&self) -> String {
        self.session
            .canvas_mol
            .name
            .as_deref()
            .map(|s| {
                s.chars()
                    .map(|c| {
                        if c.is_ascii_alphanumeric() || c == '-' || c == '_' {
                            c
                        } else {
                            '_'
                        }
                    })
                    .collect::<String>()
            })
            .filter(|s| !s.is_empty())
            .unwrap_or_else(|| "chemcraft_structure".to_string())
    }

    fn write_bytes_to_chosen_path(
        &mut self,
        title: &str,
        default_filename: &str,
        filter_label: &str,
        filter_glob: &str,
        bytes: &[u8],
    ) {
        let ext = filter_glob.trim_start_matches("*.").trim_start_matches('.');
        let exts: Vec<&str> = if ext.is_empty() || ext == "*" {
            Vec::new()
        } else {
            vec![ext]
        };
        let dlg = save_file(
            title,
            default_filename,
            filter_label,
            &exts,
            Some(Path::new(&self.export_dir)),
        );
        self.dialogs
            .ask(dlg, ChemDialogAction::SaveBytes(bytes.to_vec()));
    }

    fn write_bytes_to_path(&mut self, target_path: &Path, bytes: &[u8]) {
        let root = target_path.parent().unwrap_or_else(|| Path::new("."));
        let fname = target_path
            .file_name()
            .map(Path::new)
            .unwrap_or(target_path);
        let sb = match PathSandbox::new(None, Some(root)) {
            Ok(s) => s,
            Err(e) => {
                self.status_msg = format!("Export directory error: {e}");
                return;
            }
        };
        match sb.open_write(fname) {
            Ok((resolved, mut f)) => {
                use std::io::Write;
                if let Err(e) = f.write_all(bytes) {
                    self.status_msg = format!("Write error: {e}");
                } else {
                    if let Some(parent) = resolved.parent() {
                        self.export_dir = parent.display().to_string();
                    }
                    self.status_msg = format!("Exported -> {}", resolved.display());
                }
            }
            Err(e) => self.status_msg = format!("Sandbox write error: {e}"),
        }
    }

    fn export_cdxml_file(&mut self) {
        if self.session.canvas_mol.atoms.is_empty() {
            self.status_msg = "Canvas is empty".to_string();
            return;
        }
        let xml = write_cdxml(&self.session.canvas_mol);
        let default_name = format!("{}.cdxml", self.default_stem());
        self.write_bytes_to_chosen_path(
            "Export ChemDraw CDXML",
            &default_name,
            "ChemDraw XML",
            "*.cdxml",
            xml.as_bytes(),
        );
    }

    fn export_smiles_file(&mut self) {
        if self.session.canvas_mol.atoms.is_empty() && self.session.sar_table.rows.is_empty() {
            self.status_msg = "Nothing to export as SMILES".to_string();
            return;
        }
        let items = self.collect_handoff_items();
        let mut buf = String::new();
        for (smi, name) in items {
            buf.push_str(&format!("{smi} {name}\n"));
        }
        let default_name = format!("{}.smi", self.default_stem());
        self.write_bytes_to_chosen_path(
            "Export SMILES (.smi)",
            &default_name,
            "SMILES",
            "*.smi",
            buf.as_bytes(),
        );
    }

    fn export_sar_csv_file(&mut self) {
        if self.session.sar_table.rows.is_empty() {
            self.status_msg = "SAR table is empty".to_string();
            return;
        }
        let csv_cell = |s: &str| -> String {
            let safe = if matches!(s.as_bytes().first(), Some(b'=' | b'+' | b'-' | b'@'))
                && s.parse::<f64>().is_err()
            {
                format!("'{s}")
            } else {
                s.to_string()
            };
            if safe.contains([',', '"', '\n', '\r']) {
                format!("\"{}\"", safe.replace('"', "\"\""))
            } else {
                safe
            }
        };
        let mut buf =
            String::from("ID,Name,Series,Canonical SMILES,MW,pIC50,DockCraft_affinity_kcal\n");
        for r in &self.session.sar_table.rows {
            let pic50 = r.pic50.map_or(String::new(), |v| format!("{v:.3}"));
            let dock = r
                .properties
                .get("DockCraft_affinity_kcal")
                .map(String::as_str)
                .unwrap_or("");
            let line = [
                csv_cell(&r.id.to_string()),
                csv_cell(&r.name),
                csv_cell(&r.series),
                csv_cell(&r.canonical_smiles),
                format!("{:.2}", r.mw),
                pic50,
                csv_cell(dock),
            ]
            .join(",");
            buf.push_str(&line);
            buf.push('\n');
        }
        self.write_bytes_to_chosen_path(
            "Export SAR Table CSV",
            "chemcraft_sar_table.csv",
            "CSV Table",
            "*.csv",
            buf.as_bytes(),
        );
    }

    fn run_cmd(&mut self, cmd: Command) {
        let activates_sar = matches!(
            cmd,
            Command::LoadDemo50Sar
                | Command::AddCanvasToSar { .. }
                | Command::FilterSarByCanvasSubstructure
        );
        match CommandRegistry::execute(&mut self.session, cmd) {
            Ok(out) => {
                if activates_sar {
                    self.workspace.open(MiniAppId::LibrarySarProtac);
                    self.workspace
                        .set_sub_tab(MiniAppId::LibrarySarProtac, CompoundSubTab::Primary);
                }
                self.dirty = true;
                self.status_msg = out.message;
                self.smiles_bar = self.session.canvas_mol.canonical_smiles();
            }
            Err(e) => {
                self.status_msg = format!("Error: {e}");
            }
        }
    }

    fn export_canvas(&mut self, fmt: &str) {
        if self.session.canvas_mol.atoms.is_empty() {
            self.status_msg = "Canvas is empty".to_string();
            return;
        }
        let mut style = SceneStyle::from_preset(self.depiction_preset, 640.0, 460.0, false);
        style.show_terminal_methyls = self.show_terminal_methyls;
        let scene = build_mol_scene(
            &self.session.canvas_mol,
            &style,
            &self.session.reaction_arrows,
        );
        let bytes: Vec<u8> = match fmt {
            "svg" => scene_to_svg(&scene).into_bytes(),
            "eps" => scene_to_eps(&scene).into_bytes(),
            "pdf" => match scene_to_pdf_bytes(&scene) {
                Ok(b) => b,
                Err(e) => {
                    self.status_msg = format!("PDF error: {e}");
                    return;
                }
            },
            "png" => match scene_to_png_bytes(&scene, 2.0) {
                Ok(b) => b,
                Err(e) => {
                    self.status_msg = format!("PNG error: {e}");
                    return;
                }
            },
            _ => return,
        };
        let default_name = format!("{}.{fmt}", self.default_stem());
        let (title, flabel, fglob) = match fmt {
            "svg" => ("Export 2D Vector SVG", "SVG Vector", "*.svg"),
            "eps" => ("Export 2D Vector EPS", "EPS Vector", "*.eps"),
            "pdf" => ("Export 2D Vector PDF", "PDF Document", "*.pdf"),
            _ => ("Export PNG Image", "PNG Image", "*.png"),
        };
        self.write_bytes_to_chosen_path(title, &default_name, flabel, fglob, &bytes);
    }

    fn export_structure_3d(&mut self, fmt: &str) {
        if self.session.canvas_mol.atoms.is_empty() {
            self.status_msg = "Canvas is empty".to_string();
            return;
        }
        let (title, default_name, flabel, fglob, content) = match fmt {
            "pdb" => (
                "Export 3D Structure PDB",
                format!("{}.pdb", self.default_stem()),
                "Protein Data Bank",
                "*.pdb",
                write_pdb(&self.session.canvas_mol),
            ),
            "sdf3d" => (
                "Export 3D Conformer SDF",
                format!("{}_3d.sdf", self.default_stem()),
                "MDL 3D SDF",
                "*.sdf",
                write_sdf_v2000_3d(&[&self.session.canvas_mol]),
            ),
            _ => return,
        };
        self.write_bytes_to_chosen_path(title, &default_name, flabel, fglob, content.as_bytes());
    }

    fn open_in_chimerax(&mut self) {
        match launch_chimerax_for_mol(&self.session.canvas_mol, None) {
            Ok(bundle) => {
                self.status_msg = format!("Launched ChimeraX ({})", bundle.script_path.display());
            }
            Err(e) => {
                self.status_msg = format!("ChimeraX: {e}");
            }
        }
    }

    fn send_to_dockcraft(&mut self, include_sar_library: bool) {
        let mut owned: Vec<mol_graph::MolGraph> = Vec::new();
        if include_sar_library && !self.session.sar_table.rows.is_empty() {
            for r in &self.session.sar_table.rows {
                let mut m = r.mol.clone();
                m.properties
                    .insert("CHEMCRAFT_ROW_ID".to_string(), r.id.to_string());
                m.properties.insert(
                    "CHEMCRAFT_PARENT_SMILES".to_string(),
                    r.canonical_smiles.clone(),
                );
                if let Some(p) = r.pic50 {
                    m.properties.insert("pIC50".to_string(), format!("{p:.3}"));
                }
                mol_graph::ensure_3d_coords(&mut m);
                owned.push(m);
            }
        } else if !self.session.canvas_mol.atoms.is_empty() {
            let mut m = (*self.session.canvas_mol).clone();
            mol_graph::ensure_3d_coords(&mut m);
            owned.push(m);
        }
        if owned.is_empty() {
            self.status_msg =
                "No molecule on canvas or SAR library to send to DockCraft".to_string();
            return;
        }
        let refs: Vec<&mol_graph::MolGraph> = owned.iter().collect();
        let stem = if include_sar_library && owned.len() > 1 {
            Some("chemcraft_screening_library")
        } else {
            None
        };
        match launch_dockcraft_for_mols(&refs, None, stem) {
            Ok(bundle) => {
                self.status_msg = format!(
                    "Sent {} 3D molecule(s) to DockCraft -> {}",
                    bundle.molecule_count,
                    bundle.sdf_path.display()
                );
            }
            Err(e) => {
                self.status_msg = format!("DockCraft: {e}");
            }
        }
    }

    fn collect_handoff_items(&self) -> Vec<(String, String)> {
        let mut items = Vec::new();
        if !self.session.canvas_mol.atoms.is_empty() {
            let smi = self.session.canvas_mol.canonical_smiles();
            if !smi.is_empty() {
                items.push((
                    smi,
                    self.session
                        .canvas_mol
                        .name
                        .clone()
                        .unwrap_or_else(|| "Canvas_Mol".to_string()),
                ));
            }
        }
        for row in &self.session.sar_table.rows {
            if !row.canonical_smiles.is_empty() {
                items.push((row.canonical_smiles.clone(), row.name.clone()));
            }
        }
        items
    }

    fn open_in_admetcraft(&mut self) {
        let items = self.collect_handoff_items();
        if items.is_empty() {
            self.status_msg = "Nothing to send to AdmetCraft".to_string();
            return;
        }
        let dir = PathBuf::from(&self.export_dir);
        match launch_admetcraft_gui(&items, &dir) {
            Ok(smi_path) => {
                self.status_msg = format!(
                    "Launched AdmetCraft GUI (handoff SMILES saved at {})",
                    smi_path.display()
                );
            }
            Err(e) => {
                self.status_msg = format!("AdmetCraft launch: {e}");
            }
        }
    }

    fn export_admet_boiled_egg(&mut self) {
        let items = self.collect_handoff_items();
        if items.is_empty() {
            self.status_msg = "No compounds available for BOILED-Egg plot".to_string();
            return;
        }
        let work_dir = PathBuf::from(&self.export_dir);
        let out_path = work_dir.join("admetcraft_boiled_egg.svg");
        match self.session.admet_provider().run_full_report_with_handoff(
            &items,
            &work_dir,
            &out_path,
            AdmetReportFormat::SvgBoiledEgg,
        ) {
            Ok(msg) => self.status_msg = msg,
            Err(e) => self.status_msg = format!("AdmetCraft BOILED-Egg error: {e}"),
        }
    }

    fn export_admet_report(&mut self, fmt: AdmetReportFormat, filename: &str) {
        let evals: Vec<_> = self
            .session
            .sar_table
            .rows
            .iter()
            .filter_map(|r| r.admet_full.clone())
            .collect();
        if evals.is_empty() {
            self.status_msg = "No ADMET profiles to export".to_string();
            return;
        }
        let root = PathBuf::from(&self.export_dir);
        let sb = match PathSandbox::new(None, Some(&root)) {
            Ok(s) => s,
            Err(e) => {
                self.status_msg = format!("Export root error: {e}");
                return;
            }
        };
        match sb.open_write(std::path::Path::new(filename)) {
            Ok((resolved, mut f)) => {
                match self
                    .session
                    .admet_provider()
                    .export_report_to_writer(&evals, &mut f, fmt)
                {
                    Ok(()) => {
                        self.status_msg = format!("Exported ADMET report -> {}", resolved.display())
                    }
                    Err(e) => self.status_msg = format!("ADMET export error: {e}"),
                }
            }
            Err(e) => self.status_msg = format!("Sandbox error: {e}"),
        }
    }

    fn draw_2d_canvas_at_rect(&mut self, ui: &mut egui::Ui, rect: egui::Rect) {
        let resp = ui.interact(
            rect,
            ui.scope_id().with("chemcraft_2d_canvas"),
            Sense::click_and_drag(),
        );
        let painter = ui.painter_at(rect);
        let bg = if self.dark_mode {
            Color32::from_rgb(22, 24, 30)
        } else {
            Color32::from_rgb(255, 255, 255)
        };
        painter.rect_filled(rect, 4.0, bg);
        painter.rect_stroke(
            rect,
            4.0,
            Stroke::new(1.0_f32, Color32::from_gray(165)),
            egui::StrokeKind::Inside,
        );

        let mut style = SceneStyle::from_preset(
            self.depiction_preset,
            f64::from(rect.width()),
            f64::from(rect.height()),
            self.dark_mode,
        );
        style.padding = 30.0;
        style.show_terminal_methyls = self.show_terminal_methyls;
        let mut scene = build_mol_scene(
            &self.session.canvas_mol,
            &style,
            &self.session.reaction_arrows,
        );
        scene.background = None;
        paint_vector_scene(&painter, &scene, rect.min);

        self.hover_atom = None;
        let mut hover_bond = None;
        if let Some(hp) = resp.hover_pos() {
            let lx = f64::from(hp.x - rect.min.x);
            let ly = f64::from(hp.y - rect.min.y);
            let mut best_d = 16.0_f64;
            for (aid, sp) in scene.atom_screen_pos.iter().enumerate() {
                let d = (sp[0] - lx).hypot(sp[1] - ly);
                if d < best_d {
                    best_d = d;
                    self.hover_atom = Some(aid);
                }
            }
            if self.hover_atom.is_none() {
                let mut best_bd = 11.0_f64;
                for (bid, b) in self.session.canvas_mol.bonds.iter().enumerate() {
                    if let (Some(p1), Some(p2)) = (
                        scene.atom_screen_pos.get(b.u),
                        scene.atom_screen_pos.get(b.v),
                    ) {
                        let d = point_seg_dist([lx, ly], *p1, *p2);
                        if d < best_bd {
                            best_bd = d;
                            hover_bond = Some(bid);
                        }
                    }
                }
            }
        }

        if let Some(aid) = self.hover_atom
            && let Some(sp) = scene.atom_screen_pos.get(aid)
        {
            painter.circle_stroke(
                Pos2::new(rect.min.x + sp[0] as f32, rect.min.y + sp[1] as f32),
                11.0,
                Stroke::new(1.8_f32, Color32::from_rgb(34, 139, 230)),
            );
        } else if let Some(bid) = hover_bond
            && let (Some(p1), Some(p2)) = (
                scene
                    .atom_screen_pos
                    .get(self.session.canvas_mol.bonds[bid].u),
                scene
                    .atom_screen_pos
                    .get(self.session.canvas_mol.bonds[bid].v),
            )
        {
            painter.line_segment(
                [
                    Pos2::new(rect.min.x + p1[0] as f32, rect.min.y + p1[1] as f32),
                    Pos2::new(rect.min.x + p2[0] as f32, rect.min.y + p2[1] as f32),
                ],
                Stroke::new(4.5_f32, Color32::from_rgba_unmultiplied(34, 139, 230, 110)),
            );
        }

        if self.tool == DrawTool::SelectMove {
            if resp.drag_started() {
                self.drag_start_atom = self.hover_atom;
                if self.drag_start_atom.is_some() {
                    self.session.push_undo_checkpoint();
                }
            }
            if resp.dragged()
                && let Some(aid) = self.drag_start_atom
            {
                let delta = resp.drag_delta();
                let scale = scene.world_scale.max(1e-6);
                let mut mol = (*self.session.canvas_mol).clone();
                if let Some(Some(coords)) = mol.atoms.get_mut(aid).map(|a| &mut a.coords_2d) {
                    coords[0] += f64::from(delta.x) / scale;
                    coords[1] -= f64::from(delta.y) / scale;
                }
                self.session.canvas_mol = std::sync::Arc::new(mol);
            }
            if resp.drag_stopped() {
                self.drag_start_atom = None;
            }
        } else if matches!(self.tool, DrawTool::BondOrder(_) | DrawTool::BondStereo(_)) {
            // ChemDraw-style click-and-drag bond creation with 30° angle snapping & ring closure!
            if resp.drag_started() {
                self.bond_drag_from = self.hover_atom;
            }
            if let Some(from_id) = self.bond_drag_from
                && let Some(start_sp) = scene.atom_screen_pos.get(from_id)
                && let Some(ptr) = resp.interact_pointer_pos()
            {
                let start_pos = Pos2::new(
                    rect.min.x + start_sp[0] as f32,
                    rect.min.y + start_sp[1] as f32,
                );
                let target_atom = self.hover_atom.filter(|&a| a != from_id);
                let preview_end = if let Some(tid) = target_atom
                    && let Some(tsp) = scene.atom_screen_pos.get(tid)
                {
                    Pos2::new(rect.min.x + tsp[0] as f32, rect.min.y + tsp[1] as f32)
                } else {
                    let dx = f64::from(ptr.x - start_pos.x);
                    let dy = f64::from(ptr.y - start_pos.y);
                    let snap_step = std::f64::consts::PI / 6.0; // 30 degrees
                    let ang = (dy.atan2(dx) / snap_step).round() * snap_step;
                    let scr_len = (STANDARD_BOND_LEN * scene.world_scale).clamp(24.0, 62.0);
                    Pos2::new(
                        start_pos.x + (scr_len * ang.cos()) as f32,
                        start_pos.y + (scr_len * ang.sin()) as f32,
                    )
                };
                painter.line_segment(
                    [start_pos, preview_end],
                    Stroke::new(2.2_f32, Color32::from_rgb(34, 139, 230)),
                );
                painter.circle_filled(preview_end, 4.0, Color32::from_rgb(34, 139, 230));
            }
            if resp.drag_stopped()
                && let Some(from_id) = self.bond_drag_from.take()
                && let Some(ptr) = resp.interact_pointer_pos()
            {
                let (ord, st) = match self.tool {
                    DrawTool::BondOrder(o) => (o, BondStereo::None),
                    DrawTool::BondStereo(s) => (BondOrder::Single, s),
                    _ => (BondOrder::Single, BondStereo::None),
                };
                if let Some(target_id) = self.hover_atom.filter(|&a| a != from_id) {
                    self.run_cmd(Command::ConnectAtoms {
                        u: from_id,
                        v: target_id,
                        order: ord,
                        stereo: st,
                    });
                } else if let Some(start_sp) = scene.atom_screen_pos.get(from_id) {
                    let start_pos = Pos2::new(
                        rect.min.x + start_sp[0] as f32,
                        rect.min.y + start_sp[1] as f32,
                    );
                    let dx = f64::from(ptr.x - start_pos.x);
                    let dy = f64::from(ptr.y - start_pos.y);
                    if dx.hypot(dy) > 10.0 {
                        let snap_step = std::f64::consts::PI / 6.0;
                        let ang = (dy.atan2(dx) / snap_step).round() * snap_step;
                        let scr_len = (STANDARD_BOND_LEN * scene.world_scale).clamp(24.0, 62.0);
                        let ex = start_sp[0] + scr_len * ang.cos();
                        let ey = start_sp[1] + scr_len * ang.sin();
                        let [wx, wy] = scene.screen_to_world(ex, ey);
                        self.run_cmd(Command::AddAtomAndBond {
                            from_atom: Some(from_id),
                            symbol: self.active_element.to_string(),
                            x: wx,
                            y: wy,
                            order: ord,
                            stereo: st,
                        });
                    }
                }
            }
        }

        if resp.clicked()
            && let Some(cp) = resp.interact_pointer_pos()
        {
            let [wx, wy] =
                scene.screen_to_world(f64::from(cp.x - rect.min.x), f64::from(cp.y - rect.min.y));
            match self.tool {
                DrawTool::Ring(r) => {
                    self.run_cmd(Command::StampRing {
                        template_name: r.label().to_string(),
                        x: wx,
                        y: wy,
                        snap_atom: self.hover_atom,
                        snap_bond: hover_bond,
                    });
                }
                DrawTool::BondOrder(order) => {
                    // ChemDraw UX: clicking an existing bond cycles Single -> Double -> Triple!
                    if self.hover_atom.is_none()
                        && let Some(bid) = hover_bond
                    {
                        self.run_cmd(Command::CycleBondOrder { bond_idx: bid });
                    } else {
                        let new_pos = [wx, wy];
                        self.run_cmd(Command::AddAtomAndBond {
                            from_atom: self.hover_atom.or_else(|| {
                                if self.session.canvas_mol.atoms.is_empty() {
                                    None
                                } else {
                                    Some(self.session.canvas_mol.atoms.len() - 1)
                                }
                            }),
                            symbol: self.active_element.to_string(),
                            x: new_pos[0],
                            y: new_pos[1],
                            order,
                            stereo: BondStereo::None,
                        });
                    }
                }
                DrawTool::BondStereo(stereo) => {
                    if let Some(bid) = hover_bond {
                        self.run_cmd(Command::SetBondStereo {
                            bond_idx: bid,
                            stereo,
                        });
                    } else if let Some(aid) = self.hover_atom {
                        self.run_cmd(Command::AddAtomAndBond {
                            from_atom: Some(aid),
                            symbol: self.active_element.to_string(),
                            x: wx,
                            y: wy,
                            order: BondOrder::Single,
                            stereo,
                        });
                    }
                }
                DrawTool::ChargePlus => {
                    if let Some(aid) = self.hover_atom {
                        self.run_cmd(Command::AdjustAtomCharge {
                            atom_idx: aid,
                            delta: 1,
                        });
                    }
                }
                DrawTool::ChargeMinus => {
                    if let Some(aid) = self.hover_atom {
                        self.run_cmd(Command::AdjustAtomCharge {
                            atom_idx: aid,
                            delta: -1,
                        });
                    }
                }
                DrawTool::Superatom(abbr) => {
                    self.run_cmd(Command::AttachSuperatomGroup {
                        atom_idx: self.hover_atom,
                        abbr: abbr.to_string(),
                    });
                }
                DrawTool::Eraser => {
                    if let Some(aid) = self.hover_atom {
                        self.run_cmd(Command::EraseAtom { atom_idx: aid });
                    } else if let Some(bid) = hover_bond {
                        self.run_cmd(Command::EraseBond { bond_idx: bid });
                    }
                }
                DrawTool::ReactionArrow => {
                    let sx = f64::from(cp.x - rect.min.x);
                    let sy = f64::from(cp.y - rect.min.y);
                    if let Some(start) = self.arrow_start.take() {
                        let mut arrows = (*self.session.reaction_arrows).clone();
                        arrows.push(ReactionArrow {
                            start,
                            end: [sx, sy],
                            label_above: "Reagent".to_string(),
                            label_below: String::new(),
                        });
                        self.session.reaction_arrows = std::sync::Arc::new(arrows);
                    } else {
                        self.arrow_start = Some([sx, sy]);
                    }
                }
                DrawTool::SelectMove => {}
            }
        }
        self.last_canvas_scene = Some(scene);
    }

    fn save_project_to_path(&mut self, path: &std::path::Path) {
        let mut final_path = path.to_path_buf();
        if final_path.extension().is_none() {
            final_path.set_extension("chemcraft");
        }
        let root = final_path
            .parent()
            .unwrap_or_else(|| std::path::Path::new("."));
        let fname = final_path
            .file_name()
            .map(std::path::Path::new)
            .unwrap_or(final_path.as_path());
        let sb = match PathSandbox::new(None, Some(root)) {
            Ok(s) => s,
            Err(e) => {
                self.status_msg = format!("Save directory error: {e}");
                return;
            }
        };
        let target = match sb.validate_write_path(fname) {
            Ok(p) => p,
            Err(e) => {
                self.status_msg = format!("Save sandbox error: {e}");
                return;
            }
        };
        let doc = ChemCraftDocument {
            format: "chemcraft".to_string(),
            version: env!("CARGO_PKG_VERSION").to_string(),
            canvas_mol: (*self.session.canvas_mol).clone(),
            reaction_arrows: (*self.session.reaction_arrows).clone(),
            sar_table: (*self.session.sar_table).clone(),
            sar_filter: self.session.sar_filter.clone(),
        };
        match save_chemcraft_atomic(&target, &doc) {
            Ok(()) => {
                self.current_project_path = Some(target.clone());
                if let Some(parent) = target.parent() {
                    self.export_dir = parent.display().to_string();
                }
                self.dirty = false;
                self.status_msg = format!("Saved project -> {}", target.display());
            }
            Err(e) => self.status_msg = format!("Save error: {e}"),
        }
    }

    fn save_project(&mut self) {
        if let Some(existing) = self.current_project_path.clone() {
            self.save_project_to_path(&existing);
        } else {
            self.save_project_as();
        }
    }

    fn save_project_as(&mut self) {
        let default_name = self
            .current_project_path
            .as_ref()
            .and_then(|p| p.file_name())
            .and_then(|n| n.to_str())
            .map(String::from)
            .unwrap_or_else(|| format!("{}.chemcraft", self.default_stem()));
        let start_dir = self
            .current_project_path
            .as_ref()
            .and_then(|p| p.parent())
            .unwrap_or_else(|| Path::new(&self.export_dir));
        self.dialogs.ask(
            save_file(
                "Save ChemCraft Project (.chemcraft)",
                &default_name,
                "ChemCraft Project",
                &["chemcraft"],
                Some(start_dir),
            ),
            ChemDialogAction::SaveProjectAs,
        );
    }

    fn open_file_dialog(&mut self) {
        let mut dlg = FileDialog::open()
            .title("Open Chemical Structure, Project, or Library")
            .filter(
                "Chemistry files",
                &[
                    "chemcraft",
                    "cdxml",
                    "rxn",
                    "sdf",
                    "mol",
                    "smi",
                    "dwar",
                    "csv",
                ],
            );
        let start = Path::new(&self.export_dir);
        if start.is_dir() {
            dlg = dlg.directory(start);
        }
        self.dialogs.ask(dlg, ChemDialogAction::OpenFile);
    }

    fn import_file(&mut self, path_str: &str) {
        let p = PathBuf::from(path_str.trim());
        let root = p.parent().unwrap_or_else(|| std::path::Path::new("."));
        let sb = match PathSandbox::new(Some(root), None) {
            Ok(s) => s,
            Err(e) => {
                self.status_msg = format!("Import root error: {e}");
                return;
            }
        };
        let rel = p
            .file_name()
            .map(std::path::Path::new)
            .unwrap_or(p.as_path());
        let ext = p
            .extension()
            .and_then(|s| s.to_str())
            .unwrap_or("")
            .to_ascii_lowercase();

        if ext == "chemcraft" {
            match sb.read_bytes_limited(rel, MAX_READ_BYTES).and_then(|b| {
                decode_chemcraft_zip(&b).map_err(|e| {
                    mol_engine::SandboxError::Io(p.display().to_string(), std::io::Error::other(e))
                })
            }) {
                Ok(doc) => {
                    self.session.canvas_mol = std::sync::Arc::new(doc.canvas_mol);
                    self.session.reaction_arrows = std::sync::Arc::new(doc.reaction_arrows);
                    self.session.sar_table = std::sync::Arc::new(doc.sar_table);
                    self.session.sar_filter = doc.sar_filter;
                    self.smiles_bar = self.session.canvas_mol.canonical_smiles();
                    self.session.on_canvas_changed();
                    self.current_project_path = Some(p.clone());
                    self.dirty = false;
                    self.status_msg = format!("Loaded {}", p.display());
                }
                Err(e) => self.status_msg = format!("Load .chemcraft error: {e}"),
            }
            return;
        }

        let (_, text) = match sb.read_to_string_limited(rel, MAX_READ_BYTES) {
            Ok(r) => r,
            Err(e) => {
                self.status_msg = format!("Import sandbox error: {e}");
                return;
            }
        };
        if ext == "cdxml" {
            self.run_cmd(Command::LoadCdxml { xml: text });
            return;
        }
        if ext == "rxn" {
            match parse_rxn_v2000(&text) {
                Ok(rxn) => {
                    self.reaction_state.rsmiles = rxn.reaction_smiles();
                    self.inspector_tab = InspectorTab::Reaction;
                    self.workspace.open(MiniAppId::SynthesisRetro);
                    self.workspace
                        .set_sub_tab(MiniAppId::SynthesisRetro, CompoundSubTab::Secondary);
                    self.status_msg = format!(
                        "Loaded MDL .rxn ({} reactants -> {} products)",
                        rxn.reactants.len(),
                        rxn.products.len()
                    );
                }
                Err(e) => self.status_msg = format!("Parse .rxn error: {e}"),
            }
            return;
        }
        if ext == "smi" || ext == "smiles" {
            let mut lines = text.lines().map(str::trim).filter(|l| !l.is_empty());
            if let Some(first) = lines.next() {
                let mut parts = first.splitn(2, char::is_whitespace);
                let smi = parts.next().unwrap_or("").to_string();
                let nm = parts
                    .next()
                    .map(|s| s.trim().to_string())
                    .filter(|s| !s.is_empty());
                if !smi.is_empty() {
                    self.run_cmd(Command::LoadCanvasSmiles {
                        smiles: smi,
                        name: nm,
                    });
                }
            }
            return;
        }
        let mut table = mol_sar::SarTable::new();
        let res = if ext == "dwar" {
            parse_dwar(&text).map(|recs| {
                for r in recs.into_iter().flatten() {
                    table.add_molecule(r.mol, None, None);
                }
            })
        } else if ext == "csv" {
            parse_csv_records(&text).map(|recs| {
                for r in recs.into_iter().flatten() {
                    table.add_molecule(r.mol, None, None);
                }
            })
        } else {
            for r in parse_sdf(&text).into_iter().flatten() {
                table.add_molecule(r, None, None);
            }
            Ok(())
        };
        match res {
            Ok(()) if !table.rows.is_empty() => {
                self.session.sar_table = std::sync::Arc::new(table);
                self.session.recompute_sar_analytics();
                self.workspace.open(MiniAppId::LibrarySarProtac);
                self.workspace
                    .set_sub_tab(MiniAppId::LibrarySarProtac, CompoundSubTab::Primary);
                self.status_msg = format!(
                    "Imported {} compounds from {}",
                    self.session.sar_table.rows.len(),
                    p.display()
                );
            }
            Ok(()) => self.status_msg = "No valid molecules found in file".to_string(),
            Err(e) => self.status_msg = format!("Parse error: {e}"),
        }
    }

    fn handle_shortcuts(&mut self, ctx: &egui::Context) {
        let wants_kb = ctx.egui_wants_keyboard_input();
        let mut toggle_palette = false;
        let mut do_open = false;
        let mut do_save = false;
        let mut do_save_as = false;
        let mut do_new = false;
        let mut do_undo = false;
        let mut do_redo = false;
        let mut do_clean = false;
        let mut do_demo = false;
        let mut atom_hotkey: Option<&'static str> = None;

        ctx.input(|i| {
            if i.modifiers.ctrl && i.key_pressed(egui::Key::K) {
                toggle_palette = true;
            }
            if i.modifiers.ctrl && i.key_pressed(egui::Key::O) {
                do_open = true;
            }
            if i.modifiers.ctrl && i.modifiers.shift && i.key_pressed(egui::Key::S) {
                do_save_as = true;
            } else if i.modifiers.ctrl && i.key_pressed(egui::Key::S) {
                do_save = true;
            }
            if i.modifiers.ctrl && i.key_pressed(egui::Key::N) {
                do_new = true;
            }
            if !wants_kb && i.modifiers.ctrl && i.key_pressed(egui::Key::Z) {
                do_undo = true;
            }
            if !wants_kb && i.modifiers.ctrl && i.key_pressed(egui::Key::Y) {
                do_redo = true;
            }
            if !wants_kb && i.modifiers.ctrl && i.key_pressed(egui::Key::L) {
                do_clean = true;
            }
            if !wants_kb && i.modifiers.ctrl && i.key_pressed(egui::Key::E) {
                do_demo = true;
            }
            if !wants_kb && !i.modifiers.ctrl {
                let keys: [(egui::Key, &'static str); 8] = [
                    (egui::Key::C, "C"),
                    (egui::Key::N, "N"),
                    (egui::Key::O, "O"),
                    (egui::Key::S, "S"),
                    (egui::Key::P, "P"),
                    (egui::Key::F, "F"),
                    (egui::Key::B, "Br"),
                    (egui::Key::I, "I"),
                ];
                for (k, sym) in keys {
                    if i.key_pressed(k) {
                        atom_hotkey = Some(sym);
                    }
                }
            }
        });

        if toggle_palette {
            self.palette_open = !self.palette_open;
            self.palette_query.clear();
        }
        if do_open {
            self.open_file_dialog();
        }
        if do_save_as {
            self.save_project_as();
        } else if do_save {
            self.save_project();
        }
        if do_new {
            self.current_project_path = None;
            self.run_cmd(Command::ClearCanvas);
        }
        if do_undo {
            self.run_cmd(Command::Undo);
        }
        if do_redo {
            self.run_cmd(Command::Redo);
        }
        if do_clean {
            self.run_cmd(Command::CleanStructure2D);
        }
        if do_demo {
            self.run_cmd(Command::LoadDemo50Sar);
        }
        if let Some(sym) = atom_hotkey {
            self.active_element = sym;
            if let Some(aid) = self.hover_atom {
                self.run_cmd(Command::SetAtomElement {
                    atom_idx: aid,
                    symbol: sym.to_string(),
                });
            }
        }
    }

    fn render_palette(&mut self, ctx: &egui::Context) {
        if !self.palette_open {
            return;
        }
        let mut close = false;
        let mut chosen = None;
        egui::Window::new("ChemCraft Command Palette (Ctrl+K)")
            .collapsible(false)
            .resizable(false)
            .anchor(egui::Align2::CENTER_TOP, [0.0, 64.0])
            .fixed_size([440.0, 290.0])
            .show(ctx, |ui| {
                let resp = ui.add(
                    egui::TextEdit::singleline(&mut self.palette_query)
                        .hint_text("Type a command..."),
                );
                resp.request_focus();
                let q = self.palette_query.trim().to_ascii_lowercase();
                ui.separator();
                egui::ScrollArea::vertical().show(ui, |ui| {
                    for info in CommandRegistry::catalog() {
                        if !q.is_empty() && !info.label.to_ascii_lowercase().contains(&q) {
                            continue;
                        }
                        let lbl = match info.shortcut {
                            Some(sc) => format!("{}   [{sc}]", info.label),
                            None => info.label.to_string(),
                        };
                        if ui.button(lbl).clicked() {
                            chosen = Some(info.id);
                            close = true;
                        }
                    }
                });
                if ui.input(|i| i.key_pressed(egui::Key::Escape)) {
                    close = true;
                }
            });
        if let Some(id) = chosen {
            match id {
                "clean_2d" => self.run_cmd(Command::CleanStructure2D),
                "load_demo_50" => self.run_cmd(Command::LoadDemo50Sar),
                "swap_tetrazole" => self.run_cmd(Command::MutateBioisostere {
                    atom_index: None,
                    group: "tetrazole".to_string(),
                }),
                "add_to_sar" => self.run_cmd(Command::AddCanvasToSar { pic50: None }),
                "filter_canvas" => self.run_cmd(Command::FilterSarByCanvasSubstructure),
                "clear_filter" => self.run_cmd(Command::ClearSarFilter),
                "clear_canvas" => self.run_cmd(Command::ClearCanvas),
                "import_dockcraft_results" => {
                    self.run_cmd(Command::ImportDockCraftResults { path: None });
                }
                "undo" => self.run_cmd(Command::Undo),
                "redo" => self.run_cmd(Command::Redo),
                _ => {}
            }
        }
        if close {
            self.palette_open = false;
        }
    }

    fn handle_library_action(&mut self, action: LibraryAction) {
        match action {
            LibraryAction::None => {}
            LibraryAction::LoadHitOnCanvas(mut mol) => {
                mol_vector::ensure_2d_coords(&mut mol);
                mol_graph::ensure_3d_coords(&mut mol);
                self.session.canvas_mol = std::sync::Arc::new(mol);
                self.smiles_bar = self.session.canvas_mol.canonical_smiles();
                self.session.on_canvas_changed();
                self.dirty = true;
                self.viewer_3d.reset();
                self.status_msg = format!(
                    "Loaded library hit {} (2D + 3D)",
                    self.session.canvas_mol.name.as_deref().unwrap_or("Mol")
                );
            }
            LibraryAction::SendResultToSar(table) => {
                self.session.sar_table = std::sync::Arc::new(table);
                self.session.recompute_sar_analytics();
                self.workspace.open(MiniAppId::LibrarySarProtac);
                self.workspace
                    .set_sub_tab(MiniAppId::LibrarySarProtac, CompoundSubTab::Primary);
                self.dirty = true;
                self.status_msg = format!(
                    "Loaded {} library hits into SAR Explorer",
                    self.session.sar_table.rows.len()
                );
            }
            LibraryAction::ImportLatestDockCraft => {
                if self.session.sar_table.rows.is_empty()
                    && let Some(ref res) = self.library_win.last_result
                {
                    let mut table = mol_sar::SarTable::new();
                    for h in &res.hits {
                        table.add_molecule(h.mol.clone(), None, Some(res.source_label.clone()));
                    }
                    self.session.sar_table = std::sync::Arc::new(table);
                    self.session.recompute_sar_analytics();
                }
                self.workspace.open(MiniAppId::LibrarySarProtac);
                self.workspace
                    .set_sub_tab(MiniAppId::LibrarySarProtac, CompoundSubTab::Primary);
                self.run_cmd(Command::ImportDockCraftResults { path: None });
            }
            LibraryAction::Status(msg) => {
                self.status_msg = msg;
            }
        }
    }

    fn handle_synth_action(&mut self, action: SynthAction) {
        match action {
            SynthAction::None => {}
            SynthAction::LoadSmiles { smiles, name } => {
                self.run_cmd(Command::LoadCanvasSmiles {
                    smiles,
                    name: Some(name.clone()),
                });
                self.viewer_3d.reset();
                self.status_msg = format!("Synth Module: loaded {name} into 2D/3D canvas");
            }
            SynthAction::OpenReactionSmiles { reaction_smiles } => {
                self.reaction_state.rsmiles = reaction_smiles;
                self.inspector_tab = InspectorTab::Reaction;
                self.workspace.open(MiniAppId::SynthesisRetro);
                self.workspace
                    .set_sub_tab(MiniAppId::SynthesisRetro, CompoundSubTab::Secondary);
                self.status_msg =
                    "Synth Module: sent reaction step to Reaction & Stoichiometry tab".to_string();
            }
            SynthAction::Status(msg) => {
                self.status_msg = msg;
            }
        }
    }

    pub(crate) fn run_tutorial_action(
        &mut self,
        track_idx: usize,
        step_idx: usize,
        action: TutorialAction,
    ) {
        let feedback = match action {
            TutorialAction::LoadMolecule2DTemplate {
                smiles,
                name,
                view_mode,
            } => {
                self.run_cmd(Command::LoadCanvasSmiles {
                    smiles: smiles.to_string(),
                    name: Some(name.to_string()),
                });
                self.run_cmd(Command::CleanStructure2D);
                self.view_mode = match view_mode {
                    TutorialViewMode::View2D => ViewMode::View2D,
                    TutorialViewMode::Split2D3D => ViewMode::Split2D3D,
                };
                self.tool = DrawTool::Ring(RingTemplate::Benzene);
                format!("✓ Cargado {name} en lienzo 2D limpio y seleccionada plantilla de Benceno.")
            }
            TutorialAction::RunBioisostereAndInspect => {
                let before_smi = self.session.canvas_mol.canonical_smiles();
                if !before_smi.contains("C(=O)O") && !before_smi.contains("O=C(O)") {
                    self.run_cmd(Command::LoadCanvasSmiles {
                        smiles: "CC(C)Cc1ccc(C(C)C(=O)O)cc1".to_string(),
                        name: Some("Ibuprofen".to_string()),
                    });
                }
                self.run_cmd(Command::MutateBioisostere {
                    atom_index: None,
                    group: "tetrazole".to_string(),
                });
                let _ = self.session.refresh_live_admet();
                self.workspace.open(MiniAppId::PropertiesAdmetCraftBio);
                self.workspace
                    .set_sub_tab(MiniAppId::PropertiesAdmetCraftBio, CompoundSubTab::Tertiary);
                if let Some(ref delta) = self.session.last_mutation_delta {
                    format!(
                        "✓ Bioisóstero aplicado (-COOH → Tetrazol): ΔMW {:+.1}, ΔLogP {:+.2}, ΔTPSA {:+.1}.",
                        delta.delta_mw, delta.delta_logp, delta.delta_tpsa
                    )
                } else {
                    "✓ Bioisóstero aplicado y perfil ADMET actualizado.".to_string()
                }
            }
            TutorialAction::BuildConformerAndOpen3D {
                smiles,
                name,
                view_mode,
                enable_spin,
            } => {
                self.run_cmd(Command::LoadCanvasSmiles {
                    smiles: smiles.to_string(),
                    name: Some(name.to_string()),
                });
                let mut mol = (*self.session.canvas_mol).clone();
                generate_3d_coords(&mut mol);
                self.session.canvas_mol = std::sync::Arc::new(mol);
                self.viewer_3d.reset();
                self.viewer_3d.show_hydrogens = true;
                self.viewer_3d.show_labels = true;
                self.viewer_3d.auto_spin = enable_spin;
                self.view_mode = match view_mode {
                    TutorialViewMode::View2D => ViewMode::View2D,
                    TutorialViewMode::Split2D3D => ViewMode::Split2D3D,
                };
                self.workspace.open(MiniAppId::Viewer3dConformers);
                format!(
                    "✓ Confórmero 3D generado para {name} ({} átomos pesados) y abierto Visor 3D.",
                    self.session.canvas_mol.heavy_atom_count()
                )
            }
            TutorialAction::ComputeAdmetAndOpen => {
                let res = self.session.refresh_live_admet();
                self.workspace.open(MiniAppId::PropertiesAdmetCraftBio);
                self.workspace
                    .set_sub_tab(MiniAppId::PropertiesAdmetCraftBio, CompoundSubTab::Tertiary);
                match res {
                    Ok(()) => "✓ Perfil ADMET bajo demanda calculado y puentes CraftBio listos."
                        .to_string(),
                    Err(e) => format!("Perfil ADMET abierto (nota: {e})."),
                }
            }
            TutorialAction::OpenMiniApp { app, sub_tab } => {
                self.workspace.open(app);
                if let Some(st) = sub_tab {
                    self.workspace.set_sub_tab(app, st);
                }
                format!("✓ Abierta mini-app {}.", mini_app_meta(app).short_label)
            }
            TutorialAction::LoadDemoSarAndShowCliffs => {
                self.run_cmd(Command::LoadDemo50Sar);
                self.sar_tab = SarTab::MmpCliffs;
                self.workspace.open(MiniAppId::LibrarySarProtac);
                self.workspace
                    .set_sub_tab(MiniAppId::LibrarySarProtac, CompoundSubTab::Primary);
                format!(
                    "✓ Cargados {} compuestos SAR: {} Activity Cliffs (SALI) y {} pares MMP detectados.",
                    self.session.sar_table.rows.len(),
                    self.session.cliffs.len(),
                    self.session.mmps.len()
                )
            }
            TutorialAction::RunInSilicoLibraryDemo => {
                self.library_win.use_custom_path = false;
                self.library_win.selected_preset = BuiltinLibraryPreset::FdaApprovedMedicines;
                self.library_win.filter.funnel = ScreeningFunnel::RuleOf5;
                self.library_win.filter.generate_3d = true;
                self.library_win.diversity_count = 12;
                let msg = match self.library_win.run_filter_and_prepare() {
                    Ok(m) => format!("✓ {m}"),
                    Err(e) => format!("Aviso en cribado: {e}"),
                };
                self.workspace.open(MiniAppId::LibrarySarProtac);
                self.workspace
                    .set_sub_tab(MiniAppId::LibrarySarProtac, CompoundSubTab::Secondary);
                msg
            }
            TutorialAction::RunProtacEnumerationDemo => {
                self.enumerator.mode = EnumerateMode::Protac;
                self.enumerator.warhead_smiles = "Cc1ccc(-c2ccnc(N[*])n2)cc1".to_string();
                self.enumerator.e3_preset = E3RecruiterPreset::All;
                self.enumerator.max_compounds = 24;
                let n = self.enumerator.run_enumeration();
                if let Some(first_mol) = self
                    .enumerator
                    .last_library
                    .as_ref()
                    .and_then(|lib| lib.compounds.first())
                    .and_then(|c| c.mol.clone())
                {
                    let mut m = first_mol;
                    mol_vector::ensure_2d_coords(&mut m);
                    mol_graph::ensure_3d_coords(&mut m);
                    self.session.canvas_mol = std::sync::Arc::new(m);
                    self.smiles_bar = self.session.canvas_mol.canonical_smiles();
                    self.session.on_canvas_changed();
                    self.viewer_3d.reset();
                }
                self.workspace.open(MiniAppId::LibrarySarProtac);
                self.workspace
                    .set_sub_tab(MiniAppId::LibrarySarProtac, CompoundSubTab::Tertiary);
                format!(
                    "✓ Ensamblados {n} degradadores PROTAC válidos (CRBN/VHL/IAP) y cargado el primero en 2D+3D."
                )
            }
            TutorialAction::RunRetrosynthesisDemo { smiles, name } => {
                self.run_cmd(Command::LoadCanvasSmiles {
                    smiles: smiles.to_string(),
                    name: Some(name.to_string()),
                });
                let mol_ref = std::sync::Arc::clone(&self.session.canvas_mol);
                self.synth_module.run_for_molecule(&mol_ref);
                self.workspace.open(MiniAppId::SynthesisRetro);
                self.workspace
                    .set_sub_tab(MiniAppId::SynthesisRetro, CompoundSubTab::Primary);
                let routes = self
                    .synth_module
                    .last_plan
                    .as_ref()
                    .map_or(0, |p| p.routes.len());
                format!(
                    "✓ Retro* completado para {name}: {routes} ruta(s) resueltas con bloques comerciales."
                )
            }
            TutorialAction::SendSynthStepToReactionCalc => {
                let rxn_smi = self
                    .synth_module
                    .last_plan
                    .as_ref()
                    .and_then(|p| p.routes.first())
                    .and_then(|r| r.steps.first())
                    .map(|s| s.reaction_smiles.clone())
                    .unwrap_or_else(|| {
                        "O=C(O)c1ccccc1O.CC(=O)OC(=O)C>>CC(=O)Oc1ccccc1C(=O)O.CC(=O)O".to_string()
                    });
                self.reaction_state.rsmiles = rxn_smi;
                self.reaction_state.limiting_mmol = "2.5".to_string();
                self.reaction_state.actual_mg = "420.0".to_string();
                self.inspector_tab = InspectorTab::Reaction;
                self.workspace.open(MiniAppId::SynthesisRetro);
                self.workspace
                    .set_sub_tab(MiniAppId::SynthesisRetro, CompoundSubTab::Secondary);
                "✓ Reacción cargada en la calculadora de Estequiometría (escala 2.5 mmol) lista para exportar .rxn."
                    .to_string()
            }
        };

        self.status_msg = feedback.clone();
        self.tutorial_state
            .mark_completed(track_idx, step_idx, feedback);
    }

    fn ui_left_tools_panel(&mut self, ui: &mut egui::Ui) {
        egui::ScrollArea::vertical().show(ui, |ui| {
            ui.strong("Drawing Tools");
            ui.separator();
            let bond_tools = [
                ("Select / Move", DrawTool::SelectMove),
                (
                    "Single (-) [Drag/Cycle]",
                    DrawTool::BondOrder(BondOrder::Single),
                ),
                ("Double (=)", DrawTool::BondOrder(BondOrder::Double)),
                ("Triple (#)", DrawTool::BondOrder(BondOrder::Triple)),
                ("Aromatic (:)", DrawTool::BondOrder(BondOrder::Aromatic)),
                ("Dative (->)", DrawTool::BondOrder(BondOrder::Dative)),
                ("Up Wedge", DrawTool::BondStereo(BondStereo::Up)),
                ("Down Hash", DrawTool::BondStereo(BondStereo::Down)),
                ("Wavy Stereo", DrawTool::BondStereo(BondStereo::Wavy)),
                ("Formal Charge (+)", DrawTool::ChargePlus),
                ("Formal Charge (-)", DrawTool::ChargeMinus),
                ("Reaction Arrow", DrawTool::ReactionArrow),
                ("Eraser", DrawTool::Eraser),
            ];
            for (lbl, t) in bond_tools {
                if ui.selectable_label(self.tool == t, lbl).clicked() {
                    self.tool = t;
                }
            }

            ui.add_space(6.0);
            ui.strong("Ring Templates (Snap)");
            ui.separator();
            let rings = [
                RingTemplate::Benzene,
                RingTemplate::Cyclohexane,
                RingTemplate::Cyclopentane,
                RingTemplate::Pyridine,
                RingTemplate::Indole,
                RingTemplate::Naphthalene,
                RingTemplate::Tetrazole,
                RingTemplate::Cyclobutane,
                RingTemplate::Cyclopropane,
            ];
            for r in rings {
                if ui
                    .selectable_label(self.tool == DrawTool::Ring(r), r.label())
                    .clicked()
                {
                    self.tool = DrawTool::Ring(r);
                }
            }

            ui.add_space(6.0);
            ui.strong("Superatoms / Groups");
            ui.separator();
            ui.horizontal_wrapped(|ui| {
                for grp in SUPERATOM_GROUPS.iter().take(24) {
                    let active = self.tool == DrawTool::Superatom(grp.abbr);
                    if ui
                        .selectable_label(active, grp.abbr)
                        .on_hover_text(format!("{} — click atom to attach", grp.label))
                        .clicked()
                    {
                        self.tool = DrawTool::Superatom(grp.abbr);
                    }
                }
            });

            ui.add_space(6.0);
            ui.strong("Element Hotkeys");
            ui.separator();
            ui.horizontal_wrapped(|ui| {
                for sym in ["C", "N", "O", "S", "P", "F", "Cl", "Br", "I"] {
                    if ui
                        .selectable_label(self.active_element == sym, sym)
                        .clicked()
                    {
                        self.active_element = sym;
                    }
                }
            });
            ui.label(
                RichText::new("Hover atom + key changes element")
                    .small()
                    .weak(),
            );

            ui.add_space(8.0);
            ui.strong("Import (.cdxml/.rxn/.sdf/.dwar/.csv)");
            ui.add(
                egui::TextEdit::singleline(&mut self.file_path_bar)
                    .hint_text("/path/to/file.cdxml"),
            );
            ui.horizontal(|ui| {
                if ui.button("Browse...").clicked() {
                    self.open_file_dialog();
                }
                if ui.button("Import").clicked() && !self.file_path_bar.trim().is_empty() {
                    let p = self.file_path_bar.clone();
                    self.import_file(&p);
                }
            });
        });
    }

    fn ui_admet_and_craftbio_contents(&mut self, ui: &mut egui::Ui, id_prefix: &str) {
        ui.strong("AdmetCraft & ADMET (On-Demand)");
        ui.label(
            RichText::new(
                "ADMET no se ejecuta automáticamente al dibujar; abre AdmetCraft para el estudio completo (QED, CNS MPO, LogD7.4, hERG/AMES/DILI, BOILED-Egg) o calcula una vista rápida aquí bajo demanda.",
            )
            .small()
            .weak(),
        );
        ui.add_space(6.0);
        ui.horizontal_wrapped(|ui| {
            if ui
                .button(RichText::new("Abrir en AdmetCraft").strong())
                .on_hover_text(
                    "Guarda el lote/molécula actual y abre la aplicación dedicada AdmetCraft con los datos cargados",
                )
                .clicked()
            {
                self.open_in_admetcraft();
            }
            if ui.button("Calcular ADMET aquí").clicked() {
                match self.session.refresh_live_admet() {
                    Ok(()) => {
                        self.status_msg = "Computed on-demand ADMET summary".to_string();
                    }
                    Err(e) => {
                        self.status_msg = format!("ADMET error: {e}");
                    }
                }
            }
            if ui
                .button("Bioisóstero: -COOH → Tetrazol")
                .on_hover_text("Reemplaza un grupo carboxilo por tetrazol y calcula el delta ADMET")
                .clicked()
            {
                self.run_cmd(Command::MutateBioisostere {
                    atom_index: None,
                    group: "tetrazole".to_string(),
                });
                let _ = self.session.refresh_live_admet();
            }
        });
        ui.horizontal_wrapped(|ui| {
            if self.admet_cli_available {
                ui.checkbox(
                    &mut self.session.prefer_external_admet,
                    "Usar `admetcraft-cli`",
                );
            }
            ui.checkbox(
                &mut self.session.auto_refresh_admet,
                "Auto-actualizar al editar",
            );
        });
        ui.separator();

        if let Some(ref eval) = self.session.live_admet {
            draw_admet_radar(ui, eval, self.dark_mode);
            ui.separator();
            let p = &eval.physicochemical;
            let l = &eval.lipophilicity;
            let pk = &eval.pharmacokinetics;
            let d = &eval.druglikeness;
            let mc = &eval.medicinal_chemistry;

            egui::Grid::new((id_prefix, "admet_kv"))
                .spacing([10.0, 4.0])
                .striped(true)
                .show(ui, |ui| {
                    ui.label("Formula");
                    ui.monospace(&p.formula);
                    ui.end_row();
                    ui.label("MW");
                    ui.label(format!("{:.2} g/mol", p.molecular_weight));
                    ui.end_row();
                    ui.label("Consensus LogP");
                    ui.label(format!("{:.2}", l.consensus_log_p));
                    ui.end_row();
                    ui.label("TPSA");
                    ui.label(format!("{:.2} A^2", p.tpsa));
                    ui.end_row();
                    ui.label("HBD / HBA");
                    ui.label(format!("{} / {}", p.h_bond_donors, p.h_bond_acceptors));
                    ui.end_row();
                    ui.label("ESOL LogS");
                    ui.label(format!(
                        "{:.2} ({})",
                        eval.water_solubility.esol.log_s,
                        eval.water_solubility.esol.class.as_str()
                    ));
                    ui.end_row();
                    ui.label("GI Absorption");
                    ui.strong(pk.gi_absorption.as_str());
                    ui.end_row();
                    ui.label("BBB Permeant");
                    ui.label(if pk.bbb_permeant { "Yes" } else { "No" });
                    ui.end_row();
                    ui.label("Lipinski");
                    ui.label(format!("{} violation(s)", d.lipinski.violations));
                    ui.end_row();
                    ui.label("PAINS / Brenk");
                    ui.label(format!(
                        "{} / {} alert(s)",
                        mc.pains_alerts.len(),
                        mc.brenk_alerts.len()
                    ));
                    ui.end_row();
                });

            if let Some(ref delta) = self.session.last_mutation_delta {
                ui.add_space(6.0);
                ui.group(|ui| {
                    ui.strong("Bioisostere ADMET Delta");
                    ui.label(format!(
                        "ΔMW: {:+.2} | ΔLogP: {:+.2} | ΔTPSA: {:+.2}",
                        delta.delta_mw, delta.delta_logp, delta.delta_tpsa
                    ));
                    ui.label(format!("GI: {} -> {}", delta.old_gi, delta.new_gi));
                });
            }
        } else {
            ui.label(
                RichText::new(
                    "Sin cálculo ADMET activo. Pulsa 'Abrir en AdmetCraft' o 'Calcular ADMET aquí' cuando lo necesites.",
                )
                .weak(),
            );
        }

        ui.add_space(10.0);
        ui.separator();
        ui.strong("Exportar informes ADMET");
        ui.horizontal_wrapped(|ui| {
            if ui.button("ADMET CSV").clicked() {
                self.export_admet_report(AdmetReportFormat::Csv, "chemcraft_admet.csv");
            }
            if ui.button("ADMET PDF").clicked() {
                self.export_admet_report(AdmetReportFormat::Pdf, "chemcraft_admet.pdf");
            }
            if ui
                .add_enabled(
                    self.admet_cli_available,
                    egui::Button::new("BOILED-Egg SVG"),
                )
                .clicked()
            {
                self.export_admet_boiled_egg();
            }
        });

        ui.add_space(10.0);
        ui.separator();
        ui.strong("Puentes CraftBio y Exportación Directa");
        ui.small(format!(
            "Detectados: AdmetCraft GUI={} | DockCraft={} | ChimeraX={}",
            if self.admet_gui_available {
                "sí"
            } else {
                "bundle"
            },
            if self.dockcraft_available {
                "sí"
            } else {
                "bundle"
            },
            if self.chimerax_available {
                "sí"
            } else {
                "bundle"
            }
        ));
        ui.horizontal_wrapped(|ui| {
            if ui.button("Enviar molécula a DockCraft (3D)").clicked() {
                self.send_to_dockcraft(false);
            }
            if !self.session.sar_table.rows.is_empty()
                && ui.button("Enviar biblioteca SAR a DockCraft").clicked()
            {
                self.send_to_dockcraft(true);
            }
            if ui.button("Abrir en ChimeraX").clicked() {
                self.open_in_chimerax();
            }
            if ui.button("Exportar CDXML").clicked() {
                self.export_cdxml_file();
            }
            if ui.button("Exportar 3D SDF").clicked() {
                self.export_structure_3d("sdf3d");
            }
            if ui.button("Exportar 3D PDB").clicked() {
                self.export_structure_3d("pdb");
            }
        });
    }

    fn ui_sar_explorer_contents(&mut self, ui: &mut egui::Ui, id_prefix: &str) {
        let matching = self
            .session
            .sar_table
            .matching_indices(&self.session.sar_filter)
            .unwrap_or_else(|_| (0..self.session.sar_table.rows.len()).collect());

        ui.horizontal_wrapped(|ui| {
            ui.strong(format!(
                "SAR Explorer ({}/{} shown)",
                matching.len(),
                self.session.sar_table.rows.len()
            ));
            if ui.button("Demo 50 SAR").clicked() {
                self.run_cmd(Command::LoadDemo50Sar);
            }
            if ui.button("+ Añadir lienzo al SAR").clicked() {
                self.run_cmd(Command::AddCanvasToSar { pic50: None });
            }
            if !self.session.sar_table.rows.is_empty()
                && ui.button("Enviar biblioteca a DockCraft").clicked()
            {
                self.send_to_dockcraft(true);
            }
            if !self.session.sar_table.rows.is_empty() {
                let import_btn = ui
                    .button("Importar resultados de DockCraft")
                    .on_hover_text(
                        "Trae el último *.dockcraft.json que DockCraft haya dejado y añade sus columnas de acoplamiento junto a las medidas. Se puede deshacer.",
                    );
                if import_btn.clicked() {
                    self.run_cmd(Command::ImportDockCraftResults { path: None });
                }
                if ui.small_button("Exportar CSV").clicked() {
                    self.export_sar_csv_file();
                }
            }
            ui.separator();
            ui.selectable_value(&mut self.sar_tab, SarTab::Table, "SAR Table");
            ui.selectable_value(&mut self.sar_tab, SarTab::Scatter, "Scatter Plot");
            ui.selectable_value(
                &mut self.sar_tab,
                SarTab::MmpCliffs,
                "MMP & Activity Cliffs",
            );
            ui.selectable_value(&mut self.sar_tab, SarTab::Filters, "Dynamic Filters");
        });
        ui.separator();

        match self.sar_tab {
            SarTab::Table => {
                egui::ScrollArea::both()
                    .id_salt((id_prefix, "sar_table_scroll"))
                    .show(ui, |ui| {
                        let mut select_row = None;
                        egui::Grid::new((id_prefix, "sar_grid"))
                            .striped(true)
                            .show(ui, |ui| {
                                for h in [
                                    "2D Structure",
                                    "Compound",
                                    "Series",
                                    "pIC50",
                                    "Dock",
                                    "LE",
                                    "MW",
                                    "cLogP",
                                    "TPSA",
                                    "GI",
                                    "PAINS",
                                ] {
                                    ui.strong(h);
                                }
                                ui.end_row();

                                for &idx in matching.iter().take(100) {
                                    let row = &self.session.sar_table.rows[idx];
                                    let _ = paint_mol_thumbnail(
                                        ui,
                                        &row.mol,
                                        &[],
                                        Vec2::new(84.0, 52.0),
                                        self.dark_mode,
                                    );
                                    let is_sel = self.session.selected_sar_idx == Some(idx);
                                    if ui.selectable_label(is_sel, &row.name).clicked() {
                                        select_row = Some(idx);
                                    }
                                    ui.label(&row.series);
                                    ui.label(
                                        row.pic50.map_or("-".to_string(), |v| format!("{v:.2}")),
                                    );
                                    match row.properties.get("DockCraft_affinity_kcal") {
                                        Some(v) => {
                                            let label = ui.label(v);
                                            if let Some(pose) =
                                                row.properties.get("DockCraft_pose_sdf")
                                            {
                                                label.on_hover_text(format!(
                                                    "{} kcal/mol{}\nPose: {pose}",
                                                    v,
                                                    row.properties
                                                        .get("DockCraft_receptor")
                                                        .map_or_else(String::new, |r| format!(
                                                            " vs {r}"
                                                        ))
                                                ));
                                            }
                                        }
                                        None => {
                                            ui.label("-");
                                        }
                                    }
                                    ui.label(
                                        row.ligand_efficiency()
                                            .map_or("-".to_string(), |v| format!("{v:.2}")),
                                    );
                                    ui.label(format!("{:.1}", row.mw));
                                    if let Some(ref a) = row.admet {
                                        ui.label(format!("{:.2}", a.consensus_logp));
                                        ui.label(format!("{:.1}", a.tpsa));
                                        ui.label(&a.gi_absorption);
                                        ui.label(a.pains_count.to_string());
                                    } else {
                                        ui.label("-");
                                        ui.label("-");
                                        ui.label("-");
                                        ui.label("-");
                                    }
                                    ui.end_row();
                                }
                            });
                        if let Some(idx) = select_row {
                            self.run_cmd(Command::SelectSarRow { idx });
                        }
                    });
            }
            SarTab::Scatter => {
                ui.horizontal(|ui| {
                    ui.label("X Axis:");
                    for ax in [
                        ScatterXAxis::ConsensusLogP,
                        ScatterXAxis::Tpsa,
                        ScatterXAxis::MolecularWeight,
                        ScatterXAxis::LigandEfficiency,
                    ] {
                        ui.selectable_value(&mut self.scatter_x, ax, ax.label());
                    }
                });
                if let Some(clicked) = draw_sar_scatter(
                    ui,
                    &self.session.sar_table.rows,
                    &matching,
                    self.scatter_x,
                    self.session.selected_sar_idx,
                    self.dark_mode,
                ) {
                    self.run_cmd(Command::SelectSarRow { idx: clicked });
                }
            }
            SarTab::MmpCliffs => {
                egui::ScrollArea::vertical()
                    .id_salt((id_prefix, "sar_cliffs_scroll"))
                    .show(ui, |ui| {
                        ui.strong(format!(
                            "Activity Cliffs (SALI) — {} pairs",
                            self.session.cliffs.len()
                        ));
                        egui::Grid::new((id_prefix, "cliffs_grid"))
                            .striped(true)
                            .show(ui, |ui| {
                                ui.strong("Mol A");
                                ui.strong("Mol B");
                                ui.strong("Tanimoto");
                                ui.strong("|ΔpIC50|");
                                ui.strong("SALI");
                                ui.end_row();
                                for c in self.session.cliffs.iter().take(25) {
                                    ui.label(&c.mol_a_name);
                                    ui.label(&c.mol_b_name);
                                    ui.label(format!("{:.2}", c.tanimoto));
                                    ui.label(format!("{:.2}", c.delta_pic50));
                                    ui.strong(format!("{:.1}", c.sali));
                                    ui.end_row();
                                }
                            });
                        ui.add_space(8.0);
                        ui.strong(format!(
                            "Matched Molecular Pairs — {} transforms",
                            self.session.mmps.len()
                        ));
                        egui::Grid::new((id_prefix, "mmp_grid"))
                            .striped(true)
                            .show(ui, |ui| {
                                ui.strong("Mol A -> Mol B");
                                ui.strong("R1 -> R2");
                                ui.strong("ΔpIC50");
                                ui.strong("ΔcLogP");
                                ui.strong("ΔTPSA");
                                ui.end_row();
                                for m in self.session.mmps.iter().take(25) {
                                    ui.label(format!("{} -> {}", m.mol_a_name, m.mol_b_name));
                                    ui.monospace(format!("{} -> {}", m.r_group_a, m.r_group_b));
                                    ui.label(
                                        m.delta_pic50
                                            .map_or("-".to_string(), |v| format!("{v:+.2}")),
                                    );
                                    ui.label(
                                        m.delta_clogp
                                            .map_or("-".to_string(), |v| format!("{v:+.2}")),
                                    );
                                    ui.label(
                                        m.delta_tpsa
                                            .map_or("-".to_string(), |v| format!("{v:+.1}")),
                                    );
                                    ui.end_row();
                                }
                            });
                    });
            }
            SarTab::Filters => {
                ui.horizontal_wrapped(|ui| {
                    ui.label("SMARTS Substructure:");
                    let mut sm = self
                        .session
                        .sar_filter
                        .smarts_substructure
                        .clone()
                        .unwrap_or_default();
                    if ui.text_edit_singleline(&mut sm).changed() {
                        self.session.sar_filter.smarts_substructure =
                            if sm.trim().is_empty() { None } else { Some(sm) };
                    }
                    if ui.button("Use Canvas").clicked() {
                        self.run_cmd(Command::FilterSarByCanvasSubstructure);
                    }
                    if ui.button("Clear All").clicked() {
                        self.run_cmd(Command::ClearSarFilter);
                    }
                });
                ui.horizontal_wrapped(|ui| {
                    ui.label("Tanimoto Reference SMILES:");
                    let mut rf = self
                        .session
                        .sar_filter
                        .reference_smiles
                        .clone()
                        .unwrap_or_else(|| self.session.canvas_mol.canonical_smiles());
                    if ui.text_edit_singleline(&mut rf).changed() {
                        self.session.sar_filter.reference_smiles = Some(rf);
                    }
                    ui.add(
                        egui::Slider::new(&mut self.session.sar_filter.min_tanimoto, 0.0..=1.0)
                            .text("Min Tanimoto"),
                    );
                });
                ui.horizontal_wrapped(|ui| {
                    ui.checkbox(
                        &mut self.session.sar_filter.lipinski_pass_only,
                        "Lipinski 0 violations",
                    );
                    ui.checkbox(
                        &mut self.session.sar_filter.no_pains_only,
                        "No PAINS alerts",
                    );
                });
            }
        }
    }

    fn ui_enumerator_contents(&mut self, ui: &mut egui::Ui) {
        ui.strong("Enumerador Combinatorio R-Group / Markush y Ensamblador PROTAC");
        ui.small("Genera bibliotecas deduplicadas y validadas con puntos de unión [*] / [R1] / [R2] o ensambla degradadores bifuncionales (Warhead + Linker + CRBN/VHL/IAP).");
        ui.add_space(4.0);

        ui.horizontal_wrapped(|ui| {
            ui.label("Modo:");
            ui.selectable_value(
                &mut self.enumerator.mode,
                EnumerateMode::RGroup,
                "R-Group / Markush",
            );
            ui.selectable_value(
                &mut self.enumerator.mode,
                EnumerateMode::Protac,
                "PROTAC / TPD Bifuncional",
            );
            ui.separator();
            ui.label("Máx. compuestos:");
            ui.add(egui::DragValue::new(&mut self.enumerator.max_compounds).range(1..=200));
        });

        match self.enumerator.mode {
            EnumerateMode::RGroup => {
                ui.horizontal_wrapped(|ui| {
                    ui.label("Scaffold SMILES:");
                    ui.add(
                        egui::TextEdit::singleline(&mut self.enumerator.scaffold_smiles)
                            .desired_width(300.0)
                            .font(egui::TextStyle::Monospace)
                            .hint_text("c1ccc([*])cc1 o c1cc([R1])ccc1[R2]"),
                    );
                });
                ui.horizontal_wrapped(|ui| {
                    if ui.small_button("Preset: Fenilo [*]").clicked() {
                        self.enumerator.scaffold_smiles = "c1ccc([*])cc1".to_string();
                    }
                    if ui
                        .small_button("Preset: Quinasa 2-Aminopirimidina")
                        .clicked()
                    {
                        self.enumerator.scaffold_smiles = "Cc1ccc(-c2ccnc(N[*])n2)cc1".to_string();
                    }
                    if ui.small_button("Preset: Doble [R1]/[R2]").clicked() {
                        self.enumerator.scaffold_smiles = "c1cc([R1])ccc1C(=O)N[R2]".to_string();
                    }
                    if ui.small_button("Usar lienzo + [*]").clicked()
                        && !self.session.canvas_mol.atoms.is_empty()
                    {
                        self.enumerator.scaffold_smiles =
                            format!("({})[*]", self.session.canvas_mol.canonical_smiles());
                    }
                });
            }
            EnumerateMode::Protac => {
                ui.horizontal_wrapped(|ui| {
                    ui.label("Warhead POI SMILES:");
                    ui.add(
                        egui::TextEdit::singleline(&mut self.enumerator.warhead_smiles)
                            .desired_width(300.0)
                            .font(egui::TextStyle::Monospace)
                            .hint_text("Cc1ccc(-c2ccnc(N[*])n2)cc1"),
                    );
                });
                ui.horizontal_wrapped(|ui| {
                    ui.label("Reclutador E3:");
                    for (preset, lbl) in [
                        (E3RecruiterPreset::All, "Todos (CRBN+VHL+IAP)"),
                        (E3RecruiterPreset::Crbn, "CRBN (Poma/Lena)"),
                        (E3RecruiterPreset::Vhl, "VHL (VH032)"),
                        (E3RecruiterPreset::Iap, "cIAP"),
                    ] {
                        ui.selectable_value(&mut self.enumerator.e3_preset, preset, lbl);
                    }
                });
                ui.horizontal_wrapped(|ui| {
                    if ui.small_button("Preset: Warhead Quinasa").clicked() {
                        self.enumerator.warhead_smiles = "Cc1ccc(-c2ccnc(N[*])n2)cc1".to_string();
                    }
                    if ui.small_button("Preset: Warhead Aril-amida").clicked() {
                        self.enumerator.warhead_smiles = "c1ccc(NC(=O)c2ccncc2[*])cc1".to_string();
                    }
                    if ui.small_button("Usar lienzo").clicked()
                        && !self.session.canvas_mol.atoms.is_empty()
                    {
                        self.enumerator.warhead_smiles = self.session.canvas_mol.canonical_smiles();
                    }
                });
            }
        }

        ui.add_space(4.0);
        let mut send_all_to_sar = false;
        let mut load_enum_mol: Option<mol_graph::MolGraph> = None;
        let mut send_enum_to_dockcraft = false;

        ui.horizontal_wrapped(|ui| {
            let run_btn = egui::Button::new(
                RichText::new("▶ Enumerar Biblioteca")
                    .strong()
                    .color(Color32::WHITE),
            )
            .fill(Color32::from_rgb(79, 70, 229));
            if ui.add(run_btn).clicked() {
                let n = self.enumerator.run_enumeration();
                self.status_msg = format!("Enumerador: {n} compuestos válidos generados.");
            }
            if self
                .enumerator
                .last_library
                .as_ref()
                .is_some_and(|l| !l.compounds.is_empty())
            {
                if ui.button("Enviar serie al Explorador SAR").clicked() {
                    send_all_to_sar = true;
                }
                if ui.button("Enviar serie a DockCraft (3D)").clicked() {
                    send_enum_to_dockcraft = true;
                }
            }
        });
        ui.small(&self.enumerator.status_msg);
        ui.separator();

        if let Some(ref lib) = self.enumerator.last_library {
            ui.strong(format!(
                "Resultado: {} compuestos únicos ({} intentos)",
                lib.valid_unique, lib.attempted
            ));
            egui::ScrollArea::horizontal()
                .id_salt("enum_lib_hscroll")
                .show(ui, |ui| {
                    egui::Grid::new("enum_lib_grid")
                        .striped(true)
                        .show(ui, |ui| {
                            for h in [
                                "Acción",
                                "ID",
                                "Fórmula",
                                "MW",
                                "R1 / Linker",
                                "R2 / E3",
                                "Regla",
                            ] {
                                ui.strong(h);
                            }
                            ui.end_row();
                            for c in lib.compounds.iter().take(120) {
                                if ui.small_button("Ver 2D+3D").clicked()
                                    && let Some(ref m) = c.mol
                                {
                                    load_enum_mol = Some(m.clone());
                                }
                                ui.label(&c.id);
                                ui.monospace(&c.formula);
                                ui.label(format!("{:.1}", c.mw));
                                ui.label(&c.linker_or_r1);
                                ui.label(&c.recruiter_or_r2);
                                let badge = if c.ro5_mw_ok {
                                    "Ro5"
                                } else if c.bro5_ok {
                                    "bRo5"
                                } else {
                                    ">1050 Da"
                                };
                                ui.label(badge);
                                ui.end_row();
                            }
                        });
                });
        }

        if let Some(mol) = load_enum_mol {
            self.handle_library_action(LibraryAction::LoadHitOnCanvas(mol));
        }
        if send_all_to_sar && let Some(ref lib) = self.enumerator.last_library {
            let mut tbl = mol_sar::SarTable::new();
            let series_lbl = match lib.mode {
                EnumerateMode::RGroup => "RGroup_Enum".to_string(),
                EnumerateMode::Protac => "PROTAC_Enum".to_string(),
            };
            for c in &lib.compounds {
                if let Some(ref m) = c.mol {
                    tbl.add_molecule(m.clone(), None, Some(series_lbl.clone()));
                }
            }
            self.handle_library_action(LibraryAction::SendResultToSar(tbl));
        }
        if send_enum_to_dockcraft && let Some(ref lib) = self.enumerator.last_library {
            let mut mols: Vec<mol_graph::MolGraph> = Vec::new();
            for c in &lib.compounds {
                if let Some(ref m) = c.mol {
                    let mut m3 = m.clone();
                    mol_graph::ensure_3d_coords(&mut m3);
                    mols.push(m3);
                }
            }
            if !mols.is_empty() {
                let refs: Vec<&mol_graph::MolGraph> = mols.iter().collect();
                match launch_dockcraft_for_mols(&refs, None, Some("chemcraft_enumerated")) {
                    Ok(b) => {
                        self.status_msg = format!(
                            "Sent {} enumerated 3D molecule(s) to DockCraft -> {}",
                            b.molecule_count,
                            b.sdf_path.display()
                        );
                    }
                    Err(e) => {
                        self.status_msg = format!("DockCraft: {e}");
                    }
                }
            }
        }
    }

    fn ui_viewer3d_miniapp(&mut self, ui: &mut egui::Ui) {
        ui.horizontal_wrapped(|ui| {
            ui.checkbox(&mut self.viewer_3d.show_hydrogens, "Mostrar +H");
            ui.checkbox(&mut self.viewer_3d.show_labels, "Etiquetas");
            ui.checkbox(&mut self.viewer_3d.show_ribbon, "Cinta peptídica");
            ui.checkbox(&mut self.viewer_3d.auto_spin, "Auto-rotar");
            if ui.small_button("Reconstruir 3D (EtKDG)").clicked() {
                let mut mol = (*self.session.canvas_mol).clone();
                generate_3d_coords(&mut mol);
                self.session.canvas_mol = std::sync::Arc::new(mol);
                self.status_msg = "Confórmero 3D EtKDG regenerado.".to_string();
            }
            if ui.small_button("Centrar cámara").clicked() {
                self.viewer_3d.reset();
            }
        });

        ui.horizontal_wrapped(|ui| {
            ui.label(RichText::new("Lienzo principal:").small().weak());
            ui.selectable_value(&mut self.view_mode, ViewMode::View2D, "Solo 2D");
            ui.selectable_value(&mut self.view_mode, ViewMode::Split2D3D, "Dividido 2D+3D");
            ui.selectable_value(&mut self.view_mode, ViewMode::View3D, "Solo 3D");
        });

        let view_w = ui.available_width().max(320.0);
        let view_h = 285.0_f32;
        let (rect, _) = ui.allocate_exact_size(Vec2::new(view_w, view_h), Sense::hover());
        let mol_ref = std::sync::Arc::clone(&self.session.canvas_mol);
        render_3d_viewport(ui, &mol_ref, &mut self.viewer_3d, rect, self.dark_mode);

        ui.add_space(6.0);
        if !mol_ref.atoms.is_empty() {
            let analysis = analyze_molecule(&mol_ref);
            let hbd = mol_ref
                .atoms
                .iter()
                .enumerate()
                .filter(|(idx, a)| {
                    matches!(a.symbol.as_str(), "N" | "O")
                        && (a.implicit_h > 0
                            || mol_ref.adj[*idx]
                                .iter()
                                .any(|&(nb, _)| mol_ref.atoms[nb].symbol == "H"))
                })
                .count();
            let hba = mol_ref
                .atoms
                .iter()
                .filter(|a| matches!(a.symbol.as_str(), "N" | "O"))
                .count();
            let arom_atoms = mol_ref.atoms.iter().filter(|a| a.aromatic).count();
            let coords = compute_3d_coords(&mol_ref);
            let (mut min_c, mut max_c) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
            for p in &coords {
                for k in 0..3 {
                    min_c[k] = min_c[k].min(p[k]);
                    max_c[k] = max_c[k].max(p[k]);
                }
            }
            let span = if coords.is_empty() {
                [0.0, 0.0, 0.0]
            } else {
                [
                    (max_c[0] - min_c[0]).max(0.0),
                    (max_c[1] - min_c[1]).max(0.0),
                    (max_c[2] - min_c[2]).max(0.0),
                ]
            };

            egui::Frame::group(ui.style())
                .corner_radius(egui::CornerRadius::same(6))
                .inner_margin(egui::Margin::symmetric(8, 6))
                .show(ui, |ui| {
                    ui.strong("Perfil 3D y Rasgos Farmacofóricos");
                    ui.label(format!(
                        "Fórmula: {} | MW: {:.2} Da | Átomos pesados: {} | Anillos: {} | Quirales: {}",
                        analysis.formula_with_charge,
                        analysis.molecular_weight,
                        analysis.heavy_atoms,
                        analysis.ring_count,
                        analysis.chiral_centers
                    ));
                    ui.label(format!(
                        "Donadores H (HBD): {hbd} | Aceptores H (HBA): {hba} | Átomos aromáticos: {arom_atoms} | Carga neta: {:+}",
                        analysis.net_charge
                    ));
                    ui.small(format!(
                        "Extensión 3D (caja Å): {:.1} × {:.1} × {:.1} Å — Arrastra para rotar, rueda para zoom, doble clic para centrar.",
                        span[0], span[1], span[2]
                    ));
                });

            ui.add_space(4.0);
            ui.horizontal_wrapped(|ui| {
                if ui.button("Exportar SDF 3D").clicked() {
                    self.export_structure_3d("sdf3d");
                }
                if ui.button("Exportar PDB 3D").clicked() {
                    self.export_structure_3d("pdb");
                }
                if ui.button("Enviar a DockCraft").clicked() {
                    self.send_to_dockcraft(false);
                }
                if ui.button("Abrir en ChimeraX").clicked() {
                    self.open_in_chimerax();
                }
            });
        }
    }

    fn ui_tutorials_miniapp(&mut self, ui: &mut egui::Ui) {
        let tracks = built_in_tutorials();
        if tracks.is_empty() {
            return;
        }
        if self.tutorial_state.active_track_idx >= tracks.len() {
            self.tutorial_state.active_track_idx = 0;
        }
        let t_idx = self.tutorial_state.active_track_idx;

        ui.strong("Elige un tutorial interactivo:");
        ui.horizontal_wrapped(|ui| {
            for (idx, tr) in tracks.iter().enumerate() {
                let done = self.tutorial_state.completed_in_track(idx, tr.steps.len());
                let lbl = format!("{} ({}/{})", tr.title, done, tr.steps.len());
                if ui
                    .selectable_label(idx == t_idx, lbl)
                    .on_hover_text(format!("Tutorial [{}]: {}", tr.id, tr.summary))
                    .clicked()
                {
                    self.tutorial_state.active_track_idx = idx;
                    self.tutorial_state.active_step_idx = 0;
                }
            }
        });

        let track = &tracks[self.tutorial_state.active_track_idx];
        let total_steps = track.steps.len();
        if self.tutorial_state.active_step_idx >= total_steps {
            self.tutorial_state.active_step_idx = 0;
        }
        let s_idx = self.tutorial_state.active_step_idx;
        let done_count = self
            .tutorial_state
            .completed_in_track(self.tutorial_state.active_track_idx, total_steps);
        let progress = (done_count as f32) / (total_steps.max(1) as f32);

        ui.separator();
        ui.horizontal(|ui| {
            ui.strong(track.title);
            ui.small(format!("(~{} min)", track.duration_min));
        });
        ui.small(track.summary);
        ui.add(
            egui::ProgressBar::new(progress)
                .text(format!("Progreso: {done_count} / {total_steps} pasos")),
        );

        ui.add_space(6.0);
        ui.horizontal_wrapped(|ui| {
            for (idx, st) in track.steps.iter().enumerate() {
                let is_done = self
                    .tutorial_state
                    .completed_steps
                    .contains(&(self.tutorial_state.active_track_idx, idx));
                let tag = if is_done { "✓" } else { "•" };
                if ui
                    .selectable_label(
                        idx == s_idx,
                        format!("{tag} Paso {}: {}", idx + 1, st.title),
                    )
                    .clicked()
                {
                    self.tutorial_state.active_step_idx = idx;
                }
            }
        });

        let step = track.steps[s_idx].clone();
        ui.add_space(6.0);
        egui::Frame::group(ui.style())
            .corner_radius(egui::CornerRadius::same(8))
            .inner_margin(egui::Margin::symmetric(12, 10))
            .show(ui, |ui| {
                ui.heading(format!(
                    "Paso {} de {}: {}",
                    s_idx + 1,
                    total_steps,
                    step.title
                ));
                ui.add_space(4.0);
                ui.label(RichText::new("Objetivo:").strong());
                ui.label(step.objective);
                ui.add_space(4.0);
                ui.label(RichText::new("Qué observar en pantalla:").strong());
                ui.label(step.what_to_look_at);
                ui.add_space(8.0);
                ui.horizontal_wrapped(|ui| {
                    let run_btn = egui::Button::new(
                        RichText::new(step.action_label)
                            .strong()
                            .color(Color32::WHITE),
                    )
                    .fill(Color32::from_rgb(4, 120, 87));
                    if ui.add(run_btn).clicked() {
                        self.run_tutorial_action(
                            self.tutorial_state.active_track_idx,
                            s_idx,
                            step.action.clone(),
                        );
                    }
                    if let Some(app_id) = step.related_app
                        && ui
                            .button(format!("Abrir {}", mini_app_meta(app_id).short_label))
                            .clicked()
                    {
                        self.workspace.open(app_id);
                    }
                });
            });

        ui.add_space(6.0);
        ui.horizontal(|ui| {
            if ui
                .add_enabled(s_idx > 0, egui::Button::new("← Paso anterior"))
                .clicked()
            {
                self.tutorial_state.active_step_idx = s_idx.saturating_sub(1);
            }
            if ui
                .add_enabled(
                    s_idx + 1 < total_steps,
                    egui::Button::new("Siguiente paso →"),
                )
                .clicked()
            {
                self.tutorial_state.active_step_idx = (s_idx + 1).min(total_steps - 1);
            }
        });
        ui.add_space(4.0);
        let fb_col = if self.dark_mode {
            Color32::from_rgb(52, 211, 153)
        } else {
            Color32::from_rgb(6, 120, 82)
        };
        ui.colored_label(fb_col, &self.tutorial_state.last_feedback);
    }

    fn render_floating_mini_apps(&mut self, ctx: &egui::Context) {
        let open_list: Vec<MiniAppId> = self.workspace.open_apps.iter().copied().collect();
        let screen = ctx.content_rect();

        for id in open_list {
            let meta = mini_app_meta(id);
            let mut is_open = true;
            let is_expanded = self.workspace.expanded_apps.contains(&id);
            let mut toggle_expand = false;

            let default_x = (screen.left() + meta.default_offset[0]).clamp(
                screen.left() + 24.0,
                (screen.right() - 420.0).max(screen.left() + 24.0),
            );
            let default_y = (screen.top() + meta.default_offset[1]).clamp(
                screen.top() + 86.0,
                (screen.bottom() - 340.0).max(screen.top() + 86.0),
            );

            let mut win = egui::Window::new(meta.window_title)
                .id(egui::Id::unique(("chemcraft_miniapp_win", id)))
                .open(&mut is_open)
                .movable(true)
                .resizable(true)
                .collapsible(true)
                .min_size(egui::vec2(380.0, 320.0))
                .default_pos(egui::pos2(default_x, default_y))
                .default_size(egui::vec2(meta.default_size[0], meta.default_size[1]));

            if is_expanded {
                win = win.fixed_size(egui::vec2(
                    (screen.width() * 0.84).max(640.0),
                    (screen.height() * 0.80).max(520.0),
                ));
            }

            win.show(ctx, |ui| {
                if render_mini_app_header(ui, meta, is_expanded) {
                    toggle_expand = true;
                }
                egui::ScrollArea::vertical()
                    .id_salt(("chemcraft_miniapp_scroll", id))
                    .show(ui, |ui| match id {
                        MiniAppId::Tutorials => self.ui_tutorials_miniapp(ui),
                        MiniAppId::Viewer3dConformers => self.ui_viewer3d_miniapp(ui),
                        MiniAppId::LibrarySarProtac => {
                            let mut sub = self.workspace.sub_tab(id);
                            ui.horizontal_wrapped(|ui| {
                                ui.selectable_value(
                                    &mut sub,
                                    CompoundSubTab::Primary,
                                    "Explorador SAR y Pares MMP",
                                );
                                ui.selectable_value(
                                    &mut sub,
                                    CompoundSubTab::Secondary,
                                    "Bibliotecas (In Silico)",
                                );
                                ui.selectable_value(
                                    &mut sub,
                                    CompoundSubTab::Tertiary,
                                    "Enumerador R-Group / PROTAC",
                                );
                            });
                            self.workspace.set_sub_tab(id, sub);
                            ui.separator();
                            match sub {
                                CompoundSubTab::Primary => {
                                    self.ui_sar_explorer_contents(ui, "miniapp");
                                }
                                CompoundSubTab::Secondary => {
                                    let exp_dir = self.export_dir.clone();
                                    let act = render_library_contents(
                                        ctx,
                                        ui,
                                        &mut self.library_win,
                                        &exp_dir,
                                    );
                                    self.handle_library_action(act);
                                }
                                CompoundSubTab::Tertiary => {
                                    self.ui_enumerator_contents(ui);
                                }
                            }
                        }
                        MiniAppId::SynthesisRetro => {
                            let mut sub = self.workspace.sub_tab(id);
                            ui.horizontal_wrapped(|ui| {
                                ui.selectable_value(
                                    &mut sub,
                                    CompoundSubTab::Primary,
                                    "Retrosíntesis Retro* (Buy-Before-Make)",
                                );
                                ui.selectable_value(
                                    &mut sub,
                                    CompoundSubTab::Secondary,
                                    "Reacción y Estequiometría (.rxn)",
                                );
                            });
                            self.workspace.set_sub_tab(id, sub);
                            ui.separator();
                            match sub {
                                CompoundSubTab::Primary => {
                                    let mol_ref = std::sync::Arc::clone(&self.session.canvas_mol);
                                    let exp_dir = self.export_dir.clone();
                                    let act = render_synth_contents(
                                        ui,
                                        &mut self.synth_module,
                                        &mol_ref,
                                        &exp_dir,
                                        self.dark_mode,
                                    );
                                    self.handle_synth_action(act);
                                }
                                CompoundSubTab::Secondary | CompoundSubTab::Tertiary => {
                                    render_reaction_tab(
                                        ui,
                                        &mut self.reaction_state,
                                        &self.session.canvas_mol,
                                        &self.export_dir,
                                    );
                                }
                            }
                        }
                        MiniAppId::PropertiesAdmetCraftBio => {
                            let mut sub = self.workspace.sub_tab(id);
                            ui.horizontal_wrapped(|ui| {
                                ui.selectable_value(
                                    &mut sub,
                                    CompoundSubTab::Primary,
                                    "Análisis Molecular y ESI-MS",
                                );
                                ui.selectable_value(
                                    &mut sub,
                                    CompoundSubTab::Secondary,
                                    "Espectros 1H / 13C NMR",
                                );
                                ui.selectable_value(
                                    &mut sub,
                                    CompoundSubTab::Tertiary,
                                    "ADMET, Bioisósteros y CraftBio",
                                );
                            });
                            self.workspace.set_sub_tab(id, sub);
                            ui.separator();
                            match sub {
                                CompoundSubTab::Primary => {
                                    render_analysis_tab(ui, &self.session.canvas_mol);
                                }
                                CompoundSubTab::Secondary => {
                                    render_nmr_tab(ui, &self.session.canvas_mol, self.dark_mode);
                                }
                                CompoundSubTab::Tertiary => {
                                    self.ui_admet_and_craftbio_contents(ui, "miniapp");
                                }
                            }
                        }
                    });
            });

            if !is_open {
                self.workspace.open_apps.remove(&id);
                self.workspace.expanded_apps.remove(&id);
            } else if toggle_expand {
                self.workspace.toggle_expand(id);
            }
        }
    }
}

impl eframe::App for ChemCraftApp {
    fn ui(&mut self, ui: &mut egui::Ui, frame: &mut eframe::Frame) {
        let ctx = ui.ctx().clone();
        self.process_dialog_answers();
        self.sync_theme_and_prefs(&ctx);
        // Honor any external toggle of legacy module `open` flags before handling UI clicks
        if self.library_win.open {
            self.workspace.open(MiniAppId::LibrarySarProtac);
        }
        if self.synth_module.open {
            self.workspace.open(MiniAppId::SynthesisRetro);
        }
        self.handle_shortcuts(&ctx);
        self.render_palette(&ctx);

        // Keep polling any background headless DockCraft screening job even if window is closed
        let bg_lib_act = poll_library_docking_job(&ctx, &mut self.library_win);
        self.handle_library_action(bg_lib_act);

        // 1. Craft-Standard Top Menu Bar
        egui::Panel::top("craft_menu_bar").show(ui, |ui| {
            egui::MenuBar::new().ui(ui, |ui| {
                ui.label(
                    RichText::new("ChemCraft")
                        .strong()
                        .color(AppAccent::CHEM_AMBER.color32()),
                );
                ui.label(
                    RichText::new(format!("v{}", env!("CARGO_PKG_VERSION")))
                        .small()
                        .weak(),
                );
                ui.separator();

                ui.menu_button("File", |ui| {
                    if ui.button("New Workspace (Ctrl+N)").clicked() {
                        self.current_project_path = None;
                        self.dirty = false;
                        self.run_cmd(Command::ClearCanvas);
                        ui.close();
                    }
                    if ui.button("Open... (Ctrl+O)").clicked() {
                        self.open_file_dialog();
                        ui.close();
                    }
                    ui.separator();
                    if ui.button("Save (Ctrl+S)").clicked() {
                        self.save_project();
                        ui.close();
                    }
                    if ui.button("Save As... (Ctrl+Shift+S)").clicked() {
                        self.save_project_as();
                        ui.close();
                    }
                    ui.separator();
                    ui.menu_button("Export Structure / Data", |ui| {
                        if ui.button("ChemDraw (.cdxml)...").clicked() {
                            self.export_cdxml_file();
                            ui.close();
                        }
                        if ui.button("SMILES (.smi)...").clicked() {
                            self.export_smiles_file();
                            ui.close();
                        }
                        ui.separator();
                        if ui.button("2D Vector SVG (.svg)...").clicked() {
                            self.export_canvas("svg");
                            ui.close();
                        }
                        if ui.button("2D Vector PDF (.pdf)...").clicked() {
                            self.export_canvas("pdf");
                            ui.close();
                        }
                        if ui.button("2D Vector EPS (.eps)...").clicked() {
                            self.export_canvas("eps");
                            ui.close();
                        }
                        if ui.button("High-DPI PNG (.png)...").clicked() {
                            self.export_canvas("png");
                            ui.close();
                        }
                        ui.separator();
                        if ui.button("3D Conformer SDF (.sdf)...").clicked() {
                            self.export_structure_3d("sdf3d");
                            ui.close();
                        }
                        if ui.button("3D Structure PDB (.pdb)...").clicked() {
                            self.export_structure_3d("pdb");
                            ui.close();
                        }
                        if !self.session.sar_table.rows.is_empty() {
                            ui.separator();
                            if ui.button("SAR Table CSV (.csv)...").clicked() {
                                self.export_sar_csv_file();
                                ui.close();
                            }
                        }
                    });
                    ui.separator();
                    if ui.button("Set Default Export Folder...").clicked() {
                        self.dialogs.ask(
                            choose_folder(
                                "Select Default Export Folder",
                                Some(Path::new(&self.export_dir)),
                            ),
                            ChemDialogAction::ChooseExportDir,
                        );
                        ui.close();
                    }
                });

                ui.menu_button("Edit", |ui| {
                    if ui
                        .add_enabled(self.session.can_undo(), egui::Button::new("Undo (Ctrl+Z)"))
                        .clicked()
                    {
                        self.run_cmd(Command::Undo);
                        ui.close();
                    }
                    if ui
                        .add_enabled(self.session.can_redo(), egui::Button::new("Redo (Ctrl+Y)"))
                        .clicked()
                    {
                        self.run_cmd(Command::Redo);
                        ui.close();
                    }
                    ui.separator();
                    if ui.button("Clean 2D Structure (Ctrl+L)").clicked() {
                        self.run_cmd(Command::CleanStructure2D);
                        ui.close();
                    }
                    if ui.button("Flip Horizontal").clicked() {
                        self.run_cmd(Command::FlipCanvas2D { horizontal: true });
                        ui.close();
                    }
                    if ui.button("Flip Vertical").clicked() {
                        self.run_cmd(Command::FlipCanvas2D { horizontal: false });
                        ui.close();
                    }
                    if ui.button("Rotate 30°").clicked() {
                        self.run_cmd(Command::RotateCanvas2D { degrees: 30.0 });
                        ui.close();
                    }
                    ui.separator();
                    if ui.button("Bioisostere: -COOH -> Tetrazole").clicked() {
                        self.run_cmd(Command::MutateBioisostere {
                            atom_index: None,
                            group: "tetrazole".to_string(),
                        });
                        ui.close();
                    }
                    if ui.button("Copy Canonical SMILES").clicked() {
                        ui.ctx()
                            .copy_text(self.session.canvas_mol.canonical_smiles());
                        self.status_msg = "Copied canonical SMILES to clipboard".to_string();
                        ui.close();
                    }
                    ui.separator();
                    if ui.button("Clear Canvas").clicked() {
                        self.run_cmd(Command::ClearCanvas);
                        ui.close();
                    }
                });

                ui.menu_button("View", |ui| {
                    ui.selectable_value(&mut self.view_mode, ViewMode::View2D, "2D Vector");
                    ui.selectable_value(&mut self.view_mode, ViewMode::View3D, "3D Interactive");
                    ui.selectable_value(&mut self.view_mode, ViewMode::Split2D3D, "Split 2D + 3D");
                    ui.separator();
                    ui.menu_button("2D Depiction Style", |ui| {
                        ui.selectable_value(
                            &mut self.depiction_preset,
                            DepictionPreset::ModernColor,
                            "Modern Color",
                        );
                        ui.selectable_value(
                            &mut self.depiction_preset,
                            DepictionPreset::Acs1996Bw,
                            "ACS Document 1996 (B&W)",
                        );
                        ui.selectable_value(
                            &mut self.depiction_preset,
                            DepictionPreset::Acs1996Color,
                            "ACS Document 1996 (Color)",
                        );
                    });
                    ui.checkbox(&mut self.show_terminal_methyls, "Show Terminal CH3");
                    ui.separator();
                    ui.checkbox(
                        &mut self.viewer_3d.show_hydrogens,
                        "3D: Show Hydrogens (+H)",
                    );
                    ui.checkbox(&mut self.viewer_3d.show_labels, "3D: Atom Labels");
                    ui.checkbox(&mut self.viewer_3d.show_ribbon, "3D: Backbone Ribbon");
                    ui.checkbox(&mut self.viewer_3d.auto_spin, "3D: Auto-Spin");
                    if ui.button("3D: Rebuild Conformer").clicked() {
                        let mut mol = (*self.session.canvas_mol).clone();
                        generate_3d_coords(&mut mol);
                        self.session.canvas_mol = std::sync::Arc::new(mol);
                        self.status_msg = "Regenerated 3D conformer coordinates".to_string();
                        ui.close();
                    }
                    if ui.button("3D: Reset Camera").clicked() {
                        self.viewer_3d.reset();
                        ui.close();
                    }
                    ui.separator();
                    ui.checkbox(
                        &mut self.workspace.show_classic_inspector,
                        "Inspector lateral clásico",
                    );
                    ui.checkbox(&mut self.sar_panel_open, "Panel SAR inferior acoplado");
                    ui.checkbox(&mut self.dark_mode, "Dark Theme");
                    if ui
                        .checkbox(
                            &mut self.use_system_titlebar,
                            "Use System Titlebar (restart)",
                        )
                        .changed()
                    {
                        self.status_msg =
                            "Window frame preference saved (restart to apply).".to_string();
                    }
                });

                ui.menu_button("Tools", |ui| {
                    if ui.button("Screening Libraries (In Silico)...").clicked() {
                        self.workspace.open(MiniAppId::LibrarySarProtac);
                        self.workspace
                            .set_sub_tab(MiniAppId::LibrarySarProtac, CompoundSubTab::Secondary);
                        ui.close();
                    }
                    if ui.button("R-Group & PROTAC Enumerator...").clicked() {
                        self.workspace.open(MiniAppId::LibrarySarProtac);
                        self.workspace
                            .set_sub_tab(MiniAppId::LibrarySarProtac, CompoundSubTab::Tertiary);
                        ui.close();
                    }
                    if ui
                        .button("Synth Module (Retrosynthesis & Buy-Before-Make)...")
                        .clicked()
                    {
                        self.workspace.open(MiniAppId::SynthesisRetro);
                        self.workspace
                            .set_sub_tab(MiniAppId::SynthesisRetro, CompoundSubTab::Primary);
                        ui.close();
                    }
                    ui.separator();
                    if ui.button("Add Canvas Molecule to SAR Table").clicked() {
                        self.run_cmd(Command::AddCanvasToSar { pic50: None });
                        ui.close();
                    }
                    if ui.button("Filter SAR by Canvas Substructure").clicked() {
                        self.run_cmd(Command::FilterSarByCanvasSubstructure);
                        ui.close();
                    }
                    ui.separator();
                    if ui.button("Command Palette (Ctrl+K)").clicked() {
                        self.palette_open = true;
                        ui.close();
                    }
                });

                ui.menu_button("Mini-Apps", |ui| {
                    for meta in ALL_MINI_APPS {
                        let open = self.workspace.open_apps.contains(&meta.id);
                        let lbl = if open {
                            format!("✓ {}", meta.short_label)
                        } else {
                            meta.short_label.to_string()
                        };
                        if ui.button(lbl).on_hover_text(meta.what_it_does).clicked() {
                            self.workspace.toggle(meta.id);
                            ui.close();
                        }
                    }
                    ui.separator();
                    if ui
                        .checkbox(
                            &mut self.workspace.show_classic_inspector,
                            "Mostrar inspector lateral derecho clásico",
                        )
                        .changed()
                    {
                        ui.close();
                    }
                    if ui
                        .checkbox(
                            &mut self.sar_panel_open,
                            "Mostrar panel SAR inferior acoplado",
                        )
                        .changed()
                    {
                        ui.close();
                    }
                    if ui
                        .checkbox(
                            &mut self.workspace.left_panel_floating,
                            "Paleta de dibujo 2D como ventana flotante móvil",
                        )
                        .changed()
                    {
                        self.workspace.show_left_panel = true;
                        ui.close();
                    }
                });

                ui.menu_button("Apps", |ui| {
                    if ui.button("Abrir en AdmetCraft (ADMET & PK)").clicked() {
                        self.open_in_admetcraft();
                        ui.close();
                    }
                    if ui.button("Enviar a DockCraft (3D Docking)").clicked() {
                        self.send_to_dockcraft(false);
                        ui.close();
                    }
                    if ui.button("Visualizar en ChimeraX (3D Avanzado)").clicked() {
                        self.open_in_chimerax();
                        ui.close();
                    }
                });

                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    if ui
                        .small_button("Ctrl+K")
                        .on_hover_text("Command Palette")
                        .clicked()
                    {
                        self.palette_open = true;
                    }
                    let mut mode = ThemeMode::from_dark(self.dark_mode);
                    if theme_toggle_button(ui, &mut mode) {
                        self.dark_mode = mode.is_dark();
                    }
                });
            });
        });

        // 1b. Quick Action Toolbar
        egui::Panel::top("craft_quick_toolbar").show(ui, |ui| {
            ui.horizontal(|ui| {
                if ui
                    .button("Open")
                    .on_hover_text("Open project or structure (Ctrl+O)")
                    .clicked()
                {
                    self.open_file_dialog();
                }
                if ui
                    .button("Save")
                    .on_hover_text("Save .chemcraft project (Ctrl+S)")
                    .clicked()
                {
                    self.save_project();
                }
                ui.menu_button("Export ▾", |ui| {
                    if ui.button("ChemDraw (.cdxml)...").clicked() {
                        self.export_cdxml_file();
                        ui.close();
                    }
                    if ui.button("SVG (.svg)...").clicked() {
                        self.export_canvas("svg");
                        ui.close();
                    }
                    if ui.button("PDF (.pdf)...").clicked() {
                        self.export_canvas("pdf");
                        ui.close();
                    }
                    if ui.button("EPS (.eps)...").clicked() {
                        self.export_canvas("eps");
                        ui.close();
                    }
                    if ui.button("PNG (.png)...").clicked() {
                        self.export_canvas("png");
                        ui.close();
                    }
                    ui.separator();
                    if ui.button("3D SDF (.sdf)...").clicked() {
                        self.export_structure_3d("sdf3d");
                        ui.close();
                    }
                    if ui.button("3D PDB (.pdb)...").clicked() {
                        self.export_structure_3d("pdb");
                        ui.close();
                    }
                    if ui.button("SMILES (.smi)...").clicked() {
                        self.export_smiles_file();
                        ui.close();
                    }
                });
                ui.separator();
                if ui
                    .add_enabled(self.session.can_undo(), egui::Button::new("Undo"))
                    .clicked()
                {
                    self.run_cmd(Command::Undo);
                }
                if ui
                    .add_enabled(self.session.can_redo(), egui::Button::new("Redo"))
                    .clicked()
                {
                    self.run_cmd(Command::Redo);
                }
                if ui.button("Clean 2D").clicked() {
                    self.run_cmd(Command::CleanStructure2D);
                }
                ui.separator();
                ui.selectable_value(&mut self.view_mode, ViewMode::View2D, "2D");
                ui.selectable_value(&mut self.view_mode, ViewMode::View3D, "3D");
                ui.selectable_value(&mut self.view_mode, ViewMode::Split2D3D, "Split");

                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    if ui
                        .button("ChimeraX")
                        .on_hover_text(if self.chimerax_available {
                            "Exporta .pdb/.sdf + .cxc y abre en UCSF ChimeraX"
                        } else {
                            "Exporta bundle .pdb/.sdf + .cxc para UCSF ChimeraX"
                        })
                        .clicked()
                    {
                        self.open_in_chimerax();
                    }
                    if ui
                        .button("DockCraft")
                        .on_hover_text(if self.dockcraft_available {
                            "Exporta SDF 3D + SMILES y abre DockCraft"
                        } else {
                            "Exporta bundle SDF 3D + SMILES listo para DockCraft"
                        })
                        .clicked()
                    {
                        self.send_to_dockcraft(false);
                    }
                    if ui
                        .button(RichText::new("AdmetCraft").strong())
                        .on_hover_text(if self.admet_gui_available {
                            "Envía la molécula/lote actual y abre AdmetCraft para perfil ADMET completo"
                        } else {
                            "Guarda SMILES handoff para abrir en AdmetCraft"
                        })
                        .clicked()
                    {
                        self.open_in_admetcraft();
                    }
                    ui.label(RichText::new("Abrir en:").small().weak());
                });
            });
        });

        // 1c. Floating Mini-Applications Launcher Bar (1-click movable, resizable & expandable windows)
        egui::Panel::top("chemcraft_miniapp_bar").show(ui, |ui| {
            ui.horizontal_wrapped(|ui| {
                ui.strong("Mini-Apps:");
                for meta in ALL_MINI_APPS {
                    let is_open = self.workspace.open_apps.contains(&meta.id);
                    let [r, g, b] = meta.accent_rgb;
                    let text = if is_open {
                        RichText::new(meta.short_label)
                            .strong()
                            .color(Color32::WHITE)
                    } else {
                        RichText::new(meta.short_label)
                    };
                    let btn = egui::Button::new(text).fill(if is_open {
                        Color32::from_rgb(r, g, b)
                    } else {
                        ui.visuals().widgets.inactive.weak_bg_fill
                    });
                    if ui
                        .add(btn)
                        .on_hover_text(format!(
                            "{}\n\nCuándo usarla: {}\nFlujo: {}",
                            meta.what_it_does, meta.when_to_use, meta.io_summary
                        ))
                        .clicked()
                    {
                        self.workspace.toggle(meta.id);
                    }
                }
                ui.separator();
                if ui
                    .small_button(if self.workspace.left_panel_floating {
                        "📌 Herramientas"
                    } else {
                        "🪟 Flotar Herramientas"
                    })
                    .on_hover_text(
                        "Alterna la paleta de dibujo 2D entre panel lateral acoplado o ventana flotante móvil",
                    )
                    .clicked()
                {
                    self.workspace.left_panel_floating = !self.workspace.left_panel_floating;
                    self.workspace.show_left_panel = true;
                }
                if ui
                    .small_button(if self.workspace.show_classic_inspector {
                        "📋 Inspector: ON"
                    } else {
                        "📋 Inspector clásico"
                    })
                    .on_hover_text(
                        "Muestra/oculta el panel lateral derecho clásico (Analysis / NMR / Reaction / ADMET)",
                    )
                    .clicked()
                {
                    self.workspace.show_classic_inspector = !self.workspace.show_classic_inspector;
                }
            });
        });

        // Bottom Status Bar (with version and active project file)
        status_bar(
            ui,
            "status_bar",
            |ui| {
                ui.label(&self.status_msg);
            },
            |ui| {
                ui.label(
                    RichText::new(format!("ChemCraft v{}", env!("CARGO_PKG_VERSION")))
                        .small()
                        .strong(),
                );
                ui.separator();
                let proj_lbl = self
                    .current_project_path
                    .as_ref()
                    .and_then(|p| p.file_name())
                    .and_then(|s| s.to_str())
                    .unwrap_or("Unsaved");
                let dirty_mark = if self.dirty { "*" } else { "" };
                ui.label(
                    RichText::new(format!("Project: {proj_lbl}{dirty_mark}"))
                        .small()
                        .weak(),
                );
            },
        );

        // 2. Left Panel (or Detachable Movable Window): Vector Drawing Tools, Superatoms & Ring Templates
        if self.workspace.show_left_panel {
            if self.workspace.left_panel_floating {
                let mut open = self.workspace.show_left_panel;
                let screen = ctx.content_rect();
                egui::Window::new("✏️ Herramientas de Dibujo 2D y Plantillas")
                    .open(&mut open)
                    .movable(true)
                    .resizable(true)
                    .collapsible(true)
                    .default_pos(egui::pos2(screen.left() + 18.0, screen.top() + 96.0))
                    .default_size([220.0, 590.0])
                    .min_size([185.0, 280.0])
                    .show(&ctx, |ui| self.ui_left_tools_panel(ui));
                self.workspace.show_left_panel = open;
            } else {
                egui::Panel::left("tools_panel")
                    .resizable(true)
                    .default_size(195.0)
                    .show(ui, |ui| self.ui_left_tools_panel(ui));
            }
        }

        // 3. Optional Classic Right Dockable Panel (hidden by default so central 2D/3D canvas stays spacious)
        if self.workspace.show_classic_inspector {
            egui::Panel::right("admet_inspector")
                .resizable(true)
                .default_size(355.0)
                .min_size(310.0)
                .show(ui, |ui| {
                    ui.horizontal_wrapped(|ui| {
                        ui.selectable_value(
                            &mut self.inspector_tab,
                            InspectorTab::Analysis,
                            "Analysis",
                        );
                        ui.selectable_value(
                            &mut self.inspector_tab,
                            InspectorTab::Nmr,
                            "1H/13C NMR",
                        );
                        ui.selectable_value(
                            &mut self.inspector_tab,
                            InspectorTab::Reaction,
                            "Reaction",
                        );
                        ui.selectable_value(&mut self.inspector_tab, InspectorTab::Admet, "ADMET");
                    });
                    ui.separator();
                    egui::ScrollArea::vertical().show(ui, |ui| match self.inspector_tab {
                        InspectorTab::Analysis => {
                            render_analysis_tab(ui, &self.session.canvas_mol);
                        }
                        InspectorTab::Nmr => {
                            render_nmr_tab(ui, &self.session.canvas_mol, self.dark_mode);
                        }
                        InspectorTab::Reaction => {
                            render_reaction_tab(
                                ui,
                                &mut self.reaction_state,
                                &self.session.canvas_mol,
                                &self.export_dir,
                            );
                        }
                        InspectorTab::Admet => {
                            self.ui_admet_and_craftbio_contents(ui, "classic");
                        }
                    });
                });
        }

        // 4. Central Panel: Spacious 2D/3D Drawing Viewport + Optional Classic Docked Bottom SAR Explorer
        egui::CentralPanel::default().show(ui, |ui| {
            // Single compact structure input bar (SMILES + Name-to-Structure)
            ui.horizontal_wrapped(|ui| {
                ui.label("SMILES:");
                let smi_edit = ui.add(
                    egui::TextEdit::singleline(&mut self.smiles_bar)
                        .desired_width(240.0)
                        .font(egui::TextStyle::Monospace),
                );
                if (smi_edit.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter)))
                    || ui.button("Load").clicked()
                {
                    let smi = self.smiles_bar.clone();
                    self.run_cmd(Command::LoadCanvasSmiles {
                        smiles: smi,
                        name: None,
                    });
                }
                if ui
                    .small_button("Copy")
                    .on_hover_text("Copy canonical SMILES")
                    .clicked()
                {
                    ui.ctx()
                        .copy_text(self.session.canvas_mol.canonical_smiles());
                    self.status_msg = "Copied canonical SMILES to clipboard".to_string();
                }
                ui.separator();
                ui.label("Name:");
                let name_edit = ui.add(
                    egui::TextEdit::singleline(&mut self.name_bar)
                        .desired_width(150.0)
                        .hint_text("aspirin, caffeine..."),
                );
                if (name_edit.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter)))
                    || ui.button("Resolve").clicked()
                {
                    let nm = self.name_bar.clone();
                    if !nm.trim().is_empty() {
                        self.run_cmd(Command::LoadChemicalName { name: nm });
                    }
                }
                ui.separator();
                if self.view_mode != ViewMode::View2D {
                    ui.checkbox(&mut self.viewer_3d.show_hydrogens, "Show H");
                    ui.checkbox(&mut self.viewer_3d.show_labels, "Labels");
                    ui.checkbox(&mut self.viewer_3d.auto_spin, "Spin");
                    if ui.small_button("Rebuild 3D").clicked() {
                        let mut mol = (*self.session.canvas_mol).clone();
                        generate_3d_coords(&mut mol);
                        self.session.canvas_mol = std::sync::Arc::new(mol);
                        self.status_msg = "Regenerated 3D conformer coordinates".to_string();
                    }
                    if ui.small_button("Reset Cam").clicked() {
                        self.viewer_3d.reset();
                    }
                }
            });

            if !self.session.canvas_mol.warnings.is_empty() {
                ui.colored_label(
                    Color32::from_rgb(220, 48, 48),
                    format!("Warnings: {}", self.session.canvas_mol.warnings.join(" | ")),
                );
            }

            // Center-top interactive viewport (2D Vector / 3D Conformer / Split 2D + 3D)
            let canvas_h = if self.sar_panel_open {
                (ui.available_height() * 0.44).max(220.0)
            } else {
                ui.available_height().max(320.0)
            };
            let canvas_w = ui.available_width().max(320.0);
            let (full_rect, _) =
                ui.allocate_exact_size(Vec2::new(canvas_w, canvas_h), Sense::hover());

            match self.view_mode {
                ViewMode::View2D => {
                    self.draw_2d_canvas_at_rect(ui, full_rect);
                }
                ViewMode::View3D => {
                    let mol_ref = std::sync::Arc::clone(&self.session.canvas_mol);
                    render_3d_viewport(
                        ui,
                        &mol_ref,
                        &mut self.viewer_3d,
                        full_rect,
                        self.dark_mode,
                    );
                }
                ViewMode::Split2D3D => {
                    let gap = 6.0_f32;
                    let half_w = ((full_rect.width() - gap) * 0.5).max(120.0);
                    let left_rect = egui::Rect::from_min_size(
                        full_rect.min,
                        Vec2::new(half_w, full_rect.height()),
                    );
                    let right_rect = egui::Rect::from_min_size(
                        Pos2::new(full_rect.min.x + half_w + gap, full_rect.min.y),
                        Vec2::new(half_w, full_rect.height()),
                    );
                    self.draw_2d_canvas_at_rect(ui, left_rect);
                    let mol_ref = std::sync::Arc::clone(&self.session.canvas_mol);
                    render_3d_viewport(
                        ui,
                        &mol_ref,
                        &mut self.viewer_3d,
                        right_rect,
                        self.dark_mode,
                    );
                }
            }

            if self.sar_panel_open {
                ui.separator();
                self.ui_sar_explorer_contents(ui, "docked");
            }
        });

        // 5. Floating Movable, Resizable & Expandable Mini-Applications above the spacious canvas!
        self.render_floating_mini_apps(&ctx);

        self.library_win.open = self
            .workspace
            .open_apps
            .contains(&MiniAppId::LibrarySarProtac);
        self.synth_module.open = self
            .workspace
            .open_apps
            .contains(&MiniAppId::SynthesisRetro);

        self.dialogs.show_pending(Some(frame));
    }

    fn on_exit(&mut self) {
        let _ = self.settings.flush();
    }
}

fn point_seg_dist(p: [f64; 2], a: [f64; 2], b: [f64; 2]) -> f64 {
    let dx = b[0] - a[0];
    let dy = b[1] - a[1];
    let len2 = dx * dx + dy * dy;
    if len2 < 1e-6 {
        return (p[0] - a[0]).hypot(p[1] - a[1]);
    }
    let t = (((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / len2).clamp(0.0, 1.0);
    let px = a[0] + t * dx;
    let py = a[1] + t * dy;
    (p[0] - px).hypot(p[1] - py)
}

fn main() -> eframe::Result<()> {
    let (settings, _load_report) = open_settings::<ChemPrefs>();
    let use_system_titlebar = settings.get().ui.use_system_titlebar;
    let title = format!(
        "ChemCraft v{} — 2D/3D Molecular Design & SAR Studio",
        env!("CARGO_PKG_VERSION")
    );
    let win_cfg = WindowConfig::new(&title, "ChemCraft", "CHEMCRAFT")
        .with_size([1380.0, 860.0], [980.0, 660.0])
        .with_titlebar_pref(titlebar_pref(use_system_titlebar));
    let (viewport, _decision) = build_viewport(&win_cfg);
    let options = eframe::NativeOptions {
        viewport,
        ..Default::default()
    };
    eframe::run_native(
        "ChemCraft",
        options,
        Box::new(move |cc| Ok(Box::new(ChemCraftApp::with_context(&cc.egui_ctx, settings)))),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn app_starts_clean_with_optional_sar_and_runs_bioisostere_and_undo() {
        let mut app = ChemCraftApp::new();
        assert!(!app.sar_panel_open);
        assert!(!app.synth_module.open);
        assert!(app.synth_module.last_plan.is_none());
        app.synth_module.run_for_molecule(&app.session.canvas_mol);
        assert!(
            app.synth_module
                .last_plan
                .as_ref()
                .is_some_and(|p| !p.routes.is_empty()),
            "expected Synth Module to solve default canvas molecule on demand"
        );
        assert!(!app.synth_module.thumb_cache.is_empty());
        assert!(!app.synth_module.open);
        assert!(app.session.sar_table.rows.is_empty());
        // ADMET is strictly on-demand at startup (not auto-computed)
        assert!(app.session.live_admet.is_none());
        assert_eq!(app.inspector_tab, InspectorTab::Analysis);
        assert!(app.session.refresh_live_admet().is_ok());
        assert!(app.session.live_admet.is_some());
        app.run_cmd(Command::LoadDemo50Sar);
        assert_eq!(app.session.sar_table.rows.len(), 50);
        app.run_cmd(Command::MutateBioisostere {
            atom_index: None,
            group: "tetrazole".to_string(),
        });
        assert!(app.session.last_mutation_delta.is_some());
        assert!(app.session.can_undo());
        app.run_cmd(Command::Undo);
        assert!(app.session.can_redo());

        // Topological 2D edits keep cached 3D coords populated and support Undo
        let before_atoms = app.session.canvas_mol.atoms.len();
        app.run_cmd(Command::AddAtomAndBond {
            from_atom: Some(0),
            symbol: "F".to_string(),
            x: 2.0,
            y: 1.0,
            order: BondOrder::Single,
            stereo: BondStereo::None,
        });
        assert_eq!(app.session.canvas_mol.atoms.len(), before_atoms + 1);
        assert!(
            app.session
                .canvas_mol
                .atoms
                .iter()
                .all(|a| a.coords_3d.is_some())
        );
        app.run_cmd(Command::EraseAtom {
            atom_idx: before_atoms,
        });
        assert_eq!(app.session.canvas_mol.atoms.len(), before_atoms);
        app.run_cmd(Command::Undo);
        assert_eq!(app.session.canvas_mol.atoms.len(), before_atoms + 1);
    }
}
