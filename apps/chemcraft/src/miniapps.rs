#![forbid(unsafe_code)]

//! Floating, Movable, Resizable & Expandable Mini-Application Window Manager for ChemCraft.
//!
//! Organizes every ChemCraft capability (2D/3D conformers, SAR / Libraries / R-Group & PROTAC
//! enumeration, Retro* synthesis & stoichiometry, Molecular Analysis / 1H-13C NMR / ADMET /
//! CraftBio bridges, and Interactive Tutorials) as self-contained floating windows with a crisp
//! Spanish educational banner (*Qué hace*, *Cuándo usarla*, *Flujo I/O*), drag-to-move, free
//! resizing, and 1-click Expand/Restore.

use std::collections::{BTreeMap, BTreeSet};

use eframe::egui;
use mol_synth::{
    E3RecruiterPreset, EnumerateMode, EnumeratedLibrary, enumerate_protac_library,
    enumerate_rgroup_library,
};

/// Every standalone floating mini-application inside ChemCraft.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum MiniAppId {
    Tutorials,
    Viewer3dConformers,
    LibrarySarProtac,
    SynthesisRetro,
    PropertiesAdmetCraftBio,
}

/// Descriptive metadata and initial geometry for one floating mini-application.
#[derive(Debug, Clone, Copy)]
pub struct MiniAppMeta {
    pub id: MiniAppId,
    pub short_label: &'static str,
    pub window_title: &'static str,
    pub what_it_does: &'static str,
    pub when_to_use: &'static str,
    pub io_summary: &'static str,
    pub default_size: [f32; 2],
    pub default_offset: [f32; 2],
    pub accent_rgb: [u8; 3],
}

pub const ALL_MINI_APPS: &[MiniAppMeta] = &[
    MiniAppMeta {
        id: MiniAppId::Tutorials,
        short_label: "📚 Tutoriales",
        window_title: "📚 Tutoriales Interactivos · Aprende ChemCraft Paso a Paso",
        what_it_does: "Guías prácticas paso a paso con botones que ejecutan flujos reales sobre el lienzo 2D/3D, SAR, enumeración PROTAC y retrosíntesis.",
        when_to_use: "Ideal al empezar o para dominar el editor 2D, confórmeros 3D, enumeración R-Group/PROTAC, Retro* y puentes CraftBio.",
        io_summary: "Paso interactivo → muta el lienzo real, lanza análisis y abre la mini-app correspondiente.",
        default_size: [540.0, 560.0],
        default_offset: [220.0, 96.0],
        accent_rgb: [16, 185, 129],
    },
    MiniAppMeta {
        id: MiniAppId::Viewer3dConformers,
        short_label: "🧊 Visor 3D y Confórmeros",
        window_title: "🧊 Visor 3D Interactivo, Confórmeros y Perfil Farmacofórico",
        what_it_does: "Muestra y rota la molécula activa en 3D (+H, etiquetas, cinta de esqueleto para PDB), regenera confórmeros EtKDG y resume rasgos farmacofóricos.",
        when_to_use: "Para inspeccionar geometría espacial, centros quirales, donadores/aceptores H y preparar ligandos 3D para DockCraft o ChimeraX.",
        io_summary: "Molécula 2D/SMILES → Confórmero 3D interactivo + resumen farmacofórico → SDF 3D / PDB / DockCraft / ChimeraX.",
        default_size: [560.0, 600.0],
        default_offset: [264.0, 102.0],
        accent_rgb: [6, 182, 212],
    },
    MiniAppMeta {
        id: MiniAppId::LibrarySarProtac,
        short_label: "📊 Biblioteca · SAR · PROTAC",
        window_title: "📊 Biblioteca Química, SAR, Pares MMP y Enumerador R-Group / PROTAC",
        what_it_does: "Reúne el Explorador SAR (Tabla, Scatter, Activity Cliffs SALI, Pares MMP y Filtros), cribado de bibliotecas in silico y enumerador combinatorio R-Group / PROTAC.",
        when_to_use: "Para analizar series químicas, detectar acantilados SAR, preparar librerías 3D o ensamblar análogos R-Group y degradadores PROTAC (CRBN/VHL/IAP).",
        io_summary: "SDF/CSV/DWAR o scaffold [*] → Tabla SAR + Cliffs/MMP + biblioteca enumerada → Lienzo / DockCraft / CSV.",
        default_size: [740.0, 620.0],
        default_offset: [308.0, 108.0],
        accent_rgb: [99, 102, 241],
    },
    MiniAppMeta {
        id: MiniAppId::SynthesisRetro,
        short_label: "⚗️ Síntesis y Retro*",
        window_title: "⚗️ Síntesis, Retrosíntesis Retro* (Buy-Before-Make) y Reacciones (.rxn)",
        what_it_does: "Planifica rutas retrosintéticas multipaso Retro* contra el catálogo comercial Buy-Before-Make y calcula estequiometría, economía atómica y rendimiento (.rxn).",
        when_to_use: "Cuando quieras saber si un hit es comprable, diseñar su síntesis con protocolos ELN o balancear reacciones SMILES / MDL .rxn.",
        io_summary: "Molécula del lienzo / Reaction SMILES → Rutas Retro*, bloques en stock y tabla estequiométrica → ELN .md / JSON / .rxn.",
        default_size: [780.0, 640.0],
        default_offset: [348.0, 114.0],
        accent_rgb: [217, 119, 6],
    },
    MiniAppMeta {
        id: MiniAppId::PropertiesAdmetCraftBio,
        short_label: "🔬 Propiedades · NMR · CraftBio",
        window_title: "🔬 Propiedades, 1H/13C NMR, ADMET, Bioisósteros y Puentes CraftBio",
        what_it_does: "Calcula masa exacta, análisis elemental, aductos ESI-MS, espectros 1H/13C NMR simulados, perfil ADMET bajo demanda (radar + alertas), bioisósteros y puentes CraftBio.",
        when_to_use: "Para caracterizar analíticamente la estructura activa, comparar deltas bioisostéricos y transferir compuestos a AdmetCraft, DockCraft o ChimeraX.",
        io_summary: "Molécula activa → Masa/ESI-MS, 1H/13C NMR, Radar ADMET → CSV / PDF / BOILED-Egg / AdmetCraft / DockCraft / ChimeraX.",
        default_size: [600.0, 640.0],
        default_offset: [428.0, 126.0],
        accent_rgb: [124, 58, 237],
    },
];

#[must_use]
pub fn mini_app_meta(id: MiniAppId) -> &'static MiniAppMeta {
    ALL_MINI_APPS
        .iter()
        .find(|m| m.id == id)
        .unwrap_or(&ALL_MINI_APPS[0])
}

/// Active sub-tab selection for compound mini-apps that unite related tools cleanly.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CompoundSubTab {
    Primary,
    Secondary,
    Tertiary,
}

/// Manages open/expanded floating mini-app windows and detachable workspace panels.
#[derive(Debug, Clone)]
pub struct WorkspaceWindows {
    pub open_apps: BTreeSet<MiniAppId>,
    pub expanded_apps: BTreeSet<MiniAppId>,
    pub sub_tabs: BTreeMap<MiniAppId, CompoundSubTab>,
    /// Keep classic right inspector sidebar hidden by default so the central 2D/3D canvas stays spacious.
    pub show_classic_inspector: bool,
    /// Allow detaching the left Drawing Tools palette as a movable/resizable floating window.
    pub left_panel_floating: bool,
    pub show_left_panel: bool,
}

impl Default for WorkspaceWindows {
    fn default() -> Self {
        let mut open_apps = BTreeSet::new();
        open_apps.insert(MiniAppId::Tutorials);
        Self {
            open_apps,
            expanded_apps: BTreeSet::new(),
            sub_tabs: BTreeMap::new(),
            show_classic_inspector: false,
            left_panel_floating: false,
            show_left_panel: true,
        }
    }
}

impl WorkspaceWindows {
    pub fn open(&mut self, id: MiniAppId) {
        self.open_apps.insert(id);
    }

    pub fn toggle(&mut self, id: MiniAppId) {
        if !self.open_apps.remove(&id) {
            self.open_apps.insert(id);
        } else {
            self.expanded_apps.remove(&id);
        }
    }

    pub fn toggle_expand(&mut self, id: MiniAppId) {
        if !self.expanded_apps.remove(&id) {
            self.expanded_apps.insert(id);
        }
    }

    #[must_use]
    pub fn sub_tab(&self, id: MiniAppId) -> CompoundSubTab {
        self.sub_tabs
            .get(&id)
            .copied()
            .unwrap_or(CompoundSubTab::Primary)
    }

    pub fn set_sub_tab(&mut self, id: MiniAppId, tab: CompoundSubTab) {
        self.sub_tabs.insert(id, tab);
    }
}

/// Interactive state for the R-Group / PROTAC combinatorial enumerator inside `LibrarySarProtac`.
#[derive(Debug, Clone)]
pub struct EnumeratorState {
    pub mode: EnumerateMode,
    pub scaffold_smiles: String,
    pub warhead_smiles: String,
    pub e3_preset: E3RecruiterPreset,
    pub max_compounds: usize,
    pub last_library: Option<EnumeratedLibrary>,
    pub status_msg: String,
}

impl Default for EnumeratorState {
    fn default() -> Self {
        Self {
            mode: EnumerateMode::RGroup,
            scaffold_smiles: "Cc1ccc(-c2ccnc(N[*])n2)cc1".to_string(),
            warhead_smiles: "Cc1ccc(-c2ccnc(N[*])n2)cc1".to_string(),
            e3_preset: E3RecruiterPreset::All,
            max_compounds: 32,
            last_library: None,
            status_msg: "Introduce un scaffold con [*] / [R1] o un warhead PROTAC y pulsa 'Enumerar Biblioteca'."
                .to_string(),
        }
    }
}

impl EnumeratorState {
    pub fn run_enumeration(&mut self) -> usize {
        let max_n = self.max_compounds.clamp(1, 200);
        let lib = match self.mode {
            EnumerateMode::RGroup => {
                enumerate_rgroup_library(self.scaffold_smiles.trim(), None, None, max_n)
            }
            EnumerateMode::Protac => {
                enumerate_protac_library(self.warhead_smiles.trim(), self.e3_preset, None, max_n)
            }
        };
        let count = lib.valid_unique;
        self.status_msg = format!(
            "Enumerados {} compuestos válidos y únicos ({} intentos, {} descartados/duplicados).",
            lib.valid_unique, lib.attempted, lib.invalid_or_duplicate
        );
        self.last_library = Some(lib);
        count
    }
}

/// Paint a crisp, plain-Spanish educational banner at the top of every floating mini-application.
pub fn render_mini_app_header(ui: &mut egui::Ui, meta: &MiniAppMeta, is_expanded: bool) -> bool {
    let mut toggle_expand_clicked = false;
    let [r, g, b] = meta.accent_rgb;
    let accent = egui::Color32::from_rgb(r, g, b);
    let dark = ui.visuals().dark_mode;
    let title_col = if dark {
        accent
    } else {
        egui::Color32::from_rgb(
            (f32::from(r) * 0.72) as u8,
            (f32::from(g) * 0.72) as u8,
            (f32::from(b) * 0.72) as u8,
        )
    };
    let bg = egui::Color32::from_rgba_unmultiplied(r, g, b, if dark { 24 } else { 18 });

    egui::Frame::group(ui.style())
        .fill(bg)
        .stroke(egui::Stroke::new(1.2_f32, accent))
        .corner_radius(egui::CornerRadius::same(8))
        .inner_margin(egui::Margin::symmetric(10, 7))
        .show(ui, |ui| {
            ui.horizontal(|ui| {
                ui.colored_label(
                    title_col,
                    egui::RichText::new(meta.short_label).strong().size(13.5),
                );
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    let lbl = if is_expanded {
                        "🗗 Tamaño normal"
                    } else {
                        "🗖 Ampliar ventana"
                    };
                    if ui
                        .small_button(lbl)
                        .on_hover_text(
                            "Alterna entre tamaño normal y ventana ampliada (también puedes arrastrar y estirar cualquier borde)",
                        )
                        .clicked()
                    {
                        toggle_expand_clicked = true;
                    }
                });
            });
            ui.label(egui::RichText::new(meta.what_it_does).size(12.0));
            ui.small(format!("Cuándo usarla: {}", meta.when_to_use));
            ui.small(format!("Flujo I/O: {}", meta.io_summary));
        });
    ui.add_space(4.0);
    toggle_expand_clicked
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn workspace_windows_tracks_open_expand_and_subtab_state_cleanly() {
        let mut ws = WorkspaceWindows::default();
        assert!(ws.open_apps.contains(&MiniAppId::Tutorials));
        assert!(!ws.show_classic_inspector);
        assert!(ws.show_left_panel);
        assert!(!ws.left_panel_floating);

        ws.open(MiniAppId::LibrarySarProtac);
        assert!(ws.open_apps.contains(&MiniAppId::LibrarySarProtac));
        ws.set_sub_tab(MiniAppId::LibrarySarProtac, CompoundSubTab::Tertiary);
        assert_eq!(
            ws.sub_tab(MiniAppId::LibrarySarProtac),
            CompoundSubTab::Tertiary
        );

        ws.toggle_expand(MiniAppId::LibrarySarProtac);
        assert!(ws.expanded_apps.contains(&MiniAppId::LibrarySarProtac));

        // Toggling off closes the window and clears its expanded flag
        ws.toggle(MiniAppId::LibrarySarProtac);
        assert!(!ws.open_apps.contains(&MiniAppId::LibrarySarProtac));
        assert!(!ws.expanded_apps.contains(&MiniAppId::LibrarySarProtac));
    }

    #[test]
    fn all_mini_apps_have_distinct_ids_and_spanish_metadata() {
        let mut seen = BTreeSet::new();
        for meta in ALL_MINI_APPS {
            assert!(seen.insert(meta.id), "duplicate MiniAppId: {:?}", meta.id);
            assert!(!meta.short_label.is_empty());
            assert!(!meta.window_title.is_empty());
            assert!(!meta.what_it_does.is_empty());
            assert!(!meta.when_to_use.is_empty());
            assert!(!meta.io_summary.is_empty());
            assert!(meta.default_size[0] >= 380.0 && meta.default_size[1] >= 320.0);
            assert_eq!(mini_app_meta(meta.id).id, meta.id);
        }
        assert_eq!(seen.len(), 5);
    }

    #[test]
    fn enumerator_state_generates_valid_rgroup_and_protac_libraries() {
        let mut st = EnumeratorState::default();
        let n_rg = st.run_enumeration();
        assert!(n_rg >= 20, "expected >=20 R-group analogs, got {n_rg}");

        st.mode = EnumerateMode::Protac;
        st.max_compounds = 24;
        let n_p = st.run_enumeration();
        assert!(n_p >= 20, "expected >=20 PROTAC degraders, got {n_p}");
        assert!(
            st.last_library
                .as_ref()
                .is_some_and(|l| l.compounds.iter().all(|c| c.mol.is_some()))
        );
    }
}
