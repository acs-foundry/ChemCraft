#![forbid(unsafe_code)]

//! Interactive In-App Guided Tutorials for ChemCraft.
//!
//! Each tutorial walks the user step-by-step through a real medicinal-chemistry small-molecule
//! workflow inside ChemCraft, explaining *what* the tool does, *what to inspect*, and providing
//! a 1-click "▶ Ejecutar este paso ahora" button that mutates live `ChemCraftApp` state and opens
//! the corresponding floating mini-application.

use std::collections::BTreeSet;

use crate::miniapps::{CompoundSubTab, MiniAppId};

/// Target canvas mode requested by a tutorial step.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TutorialViewMode {
    View2D,
    Split2D3D,
}

/// Concrete live action triggered when the user clicks "▶ Ejecutar este paso ahora".
#[derive(Debug, Clone)]
pub enum TutorialAction {
    /// Load a molecule onto the 2D canvas, optionally clean 2D coordinates, switch view mode, and select Benzene ring tool.
    LoadMolecule2DTemplate {
        smiles: &'static str,
        name: &'static str,
        view_mode: TutorialViewMode,
    },
    /// Execute a real bioisosteric replacement (`-COOH -> Tetrazole`), refresh ADMET, and open Properties/ADMET.
    RunBioisostereAndInspect,
    /// Load a 3D drug molecule, regenerate EtKDG 3D coordinates, enable 3D view/spin, and open the 3D Viewer mini-app.
    BuildConformerAndOpen3D {
        smiles: &'static str,
        name: &'static str,
        view_mode: TutorialViewMode,
        enable_spin: bool,
    },
    /// Refresh on-demand ADMET summary for the active molecule and open the Properties/ADMET/CraftBio mini-app.
    ComputeAdmetAndOpen,
    /// Open a specific floating mini-application (and optional compound sub-tab).
    OpenMiniApp {
        app: MiniAppId,
        sub_tab: Option<CompoundSubTab>,
    },
    /// Load the 50-compound demo SAR series and show MMP & Activity Cliffs (SALI) inside Library/SAR/PROTAC.
    LoadDemoSarAndShowCliffs,
    /// Run a real in-silico virtual screening funnel (FDA Approved + Rule of 5 + 3D) in the Library mini-app.
    RunInSilicoLibraryDemo,
    /// Assemble a real bifunctional PROTAC library (CRBN/VHL/IAP) from a kinase warhead and display it.
    RunProtacEnumerationDemo,
    /// Load a kinase biaryl inhibitor and run Buy-Before-Make + Retro* retrosynthesis planning.
    RunRetrosynthesisDemo {
        smiles: &'static str,
        name: &'static str,
    },
    /// Send the planned synthesis step (or Aspirin preset) to the Stoichiometry & `.rxn` calculator.
    SendSynthStepToReactionCalc,
}

/// One guided step inside an interactive tutorial track.
#[derive(Debug, Clone)]
pub struct TutorialStep {
    pub title: &'static str,
    pub objective: &'static str,
    pub what_to_look_at: &'static str,
    pub action_label: &'static str,
    pub action: TutorialAction,
    pub related_app: Option<MiniAppId>,
}

/// Complete multi-step interactive tutorial track.
#[derive(Debug, Clone)]
pub struct TutorialTrack {
    pub id: &'static str,
    pub title: &'static str,
    pub duration_min: u8,
    pub summary: &'static str,
    pub steps: Vec<TutorialStep>,
}

/// Runtime state of the Interactive Tutorials mini-app.
#[derive(Debug, Clone)]
pub struct TutorialState {
    pub active_track_idx: usize,
    pub active_step_idx: usize,
    pub completed_steps: BTreeSet<(usize, usize)>,
    pub last_feedback: String,
}

impl Default for TutorialState {
    fn default() -> Self {
        Self {
            active_track_idx: 0,
            active_step_idx: 0,
            completed_steps: BTreeSet::new(),
            last_feedback: "Selecciona un tutorial y pulsa '▶ Ejecutar este paso ahora' para aprender practicando sobre ChemCraft."
                .to_string(),
        }
    }
}

impl TutorialState {
    pub fn mark_completed(&mut self, track_idx: usize, step_idx: usize, feedback: String) {
        self.completed_steps.insert((track_idx, step_idx));
        self.last_feedback = feedback;
    }

    #[must_use]
    pub fn completed_in_track(&self, track_idx: usize, total_steps: usize) -> usize {
        (0..total_steps)
            .filter(|&s| self.completed_steps.contains(&(track_idx, s)))
            .count()
    }
}

/// Build the 4 curated Spanish interactive tutorials covering every small-molecule ChemCraft workflow.
#[must_use]
pub fn built_in_tutorials() -> Vec<TutorialTrack> {
    vec![
        TutorialTrack {
            id: "sketcher_2d_bioisosteres",
            title: "1. Editor 2D, Plantillas y Bioisósteros",
            duration_min: 3,
            summary: "Domina el lienzo 2D vectorial, las plantillas de anillos y superátomos, la limpieza geométrica ACS, el reemplazo bioisostérico (-COOH → Tetrazol) y los espectros 1H/13C NMR.",
            steps: vec![
                TutorialStep {
                    title: "Cargar Ibuprofeno en lienzo 2D limpio y activar plantilla de Benceno",
                    objective: "Carga Ibuprofeno desde SMILES, limpia la geometría 2D estilo ACS y selecciona la plantilla de anillo Benceno lista para acoplar sobre cualquier átomo o enlace.",
                    what_to_look_at: "Observa el lienzo central en modo 2D puro: pasa el ratón sobre un átomo o enlace y haz clic para añadir anillos, o pulsa teclas C/N/O/F sobre un átomo para cambiar su elemento.",
                    action_label: "▶ Ejecutar este paso ahora (Cargar 2D + Plantilla)",
                    action: TutorialAction::LoadMolecule2DTemplate {
                        smiles: "CC(C)Cc1ccc(C(C)C(=O)O)cc1",
                        name: "Ibuprofen",
                        view_mode: TutorialViewMode::View2D,
                    },
                    related_app: None,
                },
                TutorialStep {
                    title: "Reemplazo bioisostérico (-COOH → Tetrazol) con delta ADMET",
                    objective: "Transforma automáticamente el ácido carboxílico en un anillo de tetrazol bioisostérico, calcula el perfil ADMET bajo demanda y muestra ΔMW, ΔLogP y ΔTPSA.",
                    what_to_look_at: "En la mini-app flotante 'Propiedades · NMR · CraftBio' verás el radar de biodisponibilidad de 6 ejes y la caja 'Bioisostere ADMET Delta'. Puedes deshacer con Ctrl+Z.",
                    action_label: "▶ Ejecutar este paso ahora (Mutar a Tetrazol + ADMET)",
                    action: TutorialAction::RunBioisostereAndInspect,
                    related_app: Some(MiniAppId::PropertiesAdmetCraftBio),
                },
                TutorialStep {
                    title: "Inspeccionar espectros simulados 1H y 13C NMR",
                    objective: "Cambia a la pestaña de RMN dentro de la mini-app de Propiedades para revisar desplazamientos químicos δ (ppm), multiplicidad e integración de la molécula modificada.",
                    what_to_look_at: "Comprueba los espectros de barras interactivos (ppm altos a la izquierda, estándar NMR) y las tablas de asignación de protones y carbonos.",
                    action_label: "▶ Ejecutar este paso ahora (Ver 1H / 13C NMR)",
                    action: TutorialAction::OpenMiniApp {
                        app: MiniAppId::PropertiesAdmetCraftBio,
                        sub_tab: Some(CompoundSubTab::Secondary),
                    },
                    related_app: Some(MiniAppId::PropertiesAdmetCraftBio),
                },
            ],
        },
        TutorialTrack {
            id: "conformers_3d_pharmacophore",
            title: "2. Confórmeros 3D y Perfil Farmacofórico",
            duration_min: 3,
            summary: "Genera coordenadas 3D EtKDG para un fármaco sulfonamida (Celecoxib), explora su conformación espacial interactiva, rasgos farmacofóricos (HBD/HBA/anillos) y exportación 3D.",
            steps: vec![
                TutorialStep {
                    title: "Generar confórmero 3D de Celecoxib y abrir Visor 3D flotante",
                    objective: "Carga Celecoxib en vista dividida 2D + 3D, regenera el confórmero 3D tetraédrico/aromático con hidrógenos explícitos visuales y abre la mini-app flotante 'Visor 3D y Confórmeros'.",
                    what_to_look_at: "Arrastra con el ratón sobre el visor 3D para rotar la sulfonamida y el grupo -CF3, usa la rueda para zoom y observa el resumen farmacofórico (HBD, HBA, anillos, extensión 3D en Å).",
                    action_label: "▶ Ejecutar este paso ahora (Celecoxib 3D + Farmacóforo)",
                    action: TutorialAction::BuildConformerAndOpen3D {
                        smiles: "Cc1ccc(-c2cc(C(F)(F)F)nn2-c2ccc(S(N)(=O)=O)cc2)cc1",
                        name: "Celecoxib",
                        view_mode: TutorialViewMode::Split2D3D,
                        enable_spin: true,
                    },
                    related_app: Some(MiniAppId::Viewer3dConformers),
                },
                TutorialStep {
                    title: "Verificar fórmula, masa exacta y aductos ESI+/ESI- MS",
                    objective: "Abre la pestaña de Análisis Molecular para auditar masa monoisotópica, DBE (insaturaciones), centros quirales, análisis elemental y picos ESI-MS.",
                    what_to_look_at: "Puedes arrastrar ambas ventanas flotantes lado a lado o pulsar '🗖 Ampliar ventana' para inspeccionar todos los aductos [M+H]+ / [M+Na]+ / [M-H]-.",
                    action_label: "▶ Ejecutar este paso ahora (Análisis Molecular y MS)",
                    action: TutorialAction::OpenMiniApp {
                        app: MiniAppId::PropertiesAdmetCraftBio,
                        sub_tab: Some(CompoundSubTab::Primary),
                    },
                    related_app: Some(MiniAppId::PropertiesAdmetCraftBio),
                },
                TutorialStep {
                    title: "Calcular perfil ADMET y revisar puentes a DockCraft / ChimeraX",
                    objective: "Evalúa Lipinski, LogP consenso, ESOL LogS, absorción GI y permeabilidad BBB bajo demanda, dejándolo listo para enviar en 1 clic a DockCraft o AdmetCraft.",
                    what_to_look_at: "Fíjate en la sección 'Puentes CraftBio' al final de la mini-app: permite lanzar AdmetCraft, DockCraft o ChimeraX con la geometría 3D ya generada.",
                    action_label: "▶ Ejecutar este paso ahora (Calcular ADMET + Puentes)",
                    action: TutorialAction::ComputeAdmetAndOpen,
                    related_app: Some(MiniAppId::PropertiesAdmetCraftBio),
                },
            ],
        },
        TutorialTrack {
            id: "library_sar_rgroup_protac",
            title: "3. Bibliotecas, SAR, Pares MMP y PROTAC",
            duration_min: 4,
            summary: "Explora series químicas con Activity Cliffs (SALI) y pares MMP, filtra bibliotecas comerciales in silico con embudo Rule-of-5 y ensambla degradadores PROTAC (CRBN/VHL).",
            steps: vec![
                TutorialStep {
                    title: "Cargar serie SAR demo (50 compuestos) y analizar Activity Cliffs / MMP",
                    objective: "Importa la serie de 50 inhibidores con pIC50, calcula descriptores ADMET y abre la pestaña 'MMP & Activity Cliffs' en la mini-app flotante de Biblioteca/SAR.",
                    what_to_look_at: "Identifica los pares de mayor SALI (pequeño cambio estructural → gran salto de potencia) y las transformaciones R1 → R2 con ΔpIC50, ΔcLogP y ΔTPSA.",
                    action_label: "▶ Ejecutar este paso ahora (Demo 50 SAR + Cliffs/MMP)",
                    action: TutorialAction::LoadDemoSarAndShowCliffs,
                    related_app: Some(MiniAppId::LibrarySarProtac),
                },
                TutorialStep {
                    title: "Filtrar biblioteca in silico FDA (Rule of 5 + diversidad MaxMin + 3D)",
                    objective: "Ejecuta el embudo de cribado virtual sobre el catálogo interno de fármacos aprobados por la FDA con ionización pH 7.4, selección diversa MaxMin y generación 3D.",
                    what_to_look_at: "En la sub-pestaña 'Bibliotecas (In Silico)' verás la tabla de hits filtrados listos para cargar en el lienzo 2D+3D, enviar al SAR o exportar como SDF 3D a DockCraft.",
                    action_label: "▶ Ejecutar este paso ahora (Cribado Biblioteca FDA 3D)",
                    action: TutorialAction::RunInSilicoLibraryDemo,
                    related_app: Some(MiniAppId::LibrarySarProtac),
                },
                TutorialStep {
                    title: "Ensamblar biblioteca combinatoria PROTAC (CRBN / VHL / IAP)",
                    objective: "Combina un warhead quinasa 2-aminopirimidina con punto de salida [*] frente a linkers alquilo/PEG/rígidos y reclutadores E3 (Pomalidomida, Lenalidomida, VH032, IAP).",
                    what_to_look_at: "En la sub-pestaña 'Enumerador R-Group / PROTAC' revisa los degradadores validados y deduplicados (MW, fórmula, bRo5) y pulsa 'Ver 2D+3D' o 'Enviar serie al SAR'.",
                    action_label: "▶ Ejecutar este paso ahora (Ensamblar PROTACs)",
                    action: TutorialAction::RunProtacEnumerationDemo,
                    related_app: Some(MiniAppId::LibrarySarProtac),
                },
            ],
        },
        TutorialTrack {
            id: "retrosynthesis_reactions",
            title: "4. Síntesis, Retrosíntesis Retro* y Reacciones",
            duration_min: 4,
            summary: "Evalúa disponibilidad comercial (Buy-Before-Make), planifica rutas sintéticas multipaso con Retro* y calcula estequiometría, economía atómica y rendimiento en MDL .rxn.",
            steps: vec![
                TutorialStep {
                    title: "Planificar ruta Retro* multipaso para un inhibidor quinasa biárilo",
                    objective: "Carga un inhibidor 2-aminopirimidina biárilo y ejecuta el motor Retro* + Buy-Before-Make para desconectarlo hasta bloques comerciales en stock.",
                    what_to_look_at: "En la mini-app 'Síntesis y Retro*' inspecciona los chips de rutas, la cobertura de stock (100% comprable), las miniaturas 2D paso a paso y la lista de compra con $/g.",
                    action_label: "▶ Ejecutar este paso ahora (Planificar Retro*)",
                    action: TutorialAction::RunRetrosynthesisDemo {
                        smiles: "Cc1ccc(-c2ccnc(Nc3ccc(S(N)(=O)=O)cc3)n2)cc1",
                        name: "Kinase_Biaryl_Inhibitor",
                    },
                    related_app: Some(MiniAppId::SynthesisRetro),
                },
                TutorialStep {
                    title: "Calcular estequiometría, economía atómica y rendimiento (.rxn)",
                    objective: "Transfiere el primer paso retrosintético a la calculadora de reacciones SMILES / MDL .rxn en escala de 2.5 mmol para estimar masas exactas y rendimiento aislado.",
                    what_to_look_at: "En la sub-pestaña 'Reacción y Estequiometría (.rxn)' verás Atom Economy (%), rendimiento teórico/aislado en mg y el botón 'Export MDL .rxn (V2000)'.",
                    action_label: "▶ Ejecutar este paso ahora (Estequiometría y .rxn)",
                    action: TutorialAction::SendSynthStepToReactionCalc,
                    related_app: Some(MiniAppId::SynthesisRetro),
                },
            ],
        },
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ChemCraftApp;

    #[test]
    fn built_in_tutorials_cover_four_tracks_with_non_empty_steps() {
        let tracks = built_in_tutorials();
        assert_eq!(tracks.len(), 4);
        let mut ids = BTreeSet::new();
        for t in &tracks {
            assert!(ids.insert(t.id), "duplicate tutorial id: {}", t.id);
            assert!(!t.steps.is_empty(), "track {} has no steps", t.id);
            assert!(t.duration_min > 0);
            for st in &t.steps {
                assert!(!st.title.is_empty());
                assert!(!st.objective.is_empty());
                assert!(!st.what_to_look_at.is_empty());
                assert!(!st.action_label.is_empty());
            }
        }
        let mut state = TutorialState::default();
        assert_eq!(state.completed_in_track(0, tracks[0].steps.len()), 0);
        state.mark_completed(0, 0, "ok".to_string());
        assert_eq!(state.completed_in_track(0, tracks[0].steps.len()), 1);
    }

    #[test]
    fn every_tutorial_step_executes_real_app_state_cleanly() {
        let mut app = ChemCraftApp::new();
        let tracks = built_in_tutorials();
        let mut saw_live_admet = false;
        for (t_idx, track) in tracks.iter().enumerate() {
            for (s_idx, step) in track.steps.iter().enumerate() {
                app.run_tutorial_action(t_idx, s_idx, step.action.clone());
                if app.session.live_admet.is_some() {
                    saw_live_admet = true;
                }
                assert!(
                    app.tutorial_state.completed_steps.contains(&(t_idx, s_idx)),
                    "step ({t_idx}, {s_idx}) '{}' was not marked completed",
                    step.title
                );
                assert!(
                    !app.tutorial_state.last_feedback.is_empty(),
                    "step ({t_idx}, {s_idx}) produced empty feedback"
                );
                assert!(
                    !app.session.canvas_mol.atoms.is_empty(),
                    "step ({t_idx}, {s_idx}) left canvas molecule empty"
                );
            }
        }
        // Verify side effects across all 4 small-molecule workflows actually ran real engines:
        assert!(saw_live_admet);
        assert!(app.library_win.last_result.is_some());
        assert!(app.enumerator.last_library.is_some());
        assert!(
            app.synth_module
                .last_plan
                .as_ref()
                .is_some_and(|p| !p.routes.is_empty())
        );
        assert!(!app.session.sar_table.rows.is_empty());
    }
}
