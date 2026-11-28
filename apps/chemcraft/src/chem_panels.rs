//! ChemDraw-grade Inspector Tabs (`Analysis`, `1H/13C NMR`, `Reaction & Stoichiometry`)
//! and lightweight native OS file picker helper (`zenity` / `kdialog`).

use std::path::Path;
use std::process::Command as ProcessCommand;

use eframe::egui;
use egui::{Color32, Pos2, Rect, RichText, Sense, Stroke, Vec2};
use mol_admet_bridge::{NmrPeak, predict_nmr};
use mol_engine::PathSandbox;
use mol_graph::{
    MolGraph, ReactionRole, analyze_molecule, compute_stoichiometry, parse_reaction_smiles,
    write_rxn_v2000,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InspectorTab {
    Admet,
    Analysis,
    Nmr,
    Reaction,
}

#[derive(Debug, Clone)]
pub struct ReactionPanelState {
    pub rsmiles: String,
    pub limiting_mmol: String,
    pub actual_mg: String,
    pub status: String,
}

impl Default for ReactionPanelState {
    fn default() -> Self {
        Self {
            rsmiles: "CC(=O)O.CCO>>CC(=O)OCC.O".to_string(),
            limiting_mmol: "1.0".to_string(),
            actual_mg: "75.0".to_string(),
            status: String::new(),
        }
    }
}

/// Launch native Linux Save As dialog (`zenity` or `kdialog`).
pub fn pick_save_file_native(
    title: &str,
    default_path: &str,
    filter_label: &str,
    filter_glob: &str,
) -> Option<String> {
    let title_arg = format!("--title={title}");
    let fname_arg = format!("--filename={default_path}");
    let filter_arg = format!("--file-filter={filter_label} ({filter_glob}) | {filter_glob}");
    if let Ok(out) = ProcessCommand::new("zenity")
        .args([
            "--file-selection",
            "--save",
            "--confirm-overwrite",
            &title_arg,
            &fname_arg,
            &filter_arg,
            "--file-filter=All files | *",
        ])
        .output()
        && out.status.success()
    {
        let s = String::from_utf8_lossy(&out.stdout).trim().to_string();
        if !s.is_empty() {
            return Some(s);
        }
    }
    let kd_filter = format!("{filter_glob}|{filter_label}");
    if let Ok(out) = ProcessCommand::new("kdialog")
        .args(["--getsavefilename", default_path, &kd_filter])
        .output()
        && out.status.success()
    {
        let s = String::from_utf8_lossy(&out.stdout).trim().to_string();
        if !s.is_empty() {
            return Some(s);
        }
    }
    None
}

/// Render ChemDraw-equivalent Analysis Window (Exact Mass, MW, DBE, Elemental Analysis, MS Adducts).
pub fn render_analysis_tab(ui: &mut egui::Ui, mol: &MolGraph) {
    if mol.atoms.is_empty() {
        ui.label("Draw or load a molecule to inspect exact mass & elemental analysis.");
        return;
    }
    let a = analyze_molecule(mol);
    ui.strong("Molecular Analysis");
    ui.separator();

    egui::Grid::new("analysis_core_grid")
        .num_columns(2)
        .spacing([10.0, 4.0])
        .striped(true)
        .show(ui, |ui| {
            ui.label("Formula:");
            ui.strong(&a.formula_with_charge);
            ui.end_row();
            ui.label("Exact Mass:");
            ui.monospace(format!("{:.4} Da", a.exact_mass));
            ui.end_row();
            ui.label("Mol. Weight:");
            ui.monospace(format!("{:.2} Da", a.molecular_weight));
            ui.end_row();
            ui.label("DBE (RDB):");
            ui.monospace(format!("{:.1}", a.dbe));
            ui.end_row();
            ui.label("Net Charge:");
            ui.monospace(format!("{:+}", a.net_charge));
            ui.end_row();
            ui.label("Atoms / Rings:");
            ui.monospace(format!(
                "{} ({} heavy) / {}R (Chiral: {})",
                a.total_atoms, a.heavy_atoms, a.ring_count, a.chiral_centers
            ));
            ui.end_row();
        });

    ui.add_space(6.0);
    ui.strong("Elemental Analysis (%)");
    ui.add(
        egui::Label::new(
            RichText::new(&a.elemental_analysis_string)
                .monospace()
                .small(),
        )
        .wrap(),
    );
    if ui.small_button("Copy Elemental Analysis").clicked() {
        ui.ctx().copy_text(a.elemental_analysis_string.clone());
    }

    ui.add_space(6.0);
    ui.strong("ESI+/ESI- Mass Spec Adducts (m/z)");
    egui::Grid::new("analysis_adducts_grid")
        .num_columns(3)
        .spacing([10.0, 3.0])
        .striped(true)
        .show(ui, |ui| {
            ui.strong("Ion");
            ui.strong("Mode");
            ui.strong("m/z");
            ui.end_row();
            for ad in &a.adducts {
                ui.monospace(&ad.label);
                ui.label(RichText::new(&ad.mode).small().color(if ad.mode == "ESI+" {
                    Color32::from_rgb(35, 130, 210)
                } else {
                    Color32::from_rgb(200, 75, 75)
                }));
                ui.monospace(format!("{:.4}", ad.mz));
                ui.end_row();
            }
        });
}

/// Render ChemDraw-equivalent simulated 1H & 13C NMR stick spectra + peak tables.
pub fn render_nmr_tab(ui: &mut egui::Ui, mol: &MolGraph, dark_mode: bool) {
    if mol.atoms.is_empty() {
        ui.label("Draw or load a molecule to predict 1H & 13C NMR spectra.");
        return;
    }
    let spec = predict_nmr(mol);

    ui.strong(format!("Predicted 1H NMR ({} groups)", spec.h1_peaks.len()));
    paint_nmr_stick_plot(ui, &spec.h1_peaks, 12.5, 0.0, "1H δ (ppm)", dark_mode, true);
    egui::Grid::new("nmr_h1_grid")
        .num_columns(4)
        .striped(true)
        .show(ui, |ui| {
            ui.strong("δ (ppm)");
            ui.strong("Mult");
            ui.strong("Int");
            ui.strong("Environment");
            ui.end_row();
            for pk in &spec.h1_peaks {
                ui.monospace(format!("{:.2}", pk.shift_ppm));
                ui.monospace(&pk.multiplicity);
                ui.monospace(format!("{}H", pk.integration));
                ui.label(RichText::new(&pk.environment).small());
                ui.end_row();
            }
        });

    ui.add_space(8.0);
    ui.strong(format!(
        "Predicted 13C NMR ({} signals)",
        spec.c13_peaks.len()
    ));
    paint_nmr_stick_plot(
        ui,
        &spec.c13_peaks,
        220.0,
        0.0,
        "13C δ (ppm)",
        dark_mode,
        false,
    );
    egui::Grid::new("nmr_c13_grid")
        .num_columns(3)
        .striped(true)
        .show(ui, |ui| {
            ui.strong("δ (ppm)");
            ui.strong("C");
            ui.strong("Environment");
            ui.end_row();
            for pk in &spec.c13_peaks {
                ui.monospace(format!("{:.1}", pk.shift_ppm));
                ui.monospace(format!("{}", pk.integration));
                ui.label(RichText::new(&pk.environment).small());
                ui.end_row();
            }
        });
}

fn paint_nmr_stick_plot(
    ui: &mut egui::Ui,
    peaks: &[NmrPeak],
    ppm_max: f64,
    ppm_min: f64,
    axis_label: &str,
    dark_mode: bool,
    is_h1: bool,
) {
    let (rect, _) = ui.allocate_exact_size(Vec2::new(280.0, 74.0), Sense::hover());
    let painter = ui.painter_at(rect);
    let bg = if dark_mode {
        Color32::from_rgb(26, 28, 35)
    } else {
        Color32::from_rgb(248, 250, 253)
    };
    let axis_col = if dark_mode {
        Color32::from_gray(160)
    } else {
        Color32::from_gray(90)
    };
    let peak_col = if is_h1 {
        Color32::from_rgb(35, 130, 215)
    } else {
        Color32::from_rgb(200, 72, 62)
    };
    painter.rect_filled(rect, 3.0, bg);
    painter.rect_stroke(
        rect,
        3.0,
        Stroke::new(0.8_f32, Color32::from_gray(130)),
        egui::StrokeKind::Inside,
    );

    let baseline_y = rect.max.y - 16.0;
    let left_x = rect.min.x + 8.0;
    let right_x = rect.max.x - 8.0;
    painter.line_segment(
        [
            Pos2::new(left_x, baseline_y),
            Pos2::new(right_x, baseline_y),
        ],
        Stroke::new(1.0_f32, axis_col),
    );

    let max_int = peaks
        .iter()
        .map(|p| p.integration)
        .max()
        .unwrap_or(1)
        .max(1) as f64;
    let span = (ppm_max - ppm_min).max(1.0);
    for pk in peaks {
        // NMR convention: high ppm on LEFT, 0 ppm on RIGHT!
        let frac = ((ppm_max - pk.shift_ppm) / span).clamp(0.0, 1.0);
        let x = left_x + (frac as f32) * (right_x - left_x);
        let h_frac = if is_h1 {
            (0.35 + 0.65 * (pk.integration as f64 / max_int)) as f32
        } else {
            0.78_f32
        };
        let top_y = baseline_y - h_frac * (rect.height() - 26.0);
        painter.line_segment(
            [Pos2::new(x, baseline_y), Pos2::new(x, top_y)],
            Stroke::new(1.5_f32, peak_col),
        );
    }

    painter.text(
        Pos2::new(left_x + 2.0, rect.max.y - 12.0),
        egui::Align2::LEFT_TOP,
        format!("{ppm_max:.0}"),
        egui::FontId::monospace(9.0),
        axis_col,
    );
    painter.text(
        Pos2::new(rect.center().x, rect.max.y - 12.0),
        egui::Align2::CENTER_TOP,
        axis_label,
        egui::FontId::proportional(9.5),
        axis_col,
    );
    painter.text(
        Pos2::new(right_x - 2.0, rect.max.y - 12.0),
        egui::Align2::RIGHT_TOP,
        format!("{ppm_min:.0}"),
        egui::FontId::monospace(9.0),
        axis_col,
    );
    let _ = Rect::NOTHING;
}

/// Render Reaction SMILES / `.rxn` V2000 & Stoichiometry / Atom Economy / Yield Calculator.
pub fn render_reaction_tab(
    ui: &mut egui::Ui,
    state: &mut ReactionPanelState,
    canvas_mol: &MolGraph,
    export_dir: &str,
) {
    ui.strong("Reaction & Stoichiometry (.rxn)");
    ui.separator();
    ui.label(RichText::new("Reaction SMILES (R1.R2>Agents>P1.P2):").small());
    ui.add(
        egui::TextEdit::singleline(&mut state.rsmiles)
            .desired_width(275.0)
            .font(egui::TextStyle::Monospace),
    );
    ui.horizontal_wrapped(|ui| {
        if ui.small_button("Use Canvas as Reactant").clicked() && !canvas_mol.atoms.is_empty() {
            let smi = canvas_mol.canonical_smiles();
            state.rsmiles = format!("{smi}>>{smi}");
        }
        if ui.small_button("Preset: Fischer Esterification").clicked() {
            state.rsmiles = "CC(=O)O.CCO>>CC(=O)OCC.O".to_string();
        }
        if ui.small_button("Preset: Aspirin Synthesis").clicked() {
            state.rsmiles =
                "O=C(O)c1ccccc1O.CC(=O)OC(=O)C>>CC(=O)Oc1ccccc1C(=O)O.CC(=O)O".to_string();
        }
    });

    ui.horizontal(|ui| {
        ui.label("Scale (mmol):");
        ui.add(egui::TextEdit::singleline(&mut state.limiting_mmol).desired_width(58.0));
        ui.label("Actual (mg):");
        ui.add(egui::TextEdit::singleline(&mut state.actual_mg).desired_width(58.0));
    });

    let mmol = state.limiting_mmol.trim().parse::<f64>().unwrap_or(1.0);
    let actual_mg = state.actual_mg.trim().parse::<f64>().ok();

    match parse_reaction_smiles(&state.rsmiles) {
        Ok(rxn) => {
            let st = compute_stoichiometry(&rxn, mmol, actual_mg);
            ui.add_space(4.0);
            ui.horizontal_wrapped(|ui| {
                ui.strong(format!("Atom Economy: {:.1}%", st.atom_economy_pct));
                ui.separator();
                ui.strong(format!("Theo. Yield: {:.2} mg", st.theoretical_yield_mg));
                if let Some(yp) = st.isolated_yield_pct {
                    ui.separator();
                    ui.colored_label(
                        Color32::from_rgb(38, 152, 72),
                        RichText::new(format!("Yield: {yp:.1}%")).strong(),
                    );
                }
            });
            ui.add_space(4.0);
            egui::Grid::new("rxn_stoich_grid")
                .num_columns(5)
                .striped(true)
                .show(ui, |ui| {
                    ui.strong("Role");
                    ui.strong("Formula");
                    ui.strong("MW");
                    ui.strong("mmol");
                    ui.strong("mg");
                    ui.end_row();
                    for row in &st.rows {
                        let role_txt = match row.role {
                            ReactionRole::Reactant => "React.",
                            ReactionRole::Agent => "Agent",
                            ReactionRole::Product => "Prod.",
                        };
                        ui.label(role_txt);
                        ui.monospace(&row.formula);
                        ui.monospace(format!("{:.1}", row.mw));
                        ui.monospace(format!("{:.2}", row.mmol));
                        ui.monospace(format!("{:.1}", row.mass_mg));
                        ui.end_row();
                    }
                });
            ui.add_space(4.0);
            if ui.button("Export MDL .rxn (V2000)").clicked() {
                let root = Path::new(export_dir);
                match PathSandbox::new(None, Some(root))
                    .and_then(|sb| sb.open_write(Path::new("chemcraft_reaction.rxn")))
                {
                    Ok((resolved, mut f)) => {
                        use std::io::Write;
                        let content = write_rxn_v2000(&rxn);
                        match f.write_all(content.as_bytes()) {
                            Ok(()) => {
                                state.status = format!("Exported {}", resolved.display());
                            }
                            Err(e) => state.status = format!("Write error: {e}"),
                        }
                    }
                    Err(e) => state.status = format!("Sandbox error: {e}"),
                }
            }
            if !state.status.is_empty() {
                ui.label(RichText::new(&state.status).small().weak());
            }
        }
        Err(e) => {
            ui.colored_label(Color32::from_rgb(210, 65, 65), format!("Invalid Rxn: {e}"));
        }
    }
}
