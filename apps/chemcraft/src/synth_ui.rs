#![forbid(unsafe_code)]

use std::collections::HashMap;
use std::path::PathBuf;

use egui::{Color32, RichText, Vec2};
use mol_graph::{MolGraph, parse_smiles};
use mol_synth::{
    RetroConfig, StockTier, SynthesisPlan, default_stock_catalog, plan_retrosynthesis,
};
use mol_vector::generate_2d_coords;

use crate::chem_panels::pick_save_file_native;
use crate::widgets::paint_mol_thumbnail;

/// Action emitted by the on-demand `Synth Module` panel back to `ChemCraftApp`.
pub enum SynthAction {
    None,
    LoadSmiles { smiles: String, name: String },
    OpenReactionSmiles { reaction_smiles: String },
    Status(String),
}

/// Persistent state of the on-demand `Synth Module` (hidden by default so the main 2D canvas stays uncluttered).
pub struct SynthModuleState {
    pub open: bool,
    pub max_steps: usize,
    pub max_routes: usize,
    pub target_scale_mg: f64,
    pub selected_route: usize,
    pub last_plan: Option<SynthesisPlan>,
    pub thumb_cache: HashMap<String, MolGraph>,
    pub status_msg: String,
}

impl Default for SynthModuleState {
    fn default() -> Self {
        Self {
            open: false,
            max_steps: 4,
            max_routes: 5,
            target_scale_mg: 100.0,
            selected_route: 0,
            last_plan: None,
            thumb_cache: HashMap::new(),
            status_msg:
                "Click 'Plan Synthesis' to run Buy-Before-Make + Retro* on the active molecule."
                    .to_string(),
        }
    }
}

impl SynthModuleState {
    pub fn run_for_molecule(&mut self, mol: &MolGraph) {
        if mol.atoms.is_empty() {
            self.status_msg = "Canvas molecule is empty — draw or load a target first.".to_string();
            return;
        }
        let cfg = RetroConfig {
            max_steps: self.max_steps.clamp(1, 6),
            max_routes: self.max_routes.clamp(1, 8),
            target_scale_mg: self.target_scale_mg.max(1.0),
            ..Default::default()
        };
        let plan = plan_retrosynthesis(mol, default_stock_catalog(), &cfg);
        self.thumb_cache.clear();
        let mut all_smis: Vec<String> = Vec::new();
        all_smis.push(plan.target_smiles.clone());
        for r in &plan.routes {
            for st in &r.steps {
                all_smis.push(st.product_smiles.clone());
                for rs in &st.reactant_smiles {
                    all_smis.push(rs.clone());
                }
            }
        }
        for smi in all_smis {
            if self.thumb_cache.contains_key(&smi) {
                continue;
            }
            if let Ok(mut g) = parse_smiles(&smi) {
                generate_2d_coords(&mut g);
                self.thumb_cache.insert(smi, g);
            }
        }
        self.status_msg = format!(
            "Found {} route(s) for {} ({} nodes, {} ms, {} catalog BBs)",
            plan.routes.len(),
            plan.target_name,
            plan.nodes_expanded,
            plan.elapsed_ms,
            plan.catalog_size
        );
        self.selected_route = 0;
        self.last_plan = Some(plan);
    }
}

/// Renders the **`Synth Module`** Buy-Before-Make + Retro* contents inside any `egui::Ui`.
pub fn render_synth_contents(
    ui: &mut egui::Ui,
    state: &mut SynthModuleState,
    canvas_mol: &MolGraph,
    export_dir: &str,
    dark_mode: bool,
) -> SynthAction {
    let mut action = SynthAction::None;

    // --- Top Controls Row ---
    ui.horizontal_wrapped(|ui| {
        let target_lbl = canvas_mol
            .name
            .as_deref()
            .filter(|s| !s.is_empty())
            .unwrap_or("Canvas Target");
        ui.strong(format!(
            "Target: {target_lbl} ({}, {:.1} Da)",
            canvas_mol.formula(),
            canvas_mol.molecular_weight()
        ));
        ui.separator();
        ui.label("Max Steps:");
        ui.add(egui::DragValue::new(&mut state.max_steps).range(1..=6));
        ui.label("Scale (mg):");
        ui.add(
            egui::DragValue::new(&mut state.target_scale_mg)
                .speed(10.0)
                .range(1.0..=50_000.0),
        );
        if ui
            .button(
                RichText::new("Plan Synthesis")
                    .strong()
                    .color(Color32::from_rgb(30, 120, 210)),
            )
            .clicked()
        {
            state.run_for_molecule(canvas_mol);
        }

        if let Some(ref plan) = state.last_plan {
            ui.separator();
            if ui.small_button("Export ELN (.md)").clicked() {
                let def = PathBuf::from(export_dir)
                    .join("synth_protocol.md")
                    .display()
                    .to_string();
                if let Some(path) =
                    pick_save_file_native("Export ELN Protocol", &def, "Markdown", "*.md")
                {
                    match std::fs::write(&path, plan.render_markdown()) {
                        Ok(()) => {
                            action = SynthAction::Status(format!(
                                "Synth Module: exported ELN protocol -> {path}"
                            ));
                        }
                        Err(e) => {
                            state.status_msg = format!("Export error: {e}");
                        }
                    }
                }
            }
            if ui.small_button("Export JSON").clicked() {
                let def = PathBuf::from(export_dir)
                    .join("synth_plan.json")
                    .display()
                    .to_string();
                if let Some(path) =
                    pick_save_file_native("Export Synthesis Plan", &def, "JSON", "*.json")
                    && let Ok(js) = serde_json::to_string_pretty(plan)
                {
                    match std::fs::write(&path, js) {
                        Ok(()) => {
                            action = SynthAction::Status(format!(
                                "Synth Module: exported JSON -> {path}"
                            ));
                        }
                        Err(e) => {
                            state.status_msg = format!("Export error: {e}");
                        }
                    }
                }
            }
        }
    });

    ui.small(&state.status_msg);
    ui.separator();

    let Some(plan) = state.last_plan.clone() else {
        ui.vertical_centered(|ui| {
            ui.add_space(32.0);
            ui.label(
                "Click 'Plan Synthesis' above to run Buy-Before-Make stock lookup + Retro* retrosynthesis.",
            );
            ui.weak(
                "Tip: Edit or load any hit from DockCraft / AdmetCraft on the canvas, then re-run.",
            );
        });
        return action;
    };

    // 1. Buy-Before-Make Status & Nearest Commercial Analogs
    let stock_col = match plan.target_stock_match.tier {
        StockTier::InCatalog => Color32::from_rgb(35, 145, 70),
        StockTier::SimpleFragmentEst => Color32::from_rgb(35, 130, 165),
        StockTier::NotInStock => Color32::from_rgb(185, 110, 20),
    };
    ui.horizontal_wrapped(|ui| {
        ui.strong("Buy-Before-Make:");
        ui.colored_label(
            stock_col,
            RichText::new(plan.target_stock_match.tier.as_str()).strong(),
        );
        if let Some(ref entry) = plan.target_stock_match.entry {
            ui.label(format!(
                "— {} ({}, ~${:.2}/g, {})",
                entry.name, entry.catalog_id, entry.est_usd_per_g, entry.supplier
            ));
        }
    });

    if !plan.nearest_commercial_analogs.is_empty() {
        egui::CollapsingHeader::new(format!(
            "Nearest Commercial Building Blocks ({})",
            plan.nearest_commercial_analogs.len()
        ))
        .default_open(false)
        .show(ui, |ui| {
            egui::Grid::new("synth_analogs_grid")
                .striped(true)
                .spacing([10.0, 4.0])
                .show(ui, |ui| {
                    for h in ["ID", "Name", "Tanimoto", "Est. $/g", "Action"] {
                        ui.strong(h);
                    }
                    ui.end_row();
                    for a in &plan.nearest_commercial_analogs {
                        ui.monospace(&a.catalog_id);
                        ui.label(&a.name);
                        ui.label(format!("{:.3}", a.tanimoto));
                        ui.label(format!("${:.2}", a.est_usd_per_g));
                        if ui.small_button("Load 2D").clicked() {
                            action = SynthAction::LoadSmiles {
                                smiles: a.canonical_smiles.clone(),
                                name: a.name.clone(),
                            };
                        }
                        ui.end_row();
                    }
                });
        });
    }

    ui.separator();

    if plan.routes.is_empty() {
        ui.colored_label(
            Color32::from_rgb(190, 90, 20),
            "No retrosynthetic disconnections found within the current step limit.",
        );
        return action;
    }

    // 2. Route Selector Chips
    ui.horizontal_wrapped(|ui| {
        ui.strong("Routes:");
        for (idx, r) in plan.routes.iter().enumerate() {
            let lbl = format!(
                "#{} ({} step{}, {:.0}% stock, {:.0}% yld)",
                r.rank,
                r.total_steps,
                if r.total_steps == 1 { "" } else { "s" },
                r.stock_coverage_pct,
                r.overall_est_yield_pct
            );
            if ui
                .selectable_label(state.selected_route == idx, lbl)
                .clicked()
            {
                state.selected_route = idx;
            }
        }
    });
    ui.separator();

    let r_idx = state.selected_route.min(plan.routes.len() - 1);
    let route = &plan.routes[r_idx];

    ui.horizontal_wrapped(|ui| {
        let cov_col = if route.fully_resolved_in_stock {
            Color32::from_rgb(35, 145, 70)
        } else {
            Color32::from_rgb(185, 115, 20)
        };
        ui.colored_label(
            cov_col,
            RichText::new(format!(
                "Stock Coverage: {:.0}% ({})",
                route.stock_coverage_pct,
                if route.fully_resolved_in_stock {
                    "All Starting Materials Buyable"
                } else {
                    "Partial"
                }
            ))
            .strong(),
        );
        ui.separator();
        ui.label(format!(
            "Steps: {} (LLS {}) | Overall Yield: {:.1}% | Atom Econ: {:.1}% | Score: {:.1}/100",
            route.total_steps,
            route.longest_linear_seq,
            route.overall_est_yield_pct,
            route.mean_atom_economy_pct,
            route.composite_score
        ));
    });
    ui.add_space(6.0);

    // 3. Forward Synthesis Steps with Live 2D Thumbnails
    for st in &route.steps {
        egui::Frame::group(ui.style()).show(ui, |ui| {
            ui.horizontal_wrapped(|ui| {
                ui.strong(format!("Step {}: {}", st.step_number, st.reaction_name));
                ui.weak(format!("[{}]", st.reaction_class));
                ui.separator();
                ui.label(format!(
                    "Est. Yield: {:.0}% | AE: {:.1}% | Scale: {:.1} mg",
                    st.est_yield_pct, st.atom_economy_pct, st.product_scale_mg
                ));
                if ui.small_button("Send to Rxn Tab").clicked() {
                    action = SynthAction::OpenReactionSmiles {
                        reaction_smiles: st.reaction_smiles.clone(),
                    };
                }
                if ui.small_button("Export .rxn").clicked()
                    && let Some(rxn_txt) = st.to_rxn_v2000()
                {
                    let def = PathBuf::from(export_dir)
                        .join(format!("step_{:02}.rxn", st.step_number))
                        .display()
                        .to_string();
                    if let Some(path) =
                        pick_save_file_native("Export MDL RXN V2000", &def, "MDL Reaction", "*.rxn")
                    {
                        match std::fs::write(&path, rxn_txt) {
                            Ok(()) => {
                                action = SynthAction::Status(format!(
                                    "Exported Step {} -> {path}",
                                    st.step_number
                                ));
                            }
                            Err(e) => {
                                state.status_msg = format!("RXN write error: {e}");
                            }
                        }
                    }
                }
            });

            ui.horizontal_wrapped(|ui| {
                ui.label(format!(
                    "Conditions: {}  |  Solvent: {}  |  Temp/Time: {}",
                    st.reagents, st.solvent, st.temperature
                ));
            });

            // 2D thumbnails: Reactant(s) -> Product
            ui.horizontal_wrapped(|ui| {
                for (ri, (rsmi, mmol, mg)) in st.reactant_masses_mg.iter().enumerate() {
                    if ri > 0 {
                        ui.strong("  +  ");
                    }
                    ui.vertical(|ui| {
                        if let Some(g) = state.thumb_cache.get(rsmi) {
                            let _ =
                                paint_mol_thumbnail(ui, g, &[], Vec2::new(150.0, 95.0), dark_mode);
                        }
                        ui.small(format!("{mg:.1} mg ({mmol:.3} mmol)"));
                        if ui.small_button("Load 2D").clicked() {
                            action = SynthAction::LoadSmiles {
                                smiles: rsmi.clone(),
                                name: format!("Step{}_R{}", st.step_number, ri + 1),
                            };
                        }
                    });
                }
                ui.strong("  ==>  ");
                ui.vertical(|ui| {
                    if let Some(g) = state.thumb_cache.get(&st.product_smiles) {
                        let _ = paint_mol_thumbnail(ui, g, &[], Vec2::new(160.0, 95.0), dark_mode);
                    }
                    ui.small(format!("Product: {:.1} mg", st.product_scale_mg));
                    if ui.small_button("Load 2D").clicked() {
                        action = SynthAction::LoadSmiles {
                            smiles: st.product_smiles.clone(),
                            name: format!("Step{}_Product", st.step_number),
                        };
                    }
                });
            });

            ui.weak(format!("Workup: {}", st.workup_notes));
            for w in &st.warnings {
                ui.colored_label(
                    Color32::from_rgb(200, 100, 15),
                    format!("Chemoselectivity / PG Alert: {w}"),
                );
            }
        });
        ui.add_space(4.0);
    }

    ui.separator();
    ui.strong(format!(
        "Buy-Before-Make Starting Materials ({})",
        route.starting_materials.len()
    ));
    egui::Grid::new("synth_shopping_grid")
        .striped(true)
        .spacing([10.0, 4.0])
        .show(ui, |ui| {
            for h in [
                "Status",
                "Building Block",
                "Formula",
                "MW",
                "Required",
                "Est. $/g",
                "Action",
            ] {
                ui.strong(h);
            }
            ui.end_row();
            for sm in &route.starting_materials {
                let sc = match sm.tier {
                    StockTier::InCatalog => Color32::from_rgb(35, 145, 70),
                    StockTier::SimpleFragmentEst => Color32::from_rgb(35, 130, 165),
                    StockTier::NotInStock => Color32::from_rgb(195, 50, 50),
                };
                ui.colored_label(sc, sm.tier.as_str());
                let nm = sm
                    .catalog_id
                    .as_deref()
                    .map_or_else(|| sm.name.clone(), |cid| format!("{} ({cid})", sm.name));
                ui.label(nm);
                ui.monospace(&sm.formula);
                ui.label(format!("{:.1}", sm.mw));
                ui.label(format!("{:.1} mg", sm.required_mg));
                ui.label(
                    sm.est_usd_per_g
                        .map_or_else(|| "—".to_string(), |p| format!("${p:.2}")),
                );
                if ui.small_button("Load 2D").clicked() {
                    action = SynthAction::LoadSmiles {
                        smiles: sm.canonical_smiles.clone(),
                        name: sm.name.clone(),
                    };
                }
                ui.end_row();
            }
        });

    action
}
