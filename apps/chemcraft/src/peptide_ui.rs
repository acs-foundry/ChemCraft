#![forbid(unsafe_code)]

//! On-demand **`Peptide & NRP Studio`** window (hidden by default so the main 2D/3D canvas
//! stays uncluttered), unifying monomer/macrocycle/HELM construction, NRPS modular BGC assembly &
//! Stachelhaus specificity prediction, bRo5 developability/proteolysis/MS analysis, SPPS planning,
//! and combinatorial optimization scans (`Ala`, `D-AA`, `N-Me`, `RCM Staple`).

use std::path::PathBuf;

use egui::{Color32, RichText};
use mol_graph::MolGraph;
use mol_peptide::{
    BuiltPeptide, MonomerCatalog, MonomerClass, NrpsBgc, NrpsCompileReport, PeptideConformation,
    PeptideDevelopabilityReport, PeptideScanBundle, ScanMode, SppsPlan,
    analyze_peptide_developability, build_peptide_or_nrp, compile_nrps_bgc,
    default_monomer_catalog, design_nrps_from_peptide, landmark_nrp_presets, parse_peptide_spec,
    plan_spps_vs_nrps, predict_a_domain_specificity, preset_nrps_bgc, propose_nrps_swap,
    render_developability_markdown, render_nrps_markdown, render_scan_csv, render_spps_markdown,
    run_peptide_scan,
};
use mol_sar::SarTable;

use crate::chem_panels::pick_save_file_native;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PeptideStudioTab {
    Builder,
    NrpsBgc,
    Developability,
    ScanOptimize,
}

/// Actions emitted by `Peptide & NRP Studio` back to `ChemCraftApp`.
pub enum PeptideStudioAction {
    None,
    LoadMolOnCanvas(MolGraph),
    LoadScanIntoSar(SarTable),
    SendToDockCraft(MolGraph),
    OpenInChimeraX(MolGraph),
    Status(String),
}

pub struct PeptideStudioState {
    pub open: bool,
    pub tab: PeptideStudioTab,
    pub catalog: MonomerCatalog,
    pub notation: String,
    pub conformation: PeptideConformation,
    pub palette_class: MonomerClass,
    pub insert_as_d: bool,
    pub insert_as_nme: bool,
    pub last_built: Option<BuiltPeptide>,
    pub last_dev: Option<PeptideDevelopabilityReport>,
    pub last_spps: Option<SppsPlan>,
    // NRPS tab state
    pub stachelhaus_query: String,
    pub active_bgc: Option<NrpsBgc>,
    pub last_nrps_report: Option<NrpsCompileReport>,
    pub swap_module_idx: usize,
    pub swap_new_monomer: String,
    pub last_swap_note: String,
    // Scan tab state
    pub scan_mode: ScanMode,
    pub last_scan: Option<PeptideScanBundle>,
    pub status_msg: String,
}

impl Default for PeptideStudioState {
    fn default() -> Self {
        Self {
            open: false,
            tab: PeptideStudioTab::Builder,
            catalog: default_monomer_catalog(),
            notation: "cyclo(Arg-Gly-Asp-D-Phe-NMe-Val)".to_string(),
            conformation: PeptideConformation::MacrocycleRing,
            palette_class: MonomerClass::Canonical,
            insert_as_d: false,
            insert_as_nme: false,
            last_built: None,
            last_dev: None,
            last_spps: None,
            stachelhaus_query: "DAWTIAAVCK".to_string(),
            active_bgc: preset_nrps_bgc("gramicidin_s"),
            last_nrps_report: None,
            swap_module_idx: 3,
            swap_new_monomer: "Ile".to_string(),
            last_swap_note: String::new(),
            scan_mode: ScanMode::NMethyl,
            last_scan: None,
            status_msg: "Select a preset or compose any L/D/NMe/NRP sequence or HELM string."
                .to_string(),
        }
    }
}

impl PeptideStudioState {
    pub fn build_current(&mut self) -> Option<MolGraph> {
        match build_peptide_or_nrp(&self.notation, self.conformation, &self.catalog) {
            Ok(built) => {
                let dev = analyze_peptide_developability(&built, &self.catalog);
                let spps = plan_spps_vs_nrps(&built.spec, &self.catalog);
                let mol = built.mol.clone();
                self.status_msg = format!(
                    "Built {} ({}, {:.4} Da, bRo5 {:.2}, {} active protease site(s))",
                    dev.name,
                    dev.formula,
                    dev.monoisotopic_mass,
                    dev.bro5_oral_score,
                    dev.active_protease_sites
                );
                self.last_built = Some(built);
                self.last_dev = Some(dev);
                self.last_spps = Some(spps);
                Some(mol)
            }
            Err(e) => {
                self.status_msg = format!("Build error: {e}");
                None
            }
        }
    }

    pub fn compile_active_bgc(&mut self) -> Option<MolGraph> {
        let bgc = self.active_bgc.clone()?;
        match compile_nrps_bgc(&bgc, &self.catalog) {
            Ok((rep, built)) => {
                let dev = analyze_peptide_developability(&built, &self.catalog);
                let spps = plan_spps_vs_nrps(&built.spec, &self.catalog);
                let mol = built.mol.clone();
                self.notation = built.notation.clone();
                self.status_msg = format!(
                    "Compiled BGC {} -> {} ({}, {:.4} Da, collinear={})",
                    rep.bgc.name,
                    rep.formula,
                    rep.predicted_notation,
                    rep.monoisotopic_mass,
                    rep.collinear_compatible
                );
                self.last_nrps_report = Some(rep);
                self.last_built = Some(built);
                self.last_dev = Some(dev);
                self.last_spps = Some(spps);
                Some(mol)
            }
            Err(e) => {
                self.status_msg = format!("NRPS compile error: {e}");
                None
            }
        }
    }

    fn append_monomer_token(&mut self, base_code: &str) {
        let mut tok = String::new();
        if self.insert_as_nme && !base_code.eq_ignore_ascii_case("Sar") {
            tok.push_str("NMe-");
        }
        if self.insert_as_d {
            tok.push_str("D-");
        }
        tok.push_str(base_code);
        if self.notation.trim().is_empty() {
            self.notation = tok;
        } else if self.notation.ends_with(')') && self.notation.starts_with("cyclo(") {
            self.notation.pop();
            self.notation.push('-');
            self.notation.push_str(&tok);
            self.notation.push(')');
        } else {
            self.notation.push('-');
            self.notation.push_str(&tok);
        }
    }
}

/// Render the on-demand **`Peptide & NRP Studio`** window when `state.open == true`.
pub fn render_peptide_studio(
    ctx: &egui::Context,
    state: &mut PeptideStudioState,
    export_dir: &str,
) -> PeptideStudioAction {
    if !state.open {
        return PeptideStudioAction::None;
    }

    let mut action = PeptideStudioAction::None;
    let mut open = state.open;

    egui::Window::new("Peptide & NRP Studio (v0.4)")
        .open(&mut open)
        .default_width(840.0)
        .default_height(640.0)
        .resizable(true)
        .show(ctx, |ui| {
            // Top tab strip + quick ecosystem handoff
            ui.horizontal_wrapped(|ui| {
                ui.selectable_value(&mut state.tab, PeptideStudioTab::Builder, "1. Builder & Palette");
                ui.selectable_value(&mut state.tab, PeptideStudioTab::NrpsBgc, "2. NRPS Assembly Line");
                ui.selectable_value(
                    &mut state.tab,
                    PeptideStudioTab::Developability,
                    "3. Developability, MS & SPPS",
                );
                ui.selectable_value(&mut state.tab, PeptideStudioTab::ScanOptimize, "4. Scan & Optimize");

                if let Some(ref built) = state.last_built {
                    ui.separator();
                    if ui.small_button("Load 2D+3D").clicked() {
                        action = PeptideStudioAction::LoadMolOnCanvas(built.mol.clone());
                    }
                    if ui.small_button("DockCraft").clicked() {
                        action = PeptideStudioAction::SendToDockCraft(built.mol.clone());
                    }
                    if ui.small_button("ChimeraX").clicked() {
                        action = PeptideStudioAction::OpenInChimeraX(built.mol.clone());
                    }
                }
            });
            ui.small(&state.status_msg);
            ui.separator();

            match state.tab {
                PeptideStudioTab::Builder => {
                    egui::ScrollArea::vertical().show(ui, |ui| {
                        ui.strong("Landmark Natural-Product NRP & Therapeutic Presets");
                        ui.horizontal_wrapped(|ui| {
                            for p in landmark_nrp_presets() {
                                if ui
                                    .button(p.name)
                                    .on_hover_text(format!("{} — {} ({})", p.organism_or_class, p.summary, p.expected_formula))
                                    .clicked()
                                {
                                    state.notation = p.notation.to_string();
                                    if let Some(mol) = state.build_current() {
                                        action = PeptideStudioAction::LoadMolOnCanvas(mol);
                                    }
                                }
                            }
                        });

                        ui.add_space(6.0);
                        ui.strong("Extended Sequence / HELM / Cross-Link Editor");
                        ui.add(
                            egui::TextEdit::multiline(&mut state.notation)
                                .desired_rows(2)
                                .desired_width(f32::INFINITY)
                                .hint_text("e.g. cyclo(Arg-Gly-Asp-D-Phe-NMe-Val) or C10-Trp-... | lactone(13->4) or PEPTIDE1{...}$$$$"),
                        );

                        ui.horizontal_wrapped(|ui| {
                            ui.label("3D Conformation:");
                            for conf in [
                                PeptideConformation::MacrocycleRing,
                                PeptideConformation::AlphaHelix,
                                PeptideConformation::BetaStrand,
                                PeptideConformation::PolyprolineII,
                                PeptideConformation::BetaHairpin,
                                PeptideConformation::Extended,
                            ] {
                                ui.selectable_value(&mut state.conformation, conf, conf.label());
                            }
                        });

                        ui.horizontal_wrapped(|ui| {
                            if ui
                                .button(
                                    RichText::new("Build & Load onto 2D+3D Canvas")
                                        .strong()
                                        .color(Color32::from_rgb(30, 120, 210)),
                                )
                                .clicked()
                                && let Some(mol) = state.build_current()
                            {
                                action = PeptideStudioAction::LoadMolOnCanvas(mol);
                            }
                            if ui.button("Wrap cyclo(...)").clicked()
                                && !state.notation.trim().starts_with("cyclo(")
                            {
                                state.notation = format!("cyclo({})", state.notation.trim());
                            }
                            if ui.button("Reverse-Design NRPS BGC").clicked()
                                && let Ok(spec) =
                                    parse_peptide_spec(&state.notation, &state.catalog, state.conformation)
                            {
                                state.active_bgc = Some(design_nrps_from_peptide(&spec));
                                state.tab = PeptideStudioTab::NrpsBgc;
                                state.status_msg = "Reverse-designed modular NRPS BGC blueprint.".to_string();
                            }
                            if ui.button("Clear").clicked() {
                                state.notation.clear();
                            }
                        });

                        if let Some(ref built) = state.last_built {
                            ui.add_space(4.0);
                            ui.horizontal_wrapped(|ui| {
                                ui.strong("HELM:");
                                ui.monospace(&built.helm);
                            });
                        }

                        ui.add_space(8.0);
                        ui.separator();
                        ui.horizontal_wrapped(|ui| {
                            ui.strong(format!(
                                "Monomer Palette ({} stereo-monomers):",
                                state.catalog.total_stereomonomer_count()
                            ));
                            ui.checkbox(&mut state.insert_as_d, "Insert as D-isomer");
                            ui.checkbox(&mut state.insert_as_nme, "Insert as N-Me");
                        });
                        ui.horizontal_wrapped(|ui| {
                            for cls in [
                                MonomerClass::Canonical,
                                MonomerClass::NrpsCationic,
                                MonomerClass::NrpsAromatic,
                                MonomerClass::NrpsAliphaticSpecial,
                                MonomerClass::BackboneModified,
                                MonomerClass::Depsipeptide,
                                MonomerClass::Heterocycle,
                                MonomerClass::StaplingClick,
                            ] {
                                ui.selectable_value(&mut state.palette_class, cls, cls.label());
                            }
                        });

                        let codes: Vec<(String, String)> = state
                            .catalog
                            .monomers
                            .iter()
                            .filter(|m| m.class == state.palette_class)
                            .map(|m| (m.code.clone(), m.name.clone()))
                            .collect();
                        ui.horizontal_wrapped(|ui| {
                            for (code, name) in codes {
                                if ui.button(&code).on_hover_text(&name).clicked() {
                                    state.append_monomer_token(&code);
                                }
                            }
                        });

                        ui.add_space(4.0);
                        ui.strong("N-Terminal Caps & Fatty-Acyl Starters:");
                        let caps: Vec<(String, String)> = state
                            .catalog
                            .n_caps
                            .iter()
                            .map(|c| (c.code.clone(), c.name.clone()))
                            .collect();
                        ui.horizontal_wrapped(|ui| {
                            for (code, name) in caps {
                                if ui.small_button(&code).on_hover_text(&name).clicked() {
                                    if state.notation.trim().is_empty() {
                                        state.notation = format!("{code}-");
                                    } else {
                                        state.notation = format!("{code}-{}", state.notation.trim());
                                    }
                                }
                            }
                        });
                    });
                }
                PeptideStudioTab::NrpsBgc => {
                    egui::ScrollArea::vertical().show(ui, |ui| {
                        ui.strong("Curated NRPS Gene Clusters");
                        ui.horizontal_wrapped(|ui| {
                            for (id, lbl) in [
                                ("gramicidin_s", "Gramicidin S (GrsA/B)"),
                                ("acv", "ACV Synthetase (PcbAB)"),
                                ("polymyxin_b1", "Polymyxin B1 (PmxA/B/E)"),
                                ("daptomycin", "Daptomycin (DptA/BC/D)"),
                            ] {
                                if ui.button(lbl).clicked() {
                                    state.active_bgc = preset_nrps_bgc(id);
                                    if let Some(mol) = state.compile_active_bgc() {
                                        action = PeptideStudioAction::LoadMolOnCanvas(mol);
                                    }
                                }
                            }
                            if ui.button("Compile Active BGC -> Molecule").clicked()
                                && let Some(mol) = state.compile_active_bgc()
                            {
                                action = PeptideStudioAction::LoadMolOnCanvas(mol);
                            }
                            if let Some(ref rep) = state.last_nrps_report
                                && ui.small_button("Export NRPS Report (.md)").clicked()
                            {
                                let def = PathBuf::from(export_dir)
                                    .join("nrps_bgc_report.md")
                                    .display()
                                    .to_string();
                                if let Some(path) = pick_save_file_native(
                                    "Export NRPS BGC Report",
                                    &def,
                                    "Markdown",
                                    "*.md",
                                ) {
                                    match std::fs::write(&path, render_nrps_markdown(rep)) {
                                        Ok(()) => {
                                            action = PeptideStudioAction::Status(format!(
                                                "Exported NRPS report -> {path}"
                                            ));
                                        }
                                        Err(e) => state.status_msg = format!("Export error: {e}"),
                                    }
                                }
                            }
                        });

                        ui.add_space(6.0);
                        ui.strong("Interactive Stachelhaus 10-aa A-Domain Specificity Lookup");
                        ui.horizontal(|ui| {
                            ui.label("10-aa Code:");
                            ui.add(
                                egui::TextEdit::singleline(&mut state.stachelhaus_query)
                                    .desired_width(140.0),
                            );
                            let pred = predict_a_domain_specificity(&state.stachelhaus_query);
                            let top_s: Vec<String> = pred
                                .candidates
                                .iter()
                                .take(3)
                                .map(|c| {
                                    format!(
                                        "{} ({:.0}% id, p={:.2})",
                                        c.monomer_code, c.identity_pct, c.probability
                                    )
                                })
                                .collect();
                            ui.label(format!("{} -> {}", pred.confidence, top_s.join(" | ")));
                        });

                        ui.add_space(6.0);
                        if let Some(ref bgc) = state.active_bgc.clone() {
                            ui.strong(format!("{} — {}", bgc.name, bgc.organism));
                            ui.small(format!("Release: {}", bgc.release_mode.label()));
                            ui.add_space(4.0);
                            egui::Grid::new("nrps_mod_grid")
                                .striped(true)
                                .min_col_width(60.0)
                                .show(ui, |ui| {
                                    ui.strong("Mod");
                                    ui.strong("Subunit");
                                    ui.strong("Domains");
                                    ui.strong("Stachelhaus");
                                    ui.strong("Monomer");
                                    ui.end_row();
                                    for m in &bgc.modules {
                                        let doms: Vec<String> = m
                                            .domains
                                            .iter()
                                            .map(|d| d.short_label())
                                            .collect();
                                        ui.label(m.module_index.to_string());
                                        ui.label(&m.subunit);
                                        ui.monospace(doms.join("-"));
                                        ui.monospace(
                                            m.stachelhaus_code.as_deref().unwrap_or("-"),
                                        );
                                        ui.strong(m.assigned_monomer.as_deref().unwrap_or("?"));
                                        ui.end_row();
                                    }
                                });

                            ui.add_space(6.0);
                            ui.horizontal_wrapped(|ui| {
                                ui.strong("Domain-Swap Engineering:");
                                ui.label("Module #");
                                ui.add(
                                    egui::DragValue::new(&mut state.swap_module_idx)
                                        .range(1..=bgc.modules.len().max(1)),
                                );
                                ui.label("New Monomer:");
                                ui.add(
                                    egui::TextEdit::singleline(&mut state.swap_new_monomer)
                                        .desired_width(90.0),
                                );
                                if ui.button("Evaluate Swap").clicked() {
                                    match propose_nrps_swap(
                                        bgc,
                                        state.swap_module_idx,
                                        &state.swap_new_monomer,
                                        &state.catalog,
                                    ) {
                                        Ok((ev, eng_built)) => {
                                            state.last_swap_note = format!(
                                                "Mod {} ({} -> {}): {} | ΔMW {:+.2}, ΔWLOGP {:+.2}, ΔbRo5 {:+.3}",
                                                ev.module_index,
                                                ev.old_monomer,
                                                ev.new_monomer,
                                                ev.interface_risk,
                                                ev.delta_mw,
                                                ev.delta_wlogp,
                                                ev.delta_bro5_score
                                            );
                                            action = PeptideStudioAction::LoadMolOnCanvas(eng_built.mol);
                                        }
                                        Err(e) => {
                                            state.last_swap_note = format!("Swap blocked: {e}");
                                        }
                                    }
                                }
                            });
                            if !state.last_swap_note.is_empty() {
                                ui.small(&state.last_swap_note);
                            }
                        }
                    });
                }
                PeptideStudioTab::Developability => {
                    egui::ScrollArea::vertical().show(ui, |ui| {
                        if state.last_dev.is_none() {
                            ui.label("Build a peptide or compile an NRPS BGC first to inspect bRo5 developability, proteolysis, ESI-MS, and SPPS strategy.");
                            return;
                        }
                        if let Some(ref dev) = state.last_dev {
                            ui.horizontal_wrapped(|ui| {
                                ui.strong(format!("{} ({})", dev.name, dev.formula));
                                if ui.small_button("Export Dev Report (.md)").clicked() {
                                    let def = PathBuf::from(export_dir)
                                        .join("peptide_developability.md")
                                        .display()
                                        .to_string();
                                    if let Some(path) = pick_save_file_native(
                                        "Export Developability Report",
                                        &def,
                                        "Markdown",
                                        "*.md",
                                    ) {
                                        let mut full_md = render_developability_markdown(dev);
                                        if let Some(ref spps) = state.last_spps {
                                            full_md.push_str("\n---\n\n");
                                            full_md.push_str(&render_spps_markdown(spps));
                                        }
                                        let _ = std::fs::write(&path, full_md);
                                        action = PeptideStudioAction::Status(format!(
                                            "Exported peptide dossier -> {path}"
                                        ));
                                    }
                                }
                            });

                            ui.label(format!(
                                "bRo5 Score: {:.3} ({}) | Exact Mass: {:.4} Da | pI: {:.2} | Charge(7.4): {:+.2} | Exposed HBD: {}/{} | WLOGP: {:.2} | TPSA: {:.1}",
                                dev.bro5_oral_score,
                                dev.bro5_verdict,
                                dev.monoisotopic_mass,
                                dev.isoelectric_point_pi,
                                dev.charge_ph74,
                                dev.hbd_exposed_est,
                                dev.hbd_total,
                                dev.wlogp,
                                dev.tpsa
                            ));
                            ui.monospace(format!(
                                "ESI-MS: [M+H]+ {:.4} | [M+2H]2+ {:.4} | [M+3H]3+ {:.4} | [M+Na]+ {:.4} | [M-H]- {:.4}",
                                dev.esi_mh_plus,
                                dev.esi_m2h_2plus,
                                dev.esi_m3h_3plus,
                                dev.esi_mna_plus,
                                dev.esi_m_minus_h
                            ));

                            ui.add_space(6.0);
                            ui.strong(format!(
                                "Proteolysis Map ({} Vulnerable / {} Shielded)",
                                dev.active_protease_sites, dev.shielded_protease_sites
                            ));
                            for ps in &dev.protease_sites {
                                let col = if ps.shielded {
                                    Color32::from_rgb(40, 150, 70)
                                } else {
                                    Color32::from_rgb(200, 60, 50)
                                };
                                ui.colored_label(
                                    col,
                                    format!(
                                        "  {} @ Pos {} ({}): {}",
                                        ps.protease,
                                        ps.position,
                                        ps.motif,
                                        ps.shield_reason
                                            .as_deref()
                                            .unwrap_or("VULNERABLE scissile bond")
                                    ),
                                );
                            }

                            if !dev.liabilities.is_empty() {
                                ui.add_space(6.0);
                                ui.strong("Stability & Safety Liabilities");
                                for l in &dev.liabilities {
                                    ui.label(format!(
                                        "  [{}] {} ({}): {}",
                                        l.severity, l.category, l.position_label, l.detail
                                    ));
                                }
                            }
                        }

                        if let Some(ref spps) = state.last_spps {
                            ui.add_space(8.0);
                            ui.separator();
                            ui.strong(format!(
                                "Synthesis Feasibility — Fmoc-SPPS vs. NRPS ({})",
                                spps.spps_difficulty
                            ));
                            ui.label(format!("Recommended Route: {}", spps.recommended_route));
                            ui.label(format!("Resin: {}", spps.recommended_resin));
                            ui.label(format!("Cyclization: {}", spps.cyclization_protocol));
                            ui.small(&spps.nrps_feasibility_note);
                        }
                    });
                }
                PeptideStudioTab::ScanOptimize => {
                    egui::ScrollArea::vertical().show(ui, |ui| {
                        ui.horizontal_wrapped(|ui| {
                            ui.label("Scan Mode:");
                            for m in [
                                ScanMode::Alanine,
                                ScanMode::DAminoAcid,
                                ScanMode::NMethyl,
                                ScanMode::RcmStaple,
                            ] {
                                ui.selectable_value(&mut state.scan_mode, m, m.label());
                            }
                            if ui.button("Run Scan").clicked() {
                                match parse_peptide_spec(
                                    &state.notation,
                                    &state.catalog,
                                    state.conformation,
                                )
                                .and_then(|sp| run_peptide_scan(&sp, state.scan_mode, &state.catalog))
                                {
                                    Ok(bundle) => {
                                        state.status_msg = format!(
                                            "Generated {} analogs via {}",
                                            bundle.variants.len(),
                                            bundle.mode.label()
                                        );
                                        state.last_scan = Some(bundle);
                                    }
                                    Err(e) => state.status_msg = format!("Scan error: {e}"),
                                }
                            }
                        });

                        if let Some(ref bundle) = state.last_scan {
                            ui.horizontal(|ui| {
                                if ui.button("Send All Variants to SAR Explorer").clicked() {
                                    let mut tbl = SarTable::new();
                                    for b in &bundle.built_mols {
                                        tbl.add_molecule(
                                            b.mol.clone(),
                                            None,
                                            Some(bundle.mode.label().to_string()),
                                        );
                                    }
                                    action = PeptideStudioAction::LoadScanIntoSar(tbl);
                                }
                                if ui.small_button("Export CSV").clicked() {
                                    let def = PathBuf::from(export_dir)
                                        .join("peptide_scan.csv")
                                        .display()
                                        .to_string();
                                    if let Some(path) = pick_save_file_native(
                                        "Export Peptide Scan CSV",
                                        &def,
                                        "CSV",
                                        "*.csv",
                                    ) {
                                        let _ = std::fs::write(&path, render_scan_csv(bundle));
                                        action = PeptideStudioAction::Status(format!(
                                            "Exported scan CSV -> {path}"
                                        ));
                                    }
                                }
                            });
                            ui.add_space(4.0);
                            egui::Grid::new("peptide_scan_grid")
                                .striped(true)
                                .min_col_width(58.0)
                                .show(ui, |ui| {
                                    ui.strong("Pos");
                                    ui.strong("Mutation");
                                    ui.strong("Mass");
                                    ui.strong("WLOGP");
                                    ui.strong("HBD_exp");
                                    ui.strong("Protease");
                                    ui.strong("bRo5 (Δ)");
                                    ui.strong("View");
                                    ui.end_row();
                                    for (vi, v) in bundle.variants.iter().enumerate() {
                                        ui.label(&v.position_label);
                                        ui.label(&v.mutation);
                                        ui.label(format!("{:.2}", v.monoisotopic_mass));
                                        ui.label(format!("{:.2}", v.wlogp));
                                        ui.label(v.hbd_exposed.to_string());
                                        ui.label(format!(
                                            "{} ({:+})",
                                            v.active_protease_sites, v.delta_protease_sites
                                        ));
                                        ui.label(format!(
                                            "{:.3} ({:+.3})",
                                            v.bro5_oral_score, v.delta_bro5_score
                                        ));
                                        if ui.small_button("Load").clicked()
                                            && let Some(bm) = bundle.built_mols.get(vi)
                                        {
                                            action =
                                                PeptideStudioAction::LoadMolOnCanvas(bm.mol.clone());
                                        }
                                        ui.end_row();
                                    }
                                });
                        }
                    });
                }
            }
        });

    state.open = open;
    action
}
