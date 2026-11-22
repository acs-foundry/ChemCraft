#![forbid(unsafe_code)]

//! Fmoc-SPPS vs. NRPS Biosynthetic Feasibility Planner and Combinatorial Peptidomimetic
//! Optimization Scans (`Alanine`, `D-Amino Acid`, `N-Methyl`, and `RCM Staple` scans).

use std::fmt::Write as _;

use serde::{Deserialize, Serialize};

use crate::builder::{BuiltPeptide, build_peptide_from_spec};
use crate::developability::{PeptideDevelopabilityReport, analyze_peptide_developability};
use crate::monomers::{CCapKind, MonomerCatalog};
use crate::notation::{CrossLink, PeptideError, PeptideSpec};

/// Difficult coupling site flagged during Fmoc-SPPS synthesis planning.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SppsCouplingAlert {
    pub junction: String,
    pub reason: String,
    pub recommended_reagent: String,
}

/// Structured Fmoc-SPPS vs. NRPS biosynthetic feasibility plan.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SppsPlan {
    pub target_name: String,
    pub spps_difficulty: String,
    pub difficulty_score: u32,
    pub recommended_route: String,
    pub recommended_resin: String,
    pub protecting_group_strategy: Vec<String>,
    pub coupling_alerts: Vec<SppsCouplingAlert>,
    pub cyclization_protocol: String,
    pub nrps_feasibility_note: String,
}

/// Plan Fmoc-SPPS synthesis strategy and compare directly against NRPS biosynthesis.
#[must_use]
pub fn plan_spps_vs_nrps(spec: &PeptideSpec, cat: &MonomerCatalog) -> SppsPlan {
    let mut score = 0u32;
    let mut alerts = Vec::new();
    let mut pgs = vec!["Standard Fmoc/tBu orthogonal protection (Nα-Fmoc; Boc for Lys/Orn/Dab/Trp; tBu for Asp/Glu/Ser/Thr/Tyr; Pbf for Arg)".to_string()];

    let mut has_exotic_nrps = false;
    let mut has_synthetic_only = false;
    let mut hydro_run = 0usize;

    for (i, r) in spec.residues.iter().enumerate() {
        let m = cat.find_monomer(&r.base_code);
        if matches!(
            r.base_code.as_str(),
            "Bmt" | "Kyn" | "mGlu" | "Hpg" | "Dpg" | "Bht"
        ) {
            has_exotic_nrps = true;
            score += 3;
        }
        if matches!(
            r.base_code.as_str(),
            "S5" | "R8" | "Pra" | "Aha" | "Aib" | "Ac5c"
        ) {
            has_synthetic_only = true;
        }
        if m.is_some_and(|ms| ms.is_hydrophobic) {
            hydro_run += 1;
            if hydro_run >= 5 {
                score += 2;
            }
        } else {
            hydro_run = 0;
        }

        if let Some(next) = spec.residues.get(i + 1) {
            let cur_bulky = m.is_some_and(|ms| ms.bulky_coupling) || r.is_n_methyl;
            let next_nme = next.is_n_methyl
                || next.base_code.eq_ignore_ascii_case("Sar")
                || next.base_code == "Pro";
            if next_nme && cur_bulky {
                score += 3;
                alerts.push(SppsCouplingAlert {
                    junction: format!("{}-{} -> {}-{}", i + 1, r.display_token(), i + 2, next.display_token()),
                    reason: "Sterically congested coupling onto N-methyl / secondary amine from bulky/β-branched residue".to_string(),
                    recommended_reagent: "HATU / HOAt / DIPEA (2× double coupling, 50°C) or Triphosgene (BTC) / Collidine".to_string(),
                });
            } else if next_nme || r.base_code == "Aib" || next.base_code == "Aib" {
                score += 1;
                alerts.push(SppsCouplingAlert {
                    junction: format!(
                        "{}-{} -> {}-{}",
                        i + 1,
                        r.display_token(),
                        i + 2,
                        next.display_token()
                    ),
                    reason: "Hindered N-Me or α,α-disubstituted coupling".to_string(),
                    recommended_reagent: "HATU / OxymaPure (double coupling, 2× 45 min)"
                        .to_string(),
                });
            }
            if m.is_some_and(|ms| ms.is_hydroxy_acid) {
                score += 2;
                alerts.push(SppsCouplingAlert {
                    junction: format!(
                        "{}-{} -> {}-{}",
                        i + 1,
                        r.display_token(),
                        i + 2,
                        next.display_token()
                    ),
                    reason: "Depsipeptide ester bond formation".to_string(),
                    recommended_reagent:
                        "DIC / DMAP (0.1 eq) or Yamaguchi esterification in DCM/DMF".to_string(),
                });
            }
        }
    }

    let recommended_resin = if spec.is_head_to_tail() || spec.has_c_term_cyclization() {
        "2-Chlorotrityl Chloride (2-CTC) resin — mild 1% TFA/DCM cleavage retains sidechain PGs for solution-phase macrocyclization"
            .to_string()
    } else {
        match spec.c_cap {
            CCapKind::Amide => {
                "Rink Amide MBHA / AM resin (releases C-terminal -NH2 upon 95% TFA cleavage)"
                    .to_string()
            }
            CCapKind::Alcohol => {
                "Pre-loaded amino-alcohol 2-CTC resin (e.g. H-Thr(tBu)-ol-2-CTC)".to_string()
            }
            CCapKind::Aldehyde => {
                "Weinreb amide resin + LiAlH4 / DIBAL-H reductive cleavage".to_string()
            }
            _ => "Wang or 2-CTC resin (releases C-terminal -COOH)".to_string(),
        }
    };

    let mut cyc_steps = Vec::new();
    for cl in &spec.cross_links {
        match *cl {
            CrossLink::HeadToTail => {
                score += 2;
                cyc_steps.push("High-dilution (0.5–1.0 mM in DMF/DCM) head-to-tail macrolactamization with PyAOP/HOAt/DIPEA prior to global TFA deprotection".to_string());
            }
            CrossLink::Macrolactone {
                from_res,
                to_oh_res,
            } => {
                score += 3;
                pgs.push(format!("Leave Mod {to_oh_res} sidechain -OH unprotected (or TBDMS-protected) for selective macrolactonization"));
                cyc_steps.push(format!("Solution-phase Yamaguchi / MNBA macrolactonization ({from_res}->{to_oh_res}) at 0.5 mM to suppress oligomerization"));
            }
            CrossLink::Lactam {
                from_acid_res,
                to_amine_res,
                from_c_terminus,
            } => {
                score += 2;
                pgs.push(format!(
                    "Orthogonal Alloc/ivDde on amine Pos {to_amine_res}{}",
                    if from_c_terminus {
                        ""
                    } else {
                        " + OAll on acid"
                    }
                ));
                cyc_steps.push(format!("Selective deprotection (Pd(PPh3)4 / 2% N2H4) + PyAOP/DIPEA lactamization ({from_acid_res}->{to_amine_res})"));
            }
            CrossLink::Disulfide { res_a, res_b } => {
                score += 1;
                pgs.push(format!(
                    "Cys(Trt) or orthogonal Cys(Acm) pair at {res_a}-{res_b}"
                ));
                cyc_steps.push(format!("Dilute oxidative folding (pH 7.8 NH4HCO3 + 10% DMSO or I2/AcOH for Acm) for Cys{res_a}–Cys{res_b}"));
            }
            CrossLink::RcmStaple { res_a, res_b } => {
                score += 2;
                cyc_steps.push(format!("On-resin Grubbs 1st-gen Ru-catalyst (20 mol% in DCE, 2× 2 h) for RCM staple ({res_a}-{res_b})"));
            }
            CrossLink::TriazoleClick {
                alkyne_res,
                azide_res,
            } => {
                score += 1;
                cyc_steps.push(format!("On-resin CuI / sodium ascorbate / DIPEA CuAAC click cyclization ({alkyne_res}-{azide_res})"));
            }
            CrossLink::Lanthionine {
                thiol_res,
                dehydro_res,
            } => {
                score += 2;
                cyc_steps.push(format!("Base-catalyzed intramolecular Michael addition Cys{thiol_res} -> Dha/Dhb{dehydro_res}"));
            }
            CrossLink::DeltaAmide { acid_res, next_res } => {
                pgs.push(format!("Use Fmoc-L-Aad(OtBu)-OH (α-OtBu protected) to couple selectively via δ-COOH ({acid_res}->{next_res})"));
            }
        }
    }
    let cyclization_protocol = if cyc_steps.is_empty() {
        "Linear cleavage in TFA/TIS/H2O (95:2.5:2.5, 2 h) + cold Et2O precipitation & RP-HPLC (C18)"
            .to_string()
    } else {
        cyc_steps.join("; ")
    };

    let spps_difficulty = match score {
        0..=2 => "Easy (Routine automated Fmoc-SPPS)".to_string(),
        3..=5 => "Moderate (Double coupling / orthogonal PG needed)".to_string(),
        6..=9 => "Challenging (Steric hindrance / complex macrocyclization)".to_string(),
        _ => "High Complexity (Specialized fragment condensation or NRPS fermentation)".to_string(),
    };

    let recommended_route = if has_synthetic_only {
        "Fmoc-SPPS (Required: contains synthetic non-biological stapling/click or Aib monomers)"
            .to_string()
    } else if has_exotic_nrps && spec.residues.len() >= 10 {
        "NRPS Fermentation + Semi-Synthesis (Complex non-proteinogenic monomers & macrocycle cheaper via engineered chassis; Fmoc-SPPS for SAR analogs)".to_string()
    } else {
        "Fmoc-SPPS for MedChem SAR (<100 mg) / Optional NRPS engineering for scale-up".to_string()
    };

    let nrps_feasibility_note = if has_synthetic_only {
        "Non-canonical synthetic monomers (S5/R8/Pra/Aha/Aib) are outside standard bacterial NRPS A-domains; keep in SPPS or use flexizyme/RaPID mRNA display.".to_string()
    } else if has_exotic_nrps {
        "Native NRPS biosynthetic precedents exist (e.g. SimA, DptA-D, PmxA-E, CepA-C); ideal candidate for BGC heterologous expression or XU domain-swap mutagenesis.".to_string()
    } else {
        "Fully expressible via engineered NRPS assembly line (A/T/C/E/nMT/TE modules) or standard Fmoc-SPPS.".to_string()
    };

    SppsPlan {
        target_name: spec.name.clone().unwrap_or_else(|| spec.to_notation()),
        spps_difficulty,
        difficulty_score: score,
        recommended_route,
        recommended_resin,
        protecting_group_strategy: pgs,
        coupling_alerts: alerts,
        cyclization_protocol,
        nrps_feasibility_note,
    }
}

/// Render Markdown ELN protocol for Fmoc-SPPS vs NRPS feasibility.
#[must_use]
pub fn render_spps_markdown(p: &SppsPlan) -> String {
    let mut md = String::new();
    let _ = writeln!(
        md,
        "# Synthetic Feasibility Plan (Fmoc-SPPS vs. NRPS): {}\n",
        p.target_name
    );
    let _ = writeln!(
        md,
        "- **SPPS Difficulty:** {} (score {})",
        p.spps_difficulty, p.difficulty_score
    );
    let _ = writeln!(md, "- **Recommended Route:** {}", p.recommended_route);
    let _ = writeln!(md, "- **Resin:** {}", p.recommended_resin);
    let _ = writeln!(
        md,
        "- **Cyclization / Cleavage:** {}",
        p.cyclization_protocol
    );
    let _ = writeln!(
        md,
        "- **NRPS Biosynthesis Note:** {}\n",
        p.nrps_feasibility_note
    );
    let _ = writeln!(md, "## Protecting-Group Strategy");
    for pg in &p.protecting_group_strategy {
        let _ = writeln!(md, "- {pg}");
    }
    if !p.coupling_alerts.is_empty() {
        let _ = writeln!(md, "\n## Difficult Coupling Junctions\n");
        let _ = writeln!(md, "| Junction | Liability | Recommended Activation |");
        let _ = writeln!(md, "|---|---|---|");
        for a in &p.coupling_alerts {
            let _ = writeln!(
                md,
                "| `{}` | {} | {} |",
                a.junction, a.reason, a.recommended_reagent
            );
        }
    }
    md
}

/// Combinatorial peptide optimization scan mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ScanMode {
    Alanine,
    DAminoAcid,
    NMethyl,
    RcmStaple,
}

impl ScanMode {
    pub fn label(self) -> &'static str {
        match self {
            Self::Alanine => "Alanine Scan (SAR Hotspots)",
            Self::DAminoAcid => "D-Amino Acid Scan (Protease Shielding)",
            Self::NMethyl => "N-Methyl Scan (bRo5 Permeability)",
            Self::RcmStaple => "i,i+4 RCM Staple Scan (Helix Stabilization)",
        }
    }

    pub fn from_str_lossy(s: &str) -> Option<Self> {
        match s.trim().to_ascii_lowercase().as_str() {
            "ala" | "alanine" | "ala-scan" => Some(Self::Alanine),
            "d" | "d-aa" | "d-scan" | "daa" | "epimer" => Some(Self::DAminoAcid),
            "nme" | "n-me" | "nmethyl" | "n-methyl" | "nme-scan" => Some(Self::NMethyl),
            "staple" | "rcm" | "staple-scan" => Some(Self::RcmStaple),
            _ => None,
        }
    }
}

/// Single analog produced by a peptide optimization scan.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ScanVariant {
    pub label: String,
    pub position_label: String,
    pub mutation: String,
    pub notation: String,
    pub helm: String,
    pub canonical_smiles: String,
    pub monoisotopic_mass: f64,
    pub wlogp: f64,
    pub tpsa: f64,
    pub hbd_exposed: usize,
    pub active_protease_sites: usize,
    pub bro5_oral_score: f64,
    pub delta_protease_sites: i32,
    pub delta_bro5_score: f64,
}

/// Complete combinatorial scan report + built 2D/3D `BuiltPeptide` objects for SDF/SAR export.
#[derive(Debug, Clone)]
pub struct PeptideScanBundle {
    pub mode: ScanMode,
    pub parent: PeptideDevelopabilityReport,
    pub variants: Vec<ScanVariant>,
    pub built_mols: Vec<BuiltPeptide>,
}

/// Run a systematic `Alanine`, `D-Amino Acid`, `N-Methyl`, or `i,i+4 RCM Staple` optimization scan.
pub fn run_peptide_scan(
    parent_spec: &PeptideSpec,
    mode: ScanMode,
    cat: &MonomerCatalog,
) -> Result<PeptideScanBundle, PeptideError> {
    let parent_built = build_peptide_from_spec(parent_spec, cat)?;
    let parent_dev = analyze_peptide_developability(&parent_built, cat);

    // Avoid mutating positions that serve as covalent cross-link anchors
    let mut locked_pos: Vec<usize> = Vec::new();
    for cl in &parent_spec.cross_links {
        match *cl {
            CrossLink::HeadToTail => {}
            CrossLink::Macrolactone {
                from_res,
                to_oh_res,
            } => {
                locked_pos.push(from_res);
                locked_pos.push(to_oh_res);
            }
            CrossLink::Lactam {
                from_acid_res,
                to_amine_res,
                ..
            } => {
                locked_pos.push(from_acid_res);
                locked_pos.push(to_amine_res);
            }
            CrossLink::Disulfide { res_a, res_b } | CrossLink::RcmStaple { res_a, res_b } => {
                locked_pos.push(res_a);
                locked_pos.push(res_b);
            }
            CrossLink::Lanthionine {
                thiol_res,
                dehydro_res,
            } => {
                locked_pos.push(thiol_res);
                locked_pos.push(dehydro_res);
            }
            CrossLink::TriazoleClick {
                alkyne_res,
                azide_res,
            } => {
                locked_pos.push(alkyne_res);
                locked_pos.push(azide_res);
            }
            CrossLink::DeltaAmide { acid_res, next_res } => {
                locked_pos.push(acid_res);
                locked_pos.push(next_res);
            }
        }
    }

    let mut variants = Vec::new();
    let mut built_mols = Vec::new();
    let n = parent_spec.residues.len();

    match mode {
        ScanMode::Alanine | ScanMode::DAminoAcid | ScanMode::NMethyl => {
            for i in 0..n {
                let pos = i + 1;
                if mode == ScanMode::Alanine && locked_pos.contains(&pos) {
                    continue;
                }
                let orig = &parent_spec.residues[i];
                let Some(mspec) = cat.find_monomer(&orig.base_code) else {
                    continue;
                };
                let mut mutated = parent_spec.clone();
                let mut_desc = match mode {
                    ScanMode::Alanine => {
                        if orig.base_code == "Ala" && !orig.is_d {
                            continue;
                        }
                        mutated.residues[i].base_code = "Ala".to_string();
                        mutated.residues[i].is_d = false;
                        format!("{} -> Ala", orig.display_token())
                    }
                    ScanMode::DAminoAcid => {
                        if !mspec.has_alpha_chiral || orig.is_d {
                            continue;
                        }
                        mutated.residues[i].is_d = true;
                        format!(
                            "{} -> {}",
                            orig.display_token(),
                            mutated.residues[i].display_token()
                        )
                    }
                    ScanMode::NMethyl => {
                        if orig.is_n_methyl || mspec.is_tertiary_n || orig.base_code == "Sar" {
                            continue;
                        }
                        if orig.base_code == "Gly" {
                            mutated.residues[i].base_code = "Sar".to_string();
                        } else {
                            mutated.residues[i].is_n_methyl = true;
                        }
                        format!(
                            "{} -> {}",
                            orig.display_token(),
                            mutated.residues[i].display_token()
                        )
                    }
                    ScanMode::RcmStaple => unreachable!(),
                };
                let lbl = format!(
                    "{}_pos{pos}_{}",
                    parent_spec
                        .name
                        .as_deref()
                        .unwrap_or("Peptide")
                        .replace(' ', "_"),
                    mutated.residues[i].display_token()
                );
                mutated.name = Some(lbl.clone());
                if let Ok(b) = build_peptide_from_spec(&mutated, cat) {
                    let dev = analyze_peptide_developability(&b, cat);
                    variants.push(ScanVariant {
                        label: lbl,
                        position_label: format!("Pos {pos}"),
                        mutation: mut_desc,
                        notation: b.notation.clone(),
                        helm: b.helm.clone(),
                        canonical_smiles: b.canonical_smiles.clone(),
                        monoisotopic_mass: dev.monoisotopic_mass,
                        wlogp: dev.wlogp,
                        tpsa: dev.tpsa,
                        hbd_exposed: dev.hbd_exposed_est,
                        active_protease_sites: dev.active_protease_sites,
                        bro5_oral_score: dev.bro5_oral_score,
                        delta_protease_sites: dev.active_protease_sites as i32
                            - parent_dev.active_protease_sites as i32,
                        delta_bro5_score: ((dev.bro5_oral_score - parent_dev.bro5_oral_score)
                            * 1000.0)
                            .round()
                            / 1000.0,
                    });
                    built_mols.push(b);
                }
            }
        }
        ScanMode::RcmStaple => {
            for i in 0..n.saturating_sub(4) {
                let p1 = i + 1;
                let p2 = i + 5;
                if locked_pos.contains(&p1) || locked_pos.contains(&p2) {
                    continue;
                }
                let mut mutated = parent_spec.clone();
                mutated.residues[i] = crate::notation::ResidueUnit {
                    base_code: "S5".to_string(),
                    is_d: false,
                    is_n_methyl: false,
                };
                mutated.residues[i + 4] = crate::notation::ResidueUnit {
                    base_code: "S5".to_string(),
                    is_d: false,
                    is_n_methyl: false,
                };
                mutated.cross_links.push(CrossLink::RcmStaple {
                    res_a: p1,
                    res_b: p2,
                });
                mutated.conformation = crate::notation::PeptideConformation::AlphaHelix;
                let lbl = format!(
                    "{}_staple_{p1}_{p2}",
                    parent_spec
                        .name
                        .as_deref()
                        .unwrap_or("Peptide")
                        .replace(' ', "_")
                );
                mutated.name = Some(lbl.clone());
                if let Ok(b) = build_peptide_from_spec(&mutated, cat) {
                    let dev = analyze_peptide_developability(&b, cat);
                    variants.push(ScanVariant {
                        label: lbl,
                        position_label: format!("Pos {p1},{p2}"),
                        mutation: format!("S5({p1})–S5({p2}) RCM i,i+4"),
                        notation: b.notation.clone(),
                        helm: b.helm.clone(),
                        canonical_smiles: b.canonical_smiles.clone(),
                        monoisotopic_mass: dev.monoisotopic_mass,
                        wlogp: dev.wlogp,
                        tpsa: dev.tpsa,
                        hbd_exposed: dev.hbd_exposed_est,
                        active_protease_sites: dev.active_protease_sites,
                        bro5_oral_score: dev.bro5_oral_score,
                        delta_protease_sites: dev.active_protease_sites as i32
                            - parent_dev.active_protease_sites as i32,
                        delta_bro5_score: ((dev.bro5_oral_score - parent_dev.bro5_oral_score)
                            * 1000.0)
                            .round()
                            / 1000.0,
                    });
                    built_mols.push(b);
                }
            }
        }
    }

    Ok(PeptideScanBundle {
        mode,
        parent: parent_dev,
        variants,
        built_mols,
    })
}

/// Export a peptide scan bundle as RFC-4180 CSV.
#[must_use]
pub fn render_scan_csv(bundle: &PeptideScanBundle) -> String {
    let mut csv = String::from(
        "Label,Position,Mutation,Notation,HELM,SMILES,ExactMass,WLOGP,TPSA,ExposedHBD,ActiveProteaseSites,DeltaProtease,bRo5Score,Delta_bRo5\n",
    );
    for v in &bundle.variants {
        let _ = writeln!(
            csv,
            "\"{}\",\"{}\",\"{}\",\"{}\",\"{}\",\"{}\",{:.4},{:.2},{:.2},{},{},{:+},{:.3},{:+.3}",
            v.label.replace('"', "\"\""),
            v.position_label,
            v.mutation.replace('"', "\"\""),
            v.notation.replace('"', "\"\""),
            v.helm.replace('"', "\"\""),
            v.canonical_smiles.replace('"', "\"\""),
            v.monoisotopic_mass,
            v.wlogp,
            v.tpsa,
            v.hbd_exposed,
            v.active_protease_sites,
            v.delta_protease_sites,
            v.bro5_oral_score,
            v.delta_bro5_score
        );
    }
    csv
}
