#![forbid(unsafe_code)]

//! Peptide & bRo5 Macrocycle Physicochemical, Permeability, Proteolysis, Stability, and ESI-MS /
//! MS-MS Characterization Engine (pure Rust, zero external dependencies).

use std::fmt::Write as _;

use mol_admet_bridge::evaluate_mol_builtin;
use mol_graph::analyze_molecule;
use serde::{Deserialize, Serialize};

use crate::builder::BuiltPeptide;
use crate::monomers::{CCapKind, MonomerCatalog};

const PROTON_MASS: f64 = 1.007_276_466;
const NA_MASS: f64 = 22.989_218;
const WATER_MASS: f64 = 18.010_565;

/// Identified proteolytic cleavage site (or shielded site) along the peptide/NRP chain.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProteaseSite {
    pub protease: String,
    /// 1-based P1 residue index (`0` for N-terminal exopeptidase).
    pub position: usize,
    pub motif: String,
    pub shielded: bool,
    pub shield_reason: Option<String>,
}

/// Chemical degradation or safety liability alert.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PeptideLiability {
    pub category: String,
    pub position_label: String,
    pub severity: String,
    pub detail: String,
}

/// Theoretical MS/MS fragment ion (`b` / `y` series).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MsFragmentIon {
    pub ion_label: String,
    pub mz: f64,
    pub sequence_fragment: String,
}

/// Comprehensive Peptide & bRo5 Macrocycle Developability Report.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PeptideDevelopabilityReport {
    pub name: String,
    pub notation: String,
    pub helm: String,
    pub canonical_smiles: String,
    pub formula: String,
    pub residue_count: usize,
    pub is_cyclic: bool,
    pub largest_ring_size: usize,
    pub monoisotopic_mass: f64,
    pub average_mw: f64,
    pub wlogp: f64,
    pub tpsa: f64,
    pub hbd_total: usize,
    pub hbd_exposed_est: usize,
    pub hba: usize,
    pub rotatable_bonds: usize,
    pub n_methyl_count: usize,
    pub d_residue_count: usize,
    pub charge_ph50: f64,
    pub charge_ph74: f64,
    pub charge_ph80: f64,
    pub isoelectric_point_pi: f64,
    pub gravy: f64,
    pub hydrophobic_moment_alpha: f64,
    pub hydrophobic_moment_beta: f64,
    pub bro5_oral_score: f64,
    pub bro5_verdict: String,
    pub active_protease_sites: usize,
    pub shielded_protease_sites: usize,
    pub protease_sites: Vec<ProteaseSite>,
    pub liabilities: Vec<PeptideLiability>,
    pub esi_mh_plus: f64,
    pub esi_m2h_2plus: f64,
    pub esi_m3h_3plus: f64,
    pub esi_mna_plus: f64,
    pub esi_m_minus_h: f64,
    pub ms_fragments: Vec<MsFragmentIon>,
}

fn net_charge_at_ph(built: &BuiltPeptide, cat: &MonomerCatalog, ph: f64) -> f64 {
    let spec = &built.spec;
    let mut q = 0.0_f64;

    // Free N-terminus (pKa ~ 8.0) when linear and uncapped (`H-`)
    if !spec.is_head_to_tail() && (spec.n_cap.is_empty() || spec.n_cap.eq_ignore_ascii_case("H")) {
        q += 1.0 / (1.0 + 10.0_f64.powf(ph - 8.0));
    }
    // Free C-terminus (pKa ~ 3.2) when linear and FreeAcid
    if !spec.has_c_term_cyclization() && spec.c_cap == CCapKind::FreeAcid {
        q -= 1.0 / (1.0 + 10.0_f64.powf(3.2 - ph));
    }

    // Track sidechain groups consumed by cross-links so they no longer ionize
    let mut consumed_sc: Vec<usize> = Vec::new();
    for cl in &spec.cross_links {
        match *cl {
            crate::notation::CrossLink::Lactam {
                from_acid_res,
                to_amine_res,
                from_c_terminus,
            } => {
                consumed_sc.push(to_amine_res);
                if !from_c_terminus {
                    consumed_sc.push(from_acid_res);
                }
            }
            crate::notation::CrossLink::Disulfide { res_a, res_b } => {
                consumed_sc.push(res_a);
                consumed_sc.push(res_b);
            }
            crate::notation::CrossLink::DeltaAmide { acid_res, .. } => {
                consumed_sc.push(acid_res);
            }
            _ => {}
        }
    }

    for (idx, r) in spec.residues.iter().enumerate() {
        if consumed_sc.contains(&(idx + 1)) {
            continue;
        }
        if let Some(m) = cat.find_monomer(&r.base_code)
            && let Some((pka, is_acid)) = m.sidechain_pka
        {
            if is_acid {
                q -= 1.0 / (1.0 + 10.0_f64.powf(pka - ph));
            } else {
                q += 1.0 / (1.0 + 10.0_f64.powf(ph - pka));
            }
        }
    }
    q
}

fn compute_isoelectric_point(built: &BuiltPeptide, cat: &MonomerCatalog) -> f64 {
    let mut lo = 0.5_f64;
    let mut hi = 13.5_f64;
    for _ in 0..42 {
        let mid = 0.5 * (lo + hi);
        if net_charge_at_ph(built, cat, mid) > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    ((0.5 * (lo + hi)) * 100.0).round() / 100.0
}

fn hydrophobic_moment(hydropathies: &[f64], angle_deg: f64) -> f64 {
    if hydropathies.is_empty() {
        return 0.0;
    }
    let rad = angle_deg.to_radians();
    let mut sx = 0.0_f64;
    let mut sy = 0.0_f64;
    for (i, &h) in hydropathies.iter().enumerate() {
        let a = (i as f64) * rad;
        sx += h * a.cos();
        sy += h * a.sin();
    }
    (((sx * sx + sy * sy).sqrt() / (hydropathies.len() as f64)) * 100.0).round() / 100.0
}

/// Compute comprehensive Peptide / bRo5 Macrocycle Developability, Proteolysis, Stability & MS.
#[must_use]
pub fn analyze_peptide_developability(
    built: &BuiltPeptide,
    cat: &MonomerCatalog,
) -> PeptideDevelopabilityReport {
    let spec = &built.spec;
    let eval = evaluate_mol_builtin(&built.mol, built.mol.name.as_deref());
    let analysis = analyze_molecule(&built.mol);
    let mono_mass = analysis.exact_mass;
    let avg_mw = eval.physicochemical.molecular_weight;
    let wlogp = eval.lipophilicity.wlogp;
    let tpsa = eval.physicochemical.tpsa;
    let hbd_total = eval.physicochemical.h_bond_donors;
    let hba = eval.physicochemical.h_bond_acceptors;
    let rot_bonds = eval.physicochemical.rotatable_bonds;
    let largest_ring_size = built.mol.rings.iter().map(Vec::len).max().unwrap_or(0);

    let n_methyl_count = spec
        .residues
        .iter()
        .filter(|r| r.is_n_methyl || r.base_code.eq_ignore_ascii_case("Sar"))
        .count();
    let d_residue_count = spec.residues.iter().filter(|r| r.is_d).count();

    // Estimate exposed HBDs (chameleonic intramolecular H-bonding in macrocycles >=14 ring atoms shields ~40% of backbone NHs)
    let backbone_nh_count = spec
        .residues
        .iter()
        .filter(|r| {
            !r.is_n_methyl
                && cat
                    .find_monomer(&r.base_code)
                    .is_some_and(|m| !m.is_tertiary_n && !m.is_hydroxy_acid)
        })
        .count();
    let shielded_nh = if spec.is_cyclic() && largest_ring_size >= 14 {
        ((backbone_nh_count as f64) * 0.45).floor() as usize
    } else {
        0
    };
    let hbd_exposed_est = hbd_total.saturating_sub(shielded_nh);

    let charge_ph50 = (net_charge_at_ph(built, cat, 5.0) * 100.0).round() / 100.0;
    let charge_ph74 = (net_charge_at_ph(built, cat, 7.4) * 100.0).round() / 100.0;
    let charge_ph80 = (net_charge_at_ph(built, cat, 8.0) * 100.0).round() / 100.0;
    let isoelectric_point_pi = compute_isoelectric_point(built, cat);

    let hydros: Vec<f64> = spec
        .residues
        .iter()
        .map(|r| cat.find_monomer(&r.base_code).map_or(0.0, |m| m.hydropathy))
        .collect();
    let gravy = if hydros.is_empty() {
        0.0
    } else {
        ((hydros.iter().sum::<f64>() / hydros.len() as f64) * 100.0).round() / 100.0
    };
    let mu_alpha = hydrophobic_moment(&hydros, 100.0);
    let mu_beta = hydrophobic_moment(&hydros, 160.0);

    // Villar / Whitty / Kihlberg bRo5 Oral Macrocycle Score (0.0..=1.0)
    let s_mw = if avg_mw <= 1000.0 {
        1.0
    } else if avg_mw <= 1300.0 {
        1.0 - (avg_mw - 1000.0) / 500.0
    } else {
        0.15
    };
    let s_hbd = if hbd_exposed_est <= 5 {
        1.0
    } else {
        (1.0 - (hbd_exposed_est as f64 - 5.0) * 0.12).clamp(0.05, 1.0)
    };
    let s_tpsa = if tpsa <= 250.0 {
        1.0
    } else {
        (1.0 - (tpsa - 250.0) / 250.0).clamp(0.05, 1.0)
    };
    let s_logp = if (0.0..=7.5).contains(&wlogp) {
        1.0
    } else if wlogp < 0.0 {
        (1.0 + wlogp / 6.0).clamp(0.10, 1.0)
    } else {
        (1.0 - (wlogp - 7.5) / 5.0).clamp(0.10, 1.0)
    };
    let s_cyc = if spec.is_cyclic() && (12..=36).contains(&largest_ring_size) {
        1.0
    } else if spec.is_cyclic() {
        0.75
    } else {
        0.35
    };
    let s_chg = if charge_ph74.abs() <= 1.2 {
        1.0
    } else {
        (1.0 - (charge_ph74.abs() - 1.2) * 0.22).clamp(0.1, 1.0)
    };
    let bro5_oral_score = (((0.20 * s_mw
        + 0.25 * s_hbd
        + 0.18 * s_tpsa
        + 0.15 * s_logp
        + 0.12 * s_cyc
        + 0.10 * s_chg)
        * 1000.0)
        .round())
        / 1000.0;
    let bro5_verdict = if bro5_oral_score >= 0.72 {
        "Oral-Feasible bRo5 (Chameleonic Macrocycle Space)".to_string()
    } else if bro5_oral_score >= 0.48 {
        "Borderline Oral (Permeation Enhancer / NMe-Optimization Recommended)".to_string()
    } else {
        "Parenteral / SC Preferred (High Polarity, Charge, or Linear Liability)".to_string()
    };

    // Protease cleavage audit (Trypsin, Chymotrypsin, Elastase, Pepsin, DPP-4, Exopeptidases)
    let mut protease_sites = Vec::new();
    let n_res = spec.residues.len();

    // N- & C-terminal exopeptidase check
    if !spec.is_head_to_tail() {
        let free_n = spec.n_cap.is_empty() || spec.n_cap.eq_ignore_ascii_case("H");
        let r1 = &spec.residues[0];
        let shield_n = if !free_n {
            Some(format!("N-capped ({})", spec.n_cap))
        } else if r1.is_d {
            Some("D-amino acid at pos 1".to_string())
        } else if r1.is_n_methyl {
            Some("N-methylated pos 1".to_string())
        } else {
            None
        };
        protease_sites.push(ProteaseSite {
            protease: "Aminopeptidase (Plasma)".to_string(),
            position: 1,
            motif: format!("H-{}", r1.display_token()),
            shielded: shield_n.is_some(),
            shield_reason: shield_n,
        });
    }
    if !spec.has_c_term_cyclization() {
        let rl = &spec.residues[n_res - 1];
        let shield_c = if spec.c_cap != CCapKind::FreeAcid {
            Some(format!("C-capped (-{})", spec.c_cap.code()))
        } else if rl.is_d {
            Some("D-amino acid at C-terminus".to_string())
        } else {
            None
        };
        protease_sites.push(ProteaseSite {
            protease: "Carboxypeptidase (Plasma)".to_string(),
            position: n_res,
            motif: format!("{}-OH", rl.display_token()),
            shielded: shield_c.is_some(),
            shield_reason: shield_c,
        });
    }

    // DPP-4 check at N-terminal Xaa1–(Pro/Ala/Aib)2 when uncapped linear
    if !spec.is_head_to_tail()
        && n_res >= 3
        && (spec.n_cap.is_empty() || spec.n_cap.eq_ignore_ascii_case("H"))
    {
        let r1 = &spec.residues[0];
        let r2 = &spec.residues[1];
        let r3 = &spec.residues[2];
        if matches!(r2.base_code.as_str(), "Pro" | "Ala" | "Ser" | "Aib") {
            let reason = if r2.base_code == "Aib" {
                Some("Blocked by α,α-disubstituted Aib2 (Semaglutide strategy)".to_string())
            } else if r2.is_d || r1.is_d {
                Some("Blocked by D-stereochemistry at P1/P2".to_string())
            } else if r3.base_code == "Pro" || r3.is_n_methyl {
                Some("Blocked by Pro/NMe at P1'".to_string())
            } else {
                None
            };
            protease_sites.push(ProteaseSite {
                protease: "DPP-4 (Incretin Peptidase)".to_string(),
                position: 2,
                motif: format!(
                    "{}-{}|{}",
                    r1.display_token(),
                    r2.display_token(),
                    r3.display_token()
                ),
                shielded: reason.is_some(),
                shield_reason: reason,
            });
        }
    }

    for i in 0..n_res.saturating_sub(1) {
        let p1 = &spec.residues[i];
        let p1p = &spec.residues[i + 1];
        let enzymes: &[(&str, bool)] = &[
            ("Trypsin", matches!(p1.base_code.as_str(), "Lys" | "Arg")),
            (
                "Chymotrypsin",
                matches!(p1.base_code.as_str(), "Phe" | "Tyr" | "Trp" | "Leu"),
            ),
            (
                "Elastase",
                matches!(p1.base_code.as_str(), "Ala" | "Val" | "Gly"),
            ),
        ];
        for &(enz, is_target) in enzymes {
            if !is_target {
                continue;
            }
            let mut reasons = Vec::new();
            if p1.is_d {
                reasons.push("D-residue at P1");
            }
            if p1p.is_d {
                reasons.push("D-residue at P1'");
            }
            if p1.is_n_methyl || p1p.is_n_methyl || p1p.base_code == "Sar" {
                reasons.push("N-methylated scissile amide");
            }
            if p1p.base_code == "Pro" || p1p.base_code == "Pip" {
                reasons.push("Pro/Pip at P1'");
            }
            if spec.is_head_to_tail() && largest_ring_size <= 24 {
                reasons.push("Tight macrocycle steric shielding");
            }
            let shielded = !reasons.is_empty();
            protease_sites.push(ProteaseSite {
                protease: enz.to_string(),
                position: i + 1,
                motif: format!("{}|{}", p1.display_token(), p1p.display_token()),
                shielded,
                shield_reason: if shielded {
                    Some(reasons.join("; "))
                } else {
                    None
                },
            });
        }
    }

    let active_protease_sites = protease_sites.iter().filter(|s| !s.shielded).count();
    let shielded_protease_sites = protease_sites.iter().filter(|s| s.shielded).count();

    // Chemical stability & safety liabilities
    let mut liabilities = Vec::new();
    let mut hydro_run = 0usize;
    for (i, r) in spec.residues.iter().enumerate() {
        let pos = i + 1;
        if let Some(next) = spec.residues.get(i + 1) {
            if r.base_code == "Asp" && matches!(next.base_code.as_str(), "Gly" | "Ser" | "Pro") {
                liabilities.push(PeptideLiability {
                    category: "Isomerization / Clipping".to_string(),
                    position_label: format!("{pos}-{}", pos + 1),
                    severity: "Moderate".to_string(),
                    detail: format!(
                        "{}-{} motif prone to succinimide isoAsp formation or acid hydrolysis",
                        r.base_code, next.base_code
                    ),
                });
            }
            if r.base_code == "Asn" && matches!(next.base_code.as_str(), "Gly" | "Ser") {
                liabilities.push(PeptideLiability {
                    category: "Deamidation".to_string(),
                    position_label: format!("{pos}-{}", pos + 1),
                    severity: "Moderate".to_string(),
                    detail: format!(
                        "{}-{} hot-spot for deamidation to Asp/isoAsp",
                        r.base_code, next.base_code
                    ),
                });
            }
        }
        if matches!(r.base_code.as_str(), "Met" | "Trp") {
            liabilities.push(PeptideLiability {
                category: "Oxidation".to_string(),
                position_label: format!("Pos {pos} ({})", r.base_code),
                severity: "Low-Mod".to_string(),
                detail: format!(
                    "{} sidechain sensitive to sulfoxide/kynurenine photo-oxidation",
                    r.base_code
                ),
            });
        }
        if i == 0
            && r.base_code == "Gln"
            && !spec.is_head_to_tail()
            && spec.n_cap.eq_ignore_ascii_case("H")
        {
            liabilities.push(PeptideLiability {
                category: "Pyroglutamate".to_string(),
                position_label: "Pos 1 (Gln)".to_string(),
                severity: "Moderate".to_string(),
                detail: "Uncapped N-terminal Gln cyclizes spontaneously to pyroglutamate (-17 Da)"
                    .to_string(),
            });
        }
        if cat
            .find_monomer(&r.base_code)
            .is_some_and(|m| m.is_hydrophobic)
        {
            hydro_run += 1;
            if hydro_run == 5 {
                liabilities.push(PeptideLiability {
                    category: "Aggregation".to_string(),
                    position_label: format!("Pos {}–{pos}", pos - 4),
                    severity: "High".to_string(),
                    detail: ">=5 consecutive hydrophobic residues: high β-sheet aggregation & SPPS coupling difficulty".to_string(),
                });
            }
        } else {
            hydro_run = 0;
        }
    }
    if charge_ph74 >= 3.5 && (mu_alpha >= 0.45 || mu_beta >= 0.45) {
        liabilities.push(PeptideLiability {
            category: "Hemolysis / Membrane Lysis".to_string(),
            position_label: "Global".to_string(),
            severity: "High".to_string(),
            detail: format!(
                "Strongly cationic amphipathic profile (charge {:+.1}, μH_β={:.2}): monitor mammalian RBC hemolysis & nephrotoxicity",
                charge_ph74, mu_beta
            ),
        });
    }

    // Theoretical ESI MS & linear b/y fragment ions
    let esi_mh_plus = ((mono_mass + PROTON_MASS) * 10_000.0).round() / 10_000.0;
    let esi_m2h_2plus = (((mono_mass + 2.0 * PROTON_MASS) / 2.0) * 10_000.0).round() / 10_000.0;
    let esi_m3h_3plus = (((mono_mass + 3.0 * PROTON_MASS) / 3.0) * 10_000.0).round() / 10_000.0;
    let esi_mna_plus = ((mono_mass + NA_MASS) * 10_000.0).round() / 10_000.0;
    let esi_m_minus_h = ((mono_mass - PROTON_MASS) * 10_000.0).round() / 10_000.0;

    let mut ms_fragments = Vec::new();
    if !spec.is_head_to_tail() && spec.residues.len() >= 2 && spec.residues.len() <= 20 {
        let res_masses: Vec<f64> = spec
            .residues
            .iter()
            .map(|r| {
                let mut m_val = cat
                    .find_monomer(&r.base_code)
                    .and_then(|ms| mol_graph::parse_smiles(&ms.smiles_template).ok())
                    .map_or(110.0, |mut g| {
                        for a in &mut g.atoms {
                            if let Some(iso) = a.isotope
                                && (101..=109).contains(&iso)
                            {
                                a.isotope = None;
                                if !a.aromatic && a.formal_charge == 0 {
                                    a.from_bracket = false;
                                    a.explicit_h = 0;
                                }
                            }
                        }
                        let _ = g.perceive();
                        analyze_molecule(&g).exact_mass - WATER_MASS
                    });
                if r.is_n_methyl {
                    m_val += 14.015_650;
                }
                m_val
            })
            .collect();
        let cap_mass = if spec.n_cap.is_empty() || spec.n_cap.eq_ignore_ascii_case("H") {
            0.0
        } else {
            cat.find_n_cap(&spec.n_cap)
                .filter(|c| !c.smiles_template.is_empty())
                .and_then(|c| mol_graph::parse_smiles(&c.smiles_template).ok())
                .map_or(0.0, |mut g| {
                    for a in &mut g.atoms {
                        if let Some(iso) = a.isotope
                            && (101..=109).contains(&iso)
                        {
                            a.isotope = None;
                            a.from_bracket = false;
                            a.explicit_h = 0;
                        }
                    }
                    let _ = g.perceive();
                    analyze_molecule(&g).exact_mass - WATER_MASS
                })
        };
        let mut cum_b = cap_mass;
        for k in 1..res_masses.len() {
            cum_b += res_masses[k - 1];
            let b_mz = ((cum_b + PROTON_MASS) * 1000.0).round() / 1000.0;
            let y_mz = (((mono_mass - cum_b) + PROTON_MASS) * 1000.0).round() / 1000.0;
            let b_seq: Vec<String> = spec.residues[..k]
                .iter()
                .map(|r| r.display_token())
                .collect();
            let y_seq: Vec<String> = spec.residues[k..]
                .iter()
                .map(|r| r.display_token())
                .collect();
            ms_fragments.push(MsFragmentIon {
                ion_label: format!("b{k}+"),
                mz: b_mz,
                sequence_fragment: b_seq.join("-"),
            });
            ms_fragments.push(MsFragmentIon {
                ion_label: format!("y{}+", res_masses.len() - k),
                mz: y_mz,
                sequence_fragment: y_seq.join("-"),
            });
        }
    }

    PeptideDevelopabilityReport {
        name: built
            .mol
            .name
            .clone()
            .unwrap_or_else(|| built.notation.clone()),
        notation: built.notation.clone(),
        helm: built.helm.clone(),
        canonical_smiles: built.canonical_smiles.clone(),
        formula: analysis.formula_with_charge,
        residue_count: spec.residues.len(),
        is_cyclic: spec.is_cyclic(),
        largest_ring_size,
        monoisotopic_mass: (mono_mass * 10_000.0).round() / 10_000.0,
        average_mw: (avg_mw * 100.0).round() / 100.0,
        wlogp: (wlogp * 100.0).round() / 100.0,
        tpsa: (tpsa * 100.0).round() / 100.0,
        hbd_total,
        hbd_exposed_est,
        hba,
        rotatable_bonds: rot_bonds,
        n_methyl_count,
        d_residue_count,
        charge_ph50,
        charge_ph74,
        charge_ph80,
        isoelectric_point_pi,
        gravy,
        hydrophobic_moment_alpha: mu_alpha,
        hydrophobic_moment_beta: mu_beta,
        bro5_oral_score,
        bro5_verdict,
        active_protease_sites,
        shielded_protease_sites,
        protease_sites,
        liabilities,
        esi_mh_plus,
        esi_m2h_2plus,
        esi_m3h_3plus,
        esi_mna_plus,
        esi_m_minus_h,
        ms_fragments,
    }
}

/// Render a publication-ready Markdown Peptide / NRP Developability & MS Dossier.
#[must_use]
pub fn render_developability_markdown(r: &PeptideDevelopabilityReport) -> String {
    let mut md = String::new();
    let _ = writeln!(
        md,
        "# Peptide / NRP Developability & Analytical Report: {}\n",
        r.name
    );
    let _ = writeln!(md, "- **Notation:** `{}`", r.notation);
    let _ = writeln!(md, "- **HELM:** `{}`", r.helm);
    let _ = writeln!(
        md,
        "- **Formula:** `{}` | **Exact Mass:** `{:.4} Da` | **Avg MW:** `{:.2} Da`",
        r.formula, r.monoisotopic_mass, r.average_mw
    );
    let _ = writeln!(
        md,
        "- **bRo5 Oral Macrocycle Score:** **`{:.3}`** — {}",
        r.bro5_oral_score, r.bro5_verdict
    );

    let _ = writeln!(md, "\n## 1. bRo5 Physicochemical & Permeability Profile\n");
    let _ = writeln!(md, "| Descriptor | Value | Descriptor | Value |");
    let _ = writeln!(md, "|---|---:|---|---:|");
    let _ = writeln!(
        md,
        "| Residues (D / NMe) | {} ({} D / {} NMe) | Largest Ring | {} atoms |",
        r.residue_count, r.d_residue_count, r.n_methyl_count, r.largest_ring_size
    );
    let _ = writeln!(
        md,
        "| WLOGP | {:.2} | TPSA (Ertl) | {:.1} Å² |",
        r.wlogp, r.tpsa
    );
    let _ = writeln!(
        md,
        "| HBD Total / Exposed | {} / {} | HBA / RotBonds | {} / {} |",
        r.hbd_total, r.hbd_exposed_est, r.hba, r.rotatable_bonds
    );
    let _ = writeln!(
        md,
        "| Net Charge (pH 5.0 / 7.4 / 8.0) | {:+.2} / {:+.2} / {:+.2} | Isoelectric Point (`pI`) | {:.2} |",
        r.charge_ph50, r.charge_ph74, r.charge_ph80, r.isoelectric_point_pi
    );
    let _ = writeln!(
        md,
        "| GRAVY Hydropathy | {:.2} | Hydrophobic Moment (α / β) | {:.2} / {:.2} |",
        r.gravy, r.hydrophobic_moment_alpha, r.hydrophobic_moment_beta
    );

    let _ = writeln!(md, "\n## 2. ESI-MS Analytical Reference Ions\n");
    let _ = writeln!(
        md,
        "`[M+H]+ = {:.4}` | `[M+2H]2+ = {:.4}` | `[M+3H]3+ = {:.4}` | `[M+Na]+ = {:.4}` | `[M-H]- = {:.4}`",
        r.esi_mh_plus, r.esi_m2h_2plus, r.esi_m3h_3plus, r.esi_mna_plus, r.esi_m_minus_h
    );

    let _ = writeln!(
        md,
        "\n## 3. Proteolysis Audit ({} Active / {} Shielded Sites)\n",
        r.active_protease_sites, r.shielded_protease_sites
    );
    if r.protease_sites.is_empty() {
        let _ = writeln!(md, "- No canonical protease cleavage sites detected.");
    } else {
        let _ = writeln!(md, "| Protease | Pos | Motif | Status | Mechanism |");
        let _ = writeln!(md, "|---|---:|---|---|---|");
        for ps in &r.protease_sites {
            let _ = writeln!(
                md,
                "| {} | {} | `{}` | {} | {} |",
                ps.protease,
                ps.position,
                ps.motif,
                if ps.shielded {
                    "Shielded"
                } else {
                    "**Vulnerable**"
                },
                ps.shield_reason
                    .as_deref()
                    .unwrap_or("Exposed scissile amide")
            );
        }
    }

    if !r.liabilities.is_empty() {
        let _ = writeln!(md, "\n## 4. Stability & Safety Alerts\n");
        for l in &r.liabilities {
            let _ = writeln!(
                md,
                "- **[{}] {} ({}):** {}",
                l.severity, l.category, l.position_label, l.detail
            );
        }
    }
    md
}
