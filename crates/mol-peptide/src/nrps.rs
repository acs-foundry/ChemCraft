#![forbid(unsafe_code)]

//! Authentic Non-Ribosomal Peptide Synthetase (NRPS) Biosynthesis Engine:
//! - Modular assembly-line representation (`C_starter`, `C_LCL`, `C_DCL`, `Cy`, `A[Stachelhaus]`,
//!   `T`, `E`, `nMT`, `Ox`, `Red`, `TE`, post-NRPS tailoring)
//! - Probabilistic 10-aa Stachelhaus A-domain specificity predictor (top-3 candidates + uncertainty)
//! - Collinearity & biosynthetic interface compatibility audit (`E -> C_DCL`, `Cy -> Cys/Ser/Thr`, etc.)
//! - Lightweight pure-Rust antiSMASH GenBank `.gbk` parser (`aSDomain` / `NRPS_PKS` / `CDS_motif`)
//! - Forward BGC -> NRP compilation, reverse Peptide -> BGC blueprinting, and Domain-Swap evaluation

use std::fmt::Write as _;

use serde::{Deserialize, Serialize};

use crate::builder::{BuiltPeptide, build_peptide_from_spec};
use crate::developability::analyze_peptide_developability;
use crate::monomers::{CCapKind, MonomerCatalog};
use crate::notation::{CrossLink, PeptideConformation, PeptideError, PeptideSpec, ResidueUnit};

/// Canonical NRPS enzymatic domain inside a biosynthetic module.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum NrpsDomain {
    CStarter,
    CLcl,
    CDcl,
    Cy,
    A,
    T,
    E,
    Nmt,
    Ox,
    Red,
    Te,
    Tailoring(String),
}

impl NrpsDomain {
    pub fn short_label(&self) -> String {
        match self {
            Self::CStarter => "C_start".to_string(),
            Self::CLcl => "C_LCL".to_string(),
            Self::CDcl => "C_DCL".to_string(),
            Self::Cy => "Cy".to_string(),
            Self::A => "A".to_string(),
            Self::T => "T".to_string(),
            Self::E => "E".to_string(),
            Self::Nmt => "nMT".to_string(),
            Self::Ox => "Ox".to_string(),
            Self::Red => "Red".to_string(),
            Self::Te => "TE".to_string(),
            Self::Tailoring(t) => format!("Tail({t})"),
        }
    }
}

/// Thioesterase / reductive product-release topology of an NRPS assembly line.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum TeReleaseMode {
    #[default]
    Hydrolysis,
    HeadToTailMacrocyclization,
    Macrolactonization {
        target_module: usize,
    },
    LariatLactam {
        target_module: usize,
    },
    Homodimerization,
    ReductiveAlcohol,
    ReductiveAldehyde,
}

impl TeReleaseMode {
    pub fn label(&self) -> String {
        match self {
            Self::Hydrolysis => "TE Hydrolysis (Linear Acid)".to_string(),
            Self::HeadToTailMacrocyclization => "TE Head-to-Tail Macrocyclization".to_string(),
            Self::Macrolactonization { target_module } => {
                format!("TE Macrolactonization (-> Mod {target_module} OH)")
            }
            Self::LariatLactam { target_module } => {
                format!("TE Lariat Lactam (-> Mod {target_module} NH2)")
            }
            Self::Homodimerization => "TE Iterative Homodimerization (2× Pentapeptide)".to_string(),
            Self::ReductiveAlcohol => "Red/TD 4e- Reduction (-CH2OH)".to_string(),
            Self::ReductiveAldehyde => "Red/TD 2e- Reduction (-CHO)".to_string(),
        }
    }
}

/// Single biosynthetic module inside an NRPS gene cluster.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NrpsModule {
    pub module_index: usize,
    pub subunit: String,
    pub domains: Vec<NrpsDomain>,
    /// Optional 10-residue Stachelhaus specificity-conferring code for the A-domain.
    pub stachelhaus_code: Option<String>,
    /// Optional explicitly assigned monomer (overridden or verified against Stachelhaus).
    pub assigned_monomer: Option<String>,
    /// Optional fatty-acyl N-starter code if module carries `CStarter`.
    pub starter_acyl: Option<String>,
    /// Non-canonical δ/γ-amide linkage to next module (e.g. ACVS Module 1 Aad).
    #[serde(default)]
    pub delta_linkage: bool,
}

/// Complete Non-Ribosomal Peptide Synthetase Biosynthetic Gene Cluster (BGC) specification.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NrpsBgc {
    pub name: String,
    pub organism: String,
    pub source_note: String,
    pub modules: Vec<NrpsModule>,
    pub release_mode: TeReleaseMode,
    pub tailoring_enzymes: Vec<String>,
}

/// Ranked Stachelhaus A-domain substrate prediction (probabilistic — never asserted as certainty).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StachelhausCandidate {
    pub monomer_code: String,
    pub reference_code_10aa: String,
    pub identity_pct: f64,
    pub similarity_pct: f64,
    pub probability: f64,
    pub reference_system: String,
}

/// Full prediction output for a single 10-aa Stachelhaus code.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StachelhausPrediction {
    pub query_code: String,
    pub confidence: String,
    pub candidates: Vec<StachelhausCandidate>,
}

struct StachelhausEntry {
    monomer: &'static str,
    code_10aa: &'static str,
    reference: &'static str,
}

const STACHELHAUS_TABLE: &[StachelhausEntry] = &[
    StachelhausEntry {
        monomer: "Phe",
        code_10aa: "DAWTIAAVCK",
        reference: "GrsA-M1 (Gramicidin S)",
    },
    StachelhausEntry {
        monomer: "Pro",
        code_10aa: "DVQFIAHVVK",
        reference: "GrsB-M1 (Gramicidin S)",
    },
    StachelhausEntry {
        monomer: "Val",
        code_10aa: "DAFWIGGTFK",
        reference: "GrsB-M2 (Gramicidin S)",
    },
    StachelhausEntry {
        monomer: "Orn",
        code_10aa: "DVGEIGSIDK",
        reference: "GrsB-M3 (Gramicidin S)",
    },
    StachelhausEntry {
        monomer: "Leu",
        code_10aa: "DAWFLGNVVK",
        reference: "GrsB-M4 (Gramicidin S)",
    },
    StachelhausEntry {
        monomer: "Ile",
        code_10aa: "DGFFLGVVYK",
        reference: "BacA-M1 (Bacitracin)",
    },
    StachelhausEntry {
        monomer: "Dab",
        code_10aa: "DVWHISLIDK",
        reference: "PmxE-M1 (Polymyxin B)",
    },
    StachelhausEntry {
        monomer: "Thr",
        code_10aa: "DFWNIGMVHK",
        reference: "DptA-M4 (Daptomycin)",
    },
    StachelhausEntry {
        monomer: "Ser",
        code_10aa: "DVWHLSLVDK",
        reference: "EntF / DptBC-M11",
    },
    StachelhausEntry {
        monomer: "Cys",
        code_10aa: "DHESDVGITY",
        reference: "PcbAB-M2 (ACVS)",
    },
    StachelhausEntry {
        monomer: "Cys",
        code_10aa: "DLYNMSLIWK",
        reference: "HMWP2 (Yersiniabactin Cy)",
    },
    StachelhausEntry {
        monomer: "Asp",
        code_10aa: "DLTKVGHVGK",
        reference: "DptA-M3 (Daptomycin)",
    },
    StachelhausEntry {
        monomer: "Asn",
        code_10aa: "DLTKLGEVGK",
        reference: "DptA-M2 (Daptomycin)",
    },
    StachelhausEntry {
        monomer: "Glu",
        code_10aa: "DAKDLGVVDK",
        reference: "SrfA-M1 (Surfactin)",
    },
    StachelhausEntry {
        monomer: "Gln",
        code_10aa: "DAQDWGVVDK",
        reference: "TycC (Tyrocidine)",
    },
    StachelhausEntry {
        monomer: "Gly",
        code_10aa: "DILQLGLIWK",
        reference: "DptA-M5 (Daptomycin)",
    },
    StachelhausEntry {
        monomer: "Ala",
        code_10aa: "DLFNNALTYK",
        reference: "SimA-M7 (Cyclosporin A)",
    },
    StachelhausEntry {
        monomer: "Trp",
        code_10aa: "DVSAIGSVAK",
        reference: "DptA-M1 (Daptomycin)",
    },
    StachelhausEntry {
        monomer: "Tyr",
        code_10aa: "DASTVAAVCK",
        reference: "TycC-M3 (Tyrocidine)",
    },
    StachelhausEntry {
        monomer: "Lys",
        code_10aa: "DAEDIGTVGK",
        reference: "BacB-M2 (Bacitracin)",
    },
    StachelhausEntry {
        monomer: "Arg",
        code_10aa: "DVEDVGAIDK",
        reference: "McyC (Microcystin)",
    },
    StachelhausEntry {
        monomer: "His",
        code_10aa: "DSALIAEVWK",
        reference: "BacC (Bacitracin)",
    },
    StachelhausEntry {
        monomer: "Met",
        code_10aa: "DGFWIGVTYK",
        reference: "NRPSpredictor2 Met",
    },
    StachelhausEntry {
        monomer: "Hpg",
        code_10aa: "DIFHLGLLCK",
        reference: "CepA/B (Vancomycin)",
    },
    StachelhausEntry {
        monomer: "Dpg",
        code_10aa: "DPYQLGMLCK",
        reference: "CepC-M7 (Vancomycin)",
    },
    StachelhausEntry {
        monomer: "Bht",
        code_10aa: "DASTLGAVCK",
        reference: "CepA-M2 (Vancomycin)",
    },
    StachelhausEntry {
        monomer: "Aad",
        code_10aa: "EPRNIVEVGK",
        reference: "PcbAB-M1 (ACVS Penicillin)",
    },
    StachelhausEntry {
        monomer: "Pip",
        code_10aa: "DVQLLAHVVK",
        reference: "RapP (Rapamycin)",
    },
    StachelhausEntry {
        monomer: "Abu",
        code_10aa: "DAWFFGLVFK",
        reference: "SimA-M2 (Cyclosporin A)",
    },
    StachelhausEntry {
        monomer: "Kyn",
        code_10aa: "DAWTVAAVCK",
        reference: "DptD-M13 (Daptomycin)",
    },
    StachelhausEntry {
        monomer: "mGlu",
        code_10aa: "DAKDLGVVYK",
        reference: "DptD-M12 (Daptomycin)",
    },
    StachelhausEntry {
        monomer: "Bmt",
        code_10aa: "DAWYLGNVVK",
        reference: "SimA-M1 (Cyclosporin A)",
    },
    StachelhausEntry {
        monomer: "Sar",
        code_10aa: "DILQIGLIWK",
        reference: "SimA-M3 (Cyclosporin A)",
    },
    StachelhausEntry {
        monomer: "Hiv",
        code_10aa: "GALWIGGTFK",
        reference: "Vlm1-M2 (Valinomycin)",
    },
    StachelhausEntry {
        monomer: "Lac",
        code_10aa: "GALFNALTYK",
        reference: "Vlm2-M2 (Valinomycin)",
    },
    StachelhausEntry {
        monomer: "Aib",
        code_10aa: "DLFYLALTYK",
        reference: "Tex1 (Alamethicin Peptaibol)",
    },
    StachelhausEntry {
        monomer: "Dap",
        code_10aa: "DVWHLSLIDK",
        reference: "VioA (Viomycin)",
    },
    StachelhausEntry {
        monomer: "Phg",
        code_10aa: "DILHLGLLCK",
        reference: "Pristina (Pristinamycin)",
    },
    StachelhausEntry {
        monomer: "bAla",
        code_10aa: "VDWVSLADK-",
        reference: "Exochelin / Bleomycin",
    },
    StachelhausEntry {
        monomer: "Dhb",
        code_10aa: "DFWNIGMVFK",
        reference: "NRPSpredictor2 Dhb",
    },
];

fn aa_chem_group(c: char) -> u8 {
    match c.to_ascii_uppercase() {
        'A' | 'V' | 'I' | 'L' | 'M' => 1,
        'F' | 'W' | 'Y' => 2,
        'S' | 'T' | 'N' | 'Q' => 3,
        'D' | 'E' => 4,
        'K' | 'R' | 'H' => 5,
        'G' | 'P' => 6,
        'C' => 7,
        _ => 0,
    }
}

/// Predict ranked substrate candidates with calibrated probabilities from a 10-aa Stachelhaus code.
#[must_use]
pub fn predict_a_domain_specificity(raw_code: &str) -> StachelhausPrediction {
    let clean: String = raw_code
        .trim()
        .to_ascii_uppercase()
        .chars()
        .filter(|c| c.is_ascii_alphabetic() || *c == '-')
        .take(10)
        .collect();
    let q_chars: Vec<char> = clean.chars().collect();
    let mut scored: Vec<(f64, f64, &StachelhausEntry)> =
        Vec::with_capacity(STACHELHAUS_TABLE.len());

    for entry in STACHELHAUS_TABLE {
        let r_chars: Vec<char> = entry.code_10aa.chars().collect();
        let len = q_chars.len().min(r_chars.len()).max(1);
        let mut exact = 0.0_f64;
        let mut sim = 0.0_f64;
        for i in 0..len {
            if q_chars[i] == r_chars[i] {
                exact += 1.0;
                sim += 1.0;
            } else if aa_chem_group(q_chars[i]) != 0
                && aa_chem_group(q_chars[i]) == aa_chem_group(r_chars[i])
            {
                sim += 0.55;
            }
        }
        let id_pct = (exact / 10.0) * 100.0;
        let sim_pct = (sim / 10.0) * 100.0;
        scored.push((id_pct, sim_pct, entry));
    }

    scored.sort_by(|a, b| b.1.total_cmp(&a.1).then_with(|| b.0.total_cmp(&a.0)));

    // Keep top 4 distinct monomer candidates and convert similarity scores via temperature softmax
    let mut top: Vec<(f64, f64, &StachelhausEntry)> = Vec::new();
    for item in scored {
        if !top.iter().any(|x| x.2.monomer == item.2.monomer) {
            top.push(item);
            if top.len() == 4 {
                break;
            }
        }
    }

    let max_s = top.first().map_or(0.0, |x| x.1);
    let weights: Vec<f64> = top.iter().map(|x| ((x.1 - max_s) / 9.5).exp()).collect();
    let w_sum: f64 = weights.iter().sum::<f64>().max(1e-9);

    let candidates: Vec<StachelhausCandidate> = top
        .iter()
        .zip(weights.iter())
        .map(|(&(id_pct, sim_pct, entry), &w)| StachelhausCandidate {
            monomer_code: entry.monomer.to_string(),
            reference_code_10aa: entry.code_10aa.to_string(),
            identity_pct: (id_pct * 10.0).round() / 10.0,
            similarity_pct: (sim_pct * 10.0).round() / 10.0,
            probability: ((w / w_sum) * 1000.0).round() / 1000.0,
            reference_system: entry.reference.to_string(),
        })
        .collect();

    let best_id = candidates.first().map_or(0.0, |c| c.identity_pct);
    let confidence = if best_id >= 80.0 {
        "High (>=80% Stachelhaus identity)".to_string()
    } else if best_id >= 60.0 {
        "Moderate (60–79%; verify alternative substrates experimentally)".to_string()
    } else {
        "Low / Ambiguous (<60%; wet-lab ATP-PPi / LC-MS verification required)".to_string()
    };

    StachelhausPrediction {
        query_code: clean,
        confidence,
        candidates,
    }
}

/// Look up canonical Stachelhaus 10-aa code for a monomer (used during reverse-BGC blueprinting).
#[must_use]
pub fn lookup_stachelhaus_for_monomer(monomer_code: &str) -> Option<&'static str> {
    STACHELHAUS_TABLE
        .iter()
        .find(|e| e.monomer.eq_ignore_ascii_case(monomer_code))
        .map(|e| e.code_10aa)
}

/// Module-level diagnostic from NRPS collinearity and interface audit.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ModuleAudit {
    pub module_index: usize,
    pub resolved_monomer: String,
    pub is_d: bool,
    pub is_n_methyl: bool,
    pub stachelhaus_prediction: Option<StachelhausPrediction>,
    pub warnings: Vec<String>,
}

/// Compiled NRPS biosynthetic analysis uniting BGC audit, uncertainty, and assembled product.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NrpsCompileReport {
    pub bgc: NrpsBgc,
    pub module_audits: Vec<ModuleAudit>,
    pub interface_warnings: Vec<String>,
    pub collinear_compatible: bool,
    pub predicted_notation: String,
    pub predicted_helm: String,
    pub canonical_smiles: String,
    pub formula: String,
    pub monoisotopic_mass: f64,
    pub validation_controls: Vec<String>,
}

/// Audit biosynthetic collinearity, stereochemical interfaces, and compile an `NrpsBgc` into its
/// predicted 2D+3D NRP product (`BuiltPeptide`) + structured uncertainty report.
pub fn compile_nrps_bgc(
    bgc: &NrpsBgc,
    cat: &MonomerCatalog,
) -> Result<(NrpsCompileReport, BuiltPeptide), PeptideError> {
    if bgc.modules.is_empty() {
        return Err(PeptideError::Nrps(
            "NRPS BGC contains zero modules".to_string(),
        ));
    }

    let mut module_audits = Vec::with_capacity(bgc.modules.len());
    let mut interface_warnings = Vec::new();
    let mut residues = Vec::new();
    let mut n_cap = "H".to_string();
    let mut delta_links = Vec::new();

    for (idx, m) in bgc.modules.iter().enumerate() {
        let mut mod_warnings = Vec::new();
        let has_a = m.domains.contains(&NrpsDomain::A);
        let has_t = m.domains.contains(&NrpsDomain::T);
        let has_e = m.domains.contains(&NrpsDomain::E);
        let has_nmt = m.domains.contains(&NrpsDomain::Nmt);
        let has_cy = m.domains.contains(&NrpsDomain::Cy);
        let has_ox = m.domains.contains(&NrpsDomain::Ox);
        let has_c_starter = m.domains.contains(&NrpsDomain::CStarter);

        if !has_a || !has_t {
            mod_warnings.push(format!(
                "Module {} missing canonical A/T core (A={has_a}, T={has_t}): non-canonical or trans-acting assembly.",
                m.module_index
            ));
        }

        if idx == 0 && has_c_starter {
            if let Some(ref acyl) = m.starter_acyl {
                n_cap = acyl.clone();
            } else {
                mod_warnings.push(
                    "Module 1 has C_starter domain but no fatty-acyl starter specified (defaulting to C10)."
                        .to_string(),
                );
                n_cap = "C10".to_string();
            }
        }

        let stach_pred = m
            .stachelhaus_code
            .as_deref()
            .filter(|s| !s.trim().is_empty())
            .map(predict_a_domain_specificity);

        let mut base_code = if let Some(ref explicit) = m.assigned_monomer {
            if let Some(ref pred) = stach_pred
                && let Some(top) = pred.candidates.first()
                && !top.monomer_code.eq_ignore_ascii_case(explicit)
            {
                mod_warnings.push(format!(
                    "Discrepancy preserved: assigned '{}' vs top Stachelhaus candidate '{}' ({:.0}% id).",
                    explicit, top.monomer_code, top.identity_pct
                ));
            }
            explicit.clone()
        } else if let Some(ref pred) = stach_pred
            && let Some(top) = pred.candidates.first()
        {
            if top.identity_pct < 80.0 {
                let alts: Vec<String> = pred
                    .candidates
                    .iter()
                    .take(3)
                    .map(|c| format!("{} ({:.0}%)", c.monomer_code, c.probability * 100.0))
                    .collect();
                mod_warnings.push(format!(
                    "Uncertain A-domain specificity: candidates [{}].",
                    alts.join(", ")
                ));
            }
            top.monomer_code.clone()
        } else {
            return Err(PeptideError::Nrps(format!(
                "Module {} lacks both assigned_monomer and stachelhaus_code",
                m.module_index
            )));
        };

        // Heterocyclization Cy (+ Ox) domain check & conversion
        if has_cy {
            let ok_cy = base_code.eq_ignore_ascii_case("Cys")
                || base_code.eq_ignore_ascii_case("Ser")
                || base_code.eq_ignore_ascii_case("Thr")
                || base_code.eq_ignore_ascii_case("Thz")
                || base_code.eq_ignore_ascii_case("Oxz")
                || base_code.eq_ignore_ascii_case("MeOxz");
            if !ok_cy {
                interface_warnings.push(format!(
                    "Module {}: Cy heterocyclization domain is incompatible with substrate '{}' (requires Cys, Ser, or Thr).",
                    m.module_index, base_code
                ));
            } else if has_ox || base_code.eq_ignore_ascii_case("Cys") {
                base_code = match base_code.to_ascii_lowercase().as_str() {
                    "cys" => "Thz".to_string(),
                    "ser" => "Oxz".to_string(),
                    "thr" => "MeOxz".to_string(),
                    _ => base_code,
                };
            }
        }

        // Check E-domain -> downstream C_DCL donor stereochemistry interface
        if has_e
            && let Some(next_mod) = bgc.modules.get(idx + 1)
            && next_mod.domains.contains(&NrpsDomain::CLcl)
            && !next_mod.domains.contains(&NrpsDomain::CDcl)
        {
            interface_warnings.push(format!(
                "Interface incompatibility: Module {} carries Epimerization (E) producing a D-donor, but Module {} has C_LCL instead of C_DCL.",
                m.module_index, next_mod.module_index
            ));
        }

        let mspec = cat
            .find_monomer(&base_code)
            .ok_or_else(|| PeptideError::UnknownMonomer(base_code.clone()))?;

        let is_d = has_e && mspec.has_alpha_chiral;
        if has_e && !mspec.has_alpha_chiral {
            interface_warnings.push(format!(
                "Module {}: E-domain present on achiral monomer '{}' (no epimerization possible).",
                m.module_index, mspec.code
            ));
        }

        let is_n_methyl =
            has_nmt && !mspec.is_tertiary_n && !mspec.code.eq_ignore_ascii_case("Sar");
        if has_nmt && mspec.is_tertiary_n && !mspec.code.eq_ignore_ascii_case("Sar") {
            interface_warnings.push(format!(
                "Module {}: nMT domain incompatible with tertiary/hydroxy monomer '{}'.",
                m.module_index, mspec.code
            ));
        }

        if m.delta_linkage && idx + 1 < bgc.modules.len() {
            delta_links.push(CrossLink::DeltaAmide {
                acid_res: idx + 1,
                next_res: idx + 2,
            });
        }

        residues.push(ResidueUnit {
            base_code: mspec.code.clone(),
            is_d,
            is_n_methyl,
        });

        module_audits.push(ModuleAudit {
            module_index: m.module_index,
            resolved_monomer: mspec.code.clone(),
            is_d,
            is_n_methyl,
            stachelhaus_prediction: stach_pred,
            warnings: mod_warnings,
        });
    }

    // Apply TE release topology (including iterative homodimerization like Gramicidin S!)
    let mut cross_links = delta_links;
    let mut c_cap = CCapKind::FreeAcid;
    let last_mod = residues.len();
    match bgc.release_mode {
        TeReleaseMode::Hydrolysis => {}
        TeReleaseMode::HeadToTailMacrocyclization => {
            cross_links.push(CrossLink::HeadToTail);
        }
        TeReleaseMode::Homodimerization => {
            let copy = residues.clone();
            residues.extend(copy);
            cross_links.push(CrossLink::HeadToTail);
        }
        TeReleaseMode::Macrolactonization { target_module } => {
            if target_module == 0 || target_module > last_mod {
                interface_warnings.push(format!(
                    "TE macrolactonization target module {target_module} out of range 1..={last_mod}."
                ));
            } else {
                cross_links.push(CrossLink::Macrolactone {
                    from_res: last_mod,
                    to_oh_res: target_module,
                });
            }
        }
        TeReleaseMode::LariatLactam { target_module } => {
            if target_module == 0 || target_module > last_mod {
                interface_warnings.push(format!(
                    "TE lariat lactam target module {target_module} out of range 1..={last_mod}."
                ));
            } else {
                cross_links.push(CrossLink::Lactam {
                    from_acid_res: last_mod,
                    to_amine_res: target_module,
                    from_c_terminus: true,
                });
            }
        }
        TeReleaseMode::ReductiveAlcohol => {
            c_cap = CCapKind::Alcohol;
        }
        TeReleaseMode::ReductiveAldehyde => {
            c_cap = CCapKind::Aldehyde;
        }
    }

    let spec = PeptideSpec {
        name: Some(bgc.name.clone()),
        n_cap,
        c_cap,
        residues,
        cross_links,
        conformation: if matches!(
            bgc.release_mode,
            TeReleaseMode::HeadToTailMacrocyclization
                | TeReleaseMode::Homodimerization
                | TeReleaseMode::Macrolactonization { .. }
                | TeReleaseMode::LariatLactam { .. }
        ) {
            PeptideConformation::MacrocycleRing
        } else {
            PeptideConformation::BetaStrand
        },
    };

    let built = build_peptide_from_spec(&spec, cat)?;
    let monoisotopic_mass = mol_graph::analyze_molecule(&built.mol).exact_mass;
    let formula = built.mol.formula();
    let collinear_compatible = interface_warnings.is_empty();

    let validation_controls = vec![
        format!(
            "HR-LC-MS confirmation: verify exact [M+H]+ = {:.4} Da and [M+2H]2+ = {:.4} Da (Δm < 5 ppm) against empty-vector / Δnrps knockout control.",
            monoisotopic_mass + 1.007276,
            (monoisotopic_mass + 2.0 * 1.007276) / 2.0
        ),
        "Marfey's FDAA derivatization (L- vs D-FDLA LC-MS) to verify E-domain epimerization stereochemistry.".to_string(),
        "ATP-PPi exchange or hydroxylamine-trapping assay for any A-domain with <80% Stachelhaus identity.".to_string(),
    ];

    let report = NrpsCompileReport {
        bgc: bgc.clone(),
        module_audits,
        interface_warnings,
        collinear_compatible,
        predicted_notation: built.notation.clone(),
        predicted_helm: built.helm.clone(),
        canonical_smiles: built.canonical_smiles.clone(),
        formula,
        monoisotopic_mass,
        validation_controls,
    };

    Ok((report, built))
}

/// Reverse-engineer a minimal collinear NRPS BGC blueprint from any `PeptideSpec`.
#[must_use]
pub fn design_nrps_from_peptide(spec: &PeptideSpec) -> NrpsBgc {
    let mut modules = Vec::with_capacity(spec.residues.len());
    let mut tailoring = Vec::new();
    let has_starter =
        !spec.is_head_to_tail() && !spec.n_cap.is_empty() && !spec.n_cap.eq_ignore_ascii_case("H");

    let delta_set: Vec<usize> = spec
        .cross_links
        .iter()
        .filter_map(|cl| match *cl {
            CrossLink::DeltaAmide { acid_res, .. } => Some(acid_res),
            _ => None,
        })
        .collect();

    for (idx, r) in spec.residues.iter().enumerate() {
        let mut domains = Vec::new();
        let prev_was_d = idx > 0 && spec.residues[idx - 1].is_d;
        if idx == 0 {
            if has_starter {
                domains.push(NrpsDomain::CStarter);
            }
        } else if matches!(r.base_code.as_str(), "Thz" | "Oxz" | "MeOxz") {
            domains.push(NrpsDomain::Cy);
        } else if prev_was_d {
            domains.push(NrpsDomain::CDcl);
        } else {
            domains.push(NrpsDomain::CLcl);
        }

        domains.push(NrpsDomain::A);
        if r.is_n_methyl || r.base_code.eq_ignore_ascii_case("Sar") {
            domains.push(NrpsDomain::Nmt);
        }
        if matches!(r.base_code.as_str(), "Thz" | "Oxz" | "MeOxz") {
            domains.push(NrpsDomain::Ox);
        }
        domains.push(NrpsDomain::T);
        if r.is_d {
            domains.push(NrpsDomain::E);
        }

        if idx + 1 == spec.residues.len() {
            match spec.c_cap {
                CCapKind::Alcohol | CCapKind::Aldehyde => domains.push(NrpsDomain::Red),
                _ => domains.push(NrpsDomain::Te),
            }
        }

        if matches!(r.base_code.as_str(), "3Cl-Tyr" | "Bht") {
            tailoring.push(format!(
                "Post-NRPS Halogenase / P450 β-hydroxylase required for Mod {} ({})",
                idx + 1,
                r.base_code
            ));
        }
        if matches!(r.base_code.as_str(), "S5" | "R8" | "Pra" | "Aha") {
            tailoring.push(format!(
                "Synthetic/non-biological monomer at Mod {} ({}): requires Fmoc-SPPS or genetic code expansion",
                idx + 1,
                r.base_code
            ));
        }

        modules.push(NrpsModule {
            module_index: idx + 1,
            subunit: format!("Nrp{}", (b'A' + ((idx / 4) as u8).min(25)) as char),
            domains,
            stachelhaus_code: lookup_stachelhaus_for_monomer(&r.base_code).map(String::from),
            assigned_monomer: Some(r.base_code.clone()),
            starter_acyl: if idx == 0 && has_starter {
                Some(spec.n_cap.clone())
            } else {
                None
            },
            delta_linkage: delta_set.contains(&(idx + 1)),
        });
    }

    let mut release_mode = match spec.c_cap {
        CCapKind::Alcohol => TeReleaseMode::ReductiveAlcohol,
        CCapKind::Aldehyde => TeReleaseMode::ReductiveAldehyde,
        _ => TeReleaseMode::Hydrolysis,
    };
    for cl in &spec.cross_links {
        match *cl {
            CrossLink::HeadToTail => release_mode = TeReleaseMode::HeadToTailMacrocyclization,
            CrossLink::Macrolactone { to_oh_res, .. } => {
                release_mode = TeReleaseMode::Macrolactonization {
                    target_module: to_oh_res,
                }
            }
            CrossLink::Lactam {
                to_amine_res,
                from_c_terminus: true,
                ..
            } => {
                release_mode = TeReleaseMode::LariatLactam {
                    target_module: to_amine_res,
                }
            }
            CrossLink::Disulfide { res_a, res_b } => {
                tailoring.push(format!(
                    "Post-NRPS oxidative disulfide isomerase (DsbA/air) for Cys{res_a}–Cys{res_b}"
                ));
            }
            CrossLink::RcmStaple { res_a, res_b } => {
                tailoring.push(format!(
                    "Chemical Grubbs Ru-catalyst RCM stapling ({res_a}-{res_b}) — non-biosynthetic"
                ));
            }
            CrossLink::TriazoleClick {
                alkyne_res,
                azide_res,
            } => {
                tailoring.push(format!(
                    "Chemical CuAAC click cyclization ({alkyne_res}-{azide_res}) — non-biosynthetic"
                ));
            }
            _ => {}
        }
    }

    NrpsBgc {
        name: spec
            .name
            .clone()
            .unwrap_or_else(|| "Reverse-Designed NRPS Blueprint".to_string()),
        organism: "In silico engineered chassis (Streptomyces / E. coli BAP1 + sfp)".to_string(),
        source_note: "Reverse-designed from peptide/NRP specification by ChemCraft mol-peptide"
            .to_string(),
        modules,
        release_mode,
        tailoring_enzymes: tailoring,
    }
}

/// Evaluation of a proposed NRPS A-domain / module swap.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NrpsSwapEvaluation {
    pub module_index: usize,
    pub old_monomer: String,
    pub new_monomer: String,
    pub new_stachelhaus_code: Option<String>,
    pub interface_risk: String,
    pub interface_notes: Vec<String>,
    pub delta_mw: f64,
    pub delta_wlogp: f64,
    pub delta_tpsa: f64,
    pub delta_charge_ph74: f64,
    pub delta_bro5_score: f64,
    pub wild_type_smiles: String,
    pub engineered_smiles: String,
    pub engineered_notation: String,
}

/// Evaluate swapping module `module_index` (1-based) in `bgc` to `new_monomer`, auditing C-domain
/// acceptor/donor interface constraints and physicochemical deltas.
pub fn propose_nrps_swap(
    bgc: &NrpsBgc,
    module_index: usize,
    new_monomer: &str,
    cat: &MonomerCatalog,
) -> Result<(NrpsSwapEvaluation, BuiltPeptide), PeptideError> {
    if module_index == 0 || module_index > bgc.modules.len() {
        return Err(PeptideError::Nrps(format!(
            "module_index {module_index} out of range 1..={}",
            bgc.modules.len()
        )));
    }
    let target_spec = cat
        .find_monomer(new_monomer)
        .ok_or_else(|| PeptideError::UnknownMonomer(new_monomer.to_string()))?;

    let (_, wt_built) = compile_nrps_bgc(bgc, cat)?;
    let wt_dev = analyze_peptide_developability(&wt_built, cat);

    let mut eng_bgc = bgc.clone();
    let old_monomer = eng_bgc.modules[module_index - 1]
        .assigned_monomer
        .clone()
        .unwrap_or_else(|| "Unknown".to_string());
    let m = &mut eng_bgc.modules[module_index - 1];
    m.assigned_monomer = Some(target_spec.code.clone());
    m.stachelhaus_code = lookup_stachelhaus_for_monomer(&target_spec.code).map(String::from);

    let mut notes = Vec::new();
    let mut risk_pts = 0u32;

    if m.domains.contains(&NrpsDomain::CLcl) || m.domains.contains(&NrpsDomain::CDcl) {
        notes.push(format!(
            "Upstream C-domain acceptor pocket in Mod {module_index} evolved for '{old_monomer}'; swapping to '{}' may require XUT/XU C-A-T unit exchange rather than isolated A-domain swap.",
            target_spec.code
        ));
        if cat
            .find_monomer(&old_monomer)
            .is_some_and(|old_m| old_m.class != target_spec.class)
        {
            risk_pts += 2;
            notes.push("Cross-class monomer swap (polarity/charge change): high C-domain proofreading risk.".to_string());
        } else {
            risk_pts += 1;
        }
    }
    if let TeReleaseMode::Macrolactonization { target_module } = bgc.release_mode
        && target_module == module_index
        && target_spec.anchor_kind != crate::monomers::SidechainAnchorKind::Hydroxyl
    {
        return Err(PeptideError::Nrps(format!(
            "cannot swap TE macrolactonization anchor Mod {module_index} to '{}' (lacks -OH nucleophile)",
            target_spec.code
        )));
    }
    if let TeReleaseMode::LariatLactam { target_module } = bgc.release_mode
        && target_module == module_index
        && target_spec.anchor_kind != crate::monomers::SidechainAnchorKind::Amine
    {
        return Err(PeptideError::Nrps(format!(
            "cannot swap TE lariat-lactam anchor Mod {module_index} to '{}' (lacks -NH2 nucleophile)",
            target_spec.code
        )));
    }

    let (eng_rep, eng_built) = compile_nrps_bgc(&eng_bgc, cat)?;
    if !eng_rep.interface_warnings.is_empty() {
        risk_pts += 2;
        notes.extend(eng_rep.interface_warnings);
    }
    let eng_dev = analyze_peptide_developability(&eng_built, cat);

    let interface_risk = match risk_pts {
        0..=1 => "Low (Conservative A-domain / XU swap)".to_string(),
        2 => "Moderate (Recommend XU/XUT C-A-T swap boundary)".to_string(),
        _ => "High (Cross-class or domain-interface clash; redesign C/E/nMT boundary)".to_string(),
    };

    let eval = NrpsSwapEvaluation {
        module_index,
        old_monomer,
        new_monomer: target_spec.code.clone(),
        new_stachelhaus_code: lookup_stachelhaus_for_monomer(&target_spec.code).map(String::from),
        interface_risk,
        interface_notes: notes,
        delta_mw: ((eng_dev.monoisotopic_mass - wt_dev.monoisotopic_mass) * 1000.0).round()
            / 1000.0,
        delta_wlogp: ((eng_dev.wlogp - wt_dev.wlogp) * 100.0).round() / 100.0,
        delta_tpsa: ((eng_dev.tpsa - wt_dev.tpsa) * 100.0).round() / 100.0,
        delta_charge_ph74: ((eng_dev.charge_ph74 - wt_dev.charge_ph74) * 100.0).round() / 100.0,
        delta_bro5_score: ((eng_dev.bro5_oral_score - wt_dev.bro5_oral_score) * 1000.0).round()
            / 1000.0,
        wild_type_smiles: wt_built.canonical_smiles,
        engineered_smiles: eng_built.canonical_smiles.clone(),
        engineered_notation: eng_built.notation.clone(),
    };
    Ok((eval, eng_built))
}

/// Lightweight pure-Rust antiSMASH GenBank (`.gbk`) NRPS annotation parser.
/// Extracts `aSDomain` features (`Condensation_LCL`, `Condensation_DCL`, `Condensation_Starter`,
/// `Heterocyclization`, `AMP-binding`, `PCP`, `Epimerization`, `nMT`, `Thioesterase`, `TD`) and
/// `/NRPS_PKS` / `/specificity` Stachelhaus qualifiers into a runnable `NrpsBgc`.
pub fn parse_antismash_gbk(gbk_text: &str) -> Result<NrpsBgc, PeptideError> {
    let mut name = "antiSMASH_NRPS_Cluster".to_string();
    let mut organism = "Unknown organism".to_string();
    let mut modules: Vec<NrpsModule> = Vec::new();
    let mut cur_domains: Vec<NrpsDomain> = Vec::new();
    let mut cur_stach: Option<String> = None;
    let mut cur_sub: Option<String> = None;
    let mut release_mode = TeReleaseMode::Hydrolysis;
    let mut subunit_idx = 0usize;

    let flush_module = |modules: &mut Vec<NrpsModule>,
                        cur_domains: &mut Vec<NrpsDomain>,
                        cur_stach: &mut Option<String>,
                        cur_sub: &mut Option<String>,
                        subunit_idx: usize| {
        if cur_domains.contains(&NrpsDomain::A) {
            let idx = modules.len() + 1;
            let starter_acyl = if idx == 1 && cur_domains.contains(&NrpsDomain::CStarter) {
                Some("C10".to_string())
            } else {
                None
            };
            modules.push(NrpsModule {
                module_index: idx,
                subunit: format!("CDS_{}", subunit_idx.max(1)),
                domains: std::mem::take(cur_domains),
                stachelhaus_code: cur_stach.take(),
                assigned_monomer: cur_sub.take(),
                starter_acyl,
                delta_linkage: false,
            });
        }
    };

    for raw_line in gbk_text.lines() {
        let line = raw_line.trim();
        if let Some(rest) = line.strip_prefix("LOCUS") {
            if let Some(tok) = rest.split_whitespace().next() {
                name = tok.to_string();
            }
        } else if let Some(rest) = line.strip_prefix("ORGANISM") {
            let org = rest.trim();
            if !org.is_empty() {
                organism = org.to_string();
            }
        } else if line.starts_with("CDS ") || line.starts_with("CDS\t") {
            subunit_idx += 1;
        } else if line.contains("/aSDomain=") || line.contains("/NRPS_PKS=\"Domain:") {
            let lo = line.to_ascii_lowercase();
            let dom = if lo.contains("condensation_starter") {
                Some(NrpsDomain::CStarter)
            } else if lo.contains("condensation_dcl") {
                Some(NrpsDomain::CDcl)
            } else if lo.contains("condensation_lcl") || lo.contains("condensation_dual") {
                Some(NrpsDomain::CLcl)
            } else if lo.contains("heterocyclization") {
                Some(NrpsDomain::Cy)
            } else if lo.contains("amp-binding") {
                Some(NrpsDomain::A)
            } else if lo.contains("pcp") || lo.contains("pp-binding") {
                Some(NrpsDomain::T)
            } else if lo.contains("epimerization") {
                Some(NrpsDomain::E)
            } else if lo.contains("nmt") || lo.contains("methyltransf") {
                Some(NrpsDomain::Nmt)
            } else if lo.contains("thioesterase") {
                release_mode = TeReleaseMode::HeadToTailMacrocyclization;
                Some(NrpsDomain::Te)
            } else if lo.contains("td") && lo.contains("domain:") {
                release_mode = TeReleaseMode::ReductiveAlcohol;
                Some(NrpsDomain::Red)
            } else {
                None
            };
            if let Some(d) = dom {
                // Start a new module whenever a new C/Cy domain or second A domain arrives after an A domain
                let starts_new = matches!(
                    d,
                    NrpsDomain::CStarter
                        | NrpsDomain::CLcl
                        | NrpsDomain::CDcl
                        | NrpsDomain::Cy
                        | NrpsDomain::A
                ) && cur_domains.contains(&NrpsDomain::A);
                if starts_new {
                    flush_module(
                        &mut modules,
                        &mut cur_domains,
                        &mut cur_stach,
                        &mut cur_sub,
                        subunit_idx,
                    );
                }
                cur_domains.push(d);
            }
        } else if let Some(pos) = line.find("Stachelhaus code:") {
            let after = line[pos + "Stachelhaus code:".len()..].trim();
            let code: String = after
                .chars()
                .take_while(|c| c.is_ascii_alphabetic() || *c == '-')
                .collect();
            if code.len() >= 8 {
                cur_stach = Some(code.to_ascii_uppercase());
            }
        } else if let Some(pos) = line.find("consensus:") {
            let after = line[pos + "consensus:".len()..].trim();
            let tok: String = after
                .chars()
                .take_while(|c| c.is_ascii_alphanumeric() || *c == '-')
                .collect();
            if !tok.is_empty()
                && !tok.eq_ignore_ascii_case("X")
                && !tok.eq_ignore_ascii_case("no_call")
            {
                cur_sub = Some(tok);
            }
        }
    }
    flush_module(
        &mut modules,
        &mut cur_domains,
        &mut cur_stach,
        &mut cur_sub,
        subunit_idx,
    );

    if modules.is_empty() {
        return Err(PeptideError::Nrps(
            "no NRPS AMP-binding modules found in antiSMASH GenBank input".to_string(),
        ));
    }

    Ok(NrpsBgc {
        name,
        organism,
        source_note: "Imported from antiSMASH GenBank (.gbk) annotations".to_string(),
        modules,
        release_mode,
        tailoring_enzymes: Vec::new(),
    })
}

/// Curated landmark NRPS Biosynthetic Gene Cluster presets (`gramicidin_s`, `acv`, `daptomycin`, `polymyxin_b1`).
#[must_use]
pub fn preset_nrps_bgc(id: &str) -> Option<NrpsBgc> {
    let m = |i: usize,
             sub: &str,
             doms: Vec<NrpsDomain>,
             stach: &str,
             mon: &str,
             acyl: Option<&str>,
             delta: bool| NrpsModule {
        module_index: i,
        subunit: sub.to_string(),
        domains: doms,
        stachelhaus_code: Some(stach.to_string()),
        assigned_monomer: Some(mon.to_string()),
        starter_acyl: acyl.map(String::from),
        delta_linkage: delta,
    };
    use NrpsDomain as D;
    match id.trim().to_ascii_lowercase().as_str() {
        "gramicidin_s" | "grs" => Some(NrpsBgc {
            name: "Gramicidin S BGC (GrsA + GrsB)".to_string(),
            organism: "Aneurinibacillus migulanus ATCC 9999".to_string(),
            source_note: "Classic 5-module homodimerizing NRPS (Kessler/Marahiel paradigm)".to_string(),
            modules: vec![
                m(1, "GrsA", vec![D::A, D::T, D::E], "DAWTIAAVCK", "Phe", None, false),
                m(2, "GrsB", vec![D::CDcl, D::A, D::T], "DVQFIAHVVK", "Pro", None, false),
                m(3, "GrsB", vec![D::CLcl, D::A, D::T], "DAFWIGGTFK", "Val", None, false),
                m(4, "GrsB", vec![D::CLcl, D::A, D::T], "DVGEIGSIDK", "Orn", None, false),
                m(5, "GrsB", vec![D::CLcl, D::A, D::T, D::Te], "DAWFLGNVVK", "Leu", None, false),
            ],
            release_mode: TeReleaseMode::Homodimerization,
            tailoring_enzymes: Vec::new(),
        }),
        "acv" | "acv_tripeptide" | "pcbab" => Some(NrpsBgc {
            name: "ACV Synthetase BGC (PcbAB)".to_string(),
            organism: "Penicillium chrysogenum".to_string(),
            source_note: "3-module β-lactam NRPS precursor with δ-amide Aad1–Cys2 and E-domain D-Val3".to_string(),
            modules: vec![
                m(1, "PcbAB", vec![D::A, D::T], "EPRNIVEVGK", "Aad", None, true),
                m(2, "PcbAB", vec![D::CLcl, D::A, D::T], "DHESDVGITY", "Cys", None, false),
                m(3, "PcbAB", vec![D::CLcl, D::A, D::T, D::E, D::Te], "DAFWIGGTFK", "Val", None, false),
            ],
            release_mode: TeReleaseMode::Hydrolysis,
            tailoring_enzymes: vec!["IPNS (Isopenicillin N synthase — Fe2+ non-heme oxidase)".to_string()],
        }),
        "polymyxin_b1" | "pmx" => Some(NrpsBgc {
            name: "Polymyxin B1 BGC (PmxA/B/E)".to_string(),
            organism: "Paenibacillus polymyxa".to_string(),
            source_note: "10-module cationic lipopeptide NRPS with 6MeC8 C_starter and 23-atom lariat TE".to_string(),
            modules: vec![
                m(1, "PmxE", vec![D::CStarter, D::A, D::T], "DVWHISLIDK", "Dab", Some("6MeC8"), false),
                m(2, "PmxE", vec![D::CLcl, D::A, D::T], "DFWNIGMVHK", "Thr", None, false),
                m(3, "PmxE", vec![D::CLcl, D::A, D::T], "DVWHISLIDK", "Dab", None, false),
                m(4, "PmxE", vec![D::CLcl, D::A, D::T], "DVWHISLIDK", "Dab", None, false),
                m(5, "PmxE", vec![D::CLcl, D::A, D::T], "DVWHISLIDK", "Dab", None, false),
                m(6, "PmxA", vec![D::CLcl, D::A, D::T, D::E], "DAWTIAAVCK", "Phe", None, false),
                m(7, "PmxA", vec![D::CDcl, D::A, D::T], "DAWFLGNVVK", "Leu", None, false),
                m(8, "PmxA", vec![D::CLcl, D::A, D::T], "DVWHISLIDK", "Dab", None, false),
                m(9, "PmxA", vec![D::CLcl, D::A, D::T], "DVWHISLIDK", "Dab", None, false),
                m(10, "PmxB", vec![D::CLcl, D::A, D::T, D::Te], "DFWNIGMVHK", "Thr", None, false),
            ],
            release_mode: TeReleaseMode::LariatLactam { target_module: 4 },
            tailoring_enzymes: vec!["EctB-like Dab transaminase".to_string()],
        }),
        "daptomycin" | "dpt" => Some(NrpsBgc {
            name: "Daptomycin BGC (DptA / DptBC / DptD)".to_string(),
            organism: "Streptomyces roseosporus NRRL 11379".to_string(),
            source_note: "13-module Ca2+-dependent lipodepsipeptide NRPS with C10 starter and Thr4-Kyn13 TE lactone".to_string(),
            modules: vec![
                m(1, "DptA", vec![D::CStarter, D::A, D::T], "DVSAIGSVAK", "Trp", Some("C10"), false),
                m(2, "DptA", vec![D::CLcl, D::A, D::T, D::E], "DLTKLGEVGK", "Asn", None, false),
                m(3, "DptA", vec![D::CDcl, D::A, D::T], "DLTKVGHVGK", "Asp", None, false),
                m(4, "DptA", vec![D::CLcl, D::A, D::T], "DFWNIGMVHK", "Thr", None, false),
                m(5, "DptA", vec![D::CLcl, D::A, D::T], "DILQLGLIWK", "Gly", None, false),
                m(6, "DptBC", vec![D::CLcl, D::A, D::T], "DVGEIGSIDK", "Orn", None, false),
                m(7, "DptBC", vec![D::CLcl, D::A, D::T], "DLTKVGHVGK", "Asp", None, false),
                m(8, "DptBC", vec![D::CLcl, D::A, D::T, D::E], "DLFNNALTYK", "Ala", None, false),
                m(9, "DptBC", vec![D::CDcl, D::A, D::T], "DLTKVGHVGK", "Asp", None, false),
                m(10, "DptBC", vec![D::CLcl, D::A, D::T], "DILQLGLIWK", "Gly", None, false),
                m(11, "DptBC", vec![D::CLcl, D::A, D::T, D::E], "DVWHLSLVDK", "Ser", None, false),
                m(12, "DptD", vec![D::CDcl, D::A, D::T], "DAKDLGVVYK", "mGlu", None, false),
                m(13, "DptD", vec![D::CLcl, D::A, D::T, D::Te], "DAWTVAAVCK", "Kyn", None, false),
            ],
            release_mode: TeReleaseMode::Macrolactonization { target_module: 4 },
            tailoring_enzymes: vec!["DptI (Glu 3-methyltransferase)".to_string(), "DptE/F (Acyl-ACP ligase)".to_string()],
        }),
        _ => None,
    }
}

/// Render a structured Markdown dossier for an NRPS BGC compilation report.
#[must_use]
pub fn render_nrps_markdown(r: &NrpsCompileReport) -> String {
    let mut md = String::new();
    let _ = writeln!(md, "# NRPS Biosynthesis & Cluster Report: {}\n", r.bgc.name);
    let _ = writeln!(md, "- **Organism:** {}", r.bgc.organism);
    let _ = writeln!(md, "- **Source / Provenance:** {}", r.bgc.source_note);
    let _ = writeln!(md, "- **Release Mode:** {}", r.bgc.release_mode.label());
    let _ = writeln!(
        md,
        "- **Predicted Product:** `{}` | **Formula:** `{}` | **Exact Mass:** `{:.4} Da`",
        r.predicted_notation, r.formula, r.monoisotopic_mass
    );
    let _ = writeln!(md, "- **HELM:** `{}`", r.predicted_helm);
    let _ = writeln!(
        md,
        "- **Collinearity / Interface Status:** {}\n",
        if r.collinear_compatible {
            "PASS (Collinear & stereochemically compatible)"
        } else {
            "WARNINGS (Review domain interfaces below)"
        }
    );

    let _ = writeln!(md, "## Module Architecture & A-Domain Specificity\n");
    let _ = writeln!(
        md,
        "| Mod | Subunit | Domains | Stachelhaus (10-aa) | Resolved Monomer | Top Candidates (Prob.) |"
    );
    let _ = writeln!(md, "|---:|---|---|---|---|---|");
    for (m, aud) in r.bgc.modules.iter().zip(r.module_audits.iter()) {
        let doms: Vec<String> = m.domains.iter().map(NrpsDomain::short_label).collect();
        let stach = m.stachelhaus_code.as_deref().unwrap_or("-");
        let mut res_lbl = String::new();
        if aud.is_n_methyl {
            res_lbl.push_str("NMe-");
        }
        if aud.is_d {
            res_lbl.push_str("D-");
        }
        res_lbl.push_str(&aud.resolved_monomer);
        let cands = aud
            .stachelhaus_prediction
            .as_ref()
            .map(|p| {
                p.candidates
                    .iter()
                    .take(3)
                    .map(|c| {
                        format!(
                            "{} ({:.0}% id, p={:.2})",
                            c.monomer_code, c.identity_pct, c.probability
                        )
                    })
                    .collect::<Vec<_>>()
                    .join("; ")
            })
            .unwrap_or_else(|| "-".to_string());
        let _ = writeln!(
            md,
            "| {} | {} | `{}` | `{}` | **{}** | {} |",
            m.module_index,
            m.subunit,
            doms.join("-"),
            stach,
            res_lbl,
            cands
        );
    }

    if !r.interface_warnings.is_empty() {
        let _ = writeln!(md, "\n## Interface & Collinearity Warnings\n");
        for w in &r.interface_warnings {
            let _ = writeln!(md, "- {w}");
        }
    }
    if !r.bgc.tailoring_enzymes.is_empty() {
        let _ = writeln!(md, "\n## Tailoring Enzymes\n");
        for t in &r.bgc.tailoring_enzymes {
            let _ = writeln!(md, "- {t}");
        }
    }
    let _ = writeln!(md, "\n## Recommended Experimental Validation Controls\n");
    for v in &r.validation_controls {
        let _ = writeln!(md, "- {v}");
    }
    md
}
