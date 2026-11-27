//! Pure-Rust Combinatorial Chemical Library Enumeration & PROTAC Bifunctional Assembler.
//!
//! Supports two complementary MedChem library generation workflows:
//! 1. **R-Group / Markush Scaffold Enumeration (`EnumerateMode::RGroup`)**: replaces `[*]` / `[R1]`
//!    (and optional `[R2]`) attachment points on a scaffold with curated MedChem substituents or
//!    user-supplied building blocks.
//! 2. **PROTAC / TPD Bifunctional Assembly (`EnumerateMode::Protac`)**: joins a POI warhead SMILES
//!    to curated E3 ligase recruiters (CRBN pomalidomide/lenalidomide, VHL VH032, cIAP) across a
//!    panel of flexible (alkyl/PEG) and rigid (piperazine/piperidine/alkyne) linkers, automatically
//!    shifting ring-closure indices so assembled SMILES never collide and validating every product
//!    through `mol_graph::parse_smiles`.

use std::collections::BTreeSet;

use mol_graph::{MolGraph, parse_smiles};
use serde::{Deserialize, Serialize};

/// Enumeration workflow mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum EnumerateMode {
    RGroup,
    Protac,
}

/// Built-in E3 ligase recruiter selection for PROTAC assembly.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum E3RecruiterPreset {
    Crbn,
    Vhl,
    Iap,
    All,
}

/// One enumerated compound validated through `mol_graph::parse_smiles`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EnumeratedCompound {
    pub id: String,
    pub canonical_smiles: String,
    pub structural_key: String,
    pub formula: String,
    pub mw: f64,
    pub heavy_atoms: usize,
    pub linker_or_r1: String,
    pub recruiter_or_r2: String,
    pub bro5_ok: bool,
    pub ro5_mw_ok: bool,
    #[serde(skip)]
    pub mol: Option<MolGraph>,
}

/// Complete enumeration report with deduplication and validation telemetry.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EnumeratedLibrary {
    pub mode: EnumerateMode,
    pub parent_input: String,
    pub attempted: usize,
    pub valid_unique: usize,
    pub invalid_or_duplicate: usize,
    pub compounds: Vec<EnumeratedCompound>,
}

#[derive(Debug, Clone, Copy)]
struct NamedFragment {
    name: &'static str,
    smiles: &'static str,
}

const DEFAULT_R_GROUPS: &[NamedFragment] = &[
    NamedFragment {
        name: "Me",
        smiles: "C",
    },
    NamedFragment {
        name: "Et",
        smiles: "CC",
    },
    NamedFragment {
        name: "iPr",
        smiles: "C(C)C",
    },
    NamedFragment {
        name: "cPr",
        smiles: "C1CC1",
    },
    NamedFragment {
        name: "F",
        smiles: "F",
    },
    NamedFragment {
        name: "Cl",
        smiles: "Cl",
    },
    NamedFragment {
        name: "Br",
        smiles: "Br",
    },
    NamedFragment {
        name: "CF3",
        smiles: "C(F)(F)F",
    },
    NamedFragment {
        name: "OMe",
        smiles: "OC",
    },
    NamedFragment {
        name: "OEt",
        smiles: "OCC",
    },
    NamedFragment {
        name: "OCF3",
        smiles: "OC(F)(F)F",
    },
    NamedFragment {
        name: "OH",
        smiles: "O",
    },
    NamedFragment {
        name: "NH2",
        smiles: "N",
    },
    NamedFragment {
        name: "NHMe",
        smiles: "NC",
    },
    NamedFragment {
        name: "NMe2",
        smiles: "N(C)C",
    },
    NamedFragment {
        name: "CN",
        smiles: "C#N",
    },
    NamedFragment {
        name: "CONH2",
        smiles: "C(=O)N",
    },
    NamedFragment {
        name: "CO2H",
        smiles: "C(=O)O",
    },
    NamedFragment {
        name: "SO2Me",
        smiles: "S(=O)(=O)C",
    },
    NamedFragment {
        name: "Morpholine",
        smiles: "N1CCOCC1",
    },
    NamedFragment {
        name: "Piperazine",
        smiles: "N1CCNCC1",
    },
    NamedFragment {
        name: "Pyrrolidine",
        smiles: "N1CCCC1",
    },
    NamedFragment {
        name: "Phenyl",
        smiles: "c1ccccc1",
    },
    NamedFragment {
        name: "Pyridyl4",
        smiles: "c1ccncc1",
    },
];

/// Linker segments written left-to-right from warhead exit-vector to E3 amide/amine (`C(=O)`).
const DEFAULT_PROTAC_LINKERS: &[NamedFragment] = &[
    NamedFragment {
        name: "Alk_C4",
        smiles: "CCCC(=O)",
    },
    NamedFragment {
        name: "Alk_C6",
        smiles: "CCCCCC(=O)",
    },
    NamedFragment {
        name: "Alk_C8",
        smiles: "CCCCCCCC(=O)",
    },
    NamedFragment {
        name: "Alk_C10",
        smiles: "CCCCCCCCCC(=O)",
    },
    NamedFragment {
        name: "PEG1",
        smiles: "CCOCC(=O)",
    },
    NamedFragment {
        name: "PEG2",
        smiles: "CCOCCOCC(=O)",
    },
    NamedFragment {
        name: "PEG3",
        smiles: "CCOCCOCCOCC(=O)",
    },
    NamedFragment {
        name: "PEG4",
        smiles: "CCOCCOCCOCCOCC(=O)",
    },
    NamedFragment {
        name: "PEG5",
        smiles: "CCOCCOCCOCCOCCOCC(=O)",
    },
    NamedFragment {
        name: "Pip_C2",
        smiles: "CCN1CCN(CC(=O))CC1",
    },
    NamedFragment {
        name: "Pip_C4",
        smiles: "CCCCN1CCN(CC(=O))CC1",
    },
    NamedFragment {
        name: "Amide_PEG2",
        smiles: "NC(=O)CCOCCOCC(=O)",
    },
    NamedFragment {
        name: "Ether_C5",
        smiles: "OCCCCC(=O)",
    },
    NamedFragment {
        name: "Alkyne_PEG2",
        smiles: "C#CCCOCCOCC(=O)",
    },
];

/// E3 ligase ligands starting at their linker-attachment NH handle.
const E3_CRBN_POMA: NamedFragment = NamedFragment {
    name: "CRBN_Pomalidomide",
    smiles: "Nc1cccc2c1C(=O)N(C3CCC(=O)NC3=O)C2=O",
};
const E3_CRBN_LENA: NamedFragment = NamedFragment {
    name: "CRBN_Lenalidomide",
    smiles: "Nc1cccc2c1CN(C3CCC(=O)NC3=O)C2=O",
};
const E3_VHL_VH032: NamedFragment = NamedFragment {
    name: "VHL_VH032",
    smiles: "N[C@@H](C(C)(C)C)C(=O)N1C[C@H](O)C[C@H]1C(=O)NCc2ccc(-c3scnc3C)cc2",
};
const E3_IAP_CORE: NamedFragment = NamedFragment {
    name: "IAP_LCL161Core",
    smiles: "N[C@@H](C1CCCCC1)C(=O)N2CCC[C@H]2C(=O)Nc3ccccc3",
};

/// Scan highest single- or two-digit (`%NN`) ring index used in a SMILES string.
fn max_ring_index(smiles: &str) -> u8 {
    let bytes = smiles.as_bytes();
    let mut i = 0usize;
    let mut in_bracket = false;
    let mut max_idx = 0u8;
    while i < bytes.len() {
        let b = bytes[i];
        if b == b'[' {
            in_bracket = true;
        } else if b == b']' {
            in_bracket = false;
        } else if !in_bracket {
            if b == b'%' && i + 2 < bytes.len() {
                if let (Some(d1), Some(d2)) = (
                    (bytes[i + 1] as char).to_digit(10),
                    (bytes[i + 2] as char).to_digit(10),
                ) {
                    let v = (d1 * 10 + d2).min(99) as u8;
                    max_idx = max_idx.max(v);
                }
                i += 2;
            } else if b.is_ascii_digit() {
                let v = b - b'0';
                max_idx = max_idx.max(v);
            }
        }
        i += 1;
    }
    max_idx
}

/// Shift all non-bracket ring closure indices in `smiles` by `offset` to avoid collisions
/// when splicing two SMILES fragments together.
#[must_use]
pub fn shift_ring_indices(smiles: &str, offset: u8) -> String {
    if offset == 0 {
        return smiles.to_string();
    }
    let bytes = smiles.as_bytes();
    let mut out = String::with_capacity(smiles.len() + 8);
    let mut i = 0usize;
    let mut in_bracket = false;
    while i < bytes.len() {
        let b = bytes[i];
        if b == b'[' {
            in_bracket = true;
            out.push('[');
        } else if b == b']' {
            in_bracket = false;
            out.push(']');
        } else if !in_bracket && b == b'%' && i + 2 < bytes.len() {
            if let (Some(d1), Some(d2)) = (
                (bytes[i + 1] as char).to_digit(10),
                (bytes[i + 2] as char).to_digit(10),
            ) {
                let shifted = ((d1 * 10 + d2) as u8).saturating_add(offset);
                out.push('%');
                out.push_str(&format!("{shifted:02}"));
                i += 3;
                continue;
            }
            out.push('%');
        } else if !in_bracket && b.is_ascii_digit() && b != b'0' {
            let shifted = (b - b'0').saturating_add(offset);
            if shifted >= 10 {
                out.push('%');
                out.push_str(&format!("{shifted:02}"));
            } else {
                out.push(char::from(b'0' + shifted));
            }
        } else {
            out.push(b as char);
        }
        i += 1;
    }
    out
}

/// Splice a substituent at a placeholder (`[*]`, `[R1]`, or `[R2]`) inside `scaffold`.
fn splice_placeholder(scaffold: &str, placeholder: &str, fragment: &str) -> String {
    let offset = max_ring_index(scaffold);
    let shifted = shift_ring_indices(fragment, offset);
    scaffold.replacen(placeholder, &shifted, 1)
}

/// Enumerate an R-group / Markush combinatorial library from a scaffold containing `[*]` / `[R1]`
/// (and optional `[R2]`).
#[must_use]
pub fn enumerate_rgroup_library(
    scaffold: &str,
    custom_r1: Option<&[(String, String)]>,
    custom_r2: Option<&[(String, String)]>,
    max_compounds: usize,
) -> EnumeratedLibrary {
    let r1_list: Vec<(String, String)> = match custom_r1 {
        Some(list) if !list.is_empty() => list.to_vec(),
        _ => DEFAULT_R_GROUPS
            .iter()
            .map(|f| (f.name.to_string(), f.smiles.to_string()))
            .collect(),
    };
    let has_r1_tag = scaffold.contains("[R1]");
    let has_star = scaffold.contains("[*]");
    let has_r2 = scaffold.contains("[R2]");

    let r1_marker = if has_r1_tag {
        "[R1]"
    } else if has_star {
        "[*]"
    } else {
        ""
    };

    let r2_list: Vec<(String, String)> = if has_r2 {
        match custom_r2 {
            Some(list) if !list.is_empty() => list.to_vec(),
            _ => DEFAULT_R_GROUPS
                .iter()
                .take(10)
                .map(|f| (f.name.to_string(), f.smiles.to_string()))
                .collect(),
        }
    } else {
        vec![("none".to_string(), String::new())]
    };

    let mut seen = BTreeSet::new();
    let mut compounds = Vec::new();
    let mut attempted = 0usize;
    let mut invalid_or_duplicate = 0usize;

    let normalized_scaffold = if r1_marker.is_empty() {
        format!("({scaffold})[*]")
    } else {
        scaffold.to_string()
    };
    let active_r1 = if r1_marker.is_empty() {
        "[*]"
    } else {
        r1_marker
    };

    'outer: for (r1_name, r1_smi) in &r1_list {
        for (r2_name, r2_smi) in &r2_list {
            if compounds.len() >= max_compounds.max(1) {
                break 'outer;
            }
            attempted += 1;
            let step1 = splice_placeholder(&normalized_scaffold, active_r1, r1_smi);
            let candidate_smi = if has_r2 {
                splice_placeholder(&step1, "[R2]", r2_smi)
            } else {
                step1
            };
            match parse_smiles(&candidate_smi) {
                Ok(mut mol) => {
                    let key = mol.structural_key();
                    if !seen.insert(key.clone()) {
                        invalid_or_duplicate += 1;
                        continue;
                    }
                    let id = if has_r2 {
                        format!("ENUM_{r1_name}_{r2_name}")
                    } else {
                        format!("ENUM_{r1_name}")
                    };
                    mol.name = Some(id.clone());
                    let mw = mol.molecular_weight();
                    let heavy = mol.heavy_atom_count();
                    let canonical = mol.canonical_smiles();
                    let formula = mol.formula();
                    compounds.push(EnumeratedCompound {
                        id,
                        canonical_smiles: canonical,
                        structural_key: key,
                        formula,
                        mw: (mw * 100.0).round() / 100.0,
                        heavy_atoms: heavy,
                        linker_or_r1: r1_name.clone(),
                        recruiter_or_r2: r2_name.clone(),
                        bro5_ok: mw <= 1000.0,
                        ro5_mw_ok: mw <= 500.0,
                        mol: Some(mol),
                    });
                }
                Err(_) => {
                    invalid_or_duplicate += 1;
                }
            }
        }
    }

    EnumeratedLibrary {
        mode: EnumerateMode::RGroup,
        parent_input: scaffold.to_string(),
        attempted,
        valid_unique: compounds.len(),
        invalid_or_duplicate,
        compounds,
    }
}

/// Assemble a bifunctional PROTAC / TPD library by linking `warhead` (with optional `[*]` exit
/// vector) across curated alkyl/PEG/rigid linkers to CRBN, VHL, and/or IAP E3 recruiters.
#[must_use]
pub fn enumerate_protac_library(
    warhead: &str,
    e3_preset: E3RecruiterPreset,
    custom_linkers: Option<&[(String, String)]>,
    max_compounds: usize,
) -> EnumeratedLibrary {
    let recruiters: Vec<NamedFragment> = match e3_preset {
        E3RecruiterPreset::Crbn => vec![E3_CRBN_POMA, E3_CRBN_LENA],
        E3RecruiterPreset::Vhl => vec![E3_VHL_VH032],
        E3RecruiterPreset::Iap => vec![E3_IAP_CORE],
        E3RecruiterPreset::All => vec![E3_CRBN_POMA, E3_CRBN_LENA, E3_VHL_VH032, E3_IAP_CORE],
    };

    let linkers: Vec<(String, String)> = match custom_linkers {
        Some(list) if !list.is_empty() => list.to_vec(),
        _ => DEFAULT_PROTAC_LINKERS
            .iter()
            .map(|f| (f.name.to_string(), f.smiles.to_string()))
            .collect(),
    };

    let mut seen = BTreeSet::new();
    let mut compounds = Vec::new();
    let mut attempted = 0usize;
    let mut invalid_or_duplicate = 0usize;

    let warhead_tpl = if warhead.contains("[*]") {
        warhead.to_string()
    } else {
        format!("{warhead}O[*]")
    };

    'outer: for rec in &recruiters {
        for (lnk_name, lnk_smi) in &linkers {
            if compounds.len() >= max_compounds.max(1) {
                break 'outer;
            }
            attempted += 1;
            // Shift linker rings above warhead, then E3 rings above both warhead + linker
            let off1 = max_ring_index(&warhead_tpl);
            let shifted_lnk = shift_ring_indices(lnk_smi, off1);
            let off2 = off1.saturating_add(max_ring_index(&shifted_lnk));
            let shifted_e3 = shift_ring_indices(rec.smiles, off2);

            // Handle linkers that embed `(=O))` inside a branch vs linear terminal `C(=O)`
            let linker_plus_e3 = if let Some(pos) = shifted_lnk.rfind("C(=O))") {
                format!(
                    "{}C(=O)({shifted_e3}){}",
                    &shifted_lnk[..pos],
                    &shifted_lnk[pos + 5..]
                )
            } else {
                format!("{shifted_lnk}{shifted_e3}")
            };

            let assembled = warhead_tpl.replacen("[*]", &linker_plus_e3, 1);
            match parse_smiles(&assembled) {
                Ok(mut mol) => {
                    let key = mol.structural_key();
                    if !seen.insert(key.clone()) {
                        invalid_or_duplicate += 1;
                        continue;
                    }
                    let id = format!("PROTAC_{}_{}", rec.name, lnk_name);
                    mol.name = Some(id.clone());
                    let mw = mol.molecular_weight();
                    let heavy = mol.heavy_atom_count();
                    let canonical = mol.canonical_smiles();
                    let formula = mol.formula();
                    compounds.push(EnumeratedCompound {
                        id,
                        canonical_smiles: canonical,
                        structural_key: key,
                        formula,
                        mw: (mw * 100.0).round() / 100.0,
                        heavy_atoms: heavy,
                        linker_or_r1: lnk_name.clone(),
                        recruiter_or_r2: rec.name.to_string(),
                        bro5_ok: mw <= 1050.0,
                        ro5_mw_ok: mw <= 500.0,
                        mol: Some(mol),
                    });
                }
                Err(_) => {
                    invalid_or_duplicate += 1;
                }
            }
        }
    }

    EnumeratedLibrary {
        mode: EnumerateMode::Protac,
        parent_input: warhead.to_string(),
        attempted,
        valid_unique: compounds.len(),
        invalid_or_duplicate,
        compounds,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ring_index_shifter_prevents_collisions_across_fused_systems() {
        assert_eq!(shift_ring_indices("c1ccccc1", 3), "c4ccccc4");
        assert_eq!(shift_ring_indices("C1CC2CCCC12", 8), "C9CC%10CCCC9%10");
    }

    #[test]
    fn rgroup_enumeration_produces_valid_deduplicated_library() {
        let lib = enumerate_rgroup_library("c1ccc([*])cc1", None, None, 100);
        assert_eq!(lib.mode, EnumerateMode::RGroup);
        assert!(
            lib.valid_unique >= 20,
            "expected >=20 valid R-group analogs, got {}",
            lib.valid_unique
        );
        assert!(lib.compounds.iter().all(|c| c.mol.is_some()));
    }

    #[test]
    fn protac_assembly_produces_valid_bifunctional_crbn_and_vhl_degraders() {
        // Kinase-like warhead with explicit exit vector `[*]` and 2 existing rings
        let warhead = "Cc1ccc(-c2ccnc(N[*])n2)cc1";
        let lib = enumerate_protac_library(warhead, E3RecruiterPreset::All, None, 100);
        assert_eq!(lib.mode, EnumerateMode::Protac);
        assert!(
            lib.valid_unique >= 40,
            "expected >=40 valid assembled PROTACs across CRBN/VHL/IAP linkers, got {}",
            lib.valid_unique
        );
        assert!(
            lib.compounds
                .iter()
                .any(|c| c.recruiter_or_r2.starts_with("CRBN"))
        );
        assert!(
            lib.compounds
                .iter()
                .any(|c| c.recruiter_or_r2.starts_with("VHL"))
        );
        assert!(lib.compounds.iter().all(|c| c.mw > 450.0));
    }
}
