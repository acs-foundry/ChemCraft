#![forbid(unsafe_code)]

//! Extended peptide / NRP sequence parser, Pistoia Alliance HELM importer/exporter,
//! covalent topology cross-links, and 8 landmark natural-product & therapeutic NRP presets.

use serde::{Deserialize, Serialize};
use thiserror::Error;

use crate::monomers::{CCapKind, MonomerCatalog};

/// Maximum residues allowed in a single interactive/MCP peptide or NRP build.
pub const MAX_NRP_RESIDUES: usize = 60;

#[derive(Debug, Error)]
pub enum PeptideError {
    #[error("empty peptide / NRP specification")]
    EmptyInput,
    #[error("invalid peptide/NRP notation: {0}")]
    Notation(String),
    #[error("unknown monomer '{0}'")]
    UnknownMonomer(String),
    #[error("incompatible chemical modification or topology: {0}")]
    Topology(String),
    #[error("NRPS biosynthesis error: {0}")]
    Nrps(String),
    #[error("molecular graph error: {0}")]
    Graph(#[from] mol_graph::MolError),
}

/// Backbone 3D conformation preset for linear or lariat regions.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
pub enum PeptideConformation {
    #[default]
    BetaStrand,
    AlphaHelix,
    Extended,
    PolyprolineII,
    BetaHairpin,
    MacrocycleRing,
}

impl PeptideConformation {
    pub fn label(self) -> &'static str {
        match self {
            Self::BetaStrand => "Beta-Strand",
            Self::AlphaHelix => "Alpha-Helix",
            Self::Extended => "Extended",
            Self::PolyprolineII => "Polyproline-II",
            Self::BetaHairpin => "Beta-Hairpin",
            Self::MacrocycleRing => "Macrocycle Ring",
        }
    }

    /// Standard L-residue `(phi, psi, omega)` torsion angles in degrees.
    pub fn phi_psi_omega_deg(self) -> (f64, f64, f64) {
        match self {
            Self::BetaStrand => (-139.0, 135.0, 180.0),
            Self::AlphaHelix => (-57.0, -47.0, 180.0),
            Self::Extended => (-180.0, 180.0, 180.0),
            Self::PolyprolineII => (-75.0, 145.0, 180.0),
            Self::BetaHairpin => (-135.0, 135.0, 180.0),
            Self::MacrocycleRing => (-120.0, 120.0, 180.0),
        }
    }

    pub fn from_str_lossy(s: &str) -> Option<Self> {
        match s.trim().to_ascii_lowercase().as_str() {
            "beta" | "beta-strand" | "beta_strand" | "strand" | "sheet" => Some(Self::BetaStrand),
            "alpha" | "alpha-helix" | "alpha_helix" | "helix" => Some(Self::AlphaHelix),
            "extended" | "ext" | "linear" => Some(Self::Extended),
            "ppii" | "polyproline" | "polyproline-ii" | "polyproline_ii" => {
                Some(Self::PolyprolineII)
            }
            "hairpin" | "beta-hairpin" | "beta_hairpin" | "turn" => Some(Self::BetaHairpin),
            "ring" | "macrocycle" | "cyclic" | "macrocycle-ring" => Some(Self::MacrocycleRing),
            _ => None,
        }
    }
}

/// Covalent cyclization or intramolecular cross-link (1-indexed residue positions).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum CrossLink {
    /// Head-to-tail macrolactam: C-terminus of last residue -> N-terminus of first residue.
    HeadToTail,
    /// Depsipeptide macrolactone: carbonyl `from_res` (usually C-term) -> sidechain OH of `to_oh_res`.
    Macrolactone { from_res: usize, to_oh_res: usize },
    /// Lactam (lariat C-term -> sidechain NH2, or sidechain Asp/Glu COOH -> sidechain Lys/Orn/Dab NH2).
    Lactam {
        from_acid_res: usize,
        to_amine_res: usize,
        from_c_terminus: bool,
    },
    /// Disulfide bridge Cys(i)-SG -- Cys(j)-SG.
    Disulfide { res_a: usize, res_b: usize },
    /// All-hydrocarbon ring-closing metathesis (RCM) staple between S5/R8 alkenes.
    RcmStaple { res_a: usize, res_b: usize },
    /// Lanthionine / methyllanthionine thioether bridge: Cys(thiol) -> Dha/Dhb(Cβ).
    Lanthionine {
        thiol_res: usize,
        dehydro_res: usize,
    },
    /// CuAAC 1,4-triazole click staple between Pra (alkyne) and Aha (azide).
    TriazoleClick { alkyne_res: usize, azide_res: usize },
    /// Sidechain δ/γ-amide backbone linkage (e.g., ACV tripeptide L-δ-Aad1 -> L-Cys2).
    DeltaAmide { acid_res: usize, next_res: usize },
}

impl CrossLink {
    pub fn label(&self) -> String {
        match self {
            Self::HeadToTail => "cyclo(head-to-tail)".to_string(),
            Self::Macrolactone {
                from_res,
                to_oh_res,
            } => {
                format!("lactone({from_res}->{to_oh_res})")
            }
            Self::Lactam {
                from_acid_res,
                to_amine_res,
                from_c_terminus,
            } => {
                if *from_c_terminus {
                    format!("lariat({from_acid_res}->{to_amine_res})")
                } else {
                    format!("lactam({from_acid_res}->{to_amine_res})")
                }
            }
            Self::Disulfide { res_a, res_b } => format!("ss({res_a}-{res_b})"),
            Self::RcmStaple { res_a, res_b } => format!("staple({res_a}-{res_b})"),
            Self::Lanthionine {
                thiol_res,
                dehydro_res,
            } => format!("lan({thiol_res}->{dehydro_res})"),
            Self::TriazoleClick {
                alkyne_res,
                azide_res,
            } => format!("click({alkyne_res}-{azide_res})"),
            Self::DeltaAmide { acid_res, next_res } => format!("delta({acid_res}->{next_res})"),
        }
    }
}

/// One parsed monomer position inside a peptide or NRP.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ResidueUnit {
    /// Canonical monomer code in `MonomerCatalog` (e.g. `"Phe"`, `"Bmt"`, `"Dab"`).
    pub base_code: String,
    /// `true` when inverted to D-stereochemistry at Cα (NRPS E-domain / D-AA).
    pub is_d: bool,
    /// `true` when backbone N is N-methylated (NRPS nMT domain).
    pub is_n_methyl: bool,
}

impl ResidueUnit {
    pub fn display_token(&self) -> String {
        let mut out = String::new();
        if self.is_n_methyl && !self.base_code.eq_ignore_ascii_case("Sar") {
            out.push_str("NMe-");
        }
        if self.is_d {
            out.push_str("D-");
        }
        out.push_str(&self.base_code);
        out
    }

    pub fn helm_token(&self, cat: &MonomerCatalog) -> String {
        if self.base_code.eq_ignore_ascii_case("Sar") {
            return "[Sar]".to_string();
        }
        let spec = cat.find_monomer(&self.base_code);
        if !self.is_d
            && !self.is_n_methyl
            && let Some(m) = spec
            && m.helm_symbol.len() == 1
        {
            return m.helm_symbol.clone();
        }
        let base_h = spec.map_or(self.base_code.as_str(), |m| m.helm_symbol.as_str());
        let mut inner = String::new();
        if self.is_n_methyl {
            inner.push_str("me");
        }
        if self.is_d {
            inner.push('d');
        }
        inner.push_str(base_h);
        format!("[{inner}]")
    }
}

/// Fully parsed peptide or Non-Ribosomal Peptide design specification.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct PeptideSpec {
    pub name: Option<String>,
    pub n_cap: String,
    pub c_cap: CCapKind,
    pub residues: Vec<ResidueUnit>,
    pub cross_links: Vec<CrossLink>,
    pub conformation: PeptideConformation,
}

impl PeptideSpec {
    #[must_use]
    pub fn is_cyclic(&self) -> bool {
        !self.cross_links.is_empty()
    }

    #[must_use]
    pub fn is_head_to_tail(&self) -> bool {
        self.cross_links
            .iter()
            .any(|cl| matches!(cl, CrossLink::HeadToTail))
    }

    #[must_use]
    pub fn has_c_term_cyclization(&self) -> bool {
        let n = self.residues.len();
        self.cross_links.iter().any(|cl| match cl {
            CrossLink::HeadToTail => true,
            CrossLink::Macrolactone { from_res, .. } => *from_res == n,
            CrossLink::Lactam {
                from_acid_res,
                from_c_terminus,
                ..
            } => *from_c_terminus && *from_acid_res == n,
            _ => false,
        })
    }

    /// Format readable extended ChemCraft peptide/NRP notation.
    #[must_use]
    pub fn to_notation(&self) -> String {
        let body: Vec<String> = self
            .residues
            .iter()
            .map(ResidueUnit::display_token)
            .collect();
        let joined = body.join("-");
        let mut s = if self.is_head_to_tail() {
            format!("cyclo({joined})")
        } else {
            let prefix = if self.n_cap.is_empty() || self.n_cap.eq_ignore_ascii_case("H") {
                String::new()
            } else {
                format!("{}-", self.n_cap)
            };
            let suffix = if self.has_c_term_cyclization() || self.c_cap == CCapKind::FreeAcid {
                String::new()
            } else {
                format!("-{}", self.c_cap.code())
            };
            format!("{prefix}{joined}{suffix}")
        };
        let extra_links: Vec<String> = self
            .cross_links
            .iter()
            .filter(|cl| !matches!(cl, CrossLink::HeadToTail))
            .map(CrossLink::label)
            .collect();
        if !extra_links.is_empty() {
            s.push_str(" | ");
            s.push_str(&extra_links.join(", "));
        }
        s
    }

    /// Export Pistoia Alliance HELM 1.0/2.0 notation (`PEPTIDE1{...}$...$$$`).
    #[must_use]
    pub fn to_helm(&self, cat: &MonomerCatalog) -> String {
        let mut tokens: Vec<String> = Vec::new();
        let mut offset = 0usize;
        if !self.is_head_to_tail()
            && !self.n_cap.is_empty()
            && !self.n_cap.eq_ignore_ascii_case("H")
            && let Some(cap) = cat.find_n_cap(&self.n_cap)
            && !cap.helm_symbol.is_empty()
        {
            tokens.push(cap.helm_symbol.clone());
            offset = 1;
        }
        for r in &self.residues {
            tokens.push(r.helm_token(cat));
        }
        if !self.has_c_term_cyclization() && !self.c_cap.helm_symbol().is_empty() {
            tokens.push(self.c_cap.helm_symbol().to_string());
        }
        let n = self.residues.len();
        let mut conns = Vec::new();
        for cl in &self.cross_links {
            match *cl {
                CrossLink::HeadToTail => {
                    conns.push(format!(
                        "PEPTIDE1,PEPTIDE1,{}:R1-{}:R2",
                        1 + offset,
                        n + offset
                    ));
                }
                CrossLink::Macrolactone {
                    from_res,
                    to_oh_res,
                } => {
                    conns.push(format!(
                        "PEPTIDE1,PEPTIDE1,{}:R2-{}:R3",
                        from_res + offset,
                        to_oh_res + offset
                    ));
                }
                CrossLink::Lactam {
                    from_acid_res,
                    to_amine_res,
                    from_c_terminus,
                } => {
                    let att = if from_c_terminus { "R2" } else { "R3" };
                    conns.push(format!(
                        "PEPTIDE1,PEPTIDE1,{}:{att}-{}:R3",
                        from_acid_res + offset,
                        to_amine_res + offset
                    ));
                }
                CrossLink::Disulfide { res_a, res_b } => {
                    conns.push(format!(
                        "PEPTIDE1,PEPTIDE1,{}:R3-{}:R3",
                        res_a + offset,
                        res_b + offset
                    ));
                }
                CrossLink::RcmStaple { res_a, res_b } => {
                    conns.push(format!(
                        "PEPTIDE1,PEPTIDE1,{}:R3-{}:R3\"RCM\"",
                        res_a + offset,
                        res_b + offset
                    ));
                }
                CrossLink::Lanthionine {
                    thiol_res,
                    dehydro_res,
                } => {
                    conns.push(format!(
                        "PEPTIDE1,PEPTIDE1,{}:R3-{}:R3\"LAN\"",
                        thiol_res + offset,
                        dehydro_res + offset
                    ));
                }
                CrossLink::TriazoleClick {
                    alkyne_res,
                    azide_res,
                } => {
                    conns.push(format!(
                        "PEPTIDE1,PEPTIDE1,{}:R3-{}:R3\"CLICK\"",
                        alkyne_res + offset,
                        azide_res + offset
                    ));
                }
                CrossLink::DeltaAmide { acid_res, next_res } => {
                    conns.push(format!(
                        "PEPTIDE1,PEPTIDE1,{}:R3-{}:R1\"DELTA\"",
                        acid_res + offset,
                        next_res + offset
                    ));
                }
            }
        }
        format!("PEPTIDE1{{{}}}${}$$$", tokens.join("."), conns.join("|"))
    }
}

/// Curated landmark natural-product NRP & therapeutic peptide preset.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NrpPreset {
    pub id: &'static str,
    pub name: &'static str,
    pub organism_or_class: &'static str,
    pub notation: &'static str,
    pub expected_formula: &'static str,
    pub expected_monoisotopic_mass: f64,
    pub summary: &'static str,
}

/// All 8 curated landmark NRP and therapeutic peptide presets.
#[must_use]
pub fn landmark_nrp_presets() -> Vec<NrpPreset> {
    vec![
        NrpPreset {
            id: "cyclosporin_a",
            name: "Cyclosporin A",
            organism_or_class: "Tolypocladium inflatum (NRPS Macrocycle, Immunosuppressant)",
            notation: "cyclo(NMe-Bmt-Abu-Sar-NMe-Leu-Val-NMe-Leu-Ala-D-Ala-NMe-Leu-NMe-Leu-NMe-Val)",
            expected_formula: "C62H111N11O12",
            expected_monoisotopic_mass: 1201.8414,
            summary: "Oral bRo5 cyclic undecapeptide with 7 N-methylated amides, D-Ala8, Abu2, and signature MeBmt1.",
        },
        NrpPreset {
            id: "gramicidin_s",
            name: "Gramicidin S",
            organism_or_class: "Aneurinibacillus migulanus (NRPS Homodimer, Antimicrobial)",
            notation: "cyclo(Val-Orn-Leu-D-Phe-Pro-Val-Orn-Leu-D-Phe-Pro)",
            expected_formula: "C60H92N12O10",
            expected_monoisotopic_mass: 1140.7059,
            summary: "Classic amphipathic β-sheet cyclic decapeptide dimerized by GrsB TE with 2× L-Orn and 2× D-Phe.",
        },
        NrpPreset {
            id: "cilengitide",
            name: "Cilengitide",
            organism_or_class: "Synthetic Cyclic Pentapeptide (αvβ3/αvβ5 Integrin Inhibitor)",
            notation: "cyclo(Arg-Gly-Asp-D-Phe-NMe-Val)",
            expected_formula: "C27H40N8O7",
            expected_monoisotopic_mass: 588.3020,
            summary: "Head-to-tail RGDf[NMe-V] cyclic pentapeptide combining D-Phe turn induction and NMe-Val stabilization.",
        },
        NrpPreset {
            id: "octreotide",
            name: "Octreotide",
            organism_or_class: "Somatostatin Analog (Disulfide Cyclic Octapeptide)",
            notation: "D-Phe-Cys-Phe-D-Trp-Lys-Thr-Cys-Thr-ol | ss(2-7)",
            expected_formula: "C49H66N10O10S2",
            expected_monoisotopic_mass: 1018.4405,
            summary: "Protease-resistant somatostatin agonist with Cys2–Cys7 disulfide, D-Phe1, D-Trp4, and C-terminal Thr-ol.",
        },
        NrpPreset {
            id: "polymyxin_b1",
            name: "Polymyxin B1",
            organism_or_class: "Paenibacillus polymyxa (NRPS Cationic Lipopeptide)",
            notation: "6MeC8-Dab-Thr-Dab-Dab-Dab-D-Phe-Leu-Dab-Dab-Thr | lariat(10->4)",
            expected_formula: "C56H98N16O13",
            expected_monoisotopic_mass: 1202.7499,
            summary: "Last-resort Gram-negative lipopeptide with 6-methyloctanoyl starter, 6× L-Dab, D-Phe6, and 23-atom lariat ring.",
        },
        NrpPreset {
            id: "daptomycin",
            name: "Daptomycin",
            organism_or_class: "Streptomyces roseosporus (NRPS Lipodepsipeptide)",
            notation: "C10-Trp-D-Asn-Asp-Thr-Gly-Orn-Asp-D-Ala-Asp-Gly-D-Ser-mGlu-Kyn | lactone(13->4)",
            expected_formula: "C72H101N17O26",
            expected_monoisotopic_mass: 1619.7104,
            summary: "Anionic Ca2+-dependent lipodepsipeptide with C10 starter, D-Asn/D-Ala/D-Ser, Orn, 3-MeGlu, Kyn, and Thr4–Kyn13 macrolactone.",
        },
        NrpPreset {
            id: "acv_tripeptide",
            name: "ACV Tripeptide",
            organism_or_class: "Penicillium / Acremonium (ACVS 3-Module NRPS, β-Lactam Precursor)",
            notation: "Aad-Cys-D-Val | delta(1->2)",
            expected_formula: "C14H25N3O6S",
            expected_monoisotopic_mass: 363.1464,
            summary: "L-δ-(α-aminoadipoyl)-L-Cys-D-Val synthesized by ACVS with non-canonical δ-amide linkage and E-domain D-Val.",
        },
        NrpPreset {
            id: "semaglutide_core",
            name: "Semaglutide GLP-1(7-16) Stabilized Core",
            organism_or_class: "DPP-4-Resistant Therapeutic GLP-1 Mimetic Core",
            notation: "His-Aib-Glu-Gly-Thr-Phe-Thr-Ser-Asp-Val-NH2",
            expected_formula: "C46H69N13O17",
            expected_monoisotopic_mass: 1075.4934,
            summary: "Demonstrates α,α-disubstituted Aib2 blocking DPP-4 cleavage compared with native GLP-1 (Ala2).",
        },
    ]
}

/// Look up a curated preset by id or partial name (case-insensitive).
#[must_use]
pub fn find_nrp_preset(query: &str) -> Option<NrpPreset> {
    let q = query.trim().to_ascii_lowercase().replace([' ', '-'], "_");
    landmark_nrp_presets()
        .into_iter()
        .find(|p| p.id == q || p.name.to_ascii_lowercase().replace([' ', '-'], "_") == q)
}

fn parse_single_residue_token(
    raw: &str,
    cat: &MonomerCatalog,
) -> Result<ResidueUnit, PeptideError> {
    let mut tok = raw.trim().trim_matches(['[', ']']).to_string();
    if tok.is_empty() {
        return Err(PeptideError::Notation("empty residue token".to_string()));
    }

    // Special direct alias: Sar == NMe-Gly
    if tok.eq_ignore_ascii_case("Sar")
        || tok.eq_ignore_ascii_case("MeGly")
        || tok.eq_ignore_ascii_case("NMe-Gly")
        || tok.eq_ignore_ascii_case("NMeGly")
    {
        return Ok(ResidueUnit {
            base_code: "Sar".to_string(),
            is_d: false,
            is_n_methyl: false, // methyl is built into Sar template
        });
    }

    let mut is_n_methyl = false;
    let mut is_d = false;

    // Peel N-methyl prefixes (`NMe-`, `NMe`, `Me-`, `Me` before residue name, avoiding `Met`/`MeOxz`)
    for pref in ["NMe-", "nme-", "N-Me-", "me-", "NMe", "nme"] {
        if let Some(rest) = tok.strip_prefix(pref)
            && !rest.is_empty()
        {
            is_n_methyl = true;
            tok = rest.to_string();
            break;
        }
    }
    if !is_n_methyl
        && (tok.starts_with("Me") || tok.starts_with("me"))
        && !tok.eq_ignore_ascii_case("Met")
        && !tok.eq_ignore_ascii_case("MeOxz")
        && !tok.eq_ignore_ascii_case("MeGlu")
        && tok.len() >= 3
    {
        let rest = &tok[2..];
        let rest_clean = rest.strip_prefix('-').unwrap_or(rest);
        if cat.find_monomer(rest_clean).is_some()
            || rest_clean
                .strip_prefix("D-")
                .or_else(|| rest_clean.strip_prefix("d-"))
                .is_some_and(|r| cat.find_monomer(r).is_some())
        {
            is_n_methyl = true;
            tok = rest_clean.to_string();
        }
    }

    // Peel D-stereochemistry prefixes (`D-`, `d-`, `d` + uppercase letter like `dF` / `dPhe`)
    if let Some(rest) = tok.strip_prefix("D-").or_else(|| tok.strip_prefix("d-")) {
        is_d = true;
        tok = rest.to_string();
    } else if let Some(rest) = tok.strip_prefix("L-").or_else(|| tok.strip_prefix("l-")) {
        tok = rest.to_string();
    } else if tok.starts_with('d') && tok.len() >= 2 {
        let rest = &tok[1..];
        if rest.chars().next().is_some_and(|c| c.is_ascii_uppercase())
            && !tok.eq_ignore_ascii_case("Dab")
            && !tok.eq_ignore_ascii_case("Dap")
            && !tok.eq_ignore_ascii_case("Deg")
            && !tok.eq_ignore_ascii_case("Dha")
            && !tok.eq_ignore_ascii_case("Dhb")
            && !tok.eq_ignore_ascii_case("Dpg")
            && cat.find_monomer(rest).is_some()
        {
            is_d = true;
            tok = rest.to_string();
        }
    }

    // Allow `NMe-` after `D-` too (e.g. `D-NMe-Phe`)
    if !is_n_methyl {
        for pref in ["NMe-", "nme-", "Me-"] {
            if let Some(rest) = tok.strip_prefix(pref) {
                is_n_methyl = true;
                tok = rest.to_string();
                break;
            }
        }
    }

    let spec = cat
        .find_monomer(&tok)
        .ok_or_else(|| PeptideError::UnknownMonomer(raw.trim().to_string()))?;

    if is_d && !spec.has_alpha_chiral {
        return Err(PeptideError::Topology(format!(
            "monomer '{}' is achiral at Cα and cannot take D-stereochemistry",
            spec.code
        )));
    }
    if is_n_methyl && spec.is_tertiary_n {
        return Err(PeptideError::Topology(format!(
            "monomer '{}' already has a tertiary backbone N or depsipeptide O and cannot be N-methylated",
            spec.code
        )));
    }

    Ok(ResidueUnit {
        base_code: spec.code.clone(),
        is_d,
        is_n_methyl,
    })
}

fn parse_pair_indices(inner: &str, sep: &str) -> Result<(usize, usize), PeptideError> {
    let parts: Vec<&str> = inner.split(sep).map(str::trim).collect();
    if parts.len() != 2 {
        return Err(PeptideError::Notation(format!(
            "expected two 1-based residue indices separated by '{sep}', got '{inner}'"
        )));
    }
    let a: usize = parts[0]
        .parse()
        .map_err(|_| PeptideError::Notation(format!("invalid residue index '{}'", parts[0])))?;
    let b: usize = parts[1]
        .parse()
        .map_err(|_| PeptideError::Notation(format!("invalid residue index '{}'", parts[1])))?;
    if a == 0 || b == 0 || a == b {
        return Err(PeptideError::Notation(format!(
            "cross-link indices must be distinct and >= 1, got ({a}, {b})"
        )));
    }
    Ok((a, b))
}

fn parse_cross_link_directive(dir: &str) -> Result<CrossLink, PeptideError> {
    let s = dir.trim();
    let lo = s.to_ascii_lowercase();
    if lo == "cyclo" || lo == "head-to-tail" || lo == "head_to_tail" {
        return Ok(CrossLink::HeadToTail);
    }
    let open = s
        .find('(')
        .ok_or_else(|| PeptideError::Notation(format!("invalid cross-link '{s}'")))?;
    let close = s
        .rfind(')')
        .ok_or_else(|| PeptideError::Notation(format!("unclosed ')' in cross-link '{s}'")))?;
    let kind = s[..open].trim().to_ascii_lowercase();
    let inner = s[open + 1..close].trim();
    let sep = if inner.contains("->") { "->" } else { "-" };
    let (a, b) = parse_pair_indices(inner, sep)?;
    match kind.as_str() {
        "ss" | "disulfide" => Ok(CrossLink::Disulfide { res_a: a, res_b: b }),
        "lactone" | "macrolactone" | "ester" => Ok(CrossLink::Macrolactone {
            from_res: a,
            to_oh_res: b,
        }),
        "lariat" | "lariat_lactam" => Ok(CrossLink::Lactam {
            from_acid_res: a,
            to_amine_res: b,
            from_c_terminus: true,
        }),
        "lactam" | "macrolactam" => Ok(CrossLink::Lactam {
            from_acid_res: a,
            to_amine_res: b,
            from_c_terminus: false,
        }),
        "staple" | "rcm" => Ok(CrossLink::RcmStaple { res_a: a, res_b: b }),
        "lan" | "lanthionine" | "thioether" => Ok(CrossLink::Lanthionine {
            thiol_res: a,
            dehydro_res: b,
        }),
        "click" | "triazole" => Ok(CrossLink::TriazoleClick {
            alkyne_res: a,
            azide_res: b,
        }),
        "delta" | "delta_amide" => Ok(CrossLink::DeltaAmide {
            acid_res: a,
            next_res: b,
        }),
        _ => Err(PeptideError::Notation(format!(
            "unknown cross-link directive '{kind}' (supported: cyclo, ss, lactone, lariat, lactam, staple, lan, click, delta)"
        ))),
    }
}

/// Tokenize dash-separated sequences while keeping `D-X`, `L-X`, `NMe-X`, `3Cl-Tyr`, `C18diacid-gGlu` intact.
fn tokenize_extended_sequence(body: &str) -> Vec<String> {
    let raw_parts: Vec<&str> = body
        .split(['-', ' ', ';'])
        .map(str::trim)
        .filter(|t| !t.is_empty())
        .collect();
    let mut merged: Vec<String> = Vec::new();
    let mut i = 0usize;
    while i < raw_parts.len() {
        let cur = raw_parts[i];
        // Compound prefixes that bind to the immediately following token
        if i + 1 < raw_parts.len() {
            if cur.eq_ignore_ascii_case("NMe")
                && i + 2 < raw_parts.len()
                && (raw_parts[i + 1].eq_ignore_ascii_case("D")
                    || raw_parts[i + 1].eq_ignore_ascii_case("L"))
            {
                merged.push(format!("{}-{}-{}", cur, raw_parts[i + 1], raw_parts[i + 2]));
                i += 3;
                continue;
            }
            if cur.eq_ignore_ascii_case("D")
                || cur.eq_ignore_ascii_case("L")
                || cur.eq_ignore_ascii_case("NMe")
                || cur.eq_ignore_ascii_case("beta")
                || cur.eq_ignore_ascii_case("3Cl")
                || cur.eq_ignore_ascii_case("C18diacid")
                || cur.eq_ignore_ascii_case("6")
            {
                merged.push(format!("{cur}-{}", raw_parts[i + 1]));
                i += 2;
                continue;
            }
        }
        merged.push(cur.to_string());
        i += 1;
    }
    merged
}

/// Parse extended peptide/NRP notation, preset name, or Pistoia HELM string into `PeptideSpec`.
pub fn parse_peptide_spec(
    raw: &str,
    cat: &MonomerCatalog,
    default_conf: PeptideConformation,
) -> Result<PeptideSpec, PeptideError> {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return Err(PeptideError::EmptyInput);
    }

    // 1. Direct HELM string
    if trimmed.starts_with("PEPTIDE1{") {
        return parse_helm(trimmed, cat, default_conf);
    }

    // 2. Landmark preset ID / name
    if let Some(preset) = find_nrp_preset(trimmed) {
        let mut spec = parse_peptide_spec(preset.notation, cat, default_conf)?;
        spec.name = Some(preset.name.to_string());
        return Ok(spec);
    }

    // 3. Split optional cross-link directives after `|`
    let (main_part, links_part) = match trimmed.split_once('|') {
        Some((a, b)) => (a.trim(), Some(b.trim())),
        None => (trimmed, None),
    };

    let mut cross_links = Vec::new();
    if let Some(lp) = links_part {
        for dir in lp.split(',').map(str::trim).filter(|s| !s.is_empty()) {
            cross_links.push(parse_cross_link_directive(dir)?);
        }
    }

    // 4. Detect head-to-tail wrapper `cyclo(...)` or `c(...)`
    let mut seq_str = main_part.to_string();
    if let Some(rest) = seq_str
        .strip_prefix("cyclo(")
        .or_else(|| seq_str.strip_prefix("CYCLO("))
        .or_else(|| seq_str.strip_prefix("c("))
        .or_else(|| seq_str.strip_prefix("C("))
        && let Some(inner) = rest.strip_suffix(')')
    {
        if !cross_links
            .iter()
            .any(|c| matches!(c, CrossLink::HeadToTail))
        {
            cross_links.push(CrossLink::HeadToTail);
        }
        seq_str = inner.trim().to_string();
    }

    let mut n_cap = "H".to_string();
    let mut c_cap = CCapKind::FreeAcid;
    let mut residues = Vec::new();

    let has_delims = seq_str.contains('-')
        || seq_str.contains(' ')
        || seq_str.contains(';')
        || seq_str.contains('[');

    if !has_delims {
        // Check if single multi-letter monomer code (e.g., "Ala", "Bmt", "Aib")
        if seq_str.len() >= 2
            && cat.find_monomer(&seq_str).is_some()
            && seq_str.chars().skip(1).any(|c| c.is_ascii_lowercase())
        {
            residues.push(parse_single_residue_token(&seq_str, cat)?);
        } else {
            // Compact 1-letter sequence (uppercase = L, lowercase = D, e.g., "RGDfV")
            for ch in seq_str.chars() {
                let is_d = ch.is_ascii_lowercase();
                let cu = ch.to_ascii_uppercase();
                let spec = cat
                    .find_monomer(&cu.to_string())
                    .ok_or_else(|| PeptideError::UnknownMonomer(ch.to_string()))?;
                if is_d && !spec.has_alpha_chiral {
                    return Err(PeptideError::Topology(format!(
                        "achiral residue '{}' cannot be lowercase D-form",
                        spec.code
                    )));
                }
                residues.push(ResidueUnit {
                    base_code: spec.code.clone(),
                    is_d,
                    is_n_methyl: false,
                });
            }
        }
    } else {
        let mut tokens = tokenize_extended_sequence(&seq_str);
        if tokens.is_empty() {
            return Err(PeptideError::EmptyInput);
        }

        // Optional leading N-cap (if first token is an N-cap and NOT a monomer like "H" = His when followed by only 1-letter?)
        if tokens.len() >= 2 {
            let first = &tokens[0];
            if first.eq_ignore_ascii_case("H") && tokens.len() >= 2 && tokens[1].len() >= 2 {
                n_cap = "H".to_string();
                tokens.remove(0);
            } else if !first.eq_ignore_ascii_case("H")
                && cat.find_monomer(first).is_none()
                && let Some(cap) = cat.find_n_cap(first)
            {
                n_cap = cap.code.clone();
                tokens.remove(0);
            }
        }

        // Optional trailing C-cap (`-NH2`, `-OH`, `-ol`, `-CHO`, `-OMe`, `-NHEt`)
        if tokens.len() >= 2
            && let Some(last) = tokens.last()
            && let Some(cc) = CCapKind::parse_suffix(last)
        {
            c_cap = cc;
            tokens.pop();
        }

        // Single compact 1-letter block inside `Ac-AGYS-NH2`
        if tokens.len() == 1
            && tokens[0].len() > 3
            && cat.find_monomer(&tokens[0]).is_none()
            && !tokens[0].contains('[')
            && !tokens[0].contains('-')
        {
            for ch in tokens[0].chars() {
                let is_d = ch.is_ascii_lowercase();
                let cu = ch.to_ascii_uppercase();
                let spec = cat
                    .find_monomer(&cu.to_string())
                    .ok_or_else(|| PeptideError::UnknownMonomer(ch.to_string()))?;
                residues.push(ResidueUnit {
                    base_code: spec.code.clone(),
                    is_d: is_d && spec.has_alpha_chiral,
                    is_n_methyl: false,
                });
            }
        } else {
            for tok in &tokens {
                residues.push(parse_single_residue_token(tok, cat)?);
            }
        }
    }

    if residues.is_empty() {
        return Err(PeptideError::EmptyInput);
    }
    if residues.len() > MAX_NRP_RESIDUES {
        return Err(PeptideError::Notation(format!(
            "residue count {} exceeds limit of {MAX_NRP_RESIDUES}",
            residues.len()
        )));
    }

    validate_cross_links(&residues, &cross_links)?;

    let conformation = if !cross_links.is_empty() && default_conf == PeptideConformation::BetaStrand
    {
        PeptideConformation::MacrocycleRing
    } else {
        default_conf
    };

    Ok(PeptideSpec {
        name: None,
        n_cap,
        c_cap,
        residues,
        cross_links,
        conformation,
    })
}

/// Parse Pistoia Alliance HELM (`PEPTIDE1{A.G.[dF].[meV]}$PEPTIDE1,PEPTIDE1,1:R1-4:R2$$$`).
pub fn parse_helm(
    helm: &str,
    cat: &MonomerCatalog,
    default_conf: PeptideConformation,
) -> Result<PeptideSpec, PeptideError> {
    let sections: Vec<&str> = helm.trim().split('$').collect();
    if sections.is_empty() {
        return Err(PeptideError::EmptyInput);
    }
    let poly_sec = sections[0].trim();
    let start = poly_sec
        .find('{')
        .ok_or_else(|| PeptideError::Notation("HELM missing '{'".to_string()))?;
    let end = poly_sec
        .rfind('}')
        .ok_or_else(|| PeptideError::Notation("HELM missing '}'".to_string()))?;
    let raw_monomers: Vec<&str> = poly_sec[start + 1..end]
        .split('.')
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .collect();
    if raw_monomers.is_empty() {
        return Err(PeptideError::EmptyInput);
    }
    // Enforce the cap here, on the raw monomer list, before anything is built. The other
    // parsing path checked it and this one did not, so HELM was a way around it: 60 kB of
    // notation expanded into a polymer large enough to exhaust several GiB.
    if raw_monomers.len() > MAX_NRP_RESIDUES {
        return Err(PeptideError::Notation(format!(
            "HELM residue count {} exceeds limit of {MAX_NRP_RESIDUES}",
            raw_monomers.len()
        )));
    }

    let mut n_cap = "H".to_string();
    let mut c_cap = CCapKind::FreeAcid;
    let mut offset = 0usize;
    let mut slice_start = 0usize;
    let mut slice_end = raw_monomers.len();

    if raw_monomers.len() >= 2 {
        let first_clean = raw_monomers[0].trim_matches(['[', ']']);
        if cat.find_monomer(first_clean).is_none()
            && let Some(cap) = cat.find_n_cap(first_clean)
        {
            n_cap = cap.code.clone();
            slice_start = 1;
            offset = 1;
        }
    }
    if slice_end > slice_start + 1 {
        let last_clean = raw_monomers[slice_end - 1].trim_matches(['[', ']']);
        if let Some(cc) = CCapKind::parse_suffix(last_clean) {
            c_cap = cc;
            slice_end -= 1;
        }
    }

    let mut residues = Vec::new();
    for tok in &raw_monomers[slice_start..slice_end] {
        residues.push(parse_single_residue_token(tok, cat)?);
    }

    let mut cross_links = Vec::new();
    if let Some(conn_sec) = sections.get(1)
        && !conn_sec.trim().is_empty()
    {
        for conn in conn_sec.split('|').map(str::trim).filter(|s| !s.is_empty()) {
            // Format: PEPTIDE1,PEPTIDE1,1:R1-5:R2["TAG"]
            let parts: Vec<&str> = conn.split(',').collect();
            if parts.len() < 3 {
                continue;
            }
            let spec_part = parts[2].trim();
            let (edge_str, tag) = match spec_part.split_once('"') {
                Some((e, rest)) => (e.trim(), rest.trim_matches('"').to_ascii_uppercase()),
                None => (spec_part, String::new()),
            };
            let Some((left, right)) = edge_str.split_once('-') else {
                continue;
            };
            let parse_end = |s: &str| -> Result<(usize, String), PeptideError> {
                let (idx_s, att) = s.split_once(':').ok_or_else(|| {
                    PeptideError::Notation(format!("invalid HELM attachment '{s}'"))
                })?;
                let raw_idx: usize = idx_s
                    .trim()
                    .parse()
                    .map_err(|_| PeptideError::Notation(format!("invalid HELM index '{idx_s}'")))?;
                let res_idx = raw_idx.checked_sub(offset).ok_or_else(|| {
                    PeptideError::Notation(format!("HELM index {raw_idx} points to terminal cap"))
                })?;
                Ok((res_idx, att.trim().to_ascii_uppercase()))
            };
            let (i1, a1) = parse_end(left)?;
            let (i2, a2) = parse_end(right)?;
            if tag == "RCM" {
                cross_links.push(CrossLink::RcmStaple {
                    res_a: i1,
                    res_b: i2,
                });
            } else if tag == "LAN" {
                cross_links.push(CrossLink::Lanthionine {
                    thiol_res: i1,
                    dehydro_res: i2,
                });
            } else if tag == "CLICK" {
                cross_links.push(CrossLink::TriazoleClick {
                    alkyne_res: i1,
                    azide_res: i2,
                });
            } else if tag == "DELTA" {
                cross_links.push(CrossLink::DeltaAmide {
                    acid_res: i1,
                    next_res: i2,
                });
            } else if (a1 == "R1" && a2 == "R2") || (a1 == "R2" && a2 == "R1") {
                cross_links.push(CrossLink::HeadToTail);
            } else if a1 == "R2" && a2 == "R3" {
                let target_spec = residues
                    .get(i2.saturating_sub(1))
                    .and_then(|r| cat.find_monomer(&r.base_code));
                if target_spec.is_some_and(|m| {
                    m.anchor_kind == crate::monomers::SidechainAnchorKind::Hydroxyl
                }) {
                    cross_links.push(CrossLink::Macrolactone {
                        from_res: i1,
                        to_oh_res: i2,
                    });
                } else {
                    cross_links.push(CrossLink::Lactam {
                        from_acid_res: i1,
                        to_amine_res: i2,
                        from_c_terminus: true,
                    });
                }
            } else if a1 == "R3" && a2 == "R3" {
                let m1 = residues
                    .get(i1.saturating_sub(1))
                    .and_then(|r| cat.find_monomer(&r.base_code));
                let m2 = residues
                    .get(i2.saturating_sub(1))
                    .and_then(|r| cat.find_monomer(&r.base_code));
                if m1.is_some_and(|m| m.anchor_kind == crate::monomers::SidechainAnchorKind::Thiol)
                    && m2.is_some_and(|m| {
                        m.anchor_kind == crate::monomers::SidechainAnchorKind::Thiol
                    })
                {
                    cross_links.push(CrossLink::Disulfide {
                        res_a: i1,
                        res_b: i2,
                    });
                } else {
                    cross_links.push(CrossLink::Lactam {
                        from_acid_res: i1,
                        to_amine_res: i2,
                        from_c_terminus: false,
                    });
                }
            }
        }
    }

    validate_cross_links(&residues, &cross_links)?;
    let conformation = if !cross_links.is_empty() {
        PeptideConformation::MacrocycleRing
    } else {
        default_conf
    };
    Ok(PeptideSpec {
        name: None,
        n_cap,
        c_cap,
        residues,
        cross_links,
        conformation,
    })
}

fn validate_cross_links(residues: &[ResidueUnit], links: &[CrossLink]) -> Result<(), PeptideError> {
    let n = residues.len();
    for cl in links {
        match *cl {
            CrossLink::HeadToTail => {
                if n < 3 {
                    return Err(PeptideError::Topology(
                        "head-to-tail macrocyclization requires at least 3 residues".to_string(),
                    ));
                }
            }
            CrossLink::Macrolactone {
                from_res,
                to_oh_res,
            } => {
                if from_res > n || to_oh_res > n {
                    return Err(PeptideError::Topology(format!(
                        "macrolactone ({from_res}->{to_oh_res}) exceeds peptide length {n}"
                    )));
                }
            }
            CrossLink::Lactam {
                from_acid_res,
                to_amine_res,
                ..
            } => {
                if from_acid_res > n || to_amine_res > n {
                    return Err(PeptideError::Topology(format!(
                        "lactam ({from_acid_res}->{to_amine_res}) exceeds peptide length {n}"
                    )));
                }
            }
            CrossLink::Disulfide { res_a, res_b } => {
                if res_a > n || res_b > n {
                    return Err(PeptideError::Topology(format!(
                        "disulfide ({res_a}-{res_b}) exceeds peptide length {n}"
                    )));
                }
            }
            CrossLink::RcmStaple { res_a, res_b } => {
                if res_a > n || res_b > n {
                    return Err(PeptideError::Topology(format!(
                        "RCM staple ({res_a}-{res_b}) exceeds peptide length {n}"
                    )));
                }
            }
            CrossLink::Lanthionine {
                thiol_res,
                dehydro_res,
            } => {
                if thiol_res > n || dehydro_res > n {
                    return Err(PeptideError::Topology(format!(
                        "lanthionine ({thiol_res}->{dehydro_res}) exceeds peptide length {n}"
                    )));
                }
            }
            CrossLink::TriazoleClick {
                alkyne_res,
                azide_res,
            } => {
                if alkyne_res > n || azide_res > n {
                    return Err(PeptideError::Topology(format!(
                        "click bridge ({alkyne_res}-{azide_res}) exceeds peptide length {n}"
                    )));
                }
            }
            CrossLink::DeltaAmide { acid_res, next_res } => {
                if acid_res > n || next_res > n || next_res != acid_res + 1 {
                    return Err(PeptideError::Topology(format!(
                        "delta-amide ({acid_res}->{next_res}) must connect adjacent residues within 1..={n}"
                    )));
                }
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod residue_cap_tests {
    use super::*;

    fn catalog() -> MonomerCatalog {
        MonomerCatalog::builtin()
    }

    /// Audit finding F3-17: `MAX_NRP_RESIDUES` was enforced on one parsing path and not on
    /// `parse_helm`, so HELM notation walked straight past it. A 5000-residue polymer was
    /// accepted and built; 30 000 residues (60 kB of input) exhausted a 6 GiB budget.
    #[test]
    fn helm_notation_obeys_the_residue_cap() {
        let cat = catalog();
        let huge = format!(
            "PEPTIDE1{{{}}}$$$$",
            vec!["A"; MAX_NRP_RESIDUES + 1].join(".")
        );
        let res = parse_helm(&huge, &cat, PeptideConformation::default());
        assert!(
            res.is_err(),
            "a {}-residue HELM polymer must be refused (cap {MAX_NRP_RESIDUES})",
            MAX_NRP_RESIDUES + 1
        );

        let massive = format!("PEPTIDE1{{{}}}$$$$", vec!["A"; 5_000].join("."));
        assert!(
            parse_helm(&massive, &cat, PeptideConformation::default()).is_err(),
            "5000 residues must be refused, not built"
        );
    }

    /// ...while a polymer inside the cap still parses, including the N/C cap forms.
    #[test]
    fn helm_within_the_cap_still_parses() {
        let cat = catalog();
        let ok = format!("PEPTIDE1{{{}}}$$$$", vec!["A"; MAX_NRP_RESIDUES].join("."));
        let spec = parse_helm(&ok, &cat, PeptideConformation::default());
        assert!(
            spec.is_ok(),
            "exactly {MAX_NRP_RESIDUES} residues must still be accepted: {:?}",
            spec.err()
        );

        let short = parse_helm("PEPTIDE1{A.G.V}$$$$", &cat, PeptideConformation::default());
        assert!(short.is_ok(), "{:?}", short.err());
    }
}
