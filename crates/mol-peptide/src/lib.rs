#![forbid(unsafe_code)]

//! Layer 2 (`mol-peptide`): Clean-room Peptide, Macrocycle, and Non-Ribosomal Peptide (NRP / NRPS)
//! Design Studio for ChemCraft.
//!
//! Unifies:
//! - 57+ curated base monomers (20 L/19 D canonical AAs, N-Me, α,α-disubstituted, β-AAs,
//!   depsipeptide α-hydroxy acids, NRPS cationic/aromatic/signature monomers including `Bmt`,
//!   `Orn`, `Dab`, `Kyn`, `Hpg`, `Dpg`, `Aad`, `mGlu`, `Pip`, `Thz`/`Oxz`, and RCM/CuAAC stapling
//!   blocks) + 14 N-terminal caps/lipids + 6 C-terminal caps + custom JSON loading.
//! - Extended sequence notation + Pistoia Alliance HELM I/O + macrocycle/depsipeptide/lariat/
//!   disulfide/staple/click/δ-amide covalent graph assembly with clean 2D & 3D geometry.
//! - Modular NRPS biosynthesis engine (probabilistic 10-aa Stachelhaus A-domain specificity,
//!   collinearity & `E -> C_DCL` interface audit, antiSMASH `.gbk` importer, forward compilation,
//!   reverse BGC blueprinting, and domain-swap evaluation).
//! - bRo5 oral macrocycle developability, exact `pI`, GRAVY/hydrophobic moment, proteolysis map
//!   (incl. DPP-4 `Aib2` blocking), ESI-MS/`b,y` fragments, Fmoc-SPPS vs. NRPS planning, and
//!   combinatorial optimization scans (`Ala`, `D-AA`, `N-Me`, `RCM Staple`).

pub mod builder;
pub mod developability;
pub mod geometry;
pub mod monomers;
pub mod notation;
pub mod nrps;
pub mod spps_and_scan;

pub use builder::{BuiltPeptide, build_peptide_from_spec, build_peptide_or_nrp};
pub use developability::{
    MsFragmentIon, PeptideDevelopabilityReport, PeptideLiability, ProteaseSite,
    analyze_peptide_developability, render_developability_markdown,
};
pub use monomers::{
    CCapKind, MonomerCatalog, MonomerClass, MonomerSpec, NCapSpec, SidechainAnchorKind,
};
pub use notation::{
    CrossLink, MAX_NRP_RESIDUES, NrpPreset, PeptideConformation, PeptideError, PeptideSpec,
    ResidueUnit, find_nrp_preset, landmark_nrp_presets, parse_helm, parse_peptide_spec,
};
pub use nrps::{
    ModuleAudit, NrpsBgc, NrpsCompileReport, NrpsDomain, NrpsModule, NrpsSwapEvaluation,
    StachelhausCandidate, StachelhausPrediction, TeReleaseMode, compile_nrps_bgc,
    design_nrps_from_peptide, lookup_stachelhaus_for_monomer, parse_antismash_gbk,
    predict_a_domain_specificity, preset_nrps_bgc, propose_nrps_swap, render_nrps_markdown,
};
pub use spps_and_scan::{
    PeptideScanBundle, ScanMode, ScanVariant, SppsCouplingAlert, SppsPlan, plan_spps_vs_nrps,
    render_scan_csv, render_spps_markdown, run_peptide_scan,
};

/// Return the default monomer catalog augmented with any JSON files in `$CRAFT_DATA_DIR/monomers/`.
#[must_use]
pub fn default_monomer_catalog() -> MonomerCatalog {
    let mut cat = MonomerCatalog::builtin();
    // `CRAFT_DATA_DIR`, else the XDG data directory. Never a path hardcoded for one machine.
    let data_dir = std::env::var("CRAFT_DATA_DIR")
        .ok()
        .filter(|v| !v.trim().is_empty())
        .map(std::path::PathBuf::from)
        .or_else(|| {
            std::env::var("XDG_DATA_HOME")
                .ok()
                .filter(|v| !v.trim().is_empty())
                .map(|x| std::path::PathBuf::from(x).join("CraftBio"))
        })
        .or_else(|| {
            std::env::var("HOME")
                .ok()
                .filter(|v| !v.trim().is_empty())
                .map(|h| std::path::PathBuf::from(h).join(".local/share/CraftBio"))
        })
        .unwrap_or_else(|| std::path::PathBuf::from("CraftBio"));
    let custom_dir = data_dir.join("monomers");
    cat.load_custom_dir(&custom_dir);
    cat
}

#[cfg(test)]
mod tests {
    use super::*;
    use mol_graph::{analyze_molecule, write_pdb, write_sdf_v2000_3d};

    #[test]
    fn all_8_landmark_nrp_presets_match_exact_literature_formula_and_mass()
    -> Result<(), PeptideError> {
        let cat = MonomerCatalog::builtin();
        assert!(cat.total_stereomonomer_count() >= 90);

        for preset in landmark_nrp_presets() {
            let built = build_peptide_or_nrp(preset.id, PeptideConformation::MacrocycleRing, &cat)?;
            assert!(
                built.mol.warnings.is_empty(),
                "preset {} emitted valence warnings: {:?}",
                preset.id,
                built.mol.warnings
            );
            assert_eq!(
                built.mol.formula(),
                preset.expected_formula,
                "formula mismatch for {} (got {}, expected {}, SMILES={})",
                preset.name,
                built.mol.formula(),
                preset.expected_formula,
                built.canonical_smiles
            );
            let mass = analyze_molecule(&built.mol).exact_mass;
            assert!(
                (mass - preset.expected_monoisotopic_mass).abs() < 0.005,
                "monoisotopic mass mismatch for {}: got {mass:.4}, expected {:.4}",
                preset.name,
                preset.expected_monoisotopic_mass
            );
            assert!(
                built
                    .mol
                    .atoms
                    .iter()
                    .all(|a| a.coords_2d.is_some() && a.coords_3d.is_some()),
                "preset {} missing 2D or 3D coords",
                preset.name
            );
            // Verify 3D SDF V2000/V3000 & PDB export cleanly
            let pdb = write_pdb(&built.mol);
            assert!(pdb.contains("CONECT") && pdb.contains("END"));
            if preset.id == "octreotide" {
                assert!(
                    pdb.contains("SSBOND "),
                    "Octreotide PDB missing SSBOND record"
                );
            }
            let sdf = write_sdf_v2000_3d(&[&built.mol]);
            assert!(sdf.contains("M  END"));
        }
        Ok(())
    }

    #[test]
    fn l_vs_d_stereochemistry_helm_roundtrip_and_rcm_click_staples() -> Result<(), PeptideError> {
        let cat = MonomerCatalog::builtin();

        // L-Ala vs D-Ala CIP & structural key distinction
        let l_ala = build_peptide_or_nrp("Ala", PeptideConformation::BetaStrand, &cat)?;
        let d_ala = build_peptide_or_nrp("D-Ala", PeptideConformation::BetaStrand, &cat)?;
        let cip_l: Vec<char> = l_ala.mol.atoms.iter().filter_map(|a| a.cip_label).collect();
        let cip_d: Vec<char> = d_ala.mol.atoms.iter().filter_map(|a| a.cip_label).collect();
        assert_eq!(cip_l, vec!['S']);
        assert_eq!(cip_d, vec!['R']);
        assert_ne!(l_ala.canonical_smiles, d_ala.canonical_smiles);

        // HELM round-trip on Cilengitide and Octreotide
        for id in ["cilengitide", "octreotide", "polymyxin_b1", "daptomycin"] {
            let b1 = build_peptide_or_nrp(id, PeptideConformation::MacrocycleRing, &cat)?;
            let spec2 = parse_helm(&b1.helm, &cat, PeptideConformation::MacrocycleRing)?;
            let b2 = build_peptide_from_spec(&spec2, &cat)?;
            assert_eq!(
                b1.mol.formula(),
                b2.mol.formula(),
                "HELM round-trip altered formula for {id}: {}",
                b1.helm
            );
            assert_eq!(
                b1.mol.structural_key(),
                b2.mol.structural_key(),
                "HELM round-trip altered structure for {id}: {}",
                b1.helm
            );
        }

        // RCM all-hydrocarbon i,i+4 staple loses exact C2H4 (28.0313 Da)
        let unstapled = build_peptide_or_nrp(
            "Ac-Ala-S5-Glu-Leu-Lys-S5-Ala-NH2",
            PeptideConformation::AlphaHelix,
            &cat,
        )?;
        let stapled = build_peptide_or_nrp(
            "Ac-Ala-S5-Glu-Leu-Lys-S5-Ala-NH2 | staple(2-6)",
            PeptideConformation::AlphaHelix,
            &cat,
        )?;
        assert!(
            stapled.mol.warnings.is_empty(),
            "{:?}",
            stapled.mol.warnings
        );
        let dm =
            analyze_molecule(&unstapled.mol).exact_mass - analyze_molecule(&stapled.mol).exact_mass;
        assert!(
            (dm - 28.0313).abs() < 0.003,
            "RCM staple should lose C2H4 (28.0313 Da), got {dm:.4}"
        );

        // CuAAC click staple preserves exact formula (isomerization alkyne+azide -> 1,2,3-triazole)
        let unclicked = build_peptide_or_nrp(
            "Ac-Ala-Pra-Gly-Ala-Aha-NH2",
            PeptideConformation::AlphaHelix,
            &cat,
        )?;
        let clicked = build_peptide_or_nrp(
            "Ac-Ala-Pra-Gly-Ala-Aha-NH2 | click(2-5)",
            PeptideConformation::AlphaHelix,
            &cat,
        )?;
        assert!(
            clicked.mol.warnings.is_empty(),
            "{:?}",
            clicked.mol.warnings
        );
        assert_eq!(unclicked.mol.formula(), clicked.mol.formula());

        // Depsipeptide ester linkage (Valinomycin subunit)
        let depsi = build_peptide_or_nrp(
            "cyclo(D-Val-Lac-Val-D-Hiv)",
            PeptideConformation::MacrocycleRing,
            &cat,
        )?;
        assert!(depsi.mol.warnings.is_empty(), "{:?}", depsi.mol.warnings);
        assert_eq!(depsi.mol.formula(), "C18H30N2O6");
        Ok(())
    }

    #[test]
    fn nrps_bgc_stachelhaus_antismash_gbk_and_domain_swap_work() -> Result<(), PeptideError> {
        let cat = MonomerCatalog::builtin();

        // 1. Stachelhaus probabilistic predictor
        let pred = predict_a_domain_specificity("DAWTIAAVCK");
        assert_eq!(pred.candidates[0].monomer_code, "Phe");
        assert!((pred.candidates[0].identity_pct - 100.0).abs() < 0.1);

        // 2. Compile Gramicidin S, ACV, Polymyxin B1, Daptomycin BGCs and verify exact natural product formulas
        for (bgc_id, expected_formula) in [
            ("gramicidin_s", "C60H92N12O10"),
            ("acv", "C14H25N3O6S"),
            ("polymyxin_b1", "C56H98N16O13"),
            ("daptomycin", "C72H101N17O26"),
        ] {
            let bgc = preset_nrps_bgc(bgc_id)
                .ok_or_else(|| PeptideError::Nrps(format!("missing preset {bgc_id}")))?;
            let (rep, built) = compile_nrps_bgc(&bgc, &cat)?;
            assert!(
                rep.collinear_compatible,
                "unexpected interface warning on {bgc_id}: {:?}",
                rep.interface_warnings
            );
            assert_eq!(built.mol.formula(), expected_formula);
        }

        // 3. Detect E-domain -> C_LCL interface mismatch cleanly
        let mut broken = preset_nrps_bgc("gramicidin_s")
            .ok_or_else(|| PeptideError::Nrps("missing grs".to_string()))?;
        broken.modules[1].domains[0] = NrpsDomain::CLcl;
        let (broken_rep, _) = compile_nrps_bgc(&broken, &cat)?;
        assert!(!broken_rep.collinear_compatible);
        assert!(
            broken_rep
                .interface_warnings
                .iter()
                .any(|w| w.contains("C_DCL"))
        );

        // 4. Domain swap (Gramicidin S Mod 3 Val -> Ile)
        let grs = preset_nrps_bgc("gramicidin_s")
            .ok_or_else(|| PeptideError::Nrps("missing grs".to_string()))?;
        let (swap, eng_built) = propose_nrps_swap(&grs, 3, "Ile", &cat)?;
        // Homodimer adds 2x CH2 (+28.0313 Da)
        assert!(
            (swap.delta_mw - 28.0313).abs() < 0.01,
            "delta_mw={}",
            swap.delta_mw
        );
        assert_eq!(eng_built.mol.formula(), "C62H96N12O10");

        // 5. Parse lightweight antiSMASH GenBank (.gbk) snippet
        let gbk = "\
LOCUS       BGC0000367_GRSA_TEST   12000 bp    DNA     linear   BCT
ORGANISM    Aneurinibacillus migulanus
FEATURES             Location/Qualifiers
     CDS             1..3000
     aSDomain        100..1200
                     /aSDomain=\"AMP-binding\"
                     /specificity=\"Stachelhaus code: DAWTIAAVCK\"
     aSDomain        1300..1500
                     /aSDomain=\"PCP\"
     aSDomain        1600..2800
                     /aSDomain=\"Epimerization\"
     CDS             3001..9000
     aSDomain        3100..4000
                     /aSDomain=\"Condensation_DCL\"
     aSDomain        4100..5200
                     /aSDomain=\"AMP-binding\"
                     /specificity=\"Stachelhaus code: DVQFIAHVVK\"
     aSDomain        5300..5500
                     /aSDomain=\"PCP\"
     aSDomain        5600..6500
                     /aSDomain=\"Condensation_LCL\"
     aSDomain        6600..7700
                     /aSDomain=\"AMP-binding\"
                     /specificity=\"Stachelhaus code: DAFWIGGTFK\"
     aSDomain        7800..8000
                     /aSDomain=\"PCP\"
     aSDomain        8100..8800
                     /aSDomain=\"Thioesterase\"
";
        let imported = parse_antismash_gbk(gbk)?;
        assert_eq!(imported.modules.len(), 3);
        let (imp_rep, imp_built) = compile_nrps_bgc(&imported, &cat)?;
        assert!(imp_rep.collinear_compatible);
        assert_eq!(imp_rep.module_audits[0].resolved_monomer, "Phe");
        assert!(imp_rep.module_audits[0].is_d);
        assert_eq!(imp_rep.module_audits[1].resolved_monomer, "Pro");
        assert_eq!(imp_rep.module_audits[2].resolved_monomer, "Val");
        assert!(imp_built.mol.warnings.is_empty());
        Ok(())
    }

    #[test]
    fn developability_dpp4_aib2_shielding_bro5_spps_and_scans_work() -> Result<(), PeptideError> {
        let cat = MonomerCatalog::builtin();

        // Native GLP-1 (Ala2) is vulnerable to DPP-4; Semaglutide core (Aib2) shields DPP-4!
        let glp1_wt = build_peptide_or_nrp(
            "His-Ala-Glu-Gly-Thr-Phe-Thr-Ser-Asp-Val-NH2",
            PeptideConformation::AlphaHelix,
            &cat,
        )?;
        let sema = build_peptide_or_nrp("semaglutide_core", PeptideConformation::AlphaHelix, &cat)?;
        let dev_wt = analyze_peptide_developability(&glp1_wt, &cat);
        let dev_sema = analyze_peptide_developability(&sema, &cat);
        assert!(
            dev_wt
                .protease_sites
                .iter()
                .any(|s| s.protease.contains("DPP-4") && !s.shielded)
        );
        assert!(
            dev_sema
                .protease_sites
                .iter()
                .any(|s| s.protease.contains("DPP-4") && s.shielded)
        );

        // Cyclosporin A scores strongly in bRo5 Oral Macrocycle space vs linear polar Polymyxin B1
        let csa = build_peptide_or_nrp("cyclosporin_a", PeptideConformation::MacrocycleRing, &cat)?;
        let pmx = build_peptide_or_nrp("polymyxin_b1", PeptideConformation::MacrocycleRing, &cat)?;
        let dev_csa = analyze_peptide_developability(&csa, &cat);
        let dev_pmx = analyze_peptide_developability(&pmx, &cat);
        assert!(
            dev_csa.bro5_oral_score > dev_pmx.bro5_oral_score + 0.30,
            "Cyclosporin A ({}) should strongly beat cationic Polymyxin B1 ({}) in bRo5 oral score",
            dev_csa.bro5_oral_score,
            dev_pmx.bro5_oral_score
        );
        assert!(dev_pmx.charge_ph74 > 4.0);
        assert!(
            dev_pmx
                .liabilities
                .iter()
                .any(|l| l.category.contains("Hemolysis"))
        );

        // SPPS plan & N-Methyl / D-AA scan on Cilengitide / linear pentapeptide
        let spps_csa = plan_spps_vs_nrps(&csa.spec, &cat);
        assert!(spps_csa.recommended_route.contains("NRPS"));
        let lin = parse_peptide_spec(
            "Ac-Arg-Gly-Asp-Phe-Val-NH2",
            &cat,
            PeptideConformation::BetaStrand,
        )?;
        let d_scan = run_peptide_scan(&lin, ScanMode::DAminoAcid, &cat)?;
        assert_eq!(d_scan.variants.len(), 4); // Gly is achiral, 4 chiral residues scanned
        let nme_scan = run_peptide_scan(&lin, ScanMode::NMethyl, &cat)?;
        assert_eq!(nme_scan.variants.len(), 5);
        assert!(render_scan_csv(&nme_scan).lines().count() == 6);
        Ok(())
    }
}
