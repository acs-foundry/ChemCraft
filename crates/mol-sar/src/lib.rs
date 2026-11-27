//! Layer 2 (`mol-sar`): High-performance in-memory DataWarrior-style SAR engine:
//! ECFP4 Morgan fingerprints (stereo-aware), parallel Tanimoto & SMARTS filters,
//! reproducible anchored Matched Molecular Pairs (MMP), deduplicated Activity Cliffs (SALI),
//! Ligand Efficiency metrics, and AdmetCraft bridge enrichment.

#![forbid(unsafe_code)]

pub mod dbtl;
pub mod fp;
pub mod mmp;
pub mod table;

pub use dbtl::{
    DbtlCompoundSummary, DbtlIngestReport, ingest_experimental_assay_text, render_dbtl_sar_csv,
    render_next_round_seeds_smi, render_qsar_recalibration_csv,
};
pub use fp::MorganFp;
pub use mmp::{
    ActivityCliff, DuplicatePair, MatchedMolecularPair, find_activity_cliffs,
    find_duplicate_compounds, find_matched_molecular_pairs,
};
pub use table::{SarFilter, SarRow, SarTable};

#[cfg(test)]
mod tests {
    use super::*;
    use mol_admet_bridge::AdmetCraftProvider;
    use mol_graph::{MolError, parse_smiles};

    #[test]
    fn ecfp4_tanimoto_substructure_mmp_and_cliffs_work_on_50_compounds() {
        let mut table = SarTable::demo_50_compounds();
        assert_eq!(table.rows.len(), 50);

        table.enrich_with_admet(&AdmetCraftProvider::default());
        assert!(table.rows.iter().all(|r| r.admet.is_some()));

        let filter = SarFilter {
            smarts_substructure: Some("c1ccccc1".to_string()),
            ..SarFilter::default()
        };
        let aromatic_hits = table.matching_indices(&filter);
        assert!(aromatic_hits.is_ok());
        let hits = aromatic_hits.unwrap_or_default();
        assert!(hits.len() >= 45);

        let cliffs = find_activity_cliffs(&table.rows, 0.35);
        assert!(!cliffs.is_empty());

        // MMP must be 100% deterministic across repeated invocations
        let mmps_1 = find_matched_molecular_pairs(&table.rows);
        let mmps_2 = find_matched_molecular_pairs(&table.rows);
        let mmps_3 = find_matched_molecular_pairs(&table.rows);
        assert!(!mmps_1.is_empty());
        assert_eq!(mmps_1, mmps_2);
        assert_eq!(mmps_2, mmps_3);

        if let (Ok(m1), Ok(m2)) = (
            parse_smiles("CC(C)Cc1ccc(C(C)C(=O)O)cc1"),
            parse_smiles("CC(C)Cc1ccc(C(C)c2nnn[nH]2)cc1"),
        ) {
            let sim = MorganFp::ecfp4(&m1).tanimoto(&MorganFp::ecfp4(&m2));
            assert!(
                sim > 0.25 && sim < 1.0,
                "expected moderate ECFP4 similarity for bioisostere pair, got {sim}"
            );
            assert!((MorganFp::ecfp4(&m1).tanimoto(&MorganFp::ecfp4(&m1)) - 1.0).abs() < 1e-12);
        }
    }

    #[test]
    fn sar_separates_enantiomers_deduplicates_kekule_and_avoids_pos_isomer_mmp()
    -> Result<(), MolError> {
        let l_ala = parse_smiles("N[C@@H](C)C(=O)O")?;
        let d_ala = parse_smiles("N[C@H](C)C(=O)O")?;
        let sim_en = MorganFp::ecfp4(&l_ala).tanimoto(&MorganFp::ecfp4(&d_ala));
        assert!(
            sim_en < 1.0,
            "enantiomers must not have Tanimoto=1.0, got {sim_en}"
        );

        // Duplicate Kekule vs Aromatic naphthalene with discordant activity MUST be detected & excluded from SALI cliffs
        let mut table = SarTable::new();
        let mut n1 = parse_smiles("c1ccc2ccccc2c1")?;
        n1.name = Some("Naph_Arom".to_string());
        let mut n2 = parse_smiles("C1=CC=C2C=CC=CC2=C1")?;
        n2.name = Some("Naph_Kekule".to_string());
        table.add_molecule(n1, Some(5.0), None);
        table.add_molecule(n2, Some(7.5), None);

        let dups = find_duplicate_compounds(&table.rows);
        assert_eq!(dups.len(), 1);
        assert!(find_activity_cliffs(&table.rows, 0.5).is_empty());

        // Position isomers (2-chlorotoluene vs 4-chlorotoluene) must NOT form a fake core MMP
        let mut iso = SarTable::new();
        let mut o_cl = parse_smiles("Cc1ccccc1Cl")?;
        o_cl.name = Some("2-Chlorotoluene".to_string());
        let mut p_cl = parse_smiles("Cc1ccc(Cl)cc1")?;
        p_cl.name = Some("4-Chlorotoluene".to_string());
        iso.add_molecule(o_cl, Some(6.0), None);
        iso.add_molecule(p_cl, Some(6.1), None);
        let iso_mmps = find_matched_molecular_pairs(&iso.rows);
        assert!(
            iso_mmps.is_empty(),
            "position isomers falsely paired: {iso_mmps:?}"
        );

        // Toluene vs Phenol uses anchored core `[*]c1ccccc1` and `[*]C` -> `[*]O` (not methane -> water)
        let mut tp = SarTable::new();
        let mut tol = parse_smiles("Cc1ccccc1")?;
        tol.name = Some("Toluene".to_string());
        let mut phe = parse_smiles("Oc1ccccc1")?;
        phe.name = Some("Phenol".to_string());
        tp.add_molecule(tol, Some(4.5), None);
        tp.add_molecule(phe, Some(5.5), None);
        let tp_mmps = find_matched_molecular_pairs(&tp.rows);
        assert_eq!(tp_mmps.len(), 1);
        assert!(tp_mmps[0].core_smiles.contains('*'));
        assert!(tp_mmps[0].r_group_a.contains('*') && tp_mmps[0].r_group_b.contains('*'));
        Ok(())
    }
}
