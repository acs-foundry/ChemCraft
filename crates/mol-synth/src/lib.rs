#![forbid(unsafe_code)]

//! `mol-synth`: Pure-Rust Buy-Before-Make building-block lookup, graph-based MedChem
//! retrosynthetic disconnection engine, Retro* AND-OR multi-step route planner, stoichiometry
//! calculator, MDL `.rxn` exporter, and ELN lab protocol generator for ChemCraft's on-demand
//! `Synth Module`.

pub mod data_dir;
pub mod enumerate;
pub mod retro_star;
pub mod stock;
pub mod templates;

use mol_graph::{MolError, parse_smiles};

pub use data_dir::{building_blocks_dir, craft_data_dir, default_stock_tsv_path};
pub use enumerate::{
    E3RecruiterPreset, EnumerateMode, EnumeratedCompound, EnumeratedLibrary,
    enumerate_protac_library, enumerate_rgroup_library, shift_ring_indices,
};
pub use retro_star::{
    RetroConfig, StartingMaterial, SynthesisPlan, SynthesisRoute, SynthesisStep,
    plan_retrosynthesis,
};
pub use stock::{
    StockAnalog, StockCatalog, StockEntry, StockMatch, StockTier, default_stock_catalog,
};
pub use templates::{Disconnection, enumerate_disconnections};

/// Convenience entry point: parse SMILES and plan multi-step retrosynthesis against the default
/// local building-block catalog (`CRAFT_DATA_DIR/building_blocks` + built-in MedChem core).
pub fn plan_retrosynthesis_smiles(
    smiles: &str,
    name: Option<&str>,
    cfg: &RetroConfig,
) -> Result<SynthesisPlan, MolError> {
    let mut mol = parse_smiles(smiles)?;
    if let Some(n) = name {
        mol.name = Some(n.to_string());
    }
    Ok(plan_retrosynthesis(&mol, default_stock_catalog(), cfg))
}

#[cfg(test)]
mod tests {
    use super::*;
    use mol_graph::parse_rxn_v2000;

    #[test]
    fn solves_real_medchem_targets_to_buyable_stock_and_exports_rxn()
    -> Result<(), Box<dyn std::error::Error>> {
        let stock = StockCatalog::builtin_medchem();
        let cfg = RetroConfig::default();

        // 1. Paracetamol (1-step amide coupling from 4-aminophenol + acetic acid, both in stock!)
        let para = parse_smiles("CC(=O)Nc1ccc(O)cc1")?;
        let plan_para = plan_retrosynthesis(&para, &stock, &cfg);
        assert!(!plan_para.routes.is_empty());
        assert!(plan_para.routes[0].fully_resolved_in_stock);
        let rxn_str = plan_para.routes[0].steps[0]
            .to_rxn_v2000()
            .ok_or("expected valid RXN V2000")?;
        let parsed_rxn = parse_rxn_v2000(&rxn_str)?;
        assert_eq!(parsed_rxn.products.len(), 1);
        assert!(parsed_rxn.reactants.len() >= 2);

        // 2. Kinase-style 2-aminopyrimidine biaryl inhibitor (multi-step SNAr + Suzuki!)
        let kinase = parse_smiles("Cc1ccc(-c2ccnc(Nc3ccc(S(N)(=O)=O)cc3)n2)cc1")?;
        let plan_kin = plan_retrosynthesis(&kinase, &stock, &cfg);
        assert!(!plan_kin.routes.is_empty());
        assert!(
            plan_kin.routes[0].total_steps >= 2 && plan_kin.routes[0].fully_resolved_in_stock,
            "expected multi-step fully resolved kinase route, got {:?}",
            plan_kin.routes[0]
        );

        // 3. Macrocyclic CDK2 Rank #1 hit (ACT_CHEMBL340259): intramolecular N-alkylation macrocyclization -> Arcyriarubin A + PEG dibromide!
        let cdk2_hit = parse_smiles("c1ccc2c(c1)c1cn2CCOCCOCCOCCn2cc(c3ccccc32)c2c1c([nH]c2=O)=O")?;
        let plan_cdk2 = plan_retrosynthesis(&cdk2_hit, &stock, &cfg);
        assert!(
            !plan_cdk2.routes.is_empty(),
            "expected retrosynthesis routes for macrocyclic CDK2 hit"
        );
        assert!(
            plan_cdk2.routes.iter().any(|r| r.fully_resolved_in_stock),
            "expected a fully stock-resolved macrocyclization route for ACT_CHEMBL340259, got {:?}",
            plan_cdk2.routes
        );

        // 4. Tetrazole bioisostere of Ibuprofen -> nitrile cycloaddition
        let ibu_tet = parse_smiles("CC(C)Cc1ccc(C(C)c2nnn[nH]2)cc1")?;
        let plan_tet = plan_retrosynthesis(&ibu_tet, &stock, &cfg);
        assert!(
            plan_tet
                .routes
                .iter()
                .any(|r| r.steps.iter().any(|s| s.template_id == "TETRAZOLE_CLICK"))
        );

        Ok(())
    }
}
