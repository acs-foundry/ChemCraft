//! Layer 0 (`mol-graph`): Clean-room editable molecular graph, iterative Hückel aromaticity/SSSR,
//! CIP R/S & E/Z stereochemistry, SMILES/SMARTS, SDF V2000/V3000, DataWarrior `.dwar`, CSV, and mutations.

#![forbid(unsafe_code)]

pub mod analysis;
pub mod cdxml;
pub mod coords3d;
pub mod dwar;
pub mod element;
pub mod error;
pub mod groups;
pub mod mol;
pub mod mutate;
pub mod names;
pub mod peptide;
pub mod rxn;
pub mod sdf;
pub mod smarts;
pub mod smiles;

pub use analysis::{
    ElementPercent, MassAdduct, MolecularAnalysis, analyze_molecule, monoisotopic_mass,
};
pub use cdxml::{parse_cdxml, write_cdxml};
pub use coords3d::{
    compute_3d_coords, covalent_radius, ensure_3d_coords, generate_3d_coords, ideal_bond_length,
    place_nerf_atom, relax_3d_coords_pinned, write_pdb, write_sdf_v2000_3d,
};
pub use dwar::{
    TabularMolRecord, parse_csv_records, parse_csv_records_with_limit, parse_dwar,
    parse_dwar_with_limit, write_dwar,
};
pub use element::{ElementMeta, element_by_number, element_by_symbol, max_neutral_valence};
pub use error::MolError;
pub use groups::{
    SUPERATOM_GROUPS, SuperatomGroup, attach_superatom_group, build_superatom_fragment,
    find_superatom_group,
};
pub use mol::{
    BondCisTrans, BondOrder, BondStereo, ChiralTag, MolAtom, MolBond, MolGraph, PdbAtomMeta,
    highest_cip_end_nbr,
};
pub use mutate::{
    SubstituentKind, attach_fragment_at, mutate_molecule, replace_carboxyl_with_tetrazole,
};
pub use names::{known_chemical_names, resolve_chemical_name};
pub use peptide::{
    MAX_PEPTIDE_RESIDUES, ParsedPeptide, SecondaryStructure, build_peptide, canonical_amino_acids,
    parse_peptide_sequence,
};
pub use rxn::{
    ReactionRole, ReactionScheme, ReactionStoichiometry, StoichRow, compute_stoichiometry,
    parse_reaction_smiles, parse_rxn_v2000, write_rxn_v2000,
};
pub use sdf::{
    parse_mol_block, parse_sdf, parse_sdf_with_limit, write_mol_v2000_block, write_sdf_v2000,
    write_sdf_v3000,
};
pub use smarts::SmartsQuery;
pub use smiles::{DEFAULT_MAX_ATOMS, parse_smiles, parse_smiles_line, parse_smiles_with_limit};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_aspirin_steroid_and_bioisostere_tetrazole() -> Result<(), MolError> {
        let aspirin = parse_smiles("CC(=O)Oc1ccccc1C(=O)O")?;
        assert_eq!(aspirin.formula(), "C9H8O4");
        assert_eq!(aspirin.aromatic_atom_count(), 6);

        let tetrazole_analogue = replace_carboxyl_with_tetrazole(&aspirin)?;
        assert_eq!(tetrazole_analogue.formula(), "C9H8N4O2");
        assert_ne!(
            tetrazole_analogue.structural_key(),
            aspirin.structural_key()
        );

        let ibuprofen = parse_smiles("CC(C)Cc1ccc(C(C)C(=O)O)cc1")?;
        let ibu_tetrazole = mutate_molecule(&ibuprofen, None, "tetrazole")?;
        let smi_q = SmartsQuery::parse("c1nnn[nH]1")?;
        assert!(smi_q.matches(&ibu_tetrazole));
        Ok(())
    }

    #[test]
    fn kekule_equivalence_fused_valence_pyrrole_acetate_ferrocene_and_macrocycle()
    -> Result<(), MolError> {
        // 1. Naphthalene & Caffeine must have 0 valence errors and Kekule == Aromatic structural_key
        let naph_arom = parse_smiles("c1ccc2ccccc2c1")?;
        let naph_kek = parse_smiles("C1=CC=C2C=CC=CC2=C1")?;
        assert!(naph_arom.warnings.is_empty(), "{:?}", naph_arom.warnings);
        assert!(naph_kek.warnings.is_empty(), "{:?}", naph_kek.warnings);
        assert_eq!(naph_arom.structural_key(), naph_kek.structural_key());

        let caffeine = parse_smiles("Cn1cnc2c1c(=O)n(C)c(=O)n2C")?;
        assert!(caffeine.warnings.is_empty(), "{:?}", caffeine.warnings);

        // 2. Pyrrole NH preserved in both aromatic and Kekule forms (C4H5N)
        let pyr_arom = parse_smiles("c1cc[nH]c1")?;
        let pyr_kek = parse_smiles("C1=CC=CN1")?;
        assert_eq!(pyr_arom.formula(), "C4H5N");
        assert_eq!(pyr_kek.formula(), "C4H5N");
        assert_eq!(pyr_arom.structural_key(), pyr_kek.structural_key());

        // 3. Carboxylate read back from SDF keeps exact 59.04 Da (no phantom H on O-)
        let acetate = parse_smiles("CC(=O)[O-]")?;
        assert!((acetate.molecular_weight() - 59.044).abs() < 0.05);
        let sdf = write_sdf_v2000(&[&acetate]);
        let back = parse_sdf(&sdf)
            .into_iter()
            .next()
            .ok_or(MolError::EmptyInput)??;
        assert!((back.molecular_weight() - 59.044).abs() < 0.05);
        assert_eq!(back.formula(), "C2H3O2");

        // 4. Ferrocene keeps Iron (Fe)
        let ferrocene = parse_smiles("[cH-]1cccc1.[cH-]1cccc1.[Fe+2]")?;
        assert_eq!(ferrocene.formula(), "C10H10Fe");

        // 5. Macrocycle >16 atoms detected by SSSR
        let c20 = parse_smiles("C1CCCCCCCCCCCCCCCCCCC1")?;
        assert_eq!(c20.rings.len(), 1);
        assert_eq!(c20.rings[0].len(), 20);
        Ok(())
    }

    #[test]
    fn cip_rs_and_ez_stereochemistry_and_canonical_distinction() -> Result<(), MolError> {
        let l_ala = parse_smiles("N[C@@H](C)C(=O)O")?;
        let l_ala_alt = parse_smiles("C[C@H](N)C(=O)O")?;
        let d_ala = parse_smiles("N[C@H](C)C(=O)O")?;

        assert_eq!(l_ala.atoms[1].cip_label, Some('S'));
        assert_eq!(l_ala_alt.atoms[1].cip_label, Some('S'));
        assert_eq!(d_ala.atoms[1].cip_label, Some('R'));
        assert_ne!(l_ala.structural_key(), d_ala.structural_key());

        let e_alkene = parse_smiles("F/C=C/Cl")?;
        let z_alkene = parse_smiles("F/C=C\\Cl")?;
        let db_e = e_alkene
            .bonds
            .iter()
            .find(|b| b.order == BondOrder::Double)
            .ok_or(MolError::EmptyInput)?;
        let db_z = z_alkene
            .bonds
            .iter()
            .find(|b| b.order == BondOrder::Double)
            .ok_or(MolError::EmptyInput)?;
        assert_eq!(db_e.cis_trans, BondCisTrans::E);
        assert_eq!(db_z.cis_trans, BondCisTrans::Z);
        assert_ne!(e_alkene.structural_key(), z_alkene.structural_key());

        // Multi-sphere CIP beyond sphere 1: L-cysteine = (R) despite S having lower degree-sum than COOH O's
        let l_cys = parse_smiles("N[C@@H](CS)C(=O)O")?;
        let d_cys = parse_smiles("N[C@H](CS)C(=O)O")?;
        let l_cys_c = l_cys
            .atoms
            .iter()
            .find(|a| a.cip_label.is_some())
            .ok_or(MolError::EmptyInput)?;
        let d_cys_c = d_cys
            .atoms
            .iter()
            .find(|a| a.cip_label.is_some())
            .ok_or(MolError::EmptyInput)?;
        assert_eq!(l_cys_c.cip_label, Some('R'));
        assert_eq!(d_cys_c.cip_label, Some('S'));

        // Idempotent canonical SMILES & DWAR round-trip across 6 chiral / E-Z cases including Penicillin G
        let pen_g = "CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)C(=O)N2[C@H]1C(=O)O";
        for smi in [
            "N[C@@H](C)C(=O)O",
            "N[C@H](C)C(=O)O",
            "N[C@@H](CS)C(=O)O",
            pen_g,
            "C[C@H](O)[C@@H](O)C(=O)O",
            "F/C(Cl)=C(/Br)I",
            "F/C(Cl)=C(\\Br)I",
        ] {
            let m1 = parse_smiles(smi)?;
            let c1 = m1.canonical_smiles();
            let m2 = parse_smiles(&c1)?;
            let c2 = m2.canonical_smiles();
            assert_eq!(c1, c2, "canonical_smiles not idempotent for {smi}");
            assert_eq!(
                m1.structural_key(),
                m2.structural_key(),
                "structural_key oscillated for {smi}"
            );
            let mut cip1: Vec<char> = m1.atoms.iter().filter_map(|a| a.cip_label).collect();
            let mut cip2: Vec<char> = m2.atoms.iter().filter_map(|a| a.cip_label).collect();
            cip1.sort_unstable();
            cip2.sort_unstable();
            assert_eq!(cip1, cip2, "CIP label multiset flipped for {smi}");
        }

        // Different input orderings of L-alanine must give IDENTICAL canonical SMILES
        assert_eq!(l_ala.canonical_smiles(), l_ala_alt.canonical_smiles());

        // Penicillin G full .dwar write -> reload -> write round-trip preserves exact stereocenters (C5=R, C6=R, C2=S)
        let mut pen_mol = parse_smiles(pen_g)?;
        pen_mol.name = Some("Penicillin_G".to_string());
        let stereo_by_role = |mol: &MolGraph| -> Vec<(&'static str, char)> {
            let mut out = Vec::new();
            for (i, a) in mol.atoms.iter().enumerate() {
                if let Some(cip) = a.cip_label {
                    let has_s = mol.adj[i]
                        .iter()
                        .any(|&(nb, _)| mol.atoms[nb].symbol == "S");
                    let has_amide_n = mol.adj[i]
                        .iter()
                        .any(|&(nb, _)| mol.atoms[nb].symbol == "N" && mol.adj[nb].len() == 2);
                    let role = if has_s {
                        "C5_bridgehead"
                    } else if has_amide_n {
                        "C6_lactam"
                    } else {
                        "C2_carboxyl"
                    };
                    out.push((role, cip));
                }
            }
            out.sort();
            out
        };
        let expected_roles = vec![
            ("C2_carboxyl", 'S'),
            ("C5_bridgehead", 'R'),
            ("C6_lactam", 'R'),
        ];
        assert_eq!(stereo_by_role(&pen_mol), expected_roles);
        let dwar1 = write_dwar(&[(&pen_mol, &pen_mol.properties)]);
        let reloaded = parse_dwar(&dwar1)?
            .into_iter()
            .next()
            .ok_or(MolError::EmptyInput)??;
        assert_eq!(stereo_by_role(&reloaded.mol), expected_roles);
        assert_eq!(reloaded.mol.structural_key(), pen_mol.structural_key());
        let dwar2 = write_dwar(&[(&reloaded.mol, &reloaded.columns)]);
        assert_eq!(dwar1, dwar2);

        // Conjugated polyenes: hexa-2,4-diene -> 3 unique isomers; unsymmetric hepta-2,4-diene -> 4 unique isomers!
        let hex_keys: std::collections::BTreeSet<String> = [
            "C/C=C/C=C/C",
            "C/C=C/C=C\\C",
            "C/C=C\\C=C/C",
            "C/C=C\\C=C\\C",
        ]
        .iter()
        .map(|s| {
            let m = parse_smiles(s)?;
            let c1 = m.canonical_smiles();
            let m2 = parse_smiles(&c1)?;
            assert_eq!(
                c1,
                m2.canonical_smiles(),
                "hexa-2,4-diene not idempotent for {s}"
            );
            Ok(m.structural_key())
        })
        .collect::<Result<_, MolError>>()?;
        assert_eq!(
            hex_keys.len(),
            3,
            "hexa-2,4-diene must produce 3 unique keys (EE, EZ=ZE, ZZ)"
        );

        let hep_keys: std::collections::BTreeSet<String> = [
            "C/C=C/C=C/CC",
            "C/C=C/C=C\\CC",
            "C/C=C\\C=C/CC",
            "C/C=C\\C=C\\CC",
        ]
        .iter()
        .map(|s| {
            let m = parse_smiles(s)?;
            let c1 = m.canonical_smiles();
            let m2 = parse_smiles(&c1)?;
            assert_eq!(
                c1,
                m2.canonical_smiles(),
                "hepta-2,4-diene not idempotent for {s}"
            );
            Ok(m.structural_key())
        })
        .collect::<Result<_, MolError>>()?;
        assert_eq!(
            hep_keys.len(),
            4,
            "hepta-2,4-diene must preserve all 4 E/Z isomers"
        );

        // Trisubstituted conjugated polyenes (Tretinoin, Isotretinoin, 11-cis-Retinal, Lycopene)
        // must preserve EVERY stereogenic double bond across canonical SMILES emission!
        let count_ez = |m: &MolGraph| -> usize {
            m.bonds
                .iter()
                .filter(|b| matches!(b.cis_trans, BondCisTrans::E | BondCisTrans::Z))
                .count()
        };
        let tretinoin = "CC1=C(C(CCC1)(C)C)/C=C/C(=C/C=C/C(=C/C(=O)O)/C)/C";
        let isotretinoin = "CC1=C(C(CCC1)(C)C)/C=C/C(=C/C=C/C(=C\\C(=O)O)/C)/C";
        let retinal_11cis = "CC1=C(C(CCC1)(C)C)/C=C/C(=C/C=C\\C(=C/C=O)/C)/C";
        let lycopene = "CC(=CCC/C(=C/C=C/C(=C/C=C/C(=C/C=C/C=C(\\C)/C=C/C=C(/C)/C=C/C=C(/C)\\CC/C=C(/C)\\C)/C)/C)/C)C";
        for (name, smi) in [
            ("Tretinoin", tretinoin),
            ("Isotretinoin", isotretinoin),
            ("11-cis-Retinal", retinal_11cis),
            ("Lycopene", lycopene),
        ] {
            let m1 = parse_smiles(smi)?;
            let n1 = count_ez(&m1);
            let can = m1.canonical_smiles();
            let m2 = parse_smiles(&can)?;
            let n2 = count_ez(&m2);
            assert_eq!(
                n1, n2,
                "{name} lost E/Z geometry in canonical SMILES: {can}"
            );
            assert_eq!(
                can,
                m2.canonical_smiles(),
                "{name} canonical SMILES not idempotent"
            );
            assert_eq!(
                m1.structural_key(),
                m2.structural_key(),
                "{name} key changed"
            );
        }
        assert_ne!(
            parse_smiles(tretinoin)?.structural_key(),
            parse_smiles(isotretinoin)?.structural_key()
        );
        Ok(())
    }

    #[test]
    fn smarts_rejects_min_charge_and_honors_precedence_and_r2() -> Result<(), MolError> {
        // Must return Err cleanly instead of panicking on -128
        assert!(matches!(
            SmartsQuery::parse("[-128]"),
            Err(MolError::InvalidSmarts { .. })
        ));

        // [C,N;H1] = (C OR N) AND H1 -> matches aniline NH2? no, matches methanol OH? no, matches dimethylamine [NH] (H1)!
        let q = SmartsQuery::parse("[C,N;H1]")?;
        let nh = parse_smiles("CNC")?;
        let pyridine = parse_smiles("c1ccncc1")?;
        assert!(q.matches(&nh));
        assert!(!q.matches(&pyridine));

        // [R2] matches naphthalene bridgehead carbons (2 SSSR rings) and NOT benzene (1 ring)
        let r2 = SmartsQuery::parse("[R2]")?;
        let naph = parse_smiles("c1ccc2ccccc2c1")?;
        let benz = parse_smiles("c1ccccc1")?;
        assert!(r2.matches(&naph));
        assert!(!r2.matches(&benz));
        Ok(())
    }

    #[test]
    fn sdf_v2000_v3000_dwar_and_headerless_rejection() -> Result<(), MolError> {
        let mut caffeine = parse_smiles("Cn1cnc2c1c(=O)n(C)c(=O)n2C")?;
        caffeine.name = Some("Caffeine".to_string());
        caffeine
            .properties
            .insert("pIC50".to_string(), "6.42".to_string());

        let v2000 = write_sdf_v2000(&[&caffeine]);
        let m2 = parse_sdf(&v2000)
            .into_iter()
            .next()
            .ok_or(MolError::EmptyInput)??;
        assert_eq!(m2.formula(), "C8H10N4O2");
        assert_eq!(m2.properties.get("pIC50").map(String::as_str), Some("6.42"));

        let v3000 = write_sdf_v3000(&caffeine);
        let m3 = parse_sdf(&v3000)
            .into_iter()
            .next()
            .ok_or(MolError::EmptyInput)??;
        assert_eq!(m3.formula(), "C8H10N4O2");

        // Uracil, pyrrole, and indole must preserve their NH protons through V2000 & V3000 SDF
        for (smi, expected_formula) in [
            ("O=c1[nH]c(=O)[nH]cc1", "C4H4N2O2"),
            ("c1cc[nH]c1", "C4H5N"),
            ("c1ccc2[nH]ccc2c1", "C8H7N"),
        ] {
            let m = parse_smiles(smi)?;
            assert_eq!(m.formula(), expected_formula);
            let sdf2 = write_sdf_v2000(&[&m]);
            let r2 = parse_sdf(&sdf2)
                .into_iter()
                .next()
                .ok_or(MolError::EmptyInput)??;
            assert_eq!(r2.formula(), expected_formula, "V2000 lost NH for {smi}");
            let sdf3 = write_sdf_v3000(&m);
            let r3 = parse_sdf(&sdf3)
                .into_iter()
                .next()
                .ok_or(MolError::EmptyInput)??;
            assert_eq!(r3.formula(), expected_formula, "V3000 lost NH for {smi}");
        }

        // External aromatic bond-type-4 V2000 SDF of uracil WITHOUT explicit H atoms must keep C4H4N2O2
        let external_uracil_sdf = "\
uracil_type4
  Ext
  
  8  8  0  0  0  0  0  0  0  0999 V2000
    0.0000    1.5000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3000   -0.7500    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3000   -2.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6000   -3.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -3.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.3000   -2.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3000   -0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  4  0  0  0  0
  3  4  4  0  0  0  0
  4  5  2  0  0  0  0
  4  6  4  0  0  0  0
  6  7  4  0  0  0  0
  7  8  4  0  0  0  0
  8  2  4  0  0  0  0
M  END
$$$$
";
        let ext_u = parse_sdf(external_uracil_sdf)
            .into_iter()
            .next()
            .ok_or(MolError::EmptyInput)??;
        assert_eq!(ext_u.formula(), "C4H4N2O2");

        // Full R/S + E/Z SDF V2000 & V3000 interchange preservation
        for smi in [
            "N[C@@H](C)C(=O)O",
            "N[C@H](C)C(=O)O",
            "N[C@@H](CS)C(=O)O",
            "CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)C(=O)N2[C@H]1C(=O)O",
            "F/C=C/Cl",
            "F/C=C\\Cl",
            "C/C=C/C=C\\CC",
        ] {
            let orig = parse_smiles(smi)?;
            let expected_key = orig.structural_key();
            let sdf2 = write_sdf_v2000(&[&orig]);
            let back2 = parse_sdf(&sdf2)
                .into_iter()
                .next()
                .ok_or(MolError::EmptyInput)??;
            assert_eq!(
                back2.structural_key(),
                expected_key,
                "SDF V2000 lost stereo for {smi}"
            );
            let sdf3 = write_sdf_v3000(&orig);
            let back3 = parse_sdf(&sdf3)
                .into_iter()
                .next()
                .ok_or(MolError::EmptyInput)??;
            assert_eq!(
                back3.structural_key(),
                expected_key,
                "SDF V3000 lost stereo for {smi}"
            );
        }

        let dwar = write_dwar(&[(&caffeine, &caffeine.properties)]);
        let rec = parse_dwar(&dwar)?
            .into_iter()
            .next()
            .ok_or(MolError::EmptyInput)??;
        assert_eq!(rec.name, "Caffeine");

        // Corrupt V2000 >999 columns rejected cleanly
        let corrupt = "bad\n  t\n\n1200  0  0  0  0  0  0  0  0999 V2000\nM  END\n";
        assert!(matches!(
            parse_sdf(corrupt).first(),
            Some(Err(MolError::InvalidSdf { .. }))
        ));

        // Headerless CSV without SMILES header rejected cleanly
        assert!(matches!(
            parse_csv_records("ID,Value\n1,42\n"),
            Err(MolError::InvalidCsv(_))
        ));
        Ok(())
    }

    #[test]
    fn valence_warning_flags_pentavalent_carbon() -> Result<(), MolError> {
        let mut mol = MolGraph::new();
        let c = mol.add_atom(MolAtom::new(0, "C")?);
        for _ in 0..5 {
            let f = mol.add_atom(MolAtom::new(0, "F")?);
            mol.add_bond(c, f, BondOrder::Single, BondStereo::None)?;
        }
        mol.perceive()?;
        assert!(mol.atoms[c].valence_exceeded);
        assert!(mol.warnings.iter().any(|w| w.contains("Valence error")));
        Ok(())
    }

    #[test]
    fn generates_valid_3d_coords_pdb_and_all_20_amino_acid_peptides() -> Result<(), MolError> {
        // 1. Small drug molecule 3D generation + PDB + 3D SDF roundtrip
        let mut sildenafil =
            parse_smiles("CCCc1nn(C)c2c(=O)[nH]c(-c3cc(S(=O)(=O)N4CCN(C)CC4)ccc3OCC)nc12")?;
        sildenafil.name = Some("Sildenafil".to_string());
        generate_3d_coords(&mut sildenafil);
        assert!(
            sildenafil
                .atoms
                .iter()
                .all(|a| { a.coords_3d.is_some_and(|p| p.iter().all(|v| v.is_finite())) })
        );
        let z_span: f64 = sildenafil
            .atoms
            .iter()
            .filter_map(|a| a.coords_3d.map(|p| p[2].abs()))
            .fold(0.0, f64::max);
        assert!(
            z_span > 1.0,
            "expected genuinely 3D conformer, z_span={z_span}"
        );

        // Verify Ibuprofen isobutyl/chiral C(sp3) bond angles and non-planar 3D span
        let mut ibu = parse_smiles("CC(C)Cc1ccc(C(C)C(=O)O)cc1")?;
        generate_3d_coords(&mut ibu);
        let ibu_z: f64 = ibu
            .atoms
            .iter()
            .filter_map(|a| a.coords_3d.map(|p| p[2].abs()))
            .fold(0.0, f64::max);
        assert!(ibu_z > 0.85, "expected non-flat Ibuprofen 3D, z={ibu_z}");
        for (ci, a) in ibu.atoms.iter().enumerate() {
            if !a.aromatic && a.symbol == "C" && ibu.adj[ci].len() >= 3 {
                let is_carbonyl = ibu.adj[ci]
                    .iter()
                    .any(|&(_, bid)| ibu.bonds[bid].order == BondOrder::Double);
                if !is_carbonyl {
                    let pc = a.coords_3d.ok_or(MolError::EmptyInput)?;
                    let p0 = ibu.atoms[ibu.adj[ci][0].0]
                        .coords_3d
                        .ok_or(MolError::EmptyInput)?;
                    let p1 = ibu.atoms[ibu.adj[ci][1].0]
                        .coords_3d
                        .ok_or(MolError::EmptyInput)?;
                    let u = [p0[0] - pc[0], p0[1] - pc[1], p0[2] - pc[2]];
                    let v = [p1[0] - pc[0], p1[1] - pc[1], p1[2] - pc[2]];
                    let du = (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt();
                    let dv = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
                    let cos_a =
                        ((u[0] * v[0] + u[1] * v[1] + u[2] * v[2]) / (du * dv)).clamp(-1.0, 1.0);
                    let deg = cos_a.acos().to_degrees();
                    assert!(
                        (96.0..=120.0).contains(&deg),
                        "unexpected sp3 angle {deg:.1} on Ibuprofen atom {ci}"
                    );
                }
            }
        }

        for b in &sildenafil.bonds {
            let p1 = sildenafil.atoms[b.u]
                .coords_3d
                .ok_or(MolError::EmptyInput)?;
            let p2 = sildenafil.atoms[b.v]
                .coords_3d
                .ok_or(MolError::EmptyInput)?;
            let d = ((p1[0] - p2[0]).powi(2) + (p1[1] - p2[1]).powi(2) + (p1[2] - p2[2]).powi(2))
                .sqrt();
            assert!(
                (0.95..=2.25).contains(&d),
                "unrealistic 3D bond distance {d} between atoms {} and {}",
                b.u,
                b.v
            );
        }

        let pdb = write_pdb(&sildenafil);
        assert!(pdb.contains("HETATM") && pdb.contains("CONECT") && pdb.contains("END"));

        let sdf3d = write_sdf_v2000_3d(&[&sildenafil]);
        let parsed3d = parse_sdf(&sdf3d)
            .into_iter()
            .next()
            .ok_or(MolError::EmptyInput)??;
        assert_eq!(parsed3d.formula(), sildenafil.formula());
        assert_eq!(parsed3d.structural_key(), sildenafil.structural_key());
        assert!(parsed3d.atoms.iter().any(|a| a.coords_3d.is_some()));

        // Aromatic [nH] preserved through 3D SDF roundtrip (Indole)
        let mut indole = parse_smiles("c1ccc2[nH]ccc2c1")?;
        generate_3d_coords(&mut indole);
        let indole_sdf = write_sdf_v2000_3d(&[&indole]);
        let indole_back = parse_sdf(&indole_sdf)
            .into_iter()
            .next()
            .ok_or(MolError::EmptyInput)??;
        assert_eq!(indole_back.formula(), "C8H7N");
        assert_eq!(indole_back.structural_key(), indole.structural_key());

        // >999 atoms automatically switches to valid CTAB V3000 3D
        let mut big_mol = MolGraph::new();
        for i in 0..1005 {
            let mut a = MolAtom::new(i, "C")?;
            a.coords_3d = Some([i as f64 * 1.5, 0.0, 0.25]);
            big_mol.add_atom(a);
            if i > 0 {
                big_mol.add_bond(i - 1, i, BondOrder::Single, BondStereo::None)?;
            }
        }
        let big_sdf = write_sdf_v2000_3d(&[&big_mol]);
        assert!(big_sdf.contains("999 V3000") && big_sdf.contains("M  V30 COUNTS 1005 1004"));

        // 2. Peptide sequence parsing (1-letter, 3-letter, caps, errors)
        let p_3 = parse_peptide_sequence("Ac-Ala-Gly-Cys-Pro-Tyr-NH2")?;
        assert_eq!(p_3.one_letter, "AGCPY");
        assert!(p_3.n_acetyl && p_3.c_amide);

        // Build full 20-amino-acid peptide in both BetaStrand and AlphaHelix
        let all20_strand = build_peptide("ACDEFGHIKLMNPQRSTVWY", SecondaryStructure::BetaStrand)?;
        let all20_helix = build_peptide(
            "Ac-ACDEFGHIKLMNPQRSTVWY-NH2",
            SecondaryStructure::AlphaHelix,
        )?;
        assert!(
            all20_strand.warnings.is_empty(),
            "{:?}",
            all20_strand.warnings
        );
        assert!(
            all20_helix.warnings.is_empty(),
            "{:?}",
            all20_helix.warnings
        );
        assert!(all20_strand.atoms.iter().all(|a| a.coords_2d.is_some()));
        assert!(all20_strand.atoms.iter().all(|a| a.coords_3d.is_some()));

        let pep_pdb = write_pdb(&all20_helix);
        assert!(
            pep_pdb.contains("ATOM  ") && pep_pdb.contains("ALA A") && pep_pdb.contains("TRP A")
        );
        assert!(pep_pdb.contains("\nTER\n"));

        // Proline closes a 5-membered pyrrolidine ring
        let pro = build_peptide("PG", SecondaryStructure::Extended)?;
        assert!(pro.rings.iter().any(|r| r.len() == 5));

        // L-Ala CA CIP should match S
        let l_ala_pep = build_peptide("Ala", SecondaryStructure::BetaStrand)?;
        let ca = l_ala_pep
            .atoms
            .iter()
            .find(|a| a.pdb_meta.as_ref().is_some_and(|m| m.atom_name == "CA"))
            .ok_or(MolError::EmptyInput)?;
        assert_eq!(ca.cip_label, Some('S'));

        // Reject invalid residue and oversized sequence cleanly
        assert!(build_peptide("AXG", SecondaryStructure::BetaStrand).is_err());
        let too_long = "A".repeat(MAX_PEPTIDE_RESIDUES + 1);
        assert!(build_peptide(&too_long, SecondaryStructure::BetaStrand).is_err());
        Ok(())
    }

    #[test]
    fn chemdraw_parity_analysis_cdxml_rxn_superatoms_and_names() -> Result<(), MolError> {
        // 1. Exact mass, DBE, charge-annotated formula, elemental analysis, ESI adducts
        let aspirin = parse_smiles("CC(=O)Oc1ccccc1C(=O)O")?;
        let an = analyze_molecule(&aspirin);
        assert_eq!(an.formula_with_charge, "C9H8O4");
        assert!(
            (an.exact_mass - 180.0423).abs() < 0.001,
            "exact_mass={}",
            an.exact_mass
        );
        assert!((an.dbe - 6.0).abs() < 0.01);
        assert!(
            an.adducts
                .iter()
                .any(|a| a.label == "[M+H]+" && (a.mz - 181.0495).abs() < 0.002)
        );
        assert!(an.elemental_analysis_string.contains("C, 60.0"));

        let acetate = parse_smiles("CC(=O)[O-]")?;
        assert_eq!(analyze_molecule(&acetate).formula_with_charge, "C2H3O2-");

        // 2. CDXML roundtrip preserving chirality and aromatic [nH]
        let ibu_r = parse_smiles("CC(C)Cc1ccc([C@@H](C)C(=O)O)cc1")?;
        let cdxml = write_cdxml(&ibu_r);
        assert!(cdxml.contains("<CDXML") && cdxml.contains("Wedge"));
        let ibu_back = parse_cdxml(&cdxml)?;
        assert_eq!(ibu_back.formula(), ibu_r.formula());
        assert_eq!(ibu_back.structural_key(), ibu_r.structural_key());

        let indole = parse_smiles("c1ccc2[nH]ccc2c1")?;
        let indole_back = parse_cdxml(&write_cdxml(&indole))?;
        assert_eq!(indole_back.formula(), "C8H7N");
        assert_eq!(indole_back.structural_key(), indole.structural_key());

        // 3. MDL .rxn V2000 + stoichiometry / atom economy / yield
        let rxn = parse_reaction_smiles("CC(=O)O.CCO>>CC(=O)OCC.O")?;
        let rxn_v2000 = write_rxn_v2000(&rxn);
        assert!(rxn_v2000.starts_with("$RXN"));
        let rxn_back = parse_rxn_v2000(&rxn_v2000)?;
        assert_eq!(rxn_back.reactants.len(), 2);
        assert_eq!(rxn_back.products.len(), 2);
        let stoich = compute_stoichiometry(&rxn_back, 1.0, Some(70.48));
        assert!((stoich.atom_economy_pct - 83.0).abs() < 1.0);
        assert!((stoich.theoretical_yield_mg - 88.11).abs() < 0.5);
        assert!(
            stoich
                .isolated_yield_pct
                .is_some_and(|y| (y - 80.0).abs() < 1.0)
        );

        // 4. Superatoms & protecting groups
        let aniline = parse_smiles("Nc1ccccc1")?;
        let boc_aniline = attach_superatom_group(&aniline, Some(0), "Boc")?;
        assert_eq!(boc_aniline.formula(), "C11H15NO2");

        // 5. Name-to-Structure (dictionary + IUPAC systematic parser)
        assert_eq!(resolve_chemical_name("Aspirin")?.formula(), "C9H8O4");
        assert_eq!(
            resolve_chemical_name("4-chlorophenol")?.formula(),
            "C6H5ClO"
        );
        assert_eq!(
            resolve_chemical_name("2-methylpropan-1-ol")?.formula(),
            "C4H10O"
        );
        assert_eq!(resolve_chemical_name("cyclopentanone")?.formula(), "C5H8O");
        assert_eq!(resolve_chemical_name("hexanoic acid")?.formula(), "C6H12O2");
        Ok(())
    }

    #[test]
    fn baricitinib_3d_sulfonamide_oxygens_are_tetrahedral_not_collinear() -> Result<(), MolError> {
        let bari = parse_smiles("CCS(=O)(=O)N1CC(C1)(CC#N)n2cc(cn2)c3ncnc4c3cc[nH]4")?;
        let coords = compute_3d_coords(&bari);
        let s_idx = bari
            .atoms
            .iter()
            .position(|a| a.symbol == "S")
            .ok_or(MolError::EmptyInput)?;
        let oxy: Vec<usize> = bari.adj[s_idx]
            .iter()
            .map(|&(nb, _)| nb)
            .filter(|&nb| bari.atoms[nb].symbol == "O")
            .collect();
        assert_eq!(oxy.len(), 2, "Baricitinib sulfonamide must have 2 oxygens");
        let ps = coords[s_idx];
        let p1 = coords[oxy[0]];
        let p2 = coords[oxy[1]];
        let oo = ((p1[0] - p2[0]).powi(2) + (p1[1] - p2[1]).powi(2) + (p1[2] - p2[2]).powi(2))
            .sqrt();
        assert!(
            oo > 2.20,
            "Sulfonamide O=S=O oxygens collapsed (d={oo:.3} Å)"
        );
        let v1 = [p1[0] - ps[0], p1[1] - ps[1], p1[2] - ps[2]];
        let v2 = [p2[0] - ps[0], p2[1] - ps[1], p2[2] - ps[2]];
        let l1 = (v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]).sqrt();
        let l2 = (v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]).sqrt();
        let cos_a = ((v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]) / (l1 * l2)).clamp(-1.0, 1.0);
        let ang_deg = cos_a.acos().to_degrees();
        assert!(
            (100.0..=135.0).contains(&ang_deg),
            "Sulfonamide O-S-O angle out of tetrahedral range: {ang_deg:.1}°"
        );
        Ok(())
    }
}
