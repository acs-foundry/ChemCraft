//! Layer 3 (`mol-engine`): Session state, Copy-on-Write Undo/Redo, `CommandRegistry`,
//! atomic `.chemcraft` ZIP+JSON project persistence, and hardened `PathSandbox`.

#![forbid(unsafe_code)]

pub mod chimerax;
pub mod commands;
pub mod dockcraft;
pub mod library;
pub mod project;
pub mod sandbox;

pub use chimerax::{
    ChimeraXBundle, ChimeraXError, default_chimerax_export_dir, is_chimerax_available,
    launch_chimerax_for_mol, write_chimerax_bundle,
};
pub use commands::{
    ChemSession, Command, CommandError, CommandInfo, CommandOutcome, CommandRegistry,
};
pub use dockcraft::{
    DockCraftBundle, DockCraftError, DockCraftResultRow, DockCraftResults, DockCraftRun,
    MergeOutcome, default_dockcraft_handoff_dir, default_dockcraft_results_dir,
    is_dockcraft_available, is_dockcraft_cli_available, latest_dockcraft_results,
    launch_dockcraft_for_mols, merge_dockcraft_results, read_dockcraft_results,
    write_dockcraft_bundle,
};
pub use library::{
    BuiltinLibraryPreset, DEFAULT_MAX_LIBRARY_RECORDS, LibraryError, LibraryFilterConfig,
    LibraryScreeningResult, MAX_LIBRARY_UNCOMPRESSED_BYTES, ScreenedLibraryHit, ScreeningFunnel,
    default_library_dir, discover_local_libraries, export_library_hits_sdf3d,
    filter_and_prepare_library, load_builtin_library, materialize_builtin_library,
    maxmin_diversity_pick, protonate_at_ph74, stream_load_library, stream_load_library_from_file,
};
pub use project::{
    ChemCraftDocument, ProjectError, decode_chemcraft_zip, encode_chemcraft_zip,
    save_chemcraft_atomic,
};
pub use sandbox::{MAX_READ_BYTES, PathSandbox, SandboxError};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn undo_redo_and_chemcraft_zip_roundtrip_work() -> Result<(), Box<dyn std::error::Error>> {
        let mut session = ChemSession::new();
        CommandRegistry::execute(&mut session, Command::LoadDemo50Sar)?;
        assert_eq!(session.sar_table.rows.len(), 50);

        // Mutate carboxylic acid -> tetrazole on Ibuprofen and check ADMET delta
        CommandRegistry::execute(
            &mut session,
            Command::MutateBioisostere {
                atom_index: None,
                group: "tetrazole".to_string(),
            },
        )?;
        assert!(session.last_mutation_delta.is_some());
        assert!(session.can_undo());

        CommandRegistry::execute(&mut session, Command::Undo)?;
        assert!(session.can_redo());

        let doc = ChemCraftDocument {
            format: "chemcraft".to_string(),
            version: "0.1.0".to_string(),
            canvas_mol: (*session.canvas_mol).clone(),
            reaction_arrows: (*session.reaction_arrows).clone(),
            sar_table: (*session.sar_table).clone(),
            sar_filter: session.sar_filter.clone(),
        };
        let bytes = encode_chemcraft_zip(&doc)?;
        let decoded = decode_chemcraft_zip(&bytes)?;
        assert_eq!(decoded.sar_table.rows.len(), 50);
        assert_eq!(
            decoded.canvas_mol.canonical_smiles(),
            session.canvas_mol.canonical_smiles()
        );
        Ok(())
    }

    #[test]
    fn sandbox_blocks_symlinks_and_path_traversal() -> Result<(), Box<dyn std::error::Error>> {
        let base = std::env::temp_dir().join("chemcraft_sandbox_test");
        let _ = std::fs::remove_dir_all(&base);
        let r = base.join("read");
        let w = base.join("write");
        let outside = base.join("outside.txt");
        std::fs::create_dir_all(&r)?;
        std::fs::create_dir_all(&w)?;
        std::fs::write(&outside, "secret")?;
        #[cfg(unix)]
        std::os::unix::fs::symlink(&outside, r.join("link.smi"))?;

        let sb = PathSandbox::new(Some(&r), Some(&w))?;
        assert!(
            sb.read_to_string_limited(std::path::Path::new("../outside.txt"), 1024)
                .is_err()
        );
        #[cfg(unix)]
        assert!(
            sb.read_to_string_limited(std::path::Path::new("link.smi"), 1024)
                .is_err()
        );
        assert!(sb.open_write(std::path::Path::new("../pwned.svg")).is_err());

        let _ = std::fs::remove_dir_all(&base);
        Ok(())
    }

    #[test]
    fn zip_decompression_bomb_and_4000_openat2_dir_races_are_blocked()
    -> Result<(), Box<dyn std::error::Error>> {
        use flate2::{Compression, write::DeflateEncoder};
        use std::io::Write;

        // 1. Decompression bomb test: ~40 MB zeros compressed to tiny archive declaring uncomp=42 bytes
        let huge = vec![b' '; 40 * 1024 * 1024];
        let mut enc = DeflateEncoder::new(Vec::new(), Compression::fast());
        enc.write_all(&huge)?;
        let comp = enc.finish()?;
        let mut zip = Vec::new();
        zip.extend_from_slice(&0x0403_4b50u32.to_le_bytes());
        zip.extend_from_slice(&20u16.to_le_bytes());
        zip.extend_from_slice(&0u16.to_le_bytes());
        zip.extend_from_slice(&8u16.to_le_bytes());
        zip.extend_from_slice(&0u16.to_le_bytes());
        zip.extend_from_slice(&0u16.to_le_bytes());
        zip.extend_from_slice(&0u32.to_le_bytes());
        zip.extend_from_slice(&(comp.len() as u32).to_le_bytes());
        zip.extend_from_slice(&42u32.to_le_bytes()); // spoofed 42-byte header!
        zip.extend_from_slice(&13u16.to_le_bytes());
        zip.extend_from_slice(&0u16.to_le_bytes());
        zip.extend_from_slice(b"manifest.json");
        zip.extend_from_slice(&comp);
        assert!(matches!(
            decode_chemcraft_zip(&zip),
            Err(ProjectError::InvalidArchive(_))
        ));

        // 2. 4,000 intermediate-directory symlink swap race test
        #[cfg(unix)]
        {
            use std::os::unix::fs::symlink;
            use std::sync::{
                Arc,
                atomic::{AtomicBool, Ordering},
            };
            let base = std::env::temp_dir().join("chemcraft_sb_race_4k");
            let _ = std::fs::remove_dir_all(&base);
            let r = base.join("read");
            let w = base.join("write");
            let outside = base.join("outside");
            std::fs::create_dir_all(r.join("sub"))?;
            std::fs::create_dir_all(w.join("sub"))?;
            std::fs::create_dir_all(&outside)?;
            std::fs::write(outside.join("leak.smi"), "SECRET_OUTSIDE")?;
            std::fs::write(r.join("sub/leak.smi"), "SAFE_INSIDE")?;

            let sb = PathSandbox::new(Some(&r), Some(&w))?;
            let stop = Arc::new(AtomicBool::new(false));
            let stop_r = Arc::clone(&stop);
            let r_link = r.join("sub");
            let w_link = w.join("sub");

            let racer = std::thread::spawn(move || {
                let mut t = false;
                while !stop_r.load(Ordering::Relaxed) {
                    t = !t;
                    let _ = std::fs::remove_file(&r_link);
                    let _ = std::fs::remove_dir_all(&r_link);
                    let _ = std::fs::remove_file(&w_link);
                    let _ = std::fs::remove_dir_all(&w_link);
                    if t {
                        let _ = symlink(&outside, &r_link);
                        let _ = symlink(&outside, &w_link);
                    } else {
                        let _ = std::fs::create_dir(&r_link);
                        let _ = std::fs::write(r_link.join("leak.smi"), "SAFE_INSIDE");
                        let _ = std::fs::create_dir(&w_link);
                    }
                }
            });

            for _ in 0..4_000 {
                if let Ok((_, txt)) =
                    sb.read_to_string_limited(std::path::Path::new("sub/leak.smi"), 1024)
                {
                    assert_ne!(txt, "SECRET_OUTSIDE", "leaked outside read-root!");
                }
                if let Ok((_, mut f)) = sb.open_write(std::path::Path::new("sub/pwn.svg")) {
                    let _ = f.write_all(b"x");
                }
            }
            stop.store(true, Ordering::Relaxed);
            let _ = racer.join();
            assert!(
                !base.join("outside/pwn.svg").exists(),
                "wrote outside write-root!"
            );
            let _ = std::fs::remove_dir_all(&base);
        }
        Ok(())
    }

    #[test]
    fn chemical_library_streaming_gzip_ph74_and_funnel_screening_work()
    -> Result<(), Box<dyn std::error::Error>> {
        let base = std::env::temp_dir().join("chemcraft_lib_screen_test");
        let _ = std::fs::remove_dir_all(&base);
        std::fs::create_dir_all(&base)?;

        // 1. Built-in presets parse cleanly with zero empty entries
        let fda = load_builtin_library(BuiltinLibraryPreset::FdaApprovedMedicines);
        let frags = load_builtin_library(BuiltinLibraryPreset::FragmentRo3);
        let kin = load_builtin_library(BuiltinLibraryPreset::KinaseAndGpcrScaffolds);
        assert!(fda.len() >= 30);
        assert!(frags.len() >= 25);
        assert!(kin.len() >= 20);

        // 2. Write a gzipped SMILES library and stream-load it back
        let smi_path =
            materialize_builtin_library(BuiltinLibraryPreset::FdaApprovedMedicines, &base)?;
        let smi_gz = base.join("fda.smi.gz");
        {
            use std::io::Write;
            let raw = std::fs::read(&smi_path)?;
            let mut enc = flate2::write::GzEncoder::new(
                std::fs::File::create(&smi_gz)?,
                flate2::Compression::fast(),
            );
            enc.write_all(&raw)?;
            enc.finish()?;
        }
        let discovered = discover_local_libraries(&base);
        assert!(discovered.len() >= 2);
        let streamed = stream_load_library(&smi_gz, 100)?;
        assert_eq!(streamed.len(), fda.len());

        // 3. Parallel screening funnel + pH 7.4 ionization + MaxMin diversity pick + 3D SDF export
        let cfg = LibraryFilterConfig {
            funnel: ScreeningFunnel::RuleOf5,
            diversity_pick_n: Some(12),
            protonate_ph74: true,
            generate_3d: true,
            ..LibraryFilterConfig::default()
        };
        let res = filter_and_prepare_library("FDA.smi.gz", streamed, &cfg)?;
        assert_eq!(res.hits.len(), 12);
        assert!(
            res.hits
                .iter()
                .all(|h| h.mol.atoms.iter().all(|a| a.coords_3d.is_some()))
        );

        // Verify Ibuprofen carboxylate ionized at pH 7.4 when unfiltered
        let mut ibu = mol_graph::parse_smiles("CC(C)Cc1ccc(C(C)C(=O)O)cc1")?;
        protonate_at_ph74(&mut ibu);
        assert!(
            ibu.atoms.iter().any(|a| a.formal_charge == -1),
            "expected deprotonated carboxylate on Ibuprofen at pH 7.4"
        );

        let sdf_out = base.join("screen_ready_3d.sdf");
        let n_written = export_library_hits_sdf3d(&res, &sdf_out)?;
        assert_eq!(n_written, 12);
        assert!(sdf_out.exists());

        let _ = std::fs::remove_dir_all(&base);
        Ok(())
    }

    #[test]
    fn chemdraw_editing_commands_charge_cycle_superatom_flip_rotate_and_name_work()
    -> Result<(), Box<dyn std::error::Error>> {
        let mut session = ChemSession::new();
        CommandRegistry::execute(
            &mut session,
            Command::LoadChemicalName {
                name: "aniline".to_string(),
            },
        )?;
        assert_eq!(session.canvas_mol.formula(), "C6H7N");

        // Attach Boc protecting group to N0, then test Undo
        CommandRegistry::execute(
            &mut session,
            Command::AttachSuperatomGroup {
                atom_idx: Some(0),
                abbr: "Boc".to_string(),
            },
        )?;
        assert_eq!(session.canvas_mol.formula(), "C11H15NO2");
        CommandRegistry::execute(&mut session, Command::Undo)?;
        assert_eq!(session.canvas_mol.formula(), "C6H7N");

        // Adjust formal charge (+1 on N0 -> anilinium C6H8N+)
        CommandRegistry::execute(
            &mut session,
            Command::AdjustAtomCharge {
                atom_idx: 0,
                delta: 1,
            },
        )?;
        assert_eq!(session.canvas_mol.atoms[0].formal_charge, 1);

        // Cycle bond order Single -> Double -> Triple -> Single on Ethane
        CommandRegistry::execute(
            &mut session,
            Command::LoadCanvasSmiles {
                smiles: "CC".to_string(),
                name: None,
            },
        )?;
        CommandRegistry::execute(&mut session, Command::CycleBondOrder { bond_idx: 0 })?;
        assert_eq!(session.canvas_mol.formula(), "C2H4");
        CommandRegistry::execute(&mut session, Command::CycleBondOrder { bond_idx: 0 })?;
        assert_eq!(session.canvas_mol.formula(), "C2H2");

        // ConnectAtoms closes propane into cyclopropane
        CommandRegistry::execute(
            &mut session,
            Command::LoadCanvasSmiles {
                smiles: "CCC".to_string(),
                name: None,
            },
        )?;
        CommandRegistry::execute(
            &mut session,
            Command::ConnectAtoms {
                u: 0,
                v: 2,
                order: mol_graph::BondOrder::Single,
                stereo: mol_graph::BondStereo::None,
            },
        )?;
        assert_eq!(session.canvas_mol.formula(), "C3H6");
        assert_eq!(session.canvas_mol.rings.len(), 1);

        // Flip H preserves stereochemical parity & Rotate 30° stays finite
        CommandRegistry::execute(&mut session, Command::FlipCanvas2D { horizontal: true })?;
        CommandRegistry::execute(&mut session, Command::RotateCanvas2D { degrees: 30.0 })?;
        assert!(
            session
                .canvas_mol
                .atoms
                .iter()
                .all(|a| a.coords_2d.is_some())
        );
        Ok(())
    }
}
