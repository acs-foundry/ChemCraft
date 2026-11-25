//! UCSF ChimeraX interop: self-contained 3D bundle export (`.pdb` / 3D `.sdf` + `.cxc`)
//! and Flatpak-sandbox-aware launcher (`chimerax`).

use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use mol_graph::{MolGraph, write_pdb, write_sdf_v2000_3d};
use thiserror::Error;

#[derive(Debug, Error)]
pub enum ChimeraXError {
    #[error("canvas molecule is empty")]
    EmptyMolecule,
    #[error("ChimeraX bundle I/O failed: {0}")]
    Io(String),
    #[error("could not launch `chimerax` (bundle exported at {bundle_dir}): {msg}")]
    Launch { bundle_dir: String, msg: String },
}

#[derive(Debug, Clone)]
pub struct ChimeraXBundle {
    pub dir: PathBuf,
    pub script_path: PathBuf,
    pub structure_path: PathBuf,
    pub is_peptide: bool,
}

/// Choose a ChimeraX export directory readable inside ChimeraX's Flatpak sandbox
/// (`$XDG_DATA_HOME/ChimeraX/chemcraft_export` or `$HOME/Projects/research/data/outputs/chemcraft_export`).
pub fn default_chimerax_export_dir() -> PathBuf {
    if let Ok(xdg_data) = std::env::var("XDG_DATA_HOME")
        && !xdg_data.trim().is_empty()
    {
        return PathBuf::from(xdg_data)
            .join("ChimeraX")
            .join("chemcraft_export");
    }
    if let Ok(home) = std::env::var("HOME")
        && !home.trim().is_empty()
    {
        let home_p = PathBuf::from(home);
        let xdg_default = home_p
            .join(".local")
            .join("share")
            .join("ChimeraX")
            .join("chemcraft_export");
        return xdg_default;
    }
    PathBuf::from("chemcraft_chimerax")
}

/// Check if a `chimerax` executable exists on `$PATH`.
pub fn is_chimerax_available() -> bool {
    if let Ok(path_var) = std::env::var("PATH") {
        for dir in std::env::split_paths(&path_var) {
            if dir.join("chimerax").is_file() {
                return true;
            }
        }
    }
    false
}

/// Sanitize a base filename prefix safely (ASCII alphanumeric + `_`/`-`).
fn sanitize_stem(raw: &str) -> String {
    let mut s: String = raw
        .chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || c == '_' || c == '-' {
                c
            } else {
                '_'
            }
        })
        .take(40)
        .collect();
    let trimmed = s.trim_matches('_').to_string();
    if trimmed.is_empty() {
        s = "chemcraft_structure".to_string();
        s
    } else {
        trimmed
    }
}

/// Write a self-contained ChimeraX visualization bundle (`<stem>.cxc` + `<stem>.pdb` / `.sdf`) into `dir`.
///
/// The generated `.cxc` script uses relative sibling filenames so ChimeraX can open it cleanly
/// without host-path leakage and inside sandboxed environments.
pub fn write_chimerax_bundle(
    mol: &MolGraph,
    dir: &Path,
    stem_override: Option<&str>,
) -> Result<ChimeraXBundle, ChimeraXError> {
    if mol.atoms.is_empty() {
        return Err(ChimeraXError::EmptyMolecule);
    }
    std::fs::create_dir_all(dir)
        .map_err(|e| ChimeraXError::Io(format!("create_dir_all {}: {e}", dir.display())))?;

    let is_peptide = mol
        .atoms
        .iter()
        .any(|a| a.pdb_meta.as_ref().is_some_and(|m| !m.is_hetatm));
    let stem = sanitize_stem(
        stem_override
            .or(mol.name.as_deref())
            .unwrap_or(if is_peptide {
                "chemcraft_biomolecule"
            } else {
                "chemcraft_molecule"
            }),
    );

    let pdb_name = format!("{stem}.pdb");
    let sdf_name = format!("{stem}.sdf");
    let cxc_name = format!("{stem}.cxc");
    let pdb_path = dir.join(&pdb_name);
    let sdf_path = dir.join(&sdf_name);
    let script_path = dir.join(&cxc_name);

    let pdb_text = write_pdb(mol);
    std::fs::write(&pdb_path, pdb_text.as_bytes())
        .map_err(|e| ChimeraXError::Io(format!("write {}: {e}", pdb_path.display())))?;

    let sdf_text = write_sdf_v2000_3d(&[mol]);
    std::fs::write(&sdf_path, sdf_text.as_bytes())
        .map_err(|e| ChimeraXError::Io(format!("write {}: {e}", sdf_path.display())))?;

    let mut cxc = String::new();
    let _ = writeln!(cxc, "# ChemCraft -> UCSF ChimeraX visualization script");
    if is_peptide {
        let _ = writeln!(cxc, "open {pdb_name}");
        let _ = writeln!(cxc, "preset \"overall look\" publication");
        let _ = writeln!(cxc, "cartoon");
        let _ = writeln!(cxc, "show atoms");
        let _ = writeln!(cxc, "style stick");
        let _ = writeln!(cxc, "color byhetero");
        let _ = writeln!(cxc, "hbonds reveal true");
    } else {
        // 3D SDF V2000 preserves exact double/triple/aromatic bond orders natively in ChimeraX
        let _ = writeln!(cxc, "open {sdf_name}");
        let _ = writeln!(cxc, "preset \"overall look\" publication");
        let _ = writeln!(cxc, "style ball");
        let _ = writeln!(cxc, "color byhetero");
    }
    let _ = writeln!(cxc, "lighting soft");
    let _ = writeln!(cxc, "set bgColor white");
    let _ = writeln!(cxc, "view");

    std::fs::write(&script_path, cxc.as_bytes())
        .map_err(|e| ChimeraXError::Io(format!("write {}: {e}", script_path.display())))?;

    let structure_path = if is_peptide { pdb_path } else { sdf_path };
    Ok(ChimeraXBundle {
        dir: dir.to_path_buf(),
        script_path,
        structure_path,
        is_peptide,
    })
}

/// Export a bundle inside `dir` (or default Flatpak-whitelisted ChimeraX export dir) and launch `chimerax`.
pub fn launch_chimerax_for_mol(
    mol: &MolGraph,
    dir: Option<&Path>,
) -> Result<ChimeraXBundle, ChimeraXError> {
    let target_dir = dir
        .map(Path::to_path_buf)
        .unwrap_or_else(default_chimerax_export_dir);
    let bundle = write_chimerax_bundle(mol, &target_dir, None)?;
    Command::new("chimerax")
        .arg(&bundle.script_path)
        .current_dir(&bundle.dir)
        .stdin(Stdio::null())
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .spawn()
        .map_err(|e| ChimeraXError::Launch {
            bundle_dir: bundle.dir.display().to_string(),
            msg: e.to_string(),
        })?;
    Ok(bundle)
}
