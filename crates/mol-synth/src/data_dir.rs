#![forbid(unsafe_code)]

use std::path::PathBuf;

/// Resolve the data directory from explicit values, so the policy can be tested without
/// touching the process environment.
///
/// Order: `CRAFT_DATA_DIR` · `$XDG_DATA_HOME/CraftBio` · `$HOME/.local/share/CraftBio` ·
/// a relative `CraftBio` as a last resort. Blank values are treated as unset.
#[must_use]
fn resolve_craft_data_dir(
    craft_data_dir: Option<&str>,
    xdg_data_home: Option<&str>,
    home: Option<&str>,
) -> PathBuf {
    fn non_blank(v: Option<&str>) -> Option<&str> {
        v.map(str::trim).filter(|s| !s.is_empty())
    }

    if let Some(dir) = non_blank(craft_data_dir) {
        return PathBuf::from(dir);
    }
    if let Some(xdg) = non_blank(xdg_data_home) {
        return PathBuf::from(xdg).join("CraftBio");
    }
    if let Some(home) = non_blank(home) {
        return PathBuf::from(home).join(".local/share/CraftBio");
    }
    PathBuf::from("CraftBio")
}

/// Returns the centralized Craft bio-data directory.
///
/// `CRAFT_DATA_DIR` if set, otherwise `$XDG_DATA_HOME/CraftBio`, otherwise
/// `$HOME/.local/share/CraftBio`. Set `CRAFT_DATA_DIR` to keep multi-gigabyte building-block
/// catalogs off the root partition.
#[must_use]
pub fn craft_data_dir() -> PathBuf {
    resolve_craft_data_dir(
        std::env::var("CRAFT_DATA_DIR").ok().as_deref(),
        std::env::var("XDG_DATA_HOME").ok().as_deref(),
        std::env::var("HOME").ok().as_deref(),
    )
}

/// Permanent building-block catalog directory (`CRAFT_DATA_DIR/building_blocks`).
#[must_use]
pub fn building_blocks_dir() -> PathBuf {
    craft_data_dir().join("building_blocks")
}

/// Default local commercial building-block TSV catalog (`CRAFT_DATA_DIR/building_blocks/medchem_stock.tsv`).
#[must_use]
pub fn default_stock_tsv_path() -> PathBuf {
    building_blocks_dir().join("medchem_stock.tsv")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_resolution_order_is_craft_dir_then_xdg_then_home() {
        assert_eq!(
            resolve_craft_data_dir(Some("  /opt/craft  "), Some("/x"), Some("/home/u")),
            PathBuf::from("/opt/craft"),
            "CRAFT_DATA_DIR wins and is trimmed"
        );
        assert_eq!(
            resolve_craft_data_dir(None, Some("/x"), Some("/home/u")),
            PathBuf::from("/x/CraftBio")
        );
        assert_eq!(
            resolve_craft_data_dir(Some("   "), Some("/x"), Some("/home/u")),
            PathBuf::from("/x/CraftBio"),
            "a blank value counts as unset"
        );
        assert_eq!(
            resolve_craft_data_dir(None, None, Some("/home/u")),
            PathBuf::from("/home/u/.local/share/CraftBio")
        );
        assert_eq!(
            resolve_craft_data_dir(None, None, None),
            PathBuf::from("CraftBio"),
            "last resort is relative, never an absolute guess"
        );
    }

    /// The default must derive from the caller's environment and nothing else. This is a
    /// public repository: a path hardcoded for one machine both leaks a username and makes
    /// the default behaviour wrong everywhere else.
    #[test]
    fn the_default_data_dir_derives_only_from_the_environment() {
        assert_eq!(
            resolve_craft_data_dir(None, None, Some("/home/someone-else")),
            PathBuf::from("/home/someone-else/.local/share/CraftBio"),
            "a foreign environment must produce a path inside that environment"
        );
        let p = craft_data_dir();
        assert!(
            !p.to_string_lossy().contains("/Datos/CraftBio"),
            "default data dir still resolves to a hardcoded developer path: {}",
            p.display()
        );
    }
}
