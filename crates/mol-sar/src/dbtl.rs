//! Design–Build–Test–Learn (DBTL) Experimental Assay Ingestion Engine.
//!
//! Closes the laboratory feedback loop: ingests experimental biochemical/cellular assay tables
//! (`IC50_nM`, `EC50_nM`, `Kd_nM`, `Ki_nM`, `pIC50`, optional `CLint`/solubility), converts potency
//! cleanly to `pIC50`, enriches with local AdmetCraft profiles, recomputes Matched Molecular Pairs
//! (MMPs) and Activity Cliffs (SALI), and emits:
//! - updated SAR table CSV (`pIC50`, `LE`, `LLE`, ADMET),
//! - QSAR calibration dataset (`smiles,y`) ready for `admetcraft-cli train-qsar`, and
//! - top high-LLE seed SMILES for Round 2 library enumeration / lead optimization.

use std::fmt::Write as _;

use mol_admet_bridge::AdmetCraftProvider;
use mol_graph::parse_smiles;
use serde::{Deserialize, Serialize};

use crate::mmp::{
    ActivityCliff, DuplicatePair, MatchedMolecularPair, find_activity_cliffs,
    find_duplicate_compounds, find_matched_molecular_pairs,
};
use crate::table::SarTable;

/// Summary metrics for one ingested experimental compound.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DbtlCompoundSummary {
    pub name: String,
    pub canonical_smiles: String,
    pub pic50: Option<f64>,
    pub le: Option<f64>,
    pub lle: Option<f64>,
    pub mw: f64,
    pub clogp: Option<f64>,
    pub tpsa: Option<f64>,
    pub lipinski_violations: Option<usize>,
    pub pains_count: Option<usize>,
    pub series: String,
}

/// Complete DBTL assay ingestion report.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DbtlIngestReport {
    pub total_rows_read: usize,
    pub valid_compounds: usize,
    pub skipped_invalid_smiles: usize,
    pub compounds_with_potency: usize,
    pub duplicates_found: usize,
    pub mmp_pairs_found: usize,
    pub activity_cliffs_found: usize,
    pub top_candidates: Vec<DbtlCompoundSummary>,
    pub top_mmp_gains: Vec<MatchedMolecularPair>,
    pub activity_cliffs: Vec<ActivityCliff>,
    pub duplicates: Vec<DuplicatePair>,
}

/// Detect delimiter (`\t` vs `,`) and split a single CSV/TSV line honoring simple double quotes.
fn split_delimited_line(line: &str, delim: char) -> Vec<String> {
    let mut fields = Vec::new();
    let mut cur = String::new();
    let mut in_quotes = false;
    for ch in line.chars() {
        if ch == '"' {
            in_quotes = !in_quotes;
        } else if ch == delim && !in_quotes {
            fields.push(cur.trim().to_string());
            cur.clear();
        } else {
            cur.push(ch);
        }
    }
    fields.push(cur.trim().to_string());
    fields
}

fn parse_pic50_from_assay_fields(headers: &[String], cols: &[String]) -> Option<f64> {
    for (idx, h) in headers.iter().enumerate() {
        let Some(raw) = cols.get(idx) else { continue };
        let clean = raw.trim().trim_start_matches(['>', '<', '~', '=']).trim();
        let Ok(val) = clean.parse::<f64>() else {
            continue;
        };
        if !val.is_finite() {
            continue;
        }
        let hl = h.to_ascii_lowercase();
        if hl == "pic50" || hl == "pec50" || hl == "pkd" || hl == "pki" {
            return Some((val * 1000.0).round() / 1000.0);
        }
        if (hl.contains("ic50") || hl.contains("ec50") || hl.contains("kd") || hl.contains("ki"))
            && val > 0.0
        {
            let p = if hl.contains("um") || hl.contains("µm") {
                6.0 - val.log10()
            } else {
                // Default biochemical assay unit: nM
                9.0 - val.log10()
            };
            return Some((p * 1000.0).round() / 1000.0);
        }
    }
    None
}

/// Ingest an experimental CSV/TSV assay file into a `SarTable`, enrich with AdmetCraft,
/// compute MMPs + Activity Cliffs, and return both the enriched table and DBTL report.
pub fn ingest_experimental_assay_text(
    assay_text: &str,
    base_table: Option<SarTable>,
    min_cliff_tanimoto: f64,
) -> Result<(SarTable, DbtlIngestReport), String> {
    let mut lines = assay_text.lines().filter(|l| !l.trim().is_empty());
    let header_line = lines
        .next()
        .ok_or_else(|| "empty experimental assay file".to_string())?;
    let delim = if header_line.contains('\t') {
        '\t'
    } else {
        ','
    };
    let headers: Vec<String> = split_delimited_line(header_line, delim)
        .into_iter()
        .map(|h| h.trim_start_matches('\u{feff}').to_string())
        .collect();

    let smiles_col = headers
        .iter()
        .position(|h| {
            let hl = h.to_ascii_lowercase();
            hl == "smiles" || hl == "canonical_smiles" || hl == "smi"
        })
        .ok_or_else(|| {
            format!("assay table missing required 'smiles' column (found: {headers:?})")
        })?;
    let name_col = headers.iter().position(|h| {
        let hl = h.to_ascii_lowercase();
        matches!(
            hl.as_str(),
            "compound_id" | "name" | "id" | "mol_id" | "title"
        )
    });
    let series_col = headers
        .iter()
        .position(|h| h.eq_ignore_ascii_case("series"));

    let mut table = base_table.unwrap_or_default();
    let mut total_rows_read = 0usize;
    let mut skipped_invalid_smiles = 0usize;

    for raw_line in lines {
        total_rows_read += 1;
        let cols = split_delimited_line(raw_line, delim);
        let Some(smi) = cols.get(smiles_col).map(|s| s.trim()) else {
            skipped_invalid_smiles += 1;
            continue;
        };
        if smi.is_empty() {
            skipped_invalid_smiles += 1;
            continue;
        }
        let Ok(mut mol) = parse_smiles(smi) else {
            skipped_invalid_smiles += 1;
            continue;
        };
        let name = name_col
            .and_then(|i| cols.get(i))
            .filter(|s| !s.trim().is_empty())
            .cloned()
            .unwrap_or_else(|| format!("ASSAY_{total_rows_read}"));
        mol.name = Some(name);
        for (idx, h) in headers.iter().enumerate() {
            if let Some(v) = cols.get(idx)
                && !v.is_empty()
            {
                mol.properties.insert(h.clone(), v.clone());
            }
        }
        let pic50 = parse_pic50_from_assay_fields(&headers, &cols);
        let series = series_col
            .and_then(|i| cols.get(i))
            .filter(|s| !s.trim().is_empty())
            .cloned();
        table.add_molecule(mol, pic50, series);
    }

    table.enrich_with_admet(&AdmetCraftProvider::default());
    let duplicates = find_duplicate_compounds(&table.rows);
    let mut mmps = find_matched_molecular_pairs(&table.rows);
    mmps.sort_by(|a, b| {
        b.delta_pic50
            .unwrap_or(0.0)
            .abs()
            .total_cmp(&a.delta_pic50.unwrap_or(0.0).abs())
    });
    let cliffs = find_activity_cliffs(&table.rows, min_cliff_tanimoto);

    let mut summaries: Vec<DbtlCompoundSummary> = table
        .rows
        .iter()
        .map(|r| DbtlCompoundSummary {
            name: r.name.clone(),
            canonical_smiles: r.canonical_smiles.clone(),
            pic50: r.pic50,
            le: r.ligand_efficiency().map(|v| (v * 1000.0).round() / 1000.0),
            lle: r
                .lipophilic_ligand_efficiency()
                .map(|v| (v * 1000.0).round() / 1000.0),
            mw: (r.mw * 100.0).round() / 100.0,
            clogp: r
                .admet
                .as_ref()
                .map(|a| (a.consensus_logp * 100.0).round() / 100.0),
            tpsa: r.admet.as_ref().map(|a| (a.tpsa * 10.0).round() / 10.0),
            lipinski_violations: r.admet.as_ref().map(|a| a.lipinski_violations),
            pains_count: r.admet.as_ref().map(|a| a.pains_count),
            series: r.series.clone(),
        })
        .collect();

    summaries.sort_by(|a, b| {
        b.pic50
            .unwrap_or(f64::NEG_INFINITY)
            .total_cmp(&a.pic50.unwrap_or(f64::NEG_INFINITY))
            .then_with(|| {
                b.lle
                    .unwrap_or(f64::NEG_INFINITY)
                    .total_cmp(&a.lle.unwrap_or(f64::NEG_INFINITY))
            })
    });

    let compounds_with_potency = table.rows.iter().filter(|r| r.pic50.is_some()).count();
    let report = DbtlIngestReport {
        total_rows_read,
        valid_compounds: table.rows.len(),
        skipped_invalid_smiles,
        compounds_with_potency,
        duplicates_found: duplicates.len(),
        mmp_pairs_found: mmps.len(),
        activity_cliffs_found: cliffs.len(),
        top_candidates: summaries.into_iter().take(20).collect(),
        top_mmp_gains: mmps.into_iter().take(20).collect(),
        activity_cliffs: cliffs,
        duplicates,
    };

    Ok((table, report))
}

/// Format the updated SAR table with measured `pIC50`, `LE`, `LLE`, and ADMET columns.
#[must_use]
pub fn render_dbtl_sar_csv(table: &SarTable) -> String {
    let mut out = String::from(
        "name,canonical_smiles,pic50,le,lle,mw,clogp,tpsa,lipinski_violations,pains_count,series\n",
    );
    for r in &table.rows {
        let p = r.pic50.map_or(String::new(), |v| format!("{v:.3}"));
        let le = r
            .ligand_efficiency()
            .map_or(String::new(), |v| format!("{v:.3}"));
        let lle = r
            .lipophilic_ligand_efficiency()
            .map_or(String::new(), |v| format!("{v:.3}"));
        let clogp = r
            .admet
            .as_ref()
            .map_or(String::new(), |a| format!("{:.2}", a.consensus_logp));
        let tpsa = r
            .admet
            .as_ref()
            .map_or(String::new(), |a| format!("{:.1}", a.tpsa));
        let lip = r
            .admet
            .as_ref()
            .map_or(String::new(), |a| a.lipinski_violations.to_string());
        let pains = r
            .admet
            .as_ref()
            .map_or(String::new(), |a| a.pains_count.to_string());
        let _ = writeln!(
            out,
            "\"{}\",\"{}\",{p},{le},{lle},{:.2},{clogp},{tpsa},{lip},{pains},\"{}\"",
            r.name.replace('"', "'"),
            r.canonical_smiles.replace('"', "'"),
            r.mw,
            r.series.replace('"', "'")
        );
    }
    out
}

/// Format a clean `smiles,y` CSV dataset for local `admetcraft-cli train-qsar` recalibration.
#[must_use]
pub fn render_qsar_recalibration_csv(table: &SarTable) -> String {
    let mut out = String::from("smiles,y,compound_id\n");
    for r in &table.rows {
        if let Some(p) = r.pic50 {
            let _ = writeln!(
                out,
                "\"{}\",{p:.4},\"{}\"",
                r.canonical_smiles.replace('"', "'"),
                r.name.replace('"', "'")
            );
        }
    }
    out
}

/// Format top measured leads as `.smi` seeds for Round 2 optimization (`enumerate` / `lead-opt`).
#[must_use]
pub fn render_next_round_seeds_smi(report: &DbtlIngestReport, top_n: usize) -> String {
    let mut out = String::new();
    for c in report.top_candidates.iter().take(top_n.max(1)) {
        let _ = writeln!(out, "{}\t{}", c.canonical_smiles, c.name);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dbtl_assay_ingest_converts_ic50_nm_computes_mmp_and_exports_qsar_csv() -> Result<(), String>
    {
        let csv = "\
compound_id,smiles,ic50_nm,series
CMPD_1,Cc1ccccc1,1000,CoreA
CMPD_2,Oc1ccccc1,10,CoreA
CMPD_3,Fc1ccccc1,250,CoreA
";
        let (table, report) = ingest_experimental_assay_text(csv, None, 0.35)?;
        assert_eq!(report.valid_compounds, 3);
        assert_eq!(report.compounds_with_potency, 3);
        assert!(report.mmp_pairs_found >= 1);
        assert_eq!(report.top_candidates[0].name, "CMPD_2");
        assert!((report.top_candidates[0].pic50.unwrap_or(0.0) - 8.0).abs() < 1e-3);

        let qsar_csv = render_qsar_recalibration_csv(&table);
        assert!(qsar_csv.starts_with("smiles,y,compound_id\n"));
        assert!(qsar_csv.contains("CMPD_2"));
        let seeds = render_next_round_seeds_smi(&report, 2);
        assert!(seeds.contains("CMPD_2"));
        Ok(())
    }
}
