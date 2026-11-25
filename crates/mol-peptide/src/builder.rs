#![forbid(unsafe_code)]

//! Atomic `MolGraph` polymerization, cyclization, and cross-link builder for linear & cyclic
//! peptides, depsipeptides, lipopeptides, stapled helices, and Non-Ribosomal Peptides (NRPs).

use std::collections::BTreeSet;

use mol_graph::{BondOrder, BondStereo, ChiralTag, MolAtom, MolGraph, PdbAtomMeta, parse_smiles};

use crate::geometry::{ResidueAtoms, layout_peptide_2d, layout_peptide_3d};
use crate::monomers::{CCapKind, MonomerCatalog, SidechainAnchorKind};
use crate::notation::{CrossLink, PeptideError, PeptideSpec};

/// Built peptide/NRP structure paired with its normalized spec, HELM string, and canonical SMILES.
#[derive(Debug, Clone)]
pub struct BuiltPeptide {
    pub spec: PeptideSpec,
    pub mol: MolGraph,
    pub notation: String,
    pub helm: String,
    pub canonical_smiles: String,
}

#[derive(Debug, Clone)]
struct StampedResidue {
    n_id: usize,
    ca_id: usize,
    c_id: usize,
    oxt_id: usize,
    o_id: Option<usize>,
    sc_anchor_id: Option<usize>,
    sc_leave_id: Option<usize>,
    all_atoms: Vec<usize>,
    anchor_kind: SidechainAnchorKind,
    is_d: bool,
}

/// Merge `frag` into `dst`, returning the index offset applied to every atom of `frag`.
fn append_fragment(
    dst: &mut MolGraph,
    frag: &MolGraph,
    res_name: &str,
    res_seq: u16,
    is_hetatm: bool,
) -> Result<usize, PeptideError> {
    let offset = dst.atoms.len();
    for (local_i, a) in frag.atoms.iter().enumerate() {
        let mut copy = a.clone();
        copy.id = offset + local_i;
        copy.chiral_seq = a
            .chiral_seq
            .iter()
            .map(|&nbr| {
                if nbr == usize::MAX {
                    usize::MAX
                } else {
                    nbr + offset
                }
            })
            .collect();
        copy.pdb_meta = Some(PdbAtomMeta {
            atom_name: format!("{}{}", a.symbol, local_i + 1),
            res_name: res_name.to_string(),
            res_seq,
            chain: 'A',
            is_hetatm,
        });
        dst.add_atom(copy);
    }
    for b in &frag.bonds {
        let eid = dst.add_bond(b.u + offset, b.v + offset, b.order, b.stereo)?;
        dst.bonds[eid].stereo_from_u = b.stereo_from_u;
        dst.bonds[eid].cis_trans = b.cis_trans;
    }
    Ok(offset)
}

fn find_marker(mol: &MolGraph, start: usize, end: usize, iso: u16) -> Option<usize> {
    (start..end).find(|&i| mol.atoms[i].isotope == Some(iso))
}

fn strip_markers(mol: &mut MolGraph, start: usize, end: usize) {
    for a in &mut mol.atoms[start..end] {
        if let Some(iso) = a.isotope
            && (101..=109).contains(&iso)
        {
            a.isotope = None;
            // Release bracket H-lock on non-aromatic, neutral marker atoms so condensation updates H cleanly
            if !a.aromatic && a.formal_charge == 0 {
                a.from_bracket = false;
                a.explicit_h = 0;
            }
        }
    }
}

fn delete_atoms_remap(mol: &mut MolGraph, to_delete: &BTreeSet<usize>) -> Vec<Option<usize>> {
    let old_atoms = std::mem::take(&mut mol.atoms);
    let old_bonds = std::mem::take(&mut mol.bonds);
    mol.adj.clear();
    let mut remap = vec![None; old_atoms.len()];
    for (old_id, mut atom) in old_atoms.into_iter().enumerate() {
        if to_delete.contains(&old_id) {
            continue;
        }
        let new_id = mol.atoms.len();
        atom.id = new_id;
        remap[old_id] = Some(new_id);
        mol.atoms.push(atom);
        mol.adj.push(Vec::new());
    }
    for atom in &mut mol.atoms {
        if !atom.chiral_seq.is_empty() {
            atom.chiral_seq = atom
                .chiral_seq
                .iter()
                .map(|&nbr| {
                    if nbr == usize::MAX {
                        usize::MAX
                    } else {
                        remap.get(nbr).copied().flatten().unwrap_or(usize::MAX)
                    }
                })
                .collect();
        }
    }
    for b in old_bonds {
        if let (Some(nu), Some(nv)) = (remap[b.u], remap[b.v])
            && let Ok(eid) = mol.add_bond(nu, nv, b.order, b.stereo)
        {
            mol.bonds[eid].stereo_from_u = b.stereo_from_u;
            mol.bonds[eid].cis_trans = b.cis_trans;
        }
    }
    remap
}

/// Assign exact PDB IUPAC sidechain heavy-atom names for all 20 canonical amino acids
/// so exported PDBs pass structural validation (BiomolCraft `info`, MolProbity, ChimeraX)
/// without false truncated-sidechain alerts.
#[allow(clippy::too_many_arguments)]
fn assign_iupac_sidechain_names(
    mol: &mut MolGraph,
    pdb_code: &str,
    n_id: usize,
    ca_id: usize,
    c_id: usize,
    o_id: Option<usize>,
    oxt_id: usize,
) {
    let set_name = |m: &mut MolGraph, aid: usize, name: &str| {
        if let Some(meta) = m.atoms.get_mut(aid).and_then(|a| a.pdb_meta.as_mut()) {
            meta.atom_name = name.to_string();
        }
    };
    let is_backbone = |aid: usize| {
        aid == n_id || aid == ca_id || aid == c_id || aid == oxt_id || o_id == Some(aid)
    };
    let sc_nbrs = |m: &MolGraph, u: usize, exclude: &[usize]| -> Vec<usize> {
        m.adj[u]
            .iter()
            .map(|&(nb, _)| nb)
            .filter(|&nb| !is_backbone(nb) && !exclude.contains(&nb))
            .collect()
    };

    let cb_list = sc_nbrs(mol, ca_id, &[]);
    let Some(&cb) = cb_list.first() else {
        return;
    };

    match pdb_code {
        "ALA" => set_name(mol, cb, "CB"),
        "SER" => {
            set_name(mol, cb, "CB");
            if let Some(&og) = sc_nbrs(mol, cb, &[]).first() {
                set_name(mol, og, "OG");
            }
        }
        "CYS" => {
            set_name(mol, cb, "CB");
            if let Some(&sg) = sc_nbrs(mol, cb, &[]).first() {
                set_name(mol, sg, "SG");
            }
        }
        "THR" => {
            set_name(mol, cb, "CB");
            for nb in sc_nbrs(mol, cb, &[]) {
                if mol.atoms[nb].symbol == "O" {
                    set_name(mol, nb, "OG1");
                } else {
                    set_name(mol, nb, "CG2");
                }
            }
        }
        "VAL" => {
            set_name(mol, cb, "CB");
            let gs = sc_nbrs(mol, cb, &[]);
            if let Some(&g1) = gs.first() {
                set_name(mol, g1, "CG1");
            }
            if let Some(&g2) = gs.get(1) {
                set_name(mol, g2, "CG2");
            }
        }
        "LEU" => {
            set_name(mol, cb, "CB");
            if let Some(&cg) = sc_nbrs(mol, cb, &[]).first() {
                set_name(mol, cg, "CG");
                let ds = sc_nbrs(mol, cg, &[cb]);
                if let Some(&d1) = ds.first() {
                    set_name(mol, d1, "CD1");
                }
                if let Some(&d2) = ds.get(1) {
                    set_name(mol, d2, "CD2");
                }
            }
        }
        "ILE" => {
            set_name(mol, cb, "CB");
            for nb in sc_nbrs(mol, cb, &[]) {
                let next = sc_nbrs(mol, nb, &[cb]);
                if let Some(&cd1) = next.first() {
                    set_name(mol, nb, "CG1");
                    set_name(mol, cd1, "CD1");
                } else {
                    set_name(mol, nb, "CG2");
                }
            }
        }
        "MET" => {
            set_name(mol, cb, "CB");
            if let Some(&cg) = sc_nbrs(mol, cb, &[]).first() {
                set_name(mol, cg, "CG");
                if let Some(&sd) = sc_nbrs(mol, cg, &[cb]).first() {
                    set_name(mol, sd, "SD");
                    if let Some(&ce) = sc_nbrs(mol, sd, &[cg]).first() {
                        set_name(mol, ce, "CE");
                    }
                }
            }
        }
        "PRO" => {
            set_name(mol, cb, "CB");
            if let Some(&cg) = sc_nbrs(mol, cb, &[]).first() {
                set_name(mol, cg, "CG");
                if let Some(&cd) = sc_nbrs(mol, cg, &[cb]).first() {
                    set_name(mol, cd, "CD");
                }
            }
        }
        "ASP" => {
            set_name(mol, cb, "CB");
            if let Some(&cg) = sc_nbrs(mol, cb, &[]).first() {
                set_name(mol, cg, "CG");
                let mut os = sc_nbrs(mol, cg, &[cb]);
                os.sort_by_key(|&o| {
                    let is_dbl = mol.adj[cg]
                        .iter()
                        .any(|&(nb, bid)| nb == o && mol.bonds[bid].order == BondOrder::Double);
                    if is_dbl { 0u8 } else { 1u8 }
                });
                if let Some(&o1) = os.first() {
                    set_name(mol, o1, "OD1");
                }
                if let Some(&o2) = os.get(1) {
                    set_name(mol, o2, "OD2");
                }
            }
        }
        "ASN" => {
            set_name(mol, cb, "CB");
            if let Some(&cg) = sc_nbrs(mol, cb, &[]).first() {
                set_name(mol, cg, "CG");
                for nb in sc_nbrs(mol, cg, &[cb]) {
                    if mol.atoms[nb].symbol == "O" {
                        set_name(mol, nb, "OD1");
                    } else if mol.atoms[nb].symbol == "N" {
                        set_name(mol, nb, "ND2");
                    }
                }
            }
        }
        "GLU" => {
            set_name(mol, cb, "CB");
            if let Some(&cg) = sc_nbrs(mol, cb, &[]).first() {
                set_name(mol, cg, "CG");
                if let Some(&cd) = sc_nbrs(mol, cg, &[cb]).first() {
                    set_name(mol, cd, "CD");
                    let mut os = sc_nbrs(mol, cd, &[cg]);
                    os.sort_by_key(|&o| {
                        let is_dbl = mol.adj[cd]
                            .iter()
                            .any(|&(nb, bid)| nb == o && mol.bonds[bid].order == BondOrder::Double);
                        if is_dbl { 0u8 } else { 1u8 }
                    });
                    if let Some(&o1) = os.first() {
                        set_name(mol, o1, "OE1");
                    }
                    if let Some(&o2) = os.get(1) {
                        set_name(mol, o2, "OE2");
                    }
                }
            }
        }
        "GLN" => {
            set_name(mol, cb, "CB");
            if let Some(&cg) = sc_nbrs(mol, cb, &[]).first() {
                set_name(mol, cg, "CG");
                if let Some(&cd) = sc_nbrs(mol, cg, &[cb]).first() {
                    set_name(mol, cd, "CD");
                    for nb in sc_nbrs(mol, cd, &[cg]) {
                        if mol.atoms[nb].symbol == "O" {
                            set_name(mol, nb, "OE1");
                        } else if mol.atoms[nb].symbol == "N" {
                            set_name(mol, nb, "NE2");
                        }
                    }
                }
            }
        }
        "LYS" => {
            set_name(mol, cb, "CB");
            if let Some(&cg) = sc_nbrs(mol, cb, &[]).first() {
                set_name(mol, cg, "CG");
                if let Some(&cd) = sc_nbrs(mol, cg, &[cb]).first() {
                    set_name(mol, cd, "CD");
                    if let Some(&ce) = sc_nbrs(mol, cd, &[cg]).first() {
                        set_name(mol, ce, "CE");
                        if let Some(&nz) = sc_nbrs(mol, ce, &[cd]).first() {
                            set_name(mol, nz, "NZ");
                        }
                    }
                }
            }
        }
        "ARG" => {
            set_name(mol, cb, "CB");
            if let Some(&cg) = sc_nbrs(mol, cb, &[]).first() {
                set_name(mol, cg, "CG");
                if let Some(&cd) = sc_nbrs(mol, cg, &[cb]).first() {
                    set_name(mol, cd, "CD");
                    if let Some(&ne) = sc_nbrs(mol, cd, &[cg]).first() {
                        set_name(mol, ne, "NE");
                        if let Some(&cz) = sc_nbrs(mol, ne, &[cd]).first() {
                            set_name(mol, cz, "CZ");
                            let nhs = sc_nbrs(mol, cz, &[ne]);
                            if let Some(&nh1) = nhs.first() {
                                set_name(mol, nh1, "NH1");
                            }
                            if let Some(&nh2) = nhs.get(1) {
                                set_name(mol, nh2, "NH2");
                            }
                        }
                    }
                }
            }
        }
        "PHE" | "TYR" => {
            set_name(mol, cb, "CB");
            if let Some(&cg) = sc_nbrs(mol, cb, &[]).first() {
                set_name(mol, cg, "CG");
                let cds = sc_nbrs(mol, cg, &[cb]);
                if let (Some(&cd1), Some(&cd2)) = (cds.first(), cds.get(1)) {
                    set_name(mol, cd1, "CD1");
                    set_name(mol, cd2, "CD2");
                    if let Some(&ce1) = sc_nbrs(mol, cd1, &[cg]).first() {
                        set_name(mol, ce1, "CE1");
                        if let Some(&ce2) = sc_nbrs(mol, cd2, &[cg]).first() {
                            set_name(mol, ce2, "CE2");
                            if let Some(&cz) = sc_nbrs(mol, ce1, &[cd1]).first() {
                                set_name(mol, cz, "CZ");
                                if pdb_code == "TYR"
                                    && let Some(&oh) = sc_nbrs(mol, cz, &[ce1, ce2]).first()
                                {
                                    set_name(mol, oh, "OH");
                                }
                            }
                        }
                    }
                }
            }
        }
        "HIS" => {
            set_name(mol, cb, "CB");
            if let Some(&cg) = sc_nbrs(mol, cb, &[]).first() {
                set_name(mol, cg, "CG");
                let mut nd1 = None;
                let mut cd2 = None;
                for nb in sc_nbrs(mol, cg, &[cb]) {
                    if mol.atoms[nb].symbol == "N" {
                        nd1 = Some(nb);
                        set_name(mol, nb, "ND1");
                    } else {
                        cd2 = Some(nb);
                        set_name(mol, nb, "CD2");
                    }
                }
                if let Some(n1) = nd1
                    && let Some(&ce1) = sc_nbrs(mol, n1, &[cg]).first()
                {
                    set_name(mol, ce1, "CE1");
                    if let Some(&ne2) = sc_nbrs(mol, ce1, &[n1]).first() {
                        set_name(mol, ne2, "NE2");
                    }
                } else if let Some(c2) = cd2
                    && let Some(&ne2) = sc_nbrs(mol, c2, &[cg]).first()
                {
                    set_name(mol, ne2, "NE2");
                }
            }
        }
        "TRP" => {
            set_name(mol, cb, "CB");
            if let Some(&cg) = sc_nbrs(mol, cb, &[]).first() {
                set_name(mol, cg, "CG");
                let mut cd1 = None;
                let mut cd2 = None;
                for nb in sc_nbrs(mol, cg, &[cb]) {
                    if mol.adj[nb].len() == 2 {
                        cd1 = Some(nb);
                        set_name(mol, nb, "CD1");
                    } else {
                        cd2 = Some(nb);
                        set_name(mol, nb, "CD2");
                    }
                }
                if let (Some(d1), Some(d2)) = (cd1, cd2)
                    && let Some(&ne1) = sc_nbrs(mol, d1, &[cg]).first()
                {
                    set_name(mol, ne1, "NE1");
                    if let Some(&ce2) = sc_nbrs(mol, ne1, &[d1]).first() {
                        set_name(mol, ce2, "CE2");
                        if let Some(&ce3) = sc_nbrs(mol, d2, &[cg, ce2]).first() {
                            set_name(mol, ce3, "CE3");
                            if let Some(&cz2) = sc_nbrs(mol, ce2, &[ne1, d2]).first() {
                                set_name(mol, cz2, "CZ2");
                                if let Some(&cz3) = sc_nbrs(mol, ce3, &[d2]).first() {
                                    set_name(mol, cz3, "CZ3");
                                    if let Some(&ch2) = sc_nbrs(mol, cz2, &[ce2]).first() {
                                        set_name(mol, ch2, "CH2");
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        _ => {}
    }
}

/// Replace old atom `old_nbr` with `new_nbr` inside `target_atom.chiral_seq` (keeps exact R/S after condensation).
fn patch_chiral_seq(mol: &mut MolGraph, target_atom: usize, old_nbr: usize, new_nbr: usize) {
    if let Some(a) = mol.atoms.get_mut(target_atom) {
        for slot in &mut a.chiral_seq {
            if *slot == old_nbr {
                *slot = new_nbr;
            }
        }
    }
}

/// Build a chemically exact 2D + 3D `MolGraph` from any parsed `PeptideSpec`.
pub fn build_peptide_from_spec(
    spec: &PeptideSpec,
    cat: &MonomerCatalog,
) -> Result<BuiltPeptide, PeptideError> {
    if spec.residues.is_empty() {
        return Err(PeptideError::EmptyInput);
    }

    let mut mol = MolGraph::new();
    let mut stamped: Vec<StampedResidue> = Vec::with_capacity(spec.residues.len());
    let mut to_delete: BTreeSet<usize> = BTreeSet::new();

    // 1. Optional N-terminal cap or fatty-acyl starter (skipped for head-to-tail macrocycles)
    let mut n_cap_c: Option<usize> = None;
    if !spec.is_head_to_tail() && !spec.n_cap.is_empty() && !spec.n_cap.eq_ignore_ascii_case("H") {
        let cap_spec = cat
            .find_n_cap(&spec.n_cap)
            .ok_or_else(|| PeptideError::UnknownMonomer(spec.n_cap.clone()))?;
        if !cap_spec.smiles_template.is_empty() {
            let frag = parse_smiles(&cap_spec.smiles_template)?;
            let start = append_fragment(&mut mol, &frag, &cap_spec.pdb_code, 0, true)?;
            let end = mol.atoms.len();
            let c_id = find_marker(&mol, start, end, 103).ok_or_else(|| {
                PeptideError::Notation(format!("N-cap '{}' missing [103C]", cap_spec.code))
            })?;
            let oxt_id = find_marker(&mol, start, end, 104).ok_or_else(|| {
                PeptideError::Notation(format!("N-cap '{}' missing [104O]", cap_spec.code))
            })?;
            strip_markers(&mut mol, start, end);
            to_delete.insert(oxt_id);
            n_cap_c = Some(c_id);
        }
    }

    // 2. Stamp every residue monomer with exact L/D stereochemistry and optional N-methylation
    for (idx, ru) in spec.residues.iter().enumerate() {
        let mspec = cat
            .find_monomer(&ru.base_code)
            .ok_or_else(|| PeptideError::UnknownMonomer(ru.base_code.clone()))?;
        let frag = parse_smiles(&mspec.smiles_template)?;
        let res_seq = (idx + 1) as u16;
        let is_hetatm =
            mspec.class != crate::monomers::MonomerClass::Canonical || ru.is_d || ru.is_n_methyl;
        let start = append_fragment(&mut mol, &frag, &mspec.pdb_code, res_seq, is_hetatm)?;
        let end = mol.atoms.len();

        let n_id = find_marker(&mol, start, end, 101).ok_or_else(|| {
            PeptideError::Notation(format!("monomer '{}' missing [101N]/[101O]", mspec.code))
        })?;
        let ca_id = find_marker(&mol, start, end, 102).ok_or_else(|| {
            PeptideError::Notation(format!("monomer '{}' missing [102C]", mspec.code))
        })?;
        let c_id = find_marker(&mol, start, end, 103).ok_or_else(|| {
            PeptideError::Notation(format!("monomer '{}' missing [103C]", mspec.code))
        })?;
        let oxt_id = find_marker(&mol, start, end, 104).ok_or_else(|| {
            PeptideError::Notation(format!("monomer '{}' missing [104O]", mspec.code))
        })?;
        let sc_anchor_id = find_marker(&mol, start, end, 105);
        let sc_leave_id = find_marker(&mol, start, end, 106);

        // Locate carbonyl =O bonded to c_id
        let o_id = mol.adj[c_id]
            .iter()
            .find(|&&(nb, bid)| {
                nb != oxt_id
                    && mol.atoms[nb].symbol == "O"
                    && mol.bonds[bid].order == BondOrder::Double
            })
            .map(|&(nb, _)| nb);

        strip_markers(&mut mol, start, end);

        // Assign standard PDB backbone atom names (N, CA, C, O, SG)
        if let Some(m) = mol.atoms[n_id].pdb_meta.as_mut() {
            m.atom_name = if mspec.is_hydroxy_acid {
                "O1".to_string()
            } else {
                "N".to_string()
            };
        }
        if let Some(m) = mol.atoms[ca_id].pdb_meta.as_mut() {
            m.atom_name = "CA".to_string();
        }
        if let Some(m) = mol.atoms[c_id].pdb_meta.as_mut() {
            m.atom_name = "C".to_string();
        }
        if let Some(oid) = o_id
            && let Some(m) = mol.atoms[oid].pdb_meta.as_mut()
        {
            m.atom_name = "O".to_string();
        }
        if let Some(sid) = sc_anchor_id
            && mspec.anchor_kind == SidechainAnchorKind::Thiol
            && let Some(m) = mol.atoms[sid].pdb_meta.as_mut()
        {
            m.atom_name = "SG".to_string();
        }
        assign_iupac_sidechain_names(&mut mol, &mspec.pdb_code, n_id, ca_id, c_id, o_id, oxt_id);

        // Apply D-isomer inversion on Cα
        if ru.is_d {
            mol.atoms[ca_id].chiral = match mol.atoms[ca_id].chiral {
                ChiralTag::Clockwise => ChiralTag::CounterClockwise,
                ChiralTag::CounterClockwise => ChiralTag::Clockwise,
                other => other,
            };
        }

        let mut all_atoms: Vec<usize> = (start..end).collect();

        // Apply composable N-methylation on backbone N
        if ru.is_n_methyl {
            if mspec.is_tertiary_n {
                return Err(PeptideError::Topology(format!(
                    "cannot N-methylate tertiary/hydroxy residue '{}'",
                    mspec.code
                )));
            }
            let mut me = MolAtom::new(0, "C")?;
            me.pdb_meta = Some(PdbAtomMeta {
                atom_name: "CN".to_string(),
                res_name: mspec.pdb_code.clone(),
                res_seq,
                chain: 'A',
                is_hetatm: true,
            });
            let me_id = mol.add_atom(me);
            mol.add_bond(n_id, me_id, BondOrder::Single, BondStereo::None)?;
            all_atoms.push(me_id);
        }

        stamped.push(StampedResidue {
            n_id,
            ca_id,
            c_id,
            oxt_id,
            o_id,
            sc_anchor_id,
            sc_leave_id,
            all_atoms,
            anchor_kind: mspec.anchor_kind,
            is_d: ru.is_d,
        });
    }

    // 3. Connect N-cap to residue 1
    if let Some(cap_c) = n_cap_c {
        mol.add_bond(cap_c, stamped[0].n_id, BondOrder::Single, BondStereo::None)?;
    }

    // Check for δ-amide linkages replacing normal α-carboxyl polymerization (e.g. ACV tripeptide)
    let delta_pairs: BTreeSet<(usize, usize)> = spec
        .cross_links
        .iter()
        .filter_map(|cl| match *cl {
            CrossLink::DeltaAmide { acid_res, next_res } => Some((acid_res - 1, next_res - 1)),
            _ => None,
        })
        .collect();

    // 4. Polymerize backbone residue i -> i+1
    for i in 0..stamped.len().saturating_sub(1) {
        let next_n = stamped[i + 1].n_id;
        if delta_pairs.contains(&(i, i + 1)) {
            let sc_c = stamped[i].sc_anchor_id.ok_or_else(|| {
                PeptideError::Topology(format!(
                    "residue {} lacks a sidechain carboxyl for delta-amide",
                    i + 1
                ))
            })?;
            let sc_o = stamped[i].sc_leave_id.ok_or_else(|| {
                PeptideError::Topology(format!(
                    "residue {} lacks a sidechain leaving OH for delta-amide",
                    i + 1
                ))
            })?;
            to_delete.insert(sc_o);
            mol.add_bond(sc_c, next_n, BondOrder::Single, BondStereo::None)?;
        } else {
            to_delete.insert(stamped[i].oxt_id);
            mol.add_bond(stamped[i].c_id, next_n, BondOrder::Single, BondStereo::None)?;
        }
    }

    // 5. Apply covalent cyclizations & cross-links
    let last_idx = stamped.len() - 1;
    let mut c_term_consumed = false;
    for cl in &spec.cross_links {
        match *cl {
            CrossLink::HeadToTail => {
                to_delete.insert(stamped[last_idx].oxt_id);
                mol.add_bond(
                    stamped[last_idx].c_id,
                    stamped[0].n_id,
                    BondOrder::Single,
                    BondStereo::None,
                )?;
                c_term_consumed = true;
            }
            CrossLink::Macrolactone {
                from_res,
                to_oh_res,
            } => {
                let d_i = from_res - 1;
                let a_i = to_oh_res - 1;
                if stamped[a_i].anchor_kind != SidechainAnchorKind::Hydroxyl {
                    return Err(PeptideError::Topology(format!(
                        "macrolactone acceptor residue {to_oh_res} ({}) must carry a sidechain -OH (Ser/Thr/Tyr/Hpg/Bmt)",
                        spec.residues[a_i].base_code
                    )));
                }
                let oh_id = stamped[a_i].sc_anchor_id.ok_or_else(|| {
                    PeptideError::Topology(format!("residue {to_oh_res} missing hydroxyl anchor"))
                })?;
                if d_i == last_idx {
                    to_delete.insert(stamped[d_i].oxt_id);
                    mol.add_bond(
                        stamped[d_i].c_id,
                        oh_id,
                        BondOrder::Single,
                        BondStereo::None,
                    )?;
                    c_term_consumed = true;
                } else {
                    let sc_c = stamped[d_i].sc_anchor_id.ok_or_else(|| {
                        PeptideError::Topology(format!(
                            "non-terminal lactone donor residue {from_res} must have sidechain -COOH"
                        ))
                    })?;
                    let sc_o = stamped[d_i].sc_leave_id.ok_or_else(|| {
                        PeptideError::Topology(format!("residue {from_res} missing sidechain OH"))
                    })?;
                    to_delete.insert(sc_o);
                    mol.add_bond(sc_c, oh_id, BondOrder::Single, BondStereo::None)?;
                }
            }
            CrossLink::Lactam {
                from_acid_res,
                to_amine_res,
                from_c_terminus,
            } => {
                let d_i = from_acid_res - 1;
                let a_i = to_amine_res - 1;
                if stamped[a_i].anchor_kind != SidechainAnchorKind::Amine {
                    return Err(PeptideError::Topology(format!(
                        "lactam acceptor residue {to_amine_res} ({}) must carry a sidechain -NH2 (Lys/Orn/Dab/Dap)",
                        spec.residues[a_i].base_code
                    )));
                }
                let nh_id = stamped[a_i].sc_anchor_id.ok_or_else(|| {
                    PeptideError::Topology(format!("residue {to_amine_res} missing amine anchor"))
                })?;
                if from_c_terminus
                    || (d_i == last_idx
                        && stamped[d_i].anchor_kind != SidechainAnchorKind::Carboxyl)
                {
                    to_delete.insert(stamped[d_i].oxt_id);
                    mol.add_bond(
                        stamped[d_i].c_id,
                        nh_id,
                        BondOrder::Single,
                        BondStereo::None,
                    )?;
                    if d_i == last_idx {
                        c_term_consumed = true;
                    }
                } else {
                    if stamped[d_i].anchor_kind != SidechainAnchorKind::Carboxyl {
                        return Err(PeptideError::Topology(format!(
                            "sidechain lactam donor residue {from_acid_res} ({}) must carry -COOH (Asp/Glu/Aad/mGlu)",
                            spec.residues[d_i].base_code
                        )));
                    }
                    let sc_c = stamped[d_i].sc_anchor_id.ok_or_else(|| {
                        PeptideError::Topology(format!("residue {from_acid_res} missing carboxyl"))
                    })?;
                    let sc_o = stamped[d_i].sc_leave_id.ok_or_else(|| {
                        PeptideError::Topology(format!("residue {from_acid_res} missing -OH"))
                    })?;
                    to_delete.insert(sc_o);
                    mol.add_bond(sc_c, nh_id, BondOrder::Single, BondStereo::None)?;
                }
            }
            CrossLink::Disulfide { res_a, res_b } => {
                let ia = res_a - 1;
                let ib = res_b - 1;
                if stamped[ia].anchor_kind != SidechainAnchorKind::Thiol
                    || stamped[ib].anchor_kind != SidechainAnchorKind::Thiol
                {
                    return Err(PeptideError::Topology(format!(
                        "disulfide ({res_a}-{res_b}) requires thiol residues (Cys) at both positions"
                    )));
                }
                let sa = stamped[ia]
                    .sc_anchor_id
                    .ok_or_else(|| PeptideError::Topology("missing SG on res_a".to_string()))?;
                let sb = stamped[ib]
                    .sc_anchor_id
                    .ok_or_else(|| PeptideError::Topology("missing SG on res_b".to_string()))?;
                mol.add_bond(sa, sb, BondOrder::Single, BondStereo::None)?;
            }
            CrossLink::RcmStaple { res_a, res_b } => {
                let ia = res_a - 1;
                let ib = res_b - 1;
                if stamped[ia].anchor_kind != SidechainAnchorKind::AlkeneStaple
                    || stamped[ib].anchor_kind != SidechainAnchorKind::AlkeneStaple
                {
                    return Err(PeptideError::Topology(format!(
                        "RCM staple ({res_a}-{res_b}) requires alkene stapling monomers (S5/R8)"
                    )));
                }
                let ca = stamped[ia].sc_anchor_id.ok_or_else(|| {
                    PeptideError::Topology("missing alkene C on res_a".to_string())
                })?;
                let la = stamped[ia].sc_leave_id.ok_or_else(|| {
                    PeptideError::Topology("missing terminal =CH2 on res_a".to_string())
                })?;
                let cb = stamped[ib].sc_anchor_id.ok_or_else(|| {
                    PeptideError::Topology("missing alkene C on res_b".to_string())
                })?;
                let lb = stamped[ib].sc_leave_id.ok_or_else(|| {
                    PeptideError::Topology("missing terminal =CH2 on res_b".to_string())
                })?;
                to_delete.insert(la);
                to_delete.insert(lb);
                mol.add_bond(ca, cb, BondOrder::Double, BondStereo::None)?;
            }
            CrossLink::Lanthionine {
                thiol_res,
                dehydro_res,
            } => {
                let it = thiol_res - 1;
                let id = dehydro_res - 1;
                if stamped[it].anchor_kind != SidechainAnchorKind::Thiol
                    || stamped[id].anchor_kind != SidechainAnchorKind::Dehydro
                {
                    return Err(PeptideError::Topology(format!(
                        "lanthionine ({thiol_res}->{dehydro_res}) requires Cys -> Dha/Dhb"
                    )));
                }
                let s_id = stamped[it]
                    .sc_anchor_id
                    .ok_or_else(|| PeptideError::Topology("missing SG".to_string()))?;
                let cb_id = stamped[id]
                    .sc_anchor_id
                    .ok_or_else(|| PeptideError::Topology("missing Dha/Dhb Cβ".to_string()))?;
                // Reduce Cα=Cβ double bond to single bond upon Michael addition
                let ca_d = stamped[id].ca_id;
                if let Some(&(_, bid)) = mol.adj[ca_d].iter().find(|&&(nb, _)| nb == cb_id) {
                    mol.bonds[bid].order = BondOrder::Single;
                }
                mol.add_bond(s_id, cb_id, BondOrder::Single, BondStereo::None)?;
            }
            CrossLink::TriazoleClick {
                alkyne_res,
                azide_res,
            } => {
                let ip = alkyne_res - 1;
                let ia = azide_res - 1;
                if stamped[ip].anchor_kind != SidechainAnchorKind::AlkyneClick
                    || stamped[ia].anchor_kind != SidechainAnchorKind::AzideClick
                {
                    return Err(PeptideError::Topology(format!(
                        "click bridge ({alkyne_res}-{azide_res}) requires Pra (alkyne) and Aha (azide)"
                    )));
                }
                let c_int = stamped[ip]
                    .sc_anchor_id
                    .ok_or_else(|| PeptideError::Topology("missing Pra C".to_string()))?;
                let c_term = stamped[ip]
                    .sc_leave_id
                    .ok_or_else(|| PeptideError::Topology("missing Pra terminal C".to_string()))?;
                let n1 = stamped[ia]
                    .sc_anchor_id
                    .ok_or_else(|| PeptideError::Topology("missing Aha N1".to_string()))?;
                // Locate N2 and N3 of azide chain
                let n2 = mol.adj[n1]
                    .iter()
                    .find(|&&(nb, _)| mol.atoms[nb].symbol == "N")
                    .map(|&(nb, _)| nb)
                    .ok_or_else(|| PeptideError::Topology("malformed azide N2".to_string()))?;
                let n3 = mol.adj[n2]
                    .iter()
                    .find(|&&(nb, _)| nb != n1 && mol.atoms[nb].symbol == "N")
                    .map(|&(nb, _)| nb)
                    .ok_or_else(|| PeptideError::Topology("malformed azide N3".to_string()))?;
                // Neutralize azide charges & form aromatic 1,4-disubstituted 1,2,3-triazole
                mol.atoms[n1].formal_charge = 0;
                mol.atoms[n1].from_bracket = false;
                mol.atoms[n2].formal_charge = 0;
                mol.atoms[n2].from_bracket = false;
                mol.atoms[n3].formal_charge = 0;
                mol.atoms[n3].from_bracket = false;
                if let Some(&(_, bid)) = mol.adj[c_int].iter().find(|&&(nb, _)| nb == c_term) {
                    mol.bonds[bid].order = BondOrder::Double;
                }
                if let Some(&(_, bid)) = mol.adj[n1].iter().find(|&&(nb, _)| nb == n2) {
                    mol.bonds[bid].order = BondOrder::Single;
                }
                if let Some(&(_, bid)) = mol.adj[n2].iter().find(|&&(nb, _)| nb == n3) {
                    mol.bonds[bid].order = BondOrder::Double;
                }
                mol.add_bond(c_term, n1, BondOrder::Single, BondStereo::None)?;
                mol.add_bond(c_int, n3, BondOrder::Single, BondStereo::None)?;
            }
            CrossLink::DeltaAmide { .. } => {}
        }
    }

    // 6. Apply C-terminal modification if C-terminus was not consumed by cyclization
    if !c_term_consumed {
        let last = &stamped[last_idx];
        let last_seq = (last_idx + 1) as u16;
        match spec.c_cap {
            CCapKind::FreeAcid => {
                if let Some(m) = mol.atoms[last.oxt_id].pdb_meta.as_mut() {
                    m.atom_name = "OXT".to_string();
                }
            }
            CCapKind::Amide => {
                to_delete.insert(last.oxt_id);
                let mut nh2 = MolAtom::new(0, "N")?;
                nh2.pdb_meta = Some(PdbAtomMeta {
                    atom_name: "N".to_string(),
                    res_name: "NH2".to_string(),
                    res_seq: last_seq + 1,
                    chain: 'A',
                    is_hetatm: false,
                });
                let nid = mol.add_atom(nh2);
                mol.add_bond(last.c_id, nid, BondOrder::Single, BondStereo::None)?;
            }
            CCapKind::Alcohol => {
                // Convert C-terminal -C(=O)OH into -CH2OH (e.g. Octreotide Thr-ol): delete carbonyl =O!
                if let Some(oid) = last.o_id {
                    to_delete.insert(oid);
                }
                if let Some(m) = mol.atoms[last.oxt_id].pdb_meta.as_mut() {
                    m.atom_name = "O".to_string();
                }
                patch_chiral_seq(&mut mol, last.ca_id, last.c_id, last.c_id);
            }
            CCapKind::Aldehyde => {
                to_delete.insert(last.oxt_id);
            }
            CCapKind::MethylEster => {
                let mut me = MolAtom::new(0, "C")?;
                me.pdb_meta = Some(PdbAtomMeta {
                    atom_name: "CM".to_string(),
                    res_name: "OME".to_string(),
                    res_seq: last_seq + 1,
                    chain: 'A',
                    is_hetatm: true,
                });
                let mid = mol.add_atom(me);
                mol.add_bond(last.oxt_id, mid, BondOrder::Single, BondStereo::None)?;
            }
            CCapKind::EthylAmide => {
                to_delete.insert(last.oxt_id);
                let mut n_e = MolAtom::new(0, "N")?;
                n_e.pdb_meta = Some(PdbAtomMeta {
                    atom_name: "N".to_string(),
                    res_name: "NHE".to_string(),
                    res_seq: last_seq + 1,
                    chain: 'A',
                    is_hetatm: true,
                });
                let nid = mol.add_atom(n_e);
                let c1 = mol.add_atom(MolAtom::new(0, "C")?);
                let c2 = mol.add_atom(MolAtom::new(0, "C")?);
                mol.add_bond(last.c_id, nid, BondOrder::Single, BondStereo::None)?;
                mol.add_bond(nid, c1, BondOrder::Single, BondStereo::None)?;
                mol.add_bond(c1, c2, BondOrder::Single, BondStereo::None)?;
            }
        }
    }

    // 7. Compact deleted leaving atoms and remap residue indices
    let remap = delete_atoms_remap(&mut mol, &to_delete);
    let res_atoms: Vec<ResidueAtoms> = stamped
        .iter()
        .filter_map(|sr| {
            Some(ResidueAtoms {
                n_id: remap[sr.n_id]?,
                ca_id: remap[sr.ca_id]?,
                c_id: remap[sr.c_id]?,
                o_id: sr.o_id.and_then(|o| remap[o]),
                all_atoms: sr.all_atoms.iter().filter_map(|&a| remap[a]).collect(),
                is_d: sr.is_d,
            })
        })
        .collect();

    // 8. Run full chemical perception (SSSR macrocycles, Hückel aromaticity, exact H, CIP R/S & E/Z)
    mol.perceive()?;

    // 9. Generate clean macrocycle/linear 2D and chiral-preserving 3D coordinates
    layout_peptide_2d(&mut mol, spec, &res_atoms);
    layout_peptide_3d(&mut mol, spec, &res_atoms);

    let notation = spec.to_notation();
    let helm = spec.to_helm(cat);
    let canonical_smiles = mol.canonical_smiles();

    let display_name = spec
        .name
        .clone()
        .unwrap_or_else(|| format!("NRP {notation}"));
    mol.name = Some(display_name);
    mol.properties
        .insert("PEPTIDE_SEQUENCE".to_string(), notation.clone());
    mol.properties
        .insert("PEPTIDE_HELM".to_string(), helm.clone());
    mol.properties.insert(
        "PEPTIDE_RESIDUES".to_string(),
        spec.residues.len().to_string(),
    );
    mol.properties.insert(
        "SECONDARY_STRUCTURE".to_string(),
        spec.conformation.label().to_string(),
    );
    mol.properties.insert(
        "PEPTIDE_CYCLIC".to_string(),
        if spec.is_cyclic() { "true" } else { "false" }.to_string(),
    );

    Ok(BuiltPeptide {
        spec: spec.clone(),
        mol,
        notation,
        helm,
        canonical_smiles,
    })
}

/// Convenience one-step builder from raw notation, preset ID, or HELM string.
pub fn build_peptide_or_nrp(
    input: &str,
    conformation: crate::notation::PeptideConformation,
    cat: &MonomerCatalog,
) -> Result<BuiltPeptide, PeptideError> {
    let spec = crate::notation::parse_peptide_spec(input, cat, conformation)?;
    build_peptide_from_spec(&spec, cat)
}
