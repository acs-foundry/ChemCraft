use std::collections::{BTreeSet, VecDeque};
use std::f64::consts::PI;

use mol_graph::{
    BondCisTrans, BondOrder, BondStereo, MolAtom, MolError, MolGraph, highest_cip_end_nbr,
};

/// Standard 2D chemistry coordinate bond length.
pub const STANDARD_BOND_LEN: f64 = 1.50;

/// Built-in ring templates for interactive ChemDraw/VectorCraft-style stamping.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RingTemplate {
    Cyclopropane,
    Cyclobutane,
    Cyclopentane,
    Cyclohexane,
    Benzene,
    Naphthalene,
    Pyridine,
    Indole,
    Tetrazole,
}

impl RingTemplate {
    pub fn label(self) -> &'static str {
        match self {
            Self::Cyclopropane => "Cyclopropane",
            Self::Cyclobutane => "Cyclobutane",
            Self::Cyclopentane => "Cyclopentane",
            Self::Cyclohexane => "Cyclohexane",
            Self::Benzene => "Benzene",
            Self::Naphthalene => "Naphthalene",
            Self::Pyridine => "Pyridine",
            Self::Indole => "Indole",
            Self::Tetrazole => "Tetrazole",
        }
    }

    pub fn smiles(self) -> &'static str {
        match self {
            Self::Cyclopropane => "C1CC1",
            Self::Cyclobutane => "C1CCC1",
            Self::Cyclopentane => "C1CCCC1",
            Self::Cyclohexane => "C1CCCCC1",
            Self::Benzene => "c1ccccc1",
            Self::Naphthalene => "c1ccc2ccccc2c1",
            Self::Pyridine => "c1ccncc1",
            Self::Indole => "c1ccc2[nH]ccc2c1",
            Self::Tetrazole => "c1nnn[nH]1",
        }
    }

    pub fn polygon_sides(self) -> usize {
        match self {
            Self::Cyclopropane => 3,
            Self::Cyclobutane => 4,
            Self::Cyclopentane | Self::Tetrazole => 5,
            Self::Cyclohexane
            | Self::Benzene
            | Self::Pyridine
            | Self::Naphthalene
            | Self::Indole => 6,
        }
    }
}

/// Magnetic snap target when stamping a ring onto the drawing canvas.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SnapTarget {
    Free {
        x: f64,
        y: f64,
    },
    Atom {
        atom_idx: usize,
    },
    Bond {
        bond_idx: usize,
        cursor_x: f64,
        cursor_y: f64,
    },
}

/// Ensure all atoms have valid, finite 2D coordinates; generate if missing.
pub fn ensure_2d_coords(mol: &mut MolGraph) {
    let missing_or_bad = mol.atoms.iter().any(|a| match a.coords_2d {
        Some([x, y]) => !x.is_finite() || !y.is_finite(),
        None => true,
    });
    if missing_or_bad {
        generate_2d_coords(mol);
    }
}

/// Clean-room 2D coordinate layout generator handling fused ring systems (e.g. steroids,
/// indoles, purines), spiro junctions, acyclic chains, and multi-component mixtures.
pub fn generate_2d_coords(mol: &mut MolGraph) {
    let n = mol.atoms.len();
    if n == 0 {
        return;
    }
    let mut coords: Vec<Option<[f64; 2]>> = vec![None; n];
    let comps = connected_components(mol);
    let mut cursor_x = 0.0;

    for comp in comps {
        if comp.len() == 1 {
            coords[comp[0]] = Some([cursor_x, 0.0]);
            cursor_x += STANDARD_BOND_LEN * 2.2;
            continue;
        }
        let comp_set: BTreeSet<usize> = comp.iter().copied().collect();
        let comp_rings: Vec<&Vec<usize>> = mol
            .rings
            .iter()
            .filter(|r| !r.is_empty() && r.iter().all(|v| comp_set.contains(v)))
            .collect();

        // Place the primary ring cluster first (rings sharing vertices); pendant rings are grown attached to their linker!
        let mut ring_placed = vec![false; comp_rings.len()];
        if !comp_rings.is_empty() {
            let first_ri = (0..comp_rings.len())
                .max_by_key(|&i| comp_rings[i].len())
                .unwrap_or(0);
            place_regular_ring(comp_rings[first_ri], [cursor_x, 0.0], 0.0, &mut coords);
            ring_placed[first_ri] = true;
            expand_connected_rings(mol, &comp_rings, &mut ring_placed, &mut coords, cursor_x);

            // Run Tutte planar cage solver on compact bridged cages (<= 10 ring atoms, all rings in 1 block:
            // covers Cubane, Prismane, Norbornane, Camphor, Adamantane, Methyl-Adamantane; excludes Morphine/Artemisinin/Drugs).
            let ring_atoms: BTreeSet<usize> = comp_rings
                .iter()
                .enumerate()
                .filter(|&(i, _)| ring_placed[i])
                .flat_map(|(_, r)| r.iter().copied())
                .collect();
            let ring_atom_list: Vec<usize> = ring_atoms.iter().copied().collect();
            let all_rings_in_block = ring_placed.iter().all(|&p| p);
            let is_bridged = (0..comp_rings.len()).any(|i| {
                let si: BTreeSet<usize> = comp_rings[i].iter().copied().collect();
                comp_rings
                    .iter()
                    .skip(i + 1)
                    .any(|rj| rj.iter().filter(|v| si.contains(v)).count() >= 3)
            });
            let eligible_cage = all_rings_in_block
                && ((ring_atom_list.len() == comp.len() && ring_atom_list.len() <= 14)
                    || (is_bridged && ring_atom_list.len() <= 10));
            let base_crossings = if eligible_cage {
                count_subset_crossings(mol, &ring_atoms, &coords)
            } else {
                0
            };
            if eligible_cage
                && base_crossings > 0
                && let Some((cage_coords, cage_crossings)) =
                    solve_bridged_tutte_cage(mol, &comp_rings, &ring_atom_list, cursor_x)
                && cage_crossings < base_crossings
            {
                for (v, pt) in cage_coords {
                    coords[v] = Some(pt);
                }
            }
        } else {
            // Acyclic component root
            let root = comp[0];
            coords[root] = Some([cursor_x, 0.0]);
        }

        // BFS place chains & pendant ring clusters outward along each attaching bond
        let mut q = VecDeque::new();
        for &v in &comp {
            if coords[v].is_some() {
                q.push_back(v);
            }
        }
        while let Some(u) = q.pop_front() {
            let pu = coords[u].unwrap_or([0.0, 0.0]);
            let placed_nbrs: Vec<usize> = mol.adj[u]
                .iter()
                .map(|&(v, _)| v)
                .filter(|&v| coords[v].is_some())
                .collect();
            let unplaced: Vec<usize> = mol.adj[u]
                .iter()
                .map(|&(v, _)| v)
                .filter(|&v| coords[v].is_none())
                .collect();
            if unplaced.is_empty() {
                continue;
            }

            let dbl_parent = placed_nbrs.iter().copied().find_map(|pn| {
                mol.adj[u].iter().find(|&&(nb, bid)| {
                    nb == pn
                        && mol.bonds[bid].order == BondOrder::Double
                        && matches!(mol.bonds[bid].cis_trans, BondCisTrans::E | BondCisTrans::Z)
                })
            });

            if let Some(&(v_dbl, dbl_bid)) = dbl_parent {
                let pv = coords[v_dbl].unwrap_or([0.0, 0.0]);
                let dx = pu[0] - pv[0];
                let dy = pu[1] - pv[1];
                let alpha = dy.atan2(dx);
                let len = dx.hypot(dy).max(1e-6);
                let nx = -dy / len;
                let ny = dx / len;
                let hi_other = highest_cip_end_nbr(mol, v_dbl, u);
                let hi_self = highest_cip_end_nbr(mol, u, v_dbl);
                let other_side = hi_other
                    .and_then(|ho| coords[ho])
                    .map_or(1.0, |pho| (pho[0] - pv[0]) * nx + (pho[1] - pv[1]) * ny);
                let want_self_hi_sign = if mol.bonds[dbl_bid].cis_trans == BondCisTrans::Z {
                    if other_side >= 0.0 { 1.0 } else { -1.0 }
                } else if other_side >= 0.0 {
                    -1.0
                } else {
                    1.0
                };
                for (i, v) in unplaced.into_iter().enumerate() {
                    let is_hi = hi_self == Some(v) || (hi_self.is_none() && i == 0);
                    let target_sign = if is_hi {
                        want_self_hi_sign
                    } else {
                        -want_self_hi_sign
                    };
                    let ang = if target_sign > 0.0 {
                        alpha + PI / 3.0
                    } else {
                        alpha - PI / 3.0
                    };
                    coords[v] = Some([
                        pu[0] + STANDARD_BOND_LEN * ang.cos(),
                        pu[1] + STANDARD_BOND_LEN * ang.sin(),
                    ]);
                    q.push_back(v);
                    expand_ring_cluster_if_needed(
                        mol,
                        &comp_rings,
                        &mut ring_placed,
                        &mut coords,
                        cursor_x,
                        &mut q,
                    );
                }
                continue;
            }

            let base_ang = if placed_nbrs.is_empty() {
                PI / 6.0
            } else if placed_nbrs.len() == 1 {
                let pp = coords[placed_nbrs[0]].unwrap_or(pu);
                let dx = pp[0] - pu[0];
                let dy = pp[1] - pu[1];
                (-dy).atan2(-dx)
            } else {
                let mut sum_x = 0.0;
                let mut sum_y = 0.0;
                let mut nbr_angles: Vec<f64> = Vec::with_capacity(placed_nbrs.len());
                for &pn in &placed_nbrs {
                    let pp = coords[pn].unwrap_or(pu);
                    let dx = pp[0] - pu[0];
                    let dy = pp[1] - pu[1];
                    let l = dx.hypot(dy).max(1e-6);
                    sum_x += dx / l;
                    sum_y += dy / l;
                    nbr_angles.push(dy.atan2(dx));
                }
                let neg_len = sum_x.hypot(sum_y);
                let neg_ang = (-sum_y).atan2(-sum_x);
                let min_sep = nbr_angles
                    .iter()
                    .map(|&a| {
                        let mut d = (neg_ang - a).abs();
                        while d > PI {
                            d = (2.0 * PI - d).abs();
                        }
                        d
                    })
                    .fold(PI, f64::min);

                let has_exocyclic_placed = placed_nbrs.iter().any(|&pn| {
                    mol.adj[u]
                        .iter()
                        .find(|&&(nb, _)| nb == pn)
                        .is_some_and(|&(_, bid)| !mol.bonds[bid].in_ring)
                });

                // When a ring atom has >= 3 placed neighbors including an exocyclic substituent
                // (e.g. gem-disubstituted 4-membered azetidine C3 in Baricitinib where 2 ring
                // edges at ±45° outweigh the exocyclic bond at 180°), `-sum` folds directly onto
                // the exocyclic bond; bisect the widest open sector instead.
                if placed_nbrs.len() >= 3 && has_exocyclic_placed && (neg_len < 1e-4 || min_sep < 0.35) {
                    nbr_angles.sort_by(f64::total_cmp);
                    let mut best_gap = -1.0_f64;
                    let mut best_mid = neg_ang;
                    let m = nbr_angles.len();
                    for idx in 0..m {
                        let a1 = nbr_angles[idx];
                        let a2 = if idx + 1 < m {
                            nbr_angles[idx + 1]
                        } else {
                            nbr_angles[0] + 2.0 * PI
                        };
                        let gap = a2 - a1;
                        if gap > best_gap {
                            best_gap = gap;
                            best_mid = a1 + 0.5 * gap;
                        }
                    }
                    best_mid
                } else {
                    neg_ang
                }
            };

            let k = unplaced.len();
            for (i, v) in unplaced.into_iter().enumerate() {
                if coords[v].is_some() {
                    continue;
                }
                let offset = if k == 1 {
                    let is_dbl_adj = mol.adj[v]
                        .iter()
                        .any(|&(_, bid)| mol.bonds[bid].order == BondOrder::Double);
                    if is_dbl_adj && placed_nbrs.len() == 1 {
                        PI / 6.0
                    } else {
                        0.0
                    }
                } else {
                    let spread = if k == 2 { PI * 2.0 / 3.0 } else { PI * 0.55 };
                    (i as f64 - (k as f64 - 1.0) * 0.5) * spread
                };
                let ang = base_ang + offset;
                coords[v] = Some([
                    pu[0] + STANDARD_BOND_LEN * ang.cos(),
                    pu[1] + STANDARD_BOND_LEN * ang.sin(),
                ]);
                q.push_back(v);
                expand_ring_cluster_if_needed(
                    mol,
                    &comp_rings,
                    &mut ring_placed,
                    &mut coords,
                    cursor_x,
                    &mut q,
                );
            }
        }

        // Track horizontal extent for side-by-side disconnected components
        let max_x = comp
            .iter()
            .filter_map(|&v| coords[v].map(|p| p[0]))
            .fold(cursor_x, f64::max);
        cursor_x = max_x + STANDARD_BOND_LEN * 2.5;
    }

    for (idx, c) in coords.into_iter().enumerate() {
        mol.atoms[idx].coords_2d = Some(c.unwrap_or([0.0, 0.0]));
    }
    clean_structure_2d(mol, 35);
}

/// Interactive geometric optimizer (`Clean Structure 2D`) that regularizes bond lengths
/// toward `STANDARD_BOND_LEN`, straightens bond angles, and separates non-bonded clashes.
pub fn clean_structure_2d(mol: &mut MolGraph, iterations: usize) {
    let n = mol.atoms.len();
    if n <= 1 {
        return;
    }
    for (idx, a) in mol.atoms.iter_mut().enumerate() {
        if a.coords_2d.is_none() {
            a.coords_2d = Some([idx as f64 * STANDARD_BOND_LEN, 0.0]);
        }
    }

    let mut pos: Vec<[f64; 2]> = mol
        .atoms
        .iter()
        .map(|a| a.coords_2d.unwrap_or([0.0, 0.0]))
        .collect();

    let iters = if n > 200 {
        iterations.min(8)
    } else {
        iterations
    };
    for _ in 0..iters {
        let mut disp = vec![[0.0_f64; 2]; n];

        // 1. Bond length spring to STANDARD_BOND_LEN
        for b in &mol.bonds {
            let dx = pos[b.v][0] - pos[b.u][0];
            let dy = pos[b.v][1] - pos[b.u][1];
            let dist = dx.hypot(dy).max(1e-4);
            let diff = (dist - STANDARD_BOND_LEN) / dist * 0.35;
            disp[b.u][0] += dx * diff;
            disp[b.u][1] += dy * diff;
            disp[b.v][0] -= dx * diff;
            disp[b.v][1] -= dy * diff;
        }

        // 2. Angular spread on degree-2/3 acyclic vertices (~120 deg)
        for u in 0..n {
            if mol.atoms[u].in_ring || mol.adj[u].len() < 2 {
                continue;
            }
            let nbrs = &mol.adj[u];
            for i in 0..nbrs.len() {
                for j in (i + 1)..nbrs.len() {
                    let a = nbrs[i].0;
                    let c = nbrs[j].0;
                    let dx = pos[c][0] - pos[a][0];
                    let dy = pos[c][1] - pos[a][1];
                    let d = dx.hypot(dy).max(1e-4);
                    let target = STANDARD_BOND_LEN * 3.0_f64.sqrt();
                    if d < target {
                        let push = (target - d) / d * 0.12;
                        disp[a][0] -= dx * push;
                        disp[a][1] -= dy * push;
                        disp[c][0] += dx * push;
                        disp[c][1] += dy * push;
                    }
                }
            }
        }

        // 3. Soft non-bonded clash avoidance (bounded to <= 400 atoms to guarantee interactive speed)
        if n <= 400 {
            let min_sep = STANDARD_BOND_LEN * 0.85;
            for i in 0..n {
                for j in (i + 1)..n {
                    let connected = mol.adj[i].iter().any(|&(v, _)| v == j);
                    if connected {
                        continue;
                    }
                    let dx = pos[j][0] - pos[i][0];
                    let dy = pos[j][1] - pos[i][1];
                    let d = dx.hypot(dy).max(1e-4);
                    if d < min_sep {
                        let push = (min_sep - d) / d * 0.25;
                        disp[i][0] -= dx * push;
                        disp[i][1] -= dy * push;
                        disp[j][0] += dx * push;
                        disp[j][1] += dy * push;
                    }
                }
            }
        }

        let prev_ring_crossings = if n <= 64 {
            count_ring_crossings(mol, &pos)
        } else {
            0
        };
        let mut next_pos = pos.clone();
        for (p, d) in next_pos.iter_mut().zip(disp.iter()) {
            let sx = d[0].clamp(-0.35, 0.35);
            let sy = d[1].clamp(-0.35, 0.35);
            if sx.is_finite() && sy.is_finite() {
                p[0] += sx;
                p[1] += sy;
            }
        }
        if n > 64 || count_ring_crossings(mol, &next_pos) <= prev_ring_crossings {
            pos = next_pos;
        } else {
            break;
        }
    }

    for (idx, p) in pos.into_iter().enumerate() {
        mol.atoms[idx].coords_2d = Some(p);
    }
}

type AtomCoordList = Vec<(usize, [f64; 2])>;

fn expand_connected_rings(
    mol: &MolGraph,
    comp_rings: &[&Vec<usize>],
    ring_placed: &mut [bool],
    coords: &mut [Option<[f64; 2]>],
    cursor_x: f64,
) -> Vec<usize> {
    let mut newly_placed = Vec::new();
    loop {
        let mut progressed = false;
        for (ri, ring) in comp_rings.iter().enumerate() {
            if ring_placed[ri] {
                continue;
            }
            let placed_verts: Vec<(usize, usize)> = ring
                .iter()
                .enumerate()
                .filter(|&(_, &v)| coords[v].is_some())
                .map(|(k, &v)| (k, v))
                .collect();
            if placed_verts.is_empty() {
                continue;
            }
            if placed_verts.len() >= 2 {
                if let Some((k1, k2)) = find_consecutive_placed(ring, coords) {
                    place_fused_ring(ring, k1, k2, coords);
                } else {
                    place_remaining_ring_verts(ring, coords);
                }
            } else {
                let (k0, v0) = placed_verts[0];
                let p0 = coords[v0].unwrap_or([cursor_x, 0.0]);
                let out_dir = avg_away_dir(mol, v0, coords);
                let radius = STANDARD_BOND_LEN / (2.0 * (PI / ring.len() as f64).sin());
                let center = [p0[0] + out_dir[0] * radius, p0[1] + out_dir[1] * radius];
                let base_angle = (p0[1] - center[1]).atan2(p0[0] - center[0]);
                let step = 2.0 * PI / ring.len() as f64;
                for (idx, &v) in ring.iter().enumerate() {
                    if coords[v].is_none() {
                        let delta = (idx as isize - k0 as isize) as f64;
                        let ang = base_angle + delta * step;
                        coords[v] = Some([
                            center[0] + radius * ang.cos(),
                            center[1] + radius * ang.sin(),
                        ]);
                    }
                }
            }
            for &v in *ring {
                newly_placed.push(v);
            }
            ring_placed[ri] = true;
            progressed = true;
        }
        if !progressed {
            break;
        }
    }
    newly_placed
}

fn expand_ring_cluster_if_needed(
    mol: &MolGraph,
    comp_rings: &[&Vec<usize>],
    ring_placed: &mut [bool],
    coords: &mut [Option<[f64; 2]>],
    cursor_x: f64,
    q: &mut VecDeque<usize>,
) {
    for v in expand_connected_rings(mol, comp_rings, ring_placed, coords, cursor_x) {
        q.push_back(v);
    }
}

fn solve_bridged_tutte_cage(
    mol: &MolGraph,
    comp_rings: &[&Vec<usize>],
    ring_atoms: &[usize],
    cursor_x: f64,
) -> Option<(AtomCoordList, usize)> {
    let core_set: BTreeSet<usize> = ring_atoms.iter().copied().collect();
    if !(4..=14).contains(&core_set.len()) {
        return None;
    }

    // Build candidate outer faces: SSSR rings + XOR simple cycles of adjacent pairs
    let mut cycles: Vec<Vec<usize>> = comp_rings
        .iter()
        .filter(|r| r.len() >= 4)
        .map(|r| (*r).clone())
        .collect();
    for i in 0..comp_rings.len() {
        let ei = ring_edge_set(mol, comp_rings[i]);
        for rj in comp_rings.iter().skip(i + 1) {
            let ej = ring_edge_set(mol, rj);
            let sym: BTreeSet<usize> = ei.symmetric_difference(&ej).copied().collect();
            if let Some(cyc) = ordered_cycle_from_edges(mol, &sym)
                && (4..=10).contains(&cyc.len())
            {
                cycles.push(cyc);
            }
        }
    }

    let mut best: Option<(usize, f64, AtomCoordList)> = None;
    for cyc in &cycles {
        let k = cyc.len();
        let radius = STANDARD_BOND_LEN * ((k as f64) / (2.0 * PI)).max(1.0) * 1.18;
        let start_ang = PI / (k as f64);
        let mut trial: Vec<[f64; 2]> = vec![[cursor_x, 0.0]; mol.atoms.len()];
        let mut pinned = vec![false; mol.atoms.len()];
        for (idx, &v) in cyc.iter().enumerate() {
            let ang = start_ang + 2.0 * PI * (idx as f64) / (k as f64);
            trial[v] = [cursor_x + radius * ang.cos(), radius * ang.sin()];
            pinned[v] = true;
        }
        // Tutte harmonic Jacobi relaxation for interior vertices
        for _ in 0..60 {
            let mut next = trial.clone();
            for &u in ring_atoms {
                if pinned[u] {
                    continue;
                }
                let mut sx = 0.0;
                let mut sy = 0.0;
                let mut cnt = 0.0;
                for &(v, _) in &mol.adj[u] {
                    if core_set.contains(&v) {
                        sx += trial[v][0];
                        sy += trial[v][1];
                        cnt += 1.0;
                    }
                }
                if cnt > 0.0 {
                    next[u] = [sx / cnt, sy / cnt];
                }
            }
            trial = next;
        }
        // Expand unpinned interior vertices mildly to balance polyhedral bond length CV
        for &u in ring_atoms {
            if !pinned[u] {
                trial[u] = [
                    cursor_x + (trial[u][0] - cursor_x) * 1.24,
                    trial[u][1] * 1.24,
                ];
            }
        }

        // Strict geometric invariants: zero bond collapse, bounded length ratio, minimum atom spacing
        let (min_bl, max_bl) = core_bond_len_bounds(mol, &core_set, &trial);
        if min_bl < STANDARD_BOND_LEN * 0.35 || max_bl / min_bl.max(1e-6) > 3.8 {
            continue;
        }
        let crossings = count_subset_crossings_slice(mol, &core_set, &trial);
        let min_dist = min_pairwise_dist(ring_atoms, &trial);
        if min_dist < STANDARD_BOND_LEN * 0.35 {
            continue;
        }
        let better = match &best {
            None => true,
            Some((bc, bd, _)) => crossings < *bc || (crossings == *bc && min_dist > *bd + 1e-4),
        };
        if better {
            let mapped = ring_atoms.iter().map(|&v| (v, trial[v])).collect();
            best = Some((crossings, min_dist, mapped));
        }
    }
    best.map(|(c, _, m)| (m, c))
}

fn core_bond_len_bounds(
    mol: &MolGraph,
    core_set: &BTreeSet<usize>,
    pos: &[[f64; 2]],
) -> (f64, f64) {
    let mut min_l = f64::INFINITY;
    let mut max_l = 0.0_f64;
    for b in &mol.bonds {
        if core_set.contains(&b.u) && core_set.contains(&b.v) {
            let d = (pos[b.u][0] - pos[b.v][0]).hypot(pos[b.u][1] - pos[b.v][1]);
            min_l = min_l.min(d);
            max_l = max_l.max(d);
        }
    }
    if min_l.is_finite() {
        (min_l, max_l)
    } else {
        (0.0, f64::INFINITY)
    }
}

fn ring_edge_set(mol: &MolGraph, ring: &[usize]) -> BTreeSet<usize> {
    let mut edges = BTreeSet::new();
    for w in 0..ring.len() {
        let u = ring[w];
        let v = ring[(w + 1) % ring.len()];
        if let Some(&(_, bid)) = mol.adj[u].iter().find(|&&(nb, _)| nb == v) {
            edges.insert(bid);
        }
    }
    edges
}

fn ordered_cycle_from_edges(mol: &MolGraph, edge_ids: &BTreeSet<usize>) -> Option<Vec<usize>> {
    if edge_ids.len() < 3 {
        return None;
    }
    let mut adj: std::collections::HashMap<usize, Vec<usize>> = std::collections::HashMap::new();
    for &bid in edge_ids {
        let b = &mol.bonds[bid];
        adj.entry(b.u).or_default().push(b.v);
        adj.entry(b.v).or_default().push(b.u);
    }
    if adj.len() != edge_ids.len() || adj.values().any(|ns| ns.len() != 2) {
        return None;
    }
    let start = *adj.keys().min()?;
    let mut cyc = vec![start];
    let mut prev = usize::MAX;
    let mut cur = start;
    for _ in 1..adj.len() {
        let next = *adj.get(&cur)?.iter().find(|&&n| n != prev)?;
        prev = cur;
        cur = next;
        cyc.push(cur);
    }
    (cyc.len() == adj.len()).then_some(cyc)
}

fn min_pairwise_dist(nodes: &[usize], pos: &[[f64; 2]]) -> f64 {
    let mut md = f64::INFINITY;
    for i in 0..nodes.len() {
        for j in (i + 1)..nodes.len() {
            let a = pos[nodes[i]];
            let b = pos[nodes[j]];
            let d = (a[0] - b[0]).hypot(a[1] - b[1]);
            if d < md {
                md = d;
            }
        }
    }
    if md.is_finite() { md } else { 0.0 }
}

fn count_subset_crossings(
    mol: &MolGraph,
    subset: &BTreeSet<usize>,
    coords: &[Option<[f64; 2]>],
) -> usize {
    let pts: Vec<[f64; 2]> = coords.iter().map(|c| c.unwrap_or([0.0, 0.0])).collect();
    count_subset_crossings_slice(mol, subset, &pts)
}

/// Count pairs of non-adjacent covalent bonds whose 2D segments intersect strictly.
pub fn count_bond_crossings(mol: &MolGraph) -> usize {
    let pos: Vec<[f64; 2]> = mol
        .atoms
        .iter()
        .map(|a| a.coords_2d.unwrap_or([0.0, 0.0]))
        .collect();
    count_all_crossings(mol, &pos)
}

fn count_all_crossings(mol: &MolGraph, pos: &[[f64; 2]]) -> usize {
    let edges: Vec<(usize, usize)> = mol.bonds.iter().map(|b| (b.u, b.v)).collect();
    count_edge_crossings(&edges, pos)
}

fn count_ring_crossings(mol: &MolGraph, pos: &[[f64; 2]]) -> usize {
    let edges: Vec<(usize, usize)> = mol
        .bonds
        .iter()
        .filter(|b| b.in_ring)
        .map(|b| (b.u, b.v))
        .collect();
    count_edge_crossings(&edges, pos)
}

fn count_subset_crossings_slice(
    mol: &MolGraph,
    subset: &BTreeSet<usize>,
    pos: &[[f64; 2]],
) -> usize {
    let edges: Vec<(usize, usize)> = mol
        .bonds
        .iter()
        .filter(|b| subset.contains(&b.u) && subset.contains(&b.v))
        .map(|b| (b.u, b.v))
        .collect();
    count_edge_crossings(&edges, pos)
}

fn count_edge_crossings(edges: &[(usize, usize)], pos: &[[f64; 2]]) -> usize {
    let mut crossings = 0usize;
    for i in 0..edges.len() {
        let (u1, v1) = edges[i];
        for &(u2, v2) in edges.iter().skip(i + 1) {
            if u1 == u2 || u1 == v2 || v1 == u2 || v1 == v2 {
                continue;
            }
            if segments_strictly_intersect(pos[u1], pos[v1], pos[u2], pos[v2]) {
                crossings += 1;
            }
        }
    }
    crossings
}

fn segments_strictly_intersect(a: [f64; 2], b: [f64; 2], c: [f64; 2], d: [f64; 2]) -> bool {
    let orient = |p: [f64; 2], q: [f64; 2], r: [f64; 2]| -> f64 {
        (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0])
    };
    let on_seg = |p: [f64; 2], q: [f64; 2], r: [f64; 2]| -> bool {
        q[0] <= p[0].max(r[0]) + 1e-5
            && q[0] + 1e-5 >= p[0].min(r[0])
            && q[1] <= p[1].max(r[1]) + 1e-5
            && q[1] + 1e-5 >= p[1].min(r[1])
    };
    let o1 = orient(a, b, c);
    let o2 = orient(a, b, d);
    let o3 = orient(c, d, a);
    let o4 = orient(c, d, b);
    let eps = 1e-7;
    if (o1 * o2 < -eps) && (o3 * o4 < -eps) {
        return true;
    }
    // Catch near-collinear overlapping segments (tolerance scales with bond length)
    let col_eps = 1e-4;
    (o1.abs() < col_eps && on_seg(a, c, b))
        || (o2.abs() < col_eps && on_seg(a, d, b))
        || (o3.abs() < col_eps && on_seg(c, a, d))
        || (o4.abs() < col_eps && on_seg(c, b, d))
}

/// Interactive ChemDraw/VectorCraft template stamper with snap-to-atom and fused-bond snapping.
pub fn stamp_ring_template(
    mol: &mut MolGraph,
    template: RingTemplate,
    target: SnapTarget,
) -> Result<(), MolError> {
    let parsed = mol_graph::parse_smiles(template.smiles())?;
    let mut template_mol = parsed;
    generate_2d_coords(&mut template_mol);

    // Simple templates (monocyclic) snap smoothly onto bond or atom
    match target {
        SnapTarget::Free { x, y } => {
            let center = centroid_of(&template_mol);
            let offset = mol.atoms.len();
            for mut a in template_mol.atoms {
                let p = a.coords_2d.unwrap_or([0.0, 0.0]);
                a.coords_2d = Some([p[0] - center[0] + x, p[1] - center[1] + y]);
                mol.add_atom(a);
            }
            for b in template_mol.bonds {
                mol.add_bond(offset + b.u, offset + b.v, b.order, b.stereo)?;
            }
        }
        SnapTarget::Atom { atom_idx } => {
            if atom_idx >= mol.atoms.len() {
                return Err(MolError::Graph("invalid snap atom".to_string()));
            }
            ensure_2d_coords(mol);
            let anchor = mol.atoms[atom_idx].coords_2d.unwrap_or([0.0, 0.0]);
            let dir = avg_away_from_mol(mol, atom_idx);
            // Spiro / single-bond attach atom 0 of template cleanly outward
            let p0 = template_mol.atoms[0].coords_2d.unwrap_or([0.0, 0.0]);
            let t_center = centroid_of(&template_mol);
            let t_ang = (t_center[1] - p0[1]).atan2(t_center[0] - p0[0]);
            let target_ang = dir[1].atan2(dir[0]);
            let rot = target_ang - t_ang;
            let attach_pos = [
                anchor[0] + dir[0] * STANDARD_BOND_LEN,
                anchor[1] + dir[1] * STANDARD_BOND_LEN,
            ];
            let offset = mol.atoms.len();
            for mut a in template_mol.atoms {
                let p = a.coords_2d.unwrap_or([0.0, 0.0]);
                let rx = (p[0] - p0[0]) * rot.cos() - (p[1] - p0[1]) * rot.sin();
                let ry = (p[0] - p0[0]) * rot.sin() + (p[1] - p0[1]) * rot.cos();
                a.coords_2d = Some([attach_pos[0] + rx, attach_pos[1] + ry]);
                mol.add_atom(a);
            }
            for b in template_mol.bonds {
                mol.add_bond(offset + b.u, offset + b.v, b.order, b.stereo)?;
            }
            mol.add_bond(atom_idx, offset, BondOrder::Single, BondStereo::None)?;
        }
        SnapTarget::Bond {
            bond_idx,
            cursor_x,
            cursor_y,
        } => {
            if bond_idx >= mol.bonds.len() {
                return Err(MolError::Graph("invalid snap bond".to_string()));
            }
            ensure_2d_coords(mol);
            let u = mol.bonds[bond_idx].u;
            let v = mol.bonds[bond_idx].v;
            let pu = mol.atoms[u].coords_2d.unwrap_or([0.0, 0.0]);
            let pv = mol.atoms[v].coords_2d.unwrap_or([STANDARD_BOND_LEN, 0.0]);
            let mid = [(pu[0] + pv[0]) * 0.5, (pu[1] + pv[1]) * 0.5];
            let dx = pv[0] - pu[0];
            let dy = pv[1] - pu[1];
            let len = dx.hypot(dy).max(1e-6);
            // Pick normal pointing toward cursor side!
            let mut nx = -dy / len;
            let mut ny = dx / len;
            let side_dot = (cursor_x - mid[0]) * nx + (cursor_y - mid[1]) * ny;
            if side_dot < 0.0 {
                nx = -nx;
                ny = -ny;
            }
            let sides = template.polygon_sides();
            let radius = len / (2.0 * (PI / sides as f64).sin());
            let apothem = (radius * radius - (len * 0.5).powi(2)).max(0.0).sqrt();
            let center = [mid[0] + nx * apothem, mid[1] + ny * apothem];

            // Build fused polygon sharing vertices u and v cleanly!
            let start_ang = (pv[1] - center[1]).atan2(pv[0] - center[0]);
            let u_ang = (pu[1] - center[1]).atan2(pu[0] - center[0]);
            let mut step = 2.0 * PI / sides as f64;
            let cross = (pv[0] - center[0]) * (pu[1] - center[1])
                - (pv[1] - center[1]) * (pu[0] - center[0]);
            if cross > 0.0 {
                // Ensure we step away from shared edge (v -> ... -> u)
                let _ = u_ang;
            } else {
                step = -step;
            }

            let is_arom = matches!(
                template,
                RingTemplate::Benzene
                    | RingTemplate::Pyridine
                    | RingTemplate::Tetrazole
                    | RingTemplate::Naphthalene
                    | RingTemplate::Indole
            );
            let border = if is_arom {
                BondOrder::Aromatic
            } else {
                BondOrder::Single
            };

            let mut prev = v;
            for k in 1..=(sides - 2) {
                let ang = start_ang + step * (k as f64);
                let sym = if template == RingTemplate::Pyridine && k == 2 {
                    "N"
                } else {
                    "C"
                };
                let mut na = MolAtom::new(0, sym)?;
                na.aromatic = is_arom;
                na.coords_2d = Some([
                    center[0] + radius * ang.cos(),
                    center[1] + radius * ang.sin(),
                ]);
                let nid = mol.add_atom(na);
                mol.add_bond(prev, nid, border, BondStereo::None)?;
                prev = nid;
            }
            mol.add_bond(prev, u, border, BondStereo::None)?;
        }
    }

    mol.perceive()?;
    Ok(())
}

fn connected_components(mol: &MolGraph) -> Vec<Vec<usize>> {
    let n = mol.atoms.len();
    let mut seen = vec![false; n];
    let mut out = Vec::new();
    for s in 0..n {
        if seen[s] {
            continue;
        }
        let mut comp = Vec::new();
        let mut q = VecDeque::new();
        seen[s] = true;
        q.push_back(s);
        while let Some(u) = q.pop_front() {
            comp.push(u);
            for &(v, _) in &mol.adj[u] {
                if !seen[v] {
                    seen[v] = true;
                    q.push_back(v);
                }
            }
        }
        out.push(comp);
    }
    out
}

fn place_regular_ring(
    ring: &[usize],
    center: [f64; 2],
    start_angle: f64,
    coords: &mut [Option<[f64; 2]>],
) {
    let k = ring.len() as f64;
    let radius = STANDARD_BOND_LEN / (2.0 * (PI / k).sin());
    for (i, &v) in ring.iter().enumerate() {
        let ang = start_angle + 2.0 * PI * (i as f64) / k;
        coords[v] = Some([
            center[0] + radius * ang.cos(),
            center[1] + radius * ang.sin(),
        ]);
    }
}

fn find_consecutive_placed(ring: &[usize], coords: &[Option<[f64; 2]>]) -> Option<(usize, usize)> {
    let m = ring.len();
    for i in 0..m {
        let j = (i + 1) % m;
        if coords[ring[i]].is_some() && coords[ring[j]].is_some() {
            return Some((i, j));
        }
    }
    None
}

fn place_fused_ring(ring: &[usize], k1: usize, k2: usize, coords: &mut [Option<[f64; 2]>]) {
    let p1 = coords[ring[k1]].unwrap_or([0.0, 0.0]);
    let p2 = coords[ring[k2]].unwrap_or([STANDARD_BOND_LEN, 0.0]);
    let mid = [(p1[0] + p2[0]) * 0.5, (p1[1] + p2[1]) * 0.5];
    let dx = p2[0] - p1[0];
    let dy = p2[1] - p1[1];
    let edge_len = dx.hypot(dy).max(1e-6);

    // Choose normal away from existing placed neighbors of p1/p2
    let n1 = [-dy / edge_len, dx / edge_len];
    let n2 = [dy / edge_len, -dx / edge_len];
    let existing_centroid = placed_centroid(coords, [p1, p2]);
    let d1 = (mid[0] + n1[0] - existing_centroid[0]).hypot(mid[1] + n1[1] - existing_centroid[1]);
    let d2 = (mid[0] + n2[0] - existing_centroid[0]).hypot(mid[1] + n2[1] - existing_centroid[1]);
    let normal = if d1 >= d2 { n1 } else { n2 };

    let sides = ring.len() as f64;
    let radius = edge_len / (2.0 * (PI / sides).sin());
    let apothem = (radius * radius - (edge_len * 0.5).powi(2)).max(0.0).sqrt();
    let center = [mid[0] + normal[0] * apothem, mid[1] + normal[1] * apothem];

    let a1 = (p1[1] - center[1]).atan2(p1[0] - center[0]);
    let a2 = (p2[1] - center[1]).atan2(p2[0] - center[0]);
    let mut step = 2.0 * PI / sides;
    let diff = normalize_angle(a2 - a1);
    if diff < 0.0 {
        step = -step;
    }

    for i in 0..ring.len() {
        let v = ring[i];
        if coords[v].is_some() {
            continue;
        }
        let dist_from_k1 = (i as isize - k1 as isize).rem_euclid(ring.len() as isize) as f64;
        let ang = a1 + dist_from_k1 * step;
        coords[v] = Some([
            center[0] + radius * ang.cos(),
            center[1] + radius * ang.sin(),
        ]);
    }
}

fn place_remaining_ring_verts(ring: &[usize], coords: &mut [Option<[f64; 2]>]) {
    let mut sum = [0.0, 0.0];
    let mut cnt = 0.0;
    for &v in ring {
        if let Some(p) = coords[v] {
            sum[0] += p[0];
            sum[1] += p[1];
            cnt += 1.0;
        }
    }
    let c = if cnt > 0.0 {
        [sum[0] / cnt, sum[1] / cnt]
    } else {
        [0.0, 0.0]
    };
    for (idx, &v) in ring.iter().enumerate() {
        if coords[v].is_none() {
            let ang = 2.0 * PI * (idx as f64) / (ring.len() as f64);
            coords[v] = Some([
                c[0] + STANDARD_BOND_LEN * ang.cos(),
                c[1] + STANDARD_BOND_LEN * ang.sin(),
            ]);
        }
    }
}

fn placed_centroid(coords: &[Option<[f64; 2]>], fallback: [[f64; 2]; 2]) -> [f64; 2] {
    let mut sx = 0.0;
    let mut sy = 0.0;
    let mut n = 0.0;
    for c in coords.iter().flatten() {
        sx += c[0];
        sy += c[1];
        n += 1.0;
    }
    if n > 2.0 {
        [sx / n, sy / n]
    } else {
        [
            (fallback[0][0] + fallback[1][0]) * 0.5,
            (fallback[0][1] + fallback[1][1]) * 0.5,
        ]
    }
}

fn avg_away_dir(mol: &MolGraph, u: usize, coords: &[Option<[f64; 2]>]) -> [f64; 2] {
    let pu = coords[u].unwrap_or([0.0, 0.0]);
    let mut sx = 0.0;
    let mut sy = 0.0;
    for &(v, _) in &mol.adj[u] {
        if let Some(pv) = coords[v] {
            let dx = pv[0] - pu[0];
            let dy = pv[1] - pu[1];
            let d = dx.hypot(dy).max(1e-6);
            sx += dx / d;
            sy += dy / d;
        }
    }
    let len = sx.hypot(sy);
    if len < 1e-4 {
        [1.0, 0.0]
    } else {
        [-sx / len, -sy / len]
    }
}

fn avg_away_from_mol(mol: &MolGraph, u: usize) -> [f64; 2] {
    let pu = mol.atoms[u].coords_2d.unwrap_or([0.0, 0.0]);
    let mut sx = 0.0;
    let mut sy = 0.0;
    for &(v, _) in &mol.adj[u] {
        if let Some(pv) = mol.atoms[v].coords_2d {
            let dx = pv[0] - pu[0];
            let dy = pv[1] - pu[1];
            let d = dx.hypot(dy).max(1e-6);
            sx += dx / d;
            sy += dy / d;
        }
    }
    let len = sx.hypot(sy);
    if len < 1e-4 {
        [1.0, 0.0]
    } else {
        [-sx / len, -sy / len]
    }
}

fn centroid_of(mol: &MolGraph) -> [f64; 2] {
    let mut sx = 0.0;
    let mut sy = 0.0;
    let mut cnt = 0.0;
    for a in &mol.atoms {
        if let Some(p) = a.coords_2d {
            sx += p[0];
            sy += p[1];
            cnt += 1.0;
        }
    }
    if cnt > 0.0 {
        [sx / cnt, sy / cnt]
    } else {
        [0.0, 0.0]
    }
}

fn normalize_angle(mut a: f64) -> f64 {
    while a > PI {
        a -= 2.0 * PI;
    }
    while a < -PI {
        a += 2.0 * PI;
    }
    a
}
