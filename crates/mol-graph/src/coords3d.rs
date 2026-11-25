//! Clean-room 3D coordinate generation, steric/geometric refinement, PDB & 3D SDF export.

use std::collections::{HashMap, VecDeque};
use std::f64::consts::PI;
use std::fmt::Write;

use crate::mol::{BondCisTrans, BondOrder, BondStereo, ChiralTag, MolGraph, highest_cip_end_nbr};

/// Covalent radius in Angstroms for common organic/medicinal elements.
pub fn covalent_radius(symbol: &str) -> f64 {
    match symbol {
        "H" => 0.31,
        "B" => 0.84,
        "C" => 0.76,
        "N" => 0.71,
        "O" => 0.66,
        "F" => 0.57,
        "Si" => 1.11,
        "P" => 1.07,
        "S" => 1.05,
        "Cl" => 1.02,
        "Se" => 1.20,
        "Br" => 1.20,
        "I" => 1.39,
        _ => 0.77,
    }
}

/// Ideal covalent bond length (Angstroms) between atoms `u` and `v`.
pub fn ideal_bond_length(mol: &MolGraph, u: usize, v: usize, order: BondOrder) -> f64 {
    let base = covalent_radius(&mol.atoms[u].symbol) + covalent_radius(&mol.atoms[v].symbol);
    let scale = match order {
        BondOrder::Single | BondOrder::Dative => 1.0,
        BondOrder::Aromatic => 0.915,
        BondOrder::Double => 0.875,
        BondOrder::Triple => 0.79,
    };
    (base * scale).clamp(0.95, 2.50)
}

/// Return `coords_3d` for every atom, generating clean 3D coordinates on the fly if missing.
pub fn compute_3d_coords(mol: &MolGraph) -> Vec<[f64; 3]> {
    if mol.atoms.is_empty() {
        return Vec::new();
    }
    let all_present = mol
        .atoms
        .iter()
        .all(|a| a.coords_3d.is_some_and(|p| p.iter().all(|c| c.is_finite())));
    if all_present {
        return mol
            .atoms
            .iter()
            .map(|a| a.coords_3d.unwrap_or([0.0, 0.0, 0.0]))
            .collect();
    }
    generate_3d_vec(mol)
}

/// Populate `coords_3d` on `mol` only if any atom lacks valid 3D coordinates.
pub fn ensure_3d_coords(mol: &mut MolGraph) {
    if mol.atoms.is_empty() {
        return;
    }
    let all_present = mol
        .atoms
        .iter()
        .all(|a| a.coords_3d.is_some_and(|p| p.iter().all(|c| c.is_finite())));
    if !all_present {
        generate_3d_coords(mol);
    }
}

/// Regenerate clean 3D conformer coordinates for all atoms in `mol`.
pub fn generate_3d_coords(mol: &mut MolGraph) {
    let coords = generate_3d_vec(mol);
    for (atom, c) in mol.atoms.iter_mut().zip(coords) {
        atom.coords_3d = Some(c);
    }
}

fn atom_is_sp(mol: &MolGraph, i: usize) -> bool {
    // Hypervalent sulfur (sulfone/sulfonamide/sulfate: S(=O)(=O)...) & phosphorus (P(=O)...)
    // with >= 3 neighbors are tetrahedral/pyramidal, never linear sp!
    if matches!(mol.atoms[i].symbol.as_str(), "S" | "P") && mol.adj[i].len() >= 3 {
        return false;
    }
    mol.adj[i]
        .iter()
        .any(|&(_, bid)| mol.bonds[bid].order == BondOrder::Triple)
        || (mol.adj[i]
            .iter()
            .filter(|&&(_, bid)| mol.bonds[bid].order == BondOrder::Double)
            .count()
            >= 2)
}

fn atom_is_sp2(mol: &MolGraph, i: usize) -> bool {
    // Hypervalent sulfur (sulfone/sulfonamide/sulfoxide) & phosphorus (phosphate) are tetrahedral/pyramidal
    if matches!(mol.atoms[i].symbol.as_str(), "S" | "P") && mol.adj[i].len() >= 3 {
        return false;
    }
    if mol.atoms[i].aromatic {
        return true;
    }
    if mol.adj[i].iter().any(|&(_, bid)| {
        matches!(
            mol.bonds[bid].order,
            BondOrder::Double | BondOrder::Aromatic
        )
    }) {
        return true;
    }
    // Amide/ester/aniline-like resonance nitrogen or oxygen conjugated to sp2 neighbor
    if matches!(mol.atoms[i].symbol.as_str(), "N" | "O" | "S")
        && mol.adj[i].len() >= 2
        && mol.adj[i].iter().any(|&(nb, _)| {
            mol.atoms[nb].aromatic
                || mol.adj[nb]
                    .iter()
                    .any(|&(_, b2)| mol.bonds[b2].order == BondOrder::Double)
        })
    {
        return true;
    }
    false
}

fn ideal_angle_rad(mol: &MolGraph, center: usize, n1: usize, n2: usize) -> f64 {
    // Check if both neighbors share a tiny ring with center (cyclopropane/cyclobutane)
    for ring in &mol.rings {
        if ring.len() <= 5 && ring.contains(&center) && ring.contains(&n1) && ring.contains(&n2) {
            return PI * (ring.len() as f64 - 2.0) / ring.len() as f64;
        }
    }
    if atom_is_sp(mol, center) {
        PI
    } else if atom_is_sp2(mol, center) {
        2.0 * PI / 3.0
    } else {
        // Tetrahedral 109.471 deg
        109.471_f64.to_radians()
    }
}

fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn norm3(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

fn unit3(a: [f64; 3], fallback: [f64; 3]) -> [f64; 3] {
    let n = norm3(a);
    if n > 1e-9 && n.is_finite() {
        [a[0] / n, a[1] / n, a[2] / n]
    } else {
        fallback
    }
}

fn ortho_pair(axis: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let ax = unit3(axis, [1.0, 0.0, 0.0]);
    let ref_v = if ax[2].abs() < 0.8 {
        [0.0, 0.0, 1.0]
    } else {
        [0.0, 1.0, 0.0]
    };
    let u = unit3(cross3(ax, ref_v), [0.0, 1.0, 0.0]);
    let v = unit3(cross3(ax, u), [0.0, 0.0, 1.0]);
    (u, v)
}

fn rodrigues_rotate(v: [f64; 3], axis: [f64; 3], theta: f64) -> [f64; 3] {
    let k = unit3(axis, [0.0, 0.0, 1.0]);
    let (s, c) = theta.sin_cos();
    let kv = dot3(k, v);
    let kxv = cross3(k, v);
    [
        v[0] * c + kxv[0] * s + k[0] * kv * (1.0 - c),
        v[1] * c + kxv[1] * s + k[1] * kv * (1.0 - c),
        v[2] * c + kxv[2] * s + k[2] * kv * (1.0 - c),
    ]
}

fn generate_3d_vec(mol: &MolGraph) -> Vec<[f64; 3]> {
    let n = mol.atoms.len();
    if n == 0 {
        return Vec::new();
    }

    // Seed positions cleanly: if good 2D coords exist, start XY from them with small non-zero Z
    // breaking symmetry so sp3 carbons/saturated rings pop smoothly into 3D instead of getting trapped at Z=0!
    let mut pos: Vec<Option<[f64; 3]>> = vec![None; n];

    // 1. Place rings in 3D first (regular polygon in XY + puckering for saturated rings)
    let mut ring_placed = vec![false; mol.rings.len()];
    // Sort largest or aromatic first
    let mut ring_order: Vec<usize> = (0..mol.rings.len()).collect();
    ring_order.sort_by_key(|&ri| {
        let r = &mol.rings[ri];
        let arom = r.iter().all(|&a| mol.atoms[a].aromatic);
        (if arom { 0usize } else { 1usize }, r.len())
    });

    while let Some(&first_ri) = ring_order.iter().find(|&&ri| !ring_placed[ri]) {
        // Place seed ring
        place_seed_ring(mol, &mol.rings[first_ri], &mut pos);
        ring_placed[first_ri] = true;

        // Expand fused/spiro rings
        loop {
            let mut progress = false;
            for &ri in &ring_order {
                if ring_placed[ri] {
                    continue;
                }
                let r = &mol.rings[ri];
                let placed_cnt = r.iter().filter(|&&a| pos[a].is_some()).count();
                if placed_cnt >= 2 {
                    place_fused_ring(mol, r, &mut pos);
                    ring_placed[ri] = true;
                    progress = true;
                } else if placed_cnt == 1 {
                    place_spiro_ring(mol, r, &mut pos);
                    ring_placed[ri] = true;
                    progress = true;
                }
            }
            if !progress {
                break;
            }
        }
    }

    // 2. BFS remaining atoms component-by-component
    let mut comp_offset_x = 0.0_f64;
    for start in 0..n {
        if pos[start].is_none() {
            // Check if this atom belongs to a component that already has some placed atoms (e.g. a ring)
            let comp_atoms = bfs_component(mol, start);
            let existing: Vec<usize> = comp_atoms
                .iter()
                .copied()
                .filter(|&a| pos[a].is_some())
                .collect();
            if existing.is_empty() {
                pos[start] = Some([comp_offset_x, 0.0, 0.0]);
            }
            expand_bfs_3d(mol, &comp_atoms, &mut pos);

            // Shift component so components never collide
            let mut min_x = f64::INFINITY;
            let mut max_x = f64::NEG_INFINITY;
            for &a in &comp_atoms {
                if let Some(p) = pos[a] {
                    min_x = min_x.min(p[0]);
                    max_x = max_x.max(p[0]);
                }
            }
            if comp_offset_x > 0.0 && min_x.is_finite() {
                let shift = comp_offset_x - min_x;
                for &a in &comp_atoms {
                    if let Some(p) = &mut pos[a] {
                        p[0] += shift;
                    }
                }
                comp_offset_x += (max_x - min_x).max(2.0) + 4.0;
            } else if max_x.is_finite() {
                comp_offset_x = max_x + 4.0;
            }
        }
    }

    let mut pts: Vec<[f64; 3]> = pos
        .into_iter()
        .enumerate()
        .map(|(i, p)| {
            p.unwrap_or_else(|| {
                mol.atoms[i]
                    .coords_2d
                    .map_or([0.0, 0.0, 0.0], |[x, y]| [x, y, 0.0])
            })
        })
        .collect();

    // Authentic out-of-plane pyramidal break ONLY on genuinely coplanar degree >= 3 sp3 atoms
    for i in 0..n {
        if !atom_is_sp(mol, i) && !atom_is_sp2(mol, i) && mol.adj[i].len() >= 3 {
            let p0 = pts[i];
            let u1 = unit3(sub3(pts[mol.adj[i][0].0], p0), [1.0, 0.0, 0.0]);
            let u2 = unit3(sub3(pts[mol.adj[i][1].0], p0), [0.0, 1.0, 0.0]);
            let u3 = unit3(sub3(pts[mol.adj[i][2].0], p0), [0.0, 0.0, 1.0]);
            let vol = dot3(cross3(u1, u2), u3).abs();
            if vol < 0.12 {
                let n_vec = unit3(cross3(u1, u2), [0.0, 0.0, 1.0]);
                let sign = if i.is_multiple_of(2) { 0.48 } else { -0.48 };
                pts[i] = add3(pts[i], scale3(n_vec, sign));
            }
        }
    }

    // For compact achiral bridged polycyclic cages (Adamantane, Cubane, etc.), push vertices
    // radially outward from the cage centroid so no bridgehead traps inside-out.
    inflate_bridged_cages(mol, &mut pts);

    let unpinned = vec![false; n];
    relax_3d_coords_pinned(mol, &mut pts, &unpinned, 120);
    enforce_chiral_signs(mol, &mut pts);
    relax_3d_coords_pinned(mol, &mut pts, &unpinned, 90);
    enforce_chiral_signs(mol, &mut pts);
    relax_3d_coords_pinned(mol, &mut pts, &unpinned, 45);

    // Center centroid at origin
    let mut cx = 0.0;
    let mut cy = 0.0;
    let mut cz = 0.0;
    for p in &pts {
        cx += p[0];
        cy += p[1];
        cz += p[2];
    }
    let inv = 1.0 / pts.len() as f64;
    for p in &mut pts {
        p[0] -= cx * inv;
        p[1] -= cy * inv;
        p[2] -= cz * inv;
    }
    pts
}

fn bfs_component(mol: &MolGraph, start: usize) -> Vec<usize> {
    let mut vis = vec![false; mol.atoms.len()];
    let mut q = VecDeque::new();
    let mut out = Vec::new();
    vis[start] = true;
    q.push_back(start);
    while let Some(u) = q.pop_front() {
        out.push(u);
        for &(v, _) in &mol.adj[u] {
            if !vis[v] {
                vis[v] = true;
                q.push_back(v);
            }
        }
    }
    out
}

fn place_seed_ring(mol: &MolGraph, ring: &[usize], pos: &mut [Option<[f64; 3]>]) {
    let k = ring.len().max(3);
    let planar = ring.iter().all(|&a| atom_is_sp2(mol, a));
    let bl = 1.44;
    let radius = bl / (2.0 * (PI / k as f64).sin()).max(0.2);
    for (idx, &a) in ring.iter().enumerate() {
        if pos[a].is_some() {
            continue;
        }
        let ang = 2.0 * PI * (idx as f64) / (k as f64);
        let z = if planar {
            0.0
        } else if k == 6 {
            if idx % 2 == 0 { 0.30 } else { -0.30 }
        } else if k == 5 && idx == 0 {
            0.34
        } else if k >= 7 {
            0.34 * (ang * 2.0).sin()
        } else {
            0.0
        };
        pos[a] = Some([radius * ang.cos(), radius * ang.sin(), z]);
    }
}

fn place_fused_ring(mol: &MolGraph, ring: &[usize], pos: &mut [Option<[f64; 3]>]) {
    let k = ring.len().max(3);
    // Find two consecutive atoms along ring already placed
    let mut anchor = None;
    for idx in 0..k {
        let a = ring[idx];
        let b = ring[(idx + 1) % k];
        if pos[a].is_some() && pos[b].is_some() {
            anchor = Some((idx, a, b));
            break;
        }
    }
    let Some((start_idx, a, b)) = anchor else {
        place_seed_ring(mol, ring, pos);
        return;
    };
    let Some(pa) = pos[a] else { return };
    let Some(pb) = pos[b] else { return };

    // Outward normal away from already-placed neighbors' centroid
    let mut cent = [0.0, 0.0, 0.0];
    let mut cnt = 0.0;
    for &(nb, _) in &mol.adj[a] {
        if nb != b
            && let Some(pn) = pos[nb]
        {
            cent = add3(cent, pn);
            cnt += 1.0;
        }
    }
    for &(nb, _) in &mol.adj[b] {
        if nb != a
            && let Some(pn) = pos[nb]
        {
            cent = add3(cent, pn);
            cnt += 1.0;
        }
    }
    let edge = sub3(pb, pa);
    let edge_len = norm3(edge).max(1.2);
    let mid = scale3(add3(pa, pb), 0.5);
    let raw_in = if cnt > 0.0 {
        sub3(mid, scale3(cent, 1.0 / cnt))
    } else {
        [0.0, 1.0, 0.0]
    };
    // Project raw_in perpendicular to edge
    let u_edge = unit3(edge, [1.0, 0.0, 0.0]);
    let perp = sub3(raw_in, scale3(u_edge, dot3(raw_in, u_edge)));
    let out_dir = unit3(perp, cross3(u_edge, [0.0, 0.0, 1.0]));
    let normal = unit3(cross3(u_edge, out_dir), [0.0, 0.0, 1.0]);

    let radius = edge_len / (2.0 * (PI / k as f64).sin()).max(0.2);
    let apothem = (radius * radius - 0.25 * edge_len * edge_len)
        .max(0.0)
        .sqrt();
    let ring_center = add3(mid, scale3(out_dir, apothem));

    let planar = ring.iter().all(|&atm| atom_is_sp2(mol, atm));
    // Angle of a and b around ring_center in plane (out_dir, u_edge)
    let va = sub3(pa, ring_center);
    let vb = sub3(pb, ring_center);
    let ang_a = dot3(va, u_edge).atan2(dot3(va, out_dir));
    let ang_b = dot3(vb, u_edge).atan2(dot3(vb, out_dir));
    let mut step = ang_b - ang_a;
    while step > PI {
        step -= 2.0 * PI;
    }
    while step < -PI {
        step += 2.0 * PI;
    }
    let sign = if step >= 0.0 { 1.0 } else { -1.0 };
    let delta = sign * 2.0 * PI / (k as f64);

    for offset in 0..k {
        let idx = (start_idx + offset) % k;
        let atm = ring[idx];
        if pos[atm].is_some() {
            continue;
        }
        let ang = ang_a + delta * (offset as f64);
        let in_plane = add3(
            ring_center,
            add3(
                scale3(out_dir, radius * ang.cos()),
                scale3(u_edge, radius * ang.sin()),
            ),
        );
        let puck = if planar {
            0.0
        } else if offset % 2 == 0 {
            0.20
        } else {
            -0.20
        };
        pos[atm] = Some(add3(in_plane, scale3(normal, puck)));
    }
}

fn place_spiro_ring(mol: &MolGraph, ring: &[usize], pos: &mut [Option<[f64; 3]>]) {
    let k = ring.len().max(3);
    let Some((anchor_idx, &anchor)) = ring.iter().enumerate().find(|(_, a)| pos[**a].is_some())
    else {
        return;
    };
    let Some(pa) = pos[anchor] else { return };
    let mut sum = [0.0, 0.0, 0.0];
    let mut cnt = 0.0;
    for &(nb, _) in &mol.adj[anchor] {
        if let Some(pn) = pos[nb] {
            sum = add3(sum, sub3(pn, pa));
            cnt += 1.0;
        }
    }
    let out_dir = if cnt > 0.0 {
        unit3(scale3(sum, -1.0), [1.0, 0.0, 0.0])
    } else {
        [1.0, 0.0, 0.0]
    };
    let (u, v) = ortho_pair(out_dir);
    let radius = 1.46 / (2.0 * (PI / k as f64).sin()).max(0.2);
    let center = add3(pa, scale3(out_dir, radius));
    for offset in 1..k {
        let idx = (anchor_idx + offset) % k;
        let atm = ring[idx];
        if pos[atm].is_some() {
            continue;
        }
        let ang = PI + 2.0 * PI * (offset as f64) / (k as f64);
        let p = add3(
            center,
            add3(
                scale3(out_dir, radius * ang.cos()),
                scale3(u, radius * ang.sin()),
            ),
        );
        let p = add3(p, scale3(v, if offset % 2 == 0 { 0.12 } else { -0.12 }));
        pos[atm] = Some(p);
    }
}

fn expand_bfs_3d(mol: &MolGraph, comp: &[usize], pos: &mut [Option<[f64; 3]>]) {
    let mut q = VecDeque::new();
    let mut in_q = vec![false; mol.atoms.len()];
    for &a in comp {
        if pos[a].is_some() {
            q.push_back(a);
            in_q[a] = true;
        }
    }

    while let Some(u) = q.pop_front() {
        let Some(pu) = pos[u] else { continue };
        let mut placed_nbrs = Vec::new();
        let mut unplaced_nbrs = Vec::new();
        for &(v, bid) in &mol.adj[u] {
            if pos[v].is_some() {
                placed_nbrs.push((v, bid));
            } else {
                unplaced_nbrs.push((v, bid));
            }
        }
        if unplaced_nbrs.is_empty() {
            continue;
        }

        // Acyclic E/Z double-bond placement if `u` sits at double bond end attached to placed `v_dbl`
        if let Some(&(v_dbl, dbl_bid)) = placed_nbrs.iter().find(|&&(_, bid)| {
            mol.bonds[bid].order == BondOrder::Double
                && matches!(mol.bonds[bid].cis_trans, BondCisTrans::E | BondCisTrans::Z)
        }) && let Some(pv) = pos[v_dbl]
        {
            let axis = unit3(sub3(pu, pv), [1.0, 0.0, 0.0]);
            // Find reference vector for `hi_other` perpendicular to axis
            let hi_other = highest_cip_end_nbr(mol, v_dbl, u);
            let ref_perp = hi_other
                .and_then(|ho| pos[ho])
                .map(|pho| {
                    let d = sub3(pho, pv);
                    let perp = sub3(d, scale3(axis, dot3(d, axis)));
                    unit3(perp, ortho_pair(axis).0)
                })
                .unwrap_or_else(|| ortho_pair(axis).0);
            let hi_self = highest_cip_end_nbr(mol, u, v_dbl);
            let want_hi_sign = if mol.bonds[dbl_bid].cis_trans == BondCisTrans::Z {
                1.0
            } else {
                -1.0
            };
            let (s_ang, c_ang) = (PI / 3.0).sin_cos();
            for (idx, &(w, bid)) in unplaced_nbrs.iter().enumerate() {
                let bl = ideal_bond_length(mol, u, w, mol.bonds[bid].order);
                let is_hi = hi_self == Some(w) || (hi_self.is_none() && idx == 0);
                let sign = if is_hi { want_hi_sign } else { -want_hi_sign };
                let dir = add3(scale3(axis, c_ang), scale3(ref_perp, sign * s_ang));
                pos[w] = Some(add3(pu, scale3(unit3(dir, [1.0, 0.0, 0.0]), bl)));
                if !in_q[w] {
                    q.push_back(w);
                    in_q[w] = true;
                }
            }
            continue;
        }

        // General hybrid placement
        let is_sp = atom_is_sp(mol, u);
        let is_sp2 = atom_is_sp2(mol, u);
        let back_dir = if placed_nbrs.is_empty() {
            [-1.0, 0.0, 0.0]
        } else {
            let mut sum = [0.0, 0.0, 0.0];
            for &(pn, _) in &placed_nbrs {
                if let Some(pp) = pos[pn] {
                    sum = add3(sum, unit3(sub3(pp, pu), [1.0, 0.0, 0.0]));
                }
            }
            unit3(sum, [-1.0, 0.0, 0.0])
        };
        let fwd_dir = scale3(back_dir, -1.0);

        // Try to match dihedral plane from grandparent to encourage trans zigzag conformations
        let ref_plane_vec = placed_nbrs.first().and_then(|&(parent, _)| {
            let pp = pos[parent]?;
            mol.adj[parent].iter().find_map(|&(gp, _)| {
                if gp != u {
                    let pgp = pos[gp]?;
                    let v = sub3(pgp, pp);
                    let perp = sub3(v, scale3(fwd_dir, dot3(v, fwd_dir)));
                    if norm3(perp) > 1e-4 {
                        return Some(unit3(perp, [0.0, 1.0, 0.0]));
                    }
                }
                None
            })
        });
        let (u_orth, v_orth) = match ref_plane_vec {
            Some(ref_v) => {
                let uo = unit3(ref_v, [0.0, 1.0, 0.0]);
                let vo = unit3(cross3(fwd_dir, uo), [0.0, 0.0, 1.0]);
                (uo, vo)
            }
            None => ortho_pair(fwd_dir),
        };

        let k = unplaced_nbrs.len();
        let tilt = if is_sp {
            0.0
        } else if is_sp2 {
            PI / 3.0 // 120 deg from back_dir -> cos(60) along fwd_dir
        } else {
            70.529_f64.to_radians() // 109.471 deg from back_dir
        };
        let (s_t, c_t) = tilt.sin_cos();

        // When 2 neighbors are already placed on a tetrahedral sp3/sulfone atom, span out of their plane!
        if !is_sp
            && !is_sp2
            && placed_nbrs.len() == 2
            && (k == 1 || k == 2)
            && let (Some(p1), Some(p2)) = (pos[placed_nbrs[0].0], pos[placed_nbrs[1].0])
        {
            let v1 = unit3(sub3(p1, pu), [1.0, 0.0, 0.0]);
            let v2 = unit3(sub3(p2, pu), [0.0, 1.0, 0.0]);
            let bis = unit3(scale3(add3(v1, v2), -1.0), fwd_dir);
            let perp_n = unit3(cross3(v1, v2), ortho_pair(bis).1);
            // Solve dot(c_in * bis, v1) = -1/3 so angle(new, v1) = angle(new, v2) = 109.47°
            let d_bv = dot3(bis, v1);
            let c_in = if d_bv.abs() > 0.15 {
                ((-1.0 / 3.0) / d_bv).clamp(0.05, 0.90)
            } else {
                0.577
            };
            let c_perp = (1.0 - c_in * c_in).max(0.0).sqrt();
            for (idx, &(w, bid)) in unplaced_nbrs.iter().enumerate() {
                let bl = ideal_bond_length(mol, u, w, mol.bonds[bid].order);
                let sign = if k == 2 {
                    if idx == 0 { 1.0 } else { -1.0 }
                } else if u.is_multiple_of(2) {
                    1.0
                } else {
                    -1.0
                };
                let dir = add3(scale3(bis, c_in), scale3(perp_n, sign * c_perp));
                pos[w] = Some(add3(pu, scale3(unit3(dir, fwd_dir), bl)));
                if !in_q[w] {
                    q.push_back(w);
                    in_q[w] = true;
                }
            }
            continue;
        }

        for (idx, &(w, bid)) in unplaced_nbrs.iter().enumerate() {
            let bl = ideal_bond_length(mol, u, w, mol.bonds[bid].order);
            let dir = if k == 1 {
                if is_sp {
                    fwd_dir
                } else {
                    // u_orth points TOWARD grandparent (syn/0 deg); PI points AWAY (trans/anti 180 deg).
                    let parent_is_sp2 = placed_nbrs
                        .first()
                        .is_some_and(|&(p, _)| atom_is_sp2(mol, p));
                    let dihedral = if is_sp2 && parent_is_sp2 {
                        // Conjugated sp2-sp2 prefers planar trans (180 deg)
                        PI
                    } else if parent_is_sp2 && !is_sp2 {
                        // Ar-CH2-R / Ar-CH(R)-R prefers ~90 deg perpendicular to aromatic plane
                        if u.is_multiple_of(2) {
                            PI * 0.5
                        } else {
                            -PI * 0.5
                        }
                    } else if is_sp2 {
                        // sp3-C(=O) slightly twisted off eclipsed
                        if u.is_multiple_of(2) {
                            PI * 0.72
                        } else {
                            -PI * 0.72
                        }
                    } else {
                        // sp3-sp3 alkyl chain prefers staggered trans/anti zigzag (180 deg)
                        PI
                    };
                    let rot_orth = add3(
                        scale3(u_orth, dihedral.cos()),
                        scale3(v_orth, dihedral.sin()),
                    );
                    add3(scale3(fwd_dir, c_t), scale3(rot_orth, s_t))
                }
            } else if k == 2 && placed_nbrs.len() == 1 && !is_sp2 {
                // Two remaining substituents on tetrahedral sp3 atom (e.g. isopropyl -CH(CH3)2 or Ar-CH(CH3)COOH):
                // Place at PI ± 60° (±120° from u_orth): exact 109.47° tetrahedral angle, maximal |sin|=0.866
                // clearance out of aromatic planes AND cleanly staggered away from alkyl grandparents!
                let phi = PI + if idx == 0 { PI / 3.0 } else { -PI / 3.0 };
                let rot_perp = add3(scale3(u_orth, phi.cos()), scale3(v_orth, phi.sin()));
                add3(scale3(fwd_dir, c_t), scale3(rot_perp, s_t))
            } else if k == 3 && placed_nbrs.len() == 1 && !is_sp2 {
                // Three remaining substituents (e.g. -C(CH3)3, -CF3): staggered anti against grandparent
                let phi = PI + 2.0 * PI * (idx as f64) / 3.0;
                let rot_perp = add3(scale3(u_orth, phi.cos()), scale3(v_orth, phi.sin()));
                add3(scale3(fwd_dir, c_t), scale3(rot_perp, s_t))
            } else if k == 2 && placed_nbrs.len() == 1 && is_sp2 {
                // Two substituents on trigonal planar sp2 atom (e.g. -C(=O)OH, -NO2):
                // Place at 0° and 180° (in plane, 120° apart)
                let phi = if idx == 0 { 0.0 } else { PI };
                let perp = add3(scale3(u_orth, phi.cos()), scale3(v_orth, phi.sin()));
                add3(scale3(fwd_dir, c_t), scale3(perp, s_t))
            } else if placed_nbrs.is_empty() && k >= 3 {
                // Seed center atom
                let theta = if k == 4 {
                    109.471_f64.to_radians()
                } else {
                    2.0 * PI / 3.0
                };
                if idx == 0 {
                    [1.0, 0.0, 0.0]
                } else {
                    let phi = 2.0 * PI * ((idx - 1) as f64) / ((k - 1) as f64);
                    let base = [
                        -theta.cos().abs(),
                        theta.sin() * phi.cos(),
                        theta.sin() * phi.sin(),
                    ];
                    unit3(base, [0.0, 1.0, 0.0])
                }
            } else {
                let phi = 2.0 * PI * (idx as f64) / (k as f64);
                let perp = add3(scale3(u_orth, phi.cos()), scale3(v_orth, phi.sin()));
                add3(scale3(fwd_dir, c_t), scale3(perp, s_t))
            };
            pos[w] = Some(add3(pu, scale3(unit3(dir, fwd_dir), bl)));
            if !in_q[w] {
                q.push_back(w);
                in_q[w] = true;
            }
        }
    }
}

fn inflate_bridged_cages(mol: &MolGraph, pts: &mut [[f64; 3]]) {
    let n = mol.atoms.len();
    if mol.rings.len() < 3 || n < 6 {
        return;
    }
    // Group ring atoms into connected ring systems via ring bonds
    let mut visited = vec![false; n];
    for start in 0..n {
        if visited[start] || !mol.atoms[start].in_ring {
            continue;
        }
        let mut comp = Vec::new();
        let mut q = VecDeque::new();
        visited[start] = true;
        q.push_back(start);
        while let Some(u) = q.pop_front() {
            comp.push(u);
            for &(v, bid) in &mol.adj[u] {
                if mol.bonds[bid].in_ring && !visited[v] {
                    visited[v] = true;
                    q.push_back(v);
                }
            }
        }
        let k = comp.len();
        if !(6..=24).contains(&k) {
            continue;
        }
        let mut in_comp = vec![false; n];
        for &a in &comp {
            in_comp[a] = true;
        }
        // Count internal ring edges and degree>=3 bridgeheads inside this ring system
        let mut internal_edges = 0usize;
        let mut bridgeheads = 0usize;
        let mut sp3_count = 0usize;
        let mut any_chiral = false;
        for &a in &comp {
            let deg_in = mol.adj[a].iter().filter(|&&(v, _)| in_comp[v]).count();
            internal_edges += deg_in;
            if deg_in >= 3 {
                bridgeheads += 1;
            }
            if !mol.atoms[a].aromatic && !atom_is_sp2(mol, a) {
                sp3_count += 1;
            }
            if matches!(
                mol.atoms[a].chiral,
                ChiralTag::Clockwise | ChiralTag::CounterClockwise
            ) {
                any_chiral = true;
            }
        }
        internal_edges /= 2;
        let cyclomatic = internal_edges + 1 - k;
        // Require a genuine 3D sp3 polycyclic cage (e.g. Adamantane cyclomatic=3, Cubane cyclomatic=5)
        if cyclomatic < 3 || bridgeheads < 2 || sp3_count * 10 < k * 7 || any_chiral {
            continue;
        }

        // All-pairs shortest-path matrix inside `comp`
        let mut idx_map = vec![usize::MAX; n];
        for (i, &a) in comp.iter().enumerate() {
            idx_map[a] = i;
        }
        let mut d2 = vec![vec![0.0_f64; k]; k];
        // Estimate average ring size to scale d=2 nicely (4-ring cubane vs 6-ring adamantane)
        let avg_ring_len = {
            let matching: Vec<usize> = mol
                .rings
                .iter()
                .filter(|r| r.iter().all(|&a| in_comp[a]))
                .map(Vec::len)
                .collect();
            if matching.is_empty() {
                6.0
            } else {
                matching.iter().sum::<usize>() as f64 / matching.len() as f64
            }
        };
        let d2_step = if avg_ring_len < 4.8 { 2.38 } else { 2.95 };
        for i in 0..k {
            let mut dist = vec![usize::MAX; k];
            dist[i] = 0;
            let mut bfs = VecDeque::new();
            bfs.push_back(i);
            while let Some(u) = bfs.pop_front() {
                for &(v_glob, _) in &mol.adj[comp[u]] {
                    let v = idx_map[v_glob];
                    if v < k && dist[v] == usize::MAX {
                        dist[v] = dist[u] + 1;
                        bfs.push_back(v);
                    }
                }
            }
            for j in 0..k {
                let sp = dist[j];
                d2[i][j] = match sp {
                    0 => 0.0,
                    1 => 1.54 * 1.54,
                    _ => d2_step * (sp as f64),
                };
            }
        }

        // Classical MDS double-centering: B = -0.5 * (D^2 - row_mean - col_mean + grand_mean)
        let mut row_mean = vec![0.0_f64; k];
        let mut grand = 0.0_f64;
        for i in 0..k {
            let s: f64 = d2[i].iter().sum();
            row_mean[i] = s / (k as f64);
            grand += s;
        }
        grand /= (k * k) as f64;
        let mut b = vec![vec![0.0_f64; k]; k];
        for i in 0..k {
            for j in 0..k {
                b[i][j] = -0.5 * (d2[i][j] - row_mean[i] - row_mean[j] + grand);
            }
        }

        // Exact symmetric Jacobi eigendecomposition (handles 3-fold degenerate Oh/Td cages cleanly)
        let mut v_mat = vec![vec![0.0_f64; k]; k];
        for (i, row) in v_mat.iter_mut().enumerate() {
            row[i] = 1.0;
        }
        for _sweep in 0..24 {
            let mut max_off = 0.0_f64;
            for p in 0..k {
                for q_idx in (p + 1)..k {
                    let apq = b[p][q_idx];
                    if apq.abs() > max_off {
                        max_off = apq.abs();
                    }
                    if apq.abs() < 1e-12 {
                        continue;
                    }
                    let app = b[p][p];
                    let aqq = b[q_idx][q_idx];
                    let tau = (aqq - app) / (2.0 * apq);
                    let t = if tau >= 0.0 {
                        1.0 / (tau + (1.0 + tau * tau).sqrt())
                    } else {
                        -1.0 / (-tau + (1.0 + tau * tau).sqrt())
                    };
                    let c = 1.0 / (1.0 + t * t).sqrt();
                    let s = t * c;
                    b[p][p] = app - t * apq;
                    b[q_idx][q_idx] = aqq + t * apq;
                    b[p][q_idx] = 0.0;
                    b[q_idx][p] = 0.0;
                    for r in 0..k {
                        if r != p && r != q_idx {
                            let arp = b[r][p];
                            let arq = b[r][q_idx];
                            b[r][p] = c * arp - s * arq;
                            b[p][r] = b[r][p];
                            b[r][q_idx] = s * arp + c * arq;
                            b[q_idx][r] = b[r][q_idx];
                        }
                        let vrp = v_mat[r][p];
                        let vrq = v_mat[r][q_idx];
                        v_mat[r][p] = c * vrp - s * vrq;
                        v_mat[r][q_idx] = s * vrp + c * vrq;
                    }
                }
            }
            if max_off < 1e-10 {
                break;
            }
        }

        let mut eig_order: Vec<usize> = (0..k).collect();
        eig_order.sort_by(|&a, &c| b[c][c].total_cmp(&b[a][a]));
        let e0 = eig_order[0];
        let e1 = eig_order[1];
        let e2 = eig_order[2];
        let s0 = b[e0][e0].max(0.0).sqrt();
        let s1 = b[e1][e1].max(0.0).sqrt();
        let s2 = b[e2][e2].max(0.0).sqrt();
        if s0 < 1e-3 || s1 < 1e-3 || s2 < 1e-3 {
            continue;
        }

        // Preserve original centroid
        let mut cent = [0.0; 3];
        for &a in &comp {
            cent = add3(cent, pts[a]);
        }
        cent = scale3(cent, 1.0 / k as f64);

        for (i, &a) in comp.iter().enumerate() {
            pts[a] = add3(
                cent,
                [s0 * v_mat[i][e0], s1 * v_mat[i][e1], s2 * v_mat[i][e2]],
            );
        }

        // Push any directly attached exocyclic substituents radially outward from the cage
        for &a in &comp {
            let out_dir = unit3(sub3(pts[a], cent), [1.0, 0.0, 0.0]);
            for &(exo, bid) in &mol.adj[a] {
                if in_comp[exo] {
                    continue;
                }
                let bl = ideal_bond_length(mol, a, exo, mol.bonds[bid].order);
                let new_exo = add3(pts[a], scale3(out_dir, bl));
                let delta = sub3(new_exo, pts[exo]);
                // Shift the connected substituent branch outside `comp`
                let mut sub_vis = vec![false; n];
                sub_vis[a] = true;
                sub_vis[exo] = true;
                let mut sq = VecDeque::new();
                sq.push_back(exo);
                while let Some(u) = sq.pop_front() {
                    pts[u] = add3(pts[u], delta);
                    for &(w, _) in &mol.adj[u] {
                        if !sub_vis[w] && !in_comp[w] {
                            sub_vis[w] = true;
                            sq.push_back(w);
                        }
                    }
                }
            }
        }
    }
}

enum ChiralRestraint {
    /// 3 heavy ligands `(i1, i2, i3)` around `ci`, with target sign `+1.0` or `-1.0` on `det(u1, u2, u3)`.
    ThreeHeavy {
        ci: usize,
        i1: usize,
        i2: usize,
        i3: usize,
        want_sign: f64,
    },
    /// 4 heavy ligands `(i0, i1, i2, i3)`, with target sign `+1.0` or `-1.0` on `det(p1-p0, p2-p0, p3-p0)`.
    FourHeavy {
        i0: usize,
        i1: usize,
        i2: usize,
        i3: usize,
        want_sign: f64,
    },
}

/// Coordinate-relaxation solver keeping `pinned[i] == true` atoms fixed while relaxing bond
/// lengths, angles, planarity, chiral volumes, E/Z geometry, and nonbonded clashes.
pub fn relax_3d_coords_pinned(mol: &MolGraph, pts: &mut [[f64; 3]], pinned: &[bool], iters: usize) {
    let n = mol.atoms.len();
    if n < 2 {
        return;
    }

    // Build 1-2 bonded set and 1-3 angle list
    let mut is_12_13 = vec![vec![false; n]; n];
    let mut bonds_ideal = Vec::with_capacity(mol.bonds.len());
    for b in &mol.bonds {
        is_12_13[b.u][b.v] = true;
        is_12_13[b.v][b.u] = true;
        bonds_ideal.push((b.u, b.v, ideal_bond_length(mol, b.u, b.v, b.order)));
    }

    let mut angles_13 = Vec::new();
    let mut sp2_triples = Vec::new();
    for c in 0..n {
        let nbrs = &mol.adj[c];
        for i in 0..nbrs.len() {
            let (a, ba) = nbrs[i];
            let d_ca = ideal_bond_length(mol, c, a, mol.bonds[ba].order);
            for &(b, bb) in nbrs.iter().skip(i + 1) {
                is_12_13[a][b] = true;
                is_12_13[b][a] = true;
                let d_cb = ideal_bond_length(mol, c, b, mol.bonds[bb].order);
                let theta = ideal_angle_rad(mol, c, a, b);
                let target_13 = (d_ca * d_ca + d_cb * d_cb - 2.0 * d_ca * d_cb * theta.cos())
                    .max(0.5)
                    .sqrt();
                angles_13.push((a, b, target_13));
            }
        }
        if atom_is_sp2(mol, c) && nbrs.len() == 3 {
            sp2_triples.push((c, nbrs[0].0, nbrs[1].0, nbrs[2].0));
        }
    }

    // Collect planar aromatic / fully-sp2 rings + their directly attached exocyclic atoms
    // so aromatic rings (benzene, pyridine, indole, naphthalene, caffeine) stay razor-flat!
    let mut planar_ring_groups: Vec<(Vec<usize>, Vec<usize>)> = Vec::new();
    let planar_rings: Vec<&Vec<usize>> = mol
        .rings
        .iter()
        .filter(|r| {
            r.len() <= 7
                && r.iter()
                    .all(|&a| mol.atoms[a].aromatic || atom_is_sp2(mol, a))
        })
        .collect();
    if !planar_rings.is_empty() {
        // Merge fused planar rings into shared flat planes (e.g. naphthalene, indole, purine)
        let mut ring_used = vec![false; planar_rings.len()];
        for i in 0..planar_rings.len() {
            if ring_used[i] {
                continue;
            }
            ring_used[i] = true;
            let mut ring_atoms: Vec<usize> = planar_rings[i].clone();
            loop {
                let mut added = false;
                for j in 0..planar_rings.len() {
                    if ring_used[j] {
                        continue;
                    }
                    let shared = planar_rings[j]
                        .iter()
                        .filter(|a| ring_atoms.contains(a))
                        .count();
                    if shared >= 2 {
                        ring_used[j] = true;
                        for &a in planar_rings[j] {
                            if !ring_atoms.contains(&a) {
                                ring_atoms.push(a);
                            }
                        }
                        added = true;
                    }
                }
                if !added {
                    break;
                }
            }
            // Directly attached neighbors of aromatic/sp2 ring atoms also lie in the ring plane
            let mut exo_atoms: Vec<usize> = Vec::new();
            for &ra in &ring_atoms {
                for &(nb, _) in &mol.adj[ra] {
                    if !ring_atoms.contains(&nb) && !exo_atoms.contains(&nb) {
                        exo_atoms.push(nb);
                    }
                }
            }
            // Keep first ring's ordered vertices to compute a consistent Newell normal
            planar_ring_groups.push((ring_atoms, exo_atoms));
        }
    }

    // Collect tetrahedral chiral volume restraints directly from OpenSMILES `chiral_seq`
    let mut chiral_restraints = Vec::new();
    for (ci, atom) in mol.atoms.iter().enumerate() {
        if !matches!(
            atom.chiral,
            ChiralTag::Clockwise | ChiralTag::CounterClockwise
        ) || atom.chiral_seq.len() != 4
        {
            continue;
        }
        let base_sign = if atom.chiral == ChiralTag::Clockwise {
            1.0
        } else {
            -1.0
        };
        let h_pos = atom.chiral_seq.iter().position(|&x| x == usize::MAX);
        match h_pos {
            Some(hp) => {
                let heavy: Vec<usize> = atom
                    .chiral_seq
                    .iter()
                    .copied()
                    .filter(|&x| x != usize::MAX && x < n)
                    .collect();
                if heavy.len() == 3 {
                    let parity = if hp.is_multiple_of(2) { 1.0 } else { -1.0 };
                    chiral_restraints.push(ChiralRestraint::ThreeHeavy {
                        ci,
                        i1: heavy[0],
                        i2: heavy[1],
                        i3: heavy[2],
                        want_sign: base_sign * parity,
                    });
                }
            }
            None => {
                let s = &atom.chiral_seq;
                if s.iter().all(|&x| x < n) {
                    chiral_restraints.push(ChiralRestraint::FourHeavy {
                        i0: s[0],
                        i1: s[1],
                        i2: s[2],
                        i3: s[3],
                        want_sign: base_sign,
                    });
                }
            }
        }
    }

    // Collect acyclic E/Z double-bond 1-4 restraints `(a, u, v, b, target_14)`
    let mut ez_restraints = Vec::new();
    for b in &mol.bonds {
        if b.order == BondOrder::Double
            && !b.in_ring
            && matches!(b.cis_trans, BondCisTrans::E | BondCisTrans::Z)
            && let (Some(a_hi), Some(b_hi)) = (
                highest_cip_end_nbr(mol, b.u, b.v),
                highest_cip_end_nbr(mol, b.v, b.u),
            )
        {
            let target_14 = if b.cis_trans == BondCisTrans::Z {
                2.82
            } else {
                3.82
            };
            ez_restraints.push((a_hi, b.u, b.v, b_hi, target_14));
        }
    }

    let mut disp = vec![[0.0_f64; 3]; n];
    for _ in 0..iters {
        disp.fill([0.0, 0.0, 0.0]);

        // 1. Bond stretch
        for &(u, v, target) in &bonds_ideal {
            let d = sub3(pts[v], pts[u]);
            let len = norm3(d);
            if len > 1e-6 {
                let f = 0.42 * (len - target) / len;
                let df = scale3(d, f);
                disp[u] = add3(disp[u], df);
                disp[v] = sub3(disp[v], df);
            }
        }

        // 2. 1-3 angle spring
        for &(a, b, target) in &angles_13 {
            let d = sub3(pts[b], pts[a]);
            let len = norm3(d);
            if len > 1e-6 {
                let f = 0.24 * (len - target) / len;
                let df = scale3(d, f);
                disp[a] = add3(disp[a], df);
                disp[b] = sub3(disp[b], df);
            }
        }

        // 3. Out-of-plane planarity for degree-3 sp2 centers (symmetric distribution)
        for &(c, a, b, d_atm) in &sp2_triples {
            let plane_n = cross3(sub3(pts[b], pts[a]), sub3(pts[d_atm], pts[a]));
            let n_len = norm3(plane_n);
            if n_len > 1e-4 {
                let u_n = scale3(plane_n, 1.0 / n_len);
                let h = dot3(sub3(pts[c], pts[a]), u_n);
                let corr_c = scale3(u_n, -0.30 * h);
                let corr_n = scale3(u_n, 0.10 * h);
                disp[c] = add3(disp[c], corr_c);
                disp[a] = add3(disp[a], corr_n);
                disp[b] = add3(disp[b], corr_n);
                disp[d_atm] = add3(disp[d_atm], corr_n);
            }
        }

        // 3b. Exact aromatic & fused planar ring flattening (momentum-conserving around ring centroid;
        // exocyclic attachments are kept coplanar with zero net force by `sp2_triples` above).
        for (ring_atoms, _exo_atoms) in &planar_ring_groups {
            if ring_atoms.len() < 3 {
                continue;
            }
            let mut cent = [0.0; 3];
            for &a in ring_atoms {
                cent = add3(cent, pts[a]);
            }
            cent = scale3(cent, 1.0 / ring_atoms.len() as f64);
            let mut normal = [0.0; 3];
            for &u in ring_atoms {
                let ru = sub3(pts[u], cent);
                for &(v, bid) in &mol.adj[u] {
                    if mol.bonds[bid].in_ring && ring_atoms.contains(&v) {
                        let rv = sub3(pts[v], cent);
                        let cr = cross3(ru, rv);
                        if dot3(cr, normal) >= 0.0 || norm3(normal) < 1e-9 {
                            normal = add3(normal, cr);
                        } else {
                            normal = sub3(normal, cr);
                        }
                    }
                }
            }
            let n_len = norm3(normal);
            if n_len > 1e-5 {
                let u_n = scale3(normal, 1.0 / n_len);
                for &a in ring_atoms {
                    let h = dot3(sub3(pts[a], cent), u_n);
                    disp[a] = add3(disp[a], scale3(u_n, -0.42 * h));
                }
            }
        }

        // 4. Tetrahedral signed chiral-volume restoring force (works on acyclic & ring bridgeheads!)
        for cr in &chiral_restraints {
            match *cr {
                ChiralRestraint::ThreeHeavy {
                    ci,
                    i1,
                    i2,
                    i3,
                    want_sign,
                } => {
                    let u1 = sub3(pts[i1], pts[ci]);
                    let u2 = sub3(pts[i2], pts[ci]);
                    let u3 = sub3(pts[i3], pts[ci]);
                    let vol = dot3(cross3(u1, u2), u3);
                    let target_v = 1.95;
                    if want_sign * vol < target_v {
                        let k = (0.24 * (target_v - want_sign * vol)).min(0.52) * want_sign;
                        let g1 = cross3(u2, u3);
                        let g2 = cross3(u3, u1);
                        let g3 = cross3(u1, u2);
                        let gc = scale3(add3(add3(g1, g2), g3), -1.0);
                        disp[i1] = add3(disp[i1], scale3(g1, k));
                        disp[i2] = add3(disp[i2], scale3(g2, k));
                        disp[i3] = add3(disp[i3], scale3(g3, k));
                        disp[ci] = add3(disp[ci], scale3(gc, k));
                    }
                }
                ChiralRestraint::FourHeavy {
                    i0,
                    i1,
                    i2,
                    i3,
                    want_sign,
                } => {
                    let u1 = sub3(pts[i1], pts[i0]);
                    let u2 = sub3(pts[i2], pts[i0]);
                    let u3 = sub3(pts[i3], pts[i0]);
                    let vol = dot3(cross3(u1, u2), u3);
                    let target_v = 5.2;
                    if want_sign * vol < target_v {
                        let k = (0.065 * (target_v - want_sign * vol)).min(0.42) * want_sign;
                        let g1 = cross3(u2, u3);
                        let g2 = cross3(u3, u1);
                        let g3 = cross3(u1, u2);
                        let g0 = scale3(add3(add3(g1, g2), g3), -1.0);
                        disp[i1] = add3(disp[i1], scale3(g1, k));
                        disp[i2] = add3(disp[i2], scale3(g2, k));
                        disp[i3] = add3(disp[i3], scale3(g3, k));
                        disp[i0] = add3(disp[i0], scale3(g0, k));
                    }
                }
            }
        }

        // 5. Acyclic E/Z double-bond cis/trans 1-4 spring & coplanarity
        for &(a, u, v, b, target_14) in &ez_restraints {
            let d = sub3(pts[b], pts[a]);
            let len = norm3(d);
            if len > 1e-5 {
                let f = 0.22 * (len - target_14) / len;
                let df = scale3(d, f);
                disp[a] = add3(disp[a], df);
                disp[b] = sub3(disp[b], df);
            }
            // Pull `b` onto plane of `(a, u, v)`
            let pn = cross3(sub3(pts[u], pts[a]), sub3(pts[v], pts[u]));
            let pnl = norm3(pn);
            if pnl > 1e-4 {
                let un = scale3(pn, 1.0 / pnl);
                let h = dot3(sub3(pts[b], pts[v]), un);
                disp[b] = add3(disp[b], scale3(un, -0.22 * h));
            }
        }

        // 6. Nonbonded VdW soft repulsion (avoid self-clashes)
        if n <= 400 {
            for i in 0..n {
                for j in (i + 1)..n {
                    if is_12_13[i][j] {
                        continue;
                    }
                    let d = sub3(pts[j], pts[i]);
                    let d2 = dot3(d, d);
                    let r_min = 2.45;
                    if d2 > 1e-6 && d2 < r_min * r_min {
                        let dist = d2.sqrt();
                        let push = 0.22 * (r_min - dist) / dist;
                        let df = scale3(d, push);
                        disp[i] = sub3(disp[i], df);
                        disp[j] = add3(disp[j], df);
                    }
                }
            }
        }

        for i in 0..n {
            if pinned.get(i).copied().unwrap_or(false) {
                continue;
            }
            let mag = norm3(disp[i]);
            let step = if mag > 0.28 {
                scale3(disp[i], 0.28 / mag)
            } else {
                disp[i]
            };
            if step.iter().all(|x| x.is_finite()) {
                pts[i] = add3(pts[i], step);
            }
        }
    }
}

/// Adjust tetrahedral centers whose signed 3D volume disagrees with OpenSMILES `@` / `@@`.
fn enforce_chiral_signs(mol: &MolGraph, pts: &mut [[f64; 3]]) {
    for (ci, atom) in mol.atoms.iter().enumerate() {
        if !matches!(
            atom.chiral,
            ChiralTag::Clockwise | ChiralTag::CounterClockwise
        ) {
            continue;
        }
        if atom.chiral_seq.len() != 4 {
            continue;
        }
        let real_nbrs: Vec<usize> = atom
            .chiral_seq
            .iter()
            .copied()
            .filter(|&x| x != usize::MAX && x < mol.atoms.len())
            .collect();
        if real_nbrs.len() < 3 {
            continue;
        }
        let pc = pts[ci];
        let p_vec = |id: usize| -> [f64; 3] {
            if id == usize::MAX || id >= mol.atoms.len() {
                // Suppressed/implicit H points opposite sum of the 3 heavy ligands
                let mut sum = [0.0, 0.0, 0.0];
                for &rn in &real_nbrs {
                    sum = add3(sum, unit3(sub3(pts[rn], pc), [1.0, 0.0, 0.0]));
                }
                scale3(unit3(sum, [0.0, 0.0, 1.0]), -1.09)
            } else {
                sub3(pts[id], pc)
            }
        };
        let v0 = p_vec(atom.chiral_seq[0]);
        let v1 = p_vec(atom.chiral_seq[1]);
        let v2 = p_vec(atom.chiral_seq[2]);
        let v3 = p_vec(atom.chiral_seq[3]);
        // Signed volume viewing from v0 down to center: ((v1-v0) x (v2-v0)) . (v3-v0)
        // Right-hand rule: looking from v0 toward center (along -v0), clockwise (v1->v2->v3) gives POSITIVE det!
        let det = dot3(cross3(sub3(v1, v0), sub3(v2, v0)), sub3(v3, v0));
        let want_cw = atom.chiral == ChiralTag::Clockwise;
        let current_cw = det > 0.0;
        if current_cw != want_cw {
            let leaves: Vec<usize> = real_nbrs
                .iter()
                .copied()
                .filter(|&nb| mol.adj[nb].len() == 1)
                .collect();
            if leaves.len() >= 2 {
                let a = leaves[0];
                let b = leaves[1];
                let la = norm3(sub3(pts[a], pc)).max(1.0);
                let lb = norm3(sub3(pts[b], pc)).max(1.0);
                let da = unit3(sub3(pts[a], pc), [1.0, 0.0, 0.0]);
                let db = unit3(sub3(pts[b], pc), [0.0, 1.0, 0.0]);
                pts[a] = add3(pc, scale3(db, la));
                pts[b] = add3(pc, scale3(da, lb));
            } else if let Some(&leaf) = leaves.first() {
                // Reflect single leaf across plane formed by center and two other neighbors
                let others: Vec<usize> = real_nbrs.iter().copied().filter(|&x| x != leaf).collect();
                if others.len() >= 2 {
                    let n = unit3(
                        cross3(sub3(pts[others[0]], pc), sub3(pts[others[1]], pc)),
                        [0.0, 0.0, 1.0],
                    );
                    let dl = sub3(pts[leaf], pc);
                    let proj = dot3(dl, n);
                    pts[leaf] = add3(pc, sub3(dl, scale3(n, 2.0 * proj)));
                }
            } else if real_nbrs.len() == 3 {
                // Ring stereocenter with implicit H: invert central atom `ci` across plane of its 3 heavy neighbors!
                let p0 = pts[real_nbrs[0]];
                let p1 = pts[real_nbrs[1]];
                let p2 = pts[real_nbrs[2]];
                let n = unit3(cross3(sub3(p1, p0), sub3(p2, p0)), [0.0, 0.0, 1.0]);
                let h = dot3(sub3(pc, p0), n);
                let shift = if h.abs() < 0.25 {
                    if h >= 0.0 { -0.75 } else { 0.75 }
                } else {
                    -2.2 * h
                };
                pts[ci] = add3(pc, scale3(n, shift));
            } else if real_nbrs.len() == 4 {
                // Quaternary ring stereocenter: reflect first exocyclic substituent (if any) across other 3
                if let Some(&exo) = real_nbrs.iter().find(|&&nb| !mol.atoms[nb].in_ring) {
                    let others: Vec<usize> =
                        real_nbrs.iter().copied().filter(|&x| x != exo).collect();
                    if others.len() >= 2 {
                        let n = unit3(
                            cross3(sub3(pts[others[0]], pc), sub3(pts[others[1]], pc)),
                            [0.0, 0.0, 1.0],
                        );
                        let dl = sub3(pts[exo], pc);
                        let proj = dot3(dl, n);
                        pts[exo] = add3(pc, sub3(dl, scale3(n, 2.0 * proj)));
                    }
                }
            }
        }
    }
}

/// Write a `MolGraph` as a standard PDB file (with `ATOM` for peptides and `HETATM` + `CONECT` for ligands).
pub fn write_pdb(mol: &MolGraph) -> String {
    let mut out = String::new();
    let title = mol.name.as_deref().unwrap_or("CHEMCRAFT 3D STRUCTURE");
    let _ = writeln!(out, "HEADER    CHEMCRAFT MOLECULE");
    let _ = writeln!(out, "TITLE     {}", title.to_ascii_uppercase());
    if let Some(seq) = mol.properties.get("PEPTIDE_SEQUENCE") {
        let _ = writeln!(out, "REMARK   1 PEPTIDE SEQUENCE: {seq}");
    }
    if let Some(ss) = mol.properties.get("SECONDARY_STRUCTURE") {
        let _ = writeln!(out, "REMARK   2 SECONDARY STRUCTURE: {ss}");
    }
    if let Some(helm) = mol.properties.get("PEPTIDE_HELM") {
        let _ = writeln!(out, "REMARK   3 HELM: {helm}");
    }

    let coords = compute_3d_coords(mol);

    // Emit standard PDB SSBOND records for covalent S-S disulfide bridges between residues.
    let mut ss_serial = 0usize;
    for b in &mol.bonds {
        let au = &mol.atoms[b.u];
        let av = &mol.atoms[b.v];
        if au.symbol == "S"
            && av.symbol == "S"
            && let (Some(mu), Some(mv)) = (&au.pdb_meta, &av.pdb_meta)
        {
            ss_serial += 1;
            let pu = coords[b.u];
            let pv = coords[b.v];
            let dist =
                ((pu[0] - pv[0]).powi(2) + (pu[1] - pv[1]).powi(2) + (pu[2] - pv[2]).powi(2))
                    .sqrt();
            let _ = writeln!(
                out,
                "SSBOND {:>3} {:>3.3} {} {:>4}    {:>3.3} {} {:>4}                       1555   1555 {:>5.2}",
                ss_serial,
                mu.res_name,
                mu.chain,
                mu.res_seq,
                mv.res_name,
                mv.chain,
                mv.res_seq,
                dist
            );
        }
    }

    let mut elem_counts: HashMap<&str, usize> = HashMap::new();
    let mut any_polymer = false;

    for (idx, a) in mol.atoms.iter().enumerate() {
        let [x, y, z] = coords[idx];
        let c = elem_counts.entry(a.symbol.as_str()).or_insert(0);
        *c += 1;
        let serial = idx + 1;
        let (rec_kind, raw_aname, res_name, chain, res_seq) = if let Some(meta) = &a.pdb_meta {
            if !meta.is_hetatm {
                any_polymer = true;
            }
            (
                if meta.is_hetatm { "HETATM" } else { "ATOM  " },
                meta.atom_name.clone(),
                meta.res_name.clone(),
                meta.chain,
                meta.res_seq,
            )
        } else {
            let fallback_name = format!("{}{}", a.symbol, *c);
            ("HETATM", fallback_name, "UNL".to_string(), 'A', 1)
        };
        let aname_fmt = if raw_aname.len() >= 4 {
            format!("{raw_aname:<4.4}")
        } else {
            format!(" {raw_aname:<3}")
        };
        let _ = writeln!(
            out,
            "{rec_kind}{serial:>5} {aname_fmt} {res_name:>3.3} {chain}{res_seq:>4}    {x:>8.3}{y:>8.3}{z:>8.3}  1.00  0.00          {:>2}",
            a.symbol.to_ascii_uppercase()
        );
    }

    if any_polymer {
        let _ = writeln!(out, "TER");
    }

    // Write CONECT records for every atom's bonded neighbors (chunks of 4 per line per PDB spec)
    for (u, nbrs) in mol.adj.iter().enumerate() {
        if nbrs.is_empty() {
            continue;
        }
        let mut sorted_nbrs: Vec<usize> = nbrs.iter().map(|&(v, _)| v + 1).collect();
        sorted_nbrs.sort_unstable();
        for chunk in sorted_nbrs.chunks(4) {
            let _ = write!(out, "CONECT{:>5}", u + 1);
            for v1 in chunk {
                let _ = write!(out, "{v1:>5}");
            }
            let _ = writeln!(out);
        }
    }
    let _ = writeln!(out, "END");
    out
}

/// Compute explicit 3D hydrogens for aromatic `[nH]` atoms (pyrrole, indole, uracil, sildenafil,
/// tetrazole) so external 3D consumers (DockCraft, ChimeraX, RDKit) preserve exact donor valences.
fn aromatic_nh_hydrogens_3d(mol: &MolGraph, coords: &[[f64; 3]]) -> Vec<(usize, [f64; 3])> {
    let mut out = Vec::new();
    for (u, a) in mol.atoms.iter().enumerate() {
        if a.symbol == "N" && a.aromatic && a.formal_charge == 0 && a.total_h() == 1 {
            let pu = coords[u];
            let mut sum = [0.0, 0.0, 0.0];
            for &(nb, _) in &mol.adj[u] {
                sum = add3(sum, unit3(sub3(coords[nb], pu), [1.0, 0.0, 0.0]));
            }
            let dir = unit3(scale3(sum, -1.0), [0.0, 1.0, 0.0]);
            out.push((u, add3(pu, scale3(dir, 1.01))));
        }
    }
    out
}

fn write_sdf_v3000_3d_record(
    out: &mut String,
    mol: &MolGraph,
    coords: &[[f64; 3]],
    extra_hs: &[(usize, [f64; 3])],
) {
    let name = mol.name.as_deref().unwrap_or("ChemCraft-3D");
    let _ = writeln!(out, "{name}");
    // CTAB spec: dimension code ("2D"/"3D") occupies columns 21-22 (0-indexed 20..22)
    let _ = writeln!(out, "  ChemCraft 111926  3D");
    let _ = writeln!(out);
    let _ = writeln!(out, "  0  0  0     0  0            999 V3000");
    let _ = writeln!(out, "M  V30 BEGIN CTAB");
    let _ = writeln!(
        out,
        "M  V30 COUNTS {} {} 0 0 0",
        mol.atoms.len() + extra_hs.len(),
        mol.bonds.len() + extra_hs.len()
    );
    let _ = writeln!(out, "M  V30 BEGIN ATOM");
    for (idx, a) in mol.atoms.iter().enumerate() {
        let [x, y, z] = coords[idx];
        if a.formal_charge != 0 {
            let _ = writeln!(
                out,
                "M  V30 {} {} {:.4} {:.4} {:.4} 0 CHG={}",
                a.id + 1,
                a.symbol,
                x,
                y,
                z,
                a.formal_charge
            );
        } else {
            let _ = writeln!(
                out,
                "M  V30 {} {} {:.4} {:.4} {:.4} 0",
                a.id + 1,
                a.symbol,
                x,
                y,
                z
            );
        }
    }
    for (h_i, &(_, [hx, hy, hz])) in extra_hs.iter().enumerate() {
        let _ = writeln!(
            out,
            "M  V30 {} H {:.4} {:.4} {:.4} 0",
            mol.atoms.len() + h_i + 1,
            hx,
            hy,
            hz
        );
    }
    let _ = writeln!(out, "M  V30 END ATOM");
    let _ = writeln!(out, "M  V30 BEGIN BOND");
    for (idx, b) in mol.bonds.iter().enumerate() {
        let t = match b.order {
            BondOrder::Single => 1,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            BondOrder::Aromatic => 4,
            BondOrder::Dative => 9,
        };
        let cfg = match b.stereo {
            BondStereo::Up => " CFG=1",
            BondStereo::Down => " CFG=3",
            BondStereo::Wavy => " CFG=2",
            _ => "",
        };
        let _ = writeln!(
            out,
            "M  V30 {} {} {} {}{}",
            idx + 1,
            t,
            b.u + 1,
            b.v + 1,
            cfg
        );
    }
    for (h_i, &(parent, _)) in extra_hs.iter().enumerate() {
        let _ = writeln!(
            out,
            "M  V30 {} 1 {} {}",
            mol.bonds.len() + h_i + 1,
            parent + 1,
            mol.atoms.len() + h_i + 1
        );
    }
    let _ = writeln!(out, "M  V30 END BOND");
    let _ = writeln!(out, "M  V30 END CTAB");
    let _ = writeln!(out, "M  END");
    for (k, v) in &mol.properties {
        let _ = writeln!(out, "> <{k}>\n{v}\n");
    }
    let _ = writeln!(out, "$$$$");
}

/// Write `MolGraph` list as a 3D SDF block (preserving/generating 3D `[x, y, z]` coordinates,
/// explicit aromatic `[nH]` hydrogens, and automatically switching to CTAB V3000 when >999 atoms/bonds).
pub fn write_sdf_v2000_3d(mols: &[&MolGraph]) -> String {
    let mut out = String::new();
    for mol in mols {
        let coords = compute_3d_coords(mol);
        let extra_hs = aromatic_nh_hydrogens_3d(mol, &coords);
        let total_atoms = mol.atoms.len() + extra_hs.len();
        let total_bonds = mol.bonds.len() + extra_hs.len();
        if total_atoms > 999 || total_bonds > 999 {
            write_sdf_v3000_3d_record(&mut out, mol, &coords, &extra_hs);
            continue;
        }
        let name = mol.name.as_deref().unwrap_or("ChemCraft-3D");
        let _ = writeln!(out, "{name}");
        // CTAB spec: dimension code ("2D"/"3D") occupies columns 21-22 (0-indexed 20..22)
        let _ = writeln!(out, "  ChemCraft 111926  3D");
        let _ = writeln!(out);
        let _ = writeln!(
            out,
            "{:>3}{:>3}  0  0  0  0  0  0  0  0999 V2000",
            total_atoms, total_bonds
        );
        for (idx, a) in mol.atoms.iter().enumerate() {
            let [x, y, z] = coords[idx];
            let _ = writeln!(
                out,
                "{:>10.4}{:>10.4}{:>10.4} {:<3} 0  0  0  0  0  0  0  0  0  0  0  0",
                x, y, z, a.symbol
            );
        }
        for &(_, [hx, hy, hz]) in &extra_hs {
            let _ = writeln!(
                out,
                "{:>10.4}{:>10.4}{:>10.4} {:<3} 0  0  0  0  0  0  0  0  0  0  0  0",
                hx, hy, hz, "H"
            );
        }
        for b in &mol.bonds {
            let t = match b.order {
                BondOrder::Single => 1,
                BondOrder::Double => 2,
                BondOrder::Triple => 3,
                BondOrder::Aromatic => 4,
                BondOrder::Dative => 9,
            };
            let s = match b.stereo {
                BondStereo::Up => 1,
                BondStereo::Down => 6,
                BondStereo::Wavy => 4,
                _ => 0,
            };
            let _ = writeln!(out, "{:>3}{:>3}{:>3}{:>3}  0  0  0", b.u + 1, b.v + 1, t, s);
        }
        for (h_i, &(parent, _)) in extra_hs.iter().enumerate() {
            let h_id = mol.atoms.len() + h_i + 1;
            let _ = writeln!(out, "{:>3}{:>3}{:>3}{:>3}  0  0  0", parent + 1, h_id, 1, 0);
        }
        let charged: Vec<(usize, i8)> = mol
            .atoms
            .iter()
            .filter(|a| a.formal_charge != 0)
            .map(|a| (a.id + 1, a.formal_charge))
            .collect();
        for chunk in charged.chunks(8) {
            let _ = write!(out, "M  CHG{:>3}", chunk.len());
            for (aid, chg) in chunk {
                let _ = write!(out, "{aid:>4}{chg:>4}");
            }
            let _ = writeln!(out);
        }
        let _ = writeln!(out, "M  END");
        for (k, v) in &mol.properties {
            let _ = writeln!(out, "> <{k}>\n{v}\n");
        }
        let _ = writeln!(out, "$$$$");
    }
    out
}

/// Quick NeRF atom placement helper shared with peptide and macrocycle builders.
pub fn place_nerf_atom(
    p_a: [f64; 3],
    p_b: [f64; 3],
    p_c: [f64; 3],
    bond_len: f64,
    bond_angle_rad: f64,
    torsion_rad: f64,
) -> [f64; 3] {
    let bc = unit3(sub3(p_c, p_b), [1.0, 0.0, 0.0]);
    let ab = sub3(p_b, p_a);
    let n = unit3(cross3(ab, bc), ortho_pair(bc).0);
    let m = unit3(cross3(n, bc), ortho_pair(bc).1);
    let theta = PI - bond_angle_rad;
    let d_local = [
        bond_len * theta.cos(),
        bond_len * theta.sin() * torsion_rad.cos(),
        bond_len * theta.sin() * torsion_rad.sin(),
    ];
    add3(
        p_c,
        add3(
            scale3(bc, d_local[0]),
            add3(scale3(m, d_local[1]), scale3(n, d_local[2])),
        ),
    )
}

#[allow(dead_code)]
fn _keep_rodrigues_used() -> [f64; 3] {
    rodrigues_rotate([1.0, 0.0, 0.0], [0.0, 0.0, 1.0], 0.0)
}
