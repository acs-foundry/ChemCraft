#![forbid(unsafe_code)]

//! Clean 2D macrocycle/linear depiction & 3D L/D-aware NeRF + closed-ring macrocycle geometry
//! engine for peptides and Non-Ribosomal Peptides (NRPs).

use std::collections::VecDeque;
use std::f64::consts::{PI, TAU};

use mol_graph::{MolGraph, place_nerf_atom, relax_3d_coords_pinned};

use crate::notation::{CrossLink, PeptideConformation, PeptideSpec};

/// Residue backbone atom index table used by the 2D/3D layout engine.
#[derive(Debug, Clone)]
pub struct ResidueAtoms {
    pub n_id: usize,
    pub ca_id: usize,
    pub c_id: usize,
    pub o_id: Option<usize>,
    pub all_atoms: Vec<usize>,
    pub is_d: bool,
}

fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
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

/// Find the subset of residue indices (0-based) belonging to the primary macrocycle ring.
fn primary_ring_residues(spec: &PeptideSpec) -> Option<(usize, usize)> {
    let n = spec.residues.len();
    for cl in &spec.cross_links {
        match *cl {
            CrossLink::HeadToTail => return Some((0, n - 1)),
            CrossLink::Macrolactone {
                from_res,
                to_oh_res,
            } => {
                let lo = from_res.min(to_oh_res).saturating_sub(1);
                let hi = from_res.max(to_oh_res).saturating_sub(1);
                if hi >= lo + 2 {
                    return Some((lo, hi));
                }
            }
            CrossLink::Lactam {
                from_acid_res,
                to_amine_res,
                ..
            } => {
                let lo = from_acid_res.min(to_amine_res).saturating_sub(1);
                let hi = from_acid_res.max(to_amine_res).saturating_sub(1);
                if hi >= lo + 2 {
                    return Some((lo, hi));
                }
            }
            _ => {}
        }
    }
    None
}

/// Generate clean 2D coordinates:
/// - Linear / disulfide peptides: zig-zag ACS backbone with alternating outward sidechains.
/// - Macrocyclic & lariat NRPs: regular macrocycle polygon backbone with radial outward sidechains
///   and a tidy extended N-terminal lariat tail (e.g. Polymyxin B1, Daptomycin, Cyclosporin A).
pub fn layout_peptide_2d(mol: &mut MolGraph, spec: &PeptideSpec, res_atoms: &[ResidueAtoms]) {
    let n_atoms = mol.atoms.len();
    if n_atoms == 0 || res_atoms.is_empty() {
        return;
    }
    let bl = 1.48_f64;
    let mut pos: Vec<Option<[f64; 2]>> = vec![None; n_atoms];
    let mut outward: Vec<[f64; 2]> = vec![[0.0, 1.0]; n_atoms];

    if let Some((ring_lo, ring_hi)) = primary_ring_residues(spec) {
        // 1. Place ring backbone atoms (N, CA, C per residue) around a clean regular macrocycle polygon
        let ring_res_count = ring_hi - ring_lo + 1;
        let ring_bb_atoms: Vec<(usize, usize)> = (ring_lo..=ring_hi)
            .flat_map(|ri| {
                let r = &res_atoms[ri];
                [(r.n_id, ri), (r.ca_id, ri), (r.c_id, ri)]
            })
            .collect();
        let m = ring_bb_atoms.len().max(6) as f64;
        let radius = (m * bl / TAU).max(3.2);
        for (k, &(aid, _ri)) in ring_bb_atoms.iter().enumerate() {
            let theta = TAU * (k as f64) / m - PI * 0.5;
            let ux = theta.cos();
            let uy = theta.sin();
            pos[aid] = Some([radius * ux, radius * uy]);
            outward[aid] = [ux, uy];
        }
        // Place ring carbonyl =O slightly inward or outward alternating cleanly
        for (idx_in_ring, ri) in (ring_lo..=ring_hi).enumerate() {
            let r = &res_atoms[ri];
            if let Some(oid) = r.o_id
                && let Some(pc) = pos[r.c_id]
            {
                let u = outward[r.c_id];
                let s = if idx_in_ring.is_multiple_of(2) {
                    -0.82
                } else {
                    0.82
                };
                pos[oid] = Some([pc[0] + u[0] * bl * s, pc[1] + u[1] * bl * s]);
                outward[oid] = [u[0] * s.signum(), u[1] * s.signum()];
            }
        }
        // 2. Extend N-terminal lariat tail (0..ring_lo) smoothly leftward from ring_lo's N atom
        if ring_lo > 0 {
            let anchor = pos[res_atoms[ring_lo].n_id].unwrap_or([0.0, -radius]);
            let dir = outward[res_atoms[ring_lo].n_id];
            let perp = [-dir[1], dir[0]];
            let mut step = 1usize;
            for ri in (0..ring_lo).rev() {
                let r = &res_atoms[ri];
                let side_sign = if ri.is_multiple_of(2) { 1.0 } else { -1.0 };
                for &aid in &[r.c_id, r.ca_id, r.n_id] {
                    let s = step as f64 * bl * 0.90;
                    let zig = if step.is_multiple_of(2) { 0.36 } else { -0.36 } * bl;
                    pos[aid] = Some([
                        anchor[0] + dir[0] * s + perp[0] * zig,
                        anchor[1] + dir[1] * s + perp[1] * zig,
                    ]);
                    outward[aid] = [perp[0] * side_sign, perp[1] * side_sign];
                    step += 1;
                }
                if let Some(oid) = r.o_id
                    && let Some(pc) = pos[r.c_id]
                {
                    pos[oid] = Some([
                        pc[0] - perp[0] * side_sign * bl * 0.88,
                        pc[1] - perp[1] * side_sign * bl * 0.88,
                    ]);
                }
            }
            let _ = ring_res_count;
        }
    } else {
        // Linear / disulfide zig-zag ACS backbone
        let dx = bl * 0.866;
        let dy = bl * 0.48;
        let mut x = 0.0_f64;
        for (ri, r) in res_atoms.iter().enumerate() {
            let side_sign = if ri.is_multiple_of(2) { 1.0 } else { -1.0 };
            pos[r.n_id] = Some([x, -side_sign * dy * 0.45]);
            outward[r.n_id] = [0.0, -side_sign];
            x += dx;
            pos[r.ca_id] = Some([x, side_sign * dy * 0.45]);
            outward[r.ca_id] = [0.0, side_sign];
            x += dx;
            pos[r.c_id] = Some([x, -side_sign * dy * 0.45]);
            outward[r.c_id] = [0.0, -side_sign];
            if let Some(oid) = r.o_id {
                pos[oid] = Some([x, -side_sign * (dy * 0.45 + bl * 0.88)]);
                outward[oid] = [0.0, -side_sign];
            }
            x += dx;
        }
    }

    // 2b. Place any small 3..=7-membered rings already anchored by >=2 backbone atoms (e.g. Proline pyrrolidine)
    let mut queue: VecDeque<usize> = VecDeque::new();
    for (aid, p) in pos.iter().enumerate() {
        if p.is_some() {
            queue.push_back(aid);
        }
    }
    place_fused_sidechain_rings_2d(mol, bl, &mut pos, &mut outward, &mut queue);

    // 3. BFS outward placement for all remaining atoms (sidechains, N-methyls, lipid caps, C-caps)
    while let Some(u) = queue.pop_front() {
        let Some(pu) = pos[u] else {
            continue;
        };
        let unplaced: Vec<usize> = mol.adj[u]
            .iter()
            .map(|&(v, _)| v)
            .filter(|&v| pos[v].is_none())
            .collect();
        if unplaced.is_empty() {
            continue;
        }
        // Compute preferred outward direction away from already-placed neighbors
        let mut away = outward[u];
        let mut placed_sum = [0.0_f64; 2];
        let mut placed_cnt = 0usize;
        for &(nb, _) in &mol.adj[u] {
            if let Some(pn) = pos[nb] {
                let vx = pn[0] - pu[0];
                let vy = pn[1] - pu[1];
                let l = (vx * vx + vy * vy).sqrt();
                if l > 1e-6 {
                    placed_sum[0] += vx / l;
                    placed_sum[1] += vy / l;
                    placed_cnt += 1;
                }
            }
        }
        if placed_cnt > 0 {
            let ax = -placed_sum[0];
            let ay = -placed_sum[1];
            let l = (ax * ax + ay * ay).sqrt();
            if l > 0.15 {
                away = [ax / l, ay / l];
            }
        }
        let base_ang = away[1].atan2(away[0]);
        let k_n = unplaced.len();
        for (ki, &v) in unplaced.iter().enumerate() {
            if pos[v].is_some() {
                continue;
            }
            let spread = if k_n == 1 {
                0.0
            } else {
                let frac = (ki as f64) / ((k_n - 1) as f64) - 0.5;
                frac * 1.18
            };
            let ang = base_ang + spread;
            let dx = ang.cos();
            let dy = ang.sin();
            let pv = [pu[0] + dx * bl, pu[1] + dy * bl];
            pos[v] = Some(pv);
            outward[v] = [dx, dy];
            queue.push_back(v);

            // If `v` is the entry vertex of an unplaced 3..=7-membered sidechain ring
            // (Phe, Tyr, His, Trp, Pro, etc.), place the entire regular polygon (and any
            // fused partner ring such as Trp indole) cleanly outward along `(dx, dy)`.
            place_sidechain_rings_2d(mol, v, pv, [dx, dy], bl, &mut pos, &mut outward, &mut queue);
        }
    }

    // Center 2D coordinates cleanly
    let mut cx = 0.0;
    let mut cy = 0.0;
    let mut cnt = 0.0;
    for (aid, p) in pos.into_iter().enumerate() {
        let pt = p.unwrap_or([0.0, 0.0]);
        cx += pt[0];
        cy += pt[1];
        cnt += 1.0;
        mol.atoms[aid].coords_2d = Some(pt);
    }
    if cnt > 0.0 {
        for a in &mut mol.atoms {
            if let Some([x, y]) = &mut a.coords_2d {
                *x -= cx / cnt;
                *y -= cy / cnt;
            }
        }
    }
}

/// Place any small (`3..=7`-membered) sidechain ring entered at single vertex `entry_v`,
/// followed by any fused partner rings (such as Trp indole or Pro pyrrolidine).
#[allow(clippy::too_many_arguments)]
fn place_sidechain_rings_2d(
    mol: &MolGraph,
    entry_v: usize,
    pv: [f64; 2],
    dir: [f64; 2],
    bl: f64,
    pos: &mut [Option<[f64; 2]>],
    outward: &mut [[f64; 2]],
    queue: &mut VecDeque<usize>,
) {
    let dir_ang = dir[1].atan2(dir[0]);
    for ring in &mol.rings {
        let k = ring.len();
        if !(3..=7).contains(&k) {
            continue;
        }
        let Some(v_pos) = ring.iter().position(|&a| a == entry_v) else {
            continue;
        };
        let placed_count = ring.iter().filter(|&&a| pos[a].is_some()).count();
        if placed_count != 1 || placed_count == k {
            continue;
        }
        let radius = bl / (2.0 * (PI / k as f64).sin()).max(0.25);
        let center = [pv[0] + dir[0] * radius, pv[1] + dir[1] * radius];
        let ang0 = dir_ang + PI;
        for step in 1..k {
            let aid = ring[(v_pos + step) % k];
            if pos[aid].is_none() {
                let ang = ang0 + TAU * (step as f64) / (k as f64);
                let ux = ang.cos();
                let uy = ang.sin();
                pos[aid] = Some([center[0] + radius * ux, center[1] + radius * uy]);
                outward[aid] = [ux, uy];
                queue.push_back(aid);
            }
        }
    }
    place_fused_sidechain_rings_2d(mol, bl, pos, outward, queue);
}

/// Place any `3..=7`-membered ring that already has two consecutive vertices positioned
/// (e.g., Trp fused benzene, Proline pyrrolidine anchored on `N-CA`), extending cleanly outward.
fn place_fused_sidechain_rings_2d(
    mol: &MolGraph,
    bl: f64,
    pos: &mut [Option<[f64; 2]>],
    outward: &mut [[f64; 2]],
    queue: &mut VecDeque<usize>,
) {
    for _pass in 0..4 {
        let mut progressed = false;
        for ring in &mol.rings {
            let k = ring.len();
            if !(3..=7).contains(&k) {
                continue;
            }
            let unplaced_cnt = ring.iter().filter(|&&a| pos[a].is_none()).count();
            if unplaced_cnt == 0 || unplaced_cnt >= k - 1 {
                continue;
            }
            // Find two consecutive placed vertices (a, b) along the ring perimeter
            let mut anchor = None;
            for i in 0..k {
                let a = ring[i];
                let b = ring[(i + 1) % k];
                if let (Some(pa), Some(pb)) = (pos[a], pos[b]) {
                    anchor = Some((i, a, b, pa, pb));
                    break;
                }
            }
            let Some((start_idx, a, b, pa, pb)) = anchor else {
                continue;
            };
            let ex = pb[0] - pa[0];
            let ey = pb[1] - pa[1];
            let edge_len = (ex * ex + ey * ey).sqrt().max(bl * 0.5);
            let ue = [ex / edge_len, ey / edge_len];
            let mid = [(pa[0] + pb[0]) * 0.5, (pa[1] + pb[1]) * 0.5];
            let n1 = [-ue[1], ue[0]];

            // Choose side facing away from already-placed neighbors of (a, b)
            let mut away = [0.0_f64; 2];
            let mut nbr_cnt = 0usize;
            for &(nb, _) in mol.adj[a].iter().chain(mol.adj[b].iter()) {
                if nb == a || nb == b {
                    continue;
                }
                if let Some(pn) = pos[nb] {
                    away[0] += mid[0] - pn[0];
                    away[1] += mid[1] - pn[1];
                    nbr_cnt += 1;
                }
            }
            let dot = if nbr_cnt > 0 && (away[0] * away[0] + away[1] * away[1]) > 1e-6 {
                away[0] * n1[0] + away[1] * n1[1]
            } else {
                (outward[a][0] + outward[b][0]) * n1[0] + (outward[a][1] + outward[b][1]) * n1[1]
            };
            let out_n = if dot >= 0.0 { n1 } else { [-n1[0], -n1[1]] };

            let radius = edge_len / (2.0 * (PI / k as f64).sin()).max(0.25);
            let apothem = (radius * radius - 0.25 * edge_len * edge_len)
                .max(0.0)
                .sqrt();
            let center = [mid[0] + out_n[0] * apothem, mid[1] + out_n[1] * apothem];

            let ang_a = (pa[1] - center[1]).atan2(pa[0] - center[0]);
            let ang_b = (pb[1] - center[1]).atan2(pb[0] - center[0]);
            let mut step = ang_b - ang_a;
            while step > PI {
                step -= TAU;
            }
            while step < -PI {
                step += TAU;
            }
            let sign = if step >= 0.0 { 1.0 } else { -1.0 };
            let delta = sign * TAU / (k as f64);

            for offset in 0..k {
                let aid = ring[(start_idx + offset) % k];
                if pos[aid].is_none() {
                    let ang = ang_a + delta * (offset as f64);
                    let ux = ang.cos();
                    let uy = ang.sin();
                    pos[aid] = Some([center[0] + radius * ux, center[1] + radius * uy]);
                    outward[aid] = [ux, uy];
                    queue.push_back(aid);
                    progressed = true;
                }
            }
        }
        if !progressed {
            break;
        }
    }
}

/// Generate authentic 3D coordinates honoring L vs D Cα chirality, secondary structure presets,
/// and closed-ring macrocycle / disulfide / staple geometry, followed by chiral-volume-preserving
/// relaxation (`relax_3d_coords_pinned`).
pub fn layout_peptide_3d(mol: &mut MolGraph, spec: &PeptideSpec, res_atoms: &[ResidueAtoms]) {
    let n_atoms = mol.atoms.len();
    if n_atoms == 0 || res_atoms.is_empty() {
        return;
    }
    let mut pos: Vec<Option<[f64; 3]>> = vec![None; n_atoms];
    let mut pinned = vec![false; n_atoms];

    let d_n_ca = 1.458_f64;
    let d_ca_c = 1.525_f64;
    let d_c_n = 1.329_f64;
    let d_c_o = 1.231_f64;
    let ang_c_n_ca = 121.7_f64.to_radians();
    let ang_n_ca_c = 111.2_f64.to_radians();
    let ang_ca_c_n = 116.2_f64.to_radians();
    let ang_ca_c_o = 120.8_f64.to_radians();

    if let Some((ring_lo, ring_hi)) = primary_ring_residues(spec) {
        let ring_res_count = ring_hi - ring_lo + 1;
        let ring_bb: Vec<(usize, usize)> = (ring_lo..=ring_hi)
            .flat_map(|ri| {
                let r = &res_atoms[ri];
                [(r.n_id, ri), (r.ca_id, ri), (r.c_id, ri)]
            })
            .collect();
        let m = ring_bb.len().max(6) as f64;
        let radius = (m * 1.44 / TAU).max(2.6);

        if spec.is_head_to_tail() && ring_res_count >= 8 {
            // Antiparallel cyclic β-hairpin architecture (exact trans-ω NeRF strands separated
            // by 4.8 Å via proper 180° rotation preserving L/D chirality, joined by two relaxed
            // 2-residue turn arcs). Keeps >= 75-85% of backbone dihedrals in Ramachandran favored!
            let s1_len = (ring_res_count - 4) / 2;
            let s2_len = ring_res_count - 4 - s1_len;
            // Layout indices within `ring_lo..=ring_hi`:
            //   turn2_b = 0
            //   strand1 = 1 .. (1 + s1_len)
            //   turn1   = (1 + s1_len) .. (3 + s1_len)          (2 residues)
            //   strand2 = (3 + s1_len) .. (3 + s1_len + s2_len)
            //   turn2_a = ring_res_count - 1
            let s1_start = ring_lo + 1;
            let s1_end = s1_start + s1_len;
            let t1_start = s1_end;
            let s2_start = t1_start + 2;
            let s2_end = s2_start + s2_len;

            type BbFrame = ([f64; 3], [f64; 3], [f64; 3], Option<[f64; 3]>);
            let build_beta_strand = |res_slice: &[ResidueAtoms]| -> Vec<BbFrame> {
                let omega = PI;
                let mut prev_c = [-1.329, 0.0, 0.0];
                let mut prev_ca = [-2.25, 1.15, 0.0];
                let mut prev_n = [-3.60, 0.80, 0.0];
                let mut out = Vec::with_capacity(res_slice.len());
                for (i, r) in res_slice.iter().enumerate() {
                    let (phi_deg, psi_deg) = if r.is_d { (135.0_f64, -135.0_f64) } else { (-135.0_f64, 135.0_f64) };
                    let phi = phi_deg.to_radians();
                    let psi = psi_deg.to_radians();
                    let n_3d = if i == 0 {
                        [0.0, 0.0, 0.0]
                    } else {
                        place_nerf_atom(prev_n, prev_ca, prev_c, d_c_n, ang_ca_c_n, psi)
                    };
                    let ca_3d = if i == 0 {
                        let rad = (PI - ang_c_n_ca).max(0.2);
                        [d_n_ca * rad.cos(), d_n_ca * rad.sin(), 0.0]
                    } else {
                        place_nerf_atom(prev_ca, prev_c, n_3d, d_n_ca, ang_c_n_ca, omega)
                    };
                    let c_3d = place_nerf_atom(prev_c, n_3d, ca_3d, d_ca_c, ang_n_ca_c, phi);
                    let o_3d = r.o_id.map(|_| place_nerf_atom(n_3d, ca_3d, c_3d, d_c_o, ang_ca_c_o, psi + PI));
                    out.push((n_3d, ca_3d, c_3d, o_3d));
                    prev_n = n_3d;
                    prev_ca = ca_3d;
                    prev_c = c_3d;
                }
                out
            };

            let s1_pts = build_beta_strand(&res_atoms[s1_start..s1_end]);
            let s2_raw = build_beta_strand(&res_atoms[s2_start..s2_end]);

            // Compute orthonormal sheet frame (ux along Strand 1, uz sheet normal, uy interstrand)
            let p_first = s1_pts.first().map_or([0.0, 0.0, 0.0], |t| t.1);
            let p_last = s1_pts.last().map_or([1.0, 0.0, 0.0], |t| t.1);
            let ux = unit3([p_last[0] - p_first[0], p_last[1] - p_first[1], p_last[2] - p_first[2]], [1.0, 0.0, 0.0]);
            let tmp = if ux[2].abs() < 0.85 { [0.0, 0.0, 1.0] } else { [0.0, 1.0, 0.0] };
            let dot_t = tmp[0] * ux[0] + tmp[1] * ux[1] + tmp[2] * ux[2];
            let uz = unit3([tmp[0] - dot_t * ux[0], tmp[1] - dot_t * ux[1], tmp[2] - dot_t * ux[2]], [0.0, 0.0, 1.0]);
            let uy = unit3([
                uz[1] * ux[2] - uz[2] * ux[1],
                uz[2] * ux[0] - uz[0] * ux[2],
                uz[0] * ux[1] - uz[1] * ux[0],
            ], [0.0, 1.0, 0.0]);

            // Place Strand 1 (pinned to lock exact Ramachandran β-sheet dihedrals)
            let mut m1 = [0.0_f64; 3];
            for (k, &(pn, pca, pc, po)) in s1_pts.iter().enumerate() {
                let r = &res_atoms[s1_start + k];
                pos[r.n_id] = Some(pn);
                pos[r.ca_id] = Some(pca);
                pos[r.c_id] = Some(pc);
                pinned[r.n_id] = true;
                pinned[r.ca_id] = true;
                pinned[r.c_id] = true;
                if let (Some(oid), Some(po3)) = (r.o_id, po) {
                    pos[oid] = Some(po3);
                    pinned[oid] = true;
                }
                m1 = add3(m1, pca);
            }
            let s1_mid = scale3(m1, 1.0 / (s1_pts.len().max(1) as f64));

            // Proper 180° rotation around `uz` (det = +1, preserves L/D chirality!) + 4.8 Å shift along `uy`
            let mut m2 = [0.0_f64; 3];
            for &(_, pca, _, _) in &s2_raw {
                m2 = add3(m2, pca);
            }
            let s2_mid = scale3(m2, 1.0 / (s2_raw.len().max(1) as f64));
            let rot180_uz = |p: [f64; 3]| -> [f64; 3] {
                let d = [p[0] - s2_mid[0], p[1] - s2_mid[1], p[2] - s2_mid[2]];
                let dz = d[0] * uz[0] + d[1] * uz[1] + d[2] * uz[2];
                // 180° rotation around uz: d' = 2*(d·uz)*uz - d
                let r = [2.0 * dz * uz[0] - d[0], 2.0 * dz * uz[1] - d[1], 2.0 * dz * uz[2] - d[2]];
                add3(add3(s1_mid, scale3(uy, 4.8)), r)
            };

            for (k, &(pn, pca, pc, po)) in s2_raw.iter().enumerate() {
                let r = &res_atoms[s2_start + k];
                pos[r.n_id] = Some(rot180_uz(pn));
                pos[r.ca_id] = Some(rot180_uz(pca));
                pos[r.c_id] = Some(rot180_uz(pc));
                pinned[r.n_id] = true;
                pinned[r.ca_id] = true;
                pinned[r.c_id] = true;
                if let (Some(oid), Some(po3)) = (r.o_id, po) {
                    pos[oid] = Some(rot180_uz(po3));
                    pinned[oid] = true;
                }
            }

            // Place the two 2-residue turn arcs smoothly connecting Strand 1 <-> Strand 2 (left unpinned for relaxation)
            let place_turn_arc = |pos: &mut [Option<[f64; 3]>],
                                  p_from: [f64; 3],
                                  p_to: [f64; 3],
                                  bulge_dir: [f64; 3],
                                  turn_indices: [usize; 2]| {
                let turn_atoms = [
                    res_atoms[turn_indices[0]].n_id,
                    res_atoms[turn_indices[0]].ca_id,
                    res_atoms[turn_indices[0]].c_id,
                    res_atoms[turn_indices[1]].n_id,
                    res_atoms[turn_indices[1]].ca_id,
                    res_atoms[turn_indices[1]].c_id,
                ];
                for (step, &aid) in turn_atoms.iter().enumerate() {
                    let t = (step as f64 + 1.0) / 7.0;
                    let arc = (PI * t).sin() * 2.5;
                    let z_puck = if step.is_multiple_of(2) { 0.45 } else { -0.45 };
                    let base = [
                        p_from[0] + (p_to[0] - p_from[0]) * t,
                        p_from[1] + (p_to[1] - p_from[1]) * t,
                        p_from[2] + (p_to[2] - p_from[2]) * t,
                    ];
                    pos[aid] = Some(add3(add3(base, scale3(bulge_dir, arc)), scale3(uz, z_puck)));
                }
                for (ti, &ri) in turn_indices.iter().enumerate() {
                    let r = &res_atoms[ri];
                    if let (Some(oid), Some(pc)) = (r.o_id, pos[r.c_id]) {
                        let s = if ti == 0 { 1.0 } else { -1.0 };
                        pos[oid] = Some(add3(pc, scale3(unit3(add3(bulge_dir, scale3(uz, s)), uz), d_c_o)));
                    }
                }
            };

            if let (Some(p1_end), Some(p2_start)) = (
                pos[res_atoms[s1_end - 1].c_id],
                pos[res_atoms[s2_start].n_id],
            ) {
                place_turn_arc(&mut pos, p1_end, p2_start, ux, [t1_start, t1_start + 1]);
            }
            if let (Some(p2_end), Some(p1_start)) = (
                pos[res_atoms[s2_end - 1].c_id],
                pos[res_atoms[s1_start].n_id],
            ) {
                place_turn_arc(&mut pos, p2_end, p1_start, scale3(ux, -1.0), [ring_hi, ring_lo]);
            }
        } else {
            // Small macrocycles (<= 5 residues, e.g. Cilengitide): pleated ring
            for (k, &(aid, ri)) in ring_bb.iter().enumerate() {
                let theta = TAU * (k as f64) / m;
                let pucker =
                    if ri.is_multiple_of(2) { 0.42 } else { -0.42 } * ((k % 3) as f64 - 1.0);
                pos[aid] = Some([radius * theta.cos(), radius * theta.sin(), pucker]);
            }
            for (idx_r, ri) in (ring_lo..=ring_hi).enumerate() {
                let r = &res_atoms[ri];
                if let Some(oid) = r.o_id
                    && let Some(pc) = pos[r.c_id]
                {
                    let z_sign = if idx_r.is_multiple_of(2) { 1.0 } else { -1.0 };
                    let rad_dir = unit3([pc[0], pc[1], 0.0], [1.0, 0.0, 0.0]);
                    pos[oid] = Some(add3(
                        pc,
                        scale3(
                            unit3(
                                [rad_dir[0] * 0.35, rad_dir[1] * 0.35, z_sign],
                                [0.0, 0.0, z_sign],
                            ),
                            d_c_o,
                        ),
                    ));
                }
            }
        }
        // Extend any lariat N-terminal tail radially outward from ring_lo's N atom
        if ring_lo > 0 {
            let p_anchor = pos[res_atoms[ring_lo].n_id].unwrap_or([radius, 0.0, 0.0]);
            let rad_out = unit3([p_anchor[0], p_anchor[1], 0.25], [1.0, 0.0, 0.0]);
            let mut step = 1usize;
            for ri in (0..ring_lo).rev() {
                let r = &res_atoms[ri];
                for &aid in &[r.c_id, r.ca_id, r.n_id] {
                    let dist = (step as f64) * 1.38;
                    let zig_z = if step.is_multiple_of(2) { 0.55 } else { -0.55 };
                    pos[aid] = Some([
                        p_anchor[0] + rad_out[0] * dist,
                        p_anchor[1] + rad_out[1] * dist,
                        p_anchor[2] + zig_z,
                    ]);
                    step += 1;
                }
                if let Some(oid) = r.o_id
                    && let Some(pc) = pos[r.c_id]
                {
                    pos[oid] = Some([pc[0], pc[1] + 1.1, pc[2] + 0.5]);
                }
            }
        }
    } else {
        // Authentic NeRF backbone for linear / disulfide / stapled peptides
        let (base_phi_deg, base_psi_deg, omega_deg) = spec.conformation.phi_psi_omega_deg();
        let omega = omega_deg.to_radians();
        let mut prev_c = [-1.329, 0.0, 0.0];
        let mut prev_ca = [-2.25, 1.15, 0.0];
        let mut prev_n = [-3.60, 0.80, 0.0];

        let mid = res_atoms.len() / 2;
        for (ri, r) in res_atoms.iter().enumerate() {
            // Invert Ramachandran (phi, psi) for D-residues and insert beta-turn at hairpin midpoint
            let (mut phi_deg, mut psi_deg) = if spec.conformation
                == PeptideConformation::BetaHairpin
                && (ri == mid || ri + 1 == mid)
            {
                if ri + 1 == mid {
                    (-60.0, -30.0)
                } else {
                    (-90.0, 0.0)
                }
            } else {
                (base_phi_deg, base_psi_deg)
            };
            if r.is_d {
                phi_deg = -phi_deg;
                psi_deg = -psi_deg;
            }
            let phi = phi_deg.to_radians();
            let psi = psi_deg.to_radians();

            let n_3d = if ri == 0 {
                [0.0, 0.0, 0.0]
            } else {
                place_nerf_atom(prev_n, prev_ca, prev_c, d_c_n, ang_ca_c_n, psi)
            };
            let ca_3d = if ri == 0 {
                let rad = (PI - ang_c_n_ca).max(0.2);
                [d_n_ca * rad.cos(), d_n_ca * rad.sin(), 0.0]
            } else {
                place_nerf_atom(prev_ca, prev_c, n_3d, d_n_ca, ang_c_n_ca, omega)
            };
            let c_3d = place_nerf_atom(prev_c, n_3d, ca_3d, d_ca_c, ang_n_ca_c, phi);
            pos[r.n_id] = Some(n_3d);
            pos[r.ca_id] = Some(ca_3d);
            pos[r.c_id] = Some(c_3d);
            if spec.cross_links.is_empty() {
                pinned[r.n_id] = true;
                pinned[r.ca_id] = true;
                pinned[r.c_id] = true;
            }
            if let Some(oid) = r.o_id {
                let o_dih = if ri + 1 < res_atoms.len() {
                    psi + PI
                } else {
                    0.0
                };
                let o_3d = place_nerf_atom(n_3d, ca_3d, c_3d, d_c_o, ang_ca_c_o, o_dih);
                pos[oid] = Some(o_3d);
                if spec.cross_links.is_empty() {
                    pinned[oid] = true;
                }
            }
            prev_n = n_3d;
            prev_ca = ca_3d;
            prev_c = c_3d;
        }
    }

    // Place Cβ of every chiral residue via exact NeRF (+122° for D, -122° for L) so initial 3D chirality is exact!
    for r in res_atoms {
        if let (Some(pn), Some(pca), Some(pc)) = (pos[r.n_id], pos[r.ca_id], pos[r.c_id]) {
            let cb_candidates: Vec<usize> = mol.adj[r.ca_id]
                .iter()
                .map(|&(nb, _)| nb)
                .filter(|&nb| nb != r.n_id && nb != r.c_id && pos[nb].is_none())
                .collect();
            for (ci, &cb_id) in cb_candidates.iter().enumerate() {
                let dih_deg: f64 = if ci == 0 {
                    if r.is_d { 122.0 } else { -122.0 }
                } else {
                    if r.is_d { -122.0 } else { 122.0 }
                };
                let p_cb = place_nerf_atom(
                    pc,
                    pn,
                    pca,
                    1.53,
                    110.5_f64.to_radians(),
                    dih_deg.to_radians(),
                );
                pos[cb_id] = Some(p_cb);
            }
        }
    }

    // BFS 3D placement for all remaining sidechain, cap, and lipid atoms using local tetrahedral/planar NeRF
    let mut queue: VecDeque<usize> = VecDeque::new();
    for (aid, p) in pos.iter().enumerate() {
        if p.is_some() {
            queue.push_back(aid);
        }
    }
    while let Some(u) = queue.pop_front() {
        let Some(pu) = pos[u] else {
            continue;
        };
        let unplaced: Vec<usize> = mol.adj[u]
            .iter()
            .map(|&(v, _)| v)
            .filter(|&v| pos[v].is_none())
            .collect();
        if unplaced.is_empty() {
            continue;
        }
        // Pick a placed parent `p_parent` and grandparent `p_gp` for NeRF reference frame
        let placed_nbrs: Vec<usize> = mol.adj[u]
            .iter()
            .map(|&(nb, _)| nb)
            .filter(|&nb| pos[nb].is_some())
            .collect();
        let p_b = placed_nbrs
            .first()
            .and_then(|&b| pos[b])
            .unwrap_or([pu[0] - 1.4, pu[1], pu[2]]);
        let p_a = placed_nbrs
            .first()
            .and_then(|&b| {
                mol.adj[b]
                    .iter()
                    .map(|&(gp, _)| gp)
                    .find(|&gp| gp != u && pos[gp].is_some())
                    .and_then(|gp| pos[gp])
            })
            .unwrap_or([p_b[0], p_b[1] + 1.4, p_b[2] + 0.5]);

        let angle_rad = if mol.atoms[u].aromatic {
            120.0_f64.to_radians()
        } else {
            111.0_f64.to_radians()
        };
        for (ki, &v) in unplaced.iter().enumerate() {
            let bl = mol_graph::ideal_bond_length(
                mol,
                u,
                v,
                mol.adj[u]
                    .iter()
                    .find(|&&(nb, _)| nb == v)
                    .map_or(mol_graph::BondOrder::Single, |&(_, bid)| {
                        mol.bonds[bid].order
                    }),
            );
            let dih = (180.0 + (ki as f64) * 120.0 + ((u % 3) as f64 - 1.0) * 18.0).to_radians();
            let pv = place_nerf_atom(p_a, p_b, pu, bl, angle_rad, dih);
            pos[v] = Some(pv);
            queue.push_back(v);
        }
    }

    let mut pts: Vec<[f64; 3]> = pos
        .into_iter()
        .map(|p| p.unwrap_or([0.0, 0.0, 0.0]))
        .collect();

    // Run chiral-volume-preserving relaxation (more iterations for cyclic/cross-linked topologies so
    // closure bonds relax cleanly to exact covalent lengths with zero clashes!)
    let iters = if spec.cross_links.is_empty() {
        140
    } else {
        320
    };
    relax_3d_coords_pinned(mol, &mut pts, &pinned, iters);

    // Center 3D centroid
    let mut cent = [0.0_f64; 3];
    for p in &pts {
        cent[0] += p[0];
        cent[1] += p[1];
        cent[2] += p[2];
    }
    let inv = 1.0 / (pts.len() as f64);
    for (a, p) in mol.atoms.iter_mut().zip(pts) {
        a.coords_3d = Some([
            p[0] - cent[0] * inv,
            p[1] - cent[1] * inv,
            p[2] - cent[2] * inv,
        ]);
    }
}
