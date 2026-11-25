//! Pure-Rust interactive 3D molecular ball-and-stick + biomolecule backbone ribbon viewer
//! rendered directly on `egui::Painter` (zero extra GPU/shader dependencies).

use std::f64::consts::PI;

use eframe::egui;
use egui::{Color32, Pos2, Rect, Response, Sense, Stroke, Vec2};
use mol_graph::{BondOrder, MolGraph, compute_3d_coords};

#[derive(Debug, Clone)]
pub struct Viewer3dState {
    pub yaw: f32,
    pub pitch: f32,
    pub zoom: f32,
    pub show_ribbon: bool,
    pub show_labels: bool,
    pub show_hydrogens: bool,
    pub auto_spin: bool,
}

impl Default for Viewer3dState {
    fn default() -> Self {
        Self {
            yaw: 0.45,
            pitch: 0.35,
            zoom: 1.0,
            show_ribbon: true,
            show_labels: true,
            show_hydrogens: true,
            auto_spin: false,
        }
    }
}

impl Viewer3dState {
    pub fn reset(&mut self) {
        self.yaw = 0.45;
        self.pitch = 0.35;
        self.zoom = 1.0;
    }
}

fn elem_color(symbol: &str, dark_mode: bool) -> Color32 {
    match symbol {
        "C" => {
            if dark_mode {
                Color32::from_rgb(155, 165, 178)
            } else {
                Color32::from_rgb(76, 82, 92)
            }
        }
        "N" => Color32::from_rgb(52, 122, 240),
        "O" => Color32::from_rgb(225, 58, 52),
        "S" => Color32::from_rgb(228, 184, 38),
        "P" => Color32::from_rgb(232, 128, 30),
        "F" => Color32::from_rgb(96, 208, 108),
        "Cl" => Color32::from_rgb(40, 190, 82),
        "Br" => Color32::from_rgb(168, 48, 42),
        "I" => Color32::from_rgb(148, 48, 178),
        "H" => {
            if dark_mode {
                Color32::from_rgb(230, 234, 240)
            } else {
                Color32::from_rgb(242, 245, 250)
            }
        }
        "Fe" => Color32::from_rgb(208, 108, 52),
        _ => Color32::from_rgb(190, 120, 215),
    }
}

fn elem_radius(symbol: &str) -> f32 {
    match symbol {
        "H" => 0.25,
        "F" | "O" | "N" => 0.40,
        "C" => 0.43,
        "S" | "P" | "Cl" => 0.50,
        "Br" | "I" | "Fe" => 0.56,
        _ => 0.44,
    }
}

fn scale_rgb(c: Color32, f: f32) -> Color32 {
    let f = f.clamp(0.0, 2.2);
    Color32::from_rgb(
        ((c.r() as f32 * f).clamp(0.0, 255.0)) as u8,
        ((c.g() as f32 * f).clamp(0.0, 255.0)) as u8,
        ((c.b() as f32 * f).clamp(0.0, 255.0)) as u8,
    )
}

fn mix_white(c: Color32, t: f32) -> Color32 {
    let t = t.clamp(0.0, 1.0);
    Color32::from_rgb(
        (c.r() as f32 * (1.0 - t) + 255.0 * t) as u8,
        (c.g() as f32 * (1.0 - t) + 255.0 * t) as u8,
        (c.b() as f32 * (1.0 - t) + 255.0 * t) as u8,
    )
}

fn rainbow_color(t: f32) -> Color32 {
    let t = t.clamp(0.0, 1.0);
    let (r, g, b) = if t < 0.25 {
        let u = t / 0.25;
        (45.0, 110.0 + 125.0 * u, 240.0)
    } else if t < 0.50 {
        let u = (t - 0.25) / 0.25;
        (45.0, 235.0, 240.0 - 140.0 * u)
    } else if t < 0.75 {
        let u = (t - 0.50) / 0.25;
        (45.0 + 200.0 * u, 235.0 - 25.0 * u, 100.0 - 65.0 * u)
    } else {
        let u = (t - 0.75) / 0.25;
        (245.0, 210.0 - 155.0 * u, 35.0)
    };
    Color32::from_rgb(r as u8, g as u8, b as u8)
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

fn is_atom_sp(mol: &MolGraph, i: usize) -> bool {
    mol.adj[i]
        .iter()
        .any(|&(_, bid)| mol.bonds[bid].order == BondOrder::Triple)
        || mol.adj[i]
            .iter()
            .filter(|&&(_, bid)| mol.bonds[bid].order == BondOrder::Double)
            .count()
            >= 2
}

fn is_atom_sp2(mol: &MolGraph, i: usize) -> bool {
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

/// Generate clean visual 3D positions `(parent_heavy_idx, xyz)` for implicit hydrogens
/// so tetrahedral, trigonal-planar, and `-OH`/`-NH` groups display accurately in 3D.
fn build_visual_hydrogens(mol: &MolGraph, heavy_coords: &[[f64; 3]]) -> Vec<(usize, [f64; 3])> {
    let mut out = Vec::new();
    for (i, a) in mol.atoms.iter().enumerate() {
        let nh = usize::from(a.implicit_h);
        if nh == 0 || a.symbol == "H" {
            continue;
        }
        let bl = match a.symbol.as_str() {
            "O" => 0.96,
            "N" => 1.01,
            "S" => 1.34,
            "P" => 1.42,
            _ => 1.09,
        };
        let p = heavy_coords[i];
        let nbrs: Vec<usize> = mol.adj[i].iter().map(|&(nb, _)| nb).collect();
        let sp = is_atom_sp(mol, i);
        let sp2 = is_atom_sp2(mol, i);

        match nbrs.len() {
            0 => {
                // Isolated atom (CH4, NH3, H2O)
                let dirs: &[[f64; 3]] = match nh {
                    1 => &[[1.0, 0.0, 0.0]],
                    2 => &[[0.8165, 0.5774, 0.0], [-0.8165, 0.5774, 0.0]],
                    3 => &[
                        [0.9428, 0.0, -0.3333],
                        [-0.4714, 0.8165, -0.3333],
                        [-0.4714, -0.8165, -0.3333],
                    ],
                    _ => &[
                        [0.5774, 0.5774, 0.5774],
                        [0.5774, -0.5774, -0.5774],
                        [-0.5774, 0.5774, -0.5774],
                        [-0.5774, -0.5774, 0.5774],
                    ],
                };
                for d in dirs.iter().take(nh) {
                    out.push((i, add3(p, scale3(unit3(*d, [1.0, 0.0, 0.0]), bl))));
                }
            }
            1 => {
                let par = nbrs[0];
                let fwd = unit3(sub3(p, heavy_coords[par]), [1.0, 0.0, 0.0]);
                // Reference vector from grandparent for staggered placement
                let ref_v = mol.adj[par].iter().find_map(|&(gp, _)| {
                    if gp != i {
                        let v = sub3(heavy_coords[gp], heavy_coords[par]);
                        let perp = sub3(v, scale3(fwd, dot3(v, fwd)));
                        if norm3(perp) > 1e-4 {
                            return Some(unit3(perp, [0.0, 1.0, 0.0]));
                        }
                    }
                    None
                });
                let (u_o, v_o) = match ref_v {
                    Some(rv) => {
                        let uo = unit3(rv, [0.0, 1.0, 0.0]);
                        let vo = unit3(cross3(fwd, uo), [0.0, 0.0, 1.0]);
                        (uo, vo)
                    }
                    None => ortho_pair(fwd),
                };
                if sp && nh == 1 {
                    out.push((i, add3(p, scale3(fwd, bl))));
                } else if sp2 {
                    let (st, ct) = (PI / 3.0).sin_cos();
                    for k in 0..nh {
                        let phi = if nh == 1 {
                            PI
                        } else if k == 0 {
                            0.0
                        } else {
                            PI
                        };
                        let orth = add3(scale3(u_o, phi.cos()), scale3(v_o, phi.sin()));
                        let dir = unit3(add3(scale3(fwd, ct), scale3(orth, st)), fwd);
                        out.push((i, add3(p, scale3(dir, bl))));
                    }
                } else {
                    // sp3: -OH (1H), -NH2 (2H), -CH3 (3H), staggered anti to grandparent!
                    let tilt = 70.529_f64.to_radians();
                    let (st, ct) = tilt.sin_cos();
                    for k in 0..nh {
                        let phi = match nh {
                            1 => PI,
                            2 => PI + if k == 0 { PI / 3.0 } else { -PI / 3.0 },
                            _ => PI + 2.0 * PI * (k as f64) / (nh as f64),
                        };
                        let orth = add3(scale3(u_o, phi.cos()), scale3(v_o, phi.sin()));
                        let dir = unit3(add3(scale3(fwd, ct), scale3(orth, st)), fwd);
                        out.push((i, add3(p, scale3(dir, bl))));
                    }
                }
            }
            2 => {
                let v1 = unit3(sub3(heavy_coords[nbrs[0]], p), [1.0, 0.0, 0.0]);
                let v2 = unit3(sub3(heavy_coords[nbrs[1]], p), [0.0, 1.0, 0.0]);
                let bis = unit3(scale3(add3(v1, v2), -1.0), ortho_pair(v1).0);
                if sp2 || nh == 1 && is_atom_sp2(mol, i) {
                    out.push((i, add3(p, scale3(bis, bl))));
                } else {
                    let perp = unit3(cross3(v1, v2), ortho_pair(bis).1);
                    let d_bv = dot3(bis, v1);
                    let c_in = if d_bv.abs() > 0.15 {
                        ((-1.0 / 3.0) / d_bv).clamp(0.15, 0.88)
                    } else {
                        0.577
                    };
                    let c_p = (1.0 - c_in * c_in).max(0.0).sqrt();
                    for k in 0..nh.min(2) {
                        let sign = if k == 0 { 1.0 } else { -1.0 };
                        let dir = unit3(add3(scale3(bis, c_in), scale3(perp, sign * c_p)), bis);
                        out.push((i, add3(p, scale3(dir, bl))));
                    }
                }
            }
            _ => {
                // Degree 3 heavy atom with 1 implicit H (e.g. chiral C-H or tertiary C-H)
                let mut sum = [0.0, 0.0, 0.0];
                for &nb in &nbrs {
                    sum = add3(sum, unit3(sub3(heavy_coords[nb], p), [1.0, 0.0, 0.0]));
                }
                let fallback = if nbrs.len() >= 2 {
                    let v1 = sub3(heavy_coords[nbrs[0]], p);
                    let v2 = sub3(heavy_coords[nbrs[1]], p);
                    unit3(cross3(v1, v2), [0.0, 0.0, 1.0])
                } else {
                    [0.0, 0.0, 1.0]
                };
                let dir = unit3(scale3(sum, -1.0), fallback);
                out.push((i, add3(p, scale3(dir, bl))));
            }
        }
    }
    out
}

#[derive(Clone)]
enum DrawItem {
    RibbonSeg {
        z: f32,
        p1: Pos2,
        p2: Pos2,
        width: f32,
        color: Color32,
    },
    Bond {
        z: f32,
        p1: Pos2,
        p2: Pos2,
        r1: f32,
        r2: f32,
        order: BondOrder,
        width: f32,
        c1: Color32,
        c2: Color32,
        inward_dir: Option<Vec2>,
    },
    Atom {
        z: f32,
        pos: Pos2,
        radius: f32,
        color: Color32,
        label: Option<String>,
        hovered: bool,
        is_hydrogen: bool,
    },
}

impl DrawItem {
    fn depth(&self) -> f32 {
        match self {
            Self::RibbonSeg { z, .. } | Self::Bond { z, .. } | Self::Atom { z, .. } => *z,
        }
    }
}

fn draw_dashed_segment(painter: &egui::Painter, a: Pos2, b: Pos2, stroke: Stroke) {
    let d = b - a;
    let len = d.length();
    if len < 4.0 {
        return;
    }
    let dir = d / len;
    let dash = 4.2_f32;
    let gap = 3.2_f32;
    let mut t = 1.5_f32;
    while t < len - 1.5 {
        let t_end = (t + dash).min(len - 1.5);
        painter.line_segment([a + dir * t, a + dir * t_end], stroke);
        t += dash + gap;
    }
}

/// Draw an interactive 3D viewport into `rect`, returning the `egui::Response` for drag/zoom handling.
pub fn render_3d_viewport(
    ui: &mut egui::Ui,
    mol: &MolGraph,
    state: &mut Viewer3dState,
    rect: Rect,
    dark_mode: bool,
) -> Response {
    let resp = ui.interact(
        rect,
        ui.id().with("chemcraft_3d_view"),
        Sense::click_and_drag(),
    );

    if resp.dragged() {
        let d = resp.drag_delta();
        state.yaw += d.x * 0.011;
        state.pitch = (state.pitch + d.y * 0.011).clamp(-1.52, 1.52);
    }
    if resp.double_clicked() {
        state.reset();
    }
    if resp.hovered() {
        let scroll = ui.input(|i| i.smooth_scroll_delta.y);
        if scroll.abs() > 0.1 {
            state.zoom = (state.zoom * (1.0 + scroll * 0.0022)).clamp(0.35, 4.5);
        }
    }
    if state.auto_spin {
        state.yaw += 0.012;
        ui.ctx().request_repaint();
    }

    let painter = ui.painter_at(rect);
    let bg = if dark_mode {
        Color32::from_rgb(16, 19, 26)
    } else {
        Color32::from_rgb(248, 250, 253)
    };
    painter.rect_filled(rect, 4.0, bg);
    painter.rect_stroke(
        rect,
        4.0,
        Stroke::new(
            1.0_f32,
            Color32::from_gray(if dark_mode { 75 } else { 170 }),
        ),
        egui::StrokeKind::Inside,
    );

    if mol.atoms.is_empty() {
        painter.text(
            rect.center(),
            egui::Align2::CENTER_CENTER,
            "3D Viewer — Load SMILES or Draw a Molecule",
            egui::FontId::proportional(14.0),
            Color32::GRAY,
        );
        return resp;
    }

    let coords = compute_3d_coords(mol);
    let hydrogens = if state.show_hydrogens && mol.atoms.len() <= 140 {
        build_visual_hydrogens(mol, &coords)
    } else {
        Vec::new()
    };

    // Compute centered radius over heavy atoms (plus hydrogens for framing)
    let mut cent = [0.0_f64; 3];
    for p in &coords {
        cent[0] += p[0];
        cent[1] += p[1];
        cent[2] += p[2];
    }
    let inv_n = 1.0 / coords.len() as f64;
    cent[0] *= inv_n;
    cent[1] *= inv_n;
    cent[2] *= inv_n;

    let mut max_r = 2.0_f64;
    for p in &coords {
        let d = sub3(*p, cent);
        max_r = max_r.max(norm3(d));
    }
    for &(_, hp) in &hydrogens {
        let d = sub3(hp, cent);
        max_r = max_r.max(norm3(d));
    }

    let (sy, cy) = state.yaw.sin_cos();
    let (sp, cp) = state.pitch.sin_cos();
    let cam_dist = (max_r as f32 * 2.8 + 8.5) / state.zoom;
    let base_scale = rect.width().min(rect.height()) * 0.40 * state.zoom / (max_r as f32);
    let focal = cam_dist;

    let project_pt = |p: [f64; 3]| -> (Pos2, f32, f32) {
        let x0 = (p[0] - cent[0]) as f32;
        let y0 = (p[1] - cent[1]) as f32;
        let z0 = (p[2] - cent[2]) as f32;

        let x1 = x0 * cy + z0 * sy;
        let z1 = -x0 * sy + z0 * cy;
        let y2 = y0 * cp - z1 * sp;
        let z2 = y0 * sp + z1 * cp;

        let depth = (cam_dist - z2).max(1.5);
        let persp = focal / depth;
        let sx = rect.center().x + x1 * base_scale * persp;
        let sy_scr = rect.center().y - y2 * base_scale * persp;
        (Pos2::new(sx, sy_scr), z2, persp)
    };

    let proj: Vec<(Pos2, f32, f32)> = coords.iter().map(|&p| project_pt(p)).collect();
    let heavy_radii: Vec<f32> = mol
        .atoms
        .iter()
        .enumerate()
        .map(|(i, a)| (elem_radius(&a.symbol) * base_scale * 0.50 * proj[i].2).clamp(5.5, 26.0))
        .collect();

    // Hover detection (closest frontmost projected heavy atom)
    let mut hovered_atom: Option<(usize, f32)> = None;
    if let Some(hp) = resp.hover_pos() {
        for (i, (pos, z, _)) in proj.iter().enumerate() {
            let r = heavy_radii[i].max(8.0);
            if pos.distance(hp) <= r + 3.0 && hovered_atom.is_none_or(|(_, best_z)| *z > best_z) {
                hovered_atom = Some((i, *z));
            }
        }
    }

    // Precompute ring centroids in 3D so aromatic inner dashed lines always sit inside their ring
    let ring_centroids: Vec<[f64; 3]> = mol
        .rings
        .iter()
        .map(|r| {
            let mut c = [0.0; 3];
            for &a in r {
                c = add3(c, coords[a]);
            }
            scale3(c, 1.0 / r.len().max(1) as f64)
        })
        .collect();

    let mut items: Vec<DrawItem> =
        Vec::with_capacity(mol.atoms.len() + mol.bonds.len() + hydrogens.len() * 2 + 32);

    // Optional C-alpha backbone spline ribbon for PDB-loaded biomolecules
    if state.show_ribbon {
        let ca_indices: Vec<usize> = mol
            .atoms
            .iter()
            .enumerate()
            .filter_map(|(idx, a)| {
                a.pdb_meta
                    .as_ref()
                    .filter(|m| !m.is_hetatm && m.atom_name == "CA")
                    .map(|_| idx)
            })
            .collect();
        if ca_indices.len() >= 2 {
            for w in 0..(ca_indices.len() - 1) {
                let a = ca_indices[w];
                let b = ca_indices[w + 1];
                let (p1, z1, per1) = proj[a];
                let (p2, z2, per2) = proj[b];
                let t = w as f32 / (ca_indices.len() - 1).max(1) as f32;
                let width = (5.2 * (per1 + per2) * 0.5).clamp(2.5, 11.0);
                items.push(DrawItem::RibbonSeg {
                    z: (z1 + z2) * 0.5 - 0.15,
                    p1,
                    p2,
                    width,
                    color: rainbow_color(t),
                });
            }
        }
    }

    // Heavy-atom bonds
    for b in &mol.bonds {
        let (p1, z1, per1) = proj[b.u];
        let (p2, z2, per2) = proj[b.v];
        let avg_z = (z1 + z2) * 0.5;
        let depth_shade =
            (0.75 + (avg_z / (max_r as f32 * 2.0)).clamp(-0.26, 0.24)).clamp(0.45, 1.0);
        let c1 = scale_rgb(elem_color(&mol.atoms[b.u].symbol, dark_mode), depth_shade);
        let c2 = scale_rgb(elem_color(&mol.atoms[b.v].symbol, dark_mode), depth_shade);
        let w = (0.14 * base_scale * (per1 + per2) * 0.5).clamp(2.0, 6.2);

        // Compute inward screen direction for aromatic rings or double-bond molecular plane
        let inward_dir = if b.order == BondOrder::Aromatic {
            mol.rings
                .iter()
                .enumerate()
                .find(|(_, r)| r.contains(&b.u) && r.contains(&b.v))
                .map(|(ri, _)| {
                    let (rc_scr, _, _) = project_pt(ring_centroids[ri]);
                    let mid = p1 + (p2 - p1) * 0.5;
                    let d = p2 - p1;
                    let len = d.length().max(1e-3);
                    let perp = Vec2::new(-d.y / len, d.x / len);
                    let to_c = rc_scr - mid;
                    if to_c.x * perp.x + to_c.y * perp.y >= 0.0 {
                        perp
                    } else {
                        -perp
                    }
                })
        } else {
            None
        };

        items.push(DrawItem::Bond {
            z: avg_z,
            p1,
            p2,
            r1: heavy_radii[b.u],
            r2: heavy_radii[b.v],
            order: b.order,
            width: w,
            c1,
            c2,
            inward_dir,
        });
    }

    // Heavy atoms
    for (i, a) in mol.atoms.iter().enumerate() {
        let (pos, z, _) = proj[i];
        let r = heavy_radii[i];
        let depth_shade = (0.88 + (z / (max_r as f32 * 2.2)).clamp(-0.28, 0.22)).clamp(0.48, 1.05);
        let col = scale_rgb(elem_color(&a.symbol, dark_mode), depth_shade);
        let label = if state.show_labels && a.symbol != "C" && a.symbol != "H" {
            Some(a.symbol.clone())
        } else {
            None
        };
        items.push(DrawItem::Atom {
            z,
            pos,
            radius: r,
            color: col,
            label,
            hovered: hovered_atom.is_some_and(|(hi, _)| hi == i),
            is_hydrogen: a.symbol == "H",
        });
    }

    // Visual implicit hydrogens + X-H bonds
    for &(parent_idx, h_xyz) in &hydrogens {
        let (p1, z1, per1) = proj[parent_idx];
        let (p2, z2, per2) = project_pt(h_xyz);
        let h_rad = (elem_radius("H") * base_scale * 0.44 * per2).clamp(3.2, 13.0);
        let avg_z = (z1 + z2) * 0.5;
        let depth_shade =
            (0.78 + (avg_z / (max_r as f32 * 2.0)).clamp(-0.24, 0.22)).clamp(0.50, 1.0);
        let c1 = scale_rgb(
            elem_color(&mol.atoms[parent_idx].symbol, dark_mode),
            depth_shade,
        );
        let c2 = scale_rgb(elem_color("H", dark_mode), depth_shade);
        let w = (0.095 * base_scale * (per1 + per2) * 0.5).clamp(1.3, 3.8);

        items.push(DrawItem::Bond {
            z: avg_z,
            p1,
            p2,
            r1: heavy_radii[parent_idx],
            r2: h_rad,
            order: BondOrder::Single,
            width: w,
            c1,
            c2,
            inward_dir: None,
        });
        items.push(DrawItem::Atom {
            z: z2,
            pos: p2,
            radius: h_rad,
            color: c2,
            label: None,
            hovered: false,
            is_hydrogen: true,
        });
    }

    // Back-to-front painter's sort (smaller z = further away; larger z = closer to camera)
    items.sort_by(|a, b| a.depth().total_cmp(&b.depth()));

    let bond_outline = if dark_mode {
        Color32::from_rgba_unmultiplied(245, 245, 255, 50)
    } else {
        Color32::from_rgba_unmultiplied(25, 30, 42, 65)
    };

    for item in items {
        match item {
            DrawItem::RibbonSeg {
                p1,
                p2,
                width,
                color,
                ..
            } => {
                painter.line_segment(
                    [p1, p2],
                    Stroke::new(width + 2.0, color.linear_multiply(0.35)),
                );
                painter.line_segment([p1, p2], Stroke::new(width, color));
            }
            DrawItem::Bond {
                p1,
                p2,
                r1,
                r2,
                order,
                width,
                c1,
                c2,
                inward_dir,
                ..
            } => {
                let d = p2 - p1;
                let len = d.length();
                if len < 1.5 {
                    continue;
                }
                let unit = d / len;
                // Clip bond sticks cleanly to sphere surfaces so bonds never pierce over atoms
                let trim1 = (r1 * 0.68).min(len * 0.42);
                let trim2 = (r2 * 0.68).min(len * 0.42);
                let q1 = p1 + unit * trim1;
                let q2 = p2 - unit * trim2;
                if (q2 - q1).length() < 1.0 {
                    continue;
                }
                let perp = Vec2::new(-unit.y, unit.x);
                let mid = q1 + (q2 - q1) * 0.5;
                match order {
                    BondOrder::Aromatic => {
                        // Standard 3D aromatic bond: 1 solid primary stick + 1 dashed inner resonance line
                        painter.line_segment([q1, q2], Stroke::new(width + 1.0, bond_outline));
                        painter.line_segment([q1, mid], Stroke::new(width * 0.92, c1));
                        painter.line_segment([mid, q2], Stroke::new(width * 0.92, c2));
                        let in_v = inward_dir.unwrap_or(perp);
                        let off = in_v * (width * 1.05).max(3.2);
                        let inner_col = if dark_mode {
                            Color32::from_rgba_unmultiplied(210, 220, 235, 195)
                        } else {
                            Color32::from_rgba_unmultiplied(70, 80, 96, 185)
                        };
                        draw_dashed_segment(
                            &painter,
                            q1 + off + unit * 1.5,
                            q2 + off - unit * 1.5,
                            Stroke::new((width * 0.48).max(1.2), inner_col),
                        );
                    }
                    BondOrder::Double => {
                        let off = perp * (width * 0.65).max(2.0);
                        painter.line_segment(
                            [q1 + off, q2 + off],
                            Stroke::new(width * 0.72 + 0.8, bond_outline),
                        );
                        painter.line_segment(
                            [q1 - off, q2 - off],
                            Stroke::new(width * 0.72 + 0.8, bond_outline),
                        );
                        painter.line_segment([q1 + off, mid + off], Stroke::new(width * 0.66, c1));
                        painter.line_segment([mid + off, q2 + off], Stroke::new(width * 0.66, c2));
                        painter.line_segment([q1 - off, mid - off], Stroke::new(width * 0.66, c1));
                        painter.line_segment([mid - off, q2 - off], Stroke::new(width * 0.66, c2));
                    }
                    BondOrder::Triple => {
                        let off = perp * (width * 0.95).max(2.6);
                        painter.line_segment([q1, mid], Stroke::new(width * 0.58, c1));
                        painter.line_segment([mid, q2], Stroke::new(width * 0.58, c2));
                        painter.line_segment([q1 + off, mid + off], Stroke::new(width * 0.44, c1));
                        painter.line_segment([mid + off, q2 + off], Stroke::new(width * 0.44, c2));
                        painter.line_segment([q1 - off, mid - off], Stroke::new(width * 0.44, c1));
                        painter.line_segment([mid - off, q2 - off], Stroke::new(width * 0.44, c2));
                    }
                    _ => {
                        painter.line_segment([q1, q2], Stroke::new(width + 1.1, bond_outline));
                        painter.line_segment([q1, mid], Stroke::new(width, c1));
                        painter.line_segment([mid, q2], Stroke::new(width, c2));
                    }
                }
            }
            DrawItem::Atom {
                pos,
                radius,
                color,
                label,
                hovered,
                is_hydrogen,
                ..
            } => {
                // Shaded 3D sphere (ambient base + diffuse offset + specular highlight)
                painter.circle_filled(pos, radius, scale_rgb(color, 0.72));
                painter.circle_filled(
                    pos + Vec2::new(-radius * 0.16, -radius * 0.16),
                    radius * 0.78,
                    color,
                );
                painter.circle_filled(
                    pos + Vec2::new(-radius * 0.28, -radius * 0.28),
                    radius * 0.42,
                    mix_white(color, 0.45),
                );
                painter.circle_filled(
                    pos + Vec2::new(-radius * 0.34, -radius * 0.34),
                    radius * 0.18,
                    mix_white(color, 0.85),
                );
                painter.circle_stroke(
                    pos,
                    radius,
                    Stroke::new(
                        if hovered {
                            2.2_f32
                        } else if is_hydrogen {
                            0.75_f32
                        } else {
                            0.9_f32
                        },
                        if hovered {
                            Color32::from_rgb(34, 139, 230)
                        } else if is_hydrogen {
                            Color32::from_rgba_unmultiplied(80, 90, 108, 130)
                        } else {
                            Color32::from_rgba_unmultiplied(15, 18, 25, 120)
                        },
                    ),
                );
                if let Some(lbl) = label {
                    let font = egui::FontId::proportional((radius * 1.02).clamp(9.0, 13.5));
                    // Subtle dark shadow for crisp readability on any atom sphere color
                    painter.text(
                        pos + Vec2::new(0.8, 0.8),
                        egui::Align2::CENTER_CENTER,
                        &lbl,
                        font.clone(),
                        Color32::from_rgba_unmultiplied(10, 12, 18, 190),
                    );
                    painter.text(pos, egui::Align2::CENTER_CENTER, lbl, font, Color32::WHITE);
                }
            }
        }
    }

    // Corner HUD hints + hover atom inspector badge
    let hud_text = if let Some((hi, _)) = hovered_atom {
        let a = &mol.atoms[hi];
        if let Some(meta) = &a.pdb_meta {
            format!(
                "Atom #{}: {} ({} {} {})",
                hi + 1,
                a.symbol,
                meta.atom_name,
                meta.res_name,
                meta.res_seq
            )
        } else {
            format!(
                "Atom #{}: {}{}",
                hi + 1,
                a.symbol,
                a.cip_label.map_or(String::new(), |c| format!(" ({c})"))
            )
        }
    } else {
        let h_note = if !hydrogens.is_empty() { " | +H" } else { "" };
        format!("Drag: Orbit | Scroll: Zoom | Double-Click: Reset | 3D Conformer{h_note}")
    };
    painter.text(
        rect.left_bottom() + Vec2::new(8.0, -6.0),
        egui::Align2::LEFT_BOTTOM,
        hud_text,
        egui::FontId::proportional(11.0),
        Color32::from_gray(if dark_mode { 170 } else { 100 }),
    );

    // Bottom-right 3D axis orientation triad (X red, Y green, Z blue)
    let origin = rect.right_bottom() + Vec2::new(-28.0, -26.0);
    let axes = [
        ("X", [1.0_f32, 0.0, 0.0], Color32::from_rgb(220, 72, 72)),
        ("Y", [0.0, 1.0, 0.0], Color32::from_rgb(64, 185, 86)),
        ("Z", [0.0, 0.0, 1.0], Color32::from_rgb(66, 138, 240)),
    ];
    for (lbl, v, col) in axes {
        let rx = v[0] * cy + v[2] * sy;
        let rz1 = -v[0] * sy + v[2] * cy;
        let ry = v[1] * cp - rz1 * sp;
        let end = origin + Vec2::new(rx * 16.0, -ry * 16.0);
        painter.line_segment([origin, end], Stroke::new(1.8_f32, col));
        painter.text(
            end + Vec2::new(rx * 4.5, -ry * 4.5),
            egui::Align2::CENTER_CENTER,
            lbl,
            egui::FontId::proportional(9.5),
            col,
        );
    }

    resp
}
