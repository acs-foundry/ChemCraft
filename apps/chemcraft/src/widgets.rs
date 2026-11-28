use std::f32::consts::TAU;

use eframe::egui;
use egui::{Color32, Pos2, Rect, Sense, Stroke, Vec2};
use mol_admet_bridge::MoleculeEvaluation;
use mol_sar::SarRow;
use mol_vector::{ReactionArrow, Rgba, SceneStyle, VectorElement, VectorScene, build_mol_scene};

/// Paint any `VectorScene` directly onto an `egui::Painter` inside `target_rect`.
/// Ensures a single unified 2D molecular drawing pipeline shared across thumbnails,
/// tooltips, and SVG/EPS/PDF/PNG exports.
pub fn paint_vector_scene(painter: &egui::Painter, scene: &VectorScene, origin: Pos2) {
    for el in &scene.elements {
        match el {
            VectorElement::Line {
                x1,
                y1,
                x2,
                y2,
                width,
                color,
                dashed,
            } => {
                let p1 = Pos2::new(origin.x + *x1 as f32, origin.y + *y1 as f32);
                let p2 = Pos2::new(origin.x + *x2 as f32, origin.y + *y2 as f32);
                let stroke = Stroke::new(*width as f32, to_color32(*color));
                if *dashed {
                    paint_dashed_line(painter, p1, p2, stroke);
                } else {
                    painter.line_segment([p1, p2], stroke);
                }
            }
            VectorElement::Polyline {
                points,
                width,
                color,
                closed,
            } => {
                let pts: Vec<Pos2> = points
                    .iter()
                    .map(|p| Pos2::new(origin.x + p[0] as f32, origin.y + p[1] as f32))
                    .collect();
                let st = Stroke::new(*width as f32, to_color32(*color));
                if *closed {
                    painter.add(egui::Shape::closed_line(pts, st));
                } else {
                    painter.add(egui::Shape::line(pts, st));
                }
            }
            VectorElement::Polygon { points, fill } => {
                let pts: Vec<Pos2> = points
                    .iter()
                    .map(|p| Pos2::new(origin.x + p[0] as f32, origin.y + p[1] as f32))
                    .collect();
                painter.add(egui::Shape::convex_polygon(
                    pts,
                    to_color32(*fill),
                    Stroke::NONE,
                ));
            }
            VectorElement::Circle {
                cx,
                cy,
                r,
                fill,
                stroke,
            } => {
                let c = Pos2::new(origin.x + *cx as f32, origin.y + *cy as f32);
                if let Some(f) = fill {
                    painter.circle_filled(c, *r as f32, to_color32(*f));
                }
                if let Some((col, w)) = stroke {
                    painter.circle_stroke(c, *r as f32, Stroke::new(*w as f32, to_color32(*col)));
                }
            }
            VectorElement::CubicCurve {
                p0,
                p1,
                p2,
                p3,
                width,
                color,
            } => {
                let bezier = egui::epaint::CubicBezierShape::from_points_stroke(
                    [
                        Pos2::new(origin.x + p0[0] as f32, origin.y + p0[1] as f32),
                        Pos2::new(origin.x + p1[0] as f32, origin.y + p1[1] as f32),
                        Pos2::new(origin.x + p2[0] as f32, origin.y + p2[1] as f32),
                        Pos2::new(origin.x + p3[0] as f32, origin.y + p3[1] as f32),
                    ],
                    false,
                    Color32::TRANSPARENT,
                    Stroke::new(*width as f32, to_color32(*color)),
                );
                painter.add(bezier);
            }
            VectorElement::Label {
                text,
                cx,
                cy,
                font_size,
                color,
                ..
            } => {
                painter.text(
                    Pos2::new(origin.x + *cx as f32, origin.y + *cy as f32),
                    egui::Align2::CENTER_CENTER,
                    text,
                    egui::FontId::proportional(*font_size as f32),
                    to_color32(*color),
                );
            }
        }
    }
}

/// Paint a compact 2D structure thumbnail for SAR table cells or hover tooltips.
pub fn paint_mol_thumbnail(
    ui: &mut egui::Ui,
    mol: &mol_graph::MolGraph,
    arrows: &[ReactionArrow],
    size: Vec2,
    dark_mode: bool,
) -> VectorScene {
    let (rect, _) = ui.allocate_exact_size(size, Sense::hover());
    let bg = if dark_mode {
        Color32::from_rgb(24, 26, 32)
    } else {
        Color32::from_rgb(252, 252, 254)
    };
    ui.painter().rect_filled(rect, 4.0, bg);
    let style = SceneStyle {
        width: f64::from(rect.width()),
        height: f64::from(rect.height()),
        padding: 10.0,
        stroke_width: 1.35,
        font_size: 10.0,
        dark_mode,
        show_stereo_labels: true,
        ..SceneStyle::default()
    };
    let mut scene = build_mol_scene(mol, &style, arrows);
    scene.background = None;
    paint_vector_scene(&ui.painter_at(rect), &scene, rect.min);
    scene
}

/// Draw the 6-axis AdmetCraft Bioavailability Radar (`LIPO, SIZE, POLAR, INSOLU, INSATU, FLEX`).
pub fn draw_admet_radar(ui: &mut egui::Ui, eval: &MoleculeEvaluation, dark_mode: bool) {
    let size = Vec2::new(220.0, 175.0);
    let (rect, _) = ui.allocate_exact_size(size, Sense::hover());
    let painter = ui.painter_at(rect);
    let center = rect.center();
    let max_r = 56.0_f32;

    let labels = ["LIPO", "SIZE", "POLAR", "INSOLU", "INSATU", "FLEX"];
    let grid_col = if dark_mode {
        Color32::from_gray(75)
    } else {
        Color32::from_gray(205)
    };
    let text_col = if dark_mode {
        Color32::from_gray(225)
    } else {
        Color32::from_gray(45)
    };

    for ring in [0.25_f32, 0.50, 0.75, 1.0] {
        let pts: Vec<Pos2> = (0..6)
            .map(|k| {
                let ang = -std::f32::consts::FRAC_PI_2 + (k as f32) * TAU / 6.0;
                Pos2::new(
                    center.x + max_r * ring * ang.cos(),
                    center.y + max_r * ring * ang.sin(),
                )
            })
            .collect();
        painter.add(egui::Shape::closed_line(
            pts,
            Stroke::new(1.0_f32, grid_col),
        ));
    }

    // Optimal bioavailability zone polygon
    let opt_pts: Vec<Pos2> = (0..6)
        .map(|k| {
            let r = max_r * 0.68;
            let ang = -std::f32::consts::FRAC_PI_2 + (k as f32) * TAU / 6.0;
            Pos2::new(center.x + r * ang.cos(), center.y + r * ang.sin())
        })
        .collect();
    painter.add(egui::Shape::convex_polygon(
        opt_pts,
        Color32::from_rgba_unmultiplied(230, 80, 95, 36),
        Stroke::NONE,
    ));

    let axes = eval.bioavailability_radar.normalized_axes;
    let mol_pts: Vec<Pos2> = (0..6)
        .map(|k| {
            let v = (axes[k] as f32).clamp(0.05, 1.0);
            let ang = -std::f32::consts::FRAC_PI_2 + (k as f32) * TAU / 6.0;
            Pos2::new(
                center.x + max_r * v * ang.cos(),
                center.y + max_r * v * ang.sin(),
            )
        })
        .collect();
    painter.add(egui::Shape::convex_polygon(
        mol_pts.clone(),
        Color32::from_rgba_unmultiplied(224, 49, 49, 75),
        Stroke::new(2.0_f32, Color32::from_rgb(215, 42, 42)),
    ));
    for p in mol_pts {
        painter.circle_filled(p, 2.8, Color32::from_rgb(200, 30, 30));
    }

    for (k, lbl) in labels.iter().enumerate() {
        let ang = -std::f32::consts::FRAC_PI_2 + (k as f32) * TAU / 6.0;
        let lp = Pos2::new(
            center.x + (max_r + 17.0) * ang.cos(),
            center.y + (max_r + 15.0) * ang.sin(),
        );
        painter.text(
            lp,
            egui::Align2::CENTER_CENTER,
            *lbl,
            egui::FontId::proportional(9.5),
            text_col,
        );
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScatterXAxis {
    ConsensusLogP,
    Tpsa,
    MolecularWeight,
    LigandEfficiency,
}

impl ScatterXAxis {
    pub fn label(self) -> &'static str {
        match self {
            Self::ConsensusLogP => "Consensus LogP (AdmetCraft)",
            Self::Tpsa => "TPSA (A^2)",
            Self::MolecularWeight => "MW (g/mol)",
            Self::LigandEfficiency => "Ligand Efficiency (LE)",
        }
    }

    fn extract(self, row: &SarRow) -> Option<f64> {
        match self {
            Self::ConsensusLogP => row.admet.as_ref().map(|a| a.consensus_logp),
            Self::Tpsa => row.admet.as_ref().map(|a| a.tpsa),
            Self::MolecularWeight => Some(row.mw),
            Self::LigandEfficiency => row.ligand_efficiency(),
        }
    }
}

/// Interactive SAR Scatter Plot (`pIC50` vs selectable X descriptor) with 2D structure tooltips.
pub fn draw_sar_scatter(
    ui: &mut egui::Ui,
    rows: &[SarRow],
    visible_indices: &[usize],
    x_axis: ScatterXAxis,
    selected_idx: Option<usize>,
    dark_mode: bool,
) -> Option<usize> {
    let avail = ui.available_size();
    let size = Vec2::new(avail.x.max(260.0), avail.y.max(170.0));
    let (rect, response) = ui.allocate_exact_size(size, Sense::click_and_drag());
    let painter = ui.painter_at(rect);

    let bg = if dark_mode {
        Color32::from_rgb(24, 27, 34)
    } else {
        Color32::from_rgb(250, 251, 253)
    };
    let fg = if dark_mode {
        Color32::from_gray(220)
    } else {
        Color32::from_gray(40)
    };
    painter.rect_filled(rect, 4.0, bg);

    let mut pts: Vec<(usize, f64, f64)> = Vec::new();
    for &idx in visible_indices {
        let r = &rows[idx];
        if let (Some(x), Some(y)) = (x_axis.extract(r), r.pic50) {
            pts.push((idx, x, y));
        }
    }
    if pts.is_empty() {
        painter.text(
            rect.center(),
            egui::Align2::CENTER_CENTER,
            "No points with both pIC50 and ADMET descriptor (click 'Load SAR Demo')",
            egui::FontId::proportional(12.0),
            fg,
        );
        return None;
    }

    let (mut min_x, mut max_x) = (f64::INFINITY, f64::NEG_INFINITY);
    let (mut min_y, mut max_y) = (f64::INFINITY, f64::NEG_INFINITY);
    for &(_, x, y) in &pts {
        min_x = min_x.min(x);
        max_x = max_x.max(x);
        min_y = min_y.min(y);
        max_y = max_y.max(y);
    }
    let span_x = (max_x - min_x).max(0.5);
    let span_y = (max_y - min_y).max(0.5);
    min_x -= span_x * 0.08;
    max_x += span_x * 0.08;
    min_y -= span_y * 0.1;
    max_y += span_y * 0.1;

    let pad_l = 48.0_f32;
    let pad_b = 28.0_f32;
    let pad_r = 16.0_f32;
    let pad_t = 14.0_f32;
    let plot = Rect::from_min_max(
        Pos2::new(rect.min.x + pad_l, rect.min.y + pad_t),
        Pos2::new(rect.max.x - pad_r, rect.max.y - pad_b),
    );
    painter.rect_stroke(
        plot,
        0.0,
        Stroke::new(1.0_f32, Color32::from_gray(130)),
        egui::StrokeKind::Inside,
    );

    let map_x = |x: f64| plot.min.x + ((x - min_x) / (max_x - min_x)) as f32 * plot.width();
    let map_y = |y: f64| plot.max.y - ((y - min_y) / (max_y - min_y)) as f32 * plot.height();

    painter.text(
        Pos2::new(plot.center().x, rect.max.y - 11.0),
        egui::Align2::CENTER_CENTER,
        x_axis.label(),
        egui::FontId::proportional(10.5),
        fg,
    );
    painter.text(
        Pos2::new(rect.min.x + 22.0, plot.center().y),
        egui::Align2::CENTER_CENTER,
        "pIC50",
        egui::FontId::proportional(10.5),
        fg,
    );

    let hover = response.hover_pos();
    let mut hovered_idx = None;
    let mut clicked_idx = None;

    for &(idx, x, y) in &pts {
        let p = Pos2::new(map_x(x), map_y(y));
        let r = &rows[idx];
        let pass_lipinski = r.admet.as_ref().is_some_and(|a| a.lipinski_violations == 0);
        let col = if pass_lipinski {
            Color32::from_rgb(34, 139, 230)
        } else {
            Color32::from_rgb(240, 140, 0)
        };
        let is_sel = selected_idx == Some(idx);
        let rad = if is_sel { 6.5 } else { 4.8 };
        painter.circle_filled(p, rad, col);
        painter.circle_stroke(
            p,
            rad,
            Stroke::new(if is_sel { 2.0_f32 } else { 1.0_f32 }, Color32::BLACK),
        );

        if let Some(hp) = hover
            && hp.distance(p) <= 9.0
        {
            hovered_idx = Some(idx);
            if response.clicked() {
                clicked_idx = Some(idx);
            }
        }
    }

    if let Some(hi) = hovered_idx {
        let r = &rows[hi];
        egui::Tooltip::for_widget(&response)
            .at_pointer()
            .show(|ui| {
                ui.set_min_width(180.0);
                ui.strong(&r.name);
                ui.label(format!(
                    "pIC50: {}  |  MW: {:.1}",
                    r.pic50.map_or("-".to_string(), |v| format!("{v:.2}")),
                    r.mw
                ));
                if let Some(ref adm) = r.admet {
                    ui.label(format!(
                        "LogP: {:.2} | TPSA: {:.1} | GI: {}",
                        adm.consensus_logp, adm.tpsa, adm.gi_absorption
                    ));
                }
                let _ = paint_mol_thumbnail(ui, &r.mol, &[], Vec2::new(170.0, 110.0), dark_mode);
            });
    }

    clicked_idx
}

fn paint_dashed_line(painter: &egui::Painter, p1: Pos2, p2: Pos2, stroke: Stroke) {
    let delta = p2 - p1;
    let len = delta.length();
    if len < 2.0 {
        return;
    }
    let steps = ((len / 6.0) as usize).max(2);
    for s in (0..steps).step_by(2) {
        let t0 = s as f32 / steps as f32;
        let t1 = ((s + 1) as f32 / steps as f32).min(1.0);
        painter.line_segment([p1 + delta * t0, p1 + delta * t1], stroke);
    }
}

fn to_color32(c: Rgba) -> Color32 {
    Color32::from_rgba_unmultiplied(c.0, c.1, c.2, c.3)
}
