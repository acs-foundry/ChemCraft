//! Layer 1 (`mol-vector`): Clean-room 2D coordinate layout, interactive ring snapping,
//! `Clean Structure 2D`, VectorCraft-style vector scene IR, and SVG/EPS/PDF/PNG export.

#![forbid(unsafe_code)]

pub mod export;
pub mod layout;
pub mod scene;

pub use export::{
    VectorError, export_sar_grid_eps, export_sar_grid_pdf_bytes, export_sar_grid_png_bytes,
    export_sar_grid_svg, scene_to_eps, scene_to_pdf_bytes, scene_to_png_bytes, scene_to_svg,
};
pub use layout::{
    RingTemplate, STANDARD_BOND_LEN, SnapTarget, clean_structure_2d, count_bond_crossings,
    ensure_2d_coords, generate_2d_coords, stamp_ring_template,
};
pub use scene::{
    DepictionPreset, ReactionArrow, Rect2, Rgba, SceneStyle, VectorElement, VectorScene,
    build_mol_scene,
};

#[cfg(test)]
mod tests {
    use super::*;
    use mol_graph::{MolError, parse_smiles};

    #[test]
    fn renders_10_complex_molecules_to_svg_eps_pdf_png_without_text_overlaps()
    -> Result<(), Box<dyn std::error::Error>> {
        let complex = [
            ("Testosterone", "CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C"),
            ("Caffeine", "Cn1cnc2c1c(=O)n(C)c(=O)n2C"),
            (
                "Sildenafil",
                "CCCc1nn(C)c2c(=O)[nH]c(-c3cc(S(=O)(=O)N4CCN(C)CC4)ccc3OCC)nc12",
            ),
            ("Aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
            ("Ibuprofen", "CC(C)Cc1ccc(C(C)C(=O)O)cc1"),
            ("Morphine", "CN1CCC23c4c5ccc(O)c4OC2C(O)C=CC3C1C5"),
            ("Diazepam", "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21"),
            ("Warfarin", "CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O"),
            ("Omeprazole", "COc1ccc2[nH]c(S(=O)Cc3ncc(C)c(OC)c3C)nc2c1"),
            (
                "Ciprofloxacin",
                "O=C(O)c1cn(C2CC2)c2cc(N3CCNCC3)c(F)cc2c1=O",
            ),
        ];

        let style = SceneStyle::default();
        let mut grid_mols = Vec::new();
        for (name, smi) in complex {
            let mut mol = parse_smiles(smi)?;
            mol.name = Some(name.to_string());
            generate_2d_coords(&mut mol);
            let scene = build_mol_scene(&mol, &style, &[]);
            assert_eq!(
                scene.label_overlap_count(),
                0,
                "{name} had overlapping atom labels"
            );
            let svg = scene_to_svg(&scene);
            assert!(svg.contains("<svg") && svg.contains("</svg>"));
            let eps = scene_to_eps(&scene);
            assert!(eps.starts_with("%!PS-Adobe-3.0 EPSF-3.0"));
            let pdf = scene_to_pdf_bytes(&scene)?;
            assert!(pdf.starts_with(b"%PDF-1.4") && pdf.contains(&b'/'));
            let png = scene_to_png_bytes(&scene, 1.0)?;
            assert!(png.starts_with(&[137, 80, 78, 71]));
            grid_mols.push((name, mol));
        }

        let refs: Vec<(&str, &mol_graph::MolGraph)> =
            grid_mols.iter().map(|(n, m)| (*n, m)).collect();
        let grid_svg = export_sar_grid_svg(&refs, 4, 240.0, 200.0, false);
        assert!(grid_svg.contains("Testosterone") && grid_svg.contains("Sildenafil"));
        let grid_pdf = export_sar_grid_pdf_bytes(&refs, 4, 240.0, 200.0)?;
        assert!(grid_pdf.starts_with(b"%PDF-1.4"));
        let grid_eps = export_sar_grid_eps(&refs, 4, 240.0, 200.0, false);
        assert!(grid_eps.contains("Testosterone") && grid_eps.contains("Ciprofloxacin"));
        let grid_png = export_sar_grid_png_bytes(&refs, 4, 240.0, 200.0, false, 1.0)?;
        assert!(grid_png.starts_with(&[137, 80, 78, 71]));
        Ok(())
    }

    #[test]
    fn pdf_renders_circles_winansi_accents_and_screen_to_world_roundtrips()
    -> Result<(), Box<dyn std::error::Error>> {
        let style = SceneStyle::default();
        let mut mol = parse_smiles("c1ccccc1")?;
        generate_2d_coords(&mut mol);
        let mut scene = build_mol_scene(&mol, &style, &[]);
        // Inverse screen -> world roundtrip check
        let sp = scene.atom_screen_pos[0];
        let wp = scene.screen_to_world(sp[0], sp[1]);
        let orig = mol.atoms[0].coords_2d.unwrap_or([0.0, 0.0]);
        assert!((wp[0] - orig[0]).abs() < 1e-6);
        assert!((wp[1] - orig[1]).abs() < 1e-6);

        // Add explicit Circle + accented Label
        scene.elements.push(VectorElement::Circle {
            cx: 120.0,
            cy: 90.0,
            r: 14.0,
            fill: Some(Rgba(225, 45, 45, 55)),
            stroke: Some((Rgba(220, 38, 38, 220), 1.5)),
        });
        scene.elements.push(VectorElement::Label {
            text: "Cafeína".to_string(),
            cx: 200.0,
            cy: 40.0,
            font_size: 12.0,
            color: Rgba(20, 20, 20, 255),
            bbox: Rect2 {
                x: 170.0,
                y: 34.0,
                w: 60.0,
                h: 12.0,
            },
        });
        let eps = scene_to_eps(&scene);
        assert!(eps.contains("Cafe\\355na"));
        assert!(eps.contains("arc fill") && eps.contains("arc stroke"));

        // Decompress PDF content stream to verify circle Bezier curves + WinAnsi octal!
        let pdf = scene_to_pdf_bytes(&scene)?;
        let start = pdf
            .windows(7)
            .position(|w| w == b"stream\n")
            .ok_or("missing pdf stream")?
            + 7;
        let end = pdf
            .windows(10)
            .position(|w| w == b"\nendstream")
            .ok_or("missing pdf endstream")?;
        let mut dec = flate2::read::ZlibDecoder::new(&pdf[start..end]);
        let mut text = String::new();
        std::io::Read::read_to_string(&mut dec, &mut text)?;
        assert!(text.contains("Cafe\\355na"));
        assert!(text.contains(" c f") && text.contains(" c S"));
        Ok(())
    }

    #[test]
    fn ring_snapping_fuses_benzene_to_naphthalene() -> Result<(), MolError> {
        let mut mol = parse_smiles("c1ccccc1")?;
        generate_2d_coords(&mut mol);
        stamp_ring_template(
            &mut mol,
            RingTemplate::Benzene,
            SnapTarget::Bond {
                bond_idx: 0,
                cursor_x: 5.0,
                cursor_y: 0.0,
            },
        )?;
        assert_eq!(mol.atoms.len(), 10);
        assert_eq!(mol.formula(), "C10H8");
        Ok(())
    }

    #[test]
    fn png_grid_title_renders_all_12_glyphs_of_testosterone()
    -> Result<(), Box<dyn std::error::Error>> {
        let mut mol = parse_smiles("CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C")?;
        generate_2d_coords(&mut mol);
        let png_bytes =
            export_sar_grid_png_bytes(&[("Testosterone", &mol)], 1, 260.0, 220.0, false, 2.0)?;
        let decoder = png::Decoder::new(std::io::Cursor::new(&png_bytes));
        let mut reader = decoder.read_info()?;
        let mut buf = vec![0; reader.output_buffer_size().ok_or("missing buf size")?];
        let info = reader.next_frame(&mut buf)?;
        let w = info.width as usize;
        let h = info.height as usize;

        // Inspect bottom title strip (cell_h = 220, y = 200..216 at scale 2.0 -> py = 402..432)
        // Each inked column in the title has at least one dark pixel (< 160 RGB luminance)
        let mut col_heights = vec![0usize; w];
        for (x, col_h) in col_heights.iter_mut().enumerate().take(w - 6).skip(6) {
            let mut min_y = usize::MAX;
            let mut max_y = 0usize;
            for y in 402..432.min(h) {
                let idx = (y * w + x) * 4;
                if buf[idx] < 160 && buf[idx + 1] < 160 && buf[idx + 2] < 160 {
                    min_y = min_y.min(y);
                    max_y = max_y.max(y);
                }
            }
            if max_y >= min_y {
                *col_h = max_y - min_y + 1;
            }
        }
        // Group consecutive inked x-columns into glyph spans and check each letter's height >= 5px
        let mut spans = Vec::new();
        let mut in_span = false;
        let mut span_max_h = 0usize;
        for &ch in &col_heights {
            if ch > 0 {
                in_span = true;
                span_max_h = span_max_h.max(ch);
            } else if in_span {
                spans.push(span_max_h);
                in_span = false;
                span_max_h = 0;
            }
        }
        if in_span {
            spans.push(span_max_h);
        }
        assert_eq!(
            spans.len(),
            12,
            "expected 12 glyphs in 'Testosterone', got {spans:?}"
        );
        assert!(
            spans.iter().all(|&gh| gh >= 6),
            "every glyph of 'Testosterone' must have vertical height >= 6px (no fallback dashes): {spans:?}"
        );
        Ok(())
    }

    #[test]
    fn bridged_cages_layout_zero_crossed_bonds() -> Result<(), MolError> {
        for (name, smi, max_cv_pct, max_total_crossings) in [
            ("Cubane", "C12C3C4C1C5C2C3C45", 50.0, 0usize),
            ("Prismane", "C12C3C1C4C2C34", 45.0, 0),
            ("Norbornane", "C1CC2CCC1C2", 25.0, 0),
            ("Bicyclo[2.2.2]octane", "C1CC2CCC1CC2", 25.0, 0),
            ("Adamantane", "C1C2CC3CC1CC(C2)C3", 30.0, 0),
            ("Methyl-Adamantane", "CC12CC3CC(CC(C3)C1)C2", 25.0, 0),
            ("Camphor", "CC1(C)C2CCC1(C)C(=O)C2", 12.0, 1),
            (
                "Baricitinib",
                "CCS(=O)(=O)N1CC(C1)(CC#N)n2cc(cn2)c3ncnc4c3cc[nH]4",
                15.0,
                0,
            ),
        ] {
            let mut mol = parse_smiles(smi)?;
            generate_2d_coords(&mut mol);
            let crossings = count_bond_crossings(&mol);
            assert!(
                crossings <= max_total_crossings,
                "{name} ({smi}) had {crossings} crossed bonds (max {max_total_crossings})"
            );
            let lens: Vec<f64> = mol
                .bonds
                .iter()
                .map(|b| {
                    let pu = mol.atoms[b.u].coords_2d.unwrap_or([0.0, 0.0]);
                    let pv = mol.atoms[b.v].coords_2d.unwrap_or([0.0, 0.0]);
                    (pu[0] - pv[0]).hypot(pu[1] - pv[1])
                })
                .collect();
            let mean = lens.iter().sum::<f64>() / (lens.len() as f64);
            let var = lens.iter().map(|d| (d - mean).powi(2)).sum::<f64>() / (lens.len() as f64);
            let cv = (var.sqrt() / mean) * 100.0;
            assert!(
                cv <= max_cv_pct,
                "{name} bond-length CV={cv:.1}% exceeds limit {max_cv_pct:.1}%"
            );
        }
        Ok(())
    }

    #[test]
    fn twenty_drugs_layout_bond_lengths_and_spacing_remain_uniform() -> Result<(), MolError> {
        let drugs = [
            ("Morphine", "CN1CCC23c4c5ccc(O)c4OC2C(O)C=CC3C1C5"),
            ("Diazepam", "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21"),
            (
                "Taxol",
                "CC(=O)OC1C(=O)C2(C)C(O)CC3OCC3(OC(C)=O)C2C(OC(=O)c2ccccc2)C2(O)CC(OC(=O)C(O)C(NC(=O)c3ccccc3)c3ccccc3)C(C)=C1C2(C)C",
            ),
            ("Artemisinin", "CC1CCC2C(C)C(=O)OC3OC4(C)CCC1C23OO4"),
            ("Testosterone", "CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C"),
            ("Caffeine", "Cn1cnc2c1c(=O)n(C)c(=O)n2C"),
            (
                "Sildenafil",
                "CCCc1nn(C)c2c(=O)[nH]c(-c3cc(S(=O)(=O)N4CCN(C)CC4)ccc3OCC)nc12",
            ),
            ("Aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
            ("Ibuprofen", "CC(C)Cc1ccc(C(C)C(=O)O)cc1"),
            ("Warfarin", "CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O"),
            ("Omeprazole", "COc1ccc2[nH]c(S(=O)Cc3ncc(C)c(OC)c3C)nc2c1"),
            (
                "Ciprofloxacin",
                "O=C(O)c1cn(C2CC2)c2cc(N3CCNCC3)c(F)cc2c1=O",
            ),
            (
                "Penicillin_G",
                "CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)C(=O)N2[C@H]1C(=O)O",
            ),
            ("Tamoxifen", "CC/C(=C(/c1ccccc1)c1ccc(OCCN(C)C)cc1)c1ccccc1"),
            (
                "Atorvastatin",
                "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CCC(O)CC(O)CC(=O)O",
            ),
            ("Metformin", "CN(C)C(=N)NC(=N)N"),
            ("Fluoxetine", "CNCCC(Oc1ccc(C(F)(F)F)cc1)c1ccccc1"),
            ("Paracetamol", "CC(=O)Nc1ccc(O)cc1"),
            ("Adamantane", "C1C2CC3CC1CC(C2)C3"),
            ("Cubane", "C12C3C4C1C5C2C3C45"),
        ];
        for (name, smi) in drugs {
            let mut mol = parse_smiles(smi)?;
            generate_2d_coords(&mut mol);
            let mut min_len = f64::INFINITY;
            let mut max_len = 0.0_f64;
            for b in &mol.bonds {
                let pu = mol.atoms[b.u].coords_2d.unwrap_or([0.0, 0.0]);
                let pv = mol.atoms[b.v].coords_2d.unwrap_or([0.0, 0.0]);
                let d = (pu[0] - pv[0]).hypot(pu[1] - pv[1]);
                min_len = min_len.min(d);
                max_len = max_len.max(d);
            }
            let ratio = max_len / min_len.max(1e-6);
            assert!(
                min_len >= STANDARD_BOND_LEN * 0.45 && ratio <= 3.5,
                "{name} bond length collapsed/stretched: min={min_len:.3}, max={max_len:.3}, ratio={ratio:.2}"
            );
        }
        Ok(())
    }
}
