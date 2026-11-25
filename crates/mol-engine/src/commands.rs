use std::path::{Path, PathBuf};
use std::sync::Arc;

use serde::{Deserialize, Serialize};
use thiserror::Error;

use mol_admet_bridge::{
    AdmetCraftProvider, AdmetDelta, AdmetProvider, MoleculeEvaluation, compute_admet_delta,
};
use mol_graph::{
    BondOrder, BondStereo, MolAtom, MolGraph, attach_superatom_group, mutate_molecule, parse_cdxml,
    parse_smiles, resolve_chemical_name,
};
use mol_sar::{
    ActivityCliff, MatchedMolecularPair, SarFilter, SarTable, find_activity_cliffs,
    find_matched_molecular_pairs,
};
use mol_vector::{
    ReactionArrow, RingTemplate, SnapTarget, clean_structure_2d, generate_2d_coords,
    stamp_ring_template,
};

/// Copy-on-Write snapshot for fast structural Undo/Redo (`Arc`-shared when untouched).
#[derive(Debug, Clone)]
pub struct SessionSnapshot {
    pub canvas_mol: Arc<MolGraph>,
    pub reaction_arrows: Arc<Vec<ReactionArrow>>,
    pub sar_table: Arc<SarTable>,
    pub selected_sar_idx: Option<usize>,
}

/// Interactive & scriptable `ChemCraft` session state.
#[derive(Debug, Clone)]
pub struct ChemSession {
    pub canvas_mol: Arc<MolGraph>,
    pub reaction_arrows: Arc<Vec<ReactionArrow>>,
    pub sar_table: Arc<SarTable>,
    pub sar_filter: SarFilter,
    pub selected_sar_idx: Option<usize>,
    pub live_admet: Option<MoleculeEvaluation>,
    pub last_mutation_delta: Option<AdmetDelta>,
    pub cliffs: Vec<ActivityCliff>,
    pub mmps: Vec<MatchedMolecularPair>,
    pub prefer_external_admet: bool,
    /// When `false` (default), ADMET is strictly on-demand instead of auto-computing on every edit.
    pub auto_refresh_admet: bool,
    undo_stack: Vec<SessionSnapshot>,
    redo_stack: Vec<SessionSnapshot>,
}

impl Default for ChemSession {
    fn default() -> Self {
        Self::new()
    }
}

impl ChemSession {
    pub fn new() -> Self {
        Self {
            canvas_mol: Arc::new(MolGraph::new()),
            reaction_arrows: Arc::new(Vec::new()),
            sar_table: Arc::new(SarTable::new()),
            sar_filter: SarFilter::default(),
            selected_sar_idx: None,
            live_admet: None,
            last_mutation_delta: None,
            cliffs: Vec::new(),
            mmps: Vec::new(),
            prefer_external_admet: false,
            auto_refresh_admet: false,
            undo_stack: Vec::new(),
            redo_stack: Vec::new(),
        }
    }

    pub fn can_undo(&self) -> bool {
        !self.undo_stack.is_empty()
    }

    pub fn can_redo(&self) -> bool {
        !self.redo_stack.is_empty()
    }

    pub fn push_undo_checkpoint(&mut self) {
        self.push_undo();
    }

    fn push_undo(&mut self) {
        self.undo_stack.push(SessionSnapshot {
            canvas_mol: Arc::clone(&self.canvas_mol),
            reaction_arrows: Arc::clone(&self.reaction_arrows),
            sar_table: Arc::clone(&self.sar_table),
            selected_sar_idx: self.selected_sar_idx,
        });
        if self.undo_stack.len() > 64 {
            self.undo_stack.remove(0);
        }
        self.redo_stack.clear();
    }

    fn restore_snapshot(&mut self, snap: SessionSnapshot) {
        self.canvas_mol = snap.canvas_mol;
        self.reaction_arrows = snap.reaction_arrows;
        self.sar_table = snap.sar_table;
        self.selected_sar_idx = snap.selected_sar_idx;
        self.on_canvas_changed();
    }

    /// Clear cached ADMET after a canvas change unless live auto-refresh was explicitly opted into.
    pub fn on_canvas_changed(&mut self) {
        if self.auto_refresh_admet {
            let _ = self.refresh_live_admet();
        } else {
            self.live_admet = None;
        }
    }

    pub fn admet_provider(&self) -> AdmetCraftProvider {
        AdmetCraftProvider {
            prefer_external_cli: self.prefer_external_admet,
        }
    }

    pub fn refresh_live_admet(&mut self) -> Result<(), String> {
        if self.canvas_mol.atoms.is_empty() {
            self.live_admet = None;
            return Ok(());
        }
        let smi = self.canvas_mol.canonical_smiles();
        if smi.is_empty() {
            self.live_admet = None;
            return Ok(());
        }
        let eval = self
            .admet_provider()
            .evaluate_smiles(&smi, self.canvas_mol.name.as_deref())
            .map_err(|e| e.to_string())?;
        self.live_admet = Some(eval);
        Ok(())
    }

    pub fn recompute_sar_analytics(&mut self) {
        let mut table = (*self.sar_table).clone();
        table.enrich_with_admet(&self.admet_provider());
        self.cliffs = find_activity_cliffs(&table.rows, 0.40);
        self.mmps = find_matched_molecular_pairs(&table.rows);
        self.sar_table = Arc::new(table);
    }
}

/// Deterministic command enum shared across GUI palette, CLI, and MCP.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum Command {
    LoadCanvasSmiles {
        smiles: String,
        name: Option<String>,
    },
    LoadChemicalName {
        name: String,
    },
    LoadCdxml {
        xml: String,
    },
    ClearCanvas,
    CleanStructure2D,
    StampRing {
        template_name: String,
        x: f64,
        y: f64,
        snap_atom: Option<usize>,
        snap_bond: Option<usize>,
    },
    AddAtomAndBond {
        from_atom: Option<usize>,
        symbol: String,
        x: f64,
        y: f64,
        order: BondOrder,
        stereo: BondStereo,
    },
    ConnectAtoms {
        u: usize,
        v: usize,
        order: BondOrder,
        stereo: BondStereo,
    },
    SetAtomElement {
        atom_idx: usize,
        symbol: String,
    },
    AdjustAtomCharge {
        atom_idx: usize,
        delta: i8,
    },
    CycleBondOrder {
        bond_idx: usize,
    },
    SetBondStereo {
        bond_idx: usize,
        stereo: BondStereo,
    },
    AttachSuperatomGroup {
        atom_idx: Option<usize>,
        abbr: String,
    },
    FlipCanvas2D {
        horizontal: bool,
    },
    RotateCanvas2D {
        degrees: f64,
    },
    EraseAtom {
        atom_idx: usize,
    },
    EraseBond {
        bond_idx: usize,
    },
    MutateBioisostere {
        atom_index: Option<usize>,
        group: String,
    },
    AddCanvasToSar {
        pic50: Option<f64>,
    },
    LoadDemo50Sar,
    FilterSarByCanvasSubstructure,
    ClearSarFilter,
    SelectSarRow {
        idx: usize,
    },
    /// Fold a DockCraft `*.dockcraft.json` results sidecar into the SAR table. `path` picks a
    /// specific file; `None` takes the most recent one DockCraft left in its results directory.
    ImportDockCraftResults {
        path: Option<String>,
    },
    Undo,
    Redo,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CommandInfo {
    pub id: &'static str,
    pub label: &'static str,
    pub shortcut: Option<&'static str>,
}

#[derive(Debug, Error)]
pub enum CommandError {
    #[error("nothing to undo")]
    NothingToUndo,
    #[error("nothing to redo")]
    NothingToRedo,
    #[error("{0}")]
    Failed(String),
}

pub struct CommandOutcome {
    pub message: String,
}

pub struct CommandRegistry;

impl CommandRegistry {
    pub fn catalog() -> Vec<CommandInfo> {
        vec![
            CommandInfo {
                id: "clean_2d",
                label: "Clean Structure 2D",
                shortcut: Some("Ctrl+L"),
            },
            CommandInfo {
                id: "load_demo_50",
                label: "Load SAR Demo (50 Compounds)",
                shortcut: Some("Ctrl+E"),
            },
            CommandInfo {
                id: "swap_tetrazole",
                label: "Bioisostere: -COOH -> Tetrazole",
                shortcut: None,
            },
            CommandInfo {
                id: "add_to_sar",
                label: "Add Canvas Molecule to SAR Table",
                shortcut: None,
            },
            CommandInfo {
                id: "filter_canvas",
                label: "Filter SAR by Canvas Substructure",
                shortcut: None,
            },
            CommandInfo {
                id: "clear_filter",
                label: "Clear SAR Filter",
                shortcut: None,
            },
            CommandInfo {
                id: "clear_canvas",
                label: "Clear Drawing Canvas",
                shortcut: None,
            },
            CommandInfo {
                id: "import_dockcraft_results",
                label: "Import DockCraft Docking Results",
                shortcut: None,
            },
            CommandInfo {
                id: "undo",
                label: "Undo",
                shortcut: Some("Ctrl+Z"),
            },
            CommandInfo {
                id: "redo",
                label: "Redo",
                shortcut: Some("Ctrl+Y"),
            },
        ]
    }

    pub fn execute(
        session: &mut ChemSession,
        cmd: Command,
    ) -> Result<CommandOutcome, CommandError> {
        match cmd {
            Command::Undo => {
                let prev = session
                    .undo_stack
                    .pop()
                    .ok_or(CommandError::NothingToUndo)?;
                session.redo_stack.push(SessionSnapshot {
                    canvas_mol: Arc::clone(&session.canvas_mol),
                    reaction_arrows: Arc::clone(&session.reaction_arrows),
                    sar_table: Arc::clone(&session.sar_table),
                    selected_sar_idx: session.selected_sar_idx,
                });
                session.restore_snapshot(prev);
                Ok(CommandOutcome {
                    message: "Undid last structural/SAR edit".to_string(),
                })
            }
            Command::Redo => {
                let next = session
                    .redo_stack
                    .pop()
                    .ok_or(CommandError::NothingToRedo)?;
                session.undo_stack.push(SessionSnapshot {
                    canvas_mol: Arc::clone(&session.canvas_mol),
                    reaction_arrows: Arc::clone(&session.reaction_arrows),
                    sar_table: Arc::clone(&session.sar_table),
                    selected_sar_idx: session.selected_sar_idx,
                });
                session.restore_snapshot(next);
                Ok(CommandOutcome {
                    message: "Redid edit".to_string(),
                })
            }
            Command::LoadCanvasSmiles { smiles, name } => {
                session.push_undo();
                let mut mol =
                    parse_smiles(&smiles).map_err(|e| CommandError::Failed(e.to_string()))?;
                if let Some(n) = name {
                    mol.name = Some(n);
                }
                generate_2d_coords(&mut mol);
                mol_graph::generate_3d_coords(&mut mol);
                session.canvas_mol = Arc::new(mol);
                session.last_mutation_delta = None;
                session.on_canvas_changed();
                Ok(CommandOutcome {
                    message: format!("Loaded SMILES {}", session.canvas_mol.canonical_smiles()),
                })
            }
            Command::ClearCanvas => {
                session.push_undo();
                session.canvas_mol = Arc::new(MolGraph::new());
                session.reaction_arrows = Arc::new(Vec::new());
                session.live_admet = None;
                session.last_mutation_delta = None;
                Ok(CommandOutcome {
                    message: "Cleared drawing canvas".to_string(),
                })
            }
            Command::CleanStructure2D => {
                session.push_undo();
                let mut mol = (*session.canvas_mol).clone();
                clean_structure_2d(&mut mol, 45);
                session.canvas_mol = Arc::new(mol);
                Ok(CommandOutcome {
                    message: "Cleaned 2D bond lengths & angles".to_string(),
                })
            }
            Command::StampRing {
                template_name,
                x,
                y,
                snap_atom,
                snap_bond,
            } => {
                let tmpl = parse_template(&template_name).ok_or_else(|| {
                    CommandError::Failed(format!("unknown ring template '{template_name}'"))
                })?;
                session.push_undo();
                let mut mol = (*session.canvas_mol).clone();
                let target = if let Some(bid) = snap_bond {
                    SnapTarget::Bond {
                        bond_idx: bid,
                        cursor_x: x,
                        cursor_y: y,
                    }
                } else if let Some(aid) = snap_atom {
                    SnapTarget::Atom { atom_idx: aid }
                } else {
                    SnapTarget::Free { x, y }
                };
                stamp_ring_template(&mut mol, tmpl, target)
                    .map_err(|e| CommandError::Failed(e.to_string()))?;
                mol_graph::generate_3d_coords(&mut mol);
                session.canvas_mol = Arc::new(mol);
                session.on_canvas_changed();
                Ok(CommandOutcome {
                    message: format!("Stamped {}", tmpl.label()),
                })
            }
            Command::AddAtomAndBond {
                from_atom,
                symbol,
                x,
                y,
                order,
                stereo,
            } => {
                session.push_undo();
                let mut mol = (*session.canvas_mol).clone();
                let mut atom = MolAtom::new(mol.atoms.len(), &symbol)
                    .map_err(|e| CommandError::Failed(e.to_string()))?;
                atom.coords_2d = Some([x, y]);
                let nid = mol.add_atom(atom);
                if let Some(prev) = from_atom {
                    mol.add_bond(prev, nid, order, stereo)
                        .map_err(|e| CommandError::Failed(e.to_string()))?;
                }
                mol.perceive()
                    .map_err(|e| CommandError::Failed(e.to_string()))?;
                mol_graph::generate_3d_coords(&mut mol);
                session.canvas_mol = Arc::new(mol);
                session.on_canvas_changed();
                Ok(CommandOutcome {
                    message: format!("Added atom {symbol}"),
                })
            }
            Command::SetAtomElement { atom_idx, symbol } => {
                session.push_undo();
                let mut mol = (*session.canvas_mol).clone();
                mol.set_atom_element(atom_idx, &symbol)
                    .map_err(|e| CommandError::Failed(e.to_string()))?;
                mol_graph::generate_3d_coords(&mut mol);
                session.canvas_mol = Arc::new(mol);
                session.on_canvas_changed();
                Ok(CommandOutcome {
                    message: format!("Changed atom #{} to {symbol}", atom_idx + 1),
                })
            }
            Command::SetBondStereo { bond_idx, stereo } => {
                let mut mol = (*session.canvas_mol).clone();
                let b = mol.bonds.get_mut(bond_idx).ok_or_else(|| {
                    CommandError::Failed(format!("invalid bond index {bond_idx}"))
                })?;
                session.push_undo();
                b.stereo = stereo;
                mol.perceive()
                    .map_err(|e| CommandError::Failed(e.to_string()))?;
                mol_graph::generate_3d_coords(&mut mol);
                session.canvas_mol = Arc::new(mol);
                session.on_canvas_changed();
                Ok(CommandOutcome {
                    message: format!("Updated stereo on bond #{}", bond_idx + 1),
                })
            }
            Command::LoadChemicalName { name } => {
                let mut mol = resolve_chemical_name(&name)
                    .map_err(|e| CommandError::Failed(e.to_string()))?;
                generate_2d_coords(&mut mol);
                mol_graph::generate_3d_coords(&mut mol);
                session.push_undo();
                let label = mol.name.clone().unwrap_or_else(|| mol.formula());
                let smi = mol.canonical_smiles();
                session.canvas_mol = Arc::new(mol);
                session.last_mutation_delta = None;
                session.on_canvas_changed();
                Ok(CommandOutcome {
                    message: format!("Resolved '{label}' -> {smi}"),
                })
            }
            Command::LoadCdxml { xml } => {
                let mut mol = parse_cdxml(&xml).map_err(|e| CommandError::Failed(e.to_string()))?;
                if mol.atoms.iter().any(|a| a.coords_2d.is_none()) {
                    generate_2d_coords(&mut mol);
                }
                mol_graph::generate_3d_coords(&mut mol);
                session.push_undo();
                let f = mol.formula();
                let smi = mol.canonical_smiles();
                session.canvas_mol = Arc::new(mol);
                session.last_mutation_delta = None;
                session.on_canvas_changed();
                Ok(CommandOutcome {
                    message: format!("Imported CDXML ({f}, {smi})"),
                })
            }
            Command::ConnectAtoms {
                u,
                v,
                order,
                stereo,
            } => {
                if u == v {
                    return Err(CommandError::Failed(
                        "cannot connect an atom to itself".to_string(),
                    ));
                }
                let mut mol = (*session.canvas_mol).clone();
                if u >= mol.atoms.len() || v >= mol.atoms.len() {
                    return Err(CommandError::Failed("atom index out of range".to_string()));
                }
                session.push_undo();
                if let Some((_, bid)) = mol.adj[u].iter().copied().find(|&(nb, _)| nb == v) {
                    mol.bonds[bid].order = order;
                    mol.bonds[bid].stereo = stereo;
                } else {
                    mol.add_bond(u, v, order, stereo)
                        .map_err(|e| CommandError::Failed(e.to_string()))?;
                }
                mol.perceive()
                    .map_err(|e| CommandError::Failed(e.to_string()))?;
                mol_graph::generate_3d_coords(&mut mol);
                session.canvas_mol = Arc::new(mol);
                session.on_canvas_changed();
                Ok(CommandOutcome {
                    message: format!("Connected atoms #{} and #{}", u + 1, v + 1),
                })
            }
            Command::AdjustAtomCharge { atom_idx, delta } => {
                let mut mol = (*session.canvas_mol).clone();
                let a = mol.atoms.get_mut(atom_idx).ok_or_else(|| {
                    CommandError::Failed(format!("invalid atom index {atom_idx}"))
                })?;
                session.push_undo();
                a.formal_charge =
                    (i16::from(a.formal_charge) + i16::from(delta)).clamp(-4, 4) as i8;
                a.from_bracket = false;
                a.explicit_h = 0;
                let q = a.formal_charge;
                mol.perceive()
                    .map_err(|e| CommandError::Failed(e.to_string()))?;
                mol_graph::generate_3d_coords(&mut mol);
                session.canvas_mol = Arc::new(mol);
                session.on_canvas_changed();
                Ok(CommandOutcome {
                    message: format!("Atom #{} charge set to {q:+}", atom_idx + 1),
                })
            }
            Command::CycleBondOrder { bond_idx } => {
                let mut mol = (*session.canvas_mol).clone();
                let b = mol.bonds.get_mut(bond_idx).ok_or_else(|| {
                    CommandError::Failed(format!("invalid bond index {bond_idx}"))
                })?;
                session.push_undo();
                b.order = match b.order {
                    BondOrder::Single => BondOrder::Double,
                    BondOrder::Double => BondOrder::Triple,
                    BondOrder::Triple | BondOrder::Aromatic | BondOrder::Dative => {
                        BondOrder::Single
                    }
                };
                if b.order != BondOrder::Single {
                    b.stereo = BondStereo::None;
                }
                let ord_label = format!("{:?}", b.order);
                mol.perceive()
                    .map_err(|e| CommandError::Failed(e.to_string()))?;
                mol_graph::generate_3d_coords(&mut mol);
                session.canvas_mol = Arc::new(mol);
                session.on_canvas_changed();
                Ok(CommandOutcome {
                    message: format!("Bond #{} -> {ord_label}", bond_idx + 1),
                })
            }
            Command::AttachSuperatomGroup { atom_idx, abbr } => {
                session.push_undo();
                let mut updated = attach_superatom_group(&session.canvas_mol, atom_idx, &abbr)
                    .map_err(|e| CommandError::Failed(e.to_string()))?;
                generate_2d_coords(&mut updated);
                mol_graph::generate_3d_coords(&mut updated);
                let smi = updated.canonical_smiles();
                session.canvas_mol = Arc::new(updated);
                session.last_mutation_delta = None;
                session.on_canvas_changed();
                Ok(CommandOutcome {
                    message: format!("Attached {abbr} -> {smi}"),
                })
            }
            Command::FlipCanvas2D { horizontal } => {
                if session.canvas_mol.atoms.is_empty() {
                    return Err(CommandError::Failed("canvas is empty".to_string()));
                }
                session.push_undo();
                let mut mol = (*session.canvas_mol).clone();
                let mut cx = 0.0_f64;
                let mut cy = 0.0_f64;
                let mut cnt = 0usize;
                for a in &mol.atoms {
                    if let Some([x, y]) = a.coords_2d {
                        cx += x;
                        cy += y;
                        cnt += 1;
                    }
                }
                if cnt > 0 {
                    cx /= cnt as f64;
                    cy /= cnt as f64;
                    for a in &mut mol.atoms {
                        if let Some([x, y]) = &mut a.coords_2d {
                            if horizontal {
                                *x = 2.0 * cx - *x;
                            } else {
                                *y = 2.0 * cy - *y;
                            }
                        }
                    }
                    // Flip Up <-> Down stereo wedges so physical 3D chirality is preserved on 2D reflection!
                    for b in &mut mol.bonds {
                        b.stereo = match b.stereo {
                            BondStereo::Up => BondStereo::Down,
                            BondStereo::Down => BondStereo::Up,
                            other => other,
                        };
                    }
                }
                session.canvas_mol = Arc::new(mol);
                Ok(CommandOutcome {
                    message: if horizontal {
                        "Flipped structure horizontally (chirality preserved)".to_string()
                    } else {
                        "Flipped structure vertically (chirality preserved)".to_string()
                    },
                })
            }
            Command::RotateCanvas2D { degrees } => {
                if session.canvas_mol.atoms.is_empty() {
                    return Err(CommandError::Failed("canvas is empty".to_string()));
                }
                session.push_undo();
                let mut mol = (*session.canvas_mol).clone();
                let rad = degrees.to_radians();
                let (s, c) = rad.sin_cos();
                let mut cx = 0.0_f64;
                let mut cy = 0.0_f64;
                let mut cnt = 0usize;
                for a in &mol.atoms {
                    if let Some([x, y]) = a.coords_2d {
                        cx += x;
                        cy += y;
                        cnt += 1;
                    }
                }
                if cnt > 0 {
                    cx /= cnt as f64;
                    cy /= cnt as f64;
                    for a in &mut mol.atoms {
                        if let Some([x, y]) = &mut a.coords_2d {
                            let dx = *x - cx;
                            let dy = *y - cy;
                            *x = cx + dx * c - dy * s;
                            *y = cy + dx * s + dy * c;
                        }
                    }
                }
                session.canvas_mol = Arc::new(mol);
                Ok(CommandOutcome {
                    message: format!("Rotated structure by {degrees:+.0}°"),
                })
            }
            Command::EraseAtom { atom_idx } => {
                session.push_undo();
                let mut mol = (*session.canvas_mol).clone();
                mol.remove_atom(atom_idx)
                    .map_err(|e| CommandError::Failed(e.to_string()))?;
                mol_graph::generate_3d_coords(&mut mol);
                session.canvas_mol = Arc::new(mol);
                session.on_canvas_changed();
                Ok(CommandOutcome {
                    message: format!("Erased atom #{}", atom_idx + 1),
                })
            }
            Command::EraseBond { bond_idx } => {
                session.push_undo();
                let mut mol = (*session.canvas_mol).clone();
                mol.remove_bond(bond_idx)
                    .map_err(|e| CommandError::Failed(e.to_string()))?;
                mol_graph::generate_3d_coords(&mut mol);
                session.canvas_mol = Arc::new(mol);
                session.on_canvas_changed();
                Ok(CommandOutcome {
                    message: format!("Erased bond #{}", bond_idx + 1),
                })
            }
            Command::MutateBioisostere { atom_index, group } => {
                if session.canvas_mol.atoms.is_empty() {
                    return Err(CommandError::Failed("canvas is empty".to_string()));
                }
                session.push_undo();
                let old_smi = session.canvas_mol.canonical_smiles();
                let mut mutated = mutate_molecule(&session.canvas_mol, atom_index, &group)
                    .map_err(|e| CommandError::Failed(e.to_string()))?;
                generate_2d_coords(&mut mutated);
                mol_graph::generate_3d_coords(&mut mutated);
                let new_smi = mutated.canonical_smiles();
                let (delta, eval) = compute_admet_delta(
                    &session.admet_provider(),
                    &old_smi,
                    &new_smi,
                    mutated.name.as_deref(),
                )
                .map_err(|e| CommandError::Failed(e.to_string()))?;
                session.canvas_mol = Arc::new(mutated);
                session.live_admet = Some(eval);
                session.last_mutation_delta = Some(delta.clone());
                Ok(CommandOutcome {
                    message: format!(
                        "Mutated to {new_smi} (ΔLogP={:+.2}, ΔTPSA={:+.2}, GI={})",
                        delta.delta_logp, delta.delta_tpsa, delta.new_gi
                    ),
                })
            }
            Command::AddCanvasToSar { pic50 } => {
                if session.canvas_mol.atoms.is_empty() {
                    return Err(CommandError::Failed("canvas is empty".to_string()));
                }
                session.push_undo();
                let mut table = (*session.sar_table).clone();
                let mut mol = (*session.canvas_mol).clone();
                if mol.name.is_none() {
                    mol.name = Some(format!("Drawn_{}", table.rows.len() + 1));
                }
                let row_id = table.add_molecule(mol, pic50, Some("Canvas".to_string()));
                table.enrich_with_admet(&session.admet_provider());
                session.selected_sar_idx = Some(row_id);
                session.sar_table = Arc::new(table);
                session.recompute_sar_analytics();
                Ok(CommandOutcome {
                    message: format!("Added compound #{} to SAR table", row_id + 1),
                })
            }
            Command::ImportDockCraftResults { path } => {
                if session.sar_table.rows.is_empty() {
                    return Err(CommandError::Failed(
                        "the SAR table is empty; load or build the series before importing results"
                            .to_string(),
                    ));
                }
                let sidecar_path = match path {
                    Some(p) => PathBuf::from(p),
                    None => {
                        let dir = crate::dockcraft::default_dockcraft_results_dir();
                        crate::dockcraft::latest_dockcraft_results(&dir).ok_or_else(|| {
                            CommandError::Failed(format!(
                                "no *.dockcraft.json found in {}. Run a screening with \
                                 `dockcraft-cli screen --sidecar` first",
                                dir.display()
                            ))
                        })?
                    }
                };
                let results = crate::dockcraft::read_dockcraft_results(&sidecar_path)
                    .map_err(|e| CommandError::Failed(e.to_string()))?;
                let pose_dir = sidecar_path
                    .parent()
                    .map_or_else(|| PathBuf::from("."), Path::to_path_buf);

                session.push_undo();
                let mut table = (*session.sar_table).clone();
                let outcome =
                    crate::dockcraft::merge_dockcraft_results(&mut table, &results, &pose_dir);
                session.sar_table = Arc::new(table);
                session.recompute_sar_analytics();

                let mut message = format!(
                    "Imported docking results for {} of {} compounds ({} vs {})",
                    outcome.merged,
                    results.rows.len(),
                    if results.run.receptor_name.is_empty() {
                        "receptor"
                    } else {
                        &results.run.receptor_name
                    },
                    if results.run.scoring.is_empty() {
                        "scoring"
                    } else {
                        &results.run.scoring
                    }
                );
                if !outcome.unmatched.is_empty() {
                    // Say which ones did not land. A merge that quietly drops rows is worse than
                    // one that refuses, because the table still looks complete.
                    let shown: Vec<&str> = outcome
                        .unmatched
                        .iter()
                        .take(3)
                        .map(String::as_str)
                        .collect();
                    let more = outcome.unmatched.len().saturating_sub(shown.len());
                    message.push_str(&format!("; no match for {}", shown.join(", ")));
                    if more > 0 {
                        message.push_str(&format!(" and {more} more"));
                    }
                }
                Ok(CommandOutcome { message })
            }
            Command::LoadDemo50Sar => {
                session.push_undo();
                let table = SarTable::demo_50_compounds();
                session.sar_table = Arc::new(table);
                session.recompute_sar_analytics();
                if let Some(first) = session.sar_table.rows.first() {
                    let mut m = first.mol.clone();
                    generate_2d_coords(&mut m);
                    mol_graph::ensure_3d_coords(&mut m);
                    session.canvas_mol = Arc::new(m);
                    session.selected_sar_idx = Some(0);
                    session.on_canvas_changed();
                }
                Ok(CommandOutcome {
                    message: format!(
                        "Loaded {} SAR compounds + AdmetCraft enrichment",
                        session.sar_table.rows.len()
                    ),
                })
            }
            Command::FilterSarByCanvasSubstructure => {
                let smi = session.canvas_mol.canonical_smiles();
                if smi.is_empty() {
                    return Err(CommandError::Failed(
                        "draw a query fragment on canvas first".to_string(),
                    ));
                }
                session.sar_filter.smarts_substructure = Some(smi.clone());
                Ok(CommandOutcome {
                    message: format!("Filtering SAR by substructure '{smi}'"),
                })
            }
            Command::ClearSarFilter => {
                session.sar_filter = SarFilter::default();
                Ok(CommandOutcome {
                    message: "Cleared SAR filters".to_string(),
                })
            }
            Command::SelectSarRow { idx } => {
                let (mut m, row_name) = {
                    let row =
                        session.sar_table.rows.get(idx).ok_or_else(|| {
                            CommandError::Failed(format!("invalid SAR row {idx}"))
                        })?;
                    (row.mol.clone(), row.name.clone())
                };
                generate_2d_coords(&mut m);
                mol_graph::ensure_3d_coords(&mut m);
                session.canvas_mol = Arc::new(m);
                session.selected_sar_idx = Some(idx);
                session.last_mutation_delta = None;
                session.on_canvas_changed();
                Ok(CommandOutcome {
                    message: format!("Selected {row_name}"),
                })
            }
        }
    }
}

fn parse_template(name: &str) -> Option<RingTemplate> {
    match name.trim().to_ascii_lowercase().as_str() {
        "cyclopropane" => Some(RingTemplate::Cyclopropane),
        "cyclobutane" => Some(RingTemplate::Cyclobutane),
        "cyclopentane" => Some(RingTemplate::Cyclopentane),
        "cyclohexane" => Some(RingTemplate::Cyclohexane),
        "benzene" => Some(RingTemplate::Benzene),
        "naphthalene" => Some(RingTemplate::Naphthalene),
        "pyridine" => Some(RingTemplate::Pyridine),
        "indole" => Some(RingTemplate::Indole),
        "tetrazole" => Some(RingTemplate::Tetrazole),
        _ => None,
    }
}
