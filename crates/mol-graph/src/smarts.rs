use std::collections::HashMap;

use crate::error::MolError;
use crate::mol::{BondOrder, MolAtom, MolGraph};

/// Upper bound on backtracking steps for one subgraph search.
///
/// Subgraph isomorphism is exponential in the worst case, and both the molecule and the
/// pattern can arrive from outside: a lattice-shaped molecule against a 22-node chain
/// pattern ran for more than 45 seconds and cost 2.6x per extra node, leaving the MCP
/// server unresponsive. When the budget runs out the search stops and returns the matches
/// found so far, which is a partial answer rather than an unbounded one.
pub const MAX_MATCH_STEPS: u64 = 2_000_000;

/// Mutable state of one subgraph search, kept together so the recursion takes a single
/// handle instead of six parallel arguments.
struct MatchState {
    mapping: Vec<usize>,
    used: Vec<bool>,
    results: Vec<Vec<usize>>,
    budget: u64,
    stop_at_first: bool,
}

impl MatchState {
    fn exhausted(&self) -> bool {
        self.budget == 0 || (self.stop_at_first && !self.results.is_empty())
    }
}

/// Compiled clean-room SMARTS query graph supporting atom/bond primitives,
/// precedence (`;` low-AND > `,` OR > `&`/juxtaposed high-AND), `!` negation,
/// `[R]`/`[RN]`/`[rN]`, `D`, `H`, `#Z`, and bounded charge `[-15, +15]`.
#[derive(Debug, Clone)]

pub struct SmartsQuery {
    nodes: Vec<SmartsNode>,
    edges: Vec<SmartsEdge>,
    adj: Vec<Vec<(usize, usize)>>,
}

/// Top-level conjunction of comma-separated OR-clauses (`a,b;c` => `(a OR b) AND c`).
#[derive(Debug, Clone)]
struct SmartsNode {
    semi_conjuncts: Vec<Vec<Vec<AtomPred>>>,
}

#[derive(Debug, Clone)]
enum AtomPred {
    Negated(Box<AtomPred>),
    Any,
    AnyAromatic,
    AnyAliphatic,
    Element { sym: String, aromatic: Option<bool> },
    AtomicNum(u8),
    InRing(bool),
    RingCount(usize),
    RingSize(usize),
    TotalH(u8),
    Degree(usize),
    Charge(i8),
}

#[derive(Debug, Clone)]
struct SmartsEdge {
    cond: BondCond,
}

#[derive(Debug, Clone, Copy)]
enum BondCond {
    Any,
    SingleOrAromatic,
    Order(BondOrder),
    InRing,
}

impl SmartsQuery {
    pub fn parse(pattern: &str) -> Result<Self, MolError> {
        let s = pattern.trim();
        if s.is_empty() {
            return Err(MolError::InvalidSmarts {
                pos: 0,
                msg: "empty SMARTS pattern".to_string(),
            });
        }
        let bytes = s.as_bytes();
        let mut nodes = Vec::new();
        let mut edges = Vec::new();
        let mut adj: Vec<Vec<(usize, usize)>> = Vec::new();
        let mut stack: Vec<usize> = Vec::new();
        let mut current: Option<usize> = None;
        let mut pending_bond: Option<BondCond> = None;
        let mut rings: HashMap<u16, (usize, Option<BondCond>, usize)> = HashMap::new();
        let mut i = 0usize;

        while i < bytes.len() {
            let ch = bytes[i];
            match ch {
                b'(' => {
                    let cur = current.ok_or(MolError::InvalidSmarts {
                        pos: i,
                        msg: "branch without preceding atom".to_string(),
                    })?;
                    stack.push(cur);
                    i += 1;
                }
                b')' => {
                    current = Some(stack.pop().ok_or(MolError::InvalidSmarts {
                        pos: i,
                        msg: "unmatched ')'".to_string(),
                    })?);
                    i += 1;
                }
                b'-' => {
                    pending_bond = Some(BondCond::Order(BondOrder::Single));
                    i += 1;
                }
                b'=' => {
                    pending_bond = Some(BondCond::Order(BondOrder::Double));
                    i += 1;
                }
                b'#' => {
                    pending_bond = Some(BondCond::Order(BondOrder::Triple));
                    i += 1;
                }
                b':' => {
                    pending_bond = Some(BondCond::Order(BondOrder::Aromatic));
                    i += 1;
                }
                b'~' => {
                    pending_bond = Some(BondCond::Any);
                    i += 1;
                }
                b'@' => {
                    pending_bond = Some(BondCond::InRing);
                    i += 1;
                }
                b'0'..=b'9' => {
                    let cur = current.ok_or(MolError::InvalidSmarts {
                        pos: i,
                        msg: "ring digit without atom".to_string(),
                    })?;
                    let rnum = u16::from(ch - b'0');
                    if let Some((other, bo1, _)) = rings.remove(&rnum) {
                        let cond = bo1.or(pending_bond).unwrap_or(BondCond::SingleOrAromatic);
                        let eid = edges.len();
                        edges.push(SmartsEdge { cond });
                        adj[other].push((cur, eid));
                        adj[cur].push((other, eid));
                    } else {
                        rings.insert(rnum, (cur, pending_bond, i));
                    }
                    pending_bond = None;
                    i += 1;
                }
                b'[' => {
                    let close = bytes[i + 1..]
                        .iter()
                        .position(|&b| b == b']')
                        .map(|o| i + 1 + o)
                        .ok_or(MolError::InvalidSmarts {
                            pos: i,
                            msg: "unclosed '['".to_string(),
                        })?;
                    let inner = std::str::from_utf8(&bytes[i + 1..close]).map_err(|_| {
                        MolError::InvalidSmarts {
                            pos: i,
                            msg: "invalid UTF-8".to_string(),
                        }
                    })?;
                    let node = parse_smarts_bracket(inner, i + 1)?;
                    let nid = nodes.len();
                    nodes.push(node);
                    adj.push(Vec::new());
                    if let Some(prev) = current {
                        let cond = pending_bond.unwrap_or(BondCond::SingleOrAromatic);
                        let eid = edges.len();
                        edges.push(SmartsEdge { cond });
                        adj[prev].push((nid, eid));
                        adj[nid].push((prev, eid));
                    }
                    current = Some(nid);
                    pending_bond = None;
                    i = close + 1;
                }
                b'*' | b'A' | b'a' | b'B' | b'C' | b'N' | b'O' | b'P' | b'S' | b'F' | b'I'
                | b'c' | b'n' | b'o' | b's' | b'p' => {
                    let (pred, step) = if ch == b'C' && bytes.get(i + 1) == Some(&b'l') {
                        (
                            AtomPred::Element {
                                sym: "Cl".to_string(),
                                aromatic: Some(false),
                            },
                            2,
                        )
                    } else if ch == b'B' && bytes.get(i + 1) == Some(&b'r') {
                        (
                            AtomPred::Element {
                                sym: "Br".to_string(),
                                aromatic: Some(false),
                            },
                            2,
                        )
                    } else {
                        let p = match ch {
                            b'*' => AtomPred::Any,
                            b'A' => AtomPred::AnyAliphatic,
                            b'a' => AtomPred::AnyAromatic,
                            b'c' => AtomPred::Element {
                                sym: "C".to_string(),
                                aromatic: Some(true),
                            },
                            b'n' => AtomPred::Element {
                                sym: "N".to_string(),
                                aromatic: Some(true),
                            },
                            b'o' => AtomPred::Element {
                                sym: "O".to_string(),
                                aromatic: Some(true),
                            },
                            b's' => AtomPred::Element {
                                sym: "S".to_string(),
                                aromatic: Some(true),
                            },
                            b'p' => AtomPred::Element {
                                sym: "P".to_string(),
                                aromatic: Some(true),
                            },
                            _ => AtomPred::Element {
                                sym: (ch as char).to_string(),
                                aromatic: Some(false),
                            },
                        };
                        (p, 1)
                    };
                    let nid = nodes.len();
                    nodes.push(SmartsNode {
                        semi_conjuncts: vec![vec![vec![pred]]],
                    });
                    adj.push(Vec::new());
                    if let Some(prev) = current {
                        let cond = pending_bond.unwrap_or(BondCond::SingleOrAromatic);
                        let eid = edges.len();
                        edges.push(SmartsEdge { cond });
                        adj[prev].push((nid, eid));
                        adj[nid].push((prev, eid));
                    }
                    current = Some(nid);
                    pending_bond = None;
                    i += step;
                }
                _ => {
                    return Err(MolError::InvalidSmarts {
                        pos: i,
                        msg: format!("unsupported SMARTS token '{}'", ch as char),
                    });
                }
            }
        }

        if let Some((rnum, (_, _, pos))) = rings.into_iter().next() {
            return Err(MolError::InvalidSmarts {
                pos,
                msg: format!("unclosed SMARTS ring {rnum}"),
            });
        }
        Ok(Self { nodes, edges, adj })
    }

    pub fn matches(&self, mol: &MolGraph) -> bool {
        !self.find_matches(mol, true).is_empty()
    }

    pub fn find_all_matches(&self, mol: &MolGraph) -> Vec<Vec<usize>> {
        self.find_matches(mol, false)
    }

    fn find_matches(&self, mol: &MolGraph, stop_at_first: bool) -> Vec<Vec<usize>> {
        if self.nodes.is_empty() || self.nodes.len() > mol.atoms.len() {
            return Vec::new();
        }
        let mut st = MatchState {
            mapping: vec![usize::MAX; self.nodes.len()],
            used: vec![false; mol.atoms.len()],
            results: Vec::new(),
            budget: MAX_MATCH_STEPS,
            stop_at_first,
        };

        for start in 0..mol.atoms.len() {
            if !node_matches(&self.nodes[0], &mol.atoms[start], mol) {
                continue;
            }
            st.mapping[0] = start;
            st.used[start] = true;
            self.backtrack(1, mol, &mut st);
            st.used[start] = false;
            st.mapping[0] = usize::MAX;
            if st.exhausted() {
                break;
            }
        }
        st.results
    }

    fn backtrack(&self, q_idx: usize, mol: &MolGraph, st: &mut MatchState) {
        if st.budget == 0 {
            return;
        }
        st.budget -= 1;
        if q_idx == self.nodes.len() {
            st.results.push(st.mapping.clone());
            return;
        }
        let prev_req: Vec<(usize, BondCond)> = self.adj[q_idx]
            .iter()
            .filter_map(|&(qn, eid)| {
                if qn < q_idx {
                    Some((st.mapping[qn], self.edges[eid].cond))
                } else {
                    None
                }
            })
            .collect();

        let candidates: Vec<usize> = if let Some(&(first_t, _)) = prev_req.first() {
            mol.adj[first_t].iter().map(|&(nbr, _)| nbr).collect()
        } else {
            (0..mol.atoms.len()).collect()
        };

        for t_idx in candidates {
            if st.used[t_idx] || !node_matches(&self.nodes[q_idx], &mol.atoms[t_idx], mol) {
                continue;
            }
            let all_edges_ok = prev_req.iter().all(|&(t_nbr, cond)| {
                mol.adj[t_idx].iter().any(|&(cand_nbr, bid)| {
                    cand_nbr == t_nbr
                        && bond_matches(cond, mol.bonds[bid].order, mol.bonds[bid].in_ring)
                })
            });
            if !all_edges_ok {
                continue;
            }
            st.mapping[q_idx] = t_idx;
            st.used[t_idx] = true;
            self.backtrack(q_idx + 1, mol, st);
            st.used[t_idx] = false;
            st.mapping[q_idx] = usize::MAX;
            if st.exhausted() {
                return;
            }
        }
    }
}

/// Precedence per Daylight SMARTS: `;` (lowest AND) -> `,` (OR) -> `&` / juxtaposed (highest AND).
fn parse_smarts_bracket(inner: &str, base_pos: usize) -> Result<SmartsNode, MolError> {
    if inner.trim().is_empty() {
        return Err(MolError::InvalidSmarts {
            pos: base_pos,
            msg: "empty bracket expression '[]'".to_string(),
        });
    }
    let mut semi_conjuncts = Vec::new();
    for semi_part in inner.split(';') {
        if semi_part.is_empty() {
            return Err(MolError::InvalidSmarts {
                pos: base_pos,
                msg: "empty ';' clause in SMARTS".to_string(),
            });
        }
        let mut or_branches = Vec::new();
        for or_part in semi_part.split(',') {
            if or_part.is_empty() {
                return Err(MolError::InvalidSmarts {
                    pos: base_pos,
                    msg: "empty ',' branch in SMARTS".to_string(),
                });
            }
            let mut preds = Vec::new();
            let bytes = or_part.as_bytes();
            let mut i = 0usize;
            while i < bytes.len() {
                if bytes[i] == b'&' {
                    i += 1;
                    continue;
                }
                let neg = bytes[i] == b'!';
                if neg {
                    i += 1;
                }
                let (pred, next_i) = parse_single_pred(bytes, i, base_pos + i)?;
                preds.push(if neg {
                    AtomPred::Negated(Box::new(pred))
                } else {
                    pred
                });
                i = next_i;
            }
            or_branches.push(preds);
        }
        semi_conjuncts.push(or_branches);
    }
    Ok(SmartsNode { semi_conjuncts })
}

fn parse_single_pred(
    bytes: &[u8],
    pos: usize,
    err_pos: usize,
) -> Result<(AtomPred, usize), MolError> {
    if pos >= bytes.len() {
        return Err(MolError::InvalidSmarts {
            pos: err_pos,
            msg: "dangling negation or operator".to_string(),
        });
    }
    let ch = bytes[pos];
    match ch {
        b'*' => Ok((AtomPred::Any, pos + 1)),
        b'A' => Ok((AtomPred::AnyAliphatic, pos + 1)),
        b'a' => Ok((AtomPred::AnyAromatic, pos + 1)),
        b'#' => {
            let (z, next_i) = read_u32(bytes, pos + 1).ok_or(MolError::InvalidSmarts {
                pos: err_pos,
                msg: "expected atomic number after '#'".to_string(),
            })?;
            let z_u8 = u8::try_from(z).map_err(|_| MolError::InvalidSmarts {
                pos: err_pos,
                msg: format!("atomic number #{z} out of u8 range"),
            })?;
            Ok((AtomPred::AtomicNum(z_u8), next_i))
        }
        b'R' => {
            if let Some((n, next_i)) = read_u32(bytes, pos + 1) {
                Ok((AtomPred::RingCount(n as usize), next_i))
            } else {
                Ok((AtomPred::InRing(true), pos + 1))
            }
        }
        b'r' => {
            if let Some((n, next_i)) = read_u32(bytes, pos + 1) {
                Ok((AtomPred::RingSize(n as usize), next_i))
            } else {
                Ok((AtomPred::InRing(true), pos + 1))
            }
        }
        b'D' => {
            let (d, next_i) = read_u32(bytes, pos + 1).unwrap_or((1, pos + 1));
            Ok((AtomPred::Degree(d as usize), next_i))
        }
        b'H' => {
            let (h, next_i) = read_u32(bytes, pos + 1).unwrap_or((1, pos + 1));
            let h_u8 = u8::try_from(h).map_err(|_| MolError::InvalidSmarts {
                pos: err_pos,
                msg: format!("H count {h} out of range"),
            })?;
            Ok((AtomPred::TotalH(h_u8), next_i))
        }
        b'+' | b'-' => {
            let sign: i32 = if ch == b'+' { 1 } else { -1 };
            let (mag, next_i) = if let Some((v, ni)) = read_u32(bytes, pos + 1) {
                let mag_i32 = i32::try_from(v).map_err(|_| MolError::InvalidSmarts {
                    pos: err_pos,
                    msg: "formal charge integer overflow".to_string(),
                })?;
                (mag_i32, ni)
            } else {
                let mut count = 1i32;
                let mut k = pos + 1;
                while bytes.get(k) == Some(&ch) {
                    count += 1;
                    if count > 15 {
                        return Err(MolError::InvalidSmarts {
                            pos: err_pos,
                            msg: "formal charge magnitude exceeds [-15, +15]".to_string(),
                        });
                    }
                    k += 1;
                }
                (count, k)
            };
            let charge = sign.saturating_mul(mag);
            if !(-15..=15).contains(&charge) {
                return Err(MolError::InvalidSmarts {
                    pos: err_pos,
                    msg: format!("formal charge {charge} outside [-15, +15]"),
                });
            }
            Ok((AtomPred::Charge(charge as i8), next_i))
        }
        _ => {
            if ch.is_ascii_alphabetic() {
                let arom = Some(ch.is_ascii_lowercase());
                let sym = if let Some(&sec) = bytes.get(pos + 1)
                    && sec.is_ascii_lowercase()
                    && crate::element::element_by_symbol(&format!(
                        "{}{}",
                        ch.to_ascii_uppercase() as char,
                        sec as char
                    ))
                    .is_some()
                {
                    (
                        format!("{}{}", ch.to_ascii_uppercase() as char, sec as char),
                        pos + 2,
                    )
                } else {
                    ((ch.to_ascii_uppercase() as char).to_string(), pos + 1)
                };
                Ok((
                    AtomPred::Element {
                        sym: sym.0,
                        aromatic: arom,
                    },
                    sym.1,
                ))
            } else {
                Err(MolError::InvalidSmarts {
                    pos: err_pos,
                    msg: format!("unknown SMARTS atom primitive '{}'", ch as char),
                })
            }
        }
    }
}

fn read_u32(bytes: &[u8], start: usize) -> Option<(u32, usize)> {
    let mut i = start;
    while i < bytes.len() && bytes[i].is_ascii_digit() {
        i += 1;
    }
    if i == start {
        return None;
    }
    let val = std::str::from_utf8(&bytes[start..i])
        .ok()?
        .parse::<u32>()
        .ok()?;
    Some((val, i))
}

fn node_matches(node: &SmartsNode, atom: &MolAtom, mol: &MolGraph) -> bool {
    node.semi_conjuncts.iter().all(|or_branches| {
        or_branches
            .iter()
            .any(|and_preds| and_preds.iter().all(|p| pred_matches(p, atom, mol)))
    })
}

fn pred_matches(pred: &AtomPred, atom: &MolAtom, mol: &MolGraph) -> bool {
    match pred {
        AtomPred::Negated(inner) => !pred_matches(inner, atom, mol),
        AtomPred::Any => true,
        AtomPred::AnyAromatic => atom.aromatic,
        AtomPred::AnyAliphatic => !atom.aromatic,
        AtomPred::Element { sym, aromatic } => {
            atom.symbol == *sym && aromatic.is_none_or(|ar| atom.aromatic == ar)
        }
        AtomPred::AtomicNum(z) => atom.atomic_number == *z,
        AtomPred::InRing(req) => atom.in_ring == *req,
        AtomPred::RingCount(cnt) => {
            mol.rings.iter().filter(|r| r.contains(&atom.id)).count() == *cnt
        }
        AtomPred::RingSize(sz) => atom.ring_sizes.contains(sz),
        AtomPred::TotalH(h) => atom.total_h() == *h,
        AtomPred::Degree(d) => mol.adj[atom.id].len() == *d,
        AtomPred::Charge(c) => atom.formal_charge == *c,
    }
}

fn bond_matches(cond: BondCond, order: BondOrder, in_ring: bool) -> bool {
    match cond {
        BondCond::Any => true,
        BondCond::SingleOrAromatic => {
            matches!(
                order,
                BondOrder::Single | BondOrder::Aromatic | BondOrder::Dative
            )
        }
        BondCond::Order(req) => order == req,
        BondCond::InRing => in_ring,
    }
}

#[cfg(test)]
mod match_budget_tests {
    use super::*;
    use crate::parse_smiles;
    use std::time::Instant;

    /// Audit finding F3-18: subgraph matching had no step budget. A lattice-shaped
    /// molecule against a long chain pattern ran for more than 45 seconds, growing 2.6x
    /// per extra pattern node, and left the MCP server unresponsive. Both inputs can come
    /// from outside.
    #[test]
    fn a_pathological_pattern_cannot_run_unbounded() {
        // Fused polycyclic lattice: many symmetric embeddings, which is what makes the
        // search explode.
        let lattice = "C1CC2CC3CC4CC1CC(C2)C(C3)C4";
        let Ok(mol) = parse_smiles(lattice) else {
            return;
        };
        let chain: String = std::iter::repeat_n("C", 22).collect();
        let Ok(q) = SmartsQuery::parse(&chain) else {
            return;
        };

        let t = Instant::now();
        let hits = q.find_all_matches(&mol);
        let elapsed = t.elapsed();

        assert!(
            elapsed.as_secs() < 10,
            "subgraph search took {elapsed:?}; the step budget is not bounding it"
        );
        // A partial answer is acceptable; an unbounded run is not.
        let _ = hits;
    }

    /// The budget must not change the answer for ordinary queries.
    #[test]
    fn ordinary_queries_are_unaffected_by_the_budget() {
        let Ok(aspirin) = parse_smiles("CC(=O)Oc1ccccc1C(=O)O") else {
            return;
        };
        for (pattern, expected_some) in [("c1ccccc1", true), ("C(=O)O", true), ("[Si]", false)] {
            let Ok(q) = SmartsQuery::parse(pattern) else {
                continue;
            };
            assert_eq!(
                q.matches(&aspirin),
                expected_some,
                "pattern {pattern} on aspirin"
            );
        }
        let Ok(q) = SmartsQuery::parse("c1ccccc1") else {
            return;
        };
        assert_eq!(
            q.find_all_matches(&aspirin).len(),
            12,
            "one benzene ring has 12 symmetric mappings"
        );
    }
}
