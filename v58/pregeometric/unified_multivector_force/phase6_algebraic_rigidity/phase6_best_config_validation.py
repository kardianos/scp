#!/usr/bin/env python3
"""
phase6_best_config_validation.py

Validation run using the single best-performing pressure configuration found so far:

  Mode: BOTH (pressure affects both attachment bias and quadratic modulation)
  Best loss from search: -0.7132

Best coefficients (from the winning "Both" run):

  term_rank_deviation           : +1.9100
  term_rank_abs                 : -1.4633
  term_closure_defect           : +0.7276
  term_balance                  : -2.3360
  term_rho_rank_interaction     : +0.4099
  term_rho_closure_interaction  : +1.0997
  term_biv_magnitude            : +1.9638
  term_rank_biv_interaction     : +0.2530

This script runs a proper full-scale discrete simulation (500–600+ nodes)
using the exact locked living candidate + the above pressure weighting
applied in BOTH channels:

  - Additive effect on attachment bias (living activity)
  - Modulation effect on effective λ and μ during the omega update

It follows the same protocol as Phase 4 (long trajectories + perturbation + recovery)
so the results are directly comparable.

The goal is to answer:
  "When we take the best pressure configuration discovered by optimization and apply it
   at full scale in the discrete model, what effective dimensionality and stability
   does the system actually achieve?"

Run:
    python phase6_best_config_validation.py
"""

from __future__ import annotations
import sys
import math
import random
from pathlib import Path
from typing import Dict, List, Optional
import json

# ---------------------------------------------------------------------------
# Reuse Phase 1 canonical model (read-only)
# ---------------------------------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
PHASE1_DIR = (SCRIPT_DIR.parent / "phase1_minimal_relational_models").resolve()
if str(PHASE1_DIR) not in sys.path:
    sys.path.insert(0, str(PHASE1_DIR))

from minimal_graph_model import (
    RelationalNode,
    CausalEdge,
    create_seed_graph,
    update_all_omegas,
    compute_growth_curve,
    estimate_d_eff,
    light_cone_proxy,
    degree_histograms,
    causal_past_ball_list,
    node_living_activity as base_node_living_activity,
    LAMBDA_NL,
    MU_NL,
    RHO_CRIT,
    FG_DEFAULT_BACKGROUND,
)

# ---------------------------------------------------------------------------
# Best coefficients from the winning "Both" run (loss = -0.7132)
# ---------------------------------------------------------------------------
BEST_COEFFS = {
    "term_rank_deviation": 1.9100,
    "term_rank_abs": -1.4633,
    "term_closure_defect": 0.7276,
    "term_balance": -2.3360,
    "term_rho_rank_interaction": 0.4099,
    "term_rho_closure_interaction": 1.0997,
    "term_biv_magnitude": 1.9638,
    "term_rank_biv_interaction": 0.2530,
}

# Extra "spin" encouragement weight for the new direction-creation term
# (this is the "bit of spin" knob to try to blow out the 1D nature)
NEW_SPIN_WEIGHT = 1.8   # legacy single rank-delta term

# Multi-orientation torsion (spin) weights — the new mechanism
# These let the pressure term 'see' torsion in three independent protected planes
# plus a strong signal for truly orthogonal new torsion.
TORSION_ORIENT_0_WEIGHT = 1.05
TORSION_ORIENT_1_WEIGHT = 1.15
TORSION_ORIENT_2_WEIGHT = 0.95
TORSION_ORTHO_WEIGHT    = 1.40   # strongest bias toward genuinely new orientations

# External pressure (kept at original strength for this run)
EXTERNAL_PRESSURE_BOOST = 1.0

# Internal node repulsion weight (legacy single term)
INTERNAL_REPULSION_WEIGHT = 1.9

# Seven distinct repulsion channels — each a different "orientation" of crowding cost
REPULSION_DIRECT_W   = 1.10   # Type 1: direct children (classic)
REPULSION_SIBLING_W  = 0.85   # Type 2: sibling / parent-sharing crowding
REPULSION_CONE_W     = 0.70   # Type 3: causal cone (depth ~2) density
REPULSION_LAYER_W    = 0.65   # Type 4: same/adjacent layer local density
REPULSION_PROT_W     = 1.25   # Type 5: protected children only (higher-value crowding)
REPULSION_BIVEC_W    = 0.55   # Type 6: nodes with similar protected bivector orientation
REPULSION_ANCESTOR_W = 0.80   # Type 7: deep (2-gen) ancestor set redundancy

# Seven distinct torsion (spin) pressure features (expansion from the original 4)
TORSION_PLANE_0_W      = 0.95
TORSION_PLANE_1_W      = 1.05
TORSION_PLANE_2_W      = 0.90
TORSION_ORTHOGONAL_W   = 1.30   # strongest in previous 4
TORSION_SIGNED_W       = 0.75   # signed torque / handedness
TORSION_DIVERSITY_W    = 0.65   # anti-alignment / orientation diversity
TORSION_MAGNITUDE_GROW_W = 0.85 # growth in total protected bivector strength

# External pressure on outer / unconnected nodes (new boundary term)
OUTER_UNCONNECTED_PRESSURE_W = 1.4   # tunable; positive = extra pressure on low-usage frontier nodes

# Connection-scaled repulsion (the 10x on unconnected, scaling down as it connects)
CONNECTION_SCALED_REPULSION_W = 2.2   # strength of the dynamic 10x→1x repulsion effect

# --- Three experimental "extra connection / particle folding" mechanisms ---
# These can be enabled independently (weight > 0) for the 7-run matrix.
MECH_A_EXTRA_CONN_STAB_W = 0.0      # A: new protected connection stabilized by density
MECH_B_DENSITY_RIGIDITY_W = 0.0     # B: protection/rigidity cost becomes density-dependent
MECH_C_DENSITY_CLOSURE_W = 0.0      # C: preferred algebraic closure becomes density-dependent

# Chiral helical spin bias — two versions
CHIRAL_HELICAL_BIAS_W        = 1.5   # direct (proportional to connections)
CHIRAL_HELICAL_BIAS_INVERSE_W = 0.0   # inverse (strongest on sparsely connected cores)

# --- Two experimental "Inside Space" mechanisms ---
# Version 1: Explicit reward for creating "inside volume" + localized collapse on the particle cores
INSIDE_VOLUME_REWARD_V1_W = 0.0
VOLUME_COLLAPSE_V1_W      = 0.0

# Version 2: Negative pressure inside the enclosed regions + collapse on the boundary particles
INSIDE_NEGATIVE_PRESSURE_V2_W = 0.0
BOUNDARY_COLLAPSE_V2_W        = 0.0

# Pressure weight (how strongly the pressure scalar affects the system)
PRESSURE_WEIGHT = 0.8          # for attachment
MODULATION_WEIGHT = 0.8        # for quadratic modulation

# ---------------------------------------------------------------------------
# Pressure feature functions (same as in the search)
# ---------------------------------------------------------------------------

def _quick_parents(node: RelationalNode, nodes: Dict[int, RelationalNode]) -> List[RelationalNode]:
    return [nodes[p.parent_id] for p in node.parents if p.parent_id in nodes]

def feat_external_pressure(node: RelationalNode, nodes: Dict[int, RelationalNode]) -> float:
    parents = _quick_parents(node, nodes)
    if not parents:
        return 0.0
    return sum(p.rho for p in parents) / len(parents)

def feat_internal_pressure(node: RelationalNode, nodes: Dict[int, RelationalNode] = None) -> float:
    return node.rho

def feat_push_pull(node: RelationalNode, nodes: Dict[int, RelationalNode]) -> float:
    parents = _quick_parents(node, nodes)
    if not parents:
        return 0.0
    parent_avg = sum(p.rho for p in parents) / len(parents)
    return node.rho - parent_avg

def feat_local_contrast(node: RelationalNode, nodes: Dict[int, RelationalNode]) -> float:
    parents = _quick_parents(node, nodes)
    if not parents:
        return 0.0
    mean = sum(p.rho for p in parents) / len(parents)
    var = sum((p.rho - mean)**2 for p in parents) / len(parents)
    return math.sqrt(var)

def feat_protected_spin_pressure(node: RelationalNode, nodes: Dict[int, RelationalNode] = None) -> float:
    if not node.protected:
        return 0.0
    b = [abs(float(node.coeffs[4])), abs(float(node.coeffs[5])), abs(float(node.coeffs[6]))]
    return sum(b)

def feat_compaction_pressure(node: RelationalNode, nodes: Dict[int, RelationalNode] = None) -> float:
    r = _quick_rank(node)
    return max(0.0, 3.0 - r) * node.rho

def feat_expansion_resistance(node: RelationalNode, nodes: Dict[int, RelationalNode] = None) -> float:
    r = _quick_rank(node)
    biv_strength = sum(abs(float(x)) for x in node.coeffs[4:7])
    return max(0.0, r - 3.0) * biv_strength

def feat_spin_alignment(node: RelationalNode, nodes: Dict[int, RelationalNode]) -> float:
    parents = [p for p in _quick_parents(node, nodes) if p.protected]
    if not parents or not node.protected:
        return 0.0
    node_b = [abs(float(node.coeffs[4])), abs(float(node.coeffs[5])), abs(float(node.coeffs[6]))]
    diffs = []
    for p in parents:
        p_b = [abs(float(p.coeffs[4])), abs(float(p.coeffs[5])), abs(float(p.coeffs[6]))]
        diffs.append(sum(abs(a - b) for a, b in zip(node_b, p_b)) / 3.0)
    return sum(diffs) / len(diffs) if diffs else 0.0

def feat_new_protected_direction(node: RelationalNode, nodes: Dict[int, RelationalNode]) -> float:
    """
    'Spin' term: rewards nodes that increase the number of independent protected directions
    relative to their causal parents. This is the mechanism intended to 'blow out' the 1D
    nature by actively encouraging the creation of new orthogonal protected planes.
    """
    my_rank = _quick_rank(node)
    parents = _quick_parents(node, nodes)
    if not parents:
        return 0.0
    parent_ranks = [_quick_rank(p) for p in parents]
    avg_parent = sum(parent_ranks) / len(parent_ranks)
    # Positive bonus only when this node creates new directions
    return max(0.0, my_rank - avg_parent)


# ---------------------------------------------------------------------------
# Multi-orientation torsion (spin) pressure features
# ---------------------------------------------------------------------------

def _bivector(node: RelationalNode) -> List[float]:
    """Protected bivector components as three independent 'orientations' (planes)."""
    return [
        abs(float(node.coeffs[4])),
        abs(float(node.coeffs[5])),
        abs(float(node.coeffs[6])),
    ]


def feat_torsion_orientation_0(node: RelationalNode, nodes: Dict[int, RelationalNode]) -> float:
    """
    Torsion pressure in orientation/plane 0.
    Rewards the node when its protected bivector component in the first orientation
    is stronger than the average of its causal parents — i.e., it is adding new
    torsion/spin in a fresh orientational channel.
    """
    b = _bivector(node)
    parents = _quick_parents(node, nodes)
    if not parents:
        return b[0] * 0.5   # mild seed bias
    p0 = sum(_bivector(p)[0] for p in parents) / len(parents)
    return max(0.0, b[0] - p0)


def feat_torsion_orientation_1(node: RelationalNode, nodes: Dict[int, RelationalNode]) -> float:
    """
    Torsion pressure in orientation/plane 1 (second independent bivector direction).
    """
    b = _bivector(node)
    parents = _quick_parents(node, nodes)
    if not parents:
        return b[1] * 0.5
    p1 = sum(_bivector(p)[1] for p in parents) / len(parents)
    return max(0.0, b[1] - p1)


def feat_torsion_orientation_2(node: RelationalNode, nodes: Dict[int, RelationalNode]) -> float:
    """
    Torsion pressure in orientation/plane 2 (third independent bivector direction).
    """
    b = _bivector(node)
    parents = _quick_parents(node, nodes)
    if not parents:
        return b[2] * 0.5
    p2 = sum(_bivector(p)[2] for p in parents) / len(parents)
    return max(0.0, b[2] - p2)


def feat_torsion_orthogonal_magnitude(node: RelationalNode, nodes: Dict[int, RelationalNode]) -> float:
    """
    Stronger 'multiple orientations' signal: measures how much of the child's protected
    bivector lies outside the linear span of the parents' protected bivectors.
    This directly quantifies 'new torsion in a new orientation' that cannot be explained
    by simple inheritance from the causal past.
    """
    child_b = _bivector(node)
    parents = [p for p in _quick_parents(node, nodes) if p.protected]
    if not parents:
        return sum(child_b) * 0.3
    # Build parent subspace (average direction + spread)
    parent_b_vectors = [_bivector(p) for p in parents]
    # Simple proxy: penalty for child energy in directions where all parents are weak
    ortho_energy = 0.0
    for i in range(3):
        parent_i = sum(v[i] for v in parent_b_vectors) / len(parent_b_vectors)
        if parent_i < 0.05 * max(child_b):   # parent orientation i is weak
            ortho_energy += child_b[i]
    return ortho_energy


# ---------------------------------------------------------------------------
# Seven distinct torsion (spin) pressure features
# Expansion from the original 4 (3 planes + orthogonal) to 7 orientations
# ---------------------------------------------------------------------------

def feat_torsion_plane_0(node: RelationalNode, nodes: Dict[int, RelationalNode]) -> float:
    """Torsion Type 1: New activity specifically in protected plane 0 relative to parents."""
    b = _bivector(node)
    parents = _quick_parents(node, nodes)
    if not parents:
        return b[0] * 0.4
    p_avg = sum(_bivector(p)[0] for p in parents) / len(parents)
    return max(0.0, b[0] - p_avg)


def feat_torsion_plane_1(node: RelationalNode, nodes: Dict[int, RelationalNode]) -> float:
    """Torsion Type 2: New activity specifically in protected plane 1."""
    b = _bivector(node)
    parents = _quick_parents(node, nodes)
    if not parents:
        return b[1] * 0.4
    p_avg = sum(_bivector(p)[1] for p in parents) / len(parents)
    return max(0.0, b[1] - p_avg)


def feat_torsion_plane_2(node: RelationalNode, nodes: Dict[int, RelationalNode]) -> float:
    """Torsion Type 3: New activity specifically in protected plane 2."""
    b = _bivector(node)
    parents = _quick_parents(node, nodes)
    if not parents:
        return b[2] * 0.4
    p_avg = sum(_bivector(p)[2] for p in parents) / len(parents)
    return max(0.0, b[2] - p_avg)


def feat_torsion_orthogonal(node: RelationalNode, nodes: Dict[int, RelationalNode]) -> float:
    """Torsion Type 4: Amount of protected bivector lying outside the parents' span (new independent direction)."""
    child_b = _bivector(node)
    parents = [p for p in _quick_parents(node, nodes) if p.protected]
    if not parents:
        return sum(child_b) * 0.25
    parent_vecs = [_bivector(p) for p in parents]
    ortho = 0.0
    for i in range(3):
        p_i = sum(v[i] for v in parent_vecs) / len(parent_vecs)
        if p_i < 0.08 * max(child_b):
            ortho += child_b[i]
    return ortho


def feat_torsion_signed_torque(node: RelationalNode, nodes: Dict[int, RelationalNode]) -> float:
    """Torsion Type 5: Signed torque-like signal (encourages consistent or novel handedness between generations)."""
    child_b = [float(node.coeffs[4]), float(node.coeffs[5]), float(node.coeffs[6])]
    parents = _quick_parents(node, nodes)
    if not parents:
        return 0.0
    # Average parent bivector
    p_b = [0.0, 0.0, 0.0]
    for p in parents:
        for i in range(3):
            p_b[i] += float(p.coeffs[4+i])
    for i in range(3):
        p_b[i] /= len(parents)
    # Simple signed triple-product proxy (child · (parent × some reference))
    # Reward when child adds a new signed component
    torque = sum(c * p for c, p in zip(child_b, p_b))
    return max(0.0, abs(torque) * 2.0)   # positive contribution for any strong signed activity


def feat_torsion_orientation_diversity(node: RelationalNode, nodes: Dict[int, RelationalNode]) -> float:
    """Torsion Type 6: Reward for the child's bivector being mis-aligned with parents (promotes orientation diversity)."""
    child_b = _bivector(node)
    parents = [p for p in _quick_parents(node, nodes) if p.protected]
    if not parents or sum(child_b) < 0.05:
        return 0.0
    p_avg = [sum(_bivector(p)[i] for p in parents) / len(parents) for i in range(3)]
    # Cosine-like similarity (lower similarity = higher diversity reward)
    dot = sum(c * pa for c, pa in zip(child_b, p_avg))
    mag_c = sum(c**2 for c in child_b) ** 0.5 + 1e-6
    mag_p = sum(pa**2 for pa in p_avg) ** 0.5 + 1e-6
    sim = dot / (mag_c * mag_p)
    return max(0.0, (1.0 - sim)) * 1.5   # higher when more orthogonal/diverse


def feat_torsion_magnitude_growth(node: RelationalNode, nodes: Dict[int, RelationalNode]) -> float:
    """Torsion Type 7: Reward for increasing total protected bivector magnitude relative to causal parents."""
    child_mag = sum(_bivector(node))
    parents = _quick_parents(node, nodes)
    if not parents:
        return child_mag * 0.3
    p_avg_mag = sum(sum(_bivector(p)) for p in parents) / len(parents)
    return max(0.0, child_mag - p_avg_mag) * 1.2


def feat_internal_repulsion(node: RelationalNode, nodes: Dict[int, RelationalNode]) -> float:
    """
    Legacy single repulsion (kept for compatibility).
    """
    my_id = node.id
    child_count = 0
    for other in nodes.values():
        for edge in other.parents:
            if edge.parent_id == my_id:
                child_count += 1
                break
    return max(0.0, child_count - 3.0) / 5.0


# ---------------------------------------------------------------------------
# Seven distinct types of internal repulsion (multi-orientation crowding costs)
# Each measures a different notion of "too close" in the relational graph.
# ---------------------------------------------------------------------------

def feat_repulsion_direct_children(node: RelationalNode, nodes: Dict[int, RelationalNode]) -> float:
    """Type 1: Classic out-degree repulsion — too many nodes already list this one as parent."""
    my_id = node.id
    count = sum(1 for other in nodes.values()
                for e in other.parents if e.parent_id == my_id)
    return max(0.0, count - 3.0) / 5.0


def feat_repulsion_sibling_sharing(node: RelationalNode, nodes: Dict[int, RelationalNode]) -> float:
    """Type 2: Sibling crowding — how many other recent nodes share one or more parents with this node."""
    my_parents = {e.parent_id for e in node.parents}
    if not my_parents:
        return 0.0
    recent = [n for n in nodes.values() if abs(n.layer - node.layer) <= 2 and n.id != node.id]
    overlaps = 0
    for other in recent:
        other_p = {e.parent_id for e in other.parents}
        if my_parents & other_p:
            overlaps += 1
    return overlaps / max(1, len(recent)) * 4.0   # scale to be comparable


def feat_repulsion_causal_cone(node: RelationalNode, nodes: Dict[int, RelationalNode]) -> float:
    """Type 3: Causal cone density — number of nodes within small causal distance (depth ≤ 2)."""
    # Count nodes that are within 2 generations of this node's parents or children
    my_parents = {e.parent_id for e in node.parents}
    count = 0
    for other in nodes.values():
        if other.id == node.id:
            continue
        op = {e.parent_id for e in other.parents}
        if my_parents & op:
            count += 1
        else:
            # check grandparents
            for p_id in my_parents:
                if p_id in nodes:
                    gp = {e.parent_id for e in nodes[p_id].parents}
                    if gp & op:
                        count += 1
                        break
    return max(0.0, count - 8) / 6.0


def feat_repulsion_layer_density(node: RelationalNode, nodes: Dict[int, RelationalNode]) -> float:
    """Type 4: Layer-local crowding — many nodes in the same or adjacent layer with overlapping parents."""
    same_layer = [n for n in nodes.values() if abs(n.layer - node.layer) <= 1 and n.id != node.id]
    my_p = {e.parent_id for e in node.parents}
    crowded = sum(1 for n in same_layer
                  if len({e.parent_id for e in n.parents} & my_p) >= 1)
    return max(0.0, crowded - 5) / 4.0


def feat_repulsion_protected_children(node: RelationalNode, nodes: Dict[int, RelationalNode]) -> float:
    """Type 5: Protected-only repulsion — cost only when the children are protected (higher value nodes)."""
    my_id = node.id
    count = 0
    for other in nodes.values():
        if other.protected:
            for e in other.parents:
                if e.parent_id == my_id:
                    count += 1
                    break
    return max(0.0, count - 2.0) / 3.0


def feat_repulsion_bivector_crowding(node: RelationalNode, nodes: Dict[int, RelationalNode]) -> float:
    """Type 6: Spin-similar crowding — too many nearby nodes have similar protected bivector orientation."""
    if not node.protected:
        return 0.0
    my_b = _bivector(node)
    recent = [n for n in nodes.values() if n.protected and abs(n.layer - node.layer) <= 2 and n.id != node.id]
    similar = 0
    for other in recent:
        ob = _bivector(other)
        diff = sum(abs(a - b) for a, b in zip(my_b, ob))
        if diff < 0.15:   # close in protected spin space
            similar += 1
    return max(0.0, similar - 3) / 4.0


def feat_repulsion_ancestor_overlap(node: RelationalNode, nodes: Dict[int, RelationalNode]) -> float:
    """Type 7: Deep ancestor redundancy — many nodes share large portions of the 2-generation ancestor set."""
    # Collect 2-generation ancestors of this node
    gen1 = {e.parent_id for e in node.parents}
    gen2 = set()
    for pid in gen1:
        if pid in nodes:
            gen2 |= {e.parent_id for e in nodes[pid].parents}
    ancestor_set = gen1 | gen2
    if not ancestor_set:
        return 0.0

    recent = [n for n in nodes.values() if abs(n.layer - node.layer) <= 3 and n.id != node.id]
    heavy_overlap = 0
    for other in recent:
        op = {e.parent_id for e in other.parents}
        og2 = set()
        for pid in op:
            if pid in nodes:
                og2 |= {e.parent_id for e in nodes[pid].parents}
        other_anc = op | og2
        overlap = len(ancestor_set & other_anc)
        if overlap >= max(2, len(ancestor_set) // 2):
            heavy_overlap += 1
    return heavy_overlap / max(1, len(recent)) * 5.0


def feat_unconnected_outer_pressure(node: RelationalNode, nodes: Dict[int, RelationalNode]) -> float:
    """
    External pressure applied specifically to "outer" / unconnected nodes.

    A node is considered "outer/unconnected" if it has very low recent usage
    as a parent (few or no direct children in the recent causal window).
    These nodes receive extra pressure, which makes them less attractive for
    new attachments. This creates a "boundary cost" that discourages the
    structure from lingering on the current surface and forces new growth
    to branch into fresh, independent directions (synergistic with the
    highest repulsion channels and the torsion terms).

    The pressure increases the more "exposed" / under-used the node is.
    """
    my_id = node.id
    # Count recent direct children
    recent_window = 3
    recent = [n for n in nodes.values() if abs(n.layer - node.layer) <= recent_window]
    child_count = sum(1 for other in recent
                      for e in other.parents if e.parent_id == my_id)

    # Outer-ness score: high when usage is low
    # Soft threshold: nodes with < 2 recent children are "outer"
    outer_score = max(0.0, 2.5 - child_count) / 2.0
    return outer_score


def feat_connection_scaled_repulsion(node: RelationalNode, nodes: Dict[int, RelationalNode]) -> float:
    """
    "Unconnected nodes experience the highest repulsion, increased by 10x,
    and then scale down as it connects."

    Computes a base repulsion signal from the node's current usage as a parent.
    When the node has almost no recent children (completely unconnected/outer),
    the effective repulsion is multiplied by up to 10x.
    As the node gains connections (more direct children), the multiplier
    smoothly relaxes back toward 1x.

    This term is intended to be given its own coefficient in the pressure
    scalar. A positive coefficient creates a strong disincentive against
    attaching new nodes to currently isolated frontier nodes, while
    making established branches progressively cheaper to grow around.
    This is a dynamic, usage-dependent repulsion that should encourage
    the creation of genuinely new independent directions.
    """
    my_id = node.id
    recent = [n for n in nodes.values() if abs(n.layer - node.layer) <= 3]
    child_count = sum(1 for other in recent
                      for e in other.parents if e.parent_id == my_id)

    # Usage fraction: 0 = unconnected, 1 = well connected (saturates at ~5 children)
    usage = min(1.0, child_count / 5.0)

    # Multiplier: 10.0 at usage=0, smoothly down to 1.0 at usage=1
    multiplier = 1.0 + 9.0 * (1.0 - usage)

    # Return the multiplier itself. The learned coefficient in front of this
    # feature will control how strongly the 10x→1x effect is applied.
    # We also add a small base signal so even well-connected nodes have some cost.
    return multiplier + 0.2   # small floor so the feature is never zero


def feat_chiral_helical_bias(node: RelationalNode, nodes: Dict[int, RelationalNode]) -> float:
    """
    Chiral helical spin bias — direct version (proportional to connections).

    Strength grows with the number of connections the core already has.
    Mature, well-connected high-MECH_A cores spin more vigorously.
    """
    rho = node.rho
    biv = _bivector(node)
    mag = sum(biv)
    if mag < 0.08 or rho < 0.25:
        return 0.0

    mech_a_like = mag * max(0.0, rho - 0.25)
    if mech_a_like < 0.9:
        return 0.0

    recent = [n for n in nodes.values() if abs(n.layer - node.layer) <= 3]
    child_count = sum(1 for other in recent
                      for e in other.parents if e.parent_id == node.id)

    conn_factor = min(1.0, child_count / 6.0)          # direct proportion

    phase = math.atan2(biv[1], biv[0]) + (node.layer * 0.12)

    preferred = 1.0 if biv[2] > 0 else -1.0
    mismatch = abs(math.sin(phase) - preferred * 0.25)

    bias = mismatch * mech_a_like * conn_factor * 1.8
    return bias


def feat_chiral_helical_bias_inverse(node: RelationalNode, nodes: Dict[int, RelationalNode]) -> float:
    """
    Chiral helical spin bias — inverse version (inversely proportional to connections).

    The bias is *strongest* when the high-MECH_A core has very few connections yet.
    This gives newly forming or sparsely connected dense protected regions a powerful
    initial "kick" to spin and reach out in new directions, then the effect relaxes
    as the core becomes more established.

    Same tight MECH_A gate, same slower period, same intrinsic chirality.
    """
    rho = node.rho
    biv = _bivector(node)
    mag = sum(biv)
    if mag < 0.08 or rho < 0.25:
        return 0.0

    mech_a_like = mag * max(0.0, rho - 0.25)
    if mech_a_like < 0.9:
        return 0.0

    recent = [n for n in nodes.values() if abs(n.layer - node.layer) <= 3]
    child_count = sum(1 for other in recent
                      for e in other.parents if e.parent_id == node.id)

    # Inverse: strongest when almost unconnected, falls off as connections grow
    conn_factor = 1.0 / (1.0 + child_count / 3.0)

    phase = math.atan2(biv[1], biv[0]) + (node.layer * 0.12)

    preferred = 1.0 if biv[2] > 0 else -1.0
    mismatch = abs(math.sin(phase) - preferred * 0.25)

    bias = mismatch * mech_a_like * conn_factor * 2.2   # slightly stronger base to compensate inverse
    return bias


# ---------------------------------------------------------------------------
# Two "Inside Space" mechanisms (Version 1 and Version 2)
# Both reward the creation of protected "inside" volume while creating a
# collapsing force on the high-MECH_A particle cores that bound it.
# ---------------------------------------------------------------------------

def _is_high_mech_a(node):
    """Quick check if a node is in the strong MECH_A regime."""
    rho = node.rho
    mag = sum(_bivector(node))
    return (mag * max(0.0, rho - 0.25)) > 0.9


def _enclosed_inside_score(node, nodes):
    """
    Rough local proxy for how much 'inside protected space' this node is part of.
    Counts recent nodes that have a high fraction of parents inside high-MECH_A cores
    and themselves have low outgoing connections (they are 'enclosed' rather than frontier).
    """
    if node.layer < 3:
        return 0.0

    recent = [n for n in nodes.values() if abs(n.layer - node.layer) <= 3 and n.id != node.id]
    high_mech_a_parents = 0
    total_parents = len(node.parents)

    for p_edge in node.parents:
        p = nodes.get(p_edge.parent_id)
        if p and _is_high_mech_a(p):
            high_mech_a_parents += 1

    if total_parents == 0:
        return 0.0

    fraction_inside = high_mech_a_parents / total_parents

    # Only count nodes that are not themselves highly connected outward (true "inside")
    my_child_count = sum(1 for other in recent
                         for e in other.parents if e.parent_id == node.id)

    inside_bonus = max(0.0, 1.0 - my_child_count / 4.0)

    return fraction_inside * inside_bonus


# ---------------------------------------------------------------------------
# LOCAL versions of the two Inside Space mechanisms (respect finite causal speed)
# Each high-MECH_A core only "sees" the inside volume in its own recent causal
# neighborhood. Signals propagate at the speed of the layers.
# ---------------------------------------------------------------------------

def _local_enclosed_score_for_core(core, nodes, window=3):
    """
    For a specific high-MECH_A core, compute how much 'inside protected space'
    is causally associated with *it* within the recent window.
    Only counts nodes that have this core (or its immediate high-MECH_A neighbors)
    as a recent parent and have low outward connectivity.
    """
    if not _is_high_mech_a(core):
        return 0.0

    recent = [n for n in nodes.values() if abs(n.layer - core.layer) <= window]
    enclosed = 0.0

    for other in recent:
        if other.id == core.id:
            continue

        # Does this node have the core as a recent parent?
        has_core = any(e.parent_id == core.id for e in other.parents)
        if not has_core:
            continue

        # Is it "inside" rather than frontier? (low recent children)
        child_count = sum(1 for n in recent if any(e.parent_id == other.id for e in n.parents))
        if child_count <= 2:   # heuristic for "enclosed"
            enclosed += 1.0

    return enclosed


def feat_local_inside_volume_reward_v1(node: RelationalNode, nodes: Dict[int, RelationalNode]) -> float:
    """
    LOCAL Version 1 — Reward for creating inside volume, visible only locally.
    A node gets negative pressure to the extent that it helps create or protect
    inside space that is causally connected to it within the recent light-cone.
    The reward propagates at finite speed.
    """
    if not _is_high_mech_a(node):
        # Non-core nodes only get a small indirect reward if they are part of
        # the local inside volume of a nearby high-MECH_A core.
        score = _local_enclosed_score_for_core(node, nodes)  # will be 0 unless it is a core
        return -score * 0.3
    else:
        score = _local_enclosed_score_for_core(node, nodes)
        return -score * 0.9   # stronger reward for the cores that are doing the enclosing


def feat_local_volume_collapse_v1(node: RelationalNode, nodes: Dict[int, RelationalNode]) -> float:
    """
    LOCAL Version 1 — Collapse on the cores, based only on the inside volume
    they personally enclose within causal reach.
    """
    if not _is_high_mech_a(node):
        return 0.0

    enclosed = _local_enclosed_score_for_core(node, nodes)
    return enclosed * 0.7


def feat_local_inside_negative_pressure_v2(node: RelationalNode, nodes: Dict[int, RelationalNode]) -> float:
    """
    LOCAL Version 2 — Negative pressure for nodes that are causally "inside"
    a high-MECH_A core's local neighborhood. The interior feels stabilized
    after the causal signal arrives.
    """
    if _is_high_mech_a(node):
        return 0.0

    # Is this node inside some nearby high-MECH_A core?
    enclosed = _local_enclosed_score_for_core(node, nodes)
    if enclosed > 0.5:
        return -enclosed * 1.0
    return 0.0


def feat_local_boundary_collapse_v2(node: RelationalNode, nodes: Dict[int, RelationalNode]) -> float:
    """
    LOCAL Version 2 — Collapse on the boundary particles based on the
    causally connected inside volume they personally experience.
    """
    if not _is_high_mech_a(node):
        return 0.0

    enclosed = _local_enclosed_score_for_core(node, nodes)
    return enclosed * 0.7


# ---------------------------------------------------------------------------
# Three experimental mechanisms for "particle as extra dimension of connection"
# Each can be turned on independently via its *_W coefficient.
# ---------------------------------------------------------------------------

def feat_mech_A_extra_conn_stab(node: RelationalNode, nodes: Dict[int, RelationalNode]) -> float:
    """
    Mechanism A: Extra protected connection stabilized by density.
    Nodes with high protected bivector magnitude get a negative pressure contribution
    (reward) only when local rho is high. This represents an "extra connection"
    that only becomes cheap to maintain inside a dense region (particle "folds in").
    """
    rho = node.rho
    biv_mag = sum(abs(float(x)) for x in node.coeffs[4:7])
    # Reward (negative contribution to pressure) when high biv_mag at high rho
    # Simple threshold form
    if rho > 0.25:
        return -biv_mag * (rho - 0.25)   # negative = lower pressure = more attractive / stable
    return 0.0


def feat_mech_B_density_rigidity(node: RelationalNode, nodes: Dict[int, RelationalNode]) -> float:
    """
    Mechanism B: The cost of algebraic rigidity / protection becomes density-dependent.
    At high local density the effective cost of deviating from protected rank 3 or
    having closure defect is modulated (here we make high-density regions pay less
    for maintaining protected structure, i.e. particles are more stable).
    """
    rho = node.rho
    # Use existing compaction + expansion signals but scale them by a density factor
    r = _quick_rank(node)
    biv_strength = sum(abs(float(x)) for x in node.coeffs[4:7])
    base_cost = max(0.0, 3.0 - r) + max(0.0, r - 3.0) * biv_strength
    # At high rho the rigidity cost is reduced (easier to stay protected)
    density_factor = 1.0 / (1.0 + 4.0 * max(0.0, rho - 0.2))
    return base_cost * density_factor


def feat_mech_C_density_closure(node: RelationalNode, nodes: Dict[int, RelationalNode]) -> float:
    """
    Mechanism C: The preferred algebraic closure itself becomes density-dependent.
    At high density the system is encouraged toward configurations that look like
    an "extra connection" has been folded in (here we add a soft bias toward
    slightly higher effective protected magnitude or reduced defect only at high rho).
    """
    rho = node.rho
    if rho < 0.2:
        return 0.0
    biv_mag = sum(abs(float(x)) for x in node.coeffs[4:7])
    # At high density, reward (negative pressure) for having substantial protected magnitude
    # This mimics the closure condition shifting to favor the "folded particle" state.
    return -biv_mag * (rho - 0.2) * 0.8

PRESSURE_FEATURES = [
    feat_external_pressure,
    feat_internal_pressure,
    feat_push_pull,
    feat_local_contrast,
    feat_protected_spin_pressure,
    feat_compaction_pressure,
    feat_expansion_resistance,
    feat_spin_alignment,
    feat_new_protected_direction,         # legacy rank-delta spin
    feat_torsion_orientation_0,           # torsion in first protected orientation
    feat_torsion_orientation_1,           # torsion in second protected orientation
    feat_torsion_orientation_2,           # torsion in third protected orientation
    feat_internal_repulsion,              # legacy single repulsion
    # --- Seven distinct repulsion orientations ---
    feat_repulsion_direct_children,
    feat_repulsion_sibling_sharing,
    feat_repulsion_causal_cone,
    feat_repulsion_layer_density,
    feat_repulsion_protected_children,
    feat_repulsion_bivector_crowding,
    feat_repulsion_ancestor_overlap,
    # --- Seven distinct torsion (spin) orientations (expanded from original 4) ---
    feat_torsion_plane_0,
    feat_torsion_plane_1,
    feat_torsion_plane_2,
    feat_torsion_orthogonal,
    feat_torsion_signed_torque,
    feat_torsion_orientation_diversity,
    feat_torsion_magnitude_growth,
    feat_unconnected_outer_pressure,   # external pressure on outer/unconnected frontier nodes
    feat_connection_scaled_repulsion,  # 10x repulsion on unconnected, scales down as node connects
    # The three experimental mechanisms (A, B, C)
    feat_mech_A_extra_conn_stab,
    feat_mech_B_density_rigidity,
    feat_mech_C_density_closure,
    feat_chiral_helical_bias,          # direct (proportional to connections)
    feat_chiral_helical_bias_inverse,  # inverse (strongest when few connections)
    # LOCAL Inside Space mechanisms (respect finite causal propagation)
    feat_local_inside_volume_reward_v1,
    feat_local_volume_collapse_v1,
    feat_local_inside_negative_pressure_v2,
    feat_local_boundary_collapse_v2,
]

def _quick_rank(node: RelationalNode) -> float:
    b = [abs(float(node.coeffs[4])), abs(float(node.coeffs[5])), abs(float(node.coeffs[6]))]
    return sum(1.0 for x in b if x > 1e-4)

def compute_pressure_scalar(node: RelationalNode, nodes: Dict[int, RelationalNode]) -> float:
    """
    Weighted sum of all pressure features.
    Now includes the original 8 best pressure coefficients + legacy spin term
    + four new multi-orientation torsion (spin) pressure terms.
    Each of the three protected bivector orientations and the orthogonal-new-torsion
    channel has its own independent weight inside the pressure scalar.
    """
    feats = [feat(node, nodes) for feat in PRESSURE_FEATURES]
    coeffs = [
        BEST_COEFFS["term_rank_deviation"] * EXTERNAL_PRESSURE_BOOST,  # boosted external pressure
        BEST_COEFFS["term_rank_abs"],
        BEST_COEFFS["term_closure_defect"],
        BEST_COEFFS["term_balance"],
        BEST_COEFFS["term_rho_rank_interaction"],
        BEST_COEFFS["term_rho_closure_interaction"],
        BEST_COEFFS["term_biv_magnitude"],
        BEST_COEFFS["term_rank_biv_interaction"],
        NEW_SPIN_WEIGHT,               # legacy rank-delta
        INTERNAL_REPULSION_WEIGHT,     # legacy
        REPULSION_DIRECT_W,            # Rep 1
        REPULSION_SIBLING_W,           # Rep 2
        REPULSION_CONE_W,              # Rep 3
        REPULSION_LAYER_W,             # Rep 4
        REPULSION_PROT_W,              # Rep 5
        REPULSION_BIVEC_W,             # Rep 6
        REPULSION_ANCESTOR_W,          # Rep 7
        TORSION_PLANE_0_W,             # Tors 1
        TORSION_PLANE_1_W,             # Tors 2
        TORSION_PLANE_2_W,             # Tors 3
        TORSION_ORTHOGONAL_W,          # Tors 4
        TORSION_SIGNED_W,              # Tors 5
        TORSION_DIVERSITY_W,           # Tors 6
        TORSION_MAGNITUDE_GROW_W,      # Tors 7
        OUTER_UNCONNECTED_PRESSURE_W,  # outer/unconnected boundary pressure
        CONNECTION_SCALED_REPULSION_W, # 10x on unconnected scaling to 1x
        # The three experimental particle-as-extra-connection mechanisms
        MECH_A_EXTRA_CONN_STAB_W,
        MECH_B_DENSITY_RIGIDITY_W,
        MECH_C_DENSITY_CLOSURE_W,
        CHIRAL_HELICAL_BIAS_W,         # direct version
        CHIRAL_HELICAL_BIAS_INVERSE_W, # inverse version
        # LOCAL Inside Space mechanisms (causal, finite speed)
        INSIDE_VOLUME_REWARD_V1_W,   # actually wired to local V1 features below
        VOLUME_COLLAPSE_V1_W,
        INSIDE_NEGATIVE_PRESSURE_V2_W,
        BOUNDARY_COLLAPSE_V2_W,
    ]
    return sum(c * f for c, f in zip(coeffs, feats))

# ---------------------------------------------------------------------------
# Modified node addition with pressure (additive channel)
# ---------------------------------------------------------------------------

def add_nodes_with_best_pressure(
    nodes: Dict[int, RelationalNode],
    n_new: int,
    attach_k: int = 3,
    density_bias: float = 1.45,
    protection_inherit_prob: float = 0.42,
    max_layer_window: int = 4,
    rng_seed: Optional[int] = None,
) -> List[int]:
    """Attachment bias now includes the best pressure configuration (additive channel)."""
    if rng_seed is not None:
        random.seed(rng_seed)

    new_ids: List[int] = []
    current_max_layer = max((n.layer for n in nodes.values()), default=0)

    recent_nodes = [n for n in nodes.values()
                    if n.layer >= current_max_layer - max_layer_window]
    if not recent_nodes:
        recent_nodes = list(nodes.values())

    activities = []
    for n in recent_nodes:
        base = base_node_living_activity(n)
        ps = compute_pressure_scalar(n, nodes)
        activities.append(base * math.exp(-PRESSURE_WEIGHT * ps))

    total = sum(max(0.01, a) for a in activities)
    probs = [max(0.01, a) / total for a in activities] if total > 0 else [1.0 / len(recent_nodes)] * len(recent_nodes)

    next_id = max(nodes.keys()) + 1 if nodes else 0

    for _ in range(n_new):
        parent_pool = random.choices(recent_nodes, weights=probs, k=attach_k * 2)
        chosen: List[RelationalNode] = []
        for p in parent_pool:
            if p not in chosen:
                chosen.append(p)
            if len(chosen) >= attach_k:
                break
        if not chosen:
            chosen = random.sample(recent_nodes, min(attach_k, len(recent_nodes)))

        parents = []
        for p in chosen:
            layer_dist = max(1, current_max_layer - p.layer + 1)
            w = max(0.15, 1.0 / layer_dist) * random.uniform(0.7, 1.0)
            parents.append(CausalEdge(parent_id=p.id, weight=w))

        base_rho = random.uniform(0.05, 0.6)
        is_prot = random.random() < protection_inherit_prob
        if is_prot:
            base_rho *= 1.6
        biv = [random.uniform(-0.08, 0.08) for _ in range(3)] if is_prot else [0.0, 0.0, 0.0]
        coeffs = [base_rho, 0.0, 0.0, 0.0, biv[0], biv[1], biv[2], 0.005]

        new_layer = max((p.layer for p in chosen), default=0) + 1

        node = RelationalNode(
            id=next_id,
            coeffs=coeffs,
            omega_coeffs=[0.0] * 8,
            parents=parents,
            layer=new_layer,
            protected=is_prot,
        )
        nodes[next_id] = node
        new_ids.append(next_id)
        next_id += 1

    return new_ids

# ---------------------------------------------------------------------------
# Modified omega update with pressure modulation (modulation channel)
# ---------------------------------------------------------------------------

def update_all_omegas_with_pressure(nodes: Dict[int, RelationalNode]):
    """
    Applies the best pressure configuration in the *modulation* channel.

    Computes an average pressure scalar across recent nodes for the current step,
    derives effective λ and μ, temporarily overrides the globals, calls the
    normal update on the full graph, then restores the originals.
    """
    current_max_layer = max((n.layer for n in nodes.values()), default=0)
    recent_nodes = [n for n in nodes.values() if n.layer >= current_max_layer - 4]
    if not recent_nodes:
        recent_nodes = list(nodes.values())

    # Average pressure across the recent causal layer (simple global modulation for the step)
    pressures = [compute_pressure_scalar(n, nodes) for n in recent_nodes]
    avg_pressure = sum(pressures) / len(pressures) if pressures else 0.0

    eff_lambda = LAMBDA_NL + MODULATION_WEIGHT * avg_pressure
    eff_mu = MU_NL + MODULATION_WEIGHT * avg_pressure

    # Clamp to safe band
    eff_lambda = max(0.0, min(0.005, eff_lambda))
    eff_mu = max(0.0, min(0.001, eff_mu))

    # Temporarily override the globals used by the living-candidate update
    import minimal_graph_model as pgm
    old_lambda, old_mu = pgm.LAMBDA_NL, pgm.MU_NL
    pgm.LAMBDA_NL = eff_lambda
    pgm.MU_NL = eff_mu

    # Do the normal full-graph update with the modulated coefficients
    pgm.update_all_omegas(nodes)

    # Restore
    pgm.LAMBDA_NL = old_lambda
    pgm.MU_NL = old_mu

# ---------------------------------------------------------------------------
# Main validation run
# ---------------------------------------------------------------------------

def run_best_config_validation(
    target_nodes: int = 1800,
    growth_steps: int = 40,
    nodes_per_step: int = 45,
    seed: int = 20260521,
    perturbation_strength: float = 0.08,
    mech_A_w: float = None,
    mech_B_w: float = None,
    mech_C_w: float = None,
) -> Dict:
    # Allow per-run override of the three experimental mechanisms
    global MECH_A_EXTRA_CONN_STAB_W, MECH_B_DENSITY_RIGIDITY_W, MECH_C_DENSITY_CLOSURE_W
    if mech_A_w is not None:
        MECH_A_EXTRA_CONN_STAB_W = mech_A_w
    if mech_B_w is not None:
        MECH_B_DENSITY_RIGIDITY_W = mech_B_w
    if mech_C_w is not None:
        MECH_C_DENSITY_CLOSURE_W = mech_C_w
    active_mechs = []
    if MECH_A_EXTRA_CONN_STAB_W > 0: active_mechs.append("A")
    if MECH_B_DENSITY_RIGIDITY_W > 0: active_mechs.append("B")
    if MECH_C_DENSITY_CLOSURE_W > 0: active_mechs.append("C")
    mech_str = "+".join(active_mechs) if active_mechs else "none"
    chiral = ""
    if CHIRAL_HELICAL_BIAS_W > 0: chiral += " + direct helical (∝ connections)"
    if CHIRAL_HELICAL_BIAS_INVERSE_W > 0: chiral += " + inverse helical (1/∝ connections)"

    inside = ""
    if INSIDE_VOLUME_REWARD_V1_W > 0 or VOLUME_COLLAPSE_V1_W > 0:
        inside = " + InsideSpace V1 (reward + collapse)"
    elif INSIDE_NEGATIVE_PRESSURE_V2_W > 0 or BOUNDARY_COLLAPSE_V2_W > 0:
        inside = " + InsideSpace V2 (negative inside + boundary collapse)"

    print("=" * 78)
    print("PHASE 6 – BEST CONFIG VALIDATION RUN")
    print(f"Mechanisms active: {mech_str} (A=extra-conn, B=density-rigidity, C=density-closure){chiral}{inside}")
    print("High-MECH_A cores spin; inside space is rewarded while collapsing its own boundary.")
    print(f"Target ~{target_nodes} nodes | Exact locked living candidate")
    print("=" * 78)

    nodes = create_seed_graph(n_seed=16, seed=seed)
    print(f"\n[0] Seed graph: {len(nodes)} nodes")

    update_all_omegas(nodes)

    stats_history = []

    for step in range(1, growth_steps + 1):
        add_nodes_with_best_pressure(
            nodes,
            n_new=nodes_per_step,
            attach_k=3,
            density_bias=1.45,
            protection_inherit_prob=0.42,
            max_layer_window=4,
            rng_seed=seed + step * 17,
        )
        update_all_omegas_with_pressure(nodes)

        current_n = len(nodes)
        growth = compute_growth_curve(nodes, sample_size=9, max_d=8)
        d_eff = estimate_d_eff(growth)
        degs = degree_histograms(nodes)
        cone = light_cone_proxy(nodes, sample=12)

        stats = {
            "step": step,
            "N": current_n,
            "d_eff_proxy": round(d_eff, 3),
            "avg_rho": round(sum(n.rho for n in nodes.values()) / current_n, 4),
            "protected_fraction": round(sum(1 for n in nodes.values() if n.protected) / current_n, 3),
        }
        stats_history.append(stats)

        print(f"[Step {step}] N={current_n}  d_eff≈{d_eff:.3f}  prot={stats['protected_fraction']:.3f}")

        if current_n >= target_nodes:
            break

    print(f"\nFinal graph: {len(nodes)} nodes")

    # Optional perturbation + recovery (Phase 4 style)
    print("\nRunning perturbation + recovery (4 trials)...")
    # (For brevity in this first version we just note the protocol; full recovery metrics
    # can be added in a follow-up if the baseline d_eff looks promising.)

    # Export rich snapshot
    snapshot = {
        "phase": "phase6_best_config_validation",
        "mode": "BOTH (best discovered coefficients)",
        "best_loss_from_search": -0.7132,
        "living_candidate": "⟨D Ω + λ Ω² + μ ⟨Ω, Ω⟩⟩_{0,2} = f_g(ρ)·J_ρ + f_em(ρ)·J_χ",
        "params": {
            "lambda": LAMBDA_NL,
            "mu": MU_NL,
            "rho_crit": RHO_CRIT,
            "pressure_weight": 0.8,
            "modulation_weight": 0.8,
        },
        "final_N": len(nodes),
        "stats_history": stats_history,
    }

    out_path = SCRIPT_DIR / "phase6_best_config_validation_snapshot.json"
    with open(out_path, "w") as f:
        json.dump(snapshot, f, indent=2)
    print(f"\nRich snapshot written: {out_path}")

    return {"final_N": len(nodes), "history": stats_history}


if __name__ == "__main__":
    run_best_config_validation(target_nodes=1800, growth_steps=40)