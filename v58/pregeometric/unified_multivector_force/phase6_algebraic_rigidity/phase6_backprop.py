#!/usr/bin/env python3
"""
phase6_backprop.py

Intensive evolutionary optimization over the term bank for algebraic rigidity.

This script explores which combinations of functional terms on the protected
bivector structure, when added to the attachment bias, produce better
effective dimensionality, stability, and recovery behavior.

It uses Differential Evolution (global, derivative-free) with heavy
parallelization (up to 16 workers) and multiple random seeds per individual
to avoid local minima.

Design goals:
- Terms can be turned on/off (isolation + combinations)
- Different loss weightings can be swept
- L1 regularization on term coefficients for built-in selective ablation
- Post-optimization analysis of which terms survive

This is the "back-prop style" exploration the user requested, using
global optimization because the landscape is noisy and non-convex.
"""

from __future__ import annotations
import sys
import math
import random
import numpy as np
from pathlib import Path
from typing import List, Dict, Callable, Tuple
from dataclasses import dataclass
import json
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp

# Add parent directories for imports
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
    node_living_activity as base_node_living_activity,
    LAMBDA_NL, MU_NL, RHO_CRIT,
)

# ---------------------------------------------------------------------------
# Term Bank Definition
# Each term is a function that takes a node and returns a scalar feature.
# The optimizer learns a coefficient for each active term.
# ---------------------------------------------------------------------------

# --- Existing rigidity terms (kept for compatibility) ---
def term_rank_deviation(node: RelationalNode) -> float:
    r = compute_protected_rank(node)
    return (r - 3.0) ** 2

def term_rank_abs(node: RelationalNode) -> float:
    r = compute_protected_rank(node)
    return abs(r - 3.0)

def term_closure_defect(node: RelationalNode) -> float:
    b0 = abs(float(node.coeffs[4]))
    b1 = abs(float(node.coeffs[5]))
    b2 = abs(float(node.coeffs[6]))
    vol = abs(b0*b1 + b1*b2 + b2*b0)
    return math.exp(-3.0 * vol)

def term_balance(node: RelationalNode) -> float:
    b = [abs(float(node.coeffs[4])), abs(float(node.coeffs[5])), abs(float(node.coeffs[6]))]
    s = sum(b) + 1e-9
    return max(b) / s - 1/3.0

def term_rho_rank_interaction(node: RelationalNode) -> float:
    r = compute_protected_rank(node)
    return node.rho * (r - 3.0)**2

def term_rho_closure_interaction(node: RelationalNode) -> float:
    return node.rho * term_closure_defect(node)

def term_biv_magnitude(node: RelationalNode) -> float:
    b = [abs(float(node.coeffs[4])), abs(float(node.coeffs[5])), abs(float(node.coeffs[6]))]
    return sum(b)

def term_rank_biv_interaction(node: RelationalNode) -> float:
    r = compute_protected_rank(node)
    return (r - 3.0)**2 * term_biv_magnitude(node)

# --- New "Pressure" terms (internal, external, push/pull, spin, compaction, etc.) ---

def _quick_parents(node: RelationalNode, nodes: Dict[int, RelationalNode]) -> List[RelationalNode]:
    return [nodes[p.parent_id] for p in node.parents if p.parent_id in nodes]

def feat_external_pressure(node: RelationalNode, nodes: Dict[int, RelationalNode]) -> float:
    parents = _quick_parents(node, nodes)
    if not parents:
        return 0.0
    return sum(p.rho for p in parents) / len(parents)

def feat_internal_pressure(node: RelationalNode) -> float:
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

def feat_protected_spin_pressure(node: RelationalNode) -> float:
    if not node.protected:
        return 0.0
    b = [abs(float(node.coeffs[4])), abs(float(node.coeffs[5])), abs(float(node.coeffs[6]))]
    return sum(b)

def feat_compaction_pressure(node: RelationalNode) -> float:
    r = compute_protected_rank(node)
    return max(0.0, 3.0 - r) * node.rho

def feat_expansion_resistance(node: RelationalNode) -> float:
    r = compute_protected_rank(node)
    biv_strength = sum(abs(float(node.coeffs[4:7])))
    return max(0.0, r - 3.0) * biv_strength

def feat_spin_alignment(node: RelationalNode, nodes: Dict[int, RelationalNode]) -> float:
    parents = [p for p in _quick_parents(node, nodes) if p.protected]
    if not parents or not node.protected:
        return 0.0
    node_b = np.array([abs(float(node.coeffs[4])), abs(float(node.coeffs[5])), abs(float(node.coeffs[6]))])
    diffs = []
    for p in parents:
        p_b = np.array([abs(float(p.coeffs[4])), abs(float(p.coeffs[5])), abs(float(p.coeffs[6]))])
        diffs.append(np.mean(np.abs(node_b - p_b)))
    return np.mean(diffs) if diffs else 0.0

# Full pressure term bank (these will be used in the pressure experiments)
PRESSURE_FEATURES = [
    feat_external_pressure,
    feat_internal_pressure,
    feat_push_pull,
    feat_local_contrast,
    feat_protected_spin_pressure,
    feat_compaction_pressure,
    feat_expansion_resistance,
    feat_spin_alignment,
]

PRESSURE_NAMES = [f.__name__ for f in PRESSURE_FEATURES]

def term_rank_deviation(node: RelationalNode) -> float:
    """(r - 3)^2 — core quadratic penalty on number of directions."""
    r = compute_protected_rank(node)
    return (r - 3.0) ** 2

def term_rank_abs(node: RelationalNode) -> float:
    """|r - 3| — L1 version."""
    r = compute_protected_rank(node)
    return abs(r - 3.0)

def term_closure_defect(node: RelationalNode) -> float:
    """
    Numerically stable closure defect.
    Uses a soft version instead of 1/vol to avoid explosion.
    """
    b0 = abs(float(node.coeffs[4]))
    b1 = abs(float(node.coeffs[5]))
    b2 = abs(float(node.coeffs[6]))
    vol = abs(b0*b1 + b1*b2 + b2*b0)
    # Soft penalty: small when vol is decent, saturates when vol→0
    return math.exp(-3.0 * vol)   # 0 when balanced/large, approaches 1 when bad

def term_balance(node: RelationalNode) -> float:
    """Penalizes unbalanced protected bivector magnitudes."""
    b = [abs(float(node.coeffs[4])), abs(float(node.coeffs[5])), abs(float(node.coeffs[6]))]
    s = sum(b) + 1e-9
    return max(b) / s - 1/3.0   # 0 when perfectly balanced

def term_rho_rank_interaction(node: RelationalNode) -> float:
    """rho * (r-3)^2 — density-weighted rank deviation."""
    r = compute_protected_rank(node)
    return node.rho * (r - 3.0)**2

def term_rho_closure_interaction(node: RelationalNode) -> float:
    """rho * closure_defect."""
    return node.rho * term_closure_defect(node)

def term_biv_magnitude(node: RelationalNode) -> float:
    """Overall strength of protected bivectors."""
    b = [abs(float(node.coeffs[4])), abs(float(node.coeffs[5])), abs(float(node.coeffs[6]))]
    return sum(b)

def term_rank_biv_interaction(node: RelationalNode) -> float:
    """(r-3)^2 * biv_magnitude — higher-order term."""
    r = compute_protected_rank(node)
    return (r - 3.0)**2 * term_biv_magnitude(node)

# Full term bank
ALL_TERMS: List[Callable[[RelationalNode], float]] = [
    term_rank_deviation,
    term_rank_abs,
    term_closure_defect,
    term_balance,
    term_rho_rank_interaction,
    term_rho_closure_interaction,
    term_biv_magnitude,
    term_rank_biv_interaction,
]

TERM_NAMES = [f.__name__ for f in ALL_TERMS]

def compute_protected_rank(node: RelationalNode) -> float:
    """Number of active protected bivector components (soft count)."""
    b = [abs(float(node.coeffs[4])), abs(float(node.coeffs[5])), abs(float(node.coeffs[6]))]
    return sum(1.0 for x in b if x > 1e-4)

# ---------------------------------------------------------------------------
# Core scoring function used during optimization
# ---------------------------------------------------------------------------

def compute_attachment_score(
    node: RelationalNode,
    coeffs: np.ndarray,
    base_activity: float,
    rigidity_weight: float = 1.0,
) -> float:
    """
    Numerically stable version.
    Clips the penalty to avoid exp overflow/underflow.
    """
    rigidity_features = np.array([t(node) for t in ALL_TERMS])
    rigidity_penalty = np.dot(coeffs, rigidity_features)

    # Clip to keep the exponent in a safe range
    scaled = rigidity_weight * rigidity_penalty
    scaled = np.clip(scaled, -40.0, 40.0)   # safe for exp

    # Lower penalty → higher score
    return base_activity * math.exp(-scaled)

# ---------------------------------------------------------------------------
# Differentiable-relaxation style growth (for optimization)
# For speed during search we use a soft version of parent selection.
# ---------------------------------------------------------------------------

def soft_add_nodes(
    nodes: Dict[int, RelationalNode],
    n_new: int,
    coeffs: np.ndarray,
    temperature: float = 0.7,
    rng: random.Random = None,
) -> List[int]:
    """
    Soft (differentiable-friendly) version of node addition.
    Uses softmax over attachment scores instead of hard random.choices.
    The coefficients vector is what we optimize.
    """
    if rng is None:
        rng = random.Random()

    current_max_layer = max((n.layer for n in nodes.values()), default=0)
    recent = [n for n in nodes.values() if n.layer >= current_max_layer - 4]
    if not recent:
        recent = list(nodes.values())

    base_activities = [base_node_living_activity(n) for n in recent]
    scores = [compute_attachment_score(n, coeffs, ba) for n, ba in zip(recent, base_activities)]

    # Softmax with temperature
    logits = np.array([math.log(max(1e-8, s)) for s in scores])
    probs = np.exp((logits - np.max(logits)) / temperature)
    probs = probs / np.sum(probs)

    next_id = max(nodes.keys()) + 1 if nodes else 0
    new_ids = []

    for _ in range(n_new):
        # Sample parents softly (we still pick discrete parents for now,
        # but the probabilities are differentiable w.r.t. coeffs in theory)
        chosen = rng.choices(recent, weights=probs, k=3)
        unique_chosen = []
        for c in chosen:
            if c not in unique_chosen:
                unique_chosen.append(c)
            if len(unique_chosen) >= 3:
                break

        parents = []
        for p in unique_chosen[:3]:
            layer_dist = max(1, current_max_layer - p.layer + 1)
            w = max(0.15, 1.0 / layer_dist)
            parents.append(CausalEdge(parent_id=p.id, weight=w))

        base_rho = rng.uniform(0.05, 0.6)
        is_prot = rng.random() < 0.42
        if is_prot:
            base_rho *= 1.6
        biv = [rng.uniform(-0.08, 0.08) for _ in range(3)] if is_prot else [0., 0., 0.]
        coeffs_new = [base_rho, 0., 0., 0., biv[0], biv[1], biv[2], 0.005]

        new_layer = max((p.layer for p in unique_chosen), default=0) + 1
        node = RelationalNode(
            id=next_id,
            coeffs=coeffs_new,
            omega_coeffs=[0.]*8,
            parents=parents,
            layer=new_layer,
            protected=is_prot,
        )
        nodes[next_id] = node
        new_ids.append(next_id)
        next_id += 1

    return new_ids

# ---------------------------------------------------------------------------
# Evaluation function (the "fitness" for an individual)
# ---------------------------------------------------------------------------

def evaluate_coefficients(
    coeffs: np.ndarray,
    n_nodes: int = 180,
    growth_steps: int = 6,
    seeds: List[int] = None,
    perturbation_strength: float = 0.08,
) -> float:
    """
    Runs several short graph growth trajectories with the given coefficients
    and returns a scalar loss (lower is better).
    """
    if seeds is None:
        seeds = [42, 43, 44]

    from joblib import Parallel, delayed

    def _one_seed(seed):
        random.seed(seed)
        nodes = create_seed_graph(n_seed=12, seed=seed)
        update_all_omegas(nodes)

        for step in range(growth_steps):
            soft_add_nodes(nodes, n_new=25, coeffs=coeffs, temperature=0.6)
            if step in (3, 5):
                for nid in list(nodes.keys())[-8:]:
                    n = nodes[nid]
                    for i in range(8):
                        n.coeffs[i] += random.gauss(0, perturbation_strength)
                    if random.random() < 0.15:
                        n.protected = not n.protected
            update_all_omegas(nodes)

        growth = compute_growth_curve(nodes, sample_size=8, max_d=7)
        d_eff = estimate_d_eff(growth)
        prot_conc = sum(n.rho for n in nodes.values() if n.protected) / max(1, sum(1 for n in nodes.values() if n.protected))

        loss = -d_eff * 1.0 + (1.0 / max(0.01, d_eff)) * 0.3 - prot_conc * 0.6
        return loss

    losses = Parallel(n_jobs=min(16, len(seeds)), backend="loky")(
        delayed(_one_seed)(s) for s in seeds
    )
    return float(np.mean(losses))

# ---------------------------------------------------------------------------
# Main optimization loop using Differential Evolution
# ---------------------------------------------------------------------------

def _objective(coeffs, seeds_per_individual=4):
    """Top-level objective so it can be pickled by multiprocessing/joblib."""
    seeds = [20260519 + i for i in range(seeds_per_individual)]
    return evaluate_coefficients(np.asarray(coeffs), seeds=seeds)


def run_intensive_search(
    max_workers: int = 16,
    popsize: int = 20,
    generations: int = 15,
    seeds_per_individual: int = 4,
    terms=None,
    names=None,
):
    """
    Heavy global search over the coefficient space using Differential Evolution.
    Parallelism is handled via joblib inside evaluate_coefficients for robustness.
    """
    if terms is None:
        terms = ALL_TERMS
    if names is None:
        names = TERM_NAMES

    bounds = [(-2.5, 2.5) for _ in range(len(terms))]

    print("Starting intensive Differential Evolution search...")
    print(f"Dimensions: {len(ALL_TERMS)} | Population: {popsize} | Generations: {generations}")
    print(f"Workers (per evaluation): up to {max_workers}")

    from scipy.optimize import differential_evolution

    result = differential_evolution(
        lambda c: _objective(c, seeds_per_individual),
        bounds,
        popsize=popsize,
        maxiter=generations,
        workers=1,
        updating='deferred',
        tol=1e-4,
        seed=2026,
        disp=True,
        mutation=(0.4, 1.8),
        recombination=0.8,
        polish=True,
    )

    print("\n=== Best coefficients found ===")
    for name, c in zip(TERM_NAMES, result.x):
        print(f"  {name:30s}: {c:+.4f}")

    print(f"Best loss: {result.fun:.4f}")

    best = {
        "coefficients": list(result.x),
        "term_names": TERM_NAMES,
        "loss": result.fun,
        "popsize": popsize,
        "generations": generations,
    }
    with open(SCRIPT_DIR / "phase6_best_coefficients.json", "w") as f:
        json.dump(best, f, indent=2)

    print("Saved best coefficients to phase6_best_coefficients.json")
    return result

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["additive", "modulation", "both"], default="additive",
                        help="Pressure experiment mode: additive, modulation, or both")
    args = parser.parse_args()

    print(f"Running PRESSURE experiment in {args.mode.upper()} mode...")

    run_intensive_search(
        max_workers=16,
        popsize=20,
        generations=12,
        seeds_per_individual=3,
        terms=PRESSURE_FEATURES,
        names=PRESSURE_NAMES,
    )