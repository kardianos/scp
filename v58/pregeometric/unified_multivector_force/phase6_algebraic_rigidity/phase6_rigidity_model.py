#!/usr/bin/env python3
"""
phase6_rigidity_model.py

Phase 6 – Algebraic Rigidity and Multivector Closure Test Harness

This file implements the first concrete test of the hypothesis that
an algebraic rigidity (closure) condition on the protected bivector
subalgebra favors exactly three independent directions as the stable,
self-reinforcing configuration for high-density protected modes.

We extend the Phase 1 / Phase 4 relational graph model with a simple
but non-trivial "rigidity cost" and bias node attachment/growth toward
nodes whose local protected structure is close to the rigid-at-3
configuration.

The goal is to observe whether including this term produces:
- Higher and more stable effective dimensionality
- Stronger self-stabilizing attractor behavior
- Better recovery after perturbation

All evolution uses the *exact* locked living candidate with no changes
to parameters.

Run:
    python phase6_rigidity_model.py

This is the "burn down the next 3 steps" implementation:
1. Implement rigidity cost + closure defect
2. Integrate into growth/attachment rule
3. Run experiment + export data for Lean
"""

from __future__ import annotations
import sys
import math
import random
from pathlib import Path
from typing import Dict, List, Optional
import json

# ---------------------------------------------------------------------------
# Reuse the Phase 1 canonical implementation (read-only)
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
    node_living_activity as base_node_living_activity,
    add_nodes_density_biased as base_add_nodes_density_biased,
    compute_growth_curve,
    estimate_d_eff,
    light_cone_proxy,
    degree_histograms,
    causal_past_ball_list,
    LAMBDA_NL,
    MU_NL,
    RHO_CRIT,
    FG_DEFAULT_BACKGROUND,
)

# ---------------------------------------------------------------------------
# Phase 6: Algebraic Rigidity / Closure Cost
# ---------------------------------------------------------------------------

def compute_protected_rigidity(node: RelationalNode) -> float:
    """
    First-cut algebraic rigidity (closure) measure on the protected bivector structure.

    We look at the three protected bivector components (the ones that carry
    independent directional information in the model).

    Returns a cost: lower is better (closer to "rigid at exactly 3 closed directions").

    Current heuristic (to be refined with real multivector closure identities later):
    - Count how many of the three bivector blades are "active".
    - Penalize strongly for > 3 or < 2 active directions.
    - Add a simple "closure quality" term that rewards configurations where
      the three directions have comparable magnitude (balanced, "closed" truss).
    """
    # The three protected bivector coefficients live in indices 4,5,6 of coeffs
    # (e12, e13, e23) in the 8-component layout.
    b0 = abs(float(node.coeffs[4]))   # e12
    b1 = abs(float(node.coeffs[5]))   # e13
    b2 = abs(float(node.coeffs[6]))   # e23

    blades = [b0, b1, b2]
    active = [b for b in blades if b > 1e-4]

    r = len(active)  # naive rank

    # Base penalty for wrong number of directions
    if r == 3:
        base = 0.0
    elif r == 2:
        base = 1.5
    elif r == 1:
        base = 4.0
    else:
        base = 6.0

    # Closure quality bonus (balanced magnitudes are "closed" and rigid)
    if r >= 2:
        s = sum(active)
        if s > 0:
            balance = 1.0 - (max(active) - min(active)) / (s + 1e-9)
            closure_bonus = 1.0 - balance   # 0 when perfectly balanced
        else:
            closure_bonus = 1.0
    else:
        closure_bonus = 2.0

    cost = base + 1.2 * closure_bonus

    # Extra penalty if the node is protected but has bad structure
    if node.protected and r != 3:
        cost += 1.5

    return cost


def node_living_activity_with_rigidity(node: RelationalNode, alpha: float = 0.8) -> float:
    """
    Phase 6 living activity = original Phase 1 activity * rigidity bonus.

    Lower rigidity cost → higher activity → more likely to attract new nodes.
    This is the concrete embodiment of "extremize protected density subject to
    algebraic rigidity at 3".
    """
    base = base_node_living_activity(node)
    cost = compute_protected_rigidity(node)
    # Convert cost to a bonus (exponential decay is stable)
    rigidity_bonus = math.exp(-alpha * cost)
    return base * rigidity_bonus


def add_nodes_with_rigidity(
    nodes: Dict[int, RelationalNode],
    n_new: int,
    attach_k: int = 3,
    density_bias: float = 1.45,
    protection_inherit_prob: float = 0.42,
    max_layer_window: int = 4,
    rigidity_weight: float = 0.8,
    rng_seed: Optional[int] = None,
) -> List[int]:
    """
    Phase 6 version of add_nodes that includes the algebraic rigidity cost
    in the attachment bias (in addition to ρ + Ω activity).
    """
    if rng_seed is not None:
        random.seed(rng_seed)

    new_ids: List[int] = []
    current_max_layer = max((n.layer for n in nodes.values()), default=0)

    recent_nodes = [n for n in nodes.values()
                    if n.layer >= current_max_layer - max_layer_window]
    if not recent_nodes:
        recent_nodes = list(nodes.values())

    # Combined activity: original living activity * rigidity factor
    activities = []
    for n in recent_nodes:
        act = node_living_activity_with_rigidity(n, alpha=rigidity_weight)
        activities.append(max(0.01, act))

    total = sum(a ** density_bias for a in activities)
    if total <= 0:
        probs = [1.0 / len(recent_nodes)] * len(recent_nodes)
    else:
        probs = [(a ** density_bias) / total for a in activities]

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
# Phase 6 Demonstration Runner
# ---------------------------------------------------------------------------

def run_phase6_rigidity_demo(
    target_nodes: int = 420,
    growth_steps: int = 9,
    nodes_per_step: int = 45,
    rigidity_weight: float = 0.9,
    seed: int = 20260520,
) -> Dict:
    """
    Run a Phase 6 experiment with the algebraic rigidity bias included in growth.
    Compare (qualitatively in this first run) against the Phase 4 baseline behavior.
    """
    print("=" * 75)
    print("PHASE 6 – ALGEBRAIC RIGIDITY TEST (first implementation)")
    print("Living candidate (exact locked) + protected rigidity cost in attachment")
    print(f"Target ~{target_nodes} nodes | rigidity_weight={rigidity_weight}")
    print("=" * 75)

    nodes = create_seed_graph(n_seed=16, seed=seed)
    print(f"\n[0] Seed graph: {len(nodes)} nodes")

    update_all_omegas(nodes)

    stats_history = []

    for step in range(1, growth_steps + 1):
        added = add_nodes_with_rigidity(
            nodes,
            n_new=nodes_per_step,
            attach_k=3,
            density_bias=1.45,
            protection_inherit_prob=0.42,
            max_layer_window=4,
            rigidity_weight=rigidity_weight,
            rng_seed=seed + step * 31,
        )
        update_all_omegas(nodes)

        current_n = len(nodes)
        growth = compute_growth_curve(nodes, sample_size=8, max_d=8)
        d_eff = estimate_d_eff(growth)
        degs = degree_histograms(nodes)
        cone = light_cone_proxy(nodes, sample=12)

        # Compute average rigidity cost for diagnostics
        rig_costs = [compute_protected_rigidity(n) for n in nodes.values()]
        avg_rig = sum(rig_costs) / len(rig_costs) if rig_costs else 0.0

        stats = {
            "step": step,
            "N": current_n,
            "d_eff_proxy": round(d_eff, 3),
            "avg_rho": round(sum(n.rho for n in nodes.values()) / current_n, 4),
            "protected_fraction": round(sum(1 for n in nodes.values() if n.protected) / current_n, 3),
            "avg_rigidity_cost": round(avg_rig, 3),
            "light_cone_mean_branch": round(cone["mean_branch"], 3),
        }
        stats_history.append(stats)

        print(f"[Step {step}] N={current_n}  d_eff≈{d_eff:.3f}  "
              f"prot={stats['protected_fraction']:.3f}  "
              f"avg_rig_cost={avg_rig:.3f}")

        if current_n >= target_nodes:
            break

    print("\nFinal graph:", len(nodes), "nodes")
    print("Average final rigidity cost:", round(sum(compute_protected_rigidity(n) for n in nodes.values()) / len(nodes), 3))

    # Export a Lean-friendly snapshot
    snapshot = {
        "phase": "phase6_rigidity_first_run",
        "living_candidate": "⟨D Ω + λ Ω² + μ ⟨Ω, Ω⟩⟩_{0,2} = f_g(ρ)·J_ρ + f_em(ρ)·J_χ",
        "params": {"lambda": LAMBDA_NL, "mu": MU_NL, "rho_crit": RHO_CRIT, "rigidity_weight": rigidity_weight},
        "final_N": len(nodes),
        "stats_history": stats_history,
    }

    out_path = SCRIPT_DIR / "phase6_first_run_snapshot.json"
    with open(out_path, "w") as f:
        json.dump(snapshot, f, indent=2)
    print(f"\nLean-ingestible snapshot written: {out_path}")

    return {"final_N": len(nodes), "history": stats_history}


if __name__ == "__main__":
    run_phase6_rigidity_demo(target_nodes=420, growth_steps=9)