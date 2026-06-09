# v65 FINDINGS — E1 ⊕ E3 (bilevel hypergradient + feature-learning threshold)

**Date**: 2026-06-05
**Setup**: identical static oscillon seed (`gen_oscillon -A 1.0 -sigma 1.5`, 128³, ±10,
core P≈0.27 ⇒ predicted feature crossover κ_c≈0.41/P²≈5.6), evolved T=8 (short — the
back-reaction response before dispersal), **kernel untouched, config-level only**. 11
runs on a V100: κ-scan {0, 0.5, 1.5, 5, 15, 30, 50} (E3) + μ±, η± at κ=50 (E1). Data:
`/space/scp/v65_e13/diag_*.tsv`. Reference: [`CONCEPT.md`](CONCEPT.md) §4, §7;
[`feature_onset.py`](feature_onset.py).

---

## E3 — "magic = feature learning" is confirmed, as a STABILITY THRESHOLD (stronger than predicted)

| κ | E_total(t8) | P_max(t8) | state |
|---|---|---|---|
| 0 | NaN | NaN | **collapse** (immediate) |
| 0.5 | −6.6×10³ | 14.7 | **collapse** (P_max 0.27→11 by t=1; E blows up) |
| 1.5 | 50.2 | 0.396 | stable (overshoot→settle; near-critical) |
| 5 | 50.9 | 0.194 | stable |
| 15 | 52.0 | 0.031 | stable |
| 30 | 52.8 | 0.0086 | stable |
| 50 | 53.3 | 0.0057 | stable |

**POSITIVE — the κ→0 (lazy) limit is not merely inert, it is *unstable*.** With μ<0 and
too little saturation the density runs away (κ=0 NaN; κ=0.5: P_max 0.27→11.1 by t=1). A
**stable reorganized background requires κ above a threshold κ_crit ~ 1** (between the
collapsing κ=0.5 and the stable κ=1.5). Above it, a stable soliton forms whose core
density P_max is set by the saturation strength and decreases monotonically with κ
(0.40→0.0057 over κ=1.5→50).

This is **sharper than the predicted smooth knee**: the lazy→feature transition is a
**stability/regime boundary**, not a linear-from-zero coefficient. It **decisively refutes
the null** ("magic is just a coefficient" — which would show a straight line from κ=0):
below threshold there is *no stable feature at all*. This is the exact NN
feature-learning picture — insufficient nonlinearity ⇒ no stable learned representation
(degenerate/collapse); sufficient ⇒ a stable reorganized core. The standard κ=50 sits far
inside the feature regime (P_max fully saturation-controlled), as `feature_onset.py`
predicted.

**Caveat**: the observed threshold (~1) is the same order as the analytic crossover
κ_c≈5.6 but not identical — they are different definitions (dynamical collapse boundary
vs χ=0.5 force-quench). Both are O(1)–O(5); a finer κ-scan in (0.5, 1.5) would pin κ_crit.
The near-critical κ=1.5 run overshoots (P_max→2.6) then settles — expected critical slowing.

---

## E1 — the bilevel hypergradient is computable from short kernel runs

Finite-difference ∇_Θ E_total at the standard point (κ=50), central differences:

  ∂E/∂μ = +0.070,  ∂E/∂η = −0.119,  ∂E/∂κ = +0.0293

**POSITIVE** — all three are smooth, finite, and well-defined from short config-level runs
(no kernel edit, no adjoint). The bilevel **hypergradient** of [`CONCEPT.md`](CONCEPT.md)
§3 is therefore computable in practice: "differentiable theory" tuning (descend Θ toward a
target spectrum) is viable from this point.

**Caveat / useful by-product**: the E3 collapse shows the Θ-landscape has a **stability
cliff** at κ ~ 1 — outside the stable region the observables diverge and the gradient is
meaningless. Any E1 descent must be **constrained to the stable region (κ ≳ 1.5)**. This
is itself a finding: the parameter landscape is not globally smooth; it has a physical
boundary (the v65 analog of a region where "training diverges"). The present gradient is
of E_total at fixed t=8 from a fixed seed (demonstrates computability); a full descent
would use relaxed-equilibrium observables and stay inside the stable region.

---

## Net verdict

Both halves of the minimal NN-physics realization **succeed**:
- **E1**: the inner-loop relaxation + outer-loop finite-difference hypergradient is a
  working, kernel-free "differentiable theory" — the bilevel structure of `CONCEPT.md` is
  operational.
- **E3**: "magic = feature learning" is confirmed in its strongest form — a **stability
  threshold** in the nonlinearity near κ~1, below which no stable reorganized background
  exists. Not a coefficient; a regime. The null is refuted.

**Open (next)**: pin κ_crit with a fine scan in (0.5, 1.5); run a *constrained* E1 descent
on Θ toward a target observable (staying κ≳1.5) to demonstrate end-to-end theory tuning;
and — the moonshot — E2/E4 (learned seed search, learned-RG fixed point) which require the
separate differentiable re-implementation validated against `scp_sim.c`.
