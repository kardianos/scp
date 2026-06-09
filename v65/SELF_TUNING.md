# v65 — Intrinsic Self-Tuning: the Optimizer Inside the Action

**Date**: 2026-06-05
**Supersedes the E1 framing.** E1 treated the Lagrangian parameters `Θ` as knobs tuned by
an *external* finite-difference loop. The goal here is the opposite: a mechanism **inside
the physics** by which the theory relaxes itself to a point that is simultaneously
**dynamically stable** and **stationary in the action** — no outside optimizer.
**Parents**: [`CONCEPT.md`](CONCEPT.md), [`FINDINGS_E13.md`](FINDINGS_E13.md) (the measured
stability threshold κ_crit ~ 1), MEMORY `V6 — Conserved Torsion Field`.

---

## 1. The goal, precisely

Promote (a subset of) the couplings `Θ = {κ, μ, η, …}` from constants to **slow dynamical
degrees of freedom** that flow on the *same* action as the fields, with a large timescale
separation:

  fast:  ∂Φ/∂τ = −δE/δΦ        (the existing soliton relaxation)
  slow:  Θ̇    = −ε ∂E/∂Θ        (ε ≪ 1; the couplings relax too)

The desired fixed point is **self-consistent stationarity** — `δE/δΦ = 0` (a soliton)
*and* `∂E/∂Θ = 0` (self-tuned couplings) — reached autonomously. This is the
"optimizer-in-the-physics": one descent on one energy, fast and slow weights together
(the equilibrium-propagation / energy-based-learning picture, made physical).

---

## 2. The obstruction (verified) — action alone selects collapse

Naive slow flow `Θ̇ = −ε ∂E/∂Θ` **does not work for the stabilizing coupling.** For
`V(P) = (μ/2)P²/(1+κP²)`,

  ∂V/∂κ = −(μ/2) P⁴/(1+κP²)²  > 0   for all P, since μ<0,

so `∂E/∂κ = ∫∂V/∂κ > 0` **everywhere**: minimizing the energy drives **κ → 0**, and κ→0 is
exactly the **collapse** region (E3: κ=0 NaN, κ=0.5 runaway, E→−∞). *Self-tuning by the
action alone selects instability.* This is a Weinberg-flavored no-go: the action has no
interior minimum at the "good" value — it is monotone toward the degenerate one.

**Consequence (the design principle):** stability and stationarity are in **tension** — the
action pulls κ down (toward collapse), stability requires κ ≳ κ_crit. A self-tuning system
cannot optimize one and get the other for free. The only self-consistent resting point is
their **balance: the critical edge κ_crit itself** — marginal stability, action as low as
stability permits. The system must **self-organize to criticality**, not to a minimum.

---

## 3. Two principled mechanisms that land on the critical edge

### (A) Conservation-constrained extremum (action-principled — the Q-ball route)
Free minimization collapses; minimization **at fixed conserved charge** need not. The
project already has the needed conservation law: the **v6 conserved density** (`∂_t ρ +
∇·j = 0`, MEMORY: density is conserved, no potential → a protected charge), plus
topological charge. Extremize the *constrained* action

  δ( E − λ Q ) = 0,   Q = ∫ρ d³x  (or the topological charge), λ = Lagrange multiplier,

jointly over `Φ` **and** `κ`. At fixed Q the energy is bounded below (the charge cannot be
radiated away), so a **nontrivial stable soliton exists in a band of κ** (the Q-ball
mechanism), and `∂E/∂κ|_Q = 0` selects a *finite* κ* — the value at which the charge-Q
soliton is the energy extremum. This is fully action-principled: the conservation law is
what converts the runaway into a self-consistent fixed point. **The self-tuned κ* is an
output of the conserved charge, not an input** — the only number-type (a dynamical
attractor) that v62/v63 leave open for a non-algebraic value.

### (B) Homeostatic stability feedback (dynamical-systems route)
Replace the slow flow's blind action-gradient with one that **senses the stability margin**
and opposes collapse:

  κ̇ = −ε ∂E/∂κ + γ · 𝒮(Φ),   𝒮 = a collapse indicator (e.g. P_max − P*, or the
                                   growth rate of the core density / curvature).

The action term drives κ down; the feedback `𝒮` pushes κ **up** the moment the core starts
to run away. The balance is a stable fixed point **at the onset of collapse** = κ_crit.
This is **self-organized criticality** (drive toward instability + dissipation at
instability → the system parks at the critical point), the Bak–Tang–Wiesenfeld template.
Less elegant than (A) (the feedback is added, not derived), but robust and directly
runnable as a controller around the existing kernel.

**(A) and (B) are the same fixed point** — the conserved charge in (A) *is* the homeostatic
sensor in (B), because fixing Q is what prevents the density from running away. (A) derives
the feedback from a conservation law; (B) imposes it. The target is identical: **κ → κ_crit
≈ 1** (the E3-measured edge).

---

## 4. Why this is the right shape (precedent, not just analogy)

- **Self-organized criticality** (Bak–Tang–Wiesenfeld): open dissipative systems tune
  *themselves* to the critical point without parameter fine-tuning. Exactly §3(B).
- **Multiple-Point Criticality Principle** (Froggatt–Nielsen): the SM couplings sit where
  multiple vacua are degenerate — a *criticality* selection of constants, predicting the
  top/Higgs masses. The SCP version: couplings self-tune to the stable/collapse boundary.
- **Higgs near-criticality**: the measured Higgs/top put the SM vacuum at the edge of
  (meta)stability — an empirical hint that nature sits at a critical edge, not a minimum.
- **Self-tuning cosmological constant** (scale-invariance / dilaton literature): a field
  that relaxes a coupling to a special value — and Weinberg's no-go for why the naive
  version fails (the same monotone-action obstruction as §2).
- **NN equilibrium propagation / energy-based models**: one energy whose relaxation does
  both inference (fast) and learning (slow). §1 is its physical instance; §3(A) is the
  conservation-constrained version (the physics analog of a normalization/Lagrangian that
  keeps learning from collapsing to a trivial solution).

The unifying statement: **a stable particle is a self-organized critical state of a theory
whose couplings are dynamical and charge-constrained.** The constant κ=50 we use by hand is,
on this view, an approximation to a value the dynamics would *select* as the edge of
stability at the soliton's conserved charge.

---

## 5. Concrete realizations (and the kernel-policy fork)

**X1 — homeostatic controller around the existing kernel (NO kernel edit, runnable now).**
Wrap the kernel in an outer loop that is itself a *dynamical law*, not hand-tuning: after
each short relaxation, update κ_{n+1} = κ_n + γ(P_max(κ_n) − P*) (or a collapse-rate
sensor from `diag`). Hypothesis: κ_n converges autonomously to κ_crit ≈ 1 from *both* sides
(from above: stable, slowly drawn down by the action; from below: collapse sensed, pushed
up). *Confirm*: convergence to a κ* near the E3 threshold, independent of the start.
*Refute*: runaway/oscillation with no fixed point. This demonstrates the SOC mechanism
§3(B) with zero kernel changes — the controller is the "intrinsic law" in discrete time.

**X2 — couplings as genuine dynamical fields (INTRINSIC; requires kernel authorization).**
The real §1: add to `scp_sim` a slow field κ(t) (or κ(x,t)) with a kinetic term and the
v6 conservation constraint (§3A), and evolve the *coupled* EOM. The self-tuned soliton is a
true solution of the extended action. **This modifies the kernel EOM and therefore needs
explicit authorization** (CLAUDE.md kernel policy) — or a *separate* differentiable
re-implementation validated against `scp_sim.c`. This is the version that is "integrated
into the physics model" in the full sense the request asks for.

**Recommended path**: X1 first (cheap, no-edit, tests whether the SOC fixed point at κ_crit
is real and basin-robust), then — only if X1 converges — X2 as the authorized intrinsic
mechanism, with the conserved charge §3(A) as the principled stabilizer.

---

## 6. Honest risks

- **The no-go is load-bearing.** Without a conservation law (A) or a sensor (B), there is
  *no* stable self-tuned point — the system collapses. Anyone who "just adds a kinetic term
  for κ" will get runaway. §2 must be respected.
- **Criticality ≠ a sharp value.** SOC parks the system *near* κ_crit with fluctuations,
  not at a fixed real number to many digits. Whether the self-tuned κ* is sharp enough to
  be a *prediction* (vs a fuzzy critical band) is the empirical question X1 answers.
- **κ_crit depends on the conserved charge Q.** The self-tuned value is `κ_crit(Q)` — so
  this predicts a *relation* between the coupling and the soliton's charge, not a universal
  constant. That relation (κ* as a function of Q) is the falsifiable content, and it is the
  intrinsic-dynamics analog of an RG trajectory (CONCEPT §5 / E4).
- This selects a *stable coupling*, an algebraic-adjacent quantity — it does **not** by
  itself produce the transcendentals (α, the Brannen phase); those still need the
  learned-RG-fixed-point route (E4). Self-tuning fixes *stability*, not the number-type
  no-go.
