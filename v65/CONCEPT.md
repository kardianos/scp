# Physics as a Bilevel Learning System — a Quantitative Formulation

**Date**: 2026-06-05
**Status**: founding concept for v65. Companion to the conceptual map in
[`../v64/NN_PHYSICS_EXPLORATION.md`](../v64/NN_PHYSICS_EXPLORATION.md); this document
supplies the *quantitative* realization — the explicit losses, gradients, order
parameters, and the one place the analogy can do real work that algebra provably cannot.
**Parents**: [`../v63/PREDICTIVITY_AUDIT.md`](../v63/PREDICTIVITY_AUDIT.md) (constant
program closed as a fit), [`../v62/THESIS.md`](../v62/THESIS.md) (transcendentals
unreachable by any algebraic map), [`../v64/PLAN.md`](../v64/PLAN.md).

---

## 1. The thesis in one line

The SCP program is **already a learning system** — it runs steepest descent on a field
energy and calls the minima "particles." Writing it explicitly in the language of
neural-net training is not a metaphor at the level where it bites: it exposes a **bilevel
optimization** whose inner loop is the existing soliton relaxation and whose outer loop is
theory selection, and it identifies the **one legitimate route** (a learned
renormalization-group fixed point) by which the transcendental constants v62/v63 proved
*unreachable by algebra* could become *outputs* rather than inputs.

The analogy is exact where the project does gradient descent on a given energy, and
**fails precisely where NN epistemics rely on a chosen loss, a held-out set, or generic
overparametrization** — physics has none of these. The value is in the former; the
discipline is in not importing the latter (§6).

---

## 2. Two levels of weights (the quantitative skeleton)

Let `Φ(x) ∈ ℝ⁶` be the 6-field configuration (3φ ⊕ 3θ) on the grid, and
`Θ = {m², μ, κ, η, δ_i}` the Lagrangian parameters. The energy functional is

  E_Θ[Φ] = ∫ [ ½|∇φ|² + ½ m²|φ|² + V(P) + ½|∇θ|² + ½ m_θ²|θ|² − η θ·(∇×φ) ] d³x,
  V(P) = (μ/2) P² / (1 + κ P²),  P = φ₀φ₁φ₂.

**Inner weights = the field Φ.** The loss is `E_Θ[Φ]` itself (a *given* objective, not a
chosen one), and training is the gradient flow the project already runs,

  ∂Φ/∂τ = − δE_Θ/δΦ          (steepest descent; `v2/.../gradient_flow_step`, the kernel's relaxation)

whose attractors `Φ*` are solitons — protons, oscillons, braids. **A soliton is a learned
feature: the field stores, in its own configuration, the quantum numbers the architecture
admits.** (This is the SCP image of Part-1's attention⇄fast-weight equivalence: a trained
associative memory is a dynamical system whose attractors are its content.)

**Outer weights = the parameters Θ.** Here the loss is *chosen* — an observable-matching
error against data:

  L(Θ) = Σ_k w_k ( O_k[Φ*(Θ)] − O_k^target )²,   Φ*(Θ) = argmin_Φ E_Θ[Φ],

with observables `O_k` = mass, charge radius, binding energy, θ-cloud size, spectral gaps.
This is a **bilevel / meta-learning** problem: inner loop relaxes a field, outer loop tunes
the theory. *This is the layer where autodiff buys the most, because Θ is only 5–6 numbers.*

---

## 3. The hypergradient reuses machinery the project already has

The outer gradient `dL/dΘ` needs `dΦ*/dΘ`. Differentiate the inner stationarity condition
`δE_Θ/δΦ |_{Φ*} = 0`:

  H(Φ*) · (dΦ*/dΘ) + ∂²E/∂Φ∂Θ = 0   ⟹   dΦ*/dΘ = − H⁻¹ ∂²E/∂Φ∂Θ,

where **H = δ²E/δΦ²** is the fluctuation/Hessian operator about the soliton. Then

  dL/dΘ = Σ_k 2 w_k (O_k − O_k^t) [ ∂O_k/∂Φ · (−H⁻¹ ∂²E/∂Φ∂Θ) + ∂O_k/∂Θ ].

The point: **H is exactly the operator the project already computes.** The normal-mode
solver `−(P g')' + W g = (ω²/c²) m g` and the BLV `P/m` fluctuation analysis *are* the
linearized-flow operator — the SCP "neural tangent kernel." So the hypergradient is not new
infrastructure; it is the existing normal-mode operator applied to one extra right-hand
side (`∂²E/∂Φ∂Θ`, an analytic derivative of the Lagrangian). The flat directions H leaves
(the v62 "flat-direction law": Z₃-symmetric invariants are φ-independent) are the soliton's
**moduli** — the genuine flat-minimum analog, and the only legitimate "overparametrization"
here (§6).

**Cheapest realization (no adjoint, no kernel edit):** since Θ is 5–6D, skip H⁻¹ entirely
and take a finite-difference hypergradient,

  ∂L/∂Θ_a ≈ [ L(Θ + ε e_a) − L(Θ − ε e_a) ] / 2ε,

each evaluation = one relax-and-measure. 12 relaxations per outer step × ~10 steps ≈ a few
GPU-hours — the v65 **E1** experiment, runnable now with zero changes to `scp_sim.c`.

---

## 4. "Magic = feature learning," made quantitative

NN theory separates the **lazy/NTK regime** (features fixed, only the readout moves) from
the **feature-learning regime** (the representation reorganizes). v64's CONCEPT already
posits the physics version: a free/linear truncation (κ→0, η→0) is *inert*; only the
nonlinear κ-saturation and η·curl coupling let the background **reorganize** (the lapse
well, the θ-cloud). That *is* the lazy-vs-feature dichotomy, and it has a sharp order
parameter and a computed threshold.

Define the relative nonlinearity of the back-reaction force (the kernel carries
`V'(P) = μP/(1+κP²)²`):

  χ(κ,P) = 1 − V'(P;κ)/V'(P;0) = 1 − 1/(1 + κP²)².

- **Lazy** (κP² ≪ 1): χ ≈ 2κP² — *linear in κ*; the response is a fixed-kernel coefficient.
- **Feature** (κP² ≫ 1): χ → 1 — saturated; the nonlinearity dominates and the background
  reorganizes.
- **Crossover** at κP² ~ O(1), i.e. **κ_c ≈ 0.41/P²** (χ=0.5; `feature_onset.py`).

At standard κ=50, P~1 this gives κP² = 50 ⟹ χ = 0.9996: **the SCP vacuum sits deep in the
feature-learning regime** (force ~100% quenched), far past the κ_c ≈ 0.41 knee. This is the
quantitative content of "magic is on": the theory is not a perturbed free field, it is in
the reorganizing regime, by a factor ~120 in κP².

**Falsifiable E3 prediction (v64 Thread C, upgraded to an order parameter):** scan κ from 0
and measure a *background-reorganization* order parameter `M(κ)` (e.g. the lapse-well depth
‖Φ_bg(κ) − Φ_free‖, or the equilibrium θ-cloud amplitude) in a real relaxation. It must show
a **knee near κ ~ 1/P² ~ O(1)** and saturate — not rise linearly from zero. *A straight line
from κ=0 refutes "magic = a regime"* (it would make magic merely a coefficient) — a clean
null in the project's ethos.

---

## 5. The learned RG fixed point — the only legitimate home for the transcendentals

v62's no-go is exact: an algebraic construction (characters, Casimirs, dimensions,
eigenvalues of rational matrices) outputs only **algebraic numbers**, so `α` and the
Brannen phase `cos(2/3)` — *transcendental* — are unreachable by any such map. v63 names the
sole surviving route: a **dynamical/RG fixed point**, because the RG `β`-function is built
from **loop integrals (logarithms)** and its fixed points are generically transcendental.

Quantitatively: let `R_w: Φ_fine ↦ Φ_coarse` be a learned coarse-graining map (block-spin +
a small equivariant correction network, weights `w`), trained so that coarse-grained
kernel snapshots evolve under `R_w` as the kernel does at the coarse scale. Couplings flow,

  g_{n+1} = β_w(g_n),   fixed point g* = β_w(g*),   exponents = eig(∂β_w/∂g |_{g*}).

A transcendental constant can appear *only* as such a `g*` or eigenvalue — never as the
output of the algebraic skeleton. Concretely, the **v62 phase loop-ratio `−c₃/(4c₆)`** or an
**effective coupling flowing toward `α`** is the target.

- **Win condition (the program's biggest possible):** a fixed-point ratio matching a target
  transcendental to tight error, **robust to the parameterization of `R_w`** — this converts
  a v63 *input* into a dynamical *output*.
- **Decisive null:** parameterization-dependent fixed points ⟹ the dynamical-origin hope is
  refuted, closing the last live route on the constant side. Either way it is publishable.

This is **E4** — high cost (a differentiable re-implementation + a neural-RG pipeline), but
it is the *only* experiment that could legitimately reopen the v59–v63 constant program,
because it is the only one whose number-type (transcendental fixed point) can clear the v62
no-go.

---

## 6. The honest disanalogies (quantitative discipline)

What transfers, and what is metaphor that must be dropped:

| NN notion | transfers? | why / why not |
|---|---|---|
| gradient flow = training | **yes, literal** | the project already runs `∂Φ/∂τ = −δE/δΦ` |
| soliton = learned attractor | **yes** | minima of the given energy |
| Hessian = NTK | **yes** | = the existing normal-mode / BLV operator (§3) |
| differentiable parameter tuning | **yes** | Θ is 5–6D; finite-difference hypergradient (§3) |
| equivariant architecture (octonionic G₂/Spin(7)) | **yes, surgically** | hard constraint on the *algebraic* observables only — v62 forbids it producing transcendentals |
| learned RG fixed point | **yes, the moonshot** | only number-type that clears the v62 no-go (§5) |
| "improve the loss" at field level | **NO** | the action is *given*; changing it changes the physics |
| double descent / grokking / lottery ticket | **NO** | no train/test split; epistemics centered on attractor existence, not generalization |
| generic overparametrization | **NO** | the grid "excess" is mostly **gauge** (the degenerate (J,P) sector is non-dynamical) or moduli — structured, not random slack |
| charge as a soft penalty | **only numerically** | topological charge is *exactly* conserved; a penalty is a search convenience that must not leak into the ontology |

The single sharpest disanalogy: **at the inner (field) level the loss is the action, which
is given.** All NN practice that presumes a designable objective is illegitimate there. The
design freedom lives only at the outer (parameter) level, and even there the targets are
experimental data. This caps the transfer — and is exactly why E1/E3/E4 are framed as
*tuning given dynamics* and *measuring regimes*, never as "search the objective."

---

## 7. Minimal first realization (E1 ⊕ E3, runnable now)

A single combined experiment that needs **no kernel modification** (config-level only):

1. Fix a relaxed soliton (template seed, short kernel relaxation).
2. **E1 axis** — finite-difference hypergradient on `Θ = {μ, κ, η}` toward a target
   observable (mass / charge radius): `∂L/∂Θ_a ≈ (L(Θ+εe_a) − L(Θ−εe_a))/2ε`. Confirm L
   decreases in ≪ the cost of a grid scan; refute if fragmentation makes the gradient
   non-smooth (and then *report the landscape roughness* — itself a result).
3. **E3 axis** — along the κ direction, measure the background-reorganization order
   parameter `M(κ)` and test for the predicted knee at κ_c ≈ 0.41/P² (`feature_onset.py`)
   vs a linear-from-zero null.

Deliverable: either (a) a working "differentiable theory" that tunes the spectrum by
descent + a confirmed feature-learning threshold (magic is a *regime*), or (b) two clean
nulls (rough landscape / linear onset) that bound the analogy. Both advance the program;
neither touches the protected kernel.

**Sequencing** (from `NN_PHYSICS_EXPLORATION.md` §2.4): **E3 → E1 → E2 (learned seed
search) → E4 (learned-RG moonshot)**. E2/E4 require a *separate* differentiable
re-implementation validated bit-for-bit against `scp_sim.c` — not edits to it.
