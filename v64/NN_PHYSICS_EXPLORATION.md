# v64 — Neural-Net-Like Methods for the SCP Soliton Program

**Date**: 2026-06-05
**Task**: (Part 1) digest `arXiv:2606.04032v2`; (Part 2) develop a specific,
technically grounded exploration of treating this project's physics "more like a
neural net," grounded in both the paper and the project's actual tooling/state.
**Parents**: [`PLAN.md`](PLAN.md), [`FINDINGS_H1.md`](FINDINGS_H1.md),
[`../v62/THESIS.md`](../v62/THESIS.md), [`../v63/PREDICTIVITY_AUDIT.md`](../v63/PREDICTIVITY_AUDIT.md).

---

## Part 1 — The paper: *Do Transformers Need Three Projections?*

### 1.1 What it actually is
`arXiv:2606.04032v2` is **not** a physics paper. It is an **ICML 2026 machine-learning
systems paper** — Kayyam, Madan Gopal & Lewis (BrainChip Inc.),
*"Do Transformers Need Three Projections? Systematic Study of QKV Variants."* Domain:
efficient transformer architectures / on-device LLM inference. The source is a standard
`icml2026.tex` two-column submission with vision/NLP benchmark tables, not equations of
motion.

### 1.2 Core thesis
Standard self-attention learns **three** separate linear projections per token —
query `Q=XW_q`, key `K=XW_k`, value `V=XW_v` (Sec. 2.1, Eq. for `A_h =
Softmax(αQ_hK_h^T)V_h`). The paper asks whether all three are necessary and evaluates
three **weight-tying** constraints (Sec. 3.1):

- **Q=K−V** (tie query=key, separate value): `A = Softmax(αKK^T)V` — *symmetric* attention.
- **Q−K=V** (separate query, tie key=value): `A = Softmax(αQK^T)K` — *asymmetric* attention.
- **Q=K=V** (single projection): `A = Softmax(αKK^T)K`.

### 1.3 Key results (precise)
- **Q−K=V is the winner.** On a 300M-param GPT (10B SlimPajama tokens) it costs only
  **+3.1% perplexity** vs. the QKV baseline (val PPL 5.27 vs 5.11) while giving **50% KV
  cache reduction** — because with `V=K` only `K` need be cached (Tab. 4, Tab. 7,
  Sec. 4.3.3). Q=K−V *trains* nearly as well (+4.9%) but gives **zero** cache benefit
  (must still cache both K and V), so it is useless for deployment. Q=K=V collapses
  catastrophically (**+25.4%**).
- **Orthogonality / compounding.** Projection sharing (a *projection-matrix* constraint)
  is orthogonal to head sharing GQA/MQA (a *number-of-KV-heads* constraint). Combined:
  **Q-GQA-4 → 87.5%**, **Q-MQA → 96.9%** cache reduction (Sec. 3.2, abstract).
- **Why Q−K=V works.** Empirical mechanistic claim (Sec. 5 / App.): in trained QKV
  models **K and V are nearly redundant** — cosine similarity 0.73, effective rank
  687 vs 702 of 1024 — whereas **Q is distinct** (cos sim 0.42 with K, 0.31 with V).
  Attention operates in a **low-rank regime**, so tying K=V loses little, but tying
  Q=K destroys the *directionality* (`QK^T` asymmetry) that causal LM needs.
- **Symmetry diagnosis.** Q=K forces a **symmetric** attention map `KK^T` (mirror about
  `y=x`); fine for non-causal vision/set tasks, fatal for causal sequence tasks. They
  patch it with fixed 2D sinusoidal positional encodings ("(X)$^+$" variants) to
  re-inject asymmetry (Sec. 3.1, App. B).

### 1.4 The one genuinely structural result (relevant to Part 2)
**Appendix A (`app:qkv_ssm`)** is the only part with conceptual reach beyond systems
engineering. Under full collapse `q_t=k_t=v_t=z_t`, *linear/kernelized* attention
becomes a pure recurrence:
```
S_t = S_{t-1} + φ(z_t) z_t^T        (running state = sum of self outer-products)
y_t = φ(z_t)^T S_t / (φ(z_t)^T Σφ(z_i))
```
The authors note this is a **state-space model (SSM)** `h_t = A h_{t-1} + B x_t,
y_t = C h_t` with the special twist that the **readout `C` is input-conditioned**
(`y_t = φ(z_t)^T S_t`) rather than a fixed matrix — i.e. linear attention = an SSM with
*adaptive observation*, and equivalently a **fast-weight / Hebbian associative-memory**
update. *"The QKV collapse highlights a continuum between programmable memory
(attention) and dynamical systems (SSMs) … representational structure, not scale alone,
determines qualitative behavior"* (App. A, final paragraph).

### 1.5 Connections to the requested themes
- **Machine learning / optimization**: the whole paper (its native domain).
- **Effective field theory / renormalization / emergent gravity / solitons**:
  **none, explicitly.** Zero overlap in content.
- **Geometry/algebra**: only thin — symmetric vs asymmetric attention maps,
  low-rank/effective-rank structure, weight tying as a hard equality constraint.
- The **load-bearing transferable idea** is in §1.4: a learned attention layer is
  *secretly a dynamical system*, and **hard parameter-tying constraints (`V=K`) often
  cost almost nothing** because the unconstrained optimum already lives on (or near)
  the tied subspace — the trained model is redundant relative to its parameterization.

> **Bottom line for Part 1.** The paper's relevance to SCP is *methodological by
> analogy only*: (a) hard weight-tying as a near-free inductive bias when the solution
> is intrinsically low-rank/redundant; (b) the attention⇄SSM⇄fast-weight equivalence,
> i.e. *trained associative memory is a dynamical system whose attractors are its
> learned content*. Both motivate the Part-2 mapping. There is **no** direct physics
> content to import.

---

## Part 2 — Fundamental physics "more like a neural net"

The framing question: can methods built to *train and analyze neural nets* be applied
to *building this physics model*? The SCP program has two halves with opposite verdicts:

- **The constant-derivation half (v59–v63) is CLOSED as a fit** (`v63` honest
  predictivity ratio **0.36×**; `v62` no-go: `α`, `cos(2/3)` are *provably*
  transcendental ⟹ unreachable by any algebraic construction). Both documents name the
  **only live routes: dynamical/RG, or back to the soliton/simulation program.**
- **The soliton/simulation half is alive** — the 6-field Cosserat kernel
  (`sfa/sim/scp_sim.c/.cu`), gradient-flow soliton solvers (`v2/.../main.c`,
  `hopfion_composition/`), and the SFA pipeline.

This is exactly where the NN analogy can be made *operational rather than metaphorical*,
because the project **already does gradient descent on a field energy** — that is
literally training.

### 2.1 The mapping (which correspondences are real, which are superficial)

| NN concept | SCP analog | Verdict |
|---|---|---|
| **weights** `θ` | (a) field configuration `Φ(x)` on the grid; (b) Lagrangian params `{m², μ, κ, η, δ_i}` | **Real but split.** Two distinct "weight" levels — see below. |
| **loss** `L(θ)` | the action / energy `E[Φ]` (for configs); an **observable-matching error** (for params) | **Real for configs** (`E` is *given*, like a fixed loss). **Chosen for params** (we pick the target spectrum). |
| **training = SGD** | **gradient flow** `∂Φ/∂τ = −δE/δΦ` (already implemented: `v2/.../main.c:56` `gradient_flow_step`) | **Real and literal.** The project already runs steepest descent on a field energy. |
| **learned feature / attractor** | a **soliton** (proton/oscillon/braid) = local minimum the flow converges to | **Real.** Solitons *are* the attractors of the training dynamics. |
| **architecture / inductive bias** | the **octonionic G₂⊂Spin(7)⊂Spin(8) structure** + the 6-field Cosserat form (3φ⊕3θ, `η·curl` coupling) | **Real.** This is hard equivariance/weight-sharing — see §2.2(d). |
| **training dynamics across scales** | **RG flow** of `{m², μ, κ, η}`; v62 P3 "Shulga kernel" = a coarse-graining map | **Real** — and *the* live route for the transcendentals. |
| **overparametrization** | grid DOF ≫ soliton moduli; also **gauge redundancy** | **Mixed.** Some is real overparam (moduli ≪ grid); much is *unphysical gauge* — see §2.3. |
| **loss-landscape geometry** | the energy landscape; Derrick scaling; saddle/min structure (MEMORY: "Skyrmion is a lattice saddle point") | **Real and already observed.** |
| **lottery ticket** | does a sub-grid / low-res seed already contain the winning soliton basin? | **Partly real** (multigrid/seed-transfer), partly metaphor. |
| **neural tangent kernel** | linearized flow operator = the **Hessian / fluctuation operator** of `E` about a soliton (the BLV `P/m`, normal-mode `−(Pg')'+Wg=ω²mg` solvers already compute this) | **Real.** The "NTK" here is literally the project's existing normal-mode / dispersion operator. |
| **feature learning vs lazy** | does the *background field* reorganize (φ-density well, θ-cloud) or just ring linearly? v64 Thread C: linear (κ→0) is **inert** | **Real and central** — see §2.2(c). The "magic/nonlinear" sector = "feature learning"; linear = "lazy/kernel regime." |
| **double descent** | over-resolving the grid past the soliton core scale | **Superficial.** No train/test split; no evidence of a double-descent curve. Drop it. |

**The crucial split: two levels of "weights."**
1. **Inner / fast weights = the field `Φ(x)`.** Loss = action `E[Φ]`; training = gradient
   flow; the optimum = a soliton. This is *given* dynamics — the genuine, literal NN
   analogy, and it already runs.
2. **Outer / slow weights = the Lagrangian parameters `{m², μ, κ, η, δ}`.** Loss = an
   *observable-matching error* (does the relaxed soliton have the target mass / charge /
   spectrum?). This is a **bilevel / meta-learning** problem: inner loop relaxes a field,
   outer loop tunes the theory. This is where autodiff buys the most (§2.2).

Note the appealing resonance with Part 1 §1.4: the soliton is a **fast-weight associative
memory** — the field stores, in its own configuration, the "content" (quantum numbers)
selected by the architecture (the algebra/Cosserat form). The attractor *is* the learned
feature.

### 2.2 What NN methods could actually do for THIS project (concrete, falsifiable)

**(a) Differentiable simulation + autodiff through the 6-field kernel → optimize seeds
and parameters.**
The kernel is an explicit time-stepper (`compute_forces` at `scp_sim.c:348`, leapfrog
update); the v2 solvers are explicit gradient-flow loops. Explicit integrators are
**differentiable by construction** (every step is `+,×,∇` stencils). Two routes:
- *Forward-mode / finite-difference adjoint* on the **outer** params (only ~5–6 numbers:
  `m², μ, κ, η, δ_i`) — cheap, no kernel rewrite, just `(E(θ+ε)−E(θ−ε))/2ε` around a
  relaxed soliton. Immediately usable.
- *Full reverse-mode (adjoint) through the relaxation* — needs the kernel re-expressed
  in an autodiff framework (JAX/PyTorch) or a hand-written adjoint. **High cost**, and
  the kernel-modification policy forbids touching `scp_sim.c/.cu`; this would be a
  *separate* differentiable re-implementation, validated against the C kernel.
  *Falsifiable target*: optimize a Gaussian-bump seed so the relaxed object hits a
  target charge/mass — confirm it converges to a known soliton basin faster than
  hand-built seeds. Refute if autodiff seeds give no speedup over `gen_*` generators.

**(b) The action as a differentiable loss to *discover* stable configurations instead of
hand-building seeds.**
The project hand-builds seeds (`gen_braid.c`, `gen_proton_analytical.c`, …). Instead:
parameterize a seed by a small latent vector `z` → decoder → initial field; relax;
backprop a loss = (final energy) + (penalty for wrong topological charge) + (penalty for
fragmentation, which the kernel already detects). This is **"learn the seed generator"**
and directly attacks the v63-named live criterion: *does a stable particle with the right
quantum numbers emerge?* — a search problem, not a numerology problem.
*Falsifiable*: find a stable bound configuration the hand-built generators miss, or prove
(by exhaustive latent search) none exists in a given topological sector.

**(c) The "magic/nonlinear = feature learning" identification is directly testable
(this is v64 Thread C, recast).**
NN theory distinguishes the **lazy/kernel (NTK) regime** (features fixed, only readout
moves) from the **feature-learning regime** (representation reorganizes). v64's CONCEPT
already posits: a *free/linear* truncation (`κ→0`, `η→0`) is **flat and inert**; only the
nonlinear `κ` saturation and `η·curl θ` coupling let the density back-react and curve
things. **That is exactly the lazy-vs-feature-learning dichotomy.** Thread C's planned
control (zero drift at κ→0; drift scaling with κ, η) *is* a feature-learning-onset
measurement. The NN framing sharpens the prediction: there should be an **order parameter**
(e.g. background-field reorganization amplitude) that turns on with `κ` like feature
learning turns on with width/learning-rate scale. *Refute* if the response is linear in
κ from κ=0 (no threshold) — then "magic" is not a regime, just a coefficient.

**(d) Symmetry-constrained architecture = the octonionic structure as hard equivariance
(geometric deep learning).**
The Part-1 lesson — *hard weight-tying (`V=K`) is near-free when the solution is
intrinsically low-rank* — maps onto: **impose the G₂/Spin(7) algebraic constraints as a
hard equivariance/weight-tying on the field ansatz**, rather than hoping the dynamics
discover them. This is the geometric-deep-learning move (equivariant architectures). The
honest caveat from v62: the algebra can only ever output *algebraic* numbers, so this
**cannot** produce `α` or the Brannen phase — but it *can* enforce the algebraic skeleton
(Koide `Q=2/3`, the selection rule) as an exact constraint on configurations, reducing the
search space. Use it as **inductive bias on the soliton search (b)**, not as a constant
generator.

**(e) Learning a coarse-grained kernel = a learned RG step (the live route for the
transcendentals).**
v62 P3 explicitly identifies the "Shulga kernel" as a coarse-graining/integrating-out map
— *literally an RG step* — and v62/v63 say **the transcendentals can only be fixed-point
outputs of such a dynamical map**, never algebraic targets. NN analog: **train a neural
coarse-graining operator** `R: Φ_fine → Φ_coarse` (a learned RG transform; cf. RG-flow /
neural-RG literature) on field snapshots from the kernel, then study its **fixed points**.
If a transcendental (a loop-ratio like the phase `−c₃/(4c₆)`, or an effective coupling
flowing toward `α`) appears as a **fixed-point eigenvalue of the learned RG map**, that is
the *only* model-internal way such a number can be an output rather than an input.
*Falsifiable & high-value*: a measured fixed-point ratio that matches a target
transcendental within tight error, *robust to the RG-map parameterization*, would convert
a v63 "input" into an "output." Failure (parameterization-dependent fixed points) refutes
the dynamical-origin hope cleanly — itself a publishable null, in the project's ethos.

**(f) The Hessian/NTK already exists — use it for landscape analysis.**
The normal-mode solvers (`−(Pg')'+Wg=(ω²/c²)mg`) and BLV `P/m` machinery *are* the
linearized-flow (NTK-like) operator about a soliton. NN loss-landscape tools (spectral
density of the Hessian, mode connectivity between solitons, flatness/sharpness of minima)
transfer directly. E.g.: is the proton basin **flat** (many near-degenerate moduli =
the v62 "flat-direction law" reappearing dynamically) or **sharp**? Mode connectivity
between two soliton minima = a *physical* reaction path / binding channel.

### 2.3 The honest critique (where the analogy misleads)

1. **The loss is given, not chosen — at the inner level.** In NN training you *design*
   the loss; here the action `E[Φ]` is fixed by the physics. So "improve the loss"
   (architecture search on the objective) is **illegitimate at the field level** — it
   would be changing the physics. The freedom is only at the *outer* (parameter) level,
   and even there the targets are experimental data, not design choices. This is the
   sharpest disanalogy and it caps how much of NN practice transfers.

2. **Exact symmetries and conservation laws, not learned approximate ones.** A NN learns
   *approximate, soft* invariances; physics has *exact* gauge symmetry, Noether currents,
   topological charge. Treating charge conservation as a soft penalty (as in §2.2(b)) is a
   numerical convenience that **must not leak into the ontology** — the charge is exactly
   conserved, not "mostly learned." Equivariant *architecture* (§2.2(d)) respects this;
   *penalty-based* enforcement does not and can give unphysical near-misses.

3. **Overparametrization is largely unphysical gauge redundancy.** In NNs,
   overparametrization is a *feature* (smooths the landscape, enables lottery
   tickets/double descent). Here, much of the "excess" grid DOF is **gauge** (e.g. the
   degenerate `(J,P)` sector is non-dynamical — MEMORY: "forced to zero algebraically")
   or zero-mode moduli. So lottery-ticket / double-descent intuitions, which rely on
   *generic* overparametrization, **do not transfer** — the redundancy here has structure
   (gauge orbits, moduli), not random slack. Only the *moduli* part is a genuine
   "flat-minimum" analog.

4. **No generalization gap.** NN epistemics centers on train/test generalization. Solitons
   have no held-out set; the relevant question is *existence/stability of an attractor*,
   not generalization. So "double descent," "grokking," and most regularization theory are
   **metaphor, not method**, here. (Keep: optimization, landscape geometry, equivariance,
   differentiable simulation. Drop: generalization-centric tools.)

5. **The Part-1 "tying is free" lesson is conditional.** `V=K` was near-free *because the
   trained K,V were already redundant* (cos 0.73). The octonionic tying is only "free" if
   the soliton solution genuinely lives on the algebraic subspace — which v62/v63 show it
   does for the *algebraic* observables and **does not** for the transcendentals. So the
   inductive bias is free where the physics is already algebraic, and *wrong* (a forced
   near-miss) where it isn't. Use it surgically.

**Genuine methodological transfer (keep)**: differentiable simulation / autodiff on the
outer params; learned seed generators framed as soliton *search*; the Hessian/NTK
landscape tools (already half-built); the learned-RG-fixed-point route for transcendentals;
equivariant inductive bias on the search. **Metaphor only (drop)**: double descent,
grokking, generalization gap, lottery-ticket-as-stated, "improve the loss" at field level.

### 2.4 Concrete next experiments for v65+

Sized to existing tooling (CUDA 6-field kernel, v2 gradient-flow solvers, SFA pipeline),
ordered by cost-to-value.

**E1 — Finite-difference outer-parameter optimization ("differentiable theory," cheap).**
- *Hypothesis*: the relaxed-soliton observables (mass, charge radius, θ-cloud size) are
  smooth in `{μ, κ, η}`, so a 5-D finite-difference gradient can tune the theory toward a
  target without hand-scanning.
- *Method*: wrap the existing relax-and-measure pipeline (v2 gradient flow or a short
  kernel run) in a finite-difference loop over `{μ, κ, η, m², δ}`; loss = Σ(observable −
  target)². No kernel edit (pure config-level perturbation, respecting the kernel policy).
- *Confirm/refute*: converges to a target spectrum in ≪ the cost of a grid scan ⇒ confirm;
  non-smooth/noisy gradients (fragmentation discontinuities) ⇒ refute, document the
  landscape roughness.
- *Cost*: **low.** ~10–20 short CPU/GPU relaxations per gradient step; days, no new kernel.

**E2 — Learned seed generator for soliton discovery (medium).**
- *Hypothesis*: a latent-parameterized seed + relaxation, with topological-charge and
  fragmentation penalties, finds stable configurations the hand-built `gen_*` generators
  miss (or proves a sector empty).
- *Method*: separate **differentiable JAX/PyTorch re-implementation** of the 6-field
  gradient flow (validated bit-for-bit against `scp_sim.c` energy on fixed inputs — does
  **not** modify the C kernel); backprop seed-latent → relaxed energy. Use the v2
  `verify3d`/`gradient_flow_step` logic as the spec.
- *Confirm/refute*: a new stable bound state with correct quantum numbers ⇒ strong
  confirm; exhaustive latent search finds only known solitons ⇒ useful null (basin census).
- *Cost*: **medium-high.** Requires the differentiable re-implementation + validation;
  weeks. The validation against the C kernel is the gate.

**E3 — Feature-learning order parameter for the "magic" sector (cheap, high-insight;
this is v64 Thread C upgraded).**
- *Hypothesis*: background-field reorganization (the lapse well + θ-cloud amplitude) turns
  on with a **threshold** in `κ` (and `η`), like NN feature-learning onset — not linearly
  from zero.
- *Method*: the planned Thread-C runs (drift/well amplitude vs `κ`, `η`, with κ→0 and η→0
  controls), but **measured as an order parameter** (background reorganization, not just
  test-particle drift), scanning fine near κ=0 to resolve threshold vs linear onset.
- *Confirm/refute*: clear onset/threshold ⇒ "magic = feature-learning regime" is a real
  structural statement; smooth-linear-from-zero ⇒ "magic" is just a coefficient (refutes
  the regime framing), a clean null for CONCEPT.md.
- *Cost*: **low.** Already-planned runs; only the analysis/order-parameter is new.

**E4 — Learned RG fixed point as a transcendental's home (high-risk, highest-value).**
- *Hypothesis*: a trained neural coarse-graining map `R` on kernel snapshots has a
  **fixed point whose eigenvalue/loop-ratio matches a target transcendental** (the phase
  `−c₃/(4c₆)`, or an effective coupling flowing toward `α`), *robustly to `R`'s
  parameterization* — the only model-internal route v62/v63 leave open.
- *Method*: generate fine/coarse field-snapshot pairs from the kernel at two resolutions;
  train `R` (small CNN, equivariant if possible) to reproduce the coarse dynamics;
  analyze its fixed points and linearization eigenvalues; test robustness across `R`
  architectures.
- *Confirm/refute*: a parameterization-robust fixed-point ratio matching a target to tight
  error ⇒ converts a v63 "input" into a dynamical "output" (the program's biggest possible
  win). Parameterization-dependent fixed points ⇒ refutes the dynamical-origin hope, a
  decisive null that closes the last live route on the constant side.
- *Cost*: **high.** New ML pipeline + snapshot generation; weeks–months. But it is the
  *only* experiment that could reopen the v59–v63 constant program on legitimate grounds.

**Recommended order**: E3 (already-planned, cheap, sharpens v64's central claim) → E1
(cheap differentiable-theory proof of concept) → E2 (soliton search, the v63 live
criterion) → E4 (the moonshot, only if E1/E2 build the differentiable infrastructure).

---

## Summary verdict

Part 1's paper is an ML-systems result with **no physics content**; its one transferable
idea is that *hard weight-tying is near-free when the solution is intrinsically low-rank/
redundant*, plus the attention⇄SSM⇄fast-weight equivalence (trained associative memory =
a dynamical system whose attractors are its content). Part 2 finds the NN↔physics mapping
is **genuine precisely where the project already does gradient descent on a field energy**:
field=fast-weights, action=loss, gradient-flow=training, **solitons=learned attractors**,
the octonionic algebra=equivariant inductive bias, RG=cross-scale training, and the
existing normal-mode/BLV operator=the NTK/Hessian. It is **metaphor** for
generalization-centric tools (double descent, grokking, lottery ticket) because physics
has no held-out set, exact (not learned) symmetries, and gauge (not generic)
overparametrization. The highest-value transfers are differentiable-simulation parameter
tuning, learned seed generators reframing the v63-live "does a stable particle emerge"
criterion as search, the "magic=feature-learning" order-parameter test (v64 Thread C
upgraded), and — the moonshot — a **learned RG fixed point** as the only legitimate home
for the transcendental constants v62/v63 proved unreachable by algebra.
