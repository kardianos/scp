# SCP Future Directions and Open Questions

**Purpose**: Tracks future research directions, unresolved questions, and
exploratory paths. Maintain with CONCEPT.md and EM_THEORY.md.

**Current equation set (2026-06): kernel-v3** — complexified 12-field Cosserat
with gauged diagonal U(1) (`complex_phi=1`, `complex_gauge=1`, g=0.05,
m_θ=1.6). Sections below this first one predate the U(1) era (≤v65 real-field
equations) and are kept for the record; most are superseded — the real theory
has no stable particles (CONCEPT.md §2/§9).

---

## Active goal (v75 multi-fabric → C₁₂) — 2026-07-15

**Doc:** `v75/C12_ATOM_GOAL.md` · **Status:** `v75/STATUS.md`

Pursue multi-fabric **Phases 1–3** on B1 (private bags + shared \(A\)):

1. **P1** — Bound hydrogenoid (retuned multi-rev C+L on B4 θ)  
2. **P2** — Parked multi-Z nucleus + light cloud, long-T stable  
3. **P3** — A≈12-class core + multi-L → **time-stable visual C₁₂**  
   (no merge / fission / disperse)

**Stretch:** radioactive C₁₂ variant with **calculated and simulated** decay rates.

Baseline closed: F11–F16 (isolation, packaging, pair Coulomb scan, Z6+L6 PASS_park).

---

## Current Open Questions (post-v71, 2026-06-11)

1. **ℓ=2-perturbed cold fission test** — super-critical balls survive the
   symmetric saddle and hot mergers evaporate rather than split; does a cold
   deformed Q > Q_max ball undergo binary fission? (v69 queue #1)
2. **Clean big-box force exponent** — pair-dynamics force laws are
   boundary-contaminated beyond D ≈ L/2 (v70); a L ≳ 2D (N=384) run would give
   the direct n=2.00 exponent if ever needed beyond the halo measurement.
3. **Flavored branch extent** — the BVP continuation reached Δω = 0.04 easily;
   where does the flavored branch end, and does deeper detuning produce
   visually stratified (concentric-shell) baryons?
4. **η ≠ 0 in the complex theory** — the curl coupling ties the component
   (color) index to spatial direction; its effect on baryon substructure,
   flavor, and the GDR mode is unexplored. Includes the in-medium η_eff ~
   g|Φ_bg| re-emergence of the old polariton phenomenology (v68 §1.2).
5. **Quantized-orbit / positronium dynamics** — clean opposite-charge orbit at
   D ≈ 16 with tangential velocity (per-ball velocity seeding: small
   gen_qball_multi extension); phase-closure quantization test (DEBROGLIE).
6. **AB / flux-line winding** — holonomy 2πqn is now topological (gauged);
   needs a flux-tube seeder (v69 O4).
7. **Charged-bath equilibrium** — replenished co-rotating bath for a true
   detailed-balance steady state (v68 #5); star-interior/accretion program.
8. **GDR systematics** — mode frequency vs flavor asymmetry and nucleus size;
   damping mechanism (radiation vs internal mixing).
9. **EM_THEORY.md rewrite** — map the EM sector onto the gauge field A
   (photon analog); θ is now massive charged torsion matter (v69 queue #6).
10. **Gravity** — still absent from the framework; the clock-gradient force
    a = −c²(1−v²)∇(ln ω) (v67 DEBROGLIE §3) is the surviving candidate
    mechanism, untested in the gauged theory.
11. **Gauged interlock molecule** — phases alone cannot hold a static two-ball
    standoff (v71 COLLIDE); Coulomb repulsion + aligned-channel attraction
    might. One run at g=0.05 with the two-low/one-high pair.
12. **Collision systematics** — velocity and impact-parameter sweeps of the
    flavored taxonomy; the merge/bounce boundary; fc3's charge-transfer
    direction vs phase convention.
13. **Per-cluster per-flavor charge matrix** — small sfa_qcomp extension
    (cluster-masked integrals) to measure flavor exchange exactly.

---

## The Classical→Quantum Bridge — a QFI Entanglement Witness (2026-06-26)

**Goal (CONCEPT.md §9).** Show that the *classical* SCP field reproduces a
genuinely *quantum* many-body signature — a multipartite-entanglement witness —
computed from the dynamical structure factor S(q,ω) the sim already produces,
using ℏ_eff = Q (the U(1)-soliton action quantum, §7) as the quantum of action.
Target: the quantum Fisher information of Mazza et al. (*Quantum Fisher
information in a strange metal*, Nat. Phys. 2026, doi:10.1038/s41567-026-03298-0).
This is the operational test of the project's central claim — that a particle is
not an object on spacetime but a collective configuration of the fabric.

**The pipeline (read-only analysis; no kernel change).**
1. Local observable O(x,t): gauge-invariant charge density ρ_Q(x,t) (whole-ball
   witness) and per-component |Φ_a|² (flavor/quark witness). Both already derive
   from existing 30-column SFA frames.
2. S(q,ω) = |FT_{x,t} δO|² at a chosen wavevector q — pick a q away from any
   "ordering" feature (the analog of Mazza's q = (0 1̄ 0)). Needs a frame cadence
   in t fine enough to Nyquist-resolve a few × the clock ω.
3. χ″ from S by FDT. Run both kernels and compare:
   - **bridge measurement** (quantum kernel, ℏ_eff = Q):
     f_Q = 4 ∫ tanh(Qω/2T_eff)(1 − e^{−Qω/T_eff}) S(q,ω) dω;
   - **control** (classical FDT χ″ = (ω/2T_eff) S): should *not* witness
     entanglement — the quantum content lives entirely in the tanh kernel.
4. nQFI = f_Q/(h_max − h_min)². The operator range (h_max − h_min) for a
   continuous classical density is an open modeling choice — candidates: per-cell
   charge range, or the Q/3 quantization unit (§4). Also report raw f_Q growth
   and its scaling exponent, which are normalization-independent.

**Open knobs / honesty.**
- T_eff: the sim is not manifestly thermal. Extract an effective temperature from
  the small-amplitude fluctuation spectrum (equipartition), or scan T_eff and
  report f_Q(T_eff).
- Success ≠ proof the field is quantum. It tests whether ℏ_eff = Q is strong
  enough to carry entanglement, not just kinematics. A null (nQFI < 1 everywhere)
  bounds the reach of the soliton-ℏ identity and is equally publishable.

**Concrete experiments (guesses, priority order).**
1. **Does binding raise the witness?** f_Q of an isolated stable Q-ball vs. a
   single dispersing component lump (s ≡ 0, no binding — §4). Prediction: the
   color-complete bound ball registers higher nQFI than free waves. Cheapest
   first test; reuses gen_qball_boost (v=0) and gen_qball_quark seeds.
2. **Scale-free criticality at the Affleck–Dine threshold.** Sweep condensate
   amplitude through κA⁶ = ½ (A* = 0.464, §3) and attempt Mazza's data collapse
   S(q,ω,T_eff)·T_eff^α = W(ω/T_eff). A *fractional* α at fragmentation would be
   this theory's own beyond-order-parameter critical point — the direct analog of
   the strange-metal QCP, where condensate→particle is the transition.
3. **Does fusion increase multipartite content?** Track f_Q across a two-nucleon
   fusion (§6): does the bound droplet witness more entanglement than the
   separated pair? Reuses gen_qball_pair / gen_qball_multi.
4. **Flavor in the witness.** Symmetric vs flavored baryon (Δω = 0.04, §4): do
   the distinct internal clocks appear as structure in per-component S(q,ω) and a
   different nQFI?

**Tooling.** New analysis binary `sfa/analysis/sfa_qfi` (reads SFA, builds
S(q,ω) for ρ_Q and per-component densities, applies both FDT kernels, reports
f_Q, nQFI, and the scaling collapse). Shares the FT infrastructure with the
existing clock-harmonic and GDR spectral tools.

### Status — tool built and validated (2026-06-26)

`sfa/analysis/sfa_qfi` is built (Makefile ANALYSIS_TOOLS) and validated on a
fine-cadence single static ball (N=64, ω=1.45, Q=118, snap_dt=0.25, 161 frames)
and on the 30-column gauged flavored ball fl14. It recovers the internal clock
(w_peak=1.405 ≈ ω) and is symmetric across the three components, as it must be.
The first run pins down exactly what remains for a *positive* witness — the
"remainder required for realization." Each item is now concrete:

1. **Observable must be gauge-invariant.** The carrier field Φ_a carries a
   coherent clock line at ω whose f_Q≈27 is large but is a *coherent* (classical,
   zero-entanglement) oscillation and is gauge-dependent — it must NOT be counted.
   The physical witness uses ρ_Q (or |Φ_a|²). DONE: the tool reports both; ρ_Q is
   the witness, fieldN is kept only as a clock/validation diagnostic.
2. **A single static ball is a correct NULL.** Its gauge-invariant nQFI ≈ 1.5e-5
   — the f32 snapshot noise floor, i.e. ≈ 0. There is no many-body entanglement
   in one classical soliton, and the witness correctly reports none. A positive
   witness therefore *requires dynamical fluctuation content*: the GDR collective
   mode (flavored nucleus, ω≈0.14), two interacting balls (clock-interference
   beat), or the Affleck–Dine fragmentation fan. These are experiments 1–4 above.
3. **Operating point: T_eff ~ ℏ_eff·ω.** Because ℏ_eff = Q ≈ 118 is large,
   ℏ_eff·ω ≈ 170 ≫ 2T_eff for any modest T, so the tanh·(1−e^…) kernel
   *saturates to 1* and f_Q just integrates S(q,ω) — the quantum/thermal
   discrimination that makes QFI a witness is lost (confirmed: f_Q flat to 5
   digits over T_eff = 0.25→8). To use the kernel as Mazza does, the run must sit
   where ℏ_eff·ω ≈ T_eff: either **small-Q balls near the existence floor**
   (Q ~ O(1)) or **T_eff ~ Q·ω**. Fixing this operating point is the single
   sharpest gap.
4. **T_eff operational definition** — extract from the small-amplitude
   fluctuation-equipartition spectrum (open).
5. **(h_max − h_min) normalization** to convert f_Q → a bound-crossing nQFI
   (modeling choice: per-cell charge range or Q/3 unit; open).
6. **Cadence** — snap_dt ≲ 0.25 resolves the clock (validated); archived
   production SFAs at dt≈10 alias badly, so each witness experiment needs a
   purpose-run fine-cadence output. DONE as a requirement.

**Remainder summary.** Infrastructure (items 1, 6) is complete. Realization is
gated on items 2 + 3: run a *dynamical, gauge-invariant* configuration at the
operating point ℏ_eff·ω ~ T_eff and show ρ_Q-nQFI crosses 1. Items 4, 5 are
modeling choices that scale the bound but do not gate the qualitative crossing.
Cheapest first realization: the flavored-nucleus GDR (experiment 1+4) sampled at
snap_dt≈1 with T_eff set to Q·ω_GDR.

### QFI as a selection principle — goal function + parameter-space limiter

A single matter soliton factorizes, so its gauge-invariant witness is null
(item 2). Multipartite entanglement is not in the matter field alone but in the
**constraint-induced correlations** between sectors, which the current kernel
already contains:

- **gauge constraint** — the discrete Gauss law div E = g·ρ_Q ties the electric
  field E (columns E_x/E_y/E_z) nonlocally to the matter charge; A is th_x/y/z.
- **geometric constraint** — the Cosserat curl coupling η·Re[Θ̄·(∇×Φ)] locks the
  torsion partner Θ (columns theta_*/thetaim_*) to the matter vorticity ∇×Φ.

These couplings make the sectors *not independent*; their cross-coherence
Γ(q) = Σ|S_{ρ_Q,X}|/√(ΣS ΣS) (now reported by `sfa_qfi`) is the correlation an
entanglement witness reads. So the program is to drive these constraints and use
QFI in two roles:

1. **Goal function** J(θ) = nQFI[ρ_Q](θ) (optionally the joint matter+gauge+Θ
   witness), maximized over the theory's otherwise-free parameters
   θ = (g, η, m_θ, ω/Q). CONCEPT.md §8 lists g, m_θ, η as *free inputs*; this
   turns the entanglement witness into the principle that fixes them — the
   classical→quantum bridge becomes a parameter-selection law, not just a check.
2. **Parameter-space limiter** — the feasible set is
   {θ : operating point ℏ_eff·ω ≈ T_eff (item 3, use `--auto-T`);
        Gauss residual at the 1e-13 floor; VK-stable (dQ/dω<0);
        θ-drain closed (m_θ>ω); charge conserved to the integrator floor}.
   QFI>1 prunes this set to the regime that actually realizes a witness. The
   limiter keeps the search on physical, conserved, stable configurations.

**Method (small optimization loop, all config-only + read-only analysis):**
template configs over a θ-grid → short fine-cadence local runs → `sfa_qfi
--auto-T` → read nQFI[ρ_Q] and Γ from the TSV as the objective → coordinate-
ascend / prune by the feasibility limiter. Starts as a coarse (η, g) sweep at
fixed ω, then refine. No kernel change — the constraints are already in the
kernel; only η (run at 0 for particle experiments so far) and g are turned on.

**First test result (2026-06-26).** Same ball (N=64, ω=1.45, Q=118), gauged
g=0.05, m_θ=1.6, snap_dt=0.25, η scanned 0→0.5, witness via `sfa_qfi --auto-T`
(T_eff=ℏ_eff·w_peak/2=12.3, kernel un-saturated):

    eta    E-drift   nQFI[rhoQ]   nQFI[Theta]   Gamma(rhoQ<->Theta)
    0.0    -0.06%    8.3e-8       0  (Theta=0)   0
    0.1    -0.11%    2.8e-5       1.2e-5         0.996
    0.2    -0.30%    4.0e-4       1.7e-4         0.997
    0.3    -0.69%    1.7e-3       7.3e-4         0.996
    0.5    -2.08%    8.2e-3       3.4e-3         0.992

Findings:
- **The geometric constraint lifts the witness off the null.** η:0→0.5 raises the
  gauge-invariant ρ_Q witness by ~5 orders (8e-8→8e-3); the matter density gains
  real dynamical weight only because Θ is sourced and back-reacts. The gauge
  sector alone does not: Γ(ρ_Q↔E)≈0.97 at all η (E is statically Gauss-slaved to
  ρ_Q), but its f_Q stays ~1e-8 — a high coherence on a *static* signal, not
  entanglement. The geometric (Cosserat) coupling is the dynamical one.
- **Constraint-locking is near-perfect and η-independent:** Γ(ρ_Q↔Θ)≈0.996 for
  every η>0. Once the curl coupling is on, Θ is algebraically slaved to ∇×Φ, so
  the sectors are maximally correlated regardless of strength — the cleanest
  available signature of constraint-induced cross-sector structure.
- **Goal-function vs limiter is now quantitative.** nQFI[ρ_Q] ~ η^3–4 (steeply
  increasing objective) while E-drift ~ η^2 (the feasibility cost: η drives the
  θ-drain/radiation). Requiring drift < 0.5% puts the feasible optimum at
  **η ≈ 0.25–0.3**, nQFI ~ 1e-3. This is the selection principle operating: QFI
  picks the coupling, stability limits it.

**Remaining gap to nQFI>1 (genuine realization).** Even at η=0.5 the bare witness
(hrange=1) is ~8e-3 ≪ 1 — the η-induced density modulation is only ~3% (intS~2.6e-3).
Two known levers remain, both already scoped: (5) the (h_max−h_min) normalization
sets the absolute bound-crossing — with the physical per-cell/Q-quantum range this
can be ≫1; and a higher dynamical-weight configuration (large-amplitude collective
mode, or two interacting balls so Γ reflects inter-particle rather than intra-ball
locking). The single static ball with η>0 is now a *non-null* but *sub-threshold*
witness — the method is validated; the realization is reduced to fixing the
normalization and raising the dynamical amplitude.

**Feasible optimum + two-ball (2026-06-26).**
- η-scan refined: under the limiter drift<0.5%, the feasible optimum is
  **η\* ≈ 0.25** (drift −0.465%, nQFI[ρ_Q]=9.1e-4); η=0.3 (−0.69%) is infeasible.
- **Two balls (in-phase, D=8, η=0.3) cross the bound:** nQFI[ρ_Q] = 5.9e-3 (q=0),
  0.124 (q=1), **1.12 (q=2)** — > 1 at the inter-particle wavevector, witnessing
  ≥2-partite entanglement that the single ball never reaches (8e-8). Inter-ball
  coherence Γ(ρ_Q↔φMod)=0.9999 at q=2. The witness reads inter-particle
  correlation (the §5 clock-interference force), as hypothesized.
- Caveat: both still drift (η=0.3 two-ball −1.1%) — dynamical/ringing, not
  stationary. The bound-crossing is real but rides a transient, not a clean
  eigenstate. This is the geometry gap below.

### The geometry requirement — what "satisfies the η-coupled equation" means

The diag shows the η drift is dominated by an **inconsistent initial condition**,
not a runaway: the seed has Θ=0, so Θ is sourced from zero, overshoots
(Q_θ:0→18 by t=10), then rings and exchanges energy with Φ (Q_θ→9→10, E sloshes
85→71→78, bounded). The symmetric ball with Θ=0 is simply **not a stationary
solution** when η≠0. A stationary, 3D-stable η-coupled gauged Q-ball must satisfy
THREE conditions at once:

  (i)   **binding:** s = Π_a|Φ_a|² ≠ 0 throughout the core;
  (ii)  **Cosserat stationarity:** Θ = the consistent partner of ∇×Φ (its
        equilibrium under ∇²Θ − m_θ²Θ + (η/2)(∇×Φ)-source = 0), NOT Θ=0;
  (iii) **Derrick/3D stability** (no collapse, no radiative mode).

The available geometries each miss one:
- **symmetric ball** (Φ_a=f(r)û_a e^{iωt}): satisfies (i), has a stationary
  curl (r̂×û)f′, but is seeded with Θ=0 → violates (ii) → the observed ringing.
- **hedgehog** (Φ_a=f(r)x_a/r·e^{iωt}): clean stationary curl (ii) but
  s ∝ x²y²z²/r⁶ vanishes on the coordinate planes → kills (i).

**Ansatz-first structures rejected (user assessment, 2026-06-26).** Lagrangian-
first geometric candidates — Skyrme–Faddeev/Hopfion, gauged CP², non-abelian
SU(2) color — are judged likely dead ends (Hopfion almost certainly; SU(2) not a
sustainable particle). The hedgehog in particular ALWAYS kills binding (the
product-zero set is codimension-1 in any d), unfixable by any profile. You cannot
patch binding on: the geometry/mathematical structure *is* the binding, so it must
be discovered, not imposed.

### Inverted methodology — fitness-first soliton discovery (ADOPTED 2026-06-26)

Do NOT start from a Lagrangian/ansatz. Start from a **fitness viability on particle
behavior**:

1. represent a particle as a **free-form, gauge-connected field** on the fabric
   that can deform into knots or any shape (blobs with arbitrary per-component
   complex amplitudes; the gauge field threads them via Gauss's law);
2. score each shape by **viability first** — does it retain its shape under
   evolution (alive fraction, core triple-product retention S_mean, low drift) —
   then by the **QFI fitness** we already run;
3. optimize the shape (point by point, in the limit) to maximize fitness, seeking
   **stability first**;
4. only if a viable shape is found, **fit a functional form to it**; the effective
   Lagrangian is the LAST step, inferred from the discovered form — if ever.

The kernel is used only as the behavior oracle (does this shape persist?), so this
needs **no kernel change** — it is seed-generation + analysis + an optimizer loop.

**Machinery built:** `sfa/seed/gen_blob_field.c` (free-form gauge-connected
rotating-blob seed from a text "genome": omega + a list of blobs, each a center,
width, and per-component complex amplitude) and `scratchpad/fitness_search.sh`
(generate genome → kernel-evolve → score viability via `analyze_sfa` + QFI via
`sfa_qfi --auto-T`). QFI is the goal function; viability (S_mean, alive, drift)
is the stability limiter.

**First search (N=48, L=12, T=12, gauged g=0.05, η=0.25; 6 hand-seeded shapes):**

    shape                 alive   S_mean   nQFI[rhoQ]
    c0 symmetric ball     19/25   0.260    4.2e-3
    c1 flavor-asym        23/25   0.208    1.0e-2
    c2 two-lobe split     20/25   0.221    6.1e-2
    c3 tri-blob           21/25   0.242    9.3e-3
    c4 phase-twist        25/25   0.311    4.2e-3
    c5 phase-conn. pair   24/25   0.354    3.1e-1   <-- best on BOTH axes

Result: with no ansatz, the fitness already prefers **phase-structured shapes**
over the plain Q-ball. c5 (two lobes carrying a relative internal phase, tied
through the gauge) has the best shape-retention AND ~73× the ball's QFI (0.306,
approaching the nQFI=1 bound); c4 (a phase-twisted single blob, components at
0/120/240°) survives every frame. The symmetric ball is mediocre on both — the
search is pointing away from it, exactly as the inverted method intends.

Caveats: this is a 6-point proof-of-machinery, not a converged optimization;
T=12/N=48 viability is preliminary (needs longer-T death check); nQFI<1 still (no
bound-crossing yet for a single object). 

### Remote CMA-ES campaign result (2026-06-26)

Ran the full optimizer on a Tesla V100 (Vast.ai via `scp-runner`): separable
CMA-ES (`v71/fitness/cma_search.py`, pure-stdlib, parallel population) over a
K=4-blob genome (40-D), 45 generations × pop 14 (630 evaluations), each a short
gauged run (N=48, L=12, T=16, g=0.05, η=0.25) scored by `eval_genome.sh`
(viability gate × QFI). Artifacts in `v71/fitness/results/` (best.gen, cma_log.tsv,
cma.out).

Engineering notes (for re-runs): per-eval SFAs must be deleted after scoring
(eval_genome.sh does this) — without it the disk fills in ~1 generation. f32 is
required for the fitness; f16 gate-noise stalled the search. Cold-start from random
wanders (all random shapes radiate, drift≫1%, gated to ~0); **warm-start the CMA
mean at the symmetric ball with a small σ0≈0.3** so the search climbs QFI from a
viable point.

**Result — fitness-first discovery succeeded.** From the warm-started ball
(J=0.15, nQFI≈0) the search climbed monotonically to **J=1.02**:

    J:  0.15 (ball) -> 0.27 -> 0.42 -> 0.58 -> 0.68 -> 0.77 -> 1.02
    best shape: alive 16/16 frames, S_mean=1.80 (7x the ball's core retention),
                drift 1.12%, nQFI[rho_Q]=2.02 (past the >=2-partite bound)

The discovered object (regenerated, characterized with `sfa_qcomp`) is NOT the
ball: a **compact 4-lobe cluster** (all blobs within r~2.5, r_rms~3.2) that is
**charge-asymmetric / flavored** — per-component charge partition
**Q ≈ 107 : 285 : 160 (≈1 : 2.7 : 1.5)**, with the three color centroids spatially
displaced (internal flavor texture), dense core (s_max=1.06 vs 0.05), tied
together through the gauge field. Structurally a **flavored baryon/nucleus found
by fitness, not derived from a Lagrangian** — the inverted method delivered a
stable, entanglement-rich shape with no ansatz.

Honest caveats: viability is over T=16 only (drift 1.1% could grow at long T — not
a proven eternal soliton); nQFI=2.0 used `--auto-T` at coarse cadence; it is a
local optimum of this fitness/genome. 

### Scale-up campaign + the η-drain drift floor (2026-06-26) — SUPERSEDED by v72 (see RESOLVED below: the floor was tooling, not physics)

Scaled the search: **persistence fitness** (viability = √(S_mean·S_final)·alive,
so the core must survive to the END of the run, not just on average),
**warm-from-genome** (`--warmgen`, inverse-encode any genome as the CMA mean),
**60-D genome** (K=6 blobs), N=96 / L=18 / **T=30** (larger grid + ~2× time), on
an A100. Two findings, one of them limiting:

1. **The N=48 winner is not stable at scale.** Run at N=96/T=40 it drifts 7.7%
   (vs 1.1% at the short N=48/T=16 it was found at); the persistence fitness
   scores it 0.65 not 1.02 — S_final 0.77 < S_mean 1.83, i.e. it decays. Short-T /
   small-grid fitness over-rewards dense-but-transient shapes.
2. **A drift floor caps the scaled fitness, and it is intrinsic η-drain.** At
   N=96/T=30 even the warm-started ball sits at ~1.5–2% energy drift, so the
   stability gate floors J (the search plateaued at J≈0.08 from gen 5, while
   finding persistent nQFI≈2.0 shapes that are all gated). Tested the obvious
   suspect — consistent **Θ=(∇×Φ)/2 seeding** (new `gen_blob_field … curl` mode,
   finite-difference curl into the torsion columns): it removes the *early*
   transient (t=12 drift 0.05% vs growing) but the *final* drift is unchanged
   (Θ=0: 1.47%, Θ=curl: 1.64% at t=30). **So the floor is ongoing radiation, not
   an initial-condition transient.**

**The exposed tension (the real physics).** Entanglement (QFI) is produced by the
η coupling — but η is also what radiates (the θ-drain) and drives the drift that
the stability gate punishes. Low η → stable but low-QFI; high η → QFI-rich but
drifting. A blob search at fixed η cannot escape this; the genuinely stable,
entanglement-rich object is the **stationary η-soliton** — a shape whose η-curl
structure carries NO net radiation — which needs the consistent geometry / coupled
Φ–Θ(–E) BVP, not a free blob superposition. The scaled campaign makes this
tradeoff quantitative.

### Stationary η-soliton BVP — first component built (2026-06-26)

Built `sfa/seed/radial_eta_soliton.c`: the symmetric ball Φ_a=f(r)û_a (û=(1,1,1)/√3,
binding s=f⁶/27≠0) has curl f′(r̂×û), which sources a transverse torsion partner
Θ=g(r)(r̂×û). The exact stationary g solves the linear radial BVP

    g'' + (2/r)g' − (2/r²)g − (m_θ²−ω²)g = −η f'      (l=1 transverse ⇒ −2/r² term)

solved by Thomas; the tool writes the consistent (Φ,Θ) seed. Key number: the
m_θ=1.6 mass SUPPRESSES the consistent partner to g_max=0.042 vs the naive
½curl value 0.109 (2.6× smaller) — which is why the earlier ½curl seed over-shot
and didn't help.

Drift test (ball, η=0.25, T=30): **Θ=0 → −2.93%, Θ=½curl → −2.66%,
Θ=BVP-exact → −2.27%.** Monotonic improvement confirms the consistent-Θ mechanism
is real and correctly signed — but it removes only ~25% of the drain. The drift is
flat to t≈20 then *accelerates*: the residual is the **Φ back-reaction**
(η·curl Θ has a ∝r̂ component) forcing the ball into an axisymmetric f(r,θ) the
spherical ansatz cannot hold. So the radiation-free soliton is genuinely
**2D axisymmetric**, not 1D radial.

**Conclusion / next:** the transverse-Θ BVP is component one (built, verified
directional). The remaining drain is dominated by the Φ-sector 2D deformation, so
the full radiation-free η-soliton needs a coupled Φ–Θ relaxer.

### 3D coupled relaxer with external-pressure convergence (2026-06-26)

Built `sfa/seed/eta_relaxer.c`: a **3D gradient-flow relaxer** of the effective
static energy of the rotating ansatz (Laplacian + (m²−ω²)/(m_θ²−ω²) effective
masses + saturating Vt(s) + the η-curl coupling, where φ and θ each feel η×curl of
the other). Two numerical ingredients make it converge to the soliton rather than
a trivial state:

1. **Norm preservation** — the η-coupled Q-ball is a *saddle* of the fixed-ω
   energy, so plain gradient flow runs to the vacuum (φ→0). Rescaling to a fixed
   matter norm M₀ each step makes the soliton the constrained minimum.
2. **External pressure (the user's prescription)** — V_ext = P·max(0,r−R)² with a
   TIGHT cage (R≈ball radius) and LARGE P confines the field so a definite compact
   configuration forms; P is then annealed to 0. Without it, at fixed norm the
   field still spreads thin (gradient term wins). With extreme P (P₀=20, R=3) the
   field forms a compact BOUND state (s_max=0.056, E_eff=−2.14<0) and survives
   pressure removal (s_max settles ~4× the bare seed) instead of dispersing.

Bugs fixed during bring-up (recorded so they aren't re-hit): (i) in-place update
with parallel neighbour reads = data race → double-buffer (Jacobi); (ii) anneal
`P=P0(1−lev/(nlevel−1))` divides by zero at nlevel=1 → guard; (iii) dtau=0.15dx²
exceeds the reaction-stable limit → default 0.04dx².

Status (HONEST — mechanism works, convergence does NOT yet): under extreme
pressure the field forms a compact BOUND state (E_eff=−2.14<0, s_max=0.056), but
when the pressure is annealed to 0 it EXPANDS into a worse configuration — the
released seed run in the kernel drifts **−8.1%** and dies (0/7 alive), vs the bare
ball's −2.93%. So the pressure trick converges *under load* but the free
(P=0) soliton is not reached. Two concrete causes to fix next:
- **Norm not held through the anneal:** Q drops 118→39 by P=0, so the
  rescale-to-M₀ is being overwhelmed by the pressure work / expansion. The frequency
  ω should be solved as the *eigenvalue* of the norm constraint (update ω each step
  from the Rayleigh quotient) rather than held fixed in (m²−ω²); fixing both ω and
  norm over-constrains and the flow leaks charge.
- **P=0 tail doesn't converge** (maxF stuck ~0.15): needs the eigenvalue update
  above plus a much longer final relaxation, or acceleration (Nesterov / Anderson).
A third possibility the data leaves open: at η=0.25 there may be NO nearby
radiation-free stationary state (the drain is fundamental at this coupling), in
which case the relaxer should be run at lower η or the ansatz changed. The tool
(`sfa/seed/eta_relaxer.c`) and the pressure-convergence method are in place; the
soliton-finder is not yet delivering a better object than the bare ball.

### Higgs "fabric-pull" self-compression — bag sector added (2026-06-26)

User direction: the compression should be INTRINSIC, not an external pressure —
a **Higgs condensate** with a VEV (the fabric "pulled taut" to a nonzero vacuum)
and a mass **gap**, whose vacuum pressure squeezes a cavity the matter digs. This
is the MIT-bag mechanism: added to `eta_relaxer` as a real scalar H with
V_H=(λ/4)(H²−v²)² and a matter coupling +(κ/2)·s·H² (matter expels the condensate
→ digs a void → the bag constant B=(λ/4)v⁴ compresses it). Maps onto the project's
own "Higgs as vacuum void" speculation (FUTURE.md, real-field era).

Result — mechanism engages, does NOT self-sustain:
- Higgs alone (no external pressure): matter disperses before it is dense enough
  to dig the void (Hmin stays at v) — a **bootstrap trap**.
- External pressure to start + Higgs to sustain: under pressure a void forms and
  E_eff goes bound, but at P=0 it disperses again. Weak coupling → shallow void
  (B≈0.1 ≪ matter pressure ~0.5–1 → too weak); strong coupling (κ=100) → matter
  expelled from the condensate *everywhere* → disperses faster. No self-bound
  bag found in the (v, λ, κ) ranges tried.

The mechanism is correct and physically grounded; the failure is convergence.
Pre-digging the void (init H with a cavity, `Rvoid` arg) was tried too and ALSO
disperses — the matter spreads out of the cavity faster than the bag holds it, and
the void then fills back to the VEV. Across every regime/protocol tried
(Higgs-alone, external-pressure-start + Higgs-sustain, strong κ, pre-dug void) the
fixed-ω/fixed-norm relaxer does NOT find a self-bound bag.

**Assessment (honest).** The relaxer-based search for a self-compressing
stationary soliton has hit a wall: at fixed ω and fixed norm the matter's gradient
energy makes it spread, and no achievable bag/η/pressure term in the relaxer holds
it without an external cage. Two readings, both worth stating: (i) the *method* is
wrong for this — a self-bound soliton likely needs the full DYNAMICAL kernel (with
the new sectors) where charge/frequency self-regulate, not a fixed-(ω,norm)
gradient flow; or (ii) the *object* may not exist at these couplings (no
radiation-free self-compressing state near the η=0.25 ball). Distinguishing them
needs the kernel to carry the Higgs (and/or 2nd-gauge) field — a scp_sim change
requiring explicit authorization. The relaxer (`sfa/seed/eta_relaxer.c`) now
contains the full coupled Φ–Θ + Higgs-bag energy and the pressure/norm/pre-dig
machinery; it is a correct testbed, but it is not producing a bound
self-compressing soliton.

### RESOLVED (2026-07-08, v72): the stationary η-soliton exists — fixed-Q flow closes the drain

The whole 2026-06-26 stationary-η-soliton arc above is superseded by v72/FINDINGS.md.
The "intrinsic η-drain drift floor," the "BVP removes only ~25%," and the relaxer's
cage-or-die behavior were artifacts of (a) the over-constrained fixed-ω+fixed-norm
flow, (b) a √3 seed mis-normalization (radial_qball's f is PER-COMPONENT; the
relaxer/BVP tools seeded f/√3 — off-shell 27× in s), and (c) an SFA semantic bug
({1,3,2}: every eta_relaxer seed loaded scrambled). The correct variational principle
— gradient flow on E_Q = E + Q²/(2N) at fixed Q, ω = Q/N the Lagrange multiplier
(`sfa/seed/eta_qflow.c`) — converges to machine floor with NO external pressure, and
its kernel drift at η=0.25 is −0.038% over T=60, BELOW the η=0 floor (−0.072%);
gauged g=0.05: −0.072% vs −1.427% for the Θ=0 seed; at η=0.5 the stationary state
drifts −0.014% vs −5.63% for Θ=0 (400×). Branch VK-stable (dQ/dω<0, Q=88–210).
Consequence for the QFI program: the drift<0.5% limiter that pinned η*≈0.25 no
longer binds for stationary seeds. **Executed same day (v72/FINDINGS.md §6):**
the single stationary ball is a proper NULL (nQFI ≤ 4e-5 at η=0.25 AND 0.5 —
so the June single-ball η-scan signal, nQFI ~ η^3–4, was the Θ-sourcing
transient, not steady entanglement); the stationary PAIR (in-phase, D≈7.6,
gauged) crosses the bound at nQFI = 3.9 (η=0.25) / 3.5 (η=0.5), peaked at the
inter-ball wavevector — 3.5× the June transient-pair value. The witness lives in
inter-particle correlation. Open: hrange normalization; Δφ pair taxonomy
(gen_sfa_pair); f_Q across fusion; flavored-GDR witness (needs multi-frequency
eta_qflow). Additional v72 branch properties: dE/dQ = ω (Legendre, 0.1–0.3%);
dressing binding −0.3% of M (same-grid control), E=mQ crossing Q* 139.7→132.6;
branch extends BELOW the η=0 window (Q_min(0.25) ∈ (80,84) vs 86.7); prolate
matter core (aspect ∝ η², 1.03–1.10) + oblate toroidal Θ belt (0.71), exactly
û-axisymmetric (ẑ-control Q20 = 0 at the P2 magic angle).

---

## Algebraic / Dynamical-Lagrangian Track (v59–v61)

*Separate from the V50/C4 simulation work below. Authoritative sources:
`v60/lagrangian/CLOSEOUT.md` (the dynamical Lagrangian) and `v61/CLOSEOUT.md`
(curved gravity + EW-vev home). All results verified with SymPy/Maxima/Lean.*

**Done:**
- v60: `ℒ_v60` on `Cl(3,1)⊗Cl(7)_even` — OBE as a connection-eliminated sector,
  2 ghost-free TT gravitons, Koide cone Q=2/3 from EL minimization, EP-exact
  coupling, stable spectrum, nonlinear time evolution.
- v61: curved-space GR (Schwarzschild from `ρ_grav`, backreaction `m'=4πr²ρ`), the
  LIGO closure (2 TT = `h₊,h×`, GW quadrupole, speed `c`), all three classic GR
  tests (deflection `4GM/b`, GW, perihelion `6πGM/(c²a(1−e²))`, Mercury 43″/cy), and
  a dynamical *home* for the EW vev `v=784a²` (Frobenius `End(L)` Higgs).

**Open (v62 candidates):**
1. Numerical self-gravitating Koide / boson star (TOV with the GEN3 potential),
   closing v61 GEN2's interior solution quantitatively.
2. R1 equipartition/democracy selection — the extra term that picks the democratic
   `End(L)` vacuum (`784a²`), or a proof it cannot be selected.
3. FRW cosmology sourced by the matter sector.
4. The four residual value-conjectures (`α`, `v=784a²`, `φ=2/9`, `f_g~α^{21/2}`) —
   inputs, not dynamical gaps.

---

## Critical Priority — V52 Quantitative Measurements

These are detailed in `v52/PLAN.md`. Summary:

### F1: Force Power Law — F(D) exponent [V44-ERA, needs retest]

V33 measured F ∝ 1/D^1.8 for braids (mixed gravity+EM). Newton requires n=2.
Need proton-proton measurement with C4 equations at D = 15..50.
**V52 Test 1**: 6 runs × 10 min GPU.

### F2: Depletion Profile δρ(r) [V44-ERA, needs retest]

V34 measured δρ ∝ 1/r^1.2 for a z-aligned braid. Need isotropic proton
measurement with C4. If n→2, the force is automatically inverse-square.
**V52 Test 2**: 1 run × 35 min GPU.

### F25: Gravitational Coupling Constant C_proton [PARTIALLY CONFIRMED — V51]

V51 gradient test (C4, steep gradient 0.15/0.05) shows proton drifts +6
code units toward low ρ over T=400. Drift rate ~0.015/t. Need multiple
gradient strengths to confirm F ∝ ∇ρ linearity and extract C precisely.
**V52 Test 3**: 3 runs × 35 min GPU (or use V51 data).

### F6: Charge-Dependent Force [V44-ERA, needs retest]

V34 measured same-winding +27%, opposite -57% for braids. Need UUD+UUD vs
UUD+UDD comparison with C4.
**V52 Test 4**: 2 runs × 20 min GPU. Requires UDD template generation.

### F26: Equivalence Principle [UNTESTED]

Test whether different-mass particles (braid, proton, deuterium) experience
the same gravitational acceleration in the same gradient. If a = const
regardless of mass, equivalence holds.
**V52 Test 6**: 3 runs × 35 min GPU.

### F_NEW: Emergent Gauss' Law [UNTESTED]

Test whether ∮ S·dA is radius-independent around a proton (S = Poynting
flux from θ radiation). Equivalently, measure θ_rms²(r) on shells —
if ∝ 1/r², Gauss' law holds for the radiation field. This is the
foundation for 1/r² Coulomb force.
**V52 Test 7**: Reuses Test 2 data + analysis tool.

---

## V51 Results (March 2026)

### Proton-Proton Collision [CONFIRMED — V51]

Two pre-converged UUD protons at D=25, v=±0.1c. V50/C4 equations.
- Protons condense by t=50, maintain distinct cores
- Closest approach D≈19 at t≈125
- **Bound state formation**: protons oscillate in separation (8-25 code units)
- Slow inspiral from θ radiation losses (25→17 over T=300)
- Cosserat hardening prevents merger even at closest approach
- Details: `v51/RESULTS.md`

**Open**: A higher-energy collision (v=0.3c or 0.5c) should be tested to
find the threshold for proton disruption vs bound-state formation.

### Proton Gravitational Drift [CONFIRMED — V51]

Single proton in steep gradient (A_high=0.15, A_low=0.05). V50/C4 equations.
- Proton drifts +6 code units toward low ρ over T=400
- Energy conservation: -0.032% over T=500
- Gradient persists (verified by slab measurement)
- Consistent with V43 results under different equations
- Details: `v51/RESULTS.md`

### Analytical Seed Failure [CONFIRMED — V51]

The V44 analytical seed (gen_proton_analytical, 512³) does NOT condense
into a proton under C4 equations. The Cosserat strain and hardening terms
prevent the six raw braids from merging. All production runs must use
pre-converged templates from V43 (evolved without C4 terms).

**Open**: Can protons form spontaneously under C4? If not, what formation
pathway works? Lower α/β during condensation, then raise to production
values? This is a fundamental question about whether the C4 equations
support particle formation or only particle stability.

---

## V50 Results (March 2026)

### Polariton EM Wave [CONFIRMED — V50]

The Cosserat curl coupling produces a hybridized θ-δφ wave (polariton).
Dispersion relation proven in Lean (zero sorrys):
- v_phase = c × √(1 − η²/m²) = 0.9428c (predicted), 0.9377c (measured)
- Mixing ratio: δφ/δθ = −ηk/m² (~18% φ admixture at λ=8)
- Propagation mechanism: θ → curl(θ) → δφ → curl(δφ) → θ (E↔B cycle)
- True eigenmode is a plane wave (infinite transverse extent)
- Background A_bg must be 0 for clean propagation tests
- Details: `v50/EM_WAVE_RESULTS.md`, `lean/Polariton.lean`

### Two-Proton Bound State [CONFIRMED — V50]

V50/C4 produces stable two-body bound state with distinct cores (V44
allowed merger into a blob). Shell structure confirmed: core → hardened
shell → exterior. Details: `v50/RESULTS.md`

### Single-Pass Force Expansion [VERIFIED — V50]

The two-pass Cosserat+hardening force computation can be algebraically
expanded to a single GPU pass, eliminating 6×N³ intermediate arrays.
Verified in Maxima (10⁻¹⁸ precision) and Lean (zero sorrys).
Details: `v50/em_wave/single_pass.mac`, `lean/SinglePass.lean`
**Status**: Not yet implemented in the CUDA kernel.

---

## High Priority

### F17: Nuclear Binding Energy — ³He and ⁴He [V44-ERA]

²H confirmed (V42, V51). Next: ³He (2 UUD + 1 UDD) and ⁴He (2 UUD + 2 UDD).
Need N=512+ grids. Requires pre-converged UDD template (currently missing
for C4 equations).

### F18: Neutron Stability Under C4 [UNTESTED]

UDD behavior will change with Cosserat constraint. The V44 result (UDD
survives T=1000) needs verification with C4.

### F19: Force Equilibration [V44-ERA, needs retest]

The V42 strong/EM ratio equilibration (259:1 → 1:1) was measured with V44
equations. C4's Cosserat constraint changes θ dynamics fundamentally.

### F_NEW: High-Energy Proton Collision [UNTESTED]

V51 tested v=0.1c (kinetic/rest ~0.5%). Test v=0.3c and v=0.5c to find:
- At what energy do protons break apart permanently?
- Is there a transition from elastic bounce → bound state → disruption?
- Do the fragments resemble known particles (pions, kaons)?

### F3: Lorentz Contraction Verification [UNTESTED]

Boost a proton at v=0.1c, 0.3c, 0.5c. Verify γ contraction and time
dilation. The C4 equations are Lorentz-invariant (Lagrangian-derived),
so boosted solutions must transform correctly.

### F4: Isotropic Background [UNTESTED]

Test braids/protons in a random-phase (non-z-aligned) background.
If the theory requires the z-preferred background, there's a hidden
preferred-frame problem.

---

## Medium Priority

### F22: Breathing Mode Spectrum [V44-ERA]

FFT of per-shell ρ(t) to reveal full mode spectrum. The V50 two-proton
data shows compound breathing (4.5-unit beat pattern).

### F24: Mass Defect Measurement [V44-ERA, redesign for C4]

Differential method: same seeds, different D. The reference equation set
is now C4. Need pre-converged templates and T=500+ equilibration.

### F_NEW: Shell Thickness vs Parameters

Sweep β (hardening) at {0.1, 0.3, 0.5, 1.0, 2.0} with α=0.1 fixed.
How does shell thickness depend on β? Does it affect force law or binding?

### F_NEW: Proton Formation Pathway Under C4

The analytical seed failure (V51) raises the question: can protons form
spontaneously under C4? Test:
1. Start with V44 equations (α=β=0), form proton
2. Slowly ramp α and β to production values over T=200
3. Does the proton survive the transition? At what ramp rate?

### F21: Group Velocity Subluminal Inside Breathing Oscillator [UNTESTED]

Verify |∂φ/∂t| > c at antinodes is phase velocity, not group velocity.
Launch a δφ perturbation at the proton edge, track arrival at far side.

### F9: Analytical Effective Potential

Derive V_eff(D) from linearized perturbation theory. Would predict the
force law without simulation.

---

## Lower Priority / Speculative

### F10: Spin and Helical Handedness

Enhanced by V50: Cosserat constraint makes θ a true geometric quantity.
Spin identification more natural.

### F11: Gravitational Waves

Accelerated braid radiation — does it have spin-2 tensor structure?

### F12: Dark Matter Profiles

Does the depletion profile match NFW/Burkert halos? Quantitative comparison.

### F13: Quantization

Classical theory → QFT. Triple-product potential well-defined quantum mechanically?

### F14/F15: Cosmological Constant / Background Origin

The Cosserat strain energy E_strain = α|M|² contributes to vacuum energy.

### Speculative: Vacuum Refractive Index and Dark Matter

The polariton speed v < c gives the vacuum a refractive index n = c/v ≈ 1.06.
Near massive objects, depletion lowers A_bg → n closer to 1. The spatial
variation of n produces lensing that mimics dark matter halos. See
previous FUTURE.md for full discussion of MOND-like behavior, Cherenkov
radiation limit, chromatic micro-lensing smoking gun, and heliospheric
implications.

### Speculative: Higgs as Vacuum Void, Goldstone as Carrier Phase

**The vacuum is not empty.** The background φ field oscillates at
amplitude A_bg with carrier wave cos(k·z + δ + Δ). This nonzero
vacuum state is the analog of the Higgs VEV in the Standard Model.

**Higgs boson as absence of field.** A point-like region where A_bg → 0
would be:
- **Massive**: the surrounding field pushes inward to fill the void.
  The restoration force (from V(P) wanting the vacuum to remain at
  A_bg) gives it an effective mass m²_H ~ V''(P_bg).
- **Scalar (spin-0)**: no angular structure, no chirality, no carrier
  phase. Just an amplitude zero.
- **Unstable**: the void fills in as surrounding field radiates inward.
  Short lifetime — analogous to the Higgs boson's ~10⁻²² s.
- **Point-like**: unlike a proton (extended composite with internal
  structure), the void is a localized deficiency.

Contrast with the proton's depletion zone, which is a PARTIAL absence
(A_bg reduced ~2%) and is long-range. The Higgs would be TOTAL absence
— much more energetic, much more localized.

**Goldstone boson as carrier phase mode.** The carrier wave phase
can shift uniformly at zero energy cost — this is the spontaneous
breaking of translational symmetry. Long-wavelength phase fluctuations
are massless excitations: Goldstone bosons.

**Goldstone eaten by gauge coupling.** The η×curl term mixes the φ
phase mode (Goldstone) with the θ field (gauge analog). The resulting
polariton has a mass gap (v = 0.94c, not c). This mirrors the Higgs
mechanism: the Goldstone is "eaten" to give mass to the gauge boson.

Mapping:
- Higgs VEV → A_bg (nonzero vacuum amplitude)
- Higgs boson → point-like void where A_bg → 0
- Goldstone boson → carrier wave phase mode
- Goldstone eaten → polariton mass gap from η coupling
- W/Z bosons → θ field (massive through φ-θ mixing)

**The V(P) potential supports this.** V(P) = (μ/2)P²/(1+κP²) with
μ < 0 has nontrivial curvature at the background. The κ parameter
breaks scale invariance (analog of the trace anomaly in QCD). Without
κ, V is scale-free; with κ, there's a preferred scale at κP² ≈ 1.

**Proton mass decomposition analog** (cf. Ji decomposition in QCD):
- Quark kinetic energy (~32%) → φ breathing kinetic energy ½|∂φ/∂t|²
- Gluon field energy (~36%) → θ field energy ½|∂θ/∂t|² + ½|∇θ|²
- Trace anomaly (~23%) → V(P) saturation (κ breaks scale symmetry)
- Quark condensate (~9%) → background amplitude A_bg

~91% of the proton mass is energy, not Higgs mechanism — in our theory,
~91% would be breathing + θ radiation + V(P) binding, with only ~9%
from the background field (A_bg coupling).

**Testable in simulation:**
1. Punch a hole in the background (set A_bg → 0 in a small region)
   and measure the decay time and products.
2. Measure the radial mode frequency — excite a uniform δA perturbation
   and FFT the response. The frequency is the Higgs mass analog.
3. Measure long-wavelength phase fluctuations to confirm they're
   massless (Goldstone).
4. Verify the polariton mass gap matches η×curl prediction (already
   done: v = 0.94c in V50).

---

## Resolved / Abandoned

### Confirmed (no further action)
- F2: Gravity mechanism — dynamic footprint, NOT energy minimization [V33/V43/V51]
- F5: Photon speed — v = 0.9377c, matches polariton prediction [V50]
- F20: Intermediate phase group in deuterium — null result, not a bond [V43]
- Composite particles: UUD stable T=500, phase confinement works [V41]
- Nuclear binding: deuterium UUD+UDD bound T=500 [V42]
- Two-proton bound state with C4 [V50/V51]

### Abandoned
- X1: S/B two-component split (artificial)
- X2: c(ρ) speed modification (not needed)
- X3: Binding-weighted gradient coupling (net repulsive)
- X4: SPH particle-based simulation (poor conservation)
- X5: Rotating FRW expansion (no braid formation)
- X6: Asymmetric drag hypothesis (refuted by V33 drag test)
- X7: Single-particle density-dependent κ collapse (theta prevents it)
- X8: Random parameter sweep (superseded by first-principles V41)

---

## v67: The Theta-Boundedness Problem (2026-06-10)

The complexified kernel's Θ is massless, potential-free, AND charged under the
diagonal U(1) — a combination nature never allows (photon: neutral; gluon:
confined; W: massive). The measured η-drain is the dynamical symptom: an
unprotected channel bleeding every ball to the marginal state. If left unbounded,
all particle masses become epoch-dependent (drain to threshold) — fatal for the
v59 Koide bridge. Resolution routes (task #10):

1. **m_θ > ω (W-route)**: Yukawa-bound dressing, drain kinematically closed.
   v67 runs th1 (m_θ=1.6) / th2 (m_θ=0.7) test this directly.
2. **Cosserat binding (geometry route)**: complexified α|∇×Φ/2 − Θ|²
   (U(1)-invariant, v66 THEORY §1.1) — θ locked to matter's torsion, gap √(2α).
3. **Gauge the diagonal U(1)**: θ → connection; amplitude unboundedness = pure
   gauge; Gauss law pins charge to matter. The η curl coupling is already
   A·J-shaped. Deepest fix; design doc before any code.
4. **sigma_cross complexified**: density-dependent θ mass (interim).

## v73: The Gauged Spinning Ring Program (2026-07-08)

The spinning Q-ring (spin = winding of real-space circulation, J = nQ =
n·ℏ_eff) survives gauging and is the theory's first charge+spin+magnetic
object (PROCESS.md §4–5; DISCOVERIES V73). The η fabric coupling, by
contrast, fissions rings (no exact rotation symmetry at η>0 — the
symmetry audit of PROCESS.md §5.3). Open experiments, in value order:

1. **Magneton / g-factor analog** — DONE 2026-07-09 (`v73/analysis/
   magneton.c`, PROCESS.md §5.4): the gauged ring self-magnetizes within
   ~60 t.u. and holds B_z(0) ≈ −0.035 and trapped hole flux ≈ −0.38 (units
   of g) steady to t = 120; g_factor = μ_z/(Q·L_z/2M) = 1.05–1.06
   (classical extended rotor with a +6% distribution anomaly; NOT Dirac
   g = 2); E_B + E_E reproduces the kernel E_em to 4 decimals. Follow-ups:
   link-sign convention audit (μ vs B relative sign), box-size control on
   the return flux, μ on the exact gauged (Q,J) state.
2. **η + gauge coexistence window** (top open experiment) — map
   lifetime(η, g) for rings and balls: is there a window where gauge hoop
   stress outruns the η fission drive (e-fold ~30 t.u. at η = 0.25, g = 0)?
   Success: L_z/Q flat and A_v/A_u ≈ 1 for ≫10² t.u. at some (η > 0,
   g > 0). Decides whether η is a viable dressed-matter coupling or only a
   transient/spin-orbit channel.
3. **Skeptical QFI suite** — hrange scan; operator suite (ρ_Q, s, |Φ|²,
   per-component ρ — map which observables null on single solitons);
   relative-phase Δφ taxonomy via gen_sfa_pair; time-window dependence;
   phase-randomized / time-shuffled surrogates. Success: single-null +
   pair-crossing stable across the suite. Falsifier: crossing only at
   hrange = 1 or only for one operator → protocol artifact, not physics.
4. **Gauged fixed-(Q,J) relaxer** — extend eta_qflow's fixJ mode
   (PROCESS.md §5.2) to complex_gauge = 1 for the EXACT stationary gauged
   spinning ring; map the spin branch and its 2×2 VK matrix
   det ∂(Q,J)/∂(ω,Ω) — the multi-charge stability criterion.
5. **n = 2 ring** (J = 2ℏ_eff) gauged; ring–ball and ring–ring collisions
   (does spin transfer in integer units?).
6. Twist-constrained relaxation for the (n,m) taxonomy (low priority
   unless chirality becomes a goal — the twist is a different quantum
   number from spin; carried over from §4.2).
