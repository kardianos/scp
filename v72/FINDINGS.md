# v72 — The Stationary η-Soliton: the η-Drain Closed by the Fixed-Q Variational Principle

**Date**: 2026-07-08
**Question (user)**: the recent versions applied the stability heuristics (CONCEPT.md §2–3,
derived in v65) with partial success — can a new approach based on those heuristics reach
FULL stability?
**Answer**: YES. Applying the heuristics' own variational principle — minimize E at fixed
charge Q with ω as the Lagrange multiplier — to the coupled Φ–Θ system produces an exact
stationary η-coupled Q-ball whose kernel drift at η=0.25 is **below the η=0 floor**
(−0.038% vs −0.072% over T=60). The June-26 "intrinsic η-drain drift floor" was not
physics: it was two tooling bugs plus an inconsistent initial condition.

## 1. The diagnosis — why the June-26 arc stalled below full stability

The 2026-06-26 campaign (FUTURE.md) concluded that an intrinsic η-drain floor capped every
η>0 configuration, that the BVP transverse-Θ seed removed only ~25% of it, and that the
residual was a "2D Φ back-reaction" needing a coupled relaxer — which then failed unless
externally caged, and collapsed on release (−8.1%, died). Three distinct defects produced
that picture:

1. **Wrong variational principle in `eta_relaxer`.** It fixed BOTH ω (in m²−ω²) and the
   matter norm (rescaling each step). The Q-ball stability heuristic is "minimize E at
   fixed Q"; fixing both over-constrains the flow, which then leaks charge and needs the
   pressure cage. (This fix was already identified in FUTURE.md but never implemented.)
2. **Seed normalization off by √3.** `radial_eta_soliton` and `eta_relaxer` seeded
   Φ_a = f/√3, but the `radial_qball` profile f is the PER-COMPONENT amplitude
   (`gen_qball_boost` convention). Their seeds were off-shell by 27× in s = Π|Φ_a|²;
   the relaxer's "Q drops 118→39" is exactly the factor 3 in norm. All June-26
   η-drift baselines (−2.93/−2.66/−2.27%) were measured on this wrong ball.
3. **SFA semantic-code bug in `eta_relaxer`.** Its column semantics were {1,3,2} =
   {ANGLE, ACCELERATION, VELOCITY} instead of {0,1,2} = {POSITION, ANGLE, VELOCITY}:
   every relaxer seed loaded into the kernel with Φ landing in Θ and Θ discarded.
   The "released seed drifts −8.1% and dies" run was scrambled on load, not unstable.

All three defects are now addressed: (1) the correct variational principle is
implemented in the new tool below (`eta_relaxer` keeps its fixed-ω flow and is
superseded for this use); (2) the √3 normalization is fixed in `eta_relaxer.c` AND
`radial_eta_soliton.c` (both now seed the per-component f, BVP source −√3ηf′ —
verified to reproduce the identical g(r) to 1e-9); (3) the SFA semantics are fixed
in `eta_relaxer.c`.

## 2. The new tool — `sfa/seed/eta_qflow.c` (fixed-Q gradient flow)

Rotating ansatz Φ_a = φ_a(x)e^{iωt}, Θ_a = ϑ_a(x)e^{iωt}. The conserved U(1) charge is
Q = ω·N with N = ∫Σ_a(|φ_a|²+|ϑ_a|²). Eliminating ω gives the fixed-Q functional

    E_Q = ∫ Σ_a [ ½|∇φ_a|² + ½m²|φ_a|² + ½|∇ϑ_a|² + ½m_θ²|ϑ_a|² ]
          + Vt(s) − η·Re[ϑ̄·(∇×φ)]  +  Q²/(2N)

whose Q²/(2N) term is literally the charge-pressure of the stability heuristic (it
contributes +ω²·field to the force with ω ≡ Q/N). Plain gradient flow on E_Q — no norm
rescaling, no pressure cage, ω recomputed as Q/N each step — is monotone descent, and its
fixed points are EXACTLY the kernel's stationary rotating states on the same stencil
(7-pt Laplacian + central-difference curl, copied verbatim from
`compute_forces_complex`). Θ is initialized from the l=1 transverse BVP
(source −√3·η·f′, the √3 restoring the correct normalization).

Convergence (N=64, L=15, ω₀=1.45 profile, η=0.25, Q=119.93): monotone E_Q descent,
maxF 2.5e-2 → 9.1e-15 by 16k iterations, 5.4e-15 at 20k (machine floor;
`results/qflow_w145_eta025_20k.log`), ~4 min on 16 CPU cores. The banked seed
`results/qflow_w145_eta025.sfa` is an 8k-iteration re-run (maxF 6.4e-10,
`results/qflow_w145_eta025.log`) regenerated after the semantics fix of §1(3);
ω and E_Q agree with the 20k run to all printed digits. Converged state:
ω=1.444635, Q_φ=118.05, Q_θ=1.89, E_Q=180.44, s_max=0.05196 (slightly denser
than the bare ball's 0.0500).

## 3. Kernel verification — drift comparison (N=64, L=15, T=60, absorbing BC, f32)

All runs identical except the seed and η; drift = ΔE_total/E_total from diag.tsv.
Ungauged (complex_phi=1) and gauged (complex_gauge=1, g=0.05) packages:

    run          η     g     seed                 drift(0→60)  dQ_total   s_max 0→60
    r_floor      0     0     ball, Θ=0            −0.072%      −0.068%    0.0500→0.0501
    r_theta0     0.25  0     ball, Θ=0            −1.361%      −1.178%    0.0500→0.0522
    r_bvp        0.25  0     ball + BVP-Θ         −0.111%      −0.097%    0.0500→0.0520
    r_qflow2     0.25  0     eta_qflow stationary **−0.038%**  −0.040%    0.0520→0.0519
    rg_theta0    0.25  0.05  ball, Θ=0            −1.427%      −1.236%    0.0500→0.0510
    rg_qflow2    0.25  0.05  eta_qflow stationary **−0.072%**  −0.072%    0.0520→0.0503
    r5_theta0    0.5   0     ball, Θ=0            −5.634%      −4.932%    0.0500→0.0551
    r5_qflow     0.5   0     eta_qflow stationary **−0.014%**  −0.015%    0.0572→0.0571

- **The η-drain is fully closed.** The exact stationary state at η=0.25 drifts LESS than
  the η=0 ball itself: the qflow seed is an exact *discrete* stationary state, whereas
  the η=0 ball is the continuum ODE profile interpolated onto the grid (small
  discretization mismatch that radiates). Gauged, the stationary seed matches the
  ungauged floor exactly (−0.072%); the 20× residual improvement is bounded by the
  absorbing-BC/f32 floor, not by η physics.
- **The BVP-Θ seed alone already removes ≈92% of the drift** once the Φ normalization
  is correct (−1.361% → −0.111%; 97% of the above-floor drain) — the June-26
  "only ~25%" number was the √3 bug.
- **The closure holds at strong coupling.** At η=0.5 (the old QFI-feasibility wall,
  naive drift ∝ η²: −5.63% for the Θ=0 seed) the stationary state drifts −0.014% —
  400× cleaner and the lowest drift of the whole campaign. The relaxed η=0.5 state
  (ω=1.430696, Q_θ=7.73 = 6.2% dressing ∝ η², s_max 0.0572) is confirmed in-kernel to
  5 decimals in ω_core (1.43075). There is no intrinsic η-drain at ANY tested coupling —
  the drain was always the off-shell (Θ=0 / mis-normalized) initial condition.
- **The kernel confirms the Lagrange multiplier**: measured omega_core = 1.44464 at t=30
  vs the relaxer's ω = 1.444635 — agreement to 5 decimals. Q_core ≈ 119.3 ≈ Q_target.

## 4. The η-coupled branch is VK-stable (heuristic 4)

`eta_qflow` at fixed Q along the branch (N=48, L=12, η=0.25, warm-started from the
closest-Q η=0 profile; `results/vk_branch.txt`):

    Q_total   ω(Q)       Q_θ (dressing)
    88        1.467912   1.56
    95        1.460697   1.63
    105       1.453145   1.73
    118       1.445586   1.86
    140       1.435820   2.07
    170       1.426119   2.34
    210       1.416710   2.68

dQ/dω < 0 across the whole scanned range — classically stable (Vakhitov–Kolokolov). The
branch is smooth; the Θ dressing carries a slowly-growing ~1.5–2% charge fraction. At
fixed Q the η=0.25 dressing lowers ω relative to η=0 (e.g. Q=118: 1.4456 vs 1.4500) —
the cross term adds binding.

## 5. Proven properties of the stationary branch (2026-07-08, second campaign)

**Legendre / canonical-family identity.** Along the η=0.25 branch (§4 data plus E_Q
from the vk logs), centered differences give dE/dQ = ω(Q) to 0.26% → 0.10% (error
falling as the sampling tightens — finite-difference truncation). The fixed points
are not merely numerically quiet: they form the thermodynamically consistent
one-parameter Q-ball family obeying the canonical soliton relation.

**Binding crossover + dressing binding energy (same-grid control).** η=0
relaxations at Q=118/140 on the identical N=48 lattice (`vk0_Q*.log`) isolate the
dressing contribution from discretization:

    Q     E(η=0)     E(η=0.25)   ΔE_dressing        ω(η=0)→ ω(η=0.25)
    118   178.1273   177.6009    −0.526 (−0.30%)    1.44943 → 1.44559
    140   209.8590   209.2530    −0.606 (−0.29%)    1.43926 → 1.43582

The Θ dressing binds ≈0.3% of the mass at fixed Q, and moves the evaporation-bound
crossing E = mQ from Q* ≈ 139.7 (η=0) to Q* ≈ 132.6 (η=0.25) on the same grid.
(The continuum η=0 branch puts Q* ≈ 145; the ~5-unit offset from 139.7 is N=48
discretization — hence the same-grid control.) E/Q falls below m between Q=118
(1.5051) and Q=140 (1.4947) on the dressed branch.

**Branch endpoint — the dressing EXTENDS the existence window.** Probes at
Q = 87, 86, 84 all converge (ω = 1.4692, 1.4707, 1.4743; Q_θ ≈ 1.5; maxF at the
flow floor, though Q=84 descends ~100× slower — the critical slowing of an
approaching fold), while at Q=80 the flow finds NO solution: ω runs to
1.499980 ≈ m, the Θ dressing evaporates (Q_θ → 0) and the field delocalizes.
So **Q_min(η=0.25) ∈ (80, 84) — below the η=0 minimum Q_min ≈ 86.7** (ω ≈ 1.4825):
the η dressing extends the branch to lower charge, consistent with the binding it
adds (at ω = 1.4743 the dressed ball holds Q_φ = 82.5 where the bare branch needs
≈ 89.5). dQ/dω < 0 on the entire found branch (§4), all VK-stable side.

**Axisymmetric geometry, measured** (`analysis/seed_aniso.c`, results
`aniso_eta{025,05}.txt`). s-weighted quadrupole about û=(1,1,1)/√3:

    seed      Q20[s]    aspect(s)   Q20[|Θ|²]  aspect(Θ)
    η=0.25    +0.0208   1.0256      −0.392     0.713
    η=0.5     +0.0830   1.1047      −0.372     0.730

The matter core is PROLATE along û — deformation scaling as η² (Q20 ratio 4.0 for
2× η) — wrapped in an OBLATE toroidal Θ belt. The ẑ-axis control gives
Q20 = 0.00000 and aspect = 1.00000 identically: û·ẑ = 1/√3 is exactly the
P2(cos α) = 0 magic angle, so a vanishing ẑ-quadrupole is a sharp proof that the
deformation is purely axisymmetric about û. The June-26 "the radiation-free
soliton is genuinely 2D axisymmetric" prediction is now quantitative.

**Long-time stability (T=1000, η=0.5, ungauged, N=64 L=15;
`results/r5_long_diag.tsv`).** The strongest-coupling stationary state meets the
project's 10³-t.u. standard: s_max 0.0572 → 0.0571 (flat to 0.2%), E −0.245% and
Q_total −0.255% over the full run (a slow LINEAR boundary-tail absorption,
−2.5e-4%/t.u. — no acceleration, no mode growth), and omega_core = 1.430773 at
t=1000 vs the relaxer's ω = 1.430696 — the internal clock still locked to the
Lagrange multiplier after ~230 rotations.

## 6. Expansion — the QFI witness on stationary states

Four gauged fine-cadence runs (g=0.05, N=64, L=15, T=40, snap_dt=0.25, `sfa_qfi
--auto-T`, hrange=1 as in the June-26 scans; outputs on /space/scp/v72/, TSVs in
`results/qfi_*.tsv`). Pair seeds: `gen_sfa_pair` superposition of the banked
stationary ball at ±D/2 along x, in-phase (D_actual = 7.62).

    run                     drift(0→40)   nQFI[ρ_Q] by q_x mode:  0        1        2       3       4
    single η=0.25 (stat.)   −0.038%       5.5e-8   4.8e-6   3.0e-5  4.2e-5  2.8e-5
    single η=0.5  (stat.)   −0.021%       5.9e-7   2.7e-6   1.5e-5  2.4e-5  1.8e-5
    pair D=7.6 η=0.25       −0.848%       9.5e-4   0.167    1.909   3.876   2.819
    pair D=7.6 η=0.5        −0.449%       6.5e-3   0.130    1.679   3.478   2.615

1. **The single stationary ball is a proper null** — nQFI ≤ 4.2e-5, ~20× below the
   June-26 transient single-ball value at the same η (9.1e-4). The η-scan's
   single-ball witness was therefore dominated by the Θ-sourcing TRANSIENT, not by
   steady-state entanglement content; the theory requirement that one classical
   soliton factorizes (witness ≈ 0) is now confirmed on the true eigenstate at
   both couplings.
2. **The two-ball bound-crossing is real and 3.5× stronger on clean states:**
   nQFI[ρ_Q] = 3.88 (η=0.25) / 3.48 (η=0.5) peaked at q=3–4 — exactly the
   inter-ball wavevector (2π/D ↔ mode n ≈ 3.9) — versus 1.12 on the June drifting
   pair. Cleaning the constituents removes the fake single-body signal and
   *raises* the genuine inter-particle witness. The pair drift (−0.4 to −0.8%) is
   the physical in-phase attraction dynamics (s_max redistributes, Q_θ grows),
   not seed shedding.
3. Operating point via --auto-T: T_eff ≈ 21 with w_peak = 0.156 (the pair beat /
   collective mode band) — un-saturated kernel, as required by the program's
   item 3. The hrange normalization remains the open modeling choice (same as
   June); all comparisons here are at hrange=1.

## 7. What this means for the QFI program

The June-26 "exposed tension" — entanglement (QFI) needs η, but η radiates, so the
feasible optimum was pinned at η*≈0.25 by the drift<0.5% limiter — is dissolved,
and §6 executes the re-run. Net reading of the program after v72:
- the witness lives in INTER-PARTICLE correlation (pairs: nQFI ≈ 3.9 ≫ 1), not in
  any single-soliton steady state (nulls at ≤ 4e-5) — the classical→quantum bridge
  is a statement about interacting configurations, exactly as CONCEPT.md §9's
  experiment 3 hypothesized;
- the single-ball η-scan scaling (nQFI ~ η^3–4) needs reinterpretation: it measured
  the amplitude of the Θ-sourcing transient, not steady entanglement content;
- open next steps: hrange normalization (the remaining modeling choice), the
  Δφ-phase taxonomy of the pair witness (gen_sfa_pair takes a relative phase),
  fusion tracking of f_Q across a merger (experiment 3), and the flavored-GDR
  witness (experiment 4) from a stationary flavored dressing (needs the
  multi-frequency extension of eta_qflow — the ω_a-partitioned constraint).

## Artifacts

- `sfa/seed/eta_qflow.c` (new tool; in Makefile SEED_TOOLS with eta_relaxer,
  radial_eta_soliton)
- `results/qflow_w145_eta025.sfa` + `.log` (banked 8k seed), `_20k.log` (machine-floor
  convergence evidence); `results/qflow_w145_eta05.sfa` + `.log` (η=0.5 state)
- `results/{r_floor,r_theta0,r_bvp,r_qflow2,rg_theta0,rg_qflow2,r5_theta0,r5_qflow}_diag.tsv`,
  `drift.py`
- `results/vk_branch.txt` (the ω(Q) branch)
- Bug fixes: `eta_relaxer.c` (SFA semantics {1,3,2}→{0,1,2}; f/√3→f),
  `radial_eta_soliton.c` (f/√3→f with −√3ηf′ source; verified g(r) identical to 1e-9)
