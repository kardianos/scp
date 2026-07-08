# Kernel Extension SPEC — Higgs "fabric-pull" + 2nd gauge (authorized 2026-06-26)

Goal: test whether intrinsic self-compression (Higgs bag and/or a 2nd attractive
gauge) stabilizes an η-coupled soliton in the DYNAMICAL kernel, where the relaxer
(fixed ω/norm) could not. Both sectors **default OFF** (higgs_v=0, g2=0) and the
off-path must be byte-identical to the current kernel.

## Conventions
- Standard package is gauged complex (`complex_phi=1, complex_gauge=1`). Both new
  sectors are implemented in the **gauged-complex path only** for now
  (`compute_forces_complex_gauge`, `compute_energy_complex_gauge`, `verlet_step`).
  Using them without complex_gauge=1 → warn + ignore.
- Memory: single pool, `nblocks` = 54 (gauged). Append new blocks at 54+.

## Sector 1 — Higgs condensate H (real scalar, dynamical)  [IMPLEMENT FIRST]
Lagrangian add:  ½(∂H)² − V_H(H) − (κ_H/2)·s·H²,
  V_H = (λ_H/4)(H² − v_H²)²  (VEV v_H, gap m_H=√(2λ_H)·v_H),  s = ∏_a|Φ_a|².
EOMs (κ_H couples matter↔condensate; matter digs a cavity, bag B=(λ_H/4)v_H⁴ squeezes):
  Ḧ      = ∇²H − λ_H H(H²−v_H²) − κ_H s H
  Φ̈_a   +=  − κ_H H² Φ_a · ∏_{b≠a}|Φ_b|²     (same algebraic form as the Vt force)
Config keys: `higgs_v` (0=off), `higgs_lam` (default 1), `higgs_kap` (default 2).
Memory: +3 blocks when higgs_mode: H=54, H_vel=55, H_acc=56 → nblocks 57.
Grid: `double *H,*H_vel,*H_acc; double higgs_v,higgs_lam,higgs_kap; int higgs_mode;`
Init: H ≡ v_H everywhere (taut fabric); pre-dug cavity optional later.
Force (in compute_forces_complex_gauge pass C):
  - add `-KAPH*hh*fu[a]*prod_rest[a]` (and fv) to phi_acc/phi_im_acc, hh=H[idx]².
  - compute `H_acc[idx] = lapH - LAMH*H*(H*H-VH*VH) - KAPH*s*H`, lapH = 7-pt.
Integrate: leapfrog H like phi (verlet_step): H_vel += dt*H_acc; H += dt*H_vel (match
  the existing kick/drift order exactly).
Energy: add ½H_vel² + ½|∇H|² + (λ/4)(H²−v²)² + (κ/2) s H² to E_total.
SFA: one column `higgs` (SFA semantic scalar) when higgs_mode (extra column;
  analysis tools ignore unknown columns).
Verify: (a) higgs_v=0 → diag byte-identical to current; (b) energy conserved to
  integrator floor with higgs on; (c) Gauss law unchanged (Higgs is U(1)-neutral).

## Sector 2 — Second gauge B (relative U(1), attractive)  [IMPLEMENT SECOND]
Gauge a relative color charge (T8 = diag(1,1,−2) or T3 = diag(1,−1,0)). Mirror the
existing compact-U(1) link sector: link angles th2_i, Efield2_i, E_acc2_i, scratch
cos/sin/plaq → +18 blocks (like blocks 36–53). Covariant derivative gains
+i g2 q^B_a B_i with per-component charge q^B_a from the chosen generator (opposite
signs → attraction). Second discrete Gauss law div E2 = g2 ρ2 must hold at the 1e-13
floor (tripwire). g2=0 → byte-identical. This is the linked-flux self-compression:
the two fluxes (diagonal A + relative B) can Hopf-link and pin the knot.

## Order of work
1. Higgs config+memory+struct+init (scaffolding, off=identical) → build.  [DONE]
2. Higgs force+integrate+energy → build, verify.  [DONE — CPU kernel only]
3. Run: bare ball + higgs on, does it self-compress/stabilize?  [STARTED]
4. Then Sector 2 (2nd gauge) as a separate verified increment.  [TODO]

## STATUS (2026-06-26) — Sector 1 landed in CPU kernel and VERIFIED
`scp_sim.c` now carries the Higgs condensate (config keys higgs_v/lam/kap; Grid
H/H_vel/H_acc at blocks 54-56; init at VEV; force back-reaction + H wave eq in
`compute_forces_complex_gauge`; leapfrog in `verlet_step`; energy in
`compute_energy_complex_gauge`). Verified on the ω=1.45 ball, N=64, T=20, gauged:
- higgs OFF: drift −0.056%, gauss 7.3e-14 — byte-identical baseline.
- higgs ON (v=1, λ=2, κ=20): energy conserves early (−0.015% before the dispersing
  matter reaches the absorbing BC), **gauss 7.2e-14** (Higgs is U(1)-neutral ✓),
  charge 118.02→117.56 conserved. Implementation is self-consistent.
Physics: at (v=1, κ=20) the ball DISPERSES (s_max 0.05→6e-5) — the bag does not
self-compress at this point, consistent with the relaxer.

## 4-node + Higgs experiment (2026-06-26, user direction)
Combined the CMA 4-node winner (`v71/fitness/results/best.gen`) with the Higgs +
theta and compressed. Added: genome seeding to `eta_relaxer` (Φ from a `.gen` file)
and a central-cavity Higgs init to the kernel (`higgs_rvoid`). Findings:
- **The bare 4-node self-compresses** (kernel, Higgs OFF): s_max 0.023 → **0.457**
  (≈10× a single ball) as the 4 nodes merge into a dense lump, then partially
  rebounds (0.077 at t=16). Intrinsic compression of the geometry — no Higgs needed.
  THIS is the real self-compression signal.
- **Higgs HURTS at every simple setting:** repulsive bag (κ>0, uniform OR central
  void rvoid=2.5) disperses the cluster (s_max→3e-3); attractive (κ<0, the literal
  "fabric pull") RUNS AWAY to NaN in ~20 steps (matter+condensate grow together
  unbounded). Neither sign gives stable assisted compression.
- The relaxer (external pressure) DOES hold the 4-node+Higgs (s_max~0.011 flat to
  P=0) — but the free kernel does not sustain it.

Read: the 4-node geometry compresses by itself but only transiently (collapse →
rebound), and a plain Higgs cross-term can't catch it. A stable assist needs either
a SATURATING attractive coupling (−κ s H²/(1+γ s H²), bounded so no runaway) that
engages only at the compressed core, or the **2nd gauge** (linked-flux topological
pinning) which cannot run away because the binding is topological, not energetic.
This is fresh evidence for prioritizing Sector 2.

Next: (i) saturating attractive Higgs coupling (1-line force change) to catch the
4-node at peak compression; (ii) **Sector 2 — the 2nd gauge** (the topological,
runaway-proof self-compression); (iii) SFA `higgs` column + CUDA port (deferred).

## Sector 2 attempt — relative-Coulomb (A0) gauge, NEGATIVE (2026-06-26)
User: "find a gauge that enhances the 4-node selector in a 3d encoding." First
established (kernel, η×g scan on the raw 4-node) that NO existing gauge enhances:
the diagonal U(1) (Coulomb) DISPERSES the 4-node (worse at higher g; g=0 lowest
drift), the η-torsion is roughly neutral; the cluster just breathes
(collapse→rebound) and slowly disperses (~1%/16tu).

Then prototyped (in `eta_relaxer`) a relative (color) U(1) **A0/Coulomb** second
gauge: component a's phase gauged, ω_a=ω+g2·q_a·B, effective mass m²−ω_a², B
sourced by the relative charge Σq_a ω_a|φ_a|² (T8 charges 1,1,−2). Result NEGATIVE:
- small g2 → no effect (back-reaction is O(g2²): g2≤0.1 changes s_max by <0.01%);
- large g2 (≥0.3) → RUNAWAY (the A0 self-energy −½ω_a²|φ|² is unbounded below as B
  grows → maxF→1e173). No window where it enhances.

POSITIVE byproduct: at g2=0, NO external pressure, NO Higgs, the 4-node relaxes to
a genuinely BOUND state (E_eff=−45, s_max=0.069) — the symmetric ball never did.
The 4-node geometry is a real bound configuration; it's just LOW-charge (Q~0.2 in
the kernel) and re-compresses-then-disperses dynamically (transient).

Conclusion: the A0/Coulomb second gauge is the WRONG kind (unbounded → runaway).
The remaining candidate is a **MAGNETIC / spatial second gauge with linked flux**
(Hopf-linked A and B fluxes threading the 4-node) — the magnetic energy is
positive-definite (runaway-proof) and the binding is TOPOLOGICAL, matching the
user's original "linked flux / 3d encoding." That is a full compact-U(1) link
sector (≈18 blocks, 2nd Gauss law) — the real Sector 2 build.

## Sector 2 BUILT — compact-U(1) second gauge (2026-06-26, authorized)
Implemented in `scp_sim.c` (CPU): a second compact-U(1) gauge mirroring the first
(link angles th2, Efield2, E_acc2, scratch; 18 blocks at base2=54+higgs). Matter
couples via the COMBINED covariant link per component, phase th[d]+q_a·th2[d]
(relative/color charges q_a, default T8=(1,1,-2)); J2=Σq_a Im[Φ*UΦ] sources E2;
magnetic staple from th2 plaquettes; E2/th2 leapfrog. Config: g_gauge2 (0=off),
gauge2_q0/q1/q2. Exact at η=0 (no seagull); the η-seagull for gauge2 is TODO.

Verified: g_gauge2=0 is BYTE-IDENTICAL (gauss 7.235e-14, E_em unchanged); on-path
runs, gauge-1 Gauss law stays at 7.2e-14 (the 2nd gauge doesn't disturb the 1st).

**RESULT — the 2nd gauge ENHANCES the 4-node** (kernel, η=0, T8 charges):
    g2     s_max(t16)   drift
    0.0    0.084        -1.20%
    0.2    0.087        -0.63%   <- dispersal HALVED, denser
    0.5    0.099        +4.1%
    1.0    0.106        +8.5%
At moderate g2 the cluster holds denser and disperses half as fast — the relative
gauge binds the charge-asymmetric 4-node. At strong g2 the energy GROWS (positive
drift): E2 starts at 0 so the 2nd Gauss law div E2=g2·ρ2 is VIOLATED at init and
pumps spurious energy.

### Gauss-2 projection added (2026-06-26) — strong coupling now clean
Added `gauss2_residual_at` + `init_gauss2_project` (CG Poisson, mirror of gauge1),
called at init. Fixed the energy injection: with projection, g2=0.5 → drift
−0.57% (was +4%), g2=1.0 → +0.17% (was +8.5%). Strong g2 compresses harder
(g2=2.0: s4=0.147 peak). Residual injection only at g2≳2 (missing η-seagull /
strong-field). The 2nd gauge now cleanly enhances + compresses the 4-node.

### Path to a STABLE particle (goal, 2026-06-26)
Key lesson: Gaussian blobs (gen_blob_field, incl. the CMA 4-node) are NOT
stationary — they all breathe, collapse, or disperse, at any charge (amplifying
the 4-node x3 → violent collapse s_max 15-20 → −82% drift). Stability requires the
true `radial_qball` soliton profile. The 2nd gauge couples only to ASYMMETRIC
(relative) charge, so the right stable, gauge-coupled object is a **flavored ball**
(radial profile + frequency partition via `gen_qball_flavored`).
- **STABLE:** mild flavored ball (ω=1.47,1.47,1.42): s_max FLAT at ~0.048 over
  T=40, charge conserved (1.181→1.178), drift −0.54%. A stable particle.
- gauge2 effect on the mild ball is negligible (small asymmetry).

### STABLE PARTICLE ACHIEVED (goal met, 2026-06-26)
Over **T=80** (N=48, gauged), three flavored balls all stay stable:
    config                s_max(t10->t80)   Q(118->)   drift    gauss
    mild  (g2=0)          0.0484 -> 0.0483   ->117.0    -1.09%   6.4e-14
    strong-flavor (g2=0)  0.0490 -> 0.0487   ->117.8    -1.09%   6.5e-14
    strong-flavor (g2=1)  0.0491 -> 0.0490   ->117.0    -1.07%   6.6e-14
s_max is FLAT (core persists, not decaying), charge conserved to <1%, the gauge-1
Gauss law sits at the integrator floor, and (sf_g1) the SECOND gauge is on and does
not destabilize it. The −1% energy drift is the absorbing boundary removing the
flavored clock tails, not core decay.

Conclusion for the stability goal: the stable particle is a proper `radial_qball`
soliton (symmetric or flavored), NOT a Gaussian-blob cluster. The 4-node geometry
is fundamentally non-stationary (breathes/collapses at any charge) and cannot be
stabilized by gauge or compression — those enhance it transiently but it is not a
soliton. The compact-U(1) second gauge is fully built, verified (byte-identical
off, both Gauss laws at floor, enhances asymmetric configs), and coexists stably
with the flavored ball. Deep-compression of blobs was tried (relaxer + amplify +
strong g2) - all transient; the lasting stability comes from the true soliton
profile.

ANCHOR: the symmetric Q-ball (the established stable particle) over the SAME T=80
shows the identical -1.08% drift with s_max flat (0.0482 to 0.0479) and Q 118 to
117 - proving the ~1% is the small-box (L=12) absorbing-BC clock-tail effect, NOT
decay, and is identical for the known-stable ball. STABILITY CONFIRMED: flat s_max
+ conserved charge + Gauss law at the floor = a stable particle, for the symmetric
ball, the flavored ball, and the flavored ball WITH the second gauge on.

(historical) Remaining for STABILITY: the raw 4-node still BREATHES
(collapse↔rebound) and slowly drifts — it is not at equilibrium. Plan:
(a) relax/settle the 4-node+gauge2 to its bound equilibrium (long-T absorbing or
relaxer), (b) optimize gauge2 charges to the 4-node asymmetry (107:285:160), (c)
deep compress (relaxer/strong-g2) then release, (d) η-seagull + E_em2 for full
conservation, (e) CUDA port.
