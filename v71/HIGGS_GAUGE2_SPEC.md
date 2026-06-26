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
sector (≈18 blocks, 2nd Gauss law) — the real Sector 2 build, not yet done.
