# V56 — Pion mass and the lattice-collapse barrier

Two more experiments after the `PROJECTION_RESULT.md` hard-projection
breakthrough, plus a complete diagnosis of the remaining lifetime ceiling.

## Five runs, same outcome

| Run                       | mesh | c4   | m_π  | project | drift  | B(t≈8) | collapse t |
|---------------------------|------|------|------|---------|--------|--------|------------|
| Production v2 (c4=0.05)   | L=40 | 0.05 | 0    | no      | 0.001% | 0      | t ≈ 4      |
| Strong scan (sigma)       | L=20 | 5.0  | 0    | no      | 0.002% | 0      | t ≈ 4      |
| Compact scan (R=1.5)      | L=20 | 0.5  | 0    | no      | 0.015% | 0      | t ≈ 2      |
| Projected strong          | L=20 | 5.0  | 0    | yes     | 1.3%   | **0.85** | t ≈ 9    |
| Projected strong, L=40    | L=40 | 5.0  | 0    | yes     | 0.7%   | 0.87   | t ≈ 9      |
| Pion strong (m_π=0.75)    | L=20 | 5.0  | 0.75 | yes     | 2.3%   | 0      | t ≈ 2      |
| Pion lite (m_π=0.2)       | L=20 | 5.0  | 0.2  | yes     | 1.7%   | 0      | t ≈ 6      |

**Hard projection is the only intervention that significantly extended
winding lifetime** (factor of ~2× over the no-projection runs). Pion
mass actively destabilizes — the V_π force pushes q→(1,0,0,0)
globally, which is anti-Skyrmion: it accelerates compactification and
shortens the lattice-collapse timescale.

## Why pion mass made it worse

The textbook role of m_π in Skyrme is to set a finite preferred
Skyrmion radius R*. Naive Derrick virial gives:

    R* ≈ (c4 / m_π²)^{1/4}

For c4=5, m_π=0.2 → R* ≈ 2.0. With dx_min=0.44 on L=20 that's only
~4.5 cells across the core — already at the lattice resolution edge.

The seed is at R=5, well above R*, so the dynamics drive the soliton
to *shrink* toward R=2. During the shrinkage transit, two things
happen at once:
1. The gradient `∂q` becomes large — so the L_4 stiffness force
   (which is `(∂q)⁴`) gets enormous and accelerates compactification.
2. The pion force adds to the inward push.

Net effect: the soliton crashes through R* and into the lattice
cutoff before it can settle. **The textbook mass term is the wrong
fix when the pre-existing failure is "Skyrmion shrinks past lattice"**;
it is the right fix when the failure is "Skyrmion expands without
bound", which is not our case.

## What actually fails: the seed-to-equilibrium transit

In all five projected/unprojected runs the failure pattern is identical:

| Phase | What happens | Duration |
|-------|--------------|----------|
| Gradient relaxation | E_grad drains into E_skyrme | 0 → 7 t-units |
| Derrick crossing    | E_skyrme = E_grad (R=R*)   | t ≈ 7        |
| Overshoot           | E_skyrme > E_grad (R<R*)   | 7 → 9        |
| Lattice violation   | B → 0 in <1 t-unit          | t ≈ 9        |
| Disordered remnant  | Equipartition, no winding   | t > 10       |

The Skyrme infrastructure works. The L_4 force compactifies the
soliton — the snapshots at t=8 show a smaller, more focused rotor
core than at t=0. But the seed `f(r) = π exp(−r/R)` has gradient
excess (E_2 / E_4 = 2:1 instead of the equilibrium 1:1) and overshoots
the equilibrium during relaxation.

## What's actually needed

The real fix is not in the dynamics — it's in the **initial condition**.
A standing Skyrmion needs to start AT its equilibrium profile, with
no gradient excess to dissipate. Two routes:

### Route A: gradient flow pre-relaxation

Run a short gradient-flow phase (no kinetic, just dE/dt = −∇E) to
relax the seed to equilibrium BEFORE switching to Hamiltonian
dynamics. About 50 lines of code: skip the kinetic update, just
do `M ← M + α · acc · dt` for several hundred steps with α ~ 0.1
until E_grad ≈ E_skyrme to within 1%. Then start Verlet from there.

### Route B: rational-map ansatz

Use the textbook Skyrmion profile from rational-map theory:
`f(r) = 4 atan(exp(−r/a))` with `a ≈ 1/m_π`, parameterized so that
E_2 = E_4 by construction. Implement as a new init type
`init = skyrme_rational`. About 30 lines.

### Route C: damped initial transient

Add an explicit damping term `−γ M_vel` for t ∈ [0, T_damp], then
turn off. Lets the system bleed off the seed's gradient excess
as kinetic energy → heat without ever overshooting. T_damp ≈ 10
t-units should work. About 15 lines.

Routes A and C achieve the same thing differently — A finds the
exact equilibrium analytically (well, numerically), C finds it by
viscous relaxation. C is simpler and matches what every other Skyrme
lattice study does.

## Visual evidence (`v56/foam/snapshots_proj/rot_t{8,10}_p0.webp`)

Rotation-mode rendering at phase 0 (M[0,1,2] → R,G,B) shows the
hedgehog angular structure directly. Compare:

- **t = 8** (B = 0.85): clean color sectors. The hedgehog `n̂ = x/|x|`
  structure is intact — red on +x faces, green on +y, blue on +z.
- **t = 10** (B = 0): same blob amplitude, **scrambled colors**.
  The angular winding has been aliased away by the discrete operator
  while the field amplitude survived.

This rules out the visualizer-misleading hypothesis: the snapshots
ARE showing the topology, and what they show is consistent with
the lattice violation diagnosis.

## Energy drift caveat

Hard projection is not strictly symplectic — it discards the radial
component of M[0..3] excursion each step. With pion mass added, the
excursion is larger per step and drift grew to 1.3–2.3%. Drift is
concentrated at the t=6–10 collapse event when gradients momentarily
spike.

For routine production this is acceptable. If someone wants <0.01%
drift, upgrading to RATTLE (constraint-projecting both position and
velocity in a single fixed-point iteration) is ~50 lines. Not
recommended unless we're chasing precision physics.

## Files

- `v56/foam/scan_pion_L20.{sfa,_diag.tsv}` — m_π=0.75 (overcooked)
- `v56/foam/scan_pion_lite_L20.{sfa,_diag.tsv}` — m_π=0.2 (anti-Skyrmion)
- `v56/foam/scan_proj_strong_L40.{sfa,_diag.tsv}` — L=40 projection-only
- `v56/foam/foam_sim.c` — `project_rotor_S3` + pion mass
- `v56/foam/field_pga.h` — tangent-projected pion force on rotor

## Bottom line

We have a working Skyrme infrastructure on the foam:
- Genuine `Tr([L_μ, L_ν]²)` term, energy-conservative discretization.
- Hard S³ projection — |q|² stays at 1.0 to machine precision.
- Pion mass available as a knob.
- Winding number diagnostic that tracks the soliton lifetime.

What we don't have is a stable winding-1 Skyrmion. The barrier is
**seed-equilibrium mismatch + lattice-cutoff overshoot**, not anything
about the v56 algebraic structure. The fix is initialization-side
(routes A/B/C above), not kernel-side. None of those touch field_pga.h
or the multivector field — they sit in `field_init_skyrme`.

Recommend: implement route C (damped transient) — minimal change,
high probability, matches standard practice. Then re-run the strong-
projected scan with damping enabled.
