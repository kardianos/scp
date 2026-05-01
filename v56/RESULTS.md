# V56 Stage A — first runs on the L=40 mesh

Four runs of the v56 multivector Higgs kernel on the 3.26M-cell foam
mesh inherited from v55. Field: 8 components per cell with Lorentz
metric `g = diag(−,−,−,−,+,+,+,+)`. Dynamics:

```
□ M[k] = − λ (|M|²_bulk − v²) g_{kk} M[k]   with |M|²_bulk = Σ g_kk M[k]²
```

Parameters: λ = 5, v² = 0.25 (vacuum hyperboloid `|M|²_bulk = v²`),
m² = 0 (pure Higgs dynamics). T = 50 (T = 20 for the vacuum sanity).

## Energy conservation summary

| Run | Wall | E_total | drift envelope | final drift |
|-----|------|---------|----------------|-------------|
| vacuum   | 9.6 min  | 0     | 0.000% | +0.000% |
| qball    | 26.5 min | 18.81 | < 0.001% | +0.000% |
| braid_spatial | 24.4 min | 6574 | < 0.005% | -0.001% |
| braid_rotor   | 23.9 min | 5511 | < 0.006% | +0.000% |

All four runs conserve energy to better than 0.01%. The same Verlet
infrastructure that gave bounded drift in v55 (after the gradient-of-2
fix) is preserved here — the Higgs quartic is just a different
self-coupling, not a different operator.

## Run 1 — `vacuum` (sanity baseline)

**Init**: `M[4] = v = 0.5`, all others zero. The chosen vacuum point on
the |M|²_bulk = v² hyperboloid (spacelike component).

**Result**: completely stationary. Every diag column zero, bulk_min =
bulk_max = 0.250 to machine precision throughout T = 20. The kernel
correctly identifies and preserves the Higgs vacuum.

## Run 2 — `qball`

**Init**: vacuum + Gaussian bump in (M[5], M[6]) plane (A = 0.3, R = 4)
with rotation rate ω = 1, intended as a Q-ball candidate.

**Trajectory**:

| t | bulk_max | comment |
|---|----------|---------|
| 0 | 0.339 | initial bump (35% above v²) |
| 4 | 0.252 | **already collapsed** to v² |
| 9–49 | 0.247–0.260 | small-amplitude radiation in periodic box |

**Conclusion**: Q-ball did not survive. A non-topological scalar with
the standard `(|M|²−v²)²` Higgs potential has no conserved charge to
localise the bump, so the energy radiates as small waves over the box.
The v55 oscillon analog (V(P) saturation) is what made localised
breathing modes possible there — replacing it with the Higgs quartic
removes that mechanism.

## Run 3 — `braid_spatial`

**Init**: v55-style braid with carriers on the three spacelike
components (M[5], M[6], M[7]), Higgs vacuum on M[4]. M[4] adjusted
per cell so |M|²_bulk = v² exactly at t = 0.

**Trajectory**:

| t | bulk range | spatial_rms | comment |
|---|------------|-------------|---------|
| 0 | [0.250, 1.596] | 0.503 | initial state — carriers boost \|M\|² massively at center |
| 4 | [0.110, 0.585] | 0.501 | rapid initial relaxation toward v² |
| 9 | [0.078, 0.405] | 0.517 | bulk range narrowing |
| 24 | [0.179, 0.332] | 0.518 | small ripples around v² |
| 49 | [0.209, 0.348] | 0.504 | persistent ~30% excursion in bulk_max |

**Conclusion**: the spatial-carrier braid does **not collapse** to
vacuum the way the qball did. After 50 t-units the bulk_max remains
40% above v² and bulk_min 16% below — the structure is slowly
radiating but not gone. With longer runs we'd expect either complete
diffusion to vacuum (no topology to protect anything in stage A) or
some metastable Goldstone-mode pattern.

The spatial RMS stays remarkably constant at ≈ 0.50 (= sqrt(v² + small
contributions)) throughout — the Higgs constraint pulls the field back
to the vacuum norm even as cell-by-cell distribution evolves.

## Run 4 — `braid_rotor`

**Init**: same braid waves but on the three timelike rotor components
(M[0], M[1], M[2]). Because g_kk = −1 for these, their contribution to
the bulk norm is negative: `|M|²_bulk = M[4]² − M[0]² − M[1]² − M[2]²`.
At t=0 we adjust M[4] = sqrt(v² + sum_carriers²) ≈ 0.530 to keep the
bulk norm at v².

**Trajectory** (showing the qualitatively different Higgs-amplitude
oscillation — *the* most interesting result of the four runs):

| t | bulk range | rotor_rms | spatial_rms | phase |
|---|------------|-----------|-------------|-------|
| 0 | [0.250, 0.250] | 0.143 | 0.520 | seed (constrained) |
| 4 | [0.015, 0.469] | 0.480 | 0.690 | tachyonic growth begins |
| 9 | [-0.055, 0.534] | 0.900 | 1.094 | bulk goes negative |
| 14 | [-0.003, 0.533] | 1.229 | 1.328 | rotor approaching peak |
| 19 | [0.043, 0.497] | 1.449 | 1.527 | near peak |
| **24** | [0.066, 0.483] | **1.538** | **1.618** | **PEAK — both 10× initial** |
| 29 | [0.044, 0.479] | 1.500 | 1.584 | starts decaying |
| 34 | [0.076, 0.440] | 1.337 | 1.428 | decay |
| 39 | [0.118, 0.380] | 1.059 | 1.170 | decay accelerates |
| 44 | [0.159, 0.347] | 0.682 | 0.846 | converging back |
| 49 | [0.177, 0.310] | 0.224 | 0.537 | near vacuum (rotor → 0, spatial → v) |

**Conclusion**: the rotor seed triggers a complete **Higgs amplitude
oscillation** — the field excursion peaks around t=24, then collapses
back to near-vacuum by t=50. This is a *radial* mode in the
8-component multivector space, swinging away from the vacuum
hyperboloid then back. The negative-metric components on rotors give
the Higgs drive an effective sign flip that drives them away from
zero, but the constraint `|M|²_bulk = v²` (energetically preferred)
eventually pulls everything back.

The peak amplitudes (rotor_rms = 1.54, spatial_rms = 1.62) are 3× the
expected scale `v = 0.5`, confirming the field made a large excursion
off the vacuum manifold before being pulled back. Energy conservation
held throughout — this is genuinely a Hamiltonian oscillation, not a
numerical instability.

## Cross-run comparison

The four runs probe four regimes of the same dynamics:

| Run | seed energy | mechanism | outcome |
|-----|------------|-----------|---------|
| vacuum | 0 | pure stationarity | completely static (sanity check) |
| qball | 19 | small bump in spacelike pair | quick collapse to v² + radiation |
| braid_spatial | 6574 | high-amplitude spacelike braid | slow relaxation, residual structure |
| braid_rotor | 5511 | high-amplitude timelike braid | huge oscillation excursion, then back |

What did NOT happen in any of the four runs:
- A persistent localised structure (no Q-ball, no braid trapping)
- A topological soliton (we have no Skyrme term yet)
- A stable non-vacuum phase (Higgs vacuum is the unique attractor)

What DID happen:
- The infrastructure works at production scale (3.26M cells, T=50 in ~25 min wall)
- Energy conservation across vastly different field excursions
- Two qualitatively different relaxation modes:
  - Spacelike carriers → slow radiative loss, persistent envelope
  - Timelike rotors → tachyonic growth → Higgs amplitude oscillation → decay

## Implications for stage B / C

The pure Higgs `(|M|²−v²)²` potential is doing exactly what the math
predicts: any deviation from the v² hyperboloid is restored. Without
a topological charge or a conservation law (like U(1) for a Q-ball
or Z for Skyrme winding), there's nothing for solitons to hold onto.

For stable particle-like states we need to either:

1. **Stage B (Dirac equation)**: replace second-order Klein-Gordon
   with first-order Dirac. The spinor structure of the relativistic
   quaternion would give half-integer spin and Pauli exclusion, plus
   a fermion number conservation that protects soliton charge.
2. **Stage C (Skyrme term)**: keep Klein-Gordon but add a fourth-
   derivative `Tr([L_μ, L_ν]²)` term that obstructs the topological
   charge from dissipating. Skyrmions in this would carry a winding
   number conserved exactly.

The braid_rotor result is the most encouraging: it shows the
multivector field admits dramatic non-vacuum excursions while
preserving total energy. That's the prerequisite for stage B/C
topological charge to take hold of one of those excursions.

## Files

```
/home/d/code/scp/v56/foam/
  v56_vacuum_L40.sfa             46 MB
  v56_qball_L40.sfa             1.4 GB
  v56_braid_spatial_L40.sfa     2.9 GB   (most data — long-lived structure
                                          → many P-frame cells deviate)
  v56_braid_rotor_L40.sfa       2.7 GB   (large excursions → many cells
                                          deviate from temporal model)
  v56_*_diag.tsv                tabular diagnostics per t-unit
```

## Known issue: rendering

`$HOME/bin/volview`'s color pipeline maps:
- R = |column[0] · column[1] · column[2]|
- G = column[0]² + column[1]² + column[2]²
- B = column[3]² + column[4]² + column[5]²

This was designed for v55's φ[3] + θ[3] layout. For v56's 8-component
multivector with rotor (cols 0–3) and spatial (cols 4–7) the mapping
shows only the M[4] = v vacuum (uniform blue cube) and misses the
braid carrier structure entirely. **Visual debugging of v56 is
currently blocked.**

### TODO: derived-quantity post-processor

Build `v56_derived.c` that reads a v56 cell-native SFA (8 cols) and
writes a 3-column voxel SFA suitable for volview:

| Channel | Quantity | Visual interpretation |
|---------|----------|------------------------|
| R | bulk_dev = \|\|M\|²_bulk − v²\| | non-vacuum activity (red where field is off the hyperboloid) |
| G | rotor_amp = sqrt(M[0]²+M[1]²+M[2]²+M[3]²) | rotor (timelike) magnitude |
| B | spatial_amp = sqrt(M[4]²+M[5]²+M[6]²+M[7]²) | spatial (spacelike) magnitude |

Volview's standard RGB pipeline then renders an intelligible image:
solid blue = pure vacuum; red highlights = where the Higgs constraint
is violated; green tinge = where the spinor (rotor) part is active.

Ballpark cost: ~300 LOC, similar to `cellsfa_to_voxel`. Should be
written before stage B / C runs to make their output viewable.

## Next

- Implement `v56_derived` for visualisation (above).
- Render all four runs (vacuum, qball, braid_spatial, braid_rotor) at
  v34 camera with the derived-quantity SFAs, save to
  `/home/d/code/scp/v56/snapshots/`.
- Decide stage B vs stage C path. Lean toward C (Skyrme term added to
  current Klein-Gordon kernel — minimal kernel surgery) over B
  (Dirac equation requires new integrator).
