# V44 — Isolated Baryon Energetics & Neutron Degradation

## Goals

1. **F24: Controlled Mass Defect** — Run isolated UUD (proton) and UDD (neutron)
   at the SAME grid (N=512, L=100, T=500) as V42 deuterium. Compare time-averaged
   energetics to get the first real binding energy data point.

2. **F18: Neutron Degradation** — Run UDD to T=1000 (double the deuterium run).
   Track phase convergence, confinement loss, and structural decay quantitatively.
   Confirm the free neutron is unstable on its own.

3. **Charge Analysis** — Measure far-field θ structure of isolated UUD and UDD.
   Confirm charge signatures (+1 vs 0/−1) on isolated baryons at production grid.

4. **OQ1: Multipole Decomposition** — Post-process existing V41 UUD data for
   angular power spectrum of θ at r=20. Zero GPU cost.

## Simulation Plan

### Run 1: Isolated UUD Proton (N=512, L=100, T=500)
- Seed: gen_proton_analytical -N 512 -L 100 -chirality UUD
- BC: absorbing sphere (damp_width=10, damp_rate=0.005) — matches V42
- Physics: m=1.5, m_theta=0, eta=0.5, mu=-41.345, kappa=50
- Output: snap_dt=100 (6 frames), diag_dt=1.0
- GPU: V100-32GB, ~2.5 hr
- Purpose: mass defect baseline, charge measurement

### Run 2: Isolated UDD Neutron (N=512, L=100, T=1000)
- Seed: gen_proton_analytical -N 512 -L 100 -chirality UDD
- BC: absorbing sphere (damp_width=10, damp_rate=0.005) — matches V42
- Physics: same as Run 1
- Output: snap_dt=100 (11 frames), diag_dt=1.0
- GPU: V100-32GB, ~5 hr
- Purpose: mass defect baseline (T≤500), degradation study (T=500→1000),
  charge measurement, template extraction if stable

### Key Comparisons (F24)
- V42 deuterium: <E_pot> = -94.7 (time-averaged)
- V44 UUD + V44 UDD: compare <E_pot> at same conditions
- If E_deut < E_UUD + E_UDD → genuine binding energy = mass defect

### Neutron Degradation Metrics (F18)
- Phase diversity: track carrier phase spread over time (should converge for UDD)
- P_int: triple product integral (should decay)
- θ confinement: outer/inner ratio (should increase → radiation escape)
- Cluster count: connected component analysis (should merge from 3 → 1)
- Compare to UUD at same time points as control

## Data Storage
- SFA files → /space/scp/v44/
- diag.tsv, analysis → /home/d/code/scp/v44/
