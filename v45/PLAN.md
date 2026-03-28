# V45 — Differential Mass Defect: Deuterium Binding Energy

## Goal

Measure the nuclear binding energy of deuterium by comparing identical
UUD+UDD systems at different inter-baryon separations. All runs use the
same seeds (gen_deuterium.c), same grid, same BC — the ONLY variable is
the separation distance D. The mass defect is the difference in E_total
between bound and unbound configurations.

## Simulation Matrix

| Run | D (separation) | Physics expected |
|-----|---------------|-----------------|
| D0  | 0  | Maximum overlap. Phase confinement (P→0) should cause violent repulsion. |
| D5  | 5  | Deep inside confinement zone. Strong repulsion expected. |
| D10 | 10 | Near the braid core boundary. Repulsion → attraction transition? |
| D15 | 15 | Interaction surface (~r=4-6 per baryon). Force balance region. |
| D20 | 20 | Nuclear range. Strong interaction expected. |
| D40 | 40 | V42 equilibrium distance. Should be the energy minimum. |
| D60 | 60 | Weak interaction. Slow drift inward expected. |
| D80 | 80 | Non-interacting baseline. Reference energy. |

## Grid Parameters (same for ALL runs)

- N=512, L=100 (domain [-100, 100]³, dx=0.39)
- Absorbing BC: damp_width=10, damp_rate=0.005
- Physics: m=1.5, m_theta=0, eta=0.5, mu=-41.345, kappa=50
- T=500 (allows equilibration and drift measurement)
- Seeds: gen_deuterium -N 512 -L 100 -sep D
- Output: COLZSTD codec, f32, snap_dt=100, diag_dt=1.0

## Measurements

For each run, from diag.tsv:
1. **E_total(t)** — time-averaged over t=200-400 (post-transient)
2. **E_pot(t)** — binding potential energy
3. **P_int(t)** — integrated triple product

From SFA snapshots:
4. **Centroid separation D(t)** — track via |P|-weighted centroid per baryon
5. **Per-particle E_pot** within r<15 of each centroid
6. **Cluster count** — fragmentation analysis

## Mass Defect Calculation

    E_bind = <E_total>(D=40) - <E_total>(D=80)

If E_bind < 0, deuterium is genuinely bound. The binding energy per nucleon
is E_bind / 2.

## Expected Results

- D=0,5: E_total >> baseline (repulsion energy), baryons fly apart
- D=10-15: E_total > baseline but decreasing, baryons may settle or separate
- D=20: E_total ≈ baseline or slightly below (nuclear interaction)
- D=40: E_total minimum (bound state, if binding exists)
- D=60: E_total slightly above D=40 (weak attraction, baryons drift inward)
- D=80: E_total = baseline (no interaction)

The force curve F(D) can be extracted from the centroid drift rates.

## GPU Requirements

- 8 runs × N=512 × T=500 × ~170ms/step ≈ 8 × 2.5 hr = 20 hr total
- V100-32GB or RTX PRO 4500 (32 GB)
- With COLZSTD: ~2 GB per run SFA → 16 GB total
- Disk: 30 GB sufficient (run sequentially, download between runs)
- Estimated cost: ~$5-8 depending on GPU

## Execution Plan

1. Generate all 8 seeds locally (gen_deuterium with -sep 0,5,10,15,20,40,60,80)
2. Provision GPU instance
3. Run sequentially (or 2 at a time if disk allows)
4. Auto-download diag + SFA for each run
5. Analyze after all runs complete
6. Extract force curve F(D) and binding energy
