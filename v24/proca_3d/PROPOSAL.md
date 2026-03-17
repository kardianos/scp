# V24-S4: 3D Verification of Proca Mediator

## Thesis

All Proca results (ME, P1-P5) are in 1D. In 3D:
- The Proca mediator decays as e^{-m_A r}/r (proper Yukawa)
- The oscillon is spherical (different profile from 1D)
- The v21 3D oscillon uses μ=-20, κ=20, m=1.0 (same as 1D tests)
- The 3D oscillon has ω≈0.95 (5% gap margin, tighter than 1D's 13%)

Key questions:
1. Does the 3D oscillon survive with pairwise coupling?
2. If yes, what is the maximum λ before instability?
3. Can we measure the Proca-mediated force in 3D?

## Method

### Phase 1: Single 3D oscillon with pairwise coupling

1. Use the v21 3D solver (triad3d.c) as base
2. Add pairwise coupling: -λ(φ_b+φ_c) to each field's acceleration
3. Start with λ = 0.5 (safe, well below m²=1.0)
4. N=128, L=20, t=1000 (same as v21 production but shorter)
5. Monitor: fc, ω, E, peak amplitude
6. If stable: increase to λ = 0.9, 0.95, 0.99

### Phase 2: Scan λ to find maximum

7. Scan λ in {0.5, 0.7, 0.8, 0.9, 0.95, 0.99}
8. At each λ: evolve t=500 (fast check)
9. Record: survived (fc>0.5), ω, E
10. Find λ_max where 3D oscillon still exists

### Phase 3: Two-oscillon 3D interaction (if Phase 2 positive)

11. At a promising λ: equilibrate single 3D oscillon (t=500)
12. Place two along z-axis at separation D
13. Track separation vs time
14. The Proca force in 3D: F = F₀·e^{-m_A D}/D (Yukawa)

Note: 3D is expensive (N=128³ ≈ 2M points × 9 arrays). Each run takes
~30 min. Keep Phase 3 minimal: one λ, two D values.

## Important Physics Note

The 3D oscillon has tighter gap margin (5% vs 13% in 1D). The pairwise
coupling INCREASES the symmetric mass (m²_S = m²+2λ), which should HELP
the 3D oscillon (deeper below gap). But it also changes the profile shape,
which could affect the 3D stability differently than 1D.

## Reference Code

- v21 3D solver: `/home/d/code/scp/v21/src/triad3d.c` (MAIN base code)
- v24/maxwell_e: `/home/d/code/scp/v24/maxwell_e/src/maxwell_e.c` (pairwise 1D)

## Output

- `src/proca_3d.c` — 3D solver with pairwise coupling
- `data/proca3d_lam{L}_ts.tsv` — time series per λ
- `RESULTS.md`

## Parameters

μ=-20, κ=20, m=1.0, A=0.8, σ=3.0
λ scan: {0.5, 0.7, 0.8, 0.9, 0.95, 0.99}
Grid: N=128, L=20 (Phase 1-2), N=128, L=30 (Phase 3)
CFL: dt = 0.25·dx

Compile: `gcc -O3 -fopenmp -Wall -o proca_3d src/proca_3d.c -lm`
