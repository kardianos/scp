# T3: Late-Time Radiation Rate — Results

## Setup
- Grid: N=192, L=60, dx=0.628, dt=0.126
- Damping: r_perp > 54 only (0.90*L)
- Dynamics: BIMODAL params, mass^2=2.25 (from phys[14]=1.50)
- Runtime: T=2000 (15917 steps), wall time 2854s (~48 min)

## Key Result: STEADY RADIATION — Braid is NOT oscillon-like

The braid radiates continuously at a roughly constant rate. It is **not** a
long-lived quasi-stable state like an oscillon.

### Energy evolution
| Time | E_core (R<8) | E_shell (8-30) | E_total | fc (core fraction) |
|------|-------------|----------------|---------|-----|
| 0    | 6513        | 39             | 6552    | 0.995 |
| 500  | 3537        | 497            | 4580    | 0.797 |
| 1000 | 3315        | 301            | 3953    | 0.827 |
| 1500 | 3114        | 205            | 3572    | 0.893 |
| 2000 | 2783        | 226            | 3228    | 0.925 |

### Radiation rate (T > 500)
- Mean |dE_core/dt| = 0.67 per unit time
- Linear trend slope = -7.6e-5 (essentially flat — no decrease)
- Power-law fit: |dE/dt| ~ t^{+0.84} (positive exponent = rate NOT decreasing)
- 83% of dE/dt measurements are negative (net outward energy flow)

### Energy budget
- **E_core lost 57%** of initial value (6513 -> 2783)
- **E_total lost 51%** (6552 -> 3228) — half the energy absorbed by damping layer
- E_shell stays ~200-300 throughout (steady-state radiation zone)

### Topology
- Winding oscillates between +/-1, 0, occasionally +/-2 or +/-3
- Final winding = 0 (no stable topological charge)
- The winding is NOT conserved — this is consistent with the braid being unstable

## Interpretation

1. **Not an oscillon**: An oscillon would show dE/dt ~ -C/t^alpha with alpha > 1,
   meaning radiation rate decreases over time. Here the radiation rate is roughly
   constant or even slightly increasing, giving power-law exponent +0.84.

2. **Two-phase decay**:
   - T=0-150: Rapid initial shedding (E_total drops 6552->5827, 11% in 150 time units)
   - T=150-2000: Slower steady radiation (~1.7 per time unit in E_total)

3. **Core concentrating but shrinking**: fc increases from 0.78 to 0.92 after T=500,
   meaning the braid contracts as it radiates. But E_core still drops steadily.

4. **Lifetime estimate**: At current rate (~1.7/t in E_total), the braid would fully
   dissipate in another ~1900 time units (total lifetime ~4000). This is NOT the
   10^4-10^5 lifetimes seen in true oscillons.

## Conclusion

The bimodal braid with mass^2=2.25 is a **transient** structure, not a quasi-stable
soliton. It radiates energy at a roughly constant rate and will fully dissipate on
a timescale of ~4000 time units. The winding number is not conserved. This rules out
oscillon-like longevity for this parameter set.

## Files
- `src/t3.c` — simulation code
- `data/t3_results.tsv` — 81-row time series (T, E_core, E_shell, E_total, fc, peakP, winding)
