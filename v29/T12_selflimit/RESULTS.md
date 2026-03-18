# T12 Step 0: Self-Limiting Mechanism Screening — Results

## Setup

- Grid: N=96, L=30, T=300
- Background: A_bg=0.1, k_bg=pi/L, traveling wave with phase offsets 2*pi*a/3
- mass^2 = 2.25 (bimodal sweet spot)
- rho0_bg = 3 * A_bg^2 * omega_bg^2 = 0.0678
- Braid: bimodal interpolation at t=0.85
- Edge damping: r > 0.85*L (except M5 which uses full periodic BC)

## Summary Table

| # | Mechanism | fc_final | drho(r=15) | d(drho)/dt | Stabilized? | Braid OK? |
|---|-----------|----------|------------|------------|-------------|-----------|
| 1 | rho-dependent sink | 0.353 | -5.3e-2 | +1.4e-5 | YES | YES (degrading) |
| 2 | Back-pressure | 0.650 | -2.1e-2 | -2.5e-6 | YES | YES |
| 3 | Speed-dependent | 0.077 | -1.5e-2 | -3.4e-4 | NO | NO |
| 4 | Saturating potential | 0.678 | -2.2e-2 | +5.3e-5 | YES | YES |
| 5 | Conservation (periodic) | 0.389 | -3.3e-2 | +6.8e-5 | MARGINAL | YES |
| 6 | Nonlinear wave speed | 0.144 | -4.6e-2 | +1.8e-5 | YES | WEAK |
| 7 | Two-component | 0.910 | -4.2e-2 | +2.6e-5 | YES | YES (strong) |
| 8 | Mexican hat | 0.052 | -6.6e-2 | -1.3e-5 | YES (trivially) | NO |

## Detailed Analysis

### Winners: M2 (Back-Pressure) and M4 (Saturating Potential)

Both M2 and M4 show the ideal combination:
- **Braid survives and strengthens** (fc rises from 0.42 to 0.65-0.68)
- **Depletion stabilizes** at a moderate level (~-0.02, about 30% of rho0_bg)
- **Rate of change near zero** at T=300 (|d(drho)/dt| < 1e-4)
- **Energy loss moderate** (E: 11910 -> 6900, about 42% loss from edge damping)

M2 and M4 give nearly identical results, which makes sense: the back-pressure
gradient force (M2) and the saturating potential force (M4) both push field
back toward depleted regions. M4 is slightly cleaner conceptually.

### Honorable Mention: M7 (Two-Component)

M7 preserves the braid structure best (fc=0.91, barely changes from 0.99).
This is because the S and B fields are separate: the braid structure (S) is
protected from the background dynamics. However:
- Depletion is larger (-0.042) and still has a slight drift
- The model is more complex (6 field arrays instead of 3)
- The physical interpretation is less clean

M7 is worth revisiting if two-field models become important.

### Failures

**M3 (Speed-Dependent)**: The most anticipated mechanism FAILS. The braid
dissolves rapidly (fc: 0.42 -> 0.08). The rho-dependent Laplacian coefficient
destabilizes the braid core: where phi^2 is large (braid core), the effective
wave speed shoots up, causing dispersive blowout. The idea is sound but the
implementation needs the Laplacian coefficient based on ENERGY density (which
is large everywhere including far field), not phi^2 which is concentrated at
the core. Even with energy-based rho, the braid-to-background ratio is too
large for stable dynamics.

**M8 (Mexican Hat)**: The Mexican hat potential drives |phi| toward rho0=0.1
everywhere. Since the braid has |phi| >> 0.1 at the core, the potential
actively destroys the braid. By T=300, both braid and background are nearly
gone (E: 13549 -> 1299). The Mexican hat would need rho0 comparable to the
braid amplitude to preserve the braid, but then it would not limit the
background depletion.

**M1 (Density-Dependent Sink)**: Stabilizes the depletion but continuously
drains energy (E: 11910 -> 2034). The sink is non-conservative by design,
so the total energy decreases monotonically. The braid survives but weakens.
Not suitable as a fundamental mechanism since it violates energy conservation.

**M6 (Nonlinear Wave Speed)**: Similar issue to M3 — the phi^2-dependent
wave speed is dominated by the braid core, not the background. The braid
weakens (fc: 0.42 -> 0.14) and energy drains away. The depletion becomes
nearly uniform (no radial structure), which defeats the purpose.

**M5 (Conservation)**: Energy is approximately conserved (11910 -> 11562,
3% loss from the non-periodic x,y boundaries — would need truly periodic
initialization to test properly). The depletion oscillates but does not
clearly stabilize. This is the null hypothesis: the standard PDE with
periodic BC naturally redistributes energy. The oscillations suggest the
depletion is NOT self-limiting in this model — it just sloshes around.

## Key Observations

1. **The braid is fragile**: Any mechanism that modifies the Laplacian or
   wave speed based on local field amplitude tends to destroy the braid.
   The braid core has |phi| ~ 0.8, while the background has |phi| ~ 0.1.
   Any nonlinear speed that scales with |phi| will preferentially affect
   the core.

2. **Restoring forces work**: Both M2 and M4 work by adding a restoring
   force that opposes depletion. The force is gentle enough to preserve the
   braid while limiting the depletion zone.

3. **Two-component (M7) protects the braid**: Separating structural and
   background fields is the most robust way to preserve the braid. The
   coupling is weak enough that the braid structure is barely affected.

4. **Initial transient**: All mechanisms show an initial transient (T=0 to
   T~50) where the braid adjusts from its analytic initialization to its
   numerical equilibrium. The depletion measurements after T=100 are more
   reliable.

## Recommendation

**Proceed with M4 (Saturating Potential)** as the primary mechanism.
Reasons:
- Clean physical interpretation: depletion has an energy cost
- Braid survives and stabilizes
- Depletion is moderate and self-limiting
- Simple implementation (one additional force term)
- Derives from a potential, so energy is well-defined

**Secondary**: M7 (Two-Component) as alternative if we need
braid-background separation for physical reasons.

## Next Steps

1. Run M4 at higher resolution (N=128, L=60) and longer time (T=1000)
   to confirm the depletion profile stabilizes and measure its radial shape
2. Vary lambda_depl to find the optimal balance between depletion depth
   and braid stability
3. Measure the radial profile of depletion: is it 1/r? 1/r^2? exponential?
4. Test with a moving braid: does the depletion trail follow?
5. Two braids: does the depletion gradient produce attraction?
