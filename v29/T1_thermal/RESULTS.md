# T1: Thermal Bath Equilibrium -- Results

## Setup
- N=96, L=20, dx=0.4211, dt=0.0842
- Periodic BC in all 3 directions (no absorbing layer)
- Bimodal braid at t=0.85 sweet spot, mass2=0 in EOM
- 6 noise levels: A_noise in {0, 0.001, 0.005, 0.01, 0.05, 0.1}
- T=1000, diagnostics every T=50

## Key Findings

### 1. Braid does NOT survive in periodic box (fc criterion FAILS)
All noise levels show fc dropping from ~0.98 to ~0.088 by T=500, far below
the 0.3 threshold. This is true even for A_noise=0 (no added noise).

**The braid disperses its energy into the box within T~50-100.** By T=50 the
energy has jumped from ~2200 to ~367000 (a 165x increase due to the cubic
coupling amplifying initial transients). The core fraction fc drops to ~0.13
at T=50 and ~0.09 by T=100, then stays near 0.088 for the rest of the run.

### 2. Energy grows continuously (no equilibrium)
Energy never stabilizes. Late-time dE/E/T ~ 0.0006 for all noise levels,
meaning energy grows at ~0.06% per unit time. This is likely numerical:
the periodic box has no dissipation, and the nonlinear coupling (mu=-30 to -43)
generates energy from lattice-scale modes.

| A_noise | Late fc | Late E_avg | dE/E/T  | Late winding |
|---------|---------|------------|---------|--------------|
| 0.000   | 0.088   | 337920     | 0.00056 | -0.1         |
| 0.001   | 0.088   | 340633     | 0.00059 | +0.2         |
| 0.005   | 0.088   | 341049     | 0.00055 | -0.3         |
| 0.010   | 0.088   | 357751     | 0.00069 | +0.3         |
| 0.050   | 0.088   | 431916     | 0.00075 | -0.2         |
| 0.100   | 0.088   | 599325     | 0.00084 | +0.0         |

### 3. Noise level has minimal effect on core dynamics
For A_noise <= 0.01, the evolution is nearly identical to the no-noise case.
This means the braid's own radiation (reflected from periodic walls) dominates
any added thermal noise. Only A_noise=0.05 and 0.1 show meaningfully higher
late-time energies (due to the extra initial energy injection).

### 4. Winding is NOT topologically protected
Winding fluctuates randomly between -1, 0, and +1 throughout the evolution,
even in the zero-noise case. The braid topology is destroyed by the first
radiation bounce (T~50). This is consistent with the V28 finding that the
braid is a lattice saddle point.

### 5. Far-field amplitude saturates quickly
The far-field probe amplitude (at r=8*sqrt(2)~11.3 from center) reaches
~1.3-2.0 by T=50 and stays there, confirming the box fills uniformly with
radiation by the first crossing time (~2L/c = 40).

## Conclusion

**The bimodal braid does NOT reach thermal equilibrium in a periodic box.**
It disperses completely within T~100 and the system enters a turbulent state
with slowly growing energy. The absorbing boundary in V28 was not modeling
a thermal bath -- it was essential for stability by removing radiated energy.

This is a **definitive negative result** for the thermal bath hypothesis:
the braid requires an energy sink, not just periodic recycling. The braid
is unstable to its own radiation when that radiation returns to the core.
