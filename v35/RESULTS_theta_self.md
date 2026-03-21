# Computation 4 Results: theta Self-Interaction

## Summary

**Definitive null result**: Adding a triple-product self-coupling V_theta to the theta
equation does NOT produce theta-solitons. Instead, the attractive potential triggers
uniform tachyonic condensation -- theta fills the entire box homogeneously, with no
localized structure. The theta-braid concept is ruled out for this coupling form.

## Phase A: Parameter Scan

Scanned mu_theta in {0, -1e3, -1e4, -1e5, -5e5, -1e6} and kap_theta in {50, 5000, 50000}.
All runs: N=64, L=20, T=100, eta=0.5, m_theta^2=0, single phi-braid.

### Results Table (t=100)

| mu_theta | kap_theta | theta_rms | Q_max    | E_pot_theta | drift%    | Status   |
|----------|-----------|-----------|----------|-------------|-----------|----------|
| 0        | any       | 3.94e-2   | 1.19e-3  | 0           | -0.4%     | STABLE   |
| -1000    | 50        | 18.0      | 3.0e+5   | -6.70e+5    | +5.0e6%   | BLOWUP   |
| -1000    | 5000      | 0.564     | 1.85     | -5.91e+3    | -13.6%    | marginal |
| -1000    | 50000     | 0.247     | 0.137    | -3.55e+2    | -0.6%     | STABLE   |
| -10000   | 50        | 780       | 1.6e+11  | -6.71e+6    | +1.9e10%  | BLOWUP   |
| -10000   | 5000      | 8.69      | 1.6e+4   | -6.71e+4    | +2.9e5%   | BLOWUP   |
| -10000   | 50000     | 0.441     | 0.337    | -6.36e+3    | -4.9%     | STABLE*  |
| -100000  | 50        | 6950      | 9.9e+12  | -6.71e+7    | +9.6e11%  | BLOWUP   |
| -100000  | 5000      | 117       | 1.7e+8   | -6.71e+5    | +2.8e8%   | BLOWUP   |
| -100000  | 50000     | 7.94      | 6.6e+4   | -6.71e+4    | +1.2e6%   | BLOWUP   |
| -500000  | 50000     | 156       | 2.2e+8   | -3.35e+5    | +6.2e8%   | BLOWUP   |
| -1000000 | 50000     | 237       | 5.9e+8   | -6.71e+5    | +1.4e9%   | BLOWUP   |

*STABLE = drift < 5%, marginal = drift 5-20%, BLOWUP = drift > 20%.

### Stability Boundary

Stability requires kap_theta > |mu_theta| / 10 (approximately). Below this, the
tachyonic instability wins: the attractive V_theta drives runaway growth, with Q
growing until limited only by numerical overflow.

### Best Stable Cases

1. **mu=-10000, kap=50000**: theta_rms grows from 0 to 0.44, then STABILIZES (equilibrium).
   E_pot_theta saturates at -6360 (theory: mu/(2*kap)*Vol = -6400). Energy drift -4.9%.

2. **mu=-1000, kap=50000**: theta_rms grows from 0 to 0.25 by t=100, still increasing.
   Would likely reach ~0.3-0.4 at longer times. Energy drift only -0.6%.

3. **mu=-1000, kap=5000**: theta_rms = 0.56. Marginal stability (drift -13.6%).

### Key Finding: Uniform Condensation, Not Soliton Formation

Spatial analysis of the stable cases shows:

| Case              | theta^2 at core (r<5) | theta^2 far (r>10) | Ratio |
|-------------------|----------------------|---------------------|-------|
| baseline (mu=0)   | 8.76e-3              | 3.72e-3             | 2.4   |
| mu=-1e3, kap=5e4  | 1.22e-2              | 2.97e-1             | 0.04  |
| mu=-1e4, kap=5e4  | 5.27e-1              | 5.52e-1             | 1.0   |
| mu=-1e3, kap=5e3  | 9.67e-1              | 8.63e-1             | 1.1   |

The theta field is UNIFORM (ratio ~1.0) or even concentrated AWAY from the braid
(ratio 0.04 for mu=-1e3, kap=5e4). There is NO localized theta structure around the
phi-braid. The self-interaction drives a homogeneous condensate.

## Phase B: theta-Braid Initialization

Initialized a helical theta-braid (same pattern as phi-braid) at x=10, with varying
amplitudes. Ran T=200 to check survival.

| Label               | A_theta | mu_t   | kap_t | theta_rms final | Localization | Status   |
|---------------------|---------|--------|-------|-----------------|--------------|----------|
| control (mu=0)      | 0.20    | 0      | -     | 0.068           | 0.69         | DISPERSED |
| mu=-1e3, kap=5e4    | 0.05    | -1000  | 50000 | 0.279           | ~1.0         | CONDENSED |
| mu=-1e3, kap=5e4    | 0.20    | -1000  | 50000 | 0.309           | ~1.0         | CONDENSED |
| mu=-1e4, kap=5e4    | 0.05    | -10000 | 50000 | 0.439           | ~1.0         | CONDENSED |
| mu=-1e4, kap=5e4    | 0.10    | -10000 | 50000 | 0.461           | ~1.0         | CONDENSED |
| mu=-1e4, kap=5e4    | 0.20    | -10000 | 50000 | 0.484           | 1.20         | CONDENSED |
| mu=-1e4, kap=5e4    | 0.50    | -10000 | 50000 | 0.559           | ~1.0         | CONDENSED |
| no-phi (eta=0)      | 0.20    | -10000 | 50000 | 0.581           | 1.06         | CONDENSED |

"Localization" = |theta|^2 at theta init site / |theta|^2 far from both braids.
1.0 = completely uniform (no structure). <1 = anti-concentrated.

### Key Observations

1. **Without self-interaction (control)**: The theta-braid DISPERSES. theta_rms drops
   from 0.018 to 0.068 (the equilibrium curl-coupling level). No survival.

2. **With self-interaction**: All initial amplitudes converge to the SAME theta_rms
   (~0.44 for mu=-1e4). The initial theta structure is irrelevant -- the condensate
   is uniform regardless of initialization.

3. **No-phi case (eta=0)**: Same result. The theta field condenses uniformly even
   without any phi coupling. The self-interaction alone drives homogeneous condensation.

4. **Localization ratio ~1.0 everywhere**: No case shows localized theta structure.
   The explicitly-initialized theta-braid at x=10 does NOT survive as a soliton.

## Physics Interpretation

The triple-product potential V_theta = (mu/2) Q^2 / (1 + kap Q^2) with negative mu
acts as a TACHYONIC mass: the uniform Q=0 state is unstable, and Q grows until
saturated by kappa. The saturation value is:

    E_pot_theta -> mu/(2*kap) * Volume  (uniform condensate)

This matches observations within 1%.

The crucial difference from the phi-braid: the phi-braid survives because it has
both attractive potential AND the mass gap (m^2 = 2.25) that prevents spreading.
The theta field with m_theta^2 = 0 has no barrier to radial dispersion, so the
attractive potential leads to uniform condensation rather than localization.

**Could m_theta^2 > 0 help?** Unlikely. The phi braid needs a SPECIFIC balance of
mass, potential, and kinetic energy to exist. The theta field would need its own
CMA-ES search to find if any stable braid exists, and the Phase A results suggest
the answer is no: the condensate is the global energy minimum.

## Conclusion

theta-braids do NOT form under the triple-product self-coupling V_theta(Q).
The attractive potential creates a uniform tachyonic condensate, not localized
solitons. This mechanism cannot provide m_eff for the Bohr ratio calculation.

To obtain an effective theta mass, one would need:
1. A different coupling form (e.g., gradient-dependent, topological)
2. Explicit theta mass m_theta^2 > 0 plus fine-tuned potential
3. A fundamentally different approach to theta self-interaction

The simple "copy phi's potential to theta" strategy fails because the theta
field lacks the structural features (mass gap, helical topology protection)
that make phi-braids stable.
