# V24-PT Results: Pairwise + Triple Product, 120-degree Oscillon

## Answer

**No.** There is no stable 120-degree oscillon with pairwise + triple product coupling
at the proposed parameters (mu=-20, kappa=20, m=1.0). The system collapses to a
**static condensate** for all lambda in [0.80, 0.95] at amplitude A=0.8. The 3omega < m
radiation threshold is never tested because the system never oscillates.

## What Happens Instead

For all lambda in the scanned range, the initial 120-degree Gaussian immediately
(within t ~ 50-100) falls into a deeply bound **static** state:

| lambda | m_A    | E_settled | E_final  | dE/dt (late) | pk    |
|--------|--------|-----------|----------|-------------|-------|
| 0.80   | 0.447  | -16.10    | -16.27   | -1.2e-05    | 1.19  |
| 0.82   | 0.424  | -18.06    | -18.26   | -1.7e-05    | 1.21  |
| 0.84   | 0.400  | -20.14    | -20.36   | -1.6e-05    | 1.21  |
| 0.86   | 0.374  | -22.28    | -22.53   | -1.5e-05    | 1.22  |
| 0.87   | 0.361  | -23.27    | -23.69   | -2.1e-05    | 1.26  |
| 0.88   | 0.346  | -24.28    | -24.70   | -2.2e-05    | 1.33  |
| 0.89   | 0.332  | -25.29    | -25.60   | -2.6e-05    | 1.31  |
| 0.90   | 0.316  | -26.46    | -26.77   | -2.9e-05    | 1.32  |
| 0.92   | 0.283  | -28.19    | -28.57   | -3.9e-05    | 1.47  |
| 0.95   | 0.224  | -32.21    | -32.82   | -5.2e-05    | 2.02  |

The "peak" values (1.2-2.0) are STATIC field amplitudes, not oscillation amplitudes.
The DFT measures omega ~ 0.01-0.02, which is DC (no oscillation). The energy slowly
drifts more negative (E_retained ~ 1.01-1.02) as absorbing boundaries drain small
fluctuations around the static state.

## Why: Static Condensation

The combined potential V = (m^2/2)|phi|^2 + lambda*sum(phi_a*phi_b) + (mu/2)P^2/(1+kappa*P^2)
has a **deep static minimum** away from phi=0 when mu < 0.

For the 120-degree configuration (phi_a = A*cos(2*pi*a/3)):
- Pairwise: V_pw = -3*lambda*A^2/4 (attractive)
- Triple: P = A^3/4, V_triple = (mu/2)*(A^6/16)/(1+kappa*A^6/16) (attractive for mu<0)
- Mass: V_mass = 3*m^2*A^2/4 (repulsive)

At A ~ 1.3, the attractive terms dominate, creating a static well with depth ~ -20 to -30.
The system falls into this well within ~50 time units and stays there permanently.

## Amplitude Scan (lambda=0.80)

Tested whether smaller amplitudes avoid collapse:

| A_init | After equil | Behavior |
|--------|-------------|----------|
| 0.1    | pk=0.023    | Dispersing wave, no binding |
| 0.2    | pk=0.046    | Dispersing wave, no binding |
| 0.3    | pk=0.072    | Dispersing wave, E drops 40% in t=2000 |
| 0.4    | pk=0.105    | Dispersing wave, E drops 50% in t=2000 |
| 0.5    | pk=0.145    | Dispersing, near collapse threshold |
| 0.8    | pk=1.22     | Collapses to static condensate |

**There is no intermediate amplitude** that produces a stable oscillating oscillon.
Below the collapse threshold, the triple product coupling is too weak (it scales as A^6)
to bind the 120-degree mode. Above it, the system collapses to a static state.

## Why 120-degree Binding is 12x Weaker Than Symmetric

For symmetric phasing (all in phase):
  P = f^3 * cos^3(wt),  <P^2> = (3/8)*f^6

For 120-degree phasing:
  P = (f^3/4) * cos(3wt),  <P^2> = (1/32)*f^6

The time-averaged triple product potential is **12x weaker** for 120-degree than
symmetric at the same amplitude. This means the 120-degree mode needs much larger
amplitude to achieve the same binding, but large amplitude causes static collapse.

## Parameters

- mu = -20.0, kappa = 20.0, m = 1.0
- A = 0.8, sigma = adapted to m_A
- Nx = 2000, xmax = 60.0, dx = 0.060
- t_equil = 1000, t_run = 5000
- Absorbing boundary in outer 25%

## Files

- `src/pt120.c` — scanner code
- `data/pt120_summary.tsv` — tabulated results
- `data/pt120_lam*.tsv` — time series per lambda
- `data/pt120_lam*_spectrum.tsv` — DFT spectra

## Conclusion

The proposal's prediction that 3*omega < m would produce a perfectly stable 120-degree
oscillon is **not realized**. The problem is not the radiation threshold but the
**absence of a bound oscillating state** in the first place. The triple product
coupling that should bind the 120-degree mode creates a static condensate instead
of an oscillating oscillon.

For stable 120-degree oscillons, one would need:
1. A binding mechanism that is NOT the triple product (which is too weak at 120-deg)
2. Or a modified potential that prevents static condensation while allowing oscillation
3. Or a fundamentally different approach (e.g., Q-ball-type conserved charge)
