# V24-S2: Mediator Field Characterization — Results

## Summary

**The antisymmetric (Proca) mode between two oscillons has fundamentally
different character at lambda=0.99 vs lambda=0.** At lambda=0.99, the A1 mode
oscillates at omega=0.15 (sub-gap for the antisymmetric mass m_A=0.1), forming
a quasi-static, long-range mediating field. At lambda=0, the A1 mode oscillates
at omega=1.008 (at the symmetric mass gap m=1), behaving as a freely propagating
massive wave that radiates away.

The symmetric mode S dominates the inter-oscillon field in both cases. The
antisymmetric mode A1 exists only when explicitly seeded (the symmetric
initialization phi_1=phi_2=phi_3 is an exact invariant of the dynamics).

**Key negative finding**: A2 remains negligible (~10^-6) in all cases —
the phi_1 <-> phi_2 exchange symmetry breaking does not generate
the (phi_1+phi_2-2*phi_3)/sqrt(6) mode.

---

## Setup

Parameters: mu=-20, kappa=20, m=1.0, A=0.8, sigma=3.0

| Parameter | lambda=0.99 | lambda=0 |
|-----------|-------------|----------|
| m_A (antisymmetric mass) | 0.100 | 1.000 |
| Proca range 1/m_A | 10.0 | 1.0 |
| Oscillon separation D | 30 | 30 |
| D / range | 3.0 | 30 |

Grid: Nx=8000, xmax=200, dx=0.05, dt~0.025
Equilibration: t=5000; Two-body evolution: t=2000

**Symmetry-breaking seed**: eps=0.01 Gaussian at midpoint, added as
+eps to phi_1, -eps to phi_2. This seeds A1 = sqrt(2)*eps*Gaussian
without perturbing S or A2.

## Equilibrated Oscillons

| lambda | omega | peak amplitude | Energy |
|--------|-------|----------------|--------|
| 0.99 | 1.602 | 0.754 | 12.32 |
| 0.00 | 0.870 | 0.503 | 1.27 |

The pairwise coupling dramatically changes the oscillon: higher frequency
(1.60 vs 0.87), larger amplitude (0.75 vs 0.50), and ~10x more energy.
At lambda=0.99 the breathing frequency is well above the symmetric mass gap
(omega=1.60 > m=1.0) but the oscillon still survives because the effective
nonlinear frequency lies below the gap for the bound component.

## Mode Decomposition

Definitions:
- S  = (phi_1 + phi_2 + phi_3) / sqrt(3)     — symmetric (mass m)
- A1 = (phi_1 - phi_2) / sqrt(2)             — antisymmetric (mass m_A)
- A2 = (phi_1 + phi_2 - 2*phi_3) / sqrt(6)  — antisymmetric (mass m_A)

### Time-averaged RMS between oscillons (second half, t > 1000)

| Mode | lambda=0.99 | lambda=0 | Ratio |
|------|-------------|----------|-------|
| S    | 4.86e-2     | 4.65e-2  | 1.05  |
| A1   | 4.45e-3     | 4.71e-3  | 0.94  |
| A2   | 1.08e-6     | 4.70e-7  | 2.3   |

The RMS amplitudes are comparable because both cases received the same
seed perturbation. The difference is NOT in amplitude — it is in frequency
and persistence.

### Frequency Content at Midpoint (DFT, second half)

| Mode | lambda=0.99 | lambda=0 | Interpretation |
|------|-------------|----------|----------------|
| S peak | 1.734 | 0.876 | Oscillon breathing (2*omega/2) |
| A1 peak | **0.150** | **1.008** | Sub-gap mediator vs mass-gap radiation |
| A2 peak | 1.602 | 1.134 | Negligible power |

**This is the central result.** The A1 frequency spectrum shows:

- **lambda=0.99**: Peak at omega=0.15, with secondary peaks near 0.35.
  This frequency is **above m_A=0.1** (so the mode can propagate in the
  antisymmetric channel) but **far below m=1.0** (so it cannot radiate
  into the symmetric vacuum). The Proca field mediates at a frequency
  set by the oscillon interaction, not by the mass gap.

- **lambda=0**: Peak at omega=1.008, right at the mass gap m=1.0.
  Without pairwise coupling, there is no light antisymmetric channel.
  The A1 perturbation simply propagates as a massive Klein-Gordon wave
  at its minimum allowed frequency.

### Spatial Structure (t=1000 profile)

The A1 field between the oscillons (-15 < x < 15):

- **lambda=0.99**: Smoothly varying, amplitude ~3.2e-3 to 3.8e-3.
  Nearly uniform across the 30-unit gap. The field connects the
  two oscillon cores with a flat "bridge."

- **lambda=0**: Bell-shaped, amplitude ~4.5e-3 at center, rising to
  ~5.9e-3 near the oscillon edges. This is a standing wave pattern
  from two massive field tails overlapping.

Both have similar amplitudes because both received the same seed. The
qualitative difference is in the frequency (quasi-static vs oscillating)
and the spatial profile (flat bridge vs standing wave).

## Interpretation

### The Proca Mediator Channel

At lambda=0.99, the pairwise coupling creates an antisymmetric channel with
mass m_A = sqrt(m^2 - lambda) = 0.1, giving range 1/m_A = 10. A perturbation
in this channel:

1. Propagates with dispersion omega^2 = k^2 + m_A^2 (light mass)
2. At the oscillon breathing frequency omega~1.6, it radiates at
   wavenumber k = sqrt(omega^2 - m_A^2) ~ 1.6 (fast-traveling wave)
3. The low-frequency component (omega~0.15) is near-static and
   persists between the oscillons

At lambda=0, the antisymmetric mass equals the symmetric mass (m_A = m = 1).
There is no light channel. The same perturbation can only propagate at
omega >= 1.0, and it radiates away as massive Klein-Gordon waves.

### Why A1 = A2 = 0 Without Seeding

The symmetric oscillon (phi_1 = phi_2 = phi_3) is an exact fixed point
of the S3 permutation symmetry. The equations of motion preserve this:
- Triple product force: symmetric under permutations
- Pairwise coupling: -lambda*(phi_b + phi_c) treats all pairs identically

To excite the antisymmetric channel, a symmetry-breaking perturbation is
required. In the two-oscillon system, this would arise from:
- Asymmetric initialization (different amplitudes)
- Finite numerical noise (but too small to observe at double precision)
- Quantum fluctuations (not modeled)

### Implications for Inter-Oscillon Force

The Proca mediator exists as a channel but is NOT spontaneously excited by
symmetric oscillon tails overlapping. The force between identical oscillons
comes from the symmetric channel (S mode) through tail overlap, not from
the antisymmetric channel.

For the antisymmetric channel to mediate a force, there must be a source
of symmetry breaking — such as different oscillon "flavors" (phi_1 != phi_2
!= phi_3) or a background field that breaks the permutation symmetry.

## Data Files

- `data/field_decomp_lam099_t{0000,0100,0500,1000,2000}.tsv` — spatial profiles
- `data/field_decomp_lam000_t{0000,0100,0500,1000,2000}.tsv` — control profiles
- `data/midpoint_spectrum_lam099.tsv` — DFT at midpoint (S, A1, A2)
- `data/midpoint_spectrum_lam000.tsv` — DFT at midpoint (control)
- `data/twobody_lam099_ts.tsv` — time series diagnostics
- `data/twobody_lam000_ts.tsv` — time series diagnostics (control)
- `data/summary.tsv` — summary table

## Code

`src/proca_field.c` — compile with `gcc -O3 -Wall -o proca_field src/proca_field.c -lm`
