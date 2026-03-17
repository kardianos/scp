# V24-D Results: Three-Frequency Non-Degenerate Oscillon

## Summary

**Definitive negative result.** Non-degenerate frequency initialization does NOT produce
lower radiation rates than the degenerate case. The triple-product coupling P = phi_1*phi_2*phi_3
acts as a powerful equalizer: all three fields converge to the SAME frequency within ~100-200
time units regardless of initial frequency splitting. The oscillon is a strong attractor in
field-amplitude space, and frequency degeneracy is a dynamical consequence, not an assumption.

Low-frequency initialization (omega << m) does NOT produce subgap oscillons -- the fields
either (a) shed excess energy and converge to the standard omega ~ 0.87-0.89 oscillon, or
(b) collapse into a static false-vacuum state at strong coupling.

---

## Setup

Lagrangian: three massive scalars with saturating triple-product coupling.

    L = sum_a [ (1/2)(dt phi_a)^2 - (1/2)(dx phi_a)^2 - (m^2/2)phi_a^2 ]
      - (mu/2) P^2 / (1 + kappa P^2),   P = phi_1 phi_2 phi_3

Base parameters: mu=-20, kappa=20, m=1.0, A=0.8, sigma=3.0, Nx=4000, xmax=100, t=10000.

Initialization: phi_a(x,0) = A*g(x), v_a(x,0) = -omega_a * A * g(x).
Each field given a different initial velocity to set its oscillation frequency.

Code: `src/threefreq1d.c`

---

## Results Table

### Standard Coupling (mu=-20, kappa=20)

| Label | omega_1 | omega_2 | omega_3 | E(0) | E(10k) | Retained | dE/dt (late) | omega_final |
|-------|---------|---------|---------|------|--------|----------|-------------|-------------|
| control | 0.87 | 0.87 | 0.87 | 7.23 | 1.27 | 17.6% | -5.86e-6 | 0.887 |
| smallsplit | 0.90 | 0.87 | 0.84 | 7.23 | 1.28 | 17.6% | -5.90e-6 | 0.893 |
| medsplit | 0.95 | 0.87 | 0.79 | 7.25 | 1.28 | 17.6% | -5.90e-6 | 0.893 |
| largesplit | 0.95 | 0.80 | 0.65 | 6.71 | 1.27 | 18.9% | -4.30e-6 | 0.881 |
| lowfreq | 0.30 | 0.30 | 0.30 | 3.82 | 1.26 | 33.1% | -1.41e-6 | 0.875 |
| mixedlow | 0.40 | 0.30 | 0.20 | 3.86 | 1.26 | 32.8% | -1.41e-6 | 0.875 |

### Strong Coupling (mu=-60, -100) — Low-Frequency Cases

| Label | omega | mu | E(0) | E(10k) | omega_final | Status |
|-------|-------|----|------|--------|-------------|--------|
| lowfreq_mu60 | 0.30 | -60 | -0.22 | -49.6 | 0.00 | COLLAPSED (static false vacuum) |
| lowfreq_mu100 | 0.30 | -100 | -4.27 | -176.3 | 0.00 | COLLAPSED (static false vacuum) |
| mixedlow_mu60 | 0.40/0.30/0.20 | -60 | -0.19 | -49.6 | 0.00 | COLLAPSED |
| mixedlow_mu100 | 0.40/0.30/0.20 | -100 | -4.23 | -176.3 | 0.00 | COLLAPSED |

---

## Key Findings

### 1. Frequency Degeneracy is Dynamically Enforced

In ALL cases with mu=-20, all three fields converge to the SAME oscillation frequency
within ~100-200 time units. The triple-product coupling P = phi_1 * phi_2 * phi_3 creates
a restoring force that equalizes the three modes.

For the largesplit case (0.95, 0.80, 0.65 initial):
- t=0: inter-field spread ~ 0.25 (fields oscillate visibly out of phase)
- t=50: spread ~ 0.07
- t=100: spread ~ 0.01
- t=200: spread ~ 0.002
- t>200: spread < 0.003 (degenerate)

This is the same mode-equalization mechanism observed in v21 Test 2 (asymmetric amplitudes).
The coupling enforces both amplitude and frequency equalization.

### 2. No Radiation Improvement from Frequency Splitting

Comparing late-time energy loss rates:

| Case | dE/dt (late) | vs control |
|------|-------------|------------|
| control (degenerate) | -5.86e-6 | baseline |
| smallsplit (delta=0.03) | -5.90e-6 | +0.6% (worse) |
| medsplit (delta=0.08) | -5.90e-6 | +0.7% (worse) |
| largesplit (asymmetric) | -4.30e-6 | -27% (better) |
| lowfreq (0.30 degen) | -1.41e-6 | -76% (better) |
| mixedlow (0.40/0.30/0.20) | -1.41e-6 | -76% (better) |

The small/medium splits show NO improvement -- the frequencies converge too fast for
any radiation reduction to occur.

The largesplit and lowfreq cases show lower dE/dt, but this is because they start with
LESS initial energy (lower initial velocities), not because of any frequency-splitting effect.
All cases converge to the SAME final state: E_final ~ 1.26-1.28, omega ~ 0.87-0.89.
The "better" dE/dt simply reflects being closer to the equilibrium already.

### 3. Final State is Universal

Regardless of initial frequency content, ALL mu=-20 cases converge to:
- E_final = 1.26-1.28
- omega_final = 0.875-0.893
- f_core > 0.97

The oscillon has a single stable attractor. Initial conditions affect only the shedding
phase (t < 2000-3000), not the long-time asymptotic state.

### 4. Low-Frequency Oscillons Do Not Exist (at these parameters)

The proposal hypothesized that frequencies satisfying omega_1 + omega_2 + omega_3 < m
(all Omega_k below the mass gap) might create radiation-free oscillons.

This fails for two reasons:
- **At mu=-20**: the binding is too weak to support oscillons at omega ~ 0.3.
  The fields shed energy until they reach the natural omega ~ 0.87.
- **At mu=-60, -100**: the potential well is too deep; the fields collapse into a
  static (omega=0) configuration with deeply negative energy. This is a false vacuum,
  not an oscillon. The condition m^2 > |mu|*(2/kappa)^(2/3)/9 is violated:
  need m > 1.62 (mu=-60) or m > 2.01 (mu=-100), but m=1.0.

### 5. The Triple-Product Sum Frequency is Irrelevant

The proposal's key idea was that Omega_1 = omega_1 + omega_2 + omega_3 drives radiation,
and splitting the frequencies could push it below the gap. In practice, the coupling
forces all three frequencies to be equal so quickly that the transient frequency content
is irrelevant. The steady-state Omega_1 = 3*omega ~ 2.6 is always above the gap,
and radiation occurs at this harmonic regardless of initialization.

---

## Why Frequency Splitting Cannot Work (Physical Argument)

The triple-product P = phi_1 * phi_2 * phi_3 couples the fields symmetrically.
Any deviation from equal amplitudes/frequencies creates an energy gradient:
the field with higher amplitude drives the others up (and vice versa).

More precisely, the force on field a is:
    F_a = -dV/d(phi_a) = -mu * P * (product of other two) / (1 + kappa*P^2)^2

If phi_1 > phi_2, then F_2 > F_1 (phi_2 sees a larger "background" from phi_1*phi_3),
so phi_2 is driven toward phi_1. This is a robust nonlinear synchronization mechanism
analogous to coupled oscillator phase-locking.

The oscillon frequency is determined by the balance between the nonlinear frequency shift
(from V) and the mass gap (from m^2). Initial conditions can only affect the transient,
not the equilibrium frequency.

---

## Data Files

- `data/threefreq_{label}_ts.tsv` — time series (phi, peaks, energies, f_core)
- `data/threefreq_{label}_spectrum{1,2,3}.tsv` — per-field DFT spectra
- `data/threefreq_{label}_Pspectrum.tsv` — triple-product P(0,t) spectrum
- `data/threefreq_{label}_eslope.tsv` — windowed dE/dt vs time
- `data/threefreq_{label}_summary.tsv` — per-case summary line
- `data/threefreq_summary.tsv` — aggregated comparison table

Labels: control, smallsplit, medsplit, largesplit, lowfreq, mixedlow,
        lowfreq_mu60, lowfreq_mu100, mixedlow_mu60, mixedlow_mu100
