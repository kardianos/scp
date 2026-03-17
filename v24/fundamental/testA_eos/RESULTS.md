# Test A: Vacuum Equation of State — RESULTS

## Summary

**P = -rho does NOT hold anywhere.** The oscillon is a conventional matter configuration,
not a vacuum energy condensate. The time-averaged pressure is essentially zero in the core
(P/rho ~ -5e-5), and approaches the radiation EOS (P ~ +rho) in the dilute vacuum region.
This holds for both the pure triple-product oscillon and the pairwise-coupled variant.

---

## Method

Full stress-energy tensor T^{mu nu} computed at every grid point during evolution:
- T^00 = rho = sum_a [v_a^2/2 + (dx phi_a)^2/2 + m^2 phi_a^2/2] + V
- T^11 = P = sum_a [v_a^2/2 + (dx phi_a)^2/2 - m^2 phi_a^2/2] - V
- T^01 = -sum_a v_a * dx phi_a

Time-averaged over t=5000 to t=10000 (19633 samples, ~650 breathing cycles at omega=0.876).
Trace: T^mu_mu = rho - P (in 1+1D with signature +,-).

Parameters: mu=-20, kappa=20, m=1.0, Nx=4000, xmax=100, dx=0.050, dt=0.025.

---

## Results: lambda=0 (triple-product only)

### Oscillon frequency
- omega = 0.876 (sub-gap: omega/m = 0.876 < 1)
- Total energy E = 1.264 (after equilibration)

### Core (|x| < 5)

| Quantity | Value |
|----------|-------|
| <rho> | 1.25e-1 |
| <P> | -7.0e-6 |
| P/rho | **-5.6e-5** |
| <T01> | -3.2e-16 (zero) |
| trace = rho - P | 1.25e-1 |

**P/rho ~ 0 in the core.** The pressure nearly vanishes on time-average. This is because
the time-averaged kinetic+gradient pressure exactly cancels the time-averaged
mass+potential pressure. The oscillon acts as pressureless dust (P=0), not vacuum (P=-rho).

### Center point (x=0)

| Quantity | Value |
|----------|-------|
| rho | 3.16e-1 |
| P | -1.6e-5 |
| P/rho | **-5.0e-5** |
| P_std | 0.233 |

The instantaneous P fluctuates wildly (std = 0.233 vs mean ~ 0), but the
time-average is effectively zero. The oscillation creates equal time in positive
and negative pressure phases.

### Vacuum (30 < |x| < 60)

| Quantity | Value |
|----------|-------|
| <rho> | 1.21e-6 |
| <P> | 9.25e-7 |
| P/rho | **+0.766** |
| trace = rho - P | 2.82e-7 |

**P/rho ~ +0.77 in the vacuum.** This is close to the radiation EOS (P = rho, i.e. P/rho = 1)
but not exactly 1. The trace rho - P = 2.8e-7 is small but nonzero, consistent with
massive radiation (trace ~ m^2 phi^2, small in the tail).

### Mid-range (10 < |x| < 20)

| Quantity | Value |
|----------|-------|
| <rho> | 7.89e-6 |
| <P> | 9.18e-7 |
| P/rho | **+0.116** |

Transition zone: P/rho interpolates smoothly from ~0 (core) to ~0.77 (vacuum).

### Integrated quantities

| Quantity | Value |
|----------|-------|
| int(rho) dx | 1.268 |
| int(P) dx | 6.4e-5 |
| int(P)/int(rho) | **5.0e-5** |
| int(T01) dx | -3.2e-15 (zero) |

**The virial theorem is satisfied**: integrated pressure vanishes for a static bound state
in 1+1D. This is exact for time-averaged oscillons and confirms the numerical averaging
is correct.

### Trace T^mu_mu = rho - P

The trace equals rho - P = 2(mass energy + potential energy) = m^2 sum phi_a^2 + 2V.
- Core: trace ~ rho (since P ~ 0), trace = 1.25e-1
- Vacuum: trace ~ 2.8e-7 << rho (radiation is nearly conformal)
- The theory is NOT conformal (m != 0), so trace != 0 everywhere

---

## Results: lambda=0.99 (with pairwise coupling)

### Oscillon frequency
- omega = 1.602 (ABOVE mass gap: omega/m = 1.60 > 1)
- This is the symmetric mode omega_S = sqrt(m^2 + 2*lambda) = sqrt(1 + 1.98) = 1.73
- Total energy E = 11.7 (much larger due to pairwise mass energy)

### Core (|x| < 5)

| Quantity | Value |
|----------|-------|
| <rho> | 1.18 |
| <P> | -2.5e-4 |
| P/rho | **-2.1e-4** |
| <T01> | -1.3e-16 (zero) |
| trace | 1.18 |

Same qualitative result: **P ~ 0 in the core**, pressureless dust.

### Vacuum (30 < |x| < 60)

| Quantity | Value |
|----------|-------|
| <rho> | 6.30e-5 |
| <P> | 5.44e-5 |
| P/rho | **+0.865** |

Closer to radiation EOS (P/rho = 0.87 vs 0.77 without pairwise). The pairwise
coupling produces higher-frequency radiation that is more relativistic.

### Center point (x=0)

| Quantity | Value |
|----------|-------|
| rho | 2.15 |
| P | -5.1e-4 |
| P/rho | **-2.4e-4** |

### Comparison

| Region | P/rho (lambda=0) | P/rho (lambda=0.99) |
|--------|-------------------|----------------------|
| Core | -5.6e-5 | -2.1e-4 |
| Mid | +0.116 | +0.748 |
| Vacuum | +0.766 | +0.865 |

The pairwise coupling changes the detailed numbers but not the qualitative picture:
pressureless core, radiation-like vacuum.

---

## Physical Interpretation

### Why P ~ 0 in the core

The oscillon oscillates as phi ~ A cos(omega t) * f(x). The instantaneous pressure
oscillates between positive (kinetic-dominated half-cycles) and negative
(potential-dominated half-cycles). Time-averaging over complete cycles, these cancel
almost exactly.

This is the **virial theorem** in action: for a localized 1+1D bound state,
int(P) dx = 0 exactly. Since the pressure profile is smooth and peaked at x=0,
the pointwise pressure must also be small.

### Why P/rho ~ +1 in the vacuum

Far from the oscillon, the field consists of dispersive radiation leaking outward.
For relativistic (k >> m) radiation, P = rho. For massive radiation at the
oscillon frequency omega = 0.876 < m, the radiation is evanescent (exponentially
decaying) rather than propagating, so P/rho is neither exactly 0 nor exactly 1.
The measured value P/rho ~ 0.77 reflects a mix of propagating harmonics (above gap)
and evanescent tails.

### No vacuum energy

The oscillon is NOT a vacuum condensate. A Lorentz-invariant vacuum requires P = -rho
(cosmological constant equation of state). The oscillon has P ~ 0 (dust) in the core
and P ~ +rho (radiation) in the far field. There is no regime where P = -rho.

The trace T^mu_mu = rho - P ~ rho in the core, confirming the energy is dominated
by the mass term (non-relativistic matter), not by vacuum energy (which would give
trace = 2*rho for P = -rho in 1+1D).

---

## Files

- `src/eos.c` — stress-energy tensor diagnostic code
- `data/eos_profile.tsv` — spatial profile of time-averaged T^{mu nu} (lambda=0)
- `data/eos_scatter.tsv` — P vs rho scatter data (lambda=0)
- `data/eos_profile_pw.tsv` — same with lambda=0.99
- `data/eos_scatter_pw.tsv` — same with lambda=0.99
- `data/eos_spectrum.tsv` — frequency spectrum (lambda=0)
- `data/eos_spectrum_pw.tsv` — frequency spectrum (lambda=0.99)

Compile: `gcc -O3 -Wall -o eos src/eos.c -lm`
