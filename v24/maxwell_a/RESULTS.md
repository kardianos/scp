# V24-MA Results: Complex Scalar + U(1) Gauge Field

## Summary

**The triad oscillon does NOT carry a conserved U(1) charge.** The triple-product coupling P = phi_1 phi_2 phi_3 breaks gauge invariance, forcing all three fields into phase-lock. Any initial charge Q is equalized to zero within a few oscillation periods. The gauge field develops oscillating structure at the oscillon core (peak_A ~ 0.5-1.0 for e=0.1) but remains a small perturbation. At e >= 0.15, the gauge field grows without bound and destroys the oscillon.

**Stability boundary: e_crit ~ 0.12.** Below this, the oscillon survives to t=10000+ with negligible charge. Above it, the e^2 A^2 |Psi|^2 feedback loop drives runaway growth of A.

---

## Model

Four fields: phi_1, phi_2 (charged under U(1)), phi_3 (neutral), A (gauge, temporal gauge A_0=0).

    L = 1/2|D_t Psi|^2 - 1/2|D_x Psi|^2 - 1/2 m^2|Psi|^2
      + 1/2(dt phi_3)^2 - 1/2(dx phi_3)^2 - 1/2 m^2 phi_3^2
      - V(P) - 1/4 F^2

    Psi = phi_1 + i phi_2,  D_mu = d_mu - ieA_mu

    V = (mu/2) P^2 / (1 + kappa P^2),  P = phi_1 phi_2 phi_3

    j = e(phi_1 dx phi_2 - phi_2 dx phi_1) - e^2 A |Psi|^2

    Q = integral (phi_1 v_2 - phi_2 v_1) dx  [temporal gauge]

**NOTE**: P = phi_1 phi_2 phi_3 is NOT gauge invariant. Under Psi -> e^{i alpha} Psi:
phi_1 -> phi_1 cos(alpha) - phi_2 sin(alpha), phi_2 -> phi_1 sin(alpha) + phi_2 cos(alpha).
This means P changes, so Q is NOT exactly conserved even classically. The gauge coupling
is treated as a perturbation on top of the known binding mechanism.

**Parameters**: mu=-20, kappa=20, m=1.0, A=0.8, sigma=3.0 (v21 optimized).
Grid: Nx=4000, xmax=100, dx=0.05, dt=0.013. Absorbing boundary at |x| > 75.

### Gauge-invariant coupling attempted

The proposal specifies V = (mu/2)(|Psi|^2 phi_3^2)^2 / (1 + kappa (|Psi|^2 phi_3^2)^2).
This IS gauge-invariant but is 8th order in phi at small amplitude (vs 6th for the
triple product). With mu=-20, kappa=20: the coupling is too weak to form an oscillon.
The field disperses like a free massive field.

With V = (mu/2)|Psi|^2 phi_3^2 / (1 + kappa |Psi|^2 phi_3^2) (4th order), the coupling
is strong enough but creates a false vacuum at phi != 0 for mu=-20, kappa=20.
Increasing kappa to 50 stabilizes the vacuum, but the resulting oscillon is weak
(peak amplitudes ~0.05, f_core < 0.5). The biquadratic coupling is fundamentally a
two-body interaction (only phi_1 and phi_3 participate when phi_2=0), losing the
three-body resonance that makes the triad oscillon robust.

**Conclusion**: Gauge-invariant coupling does not support the triad oscillon mechanism.
We proceed with the triple-product coupling to study gauge field dynamics on a known
oscillon background.

---

## Phase 1: Oscillon with gauge coupling

### Initialization

Standard triad: phi_1 = phi_2 = phi_3 = A exp(-x^2/(2 sigma^2)), A_mu = 0.
Charge seed: v_2 += twist * omega_0 * phi_1 (velocity kick on phi_2 only).
Gauss constraint: vA(x) = -integral_{-inf}^{x} rho dx' enforced at t=0.

### e = 0 (control)

| Time | E_total | E_pot | fc | omega | Q | Notes |
|------|---------|-------|-----|-------|---|-------|
| 0 | 3.375 | -2.023 | 1.00 | — | 0.259 | Initial |
| 100 | 3.187 | -1.028 | 0.88 | — | 0.001 | Q equalized |
| 500 | 2.655 | -0.030 | 0.99 | — | 0.003 | Shedding |
| 5000 | 1.375 | -0.087 | 1.00 | — | 0.000 | Stable |
| 10000 | 1.264 | -0.171 | 1.00 | 0.876 | 0.000 | Stable oscillon |

**Key**: Q drops from 0.26 to ~0 within t < 100. The triple-product coupling forces
phi_1 and phi_2 into phase-lock (same amplitude, same phase), killing the charge.
omega = 0.876 < m = 1.0 (sub-gap, radiation forbidden).

### e = 0.05

| Time | E_total | E_pot | E_EM | fc | Q | peak_A | Notes |
|------|---------|-------|------|-----|---|--------|-------|
| 0 | 3.383 | -2.023 | 0.008 | 1.00 | 0.259 | 0.00 | |
| 500 | 2.651 | -0.044 | 0.003 | 0.99 | 0.002 | 0.23 | |
| 5000 | 1.383 | -0.049 | 0.003 | 1.00 | 0.000 | 0.19 | |
| 10000 | 1.265 | -0.059 | 0.003 | 1.00 | 0.000 | 0.16 | Stable |

omega = 0.876. Identical to e=0 within errors. E_EM ~ 0.003 (0.2% of total).
Gauge field oscillates with peak_A ~ 0.15-0.25. Q = 0 throughout.

### e = 0.1

| Time | E_total | E_pot | E_EM | fc | Q | peak_A | Notes |
|------|---------|-------|------|-----|---|--------|-------|
| 0 | 3.408 | -2.023 | 0.033 | 0.99 | 0.259 | 0.00 | |
| 100 | 3.346 | -1.571 | 0.008 | 0.87 | 0.009 | 1.43 | A develops |
| 500 | 1.451 | -0.349 | 0.002 | 0.97 | -0.007 | 1.07 | |
| 5000 | 1.291 | -0.058 | 0.006 | 1.00 | -0.004 | 0.62 | |
| 10000 | 1.229 | 0.000 | 0.003 | 1.00 | -0.003 | 0.82 | Stable |

omega = 0.888. E_EM ~ 0.003-0.033 (0.2-1% of total). Gauge field peak_A ~ 0.5-1.4
(significant structure). Q oscillates around 0 with amplitude ~0.005.
**Gauss violation grows to ~0.007** but stabilizes (not catastrophic).

### e = 0.15 (UNSTABLE)

Gauge field A grows to 2-6 within t < 600. Energy increases (pumped by scalar→gauge
transfer). Oscillon destroyed by t ~ 600. NaN at t ~ 700.

### e >= 0.3 (ALL UNSTABLE)

Immediate blowup. Even with twist=0.01, NaN within t < 200.
The e^2 A^2 |Psi|^2 term creates positive feedback: A grows → effective mass of
Psi increases → stronger current → A grows faster.

**With twist=0 (no charge seed)**: All e values are stable. The symmetric
initial condition (phi_1 = phi_2) keeps Q = 0 exactly, and A = 0 exactly (no source).
The instability is triggered by ANY asymmetry between phi_1 and phi_2.

---

## Phase 3: Charge conservation

**Q is NOT conserved** — not even approximately. The initial Q = 0.26 drops to ~0
within ~100 time units (a few oscillation periods). This is NOT a gauge-coupling
effect: it happens identically at e = 0. The mechanism:

The triple-product coupling P = phi_1 phi_2 phi_3 acts as a strong restoring force
that equalizes all three fields. When phi_2 has a different velocity from phi_1,
the coupling drives them back into phase-lock. Since Q = integral(phi_1 v_2 - phi_2 v_1)dx,
phase-locked fields (v_1/phi_1 = v_2/phi_2) give Q = 0.

**Implication**: The triad oscillon is intrinsically uncharged. The three-body
resonance mechanism REQUIRES all three fields to oscillate in phase, which is
incompatible with carrying U(1) charge (which requires a phase difference between
phi_1 and phi_2).

---

## Gauge field structure

At e = 0.1, the gauge field A develops non-trivial structure despite Q = 0:

- peak_A oscillates between 0.1 and 1.4, correlated with the oscillon breathing
- E_EM ~ 0.003-0.033 (small fraction of total energy)
- The A profile tracks the oscillon core (extends ~3 sigma)
- Gauss violation ~0.007 (manageable but growing slowly)

The gauge field acts as a **slave mode**: driven by the oscillating charged current
j = e(phi_1 dx phi_2 - phi_2 dx phi_1), it follows the scalar dynamics without
significantly back-reacting. In 1D, there is no magnetic energy (B = 0), so A is
free to grow; only the e^2 A^2 |Psi|^2 screening term provides a restoring force.

---

## Key findings

1. **The triad oscillon cannot carry U(1) charge.** The triple-product coupling
   forces phase-locking, which is incompatible with charge.

2. **Gauge-invariant coupling (|Psi|^2 phi_3^2) does not support oscillons.**
   The biquadratic form is either too weak (S^2) or creates false vacua (S).

3. **Gauge field develops structure** at e = 0.05-0.1 but remains a small
   perturbation (E_EM < 1% of total).

4. **Stability boundary at e_crit ~ 0.12.** Above this, the e^2 A^2 |Psi|^2
   feedback loop drives gauge field runaway. Any phi_1 != phi_2 asymmetry
   triggers the instability.

5. **In the uncharged sector (twist=0), all e values are stable** because the
   symmetric initial condition keeps A = 0 exactly.

6. **Omega shifts slightly with e**: 0.876 (e=0) → 0.876 (e=0.05) → 0.888 (e=0.1).
   The gauge coupling raises the effective oscillation frequency but keeps it below
   the mass gap.

---

## Implications for the SCP program

The triad oscillon is an intrinsically neutral object. To obtain charged oscillons,
one needs either:

- A **Q-ball**-type mechanism: single complex field with V(|Psi|^2) that has a
  global minimum at |Psi| != 0, giving Psi = f(r) e^{i omega t}. This naturally
  carries charge Q = omega integral |f|^2 dx.

- A coupling that does NOT force phase-locking: e.g., separate potentials for
  each field, with the gauge coupling providing the ONLY interaction between phi_1
  and phi_2.

The fundamental tension: three-body resonance (P = phi_1 phi_2 phi_3) requires
all three in-phase. U(1) charge requires phi_1 and phi_2 out-of-phase. These are
incompatible.

---

## Files

| File | Description |
|------|-------------|
| `src/maxwell_a.c` | Solver: triad + U(1) gauge field |
| `data/maxwell_e{X}_ts.tsv` | Time series for each e value |
| `data/maxwell_e{X}_spectrum.tsv` | DFT power spectrum |
| `data/maxwell_e{X}_profile.tsv` | Final field profiles |
| `PROPOSAL.md` | Original specification |
