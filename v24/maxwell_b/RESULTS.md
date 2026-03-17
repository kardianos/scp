# V24-MB Results: Three Complex Fields with Quark Charges

## Summary

The U(1) gauge coupling **destroys the oscillon** at all nonzero e_scale tested.
Even e_scale=0.1 is sufficient to destabilize the bound state. The gauge field
acts as a radiation channel that drains energy from the scalar sector. The baryon
does not carry a persistent EM charge because the oscillon disperses before any
charge structure can establish.

## Phase 1: Control (e_scale=0)

| Quantity       | Value  |
|----------------|--------|
| omega_osc      | 0.900  |
| mass gap m     | 1.000  |
| Oscillon?      | YES    |
| E_total(t=0)   | 3.365  |
| E_total(t=10k) | 2.769  |
| f_core(t=10k)  | 0.999  |
| phi_0 amplitude| ~0.67  |
| E_em           | 0      |

The control case reproduces v21: a long-lived oscillon with omega=0.90 below the
mass gap, energy slowly leaking via radiation. Core fraction stays >0.99.

## Phase 2: e_scale Scan (equal masses m=1)

| e_scale | charges (e1,e2,e3) | Q_total | E_final | E_em    | fc    | omega_peak | Oscillon? |
|---------|-------------------|---------|---------|---------|-------|------------|-----------|
| 0.0     | (0, 0, 0)         | 0       | 2.769   | 0       | 0.999 | 0.900      | YES       |
| 0.1     | (+0.067,+0.067,-0.033) | +0.1 | 1.974 | 0.005  | 0.117 | 1.002      | NO        |
| 0.3     | (+0.20,+0.20,-0.10)   | +0.3 | 0.666 | 0.013  | 0.180 | 1.002      | NO        |
| 1.0     | (+0.67,+0.67,-0.33)   | +1.0 | 0.025 | 0.024  | 0.026 | 1.020      | NO        |

### Key findings:

1. **Oscillon destroyed by gauge coupling**: At e_scale=0.1, the oscillation
   frequency shifts above the mass gap (omega=1.002 > m=1.0). The scalar field
   can now radiate freely. By t=10000, the core fraction drops from 1.0 to 0.12.

2. **Complete dispersal at e_scale=1.0**: At full quark charges, E_total drops
   from 3.37 to 0.025 by t=10000, of which 0.024 is EM field energy. The scalar
   fields are essentially zero (phi < 10^{-4}). The system becomes pure EM radiation.

3. **EM energy fraction grows with e_scale**: The EM field absorbs energy from
   the scalars. At e=1.0, 98% of remaining energy is electromagnetic.

4. **Gauge field A develops O(1) amplitude**: Even at e=0.1, A reaches ~0.5-1.0
   at the center, which is comparable to the scalar field values. This is because
   the current j ~ e*phi^2*A feeds back strongly.

5. **Chi remains tiny**: The imaginary parts chi_a stay at O(10^{-4}) or smaller
   throughout. The seed perturbation does not grow — instead the real parts phi_a
   simply disperse.

6. **No persistent charge**: Q_em (from Gauss law) is O(10^{-280}), i.e., zero
   to machine precision. The chi seed is antisymmetric (odd in x) so there is no
   net charge. The "baryon charge +1" from e1+e2+e3 does not manifest as a
   measurable EM charge because the oscillon disperses before establishing any
   charge asymmetry.

## Phase 3: UUD vs UDD Mass Comparison

Masses: m_U=1.0, m_D=0.95 (5% quark mass splitting).

### e_scale=0.3

| Config | Charges                | Q_total | E_final | E_em    | fc    |
|--------|------------------------|---------|---------|---------|-------|
| UUD    | (+0.20,+0.20,-0.10)   | +0.3    | 0.626   | 0.0008  | 0.080 |
| UDD    | (+0.20,-0.10,-0.10)   | 0.0     | 0.623   | 0.0014  | 0.042 |

### e_scale=1.0

| Config | Charges                | Q_total | E_final | E_em    | fc    |
|--------|------------------------|---------|---------|---------|-------|
| UUD    | (+0.67,+0.67,-0.33)   | +1.0    | 0.024   | 0.024   | 0.011 |
| UDD    | (+0.67,-0.33,-0.33)   | 0.0     | 0.313   | 0.004   | 0.026 |

### Key findings:

1. **At e=0.3**: UUD and UDD have nearly identical final energies (0.626 vs 0.623).
   The EM self-energy difference is only delta_E_em = 0.0006, far too small to
   produce meaningful mass splitting. Both oscillons are already dispersed (fc < 0.1).

2. **At e=1.0**: The UUD (proton analog) disperses MORE than UDD (neutron analog).
   UUD retains only E=0.024 (essentially all EM) while UDD retains E=0.313.
   This is because UUD has charges (+2/3,+2/3,-1/3) with sum +1 (net charge),
   while UDD has (+2/3,-1/3,-1/3) with sum 0 (neutral). The charged UUD radiates
   more efficiently via the gauge field.

3. **Wrong mass ordering**: In nature, M_neutron > M_proton by 1.3 MeV, with
   EM self-energy making the proton lighter despite the quark mass contribution.
   Here we see the opposite effect: the charged configuration (UUD) loses MORE
   energy via EM radiation, ending up lighter — but only because both are
   completely dispersed. There is no stable bound state to compare.

## Mechanism of Destruction

The gauge coupling destroys the oscillon through a positive feedback loop:

1. Scalar oscillation at frequency omega ~ 0.9m generates a current j ~ e*phi*dchi/dx
2. The current sources the gauge field: d_t^2 A = d_x^2 A - j
3. The gauge field introduces an effective mass correction: the term -e^2*A^2*phi
   shifts the effective frequency of phi upward
4. When omega_eff > m, the mass gap protection is lost and phi radiates freely
5. As phi disperses, the current drops, but the damage is already done

The critical insight is that even a small chi seed (10^{-3}) combined with the
gauge backreaction (e^2*A^2*phi term) is enough to push the frequency above the
mass gap.

## Does the Baryon Carry Charge +1?

**Formally yes, operationally no.** The charges are assigned as e1+e2+e3 = +1
for UUD, and the Lagrangian is gauge-invariant. However:

- The oscillon is not stable with gauge coupling on
- No persistent charge density develops (chi stays at seed level ~10^{-4})
- Gauss-law charge Q is zero to machine precision
- The "baryon" disperses into free scalar + EM radiation

For the charge to manifest, the oscillon would need to survive the gauge coupling.
This might require either (a) much weaker coupling (e_scale << 0.1), or (b) a
topological stabilization mechanism absent in the 1D model.

## Does EM Self-Energy Affect Mass Ordering?

**No, because there is no stable bound state.** The EM coupling is too strong
relative to the scalar binding: even e_scale=0.1 destroys the oscillon. Without
a stable object, "mass" is not well-defined.

The trend at e=1.0 (charged UUD disperses faster than neutral UDD) is consistent
with the physical intuition that EM self-energy adds to the proton mass, but the
effect is catastrophic rather than perturbative.

## Parameters

- mu = -20, kappa = 20
- m = 1.0 (equal mass scan), m_U = 1.0, m_D = 0.95 (UUD/UDD)
- A_init = 0.8, sigma = 3.0
- chi_seed = 10^{-3} (odd perturbation x*exp(-x^2/2sigma^2))
- Nx = 4000, xmax = 100, dx = 0.05, dt = 0.0127
- tfinal = 10000

## Files

- `src/maxwell_b.c` — solver (7 field arrays, Velocity Verlet, absorbing BC)
- `data/UUD_e{0.00,0.10,0.30,1.00}_ts.tsv` — time series for e_scale scan
- `data/UUD_proton_e{0.30,1.00}_ts.tsv` — UUD proton comparison
- `data/UDD_neutron_e{0.30,1.00}_ts.tsv` — UDD neutron comparison
- `data/*_spectrum.tsv` — DFT power spectra
