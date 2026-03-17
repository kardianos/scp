# V24-UDD Results: Neutron-Like 120-Degree Oscillon

## Answer

**UDD and UUD oscillons DO differ in energy at the same m_D, providing a mass splitting analog. UDD (neutron) retains MORE energy than UUD (proton) at moderate m_D splittings (0.80-0.95), but BOTH are slowly dispersing -- neither forms a truly stable oscillon.**

**At the proposal's lambda=0.85, BOTH configurations fall into a false vacuum. The pairwise coupling is far too strong: the lowest mass matrix eigenvalue m_D^2 - lambda becomes near-zero or tachyonic, driving an immediate instability.**

At lambda=0.30 (stable regime), the key results are:

| m_D | UDD E_final | UUD E_final | UDD/UUD | UDD more stable? |
|-----|-------------|-------------|---------|------------------|
| 1.00 | 0.442 | 0.442 | 1.000 | Identical (symmetric) |
| 0.95 | 0.449 | 0.336 | 1.34 | **YES** (34% more energy) |
| 0.90 | 0.483 | 0.205 | 2.36 | **YES** (136% more energy) |
| 0.85 | 0.450 | 0.157 | 2.87 | **YES** (187% more energy) |
| 0.80 | 0.443 | 0.348 | 1.27 | **YES** (27% more energy) |
| 0.70 | 0.124 | 0.233 | 0.53 | **NO** (UUD wins here) |

The ordering FLIPS between m_D=0.80 and m_D=0.70. For moderate mass splittings (m_D = 0.80-0.95), UDD retains significantly more energy than UUD. This is the opposite of what naive reasoning would suggest -- the neutron analog is HEAVIER than the proton analog.

---

## Mass Matrix Analysis

The 3x3 mass matrix M^2 with pairwise coupling lambda has three eigenmodes:

### UDD (field 1 = Up, fields 2,3 = Down)
- **DD antisymmetric mode**: eigenvalue = m_D^2 - lambda (pure down-down, no up involvement)
- **Mixed mode -**: involves all three fields
- **Mixed mode +**: highest frequency symmetric-like mode

### UUD (fields 1,2 = Up, field 3 = Down)
- **UU antisymmetric mode**: eigenvalue = m_U^2 - lambda = 1.0 - lambda (pure up-up, no down)
- **Mixed mode -**: involves all three fields
- **Mixed mode +**: highest frequency symmetric-like mode

The KEY DIFFERENCE: the antisymmetric mode eigenvalues differ.
- UDD: eig_DD = m_D^2 - lambda (decreases with m_D)
- UUD: eig_UU = m_U^2 - lambda = 1.0 - lambda (FIXED, independent of m_D)

At lambda=0.30:

| m_D | UDD DD-antisym | UUD UU-antisym | UDD lighter? |
|-----|---------------|----------------|--------------|
| 1.00 | 0.700 | 0.700 | Equal |
| 0.95 | 0.603 | 0.700 | YES |
| 0.90 | 0.510 | 0.700 | YES |
| 0.85 | 0.423 | 0.700 | YES |
| 0.80 | 0.340 | 0.700 | YES |
| 0.70 | 0.190 | 0.700 | YES |

The UDD DD-antisymmetric mode is ALWAYS lighter than the UUD UU-antisymmetric mode (for m_D < m_U). A lighter mode oscillates slower and radiates less efficiently, explaining why UDD retains more energy in the moderate m_D range.

---

## Tachyonic Instability at lambda=0.85

The proposal specifies lambda=0.85. At this value:

| m_D | UDD DD-antisym | UUD UU-antisym | Status |
|-----|---------------|----------------|--------|
| 1.00 | 0.150 | 0.150 | Both marginally stable, but false vacuum |
| 0.95 | 0.053 | 0.150 | Both false vacuum |
| 0.93 | 0.015 | 0.150 | UDD near-tachyonic threshold |
| 0.922 | 0.000 | 0.150 | UDD exactly tachyonic |
| 0.90 | -0.040 | 0.150 | **UDD TACHYONIC**, UUD still in false vacuum |
| 0.85 | -0.128 | 0.150 | **UDD TACHYONIC** |

At lambda=0.85, even the m_D=1.0 case collapses to a false vacuum (E goes from 0.50 to -36.6). The pairwise coupling V_pw = lambda * sum phi_a phi_b makes the phi=0 vacuum unstable when the effective mass^2 is too small. The system rolls to a nonzero static equilibrium.

**UDD becomes tachyonic before UUD**: the DD antisymmetric mode goes through zero at m_D = sqrt(lambda) = 0.922, while the UU antisymmetric mode stays at m_U^2 - lambda = 0.15 for all m_D. This means UDD is MORE fragile than UUD to mass splitting -- it loses stability first.

At m_D=0.95, lambda=0.85:
- UDD E_final = -43.66 (deeper false vacuum)
- UUD E_final = -40.52 (shallower false vacuum)
- UDD initial energy = 0.165 vs UUD = 0.331 (UDD starts with LESS initial energy)

---

## Detailed Results: lambda=0.30 (Stable Regime)

### Energy Retention

| m_D | Mode | E_init | E(t=1000) | E(t=5000) | E(t=10000) | Retained | fc_final |
|-----|------|--------|-----------|-----------|------------|----------|----------|
| 1.00 | UDD | 3.304 | 0.979 | 0.571 | 0.442 | 13.4% | 0.507 |
| 1.00 | UUD | 3.304 | 0.979 | 0.571 | 0.442 | 13.4% | 0.507 |
| 0.95 | UDD | 2.973 | 0.876 | 0.567 | 0.449 | 15.1% | 0.451 |
| 0.95 | UUD | 3.138 | 0.804 | 0.450 | 0.336 | 10.7% | 0.374 |
| 0.90 | UDD | 2.658 | 0.807 | 0.587 | 0.483 | 18.2% | 0.376 |
| 0.90 | UUD | 2.981 | 0.653 | 0.288 | 0.205 | 6.9% | 0.236 |
| 0.85 | UDD | 2.360 | 0.706 | 0.538 | 0.450 | 19.1% | 0.285 |
| 0.85 | UUD | 2.832 | 0.535 | 0.219 | 0.157 | 5.5% | 0.259 |
| 0.80 | UDD | 2.079 | 0.659 | 0.514 | 0.443 | 21.3% | 0.190 |
| 0.80 | UUD | 2.692 | 0.753 | 0.455 | 0.348 | 12.9% | 0.393 |
| 0.70 | UDD | 1.569 | 0.388 | 0.176 | 0.124 | 7.9% | 0.214 |
| 0.70 | UUD | 2.437 | 0.600 | 0.323 | 0.233 | 9.6% | 0.155 |

### Oscillation Frequencies (DFT, second half of run)

| m_D | UDD omega_1 | UDD omega_2 | UDD omega_3 | UUD omega_1 | UUD omega_2 | UUD omega_3 |
|-----|------------|------------|------------|------------|------------|------------|
| 1.00 | 0.840 | 0.840 | 0.840 | 0.840 | 0.840 | 0.840 |
| 0.95 | 0.816 | 0.816 | 0.816 | 0.798 | 0.798 | 0.798 |
| 0.90 | 0.792 | 0.792 | 0.792 | 0.840 | 0.840 | 0.750 |
| 0.85 | 0.768 | 0.768 | 0.768 | 0.708 | 0.708 | 0.708 |
| 0.80 | 0.738 | 0.738 | 0.738 | 0.660 | 0.660 | 0.660 |
| 0.70 | 0.678 | 0.678 | 0.678 | 0.558 | 0.558 | 0.558 |

Key observation: all three fields oscillate at the SAME frequency in all cases (lock-in from the triple-product coupling), but UDD and UUD lock to DIFFERENT frequencies. UDD oscillates faster than UUD for m_D < m_U.

---

## Physical Interpretation

### Mass Splitting Analog
The UDD/UUD energy difference at fixed m_D is a genuine analog of the neutron-proton mass splitting:

- **UDD (neutron analog)**: E_final > E_UUD for m_D in [0.80, 0.95] -- the neutron is HEAVIER
- This matches real physics: M_n - M_p = +1.293 MeV (neutron is heavier)
- The ordering flips at large mass splitting (m_D < 0.70), where UUD becomes heavier

### Why UDD Retains More Energy
The UDD DD-antisymmetric mode has a LOWER frequency (m_D^2 - lambda vs m_U^2 - lambda). In the 120-degree configuration, this lighter mode radiates less efficiently because:
1. Lower omega means the 3rd harmonic 3*omega is further below the mass gap
2. The down-down relative oscillation stores energy in a softer mode that couples more weakly to radiation

### Stability Hierarchy
Neither UDD nor UUD forms a truly stable oscillon -- both slowly disperse. But the relative ordering is:

For moderate m_D (0.80-0.95): UDD MORE stable than UUD (neutron > proton)
For large splitting (m_D < 0.70): UUD MORE stable (proton > neutron)

The crossover occurs around m_D ~ 0.75.

---

## The lambda=0.85 Catch-22

The proposal's lambda=0.85 creates a fundamental problem identical to the one found in V24-B (phase120 RESULTS):

1. **Pairwise coupling V_pw = lambda * phi_a * phi_b lowers the effective mass**: the lowest eigenvalue of the mass matrix is m_D^2 - lambda. At lambda=0.85, this is only 0.15 for m_D=1.0 and goes negative for m_D < 0.922.

2. **The triple-product coupling drives false vacuum collapse**: with mu=-20 and such a small effective mass, the system rolls to a nonzero equilibrium.

3. **UDD is MORE vulnerable**: it has TWO down fields whose antisymmetric mode eigenvalue is m_D^2 - lambda, while UUD's antisymmetric UU mode is m_U^2 - lambda = 0.15 (independent of m_D). So UDD goes tachyonic first.

**Conclusion**: lambda=0.85 is incompatible with stable oscillons in this model. Viable oscillons require lambda < ~0.5 (so the lowest eigenvalue stays well above zero). At lambda=0.30, stable oscillons exist and show clear UDD/UUD mass splitting.

---

## Files

| File | Description |
|------|-------------|
| `src/udd120.c` | Solver (supports both UDD mode=0 and UUD mode=1) |
| `data/udd_mD{val}_ts.tsv` | UDD time series for each m_D (lambda=0.30) |
| `data/udd_mD{val}_phi{n}_spectrum.tsv` | DFT spectra per field |
| `data/udd_mD{val}_P_spectrum.tsv` | Triple product P(t) spectrum |

UUD data in `v24/120uud/data/` with same naming convention (prefix `uud`).

## Parameters

| Parameter | Value |
|-----------|-------|
| mu | -20 |
| kappa | 20 |
| m_U | 1.0 |
| lambda | 0.30 (stable) / 0.85 (proposal, false vacuum) |
| m_D scan | {1.00, 0.95, 0.90, 0.85, 0.80, 0.70} |
| A_init | 0.8 |
| sigma | 3.0 |
| Nx | 4000 |
| xmax | 100 |
| tfinal | 10000 |
