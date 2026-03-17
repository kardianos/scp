# V24-B Results: 120-Degree Phase-Separated Oscillon

## Answer

**The 120-degree phase configuration does NOT live longer than the 0-degree symmetric oscillon. It disperses catastrophically fast.**

At the only parameter set that produces a genuine sub-gap oscillon (mu=-20, kappa=20, m=1.0):
- **Control (0-phase)**: retains 98.9% of energy over t=10000, fc=1.000 throughout
- **120-degree state**: retains only 17.0% of energy, fc drops to ~0.14 by t=10000

The 120-degree state is not just slightly worse -- it is qualitatively dead. The symmetric oscillon is rock-solid while the 120-degree state disperses into free radiation.

---

## Why It Fails

The proposal's reasoning about harmonic elimination is correct in principle but misses a fatal problem: **the 120-degree state has 16x weaker binding energy** (P_max = f^3/4 vs f^3 for symmetric). This means the nonlinear frequency shift is 16x weaker, which is insufficient to push omega below the mass gap. Without sub-gap protection, the fields radiate freely.

Specifically:
- Symmetric state: P = f^3 cos^3(wt), with P_max = f^3 -> strong frequency shift -> omega = 0.876 < m = 1.0
- 120-degree state: P = (f^3/4) cos(3wt), with P_max = f^3/4 -> weak frequency shift -> omega = 1.002 >= m = 1.0

The 120-degree DFT confirms: phi_1 oscillates at omega = 1.002 (AT the mass gap, not below it), while the control phi_1 oscillates at omega = 0.876 (well below). The 120-degree state cannot sustain sub-gap oscillation and simply disperses.

The P(t) spectrum confirms the harmonic elimination works as predicted:
- Control P(t): peak at omega = 0.88 (fundamental)
- 120-degree P(t): peak at omega = 3.0 (3rd harmonic only, as predicted by cos(theta)*cos(theta+2pi/3)*cos(theta+4pi/3) = (1/4)cos(3*theta))

So the harmonic structure is exactly as predicted, but the binding is too weak to matter.

---

## Detailed Results

### Phase 1: Equilibration (mu=-20, kappa=20, m=1.0)

Equilibrated for t=5000 with symmetric initial conditions. The oscillon sheds ~63% of initial energy and settles into a breathing mode:

| Metric | Value |
|--------|-------|
| omega_equil | 0.870 |
| E_equil | 1.274 |
| f_core | 1.000 |
| Peak amplitude | 0.35-0.50 (breathing) |
| omega/m | 0.870 (13% gap margin) |

### Phase 2: 120-Degree State Evolution (t=10000)

| Metric | Control (0-phase) | 120-degree | Ratio |
|--------|-------------------|------------|-------|
| E(0) | 1.274 | 1.318 | 1.03 |
| E(10000) | 1.259 | 0.224 | 0.18 |
| E retained | **98.9%** | **17.0%** | -- |
| f_core(10000) | **1.000** | 0.14 | -- |
| omega (DFT) | **0.876** | 1.002 | -- |
| Peak amp (late) | 0.47 | 0.03 | -- |
| Ep(t) | -0.25 typ. | -0.00 | -- |
| Phase stability | N/A (0-phase) | Phases drift but roughly maintain 120 | -- |

The 120-degree state's energy drops rapidly in the first ~500 time units (1.32 -> 0.79), then continues declining slowly. The potential energy Ep = 0 throughout, confirming there is no effective binding.

### Phase 3: Strong Binding Regime

#### mu=-60, kappa=20, m=1.0 (FALSE VACUUM)

**Invalid**: m_crit = 1.198 > m = 1.0. The phi=0 vacuum is unstable. Fields settle to a nonzero static equilibrium (not an oscillon). The 120-degree state immediately collapses to the symmetric state (phases snap to ~180/0 within t~250).

| Metric | Control | 120-degree |
|--------|---------|------------|
| omega | 0 (static) | 0 (static) |
| E(10000) | -49.8 | -46.8 |
| Phase collapse | N/A | Immediate (<250 t.u.) |

#### mu=-100, kappa=20, m=1.0 (FALSE VACUUM)

Same pathology, even more extreme. E = -176 (deeply bound false vacuum). Phase collapse within t~250.

#### mu=-40, kappa=20, m=1.0 (NEAR FALSE VACUUM BOUNDARY)

m_crit = 0.98, just barely below m=1.0. The oscillon sits right at the mass gap (omega = 1.002). Both states slowly disperse:

| Metric | Control | 120-degree |
|--------|---------|------------|
| omega | 1.002 | 1.002 |
| E retained | 87.1% | 77.1% |
| f_core | 0.14 | 0.22 |
| Phase stability | N/A | **Perfectly stable** at 120 |

Neither is a true sub-gap oscillon. The 120-degree phases are perfectly maintained because the coupling is negligible at these small amplitudes.

#### mu=-100, kappa=100, m=1.0 (TRUE VACUUM, WEAK EFFECTIVE COUPLING)

m_crit = 0.77. True vacuum, but the high saturation (kappa=100) weakens the effective coupling. Again omega = 1.002 and both states disperse:

| Metric | Control | 120-degree |
|--------|---------|------------|
| omega | 1.002 | 1.002 |
| E retained | 72.7% | 71.9% |
| Phase stability | N/A | **Perfectly stable** at 120 |

Nearly identical loss rates. The 120-degree pattern is stable but both states are dying.

#### mu=-30, kappa=20, m=1.0 (MARGINAL OSCILLON)

omega = 1.002 after equilibration -- right at mass gap. Not a true sub-gap oscillon:

| Metric | Control | 120-degree |
|--------|---------|------------|
| E retained | 96.4% | 89.0% |
| Phase stability | N/A | **Perfectly stable** at 120 |

---

## The Fundamental Tradeoff

There is a **catch-22** for the 120-degree oscillon:

1. **To survive**, the oscillon needs omega < m (sub-gap). This requires strong enough nonlinear frequency shift from the triple-product coupling.

2. **The frequency shift** scales as delta(omega^2) ~ |mu| * P_max^2. For the 120-degree state, P_max = f^3/4 (vs f^3 for symmetric), so the shift is **16x weaker**.

3. **To compensate**, you need 16x stronger coupling |mu|. But the false vacuum condition m^2 > |mu|*(2/kappa)^(2/3)/9 sets an upper limit on |mu|.

4. **At the maximum allowed |mu|** (~41 for kappa=20), the symmetric state barely achieves sub-gap oscillation (omega ~ 1.0). The 120-degree state, with 16x less frequency shift, cannot get below the gap.

5. **Increasing kappa** to allow larger |mu| simultaneously weakens the effective coupling (V saturates sooner), negating the benefit.

The 120-degree harmonic elimination trick works perfectly in theory (the 3rd harmonic channel P spectrum confirms it). But the same phase configuration that eliminates harmonics also kills the binding. You cannot have both.

---

## 3*omega < m/3 Threshold

The proposal hypothesized that 3*omega < m (i.e., omega < m/3 = 0.333) would make the 120-degree state perfectly stable. This requires extremely deep binding.

At the working oscillon (mu=-20): omega = 0.87, so 3*omega = 2.61 >> m = 1.0. The 3rd harmonic radiates freely.

At the false-vacuum parameters (mu=-60, -100): omega -> 0 (fields are static), so technically 3*omega < m, but the system is not an oscillon -- it is a static field in a false vacuum.

**No parameter regime exists where:**
- phi=0 is the true vacuum (necessary for oscillon to be meaningful)
- A genuine oscillating bound state forms with omega < m
- omega < m/3 is achieved

The maximum frequency shift achievable in the true-vacuum regime gives omega ~ 0.76 (at mu=-20, 1D). To reach omega < 0.33 would require a ~5x deeper shift, which is impossible without entering the false vacuum.

---

## Phase Stability Summary

| mu | kappa | Regime | omega | 120 phases stable? | 120 lives longer? |
|----|-------|--------|-------|--------------------|--------------------|
| -20 | 20 | **Oscillon** | 0.87 | No (drift, then disperse) | **NO** (17% vs 99%) |
| -30 | 20 | Marginal | 1.00 | Yes (perfectly) | No (89% vs 96%) |
| -40 | 20 | Near-critical | 1.00 | Yes (perfectly) | No (77% vs 87%) |
| -60 | 20 | False vacuum | 0 (static) | No (collapse to 0-phase) | N/A |
| -100 | 20 | False vacuum | 0 (static) | No (collapse to 0-phase) | N/A |
| -100 | 100 | Weak coupling | 1.00 | Yes (perfectly) | No (72% vs 73%) |

The 120-degree phases are paradoxically most stable when the binding is weakest (the coupling that would synchronize the phases to 0-degree is negligible). But in that regime, neither state forms a true oscillon.

---

## Files

| File | Description |
|------|-------------|
| `src/phase120_1d.c` | 1D solver: equilibrate + fork + compare |
| `data/phase120_ts_mu20.tsv` | 120-degree time series (mu=-20, base) |
| `data/phase120_control_ts_mu20.tsv` | 0-phase control time series (mu=-20) |
| `data/phase120_summary_mu20.tsv` | Summary table (mu=-20) |
| `data/phase120_*_mu30.tsv` | mu=-30 results |
| `data/phase120_*_mu40.tsv` | mu=-40 results |
| `data/phase120_*_mu60.tsv` | mu=-60 results (false vacuum) |
| `data/phase120_*_mu100.tsv` | mu=-100 results (false vacuum) |
| `data/phase120_*_mu100k100.tsv` | mu=-100, kappa=100 results |
| `data/phase120_ctrl_phi_spectrum.tsv` | Control phi_1 DFT (last run) |
| `data/phase120_ph120_phi_spectrum.tsv` | 120-deg phi_1 DFT (last run) |
| `data/phase120_ctrl_P_spectrum.tsv` | Control P(t) DFT (last run) |
| `data/phase120_ph120_P_spectrum.tsv` | 120-deg P(t) DFT (last run) |
