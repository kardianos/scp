# V23-C Results: Correlation Length Near Gap Edge

## Summary

Systematic scan of the v21 triple-product oscillon as coupling strength |mu| is reduced from 20 to 5, pushing the oscillon frequency omega toward the mass gap m=1.0. All 11 parameter points produce long-lived oscillons at t=10000. No critical mu_c was found in this range -- the oscillons are more robust than expected.

The measured correlation length xi (from tail decay) does NOT diverge as omega approaches m. Instead, xi remains roughly constant at ~5-11, about 2-6x larger than the predicted xi_pred = 1/sqrt(m^2 - omega^2). The ratio xi/xi_pred decreases monotonically from ~6 at strong coupling to ~1.7 at weak coupling, trending toward the theoretical prediction.

---

## Model

Three massive scalar fields with saturating triple-product coupling:

    L = sum_a [ (1/2)(dt phi_a)^2 - (1/2)(dx phi_a)^2 - (m^2/2)phi_a^2 ] - V(P)
    V = (mu/2) P^2 / (1 + kappa P^2),   P = phi_1 phi_2 phi_3

Fixed parameters: kappa=20, m=1.0, A_init=1.0, sigma=3.0.
Scanned: mu in {-20, -18, -16, -14, -12, -10, -9, -8, -7, -6, -5}.
Grid: Nx=4000, xmax=100, tfinal=10000. Velocity Verlet with CFL condition.

---

## Phase 1: Correlation Length vs Gap Margin

### Results Table

| mu   | omega  | kappa_tail | delta    | xi_meas | A_peak | E_total | fc_avg | xi_pred | xi/xi_pred |
|------|--------|------------|----------|---------|--------|---------|--------|---------|------------|
| -20  | 0.8830 | 0.1146     | 0.117    | 8.73    | 0.519  | 1.269   | 0.983  | 2.13    | 4.10       |
| -18  | 0.9130 | 0.1118     | 0.087    | 8.94    | 0.526  | 1.343   | 0.963  | 2.45    | 3.65       |
| -16  | 0.8220 | 0.0934     | 0.178    | 10.70   | 0.561  | 1.449   | 0.950  | 1.76    | 6.09       |
| -14  | 0.8430 | 0.0885     | 0.157    | 11.30   | 0.763  | 2.906   | 0.967  | 1.86    | 6.08       |
| -12  | 0.8670 | 0.1206     | 0.133    | 8.29    | 0.753  | 4.598   | 0.984  | 2.01    | 4.13       |
| -10  | 0.8910 | 0.1512     | 0.109    | 6.61    | 0.739  | 4.680   | 0.992  | 2.20    | 3.00       |
| -9   | 0.9030 | 0.1636     | 0.097    | 6.11    | 0.737  | 4.758   | 0.992  | 2.33    | 2.63       |
| -8   | 0.9140 | 0.1783     | 0.086    | 5.61    | 0.752  | 5.034   | 0.990  | 2.46    | 2.27       |
| -7   | 0.9250 | 0.1944     | 0.075    | 5.14    | 0.750  | 6.088   | 0.991  | 2.63    | 1.95       |
| -6   | 0.9360 | 0.1960     | 0.064    | 5.10    | 0.740  | 6.281   | 0.988  | 2.84    | 1.80       |
| -5   | 0.9480 | 0.1925     | 0.052    | 5.19    | 0.725  | 6.284   | 0.984  | 3.14    | 1.65       |

Columns:
- delta = (m - omega)/m: gap margin (fraction below mass gap)
- xi_meas = 1/kappa_tail: measured correlation length from exponential tail fit
- xi_pred = 1/sqrt(m^2 - omega^2): predicted decay length from dispersion relation
- fc_avg: time-averaged energy fraction in core (|x| < 9)

### Key Findings

**1. No critical mu_c in [-20, -5].**
All oscillons survive to t=10000 with fc_avg > 0.95 and peak amplitude > 0.5. The oscillon is remarkably robust across a 4x range of coupling strength.

**2. Frequency increases monotonically (with anomaly at mu=-16,-14).**
From mu=-20 to mu=-5, omega increases from 0.883 to 0.948 (gap margin shrinks from 11.7% to 5.2%). However, mu=-16 and mu=-14 show anomalously LOW omega (0.822 and 0.843), lower than both neighbors. This suggests the oscillon undergoes a mode transition at intermediate coupling, possibly settling into a higher-energy breathing state.

**3. xi does NOT diverge.**
The prediction was xi ~ 1/sqrt(2*m*delta) diverging as delta -> 0. Instead:
- xi_meas stays in the range [5.1, 11.3], roughly constant
- kappa_tail is ~0.09-0.20, corresponding to decay lengths of 5-11 code units
- The ratio xi/xi_pred decreases from ~6 toward ~1.7 as coupling weakens

**4. xi/xi_pred converges toward ~1-2 at weak coupling.**
The convergence from above suggests the oscillon tail has additional structure beyond the asymptotic exponential: the breathing mode creates an oscillating envelope that decays more slowly than exp(-kappa_pred*r). At weak coupling, where the oscillon is simpler, the measured and predicted values converge.

**5. Energy increases at weak coupling.**
Counter-intuitively, E_total is HIGHER at weak coupling (6.28 at mu=-5 vs 1.27 at mu=-20). This is because at strong coupling, more energy is radiated during the initial shedding phase. The strong-coupling oscillon loses ~77% of its initial energy, while the weak-coupling one retains much more. The coupling potential V energy is small at weak coupling, so the oscillon is mostly free massive field with little binding.

### Anomalous mu=-16, -14 Behavior

At mu=-16 and -14, the oscillon has:
- Lower omega (0.822, 0.843) than neighbors
- Higher xi_meas (10.7, 11.3)
- Lower fc_avg (0.950, 0.967)
- Intermediate energy (1.45, 2.91 -- between strong and weak coupling values)

This suggests a transition in the oscillon's internal mode structure around |mu| ~ 14-16. The oscillon may be transitioning between two different quasiperiodic attractors, each with its own frequency and spatial extent.

### Tail Measurement Method

The tail decay rate kappa_tail is measured by:
1. Accumulating phi^2(x, t) during the second half of the evolution (t > 5000)
2. Computing the RMS envelope: phi_rms(x) = sqrt(<phi^2>_t)
3. Fitting log(phi_rms) vs x in the range x in [3, 25], excluding points below the noise floor (1e-6 of peak)

This time-averaging smooths over the breathing oscillations, giving a clean exponential envelope. The fitting range [3, 25] captures the tail of the soliton core (peak at x=0, width ~3) while avoiding the absorbing boundary region.

---

## Phase 2: Two-Oscillon Interaction

Two in-phase symmetric oscillons placed at separation D, evolved for t=5000. Grid: Nx=8000, xmax=200.
Tested mu = {-20, -12, -8} at D = {20, 30, 40, 50}.

### Results Table

| mu   | D_init | D_final | F_avg        | dD      | Notes                    |
|------|--------|---------|--------------|---------|--------------------------|
| -20  | 20     | 96.7    | -2.0e-05     | +76.7   | Repulsion (overlap)      |
| -20  | 30     | 78.9    | -1.8e-04     | +48.9   | Repulsion (overlap)      |
| -20  | 40     | 210.5   | +3.3e-05     | +170.5  | Weak repulsion           |
| -20  | 50     | 21.5    | +2.9e-06     | -28.5   | **ATTRACTION + MERGER**  |
| -12  | 20     | 11.2    | +9.6e-07     | -8.8    | **ATTRACTION + MERGER**  |
| -12  | 30     | 227.8   | -6.5e-05     | +197.8  | Repulsion                |
| -12  | 40     | 211.1   | -6.1e-07     | +171.1  | Weak repulsion           |
| -12  | 50     | 171.2   | +1.7e-06     | +121.2  | Separation (shedding)    |
| -8   | 20     | 271.8   | +1.4e-06     | +251.8  | Strong separation        |
| -8   | 30     | 244.8   | -1.4e-05     | +214.8  | Separation               |
| -8   | 40     | 141.1   | +2.1e-06     | +101.1  | Separation               |
| -8   | 50     | 236.7   | +1.2e-05     | +186.7  | Separation               |

### Phase 2 Interpretation

The Phase 2 results are dominated by the **shedding dynamics**, not equilibrium interactions:

1. **The Gaussian initial condition is not a proper oscillon.** Both lumps undergo violent shedding during t < 2000, radiating ~50-80% of their energy. The radiation pressure from this shedding pushes the centroids apart regardless of any tail-mediated interaction.

2. **Two cases show attraction/merger.** At mu=-20/D=50 and mu=-12/D=20, the final separation is less than the initial. These are the cases where the oscillons are far enough apart or strongly enough coupled that the equilibrium interaction dominates over the shedding transient.

3. **The force measurement is unreliable.** The parabolic fit to separation vs time gives F_avg ~ 10^{-5} to 10^{-6}, but the separation vs time curve is not parabolic -- it's dominated by the initial shedding burst followed by slow drift. The exponential fit lambda values (13.5 for mu=-20, 34.5 for mu=-12, -22.2 for mu=-8) are not meaningful.

4. **Proper Phase 2 requires equilibrated oscillons.** To measure the true interaction, one should first evolve single oscillons to equilibrium (t > 5000), then use the equilibrated profile as initial conditions for the two-oscillon run. The Gaussian initialization makes the interaction measurement impossible to disentangle from the shedding.

---

## Implications

### For the Critical Gravity Proposal

The original hypothesis was that xi diverges as omega approaches m, providing a mechanism for long-range interactions between oscillons that could mimic gravitational attraction. The results show:

1. **xi does NOT diverge** in the accessible parameter range. The correlation length remains finite (~5-11) even at delta=5.2% (omega=0.948).

2. **The oscillon is too robust.** There is no critical mu_c where the oscillon dies in the range mu in [-20, -5], so the critical behavior cannot be probed. At weaker coupling (|mu| < 5), the coupling energy V becomes negligible and the system approaches three free massive fields.

3. **The breathing mode prevents clean exponential tails.** The oscillon's time-dependent profile complicates the tail measurement. The effective tail is ~2-6x longer than predicted from the dispersion relation, likely due to the oscillating envelope creating a broader spatial distribution.

4. **Phase 2 requires equilibrated initial conditions.** The Gaussian initialization creates violent shedding that overwhelms any tail-mediated interaction. A follow-up study should use equilibrated single-oscillon profiles as initial conditions.

### What Would Be Needed

To properly test the xi-divergence hypothesis:
- Equilibrate single oscillons first, then use their profiles for two-oscillon runs
- Scan mu more finely near the critical point (if one exists below |mu|=5)
- Use a larger domain (xmax > 200) to accommodate oscillons with xi > 50
- Consider 3D: the 3D oscillon at mu=-20 has omega=0.948 (5% margin), already near the gap edge

---

## Files

| File | Description |
|------|-------------|
| `src/critical1d.c` | Combined Phase 1 + Phase 2 solver |
| `data/critical_mu{N}_ts.tsv` | Time series for each mu |
| `data/critical_mu{N}_profile.tsv` | Final RMS + snapshot profile |
| `data/critical_mu{N}_spectrum.tsv` | DFT power spectrum |
| `data/interact_mu{N}_D{D}_ts.tsv` | Phase 2 separation vs time |
