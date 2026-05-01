# Stochastic Calculus Applications to SCP

## 1. Floquet/Kapitsa Theory — Formalizing the Strobing

The V48 field-dependent theta mass m_theta^2 = lambda_theta * P^2 oscillates at
6*omega (carrier frequency omega ~ 1.5, so 6*omega ~ 9). P^2 crosses zero six
times per carrier cycle. Despite this, V48 simulations show the proton survives
T=5000 — the time-averaged mass is "sufficient."

Deterministic multi-scale averaging (Floquet/Kapitsa theory) makes this rigorous.
For the theta equation at a single point:

    d^2 theta / dt^2 = ... - lambda_theta * P^2(t) * theta + ...

Decompose P^2(t) = <P^2> + delta_P^2(t), where delta_P^2 is zero-mean and
oscillates at 6*omega. The averaging theorem (Khasminskii) says: on timescales
long compared to 1/(6*omega) ~ 0.11, the effective dynamics are governed by:

    d^2 theta / dT^2 = ... - lambda_theta * <P^2> * theta + (correction)

The correction to the effective mass squared is:

    delta_m_eff^2 = -<(delta_P^2)^2> / (2 * (6*omega)^2)

This shifts the effective Yukawa range by a few percent downward. The separation
of scales — fast strobing period ~0.7 vs slow Yukawa range ~7 — gives a 10:1
ratio, adequate for the averaging to be accurate.

**Practical use**: Compute <P^2> and var(P^2) from simulation data at the braid
core. The effective Yukawa range is 1/sqrt(lambda_theta * <P^2> + correction),
not 1/sqrt(lambda_theta * P_peak^2). This gives the correct nuclear binding
range without needing to resolve every carrier cycle.

## 2. Reduced SDE Models for Braid Collective Coordinates

Braids have thermal-like behavior (CONCEPT.md):
- Thermal equilibrium at A_noise ~ 0.02
- Steady radiation rate ~0.5 energy/time
- Temperature T > 0 always

A reduced stochastic model for braid position:

    dX = F_grav(X) dt + sigma(X) dW

where F_grav = -C * grad(rho) is the gravitational drift and sigma encodes
radiation-reaction noise. The diffusion coefficient sigma^2 is related to the
radiation rate by a fluctuation-dissipation relation.

Similarly, the breathing mode amplitude:

    dA = -k_eff * (A - A_eq) dt + sigma_A dW

This is an Ornstein-Uhlenbeck process. The equilibrium variance is
sigma_A^2 / (2 * k_eff), predicting the thermal fluctuation amplitude of
the breathing mode. The power spectrum is Lorentzian — testable against FFT
of simulation P_int time series.

**Practical use**: Fit sigma and k_eff from short simulation runs. Then predict
long-time behavior (drift rates, equilibration times, binding lifetimes) without
running expensive T=5000 PDE simulations. Standard Langevin dynamics technique.

## 3. Fokker-Planck for Binding Statistics

If braid position obeys the SDE above, the probability density p(x,t) satisfies:

    dp/dt = -d/dx[F(x) p] + (1/2) d^2/dx^2[sigma(x)^2 p]

At steady state:

    p_ss(x) ~ exp(-Phi(x) / D_eff)

where Phi is the drift potential and D_eff = sigma^2/2.

For two baryons at separation D with Yukawa potential V(D) from the
field-dependent theta mass:

    p_ss(D) ~ exp(-V_Yukawa(D) / T_eff)

where T_eff is set by the radiation temperature.

This predicts:
- Baryons preferentially sit at the Yukawa minimum (binding)
- At high temperature, baryons are unbound (thermal deconfinement)
- The binding/unbinding transition occurs at T_eff ~ V_binding
- The deuterium binding fraction as a function of lambda_theta and T_eff

**Practical use**: Given the Yukawa potential depth (measurable from E_pot in
the simulation) and the radiation temperature (measurable from theta_rms),
predict whether binding should occur WITHOUT running the two-baryon simulation.
This gives a go/no-go criterion for choosing lambda_theta before committing
GPU hours.

## 4. Ito vs Stratonovich and Noise-Induced Drift

The gravitational coupling is state-dependent: coupling is stronger in high-rho
environments (V33: kicked braid retains 1.40x momentum at A_bg=0.15 vs 0.30x
at A_bg=0.05). This means the noise amplitude sigma(X) depends on X — the
noise is multiplicative.

Stratonovich is the correct interpretation (Wong-Zakai theorem — noise from
continuous field dynamics, not discontinuous jumps). The Stratonovich SDE
dX = a_S dt + b(X) o dW is equivalent to the Ito SDE:

    dX = [a_S + (1/2) b * db/dX] dt + b(X) dW

The extra drift (1/2) b * db/dX is the noise-induced drift. Since sigma is
larger in high-rho (stronger coupling), d(sigma)/dX > 0 toward high rho,
so the noise-induced drift pushes the braid TOWARD LOW DENSITY — same
direction as gravity.

The magnitude is ~ (A_noise/A_bg)^2 ~ 4% of the gravitational drift. Small
but has the correct sign. Most significant for weakly bound systems near the
binding threshold.

**Practical use**: If measured gravitational drift exceeds the predicted
F = -C * grad(rho), the noise-induced drift may explain the discrepancy.
Worth computing once to check.

## 5. Phase Diffusion and Carrier Linewidth

The braid's carrier phase satisfies:

    dphi = omega dt + sigma_phi dW

Phase diffusion from radiation noise gives the carrier a finite linewidth.
If the phase diffusion is significant over one carrier cycle, the strobing
becomes genuinely stochastic (not periodic), and stochastic averaging
(rather than Floquet theory) is the correct framework.

**Practical use**: Measure the carrier phase autocorrelation from simulation
data. If it decorrelates within ~10 carrier cycles, stochastic averaging
applies. If it stays coherent for hundreds of cycles (likely), Floquet
theory is sufficient.

## Priority Ranking

1. **Fokker-Planck binding prediction** — highest immediate value. Can predict
   whether lambda_theta=10 at D=14 should bind before running the simulation.
   Requires: measure V_Yukawa(D) and T_eff from existing single-proton data.

2. **Floquet correction to Yukawa range** — formalizes the strobing result.
   Requires: measure <P^2> and var(P^2) at braid core from existing data.

3. **Reduced Langevin model** — enables fast parameter sweeps.
   Requires: fit drift and diffusion from short simulations.

4. **Noise-induced drift correction** — small effect, theoretical interest.
   Requires: careful comparison of predicted vs measured drift rates.
