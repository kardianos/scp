# Avenue D: Spectral Self-Consistency and Zero-Point Constraints

## Summary

This analysis investigates whether the quantized soliton spectrum constrains the
degenerate sector coupling epsilon_0^2 (and thus the gravitational coupling g_top).

**Result: The quantum spectrum does NOT constrain epsilon_0^2.** At the gravity-scale
value epsilon_0^2 ~ 10^{-14}, all spectral effects are unmeasurably small compared
to nuclear observables. The hierarchy problem (why gravity is weak) cannot be resolved
by spectral self-consistency alone. It must come from elsewhere: anomaly cancellation,
topological quantization, or renormalization group flow.

Additionally, epsilon_0^2 > 0 **widens** the Channel C near-binding gap (wrong
direction), so the constraint coupling does not create new light particles.

**UPDATE (Avenue C result):** Avenue C has shown that epsilon_0^2 != 0 at the
algebraic level changes the target space from S^3 to S^7, killing solitons
(pi_3(S^7) = 0). Therefore epsilon_0^2 = 0 EXACTLY at the fundamental level.
The constraint eigenvalue analysis (Section 1.7) was computed hypothetically
at epsilon_0^2 > 0; at epsilon_0^2 = 0 the degenerate sector has NO coupling
to the bulk soliton. See Section 11 for implications.

---

## 1. Mode Enumeration

All known modes of the B=1 hedgehog Skyrmion at physical parameters
(e=1, rho_0=1, m_pi=0.398 code^{-1}, massive equilibrium profile):

### 1.1 Rotational modes (collective coordinates)

Quantum numbers: J = I from collective SU(2) quantization.

    E_rot = hbar^2 J(J+1) / (2 Lambda)

| State     | J=I | E_rot (code)  | E_rot (MeV)  |
|-----------|-----|---------------|--------------|
| Nucleon   | 1/2 | 6.5 x 10^{-6}| 6 x 10^{-5}  |
| Delta     | 3/2 | 2.2 x 10^{-5}| 2 x 10^{-4}  |

These are 3 rotational modes (generators of SU(2)_V).
Lambda = 86.6 code (massive profile), hbar = 0.0386 code.

Note: The Delta-N splitting = 0.2 MeV, far below experimental 293.7 MeV.
This is a known deficiency of the Skyrme model at our parameter choice.

### 1.2 K=0 Breathing modes (radial perturbation)

Sturm-Liouville: -(Pg')' + Wg = (omega^2/c^2) m g

| n | lambda  | omega  | E (MeV) | Status     |
|---|---------|--------|---------|-----------|
| 0 | 0.327   | 0.572  | 5.2     | Continuum  |
| 1 | 0.602   | 0.776  | 7.1     | Continuum  |
| 2 | 0.983   | 0.992  | 9.0     | Continuum  |
| 6 | 4.511   | 2.124  | 19.3    | Continuum  |

Threshold: lambda = m_pi^2 = 0.158. No bound states found below threshold.
No instabilities. The K=0 spectrum is purely continuum.
Mode n=6 at 19.3 MeV matches the soliton light-crossing frequency.

### 1.3 K=1 Channel A — Isorotational (L=0, I=1)

Three degenerate modes (isospin generators). Zero mode at omega = 0.

    Physical interpretation: the PION emerges from collective quantization
    of this zero mode. After semiclassical quantization, each mode acquires
    an effective frequency omega_eff = m_pi.

Contribution to quantum mass: 3 x (1/2) m_pi after Goldstone boson
treatment. In a proper quantization, these are Goldstone modes that
acquire mass through chiral symmetry breaking.

### 1.4 K=1 Channel B — Translational (L=1, I=0)

Three degenerate modes (spatial translations). Zero mode at omega = 0.

    Physical interpretation: center-of-mass motion.
    Eliminated in the rest frame. No contribution to M_quantum.

Lowest continuum eigenvalue: lambda = 0.399 (omega = 0.632, 5.7 MeV),
well above threshold m_pi^2 = 0.158.

### 1.5 K=1 Channel C — Mixed (L=1, I=1, K=1)

Nine degenerate modes (3 orbital x 3 isospin).

    Potential minimum: W/m = 0.167 at the box boundary
    Continuum threshold: m_pi^2 = 0.158
    Gap: 0.009 (5.7% above threshold)

**No bound states found.** The W/m ratio approaches m_pi^2 from above
asymptotically as r -> infinity but never reaches it.

Lowest box eigenvalue: lambda = 0.250 (omega = 0.500, 4.5 MeV).

If a bound state existed here, it would have quantum numbers:
- K=1, L=1, I=1, J^P = 1^+
- Analogous to the a_1(1260) axial-vector meson in QCD
- Mass would be near m_pi (but below it) — lighter than the pion!

### 1.6 K >= 2 modes

Not computed numerically. These include:
- K=2: five degenerate modes (d-wave), expected W/m well above threshold
- Higher K: progressively weaker coupling to soliton

At high K, modes approach free-particle behavior (W/m ~ K(K+1)/r^2 -> 0).
The density of states approaches the free-space continuum.

### 1.7 Constraint (degenerate sector) modes

From the eigenvalue problem: -nabla^2 p - S_0(r) p = E_0 p

The Skyrmion strain S_0(r) = f'^2 + 2 sin^2(f)/r^2 acts as an attractive
potential for the degenerate sector fields (p, j_1, j_2, j_3).

| l  | E_0           | kappa  | Range (fm) | Degeneracy      | Status |
|----|---------------|--------|------------|-----------------|--------|
| 0  | -1.979        | 1.407  | 0.400      | 4 (p + 3 j)     | Bound  |
| 1  | > 0           | —      | —          | 4 x 3 = 12      | Unbound|
| 2  | > 0           | —      | —          | 4 x 5 = 20      | Unbound|

Only the s-wave supports a bound state. There are 4 degenerate bound modes
(one for each component of the degenerate sector: p, j_1, j_2, j_3).

These modes have omega = sqrt(|E_0|) = 1.407 code (12.8 MeV).
This is **above** the pion mass threshold (m_pi = 0.398 code = 3.6 MeV),
so the constraint mode lies in the pion continuum and can decay.

The massive profile gives |E_0| = 1.979, vs 1.487 for the sigma-model profile.
The massive soliton is more compact, deepening the S_0(r) well.

---

## 2. Quantum Mass Formula

The quantum mass of the B=1 soliton at one-loop order:

    M_quantum = M_classical + E_rot + E_zp(bulk) + E_zp(degen) + counterterms

### 2.1 Classical contribution

    M_classical = E_sol = 110.19 code = 1003 MeV

(At the massive profile; parameter-fit adjusted to match M_nucleon = 938 MeV
by tuning rho_0 and e.)

### 2.2 Rotational contribution

    E_rot(J=1/2) = 3 hbar^2 / (8 Lambda) = 6.5 x 10^{-6} code = 0.06 keV

Completely negligible at our parameters.

### 2.3 Bulk zero-point energy (UV divergent)

    E_zp(bulk) = (1/2) hbar Sum_n omega_n

Sum over all bulk modes (K=0,1,2,...). This sum diverges as Lambda_UV^4 in 4D.

Renormalization: subtract the free-space zero-point energy.

    E_zp^{ren}(bulk) = (1/2) hbar Sum_n [omega_n - omega_n^{free}]

This is finite and O(1) in code units, modifying M_quantum by ~10 MeV.
It is an important correction but is **independent of epsilon_0^2**.

### 2.4 Degenerate sector zero-point energy

At epsilon_0^2 > 0, the degenerate sector acquires dynamics through the
constraint |q|^2 + epsilon_0^2 p^2 = rho_0^2.

**Bound modes (l=0):**

    E_zp(bound) = 4 x (1/2) hbar sqrt(|E_0|) x epsilon_0^2
                = 4 x 0.5 x 0.0386 x 1.407 x epsilon_0^2
                = 0.109 x epsilon_0^2 code
                = 0.99 x epsilon_0^2 MeV

For epsilon_0^2 ~ 10^{-14}: E_zp(bound) ~ 10^{-14} MeV = **NEGLIGIBLE**.

**Continuum modes (l >= 1):**

Same UV structure as bulk (quartic divergence). After renormalization,
the epsilon_0^2-dependent finite part comes from the difference between
the constraint-modified and free continuum spectra. This is also O(epsilon_0^2).

### 2.5 Net quantum mass

    M_quantum = M_classical + O(keV) + E_zp^{ren}(bulk) + O(epsilon_0^2 MeV)

The epsilon_0^2-dependent term is negligible for epsilon_0^2 ~ 10^{-14}.
The quantum mass is effectively independent of the gravitational coupling.

---

## 3. Parameter Over-Determination Analysis

### 3.1 Parameters

| # | Parameter  | Role                          |
|---|-----------|-------------------------------|
| 1 | rho_0     | Field amplitude (VEV)         |
| 2 | lambda    | Bulk potential coupling        |
| 3 | e         | Skyrme coupling               |
| 4 | mu        | Degenerate mass scale          |
| 5 | c         | Speed of light                 |
| 6 | epsilon_0^2 | Constraint coupling (NEW)    |

### 3.2 Observables

| # | Observable   | Value           | Constrains        |
|---|-------------|-----------------|-------------------|
| 1 | M_nucleon   | 938.3 MeV       | rho_0^3/e         |
| 2 | r_proton    | 0.841 fm        | rho_0/e           |
| 3 | m_pion      | 139.6 MeV       | mu (pion mass)     |
| 4 | Delta-N     | 293.7 MeV       | Lambda ~ e^{-3}   |
| 5 | g_A         | 1.27            | Profile shape      |
| 6 | F_pi        | 92.1 MeV        | rho_0              |
| 7 | G_Newton    | 6.674e-11       | epsilon_0^2 (NEW)  |

### 3.3 Counting

Without epsilon_0^2 (5 parameters, >= 6 observables):
- M_N + r_p fix the combination rho_0^4/e^2 = 0.0259 and rho_0^3/e = 103.13/E_{FB_coeff}
- m_pi fixes mu via m_pi = mu F_pi
- Delta-N, g_A, F_pi are PREDICTIONS of the remaining freedom (e)
- The system is overconstrained: current best fit gives e_std = 3.69, F_pi = 95 MeV
  (3% from experiment), but Delta-N = 0.2 MeV (1000x too low)

With epsilon_0^2 (6 parameters, >= 7 observables):
- G_Newton adds one observable, epsilon_0^2 adds one parameter
- Net: still overconstrained (same degree)
- But epsilon_0^2 only affects G_Newton. All other observables are
  insensitive to epsilon_0^2 at the 10^{-14} level.

**Result: epsilon_0^2 effectively decouples from the nuclear spectrum.**
It is a free parameter determined solely by G_Newton. The nuclear
observables neither constrain it nor are affected by it.

---

## 4. epsilon_0^2 Effect on Channel C

### 4.1 Mechanism

At epsilon_0^2 > 0, the constraint modifies |q|:

    |q(r)|^2 = rho_0^2 - epsilon_0^2 p_0^2(r)

where p_0(r) is the constraint eigenfunction. This changes the effective
radial density profile rho(r) -> sqrt(rho_0^2 - epsilon_0^2 p_0^2).

### 4.2 Effect on Channel C potential

Channel C has W = W_cent + W_E4 + W_pi, with P = m (both proportional to Omega).

The W/m ratio decomposes as:

At the W/m minimum (r ~ 10 code):

| Term      | W_i/m      | rho-dependence |
|-----------|-----------|----------------|
| W_cent/m  | 0.020     | Independent    |
| W_E4/m    | ~0        | Independent    |
| W_pi/m    | 0.158     | ~ 1/rho^2     |
| **Total** | **0.178** |                |

Threshold: m_pi^2 = 0.158.

The rho-dependent correction:

    delta(W_pi/m) = W_pi/m x epsilon_0^2 p_0^2(r) / rho_0^2

This is **positive** because:
- W_pi/m > 0 at the minimum (cos(f) > 0 at large r)
- p_0^2 > 0

**Sign: epsilon_0^2 > 0 INCREASES W/m, WIDENING the gap.**

The constraint coupling pushes Channel C **further from binding**, not
closer. No new light particles are created.

### 4.3 Quantitative estimate

The eigenfunction p_0(r) is concentrated near the soliton core (range 0.4 fm)
while the Channel C minimum is at r ~ 5-10 code (3-6 fm). At the minimum
location, p_0^2 ~ 10^{-23} (exponentially suppressed). The shift is:

    delta(W/m) ~ 10^{-23} x epsilon_0^2

For any conceivable epsilon_0^2, this is unmeasurably small.

---

## 5. Breathing Mode (K=0) Interaction with Constraint

The constraint bound state has omega = sqrt(|E_0|) = 1.407 code.
The breathing mode continuum starts at omega = m_pi = 0.398 code.

Since 1.407 > 0.398, the constraint mode lies **inside the continuum**.
It is a resonance, not a stable bound state, and can decay:

    constraint mode -> pion(s) + radiation

The decay width is O(epsilon_0^2), so at epsilon_0^2 ~ 10^{-14}, the
resonance is astronomically narrow (effectively stable on nuclear timescales,
but not on cosmological timescales).

The constraint mode's zero-point energy contributes to M_quantum regardless
of whether it is technically stable or resonant.

For m_pi > sqrt(|E_0|) = 1.407 code, the constraint mode would be below the
continuum threshold and truly stable. This would require m_pi > 12.8 MeV
in physical units — about 11x the actual pion mass. Not physical.

---

## 6. UV Divergence Structure

### 6.1 Bulk sector (epsilon_0^2 = 0)

Zero-point sum: E_zp = (1/2) hbar Sum omega_n

Divergence: quartic in 4D (Sum ~ Lambda_UV^4).
Renormalization: subtract free-space contribution.
Renormalized result: finite, O(1) in code units.

### 6.2 Degenerate sector (epsilon_0^2 > 0)

Same UV structure as bulk (4 scalar fields with S_0(r) potential).
The epsilon_0^2-dependent divergence has the same form as the bulk divergence.

Counterterms: polynomial in epsilon_0^2, m_pi^2, Lambda_UV.

    L_ct = c_1 Lambda_UV^4 + c_2 epsilon_0^2 Lambda_UV^2 + c_3 epsilon_0^4 log(Lambda_UV) + ...

All divergences can be absorbed by the existing parameters.
**No constraint on epsilon_0^2 from UV finiteness alone.**

### 6.3 Zeta function regularization

    zeta_S(s) = Sum_n omega_n^{-s}

    E_zp^{ren} = (1/2) hbar zeta_S(-1)

The bound state contribution (4 modes at omega = 1.407):

    zeta_{bound}(-1) = 4 x 1.407 = 5.627

This is always finite (isolated eigenvalue). The continuum requires
analytic continuation but introduces no epsilon_0^2 constraint.

---

## 7. Self-Consistency Tests

### Test 1: UV finiteness
The zero-point sum diverges regardless of epsilon_0^2.
Standard renormalization absorbs the divergence.
**No constraint on epsilon_0^2.**

### Test 2: Nucleon mass
M_quantum = 938 MeV requires balancing M_classical with quantum corrections.
The epsilon_0^2 contribution is ~ epsilon_0^2 x 1 MeV ~ 10^{-14} MeV.
**No constraint on epsilon_0^2** (completely negligible).

### Test 3: Delta-N splitting
Depends on Lambda_eff = Lambda + O(epsilon_0^2).
At epsilon_0^2 ~ 10^{-14}: shift ~ 10^{-14}.
**No constraint on epsilon_0^2.**

### Test 4: Channel C binding
Gap = 0.009 (5.7% of threshold).
epsilon_0^2 > 0 widens the gap (wrong sign).
**epsilon_0^2 cannot create new light particles** in Channel C.

### Test 5: Spectral flow
As epsilon_0^2 increases from 0, degenerate modes turn on continuously.
No level crossings or topological transitions at infinitesimal epsilon_0^2.
**No discrete constraint from spectral flow.**

### Test 6: Bohr-Sommerfeld quantization
The constraint eigenvalue E_0 = -1.979 is deep below threshold.
The Bohr-Sommerfeld condition is satisfied for any epsilon_0^2 > 0.
**No discrete constraint from quantization condition.**

---

## 8. What Would Constrain epsilon_0^2?

The spectral analysis demonstrates that epsilon_0^2 enters all nuclear
observables at order epsilon_0^2 ~ 10^{-14}, which is negligible.
The quantum spectrum is essentially identical for epsilon_0^2 = 0 and
epsilon_0^2 = 10^{-14}.

Potential constraints from **non-spectral** mechanisms:

1. **Anomaly cancellation (Avenue C)**: Algebraic constraint from
   triangle diagrams involving the topological current B^mu and
   degenerate sector. This is the most promising route.

2. **WZW coefficient (Avenue B)**: Topologically quantized coupling
   from the rotating soliton. Gives g_top ~ N/(240 pi^2) x hbar/Lambda,
   but naive estimate is ~10^{10}x too large.

3. **Technical naturalness**: The Z_2 symmetry p -> -p protects
   epsilon_0^2 = 0. Any nonzero value is technically natural as long as
   it is small (no quadratic divergence drives it to O(1) values).

4. **Asymptotic safety**: The RG flow of epsilon_0^2 may have a
   non-trivial fixed point. If epsilon_0^2 is an irrelevant coupling
   near an asymptotically safe fixed point, its value is predicted
   by the UV completion.

5. **Nonperturbative effects (instantons)**: Tunneling between
   topological sectors could generate epsilon_0^2 ~ exp(-S_instanton),
   where S_instanton is the instanton action. If S_instanton ~ 30,
   then epsilon_0^2 ~ exp(-30) ~ 10^{-13}, close to the required value.

---

## 9. Numerical Code and Verification

### Tools used

- `src/angular_modes.c`: K=1 mode solver (Channels A, B, C)
- `src/normal_modes.c`: K=0 breathing mode solver
- `src/varsig.c`: Constraint eigenvalue solver
- `src/spectral_eps.c`: **New** — combined spectral analysis for Avenue D

### Key numerical results (massive profile, e=1, rho_0=1, m_pi=0.398)

| Quantity                    | Value           |
|-----------------------------|-----------------|
| E_sol (massive)             | 110.19 code     |
| Lambda (moment of inertia)  | 86.6 code       |
| Constraint E_0 (l=0)        | -1.979          |
| Constraint kappa            | 1.407           |
| Constraint range            | 0.400 fm        |
| Channel C W/m minimum       | 0.167           |
| Channel C threshold         | 0.158 (m_pi^2)  |
| Channel C gap               | 5.7%            |
| hbar                        | 0.0386 code     |
| 1 code E                    | 9.098 MeV       |
| ZPE of constraint modes     | 0.109 eps_0^2 code = 0.99 eps_0^2 MeV |

### Verification

All solvers have been cross-checked:
- Virial theorem satisfied to machine precision (E_2 - E_4 + 3E_V = 0)
- Normal modes: variational consistency check (eigenvalue vs Rayleigh quotient)
- Constraint eigenvalue: ODE residual < 4% (Numerov method)
- Angular modes: zero modes at omega=0 confirmed for Channels A and B

---

## 10. Conclusions

### Primary finding

The quantum soliton spectrum does not constrain epsilon_0^2. At the
gravitational scale epsilon_0^2 ~ 10^{-14}, the degenerate sector's
contribution to all nuclear observables is negligible (O(10^{-14}) corrections).

### Channel C near-binding

The near-binding in Channel C (gap = 5.7%) is a genuine feature of the
Skyrme model at physical parameters. However:
- epsilon_0^2 > 0 widens the gap (wrong direction for binding)
- The eigenfunction p_0(r) is exponentially suppressed at the gap location
- No new light particles emerge from the constraint coupling

### Over-determination

The system has 6 parameters and 7+ observables, so it is formally overconstrained.
But epsilon_0^2 only affects one observable (G_Newton) and is invisible to the
other six. In practice, the system splits into:
- Nuclear sector: 5 parameters, 6+ observables (overconstrained, same as before)
- Gravitational sector: 1 parameter (epsilon_0^2), 1 observable (G_Newton) (just-determined)

### Assessment

Avenue D (spectral self-consistency) provides a **null result**: no constraint
on the gravitational coupling from the spectrum. This is the expected outcome
for a hierarchy of 10^{38} between nuclear and gravitational scales — the two
sectors simply do not talk to each other at the spectral level.

**Reinforced by Avenue C:** The topological obstruction pi_3(S^7) = 0 forces
epsilon_0^2 = 0 exactly. The degenerate sector is completely inert at the
algebraic level. Any gravitational coupling must arise dynamically — through
the WZW term for rotating solitons (Avenue B), quantum condensate formation,
or nonperturbative effects. The spectral analysis confirms this is fully
self-consistent: at epsilon_0^2 = 0, the bulk spectrum stands alone.

---

## 11. Impact of Avenue C: epsilon_0^2 = 0 Exactly

### 11.1 The topological obstruction

Avenue C has established that epsilon_0^2 != 0 at the fundamental algebraic level
changes the norm constraint from |q|^2 = rho_0^2 (target S^3) to
|q|^2 + epsilon_0^2 |d|^2 = rho_0^2 (target S^7). Since pi_3(S^7) = 0,
topological solitons with B != 0 cease to exist on S^7. Therefore:

    epsilon_0^2 = 0 EXACTLY at the algebraic level.

This is a hard topological constraint, not a perturbative approximation.

### 11.2 Consequences for the spectral analysis

**Section 1.7 (constraint modes):** The eigenvalue problem
-nabla^2 p - S_0(r) p = E_0 p was computed under the assumption that
epsilon_0^2 > 0 provides a coupling between the bulk and degenerate sectors.
At epsilon_0^2 = 0, this coupling does not exist. The degenerate sector
fields (p, j_1, j_2, j_3) are completely decoupled from the soliton.
The constraint eigenvalue E_0 = -1.979 exists as a mathematical property
of the S_0(r) potential, but it does not correspond to a physical mode
because there is no mechanism to populate it.

**Sections 2.4, 3, 4, 5 (epsilon_0^2-dependent corrections):** All corrections
proportional to epsilon_0^2 are exactly zero. The quantum mass formula simplifies:

    M_quantum = M_classical + E_rot + E_zp(bulk) + counterterms

with no degenerate sector contribution at all.

**Section 4 (Channel C):** The constraint coupling mechanism for shifting
Channel C is inoperative. However, the near-binding gap (5.7%) remains
a genuine property of the bulk spectrum and could potentially be closed
by other mechanisms:
- Higher-order Skyrme terms (L_6 sextic)
- Finite-lambda effects on the profile
- Non-hedgehog deformations of the soliton

### 11.3 What the spectral analysis DOES tell us at epsilon_0^2 = 0

Even with epsilon_0^2 = 0, the mode enumeration (Sections 1.1-1.6) remains
fully valid and important:

1. **Bulk spectrum is complete at one-loop:** K=0 breathing (no bound states),
   K=1 Channels A/B/C (no bound states, near-binding in C), rotational modes.

2. **The quantum mass is determined by bulk modes only:** No degenerate
   sector participation. The 5 parameters (rho_0, lambda, e, mu, c) must
   match all nuclear observables without help from epsilon_0^2.

3. **G_Newton requires a DYNAMICAL mechanism:** Since epsilon_0^2 = 0
   algebraically, gravity cannot come from the constraint coupling. The
   remaining candidates are:
   - **WZW term (Avenue B):** The rotating soliton sources the degenerate
     sector through the Wess-Zumino-Witten term. This does NOT require
     epsilon_0^2 != 0 — it couples B^i (spatial baryon current from rotation)
     directly to degenerate sector fields through a topologically quantized
     coefficient.
   - **Dynamical condensate:** Quantum fluctuations of p around zero could
     generate an effective epsilon_0^2 through <p^2> != 0, even though the
     tree-level coupling is zero. This would be a radiatively generated
     coupling. The question is whether the soliton background enhances
     these fluctuations enough to produce a measurable effect.

4. **The Channel C gap is a prediction:** At epsilon_0^2 = 0, the gap of
   5.7% above threshold is a pure prediction of the Skyrme model. If
   a future calculation (with L_6, finite lambda, or modified profile)
   closes this gap, it would predict a new light particle with J^P = 1^+,
   I = 1 — but this is independent of the gravitational coupling.

### 11.4 Revised conclusion

The original Avenue D conclusion (spectral self-consistency does not constrain
epsilon_0^2) is STRENGTHENED by the Avenue C result: epsilon_0^2 is not merely
unconstrained by the spectrum — it is forced to zero by topology. The spectrum
confirms this is self-consistent: at epsilon_0^2 = 0, the bulk and degenerate
sectors decouple completely, and the nuclear spectrum depends only on the 5
bulk parameters.

The gravitational coupling g_top must arise from a dynamical mechanism
(WZW for rotating solitons, or quantum condensate) rather than from a
fundamental algebraic coupling. The spectral analysis provides no constraint
on the magnitude of such dynamical effects.

---

## 12. Impact of Avenue B: WZW Overcoupling and Cancellation

### 12.1 The WZW result

Avenue B has established:
- The B^i spatial coupling vanishes by baryon current conservation (div B = 0).
  This is exact at O(Omega), not a symmetry accident.
- WZW reduces to the STATIC B^0 coupling — same mechanism as Path 3.
  Rotation adds nothing new.
- For SU(2): pi_5(SU(2)) = Z_2, so the WZW coefficient N is either 0 or 1.
  - N=0: no WZW coupling to degenerate sector. Gravity is zero at tree level.
  - N=1: g_WZW = 1/(240 pi^2) = 4.2 x 10^{-4}, giving G_eff ~ 10^{26} G_Newton.
    This is 10^{13} times TOO STRONG.

### 12.2 Can spectral cancellation rescue N=1?

If N=1, the effective gravitational coupling would be:

    g_eff = g_WZW + g_loop + g_higher + ...
          = 4.2 x 10^{-4} + (loop corrections) + ...

Matching G_Newton requires g_eff = 2.8 x 10^{-17}, demanding cancellation of
g_WZW to 1 part in 10^{13}. The question is whether loop corrections from the
bulk spectrum can provide this cancellation.

**Answer: NO, for three reasons.**

**Reason 1: Topological quantization is exact.**
The WZW coefficient N is an integer (topological winding number of
pi_5(SU(2))). It cannot receive perturbative corrections. No matter how
many loop diagrams are summed, N remains exactly 0 or 1. The coefficient
g_WZW = N/(240 pi^2) is exact to all orders in perturbation theory.

**Reason 2: At epsilon_0^2 = 0, loop corrections to the B^0-p vertex are zero.**
The one-loop correction to the B^0-p coupling involves a pion loop with both
a B^0 insertion (from the topological current) and a p external leg (from the
degenerate sector). At epsilon_0^2 = 0, there is no pion-p vertex in the
Lagrangian, so the loop diagram has no way to connect the pion propagator
to the external p leg. The correction is identically zero.

More precisely: the one-loop effective action for the B^0-p vertex is:

    Gamma^{(1)}_{B0,p} = Sum_n integral d^4k / (k^2 + omega_n^2) x V_{pion-p}(k,n)

where V_{pion-p} is the pion-p coupling vertex. At epsilon_0^2 = 0, V_{pion-p} = 0
for all modes n, so the entire sum vanishes regardless of the mode spectrum
{omega_n}. The mode enumeration (breathing modes, Channel A/B/C) does not
enter the calculation at all.

**Reason 3: The WZW term operates at a different level from loop corrections.**
The WZW term is a 5-dimensional topological term (Chern-Simons form integrated
over the 5D ball whose boundary is spacetime x field space). Loop corrections
are 4-dimensional momentum-space integrals. They cannot cancel each other
because they have different analytical structures. The WZW contribution to the
effective action is non-analytic in the coupling constants (it involves epsilon
tensors, not polynomial loop integrals).

### 12.3 Spectral implications at N=0 vs N=1

**If N=0 (no WZW coupling):**
- g_eff = 0 at all orders in perturbation theory (epsilon_0^2 = 0 kills
  the constraint coupling, and N=0 kills the WZW coupling).
- Gravity would require a nonperturbative mechanism: instanton tunneling,
  condensate formation, or physics beyond Cl^+(3,0,1).
- The bulk spectrum is completely unaffected. All mode enumeration results
  (Sections 1.1-1.6) stand unchanged.

**If N=1 (WZW too strong by 10^{13}):**
- g_eff = 4.2 x 10^{-4} with no perturbative cancellation available.
- This is a FATAL problem: nucleon-nucleon "gravity" would be 10^{13}x
  stronger than observed, comparable to the Skyrme-term nuclear force.
- The spectrum cannot fix this. No choice of (rho_0, lambda, e, mu, c)
  reduces g_WZW, since it is a topological constant.
- The only escape routes:
  1. Massive degenerate sector: if mu >> m_pi, the p-mediated force
     becomes Yukawa with range 1/mu. At mu ~ 10^7 code^{-1}
     (range ~ 10^{-7} fm), the effective "gravity" would be invisible
     at nuclear and larger scales. But this introduces an extreme
     hierarchy in the mass parameter.
  2. Non-minimal kinetic term for p: if the degenerate sector has
     a kinetic coefficient kappa^2 >> 1, then G_eff ~ g_WZW^2/kappa^2
     is suppressed. Need kappa ~ 10^{6.5} to get the right G_Newton.
  3. The theory is simply inconsistent with N=1, favoring N=0.

### 12.4 The N=0 scenario and residual gravitational mechanisms

If N=0 (the self-consistent scenario), then at the perturbative level:
- epsilon_0^2 = 0 (Avenue C: topological obstruction)
- g_WZW = 0 (Avenue B: N=0)
- g_loop = 0 (no vertex to correct at epsilon_0^2 = 0)

**Gravity is zero at all orders in perturbation theory.**

This is a sharp prediction: if CHPT describes nature, the gravitational
coupling must be NONPERTURBATIVE. Possible mechanisms:

1. **Instanton tunneling:** Euclidean configurations connecting B=0 and B=1
   sectors generate an effective B^0-p coupling of order e^{-S_inst}.
   If S_inst ~ 32, then g_top ~ e^{-32} ~ 10^{-14}, in the right ballpark.
   The instanton action depends on the bulk parameters (rho_0, e, lambda),
   so this WOULD provide a connection between the nuclear spectrum and
   the gravitational coupling — but only through nonperturbative physics
   that is invisible to the one-loop mode analysis.

2. **Condensate formation:** A <p^2> condensate could form through
   nonperturbative dynamics, generating an effective epsilon_0^2_eff.
   This is analogous to chiral symmetry breaking in QCD, where <qq> != 0
   despite the perturbative vacuum having <qq> = 0.

3. **Higher topological sectors:** The gravitational coupling could
   involve tunneling between soliton sectors (B -> B+1 -> B) that
   generates an effective long-range force. The amplitude for such
   processes is ~ e^{-E_sol/hbar} ~ e^{-103/0.039} ~ e^{-2600},
   which is astronomically small — far too small for gravity.

**None of these mechanisms are accessible through the one-loop mode spectrum.**
The spectral analysis (breathing modes, angular modes, constraint eigenvalues)
probes the perturbative structure around a single soliton. Nonperturbative
effects require information about the FULL field-theory path integral,
which is beyond the scope of normal-mode analysis.

### 12.5 Updated conclusion

The combined results of Avenues B, C, and D paint a consistent picture:

1. **epsilon_0^2 = 0 exactly** (Avenue C: topological obstruction)
2. **N = 0 or 1 for WZW** (Avenue B: pi_5(SU(2)) = Z_2)
3. **N=1 is overcoupled by 10^{13}** (Avenue B: g_WZW = 4.2 x 10^{-4})
4. **No perturbative cancellation exists** (Avenue D: spectral analysis)
5. **Spectrum is blind to epsilon_0^2** (Avenue D: O(10^{-14}) effects)

The self-consistent scenario is **N=0, epsilon_0^2 = 0**: gravity is
absent at all orders in perturbation theory, and must arise from
nonperturbative dynamics (instantons or condensates).

The perturbative spectrum — which is completely characterized by the
mode enumeration in Sections 1.1-1.6 — is necessary for nuclear physics
(masses, splittings, scattering) but irrelevant for gravity.

---

## Files

- `src/spectral_eps.c` — Combined spectral analysis code for Avenue D
- `src/angular_modes.c` — K=1 angular mode solver
- `src/normal_modes.c` — K=0 breathing mode solver
- `src/varsig.c` — Constraint eigenvalue solver
- `data/profiles/profile_massive_e1_mpi0.398.dat` — Massive equilibrium profile
