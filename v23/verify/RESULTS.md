# V23 Verification Results

Independent verification of V23-A (Hessian), V23-B (Cross-gradient), V23-C (Critical).

---

## V23-A: Hessian Splitting

### 1. Hessian Derivation from Scratch

The potential is V(P) = (mu/2) P^2 / (1 + kappa P^2), with P = phi_1 phi_2 phi_3.

**First derivative:**

    dV/dP = mu P (1 + kappa P^2)^{-1} - mu kappa P^3 (1 + kappa P^2)^{-2}
          = mu P / (1 + kappa P^2)^2

    dV/dphi_a = (dV/dP)(dP/dphi_a) = mu P (dP/dphi_a) / (1 + kappa P^2)^2

**Second derivative:**

Let Q = 1 + kappa P^2. Then dV/dphi_a = mu P (dP/dphi_a) / Q^2.

    d^2V/(dphi_a dphi_b) = mu [ (dP/dphi_a)(dP/dphi_b) Q^2 + P (d^2P/(dphi_a dphi_b)) Q^2
                                - P (dP/dphi_a) 2Q (2 kappa P (dP/dphi_b)) ] / Q^4

Simplify by collecting the Q factors:

    = mu / Q^4 * [ (dP/dphi_a)(dP/dphi_b) Q^2 + P (d^2P/(dphi_a dphi_b)) Q^2
                    - 4 kappa P^2 (dP/dphi_a)(dP/dphi_b) Q ]

    = mu / Q^3 * [ (dP/dphi_a)(dP/dphi_b)(Q - 4 kappa P^2) + P (d^2P/(dphi_a dphi_b)) Q ]

    = mu / Q^3 * [ (dP/dphi_a)(dP/dphi_b)(1 + kappa P^2 - 4 kappa P^2) + P (d^2P/(dphi_a dphi_b))(1 + kappa P^2) ]

    = mu / Q^3 * [ (dP/dphi_a)(dP/dphi_b)(1 - 3 kappa P^2) + P (d^2P/(dphi_a dphi_b))(1 + kappa P^2) ]

**CHECK: This matches the PROPOSAL.md formula (line 30). CORRECT.**

### 2. Evaluation on the Symmetric Background phi_a = f

With phi_1 = phi_2 = phi_3 = f:

    P = f^3
    dP/dphi_a = phi_b phi_c = f^2   (for all a)
    d^2P/(dphi_a dphi_b) = phi_c = f   (for a != b)
    d^2P/(dphi_a dphi_a) = 0           (for a = b)

Substituting into the second derivative formula with Q = 1 + kappa f^6:

**Diagonal (a = b):**

    H_{aa} = mu / Q^3 * [ f^4 (1 - 3 kappa f^6) + 0 ]
           = mu f^4 (1 - 3 kappa f^6) / Q^3

**Off-diagonal (a != b):**

    H_{ab} = mu / Q^3 * [ f^4 (1 - 3 kappa f^6) + f^3 * f * Q ]
           = mu / Q^3 * [ f^4 (1 - 3 kappa f^6) + f^4 (1 + kappa f^6) ]
           = mu f^4 / Q^3 * [ (1 - 3 kappa f^6) + (1 + kappa f^6) ]
           = mu f^4 (2 - 2 kappa f^6) / Q^3
           = 2 mu f^4 (1 - kappa f^6) / Q^3

**CHECK: Matches PROPOSAL.md lines 54-57. CORRECT.**

### 3. Eigenvalue Verification

The Hessian is a 3x3 matrix with diagonal d and off-diagonal o:

    d = mu f^4 (1 - 3 kappa f^6) / (1 + kappa f^6)^3
    o = 2 mu f^4 (1 - kappa f^6) / (1 + kappa f^6)^3

For a matrix with equal diagonals d and equal off-diagonals o, the eigenvalues are:
- lambda_sym = d + 2o (eigenvector (1,1,1)/sqrt(3), multiplicity 1)
- lambda_anti = d - o (eigenvector perpendicular to (1,1,1), multiplicity 2)

**Symmetric eigenvalue:**

    lambda_sym = d + 2o
               = mu f^4 / Q^3 * [ (1 - 3kf^6) + 4(1 - kf^6) ]
               = mu f^4 (5 - 7 kappa f^6) / Q^3

**CHECK: Matches PROPOSAL.md line 76. CORRECT.**

**Antisymmetric eigenvalue:**

    lambda_anti = d - o
                = mu f^4 / Q^3 * [ (1 - 3kf^6) - 2(1 - kf^6) ]
                = mu f^4 / Q^3 * [ 1 - 3kf^6 - 2 + 2kf^6 ]
                = mu f^4 (-1 - kf^6) / Q^3
                = -mu f^4 (1 + kf^6) / (1 + kf^6)^3
                = -mu f^4 / (1 + kf^6)^2

**CHECK: Matches PROPOSAL.md lines 78-81. CORRECT.**

### 4. Sign Analysis with mu = -20

With mu < 0:

    lambda_anti = -mu f^4 / (1 + kf^6)^2 = |mu| f^4 / (1 + kf^6)^2 > 0  ALWAYS

So m^2_anti = m^2 + lambda_anti > m^2. The antisymmetric mode is HARDENED (stiffened), not softened.

**This means the PROPOSAL.md got the sign wrong** (lines 92-97 state m^2_anti is "ALWAYS reduced").
The RESULTS.md correctly identifies this error (lines 29-36).

For the symmetric eigenvalue:

    lambda_sym = mu f^4 (5 - 7kf^6) / Q^3

With mu < 0:
- When 7kf^6 < 5 (i.e., f < (5/(7k))^{1/6}): lambda_sym < 0, so m^2_sym is REDUCED (softened)
- When 7kf^6 > 5: lambda_sym > 0, so m^2_sym is INCREASED

With kappa = 20: f_crit = (5/140)^{1/6} = 0.0357^{1/6} = 0.574

At f(0) = 0.506 < 0.574, so 7kf^6 = 7*20*0.506^6 = 140*0.0168 = 2.35 < 5.
Therefore lambda_sym < 0 (softened) at this amplitude.

**CHECK: The RESULTS.md conclusion that the symmetric mode is softened and the antisymmetric
mode is hardened is CORRECT. The PROPOSAL.md had the wrong expectation.**

### 5. Numerical Verification at f = 0.506, kappa = 20

    f^4 = 0.506^4 = 0.0655
    f^6 = 0.506^6 = 0.01676
    kf^6 = 20 * 0.01676 = 0.3352
    Q = 1 + 0.3352 = 1.3352
    Q^2 = 1.7828
    Q^3 = 2.3803

    lambda_sym = (-20)(0.0655)(5 - 7*0.3352) / 2.3803
               = (-20)(0.0655)(5 - 2.347) / 2.3803
               = (-20)(0.0655)(2.653) / 2.3803
               = (-20)(0.1738) / 2.3803
               = -3.475 / 2.3803
               = -1.460

    lambda_anti = -(-20)(0.0655) / 1.7828
                = 20 * 0.0655 / 1.7828
                = 1.310 / 1.7828
                = 0.735

    m^2_sym(0) = 1.0 + (-1.460) = -0.460
    m^2_anti(0) = 1.0 + 0.735 = 1.735
    Delta m^2 = -0.460 - 1.735 = -2.195

**RESULTS.md reports:**
- lambda_sym(0) = -1.455  (I get -1.460, difference < 0.4%, likely from f(0) precision)
- lambda_anti(0) = +0.735  (matches)
- m^2_sym(0) = -0.455  (I get -0.460, consistent with lambda_sym rounding)
- m^2_anti(0) = +1.735  (matches)
- Delta m^2(0) = -2.190  (I get -2.195, consistent)

**CHECK: Numerical values are CORRECT to within rounding from the f(0) = 0.506 precision.**

### 6. Code Verification

The code (hessian.c) correctly implements:
- Potential force: `-mu * P * dP / denom2` (line 61) -- CORRECT
- Eigenvalue formulas: lines 289-290 match the analytic expressions -- CORRECT
- Time-averaging: numerical integration over one period (lines 327-341) rather than
  using the approximate (3/8)f^4 formula -- this is MORE ACCURATE than the proposal
  suggested, because the Hessian is a nonlinear function of f*cos(theta) so simple
  moment formulas (3/8)f^4 etc. do not suffice. CORRECT approach.
- The Velocity Verlet integrator is standard. CFL condition is reasonable.
- The envelope tracking (recording max |phi| over the last 20% of simulation) is a
  standard technique.

**One minor issue**: The comment on line 19 states "⟨f^4 cos^4⟩ = (3/8)f^4_peak" but
the code correctly does numerical time-averaging instead. The comment is misleading but
the code is correct.

### V23-A Verdict

- **Math**: CORRECT. The Hessian derivation, eigenvalue formulas, and sign analysis are all verified.
- **Code**: CORRECT. Implements the formulas accurately, with proper time-averaging.
- **PROPOSAL.md error**: Lines 92-97 incorrectly state that lambda_anti is positive and m^2_anti
  is reduced. The sign convention was confused. With mu < 0, lambda_anti = |mu|f^4/(1+kf^6)^2 > 0,
  which INCREASES m^2_anti. The RESULTS.md correctly identifies and fixes this error.
- **Conclusions**: SOUND. The splitting exists, is large (219% of m^2), and has the
  "wrong sign" for the downstream gravity proposals (shear is hardened, not softened).
  This conclusion correctly follows from the data.

---

## V23-B: Cross-Gradient Oscillon Dynamics

### 1. Modified EOM Derivation

The cross-gradient Lagrangian density term is:

    L_cross = -(1/2) eta (d_i phi_j)(d_j phi_i)

where i, j are spatial indices and a, b are field indices, with the "elastic" identification
field_index = spatial_direction.

The Euler-Lagrange equation for phi_a gets a contribution from L_cross:

    d/d(d_i phi_a) [ -(1/2) eta (d_j phi_k)(d_k phi_j) ]

Since L_cross involves (d_i phi_j)(d_j phi_i), we need:

    dL_cross / d(d_i phi_a) = -(1/2) eta [ delta_{ja} d_j phi_i  (from first factor, j->a, gives d_a phi_i)
                                          + delta_{ia} d_k phi_j  (from second factor, but this needs more care) ]

Actually, let me be more careful. Writing L_cross = -(1/2) eta sum_{i,j} (d_i phi_j)(d_j phi_i):

    dL_cross / d(d_m phi_a) = -(1/2) eta sum_{i,j} [ delta_{im} delta_{ja} (d_j phi_i)
                                                     + (d_i phi_j) delta_{jm} delta_{ia} ]
                             = -(1/2) eta [ d_a phi_m + d_m phi_a ]

Then the EOM contribution (from the spatial part) is:

    d_m [ dL_cross / d(d_m phi_a) ] = -(1/2) eta d_m [ d_a phi_m + d_m phi_a ]
                                    = -(1/2) eta [ d_a (d_m phi_m) + d_m d_m phi_a ]
                                    = -(1/2) eta [ d_a (div phi) + nabla^2 phi_a ]

Wait -- this gives TWO terms, not one. The EOM should be:

    d^2 phi_a / dt^2 = nabla^2 phi_a + (1/2) eta [ d_a(div phi) + nabla^2 phi_a ] - m^2 phi_a - dV/dphi_a
                     = (1 + eta/2) nabla^2 phi_a + (eta/2) d_a(div phi) - m^2 phi_a - dV/dphi_a

**ERROR FOUND**: The PROPOSAL.md (line 30) gives:

    d^2 phi_a / dt^2 = nabla^2 phi_a + eta d_a(div phi) - m^2 phi_a - dV/dphi_a

This is INCORRECT. The cross-gradient term -(1/2)eta (d_i phi_j)(d_j phi_i) contributes
BOTH eta/2 * d_a(div phi) AND eta/2 * nabla^2 phi_a to the EOM.

**However**, there is an alternative convention. The Lagrangian could be written as:

    L_cross = -eta (d_i phi_j)(d_j phi_i)    (without the 1/2 factor)

In that case, the EOM contribution is:

    eta d_a(div phi) + eta nabla^2 phi_a

which still gives TWO terms, not matching the PROPOSAL.

Or, a third possibility: the intended Lagrangian is actually just:

    L_cross = -eta (d_i phi_j)(d_i phi_j)    (same indices -- this is NOT cross-gradient)

which would give eta nabla^2 phi_a, a simple rescaling.

**Let me check what the 1D reduction in PROPOSAL.md gives.**

In 1D (x-direction only), the PROPOSAL claims (lines 39-50):
- For a=1 (x-component): extra term = eta * d^2 phi_1/dx^2
- For a=2,3: extra term = 0

With the correct cross-gradient L_cross = -(1/2)eta (d_i phi_j)(d_j phi_i), in 1D:

    L_cross = -(1/2) eta (d_x phi_1)(d_x phi_1)   [only i=x, j=1 contributes since only d_x exists]

Wait -- the sum over (i,j) with i,j spatial and field indices identified means in 1D only
i = x = 1 is a valid spatial index. So:

    L_cross = -(1/2) eta (d_1 phi_1)(d_1 phi_1) = -(1/2) eta (d_x phi_1)^2

This gives a contribution to the phi_1 EOM:

    d_x [ eta (d_x phi_1) ] = eta d^2 phi_1 / dx^2

And zero for phi_2, phi_3 (since there are no d_y, d_z in 1D).

**So in 1D, the result is:**
- phi_1: (1 + eta) d^2 phi_1/dx^2  (coefficient enhanced)
- phi_2: d^2 phi_2/dx^2  (unchanged)
- phi_3: d^2 phi_3/dx^2  (unchanged)

**This matches the PROPOSAL (lines 48-50) and the code (crossgrad1d.c line 168).**

The issue is only with the 3D formula. In 3D, the correct EOM from
L_cross = -(1/2)eta (d_i phi_j)(d_j phi_i) would be:

    d^2 phi_a/dt^2 = (1 + eta/2) nabla^2 phi_a + (eta/2) d_a(div phi) - ...

But the PROPOSAL writes (line 30):

    d^2 phi_a/dt^2 = nabla^2 phi_a + eta d_a(div phi) - ...

These give DIFFERENT results in 3D. However, since V23-B only uses the 1D code, let me
check whether the 1D reduction is consistent.

In 1D, the correct EOM (from the 3D version with (1+eta/2) nabla^2 + (eta/2) d_a(div)):
- For a=1: (1+eta/2) d_xx phi_1 + (eta/2) d_x(d_x phi_1) = (1+eta) d_xx phi_1
- For a=2: (1+eta/2) d_xx phi_2 + (eta/2) d_y(d_x phi_1) = (1+eta/2) d_xx phi_2
            (since d_y does not exist in 1D, the second term vanishes)

Wait, this gives (1+eta/2) for phi_2, not 1. That does NOT match the code.

Actually let me reconsider. The issue is what "d_a(div phi)" means for a=2 in 1D:

    d_a(div phi) = d_y(d_x phi_1 + d_y phi_2 + d_z phi_3) = d_y(d_x phi_1) = 0 in 1D

So the d_a(div phi) term is zero for a=2,3 in 1D. The ONLY contribution from the cross-gradient
to phi_2 in 1D comes from the (eta/2) nabla^2 phi_a term, giving (1+eta/2) d_xx phi_2.

But the code has coefficient 1.0 for phi_2, phi_3 (line 168). This means the code is using
the PROPOSAL's (incorrect) 3D formula, not the correct one. Since the PROPOSAL drops the
eta/2 * nabla^2 phi_a term, it gets coefficients (1+eta, 1, 1) in 1D, whereas the correct
coefficients should be (1+eta, 1+eta/2, 1+eta/2).

**HOWEVER**: there is another possible reading. If the intended term is NOT
-(1/2)eta (d_i phi_j)(d_j phi_i) but rather a different normalization or a different
Lagrangian term altogether, the 1D reduction could be different. Looking at the proposal
more carefully, it says the EOM is:

    d^2 phi_a/dt^2 = nabla^2 phi_a + eta d_a(nabla . phi)

This is exactly the equation for isotropic elasticity, where the standard Lame form is:

    rho d^2 u_a/dt^2 = mu nabla^2 u_a + (lambda + mu) d_a(div u)

With mu=1 and lambda+mu=eta, i.e., lambda = eta-1. In that convention, the Lagrangian is:

    L = (1/2)(d_t u)^2 - (1/2) mu (d_i u_j)(d_i u_j) - (1/2) lambda (div u)^2

    = (1/2)(d_t u)^2 - (1/2)(d_i u_j)^2 - (1/2)(eta-1)(div u)^2

This is a DIFFERENT Lagrangian from (d_i phi_j)(d_j phi_i). The "cross-gradient" term
(d_i phi_j)(d_j phi_i) is related but distinct from (div phi)^2.

In fact, (d_i phi_j)(d_j phi_i) = (d_i phi_j)(d_j phi_i), while
(div phi)^2 = (d_i phi_i)(d_j phi_j).

These are different contractions. From the EOM the PROPOSAL gives, the implied Lagrangian
is L_extra = -(eta/2)(div phi)^2, NOT -(eta/2)(d_i phi_j)(d_j phi_i).

With L_extra = -(eta/2)(div phi)^2:

    dL_extra/d(d_m phi_a) = -eta (div phi) delta_{ma}

    => EOM contribution: d_m [ eta (div phi) delta_{ma} ] = eta d_a(div phi)

This gives EXACTLY the PROPOSAL equation. And in 1D:

    L_extra = -(eta/2)(d_x phi_1)^2

which gives phi_1 extra term eta d_xx phi_1, and nothing for phi_2, phi_3.
Coefficients: (1+eta, 1, 1). MATCHES the code.

**VERDICT on the EOM**: The PROPOSAL title says "cross-gradient" eta(d_i phi_j)(d_j phi_i)
but the actual EOM implemented corresponds to the DIVERGENCE term eta(div phi)^2 / 2.
These are different Lagrangian terms! The "cross-gradient" (d_i phi_j)(d_j phi_i) would
give a different 3D equation. However:

- In 1D, both give the same phi_1 equation (both reduce to eta*(d_x phi_1)^2).
- In 1D, the divergence term gives (1+eta, 1, 1) for the Laplacian coefficients.
- The cross-gradient would give (1+eta, 1+eta/2, 1+eta/2).
- The code implements (1+eta, 1, 1), consistent with the DIVERGENCE term.

**This is a naming error, not a physics error for the 1D investigation.** The 1D results
are correct for the Lagrangian L_extra = -(eta/2)(div phi)^2. If the intent was truly
(d_i phi_j)(d_j phi_i), the 1D code is wrong for phi_2 and phi_3.

### 2. Code Implementation Check

The code (crossgrad1d.c) correctly:
- Applies Laplacian coefficient (1+eta) to phi_1 and 1.0 to phi_2, phi_3 (line 168, 177)
- Uses force_pot() for the potential force (line 177) -- correct sign and form
- Records the full 3-field dynamics independently (not assuming symmetry)
- Adds antisymmetric perturbation delta_phi1 = +eps, delta_phi2 = -eps (lines 338-339)
- Measures growth rate from log-envelope fit (lines 448-477) -- reasonable method

### 3. Vacuum Dispersion

With the implemented Lagrangian (divergence term), vacuum dispersion:
- Longitudinal (phi_1 mode, k along x): omega^2 = (1+eta)k^2 + m^2, so c_L = sqrt(1+eta)
- Transverse (phi_2, phi_3 modes): omega^2 = k^2 + m^2, so c_T = 1

**CHECK: PROPOSAL statement "eta only changes c_L, not c_T" (line 7) is CORRECT for the
divergence-term Lagrangian.**

If the true cross-gradient (d_i phi_j)(d_j phi_i) were used instead, in 3D the transverse
speed would also be modified: c_T^2 = 1 + eta/2, c_L^2 = 1 + 3eta/2.

### 4. Antisymmetric Perturbation Measurement

The perturbation method:
- Applied as delta_phi1 = +eps*G(x), delta_phi2 = -eps*G(x) with G Gaussian (line 338-339)
- This excites the (1,-1,0)/sqrt(2) mode, which is one of the antisymmetric eigenvectors
- Growth rate measured from windowed peak envelope fit to log(amplitude) vs time

**Potential issue**: The growth rate measurement uses 20 windows over the perturbation
evolution time (1000 time units). With omega ~ 0.87-0.99, this covers ~140-160 oscillation
periods. The sliding-window peak envelope is a reasonable but crude measure.

**The RESULTS.md correctly notes** that the "grew" flag is misleading (it triggers when
final > 1.5 * initial, but initial is measured at one instant, possibly catching a
zero-crossing). The growth RATE is the reliable measure, and all rates are negative.

### V23-B Verdict

- **Math**: PARTIALLY INCORRECT. The title and PROPOSAL say "cross-gradient" (d_i phi_j)(d_j phi_i)
  but the EOM and code implement the divergence coupling (div phi)^2. These are different
  in 3D but identical in 1D for phi_1. The 1D results are correct for the implemented
  Lagrangian, but the naming is misleading. If 3D follow-up (Phase 3) were done using the
  same EOM, it would NOT correspond to the cross-gradient term as stated.
- **Code**: CORRECT for what it implements (divergence coupling in 1D). The Velocity Verlet
  integrator, CFL condition, absorbing boundaries, and perturbation measurement are all sound.
- **Conclusions**: SOUND for the investigated case. The oscillon survives eta up to 2.0 with
  progressive weakening. No antisymmetric instability is found. The breathing frequency
  approaches the mass gap. These conclusions follow from the data regardless of the
  naming confusion.
- **Flag**: If the intent was truly the cross-gradient coupling (d_i phi_j)(d_j phi_i),
  the 3D investigation would need different equations. In 3D, the cross-gradient gives
  both a d_a(div phi) term AND a nabla^2 phi_a enhancement, affecting transverse modes
  as well as longitudinal. This could change the conclusions about shear stability.

---

## V23-C: Critical Correlation Length

### 1. Tail Fitting Method

The RMS envelope method (critical1d.c lines 283-286, 349-381):

1. Accumulate phi^2(x,t) during second half of evolution (t > 5000)
2. Compute phi_rms(x) = sqrt(<phi^2>_t)
3. Fit log(phi_rms) vs x in range x in [3, 25]

**Is the RMS envelope the right thing to fit?**

For an oscillon with phi(x,t) = A(x) cos(omega t + theta(x)):

    <phi^2>_t = A(x)^2 / 2

So phi_rms = A(x)/sqrt(2), which correctly captures the spatial envelope A(x).

For large x (tail region), A(x) ~ C * exp(-kappa x) where kappa = sqrt(m^2 - omega^2).

So log(phi_rms) ~ const - kappa * x, and the slope gives -kappa. **This is correct.**

**However**, there is a subtlety: the oscillon radiates, and the radiation field has
omega = m (at threshold), giving kappa_rad = 0. So at large x, the field is a
superposition of the bound tail (decaying) and radiation (constant amplitude in 1D,
since there is no geometric spreading). The fit range [3, 25] should capture the
exponential tail before the radiation floor dominates.

The noise floor cutoff (phi_rms > peak * 1e-6) is reasonable but may not be sufficient
if radiation amplitude is above this level. The RESULTS.md acknowledges that the tail
has "additional structure beyond the asymptotic exponential" and that the ratio xi/xi_pred
ranges from 1.7 to 6.

**Assessment**: The method is reasonable but has known limitations. The main issue is that
the breathing mode creates a standing-wave pattern that can extend beyond the simple
exponential tail, biasing kappa_tail toward smaller values (larger xi). This explains
why xi_meas > xi_pred consistently.

### 2. Predicted xi = 1/sqrt(m^2 - omega^2)

For a mode with frequency omega in a massive field with mass m, the spatial decay rate
outside the source is:

    kappa = sqrt(m^2 - omega^2)

and xi = 1/kappa.

**This is CORRECT for a single-frequency mode in a free massive field.** The oscillon
tail should indeed decay as exp(-kappa r) (in 1D, without the 1/r prefactor that appears
in 3D -- note: in 1D the Green's function gives pure exponential, no 1/r factor).

At the v21 parameters (omega = 0.87 from RESULTS.md):
    kappa = sqrt(1 - 0.87^2) = sqrt(1 - 0.7569) = sqrt(0.2431) = 0.493
    xi_pred = 1/0.493 = 2.03

RESULTS.md gives xi_pred = 2.13 for mu = -20, omega = 0.883. Checking:
    kappa = sqrt(1 - 0.883^2) = sqrt(1 - 0.7797) = sqrt(0.2203) = 0.469
    xi_pred = 1/0.469 = 2.13  -- MATCHES.

**CHECK: The formula and its numerical evaluation are CORRECT.**

### 3. Phase 2 Force Measurement Reliability

The Phase 2 uses Gaussian initial conditions (not equilibrated oscillons) and measures
the force from a quadratic fit to separation vs time in the second half of the run.

**Problems correctly identified in RESULTS.md:**
1. Gaussian initial conditions cause violent shedding (50-80% energy loss)
2. Radiation pressure from shedding dominates the centroid motion
3. The separation vs time is not parabolic
4. The quadratic fit gives unreliable F_avg values

**Additional issues not mentioned:**

1. **Centroid measurement with two merging oscillons**: The left/right centroid method
   (lines 413-437 of critical1d.c) splits the grid at the midpoint and computes
   phi^2-weighted centroid in each half. When shedding radiation reaches the other half,
   the centroid is contaminated. This is especially bad when the two oscillons overlap.

2. **In-phase vs out-of-phase**: The two oscillons are initialized in-phase (both
   positive Gaussians). This is not guaranteed to match the natural phase relationship
   during evolution. Phase mismatch changes the interaction from attractive to repulsive.

3. **The D=50 "attraction" at mu=-20 and D=20 "attraction" at mu=-12** could be artifacts
   of shedding dynamics, not genuine tail-mediated attraction.

**Assessment**: The Phase 2 results are NOT RELIABLE for measuring the inter-oscillon
force, as the RESULTS.md honestly acknowledges. The conclusions about interaction range
should not be taken as quantitative. The qualitative statement that "proper Phase 2
requires equilibrated oscillons" is correct.

### 4. Mu Scan Implementation

The code correctly:
- Scans mu from -20 to -5 (line 43)
- Uses the same initialization (A=1.0, sigma=3.0) for all mu values
- Measures frequency from DFT of second-half time series (lines 303-345)
- Fits tail decay from RMS profile (lines 348-381)
- Reports alive/dead status based on peak amplitude, core fraction, and frequency (line 392)

**Issue with initialization**: The same A_init=1.0 is used for all mu values. At weak
coupling (|mu| small), the nonlinear frequency shift is smaller, so the oscillon
amplitude should be SMALLER for a given gap margin. Using A=1.0 at mu=-5 means the
initial condition is further from the true oscillon than at mu=-20. This explains the
anomalous behavior at intermediate mu values: the system is evolving toward different
attractors depending on how far the initial condition is from the true oscillon.

**Issue with the anomaly at mu=-16, -14**: The RESULTS.md notes anomalously LOW omega
at these values. This could be:
1. A mode transition (as suggested)
2. The oscillon settling to a different branch (higher-order breathing mode)
3. A numerical artifact from insufficient equilibration time

The non-monotonic behavior of omega vs mu is suspicious and warrants investigation,
but the RESULTS.md appropriately flags it as an anomaly.

### V23-C Verdict

- **Math**: CORRECT. The correlation length formula xi = 1/sqrt(m^2 - omega^2) is standard
  and correctly applied. The 1D tail decay is pure exponential (no 1/r factor), correctly
  handled in the fit.
- **Code**: CORRECT in implementation, but with limitations:
  - Same initial amplitude for all mu values biases the scan
  - Phase 2 uses un-equilibrated initial conditions, making force measurements unreliable
  - The RESULTS.md honestly acknowledges these limitations
- **Conclusions**: SOUND but INCOMPLETE.
  - The finding "xi does NOT diverge in the accessible range" is correctly stated as an
    observational result. However, it does not rule out divergence at smaller |mu| (below 5).
  - The xi/xi_pred ratio trending toward 1-2 at weak coupling supports the theoretical
    formula; the discrepancy at strong coupling is explained by radiation contamination.
  - Phase 2 results are correctly flagged as unreliable.
  - The overall conclusion that "the critical gravity proposal is not supported by these
    data" is appropriately cautious.

---

## Global Summary

| Investigation | Math | Code | Conclusions |
|--------------|------|------|-------------|
| V23-A Hessian | CORRECT | CORRECT | CORRECT (proposal sign error found and fixed in results) |
| V23-B Cross-gradient | NAMING ERROR (1D code is correct for what it implements) | CORRECT (for divergence coupling) | SOUND (for the implemented Lagrangian) |
| V23-C Critical | CORRECT | CORRECT (with noted limitations) | SOUND but incomplete |

### Issues That Could Change Conclusions

1. **V23-B**: If the true "cross-gradient" (d_i phi_j)(d_j phi_i) were implemented in 3D
   instead of the divergence coupling (d_i phi_i)^2, the transverse wave speed would also
   increase (c_T^2 = 1 + eta/2), potentially changing the shear stability results. This
   matters for any 3D follow-up.

2. **V23-C**: The non-monotonic omega(mu) behavior at mu=-16,-14 suggests that the
   A_init=1.0 Gaussian is not close enough to the true oscillon at these coupling values.
   Re-running with A_init adapted to each mu (e.g., from a rough amplitude-mu relationship)
   would give cleaner results.

3. **V23-A**: No issues that could change the conclusion. The sign of the splitting is
   algebraically determined and robust.
