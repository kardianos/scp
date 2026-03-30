# V49 Proposal: Jacobian-Squared Theta Mass (m_theta^2 = lambda_theta * J^2)

## Motivation

V48 introduced field-dependent theta mass via m_theta^2 = lambda_theta * P^2, where
P = phi_0 * phi_1 * phi_2 (triple product). While this successfully gave theta
a position-dependent mass at braid cores (V48 runs confirmed stability through
T=5000 at lambda_theta=3.5), the P^2 order parameter has known deficiencies:

1. **P != 0 in the background.** The carrier wave gives P = (A^3/4) cos(3kz),
   so <P^2> = A^6/32 ~ 3e-8. This is negligible for the mass term (Yukawa range
   ~3000 in background vs ~7 at core), but becomes problematic for any coupling
   that uses P^2 as a modulator (e.g., eta(P^2) screening), where the 6*omega
   oscillation of P^2 creates parametric pumping in the vacuum.

2. **P^2 oscillates at 6x carrier frequency**, even inside the core. The mass
   term "strobes" — P^2 crosses zero six times per carrier cycle, momentarily
   making theta massless everywhere. For the current lambda_theta*P^2
   implementation this is tolerable (the time-averaged effect is what matters
   for binding), but it prevents the formation of a truly "hard" theta shell
   with a clean spectral gap.

3. **P^2 gives a soft boundary.** The transition from P!=0 (core) to P~0 (far
   field) follows the smooth decay of the field amplitudes. There is no sharp
   surface separating the gapped and gapless regions.

## The Jacobian Determinant

Replace P^2 with J^2 as the order parameter, where:

    J = det(d_i phi_j) = epsilon_ijk (d_x phi_i)(d_y phi_j)(d_z phi_k)

This is the determinant of the 3x3 Jacobian matrix of the field map
phi: R^3 -> R^3, evaluated at each grid point.

### What J fixes and what it doesn't

**Problem 1 (background P != 0): FIXED.** For the background carrier wave
phi_a = A cos(kz + delta_a), the gradient matrix is:

    d_i phi_j = | 0  0  0 |
                | 0  0  0 |
                | *  *  * |

Only the z-row is nonzero (the wave propagates in z). This matrix has **rank 1**.
Its determinant is **exactly zero, always, at all times.** This is not an
approximation — it is a theorem of linear algebra (rank < 3 => det = 0).

Therefore m_theta^2 = lambda_theta * J^2 = 0 exactly in the background. The
photon is truly massless in free space. No parametric resonance, no Yukawa
cutoff at any range.

**Caveat — crossing waves:** Where radiated waves from different particles
intersect at non-parallel angles, the gradient matrix can gain rank, giving
J != 0 in the vacuum. However, far-field wave amplitudes are ~0.01 A_bg, so
J_vacuum ~ (0.01 A k)^3 ~ 10^{-9} J_core. The induced photon mass is
negligible (~10^{-9} of core mass). This is analogous to vacuum polarization
in QED — a tiny correction, not a fatal flaw.

**Problem 2 (strobing): NOT FIXED.** For a braid with separable time dependence
phi_a(x,t) = A_a(x) cos(omega*t + delta_a), the Jacobian factors as:

    J(x,t) = det(nabla A_a) * cos(wt) cos(wt+2pi/3) cos(wt-2pi/3)
            = det(nabla A_a) * (1/4) cos(3wt)

J oscillates at 3*omega and crosses zero 6 times per carrier cycle, just like P.
The rank of the gradient matrix is time-invariant (the helix always has 3D
structure), but det goes through zero because all three rows simultaneously
pass through zero amplitude at their respective phase crossings.

Therefore J^2, like P^2, strobes at 6*omega. The theta mass blinks on and off.
This is the SAME strobing problem as P^2. The advantage of J over P is ONLY in
the background (exactly zero vs small-but-nonzero).

**Practical implication:** The V48 P^2 simulations showed that despite strobing,
the time-averaged mass is sufficient for effective Yukawa behavior — the proton
survives and breathes coherently for T=5000. The strobing is inelegant but not
fatal in practice. J^2 would behave similarly at the core.

**Problem 3 (soft boundary): PARTIALLY FIXED.** The spatial profile of J is
sharper than P because J involves derivatives (concentrated at the braid surface
where gradients are steepest), while P involves amplitudes (spread over the core
volume). This gives a somewhat harder transition, but it is NOT a true "hard
shell" because the time-oscillation still makes the boundary permeable 6x/cycle.

### Physical interpretation

J = det(d_i phi_j) is the **topological charge density** (up to normalization).
It is nonzero only where the field map phi: R^3 -> R^3 is locally a
diffeomorphism — i.e., where all three field components have independent spatial
variation. This is precisely the condition that distinguishes a 3D topological
structure (braid, hedgehog) from a 1D wave (background carrier).

Using J^2 as the theta mass means: **theta acquires mass only where there is
genuine 3D topological structure.** The "gap" is topologically induced. No
particle = no gap = massless photon. This is the cleanest realization of the
"emergent pion" concept from the V48 proposal.

### Derrick scaling

J involves 3 spatial derivatives, so J^2 involves 6 derivative powers.
The energy integral of the mass term:

    E_J = (1/2) lambda_theta integral J^2 |theta|^2 d^3x

Under x -> lambda*x (Derrick scaling): J ~ lambda^(-3) (three derivatives),
so J^2 ~ lambda^(-6), d^3x ~ lambda^3. Total scaling: lambda^(-3).

This is STRONGLY anti-Derrick (prevents collapse) — even more so than the
Skyrme E_4 term (lambda^(-1)). Combined with gradient energy (lambda^(+1)),
it creates a preferred soliton size, stabilizing the braid against both
expansion and collapse.

**Independent utility:** Even if J^2 is not used as a theta mass, adding a
term Delta_V = +lambda_J J^2 to the phi potential would provide Derrick
stabilization of the soliton. This may be worth exploring independently.

## The Equation of Motion

### Lagrangian

    L_J = -(1/2) lambda_theta J^2 |theta|^2

where J = det(d_i phi_j) and |theta|^2 = theta_0^2 + theta_1^2 + theta_2^2.

### Theta equation

Variation with respect to theta_a is straightforward (theta appears only
algebraically in L_J):

    d^2 theta_a / dt^2 = laplacian(theta_a) - (m_theta^2 + lambda_theta J^2) theta_a
                         + eta curl(phi)_a

The effective theta mass is:

    m_theta^2_eff(x) = m_theta^2 + lambda_theta * J(x)^2

This is identical in form to the V48 P^2 implementation, just with J^2 replacing P^2.

### Phi equation

Variation with respect to phi_a is more complex because J depends on the
spatial derivatives of phi. Writing J = det(M) where M_ij = d_j phi_i, and
C_aj = cofactor(M, a, j):

    dJ / d(d_j phi_a) = C_aj

The Euler-Lagrange equation gives the force on phi_a as:

    F_J_a = -d_j [ lambda_theta |theta|^2 J C_aj ]

**IMPORTANT:** The spatial derivative d_j acts on the FULL product
(lambda_theta |theta|^2 J C_aj). The |theta|^2 factor CANNOT be pulled outside
the derivative because theta varies in space. Factoring it out would break
energy conservation and violate Newton's third law between the phi and theta
fields.

### Simplification via the Piola Identity

The cofactor matrix C_aj satisfies the **Piola identity** from continuum mechanics:

    d_j C_aj = 0    (exactly, always)

This is because C_aj involves first derivatives of phi with antisymmetric
(epsilon) structure, and mixed partial derivatives commute. The divergence
of the cofactor of any gradient field vanishes identically.

Applying the product rule to the force:

    F_J_a = -d_j [ lambda_theta |theta|^2 J C_aj ]
          = -lambda_theta [ d_j(|theta|^2 J) C_aj + |theta|^2 J (d_j C_aj) ]
          = -lambda_theta [ d_j(|theta|^2 J) C_aj + 0 ]
          = -lambda_theta C_aj d_j(|theta|^2 J)

The force simplifies to a contraction of the cofactor matrix with the gradient
of the scalar field (|theta|^2 * J). This requires:

1. C_aj — the 9 cofactor components (already computed for J)
2. d_j(|theta|^2 * J) — gradient of a scalar, computed by finite-differencing
   the product |theta|^2 * J at neighboring voxels

**No cross-second-derivatives (d_xy, d_xz, d_yz) are needed.** The Piola
identity eliminates the d_j C_aj term that would have required them. This
dramatically simplifies the implementation — no face-diagonal stencils,
no 27-point neighborhoods, just the standard 7-point stencil applied to a
pre-computed scalar field.

### Explicit cofactor computation

Let M_aj = d_j phi_a (gradient matrix). The cofactors are:

    C_00 = (d_y phi_1)(d_z phi_2) - (d_z phi_1)(d_y phi_2)
    C_01 = (d_z phi_1)(d_x phi_2) - (d_x phi_1)(d_z phi_2)
    C_02 = (d_x phi_1)(d_y phi_2) - (d_y phi_1)(d_x phi_2)

    C_10 = (d_y phi_2)(d_z phi_0) - (d_z phi_2)(d_y phi_0)
    C_11 = (d_z phi_2)(d_x phi_0) - (d_x phi_2)(d_z phi_0)
    C_12 = (d_x phi_2)(d_y phi_0) - (d_y phi_2)(d_x phi_0)

    C_20 = (d_y phi_0)(d_z phi_1) - (d_z phi_0)(d_y phi_1)
    C_21 = (d_z phi_0)(d_x phi_1) - (d_x phi_0)(d_z phi_1)
    C_22 = (d_x phi_0)(d_y phi_1) - (d_y phi_0)(d_x phi_1)

Each cofactor requires 2 multiplies and 1 subtract (18 mul + 9 sub total).
The determinant is then J = phi_0 C_00 + phi_1 C_10 + phi_2 C_20... wait, no:
J = (d_x phi_0) C_00 + (d_y phi_0) C_01 + (d_z phi_0) C_02 (expansion along
row 0). Or equivalently, J = M_a0 C_a0 (sum over a, any fixed j).

## Implementation (using Piola-simplified force)

### In compute_forces, per voxel:

```
// 1. Compute 9 gradient components (6 axis-aligned already available for
//    Laplacian + curl; need 3 additional cross-gradients)
double dxp0 = (phi[0][ip] - phi[0][im]) * idx1;  // d_x phi_0
double dyp0 = (phi[0][jp] - phi[0][jm]) * idx1;  // d_y phi_0
double dzp0 = (phi[0][kp] - phi[0][km]) * idx1;  // d_z phi_0
// ... same for phi_1, phi_2 (9 components total)

// 2. Compute cofactors (18 mul + 9 sub)
double C00 = dyp1*dzp2 - dzp1*dyp2;
double C01 = dzp1*dxp2 - dxp1*dzp2;
double C02 = dxp1*dyp2 - dyp1*dxp2;
// ... C10-C22 by cyclic permutation

// 3. Compute J = det(M) (3 mul + 2 add)
double J = dxp0*C00 + dyp0*C01 + dzp0*C02;
double J2 = J * J;

// 4. Theta force: use J2 as mass
double mtheta2_eff = MTHETA2 + LAMBDA_THETA * J2;
theta_acc[a] = lapt - mtheta2_eff * theta[a] + eta * curl_phi[a];

// 5. Store J * |theta|^2 at this voxel for gradient computation
JT2[idx] = J * theta2;   // pre-computed scalar field
```

Then in a SECOND PASS (or with neighbor lookups):

```
// 6. Phi back-reaction: F = -lambda_theta * C_aj * d_j(J * |theta|^2)
double dJT2_dx = (JT2[ip] - JT2[im]) * idx1;
double dJT2_dy = (JT2[jp] - JT2[jm]) * idx1;
double dJT2_dz = (JT2[kp] - JT2[km]) * idx1;

phi_acc[0] -= LAMBDA_THETA * (C00*dJT2_dx + C01*dJT2_dy + C02*dJT2_dz);
phi_acc[1] -= LAMBDA_THETA * (C10*dJT2_dx + C11*dJT2_dy + C12*dJT2_dz);
phi_acc[2] -= LAMBDA_THETA * (C20*dJT2_dx + C21*dJT2_dy + C22*dJT2_dz);
```

### Two-pass requirement

The phi back-reaction requires d_j(J |theta|^2), which needs J*|theta|^2
at neighboring voxels. This means the force computation needs either:

(a) A two-pass approach: first pass computes J*|theta|^2 everywhere, second
    pass computes its gradient and applies the force. This doubles the memory
    bandwidth but is straightforward.

(b) Cache the J*|theta|^2 values in a temporary array (1 extra double per
    voxel = 57 MB at N=196).

Option (b) is recommended. The temporary array JT2[] is allocated once and
reused every timestep.

## Energy accounting

The J^2 energy contribution:

    E_J = (1/2) lambda_theta * sum_over_grid[ J^2 * |theta|^2 ] * dV

This fits into the existing E_tmass diagnostic slot:

    E_tmass = (1/2) * (m_theta^2 + lambda_theta * J^2) * |theta|^2

## Comparison to P^2 (corrected)

| Property               | P^2 = (phi_0 phi_1 phi_2)^2 | J^2 = det(d_i phi_j)^2      |
|------------------------|------------------------------|-------------------------------|
| Background value       | A^6/32 ~ 3e-8 (nonzero)     | **0 exactly**                 |
| Core strobing          | 6*omega (crosses zero)       | 6*omega (crosses zero) [SAME] |
| Shell sharpness        | Soft (amplitude decay)       | Somewhat harder (gradient-based) |
| Photon in free space   | Slightly massive             | **Exactly massless**          |
| Parametric resonance   | Present (e-fold ~10^5 t)     | **Absent in background**      |
| Derrick scaling        | lambda^3 (same as mass)      | **lambda^(-3) (anti-collapse)** |
| Computation cost       | 1 multiply                   | ~30 multiplies + temp array   |
| Second derivatives     | Not needed                   | **Not needed (Piola identity)** |
| Lagrangian consistent  | Yes                          | Yes                           |
| Back-reaction          | Simple (P * dP/dphi * theta2)| C_aj d_j(J theta2) — needs 2 passes |

The key advantages of J^2 over P^2 are: exactly zero background (no parametric
resonance in vacuum), anti-Derrick stabilization, and topological selectivity
(only genuine 3D structures generate mass). The key disadvantage is ~30x
computational cost.

Both share the same core strobing behavior. Neither produces a truly "hard"
shell in the time-domain sense.

## The Strobing Problem — Open Question

No purely local, instantaneous, algebraic function of the fields gives a
non-oscillating order parameter at the braid core. Every candidate (P^2, J^2,
|phi|^2, chi^2) oscillates at harmonics of the carrier frequency, because the
fields themselves oscillate.

The only approaches that avoid strobing:

1. **Time-averaged quantities** (running window <P^2> or <J^2>) — requires
   memory buffer, breaks time-locality in the Lagrangian.

2. **Topological winding numbers** (integral invariants) — nonlocal, require
   integration over surfaces, not available per-voxel.

3. **Accept the strobing** and rely on the time-averaged mass being sufficient.
   The V48 T=5000 runs confirm this IS adequate in practice: the proton
   survives, breathes coherently, and the effective Yukawa range is set by
   <m_theta^2> = lambda_theta <J^2> = (1/2) lambda_theta det(nabla A)^2.

Option 3 is the pragmatic path for near-term experiments. The strobing is
theoretically inelegant but not practically fatal.

## Suggested parameter range

At a braid core, J ~ (A/R)^3 ~ (0.3/3)^3 ~ 10^{-3}. So J^2 ~ 10^{-6}.

For m_theta_core ~ 0.15 (Yukawa range ~7, matching V48 target):

    lambda_theta = m_theta_core^2 / J_core^2 ~ 0.02 / 10^{-6} ~ 2e4

This is much larger than the P^2 case (3.5) because J^2 is much smaller than
P^2 at the core. The numbers need calibration from a test run measuring J_rms
at the braid core.

## Implementation plan

### Phase 1: Measure J at braid cores

Before modifying the kernel, run the existing proton template through an
analysis tool that computes J = det(d_i phi_j) at each voxel for a few frames.
This gives:
- J_peak at braid core
- J profile (radial decay)
- J_rms vs P_rms correlation
- Confirms J=0 in background
- Measures <J^2> to calibrate lambda_theta

### Phase 2: Implement in CPU kernel

Add to compute_forces (two-pass approach):
1. Pass 1: Compute 9 gradients, cofactors, J, and store JT2[idx] = J * theta2
2. Pass 2: Compute d_j(JT2) at each voxel, apply phi back-reaction via C_aj

New config parameter: `mass_theta_mode = 0` (P^2, default) or `1` (J^2).
Shares lambda_theta parameter.

### Phase 3: Stability + binding test

Same protocol as V48:
1. Single proton at N=196, T=200 — verify stability
2. Measure m_theta_eff at core vs background
3. Two-proton sweep at D=10,15,25,80 — test for binding

## Risks

1. **Computational cost**: ~30x more multiplies per voxel than P^2 plus a
   temporary array and two-pass force computation. May need GPU for production.

2. **Lambda_theta calibration**: J^2 is much smaller than P^2, so lambda_theta
   needs to be ~10^4. This large parameter could cause stiffness requiring
   smaller dt.

3. **Two-pass coupling**: The phi back-reaction requires a separate pass for
   the gradient of JT2. This complicates the integration loop and doubles
   memory bandwidth for the phi force. For GPU implementation, this means
   two kernel launches instead of one.

4. **Crossing-wave mass**: Between multiple particles, intersecting radiation
   gives J != 0 in the vacuum, inducing a tiny photon mass. Quantitatively
   negligible (~10^{-9} of core), but breaks the "exactly massless" claim
   in multi-particle systems.
