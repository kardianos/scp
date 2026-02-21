# V12: The Non-Linear Strong Geon

## Thesis

V8-V11 isolated three mechanisms for confinement:

| Version | Mechanism | Result | Fatal Flaw |
|---------|-----------|--------|-----------|
| V9 | Massive spin-2 Yukawa (f2(1270)) | 887 MeV = 94.5% M_p | Q9: 1/A runaway collapses to horizon |
| V10 | Decoupled g_tt/g_rr | 6% spatial repulsion | 1/A feedback independent of spatial curvature |
| V11 | BI nonlinearity + Hopf topology | 61% slower initial dispersal | No inward binding force; disperses at c |

These were tested in isolation. V12 combines them:

- **Binding force**: V9 Yukawa massive tensor (f2(1270)), provides 1 fm confinement
- **Energy cap**: Born-Infeld nonlinearity, prevents 1/A runaway and horizon formation
- **Topology**: Hopf magnetic knot, prevents leakage without artificial Rmax box

The central claim: BI nonlinearity cuts off the divergence that creates the Q9 paradox,
allowing a self-consistent equilibrium at proton-mass energies without horizon formation.

---

## 1. Why the 1/A Runaway Fails with BI

### The V9/V10 Runaway (Linear EM)

In the Regge-Wheeler standing wave with metric f(r) = A(r):

```
E_phys = omega * u(r) / sqrt(A)          (blueshifted electric field)
B_phys = u'(r) * sqrt(A)                  (spatial derivative, no blueshift)

rho_linear = (E_phys^2 + B_phys^2) / (8pi)
           = N_l/(8pi r^2) * [omega^2 u^2 / A + A (u')^2]
```

As A -> 0: the electric term omega^2 u^2 / A diverges. This sources deeper A,
which increases rho further. **Positive feedback loop**:

```
deeper A  -->  rho ~ 1/A  -->  larger M_eff  -->  deeper A  -->  ...
```

The loop has no fixed point except A = 0 (horizon).

### The BI Cutoff

Born-Infeld electrodynamics caps the electric field: E_phys < b. This means:

```
omega * u(r) / sqrt(A) < b
=>  A > omega^2 u^2(r) / b^2  =:  A_floor(r)
```

The metric **cannot** deepen below A_floor. The 1/A divergence is explicitly
cut off.

### CRITICAL: BI Caps E, but NOT the Energy Density

**The E-field formulation of the BI energy density DIVERGES as E -> b:**

```
rho_BI(E,B) = (b^2 + B^2) / sqrt(1 + (B^2 - E^2)/b^2) - b^2
```

As the gravitational well deepens (A -> 0), the blueshifted electric field
E_phys = omega u / sqrt(A) approaches b. The denominator sqrt(1 - E^2/b^2 + ...)
goes to ZERO. The energy density rho_BI -> INFINITY. If E_phys exceeds b, the
square root becomes imaginary and the code crashes with NaN.

**This is NOT a bounded energy density.** The scalar approximation
rho_BI = b^2(sqrt(1 + rho_M/b^2) - 1), which IS bounded, does not capture
this divergence because it treats rho_M as a single scalar rather than
distinguishing E and B contributions to the BI Lagrangian.

The plan's Phase 1 Step 5 ("If E > b, reduce u amplitude") is also wrong:
reducing u destroys the normalization condition int u^2 dr = 1, which
represents the quantization of the wave (a specific quantum of energy).
Artificially deleting energy to prevent a code crash is not a physical solution.

### The D-Field Formulation (Correct Approach)

The resolution is to parameterize the BI Hamiltonian by (D, B) instead of (E, B).
In BI electrodynamics, as E -> b, the displacement field D -> infinity. The
Hamiltonian expressed in terms of D is everywhere finite and well-behaved:

```
rho_BI(D, B) = b^2 (sqrt(R) - 1)

R = 1 + (D^2 + B^2)/b^2 + |D x B|^2/b^4
```

The |D x B|^2/b^4 cross-term is the nonlinear Poynting correction. It is
negligible for Phase 1 (standing waves with D and B 90 degrees out of phase)
but essential for Phase 3 (3D Hopfion with swirling, overlapping D and B).

Properties:
- D << b, B << b: rho_BI = (D^2 + B^2)/2 = rho_Maxwell (correct weak-field limit)
- D -> infinity: rho_BI ~ b|D| (grows linearly, NEVER diverges)
- B -> infinity: rho_BI ~ b|B| (grows linearly, NEVER diverges)
- R >= 1 for all (D, B) (no square root of negative)
- No NaN possible regardless of field strength

For a standing wave, D_phys = omega u / sqrt(A) (same blueshift), and the
energy density becomes:

```
rho_BI ~ b * D_phys ~ b * omega * u / sqrt(A) ~ 1/sqrt(A)
```

The feedback loop:
```
deeper A  -->  rho_BI ~ 1/sqrt(A) (not 1/A)  -->  weaker M growth  -->  A saturates
```

This is **sub-critical**: the 1/sqrt(A) source grows too slowly to sustain the
runaway. The Hartree iteration should converge to a finite A_min > 0.

### The A_floor Estimate

At the origin (r -> 0), the mode function u(0) is determined by the eigenvalue
problem. For the l=1 mode in V9:

```
u(0) ~ (Rmax/pi)^{1/2}  (normalized standing wave)
omega^2 ~ (pi/Rmax)^2     (box eigenvalue)
```

So A_floor = omega^2 u^2(0) / b^2 ~ 1/b^2. For b >> 1 (weak BI), A_floor ~ 0
and the V9 behavior is recovered. For b ~ 1 (strong BI), A_floor ~ 1 and the
well is shallow.

The physical b depends on the interpretation:
- QCD string tension: b ~ Lambda_QCD^2 ~ (200 MeV)^2 in physical units
- In V9 code units (L_0 = 1 fm, E_0 = 197.3 MeV): b ~ 1 code unit
- The Schwinger limit: b_QED ~ 10^18 V/m (irrelevant for nuclear scale)

**The critical test**: does there exist a b value where A_floor is large enough
to prevent collapse but small enough to give a deep binding well?

---

## 2. V9 Code Already Has BI — Discovery

The V9 source (`strong_geon.c:214-239`) already implements a BI-like energy cap:

```c
double rho_M = N_l / (8.0 * M_PI * r * r) *
               (omega2 * u * u / f + f * up[i] * up[i]);

if (b_BI > 0.0 && rho_M > b_BI * b_BI) {
    rho[i] = b_BI * b_BI * (sqrt(1.0 + rho_M / (b_BI*b_BI)) - 1.0);
} else {
    rho[i] = rho_M;
}
```

This applies the scalar BI formula rho_BI = b^2(sqrt(1 + rho_M/b^2) - 1) when
the Maxwell energy density exceeds b^2. The parameter is `b_BI` (default 0 = off).

**This was never tested at physical coupling.** All V9 results used b_BI = 0.

### Correction Needed to the BI Formula

The V9 code applies BI to the TOTAL energy density rho_M. This is a scalar
approximation (Level 1). It avoids the E-field divergence because it treats
rho_M as a single scalar without resolving the E vs B contributions.

The correct vector BI treatment distinguishes E and B. As shown in Section 1,
the E-field formulation DIVERGES as E -> b. The D-field formulation is safe:

```
rho_BI(D, B) = b^2 (sqrt(R) - 1)
R = 1 + (D^2 + B^2)/b^2 + |D x B|^2/b^4
```

For a standing wave (D and B are 90 degrees out of phase, D.B = 0),
the cross-term |D x B|^2 is suppressed by the phase offset. The Hamiltonian
simplifies to:

```
At time t:  D(t) = D_0 cos(wt),  B(t) = B_0 sin(wt)
|D x B|^2 = D_0^2 B_0^2 cos^2(wt) sin^2(wt) sin^2(angle)
rho_BI(t) = b^2 (sqrt(1 + D_0^2 cos^2(wt)/b^2 + B_0^2 sin^2(wt)/b^2
                       + D_0^2 B_0^2 cos^2 sin^2 sin^2(angle)/b^4) - 1)
```

For the 1D standing wave (Phase 1), D and B are parallel (angle = 0), so
|D x B| = 0 exactly. The cross-term vanishes and the simplified formula
without the Poynting term is exact for Phase 1.

The time-averaged density is:

```
<rho_BI> = (1/T) int_0^T rho_BI(t) dt
```

where D_0 and B_0 depend on (omega, u, u', A) from the mode equation.
In the weak-field limit (D ~ E), D_0 recovers the Maxwell E_0. In the
strong-field limit (E -> b), D_0 -> infinity but rho_BI remains finite.

For the V12 code, two levels of approximation:

**Level 1 (scalar BI):** Use V9's existing formula. Quick concept test. Will
correctly show that a nonlinear energy cap prevents the 1/A runaway. However,
it cannot capture the E -> b divergence/saturation physics because it doesn't
distinguish E from B.

**Level 2 (vector BI, D-field):** Compute D_0^2 and B_0^2 separately from the
mode function. In the deep-well regime, D_phys = omega u / sqrt(A) while
B_phys = u' sqrt(A). Time-average the D-field BI energy density numerically
(Gauss quadrature over one period). Feed into Yukawa BVP. This formulation
is NaN-free for arbitrarily deep wells.

Level 1 is sufficient for a quick proof-of-concept. Level 2 is required for
physically correct results in the deep-well regime where E_phys approaches b.

---

## 3. Phase 1: BI-Modified Hartree (Spherical Standing Wave)

### Goal

Run the V9 Hartree solver with b_BI > 0 at physical coupling (kappa=35, mu=6.47).
Determine if the Q9 phase transition is eliminated or smoothed.

### Equations (V9 framework, minimal modification)

**Metric:**
```
ds^2 = -f(r) dt^2 + dr^2/f(r) + r^2 dOmega^2
f(r) = 1/(1 + 2M_eff(r)/r)     (Pade, always > 0)
```

**Regge-Wheeler eigenvalue (unchanged):**
```
-f^2 v'' + W(r) v = omega^2 v
W(r) = f l(l+1)/r^2 + f f''/2 - (f')^2/4
v = sqrt(f) u
```

**BI energy density (Level 1 — scalar):**
```
rho_M = N_l/(8pi r^2) [omega^2 u^2/f + f (u')^2]

rho_BI = b^2 (sqrt(1 + rho_M/b^2) - 1)      if rho_M > 0
```

Properties:
- rho_M << b^2: rho_BI = rho_M (Maxwell limit)
- rho_M >> b^2: rho_BI = b sqrt(rho_M) (square-root cap)
- rho_BI < infinity for all rho_M

**BI energy density (Level 2 — vector, D-field formulation):**
```
D_0^2 = N_l omega^2 u^2 / (8pi r^2 f)       (displacement amplitude squared)
B_0^2 = N_l f (u')^2 / (8pi r^2)             (magnetic amplitude squared)

alpha = D_0^2 / b^2,  beta = B_0^2 / b^2

<rho_BI> = (b^2/pi) int_0^pi [sqrt(1 + alpha cos^2(t) + beta sin^2(t)) - 1] dt
```

Note: the |D x B|^2/b^4 Poynting cross-term vanishes exactly for the 1D
standing wave (D parallel to B at every point, so D x B = 0). The formula
above is therefore exact for Phase 1, not an approximation.

Properties:
- The integrand is sqrt(1 + positive + positive) - 1 >= 0 for ALL field strengths.
- No constraint on alpha or beta — the formula is valid for arbitrarily large D_0.
- In the weak-field limit (alpha, beta << 1): <rho_BI> = (D_0^2 + B_0^2)/4 = rho_Maxwell.
- In the strong-D limit (alpha >> 1): <rho_BI> ~ b D_0 / pi (linear growth, not divergent).
- No NaN possible. No need to reduce u amplitude.

**Yukawa BVP (unchanged):**
```
M_eff'' - mu^2 M_eff = -kappa r rho_BI
BC: M_eff(0) = 0,  M_eff'(Rmax) = -mu M_eff(Rmax)
```

**Hartree loop:**
```
1. Initialize: M_eff = 0 (flat space) or seeded Yukawa well
2. Compute metric: f = 1/(1 + 2M_eff/r)
3. Solve eigenvalue: omega^2 from Regge-Wheeler in potential W(f)
4. Reconstruct u(r), normalize (int u^2 dr = 1, do NOT reduce amplitude)
5. Compute rho_BI (Level 1 or Level 2 D-field formulation)
6. Solve Yukawa BVP -> M_new
7. Under-relax: M_eff = alpha M_new + (1-alpha) M_eff
8. Repeat 2-7 until delta_M < 10^{-14}
```

There is NO "reduce u amplitude" step. The D-field BI energy density is
finite for all field strengths, so the solver handles arbitrarily deep wells
without NaN. The BI saturation is SELF-CONSISTENT: as A deepens, rho_BI
grows as 1/sqrt(A) (sub-critical), the Yukawa source weakens relative to
the linear case, and A saturates at A_floor > 0.

**Note on the Regge-Wheeler approximation:** Phase 1 modifies the energy
density source rho to use BI, but leaves the wave equation itself linear
(-f^2 v'' + Wv = omega^2 v). This is not strictly correct — true BI theory
modifies the wave equation through its nonlinear constitutive relations.
However, this is an acceptable mean-field approximation for Phase 1, which
is explicitly a 5-minute proof-of-concept to test whether capping the energy
density stops the horizon collapse. Phase 3 solves the true nonlinear BI
PDEs on the 3D grid, making Phase 1's linear wave approximation moot.

### Parameter Scan

**Scan 1: b_BI at fixed kappa=35, mu=6.47, Rmax=2.85 (V9 proton match)**
```
b_BI = {0, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0, 1e10}
```
For each b_BI: record omega^2, E_total, f_min, well_depth, convergence.

Expected behavior:
- b_BI = 0: V9 result (collapsed, omega^2 ~ 0, E ~ 45 MeV)
- b_BI = 1e10: same as b_BI = 0 (BI inactive)
- b_BI ~ O(1): A_floor > 0, well depth limited, omega^2 higher, E possibly ~ M_p

**Scan 2: Rmax at fixed b_BI (optimal from Scan 1), kappa=35, mu=6.47**
```
Rmax = {1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0}
```
Check if the Q9 sharp transition is smoothed into a continuous crossover.

**Scan 3: 3-quark geon (2xl=1 + 1xl=2) at optimal b_BI**
Same as V9 Q7, with BI. Target: E_total = 938 MeV, no horizon.

### Success Criteria

**Positive (Q9 resolved)**:
- Hartree converges to f_min > 0.1 (no near-horizon)
- E_total in range [800, 1000] MeV
- Mode is localized (R_eff < 2 fm)
- No sharp phase transition between bound and unbound

**Negative (Q9 persists)**:
- Either f_min ~ 0 (collapsed) or E_total << M_p (too light)
- Phase transition remains sharp even with BI
- BI merely shifts the critical Rmax without smoothing

---

## 4. Phase 2: Static Hopfion Source (Topological Confinement)

### Goal

Replace the oscillating standing wave (l=1 mode) with a static Hopfion
magnetic field from V11. Test whether the Yukawa metric can self-consistently
confine the Hopfion without an Rmax box.

### Why This is Different from Phase 1

Phase 1 uses spherical standing waves (angular momentum l >= 1 provides a
centrifugal barrier). The Rmax boundary is still needed to reflect the wave.

Phase 2 uses the Hopfion topology as the confinement mechanism. The Hopf field
has linked magnetic field lines that cannot separate. In flat space (V11), the
linked structure dilates at c. With Yukawa gravity, the dilation is opposed by
the attractive metric well.

The question: is the Yukawa well strong enough to prevent Hopfion dilation?

### Equations

**Hopfion B field (from V11, pure magnetic, static):**
```
B_i(x,y,z) = C_0 * n . (dn/dx_j x dn/dx_k)     (area form pullback)
```
where n: R^3 -> S^2 is the Hopf map with parameter a (core size).

**Spherically-averaged BI energy density:**
```
rho_hop(r) = (1/4pi) int_0^{2pi} dphi int_0^pi sin(theta) d(theta)
             b^2 (sqrt(1 + |B(r,theta,phi)|^2 / b^2) - 1)
```

This integral is computed numerically (Gauss-Legendre on the sphere, using
the analytical B field from V11).

**Spherically-averaged BI pressures:**

The BI stress tensor for a pure magnetic field:
```
T^i_j = (B^i B_j) / (b^2 Gamma) - delta^i_j b^2(Gamma - 1)
Gamma = sqrt(1 + B^2/b^2)
```

Radial pressure (spherical average):
```
p_r(r) = (1/4pi) int [B_r^2 / (b^2 Gamma) - b^2(Gamma-1)] sin(theta) d(theta) d(phi)
```

Tangential pressure:
```
p_t(r) = (1/4pi) int [(B_theta^2 + B_phi^2)/(2 b^2 Gamma) - b^2(Gamma-1)] dOmega
```

**Yukawa-TOV with anisotropic pressure (generalization of V9):**
```
M_eff'' - mu^2 M_eff = -kappa r rho_hop(r)    (mass equation, sourced by energy)

f(r) = 1 / (1 + 2 M_eff / r)                  (Pade metric)
```

For the static case, there is no eigenvalue to solve. The metric is determined
entirely by the energy density profile of the Hopfion.

**Self-consistency requirement:**

In curved spacetime, the Hopfion field satisfies the BI field equations.
The metric modifies the field through the covariant derivative. For a static
magnetic field in a static spherical metric:

```
nabla_mu (B^{mu nu} / Gamma) = 0
```

This is a nonlinear elliptic PDE for the magnetic field in the background
metric f(r). The self-consistent solution requires iterating:

```
1. Start with flat-space Hopfion B(r,theta,phi)
2. Compute rho_hop(r), p_r(r)
3. Solve Yukawa BVP -> M_eff(r) -> f(r)
4. Solve BI field equation in curved metric f(r) -> new B
5. Repeat until convergence
```

Step 4 is the hard part: solving the 3D nonlinear BI equations in a spherical
background metric. This can be done with the V11 RK4 code modified to include
metric factors in the curl operator.

However, for a FIRST TEST, skip step 4: just compute the metric from the
flat-space Hopfion and check whether the resulting well depth and energy are
in the right ballpark.

### Concrete Implementation

**Step 2a**: Modify V11's `init_hopfion_B()` to compute the spherical average
of the BI energy density and store as rho(r) on a 1D radial grid.

**Step 2b**: Feed rho(r) into V9's `solve_yukawa_bvp()` at physical parameters
(kappa=35, mu=6.47).

**Step 2c**: Read off f_min, well_depth, E_total = 4pi int r^2 rho(r) dr.

**Step 2d**: Check self-consistency: is the metric correction small enough
that the flat-space Hopfion is a good approximation?

### Hopfion Energy and Size

From V11 (a=1, b=1, C_0 = b*a^2/16):
```
E_total = 0.584 code units (V11 code units, not V9 code units!)
```

Unit conversion needed. V11 uses natural units with c=1, no physical scale.
V9 uses fm and GeV. The Hopfion parameter a sets the core size:

```
a = r_proton = 0.88 fm          (proton charge radius)
b = Lambda_QCD^2 / (hbar c)     (BI field strength in 1/fm^2 units)
```

The energy integral scales as:
```
E_hop = b^2 a^3 * I(b*a^2)
```
where I is a dimensionless integral depending on the BI nonlinearity parameter
b*a^2. For b*a^2 >> 1 (strong BI): I ~ const. For b*a^2 << 1 (weak BI):
I ~ (b*a^2)^2 and E_hop ~ b^4 a^7 (Maxwell scaling).

---

## 5. Phase 3: Full 3D BI-Hopfion in Yukawa Metric

### Goal

Full 3D evolution of the BI Hopfion in a self-consistent Yukawa metric.
No spherical averaging. The metric is computed from the full 3D stress tensor.

This is the definitive test but is computationally expensive (3D metric + 3D BI
evolution). Only pursue if Phase 1 or Phase 2 shows positive results.

### Architecture

Combine V11 (3D BI evolution) with V9 (Yukawa metric):

```
1. Initialize: Hopfion B field on 3D grid
2. Compute T^{mu nu}_BI at each grid point
3. Solve 3D Yukawa equation: (nabla^2 - mu^2) Phi = -kappa T^00_BI
   (via FFT or multigrid on 3D grid)
4. Compute lapse: alpha = 1/sqrt(1 - 2Phi) (Pade, always > 0)
5. Evolve BI field one timestep with ADM metric factors (curl(alpha E))
6. Monitor div(B), div(D) constraints
7. Repeat 2-6
```

The 3D Yukawa equation can be solved efficiently with FFT:
```
Phi(k) = -kappa T^00(k) / (k^2 + mu^2)
```
where T^00(k) is the Fourier transform of the energy density.

**Analytical BI Constitutive Relations (D-field — NO Newton Solver):**

Because the state variables are (D, B), the constitutive fields E and H
needed for the ADM evolution are obtained by differentiating the Hamiltonian:

```
H(D, B) = b^2 (sqrt(R) - 1)
R = 1 + (D^2 + B^2)/b^2 + |D x B|^2/b^4

E = dH/dD = (D + B x (D x B)/b^2) / sqrt(R)

H = dH/dB = (B + D x (B x D)/b^2) / sqrt(R)
```

These are EXACT, ALGEBRAIC formulas — no iteration required. At each grid
point, the RHS evaluation computes R (one sqrt), two triple cross products,
and two divisions. This completely eliminates the Newton-Raphson loop that
V11 required (8-10 iterations per point per RHS evaluation).

Verification:
- Weak field (D, B << b): R -> 1, E -> D, H -> B (Maxwell limit)
- Strong D (|D| >> b): E -> D/|D| * b (caps at b, as required)
- D x B = 0 (parallel fields): E = D/sqrt(1 + (D^2+B^2)/b^2), same as
  the simplified formula without cross-term

The triple cross product expands via BAC-CAB:
```
B x (D x B) = D(B.B) - B(B.D) = D B^2 - B(B.D)
D x (B x D) = B(D.D) - D(D.B) = B D^2 - D(D.B)
```

So the constitutive relations in component form:
```
E_i = (D_i + (D_i B^2 - B_i (B.D))/b^2) / sqrt(R)
    = (D_i (1 + B^2/b^2) - B_i (B.D)/b^2) / sqrt(R)

H_i = (B_i + (B_i D^2 - D_i (D.B))/b^2) / sqrt(R)
    = (B_i (1 + D^2/b^2) - D_i (D.B)/b^2) / sqrt(R)
```

**Metric factors in BI evolution (ADM 3+1 formalism):**

In the 3+1 ADM decomposition with lapse alpha = sqrt(-g_00) and flat spatial
metric (weak-field: g_ij = delta_ij), Maxwell's equations take the form:

```
partial_t B = -curl(alpha E)           (Faraday)
partial_t D = +curl(alpha H)           (Ampere)
```

where alpha = 1/sqrt(1 - 2Phi) is the Pade lapse (strictly positive for all
Phi, see Safeguard A below), and D, H are the BI constitutive fields. The
lapse goes INSIDE the curl operator.

**Why this is critical:** By the product rule:
```
curl(alpha E) = alpha curl(E) + grad(alpha) x E
```

The second term, grad(alpha) x E = grad(Phi) x E, is the gravitational
lensing force. It bends the light trajectory toward the center of the well
(where grad(Phi) points inward). If the lapse is placed OUTSIDE the curl
(i.e., alpha * curl(E)), this term vanishes and the light experiences time
dilation but NO trajectory bending. The Hopfion would dilate outward at a
slower clock rate, completely ignoring the gravitational binding force.

The RK4 compute_rhs() function must therefore compute:
```
rhs_B[i] = -curl_at(alpha * E, i, j, k)     (NOT: -alpha * curl_at(E, ...))
rhs_D[i] = +curl_at(alpha * H, i, j, k)     (NOT: +alpha * curl_at(H, ...))
```

This requires forming the products alpha*E_x, alpha*E_y, alpha*E_z (and
similarly alpha*H) on the grid BEFORE taking the curl. The curl then acts
on these product fields.

These metric corrections are O(Phi) ~ O(well_depth). For the proton, well_depth
~ 30-50%, so these corrections are significant and must be included. The
grad(Phi) x E lensing term is the ENTIRE mechanism by which gravity confines
the Hopfion.

### Computational Cost (Revised — No Newton Solver)

N=128, 3D grid: 2M points.
Per RHS evaluation:
  - Analytical (D,B)->E,H: O(N^3) ~ 2M points x ~30 flops = 6 x 10^7 flops
  - Form alpha*E, alpha*H: O(N^3) ~ 2M x 6 = 1.2 x 10^7 flops
  - Two curls: O(N^3) ~ 2M x 2 x 30 = 1.2 x 10^8 flops
  - No Newton iteration (eliminated by analytical inversion)
Per timestep: 4 RK4 stages x above ~ 8 x 10^8 flops
FFT for metric update: O(N^3 log N) ~ 10^7 (done every M timesteps, amortized)
Steps for T=10: ~500.
Total: ~4 x 10^11 flops ~ tens of seconds on 8 cores.

Compared to the original estimate with Newton (10 iterations per point per
RHS evaluation), this is approximately **10x faster**. The Newton solver was
the dominant cost; eliminating it makes Phase 3 computationally trivial.

### Numerical Safeguards

#### Safeguard A: Positive-Definite Lapse Function

The weak-field lapse alpha = sqrt(1 + 2Phi) becomes imaginary if Phi < -0.5.
During the initial transient of Phase 3, before the Hopfion reaches equilibrium,
violent compression could cause the FFT Yukawa solver to output Phi < -0.5.
This would instantly crash the RK4 solver with NaN.

The 1D code (V9) already solved this problem: the Pade metric f = 1/(1+2M/r)
is strictly positive for all M. Apply the same strategy to the 3D lapse:

```
alpha = 1 / sqrt(1 - 2Phi)        (Pade lapse)
```

Properties:
- Weak field (|Phi| << 1): alpha ~ 1 + Phi + O(Phi^2) (same as sqrt(1+2Phi))
- Deep well (Phi -> -infinity): alpha -> 0+ (approaches zero, NEVER crosses)
- alpha > 0 for ALL Phi < 0.5 (and Phi > 0 always safe)
- Recovers Schwarzschild lapse alpha = 1/sqrt(1-2M/r) at large r

Alternative: exponential lapse alpha = exp(Phi). Also strictly positive,
same weak-field limit. Either choice is safe; Pade is preferred for
consistency with the 1D V9 code.

Implementation: replace sqrt(1.0 + 2.0*Phi) with 1.0/sqrt(1.0 - 2.0*Phi)
everywhere in compute_rhs(). Single line change, zero performance cost.

#### Safeguard B: Divergence Constraint Monitoring

The ADM equations analytically preserve div(B) = 0 and div(D) = 0:
```
partial_t div(B) = -div(curl(alpha E)) = 0    (exact)
partial_t div(D) = +div(curl(alpha H)) = 0    (exact)
```

On a collocated grid with 4th-order finite differences, truncation errors
will cause both constraints to drift. This drift is harmless if it stays
small (< 1% of field gradients), but catastrophic if it grows — it creates
magnetic monopoles that destroy the Hopf topology.

The code MUST output at every diagnostic step:
```
max_divB = max |div(B)|    over all grid points
max_divD = max |div(D)|    over all grid points
```

where div is computed with the same 4th-order stencil used for the curl.

Acceptance criterion: max_divB / max|curl(B)| < 0.01 for the full run.
If violated: reduce dt (CFL safety factor from 0.4 to 0.2), or add
hyperbolic divergence cleaning (Dedner et al.), or switch to a staggered
Yee spatial grid (preserves div exactly, but complicates the curl(alpha E)
product).

For V11 (N=128, T=10), divB stayed at ~5e-3 (4th-order truncation of
initial data). This was acceptable. Phase 3 adds metric factors that may
amplify the drift, so monitoring is essential.

---

## 6. Caveats and Potential Failure Modes

### Caveat 1: BI Does NOT Cap Magnetic Energy Density

BI caps the electric field (E < b) but NOT the magnetic field. The magnetic
energy density b^2(sqrt(1+B^2/b^2)-1) ~ b|B| for B >> b, which is unbounded.

However, the 1/A runaway in V9/V10 is driven by the ELECTRIC component
(omega^2 u^2/A), not the magnetic component (A(u')^2). The blueshift factor
1/sqrt(A) amplifies E but suppresses B. So the BI cap on E is the relevant
mechanism.

### Caveat 2: Static Hopfion Has No Eigenfrequency

The V9 Regge-Wheeler framework assigns energy via omega^2 (mode frequency).
A static Hopfion has omega = 0. Its energy comes entirely from the magnetic
field:

```
E_hop = int b^2 (sqrt(1 + B^2/b^2) - 1) d^3x
```

For the Hopfion to match the proton mass (938 MeV), we need:
```
b^2 a^3 * I(b*a^2) = 938 MeV * a / (hbar c)
```

This constrains the (a, b) parameter space.

### Caveat 3: Hopfion Topology is NOT Conserved in BI

Unlike ideal MHD, Born-Infeld electrodynamics does not exactly conserve the
Hopf invariant. The linked field lines can reconnect through BI nonlinearity.
V11 showed that the field disperses even in flat space.

However, V11 also showed that the Hopf structure persists for t ~ a/c (one
light-crossing time). If the Yukawa metric provides binding on this timescale,
the topology may be effectively conserved.

### Caveat 4: Spherical Average May Miss the Physics

The Hopfion is toroidal, not spherical. Spherical averaging (Phase 2) loses
the directional structure. Phase 3 (full 3D) is needed for a definitive answer.

However, for the Q9 question (does BI prevent the 1/A runaway?), spherical
averaging is sufficient because the runaway is a radial phenomenon.

---

## 7. Strategic Reassessment

### Phase 1 Will Likely Fail to Produce a Stable Particle

Even with the D-field BI fix preventing NaN crashes, Phase 1 still uses l=1
spherical standing waves. The centrifugal barrier for spherical waves decays
as 1/r^2, while the Yukawa well decays as e^{-mu r}/r. At large r, the wave
will ALWAYS tunnel out. BI nonlinearity alters the core physics but does
absolutely nothing to the far-field barrier. Phase 1 will yield a smoothed-out
resonance (proving the concept that BI prevents horizon collapse), but it will
still leak.

### Phase 3 is the ONLY Path to a Stable Geon

The Hopfion topology provides what spherical standing waves cannot: confinement
without a box. The linked magnetic field lines cannot separate topologically.
The Yukawa metric provides the inward binding force. The BI nonlinearity prevents
the core from collapsing to a singularity. ALL THREE are needed simultaneously.

### Phase 2 is Skippable

A static snapshot of a dynamic knot tells nothing about stability. The Hopfion
is a time-evolving configuration — its energy, pressure, and topology are all
dynamic quantities. Spherically averaging a toroidal field and feeding it into
a 1D solver loses the directional structure that makes the Hopfion a Hopfion.

---

## 8. Implementation Order (Optimized)

### Step 1: Phase 1 Level 1 (proof of concept, ~5 minutes)

Modify V9's `compute_rho()` (already has b_BI!) to run at physical coupling:
```
strong_geon -hartree -kappa 35 -mu 6.47 -l 1 -Rmax 2.85 -b_BI 1.0
```

Scan b_BI = {0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0} at Rmax = 2.85.
This takes minutes (V9 is a 1D solver).

**Key observable**: does f_min stay above 0.1 with E_total > 800 MeV?

This proves that a nonlinear energy cap stops the 1/A horizon collapse.
The result will be a leaky resonance, not a stable particle — that is expected.

### Step 2: Phase 3 (the final test)

Skip Phase 1 Level 2 and Phase 2. Go directly to the 3D BI-Hopfion in Yukawa:

1. Build the 3D FFT Yukawa solver (Section 5).
2. Update the V11 RK4 code with the CORRECT curl(alpha E) ADM operators.
3. Use the D-field BI constitutive relations with ANALYTICAL inversion
   (E = dH/dD, H = dH/dB — no Newton solver, ~10x faster than V11).
4. Include the |D x B|^2/b^4 Poynting cross-term in the energy density.
5. Inject the Hopfion into the strong gravity metric.

If the inward pull of the Yukawa metric (grad(alpha)) is stronger than the
topological dilation pressure of the Hopfion, the knot will compress. As it
compresses, E and B approach the BI limit b. The nonlinearity stiffens the
fields, preventing collapse to a singularity.

If this balances: a stable, finite-energy, 3D particle made of self-interacting
light in a self-consistent gravitational well.

**Key observables**:
- R_eff: stable (not growing, not collapsing)
- E_total: conserved to < 1%
- Hopf invariant: approximately conserved
- f_min: finite (> 0.1), no horizon
- max_divB / max|curl(B)| < 0.01 (constraint preservation)
- alpha_min > 0 (lapse stays positive, no coordinate singularity)

---

## 9. Equations Summary

### BI Energy Density — Two Formulations

**E-field formulation (from V11, for diagnostics only — DIVERGES as E -> b):**
```
u = (B^2 Y + E^2 + (E.B)^2/b^2) / (Gamma (Y + Gamma))
Y = 1 + B^2/b^2
Gamma = sqrt(1 + (B^2 - E^2)/b^2 - (E.B)^2/b^4)
```
WARNING: Gamma -> 0 as E -> b, causing u -> infinity.

**D-field formulation (for evolution — ALWAYS finite):**
```
rho_BI = b^2 (sqrt(R) - 1)
R = 1 + (D^2 + B^2)/b^2 + |D x B|^2/b^4
```
Valid for ALL field strengths. R >= 1 always. No NaN possible.
The |D x B|^2 Poynting cross-term vanishes for parallel D,B (Phase 1)
but is essential for 3D configurations with swirling fields (Phase 3).

### Analytical BI Constitutive Relations (D-field)

```
E = dH/dD = (D(1 + B^2/b^2) - B(B.D)/b^2) / sqrt(R)
H = dH/dB = (B(1 + D^2/b^2) - D(D.B)/b^2) / sqrt(R)
```

These are exact algebraic formulas — NO Newton solver needed.
Weak-field limit: E -> D, H -> B (Maxwell). Strong-field: |E| < b (capped).

### Yukawa BVP (from V9)

```
M'' - mu^2 M = -kappa r rho(r)
M(0) = 0,  M'(Rmax) = -mu M(Rmax)
f(r) = 1/(1 + 2M(r)/r)     (Pade)
```

### Regge-Wheeler Eigenvalue (from V9)

```
-f^2 v'' + W(r) v = omega^2 v
v = sqrt(f) u,  u normalized (int u^2 dr = 1)
W = f l(l+1)/r^2 + f f''/2 - (f')^2/4
```

No BI constraint on u amplitude. The D-field formulation handles
arbitrarily deep wells without requiring amplitude reduction.

### BI Standing-Wave Energy Density (Level 2, D-field)

```
D_0^2(r) = N_l omega^2 u^2 / (8 pi r^2 f(r))
B_0^2(r) = N_l f(r) (u')^2 / (8 pi r^2)
alpha(r) = D_0^2/b^2,  beta(r) = B_0^2/b^2

<rho_BI>(r) = (b^2/pi) int_0^pi [sqrt(1 + alpha cos^2(t) + beta sin^2(t)) - 1] dt
```

Properties:
- Integrand >= 0 for all alpha, beta >= 0 (no sign issues)
- Weak field: <rho_BI> = (D_0^2 + B_0^2)/4 = rho_Maxwell
- Strong D: <rho_BI> ~ b D_0 / pi (linear, not divergent)

### 3D Yukawa via FFT (Phase 3)

```
Phi(k) = -kappa rho_BI(k) / (k^2 + mu^2)

where rho_BI(k) = FFT[T^00_BI(x)]
      Phi(x) = FFT^{-1}[Phi(k)]
```

Metric: g_00 = -(1 + 2Phi(x)),  g_ij = delta_ij  (weak-field limit)
Lapse: alpha = 1/sqrt(1 - 2Phi)  (Pade, always > 0, ~ 1 + Phi for |Phi| << 1)

### ADM Evolution Equations (Phase 3)

```
partial_t B = -curl(alpha E)           (Faraday)
partial_t D = +curl(alpha H)           (Ampere)

curl(alpha E) = alpha curl(E) + grad(alpha) x E
                \_____________/   \_____________/
                time dilation     gravitational lensing
```

The grad(alpha) x E term bends light inward toward the metric well.
Implementation: form alpha*E on the grid, THEN take the curl.

---

## 10. Implementation Traceability (IMPLEMENTATION.md)

Once any equation from this plan is implemented in code, create and maintain
an `IMPLEMENTATION.md` file in the same directory as the source. This file
provides a line-by-line traceability map between the mathematical equations
and the code that implements them.

### Purpose

V8-V11 accumulated 17 critical bugs (see `critical_bugs.md`). Most were
silent sign errors, wrong normalizations, or misplaced operators that
produced plausible-looking but incorrect results. The root cause in every
case was the same: the equation was written in the plan, the code was written
separately, and the two drifted apart without anyone noticing.

IMPLEMENTATION.md eliminates this failure mode by forcing every equation to
be paired with its verbatim code, side by side, with an explicit description
of how one maps to the other.

### Format

For each equation implemented in code:

1. **State the equation** in standard math notation (the "what").
2. **Quote the code verbatim** that implements it (the "how").
3. **Describe the mapping**: variable names, index conventions, unit
   conversions, numerical method choices, and any non-obvious correspondences
   between the math and the code.

Example entry:

```
### D-field BI Energy Density

**Equation:**
  R = 1 + (D^2 + B^2)/b^2 + |D x B|^2/b^4
  rho_BI = b^2 (sqrt(R) - 1)

**Code** (bi_geon.c:142-150):
  double D2 = Dx*Dx + Dy*Dy + Dz*Dz;
  double B2 = Bx*Bx + By*By + Bz*Bz;
  double DxBx = Dy*Bz - Dz*By;
  double DxBy = Dz*Bx - Dx*Bz;
  double DxBz = Dx*By - Dy*Bx;
  double P2 = DxBx*DxBx + DxBy*DxBy + DxBz*DxBz;
  double R = 1.0 + (D2 + B2)/(b*b) + P2/(b*b*b*b);
  double rho = b*b * (sqrt(R) - 1.0);

**Mapping:**
  - D^2 = |D|^2  →  D2 = Dx*Dx + Dy*Dy + Dz*Dz
  - |D x B|^2  →  P2 (computed via explicit cross product components)
  - b is the BI field strength parameter (global `b_param`)
  - R >= 1 always (sum of 1 + non-negative terms), so rho >= 0 (no NaN)
  - Poynting cross-term P2/b^4 is negligible for parallel D,B (Phase 1)
    but essential for 3D knotted fields (Phase 3)
```

### Rules

1. **Every equation gets an entry.** If it appears in the plan, it appears
   in IMPLEMENTATION.md with its code. No exceptions.

2. **Code is quoted verbatim.** Copy-paste from the source file, with the
   file name and line numbers. Do not paraphrase or simplify the code.

3. **When fixing a bug, fix BOTH files.** Update the code, then update the
   corresponding IMPLEMENTATION.md entry. If the two are out of sync, the
   bug is not fixed.

4. **Include numerical method details.** If the equation is continuous but
   the code uses a discrete approximation (finite differences, quadrature,
   Newton iteration), state the method and its order of accuracy.

5. **Include failure modes.** If the equation has domain restrictions
   (e.g., sqrt of negative, division by zero), state how the code handles
   them and whether the implementation is safe.

---

## 11. Document Lifecycle

### PLAN.md is Frozen at Implementation Start

Once coding begins, PLAN.md is NEVER modified again. It is the permanent
record of what was intended, what equations were chosen, and why. If a bug
is found in the plan's equations during implementation, the fix goes into
IMPLEMENTATION.md and RESULTS.md — not back into PLAN.md. This preserves
the audit trail: you can always diff the plan against the implementation
to see what changed and why.

### RESULTS.md — Living Record of Execution

Create `RESULTS.md` when implementation begins. It is updated throughout
execution, phase by phase:

**When starting a new phase**, add a section that:
1. States the phase name and goal (referencing PLAN.md section numbers).
2. Lists the specific equations being implemented (by PLAN.md equation labels).
3. States the success/failure criteria (copied verbatim from PLAN.md).
4. Lists the parameter values being tested.
5. States what diagnostics will be checked and what "correct" looks like.

**When results come back**, update the same section with:
1. Raw numerical output (tables, key values).
2. Whether each success criterion was met (PASS/FAIL with numbers).
3. Any unexpected behavior or deviations from the plan.
4. Bug fixes applied during the phase (with IMPLEMENTATION.md cross-refs).
5. Conclusion for this phase: what was learned, what it means for next phase.

Example structure:

```
# V12 Results

## Phase 1 Level 1: Scalar BI Hartree (PLAN.md Section 3)

### Setup
- Equations: PLAN.md Section 3, Level 1 scalar BI
- Parameters: kappa=35, mu=6.47, Rmax=2.85
- Success criteria: f_min > 0.1, E_total in [800, 1000] MeV (PLAN.md §3)
- b_BI scan: {0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0}

### Results
| b_BI | f_min | E_total (MeV) | omega^2 | converged |
|------|-------|---------------|---------|-----------|
| ...  | ...   | ...           | ...     | ...       |

### Assessment
- [PASS/FAIL]: f_min criterion — ...
- [PASS/FAIL]: E_total criterion — ...
- Conclusion: ...

---

## Phase 3: 3D BI-Hopfion in Yukawa (PLAN.md Section 5)
...
```

### CONCLUSION.md — Final Synthesis

Written ONCE, after all phases are complete. Contains three parts:

1. **Theory**: The conceptual argument (why V9 Yukawa + V11 BI + Hopf
   topology is the only remaining path). Summarizes PLAN.md Sections 1-7
   in narrative form, not equation-heavy.

2. **Method**: What was actually implemented and run. References
   IMPLEMENTATION.md for equation-code traceability. States any deviations
   from the original plan and why.

3. **Results**: The bottom line. Did V12 produce a stable geon? Summarizes
   RESULTS.md into a single verdict with the key numbers. States what was
   learned regardless of outcome.

### File Summary

```
v12/
  PLAN.md              — FROZEN at implementation start. The blueprint.
  IMPLEMENTATION.md    — Equation-to-code traceability map. Updated with code.
  RESULTS.md           — Phase-by-phase execution log. Updated with results.
  CONCLUSION.md        — Final synthesis. Written once at the end.
  src/                 — Source code.
  data/                — Output data files.
```
