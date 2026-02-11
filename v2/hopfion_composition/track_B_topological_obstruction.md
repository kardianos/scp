# Track B: Topological Obstruction Analysis

**Question**: Does the topology of a B=1 Skyrmion FORCE the degenerate sector
d of the Cl+(3,0,1) field to be nonzero?

**Short answer**: No. There is no topological obstruction to d = 0.
The degenerate sector is not forced to be nonzero by any purely geometric or
topological argument within the current framework. The result is a definitive
negative: d = 0 is a valid, smooth, globally well-defined field configuration
on the B=1 hedgehog background.

---

## 1. B1: Is d = 0 Allowed? (Explicit Check)

### Setup

The field Psi in Cl+(3,0,1) decomposes as Psi = q + epsilon*d where:
- q = s + f_1 e_{23} + f_2 e_{31} + f_3 e_{12}  (bulk quaternion)
- d = j_1 e_{01} + j_2 e_{02} + j_3 e_{03} + tau e_{0123}  (degenerate)
- epsilon = e_0, epsilon^2 = 0

The hedgehog B=1 ansatz is:

    q = rho_0 (cos f(r) + sin f(r) * r_hat . sigma)

where f(0) = pi, f(infinity) = 0, and sigma = (e_{23}, e_{31}, e_{12}).

### Sigma-model constraint

The constraint |Psi|^2 = <Psi Psi~>_0 = rho_0^2 gives:

    <(q + epsilon*d)(q~ + epsilon*d~)>_0 = |q|^2 + epsilon*(q d~ + d q~)|_0 = rho_0^2

Since epsilon^2 = 0, we need:
1. |q|^2 = rho_0^2  (automatically satisfied by the hedgehog)
2. <q d~ + d q~>_0 = 0  (the sigma-model constraint on d)

On the hedgehog, this second condition reduces (on the z-axis by rotational
symmetry) to:

    tau = j_r * tan(f)

where j_r is the radial component of the J-vector.

### Check: d = 0

Setting d = 0 (i.e., j_1 = j_2 = j_3 = tau = 0):
- Condition 1: |q|^2 = rho_0^2. SATISFIED (hedgehog has unit norm).
- Condition 2: <q*0 + 0*q~>_0 = 0. TRIVIALLY SATISFIED.

### Euler-Lagrange equations

The current Lagrangian is:

    L = L_2 + L_4 - V - V_D

where:
- L_2 = (1/2) eta^{mu nu} <d_mu Psi d_nu Psi~>_0 = (1/2)|d_mu q|^2
  (degenerate sector drops out because e_0^2 = 0)
- L_4 = Skyrme term, depends only on q (same reason)
- V = (lambda/4)(|q|^2 - rho_0^2)^2, depends only on q
- V_D = (mu^2/2)(|J|^2 + tau^2), depends only on d

The Euler-Lagrange equation for d is:

    delta L / delta d = -mu^2 d = 0

Since mu^2 > 0, the unique solution is d = 0.

Even if mu = 0 (no degenerate mass term), the EL equation for d becomes
0 = 0 (no term in L depends on d or its derivatives). In this case d is
completely undetermined by the action — it is a flat direction. d = 0 is one
valid solution among infinitely many.

### Smoothness check

At every point in R^3:
- q is smooth and nonzero (|q| = rho_0 > 0 everywhere)
- d = 0 is obviously smooth

There is no singularity, no discontinuity, no obstruction. The configuration
Psi = q + 0 is a perfectly valid smooth Cl+(3,0,1) field.

### Completeness of Cl+(3,0,1)

One might worry that "using only 4 of 8 components" somehow breaks the algebraic
structure. It does not:

- Cl+(3,0,1) is an 8-dimensional ALGEBRA. The element q (with d = 0) is a
  perfectly valid element of this algebra — specifically, it lies in the
  subalgebra Cl+(3,0) subset Cl+(3,0,1).

- The even subalgebra Cl+(3,0) = {1, e_{23}, e_{31}, e_{12}} is isomorphic
  to the quaternions H. It is CLOSED under multiplication, addition, and
  inversion (for nonzero elements). A field valued in H subset DH = H + epsilon*H
  is algebraically self-consistent.

- Setting d = 0 is analogous to having a real number r in C (complex numbers).
  There is no topological or algebraic reason why r cannot be real — the imaginary
  part being zero is perfectly consistent with all properties of C.

### Verdict on B1

**d = 0 IS allowed.** It satisfies:
- The sigma-model constraint
- The Euler-Lagrange equations (for any mu >= 0)
- Smoothness requirements
- Algebraic consistency of Cl+(3,0,1)

The hedgehog Skyrmion with d = 0 is the standard Skyrme model solution embedded
in the larger algebra. There is nothing inconsistent about this embedding.

---

## 2. B2: Fiber Bundle Analysis

### The bundle structure

At each spatial point x in R^3, the field value Psi(x) lives in Cl+(3,0,1).
The sigma-model constraint |q|^2 = rho_0^2 restricts q to S^3 (the 3-sphere
of radius rho_0 in quaternion space). The constraint on d is:

    <q d~ + d q~>_0 = 0

This defines, for each q in S^3, a LINEAR subspace of the degenerate sector.
Specifically, if we write xi = q^{-1} d (the "body-frame" degenerate field),
then the constraint becomes:

    <xi~ + xi>_0 = 0

i.e., the scalar part of xi + xi~ must vanish. Since xi lives in Cl-(3,0)
(odd part of the spatial Clifford algebra, spanned by {e_1, e_2, e_3, e_{123}}),
its reversal is xi~ = xi for the vector parts and xi~ = -xi for the e_{123} part.
Wait — let me be more careful.

Actually, d = j_1 e_{01} + j_2 e_{02} + j_3 e_{03} + tau e_{0123}, so
q^{-1} d is an element of the form e_0 * (something in Cl-(3,0)).

Let me parametrize differently. Write Psi = q(1 + epsilon * xi) where
xi = q^{-1} d / (doesn't quite work because of the epsilon structure).

More carefully: Psi = q + epsilon * d. The reversal is:
Psi~ = q~ + epsilon * d~ (since epsilon commutes past everything for reversal
purposes — epsilon = e_0 which is grade 1, but in the even subalgebra the
products like e_{01} reverse to -e_{01}).

Actually, let us simply work out the constraint directly on the hedgehog.

### The constraint on the hedgehog (full 3D)

The hedgehog is q = rho_0(cos f + sin f * n_hat . sigma) where n_hat = r_hat.
In components at an arbitrary point (r, theta, phi):

    s = rho_0 cos f(r)
    f_1 = rho_0 sin f(r) * sin(theta) sin(phi)     [coefficient of e_{23}]
    f_2 = rho_0 sin f(r) * (-sin(theta) cos(phi))   [coefficient of e_{31}]
    f_3 = rho_0 sin f(r) * cos(theta)                [coefficient of e_{12}]

(The exact angular dependence follows from the hedgehog map identifying spatial
directions with isospin directions.)

The sigma-model constraint q d~ + d q~ = 0 (at grade 0) gives, after expanding
the Clifford products:

    2(s*tau + f_1*j_1 + f_2*j_2 + f_3*j_3) = 0

i.e.:
    s*tau + f . J = 0

where f = (f_1, f_2, f_3) and J = (j_1, j_2, j_3). On the hedgehog:

    rho_0 cos(f) * tau + rho_0 sin(f) * (r_hat . J) = 0

    tau = -tan(f) * (r_hat . J)

(Note the sign: tau = -j_r tan(f), where j_r = r_hat . J is the radial component
of J. The sign depends on conventions; the key point is that this is one
equation for four unknowns.)

### Fiber at each point

At each spatial point, the constraint removes one degree of freedom from the
degenerate sector (4D -> 3D). The fiber of allowed d values is a 3-dimensional
LINEAR subspace, depending smoothly on q.

**Key observation**: The zero section d = 0 is always in this fiber (0 satisfies
any linear constraint). Therefore there is always a smooth global section of
this bundle that is identically zero.

### Does the fiber twist?

For the fiber to FORCE d nonzero, we would need the bundle to have no global
nonvanishing section — equivalently, the zero section would need to be
"topologically isolated" in some sense. But here the fiber is a 3-dimensional
VECTOR SPACE (a linear subspace of R^4). A vector bundle over a contractible
base (R^3) is always trivial. Even over a non-contractible base, the relevant
obstruction would live in characteristic classes.

But we can check directly:

**At r = 0**: f = pi, cos f = -1, sin f = 0.
Constraint: -rho_0 * tau + 0 = 0, so tau(0) = 0.
The J-vector (j_1, j_2, j_3) is unconstrained. Fiber is R^3 (free J, tau = 0).

**At r = r_{1/2}** (where f = pi/2): cos f = 0, sin f = 1.
Constraint: 0 + rho_0 * j_r = 0, so j_r(r_{1/2}) = 0.
Fiber: j_r = 0, tau free, tangential J free. Still 3-dimensional.

**At r -> infinity**: f -> 0, sin f -> 0, cos f -> 1.
Constraint: rho_0 * tau + 0 = 0, so tau -> 0.
Fiber: tau = 0, J free. Still R^3.

The fiber is always 3-dimensional, and d = 0 is always in it. As q varies
smoothly, the constraint equation changes smoothly but never excludes d = 0.

### Explicit smooth nonzero solutions

For completeness, here are smooth d != 0 solutions:

**Example 1**: j_r(r) = A sin^2(f(r)), j_theta = j_phi = 0.
Then tau = -A sin^2(f) * tan(f) = -A sin(f) cos(f) * (sin f / cos f) * tan f...

Let me redo this carefully. With j_r(r) = A sin^2(f), the constraint gives:
    tau = -tan(f) * j_r = -A sin^2(f) tan(f) = -A sin^2(f) sin(f)/cos(f)

At f = pi/2 (where cos f = 0), j_r = A sin^2(pi/2) = A, but tau = -A/0 = divergent.

So this is NOT smooth. The requirement that tau be smooth near f = pi/2 forces
j_r to vanish there fast enough. Specifically:

    tau = -j_r * sin(f)/cos(f)

is smooth iff j_r vanishes at least as fast as cos(f) near the zero of cos(f).

**Example 2**: j_r(r) = A cos(f(r)), j_theta = j_phi = 0.
Then tau = -A cos(f) * tan(f) = -A sin(f).
This is smooth everywhere! At f = pi/2: j_r = 0, tau = -A, both finite.

**Example 3**: j_r(r) = A cos(f(r)) f'(r), tau = -A sin(f) f'(r).
Also smooth (f' is smooth by the ODE regularity).

So there exists a continuous family of smooth d != 0 solutions, parameterized
by the amplitude A and the functional form of the transverse components. But
d = 0 (A = 0) is always one of them.

### Bundle triviality

The constraint bundle is a rank-3 vector bundle over R^3. Since R^3 is
contractible, EVERY vector bundle over R^3 is trivial (product bundle).
This means there is always a global nonvanishing section — but there is also
always a global zero section. There is no topological obstruction to the zero
section.

Even if we compactify R^3 to S^3 (one-point compactification, which is the
natural domain for a Skyrmion with f -> 0 at infinity), a rank-3 vector bundle
over S^3 is classified by pi_2(GL(3,R)) = 0. So the bundle is still trivial,
and d = 0 is still a valid global section.

### Verdict on B2

**The fiber does NOT twist in a way that forces d != 0.** The constraint bundle
is a trivial rank-3 vector bundle, and d = 0 is a valid smooth global section.
There is an infinite family of smooth d != 0 solutions (parameterized by amplitude
and angular structure), but none is preferred or required by the topology of q.

---

## 3. B3: Dual Connection Curvature

### Setup

Define the body-frame degenerate field:

    xi = q^{-1} d

The left current decomposes as:

    L_mu = Psi^{-1} d_mu Psi = l_mu + epsilon * delta_l_mu

where:
- l_mu = q^{-1} d_mu q  (bulk left current, valued in Im(H) = su(2))
- delta_l_mu = q^{-1} d_mu d - q^{-1} d * l_mu = D_mu(xi)

Here D_mu is the covariant derivative with connection l_mu, acting on xi by:

    D_mu(xi) = d_mu(xi) + [l_mu, xi]    (or d_mu(xi) - l_mu * xi + xi * l_mu, etc.)

The exact form depends on the representation. For xi in the adjoint representation
of the quaternion algebra:

    D_mu(xi) = d_mu(xi) + l_mu * xi - xi * l_mu = d_mu(xi) + [l_mu, xi]

### Curvature of the dual connection

The curvature of D_mu acting on xi is:

    [D_mu, D_nu](xi) = [R_{mu nu}, xi]

where R_{mu nu} = d_mu l_nu - d_nu l_mu + [l_mu, l_nu] is the curvature
2-form of the bulk connection l.

For the hedgehog, R_{mu nu} is NONZERO. In fact, its integral is the
topological charge:

    B = -(1/(24 pi^2)) integral epsilon^{ijk} Tr(l_i l_j l_k) d^3x = 1

So the curvature R_{mu nu} is nontrivial — the bulk connection is a nontrivial
SU(2) gauge field.

### Does nonzero curvature force xi != 0?

No. The curvature tells us about the PARALLEL TRANSPORT of xi, not about xi
itself. The equation:

    [D_mu, D_nu](xi) = [R_{mu nu}, xi]

says that if xi != 0, then parallel-transporting it around a loop will ROTATE
it (because R_{mu nu} != 0). But if xi = 0, the equation is trivially
0 = [R_{mu nu}, 0] = 0. There is no contradiction.

An analogy: a curved manifold (nonzero Riemann tensor) does not force every
vector field to be nonzero. The zero vector field is always a valid section of
any vector bundle, curved or not. The curvature affects how NONZERO vectors
rotate under transport, not whether vectors must be nonzero.

### Flatness condition for xi = 0

The dual part of the Maurer-Cartan equation is (as derived in the problem statement):

    d(delta_l) + l ^ delta_l + delta_l ^ l = 0

For delta_l = 0 (which corresponds to xi = const, and if we additionally set
xi = 0 then delta_l = D(0) = 0):

    0 + 0 + 0 = 0

This is trivially satisfied. The Maurer-Cartan equation places NO constraint
that forces delta_l != 0.

### What if we require delta_l to be a FLAT connection?

One might ask: does the dual connection delta_l need to be flat (zero curvature)
for geometric consistency? The answer is:

The curvature of delta_l as a connection (if we view it as such) would be:

    F^{delta}_{mu nu} = D_mu(delta_l_nu) - D_nu(delta_l_mu) + [delta_l_mu, delta_l_nu]

For delta_l = 0, this is trivially zero. So xi = 0 gives a flat dual connection.

If xi != 0 but constant in the body frame (D_mu xi = 0 for all mu), then
delta_l = 0 still, and the dual curvature is still zero. But D_mu xi = 0 with
nonzero xi means xi is covariantly constant — a parallel section of the
adjoint bundle. Whether such a section exists depends on the holonomy of the
bulk connection l.

### Holonomy argument

The holonomy group of the bulk SU(2) connection l on the hedgehog is
nontrivial (it is generically all of SU(2), since the curvature R_{mu nu} spans
all of su(2) at generic points). A covariantly constant section xi (with
D_mu xi = 0) exists only if the holonomy group acts trivially on xi — i.e.,
xi commutes with all elements of SU(2). The only element of Im(H) that commutes
with all of SU(2) is xi = 0 (since su(2) has trivial center).

**This means**: There is no NONZERO covariantly constant xi on the B=1 hedgehog.
Any nonzero xi must have D_mu xi != 0 at some points.

But this does NOT force xi to be nonzero! It says:
- xi = 0 is the UNIQUE covariantly constant section.
- If xi != 0 anywhere, it cannot be covariantly constant everywhere.

The holonomy of the SU(2) connection forces nonzero xi to VARY, but it does
not force xi to be nonzero in the first place.

### Verdict on B3

**The dual connection curvature does NOT force d != 0.** The curvature of the
bulk connection is nonzero (reflecting the topology of the Skyrmion), but this
only constrains the TRANSPORT of d, not its existence. The zero section d = 0
is always compatible with all integrability conditions.

The one interesting result: the holonomy of the bulk connection rules out
covariantly constant nonzero d. Any d != 0 must vary spatially.

---

## 4. B4: Alternative Topological Arguments

### 4.1 Holonomy

The path-ordered exponential of the bulk left current l around a closed loop
gamma is:

    Hol(gamma) = P exp(oint_gamma l_mu dx^mu) in SU(2)

For a large sphere enclosing the soliton, the holonomy is nontrivial (related
to the winding number B = 1). The full left current L = l + epsilon * delta_l
has holonomy:

    Hol_full(gamma) = P exp(oint_gamma L_mu dx^mu)
                    = Hol(gamma) + epsilon * (correction from delta_l)

For delta_l = 0 (d = 0), the full holonomy is simply:

    Hol_full = Hol(gamma) + 0 = Hol(gamma)

This is a perfectly valid element of the dual quaternion group (it lies in the
SU(2) subgroup). There is no requirement for the holonomy to have a nonzero
epsilon component. The dual quaternion group DH* = H* + epsilon*H contains
H* as a subgroup; holonomies can live entirely in this subgroup.

For comparison: in SO(3,1) (the Lorentz group), a pure rotation (no boost) is
a valid group element. The boost components being zero is not an inconsistency.
Similarly, in the dual quaternion group, the translational (epsilon) components
being zero is allowed — it just means the holonomy is a pure rotation with no
translation.

**Verdict**: Holonomy does NOT force d != 0.

### 4.2 Index Theory

The Dirac operator on the Skyrmion background:

    D = gamma^mu D_mu = gamma^mu (d_mu + l_mu)

has zero modes counted by the Atiyah-Singer index theorem. For an SU(2)
connection on R^3 with instanton number k (which for the Skyrmion is related
to B = 1), the index gives:

    ind(D) = 2k (for SU(2) fundamental representation)

These zero modes live in the SPINOR space, not in the adjoint (where d lives).
The degenerate sector d transforms under the adjoint action of q (conjugation),
not under the fundamental. The relevant index for adjoint fields would be
different.

More importantly: the index theorem counts zero modes of a DIFFERENTIAL OPERATOR,
not components of the field itself. Even if there are B zero modes in some
fermionic sector, this does not force a bosonic field (d) to be nonzero.
The zero modes would tell us about linearized fluctuations around the background,
not about the background itself.

For the specific case of the degenerate sector:
- The "kinetic operator" for d is ABSENT from the Lagrangian (e_0^2 = 0
  kills it).
- There is no well-defined Dirac or Laplace operator for d coupled to the
  bulk connection.
- Therefore the index theorem does not apply in the standard way.

If we were to ADD a kinetic term L_{2,D} = (kappa^2/2)|D_mu d|^2 (with some
covariant derivative), then the zero modes of the resulting operator would be
relevant. But in the CURRENT Lagrangian, no such term exists.

**Verdict**: Index theory does NOT force d != 0 (the relevant operator is
absent from the Lagrangian).

### 4.3 Homotopy

The target space for the full field Psi = q + epsilon*d, subject to |q| = rho_0,
is:

    M = {(q, d) : q in S^3, sigma-model constraint on d}

The sigma-model constraint is: s*tau + f.J = 0, which defines a codimension-1
linear subspace of R^4 (the d-space) at each q. As q varies over S^3, this
subspace rotates, defining a rank-3 vector bundle over S^3.

The total space of this bundle is:

    E = {(q, d) : q in S^3, d in R^3_q}

where R^3_q is the fiber at q. Since d can be zero (the zero section), the
relevant question for whether d MUST be nonzero is whether the zero section
can be deformed away from zero — i.e., whether there exists a nonvanishing
section. But the question is the OPPOSITE: can the zero section be maintained?
Yes, trivially.

For completeness, let us compute the homotopy of the total space. As a rank-3
vector bundle over S^3, it is homotopy equivalent to S^3 (the zero section
is a deformation retract). Therefore:

    pi_3(M) = pi_3(S^3) = Z

The topological charge of the full field Psi maps down to the topological charge
of q alone. There is no ADDITIONAL topological invariant from the d sector
(because the fibers are contractible vector spaces).

If instead we required |d| = const != 0 (d lives on a sphere rather than a
vector space), then the fiber would be S^2 and we would have a sphere bundle
over S^3. The homotopy pi_3 of such a bundle could have additional structure.
But the sigma-model constraint does NOT require |d| to be constant — it is a
LINEAR constraint, and d = 0 is always allowed.

**Verdict**: The homotopy analysis confirms pi_3(M) = Z with no additional
topological charges from d. The d-sector does not contribute independent
topological invariants.

### 4.4 Energy Argument

Even if d = 0 is kinematically allowed, could d != 0 be energetically preferred?

In the current Lagrangian:
- E_2, E_4, E_V depend only on q (proven in spec/math/03_dynamics.md, Section 5.6)
- E_D = (mu^2/2) integral (|J|^2 + tau^2) d^3x >= 0

So d != 0 INCREASES E_D without affecting E_2, E_4, E_V. The energy is strictly
MINIMIZED at d = 0 (for any mu > 0). For mu = 0, the energy is independent
of d, so d = 0 is degenerate with all other d.

What about the PSEUDOSCALAR grade contribution? As computed in the constraint
gravity analysis (constraint_gravity_hypothesis.md, Step 1):

    <(d_mu Psi)(d_mu Psi~)>_4 = 2[d_mu s * d_mu tau
                                   - d_mu f_1 * d_mu j_1
                                   - d_mu f_2 * d_mu j_2
                                   - d_mu f_3 * d_mu j_3] * e_{0123}

If this term were included in the action with some coupling g_ps:

    L_ps = g_ps * <(d_mu Psi)(d_mu Psi~)>_4 * (e_{0123})^{-1}

then the Euler-Lagrange equation for d would involve derivatives of q as
source terms. Specifically, varying L_ps with respect to tau:

    EL for tau: 2 g_ps * nabla^2 s = 0  =>  nabla^2 tau = -nabla^2 s

This is a Poisson equation for tau sourced by nabla^2 s. But as shown in the
constraint gravity analysis, the integral of nabla^2 s over all space is ZERO
(no monopole moment), so the resulting tau field falls off as 1/r^3 or faster.

Furthermore, this term is a TOTAL DERIVATIVE in disguise. The integrand
d_mu s * d_mu tau - d_mu f_i * d_mu j_i can be rewritten (on the hedgehog
with J proportional to r_hat and tau constrained) in a form that, after
integration by parts, contributes only a boundary term. This is because the
cross kinetic term is actually the epsilon^1 expansion of |d_mu Psi|^2, which
is topological (like a Chern-Simons form in one lower dimension).

The bottom line: even including the pseudoscalar Lagrangian term, the energetics
do not prefer d != 0 over d = 0. The pseudoscalar term is either a boundary
term (does not affect EL equations) or produces a source with zero monopole
moment (not relevant for long-range behavior).

**Verdict**: Energy considerations do NOT force d != 0. In the current
Lagrangian, d = 0 strictly minimizes the energy. Even with the pseudoscalar
term, d = 0 remains an extremum.

---

## 5. B5: Geometric Consistency Conditions

### 5.1 Maurer-Cartan equation

As stated in the problem, the full Maurer-Cartan equation for L = Psi^{-1} dPsi
is automatically satisfied (it is an algebraic identity for any invertible Psi).
Decomposing:

    Order epsilon^0: dl + l ^ l = 0  (bulk Maurer-Cartan, automatic)
    Order epsilon^1: d(delta_l) + l ^ delta_l + delta_l ^ l = 0

For delta_l = 0 (d = 0), the second equation is 0 = 0. Trivially satisfied.

### 5.2 Bianchi identity for the dual curvature

If we define a dual curvature F^{dual}_{mu nu} = D_mu(delta_l_nu) - D_nu(delta_l_mu),
the Bianchi identity would be:

    D_{[rho} F^{dual}_{mu nu]} = 0

For F^{dual} = 0 (from delta_l = 0), this is trivially satisfied.

### 5.3 Left current singularity analysis

At r = 0, the hedgehog has q = -rho_0 (f(0) = pi). The left current:

    l_mu = q^{-1} d_mu q = -(1/rho_0) d_mu q

Near r = 0, with f(r) ~ pi - a*r:

    q ~ rho_0(-1 + a*r * r_hat . sigma)

    d_mu q ~ rho_0 * a * (d_mu r * r_hat . sigma + r * d_mu(r_hat) . sigma)

    l_mu = q^{-1} d_mu q ~ -(1/rho_0) * rho_0 * a * (...) = -a * (...)

At r = 0 exactly: d_mu q involves d_mu(r*r_hat) = delta_mu,i * e_i (unit vectors),
so l_mu ~ -a * sigma_mu. This is REGULAR at the origin — the left current
does not diverge.

More explicitly: the hedgehog left current in spherical coordinates is:

    l_r = f' * r_hat . sigma

which at r = 0 gives l_r = f'(0) * r_hat . sigma = -a * r_hat . sigma. This is
regular because the function f(r)/r is bounded (f ~ pi - ar near 0, so the
relevant combination sin(f)/r ~ a is smooth).

The angular components:

    l_theta ~ (sin f / r) * sigma_theta

At r = 0: sin(f) ~ sin(pi - ar) ~ ar, so sin(f)/r ~ a. Regular.

So the left current l_mu is SMOOTH at the origin. There is no singularity
to resolve. The claim that "d != 0 might be needed to resolve a singularity
in l" is unfounded — there is no singularity.

### 5.4 Compatibility with Dirichlet energy minimization

If we posit a "full" Dirichlet energy functional:

    E_full = (1/2) integral |d_mu Psi|^2_{full} d^3x

where |.|_{full} is some norm that DOES see the degenerate sector (e.g., the
component-wise L^2 norm rather than the Clifford norm), then:

    E_full = (1/2) integral (|d_mu q|^2 + |d_mu d|^2) d^3x

The cross-term vanishes if we use a simple component-wise inner product.
Then E_full >= E_bulk, with equality iff d = const. The sigma-model constraint
then forces d = 0 (since d = const with the constraint requires d = 0 at
infinity, hence d = 0 everywhere). If d is not constant, E_full > E_bulk,
and d = 0 is still the energy minimizer subject to the constraint.

So even with a hypothetical "full norm" energy functional, d = 0 minimizes
the energy.

### 5.5 What geometric principle COULD force d != 0?

Having exhausted all standard geometric and topological arguments, let me
identify what WOULD be needed to force d != 0:

**A. Lagrangian coupling with the right structure.**
A term in the Lagrangian that couples d to q in a way that d = 0 is NOT a
critical point. For example:

    L_coupling = g * B^0(q) * tau

where B^0 is the baryon density (cubic in derivatives of q). With this term,
the EL equation for tau becomes:

    mu^2 * tau = g * B^0

which forces tau = g * B^0 / mu^2 (nonzero wherever B^0 != 0). But this
coupling must be ADDED to the Lagrangian — it does not follow from the algebra
or topology alone.

This is exactly the WZW term discussed in Track A (constraint_gravity_hypothesis.md).
The WZW term provides the coupling B^0 * tau, and its coefficient N_c = 3 is
topologically quantized. If the WZW term is included, then d != 0 is
DYNAMICALLY forced (by the equations of motion), not topologically.

**B. A non-quadratic norm.**
If instead of |Psi|^2 = rho_0^2 we impose a constraint that involves both
q and d in a nonlinear way (e.g., some secondary norm equals a fixed value),
then d could be forced nonzero. But the standard sigma-model constraint
|q|^2 = rho_0^2 does not do this.

**C. A gauge-fixing condition.**
If the theory has a gauge symmetry under which d transforms, and we require a
specific gauge (like "axial gauge" j_3 = 0), this could constrain d but not
force it nonzero. Gauge fixing removes redundancy, not degrees of freedom.

**D. Boundary conditions from a larger topology.**
If the spatial manifold is not R^3 but (say) S^3, and if d must satisfy
specific boundary conditions that are incompatible with d = 0 globally, then
d could be forced nonzero. But on R^3 with standard asymptotic conditions
(Psi -> rho_0 at infinity, implying d -> 0), there is no such obstruction.

### Verdict on B5

**No geometric consistency condition forces d != 0.** The Maurer-Cartan
equation, Bianchi identity, regularity conditions, and energy minimization
all allow d = 0. The only mechanism that could force d != 0 is a DYNAMICAL
one (a Lagrangian coupling term like WZW), not a topological or geometric one.

---

## 6. Summary: Is d Topologically/Geometrically Forced to Be Nonzero?

### Definitive answer: NO.

The degenerate sector d of the Cl+(3,0,1) field is NOT topologically or
geometrically forced to be nonzero on a B=1 Skyrmion background. Every
argument examined reaches the same conclusion:

| Investigation | Method | Result |
|---------------|--------|--------|
| B1: d = 0 check | Direct substitution | d = 0 satisfies all constraints and EL equations |
| B2: Fiber bundle | Bundle triviality | Zero section is smooth and globally defined |
| B3: Dual curvature | Holonomy analysis | Curvature constrains transport, not existence |
| B4a: Holonomy | Group theory | DH* contains H* as subgroup; pure-q holonomy valid |
| B4b: Index theory | Operator analysis | Relevant operator absent from Lagrangian |
| B4c: Homotopy | pi_3 computation | No additional invariants from d sector |
| B4d: Energy | Variational analysis | d = 0 minimizes energy (or is degenerate at mu=0) |
| B5: Geometric | Maurer-Cartan + Dirichlet | All conditions trivially satisfied at d = 0 |

### Why the negative result is robust

The fundamental reason d = 0 is allowed is ALGEBRAIC:

1. The subalgebra Cl+(3,0) subset Cl+(3,0,1) is closed under all operations
   (multiplication, inversion, differentiation). A field valued entirely in
   this subalgebra is self-consistent.

2. The degenerate direction e_0 (with e_0^2 = 0) makes the d sector INVISIBLE
   to the standard Lagrangian. The theory literally cannot distinguish d = 0
   from d = anything.

3. The sigma-model constraint is LINEAR in d, and d = 0 always satisfies a
   linear constraint.

4. The target space for d is a VECTOR SPACE (or vector bundle), not a sphere
   or other manifold with nontrivial topology. Vector bundles over contractible
   (or even most non-contractible) bases are trivial.

### What this means for the theory

d = 0 being allowed is not a failure — it is the EXPECTED result for the
current Lagrangian. The standard Skyrme model has been studied for 60+ years
with only the quaternionic (bulk) sector. Embedding it in Cl+(3,0,1) does not
automatically activate the extra degrees of freedom.

---

## 7. What Could Force d != 0?

Since topology and geometry alone cannot force d != 0, the mechanisms must
be DYNAMICAL — they require modifying or extending the Lagrangian.

### 7.1 The WZW term (Track A)

The Wess-Zumino-Witten topological term:

    Gamma_WZW = (N_c / 240 pi^2) integral <L^5>

when expanded to order epsilon^1 in the dual number, produces a coupling:

    L_WZW superset N_c * tau * B^0

This gives tau the equation of motion:

    kappa^2 nabla^2 tau - mu^2 tau = -N_c B^0

(assuming a kinetic term kappa^2|nabla tau|^2/2 and mass term mu^2 tau^2/2
are present from the degenerate Lagrangian L_{2,D}).

The solution is:
- Massless (mu = 0): tau(r) ~ N_c B / (4 pi kappa^2 r)  (1/r gravitational)
- Massive (mu > 0): tau(r) ~ N_c B exp(-mu r/kappa) / (4 pi kappa^2 r)  (Yukawa)

**This is a DYNAMICAL forcing.** The WZW term provides a source for tau
proportional to the baryon density. The coefficient N_c is topologically
quantized (typically N_c = 3 for SU(2) ~ QCD). This is the most promising
mechanism for forcing d != 0.

### 7.2 Pseudoscalar kinetic mixing (from grade-4 Lagrangian)

The pseudoscalar grade of the quadratic Lagrangian:

    <(d_mu Psi)(d_mu Psi~)>_4 = 2(d_mu s * d_mu tau - d_mu f_i * d_mu j_i)

provides a kinetic coupling between q and d. However, as shown in the constraint
gravity analysis (Section Step 1 Result), this term:

- Sources tau with nabla^2 s, which has ZERO monopole moment
- Produces tau ~ 1/r^3 at best (dipole), not 1/r
- May be a total derivative (boundary term)

So this mechanism IS a forcing of d != 0 in principle (tau is sourced by
nabla^2 s != 0 in the soliton core), but the resulting d field is weak and
short-ranged.

### 7.3 Skyrme coupling through full algebra norm

If the Skyrme term is computed using a norm that sees all 8 components:

    L_4^{full} = (1/4e^2) sum_{mu<nu} ||[L_mu, L_nu]||^2_{full}

where L = l + epsilon * delta_l, then the commutators would mix bulk and
degenerate sectors. Expanding:

    [L_mu, L_nu] = [l_mu, l_nu] + epsilon * ([delta_l_mu, l_nu] + [l_mu, delta_l_nu])
                   + epsilon^2 * [...] = 0

    ||[L_mu, L_nu]||^2 = ||[l_mu, l_nu]||^2
                         + 2 Re <[l_mu, l_nu], [delta_l_mu, l_nu] + [l_mu, delta_l_nu]>
                         + ...

The cross-term couples delta_l (and hence d) to the bulk curvature. This would
give a nontrivial EL equation for d sourced by the bulk field strength.

However, computing this cross-term in the ||.||^2 requires a norm that sees
the epsilon-component. The standard Clifford norm kills it (e_0^2 = 0). One
would need a modified norm — which means changing the Lagrangian.

### 7.4 Component-wise kinetic term for d

The simplest fix (mentioned in spec/math/03_dynamics.md, Section 5.6) is to
add an explicit degenerate kinetic term:

    L_{2,D} = (kappa^2/2c^2)|d_t d|^2 - (kappa^2/2)|nabla d|^2

This gives d propagation dynamics, but does NOT force d != 0 in a static
soliton background (the EL equation becomes nabla^2 d - (mu^2/kappa^2) d = 0,
which has only the trivial solution d = 0 for mu > 0).

To force d != 0 STATICALLY, we additionally need a source term — which is
what the WZW coupling provides.

### 7.5 Summary of forcing mechanisms

| Mechanism | Forces d != 0? | Long-range? | From algebra? |
|-----------|---------------|-------------|---------------|
| WZW term | YES (source = B^0) | 1/r if mu=0 | Topological (N_c quantized) |
| Grade-4 kinetic | YES (source = nabla^2 s) | No (1/r^3) | Algebraic (coefficient = 2) |
| Full-norm Skyrme | YES (in principle) | Unknown | Requires modified norm |
| Degenerate kinetic alone | NO (source = 0) | N/A | Added by hand |
| Topology/geometry | NO | N/A | N/A |

---

## 8. What Determines the Profile and Normalization?

Since d is NOT topologically forced, its profile must be determined by DYNAMICS.
The most natural scenario:

### With WZW coupling + degenerate kinetic term:

The static equation for tau on the hedgehog background is:

    -kappa^2 nabla^2 tau + mu^2 tau = N_c B^0(r)

where B^0(r) = -f' sin^2(f) / (2 pi^2 r^2).

The sigma-model constraint tau = -j_r tan(f) then determines j_r from tau:

    j_r(r) = -tau(r) / tan(f(r)) = -tau(r) cos(f) / sin(f)

This is well-defined everywhere because:
- At r = 0: f = pi, sin(f) ~ ar, tau ~ tau_0 (finite), so j_r ~ -tau_0/(a*r) * cos(pi)
  Wait: cos(f(0)) = cos(pi) = -1, sin(f(0)) = 0. So j_r = -tau(0)*(-1)/0 ... divergent?

We need tau(0) = 0 for consistency (since the constraint at r = 0 requires
tau = 0, as j_r must be finite). The WZW-sourced tau has tau(0) finite (B^0 is
finite at origin), so we'd need to check self-consistency.

Actually, let's reconsider. At r = 0:
- The hedgehog is q = -rho_0 (scalar, no angular dependence)
- The sigma-model constraint is: s*tau + f.J = 0
- At r = 0: s = -rho_0, f = (0,0,0), so: -rho_0 * tau(0) = 0, giving tau(0) = 0.
- But B^0(0) is finite (B^0 = -f'sin^2f/(2pi^2 r^2) -> -(-a)(ar)^2/(2pi^2 r^2) = a^3/(2pi^2))

So the Poisson equation gives: nabla^2 tau(0) = -N_c B^0(0)/kappa^2 != 0.
With tau(0) = 0 and nabla^2 tau(0) != 0, this means tau''(0) != 0. So tau ~ C*r^2
near the origin. The constraint then gives:

    j_r(0) = -tau(0)/tan(f(0)) = 0/0

which is indeterminate but resolved by L'Hopital: as r -> 0,
tau ~ C r^2, tan(f) ~ tan(pi - ar) = -tan(ar) ~ -ar.
So j_r ~ -C r^2 / (-ar) = C r / a -> 0. Regular!

So the WZW-sourced solution IS compatible with the sigma-model constraint and
regularity at the origin. The profile would be:

    tau(r) ~ C r^2 near r = 0  (with C = -N_c B^0(0)/(6 kappa^2))
    tau(r) ~ N_c B / (4 pi kappa^2 r) for r >> R_soliton  (if mu = 0)

    j_r(r) = -tau(r) cos(f(r)) / sin(f(r))

This has already been computed numerically in degenerate.c (with manually
chosen coupling g_top in place of N_c). The results there
(E_self ~ 0.027 code units for massless, Delta M/M ~ 10^{-4}) would apply
with the appropriate identification of g_top = N_c.

### Normalization

The normalization of d is set by:
1. The WZW coefficient N_c (topologically quantized)
2. The degenerate kinetic scale kappa
3. The degenerate mass mu

These are physical parameters. N_c = 3 is fixed by topology (for SU(2) in
the standard Skyrmion framework). kappa and mu are the remaining free parameters
of the degenerate sector.

---

## 9. Conclusions

### Main result

The topology of the B=1 Skyrmion does NOT force the degenerate sector to be
nonzero. This is a definitive negative result, established through five
independent lines of investigation (sigma-model constraint, fiber bundle
analysis, dual connection curvature, holonomy/index/homotopy, and geometric
consistency). The fundamental reason is that d = 0 lies in a closed
subalgebra Cl+(3,0) of Cl+(3,0,1), and the current Lagrangian cannot
distinguish fields in this subalgebra from the full algebra.

### Implication for gravity

The constraint gravity hypothesis (that d produces a 1/r gravitational
potential) requires a DYNAMICAL mechanism to activate d, not a topological one.
The most promising such mechanism is the WZW topological term, which couples
the pseudoscalar tau to the baryon density B^0 with a quantized coefficient N_c.
This is the subject of Track A.

If Track A confirms the WZW coupling, then d IS forced nonzero — but by the
equations of motion (dynamically), not by topology. The distinction matters:
topological forcing would be parameter-free, while dynamical forcing depends
on the coupling constants (kappa, mu) and the WZW coefficient (N_c).

### What topology DOES constrain

While topology does not force d != 0, it does constrain what d CAN be if it
IS nonzero:

1. **Regularity**: tau must vanish at r = 0 (from the sigma-model constraint
   at q = -rho_0). Any smooth d on the hedgehog must have tau(0) = 0.

2. **Asymptotics**: tau -> 0 as r -> infinity (from f -> 0 and the constraint
   tau = -j_r tan f).

3. **No covariantly constant d**: The holonomy of the bulk SU(2) connection
   forbids d from being covariantly constant (unless d = 0). Any nonzero d
   must vary spatially.

4. **Consistency at f = pi/2**: At the radius where f = pi/2, the radial
   component j_r must vanish (otherwise tau diverges through tan(f)).

These constraints determine the SHAPE of d (given its amplitude), but not
whether d is nonzero in the first place.
