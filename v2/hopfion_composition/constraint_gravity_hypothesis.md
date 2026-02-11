# Constraint Gravity from Cl+(3,0,1): The Degenerate Sector as Gauss's Law

## The Problem

We need long-range (1/r) spin-2 gravity to emerge from the Cl+(3,0,1) Skyrme
framework. Four paths were investigated:

| Path | Result | Limitation |
|------|--------|------------|
| 1. Finite-lambda | P/m = 2 algebraically | Dead end |
| 2. L_6 sextic | Short-range (~0.55 fm) | Nuclear, not gravitational |
| 3. B^0 p coupling | 1/r with correct sign | Spin-0, free parameter g_top |
| 4. Hopfion metric | Rank-2 tensor, anisotropic | Core-scale only |

Path 3 gives 1/r but has two problems:
1. The coupling g_top is a free parameter (not derived from the algebra)
2. It's scalar (spin-0) gravity, not tensor (spin-2)

Path 4 gives a genuine rank-2 tensor but only at core scale.

No path combines long-range 1/r with tensor structure.

## The Core Obstacle (as understood before this insight)

The degenerate sector (j_1, j_2, j_3, p) of Psi in Cl+(3,0,1) has no dynamics
because e_0^2 = 0. The standard Lagrangian:

    L = <(d Psi)(d Psi~)>_0

extracts the scalar part (grade 0) of the product. With e_0^2 = 0, every term
involving the degenerate sector vanishes. So (J, p) are algebraically constrained
to zero -- they can't propagate, can't mediate forces.

This was seen as the fundamental problem: the degenerate sector is "dead."

## The Circuit Analogy (Key Insight)

Consider electromagnetic field propagation in a circuit:

1. Current flows in the wire (localized source)
2. The EM field exists in the vacuum AROUND the wire
3. Energy flows through the surrounding space (Poynting vector S = E x H)
4. The field extends to infinity even though the source is finite

Now consider the electrostatic potential specifically:

- In Coulomb/temporal gauge, A_0 has **no kinetic term** in the Lagrangian
- The term (1/2)(dA_0/dt)^2 is absent -- just like e_0^2 = 0 kills the
  degenerate sector's kinetic term
- Yet A_0 produces the **1/r Coulomb potential**
- How? Through a **constraint equation** (Gauss's law):

    nabla^2 A_0 = -rho / epsilon_0

- This is NOT a wave equation -- it's Poisson's equation
- The solution is instantaneous (non-retarded): A_0 = (1/4pi eps_0) int rho/r dV
- The 1/r field is a CONSTRAINT, not a propagating degree of freedom

**The parallel to our theory is exact:**

| Electromagnetism | Cl+(3,0,1) theory |
|------------------|-------------------|
| Scalar potential A_0 | Pseudoscalar field p |
| No kinetic term (Coulomb gauge) | No kinetic term (e_0^2 = 0) |
| Charge density rho | Baryon density B^0 |
| Gauss's law: nabla^2 A_0 = -rho/eps_0 | Constraint: nabla^2 p = -(coeff) * B^0 |
| Coulomb field: A_0 ~ Q/r | Gravitational field: p ~ B/r |
| Coupling from gauge symmetry (not free) | Coupling from algebra (not free?) |

## The Reinterpretation

e_0^2 = 0 is not a bug. It is the answer.

It means the degenerate sector is a **constraint sector** -- analogous to A_0 in
electromagnetism. Constraint fields don't propagate (no waves), but they DO produce
1/r potentials via Poisson's equation when coupled to a source.

The soliton's topological charge B (baryon number) is the source -- analogous to
electric charge. The pseudoscalar p is the "Coulomb field" of baryon number.

## The Missing Coupling

For this to work, we need a coupling between p and B^0 that survives e_0^2 = 0.
The current Lagrangian uses <...>_0 (scalar grade extraction), which kills all
degenerate-sector terms.

But the full geometric product (d Psi)(d Psi~) has components at ALL grades, not
just grade 0. Consider the cross-terms between bulk and degenerate sectors:

### Algebraic structure

The field decomposes as Psi = q + d, where:
- q = s + f_1 e_{23} + f_2 e_{31} + f_3 e_{12}   (bulk quaternion)
- d = j_1 e_{01} + j_2 e_{02} + j_3 e_{03} + p e_{0123}   (degenerate)

The reversal:
- q~ = s - f_1 e_{23} - f_2 e_{31} - f_3 e_{12}   (quaternion conjugate)
- d~ = -j_1 e_{01} - j_2 e_{02} - j_3 e_{03} + p e_{0123}
  (bivectors negate, pseudoscalar stays)

The product (d_mu Psi)(d_mu Psi~) contains:
1. (d_mu q)(d_mu q~) -- purely bulk, gives |dq|^2 at grade 0 (the standard L_2)
2. (d_mu q)(d_mu d~) + (d_mu d)(d_mu q~) -- cross-terms
3. (d_mu d)(d_mu d~) -- purely degenerate, vanishes (e_0^2 = 0)

**Term 2 is the key.** Consider a representative cross-term:

    (f_1 e_{23}) * (p e_{0123})

Computing the geometric product:
    e_{23} * e_{0123} = e_2 e_3 e_0 e_1 e_2 e_3
                      = e_0 e_1 (e_2)^2 (e_3)^2  (after sorting, using e_i^2 = 1)
                      = e_{01}

This is a bivector involving e_0 LINEARLY (not squared). Its coefficient is NOT
proportional to e_0^2. It SURVIVES the e_0^2 = 0 condition.

More generally, the cross-terms produce grade-2 and grade-4 components that involve
e_0 only linearly. The pseudoscalar part (grade 4):

    <(d_mu Psi)(d_mu Psi~)>_4 = (cross-terms between q and d at grade 4)

This is a topological density -- related to the Chern-Simons / Wess-Zumino-Witten
structure. It's the pseudoscalar part of the Lagrangian that the current formulation
discards by taking only <...>_0.

### What the constraint equation would look like

If the action includes a pseudoscalar term:

    S_top = integral <(d_mu Psi)(d_mu Psi~)>_4 * (e_{0123})^{-1} d^4x

then varying with respect to p gives an Euler-Lagrange equation. Since p appears
WITHOUT derivatives in some cross-terms (it multiplies d_mu q), the resulting
equation would be a constraint:

    nabla^2 p = (algebraic coefficient from Cl(3,0,1)) * B^0

where B^0 = -f' sin^2(f) / (2 pi^2 r^2) is the baryon density.

This is Gauss's law for baryon charge. The solution:

    p(r) = (coefficient) * B / (4 pi r)    for r >> R_soliton

## Why This Could Work (and Previous Paths Missed It)

### What changes vs Path 3

Path 3 (degenerate.c) already computed this scenario with manually chosen g_top.
The new insight is that g_top might not be free -- it could be determined by the
algebra through the pseudoscalar grade of (d Psi)(d Psi~).

| Aspect | Path 3 (manual) | Constraint gravity |
|--------|-----------------|-------------------|
| Coupling g_top | Free parameter | Determined by algebra |
| Lagrangian term | Added by hand: g_top B^0 p | From <(dPsi)(dPsi~)>_4 |
| Physical origin | Postulated | Topological (WZW-like) |
| e_0^2 = 0 role | Problem (kills dynamics) | Feature (makes p a constraint) |

### The density-as-geometry connection

A soliton creates a region of non-trivial field configuration. If propagation is
c * dt = ds, then any modification of the effective propagation speed c_eff(r)
changes the effective geometry of space: slower c = "larger" space, faster c =
"smaller" space.

The BLV effective metric computes this for perturbations of the field itself:
- g^{00}(r) = 1 + c_4 sum|A_i|^2 / |q|^4   (temporal inertia)
- g^{jj}(r) = 1 + c_4 sum_{i!=j}|A_i|^2 / |q|^4   (spatial stiffness)

For the hedgehog, P/m = g^{rr}/g^{00} = 2 everywhere (conformally flat -- no
differential lensing). This is algebraic for L_2 + L_4.

But the BLV metric DOES show that the soliton core is geometrically different from
vacuum (g^{00} ~ 17 at the center). The modification just happens to be conformal
(time and space stretch equally). The soliton IS a region of modified geometry --
the modification just doesn't produce 1/r behavior because it's confined to the
core (exponential falloff).

The 1/r behavior requires a field that satisfies nabla^2 phi = source with no
mass term (Gauss's law / Poisson's equation). The constraint field p, if coupled
to B^0 through the algebraic structure, provides exactly this.

### The full picture

Combining the insights:

1. The soliton creates a region of modified geometry (BLV metric, Path 4 confirmed
   this is anisotropic for non-hedgehog topologies)

2. The topological charge B acts as the source for a constraint field p (like
   electric charge sources the Coulomb field)

3. The constraint field p extends to infinity as 1/r (Poisson equation, no mass
   gap because it's a constraint, not a propagating field)

4. The 1/r field of p modifies the effective geometry for distant observers/solitons

5. The coupling is determined by the algebra (pseudoscalar grade of the kinetic
   term), not a free parameter

## What Needs to Be Computed

### Step 1: Explicit pseudoscalar grade computation

Compute <(d_mu Psi)(d_mu Psi~)>_4 explicitly in Cl(3,0,1). Identify all terms
that couple the degenerate sector (p, J) to the bulk sector (q) and survive
e_0^2 = 0.

Key question: do the cross-terms give a coupling of the form (something) * p * B^0,
or a different structure?

### Step 2: Euler-Lagrange equation for p

If the action includes the pseudoscalar term, derive the equation of motion for p
by variation. Determine:
- Is it a constraint (nabla^2 p = source) or a wave equation?
- What is the source term? Is it B^0, or something else?
- What is the algebraic coefficient? (This determines G_eff)

### Step 3: Comparison with known topological terms

The pseudoscalar part of (d Psi)(d Psi~) may be related to:
- The Wess-Zumino-Witten term (known to couple omega-meson to B^0 with coefficient
  N_c = 3 in QCD)
- The Chern-Simons form
- The Euler density (topological invariant of the field configuration)

If it matches the WZW structure, the coupling coefficient would be topologically
quantized (integer or rational) -- genuinely parameter-free.

### Step 4: Self-consistency check

With the derived p(r) ~ B/r field:
1. Does it modify the BLV effective metric at long range?
2. Does the modification have the correct sign (attractive)?
3. What is the effective G_Newton in terms of Cl(3,0,1) parameters?
4. Is the spin-0 or spin-2 character determined by whether p alone (scalar)
   or (p, J) together (vector/tensor) enter the constraint?

### Step 5: Tensor structure

The degenerate sector has FOUR components: (j_1, j_2, j_3, p). This is a
4-vector in spacetime terms. If all four satisfy constraint equations sourced by
the baryon 4-current B^mu = (B^0, B^i):

    nabla^2 j_i = (coeff) * B^i
    nabla^2 p   = (coeff) * B^0

then the constraint field (j, p) transforms as a 4-vector -- giving spin-1 (like
electromagnetism), not spin-0 or spin-2.

For spin-2, you'd need the constraint field to be a symmetric tensor. This would
require coupling to the energy-momentum tensor T^{mu nu}, not just the baryon
current B^mu. Whether the algebra provides such a coupling is an open question.

However: even spin-1 gravity (like "graviphoton" or "gravitovector") coupled to
a conserved current gives 1/r with correct sign for like charges. The distinction
from GR (spin-2) appears in:
- Light deflection (spin-1 gives half of GR, same as spin-0)
- Gravitational wave polarization (vector vs tensor)
- Coupling universality (spin-1 couples to baryon number, not mass-energy)

## Summary of the Hypothesis

**Claim**: The degenerate sector of Cl+(3,0,1) with e_0^2 = 0 naturally produces
a 1/r gravitational potential through a constraint mechanism analogous to Gauss's
law in electromagnetism.

**Key elements**:
1. e_0^2 = 0 makes p a constraint field (feature, not bug)
2. The pseudoscalar grade <(dPsi)(dPsi~)>_4 couples p to B^0 (survives e_0^2 = 0)
3. The constraint equation nabla^2 p = (coeff) * B^0 gives p ~ B/r
4. The coupling coefficient is algebraically determined (not free)
5. The mechanism is topological (WZW-like)

**What this would mean**: The gravitational hierarchy (G_Newton << G_nuclear) is
explained by the pseudoscalar coupling being parametrically smaller than the
scalar coupling. The gravitational potential is carried by the degenerate sector
that was previously thought to be inert.

**Critical test**: Compute <(d_mu Psi)(d_mu Psi~)>_4 explicitly and verify that
it contains a p * B^0 coupling term with a definite, calculable coefficient.

---

## Step 1 Result: Explicit Pseudoscalar Grade Computation

### The computation

The field decomposes as Psi = q + d where:
- q = s + f_1 e_{23} + f_2 e_{31} + f_3 e_{12}   (bulk, grade 0+2)
- d = j_1 e_{01} + j_2 e_{02} + j_3 e_{03} + tau e_{0123}   (degenerate)

We need the pseudoscalar (grade-4) part of (d_mu q)(d_mu d~) + (d_mu d)(d_mu q~).

Using reversal rules: q~ negates bivectors, d~ negates e_{0i} but keeps e_{0123}.

**Representative products** (all evaluated in Cl(3,0,1)):

1. (d_mu s)(d_mu tau): s is grade-0, tau*e_{0123} is grade-4.
   Product: (d_mu s)(d_mu tau) * e_{0123}   →   grade 4, coefficient = (d_mu s)(d_mu tau)

2. (d_mu f_i e_{ij})(d_mu j_k e_{0k}): The geometric product e_{ij} * e_{0k} gives:
   - e_{23} * e_{01} = e_{0123}   (grade 4)
   - e_{31} * e_{02} = e_{0123}   (grade 4)
   - e_{12} * e_{03} = e_{0123}   (grade 4)
   - Cross terms (i != k): produce grade-2 bivectors, NOT grade-4

3. Including both orderings (d_mu q)(d_mu d~) + (d_mu d)(d_mu q~) and accounting
   for the reversal signs:

### Result

    <(d_mu Psi)(d_mu Psi~)>_4 = 2[d_mu s * d_mu tau
                                  - d_mu f_1 * d_mu j_1
                                  - d_mu f_2 * d_mu j_2
                                  - d_mu f_3 * d_mu j_3] * e_{0123}

This expression:
- **DOES survive e_0^2 = 0** (e_0 appears only linearly in e_{0123})
- **IS a coupling between bulk (s, f_i) and degenerate (tau, j_i) sectors**
- Has the structure of a **kinetic mixing** term (derivatives of both sectors)

### The problem: monopole structure

When we add L_ps = <(d_mu Psi)(d_mu Psi~)>_4 * (e_{0123})^{-1} to the action and
vary with respect to tau, we get the Euler-Lagrange equation:

    nabla^2 tau = nabla^2 s = nabla^2 (rho_0 cos f)

on the hedgehog background (where s = rho_0 cos f(r)).

**This is NOT the baryon density B^0.** The source is nabla^2 s, which by the
divergence theorem has ZERO total monopole moment:

    integral nabla^2 s dV = integral (grad s) . dS = 0    (f -> 0 at infinity)

Compare with the baryon density:

    integral B^0 dV = B = 1    (unit topological charge)

Since the source has zero monopole, the resulting tau field falls off as 1/r^3 or
faster (dipole or higher), NOT as 1/r. **The pseudoscalar grade of the quadratic
Lagrangian does not produce gravitational 1/r behavior.**

### Why B^0 cannot appear in any quadratic Lagrangian

The baryon density is:

    B^0 = -(1/2pi^2) epsilon_{ijk} Tr(L_i L_j L_k)

where L_i = q^{-1} d_i q are the left currents. This is **cubic** in derivatives of q.

Any quadratic Lagrangian L = (d Psi)(d Psi~) contains at most **two** derivatives.
Therefore no grade extraction of the quadratic Lagrangian can produce B^0 as a source
term -- B^0 has the wrong number of derivatives.

### What IS needed: the Wess-Zumino-Witten term

The WZW topological term is fifth-order in the field:

    Gamma_WZW = -(N_c / 240 pi^2) integral_{D^5} Tr(L ^ L ^ L ^ L ^ L)

where the integral is over a 5-dimensional disk whose boundary is spacetime.
In 3+1 dimensions, this reduces to a 4D integral involving the baryon current:

    L_WZW = (N_c / 48 pi^2) epsilon^{mu nu rho sigma} Tr(L_mu L_nu L_rho) * omega_sigma

where omega_sigma is a vector meson field (the "omega meson" in QCD).

**The key identification**: In Cl+(3,0,1), the pseudoscalar tau plays the role of
the temporal component omega_0 of the omega meson. The WZW term then gives:

    L_WZW superset N_c * tau * B^0

This provides exactly the coupling we need:
- Source is B^0 (unit monopole → 1/r)
- Coupling coefficient N_c is topologically quantized (= 3 in QCD)
- The constraint equation becomes: nabla^2 tau = (N_c / kappa^2) * B^0
- Solution: tau ~ N_c * B / (4 pi kappa^2 r)

### Revised status of the hypothesis

The constraint gravity mechanism (e_0^2 = 0 → constraint field → 1/r) is
**structurally correct** but requires the WZW topological term, not the
quadratic Lagrangian. The pseudoscalar grade of (d Psi)(d Psi~) provides
kinetic mixing between sectors (potentially important for propagation), but
the monopole coupling to B^0 comes from the topological (fifth-order) term.

| Aspect | Quadratic L (grade 4) | WZW term |
|--------|----------------------|----------|
| Order in fields | 2nd | 5th |
| Source for tau | nabla^2 s (zero monopole) | B^0 (unit monopole) |
| tau(r) behavior | 1/r^3 (dipole) | 1/r (monopole) |
| Coupling coeff | 2 (from algebra) | N_c (topological) |
| Status | Computed, insufficient | Required for gravity |

### Implications

1. **The e_0^2 = 0 constraint mechanism remains valid.** tau is still a constraint
   field with no propagating degrees of freedom, solving Poisson's equation.

2. **The coupling must come from a topological term.** This is actually good news:
   the WZW coefficient N_c is quantized (integer), giving a parameter-free coupling.

3. **The question becomes**: Does Cl+(3,0,1) naturally produce a WZW-like term?
   The standard WZW term requires SU(2) → SU(2)_L x SU(2)_R chiral symmetry.
   The even subalgebra Cl+(3,0,1) contains quaternions (isomorphic to SU(2)),
   so the mathematical structure exists. The specific form needs derivation.

4. **In standard Skyrme phenomenology**, the WZW term is what gives the omega
   meson its coupling to baryon number (with N_c = 3). This is well-established
   physics. The new element here is identifying tau (the e_{0123} component)
   as the omega_0 field and recognizing that e_0^2 = 0 makes it a constraint
   rather than a propagating field.

---

## Two Parallel Investigations (Track A and Track B)

### The "field as matter" principle

The core differentiator of this theory: particles ARE field, field IS matter. A soliton
isn't a particle moving through a field — it's a topological configuration OF the field.
This changes what "massless" means:

- **Massive** = localized field structure (soliton, topologically trapped, permanent)
- **Radiation** = propagating perturbation (waves, not trapped, disperses)
- **Constraint structure** = field configuration FORCED to exist by soliton topology.
  Not a soliton (no independent charge). Not radiation (no propagation, e_0^2 = 0).
  Not nothing (determined by geometry, extends to infinity).

The constraint structure is a **geometric shadow** of the soliton. The soliton creates
the topology; the constraint structure is what the rest of the field MUST do in response.
Like a shadow isn't a thing, but it's not nothing either.

### Sigma-model constraint result

The constraint |Psi|^2 = rho_0^2 with Psi = q + epsilon*d gives:

    q*d_tilde + d*q_tilde = 0

On the hedgehog (z-axis, by symmetry j_1 = j_2 = 0):

    tau = j_r * tan(f)

This is ONE equation for TWO unknowns (tau, j_r). The sigma-model constraint
relates them but does not uniquely determine either. A SECOND equation is needed.

This second equation — the "missing constraint" — could come from:
1. A Lagrangian term (WZW) → Track A
2. A geometric/topological condition → Track B

### Track A: WZW from Full Cl+(3,0,1) Left Current

**Idea**: The WZW 5-form Gamma = c integral <L^5>, where L = Psi^{-1} d Psi,
evaluated on the FULL Cl+(3,0,1) field (not just q), automatically produces the
B^0 tau coupling through its dual-number expansion.

**Method**: The full left current is:

    L_mu = Psi^{-1} partial_mu Psi = l_mu + epsilon * delta_l_mu

where l_mu = q^{-1} partial_mu q (bulk) and
delta_l_mu = q^{-1} partial_mu d - q^{-1} d * l_mu (degenerate, = covariant derivative).

Expanding L^5 = (l + epsilon*delta_l)^5 to order epsilon^1 (epsilon^2 = 0):

    L^5|_{order epsilon} = sum over all positions of delta_l in l^4

Since l^3 ~ B^0 (the topological density is Tr(l wedge l wedge l)), the terms
l^3 * delta_l * l and permutations produce B^0 * (degenerate sector) coupling.

**Key questions**:
1. Does the order-epsilon term of <L^5> contain tau * B^0 explicitly?
2. What is the coefficient? (Should be topologically quantized)
3. Is the resulting equation for tau a Poisson equation (constraint)?

**Deliverable**: Explicit algebraic computation in Cl(3,0,1) showing whether
<(Psi^{-1} dPsi)^5>|_{epsilon^1} contains the B^0 tau coupling term.

**Output file**: `track_A_wzw_computation.md`

### Track B: Topological Obstruction / Geometric Constraint

**Idea**: The soliton's topology FORCES d to be nonzero. If q wraps S^3 once (B=1),
the degenerate sector d, constrained by tau = j_r tan(f), cannot be continuously
zero everywhere. The topology of q creates a geometric obstruction that determines
d's profile, including its long-range (1/r?) behavior.

**Method**: Three sub-investigations:

**B1. Can d = 0 on the hedgehog?**
tau = j_r tan(f). Setting j_r = 0 gives tau = 0 everywhere. This IS consistent with
the sigma-model constraint alone. But is it consistent with the FULL geometric
structure of Cl+(3,0,1)? Check whether Psi = q + 0 (d = 0) satisfies all
integrability conditions for the left current L = Psi^{-1} dPsi.

**B2. Integrability of the dual connection**
The degenerate part delta_l_mu = D_mu(q^{-1}d) is a "dual connection." Its curvature:

    F^{dual}_{mu nu} = D_mu(delta_l_nu) - D_nu(delta_l_mu)

might be constrained by q's topology. If F^{dual} is forced to be nonzero
(analogous to how a monopole's vector potential must have a Dirac string),
then d cannot be zero.

**B3. Fiber bundle topology**
At each point, q defines a point on S^3. The allowed d forms a fiber (from
tau = j_r tan f). As q wraps S^3 (going from f=pi at origin to f=0 at infinity),
does the fiber twist in a way that forces d to be nonzero? This is a question
about the characteristic class of the bundle.

The concrete test: on the hedgehog, can we find a SMOOTH, everywhere-regular
solution d(r) that is forced to be nonzero by the topology, and if so, what
is its asymptotic behavior (1/r? 1/r^2? exponential?)?

**Key questions**:
1. Is d = 0 geometrically allowed on a B=1 hedgehog? (If yes, Track B fails)
2. If d must be nonzero, what determines its normalization?
3. What is the asymptotic behavior of the forced d(r)?

**Deliverable**: Analysis of whether q's topology forces d != 0, and if so,
the profile of the geometrically determined d(r).

**Output file**: `track_B_topological_obstruction.md`

### Expected outcomes

| Outcome | Track A | Track B | Implication |
|---------|---------|---------|-------------|
| Coupling exists | WZW gives B^0 tau | d forced nonzero | Gravity from algebra |
| Coupling absent | No B^0 tau in L^5 | d = 0 allowed | Need external coupling |
| Mixed | WZW exists but different form | Obstruction exists but short-range | Partial progress |

The most exciting result would be Track B showing d is forced nonzero with 1/r tail,
AND Track A confirming the coupling coefficient. This would mean gravity is a
GEOMETRIC CONSEQUENCE of the soliton's topology in Cl+(3,0,1), with no free parameters.

---

## Results of Parallel Investigation

### Track A: WZW on Full Cl+(3,0,1) — NEGATIVE

Full computation in `track_A_wzw_computation.md` (1652 lines).

**The WZW 5-form on Cl+(3,0,1) does NOT produce a static B^0 tau monopole coupling.**

Three independent obstructions:

1. **Grade obstruction kills <L^5>_0**: The scalar extraction of L^5 at order e_0
   vanishes identically because e_0 is grade 1, and <e_0 * X>_0 = 0 for ANY X.
   The degenerate sector is invisible to the standard WZW scalar trace.

2. **Mixed WZW <L^5>_4 exists but trivial on statics**: The grade-4 extraction
   DOES see the degenerate sector and produces a coupling B^mu w_mu (baryon
   CURRENT times degenerate gauge field). But for a static hedgehog, B^i = 0
   (no baryon flow) and w_0 = 0 (no time dependence), so the coupling vanishes.

3. **Monopole mismatch is robust**: The kinetic cross-term 2(d_mu s)(d_mu tau)
   sources tau through nabla^2 s, which has zero total monopole. No quadratic
   term can produce B^0 (unit monopole) because B^0 is cubic in derivatives.

**What the WZW DOES give**: Coupling for ROTATING Skyrmions (N/P mass splitting),
topological quantization (pi_5(S^3) = Z_2), parity violation in scattering.
But these are dynamical effects — they don't produce static potentials.

**The Lie algebra of Cl+(3,0,1)**: su(2) semidirect (V + R), where V = span{e_{0i}}
and R = span{e_{0123}} is central. The algebra is NOT semisimple — the Killing form
is degenerate on V + R. This is why the standard WZW construction (which requires a
non-degenerate invariant form) cannot work in the usual way.

### Track B: Topological Obstruction — NEGATIVE

Full analysis in `track_B_topological_obstruction.md` (935 lines).

**d = 0 is allowed. No topological obstruction forces the degenerate sector to be nonzero.**

Key findings:
- The subalgebra Cl+(3,0) is closed; Psi = q (with d = 0) is self-consistent
- The fiber bundle of allowed d over S^3 is trivial (pi_2(GL(3,R)) = 0)
- The Maurer-Cartan equation at order epsilon is trivially satisfied at delta_l = 0
- d = 0 minimizes energy (strictly for mu > 0, degenerately for mu = 0)
- No holonomy, index theory, or homotopy argument forces d != 0

### Combined Assessment

Both tracks give negative results. The "row 2" outcome was realized:

| Track | Result | Status |
|-------|--------|--------|
| A: WZW coupling | No static B^0 tau | **CLOSED** |
| B: Topological forcing | d = 0 is allowed | **CLOSED** |

**The constraint gravity hypothesis, as formulated, does not work.**

The e_0^2 = 0 constraint mechanism is structurally valid (constraint fields DO give
1/r from Poisson's equation). But there is no algebraically-determined coupling
between tau and B^0 within the Cl+(3,0,1) framework:
- The quadratic Lagrangian gives kinetic mixing (zero monopole)
- The WZW term gives current coupling (vanishes for statics)
- The topology doesn't force d != 0

### Remaining possibilities

1. **The coupling IS a free parameter** (g_top in degenerate.c). This is the honest
   answer: Cl+(3,0,1) does not determine the gravitational coupling. It must be
   set to match G_Newton (g_top ~ 2.8e-17). The theory explains the MECHANISM
   (constraint field → 1/r) but not the STRENGTH of gravity.

2. **A different algebraic structure** beyond Cl+(3,0,1) might determine the coupling.
   For instance, Cl(3,1) (Minkowski signature, e_0^2 = -1 not 0) has different
   properties. Or a larger algebra with more structure.

3. **Effective metric from nonlinear dynamics** (Path 4 generalization). The BLV
   metric IS anisotropic on non-hedgehog backgrounds. Long-range effects might
   emerge from TIME-DEPENDENT processes (scattering, radiation) rather than
   static potentials. The WZW term contributes to dynamics even though it
   doesn't produce static forces.

4. **The gravity mechanism is fundamentally different** from what was hypothesized.
   Perhaps it involves the energy-momentum tensor T^{mu nu} (which does have
   the right monopole structure: integral T^{00} dV = E_sol), not the baryon
   density B^0. But coupling to T^{mu nu} requires spin-2 (tensor) gravity,
   not the spin-0 scalar gravity of the degenerate.c approach.

---

## Path 6: Variable Signature — Clifford Bundle Approach

### The GR analogy

General relativity required a conceptual shift: the metric is not a fixed background
but a dynamical object determined by matter. Space-time geometry and matter are
coupled through Einstein's equation: G_μν = (8πG/c⁴) T_μν.

In the Cl⁺(3,0,1) framework, the algebra signature plays the role of the metric.
The parameter e₀² = 0 (projective geometric algebra) is analogous to flat space.
The question: can the field configuration itself MODIFY the effective signature?

### The continuous mapping (3,0,1) → (3,1)

The algebras Cl(3,0,1) (e₀² = 0) and Cl(3,1) (e₀² = -1) are related by a
continuous deformation of the metric signature:

    e₀² = -ε²    where ε ∈ [0, 1]

At ε = 0: Cl(3,0,1) — projective, degenerate sector has no dynamics
At ε = 1: Cl(3,1) — Lorentzian, degenerate sector has full dynamics

**Key insight**: The BLV effective metric already computes HOW MUCH the local
field configuration modifies the effective propagation geometry. For the hedgehog:

    g⁰⁰_eff(r) = 1 + c₄ × S₀(r)    where S₀ = f'² + 2sin²f/r²

and c₄ = 2ρ₀²/e² is the Skyrme coupling.

At the soliton core (r → 0):
    g⁰⁰_eff(0) = 1 + c₄ × 3a² = 1 + 6a²ρ₀²/e²  ≈ 13.1 at e=1, ρ₀=1

In vacuum (r → ∞):
    g⁰⁰_eff(∞) = 1    (flat)

The soliton IS a region of modified geometry — g⁰⁰ varies from 13 to 1.

### Self-determined signature

The variable-signature hypothesis identifies:

    ε²(x) ≡ g⁰⁰_eff(x) - 1 = c₄ S₀(x) = (2ρ₀²/e²)[f'² + 2sin²f/r²]

This is NOT a free parameter but is DETERMINED by the field itself. The soliton
creates the geometry that determines the algebra that governs the soliton. This is
the analog of Einstein's self-consistency loop:

    field → metric → algebra → field equation → field

Compare with GR:

    matter → curvature → geodesics → matter motion → matter

### Modified Lagrangian at ε² ≠ 0

When e₀² = -ε²(x), the full Cl⁺(3,0,1) Lagrangian decomposes as:

    L = ⟨(∂Ψ)(∂Ψ̃)⟩₀ = |∂q|² + ε²(x)[|∂τ|² - |∂j|²] + cross-terms

The degenerate sector (j, τ) acquires kinetic/potential energy proportional to
the LOCAL value of ε²(x). Inside the soliton (ε² ≈ 12), the degenerate sector
is nearly as dynamical as the bulk. In vacuum (ε² → 0), it reverts to being
a constraint sector.

### The monopole moment: coupling to ENERGY, not baryon number

The total "degenerate metric charge" (monopole moment of ε²):

    M_ε = ∫ ε²(r) 4πr² dr = c₄ ∫ S₀(r) 4πr² dr

Since E₂ = (ρ₀²/2) ∫ S₀ dV, we have ∫S₀ dV = 2E₂/ρ₀², giving:

    M_ε = c₄ × 2E₂/ρ₀² = (2ρ₀²/e²)(2E₂/ρ₀²) = 4E₂/e²

Using the virial theorem E₂ = E_sol/2 (sigma model):

    M_ε = 2E_sol/e²

**This is remarkable.** The monopole moment of the modified signature equals the
soliton's TOTAL ENERGY (up to a factor 1/e²). Compare:

| Source | Monopole | Couples to |
|--------|----------|------------|
| B⁰ (baryon density) | B = 1 (integer) | Baryon number |
| ε²(r) (metric deformation) | 2E_sol/e² (∝ energy) | Mass-energy |

In GR, gravity couples to the stress-energy tensor T^{μν}, not to any specific
charge. The variable-signature framework naturally produces coupling to ENERGY,
which is exactly what gravity requires for universality.

### The degenerate elliptic equation

If the WZW coupling B^μ w_μ (from Track A) is modulated by the local ε²:

    L_coupling = ε²(x) × N_c × τ × B⁰

then varying τ with the position-dependent kinetic coefficient gives:

    -∇·(ε²(x) ∇τ) = ε²(x) × N_c × B⁰

which simplifies to:

    -∇²τ = N_c B⁰    (if ε² cancels)

But ε² does NOT cancel in general because it appears in the kinetic coefficient
(which involves ∇ε² through the divergence). The full equation is:

    -ε² ∇²τ - (∇ε²)·(∇τ) = ε² N_c B⁰

### Critical issue: degeneracy at infinity

Since ε²(r) → 0 exponentially at large r (the hedgehog profile decays as
f ~ Ce^{-r}/r), the elliptic operator -∇·(ε²∇) becomes DEGENERATE at infinity.

Solving the 1D radial equation:
    dτ/dr = -Q_enc(r) / [r² ε²(r)]

where Q_enc(r) = ∫₀ʳ [source] 4πr'² dr', gives:

At large r where Q_enc → N_c and ε² ~ Ce^{-2r}/r²:
    dτ/dr ~ N_c e^{2r} / (C r²) → ∞

**The solution diverges exponentially.** This is physically correct: if there is
no "restoring force" for τ at infinity (ε² = 0 means no kinetic energy), the
source pushes τ without limit.

### Physical interpretation: signature bubbles

This divergence tells us something profound: **a soliton cannot exist in pure
Cl(3,0,1) and simultaneously generate a well-behaved far field.**

Two resolutions:

**Resolution A — Global minimum signature ε₀²**:
e₀² is not exactly zero but has a small universal floor:

    e₀² = -ε₀²    everywhere (global constant)

The local effective value is ε²_eff(r) = ε₀² + c₄ S₀(r). At large r, ε² → ε₀²
(not zero), and the degenerate sector has weak but nonzero dynamics everywhere.
The constraint equation becomes:

    -∇·(ε²_eff ∇τ) = [source]

At large r: -ε₀² ∇²τ = 0, giving τ ~ const/r (1/r behavior with coefficient
determined by the source monopole and ε₀²). The effective Newton's constant:

    G_eff ~ g_source/(4π ε₀² E_sol²)

The hierarchy G_Newton ~ 10⁻³⁸ reflects ε₀² ~ 10⁻³⁸.

This is essentially Path 3 (degenerate.c) with κ_D² ↔ ε₀², but now ε₀² has a
geometric meaning: it is the global minimum signature of the Clifford bundle.

**Resolution B — Signature bubble**:
The soliton creates a "bubble" of Cl(3,1)-like algebra in a Cl(3,0,1) background.
Inside the bubble, the degenerate sector is dynamical. Outside, it freezes.
The boundary is where ε²(r) transitions from large to negligible.

In this picture, "gravity" is not a 1/r field but a DOMAIN WALL effect — the
metric modification extends exactly as far as the soliton's field profile.
This gives the core-scale (~1.5 fm) effect already computed in Path 4, not
long-range gravity.

### Numerical computation

The following quantities are computed in `proposal/hopfion_search/src/varsig.c`:

1. ε²(r) = (2ρ₀²/e²)[f'² + 2sin²f/r²] from the hedgehog profile
2. M_ε = ∫ε²(r)4πr²dr = 2E₂/e² = E_sol/e² (verified numerically)
3. Solution of -∇·(ε²_eff∇τ) = B⁰ for ε²_eff = max(ε², ε₀²) at various ε₀²
4. Effective G from far-field τ behavior
5. Demonstration that ε₀² → 0 causes τ divergence

### Assessment

The variable-signature framework is conceptually beautiful — it gives gravity
coupling to ENERGY (not baryon number), matches the GR self-consistency loop,
and explains the hierarchy through ε₀² ~ 10⁻³⁸. However:

1. **It does not eliminate the free parameter.** The global floor ε₀² replaces
   g_top from Path 3. Same physics, different language.

2. **No mechanism determines ε₀².** Within the Clifford algebra framework,
   there is no principle that sets e₀² = -10⁻³⁸ rather than 0 or -1.

3. **Resolution B (domain wall) gives short-range, not 1/r.** This reduces
   to Path 4.

4. **The Z₂ symmetry d → -d persists.** Even with variable signature, the
   Lagrangian is even in d (no linear source), so d = 0 remains consistent.
   A source requires the WZW term (odd in d), which vanishes on statics.

The honest conclusion: the mechanism for 1/r gravity (Poisson equation sourced
by baryon/energy density) is structurally sound. The coupling constant g_top
(or equivalently ε₀²) remains a free parameter that must be fit to G_Newton.
The Cl⁺(3,0,1) framework explains WHAT gravity is (constraint field) and WHY
it's weak (ε₀² ≈ 0), but not HOW WEAK (the precise value of G_Newton).

### CORRECTION: Fundamental error in self-referential mapping

The identification ε²(x) ≡ g⁰⁰_eff(x) - 1 is **wrong** for a self-referential
mapping of so(3,1) into Cl(3,0,1). Two distinct mathematical objects were conflated:

1. **BLV effective metric g⁰⁰_eff**: A derived kinematic quantity describing how
   small perturbations propagate on the soliton background. It comes from the
   second variation δ²L/δ(∂_μφ)δ(∂_νφ) and measures the local "sound speed."
   This is a property of the PERTURBATION THEORY, not the algebra.

2. **Algebraic signature e₀²**: A structural parameter of the Clifford algebra
   that determines the commutation relations [Kᵢ, Kⱼ] = -e₀² εᵢⱼₖ Jₖ (boosts).
   This is the Inönü-Wigner contraction parameter: e₀² = 0 gives the Euclidean
   group ISO(3), e₀² = -1 gives the Lorentz group SO(3,1).

For a genuine self-referential mapping, ε² should be determined by the field's
own algebraic/norm structure — not by a derived perturbative quantity. On the
sigma model where |q(x)| = ρ₀ **everywhere** (spatially uniform), the field
carries no position-dependent information that could define a position-dependent
signature. The correct self-referential mapping gives:

    ε₀² = const    (not ε²(r))

### Re-computation with correct mapping (varsig.c rewritten)

With constant ε₀², the norm constraint |q|² + ε₀²p² = ρ₀² couples the pseudoscalar
p to the Skyrmion through the E₂ energy functional. Crucially, E₄ (Skyrme term) is
ρ-independent: for Ψ = ρq̂, [L_i,L_j] = [l̂_i,l̂_j] because ∂ρ/ρ commutes with
all Lie algebra elements. So the coupling potential is W(r) = S₀(r) = f'² + 2sin²f/r²
from E₂ only.

This gives an **eigenvalue problem** (not a Poisson equation):

    -∇²p - S₀(r)p = E₀p

**Numerical results** (Numerov method, bisection on node count, ODE residual 2×10⁻⁹):
- Ground state: **E₀ = -1.4873** (one bound state)
- No excited states (n=1 has no bound eigenvalue)
- Decay: κ = √|E₀| = 1.22, range 1/κ = 0.82 code = 0.46 fm
- Instability: at any ε₀² > 0, the Skyrmion generates a localized p-halo

**This does NOT reduce to Path 3.** The two are fundamentally different:
- Path 3: -κ²∇²p = g_top B⁰ (sourced equation → 1/r far field)
- Path 6: -∇²p = S₀p + E₀p (eigenvalue problem → exponential decay)

The constraint coupling produces a localized pseudoscalar halo, analogous to
the pion cloud. This is a nuclear-scale effect (~0.46 fm), not long-range gravity.
For 1/r gravity, the linear source term g_top B⁰ is still required as a separate
coupling.

**Status: CLOSED — constraint coupling gives localized halo, not 1/r gravity.**
