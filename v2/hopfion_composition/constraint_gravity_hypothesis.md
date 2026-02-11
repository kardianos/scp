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
