# Post-V22: Nine Proposals for Emergent Gravity

## The Problem

V2 through V22 have systematically closed every conventional path to gravity
from the oscillon/Skyrmion framework. Three fundamental obstructions remain:

1. **Mass gap → Yukawa**: Every massive mediator gives e^{-mr}/r, not 1/r
2. **Depletion → negative feedback**: Energy concentrations create fast-light
   (diverging lens), not slow-light (converging lens)
3. **No spin-2**: Three real scalar fields have no symmetric traceless tensor
   excitation matching graviton polarization

The question is no longer "how to add gravity" but "what IS the field such
that gravity is already there?"

---

## Branching Proposals

```
THE FIELD IS...
│
├─── A: A MEDIUM (vacuum has structure)
│    ├── 1. Superfluid condensate (Goldstone gravity)
│    ├── 2. Elastic continuum (shear-graviton)
│    └── 3. Pre-geometric fluid (density = determinant)
│
├─── B: GEOMETRY ITSELF
│    ├── 4. Displacement coordinates (strain = curvature)
│    ├── 5. Conformal clock field (breathing = ticking)
│    └── 6. Correlation fabric (entanglement metric)
│
└─── C: EMERGENT / STATISTICAL
     ├── 7. Near-critical order parameter (critical gravity)
     ├── 8. Phonon-graviton from oscillon lattice
     └── 9. Dual massless description (gauged conservation)
```

---

## Branch A: The Field is a Medium

### 1. Superfluid Condensate — Goldstone Gravity

The current model has phi=0 as the true vacuum and oscillons as excited lumps.
Invert this: the TRUE vacuum is a uniform condensate phi_a = v (all three
fields nonzero), and the oscillon is a *depletion* — a hole in the condensate
where P = phi_1 phi_2 phi_3 drops to zero.

The potential V = (mu/2)P^2/(1+kappa P^2) with mu<0 has a local minimum at
P=0 AND (if we relax the stability condition m^2 > |mu|(2/kappa)^{2/3}/9) a
LOWER minimum at P = v^3 != 0. In this condensed phase:

- The condensate spontaneously breaks whatever continuous symmetry exists
  → **Goldstone mode = massless**
- The oscillon is a depletion (phi→0 at core) = dark soliton in the superfluid
- The Goldstone mode mediates a **1/r** interaction between depletions
- Coupling strength G ~ 1/(rho_condensate * v^4) — could be very small for
  dense condensate

**Addresses**: Obstruction 1 (mass gap bypassed by Goldstone), Obstruction 3
(coupling could be weak).

**Key obstacle**: The triple product V(P^2) has only discrete symmetry (P→-P),
not continuous. No Goldstone theorem without continuous symmetry. Would need to
embed the three fields into a larger symmetry group (e.g., treat phi_a as an
SO(3) vector, replace P=phi_1 phi_2 phi_3 with |phi|^3).

**Test**: Modify the potential to have SO(3) invariance [V = V(|phi|^2)], find
the condensate vacuum, verify a gapless Goldstone mode, compute the coupling
to depletions. Check if oscillons still form.

---

### 2. Elastic Continuum — Shear Modes as Gravitons

Reinterpret the three fields phi_1, phi_2, phi_3 as the **displacement vector**
of an elastic medium. The field index a IS the spatial direction: phi_a(x) =
displacement of the medium at point x in direction a.

The strain tensor is eps_{ij} = (d_i phi_j + d_j phi_i)/2. Decompose:

- **Compression** (trace): theta = d_i phi_i — this is the massive mode (the
  mass gap m gives compression waves a gap, like optical phonons)
- **Shear** (traceless symmetric): sigma_{ij} = eps_{ij} - theta delta_{ij}/3
  — in a **fluid** (zero shear modulus), these are **massless**

The critical insight: shear modes are a **symmetric traceless tensor** → they
have **spin-2** character, exactly like gravitons. If the vacuum is a fluid (not
a solid), shear modes propagate freely with omega = c_s k (gapless, linear
dispersion).

An oscillon = breathing compression wave trapped below the compression gap.
It generates shear through nonlinear coupling (compression-shear interaction
proportional to eps^2). The shear field falls as **1/r** → this IS gravity.

The gravitational constant:

    G_eff ~ (nonlinear coupling)^2 / (rho_medium * c_shear^4)

This could be naturally small if the nonlinear coupling is weak (second-order).

**Addresses**: ALL THREE obstructions.
- Obstruction 1: shear modes are massless (no gap)
- Obstruction 2: positive feedback — compression concentrations generate shear
  wells that attract more compression
- Obstruction 3: shear is naturally spin-2

**Key obstacle**: The current Lagrangian L = Sum[1/2(d_t phi_a)^2 -
1/2|grad phi_a|^2 - 1/2 m^2 phi_a^2] - V(P) doesn't have elastic structure.
The mass term m^2 phi_a^2 gives mass to ALL modes (compression and shear). For
shear to be massless, the potential must depend only on the STRAIN (gradients of
phi), not on phi itself. The m^2 phi^2 term would need to be reinterpreted as
arising from boundary conditions or a condensate.

**Test**: Add gradient-coupling terms (d_i phi_j * d_j phi_i) to the Lagrangian
while keeping the triple-product coupling. Compute the linearized spectrum
around an oscillon. Separate into scalar (compression), vector, and tensor
(shear) sectors. Check if the tensor sector is gapless.

---

### 3. Pre-Geometric Fluid — Density is the Determinant

The triple product P = phi_1 phi_2 phi_3 has a geometric meaning: it is
proportional to the **determinant** of a diagonal 3x3 matrix
diag(phi_1, phi_2, phi_3). In differential geometry, det(g) determines
the volume element.

Proposal: P IS sqrt(det(g)) — the field encodes the local volume of space.
The vacuum P = 0 represents "zero volume" (degenerate geometry). An oscillon
with oscillating P represents space "breathing" — alternating between finite
and zero volume.

In this view:
- The Lagrangian density L * sqrt(det(g)) naturally includes a factor of P
- The potential V(P^2) is a cosmological-constant-like term weighted by the
  volume element
- Perturbations around a nonzero P_0 would give an emergent metric with
  curvature

The field equation for P near a static background:

    d^2 P/dt^2 = c^2 nabla^2 P - V'(P)

has a kink solution P(r) interpolating between P=0 and P=P_max. The kink IS
the boundary between "space" and "no-space" — a bubble of geometry.

**Addresses**: Reframes the problem — gravity doesn't need to be added because
geometry IS the field.

**Key obstacle**: P = phi_1 phi_2 phi_3 is the product of three independent
fields, not the determinant of a metric. Making this identification rigorous
requires the fields to transform as a frame (dreibein), not as independent
scalars.

**Test**: Compute the effective metric g_eff_{mu nu} around a static oscillon
using P as the conformal factor. Check if perturbations satisfy a wave equation
with Coulombic potential.

---

## Branch B: The Field is Geometry Itself

### 4. Displacement Coordinates — Strain is Curvature

Take proposal 2 further: the three fields ARE the spatial coordinates
themselves. The vacuum is phi_a(x) = 0 (undistorted), and the physical
position of a point x is X_a = x_a + phi_a(x).

The PHYSICAL metric in the displaced frame is:

    g_{ij} = delta_{ij} + 2 eps_{ij} + eps_{ik} eps_{kj}

where eps_{ij} = (d_i phi_j + d_j phi_i)/2 is the strain. The Riemann
curvature R_{ijkl} is quadratic in the strain gradients.

For a spherically symmetric oscillon:
- phi_a = f(r) x_a/r  (radial displacement)
- eps_{rr} = f'(r),  eps_{theta theta} = f(r)/r
- The scalar curvature R ~ d^2 f/dr^2 + (df/dr)/r

An oscillon with f(r) ~ e^{-mr} gives R ~ e^{-mr}/r^2 — Yukawa again.

BUT: what if the curvature sources a MASSLESS perturbation? In GR, the Einstein
equations are nabla^2 h = -16 pi G T, where h is the metric perturbation
(massless!) and T is the stress-energy. If the strain field phi_a has a massive
sector (compression) and a massless sector (conformal deformation), the massive
oscillon could source the massless geometric mode.

The key equation:

    nabla^2 h_{ij}^{TT} = -kappa * [d_i phi_a d_j phi_a]^{TT}

The source on the RHS is a function of the massive field, but the response
h_{ij}^{TT} is massless → falls as 1/r.

**Addresses**: Obstruction 1 (massless h from massive source), Obstruction 3
(TT tensor = spin-2).

**Key obstacle**: Why should the metric perturbation h be massless? In ordinary
elasticity, both compression and shear can be massive (in a crystal).
Masslessness requires a SYMMETRY (diffeomorphism invariance). But the
displacement field phi_a explicitly breaks diffeomorphism invariance by
providing a preferred reference frame.

**Test**: Write the effective metric from the oscillon's displacement profile.
Compute the linearized Einstein tensor. Check if there's a self-consistent
solution where the curvature sources more curvature at 1/r.

---

### 5. Conformal Breathing Clock — Frequency = Time

The oscillon breathes at frequency omega ~ 0.95m. This is an internal clock.
Define LOCAL TIME by the number of breathing cycles:

    T_local(x) = omega(x) * t_coordinate

Near another oscillon, the field environment changes → omega shifts → time
flows at a different rate → gravitational redshift.

In GR, the gravitational potential relates to the clock rate:

    g_{00} = -(1 + 2 Phi/c^2),  so  Phi = c^2 (omega/omega_0 - 1)

If the breathing frequency shift is delta_omega/omega_0 = -Phi/c^2, then
we need Phi ~ 1/r.

How could omega depend on distant matter at 1/r? The breathing frequency:

    omega^2 = m^2 - Delta(mu, kappa, A^2)

where Delta is the nonlinear frequency shift. If m^2 has a 1/r tail from a
massless mediator (like Phi in V22), then:

    delta_omega/omega ~ delta_m^2/(2m^2) = beta Phi/(2m^2) ~ beta E/(8 pi m^2 r)

This gives the right form.

**Radical twist**: Don't ADD a mediator. The breathing frequency omega is ITSELF
the mediator. If local omega depends on surrounding energy through energy
conservation d_t rho + div j = 0, and the time-averaged energy density is
constrained at all distances by the conservation law, then the frequency
modulation propagates through the conserved current — which is EXACT and does
not decay exponentially.

**Addresses**: Philosophically compelling — gravity = time dilation = clock rate
modulation.

**Key obstacle**: No mechanism for 1/r clock shift from massive fields. The
conserved energy current constrains the integral of rho, not the local value.
Still need a massless channel to carry the influence.

**Test**: Measure the breathing frequency of one oscillon as a function of
distance from a second oscillon. Plot delta_omega vs separation. Check if it
falls as 1/r or e^{-mr}/r.

---

### 6. Correlation Fabric — Geometry from Field Correlations

In quantum gravity (Ryu-Takayanagi), spatial geometry EMERGES from entanglement
structure. Regions that are highly correlated are "close" in the emergent
geometry; uncorrelated regions are "far apart."

Apply this classically: two points x and y have a field correlation:

    C(x,y) = <phi_a(x,t) phi_a(y,t)>_time

In vacuum: C ~ e^{-m|x-y|}/|x-y| (massive correlator, exponential decay).

Near an oscillon: the breathing creates oscillating correlations. The
time-averaged correlation is:

    <C> = <A^2> f(x) f(y)

This is SEPARABLE (product form), defining a different geometry than vacuum.

The metric of the "correlation space":

    ds^2_corr ~ -d(log C)/d(|x-y|) * d|x-y|^2

In vacuum: ds^2_corr ~ m * d|x-y|^2 (flat with scale m).
Near an oscillon: ds^2_corr is modified by f(r), creating curvature.

**Addresses**: Reframes gravity as correlation structure, not a field.

**Key obstacle**: The correlation function is an observable, not a dynamical
field. It doesn't exert forces. Converting this into predictions requires a
rule for how "correlation geometry" affects motion — which is essentially
postulating GR.

**Test**: Compute the full two-point correlation function of the three fields
around an oscillon (with temporal averaging). Map it to a metric. Compute
geodesics in this metric and check if they resemble gravitational orbits.

---

## Branch C: Emergent / Statistical

### 7. Near-Critical Order Parameter — Critical Gravity

The oscillon exists because omega < m (subgap, margin ~5%). This is NOT a
coincidence — the oscillon CAN ONLY EXIST near the gap edge. In 1D the margin
is 18%; in 3D it shrinks to 5%. In higher dimensions, it would shrink further.

In statistical mechanics, systems near a critical point exhibit:
- Diverging correlation length xi → infinity
- Emergent massless modes (critical fluctuations)
- Universal power-law correlations (1/r^{d-2+eta})
- Scale invariance

The oscillon's 5% gap margin means the system is **near-critical**. The
correlation length of the breathing mode is:

    xi ~ 1/sqrt(m^2 - omega^2) = 1/sqrt(m^2 - 0.95^2 m^2) ~ 3.2/m

This is large but finite. At exact criticality (omega → m), xi → infinity.

**The proposal**: As the oscillon evolves and sheds energy, it approaches
criticality (omega → m from below). At the critical point:

- A massless collective mode emerges
- This mode has power-law (1/r) correlations
- It mediates a long-range interaction between oscillons
- The coupling strength G ~ (m - omega)^alpha vanishes at criticality →
  naturally WEAK

This gives gravity as an **emergent, near-critical phenomenon**. The
gravitational constant is G ~ (gap margin)^alpha, which is naturally small
for oscillons near the edge.

**Addresses**: Obstruction 1 (massless mode at criticality), Obstruction 3
(coupling naturally weak from small gap margin).

**Key obstacle**: At exact criticality, the oscillon DIES (omega=m means it can
radiate freely). The massless mode appears precisely when the oscillon ceases
to exist. The useful regime is the APPROACH to criticality, where xi is large
but finite → still Yukawa with very long range.

**Test**: Systematically vary parameters to push omega/m from 0.95 toward 1.0.
Measure the correlation length of vacuum oscillations around the oscillon.
Check if it diverges as omega→m. Compute the interaction between two oscillons
as a function of their gap margin.

---

### 8. Phonon-Graviton from Oscillon Matter

Stop thinking about one or two oscillons. Think about 10^23 of them forming
a GAS or LIQUID. This oscillon matter has collective excitations — **phonons**
— that are automatically massless (acoustic branch with omega = c_s k).

Individual oscillon: massive (E ~ 85), interacts via Yukawa tails (short-range).
Oscillon gas: has sound waves that are MASSLESS, propagating at c_s.

These acoustic phonons could BE gravity:
- They're massless (1/r propagation)
- In a 3D isotropic gas, transverse sound is a symmetric tensor → spin-2
- They couple universally to energy density → gravitational
- The coupling strength G ~ 1/(rho_gas * c_s^2) → small for dense, stiff gas

**The equation of state** determines everything:

    P(rho) = pressure vs oscillon density
    c_s^2  = dP/drho  (sound speed)
    G_eff  = 1/(rho * c_s^2)  (emergent Newton's constant)

For GW speed = c: need c_s = c. This requires relativistic matter (radiation-
dominated gas of oscillons where P = rho c^2/3).

**Addresses**: All three obstructions.
- Massless: acoustic phonons are automatically gapless
- Spin-2: transverse phonons in 3D isotropic media are naturally spin-2
- Universal coupling: all oscillons couple to the same acoustic field

**Key obstacle**: This describes gravity WITHIN an oscillon gas, not in empty
space. The real universe has gravity in vacuum (between isolated masses). Also
requires a very specific equation of state. And the inter-oscillon potential is
Yukawa → the gas might not form or might crystallize.

**Test**: Initialize a box of N oscillons (N=8-27) at random positions. Let them
equilibrate. Measure the spectrum of collective fluctuations. Check if there's
an acoustic branch (omega proportional to k at small k). Compute the effective
equation of state.

---

### 9. Dual Massless Description — Gauged Conservation

Energy-momentum conservation d_mu T^{mu nu} = 0 is an EXACT symmetry of the
system. By Noether's theorem, it is associated with spacetime translations. In
gauge theory language, making this symmetry LOCAL (gauging it) requires
introducing a GAUGE FIELD that couples to T^{mu nu}.

That gauge field IS the graviton. This is essentially how GR can be derived as
a gauge theory of the Poincare group.

The proposal: don't add gravity — DERIVE it from gauging the energy-momentum
conservation of the three-field system.

1. Start with L = L(phi_a, d_mu phi_a) — the oscillon Lagrangian
2. Compute T^{mu nu} from Noether's theorem
3. Promote global translation symmetry to LOCAL: x^mu → x^mu + xi^mu(x)
4. This requires a connection Gamma^mu_{nu rho} (Christoffel symbols)
5. The connection's dynamics follow from consistency (Bianchi identities)
6. Result: Einstein's equations with G determined by the field content

The gravitational constant emerges as G = hbar c / M_Pl^2, where M_Pl depends
on the number and coupling of field species.

**Twist**: Instead of gauging the full Poincare group, gauge only the INTERNAL
permutation symmetry S_3 of the three fields. This discrete gauging gives a
different (possibly weaker) force. Or: gauge the approximate SO(3) rotating
the three field indices → SU(2) Yang-Mills → asymptotic freedom → coupling
WEAKENS at long distances → could mimic gravity's weakness.

**Addresses**: Automatically gives spin-2, massless, universally coupled gravity.

**Key obstacle**: This is "just adding GR." The gauge procedure doesn't
determine G; it just formalizes the coupling. The gravitational constant remains
a free parameter. And the procedure treats spacetime as fundamental (background
metric exists a priori).

**Test**: Formally construct the stress-energy tensor of the oscillon. Couple it
minimally to linearized metric perturbation h_{mu nu}. Compute the resulting
force between two oscillons at lowest order. Compare G_eff with the oscillon's
internal energy scales.

---

## Ranking

| # | Proposal | Mass gap bypass | Spin-2 | Weak coupling | Testable now |
|---|----------|-----------------|--------|---------------|--------------|
| **2** | **Elastic shear-graviton** | **Yes (shear massless)** | **Yes (sym. tensor)** | **Possibly** | **Modify Lagrangian** |
| **7** | **Critical gravity** | **Near (long xi)** | **No (scalar)** | **Yes (gap margin)** | **Yes (param scan)** |
| **8** | **Phonon-graviton** | **Yes (acoustic)** | **Yes (3D transverse)** | **Possibly** | **Multi-oscillon sim** |
| 1 | Superfluid Goldstone | Yes | No (scalar) | Possibly | Modify potential |
| 5 | Breathing clock | No | No | No | Measure delta_omega(r) |
| 4 | Displacement coords | Needs argument | Yes | Unknown | Compute metric |
| 9 | Gauged conservation | Yes (by construction) | Yes | Free parameter | Formal exercise |
| 3 | Pre-geometric fluid | Unknown | Unknown | Unknown | Conceptual |
| 6 | Correlation fabric | No | N/A | N/A | Compute C(x,y) |

## Assessment

**Proposal 2 (Elastic Shear-Graviton)** is the most promising because it
simultaneously addresses all three obstructions: shear modes are naturally
massless, naturally spin-2, and the nonlinear coupling from compression to
shear can be naturally weak. The required modification is interpretive:
identify field indices with spatial indices, and add gradient-coupling terms
that distinguish compression from shear. The mass term m^2 phi^2 would need
to become a compression-only gap (depending on div(phi), not phi itself).

**Proposal 7 (Critical Gravity)** is the most testable with existing
infrastructure: just scan the parameter space toward the gap edge and measure
correlation lengths. If xi diverges, the critical gravity mechanism is viable.

**Proposal 8 (Phonon-Graviton)** is the most physically concrete but requires
the largest computational investment (many-body oscillon simulations).

The common thread: **gravity requires a massless mode**. The only paths to
masslessness in this framework are symmetry breaking (1, 2), criticality (7),
collective behavior (8), or gauge invariance (9). Pure single-oscillon physics
in a gapped vacuum cannot produce gravity — this has been proven exhaustively
from V2 through V22.

---
---

# Part II: Combinations, Branches, and the Grand Synthesis

## Key Structural Insight: Spin-2 from 3 Fields

Before combining proposals, a critical result from elasticity theory:

The transverse-traceless (TT) strain modes of a 3-component displacement
field are **automatically spin-2**, even though the displacement itself is
spin-1. For a wave propagating in the z-direction:

- Displacement: u_x, u_y (2 transverse DOF, spin-1)
- TT strain: h_+ = eps_xx - eps_yy, h_x = 2 eps_xy (2 modes, spin-2)

The derivative d_i acting on phi_j raises the spin by 1. So **3 scalar
fields → 2 TT spin-2 modes**, matching exactly the graviton's 2 polarizations.
The spin-2 obstruction is solved by the elastic interpretation without adding
fields.

The current Lagrangian gradient term (1/2) Sum_a |grad phi_a|^2 = (1/2)(d_i
phi_a)(d_i phi_a) corresponds to elastic parameters lambda=0, mu=1, mu'=0 in
Lame notation. This already gives transverse (shear) waves at speed c_T = 1
and longitudinal (compression) waves at c_L = sqrt(2). The problem is that
the mass term m^2|phi|^2 gaps BOTH sectors equally.

---

## Combinations That Multiply Power

### 2+4: Elastic Displacement IS Coordinates (Natural Merger)

Proposals 2 and 4 are not really separate — they are two perspectives on one
idea. If phi_a is the displacement of the medium in direction a, then the
strain tensor eps_{ij} = (d_i phi_j + d_j phi_i)/2 IS the metric perturbation
h_{ij}. The key result above (TT strain = spin-2) means this single unified
proposal addresses the spin-2 obstruction directly.

---

### ~~1+2: The Condensate-Elastic Synthesis~~ — RETRACTED

**STATUS: DOES NOT WORK.** The Mexican hat potential was tested in two prior
versions and failed for specific structural reasons:

**V2 failure**: Used V = lambda(|q|^2 - rho_0^2)^2 with sigma-model Skyrmion.
The BLV effective metric satisfies P/m = 2 as an **algebraic identity** — it
holds for ANY rho(r) profile, even with varying density. The Mexican hat
changes rho(r) but cannot break the metric identity.

**V6 failure**: Explicitly removed the Mexican hat because a finite-lambda
potential makes the density mode MASSIVE (m_density = sqrt(2 lambda) rho_0),
giving Yukawa range, not 1/r. V6 replaced it with massless conserved density
+ pressure.

**Additional structural problem**: The SO(3) -> SO(2) breaking from
lambda(|phi|^2 - v^2)^2 gives 2 Goldstone modes that are spin-0/1 (angular
modes on S^2), NOT spin-2. The TT strain argument was optimistic but the
Goldstone bosons are the fundamental propagating DOF, and forces are mediated
by the field itself (scalar/vector), not by its strain.

**Lesson**: Mexican hat potentials gap the radial mode (Yukawa) and give
Goldstone modes with the wrong spin. Neither sector produces 1/r spin-2
gravity. A different mechanism is needed.

See **"Proposal 2 Revisited"** below for an alternative approach.

---

### 2+7: Critical Coupling Enhancement

The elastic medium (2) has a compression-shear coupling that is nonlinear
(second-order in strain). Near criticality (7), this coupling gets **enhanced**
by the diverging correlation length.

The mechanism: the near-critical oscillon (omega ~ 0.95m) has a large "halo"
of weakly bound field extending to xi ~ 3.2/m. This halo couples more
efficiently to shear modes than a compact oscillon would, because the coupling
integral

    integral (oscillon energy density) * (shear Green's function) d^3x

gets contributions from the extended halo out to r ~ xi.

The effective gravitational constant scales as:

    G_eff ~ G_0 * (xi / r_core)^alpha ~ G_0 * (m / (m - omega))^{alpha/2}

where G_0 is the "bare" elastic coupling and alpha is a critical exponent.

This provides a **natural explanation for why G is small but nonzero**: the
oscillon must be near-critical to exist (5% margin), and the criticality
controls the coupling strength.

---

### 2+8: Elastic Medium = Oscillon Lattice (Self-Consistency)

If the universe is filled with oscillons (matter), and the elastic medium IS
the vacuum, then the collective modes of the oscillon system must be
CONSISTENT with the elastic modes of the vacuum.

This is a self-consistency requirement: the phonons of the oscillon lattice (8)
must match the shear modes of the elastic vacuum (2). If they do, the picture
is self-consistent. If they don't, something is wrong.

Test: compute the phonon spectrum of an oscillon lattice and compare with the
elastic shear spectrum. Agreement confirms that gravity is the same at
microscopic (single oscillon) and macroscopic (many oscillons) scales.

---

### 6+2: Correlation as Diagnostic for Elastic Gravity

Proposal 6 (correlation fabric) gains teeth when combined with 2 (elastic).
The field-field correlator:

    C_{ij}(x,y) = <phi_i(x,t) phi_j(y,t)>_time

IS the elastic Green's function, which IS the gravitational propagator.
Computing C_{ij} numerically around an oscillon directly measures the emergent
metric.

The TT part of C_{ij} gives the gravitational (spin-2) propagator. The trace
gives the scalar (compression) propagator. Separating these numerically would
directly test whether the spin-2 channel is gapless while the scalar is gapped.

---

### 7+8: Phase Transition in Oscillon Matter

An oscillon gas near criticality (7+8) has a genuine **phase transition**:

- **Ordered phase** (omega < m): oscillons exist, matter exists, gravity works
  (shear modes propagate)
- **Disordered phase** (omega > m): no oscillons, no matter, no gravity
- **Critical point** (omega = m): scale-invariant, massless modes everywhere

The universe is in the ordered phase, near but not at the critical point. The
distance from criticality determines both the existence of matter AND the
strength of gravity. They are linked.

Prediction: G should depend on the particle mass spectrum. Heavier particles
(further from their gap edge) couple more weakly to gravity. In reality G is
universal — so this prediction is WRONG unless all particles sit at the same
fractional gap margin. This is satisfied if the gap margin is a property of
the MEDIUM, not the individual oscillon. All oscillons in the same medium share
the same margin because the medium's properties are uniform.

---

## Further Branches

### From Proposal 2 (Elastic)

**2a. Cosserat Medium** — The elastic medium has micro-rotations (internal
angular DOF) in addition to displacement. This gives 6 DOF (3 displacement +
3 rotation), allowing spin-2 modes directly without the strain construction.
The micro-rotation field could be identified with the antisymmetric part of
d_i phi_j. Cosserat media have additional "couple stresses" that give richer
wave propagation.

**2b. Gradient Elasticity** — Add higher-gradient terms: L ~ |d^2 phi|^2
(strain gradients). This connects to the L_6 sextic term already studied. In
gradient elasticity, the dispersion becomes omega^2 = c^2 k^2 + beta k^4 —
exactly the Casimir Space paper's dispersion relation. This would unify the
elastic-gravity picture with the Casimir hydrogen spectrum.

**2c. Relativistic Elasticity** — The current model is already relativistic
(Lorentz-invariant kinetic term). But standard elasticity is non-relativistic.
The relativistic extension (Carter-Quintana theory) has the strain defined
through the 4-velocity with the Lagrangian depending on the invariants of the
material metric. This is a well-developed formalism that could be applied
directly to the three-field system.

**2d. Defects as Gauge Fields** — In elastic media, dislocations are described
by a gauge theory (the Volterra process). The dislocation density tensor
alpha_{ij} = eps_{ikl} d_k eps_{lj} acts as the field strength of a gauge
field. If the oscillon is a particular type of defect, its dislocation field
could be identified with the electromagnetic potential. This would give BOTH
gravity (shear) and electromagnetism (dislocations) from the same elastic
medium.

**2e. Topological Elasticity** — In 3D, elastic media can have topological
band structure (topological phononic crystals). Certain phonon modes are
topologically protected at surfaces/edges. Could the massless shear modes be
topologically protected, explaining why they remain exactly massless?

### From Proposal 7 (Critical)

**7a. Self-Organized Criticality** — The oscillon does not need external tuning
to sit near criticality — it NATURALLY evolves toward the gap edge by shedding
energy. This is self-organized criticality (SOC): the system drives itself to
the critical state. The oscillon's energy loss rate dE/dt ~ -10^{-5} is
exactly this SOC process. Eventually every oscillon approaches omega -> m.
The question: is the approach rate related to G?

**7b. Dynamical Critical Exponent** — Near criticality, time and space scale
differently: xi_t ~ xi^z where z is the dynamical critical exponent. For
z != 1, the emergent "speed of gravity" differs from c. For z = 1 (Lorentz
invariance preserved at criticality), the speed of gravity equals c. What
determines z in the oscillon system?

**7c. Multi-Oscillon Criticality** — Single-oscillon criticality gives
xi ~ 3/m. But multiple oscillons could create COLLECTIVE criticality with
xi -> infinity even when individual oscillons are not critical. This is like
percolation: when the oscillon density exceeds a threshold, their halos
overlap and the system undergoes a percolation transition. At the percolation
threshold, correlations become truly long-range.

### From Proposal 8 (Phonon-Graviton)

**8a. Oscillon Superfluid** — If oscillons Bose-condense (all breathing in
phase), the condensate has a macroscopic wavefunction Psi(x). The phonon of
this condensate is the Bogoliubov mode, which is massless (omega = c_s k) and
scalar (spin-0). This gives scalar gravity.

**8b. Oscillon Liquid Crystal** — Oscillons with preferred breathing
orientations could form a nematic phase. The director field fluctuations
(Frank elastic modes) include both scalar and tensor components. The tensor
components are the gravitons.

**8c. Amorphous Oscillon Solid (Glass)** — A disordered arrangement of
oscillons. Unlike a crystal, a glass has no Bragg peaks but still has shear
rigidity. The shear modes of a glass are massless and isotropic — exactly
what is needed for gravity. A glass is more natural than a crystal because
oscillons have no reason to form a periodic lattice.

---

## ~~The Grand Synthesis: 1+2+4+7~~ — RETRACTED

The Mexican-hat-based synthesis is retracted (see corrections to 1+2 above).
The replacement is the **Revised Synthesis: 2+4+7** described next.

---
---

# Part III: Proposal 2 Revisited — Shear Gravity Without Mexican Hat

## The Core Problem

Proposal 2 requires: compression gapped (oscillons exist), shear ungapped
(gravity propagates). The mass term m^2|phi|^2 gaps both equally. The Mexican
hat was supposed to fix this but fails (wrong spin for Goldstone modes, P/m=2
identity, Yukawa radial mode). What else could selectively gap compression
while leaving shear massless?

## Why This Is Structurally Hard

For a gap (omega -> m as k -> 0), the term must be NON-DERIVATIVE: it must
contribute m^2 at zero wavenumber. The only Lorentz-scalar non-derivative
quadratic in phi is |phi|^2, which gaps everything equally.

Derivative terms like (div phi)^2 change wave SPEEDS but not gaps:
- Adding eta(d_i phi_j)(d_j phi_i): c_L^2 = 1 + eta, c_T^2 = 1. Both massless.
- Adding K(div phi)^2: c_L^2 = 1 + K, c_T^2 = 1. Both massless.

So in the VACUUM, with only quadratic terms available, selective gapping is
impossible without breaking Lorentz invariance. All roads lead to m^2|phi|^2.

## The Escape: Selective Gapping Around the Oscillon

The vacuum has a uniform gap m for all modes. But the OSCILLON BACKGROUND
modifies the effective mass through the nonlinear coupling V(P). The
perturbation spectrum around the oscillon has position-dependent, mode-
dependent effective masses.

For the symmetric oscillon phi_1 = phi_2 = phi_3 = f(r)cos(omega t):

The Hessian of V(P) in field space has eigenmodes:
- **Symmetric**: delta phi_a = delta f (all equal) — "compression-like"
- **Antisymmetric**: delta phi_1 = -delta phi_2, etc. — "shear-like"

These eigenmodes see DIFFERENT effective masses:

    m^2_sym  = m^2 + d^2V/df^2|_sym   (compression)
    m^2_anti = m^2 + d^2V/df^2|_anti  (shear-like)

For mu < 0 (attractive triple product): both effective masses are reduced
inside the oscillon core. But the reduction is DIFFERENT for the two sectors
because the Hessian eigenvalues differ.

If the antisymmetric (shear-like) effective mass is reduced MORE than the
symmetric one, then the shear sector approaches m^2_anti -> 0 inside the
core — a near-massless shear mode emerges as a RESONANCE.

## Mechanism: Shear Softening + Cross-Gradient Amplification

Add a cross-gradient coupling eta(d_i phi_j)(d_j phi_i) to the Lagrangian:

    L = Sum_a [1/2(d_t phi_a)^2 - 1/2|grad phi_a|^2]
        - 1/2 eta (d_i phi_j)(d_j phi_i)
        - 1/2 m^2 |phi|^2
        - V(phi_1 phi_2 phi_3)

This term is natural in elasticity (it corresponds to a nonzero Lame
parameter lambda). Its effects:

1. In vacuum: changes compression speed to c_L^2 = 1+eta, shear stays at
   c_T^2 = 1. Both still gapped at m. No change to vacuum physics.

2. Around the oscillon: the cross-gradient couples DIFFERENTLY to the
   compression and shear perturbation sectors. It AMPLIFIES the Hessian
   splitting between symmetric and antisymmetric modes.

3. For some critical eta_c: the antisymmetric (shear-like) effective gap
   CLOSES to zero near the oscillon core. A near-massless shear resonance
   appears.

4. This resonance extends to distances xi >> 1/m (controlled by how close
   m^2_anti is to zero), mediating a quasi-1/r interaction.

5. Combined with proposal 7 (near-criticality, omega ~ m): the shear
   resonance range xi diverges as the gap margin shrinks, giving truly
   long-range gravity at the critical point.

## Why Cross-Gradients Help

The cross-gradient term (d_i phi_j)(d_j phi_i) mixes field index j with
spatial index i. For the perturbation analysis:

- Symmetric perturbation (delta phi proportional to background gradient):
  cross-gradient adds constructively -> stiffens compression
- Antisymmetric perturbation (delta phi orthogonal to background gradient):
  cross-gradient adds destructively -> softens shear

The larger eta is, the more the compression-shear splitting grows. At eta_c,
the shear gap closes.

## The Revised Synthesis: 2+4+7

1. **The field is an elastic displacement** (2): eps_{ij} = d_i phi_j
2. **The strain IS the metric** (4): g_{ij} = delta_{ij} + 2 eps_{ij}
3. **Add cross-gradient coupling** (new): eta(d_i phi_j)(d_j phi_i)
4. **Oscillons are near-critical** (7): omega ~ m, large xi

In this picture:

- In vacuum: all modes have gap m (no change from v22)
- The oscillon background modifies the effective perturbation spectrum
- Cross-gradient coupling amplifies the compression-shear splitting
- The shear-like sector approaches zero gap near the oscillon core
- The near-massless shear resonance has spin-2 character (TT strain)
- Its range xi is controlled by the gap margin (proposal 7 connection)
- G_eff ~ (coupling)^2 / (m^4 xi^2) — could be weak if xi >> r_core

**What changes from v22:**

1. Add cross-gradient term eta(d_i phi_j)(d_j phi_i) to the EOM
2. Keep m^2|phi|^2 (unchanged vacuum)
3. Keep triple-product V(P) (unchanged binding)
4. Oscillons still breathe around phi=0 (unchanged matter emergence)
5. Diagnose the perturbation spectrum for shear-compression splitting

**What to test:**

1. Compute the Hessian of V(P) around the symmetric oscillon background.
   Find eigenvalues for symmetric vs antisymmetric perturbation modes.
   Does the splitting exist even at eta=0? (Analytic + numerical)

2. Add eta cross-gradient to 1D triad code. Scan eta to find where the
   antisymmetric effective gap closes. Do oscillons survive? (1D, fast)

3. If shear gap closes: compute the antisymmetric mode's spatial extent
   (correlation length xi). Does it grow as eta -> eta_c? (1D eigenvalue)

4. Two-oscillon test: measure force vs separation with cross-gradients.
   Is it longer-range than pure Yukawa? (3D, expensive but definitive)

## Alternative Paths for Proposal 2 (If Cross-Gradient Fails)

**Path A: One massless field as mediator**

Break the permutation symmetry: give mass to fields 1,2 but not field 3.

    V_mass = 1/2 m^2 (phi_1^2 + phi_2^2)     [field 3 massless]

Fields 1,2 form the oscillon (massive, trapped). Field 3 is massless and
coupled through the triple product P = phi_1 phi_2 phi_3. The time-averaged
source <phi_1 phi_2> drives phi_3 via Box phi_3 = -dV/dphi_3, giving
phi_3 ~ 1/r at long range.

Pro: simple, 1/r guaranteed. Con: scalar (spin-0), not spin-2. Breaks S_3
symmetry. Need to check if oscillon survives with one massless component.

**Path B: Shear oscillon in incompressible medium**

Enforce div(phi) = 0 (incompressible constraint). Only shear modes exist.
The oscillon must be made of VORTICAL motion (circulating shear waves) rather
than breathing compression. Like a smoke ring — stable, divergence-free.

Pro: only massless (shear) modes exist, automatically 1/r. Con: existing
oscillons are compressive breathers. Need to find a new type of localized
solution. Vortex rings in 3 fields are unexplored.

**Path C: Emergent gaplessness from many-body effects**

Don't try to make shear massless for a SINGLE oscillon. Instead, recognize
that a GAS of oscillons (proposal 8) has collective modes whose effective
gap can be renormalized to zero. The shear gap renormalization:

    m^2_shear,eff = m^2 - n_osc * <delta m^2_anti>

where n_osc is the oscillon density and <delta m^2_anti> is the per-oscillon
shear softening. At critical density n_c = m^2 / <delta m^2_anti>, the
collective shear gap closes → massless mode → gravity.

This is a phase transition: below n_c, gravity is Yukawa; above n_c, gravity
is 1/r. The universe must be above n_c for gravity to work.

Pro: no modification to single-oscillon physics. Con: requires many oscillons,
hard to simulate. The critical density might be unphysically large.
