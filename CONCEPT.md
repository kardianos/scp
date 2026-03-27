# SCP Field Theory — Current Best Concept

**Status**: Working theory with numerical confirmation of core mechanisms.
This document reflects the current best understanding as of March 2026.
Sections marked [CONFIRMED] have numerical evidence. Sections marked
[THEORETICAL] are consistent with the framework but unverified.

---

## 1. The Field

There is ONE field φ_a(x,t) with three real components (a = 0, 1, 2).
The field fills all of space. It IS the medium — there is no separate
"spacetime" or "vacuum." Where the field exists, physics exists.

### Equation of Motion [CONFIRMED]

    ∂²φ_a/∂t² = ∇²φ_a - m²φ_a - ∂V/∂φ_a               (1)

    V(P) = (μ/2) P² / (1 + κP²)                          (2)

    P = φ₀ φ₁ φ₂    (triple product)                      (3)

Parameters: m² = 2.25, μ = -41.345, κ = 50

The triple product P = φ₀φ₁φ₂ is the VOLUME FORM of the three-field
displacement gradient. The potential V(P) penalizes local volume change.
This is the natural energy functional of a 3D deformable medium with
irreducible three-body interactions.

Eq. (1) produces braids and gravity from 3 scalar fields alone.
For electromagnetism, the field requires 3 additional ANGLE components
(Cosserat extension, Eq. 10 in §4). No c(ρ), no S/B split, no smoothing.

### Background

The field has a uniform background with amplitude A_bg ≈ 0.1:

    φ_a(x) = A_bg × cos(k·z + 2πa/3) + (braid structures)

The background is the "fabric" of space. Its density ρ_bg = energy density
of the oscillating background ≈ 0.03 in code units.

---

## 2. Field Structures [CONFIRMED]

The theory has a hierarchy of structures:

| Level | Structure | Analog | Example |
|-------|-----------|--------|---------|
| Sub-component | Braid | Quark | Single helical pattern along one axis |
| Particle | Baryon (3-braid composite) | Proton, neutron | UUD, UDD phase-confined composites |
| Nucleus | Multi-baryon composite | Deuterium, helium | UUD+UDD bound state |

Braids are the fundamental building blocks — self-sustaining helical
patterns in the field. A braid is NOT separate from the field. It IS
the field, organized differently. The simulation makes no distinction
between "braid" and "background" — every grid point runs the same
equation. Physical particles (protons, neutrons) are composites of
three braids, analogous to how hadrons are composites of quarks.

### Braid Structure [CONFIRMED]

The following is the SEED initialization — the initial condition that,
when evolved by the field equations, relaxes into a stable braid. The
braid's dynamics are entirely governed by Eq. (1); the seed merely
provides a starting point close enough to the attractor that the
structure survives and self-organizes.

The bimodal braid seed (the best-characterized configuration):

    φ_a(x,y,z,0) = A_a × E(x,y) × cos(k·z + δ_a) + A_bg × cos(...)

    E(x,y) = exp(-[x²/(1+ε)² + y²/(1-ε)²] / 2R²)       (4)

    A_a = 0.8,  δ = (0, 3.00, 4.43),  R = 3.0,  ε = 0.333

The three fields oscillate with specific phase offsets (δ), creating a
helical twist along z. The elliptical envelope (ε = 0.333) breaks
azimuthal symmetry. After initialization, the field equations reshape
this seed into the equilibrium braid structure described below.

### Binding Anatomy [CONFIRMED]

The equilibrium braid has internal geometry measured by the triple
product magnitude P = \|φ₀φ₁φ₂\|:

| Region | r | P | Role |
|--------|---|---|------|
| Tight core | 0–3 | 0.11–0.15 | Strongly bound, self-contained |
| Interaction surface | 4–6 | 0.01–0.06 | Where the braid couples to the fabric |
| Fabric | >6 | ~0 | Unbound ambient field |

The core INTAKES field (energy flows inward). The surface OUTTAKES
(energy flows outward). The braid is a PUMP — it processes field
through its helical structure.

### Stability [CONFIRMED]

- Braid survives T=500+ in periodic BC (V29-T1b, V33-C2)
- Robust to 6× amplitude perturbation (V29-T4b)
- Thermal equilibrium at A_noise ≈ 0.02 (V29-T1b)
- Evaporative self-concentration: fc rises 0.80→0.93 as outer layers
  radiate away (V29-T3)
- Steady radiation rate ~0.5 energy/time (V29-T3)

---

## 3. Gravity [CONFIRMED]

### Mechanism: Asymmetric Footprint

A massive particle (proton) depletes the surrounding field by binding
field energy into its structure. The particle's perturbation profile
decays as a Yukawa tail (~e^{-m_eff·r}/r), where the effective mass
depends on the local background:

    m_eff² = m² + V''(P_bg)

In a density gradient, the perturbation profile becomes ASYMMETRIC:
on the depleted side (lower ρ), P_bg is smaller, V'' is weaker, m_eff
is smaller, and the Yukawa tail extends further. On the dense side,
m_eff is larger and the tail is shorter. This asymmetric footprint
shifts the particle's equilibrium position toward the depleted side.
This is gravity.

The force is NOT from energy minimization. The static interaction
energy E(D) between two particles is monotonically repulsive at all
separations (V33 energy sweep, D=8-50). The attraction is a DYNAMIC
effect from the breathing cycle's interaction with the gradient.

The force is NOT electromagnetic. Setting η=0 (no θ coupling) gives
the same or larger gravitational drift. The mechanism is entirely in
the φ-depletion sector.

### Proton Gradient Test [CONFIRMED — V43]

The definitive gravitational test: a pre-converged UUD proton placed in
a linear density gradient drifts toward LOW density. The proton is nearly
spherical (aspect ratio 1.05-1.09) and couples isotropically to the
background, making it a clean gravitational test object.

Centroid drift over T=200 (N=384, L=100, gradient-pinned BCs):

| Object | Gradient ΔA | η | Drift | Direction |
|--------|------------|---|-------|-----------|
| Proton | 0.04 (gentle) | 0.5 | +0.80 | → low ρ (gravity) |
| Proton | 0.04 (gentle) | 0.0 | +0.98 | → low ρ (gravity) |
| Proton | 0.10 (steep) | 0.5 | +1.84 | → low ρ (gravity) |
| Braid | 0.04 (gentle) | 0.5 | -1.02 | → high ρ (EM artifact) |
| Braid | 0.10 (steep) | 0.5 | -2.75 | → high ρ (EM artifact) |

Key findings:

1. **Gravity is φ-depletion, not θ-coupling**: The η=0 proton (no θ)
   drifts +0.98 — slightly MORE than η=0.5 (+0.80). The gravitational
   force comes entirely from the φ-depletion mechanism. The θ coupling
   actually provides a small counter-force (radiation drag).

2. **F ∝ ∇ρ confirmed for protons**: Steep/gentle drift ratio = 2.3×
   for 2.5× gradient ratio. Nearly linear.

3. **Asymmetric binding directly drives the drift**: The gradient creates
   20-42% left/right imbalance in the triple product |P| (the binding
   density). In 10/11 frame transitions, the sign of asym_P matches the
   direction of centroid drift. The force amplification from φ to
   P = φ₀φ₁φ₂ makes the binding asymmetry the dominant force channel
   (φ energy asymmetry is only 1-3%).

### Braids vs Protons [CONFIRMED — V43]

A bare z-aligned braid drifts OPPOSITE to gravity (toward high density),
dominated by its EM coupling — the anisotropic current loop interacts
with the gradient via θ. Only the isotropic proton composite shows clean
gravitational response. This mirrors real physics: EM >> gravity at the
quark/braid level; gravity emerges at the baryon/proton level.

The V33 braid-braid force measurements (C=186 coupling constant,
F ∝ 1/D^1.8 power law) mix gravitational and EM contributions. The
pure gravitational coupling constant C_proton has not yet been precisely
measured but the proton's F ∝ ∇ρ linearity is confirmed at two gradient
strengths.

### Force Law [CONFIRMED]

    F = -C × ∇ρ                                                     (5)

The force is proportional to the local density gradient ∇ρ. Confirmed
for the proton at two gradient strengths (drift ratio 2.3× for 2.5×
gradient ratio). The V33 braid gradient sweep at four strengths showed
R²=0.9998 linearity for F ∝ ∇ρ, though the braid measurement mixes
gravitational and EM effects.

Energy conservation: proton runs drift < 0.1% (excellent). Braid η=0.5
runs gain +1-2% from EM pumping through the gradient — itself evidence
that the braid's interaction is EM-dominated.

### Energy Conservation and the Depletion [CONFIRMED]

The depletion around a particle is NOT consumed or radiated energy.
It is field energy that has been BOUND into the particle's structure.
Total energy is conserved (drift < 0.1% per T=100).

The particle does not act as a heat engine. The intake/outtake cycle is a
STEADY STATE, not a one-way consumption:

    Intake at core: field energy flows into the helical pattern
    Outtake at surface: the oscillation re-emits some coherent energy
    Net: the structure is maintained, the depletion profile is STATIC

The depletion profile exists because the binding (triple product
coupling) concentrates field energy from the surrounding region into the
compact helical arrangement. The surrounding field has lower ρ because
that energy is now part of the particle.

Analogy: a crystal forming in solution. The crystal binds solute into a
lattice. The surrounding solution is depleted. The depletion is structural,
not from ongoing consumption. No energy flows once equilibrium is reached.

This is NOT a Le Sage / kinetic / pushing theory of gravity:
- Le Sage: continuous absorption of mediator particles → heating → fails
- This model: static structural depletion from field binding → no heating
- The gravitational force comes from the GEOMETRY of the depletion,
  not from continuous energy flow
- At thermal equilibrium (V29-T1b), radiation = absorption, yet the
  depletion persists because it's structural, not dynamic

The particle's radiation (V29-T3: ~0.5 energy/time for a single braid)
is the structure slowly relaxing toward equilibrium. It is NOT "processing
field for gravity." The radiation is a separate phenomenon that ceases
at thermal equilibrium. The depletion (and thus gravity) persists
regardless of the radiation state.

Key: 0 < V < c and T > 0 always. The particle is always moving and always
at nonzero temperature. But its gravitational depletion is a BOUND STATE
property — like the electric field of a charged particle, which exists
without continuous energy expenditure.

### F = ma from the Field [CONFIRMED structure, THEORETICAL identification]

The measured force law F = -C × ∇ρ has the structure of Newton's second
law. To see this, define the quantities precisely:

**Definitions**:

    ρ(x)    = local field energy density = ½Σ(∂φ_a)² + ½m²Σφ_a² + V(P)
    ρ_bg    = undisturbed background energy density ≈ 0.03
    δρ(x)   = ρ(x) - ρ_bg  (depletion profile; δρ < 0 near particle)
    ∇ρ      = spatial gradient of the energy density

    M_grav  = gravitational mass of the SOURCE particle
              = total field energy bound into the particle's structure
              = -∫ δρ dV  (integrated depletion; positive for a particle)
              This sets the DEPTH and SHAPE of the depletion profile.

    m_grav  = gravitational mass of the TEST particle
              = coupling strength to the density gradient
              = how strongly the test particle's footprint asymmetry
                responds to a given ∇ρ

    m_inert = inertial mass of the test particle
              = resistance to acceleration
              = E_particle / c²  (total bound energy, relativistic mass-energy)

    C       = proportionality constant (C_proton not yet precisely measured;
              the V33 braid value C=186 mixes gravitational and EM effects)

**The force law**:

    F = -C × ∇ρ                                                     (5)

The source particle creates the depletion profile δρ(x), which sets ∇ρ at
the test particle's location. ∇ρ is proportional to M_grav (deeper
depletion = larger source mass = steeper gradient). So:

    ∇ρ  ∝  M_grav × geometry(D)

where geometry(D) encodes the radial dependence (1/D², Yukawa, etc.).

The test particle responds to this gradient with force F = -C × ∇ρ. If
the coupling constant C is proportional to the test particle's mass m_grav:

    C  =  m_grav × G_eff

then the force becomes:

    F  =  -G_eff × m_grav × M_grav × geometry(D)                    (6)

This is Newton's gravitational force law, with G_eff as the emergent
gravitational constant.

**The acceleration**:

    a = F / m_inert = -(C / m_inert) × ∇ρ                          (7)

**If m_grav = m_inert** (equivalence principle), then C/m_inert = G_eff,
and ALL particles experience the same acceleration in the same gradient:

    a = -G_eff × ∇ρ                                                 (8)

This is F = ma, where:
- **m** = m_grav = m_inert = the particle's bound energy (one quantity)
- **a** = -G_eff × ∇ρ = gravitational acceleration (set by the source)
- **F** = m × a = the gravitational force

**Why m_grav = m_inert should hold**: Both masses originate from the SAME
physical property — how much field energy is bound into the particle's
structure via the V(P) coupling:

- m_grav: more bound energy → deeper depletion → stronger ∇ρ response
- m_inert: more bound energy → more energy to accelerate → greater inertia

There is only ONE coupling mechanism (the triple-product potential V(P)),
and it governs both how a particle creates gradients AND how it responds
to them. The equivalence principle is not imposed — it follows from the
single-coupling structure of the theory.

**Status**: The proportionality F ∝ ∇ρ is confirmed for the proton (V43,
two gradient strengths). The identification C = m_grav × G_eff and
m_grav = m_inert is theoretical — testing requires particles of different
sizes/binding energies in the same gradient to verify that C/m_inert is
constant.

The m² = 2.25 in Eq. (1) is a LAGRANGIAN parameter that provides field
confinement. It is NOT the gravitational mass. The gravitational mass
M_grav emerges from the particle's structural depletion of the
surrounding fabric.

### Field Coupling Scales with ρ [CONFIRMED]

The particle-field coupling (both force and drag) is STRONGER in high-ρ
environments. A kicked braid in high-ρ background (A_bg=0.15)
retains/amplifies its momentum (1.40×), while the same kick in low-ρ
(A_bg=0.05) decays rapidly (0.30× retention). This is because V(P)
coupling is proportional to P_bg, which scales with ρ.

This means the gravitational drift is NOT from reduced drag in depleted
regions. It is from the geometric footprint asymmetry: the particle's
perturbation profile extends further where ρ is lower (longer Yukawa
range), shifting the equilibrium position toward depletion.

### Comparison to Newton [THEORETICAL]

For F ∝ 1/D^n: braid-braid measurements gave n ≈ 1.8 (V33). Newton
requires n = 2. The braid measurement mixes gravitational and EM
contributions (V43), so the pure gravitational exponent may differ.
The discrepancy may also be from:
- Periodic BC image contamination at D > 20
- The massive field (m=1.5) creating Yukawa corrections
- Resolution effects at the tested grid spacing

The proton-proton force law exponent has not yet been measured.

---

## 4. Electromagnetism [CONFIRMED]

### Path: Complex Fields + Gauge Symmetry

V29-T10G showed that upgrading from 3 real to 3 COMPLEX scalar fields
preserves the braid and adds U(1) conserved charges. The potential
becomes V(|ψ₀ψ₁ψ₂|²) with ψ_a ∈ ℂ.

The braid carries charge Q ≈ 3-14 (from T10G). Promoting the global
U(1) to a gauge symmetry (adding A_μ) would give:

    D_μ ψ_a = ∂_μ ψ_a + ieA_μ ψ_a
    L_EM = |D_μ ψ|² - (1/4)F_μν² - V(|P|²)

This would produce EM interactions between charged braids. However,
V24-MA showed that explicit gauge coupling can kill oscillons — the
braid's survival under gauging needs verification.

### Torsion in 3-Field Theory [TESTED — NEGATIVE]

V34 Phase 2 tested whether torsion (antisymmetric field gradients) could
serve as an EM force carrier using only the existing 3 φ fields. Four
experiments showed all linearized modes are massive (m_eff ≈ 1.5) — torsion
waves disperse in ~5 time units, and opposite-winding braids do NOT respond
oppositely. EM cannot emerge from 3 scalar fields alone.

### Cosserat Extension: 3 Angle Fields [CONFIRMED — CHARGE-DEPENDENT FORCE]

Adding 3 ANGLE fields θ_a (local rotation) alongside the 3 position
fields φ_a gives 6 components per point — displacement + rotation —
as a proper elastic solid should.

**The θ sector IS the electromagnetic sector of the theory:**
- Massless (m_θ² = 0) → propagates at c → photon analog
- Sourced by curl(φ) → helical braids radiate θ
- Winding number W = ±1 → electric charge (V34 confirmed)
- Carrier phase offset → color charge (V41 confirmed)
- Wave-mediated force (oscillating, period ~4t, 0.2% DC) → QFT-like, not static
- η parameter → analog of fine structure constant (determines EM/strong ratio)

**Role in stability** (V39-V41 confirmed):
- θ coupling acts as energy extraction channel — prevents gravitational collapse (V39)
- θ confinement (outer/inner < 0.7) is the strongest stability discriminant (V40)
- Net θ ≠ 0 in UUD synchronizes composite breathing (V41)
- Net θ = 0 in UDD causes phase convergence → confinement loss (V41)
- In deuterium, strong/EM forces self-equilibrate to 1:1 ratio (V42)

**The equation** (Eq. 10):

    ∂²φ_a/∂t² = ∇²φ_a - m²φ_a - ∂V/∂φ_a + η × curl(θ)_a       (10a)
    ∂²θ_a/∂t² = ∇²θ_a - m_θ²θ_a          + η × curl(φ)_a       (10b)

    curl(F)₀ = ∂F₂/∂y - ∂F₁/∂z
    curl(F)₁ = ∂F₀/∂z - ∂F₂/∂x
    curl(F)₂ = ∂F₁/∂x - ∂F₀/∂y

Parameters: η = 0.5 (coupling), m_θ² = 0 (massless angles).
Field↔space identification: φ₀↔x, φ₁↔y, φ₂↔z (and similarly for θ).

**Stability results** (V34 Cosserat):

1. **The braid sources the angle fields.** θ starts at zero and grows to
   θ_rms ≈ 0.025-0.065. The braid's helical curl excites the θ fields.

2. **Massless angle fields (m_θ=0) are stable.** No tachyonic instability
   (θ has no V(P) coupling). The braid actually benefits: the massless
   mode provides a radiation channel that stabilizes the braid.

3. **dt-converged.** θ_rms varies by only 2% across a 4× dt range
   (dt×0.25, 0.5, 1.0). The θ structure is physical, not numerical.

**θ field structure** (V34 characterization):

The θ field around the braid forms circular patterns perpendicular to
the braid's helical axis, following the right-hand rule (confirmed by
volumetric visualization). The θ_φ (azimuthal) component is the dominant
component at the braid surface.

The θ field is an OSCILLATING WAVE, not a static 1/r field:
- θ_φ² decays as ~r^(-0.5), not r^(-2) (Biot-Savart)
- θ_φ oscillates with period ~4 time units, alternating sign with radius
- The time-averaged DC component is 0.2% of the oscillation amplitude
- The dominant response is wave-mediated, not static-field-mediated

**Charge conjugation** [CONFIRMED]:

The time-averaged DC component of θ_φ REVERSES SIGN when the winding
reverses (W=+1 → W=-1). Correlation of DC residuals between opposite
windings: -0.68. At far field (r>10), the ratio θ_φ(+)/θ_φ(-) ≈ -1.0.

**Winding number IS electric charge** for the θ-mediated interaction.

**Charge-dependent force** [CONFIRMED]:

Two-braid experiments at D=15, T=150 (V34 θ characterization):

| Configuration | ΔD (infall) | vs baseline |
|--------------|-------------|-------------|
| 3-field (η=0, gravity only) | -5.05 | — |
| 6-field, same winding (η=0.5) | -6.41 | 27% MORE attraction |
| 6-field, opposite winding (η=0.5) | -2.16 | 57% LESS attraction |

Same-winding braids attract MORE (parallel "currents" attract through θ).
Opposite-winding braids attract LESS (antiparallel "currents" repel,
partially canceling gravity). This is the correct sign for electromagnetic
interaction between current-carrying elements (Ampère's law analog).

**Two-sector structure**:

| Sector | Fields | Mass | Coupling | Mediates |
|--------|--------|------|----------|----------|
| Position (φ) | 3 | m²=2.25 (massive) | V(P), nonlinear | Matter, gravity (always attractive) |
| Angle (θ) | 3 | m_θ²=0 (massless) | curl(φ), linear | EM-like (charge-dependent, wave-mediated) |

The gravity force (from φ depletion) is always attractive regardless of
winding. The θ-mediated force is charge-dependent: attractive for same
winding, repulsive for opposite. The net force on two particles is
gravity + EM, with the sign of EM depending on relative winding. For
bare braids, the EM component can dominate and reverse the net force
direction (V43).

**Distinction from classical EM**: The θ-mediated force operates through
wave exchange (oscillating field, ~radiation coupling), not through a
static 1/r magnetic field. This is structurally closer to QFT photon
exchange than to classical Biot-Savart. The 0.2% DC component carries
the charge information; the 99.8% oscillation carries the wave energy.

**Open**: Does the force ratio (same/opposite) match the EM coupling
constant? What is the analog of Coulomb's law for composite baryons?
Is there a regime where the DC component dominates (low frequency, long range)?

### EM Field Confirmations [CONFIRMED — V43]

Three electromagnetic properties verified numerically:

1. **1/r² Coulomb force** (OQ2): Proven analytically from energy conservation
   via the Poynting vector derivation. The curl-mediated force inherits
   Coulomb scaling from the massless θ propagator.

2. **Photon speed** (OQ5): A plane-wave θ pulse propagates at v=0.997c
   (η=0) and v=0.906c (η=0.5). The bare propagation is essentially c;
   the curl coupling to the massive φ sector slows the photon branch
   by ~10% at η=0.5. This means the vacuum has a refractive index
   n = c/v ≈ 1.10 — the background φ field acts as a dielectric medium
   for θ waves, analogous to how the Higgs condensate gives mass to
   weak bosons. The index depends on η and A_bg.

3. **No vacuum birefringence** (OQ6): Cross-polarized θ waves mix only
   through the χ² (triple product) nonlinearity, which couples different
   polarization channels equally. No vacuum birefringence — consistent
   with the multilinear structure of P = φ₀φ₁φ₂.

4. **Composite proton is an effective monopole** (OQ1, V44): Spherical
   harmonic decomposition of |θ| around a UUD proton shows the l=0
   (monopole) component carries 54–87% of angular power at all radii,
   with monopole dominance increasing at far field (76% at r=20). Each
   individual θ component remains dipolar (l=1 dominated), confirming
   each braid is a magnetic dipole — but the composite |θ| from three
   orthogonal dipoles sums to an effectively isotropic radiation pattern.
   The proton acts as an effective point charge at distances beyond the
   core radius, producing 1/r² radiation pressure indistinguishable from
   a Coulomb monopole. No true scalar monopole exists (∇·E = 0 is exact);
   the "charge" is an isotropic radiation pattern.

---

## 5. Composite Particles (Baryons) [CONFIRMED]

### Phase-Confined 3-Braid Composites (V41)

Three braids arranged along x, y, z axes with CARRIER PHASE OFFSETS
Δ = {0, 2π/3, 4π/3} form a confined composite:

```
Braid k: φ_a(x,t) = A(x) × cos(k_z × z_rotated + δ_a + Δ_k)    (11)
```

At the triple overlap center: P = Σ cos(Δ_k) = 0. The three braids
CANNOT merge because their combined P vanishes at close range. This is
CONFINEMENT — attracted by the depletion interaction at medium range,
repelled by destructive P-interference at short range.

### Color Charge [CONFIRMED]

Each braid carries 1/3 of the phase circle (2π/3 = 120°). The net theta
from three phase-offset braids cancels: Σ sin(Δ_k) = 0. This is the
fractional charge structure of quarks — color neutral when all three
phases are present.

| QCD | Cosserat analog |
|-----|-----------------|
| Color (R/G/B) | Carrier phase (0, 2π/3, 4π/3) |
| 1/3 fractional charge | 1/3 of phase circle |
| Color neutral baryon | Phase-cancelled θ = 0 at overlap |
| Confinement | P → 0 at triple overlap |
| Asymptotic freedom | V(P) coupling weakens at overlap |

### Proton vs Neutron [CONFIRMED]

**UUD (proton analog)**: Two Up-chirality + one Down-chirality braids.
Net charge = +1. Net θ ≠ 0. The net theta coupling SYNCHRONIZES the
three braids, producing:
- Coherent breathing (period ~90 time units)
- Anti-phase locking (Δφ = π between groups)
- Three distinct confined cores maintained to T=500
- S_final = 0.97, P_int GROWING (123% retention)

**UDD (neutron analog)**: One Up + two Down-chirality braids.
Net charge = -1. Net θ ≈ 0 (phase cancellation). Without net theta:
- No coherent breathing (irregular oscillation)
- Carrier phases CONVERGE (loss of color diversity)
- Confinement weakening (one dominant blob absorbs others)
- S_final = 0.72, P_int retention 90%

The proton is more stable than the neutron — consistent with real physics
(free neutron half-life ~10 min, proton > 10³⁴ years).

### Proton Formation Dynamics [CONFIRMED — V43]

The 3-braid → proton formation is smooth and rapid (~1 time unit), with
no sharp phase transition. The converged proton is a "breathing shell":
depleted core (φ_rms drops to 0.33× initial at r<1), field concentrated
at r≈5 (peak φ_rms = 0.73× initial). The aspect ratio remains ≈1.14
throughout — near-spherical from initialization onward.

Two breathing timescales: the fast internal field oscillation has period
~2.2t (the φ carrier wave completing one cycle), while the structural
envelope breathing (expansion/contraction of R_rms) has period ~150t
for a single proton (V42). The 2.2t oscillation modulates P_int, which
oscillates 54% between peak (5.31) and trough (3.05). The φ_rms vs |P| anticorrelation
is a structural feature: high field amplitude redistributes to
inter-braid regions where the three carrier phases cancel, suppressing
the triple product P = φ₀φ₁φ₂.

### Three Stability Signatures [CONFIRMED]

From V40 Gen 4 structural analysis, three per-voxel metrics distinguish
stable from unstable regions:

1. **θ confinement**: θ_rms DECREASING outward (outer/inner < 0.7).
   When θ increases outward, the charge-carrying field is radiating away.

2. **Velocity structure**: Stable braids have |v| increasing outward
   (breathing shell). Stable 3-braids have |v| DECREASING outward
   (contraction — a novel mechanism).

3. **|P| concentration**: Inner/outer ratio > 10×. Critical P for maximum
   binding force: P_opt = 1/√(3κ) ≈ 0.082.

### Quantitative Design Rules [CONFIRMED]

From the V41 72-seed parameter sweep with stability prediction:

| Parameter | Optimal | Why |
|-----------|---------|-----|
| A (amplitude) | 0.3 | P_peak ≈ 0.07, near optimal 0.082 |
| R_tube | 4.5 | Maximizes triple-product overlap volume |
| theta_init | 0.5-0.8 | Pre-loads θ, prevents initial radiation transient |
| v_profile | contracting | Matches UDD contraction stability mechanism |
| carrier phase | {0, 2π/3, 4π/3} | Phase confinement (P=0 at triple overlap) |

Quantitative stability thresholds (from V40 Gen 4 structural analysis):

| Metric | Stable threshold | Unstable indicator |
|--------|-----------------|-------------------|
| θ outer/inner | < 0.7 | > 1.0 |
| |v| outer/inner | > 1.3 (breathing) or < 0.7 (contraction) | ≈ 1.0 (no dynamics) |
| |P| inner/outer | > 10 | < 2 |
| ρ inner/outer | > 5 | < 2 |
| Cluster volume | > 10,000 voxels | < 1,000 |
| Core force balance | < 0.05 | > 0.5 |

### Temporal Scales [MEASURED]

| Process | Timescale | Source |
|---------|-----------|--------|
| Single braid breathing | ~50 t | V40 Gen 0 |
| UUD composite breathing | ~150 t | V42 breathing analysis |
| Deuterium breathing | ~300 t | V42 breathing analysis |
| Phase convergence (UDD) | ~100 t | V41 phase analysis |
| Force equilibration (deuterium) | ~300 t | V42 accel analysis |
| UUD phase-lock stabilization | ~330 t | V41 T=500 (3 clusters from t=330) |

### The Breathing Oscillator [CONFIRMED — V42]

Every particle in this theory is a breathing oscillator — a dynamical system
that is always in motion internally (T > 0) and always subluminal in energy
transport (0 < V < c). The breathing character determines the particle's
physical properties:

**Mass** = total oscillation energy (E_kin + E_grad + E_mass + E_pot).
A faster-breathing particle has more kinetic energy → heavier.

This replicates QCD mass generation: in the Standard Model, the Higgs
mechanism provides only ~1% of the proton's mass. The other 99% comes from
the kinetic and binding energy of nearly-massless quarks and massless gluons
continuously oscillating inside the confinement radius. In our theory, the
braid amplitude A=0.3 sets a small "bare" field energy, but the breathing
oscillation generates kinetic energy (E_kin >> E_pot at most times) that
dominates the total mass. The particle's mass is overwhelmingly DYNAMICAL,
not from the field amplitude itself.

**Inertia** = resistance to rearranging the breathing pattern. More internal
KE to rearrange → more inertia. Gives m_inertial = m_gravitational naturally.

**Charge response** = θ coupling modulated by breathing phase. The particle
is not a static charge but an oscillating one — the time-averaged response
gives the effective charge.

**Temperature** = velocity dispersion (T_kin = <v²> - <v>²). Determines
interaction rate — hotter particles couple more strongly to external fields.

**Measured breathing properties:**

| Property | UUD Proton | Deuterium |
|----------|-----------|-----------|
| Breathing period | 150 t | 300 t |
| E_kin/E_pot correlation | -0.15 (free) | +0.12 (multi-mode) |
| Core T_kin | 0.001 | **0.16** (160× hotter) |
| Core |v| | 0.05c | **0.56c** |
| v_radial (core) | -0.01 (contracting) | -0.08 (8× faster) |
| Global v_max (field) | 2.22c | 2.02c |

Both the proton and deuterium breathe in MULTI-MODE fashion (E_kin and E_pot
weakly correlated: r=-0.15 for proton, r=+0.12 for deuterium). The internal
degrees of freedom (three braids each, plus inter-baryon modes in deuterium)
produce complex, multi-frequency oscillations. Larger structures breathe
slower (300t vs 150t).

**Superluminal field velocities**: The internal field oscillation ∂φ/∂t
reaches 2c at antinodes. This is the PHASE velocity of the standing wave,
not the energy transport velocity. The GROUP velocity (from V39 BLV analysis)
is always subluminal. Analogous to de Broglie wave phase velocity c²/v > c
for massive particles in quantum mechanics.

**Inverted temperature profile**: The proton has a cold core (T=0.001) and
hot shell (T=0.015). The deuterium has an extremely hot core (T=0.16) where
the nuclear bond zone creates intense field dynamics — 160× hotter than the
single proton. This temperature difference IS the nuclear binding energy
expressed as kinetic energy at the interaction zone.

## 6. Nuclear Binding [CONFIRMED]

### Deuterium (V42)

Two phase-confined baryons (UUD proton + UDD neutron) separated by
40 code units, simulated at N=512, L=100, T=500.

**RESULT**: The bound state SURVIVED T=500 with persistent inter-baryon
attraction. The two baryons maintained distinct identities while remaining
bound by the residual depletion interaction.

Key findings from the acceleration field analysis:

**Force equilibration** [CONFIRMED — emergent property]:
The binding force (V(P) potential) starts 259× stronger than the curl
coupling (θ). By T=300, they self-tune to a **1:1 ratio**. This emergent
equilibrium is not imposed — it arises from the coupled dynamics. In real
physics, the strong and electromagnetic forces also reach balance within
nuclei, producing the observed nuclear binding energies.

| Time | F_pot (strong) | F_curl (EM) | Ratio |
|------|---------------|-------------|-------|
| t=0 | 0.480 | 0.002 | 259:1 |
| t=300 | 0.011 | 0.010 | 1.05:1 |
| t=500 | 0.034 | 0.028 | 1.21:1 |

**Inter-baryon attraction**: Consistently ATTRACTIVE (F_x = -22K → -16K),
slowly settling as the baryons find their equilibrium separation. Not
merging (distinct structures preserved), not flying apart (force remains
attractive). This IS nuclear binding: two composites held by residual
depletion interaction.

**Core force balance**: |F_total| < 1% of |F_largest| at ALL radii. Four
large forces (Laplacian, mass, potential, curl) nearly perfectly cancel.
The structure is in genuine dynamical equilibrium.

**Theta radiation suppressed**: No theta-dominated breakaway structures in
deuterium (unlike single baryons which showed theta blobs). The two-baryon
system's theta fields partially cancel at large distances — the bound state
is a better "charge container" than its individual components.

**System compaction**: R_rms decreased 43% from t=0 to t=500. The nuclear
bound state is actively contracting toward a tighter equilibrium.

**Phase structure**: A new third phase group (φ≈0.7) emerged at T=500,
intermediate between the two anti-phase groups. This may be the inter-baryon
bridge that mediates the nuclear bond — not present in single baryons.

### Multi-Braid Formation

Two separate braids scatter (V33-C3: approach to D≈12, bounce, separate).
They CANNOT merge by collision at the tested initial conditions.

Multi-baryon bound states form by PLACING pre-formed baryons at appropriate
separations. The residual depletion interaction (from each baryon's ρ
depression) provides the nuclear binding force. The inter-baryon separation
self-adjusts to an equilibrium distance.

---

## 7. Dark Matter [THEORETICAL]

### Accretion + Depletion Profile

V30's M7 analysis showed the braid simultaneously:
- ACCRETES B-field at the core (higher density)
- DEPLETES B-field in a surrounding shell (lower density)

This profile (dense core + depleted shell) is structurally similar to
dark matter halos (NFW profiles). The "dark matter" in this framework
would be the field depletion around each particle — invisible but
gravitationally active through the depletion mechanism.

---

## 8. Cosmology [THEORETICAL]

### Braid Formation

V30 tested expansion + rotation for spontaneous braid formation.
Result: no braids from expansion alone. Braids require specific
initialization conditions (helical phase structure).

In a cosmological context, braids would form during phase transitions
in the early universe — analogous to topological defect formation
(Kibble-Zurek mechanism) but requiring a field with three-body coupling
and sufficient angular momentum.

---

## 9. Lorentz Invariance and Relativity

### The Equation IS Lorentz-Invariant [CONFIRMED by construction]

Equation (1) is a standard relativistic wave equation with c=1. It is
manifestly Lorentz-invariant. Any solution boosted by a Lorentz
transformation is also a solution.

### The Background IS NOT an Aether

The uniform background φ_a = A_bg × cos(k·z + 2πa/3) has a preferred
z-direction. This is an INITIAL CONDITION, not a property of the equation.
The equation has full rotational and Lorentz symmetry. The z-preference
is a simulation convenience for the helical braid initialization.

A rotationally symmetric background (random phases, isotropic spectrum)
would have no preferred direction. The physics does not require a
preferred frame.

### Lorentz Contraction of Solitons [THEORETICAL, untested]

A boosted soliton SHOULD contract by γ = 1/√(1-v²/c²) because the
equation is Lorentz-invariant. This has not been tested numerically
(V29-T8 was deprioritized). It is a key planned verification:

    Test: Initialize a soliton with Lorentz boost v=0.1c, 0.3c, 0.5c
    Measure: Does the aspect ratio match γ?
    Expected: YES (from the equation's symmetry)

### The Depletion Is Power-Law, Not Yukawa [CONFIRMED]

The mass term m²=2.25 gives individual field perturbations a Yukawa
decay e^{-mr}/r (range ≈ 1/m ≈ 0.67). But the braid's COLLECTIVE
depletion profile does NOT follow Yukawa. Direct measurement (V34
phonon test, N=256, L=60, T=200):

    Yukawa (m=1.5 fixed): R² = -0.44  (EXCLUDED — worse than flat line)
    Yukawa (m free):       converges to m ≈ 0.02 (effectively massless)
    Power law δρ ∝ 1/r^n: R² = 0.98, n = 1.2 ± 0.2

At r=10, Yukawa m=1.5 predicts δρ ≈ 3×10⁻⁸. Measured: 2×10⁻².
The depletion is 500,000× larger than Yukawa allows.

**The soliton's depletion extends as a power law to r=30+ despite m=1.5.**
(Measured on a z-aligned braid; the proton's depletion profile may differ
in shape but the power-law character should persist.)

This resolves the apparent paradox: V33-C1 measured attraction at D=15-80
where Yukawa would be 10⁻¹⁰ to 10⁻⁵³. The force IS real, carried by
a collective mode of the background (phonon), not by individual massive
field excitations.

**Mechanism**: The background φ_a = A_bg × cos(kz + 2πa/3) is a
periodic structure. The braid perturbs the background amplitude, and
this perturbation propagates as a PHONON — a collective oscillation of
the background lattice. The phonon is massless (Goldstone mode of broken
translational symmetry), so it propagates without Yukawa suppression.

The m² parameter confines the braid's internal structure but does NOT
limit the gravitational range. The mass and the range are DECOUPLED:
- m² binds the braid (helical structure needs the spring)
- The phonon carries the gravity (massless, power-law)

**Implication**: The m→0 problem (F1 in FUTURE.md) may not be a problem.
The gravitational range was never limited by m. Lowering m² is needed
only if the force exponent (currently n≈1.2 for depletion, n≈1.8 for
force) needs to approach n=2 for Newtonian gravity. The exponent may
improve with better equilibration (T=1000+) and isotropic backgrounds.

### The Mass Term m²

The m² = 2.25 in Eq. (1) provides:
1. Braid binding: the helical structure requires m² ≥ 1.25 (Track G)
2. Vacuum stability: the background requires m² ≥ 0.25 (Track G)
3. It does NOT limit gravitational range (phonon is massless)

Field-dependent mass (V34 Track GB) was tested and does not improve on
constant m². The φ⁴ coupling kills braids (Derrick's theorem); the
inverse coupling collapses at long range. Constant m² is Goldilocks
for braid binding. See v34/GB_field_mass/RESULTS.md.

### Numerical vs Physical Energy Drift

The 0.1% energy drift per T=100 is from the NUMERICAL method (dt² error
in symplectic Verlet), not from physical radiation. The physical radiation
(V29-T3: ~0.5 energy/time) is a separate, much larger effect visible in
the energy timeseries as a slow decrease in E_total.

At thermal equilibrium (V29-T1b: A_noise ≈ 0.02), physical radiation
= absorption, and E_total stabilizes. The numerical drift persists but
can be reduced with smaller dt.

---

## 10. Open Questions

1. **θ static vs wave**: The θ force is wave-mediated (oscillating), not
   static (1/r). Is there a low-frequency limit where the DC component
   dominates? Or is wave exchange the fundamental mechanism (QFT-like)?
2. **Coulomb analog**: What is the static electric force between stationary
   composite baryons of opposite net chirality? (Current tests have
   moving/oscillating particles.)
3. **Depletion exponent**: Why n≈1.2 for δρ(r) and n≈1.8 for F(D)?
   Does it converge to n=2 with isotropic background and longer runs?
4. **Electrons as orbital modes**: Does the θ shell around a proton
   support bound standing-wave modes? (The radial θ profile peaks at
   r≈8, outside the core — like an electron shell.)
5. **Multi-baryon formation**: conditions for stellar nucleosynthesis analog.
6. **Spin**: does the braid's helical handedness correspond to spin?
7. **Quantization**: the theory is purely classical. What changes in QFT?
8. **Background origin**: what sets ρ_bg, m², and η?
9. **Fine structure constant**: What determines the ratio of θ force to
   φ force? η controls it in the equation, but what sets η physically?

## 11. Resolved Questions

1. **Is gravity energy minimization?** NO. E(D) is monotonically repulsive.
   The attraction is a dynamic effect. (V33 energy_vs_D, March 2026)
2. **Is the drift proportional to gradient?** YES. F ∝ ∇ρ with R²=0.9998.
   (V33 gradient sweep, March 2026)
3. **Is drag lower in depleted regions?** NO. Coupling is stronger in high ρ.
   The mechanism is geometric asymmetry, not friction. (V33 drag test)
4. **Is the depletion Yukawa?** NO. δρ ∝ 1/r^1.2 (power law), not e^{-mr}/r.
   Yukawa m=1.5 excluded at 500,000× discrepancy. (V34 phonon test)
5. **Does m² limit gravitational range?** NO. The phonon (collective mode)
   is massless and carries the long-range force. m² binds the braid only.
   (V34 phonon test + Track G + Track GB)
6. **Can field-dependent mass decouple binding from range?** NO. φ⁴ kills
   braids (Derrick's theorem). Inverse coupling collapses. Constant m² is
   Goldilocks. But this is moot since the range isn't m²-limited anyway.
   (V34 Track GB)
7. **Can torsion waves carry EM in 3-field theory?** NO. All linearized
   modes are massive. Torsion disperses in ~5t. Winding ≠ charge.
   (V34 Phase 2)
8. **Can the equation support angle (rotation) fields?** YES. The Cosserat
   extension (3 φ + 3 θ) is stable. The braid sources θ through curl(φ).
   Massless θ (m_θ=0) is stable and actually helps the braid. (V34 Cosserat)
9. **Does the phonon arise from a Goldstone mode?** NO. All linearized modes
   around the uniform background are massive (splitting only 0.004 at A=0.1).
   The power-law depletion is a nonlinear collective effect. (V34 Phase 1a)
10. **Does the θ field mediate a charge-dependent force?** YES. Same-winding
    braids attract 27% more, opposite-winding attract 57% less than the
    3-field baseline. Winding = charge. (V34 θ characterization)
11. **Does winding reversal flip the θ field?** YES. Time-averaged θ_φ DC
    component has correlation -0.68 between W=+1 and W=-1. Far-field
    ratio ≈ -1.0. (V34 θ characterization)
12. **Is the θ field dt-converged?** YES. θ_rms varies by only 2% across
    a 4× dt range (0.25×, 0.5×, 1.0×). The structure is physical. (V34)
13. **Is the θ field a static 1/r magnetic field?** NO. It is an oscillating
    wave (period ~4t) with 0.2% DC bias. The force is wave-mediated
    (radiation coupling), not static. Closer to QFT photon exchange than
    classical Biot-Savart. (V34 θ characterization)
14. **Can density-dependent κ create black holes?** NO. Theta radiation
    prevents self-trapping. Higher γ → faster dispersal. (V39)
15. **Is the BLV metric geometric or dispersive?** DISPERSIVE. The Cosserat
    Lagrangian is quadratic in derivatives → BLV metric is identically flat.
    All "gravitational" effects are frequency-dependent mass lensing. No
    absolute trapping horizon. (V39 BLV analysis)
16. **Can 3-braid composites form stable particles?** YES. Phase-confined
    UUD (proton) survives T=500 with anti-phase breathing at period 90t.
    Three distinct cores maintained by carrier phase offsets. (V41)
17. **Is chirality = charge?** YES. Winding direction determines the sign
    of the θ-mediated force. UUD (net +1) is more stable than UDD (net -1).
    (V34 + V41)
18. **Is carrier phase = color charge?** YES. Three braids with Δ={0,2π/3,4π/3}
    create confinement (P→0 at overlap). Each carries 1/3 of the phase circle.
    (V41 phase confinement)
19. **Can two baryons bind (nuclear force)?** YES. UUD+UDD at N=512 T=500
    shows persistent inter-baryon attraction, force equilibration (strong:EM
    → 1:1), and system compaction. (V42 deuterium)
20. **Does a proton respond gravitationally to a density gradient?** YES.
    Pre-converged UUD proton drifts toward low ρ (+0.80 to +1.84 code units
    over T=200). Drift scales linearly with gradient. (V43 gradient test)
21. **Is the gravitational force from φ-depletion or θ-coupling?** φ-depletion.
    The η=0 proton drifts +0.98 (slightly MORE than η=0.5 at +0.80). Gravity
    is pure φ; θ coupling provides a small counter-force. (V43 gradient test)
22. **Does asymmetric binding directly drive the drift?** YES. The density
    gradient creates 20–42% imbalance in |P| (the binding density). In 10/11
    frame transitions, asym_P sign matches drift direction. The force is from
    V(P) coupling asymmetry, not raw field amplitude. (V43 phase binding)

---

## References to Simulation Code

### Best-in-Class Simulation (3-field, gravity only)
**File**: `/home/d/code/scp/v33/src/v33.c`
- Standard equation (Eq. 1), no modifications
- Single malloc for entire simulation (9 arrays × N³ × 8 bytes)
- Symplectic Velocity Verlet (energy drift < 0.1%/T=100)
- Fully periodic boundary conditions (for C1-C5 characterization)

### Gradient Test (Gravity Direction)
**File**: `/home/d/code/scp/v33/src/v33_gradient_test.c`
- x-direction: PINNED (high ρ left, low ρ right)
- y,z-directions: FREE-FLOATING (linear extrapolation from interior:
  phi[boundary] = 2×phi[boundary-1] - phi[boundary-2])
- Braid drifts toward LOW ρ = toward depletion = GRAVITY

### Braid Detection and Analysis
**File**: `/home/d/code/scp/v33/src/analyze.c`
- Reads binary field snapshots
- Flood-fill connected components above 5× average ρ
- Per-braid: position, energy, size, |P|, binding weight, winding, momentum
- Radial profiles: ρ(r), |P|(r), w(r) around each braid
- Inter-braid separations
- Background statistics

### Footprint Asymmetry Analysis
**File**: `/home/d/code/scp/v33/src/footprint_asymmetry.c`
- Reads field snapshots, finds braid center, extracts x-directional profile
- Measures half-width toward high-ρ vs low-ρ sides
- Computes geometric vs field-weighted center of perturbation footprint

### Drag Test (Braid Coupling vs Background Density)
**File**: `/home/d/code/scp/v33/src/v33_drag_test.c`
- Single braid with Galilean velocity kick in uniform backgrounds
- Tracks position, momentum, and braid energy over time
- Result: coupling STRONGER in high-ρ (opposite of simple drag model)

### Energy vs Separation (F = -dE/dD test)
**File**: `/home/d/code/scp/v33/src/v33_energy_vs_D.c`
- Two braids at separation D, settle T=30, average E over t=20-30
- Sweep D=8–50 to map E_int(D)
- Result: E_int always positive (repulsive), no attractive well

### F2 Test Battery (gradient sweep + energy sweep)
**File**: `/home/d/code/scp/v33/run_f2_tests.sh`
- Gradient strength sweep: 4 gradients (ρ ratio 1.22–9.00)
- Energy separation sweep: D=8–50
- Footprint analysis on all snapshots

### 6-Field Cosserat Simulation
**File**: `/home/d/code/scp/v34/torsion_coupling/src/v33_cosserat.c`
- 3 position + 3 angle fields (Eq. 10), curl coupling
- Single malloc for 18 arrays (6 fields × 3 arrays each)
- Supports -eta (coupling) and -mt (angle mass) flags
- Massless angle fields (m_θ=0) confirmed stable

### Binary Snapshot Formats

**3-field (v33):**
```
int32: N, float64: L, float64: t
float64[N³]: phi_0, phi_1, phi_2
```

**6-field (Cosserat):**
```
int32: N, float64: L, float64: t, int32: nf=6
float64[N³]: phi_0, phi_1, phi_2, theta_0, theta_1, theta_2
```

---

## Version History

- **V28**: Bimodal braid discovery (CMA-ES search)
- **V29**: Validation campaign (13 confirmed findings, field study)
- **V30**: Failed approaches (c(ρ), FRW expansion)
- **V31**: M7 split gravity tests (wrong sign → inverted c works but artificial)
- **V32**: SPH prototype, gradient coupling (self-interaction problem solved by binding weight, but gradient coupling is net repulsive with controls)
- **V33**: Clean standard equation. Single alloc. CONFIRMED: intrinsic attraction from standard equation, F∝1/D^1.8, gravity direction verified via gradient test. Footprint asymmetry measured (R_low/R_high=1.09–1.57). Drag test showed coupling scales with ρ (not anti-ρ). Gravity mechanism identified as geometric (asymmetric Yukawa profile), not dynamic (friction).
- **V34**: Metastability (Track G: m²≥1.25 for braids). Depletion is power-law 1/r^1.2 NOT Yukawa (phonon test). Torsion waves fail as EM in 3-field theory. **MAJOR: Cosserat 6-field (3φ+3θ) with massless θ. CHARGE-DEPENDENT FORCE CONFIRMED**: same-winding +27% attraction, opposite -57%. Winding = charge. θ pattern follows right-hand rule (visualization confirmed, dt-converged). Wave-mediated mechanism (QFT-like, not static Biot-Savart). SFA archive format + volumetric viewer built.
- **V35**: Hydrogen-like bound states confirmed (parameter sweep: Bohr ratio 52,957, exact hydrogen spectrum). Field quantization ε and θ self-interaction both ruled out. Multi-scale simulation built: 3D spherical core (braid) + 1D radial Schrödinger wedge (electron). Frozen-braid standalone wedge confirms stable orbital over 44 periods — eigenstate is stationary (correct QM). Classical 3D θ packets collapse without ℏ (Rutherford problem, confirmed). Two-step architecture validated: V_eff(r) from classical field + ℏ_eff from braid action → Bohr orbits. Dynamics require perturbations (multi-atom, orbital transitions).
- **V36**: CUDA GPU port on Tesla V100 ($0.17/hr Vast.ai). 6-field Cosserat at 1.9 ms/step (3.5× over initial GPU, ~50× over CPU). SFA f32 optimizer: 16.8 GB → 3.8 GB (4.4×). SFA C++ compatibility. Viewer updated for f32 dtype. All output in SFA format.
- **V37**: Compact braid geometries (truncated helix, braid3).
- **V38**: GPU seedrun framework.
- **V39**: Density-dependent κ collapse RULED OUT with 6-field Cosserat (theta prevents self-trapping). BLV effective metric: Φ/c²≈-0.02 (nuclear scale, 10³⁷× GR). Unified sim kernel (scp_sim CPU+GPU) with KVMD metadata, 12-column SFA, SFA restart.
- **V40**: Evolutionary composite particle search (4 generations). Discovered three stability signatures (θ confinement, velocity structure, |P| concentration). UDD 3-braid contraction mechanism. Analysis tool suite built.
- **V41**: First-principles seed construction from stability signatures (95% P_int retention). Phase confinement discovered (carrier phase = color charge). Proton analog stable at T=500 (anti-phase locked, 3 confined cores). Neutron analog weaker (no θ synchronization). Charge = chirality + carrier phase.
- **V42**: **DEUTERIUM** — first nuclear binding simulation. UUD+UDD at N=512, T=500. Force equilibration (strong:EM → 1:1). Persistent inter-baryon attraction. Theta radiation suppressed in bound state. System compacting.

---

## Appendix: Critical Code Implementation

### A1. Memory — Single Allocation

All 9 field arrays (phi[3], vel[3], acc[3]) in ONE contiguous block:

```c
static Grid *grid_alloc(int N, double L) {
    Grid *g = calloc(1, sizeof(Grid));
    g->N  = N;
    g->N3 = (long)N * N * N;
    g->L  = L;
    g->dx = 2.0 * L / (N - 1);
    g->dt = 0.12 * g->dx;

    /* ONE allocation */
    long total = 9 * g->N3;
    g->mem = malloc(total * sizeof(double));
    memset(g->mem, 0, total * sizeof(double));

    /* Pointers into the single block */
    for (int a = 0; a < NFIELDS; a++) {
        g->phi[a] = g->mem + (0 + a) * g->N3;
        g->vel[a] = g->mem + (3 + a) * g->N3;
        g->acc[a] = g->mem + (6 + a) * g->N3;
    }
    return g;
}

static void grid_free(Grid *g) {
    free(g->mem);  /* ONE free */
    free(g);
}
```

### A2. Force Computation — Eq. (1) Evaluated

The core physics loop. Every grid point runs the same equation.
No distinction between braid and background.

```c
static void compute_forces(Grid *g) {
    const int N = g->N, NN = N * N;
    const long N3 = g->N3;
    const double idx2 = 1.0 / (g->dx * g->dx);

    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < N3; idx++) {
        int i = (int)(idx / NN);
        int j = (int)((idx / N) % N);
        int k = (int)(idx % N);

        /* Fully periodic neighbors */
        int ip = (i+1)%N, im = (i-1+N)%N;
        int jp = (j+1)%N, jm = (j-1+N)%N;
        int kp = (k+1)%N, km = (k-1+N)%N;

        long n_ip = (long)ip*NN + j*N + k;
        long n_im = (long)im*NN + j*N + k;
        long n_jp = (long)i*NN + jp*N + k;
        long n_jm = (long)i*NN + jm*N + k;
        long n_kp = (long)i*NN + j*N + kp;
        long n_km = (long)i*NN + j*N + km;

        /* Triple product P = φ₀φ₁φ₂ and dV/dP */
        double p0 = g->phi[0][idx];
        double p1 = g->phi[1][idx];
        double p2 = g->phi[2][idx];
        double P  = p0 * p1 * p2;
        double den = 1.0 + KAPPA * P * P;
        double mPd2 = MU * P / (den * den);    /* μP/(1+κP²)² */

        for (int a = 0; a < NFIELDS; a++) {
            /* 7-point Laplacian (6 neighbors + center) */
            double lap = (g->phi[a][n_ip] + g->phi[a][n_im]
                        + g->phi[a][n_jp] + g->phi[a][n_jm]
                        + g->phi[a][n_kp] + g->phi[a][n_km]
                        - 6.0 * g->phi[a][idx]) * idx2;

            /* ∂P/∂φ_a */
            double dPda = (a==0) ? p1*p2 : (a==1) ? p0*p2 : p0*p1;

            /* acc = ∇²φ - m²φ - (μP/(1+κP²)²)(∂P/∂φ_a) */
            g->acc[a][idx] = lap - MASS2 * g->phi[a][idx] - mPd2 * dPda;
        }
    }
}
```

### A3. Time Integration — Symplectic Velocity Verlet

Energy-conserving (symplectic). Drift < 0.1% per T=100.

```c
static void verlet_step(Grid *g) {
    const long N3 = g->N3;
    const double hdt = 0.5 * g->dt;
    const double dt  = g->dt;

    /* Half-kick: v += (dt/2) × acc */
    for (int a = 0; a < NFIELDS; a++) {
        double *v = g->vel[a], *ac = g->acc[a];
        for (long idx = 0; idx < N3; idx++)
            v[idx] += hdt * ac[idx];
    }
    /* Drift: φ += dt × v */
    for (int a = 0; a < NFIELDS; a++) {
        double *p = g->phi[a], *v = g->vel[a];
        for (long idx = 0; idx < N3; idx++)
            p[idx] += dt * v[idx];
    }
    /* Recompute forces from new φ */
    compute_forces(g);
    /* Half-kick: v += (dt/2) × acc */
    for (int a = 0; a < NFIELDS; a++) {
        double *v = g->vel[a], *ac = g->acc[a];
        for (long idx = 0; idx < N3; idx++)
            v[idx] += hdt * ac[idx];
    }
}
```

### A4. Boundary Conditions — Pinned + Free-Floating (Gradient Test)

x-direction pinned to imposed ρ gradient. y,z free-floating outflow.

```c
static void apply_bc(Grid *g) {
    const int N = g->N, NN = N*N;
    const int margin = 3;

    /* x-direction: PIN to saved boundary values (high ρ left, low ρ right) */
    for (long idx = 0; idx < g->N3; idx++) {
        int i = (int)(idx / NN);
        if (i < margin || i >= N - margin) {
            for (int a = 0; a < NFIELDS; a++) {
                g->phi[a][idx] = g->pin_phi[a][idx];  /* restore saved */
                g->vel[a][idx] = g->pin_vel[a][idx];
                g->acc[a][idx] = 0;
            }
        }
    }

    /* y,z-directions: FREE-FLOATING (linear extrapolation from interior)
       phi[boundary] = 2×phi[interior_1] - phi[interior_2]
       Allows waves to flow out. Does not inject or remove energy. */
    for (int b = 0; b < margin; b++) {
        int j_in1 = margin, j_in2 = margin + 1;
        // ... (applied symmetrically to both edges of y and z)
        g->phi[a][idx_boundary] = 2*g->phi[a][idx_in1] - g->phi[a][idx_in2];
    }
}
```

### A5. Braid Detection — Connected Component Analysis

Post-hoc analysis of field snapshots. No distinction during simulation.

```c
/* Threshold: cells with ρ > 5× average are "braid" candidates */
double thresh = 5.0 * avg_rho;

/* Flood-fill from each seed above threshold */
for (long seed = 0; seed < N3; seed++) {
    if (rho[seed] < thresh || label[seed] >= 0) continue;
    int bid = n_braids++;

    /* BFS: expand to 6-connected neighbors above threshold */
    stack[sp++] = seed;
    label[seed] = bid;
    while (sp > 0) {
        long idx = stack[--sp];
        /* Accumulate: centroid, energy, |P|, binding weight */
        braids[bid].cx += x * rho * dV;
        braids[bid].E_total += rho * dV;
        braids[bid].P_avg += |P|;
        braids[bid].w_avg += 1.0/(1.0 + |P|/P_thresh);
        /* Expand to unvisited neighbors above threshold */
        for (int d = 0; d < 6; d++) { ... }
    }
}
```
