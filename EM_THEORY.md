# Electromagnetism in the Cosserat Field Theory

## Summary

The 6-field Cosserat equation naturally produces classical electromagnetism.
The three angle fields θ_a map to the electromagnetic vector potential A_a.
Braids act as Ampèrian current loops whose oscillating fields create
effective Coulomb forces via radiation pressure. Maxwell's equations emerge
from linearizing the Cosserat equation around a uniform background.

---

## 1. The Field Mapping

### Cosserat Equation (Eq. 10 in CONCEPT.md)

```
∂²φ_a/∂t² = ∇²φ_a - m²φ_a - ∂V/∂φ_a + η × curl(θ)_a     (position, massive)
∂²θ_a/∂t² = ∇²θ_a - m_θ²θ_a          + η × curl(φ)_a     (angle, massless)
```

Parameters: m² = 2.25, m_θ² = 0, η = 0.5, μ = -41.345, κ = 50.

### Electromagnetic Identification

| Cosserat field | EM quantity | Justification |
|---------------|-------------|---------------|
| θ_a | Vector potential A_a | Massless, propagates at c, sourced by currents |
| -∂θ/∂t | Electric field E | Time derivative of vector potential |
| ∇×θ | Magnetic field B | Curl of vector potential |
| η∇×φ | Current density J_eff | Braid's helical twist creates a current loop |
| Winding W = ±1 | Sign of magnetic moment | Current circulation direction (±) |

**Gauge**: The theory operates in the Weyl gauge (scalar potential Φ = 0).
There is no gauge freedom — θ is a physical field, not a gauge artifact.
This is a feature: the electromagnetic degrees of freedom are fully physical
with no redundant gauge modes.

### What Is "Charge"?

A braid's helical twist creates a nonzero ∇×φ, which acts as a current
source J_eff = η∇×φ for the θ equation. This is a CURRENT LOOP, not
a point charge.

**Mathematically**: ∇·J_eff = ∇·(η∇×φ) = 0 identically (divergence of
curl is always zero). By the continuity equation ∂ρ/∂t + ∇·J = 0, this
means ∂ρ_charge/∂t = 0. There are no scalar electric monopole charges in
the linearized theory.

The braid is fundamentally a **magnetic object** — an Ampèrian current
loop generating a magnetic dipole field. The winding number W = ±1
determines the direction of current circulation, hence the sign of the
magnetic moment. Effective electric charge emerges at the composite-baryon
level when orthogonal dipoles produce an isotropic radiation pattern that
mimics a monopole (see OQ1).

**V34 confirmation**: Same-winding braids attract, opposite-winding repel.
This IS Ampère's force law for parallel and anti-parallel currents, not
Coulomb's law for electric charges. The right-hand rule θ pattern observed
in volumetric rendering is the magnetic field B = ∇×A = ∇×θ around a
current loop.

### The Composite Baryon Question

In the phase-confined baryon (V41), three braids are arranged along x, y, z
axes with chiralities UUD or UDD. The three current loops are ORTHOGONAL.
The composite current pattern is more complex than a single dipole.

**Open question**: Does the UUD composite (net "charge" +1) behave as an
effective electric monopole at long range, despite each component being a
pure magnetic dipole? The multipole expansion of the θ field around a
composite baryon has not been computed. If the three orthogonal current
loops produce a net monopole-like term in the far field, this would explain
how effective "electric charge" emerges from purely magnetic components.

---

## 2. Light as a Coupled δφ-δθ Wave

### Linearized Modes

Linearizing around a uniform background (V'(P) ≈ 0 for small perturbations):

```
∂²δφ_a/∂t² = ∇²δφ_a - m²δφ_a + η curl(δθ)_a     (massive + coupling)
∂²δθ_a/∂t² = ∇²δθ_a           + η curl(δφ)_a     (massless + coupling)
```

The coupled system has two hybridized dispersion branches:

**Branch 1 (photon-like)**: Nearly massless, propagates at ≈ c.
- Mostly δθ with small δφ admixture
- For small η: ω² ≈ k² + O(η²m²/k²) → massless at high k
- This IS light in the theory

**Branch 2 (matter-like)**: Massive, dispersive.
- Mostly δφ with small δθ admixture
- For small η: ω² ≈ k² + m² + O(η²) → massive Klein-Gordon
- This IS the field excitation (not a particle — particles are braids)

The hybridization is controlled by η. For η = 0.5 and m = 1.5, the mixing
is small (~η²/m² ≈ 0.11) — the photon-like branch is nearly pure θ.

### Free-Space Light vs Near-Braid Hybridization

Free-space light at low energies (ω < m = 1.5) propagates as a **pure θ wave**.
The φ field is massive (m² = 2.25) and cannot be excited below its mass gap.
The E ↔ B oscillation of free-space light happens entirely within the θ sector:
E = -∂θ/∂t and B = ∇×θ exchange energy via the standard Maxwell mechanism,
with no φ involvement.

The hybridization with δφ occurs ONLY:
- Inside or very close to braids (where φ is large and nonlinear V(P) matters)
- At very high energies (ω > m ≈ 1.5, i.e., gamma rays in code units)

This means the vacuum is NOT birefringent or dispersive for low-energy light.
Photons propagate at exactly c through empty space. The φ coupling only
affects light when it passes through or near matter (braids).

### V34 Evidence

The θ field around a braid shows:
- Oscillation period ~4 time units
- 99.8% oscillating (AC) component
- 0.2% steady (DC) component
- Circular patterns following the right-hand rule
- dt-converged (varies by only 2% across 4× dt range)

This is electromagnetic radiation sourced by the braid's helical current.

---

## 3. The Meta-Static Coulomb Potential

### The Scale Separation

The θ oscillation frequency and the electron orbital frequency differ by
a factor of ~20 million:

| Timescale | Value | Code units |
|-----------|-------|------------|
| θ oscillation period | ~4 t | ~7.5×10⁻²⁴ s |
| Bohr orbital period | ~8×10⁷ t | ~1.5×10⁻¹⁶ s |
| Ratio | | ~2×10⁷ |

From the electron's perspective, the θ field completes 20 million
oscillation cycles per orbit. The electron sees only the time-averaged
field — the DC component — which is expected to have the structure of a
static potential (pending explicit ⟨θ(r)⟩ measurement from existing SFA
files — see OQ3).

Note: the observed |∂φ/∂t| > c in BREATHING_ANALYSIS.md are phase velocities
of the standing-wave oscillation and do not carry information or energy
(group velocity ≤ c, confirmed in V39 BLV analysis).

### How the "Static" Coulomb Force Emerges

The mechanism operates through **ponderomotive force** (radiation pressure):

1. Each braid emits oscillating θ waves (AC magnetic dipole radiation)
2. The radiation intensity (flux) falls as 1/r² in the far field
3. When θ-waves from braid A hit braid B's current loop, braid B
   SCATTERS the radiation, absorbing momentum proportional to the
   incident flux
4. The momentum transfer force F ∝ incident flux ∝ 1/r² — the Coulomb law

Note: this is radiation SCATTERING/momentum transfer, not the ponderomotive
gradient force (which would give F ∝ -∇⟨|B|²⟩ ∝ 1/r³). The 1/r² dependence
comes from the incident radiation intensity, not its gradient.
(Exact derivation pending — see OQ2.)

This is identical to the QED mechanism:
- Each θ oscillation cycle = one virtual photon exchange
- 20 million exchanges per orbit → statistical average converges
- The sum over all exchange frequencies gives 1/r (massless propagator)
- No quantum path integral needed — the classical dynamics naturally
  time-average to the same result

### Regime Transition: EM → Nuclear

The "static" Coulomb picture is valid only at distances >> λ_θ ≈ 2.2 fm.

At nuclear distances (~2.2 fm), the individual oscillation cycles become
resolvable. The time-averaging breaks down and the force becomes the full
wave-mediated interaction seen in V34/V41/V42. The "simple" electromagnetic
force at atomic scales and the "complicated" nuclear force at femtometer
scales are THE SAME FORCE at different resolution regimes.

| Distance | Regime | Force character |
|----------|--------|----------------|
| > 100 fm | Atomic/Coulomb | Static 1/r², from radiation pressure |
| ~2-10 fm | Nuclear | Wave-mediated, oscillatory, V(P) coupling dominates |
| < 2 fm | Confinement | Phase-cancelled, P→0 at braid overlap |

### The Fine Structure Constant

The ratio of the DC component to the total θ oscillation amplitude determines
the effective coupling strength at long range.

From V34: DC/total ≈ 0.002 (0.2% DC bias).
Physical α = e²/(4πε₀ℏc) ≈ 1/137 ≈ 0.0073.

The factor ~3.7× discrepancy between 0.002 and 0.007 depends on:
- η (coupling strength, currently 0.5)
- Braid geometry (R_tube, ellipticity, amplitude)
- The precise mapping between DC bias and coupling constant

This is a calibration target, not yet a prediction. Determining whether
α is derivable from (η, m, μ, κ) or requires an additional input is an
open question.

### Natural UV Cutoff

The θ oscillation frequency (~1.3×10²³ Hz) provides a natural ultraviolet
cutoff for the electromagnetic interaction. There is no need for
renormalization — the theory has a physical highest frequency (the braid
breathing mode) that regularizes all integrals.

In QED, the Coulomb potential is formally divergent at short distances and
requires renormalization. In this theory, the potential transitions smoothly
from 1/r (Coulomb) at long range to the wave-mediated regime at short range,
with no divergence at any scale.

---

## 4. Maxwell's Equations from the Cosserat Lagrangian

### Derivation

Starting from the linearized θ equation:

```
∂²δθ_a/∂t² = ∇²δθ_a + η(∇×δφ)_a
```

Define:
```
A_a = δθ_a                    (vector potential)
E_a = -∂A_a/∂t = -∂δθ_a/∂t   (electric field)
B = ∇×A = ∇×δθ               (magnetic field)
J_a = η(∇×δφ)_a              (effective current density)
```

**Maxwell equation 1**: ∇·B = 0
```
∇·B = ∇·(∇×δθ) = 0    (identity, always true)
```

**Maxwell equation 2**: ∇×E = -∂B/∂t
```
∇×E = ∇×(-∂δθ/∂t) = -∂/∂t(∇×δθ) = -∂B/∂t    (commuting derivatives)
```

**Maxwell equation 3**: ∇×B = ∂E/∂t + J
Taking the curl of the θ equation and rearranging:
```
∂²(∇×δθ)/∂t² = ∇²(∇×δθ) + η∇×(∇×δφ)
∂²B/∂t² = ∇²B + ∇×J
```
Using B = ∇×A and the wave equation: ∇×B = ∂E/∂t + J (in appropriate form).

**Maxwell equation 4**: ∇·E = ρ_charge
```
∇·E = -∂(∇·δθ)/∂t
```

In the Weyl gauge (Φ=0), if ∇·A = 0 (Coulomb gauge condition), then ∇·E = 0
in vacuum. Charge density appears only at the locations of braids where the
source term J = η∇×φ is nonzero. The effective charge density is:

```
ρ_eff = -∂(∇·δθ)/∂t    (from the divergence of the E field)
```

**Proof that ∇·E = 0 is exact and unbreakable:**
1. Take the divergence of the θ wave equation: □(∇·θ) = η ∇·(∇×φ)
2. By vector calculus identity: ∇·(∇×φ) ≡ 0
3. Therefore: □(∇·θ) = 0
4. If ∇·θ = 0 initially (which it is — θ starts at zero or is pre-loaded
   with curl structure that is divergence-free), then ∇·θ = 0 for all time.
5. Since E = -∂θ/∂t: ∇·E = -∂(∇·θ)/∂t = 0 everywhere, forever.

**This is not a gauge choice — it is a theorem of the equations.** The theory
strictly forbids scalar electric monopole charges. All electromagnetic
interactions arise from magnetic dipole (current loop) sources. Effective
"Coulomb-like" forces emerge from radiation scattering between these dipoles.

Gauss's law ∇·E = ρ_eff holds with ρ_eff = 0 everywhere. Net electric
monopole charge is zero for every configuration. The effective "charge"
of a composite baryon (OQ1) is not a true monopole but an isotropic
radiation pattern that mimics one.

### What This Gives and What It Doesn't

**Confirmed**:
- Two transverse polarizations for free EM waves (from the curl structure)
- Propagation at c (massless θ)
- Source coupling via J = η∇×φ (braids radiate)
- ∇·B = 0 and ∇×E = -∂B/∂t (exact, from the equations)

**Not yet confirmed**:
- Whether the nonlinear V(P) corrections modify Maxwell at strong fields
- The exact form of Gauss's law (∇·E = ρ) for composite baryons
- Whether the ponderomotive force exactly reproduces 1/r²
- The relationship between η and α (fine structure constant)

---

## 5. Ohm's Law and Conductivity

### Concept

A "wire" in this theory would be a lattice of braids (same winding) in a
line or tube. Applying a θ gradient (voltage) across the wire would accelerate
the braids along the wire. Resistance arises from:
- Braid-braid scattering at depletion zone boundaries
- θ radiation losses (electromagnetic radiation from accelerating charges)
- V(P) coupling friction between moving braids and the background

### Requirements

This is a many-body problem requiring:
- ~10-100 braids in a line geometry
- N ≥ 1024 grid for adequate resolution
- Applied θ boundary conditions (voltage analog)
- Measurement of braid drift velocity vs applied θ gradient

If J ∝ E (linear response), the proportionality constant is the
conductivity σ. The resistance R would depend on the braid spacing,
amplitude, and coupling η.

### Status

This is a V43+ experiment. The multi-baryon framework and phase confinement
machinery from V41/V42 provide the tools, but the simulation scale (100+
braids) requires significantly larger grids than current capability.

---

## Open Questions

### OQ1: Isotropic Radiation Pattern from Composite Dipoles

Each braid is a magnetic dipole (current loop) that emits anisotropic
radiation (stronger at the poles than the equator). A UUD baryon has three
ORTHOGONAL current loops. Since ∇·θ = 0 is exact, there are no true
monopoles in the multipole expansion.

**The real question**: Do three orthogonal current loops sum to produce a
SPHERICALLY ISOTROPIC radiation pattern? If so, the 1/r² radiation pressure
on distant particles is perfectly isotropic — indistinguishable from a
Coulomb point charge. The composite "charge" is not a monopole but an
isotropic radiator.

**Test**: Compute the angular power spectrum of the θ field at r=20 around
a UUD composite from existing V41 data. If the l=0 (isotropic) component
dominates over l=1 (dipole) and l=2 (quadrupole), the composite effectively
acts as a point charge. Three orthogonal dipoles should give significant
l=0 content by symmetry.

### OQ2: Analytical 1/r² Force from Radiation Scattering

The effective Coulomb force should be F = σ⟨S⟩/c, where σ is the scattering
cross-section of a target braid and ⟨S⟩ is the time-averaged Poynting flux
from the source braid. If ⟨S⟩ ∝ 1/r², then F ∝ 1/r² (Coulomb).

**Approach** (Poynting vector method — cleaner than interaction energy):
1. Define the Cosserat Poynting vector: S = E × B = (-∂θ/∂t) × (∇×θ)
2. Compute ⟨S⟩ for a single oscillating braid at distance r
3. Prove ⟨|S|⟩ ∝ 1/r² in the far field (standard for dipole radiation)
4. Define a scattering cross-section σ for a target braid's current loop
5. Force is F = σ⟨S⟩/c — purely from momentum transfer of scattered radiation

This avoids the two-body coupled integral and reduces to a one-body
radiation calculation plus a cross-section. The one-body radiation pattern
of a magnetic dipole is well-known analytically.

### OQ3: Time-Averaged θ Radial Profile (DC and RMS)

The DC component of the θ field around a single braid has been measured
as 0.2% of the oscillation amplitude (V34), but its RADIAL PROFILE has
not been measured.

**Important caveat**: Since the effective Coulomb force comes from AC
radiation scattering (OQ2), the DC component may be physically irrelevant
to the force law — it could be a small numerical artifact of the discrete
grid initialization. The AC RMS intensity profile ⟨θ_rms²(r)⟩ may be the
physically meaningful quantity for force calculations.

**Test**: Measure BOTH profiles from existing V34 or V41 SFA files:
1. DC profile: ⟨θ(r)⟩ (time average of signed θ at each radius)
2. RMS intensity profile: ⟨θ_rms²(r)⟩ (time average of θ² at each radius)

If ⟨θ_rms²⟩ ∝ 1/r² at long range, this confirms the radiation intensity
falls correctly for a 1/r² Coulomb force via momentum transfer. The DC
profile is a secondary diagnostic — if it's noisy or non-physical, the
AC radiation mechanism stands on its own.

### OQ4: Fine Structure Constant from Theory Parameters

Is α derivable from (η, m, μ, κ), or does it require an additional input?

**Test**: Compute the DC/total ratio of θ at r=10 (far field) for several
values of η = {0.1, 0.3, 0.5, 0.7, 1.0}. If DC/total ∝ η^n, the fine
structure constant α ∝ η^n, and η IS the fundamental EM coupling constant.

### OQ5: Photon Dispersion

Does a pure θ wave packet propagate at exactly c without dispersion?

**Test**: Initialize a localized transverse θ perturbation (no φ, no braid)
in an empty background. Track the wave packet centroid over T=100. Measure:
group velocity (should be c), dispersion (should be zero for massless field),
polarization (should be transverse, two modes).

### OQ6: Nonlinear Maxwell Corrections

At high field strengths (near a braid core), the V(P) nonlinearity modifies
the linearized Maxwell equations. What corrections appear?

**Test**: Solve the full nonlinear Cosserat equation perturbatively to second
order. The first-order corrections give standard Maxwell. The second-order
corrections give the nonlinear EM effects (analogous to Euler-Heisenberg
effective Lagrangian in QED, where photon-photon scattering occurs at high
intensity).

---

## Next Steps (Prioritized)

### Zero GPU cost (analytical):
1. **Full Maxwell derivation** with proper definitions, all four equations, and the ∇·E = 0 proof
2. **Ponderomotive force calculation** for two oscillating current loops
3. **Multipole expansion** of θ around a composite baryon from symmetry arguments

### Low GPU cost (post-processing existing data):
4. **⟨θ(r)⟩ time-averaged profile** from V34 or V41 SFA files
5. **Multipole decomposition** of θ from V41 UUD T=500 data

### Moderate GPU cost (new simulations):
6. **Pure θ wave packet propagation test** (N=128, T=100, minutes on CPU)
7. **η-variation of DC/total ratio** (5 single-braid runs)

### Future (V43+):
8. **Wire/Ohm's law simulation** (100+ braids, N=1024+)
9. **Photon-photon scattering** (two θ wave packets crossing)
