# V27: From Braided Soliton to Emergent Physics

## Big Picture Goal

Demonstrate that three displacement fields φ₁, φ₂, φ₃ with a braided
topology, propagating through their own self-consistent metric, produce:

1. **Stable matter** — a non-dispersing, non-breathing soliton
2. **Gravity** — spin-2 (l=2) mediated attraction between solitons
3. **Electromagnetism** — spin-1 from the torsion of the frame field
4. **Mass** — from the dynamics, not from an imposed parameter
5. **Confinement** — topological protection against unwinding

This is the COMPLETE theory: matter, forces, and spacetime from one
Lagrangian with one field (three-component displacement).

## Where We Are (V26 Results)

| What | Best Result | Gap |
|------|------------|-----|
| Spin-2 content | 41.5% l=2 (static braid) | Not dominant; static braid unwinds |
| Braid survival | |P|=0.574 (DynA propagating) | fc only 0.27; needs both |P| AND fc |
| Localization | fc=0.71 (DynB rotating) | |P| dies; braid topology lost |
| Self-consistent metric | Stable g^{ij} (Phase 4) | Not tested on propagating braid |
| Causal gravity | c propagation (V24-L4) | Only tested on oscillons |
| Confinement | None | The central unsolved problem |

**The V26 discovery**: propagation maintains |P| (braid structure) while
rotation maintains fc (localization). DynAB (combining both) inherited
rotation's |P| death — the two mechanisms INTERFERE. V27 must find how
to get BOTH simultaneously.

## Progress Milestones

### Milestone 1: Stable Propagating Braid (|P| > 0.1 AND fc > 0.5)

**The problem**: DynA has |P|=0.574 but fc=0.27. DynB has fc=0.71 but
|P|=3e-5. DynAB gets neither. Why does combining them fail?

**Hypothesis**: Rotation mixes the azimuthal (θ) and axial (kz) phase
structures, disrupting the triple-product coupling that propagation
alone maintains. The fix: decouple the propagation direction from the
rotation plane.

**Tests**:
- M1a: DynA with PAIRWISE COUPLING λ_pw=0.5 (deepens the potential
  well, may increase fc without killing |P|)
- M1b: DynA with SMALLER domain (L=10 instead of 20, stronger confinement
  from periodic BC, less room for dispersal)
- M1c: DynB with AXIAL propagation added at LOW velocity (v=0.01, not
  0.30 — gentle propagation that doesn't disrupt rotation)
- M1d: Scan the propagation speed: v_g = {0.01, 0.05, 0.1, 0.2, 0.3}
  at fixed Ω=0.1. Find the optimal balance.

**Success**: fc > 0.5 AND |P| > 0.1 at t=500.

### Milestone 2: Spin-2 Gravity from the Propagating Braid

**The problem**: DynA has l=2=26.7% but the strain is a mix of l=0 (monopole)
and l=2 (quadrupole). For gravity, l=2 must DOMINATE at long range.

**Hypothesis**: The l=0 component comes from the radial envelope (Gaussian
tube profile). The l=2 comes from the three-strand azimuthal offset (2π/3).
At long range (r >> R_tube), the l=0 decays faster (r⁻¹ in 3D) than l=2
(also r⁻¹, but from a quadrupole source → r⁻³ for the potential). So l=2
NEVER dominates at long range for a single soliton.

For TWO solitons: the INTERACTION is what matters. The l=2 strain from
soliton A at soliton B's location determines the tidal force. If this
tidal force is ATTRACTIVE and QUADRUPOLAR: that's gravity.

**Tests**:
- M2a: Compute the strain field of the propagating braid (DynA) at
  distances r = 5, 10, 15, 20 from the center. Decompose into l=0,1,2
  at each r. Does l=2/l=0 INCREASE with distance?
- M2b: Place TWO propagating braids at D=30 (both DynA, same parameters).
  Measure the FORCE (centroid acceleration). Is it attractive?
- M2c: Measure the angular pattern of the force. Is it quadrupolar?
- M2d: Apply the self-consistent metric (g_{ij} = δ_{ij} + 2ε_{ij})
  to the propagating braid. Does the braid survive with metric backreaction?

**Success**: Two braids attract with measurable quadrupolar angular pattern.

### Milestone 3: Electromagnetic Structure from Torsion

**The problem**: The torsion ω_{ij} = ½(∂_iφ_j - ∂_jφ_i) is the EM
candidate. V26 Phase 4 showed that ADDING κ_T ω² term kills the braid.
But the torsion EXISTS as an observable — it just shouldn't be penalized.

**Approach**: Don't add an EM COUPLING. Instead, MEASURE the torsion field
of the propagating braid. Compute:
- Vorticity: Ω_k = ε_{ijk} ω_{ij}
- Torsion flux: Φ_T = ∫ Ω · dA through a cross-section
- Is the flux QUANTIZED (integer × 2π from braid crossing number)?
- Does the torsion field have the angular structure of a magnetic field?

**Tests**:
- M3a: Compute ω_{ij} at every grid point for the DynA propagating braid.
  Map the vorticity field Ω_k.
- M3b: Integrate Φ_T through cross-sections at z = 0, L/4, L/2, 3L/4.
  Is Φ_T constant (flux conservation)? Is it quantized?
- M3c: Two braids with SAME torsion flux → repel or attract? (like charges)
  Two braids with OPPOSITE flux → attract? (unlike charges)

**Success**: Quantized torsion flux, EM-like force between braids.

### Milestone 4: Mass from Dynamics

**The problem**: m=1 is imposed. V26 showed m=0 always collapses.
The mass should emerge from the braid's own dynamics, not from a parameter.

**Approach**: Instead of m²φ², use a term that creates an EFFECTIVE mass
from the braid's own field configuration:

    m²_eff(x) = β · |∇φ|² = β · (∂_iφ_a)(∂_iφ_a)

This is zero in vacuum (no gradients → massless propagation) and nonzero
inside the braid (large gradients → effective mass → confinement).

This is ALREADY in the Lagrangian — it's the gradient kinetic term!
The "mass" from ½(∂_iφ_a)² acts as a position-dependent mass gap:
large inside the braid (confinement), zero outside (free propagation).

The issue: the gradient term gives ω² = k² (massless dispersion), not
ω² = k² + m² (massive). For a MASS GAP, need the effective m² to appear
at k=0. This requires nonlinearity — the triple product V(P) provides it.

**Tests**:
- M4a: Measure the dispersion relation ω(k) of small perturbations
  INSIDE the propagating braid vs OUTSIDE.
  Inside: ω² = k² + m²_eff (if m_eff exists)
  Outside: ω² = k² + m² (bare mass)
- M4b: Gradually REDUCE m from 1.0 to 0 while maintaining propagation.
  At what m does |P| collapse? Is there a minimum m for braid survival?
- M4c: Try m=0 with STRONGER triple product (μ=-100, κ=100) — does
  the stronger binding replace the mass term?

**Success**: Identify the minimum mass for braid survival, or find a
parameter regime where m=0 braids survive through propagation + strong binding.

### Milestone 5: Topological Confinement

**The problem**: The braid unwinds because nothing PREVENTS it.
The crossing number is not conserved by the dynamics.

**Approach**: Find a Lagrangian term whose variation EXACTLY conserves
the braid crossing number. This is analogous to the Skyrme term L₄
that stabilizes the Skyrmion through Derrick's theorem.

For braids: the crossing number n is related to the LINKING NUMBER
of the three field tubes. The linking number is a topological invariant
computed from:

    Lk = (1/4π) ∮∮ (dr₁ × dr₂) · (r₁-r₂)/|r₁-r₂|³

This is a NONLOCAL quantity — it can't be written as a local Lagrangian
density. But it CAN be written as:

    Lk = ∫ A · B d³x

where A is the "gauge potential" of one tube and B is the "magnetic field"
of the other. For our three fields:

    A_i = ε_{ijk} φ_a ∂_j φ_b    (Chern-Simons-like)
    B_k = ε_{ijk} ∂_i φ_a ∂_j φ_b

A Chern-Simons term ∫ A · B d³x in the Lagrangian would EXACTLY
conserve the linking number. This is a TOPOLOGICAL term — it doesn't
contribute to the energy but constrains the dynamics.

**Tests**:
- M5a: Compute the linking number Lk of the propagating braid at t=0 and
  t=500. Does it change? (Without the CS term, it probably does.)
- M5b: Add a Chern-Simons-like term to the Lagrangian:
  L_CS = λ_CS · ε_{abc} ε_{ijk} φ_a (∂_i φ_b)(∂_j φ_c)
  This is a total derivative in some cases but may affect the dynamics.
- M5c: Does the CS term PREVENT unwinding (Lk conserved)?
  Does the braid survive LONGER with it?
- M5d: Does the CS term affect the mass gap or dispersion?

**Success**: Braid with conserved linking number that survives indefinitely.

## Implementation Strategy

### Phase I (Milestones 1 + 4): Stable Propagating Braid

Start from DynA (the |P|-maintaining configuration). Add pairwise coupling
to improve fc. Scan mass reduction. All in 3D, N=128, L=20.

Deliverable: a braid with |P| > 0.1, fc > 0.5, at the MINIMUM mass.

### Phase II (Milestone 2): Gravity Between Braids

Using the best braid from Phase I: two-braid interaction test.
Self-consistent metric (g_{ij} from strain). Force measurement.

Deliverable: measured inter-braid force with angular decomposition.

### Phase III (Milestone 3): EM Structure

Torsion mapping of the Phase I braid. Flux quantization test.
Two-braid EM-like interaction.

Deliverable: vorticity field map, flux measurement, charge-like interaction.

### Phase IV (Milestone 5): Topological Confinement

Add Chern-Simons term. Test linking number conservation.
Long-time stability with and without CS term.

Deliverable: indefinitely stable braid with conserved topology.

### Phase V (Integration): Full Theory

Combine all: stable braid (Phase I) with self-consistent metric (Phase II),
torsion EM (Phase III), and topological protection (Phase IV).

Deliverable: two interacting braided solitons with emergent gravity, EM,
and topological stability — the complete theory.

## The Lagrangian (V27 Target)

    L = Σ_a [½(∂_t φ_a)² - ½(∂_i φ_a)²]     kinetic + gradient
        - ½m²_min Σφ_a²                        minimal mass (from M4)
        - (μ/2)P²/(1+κP²)                      triple product (binding)
        - λ_pw(φ₁φ₂+φ₂φ₃+φ₃φ₁)               pairwise (force hierarchy)
        + λ_CS ε_{abc}ε_{ijk}φ_a(∂_iφ_b)(∂_jφ_c)  Chern-Simons (topology)

    Derived quantities (not separate fields):
        ε_{ij} = ½(∂_iφ_j + ∂_jφ_i)            strain → gravity (spin-2)
        ω_{ij} = ½(∂_iφ_j - ∂_jφ_i)            torsion → EM (spin-1)
        g_{ij} = δ_{ij} + 2ε_{ij}               metric → ds = c·dt
        Φ_T = ∫ε_{ijk}ω_{jk} dA                 torsion flux → charge

## Success Criteria for V27

| Milestone | Minimum | Target | Breakthrough |
|-----------|---------|--------|-------------|
| M1 (stable braid) | fc>0.4, |P|>0.05 | fc>0.5, |P|>0.1 | fc>0.7, |P|>0.3 |
| M2 (gravity) | Two-braid attraction | Quadrupolar pattern | Universal coupling |
| M3 (EM) | Nonzero Φ_T | Quantized flux | Charge-like force |
| M4 (mass) | m_min identified | m_min < 0.5 | m=0 braid survives |
| M5 (confinement) | Lk measured | CS slows unwinding | Lk exactly conserved |

## Parameters

Base: μ=-20, κ=20, m=1.0 (from DynA)
Pairwise: λ_pw = 0.5 (from V24-ME)
CS: λ_CS to be determined
Grid: N=128, L=20 (single braid), N=128, L=40 (two braids)

## Reference Code

- v26/dynA_propagate/src/dynA.c (propagating braid — starting point)
- v26/dynB_rotate/src/dynB.c (rotating braid)
- v26/src/v26_phase4.c (torsion + strain analysis)
- v25/src/v25.c (self-consistent metric)
- v24/fundamental/combo_1257/src/combo1257.c (self-gravitating lattice)
