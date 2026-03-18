# Phase B Field Study: Comprehensive Findings

## The Fundamental Question
Can three real scalar fields with a triple-product potential produce gravity, EM, and mass?

## Answer: NO for gravity and EM. Partial for mass.

### Why No Gravity

The effective metric around a braid is determined by the Hessian of the potential:
M²_ab(x) = m²δ_ab + ∂²V/∂φ_a∂φ_b

For gravity, we need:
1. ALL eigenvalues > vacuum everywhere → attractive for all modes
2. 1/r decay → massless mediator
3. Universal coupling → same metric for all modes

What we get:
1. Mixed signs: center is attractive, shell at r≈4 is repulsive for one mode
2. Yukawa decay (range ~6 code units) → massive mediator
3. Three DIFFERENT eigenvalues → birefringent (NOT universal)

Root cause: The triple-product Hessian ∂²(P²/(1+κP²))/∂φ_a∂φ_b has off-diagonal
elements that can be negative (anti-correlated fields reduce the potential).
The mass term (m=1.5) ensures Yukawa decay. Both are structural features of the
model that cannot be tuned away.

### Why No EM

EM requires:
1. Phase-dependent scattering (charge)
2. Gauge invariance (conservation of charge)
3. Long-range 1/r force (massless gauge boson)

What we get:
1. No measurable wave-wave scattering even at high amplitude
   - The triple-product coupling needs persistent structured overlap (braid)
   - Plane waves don't overlap enough to scatter
2. No gauge symmetry — three real scalars have discrete symmetry S₃×(Z₂)³ only
3. No massless mode — all modes are massive (m²≥0.82)

Root cause: Real scalar fields cannot carry phase information in the way complex
fields can. The "torsion" diagnostic is a geometric effect of the braid structure,
not a field-theoretic gauge interaction.

### Why Mass Doesn't Emerge

The effective mass around a condensate varies from 0.82 (at v≈0.3) to 3.88
(at the braid center). The κ-regularization prevents the mass from ever reaching
zero or going negative. No tachyonic instability → no spontaneous symmetry
breaking → mass stays as a Lagrangian parameter.

### What DOES Work

1. **Stable soliton**: fc=0.96, survives perturbation, thermal equilibrium
2. **Mode conversion**: 3.7% channel mixing at braid center
3. **Tensor structure**: 41% anisotropy in effective metric (3 different eigenvalues)
4. **Localization**: Energy concentrated in tube with R≈3

## What Would Need to Change

### For Gravity
Option 1: Add a massless scalar field ρ that couples to the soliton energy density.
Then □ρ = -αT₀₀ gives 1/r potential. But this is just adding gravity by hand.

Option 2: Derivative coupling (BLV-type): if the kinetic term becomes
(1/2)g^μν(φ)∂_μφ·∂_νφ where g^μν depends on the background, then waves
propagate on an effective metric that IS gravity. But V6-V7 already explored
this (see MEMORY.md) and found only nuclear-scale (10^37× too strong) or
scalar (no tensor GW).

Option 3: Fundamentally different Lagrangian where the field equations
themselves contain the metric as a derived quantity. This is the Einstein-Cartan
approach — torsion IS the connection. Would need to reformulate with gauge fields.

### For EM
Need complex fields with U(1) gauge symmetry. Upgrade φ_a → ψ_a ∈ C with
L = |D_μψ_a|² - V(|ψ₀||ψ₁||ψ₂|) - (1/4)F_μν²
where D_μ = ∂_μ + ieA_μ and F_μν = ∂_μA_ν - ∂_νA_μ.
This adds EM by hand. The question is whether the braid structure is preserved.

### For Mass Emergence
Need a potential with spontaneous symmetry breaking: V has minimum at φ≠0.
The current V(P) = (μ/2)P²/(1+κP²) with μ<0 would be unstable without κ.
A Mexican-hat-like V(ρ) = λ(ρ²-v²)² with ρ = √(φ₀²+φ₁²+φ₂²) would give
a Higgs-like mechanism, but this doesn't preserve the triple-product topology.

## Additional Tests (T10F-H)

### T10F: Massless Mediator — IMPOSSIBLE
Making one field massless (m₀=0) destroys the braid (fc: 0.88→0.39).
The massless field radiates away freely. Far field decays as r^(-2.5), not 1/r.
**Fundamental tension: mass confines soliton but prevents 1/r gravity.**

### T10G: Complex Field Upgrade — VIABLE
Upgrading to 3 complex scalars preserves the braid (fc=0.91) and adds U(1) charges
(Q≈3-14). The soliton is a CHARGE-CARRYING structure. This is a necessary
ingredient for EM but not sufficient.

### T10H: Two Charged Braids — CHARGES DON'T MATTER
Same-charge, opposite-charge, and uncharged braids all attract identically
(ΔD≈-5 to -6). The force is purely Yukawa (massive field overlap), not Coulomb.
**Global U(1) symmetry does NOT produce EM-like forces. A gauge field A_μ is needed.**

### T10D: Kibble-Zurek — NO SPONTANEOUS BRAIDS
Random fields cooled with damping never form localized structures. The triple-product
coupling drives energy growth, not structure formation. Braids require specific
initialization — they are NOT natural topological defects.

## Comprehensive Conclusion

The Phase B field study has definitively answered the core questions:

| Question | Answer | Evidence |
|----------|--------|----------|
| Does the braid curve the field? | Partially | 4% speed reduction at center, but birefringent and sign-changing |
| Is there 1/r gravity? | NO | Yukawa range ~6, massive mediator, no massless mode possible |
| Is there EM? | NO | No wave-wave scattering, no charge-dependent force, no gauge symmetry |
| Does mass emerge? | NO | m=1.5 is a Lagrangian parameter; no tachyonic instability to drive SSB |
| Do braids form naturally? | NO | Random fields don't cool into braids |
| Does the complex upgrade work? | YES | Braid survives, carries charge, but charges don't matter |

### The Fundamental Architecture Problem

This model has **three structural gaps** that cannot be fixed by tuning parameters:

1. **No massless mode for gravity**: All propagating modes are massive (m²≥0.82).
   Making any field massless destroys the soliton. A separate massless field
   (graviton analog) must be ADDED to the model.

2. **No gauge symmetry for EM**: The triple-product potential has only global
   symmetries. Promoting U(1)_global → U(1)_local requires adding a gauge field A_μ,
   which is equivalent to putting EM in by hand.

3. **No spontaneous structure formation**: The V(P) potential with μ<0, κ>0
   doesn't break symmetry. It has a unique vacuum (φ=0) with no degenerate minima.
   Solitons are metastable energy lumps, not topological defects.

### Possible Paths Forward

**(a) Accept the model as nuclear physics only.** The braid resembles a nucleus
(localized, massive, finite binding energy). The Yukawa force between braids
is analogous to the nuclear force. This is consistent with the Skyrme model heritage.

**(b) Add gravity and EM explicitly.** Couple the braid field to a metric tensor
(gravity) and a gauge field (EM). This gives a Skyrme-Einstein-Maxwell system.
The braids would then be charged, gravitating solitons. But gravity and EM are
external inputs, not emergent.

**(c) Change the framework entirely.** Instead of scalar fields with a potential,
use a gauge-Higgs system where:
- The gauge field IS gravity+EM (spin-2 + spin-1 from same connection)
- The Higgs field breaks the symmetry to produce mass
- Topological defects (braids, skyrmions) are guaranteed by π₃ of the vacuum manifold
This is closer to the original CHPT vision (see MEMORY.md).

**(d) Hybrid: V6 density wave + braid.** V6 showed that a conserved density field
with a massless wave equation can produce 1/r potential. If the braid modulates
a separate density field ρ, and □ρ sources from the braid energy, this gives
a 1/r static potential. The density field is the graviton analog. But V6 found
this gives G_eff ~ 10⁴⁰ × G_N (coupling uncontrolled).

### Recommendation
Path (c) is the most promising for EMERGENT gravity+EM. Path (d) is the
quickest numerical test. Path (b) gives a working model but no new physics.
Path (a) is the most honest assessment of what the model actually does.
