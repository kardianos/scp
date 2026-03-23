# V39: Gravitational Self-Trapping — Alternate Route to Stable Particles

## Motivation

V38 showed the toroidal braid fragments (6→32 clusters). Every compact
structure tried so far eventually disperses because V(P) binding (capped
at |μ/(2κ)| = 0.413 per volume) is always less than gradient + mass energy.

BUT: V34 Track GB found that **inverse coupling** creates a SELF-REINFORCING
collapse mechanism:
```
low core mass → V(P) dominates → P grows → Σφ² grows → m_eff drops further
→ V(P) grows more → collapse until κ saturation caps V(P)
```
This produced 19,000% E_pot growth — the binding EXPLODED. The collapse
was arrested by κ saturation, creating a "frozen blob" that slowly drained.

**The idea**: What if we can harness this collapse mechanism to create a
self-sustaining structure? Not a frozen blob (which drains), but a
DYNAMIC collapsed structure — a braid-like traveling wave inside a
gravitationally self-trapped region. The collapse creates the well,
the wave dynamics maintain coherence.

## Core Question

**Is there a critical ρ (field density) above which the self-reinforcing
depletion creates a self-sustaining gravitational trap?**

In GR: Schwarzschild radius r_s = 2GM/c². Above this density, nothing escapes.
In our theory: is there an analogous critical density where V(P) binding +
inverse coupling creates an inescapable well?

## Investigation Variants

### Variant A: 1D Braid with Inverse Coupling
Simplest test. Take the working 1D oscillon from V24 and add inverse coupling.
- Equation: ∂²φ/∂t² = ∂²φ/∂x² - m_eff²(Σφ²)φ - V'(P)
- m_eff² = α/(1 + β×Σφ²) (inverse coupling)
- Sweep (α, β) to find the boundary between "stable braid" and "collapse"
- At the boundary: does a STABLE intermediate state exist?
- Fast: 1D, N=1024, can run millions of steps on CPU

### Variant B: 1D Critical Density Search
Start with a uniform field of density ρ and measure:
- Below ρ_crit: field disperses (gradient wins)
- Above ρ_crit: field collapses (inverse coupling wins)
- AT ρ_crit: what happens? Marginal stability? Oscillation?
- Map ρ_crit vs (α, β, μ, κ) parameter space
- This is the "black hole threshold" in 1D

### Variant C: 2D Radial Collapse
Cylindrically symmetric 2D (r, z):
- Initialize a dense blob of radius R with density ρ > ρ_crit
- Does it collapse into a stable ring? A point? An oscillon?
- The 2D geometry allows angular momentum (rotation could stabilize)
- N=512 radial × 512 axial, still fast on CPU

### Variant D: Hybrid Coupling — Current + Inverse
The key experiment. Two mass terms:
- m_eff² = m_constant² + m_inverse²(Σφ²)
- m_constant² = 2.25 (current, provides vacuum stability)
- m_inverse² = α/(1 + β×Σφ²) (inverse, provides self-trapping)
- At low field: m_eff² ≈ 2.25 + α (normal vacuum)
- At high field: m_eff² ≈ 2.25 + α/(1+β×large) → approaches 2.25 (core softens)
- The constant term prevents full collapse; inverse provides extra binding
- Test in 1D first, then 3D on GPU

### Variant E: 3D Toroidal Braid with Inverse Coupling
If Variant D finds a viable (α, β) in 1D, apply it to the toroidal braid:
- The torus provides the traveling wave (self-reconstruction)
- The inverse coupling provides extra binding at the core
- The κ saturation prevents singularity
- The constant m² prevents vacuum collapse
- N=256 on GPU

### Variant F: Density-Dependent κ
Instead of field-dependent mass, make the SATURATION parameter κ
field-dependent: κ_eff = κ₀ / (1 + γ×Σφ²).
- At low density: κ large → V(P) saturates quickly (weak binding)
- At high density: κ small → V(P) can go deeper (stronger binding)
- This directly increases the binding ceiling in dense regions
- The V_max = |μ/(2κ_eff)| becomes field-dependent
- At high enough density: V_max → ∞ (no ceiling!)
- This IS a black hole mechanism if γ is large enough

## Priority Order

1. **Variant A** (1D braid + inverse): fastest, establishes if the mechanism works
2. **Variant B** (1D critical density): maps the phase space
3. **Variant D** (hybrid coupling): the actual proposed physics
4. **Variant F** (density-dependent κ): potentially the most powerful — removes the ceiling
5. **Variant C** (2D radial): if 1D shows promise, test geometry effects
6. **Variant E** (3D torus + inverse): final GPU test if earlier variants work

## Theoretical Framework

### The depletion-binding feedback loop
In the current theory (constant m²):
```
More binding → deeper depletion → attracts more matter → BUT
matter arriving adds gradient energy → gradient > binding → structure disperses
```
The feedback is NEGATIVE at step 4. Gradient pressure wins.

With inverse coupling:
```
More binding → deeper depletion → LOWER m_eff at core → V(P) gets STRONGER
→ even more binding → deeper depletion → ...
→ collapse arrested by κ saturation
```
The feedback is POSITIVE. But κ caps it.

With density-dependent κ (Variant F):
```
More binding → deeper depletion → LOWER κ_eff → V(P) ceiling RISES
→ even more binding → even lower κ_eff → ceiling rises more → ...
→ true runaway collapse (no cap!)
→ only arrested by: (a) all field energy consumed, or (b) gradient pressure at very small scales
```
This IS a black hole mechanism. The question: does it stabilize into a particle, or
does it just collapse to a point and radiate everything?

### Estimating ρ_crit in 1D
For a 1D blob of width W and amplitude A:
- Gradient energy: E_grad ~ A²/W
- Mass energy: E_mass ~ m² × A² × W
- Binding energy: E_pot ~ |μ/(2κ)| × W (at saturation)

Stability: E_pot > E_grad + E_mass (binding exceeds dispersive forces)
|μ/(2κ)| × W > A²/W + m² × A² × W

At the optimal W: W_opt = A / (m × sqrt(|μ/(2κ)|))
This gives: ρ_crit ~ m × sqrt(|μ/(2κ)|) / A

With our parameters: ρ_crit ~ 1.5 × sqrt(0.413) / A ~ 0.96/A

For A=0.8 (braid amplitude): ρ_crit ~ 1.2

This is the ambient field density (~A² ~ 0.64 for the braid core), suggesting
we're already NEAR the critical density. The inverse coupling would push us
over the edge.

## Tools Needed

### 1. `v39/src/oscillon_1d.c` — 1D simulation with configurable mass coupling
- Standard 3-field equation in 1D (from V24)
- Support: constant m², inverse coupling, hybrid, density-dependent κ
- Fast: N=1024-4096, millions of steps per second on CPU
- Output: timeseries + field snapshots

### 2. `v39/src/critical_density_1d.c` — 1D critical density scanner
- Initialize uniform blob of varying density
- Sweep (ρ, α, β, γ) parameter space
- Binary search for ρ_crit at each parameter point
- Output: phase diagram

### 3. `v39/src/collapse_2d.c` — 2D radial collapse simulator
- Cylindrically symmetric (r, z)
- Inverse coupling + optional rotation
- Output: SFA-like snapshots for visualization

### 4. Modified `v38/src/seedrun_cuda.cu` — 3D with inverse coupling
- Add `-inv_alpha` and `-inv_beta` command line options
- m_eff² = MASS2 + inv_alpha/(1 + inv_beta*Sigma_phi2) in force kernel
- For Variant E (3D torus with inverse coupling)

## Success Criteria

A structure that:
- Uses the self-reinforcing depletion mechanism (not just constant m²)
- Is bounded (doesn't expand forever)
- Maintains V(P) binding for T ≥ 200 in absorbing BC
- Has internal dynamics (not a frozen blob)
- Ideally: 1 cluster throughout (sfa_frag verified)

## Budget
- Variants A, B, C, D: CPU only, minutes each
- Variant E: V100 GPU, ~30 min (~$0.07)
- Variant F: CPU first, GPU if promising
