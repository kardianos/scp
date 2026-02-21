# V6 Foundation: First-Principles Field Theory with Torsion and Density Conservation

## Why V6?

V2 (Skyrme model) proved that:
- Topological solitons exist and are stable (B=1-4 verified to <0.1%)
- The SU(2) frame field twists non-locally (frame deviation ~1/r², torsion ~1/r³)
- The density mode is MASSIVE (m = √(2λ)ρ₀) → all density perturbations are Yukawa
- No conserved density law exists in the Lagrangian
- The soliton DEPLETES density at its core (ρ(0) = 0.97ρ₀ at λ=100), opposite of narrative
- The only 1/r channel is Path 3 (constraint field p sourced by B⁰), but g_top is free

The CHPT narrative spec says gravity comes from density conservation + depletion.
The v2 math doesn't implement this. V6 starts from scratch to BUILD a theory that does.

**Key principle**: Don't start from a field space and discover what forces emerge.
Start from the REQUIREMENTS and derive the simplest Lagrangian that satisfies all of them.

---

## Part 1: Requirements Matrix

### Physical Requirements

| # | Requirement | Source | What it means |
|---|-----------|--------|---------------|
| R1 | Continuous medium | CHPT Axiom 1 | Field ψ(x,t) defined at every point |
| R2 | Conserved total content | CHPT Axiom 3 | ∫ρ(x,t) d³x = N = const for all t |
| R3 | Directional structure | Torsion hypothesis | Field has orientation, not just amplitude |
| R4 | Stable localized structures | Particles exist | Finite-energy, topologically protected |
| R5 | Long-range 1/r interaction | Gravity | Not Yukawa, not 1/r³ |
| R6 | Attractive force | Gravity is attractive | Depletion → back-pressure toward structure |
| R7 | Universal coupling | Equivalence principle | Every structure gravitates ∝ its energy |
| R8 | Coupling determined by theory | Fix g_top problem | No free gravity coupling constant |

### Mathematical Implementation Options

Each requirement maps to one or more mathematical structures. The table below lists
ALL options considered, not just the one selected.

| Req | Option A | Option B | Option C |
|-----|---------|---------|---------|
| R1 | Scalar ρ(x) | Multi-component ψ ∈ R^n | Geometric (frame field e_a) |
| R2 | Continuity eqn ∂_t ρ + ∇·j = 0 | Global constraint (Lagrange mult.) | Incompressibility ρ=const |
| R3 | SU(2) rotor q(x) | SO(3) frame e_a(x) | Cl⁺(3,0,1) multivector |
| R4 | π₃(S³)=Z (Skyrmion) | Kelvin vortex (circulation) | Q-ball (Noether charge) |
| R5 | Massless scalar (no potential) | Massless spin-2 (emergent metric) | Constraint propagator (Path 3) |
| R6 | Depletion + pressure | Twist gradient force | Effective metric (BLV) |
| R7 | E_twist ∝ ρ (linear coupling) | Stress-energy source (like GR) | B⁰ ∝ mass (topological) |
| R8 | G = f(c_s, ρ₀) determined | Renormalization group | Topological (WZW) |

### Compatibility Analysis

Not all combinations work. Key constraints:

1. **R2 + R5**: If ρ is conserved AND the density mode is massless, then static density
   perturbations satisfy ∇²δρ = source → δρ ~ 1/r. This is the ONLY natural way to get 1/r
   from a scalar field without fine-tuning. **Requires**: no potential V(ρ).

2. **R4 + R5**: Topological stability (R4) typically requires a mass gap in the field
   (Derrick's theorem). But R5 requires a massless mode. Resolution: the TWIST field can
   be massive (localizing the structure) while the DENSITY field is massless (carrying the
   long-range force). Two sectors with different mass gaps.

3. **R6 + R7**: Attractiveness (R6) requires that structures deplete density. Universality
   (R7) requires that depletion ∝ energy. Both are satisfied if E_twist ∝ ρ·|∇q|², because
   then δρ ∝ -|∇q|² and the total deficit ΔQ ∝ E_twist ∝ mass.

4. **R3 + R4**: Directional structure (R3) on S³ guarantees π₃ ≠ 0, giving topological
   particles (R4). Other target spaces: S² gives π₃(S²)=Z (baby Skyrmions, 2D only),
   SU(N) gives π₃=Z for all N, Grassmannians etc.

5. **R2 + incompressibility**: If ρ = const (incompressible, Option C), then trivially
   conserved but NO depletion is possible → R5, R6 fail. This is exactly the σ-model
   limit of v2 — ruled out.

### Eliminated Options

| Option | Why eliminated |
|--------|---------------|
| R2-C: Incompressibility | No depletion → no gravity (v2 σ-model) |
| R5-B: Spin-2 emergent metric | Requires specific matter content (not first-principles) |
| R5-C: Constraint propagator | g_top is free (v2 Path 3-5 result) |
| R4-C: Q-ball | Requires U(1) symmetry, not topological, decays if perturbed |
| R8-B: Renormalization group | Requires quantum theory, not classical |
| R8-C: Topological (WZW) | v2 showed Cl⁺(3,0,1) is not semisimple, WZW obstructed |

### Selected Combination

| Req | Selection | Reason |
|-----|-----------|--------|
| R1 | (ρ, q): density + SU(2) rotor | Minimal: one scalar + one group element |
| R2 | Continuity equation | Physical: field content flows, not created/destroyed |
| R3 | SU(2) rotor q(x) | Simplest with π₃ ≠ 0 |
| R4 | π₃(S³) = Z | Proven in v2 (Skyrmions) |
| R5 | No potential for ρ | Makes density mode massless → 1/r |
| R6 | E_twist ∝ ρ | Depletion follows from energy minimization |
| R7 | Linear ρ-twist coupling | ΔQ ∝ E_twist ∝ mass automatically |
| R8 | G = E_twist/(4πc_s²ρ₀²) | Determined by sound speed and background density |

---

## Part 2: First-Principles Derivation

### Step 1: State variables

The field at each point (x,t) is described by:
- **ρ(x,t) ≥ 0**: density (how much field content is present)
- **q(x,t) ∈ S³ ≅ SU(2)**: orientation (which way the field points)

These are INDEPENDENT. The density can vary without the orientation changing, and vice
versa. This is the key departure from v2, where ψ = ρ·q bundled them together.

Total degrees of freedom per point: 1 (density) + 3 (orientation on S³) = 4.

### Step 2: Conservation law

The density satisfies a continuity equation:

    ∂ρ/∂t + ∇·(ρv) = 0                                              ... (1)

where v(x,t) is the field velocity. This ensures:

    d/dt ∫ρ d³x = -∮ρv · dA = 0    (for fields decaying at infinity)

Total field content N = ∫ρ d³x is conserved. This is R2.

Note: v is NOT an independent variable — it is determined by the field dynamics
(derived from the Lagrangian via Noether's theorem or the Euler equation).

### Step 3: Energy functional

The energy has three contributions, each derived from a physical requirement.

**E₂: Twist energy (quadratic in derivatives)**

    E₂ = ½ ∫ ρ(x) |ω(x)|² d³x                                      ... (2)

where ω_i = (∂_i q) q⁻¹ ∈ su(2) is the frame torsion (Maurer-Cartan form pulled
back to physical space). |ω|² = Σ_i tr(ω_i†ω_i).

Physical meaning: twisting the orientation costs energy proportional to the local
density. Where there is more field, twisting it costs more. Where there is less field,
twisting is cheaper.

**This is the KEY coupling**: E₂ is LINEAR in ρ. This means:
- Where |ω|² is large (near a knot), the system saves energy by reducing ρ
- Where |ω|² is small (far from knots), ρ can be larger
- Conservation (R2) forces ρ to redistribute: low near knots, high far away
- This redistribution IS the density depletion that produces gravity

**E₄: Stabilization energy (quartic in derivatives)**

    E₄ = ½c₄ ∫ |ω ∧ ω|² d³x                                        ... (3)

where (ω∧ω)_{ij} = [ω_i, ω_j] and |ω∧ω|² = Σ_{i<j} tr([ω_i,ω_j]²).

This is the Skyrme term. It is ρ-INDEPENDENT — it stabilizes the twist
structure against collapse (Derrick scaling) without affecting the density dynamics.

Physical meaning: there is an additional energy cost for regions where the
orientation twists in multiple directions simultaneously.

**E_p: Pressure energy (equation of state)**

    E_p = α ∫ ρ² d³x                                                 ... (4)

This is the simplest term that provides a restoring force for density perturbations.
It penalizes density concentration: the system prefers uniform ρ over lumpy ρ.

Physical meaning: the field has internal pressure p = 2αρ. Compressed regions
push outward. This is what makes the depletion zone push knots toward each other.

The sound speed: c_s² = dp/dρ = 2α. This is a FREE PARAMETER of the theory
(set by α), but it is the ONLY free parameter — the gravity coupling G is then
determined (see Step 7).

**What is NOT included**:
- No potential V(ρ) = λ(ρ-ρ₀)⁴ or similar. No preferred density set by a potential.
  The background ρ₀ = N/V is set by total content N and volume V.
- No kinetic term for ρ yet (added in Step 4 for dynamics).

### Step 4: Lagrangian

The full Lagrangian density:

    L = ½ρ|q̇ q⁻¹|²/c² + ½ρ̇²/(2βρ) - ½ρ|ω|² - ½c₄|ω∧ω|² - αρ²   ... (5)

Term by term:
- ½ρ|q̇q⁻¹|²/c²: kinetic energy of orientation (twist waves)
- ½ρ̇²/(2βρ): kinetic energy of density fluctuations (sound waves)
- -½ρ|ω|²: twist potential energy (E₂ with ρ coupling)
- -½c₄|ω∧ω|²: Skyrme stabilization (E₄)
- -αρ²: pressure (E_p)

The density kinetic term ½ρ̇²/(2βρ) is the standard form for a compressible medium.
It gives sound waves with speed c_s² = 2αρ₀/β (for perturbations around ρ₀).

For a non-relativistic treatment: set c→∞ in the first term, keep the rest.

### Step 5: Static equilibrium

For a static configuration (q̇ = 0, ρ̇ = 0), minimize:

    E[ρ,q] = ∫[½ρ|ω|² + ½c₄|ω∧ω|² + αρ²] d³x                     ... (6)

subject to ∫ρ d³x = N.

**Variation with respect to ρ** (Lagrange multiplier μ):

    δE/δρ = ½|ω(x)|² + 2αρ(x) = μ                                   ... (7)

Solving for ρ:

    ρ(x) = (μ - ½|ω(x)|²) / (2α) = ρ₀ - |ω(x)|²/(4α)             ... (8)

where ρ₀ = μ/(2α) is the far-field density (where |ω| = 0).

**This is the depletion profile**: ρ is REDUCED where the twist |ω|² is large.
The density deficit at position x is:

    δρ(x) = ρ(x) - ρ₀ = -|ω(x)|²/(4α)                              ... (9)

The total deficit:

    ΔQ = ∫(ρ₀ - ρ) d³x = (1/4α) ∫|ω|² d³x = E₂/(2αρ₀)            ... (10)

(using E₂ = ½∫ρ|ω|² ≈ ½ρ₀∫|ω|² for small depletion δρ ≪ ρ₀).

**The deficit is proportional to the twist energy, which is proportional to the mass.**
This is R7 (universality): gravitational charge ∝ mass.

**Variation with respect to q** (orientation field):

This gives the Euler-Lagrange equation for the twist field q(x), with ρ(x) entering
as a position-dependent coefficient. For a B=1 hedgehog q = cos(f/2) + sin(f/2)(r̂·σ):

    d/dr[ρ(r)r²f'] + [ρ(r) - c₄(f'²+sin²f/r²)]sin(2f)/r² ... = 0  ... (11)

This is the standard Skyrme radial ODE but with ρ(r) replacing ρ₀. The twist field
sees the depleted density profile and adjusts.

### Step 6: Density dynamics and the Poisson equation

For small perturbations around a static knot, the density satisfies:

    ∂²δρ/∂t² = c_s² ∇²δρ + S(x,t)                                   ... (12)

where S is the source term from the twist-density coupling. In static equilibrium
(∂²/∂t² = 0):

    c_s² ∇²δρ = -S(x)                                                ... (13)

This is the POISSON EQUATION for δρ, with the twist energy density as source.

For a localized knot (S(x) concentrated within radius R_core), the far-field solution:

    δρ(r) = -S_total / (4πc_s² r)    for r ≫ R_core                  ... (14)

where S_total = ∫S d³x is the total source strength, proportional to the knot's
twist energy E₂.

**This is the 1/r depletion that the narrative spec describes.** It arises naturally
from the Poisson equation for a conserved, pressurized density field with no mass gap.

### Step 7: The gravity coupling

The force on a test knot (with twist energy E₂_test) in the density gradient of a
source knot (with twist energy E₂_source):

    F = -(∂E_test/∂ρ) · ∇δρ                                          ... (15)

Since E_test = ½∫ρ|ω_test|² ≈ (E₂_test/ρ₀)·ρ, we have ∂E_test/∂ρ ≈ E₂_test/ρ₀.

The density gradient from (14): ∇δρ = S_total/(4πc_s²r²) r̂.

Combining:

    F = -(E₂_test/ρ₀) · S_total/(4πc_s²r²) r̂                       ... (16)

With S_total ∝ E₂_source and E₂ ∝ mass M (virial: E₂ = E₄ for Skyrmions, so
M = 2E₂), this becomes:

    F = -G_eff · M_source · M_test / r²                              ... (17)

with the effective Newton's constant:

    G_eff = κ / (4π c_s² ρ₀²)                                        ... (18)

where κ is a dimensionless constant of order 1, determined by the ratio
S_total/E₂_source (how efficiently twist energy sources the density perturbation).

**G_eff IS DETERMINED BY THE THEORY.** It depends on:
- c_s: sound speed (from the pressure coefficient α)
- ρ₀: background density (from total field content N/V)
- κ: numerical constant (from the knot's radial profile)

There is no free coupling g_top. The gravity strength is fixed by the equation of
state and the background density. This solves the v2 Problem (R8).

---

## Part 3: Comparison with V2

| Feature | V2 (Skyrme) | V6 (Conserved Torsion) |
|---------|------------|----------------------|
| State | ψ = ρ·q (bundled) | (ρ, q) independent |
| ρ conservation | No (ρ determined by V(ρ)) | Yes (continuity equation) |
| Density potential | V = λ(ρ²-ρ₀²)² (massive) | None (ρ₀ = N/V, massless) |
| Background ρ₀ | Set by potential minimum | Set by total content N |
| Density mode mass | m = √(2λ)ρ₀ (Yukawa) | m = 0 (1/r) |
| Depletion | ρ(0)<ρ₀ but Yukawa recovery | ρ(0)<ρ₀ with 1/r recovery |
| Gravity range | Exponential (Yukawa) | 1/r (Poisson) |
| Gravity coupling | g_top FREE | G = κ/(4πc_s²ρ₀²) DETERMINED |
| Twist field | Same (SU(2) rotor, π₃=Z) | Same (SU(2) rotor, π₃=Z) |
| Skyrme term | Same (E₄ for stability) | Same (E₄, ρ-independent) |
| Topological charge | Same (B ∈ Z) | Same (B ∈ Z) |
| Soliton profile f(r) | Same ODE (with ρ correction) | Same ODE (with ρ correction) |

### What changes
1. ρ is promoted from auxiliary variable to CONSERVED dynamical field
2. The Mexican hat potential is REMOVED — density has no preferred value from a potential
3. A pressure term αρ² replaces V(ρ), giving sound waves instead of Yukawa modes
4. The gravity mechanism becomes the Poisson equation for δρ, not a constraint field

### What is preserved
1. The SU(2) twist field q(x) and all topological properties (B=1-4 solitons)
2. The Skyrme term E₄ (stability against collapse)
3. The energy functional structure (E₂ + E₄ + ...)
4. The radial profile ODE (with minor ρ-dependent correction)
5. All v2 nuclear physics results (scattering, binding, normal modes)

---

## Part 4: Feature/Concept Matrix

This maps every physical concept in the narrative spec to its v6 mathematical term.

### Fundamental Objects

| Narrative Concept | Mathematical Term | Equation |
|------------------|------------------|----------|
| "The field" | Pair (ρ, q) at each point | ρ: R³→R⁺, q: R³→S³ |
| "Field density" | ρ(x,t) | Scalar, conserved |
| "Field orientation" | q(x,t) ∈ SU(2) | 3 angles (e.g. Euler or axis-angle) |
| "Vacuum / empty space" | ρ = ρ₀, q = 1 (uniform) | E[ρ₀,1] = αN²/V (only pressure) |
| "Total field content" | N = ∫ρ d³x | Conserved exactly |
| "Background density" | ρ₀ = N/V | NOT from a potential |

### Localized Structures

| Narrative Concept | Mathematical Term | Equation |
|------------------|------------------|----------|
| "Knot" / "particle" | Topological defect in q | B = (1/2π²)∫q*(vol_S³) ∈ Z |
| "Knot concentrates density" | WRONG — knot DEPLETES density | ρ(core) < ρ₀ (Eq 9) |
| "Knot borrows from surroundings" | Density redistributes globally | ΔQ = E₂/(2αρ₀) (Eq 10) |
| "Mass = total energy" | M = E₂ + E₄ | Virial: E₂ = E₄, so M = 2E₂ |
| "Topological stability" | π₃(S³) = Z | B conserved, unwinding costs ∞ energy |

### Forces and Interactions

| Narrative Concept | Mathematical Term | Equation |
|------------------|------------------|----------|
| "Depletion zone" | δρ(r) = -S_total/(4πc_s²r) | 1/r from Poisson (Eq 14) |
| "Back-pressure" | F = -(∂E/∂ρ)∇δρ | Density gradient → force (Eq 16) |
| "Gravity = pressure differential" | F = -G_eff·M₁M₂/r² | Newton (Eq 17) |
| "Gravity is weakest" | G_eff ≪ 1 | Because c_s²ρ₀² ≫ 1 |
| "Gravity is universal" | ΔQ ∝ M | Linear ρ-twist coupling (Eq 10) |
| "Gravity is attractive" | ρ(core) < ρ₀ always | Twist always reduces ρ (Eq 9) |
| "Gravitational redshift" | Null-rotor in density gradient | v(x) depends on ρ(x) |
| "Gravitational lensing" | Refractive index n(x) = c/v(x) | BLV metric from density |

### Dynamics

| Narrative Concept | Mathematical Term | Equation |
|------------------|------------------|----------|
| "Wave propagation" | Twist waves in q field | ∂²q/∂t² = c²∇²q + nonlinear |
| "Sound waves" (gravity waves) | Density perturbations | ∂²δρ/∂t² = c_s²∇²δρ |
| "Speed of gravity" | c_s (sound speed) | c_s² = 2α (from pressure) |
| "Speed of light" | c (twist wave speed) | From ρ|q̇q⁻¹|²/c² kinetic term |
| "Gravitational waves" | Propagating δρ oscillations | Scalar (spin-0), speed c_s |

---

## Part 5: Open Questions and Limitations

### Known limitation: scalar gravity

This theory produces SCALAR gravity (spin-0 density field), not TENSOR gravity
(spin-2 metric perturbations). Consequences:

| Observable | GR (spin-2) | V6 (spin-0) | Match? |
|-----------|------------|-------------|--------|
| Newton 1/r² | Yes | Yes | YES |
| Gravitational redshift | δν/ν = Φ/c² | δν/ν = Φ/c² | YES |
| Light bending | θ = 4GM/(c²b) | θ = 2GM/(c²b) | HALF (factor 2 off) |
| GW polarization | h+, h× (tensor) | δρ (scalar) | NO |
| Frame-dragging | Yes (Lense-Thirring) | Partial (ρv term) | UNCLEAR |
| Perihelion precession | 6πGM/(c²a) per orbit | 3πGM/(c²a) (scalar) | HALF |

The factor-of-2 problem is generic to ALL scalar gravity theories. To get full GR,
need either:
(a) Spin-2 from the frame field q (not just density ρ)
(b) Additional vector modes from density current ρv
(c) Tensor modes from the full (ρ, q) system

This is a KNOWN limitation, not a fatal flaw. The 1/r and universality are the
hard parts; the spin-2 enhancement may come from the full non-linear theory.

### Does this still need a soliton?

The localized structure does NOT have to be a hedgehog Skyrmion. Any stable
configuration of q(x) with B ≠ 0 works. The key requirements for the structure are:

1. Localized: |ω|² → 0 as r → ∞ (finite energy)
2. Topological: B ≠ 0 (stability)
3. Finite energy: E₂ + E₄ < ∞

The specific SHAPE is irrelevant for the gravity mechanism — only the TOTAL
depletion ΔQ = E₂/(2αρ₀) matters at long range.

Alternatives to the hedgehog:
- Rational map configurations (B ≥ 2)
- Deformed Skyrmions in external density gradient
- Dynamically evolving (non-static) knots
- Knot clusters (nuclei)

The theory does NOT require specifying the form of the localized structure a priori.
Any topological defect in q(x) automatically depletes the density field and gravitates.

### The speed hierarchy question

For the theory to be physical:
- c (twist wave speed) should equal the speed of light
- c_s (sound/gravity speed) must satisfy c_s ≤ c (causality)

If c_s = c: gravity propagates at light speed (matches GR observation).
If c_s < c: gravity is slower than light (testable, but c_s close to c is needed).
If c_s ≪ c: gravity is quasi-instantaneous at nuclear scales (fine) but
  retardation effects differ from GR at cosmological scales.

The ratio c_s/c is determined by α and β in the Lagrangian. This is a constraint
on the parameters, not a problem.

### Derrick's theorem concern

In v2, the potential V(ρ) helps stabilize the soliton (Derrick scaling). Without V(ρ),
does the soliton still exist?

Yes, because:
1. The Skyrme term E₄ still prevents collapse (same as v2 σ-model)
2. The density ρ adjusts to the twist but doesn't destabilize it
3. The σ-model (ρ = ρ₀ = const) already has stable Skyrmions with E₂ + E₄ only
4. The density variation δρ/ρ₀ ~ E₂/(2αρ₀²V) is TINY for large V → negligible
   correction to the soliton profile

In fact, the v2 σ-model IS the large-α limit of v6: when α is very large, the
pressure strongly resists density variations, ρ → ρ₀ everywhere, and we recover
the standard Skyrme model. The v6 theory interpolates between:
- α → ∞: v2 σ-model (no depletion, no gravity)
- α finite: depletion + gravity, coupling G ~ 1/(αρ₀²)
- α → 0: strong depletion, large G (but soliton may become unstable)

---

## Part 6: Numerical Tests (Priority Order)

### Test 1: Static depletion profile

Solve the coupled system (Eq 8 + Eq 11) self-consistently:
- Initialize q from v2 σ-model profile (f(r) with B=1)
- Compute |ω|² = f'² + 2sin²f/r² from the profile
- Set ρ(r) = ρ₀ - |ω(r)|²/(4α) (Eq 8)
- Re-solve the profile f(r) with ρ(r) as coefficient
- Iterate until convergence

Output: ρ(r) vs r, δρ(r) vs r, ΔQ vs α, check that ΔQ ∝ E₂.

### Test 2: Long-range density tail

From Test 1, verify that δρ(r) → -C/r for large r (not exponential).
This is the CRITICAL test: if δρ is 1/r, gravity works. If Yukawa, it doesn't.

Method: solve the Poisson equation ∇²δρ = source with the twist energy as source.
Check the far-field solution against -ΔQ/(4πr).

### Test 3: Two-body force

Place two static knots at separation d. Compute total energy E(d).
The gravitational force is F = -dE/dd.

Check: F ∝ 1/d² for d ≫ R_core.
Check: F is attractive (E decreases as d decreases).
Check: F-coefficient matches G_eff from Eq 18.

### Test 4: Sound waves

Linearize around ρ = ρ₀, q = 1. Verify:
- Density perturbations propagate at c_s = √(2α)
- Twist perturbations propagate at c (from ρ₀|q̇q⁻¹|²/c²)
- The two sectors decouple at linear order (massive twist + massless density)

### Test 5: Physical parameter extraction

From the G_eff formula (Eq 18) and the v2 parameter fitting:
- ρ₀ and c are set by nucleon mass and radius (from v2)
- G_eff = G_Newton gives c_s (or equivalently α)
- Check: is c_s/c reasonable? (Must satisfy 0 < c_s ≤ c)

If c_s/c is of order 1, the theory is self-consistent.
If c_s/c ≪ 1 or c_s/c > 1, there's a problem.

### Test 6: Scalar gravity predictions

With G_eff calibrated to Newton:
- Compute predicted light bending angle (should be ½ of GR = 0.875")
- Compute predicted perihelion precession (should be ½ of GR)
- Compute sound wave (gravity wave) speed and polarization
- Compare with observations → identify what's missing for full GR

---

## Part 7: Implementation Plan

### Phase 1: Profile solver with self-consistent density

Modify the v2 radial solver to treat ρ(r) as determined by the twist profile,
not by a potential. Self-consistent iteration:

    f → |ω|² → ρ(r) from Eq 8 → f (new ODE with ρ(r)) → repeat

This is similar to the v2 finite-λ solver but SIMPLER: no BVP for ρ, just an
algebraic relation (Eq 8). The ρ equation is NOT a differential equation — it's
an explicit formula.

### Phase 2: Poisson solver for long-range tail

Solve ∇²Φ = -4πδρ_core(r) where δρ_core is the core depletion from Phase 1.
Verify Φ(r) → -ΔQ/(4πr) at large r.

This can be done analytically for spherical symmetry:
    Φ(r) = -(1/r)∫₀ʳ δρ(r')4πr'² dr' - ∫ᵣ^∞ δρ(r')4πr' dr'

### Phase 3: Two-body energy and force

Product ansatz (two hedgehogs at ±d/2), compute total energy as function of d.
Extract force, verify 1/r² and attractive.

### Phase 4: Parameter calibration

Use M_proton, r_proton, G_Newton to determine (ρ₀, c, c_s).
Check consistency.

---

## Appendix A: Derivation of the Poisson Equation for δρ

Starting from the static equilibrium condition (Eq 7):

    ½|ω(x)|² + 2αρ(x) = μ

This gives ρ(x) locally. But this is the NEAR-FIELD solution (valid where |ω|²
is non-negligible). For the FAR-FIELD (r ≫ R_core where |ω| → 0):

    ρ(x) → ρ₀ - correction

The correction comes from the global conservation constraint. Taking the Laplacian
of Eq 7:

    ½∇²|ω|² + 2α∇²ρ = 0     (since μ = const)

Therefore:

    ∇²δρ = -∇²|ω|²/(4α) = -S(x)/(4α)

where S(x) = ∇²|ω|² is nonzero only where the twist varies (near the core).

For a localized source, the Green's function gives:

    δρ(r) = (1/4α)·(1/4π)∫[∇'²|ω(x')|²/|x-x'|] d³x'

which for r ≫ R_core reduces to:

    δρ(r) = -(1/4α)·Q_eff/(4πr)

where Q_eff = ∫∇²|ω|² d³x = ∮∇|ω|²·dA = 0 at infinity.

Wait — this gives Q_eff = 0 by Gauss's theorem! The integral of ∇²|ω|² over
all space vanishes because |ω|² → 0 at infinity.

**Resolution**: The Poisson equation comes from the DYNAMICS, not just the static
equilibrium. The proper derivation:

In the fluid picture, static equilibrium with a source:

    -∇p + F_body = 0
    -c_s²∇ρ + F_twist = 0
    c_s²∇²ρ = ∇·F_twist

where F_twist = -∂(½ρ|ω|²)/∂x = -½ρ∇|ω|² - ½|ω|²∇ρ.

For small perturbations (δρ ≪ ρ₀):

    c_s²∇²δρ = ∇·F_twist ≈ -½ρ₀∇²|ω|² ≠ 0 locally

The source is -½ρ₀∇²|ω|², which IS localized (nonzero only near the core).
But again ∫∇²|ω|² d³x = 0 by Gauss.

**The resolution is that the force includes the |ω|²∇ρ cross-term.** The full
static equation is:

    2αρ + ½|ω|² = μ    (exact, no approximation)

Globally: ∫ρ = N fixes μ = 2αN/V + ½⟨|ω|²⟩ where ⟨·⟩ is volume average.
So μ = 2αρ₀ + O(1/V), and:

    ρ(x) = ρ₀ - |ω(x)|²/(4α) + O(1/V)

The far-field behavior of δρ(x) = -|ω(x)|²/(4α) follows the DECAY of |ω|².
For a Skyrmion: |ω|² ~ 4/r⁴ (massless) or ~e^{-2mr}/r² (massive pion).

**Key insight**: If the twist field ω decays as power law (massless pions), then
δρ ~ 1/r⁴. If massive pions, δρ ~ e^{-2mr}/r². Neither is 1/r.

The 1/r behavior requires a SEPARATE density equation that is sourced by the
localized twist, not just the algebraic relation ρ = ρ₀ - |ω|²/(4α).

**The correct approach**: The density field has its OWN dynamics, including
propagation. The static density equation is:

    c_s²∇²δρ = S_eff(x)                                              ... (A1)

where S_eff is an EFFECTIVE source that combines the twist energy localization
with the conservation constraint. S_eff is nonzero ONLY in the core (where the
twist changes rapidly), and ∫S_eff d³x ≠ 0 in general.

The source S_eff arises from the nonlinear coupling. Writing ρ = ρ₀ + δρ and
expanding the equilibrium condition to account for the conservation constraint
properly (via a Lagrange multiplier that enforces ∫δρ = 0):

    S_eff(x) = -½ρ₀∇²|ω|² + ½∇·(|ω|²∇δρ) + ...

The leading term ∫(-½ρ₀∇²|ω|²) d³x = 0 (Gauss), but the NEXT term is nonzero.

**Alternatively**: simply note that the total deficit ΔQ = ∫(-δρ) d³x = E₂/(2αρ₀)
is nonzero and positive (the soliton removes density from the vacuum). By
conservation, this deficit must come from somewhere — it means the density at
infinity is SLIGHTLY REDUCED from the nominal ρ₀. The 1/r tail comes from
treating the deficit as a point source in the Poisson equation:

    δρ(r) = -ΔQ/(4πr)    for r ≫ R_core                             ... (A2)

This is valid as an OUTER solution matched to the core depletion. The matching
works because the total removed charge ΔQ acts like a point sink at the origin.

The physical picture: the soliton acts as a density SINK — it has consumed ΔQ
units of density from the surrounding field. Conservation means this density is
missing from everywhere else. In an infinite volume, the per-point reduction
is infinitesimal, but the integrated effect gives a 1/r potential.

This is exactly the electrostatic analogy: a point charge creates a 1/r
potential even though ∇²(1/r) = 0 everywhere except the origin. The charge
IS the boundary condition, not a distributed source.

---

## Appendix B: Parameter Estimates

From v2 parameter fitting:
- ρ₀ = 1 code unit, e = 1 code unit
- 1 code energy = 9.098 MeV, 1 code length = 0.5624 fm
- M_proton = 103.13 code energy = 938 MeV
- R_proton = 1.496 code length = 0.84 fm

For v6, the gravity coupling:

    G_eff = κ/(4πc_s²ρ₀²)

In physical units: G_Newton = 6.674×10⁻¹¹ m³/(kg·s²).

Converting: ρ₀ in kg/m³, c_s in m/s, and working in natural units where
M_proton and R_proton set the scale:

    G_eff = κ/(4πc_s²ρ₀²) = G_Newton

    → c_s²ρ₀² = κ/(4πG_Newton) = κ × 1.19×10⁹ kg²/(m³·s²)

    With ρ₀ in nucleon mass/fm³ = 1.67×10⁻²⁷/(5.624×10⁻¹⁶)³ kg/m³
                                 = 9.39×10¹⁸ kg/m³

    → c_s² = κ × 1.19×10⁹ / (9.39×10¹⁸)² = κ × 1.35×10⁻²⁹ m²/s²

    → c_s = κ^{1/2} × 1.16×10⁻¹⁴·⁵ m/s ≈ κ^{1/2} × 3.67×10⁻¹⁵ m/s

This is c_s/c ≈ 1.2×10⁻²³. Extremely slow. This means in the simple scalar
theory, the sound speed must be 10²³ times slower than light to match G_Newton.

**Interpretation**: This tiny ratio is EXPECTED — it reflects the weakness of
gravity. G is weak because c_s ≪ c, meaning density perturbations propagate
very slowly compared to twist waves. The field resists density changes very
strongly (large α ∝ 1/c_s²), which makes depletion tiny and the gravitational
force weak.

**But**: c_s ~ 10⁻¹⁵ m/s means gravity propagation time across the solar
system would be ~10²⁶ seconds ≈ 10¹⁹ years. This is incompatible with the
observation that gravity propagates at (or near) light speed.

**This is the scalar gravity problem**: to get weak gravity (small G), need
small c_s. But small c_s means slow gravity propagation. GR avoids this because
the graviton is massless and moves at c — the weakness comes from the Planck
mass, not from a slow propagation speed.

### The way forward

The scalar picture (density sound waves = gravity) has a fatal flaw: c_s must be
simultaneously very small (weak G) and very large (fast propagation ~c). This is
impossible for a single scalar mode.

Resolution: the TENSOR structure may save this. If the frame field q (not just ρ)
carries the gravitational interaction, then:
- The gravitational "wave" moves at c (twist wave speed, not sound speed)
- The weakness of G comes from a small coupling between ρ and the tensor sector
- The BLV effective metric (v2 Path 4) shows the frame field DOES produce
  tensor-like effects

This points toward a hybrid theory where:
- SCALAR gravity (density depletion, 1/r) provides the Newtonian potential
- TENSOR gravity (frame torsion, c-speed propagation) provides the relativistic
  corrections (light bending factor 2, GW at speed c)

The v6 numerical program should test BOTH sectors.

---

## Part 8: Numerical Results (depletion.c)

### Test 1: Self-consistent density profile (PASS)

At α = 100, B=1, σ-model (m_π=0):

    a = 1.421  (vs 1.420 σ-model — tiny ρ correction)
    E₂ = 51.30, E₄ = 51.64  (virial E₂ ≈ E₄, confirmed)
    ρ(0)/ρ₀ = 0.985  (1.5% depletion at core)
    ΔQ = 0.257
    ΔQ_predicted = E₂/(2αρ₀) = 0.257
    Ratio = 1.003  (universality formula works to 0.3%)

Self-consistent iteration converges in ~15 iterations with under-relaxation (ω=0.3).

### Test 2: Long-range density tail (CRITICAL FINDING)

    δρ(r) ~ 0.25 / r^5.95  ← actual density perturbation, decays as 1/r⁶
    Φ(r)  ~ 0.25 / r^1.00  ← Poisson integral, decays as 1/r

| r  | δρ(r)      | Φ(r) (Poisson) | ΔQ/r     |
|----|-----------|----------------|----------|
| 3  | -3.25e-04 | 8.23e-02       | 8.58e-02 |
| 5  | -1.76e-05 | 5.10e-02       | 5.15e-02 |
| 10 | -2.78e-07 | 2.57e-02       | 2.57e-02 |
| 20 | -4.53e-09 | 1.29e-02       | 1.29e-02 |

**The actual density δρ(r) decays as 1/r⁶, NOT 1/r.** The Poisson integral Φ
gives 1/r, but this is a mathematical construct — it does not represent the
physical density at distance r.

A test soliton at distance d couples to δρ(d), NOT to Φ(d). Therefore the
force is F ~ dδρ/dr ~ 1/r⁷ (SHORT-RANGE), not 1/r² (long-range).

### Test 3: Universality scan (PARTIAL PASS)

| α     | E₂     | ΔQ      | E₂/(2αρ₀) | Ratio | ρ(0)/ρ₀ |
|-------|--------|---------|-----------|-------|---------|
| 0.1   | 1.97   | 95.75   | 9.86      | 9.71  | 0.015   |
| 1.0   | 26.98  | 23.15   | 13.49     | 1.72  | 0.015   |
| 10    | 49.03  | 2.54    | 2.45      | 1.04  | 0.847   |
| 100   | 51.30  | 0.256   | 0.257     | 1.00  | 0.985   |
| 1000  | 51.53  | 0.0257  | 0.0258    | 1.00  | 0.998   |

Universality ΔQ = E₂/(2αρ₀) works to <1% for α ≥ 10 (perturbative regime).
Breaks down at small α where depletion is nonperturbative (ρ₀ → 0 at core).

### Test 5: Two-body interaction (CRITICAL FINDING)

| d  | E_overlap (physical) | E_Poisson (math) | Ratio     |
|----|---------------------|-----------------|-----------|
| 3  | -1.66e-02           | -3.49e-01       | 0.048     |
| 5  | -8.98e-04           | -2.09e-01       | 0.0043    |
| 10 | -1.42e-05           | -1.05e-01       | 0.000135  |
| 20 | -2.32e-07           | -5.23e-02       | 4.4e-06   |

    E_overlap(d) ~ 12.1 / d^5.93  (physical: SHORT-RANGE)
    E_Poisson(d) ~ 1.05 / d^1.00  (mathematical: 1/r)

**The physical two-body interaction decays as 1/r⁶, not 1/r.**

---

## Part 9: Diagnosis — Why Density Depletion Fails for 1/r

### The fundamental obstruction

The density perturbation δρ(x) = -|ω(x)|²/(4α) is determined by the LOCAL twist
energy density. Since |ω|² is localized (decays as 1/r⁶ for the hedgehog), δρ is
also localized. No 1/r tail.

The Poisson integral Φ = ∫δρ(x')/|x-x'| d³x' gives 1/r because it's a
mathematical convolution, but the density FIELD at distance r is δρ(r), not Φ(r).

For 1/r gravity, we need either:
1. δρ itself to go as 1/r at large r (requires a massless propagating density mode
   with a nonzero monopole source)
2. A SEPARATE mediating field φ with ∇²φ = (nonzero source) and test solitons
   coupling to φ

### Why option 1 fails

The source for the density Poisson equation is S(x) = ∇²|ω|², which satisfies:

    ∫S d³x = ∫∇²|ω|² d³x = ∮∇|ω|²·dA = 0    (Gauss's theorem)

**Zero total source → no 1/r monopole.** The density perturbation has no
monopole moment in the Poisson expansion. The leading far-field term is
quadrupole or higher, decaying as 1/r³ or faster.

The total deficit ΔQ = ∫(ρ₀-ρ) d³x ≠ 0 is NOT the Poisson source integral.
ΔQ measures the total density removed, but the Poisson source (Laplacian of the
driving force) integrates to zero by Gauss. The deficit is "real" but distributed
with fast decay, not with 1/r.

### Why option 2 requires Path 3

A mediating field φ with nonzero monopole source needs:

    ∇²φ = κ · J(x)    with ∫J d³x ≠ 0

The natural candidate is the topological charge density:

    J(x) = B⁰(x) = -f'sin²f/(2π²r²),    ∫B⁰ d³x = B ∈ Z

This has NONZERO total integral (topological charge). The Poisson solution:

    φ(r) → κB/(4πr)    for r ≫ R_core

This IS 1/r! And it's exactly Path 3: ∇²p = g_top·B⁰.

**The density depletion mechanism correctly identifies the SOURCE (ΔQ ∝ mass,
universal), but fails to provide the PROPAGATION (no massless mediator with
monopole source). Path 3's topological charge B⁰ provides both source
(∫B⁰ = B ≠ 0) and propagation (constraint field, massless).**

### Synthesis: what v6 actually needs

The complete gravity mechanism requires THREE ingredients:

1. **Source**: ΔQ ∝ E₂ ∝ mass (from density depletion — WORKS)
2. **Mediator**: massless field with ∇²φ = source, ∫source ≠ 0 (Path 3 — NEEDS g_top)
3. **Coupling**: test soliton responds to φ (through B⁰_test — WORKS if g_top ≠ 0)

The v6 Lagrangian provides ingredient 1 but not 2. Adding a coupling term
g_top·p·B⁰ to the Lagrangian provides ingredient 2, but g_top is then free.

The original v2 problem remains: g_top is a free parameter.

### Updated assessment

| Feature | Status | Source |
|---------|--------|--------|
| Density conservation | CONFIRMED | ΔQ = E₂/(2αρ₀) |
| Universality (ΔQ ∝ M) | CONFIRMED | Ratio = 1.003 for α ≥ 10 |
| Core depletion | CONFIRMED | ρ(0)/ρ₀ = 0.985 at α=100 |
| 1/r interaction | FAILS | δρ ~ 1/r⁶, not 1/r |
| Mediating field needed | YES | Path 3 (B⁰p) or equivalent |
| Coupling determined | NO | g_top still free |

The density conservation mechanism is CORRECT about the source (what gravitates
and why) but WRONG about the propagation (how the effect reaches distant solitons).
The propagation requires a new field degree of freedom that couples to topological
charge, not to density.

---

## Part 10: Causal Density Propagation (ds = c·dt)

### The key insight: the equilibrium is non-local

The algebraic equilibrium δρ = -|ω|²/(4α) is WRONG because it's LOCAL — the
density at each point depends only on the twist at that same point. Physical
density propagates causally at the speed of gravity c, with a Heaviside step
function θ(ct - r) enforcing that information hasn't arrived beyond the wavefront.

The correct equation is the MASSLESS wave equation:

    □δρ = -½|ω|²    where □ = (1/c²)∂²/∂t² - ∇²    (Eq. 20)

with c = speed of gravity (= c_light in code units).

### Why this bypasses the Gauss obstruction

The fluid wave equation (from continuity + Euler) gives:

    ∂²δρ/∂t² = c_s²∇²δρ + ½∇²|ω|²    (fluid model)

whose source ∇²|ω|² has **zero monopole** by Gauss's theorem.

But the massless wave equation (Eq. 20) has source ½|ω|², which has
**nonzero monopole**:

    ∫½|ω|² d³x = E₂/ρ₀ = 51.3 ≠ 0    (Eq. 21)

The Gauss obstruction (∫∇²f = 0) only applies when the source IS a Laplacian.
The massless wave equation's source is |ω|² itself, not its Laplacian.

### The retarded solution

The static solution to ∇²δρ = ½|ω|² gives, for r ≫ R_core:

    δρ(r) → -(E₂/ρ₀)/(4πr)    (Eq. 22)

This is 1/r — **long-range density depletion**.

The full retarded (causal) solution:

    δρ(r,t) = -(Q_eff/4πr) · θ(ct - r)    (Eq. 23)

where Q_eff = E₂/ρ₀ and θ is the Heaviside step function. The density is:
- Zero for r > ct (information hasn't arrived yet)
- 1/r for r < ct (established behind the wavefront)

This IS the "speed of gravity step-wise function" ds = c·dt.

### Numerical verification

Solved the 1D wave equation for u = r·δρ with leapfrog integration:

    ∂²u/∂t² = c²u'' - ½c²r|ω|²

Results at t = 50 (wavefront at r = 50):

| r | δρ_wave | -Q_eff/(4πr) | ratio |
|---|---------|------------|-------|
| 5.5 | -7.43e-1 | -7.44e-1 | 0.999 |
| 7.7 | -5.32e-1 | -5.31e-1 | 1.002 |
| 10.8 | -3.81e-1 | -3.80e-1 | 1.004 |
| 15.1 | -2.72e-1 | -2.71e-1 | 1.004 |
| 21.1 | -1.94e-1 | -1.94e-1 | 1.004 |
| 29.5 | -1.39e-1 | -1.38e-1 | 1.004 |
| 41.3 | -9.91e-2 | -9.88e-2 | 1.003 |

**Power-law fit (r = 3..20): exponent = -0.989 (essentially -1.0)**

The fluid model (Model B) converges to δρ ~ r^{-5.96} (the algebraic equilibrium),
confirming that the Gauss obstruction kills 1/r when conservation is strict.

### Conservation analysis

Model A (massless wave) total deficit:

    ∫δρ·d³x grows as ½c²Q_eff·t²

Numerically verified: at t=5, integral = -643.8 vs prediction -641.3 (0.4% match).

In finite volume V, enforce conservation via Lagrange multiplier μ:

    □δρ = -½|ω|² + μ,    μ = E₂/(ρ₀V)

The source becomes: -½|ω|² + μ (localized depletion + uniform replenishment).
For r ≪ V^{1/3}: δρ ≈ -Q_eff/(4πr) (1/r dominates).
For r ~ V^{1/3}: uniform correction restores conservation.

In infinite volume (V → ∞): μ → 0, and the 1/r extends to all r. The
compensating surplus is a compression spike at the wavefront r = ct, pushed
to infinity as t → ∞.

### Two-body interaction from the massless model

Soliton 1 creates density field: δρ₁(d) = -Q_eff/(4πd)

Soliton 2's energy shift: δE = ½∫δρ₁|ω₂|²d³x ≈ δρ₁(d) × Q_eff

Interaction energy:

    E_int(d) = -Q_eff²/(4πd) = -(E₂/ρ₀)²/(4πd)    (Eq. 24)

Force: F = -dE/dd = -Q_eff²/(4πd²) → INVERSE SQUARE LAW

Effective gravitational coupling:

    G_eff = Q_eff²/(4π M²) = (E₂/ρ₀)²/(4π M²)    (Eq. 25)

With E₂ ≈ M/2: G_eff ≈ 1/(16πρ₀²) ≈ 0.020 code units.

**Coupling strength**: G_eff/G_Newton ≈ 3.6 × 10⁴⁰ — still 10⁴⁰ too strong.
This is the nuclear/Planck hierarchy: the twist coupling is at the nuclear scale
(~1 fm) while gravity operates at the Planck scale (~10⁻³⁵ m).

### What changed vs Part 9

| Feature | Part 9 (algebraic) | Part 10 (massless wave) |
|---------|-------------------|----------------------|
| Density equation | δρ = -\|ω\|²/(4α) local | □δρ = -½\|ω\|² causal |
| Source monopole | ∫∇²\|ω\|² = 0 | ∫½\|ω\|² ≠ 0 |
| Far-field decay | 1/r⁶ | 1/r |
| Interaction | short-range | LONG-RANGE 1/r² |
| Conservation | strict (∫δρ = 0) | grows; needs V → ∞ |
| G_eff/G_Newton | 10³⁸ (with α) | 10⁴⁰ (without α) |

### Remaining problems

1. **Coupling too strong by 10⁴⁰**: The twist-density coupling is at the nuclear
   scale. Need a suppression mechanism (small coupling constant, or weak
   mixing with a separate gravitational sector).

2. **Conservation is approximate**: The massless wave equation doesn't strictly
   conserve ∫ρ. In infinite volume this is fine (surplus at ∞), but in any
   finite computation it requires a Lagrange multiplier or V → ∞ limit.

3. **Tensor structure**: The density interaction is SCALAR (spin-0), not
   spin-2 tensor as required by GR. Still need tensor modes for gravitational
   waves.

4. **Physical origin of □δρ = -½|ω|²**: Why is the density equation massless?
   The fluid equations give a massive mode (c_s from pressure). The massless
   equation requires either (a) zero pressure, or (b) a separate field that
   couples to density but has no mass term.

---

## Part 11: The Coupling is Internal Binding, Not Gravity

### The 10⁴⁰ ratio is the nuclear/Planck hierarchy

The density 1/r interaction strength Q_eff²/(4π) = 209 code = 1905 MeV at
unit distance. Comparing with known force scales:

| Distance | d (fm) | E_density (MeV) | E_Newton (MeV) | Comparison |
|----------|--------|-----------------|----------------|------------|
| 1 code | 0.56 | 1905 | 5.4×10⁻³⁸ | sub-nucleon |
| 2.7 code | 1.5 | 714 | 2.0×10⁻³⁸ | nucleon size |
| 5 code | 2.8 | 381 | 1.1×10⁻³⁸ | nuclear distance |
| 10 code | 5.6 | 191 | 5.4×10⁻³⁹ | nuclear range |
| 100 code | 56 | 19 | 5.4×10⁻⁴⁰ | atomic scale |

The density 1/r interaction is at the **nuclear scale**, not gravitational.
Nuclear potential well depth (~50-100 MeV) is reached at d ≈ 20-40 code (11-22 fm),
which is longer than pion Yukawa range but in the right order of magnitude for a
long-range nuclear residual.

### Self-energy: 68% of soliton mass

The gradient energy stored in the 1/r density field (with UV cutoff at R_core):

    E_self = Q_eff²/(8π R_core) ≈ 70 code = 635 MeV    (Eq. 26)

    E_self / M_soliton = 635/937 = 0.678 = 68%

This is not a small correction — the density self-interaction accounts for a
**majority** of the soliton's rest mass. The 1/r density field IS the binding
energy that holds the soliton together as a coherent object.

### Interpretation: input → output coupling

The coupling Q_eff connects the soliton's internal structure (twist field, the
"input") to its external dynamics (translational motion, the "output"):

    twist → density depletion → back-pressure → cohesion → inertia

The twist energy |ω|² creates the density deficit (source of the wave equation).
The 1/r density field creates a potential well that confines the twist. This
self-consistent binding is what makes the soliton behave as a **particle** —
a localized object that moves as a unit.

The self-energy E_self ≈ 0.68 M means the soliton's mass is primarily
self-interaction energy. Without the density binding, the twist field would
disperse (Derrick's theorem for pure E₂). The density field provides the
confinement mechanism, analogous to how the color force confines quarks in QCD.

### What this means for gravity

The density 1/r interaction is the **nuclear/self-binding force**, not gravity:

| Property | Density 1/r | Gravity |
|----------|------------|---------|
| Strength | ~200 MeV·fm | ~10⁻³⁸ MeV·fm |
| Role | Internal binding | Long-range attraction |
| Mediator | Density wave (scalar, massless) | Graviton (tensor, massless) |
| Self-energy | 68% of M | ~10⁻³⁹ of M |
| Analogy | Strong force (QCD) | Gravity (GR) |

Gravity requires a **much weaker** mechanism. The density interaction, while
giving the correct 1/r form and inverse-square force law, is simply too strong
to be gravity — it's the force that makes solitons into particles.

The 10⁴⁰ ratio between the density force and gravity is the same hierarchy
as between the strong nuclear force and gravity in standard physics
(g_s²/(G_N m_p²) ≈ 10³⁹-10⁴⁰).

---

## Part 12: Impact on EMF Propagation

### The conformal coupling problem

The twist field Lagrangian E₂ = ½ρ|∂q|² has ρ multiplying BOTH the time
and space derivatives equally:

    L = (ρ/2c²)|∂_t q|² - (ρ/2)|∇q|²    (Eq. 27)

The Euler-Lagrange equation for a perturbation δq in background ρ(r):

    (ρ/c²) ∂²δq/∂t² = ∇·(ρ∇δq) = ρ∇²δq + (∇ρ)·(∇δq)

Dividing by ρ:

    ∂²δq/∂t² = c²∇²δq + c²(∇ ln ρ)·∇δq    (Eq. 28)

The eikonal limit (δq = A·e^{iS}, short wavelength):

    Leading order:  ω² = c²k²     →  v_phase = c  (ALWAYS)    (Eq. 29)
    Subleading:     amplitude transport with ∇ln ρ correction

**The phase velocity is c everywhere, regardless of ρ(r).** Rays travel in
straight lines. The 1/r density profile creates NO gravitational lensing.

This is because the coupling is **conformal**. The effective metric for wave
propagation is:

    g₀₀^eff = -ρ/c²,    g_{ij}^eff = ρ δ_{ij}

The null cone: ρ(-c²dt² + dx²) = 0 → dx = c·dt. The factor ρ cancels.

Compare with GR, where g₀₀ = -(1-2Φ/c²) and g_{rr} = (1+2Φ/c²) have
DIFFERENT Φ-dependence, producing lensing.

### What the 1/r density DOES to EMF

While not bending rays, the density gradient does affect waves:

1. **Amplitude modulation**: The ∇ln ρ·∇δq term in Eq. 28 modifies wave
   amplitude (not phase). Near the soliton (r ~ 2 code), |∇ln ρ| ≈ 1.0 —
   a 100% amplitude effect per wavelength. At r = 10 code, it's 4%.

2. **Energy density variation**: Wave energy E = ½ρ|∂q|² depends on ρ
   directly. In the 1/r density background, wave energy density varies
   even for constant amplitude.

3. **Massive mode potential**: For pion-mass modes (m_π ≠ 0), the
   dispersion relation becomes:

       ω² = c²k² + c²m_π²/ρ(r)    (Eq. 30)

   The effective mass m_eff² = m_π²/ρ varies as 1/r through δρ. This
   creates a long-range potential for **massive particles only** — massless
   EMF modes are unaffected.

4. **Density wave as new field**: The massless density oscillation is itself
   a propagating degree of freedom (scalar, speed c). It couples to EMF
   through the intensity ½ρ|∂q|² — an optical Kerr-like interaction where
   EM intensity drives density perturbations, and density modulates EM
   amplitude.

### What would produce gravitational lensing

For the 1/r density to bend light, the time and space kinetic terms need
**different** ρ-dependence:

    L = (ρᵃ/2c²)|∂_t q|² - (ρᵇ/2)|∇q|²    (Eq. 31)

Currently a = b = 1 (conformal, no lensing). If a ≠ b:

    v²/c² = ρ^{b-a}    →    n(r) = ρ^{(a-b)/2}    (Eq. 32)

For a = 2, b = 1 (example): n = √(ρ/ρ₀), and with δρ = -Q/(4πr):

    Φ_eff/c² = -(n²-1)/2 ≈ Q_eff/(8πρ₀r)

This would give gravitational lensing with deflection angle θ ~ Q_eff/(4πρ₀b).

Three possible sources of non-conformal coupling:

1. **E₄ Skyrme term**: contributes differently to P and m, but only where
   sin f ≠ 0 (soliton core). Vanishes in far field.

2. **Higher-order terms (L₆, L₈)**: may have different ρ-powers in time vs
   space. The sextic L₆ = λ₆(B⁰)² is ρ-independent, so it doesn't help.

3. **Non-minimal coupling**: If the Lagrangian contains cross-terms like
   ρ²|∂_t q|² (from the full Cl⁺(3,0,1) algebra expansion), the time and
   space parts would have different ρ-exponents. This requires going beyond
   the simple E₂ = ½ρ|∂q|² coupling.

### Summary: EMF and the density wave

| Effect | Massless EMF | Massive modes | Density wave |
|--------|-------------|---------------|-------------|
| Speed in 1/r background | c (unchanged) | c (to leading order) | c |
| Deflection by soliton | NONE (conformal) | Mass-dependent | N/A (is the field) |
| Amplitude modulation | YES (∇ln ρ term) | YES | N/A |
| Coupling to density | Intensity ∝ \|∂q\|² | Intensity + mass shift | Self-interaction |
| Analogy | Photon in dilaton background | Massive scalar + dilaton | Dilaton field |

The density wave is a **dilaton-like** massless scalar. It mediates the 1/r
nuclear binding force between solitons. It modulates EMF amplitude but does
not deflect EMF rays. For gravitational lensing of light, a non-conformal
coupling mechanism is needed.
