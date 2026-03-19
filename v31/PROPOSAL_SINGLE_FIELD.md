# Single-Field c(ρ): No Split, No Smoothing

## The Physics

One field φ_a (3 components). The field IS the medium. c depends on
local energy density directly:

    ∂²φ_a/∂t² = c²(ρ(x)) × ∇²φ_a - m²φ_a - ∂V/∂φ_a

    ρ(x) = Σ_a [½v_a² + ½|∇φ_a|² + ½m²φ_a²]

    c²(x) = ρ₀ / ρ(x)    [simplest physical formula]

Dense field → slow c. Dilute → fast c. No arbitrary bounds, no smoothing.

The braid is just a concentrated region of φ. The "background" is the
same φ at lower amplitude. No distinction imposed.

## The Numerical Challenge

The braid core has ρ ≈ 100ρ₀ → c ≈ 0.1 (10× slower than ambient).
The vacuum has ρ → 0 → c → ∞ (if truly empty).

With uniform background A_bg=0.1: ρ_min ≈ 0.03 → c_max ≈ √(ρ₀/0.03) ≈ 6.
The c ratio across the domain is ~60:1 (0.1 to 6).

Standard Verlet with global dt: dt = CFL × dx / c_max = 0.15 × 0.3 / 6 = 0.008.
At T=500: 62500 steps. Each step: same as current. Total: ~60× slower.

This is the baseline cost. Sub-cycling and implicit methods reduce it.

## Proposal C: Sub-Cycling (local adaptive timestepping)

### Concept
Divide the domain into regions with similar c. Fast regions (high c,
low ρ) take many sub-steps per global step. Slow regions (low c,
high ρ) take fewer sub-steps or a single step.

### Implementation: Multi-rate Verlet

1. Compute c²(x) = ρ₀/ρ(x) at all points. No smoothing.

2. Classify each point into SPEED TIERS:
   - Tier 0: c ≤ 1 (slow, includes braid core)
   - Tier 1: 1 < c ≤ 2
   - Tier 2: 2 < c ≤ 4
   - Tier 3: c > 4 (fast, vacuum edges)

3. Set dt_base from Tier 0 (slowest): dt_base = CFL × dx / 1.0
   Each higher tier sub-cycles with dt = dt_base / 2^tier.

4. Per global step (dt = dt_base):
   - Update ALL tiers once with dt_base (safe for Tier 0)
   - Update Tiers 1-3 again with dt_base (2× sub-steps for Tier 1)
   - Update Tiers 2-3 again twice (4× sub-steps for Tier 2)
   - Update Tier 3 again four times (8× sub-steps for Tier 3)

5. At each sub-step for a tier:
   - Recompute forces ONLY for points in that tier + their neighbors
   - Update vel, phi ONLY for points in that tier
   - Boundary between tiers: use interpolated values from coarser tier

### CPU Cost
- Tier 0 (core, ~5% of domain): 1 update per step
- Tier 1 (~15%): 2 updates
- Tier 2 (~30%): 4 updates
- Tier 3 (~50%): 8 updates

Effective: 0.05×1 + 0.15×2 + 0.30×4 + 0.50×8 = 5.55 force evaluations per step.
Standard Verlet: 1 force per step but with dt 8× smaller → 8 steps.
Net: sub-cycling is ~0.7× the cost of naive small-dt Verlet. Modest savings.

BUT: the real savings come from the global dt being set by Tier 0 (c=1)
instead of Tier 3 (c=6). dt_base = CFL×dx/1 vs dt_naive = CFL×dx/6.
So 6× fewer global steps. With 5.55× cost per global step: net 6/5.55 = 1.08×.
Almost breakeven! But we get correct physics without any capping or smoothing.

### Complexity
Medium-high. Need: tier classification, partial force updates, inter-tier
boundary handling, multi-rate Verlet bookkeeping. ~500 lines of new code.

## Proposal D: Semi-Implicit with Local Gauss-Seidel

### Concept
Instead of explicit Verlet (which needs dt < dx/c), use an implicit
scheme that's unconditionally stable. Because the equation is LOCAL
(each point depends only on neighbors), we can solve it with local
Gauss-Seidel iteration — no global matrix needed.

### Implementation: Crank-Nicolson with local sweep

The wave equation discretized:
    (φ^{n+1} - 2φ^n + φ^{n-1})/dt² = c² × ∇²[½(φ^n + φ^{n+1})]
                                        - m² × ½(φ^n + φ^{n+1})
                                        - V'(φ^n)

Rearrange for φ^{n+1}:
    φ^{n+1}(x) = [2φ^n - φ^{n-1} + dt²(c²∇²φ_avg - m²φ_avg - V')]

where φ_avg = ½(φ^n + φ^{n+1}).

This is implicit in φ^{n+1} (it appears on both sides through ∇²φ_avg).

Solve by GAUSS-SEIDEL: sweep through all grid points. At each point,
compute φ^{n+1} from the CURRENT values of its neighbors' φ^{n+1}
(already updated earlier in the sweep) and φ^n (known from last step).

For a single point with 6 neighbors:
    φ^{n+1}(x) = [2φ^n(x) - φ^{n-1}(x)
                   + ½dt²c²/dx² × Σ_neighbors φ^{n+1}(neighbor)
                   + ½dt²c²/dx² × Σ_neighbors φ^n(neighbor)
                   - ½dt²(6c²/dx² + m²) × φ^n(x)
                   - dt²V'(x)]
                  / [1 + ½dt²(6c²/dx² + m²)]

Each point's update is one division + ~20 multiplications. Local only.

Convergence: 5-20 Gauss-Seidel sweeps per timestep (for wave equations,
red-black ordering converges fast).

### The c self-consistency
c²(x) depends on ρ(x) which depends on φ and v. During the GS sweep,
ρ changes → c changes. We can either:
(a) Lag c by one step: c²(x) from ρ^n (simple, first-order in time)
(b) Update c during the sweep: c²(x) from current φ^{n+1} estimate (more accurate)

Option (a) is simpler and usually sufficient.

### CPU Cost
Per GS sweep: 1× force evaluation (same as one Verlet step).
Sweeps per step: 10-20.
But dt can be 10-100× larger than explicit (unconditionally stable).

At 15 sweeps, dt 10× larger: effective cost = 15/10 = 1.5× explicit.
At 15 sweeps, dt 50× larger: effective cost = 15/50 = 0.3× explicit.

For a 60:1 c ratio: explicit needs dt = dx/6. Implicit can use dt = dx/1 (or larger).
So 6× larger dt. At 15 sweeps: 15/6 = 2.5× cost. Modest increase.

For a 1000:1 c ratio (extreme): explicit needs dt = dx/1000. Implicit: dt = dx/1.
At 15 sweeps: 15/1000 = 0.015× cost. MASSIVE savings.

### Complexity
High. Need: Gauss-Seidel sweep with red-black ordering, Crank-Nicolson
formulation, convergence check, c-update logic. ~700 lines of new code.
But the code is straightforward (no global matrix, no CG solver).

## Recommendation

**For exploration: Proposal C** (sub-cycling). It's explicit (familiar,
debuggable), handles the 60:1 c ratio, and the physics is exact. The
implementation is moderately complex but well-understood.

**For production: Proposal D** (implicit GS). Unconditionally stable,
handles ANY c ratio, and becomes more efficient as the c ratio increases.
But harder to debug and verify.

**Either way**: single field, no split, no smoothing. The c(ρ) formula
is c² = ρ₀/ρ (physical, no bounds).

## Parameters

    φ_a: 3 real fields (braid + background in ONE set of arrays)
    c²(x) = ρ₀ / ρ(x), ρ = Σ[½v² + ½|∇φ|² + ½m²φ²]
    V(P) = (μ/2)P²/(1+κP²), P = φ₀φ₁φ₂
    m = 1.50, μ = -41.3, κ = 50
    Init: braid + uniform background A_bg=0.1
    N=128, L=20, dx=0.31

## Grid & Memory

Single field: 3 × 3 arrays × 128³ × 8 bytes = 151 MB.
Plus ρ array + c² array: +34 MB. Total: ~185 MB. Trivial.

For N=256 (higher resolution): 1.2 GB. Still fine.
