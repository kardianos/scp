# Step 2: Dimensional / Naturalness Analysis for η₁

## The Question

Given parameters {m², μ, κ, η₀}, what combinations have the correct
"dimensions" (in code units where everything is dimensionless) and
magnitude for η₁?

## Parameter Inventory

| Parameter | Value | Physical role |
|-----------|-------|---------------|
| m² | 2.25 | Field mass squared |
| μ | -41.345 | Potential coupling |
| κ | 50 | Potential saturation |
| η₀ | 0.5 | Base φ-θ coupling |
| A_bg | 0.1 | Background amplitude |

All are dimensionless in code units (c=1, lattice units).

## Candidate Expressions for η₁

η₁ multiplies |P| in the coupling, so η₁|P| must have the same units
as η₀. Since |P| = |φ₀φ₁φ₂| ~ A³, we need [η₁] = [η₀]/[A³].
In code units this is just a number.

### Ratios involving two parameters:

| Expression | Value | η_eff at |P|=0.1 | Notes |
|------------|-------|---------|-------|
| |μ|/η₀ | 82.7 | 8.77 | **In binding range** |
| κ/η₀ | 100 | 10.5 | Close but large |
| |μ|·η₀ | 20.7 | 2.57 | Too small |
| m²/η₀ | 4.5 | 0.95 | Too small |
| √κ/η₀ | 14.1 | 1.91 | Too small |
| |μ|/m² | 18.4 | 2.34 | Too small |
| κ·η₀ | 25 | 3.0 | Borderline |
| |μ|/A_bg³ | 41345 | Way too large |
| η₀/(A_bg³) | 500 | Way too large |

### Ratios involving three parameters:

| Expression | Value | η_eff at |P|=0.1 |
|------------|-------|---------|
| |μ|/(η₀²) | 165.4 | **Wait — this gives η₁·η₀ = |μ|/η₀** |
| |μ|·κ/η₀ | ~4000 | Way too large |
| m²·κ/η₀ | 225 | Too large |

## Analysis

### The leading candidate: η₁ = |μ|/η₀

Value: 82.7. At |P|=0.1: η_eff = 8.77. Right in the binding range.

Physical interpretation:
- |μ| is the strength of the binding potential V(P)
- η₀ is the strength of the EM coupling
- The ratio |μ|/η₀ measures how much STRONGER binding is than EM
- η₁ = |μ|/η₀ means: "the topology amplification brings the coupling
  up to the scale of the binding potential"

This makes physical sense: the nuclear force should be as strong as
the binding potential, not as strong as the EM coupling.

### Why not κ/η₀ = 100?

This gives η_eff = 10.5 at core — slightly too strong. But it's close.
The saturation parameter κ controls WHERE the potential peaks. Using
κ instead of |μ| gives a coupling that responds to the saturation
scale rather than the potential depth. Less physically motivated.

### The virial derivation (Step 1) gives:

η₁ = -R/(2η₀J₁) where R is the virial residual.

If η₁ = |μ|/η₀, then R = 2|μ|J₁. This means the virial residual is
proportional to |μ| × (topology-weighted curl integral). The residual
is large precisely because the potential is strong — the soliton
WANTS a stronger coupling to reach Derrick equilibrium.

## Conclusion

**η₁ = |μ|/η₀ is the most natural candidate.** It connects the
topology amplification to the binding potential strength, uses only
existing parameters, and gives the correct magnitude for nuclear
binding. Numerical verification needed (Step 3).

## Alternative: η₁ derived from P_opt

The optimal binding point P_opt = 1/√(3κ) ≈ 0.082 is where V'(P) is
maximized. If we require η_eff(P_opt) to satisfy some condition:

η₀ + η₁·P_opt = target → η₁ = (target - η₀)/P_opt

For target = √(|μ|) ≈ 6.43: η₁ = 5.93/0.082 = 72.3 → η_eff = 6.43
For target = |μ|^(1/3) ≈ 3.45: η₁ = 2.95/0.082 = 36 → η_eff = 3.45

The |μ|/η₀ = 82.7 gives η_eff(P_opt) = 0.5 + 82.7×0.082 = 7.28.
This is close to √(|μ|·η₀) = √(20.67) = 4.55. Not an exact match
to any simple function of the parameters.

The numerical cross-check (Step 3) will determine whether η₁ = |μ|/η₀
is exact or approximate.
