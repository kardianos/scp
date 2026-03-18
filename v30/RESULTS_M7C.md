# V30 M7+c(ρ_B): Best Model — STABLE AT ALL α_c

## The Model

S (braid) propagates at c_eff set by B (background):
    c_eff²(x) = 1 - α_c × (1 - ρ_B(x)/ρ₀)

B propagates at c=1 always. Coupling: g×S²×B (B→S conversion).

## Result: Braid Survives at All Tested α_c

| α_c | fc_final | min_c² (depletion) | avg_c²(core) (accretion) | Stable |
|-----|----------|-------------------|-------------------------|--------|
| 0.000 | 0.951 | 1.000 | 1.000 | YES |
| 0.050 | 0.840 | 0.950 | 1.020 | YES |
| 0.100 | 0.924 | 0.900 | 1.038 | YES |
| 0.200 | 0.978 | 0.800 | 1.076 | YES |
| 0.500 | 0.950 | 0.500 | 1.192 | YES |

N=128, L=20, T=500, A_bg=0.1, g=0.01.

## Key Finding: Simultaneous Accretion + Depletion

The B-field shows TWO effects simultaneously:
- **Core (r<3)**: B ACCUMULATES (ρ_B > ρ₀, c_eff > 1)
  The S→B radiation and direct coupling enrich B near the braid.
- **Intermediate (r≈5-15)**: B DEPLETED (ρ_B < ρ₀, c_eff < 1)
  The braid consumes B from the surrounding region.

This creates a c_eff profile:
- Core: c > 1 (fast pocket — braid propagates quickly)
- Shell: c < 1 (slow zone — other waves get trapped)
- Far field: c = 1 (normal)

## Dark Matter Halo Analog

The B-density profile (excess at core, deficit at shell, normal at far field)
matches the structure of a dark matter halo:
- Excess at center → more "mass" than luminous matter alone
- Extended slow zone → flat rotation curves
- The braid IS the luminous matter; the B-field IS the dark matter

At α_c=0.5: core c²=1.19, shell c²=0.50. Ratio = 2.4×.
This is a significant speed variation that would strongly affect
test-particle geodesics.

## Why This Works (Unlike V30 Uniform c(ρ))

The M7 separation is ESSENTIAL:
- S (braid) → creates the depletion in B
- B (background) → sets c_eff for S
- The braid doesn't freeze because c_eff comes from B, not S
- B has its own dynamics (c=1 always), so it flows freely

With uniform c(ρ) on all fields (V30 earlier attempts):
- Dense → slow → frozen (no dynamics)
- OR dense → fast → blowup
Both fail. M7 avoids this by separating source from medium.

## The Complete Model Equation

    S fields (braid):
    ∂²S_a/∂t² = c_eff²(ρ_B) × ∇²S_a - m²S_a - ∂V_S/∂S_a - g×(ΣB²)×S_a

    B fields (background):
    ∂²B_a/∂t² = ∇²B_a - m²B_a - g×(ΣS²)×B_a

    c_eff²(x) = 1 - α_c × (1 - ρ_B(x)/ρ₀)
    ρ_B(x) = Σ_a [½(∂_t B_a)² + ½m²B_a²]

    Parameters: m=1.50, μ=-41.3, κ=50, g=0.01, α_c=0.1-0.5

## Next Steps

1. **Measure the c_eff radial profile** at each α_c (is depletion zone 1/r?)
2. **Two braids with c(ρ_B)** — does the c-depletion create attraction?
3. **Pinned boundaries** — does the system reach steady-state accretion/depletion?
4. **Gravitational wave test** — accelerate a braid, measure c_eff ripples

## Files
- data/m7c/summary.tsv — final metrics for all α_c
- data/m7c/ts_ac*.tsv — time series for each α_c value
