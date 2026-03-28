# V47 Stability Analysis

## Finding: η(P) = η₀ + η₁|P| is UNSTABLE for η₁ > ~7

### Test Results (N=64, T=10, proton template)
| η₁ | Stable? | E_total at T=10 |
|-----|---------|-----------------|
| 0 | YES | 3455 |
| 1 | YES | 3471 |
| 2 | YES | 3493 |
| 5 | YES | 3321 |
| 8 | **NO** | NaN |
| 50 | **NO** | NaN (even at dt=0.001) |
| 115 | **NO** | NaN |

### Root Cause
The coupling η(P)|P| creates positive feedback:
1. Large |P| → large η_eff → strong curl force
2. Strong curl force → larger φ oscillation
3. Larger oscillation → larger |P|
4. Goto 1 (runaway)

V(P) saturates P itself, but the CURL FORCE depends on ∇×θ which
is proportional to η_eff — it grows without bound.

### Solution: Saturating η(P)

Replace η(P) = η₀ + η₁|P| with:

    η(P) = η₀ + η₁ × P² / (1 + β × P²)

- At small P: η ≈ η₀ + η₁P² (quadratic onset — gentler than linear)
- At large P: η → η₀ + η₁/β (bounded — coupling saturates)
- Maximum η_eff = η₀ + η₁/β (can tune independently)

With η₁=500, β=50: η_max = η₀ + 500/50 = 10.5
At P=0.1: η_eff = 0.5 + 500×0.01/(1+50×0.01) = 0.5 + 5/1.5 = 3.8
At P_opt=0.082: η_eff = 0.5 + 500×0.0067/(1+50×0.0067) = 0.5 + 3.36/1.34 = 3.0

This is still amplified relative to η₀=0.5 but naturally bounded.

### Alternative: η(P) with existing V(P) saturation form

    η(P) = η₀ × (1 + κ_η × P²) / (1 + β_η × P²)

This uses the SAME saturating structure as V(P). At P→0: η=η₀.
At P→∞: η = η₀ × κ_η/β_η. The ratio κ_η/β_η sets the amplification.

### Recommended Implementation

Use: η(P) = η₀ + eta1 × P² / (1 + kappa_eta × P²)

Two new parameters:
- eta1: coupling amplification strength
- kappa_eta: saturation scale (reuse κ=50 from V(P) for simplicity)

With kappa_eta = κ = 50:
- P_sat = 1/√κ_eta = 0.141 (coupling saturates around |P| ≈ 0.14)
- η_max = η₀ + eta1/kappa_eta = 0.5 + eta1/50

For η_max ≈ 8: eta1 = (8-0.5)×50 = 375
For η_max ≈ 12: eta1 = (12-0.5)×50 = 575
