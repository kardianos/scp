# V50 Comparison: C4 vs C5

## The Shared Equation System

Both C4 and C5 build on V44 + Cosserat strain + curl²-hardening.
They differ in ONE term: the chiral/twist coupling in the phi sector.

### Shared terms (identical in both):

    ∂²φ_a/∂t² = ∇²φ_a - m²φ_a - ∂V/∂φ_a
                 + η curl(θ)_a                   [V44: EM coupling]
                 - α curl(M)_a                   [Cosserat: geometric constraint]
                 - 2 curl(Q)_a                   [Hardening: shell stiffness]
                 + (variant-specific term)

    ∂²θ_a/∂t² = ∇²θ_a
                 + η curl(φ)_a                   [V44: EM coupling]
                 + 2α M_a                        [Cosserat: pull θ → curl(φ)/2]
                 - β|∇×φ|² θ_a                   [Hardening: θ heavy at twist]

    M_a = curl(φ)_a/2 - θ_a         (Cosserat mismatch)
    Q_a = (β/2)|θ|² curl(φ)_a       (hardening intermediate)

### Shared parameters:

    m² = 2.25, μ = -41.345, κ = 50, η = 0.5
    α_cs = 0.1  (Cosserat strain)
    β_h = 0.5   (curl² hardening)

### The variant-specific term:

**C4 (no additional phi term):**

    κ_h = 0, λ_tw = 0
    No chiral or twist coupling. Pure V44 + Cosserat + hardening.

**C5 (twist stiffening):**

    λ_tw = 1.0
    L_twist = (λ/2) P² [φ · curl(φ)]²

    Phi force: λ × [2P dPda h² + 2P²h curl(φ)_a - ∇(P²h) cross terms]

    This rewards strong helical twist (h² = [φ·curl(φ)]²) at braid
    cores (P²≠0), regardless of handedness. Makes the braid core
    a stiffer topological object — harder to deform or untwist.

## Test Configuration (identical for both)

    N=100, L=20, dx=0.404
    Two protons at x=±6 (initial D=12)
    BC: absorbing T<100, periodic T≥100
    Burst windows: t=100-115, 200-215, 300-315

## Energy Comparison

| Time | C4 E_total | C5 E_total | Difference |
|------|-----------|-----------|------------|
| t=0 | 4293 | 4293 | 0.00% |
| t=100 (switch) | 1997 | 1997 | 0.00% |
| t=200 | 1988 | 1988 | 0.01% |
| t=400 | 1979 | 1979 | 0.01% |

Energy evolution is **identical to 0.01%**. The twist stiffening at
λ_tw=1.0 has negligible impact on total energy dynamics. Both runs
conserve energy to <1% over 300 time units of periodic BC.

## Field Evolution Comparison

| Metric (t=400) | C4 | C5 |
|----------------|-----|-----|
| E_phi_kin | 873 | 875 |
| E_theta_kin | 195 | 196 |
| E_grad | 80 | 80 |
| E_mass | 718 | 716 |
| E_pot | -34.0 | -33.6 |
| E_tgrad | 178 | 178 |
| E_coupling | -31.1 | -31.1 |
| phi_max | 0.678 | 0.677 |
| P_max | 0.292 | 0.289 |
| P_int | 46.6 | 46.3 |
| theta_rms | 0.0313 | 0.0313 |

All values match to <1%. The twist stiffening (λ_tw=1.0) produces
no measurable difference in any diagnostic quantity.

## Cross-Section Comparison (t≈100)

| Metric | C4 | C5 |
|--------|-----|-----|
| Centroid | (-0.26, -1.28, 8.89) | (-0.25, -1.31, 8.94) |
| Peak hardening radius | 2.02 | 2.02 |
| Peak hardening value | 0.000188 | 0.000192 |
| Particle boundary | 5.62 | 5.62 |

Identical shell structure. Same core-shell-exterior organization.

## Interpretation

At λ_tw=1.0, the twist stiffening term has **no measurable effect**
on the two-proton dynamics. This is because:

1. The term L_twist = (λ/2)P²h² is localized to braid cores (P²≠0)
   and proportional to h² = [φ·curl(φ)]². At the braid core, h² is
   already at its equilibrium value from V(P) binding. The additional
   reward for twist doesn't change the equilibrium — it just makes
   the minimum deeper, which doesn't affect the dynamics unless the
   twist is being perturbed.

2. The Cosserat constraint (α=0.1) already aligns θ with the geometric
   twist. The curl²-hardening (β=0.5) already creates the shell. The
   twist stiffening is redundant — it rewards a topology that's already
   stable.

3. The term would matter more if the braids were actively being
   deformed (e.g., during a close collision or merger). At D=12
   with the current parameters, the protons interact through their
   shells but don't deform each other's cores enough to activate
   the twist reward.

## Conclusion

**C4 is sufficient.** The twist stiffening (C5) adds computational
cost (P and h at 6 neighbors for the ∇P²h gradient) without
measurable benefit at the current parameter regime. The three-term
system (V44 + Cosserat + curl²-hardening) is the minimal effective
equation set.

C5 should be kept as an option for future experiments where core
deformation becomes significant (e.g., closer initial separation,
higher collision velocities, or anti-proton interactions where
opposite-chirality braids may interact differently). But for
standard nuclear binding studies, C4 is the recommended baseline.

## Recommended Equation Set (C4)

    ∂²φ_a/∂t² = ∇²φ_a - m²φ_a - ∂V/∂φ_a + η curl(θ)_a
                 - α curl(M)_a - 2 curl(Q)_a

    ∂²θ_a/∂t² = ∇²θ_a + η curl(φ)_a + 2α M_a - β|∇×φ|² θ_a

    M_a = curl(φ)_a/2 - θ_a
    Q_a = (β/2)|θ|² curl(φ)_a
    V(P) = (μ/2) P² / (1 + κP²)

    Parameters: m²=2.25, μ=-41.345, κ=50, η=0.5, α=0.1, β=0.5
    All force signs Maxima-verified.
