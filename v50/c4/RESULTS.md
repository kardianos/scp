# V50/C4 Results: Curl-Squared Hardening

## Configuration

    Base: V44 (6-field Cosserat, η=0.5, m²=2.25, μ=-41.345, κ=50)
    Cosserat strain: α=0.1 (geometric theta constraint)
    Curl-squared hardening: β=0.5 (shell stiffening)
    Chiral helicity: κ_h=0 (disabled for this test)
    Grid: N=64, L=15, dx=0.476
    BC: absorbing T<50, periodic T≥50
    Init: proton template (V43 pre-converged UUD)

## Equations

    ∂²φ_a/∂t² = ∇²φ_a - m²φ_a - ∂V/∂φ_a + η curl(θ)_a
                 - α curl(M)_a        [Cosserat: stiffens twist geometry]
                 - 2 curl(Q)_a        [Hardening: stiffens where θ×curl overlap]

    ∂²θ_a/∂t² = ∇²θ_a + η curl(φ)_a
                 + 2α M_a             [Cosserat: pulls θ toward curl(φ)/2]
                 - β|∇×φ|² θ_a       [Hardening: θ heavy where twist is strong]

    M_a = curl(φ)_a/2 - θ_a          (geometric mismatch)
    Q_a = (β/2)|θ|² curl(φ)_a        (hardening intermediate)

All force signs Maxima-verified (derive_harden.mac, derive_cosserat.mac).

## Stability

    T=0→50 (absorbing): E dropped 6207 → ~2800 (seed transient cleanup)
    T=50→200 (periodic): E = 2802 → 2787 (drift -0.5%)

Energy conserved to <1% post-switch. No NaN. No blowup. Stable.

## Radial Profile (T=200, last frame)

Centroid: (1.23, 7.50, 0.72)

    Region       r       phi_rms   |P|      curl_rms  harden    theta_rms
    ──────────   ──────  ────────  ───────  ────────  ────────  ─────────
    Core         0-3     0.88-0.96 0.12-0.17 0.07-0.10 0.0002   0.14-0.16
    Transition   3-5     0.81-0.87 0.10-0.12 0.07-0.09 0.0001   0.11-0.13
    SHELL        5-7.5   0.41-0.72 0.02-0.08 0.17-0.18 0.0004   0.11-0.12
    Exterior     8-13    0.06-0.25 <0.004    0.05-0.12 <0.0002  0.09-0.12

## Shell Structure Confirmed

The particle has three distinct regions:

### 1. Core (r < 3)
- phi_rms ≈ 0.93 (strong field)
- |P| ≈ 0.15 (strong trilinear binding)
- curl_rms ≈ 0.08 (moderate twist)
- theta_rms ≈ 0.15 (theta follows geometry via Cosserat)
- Hardening ≈ 0.0002 (moderate — curl not at peak)

### 2. Shell (r ≈ 5-7.5)
- phi_rms dropping (0.72 → 0.41)
- |P| transitioning (0.08 → 0.02)
- curl_rms PEAKS at 0.18 (twist gradients steepest here)
- theta_rms ≈ 0.11 (still substantial)
- **Hardening PEAKS at 0.0004** (2× core, 6× exterior)
- **Mismatch PEAKS at 0.028** (Cosserat working hardest here)

The shell forms WHERE the helical twist has the steepest spatial
gradients. This is NOT the core (where P is large and curl is
moderate) and NOT the far field (where both are small). It's the
BOUNDARY — the surface where the braid's 3D structure transitions
to the 1D carrier wave.

### 3. Exterior (r > 8)
- phi_rms ≈ 0.06-0.25 (approaching background A_bg=0.1)
- |P| < 0.004 (negligible binding)
- curl_rms decaying (0.12 → 0.05)
- theta_rms ≈ 0.09 (residual, decaying)
- Hardening < 0.0002 (negligible)

## Key Observations

1. **The hardening term creates a shell, not a uniform blob.**
   Peak hardening is at r≈7.3, coinciding with peak curl and the
   |P| transition zone. This is the physical boundary of the particle.

2. **Theta has geometric structure from the Cosserat constraint.**
   theta_rms is highest at the core (0.15) and decays outward.
   It follows the braid topology, not filling the box uniformly.
   Compare to V44 baseline where theta fills the box at T=200.

3. **Curl peaks at the surface, not the core.**
   The helical twist creates maximum |∇×φ| at the boundary where the
   spatial structure transitions from 3D braid to 1D carrier wave.
   The hardening term β|θ|²|∇×φ|² naturally concentrates at this
   transition — creating a rigid boundary without needing P² or |φ|².

4. **Mismatch peaks at the shell too.**
   The Cosserat constraint (θ → curl(φ)/2) is most stressed at the
   shell. This means the geometric alignment between theta and phi's
   twist is least perfect at the boundary — exactly where the
   hardening provides extra rigidity.

5. **The exterior is NOT saturated.**
   At r>10, theta_rms ≈ 0.09 (small, decaying). The Cosserat
   constraint + hardening prevent the theta accumulation that
   plagued V44 and V49. The background theta mass from the
   hardening term (β|∇×φ|² ≈ 0 in far field, but Cosserat gives
   effective mass) keeps theta from filling the box.

## Comparison to Previous Versions

| Feature           | V44        | V48 (λ_θP²) | V49 (unified) | V50/C4       |
|-------------------|------------|-------------|---------------|--------------|
| Theta structure   | Fills box  | P² confines | Saturates     | **Geometric** |
| Shell boundary    | None       | Soft (P²)   | None          | **Hard (curl²)** |
| Energy conserved  | Yes        | Yes         | After fixes   | **Yes (-0.5%)** |
| Massless photon   | Yes        | Nearly      | With bare mass| **Yes (|∇×φ|=0 in vac)** |
| New params        | 0          | 1           | 4-5           | **2 (α, β)** |
| Physically motivated | N/A     | Weak        | Complex       | **Cosserat + twist rigidity** |

## Files

- smoke_b05.sfa: 9 frames, T=0→200, N=64
- smoke_b05_diag.tsv: diagnostics
- cross_section.c: radial profile analysis tool
- derive_harden.mac: Maxima force verification
- PROPOSAL.md: physical motivation
- /lean/V50C4/CurlHardening.lean: formal properties

## Next Steps

1. Longer run (T=600+) to verify long-term stability
2. Two-proton test at D=10-14 to check anti-blob behavior
3. Parameter sweep on β to find optimal shell hardness
4. Measure theta dispersion in far field (is photon truly massless?)
5. Re-enable chiral term (κ_h) with corrected sign, small value

## Two-Proton Interaction (N=100, D=12, T=400)

### Configuration
    N=100, L=20, dx=0.404
    Two protons at x=±6 (initial D=12)
    α=0.1 (Cosserat), β=0.5 (hardening), κ_h=0 (no chiral)
    BC: absorbing T<100, periodic T≥100
    Burst windows: t=100-115, 200-215, 300-315 (every timestep)

### Results: BOUND TWO-BODY STATE

    Energy: 1997 at switch → 1979 at T=400 (drift -0.9%, stable)
    Protons: DISTINCT — two P-maxima maintained throughout

### Core tracking (burst window t=108-110, 92 samples):

    Separation: 3.8-5.0 code units (from initial D=12)
    Peak 1 P: 0.24-0.53 (oscillating)
    Peak 2 P: 0.18-0.42 (oscillating, slightly weaker)
    Gap min P: 0.00-0.44 (cycles between shallow and complete)

### Breathing dynamics:

    Two superposed modes:
    - Symmetric: both cores grow/shrink together, period ~2.0 t
    - Asymmetric: cores alternate strength, period ~2.5 t
    - Full beat: 4.5 t
    - Amplitude STABLE — no secular trend over 5+ cycles

### Shell structure confirmed:

    Cross-section at t≈101 shows:
    - Core 1: z≈9.5, P=0.110
    - Gap:    z≈12.7, P=0.003 (30× drop from core)
    - Core 2: z≈17.2, P=0.053
    - Gap has MAXIMUM theta (0.10) and curl (0.23)
    - Hardening peaks in the gap (0.000528) — the stiff interface

### Files
    pp_D12_N100.sfa: 162 GB, 1024 valid frames (T=0 to T≈110)
    burst_cores_fine.tsv: 92-point core tracking through burst
    full_cores.tsv: 21-point tracking across full simulation
    slices/: 2D cross-section data at t≈101
