# V50 Results

## Selected Equation Set: C4

V50 explored returning to V44's clean geometric equations and adding
minimal terms to solve the blob merger problem while preserving gravity,
EM, and massless photons.

### The C4 equations (V44 + Cosserat strain + curl²-hardening):

    ∂²φ_a/∂t² = ∇²φ_a - m²φ_a - ∂V/∂φ_a + η curl(θ)_a
                 - α curl(M)_a - 2 curl(Q)_a

    ∂²θ_a/∂t² = ∇²θ_a + η curl(φ)_a + 2α M_a - β|∇×φ|² θ_a

    M_a = curl(φ)_a/2 - θ_a           (Cosserat geometric mismatch)
    Q_a = (β/2)|θ|² curl(φ)_a         (hardening intermediate)
    V(P) = (μ/2) P² / (1 + κP²)       (saturating binding potential)

### Parameters:

    m² = 2.25, μ = -41.345, κ = 50, η = 0.5
    α = 0.1   (Cosserat strain coupling)
    β = 0.5   (curl-squared hardening)

### All force signs Maxima-verified:
    - Cosserat phi force: -α curl(M)  (derive_cosserat.mac)
    - Hardening phi force: -2 curl(Q) (derive_harden.mac)
    - Both theta forces verified algebraically

### Key results:

1. **Two-proton bound state confirmed** (N=100, D=12, T=400)
   - Two P-maxima maintained throughout, separated by ~3.8-8.6 code units
   - Theta-rich hardened gap between cores prevents merger
   - Compound breathing with 4.5 time-unit beat pattern
   - Energy conserved to <1% over 300 time units of periodic BC
   - Protons breathing in symmetric + asymmetric modes, amplitude stable

2. **Shell structure confirmed** (cross-section at t≈101)
   - Core (r<3): high P, moderate curl, θ follows geometry
   - Shell (r≈5-7.5): peak curl, peak hardening, θ-rich barrier
   - Exterior (r>8): decaying fields, no saturation

3. **Theta follows braid geometry** (from Cosserat constraint)
   - No box saturation (the V44/V49 problem)
   - θ ≈ curl(φ)/2 at equilibrium
   - Massless photon in vacuum (|∇×φ|=0 → no hardening mass)

## C5 Comparison

C5 adds a twist stiffening term: L_twist = (λ/2) P² [φ·curl(φ)]²
with λ_tw=1.0. This rewards strong helical twist at braid cores
regardless of handedness.

**Result: C4 and C5 are identical to 0.01% across all metrics.**
Energy, field amplitudes, shell structure, breathing dynamics — all
match. The twist stiffening has no measurable effect because the
Cosserat + hardening system already stabilizes the braid topology.

C5 is retained as an option for future experiments where core
deformation matters (close collisions, anti-proton interactions).
For standard nuclear physics, C4 is the minimal effective system.

Note: C5's handedness-neutral twist stiffening (h²) is distinct
from the chiral handedness preference (linear h, κ_h from C2/C3).
The chiral term could have cosmological implications (parity
violation, baryogenesis) but is not needed at nuclear scale.

## Files

    v50/c4/   — C4 implementation, analysis tools, results
    v50/c5/   — C5 implementation (twist stiffening variant)
    v50/c2/   — C2 (chiral only, no Cosserat — superseded)
    v50/c3/   — C3 (Cosserat + chiral — intermediate step)
    v50/view/ — Preview converter and viewer improvements
    v50/COMPARISON_C4_C5.md — detailed comparison
    v50/HARDWARE_NOTES.md — performance analysis and EPYC 9255
    v50/proposal_*.md — design proposals from multiple sources
