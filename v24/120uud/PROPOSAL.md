# V24-UUD: Proton-Like 120° Oscillon (Two Up + One Down)

## Thesis

Break the S₃ symmetry of the three fields to create a UUD-like bound state.
Two fields ("up") have mass m_U, one field ("down") has mass m_D ≠ m_U.
All three maintain 120° phase relationships through the pairwise coupling.

This models the proton (uud): two quarks of one type, one of another,
bound in a color-singlet (120° phase) configuration.

The key physics: the 120° phase structure (from pairwise coupling) provides
the "color confinement." The mass difference m_U ≠ m_D provides the "flavor."
From outside, this looks like a particle with specific mass and charge.

## EOM

    ∂²φ₁/∂t² = ∂²φ₁/∂x² - m_U²φ₁ - λ(φ₂+φ₃) - ∂V_triple/∂φ₁   [Up]
    ∂²φ₂/∂t² = ∂²φ₂/∂x² - m_U²φ₂ - λ(φ₁+φ₃) - ∂V_triple/∂φ₂   [Up]
    ∂²φ₃/∂t² = ∂²φ₃/∂x² - m_D²φ₃ - λ(φ₁+φ₂) - ∂V_triple/∂φ₃   [Down]

## Mass Matrix

With per-field masses m₁²=m₂²=m_U², m₃²=m_D², and pairwise coupling λ:

    M² = |m_U²  λ    λ  |
         |λ     m_U²  λ  |
         |λ     λ    m_D²|

Eigenvalues:
- When m_U = m_D = m: standard case (m²+2λ, m²-λ, m²-λ)
- When m_U ≠ m_D: the degeneracy of the antisymmetric modes SPLITS

For small δm = m_D - m_U:
- Mode 1 (mostly symmetric): ≈ m² + 2λ + O(δm)
- Mode 2 (UU antisymmetric): ≈ m_U² - λ (unchanged, doesn't involve D)
- Mode 3 (UD antisymmetric): ≈ m² - λ/2 + O(δm) (shifted)

This splitting means the UUD oscillon has TWO distinct internal frequencies,
creating a richer oscillation pattern than the pure 120° state.

## Initialization

The UUD 120° state:
    φ₁ = A_U·g(x)·cos(0)
    φ₂ = A_U·g(x)·cos(2π/3)
    φ₃ = A_D·g(x)·cos(4π/3)
    v₁ = 0
    v₂ = -ω_U·A_U·g(x)·sin(2π/3)
    v₃ = -ω_D·A_D·g(x)·sin(4π/3)

where ω_U = √(m_U² - λ), ω_D = √(m_D² - λ), A_U and A_D may differ.

## Parameter Scan

Fix λ = 0.85 (strong pairwise coupling, near the 3ω < m threshold).
Scan the mass splitting:

    m_U = 1.0 (fixed)
    m_D ∈ {1.0, 0.95, 0.90, 0.85, 0.80, 0.70, 0.50}

At each m_D: measure the UUD oscillon energy (the "proton mass").
Compare with the equal-mass case (m_D = m_U = 1.0).

Also try the reverse: m_D > m_U (heavier down quark).

## Key Questions

1. Does the UUD oscillon survive with m_U ≠ m_D?
2. Does the 120° phase persist despite the mass asymmetry?
3. What is the UUD oscillon energy as a function of m_D/m_U?
4. How does it compare with the UDD oscillon (V24-UDD)?
5. Does the mass splitting create a DIFFERENT oscillation pattern
   (beat frequencies from ω_U ≠ ω_D)?

## Reference Code

- v21 1D solver: `/home/d/code/scp/v21/src/triad1d.c`
- V24-PW code (once written): for the pairwise coupling implementation

## Output

- `src/uud120.c` — solver with per-field masses + pairwise coupling
- `data/uud_mD{val}_ts.tsv` — time series per m_D
- `data/uud_summary.tsv` — energy vs m_D/m_U
- `RESULTS.md`

## Parameters

μ=-20, κ=20, m_U=1.0, λ=0.85
m_D scan: {1.0, 0.95, 0.90, 0.85, 0.80, 0.70, 0.50}
A_U=A_D=0.8 (equal initial amplitudes), σ=3.0
Nx=4000, xmax=100, tfinal=10000

Compile: `gcc -O3 -Wall -o uud120 src/uud120.c -lm`
