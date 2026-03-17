# V24-UDD: Neutron-Like 120° Oscillon (One Up + Two Down)

## Thesis

The UDD counterpart of V24-UUD. One field ("up") has mass m_U, two fields
("down") have mass m_D. All three maintain 120° phase through pairwise
coupling. This models the neutron (udd).

## EOM

    ∂²φ₁/∂t² = ∂²φ₁/∂x² - m_U²φ₁ - λ(φ₂+φ₃) - ∂V_triple/∂φ₁   [Up]
    ∂²φ₂/∂t² = ∂²φ₂/∂x² - m_D²φ₂ - λ(φ₁+φ₃) - ∂V_triple/∂φ₂   [Down]
    ∂²φ₃/∂t² = ∂²φ₃/∂x² - m_D²φ₃ - λ(φ₁+φ₂) - ∂V_triple/∂φ₃   [Down]

## Mass Matrix

    M² = |m_U²  λ    λ  |
         |λ     m_D²  λ  |
         |λ     λ    m_D²|

Eigenvalues:
- Mode 1 (mostly symmetric): involves all three
- Mode 2 (DD antisymmetric): ≈ m_D² - λ (pure D-D, doesn't involve U)
- Mode 3 (UD antisymmetric): ≈ mixed

The DD antisymmetric mode has mass m_D² - λ. If m_D < m_U, this mode
is LIGHTER than the UU antisymmetric mode in the UUD case. This means
the UDD oscillon (neutron) has a lower internal frequency than UUD (proton)
for the same m_D < m_U.

## Initialization

    φ₁ = A_U·g(x)·cos(0)                        [Up]
    φ₂ = A_D·g(x)·cos(2π/3)                     [Down]
    φ₃ = A_D·g(x)·cos(4π/3)                     [Down]
    v₁ = 0
    v₂ = -ω_D·A_D·g(x)·sin(2π/3)
    v₃ = -ω_D·A_D·g(x)·sin(4π/3)

## Parameter Scan

Same as V24-UUD: fix λ=0.85, m_U=1.0, scan m_D.

## Key Questions

1. Does the UDD energy differ from UUD at same m_D?
   (Analog of proton-neutron mass splitting)
2. Is one more stable than the other?
3. At what m_D does the mass ordering flip
   (UDD heavier vs lighter than UUD)?
4. Does the DD antisymmetric mode create different radiation
   than the UU mode in the UUD case?

## Output

- `src/udd120.c` — solver (identical structure to uud120.c but φ₂,φ₃ = D)
- `data/udd_mD{val}_ts.tsv`
- `data/udd_summary.tsv`
- `RESULTS.md`

## Parameters

μ=-20, κ=20, m_U=1.0, λ=0.85
m_D scan: {1.0, 0.95, 0.90, 0.85, 0.80, 0.70, 0.50}
A_U=A_D=0.8, σ=3.0
Nx=4000, xmax=100, tfinal=10000

Compile: `gcc -O3 -Wall -o udd120 src/udd120.c -lm`
