# Variant C — Findings: U(1)_em Coupling from Algebra

**Date**: 2026-05-22
**Script**: `03_alpha_from_u1.py`

---

## Result: Algebra confirmed; α not predicted

**Tr(Q²) over one SM generation = 8/3** — an exact algebraic value, identifying:
- $\nu$: Q² = 0
- $e$: Q² = 1
- $u$ (3 colors): 3 × (2/3)² = 4/3
- $d$ (3 colors): 3 × (1/3)² = 1/3
- **Sum: 8/3**

This is the standard QED β-function coefficient per generation, confirming the algebra matches SM phenomenology in this sector.

## Tested Normalizations and Outcomes

| Normalization | Predicted α | Ratio to empirical |
|---------------|-------------|---------------------|
| $g^2 = 16\pi^2/$Tr(Q²) | 4.71 | 646× too large |
| $1/\alpha = 4\pi \cdot $Tr(Q²) | 33.5 | 4× too small |
| $1/\alpha = 4\pi \cdot $Tr(Q²) × 21 | 703.7 | 5.1× too large |
| $21 \cdot 2\pi^2$ | 414.5 | 3× too large |

**None of the simple natural normalizations gives α = 1/137** without a tuned factor.

## Honest Assessment

α prediction from the Furey algebra alone is NOT immediate. The algebra gives the correct quantum-number content (charges 0, ±1/3, ±2/3, ±1; SU(3)_c triplet) but does not by itself determine the absolute strength of the U(1)_em coupling.

This is consistent with the broader experience of the Wyler-style derivations: getting numerically close to 137 requires a specific homogeneous space volume, and the choice of that space is not uniquely fixed by the algebra alone.

## Status

Variant C: **PARTIAL** (algebra confirmed, α not predicted). Proceed to Variant D for the instanton-action approach (the $S_{\rm em} \approx \pi^2/2$ conjecture from step 11).
