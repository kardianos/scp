# Variant A — Findings: C ⊗ H ⊗ O Algebra Construction

**Date**: 2026-05-22
**Script**: `01_choh_algebra.py`

---

## Result: Algebra built and verified

- **Dimension**: 64 over ℝ (= 2 × 4 × 8).
- **Basis**: $e_{(a,b,c)}$ with $a \in \{0,1\}$ (ℂ), $b \in \{0,1,2,3\}$ (ℍ), $c \in \{0,...,7\}$ (𝕆).
- **Multiplication**: $(e_{a_1 b_1 c_1})(e_{a_2 b_2 c_2}) = (ℂ$-mul on $a$) $\otimes$ (ℍ-mul on $b$) $\otimes$ (𝕆-mul on $c$).
- **Identity**: $1 \otimes 1 \otimes e_0$.
- **Non-associative**: associator on three octonionic basis elements has norm 2 (expected).
- **Alternative law VIOLATED**: $(xy)x \neq x(yx)$ for general elements. This is correct mathematics — tensor products do not preserve alternativity in general; only 𝕆 alone is alternative.
- **Structure constants**: 4096 nonzero out of 262144 total (1.56% density).

## Key Algebraic Notes

The alternative-law violation is **not a bug** — it is a mathematical feature. Furey's program uses ℂ ⊗ ℍ ⊗ 𝕆 not because the algebra itself is "nice," but because its left-multiplication ideals decompose into the Standard Model fermion representations. The representation theory matters, not the algebraic axioms.

## Infrastructure Built

- `choh_structure.npz`: stored 64×64×64 structure tensor for use by downstream variants.
- Multiplication function `CHO_mul(x, y)` operating on 64-element real arrays.
- Conjugation function `CHO_conj(x)` (factorwise conjugation).

This is the foundation for Variants B–G.

## Status

Variant A: **COMPLETE**. Proceed to Variant B (SM idempotents via ℂ ⊗ 𝕆 ≅ Cl(6)).
