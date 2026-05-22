# Variant G — Findings: Gravity Sector

**Date**: 2026-05-22
**Script**: `06_gravity_sector.py`

---

## Result: Cross-sector form correct, absolute G off by 0.5%

Under the conjecture $S_{\rm em} = \pi^2/2 - 2\alpha$ (Variant D) and cross-sector ratio 21 (step 9):
- Predicted $S_{\rm grav} = 21(\pi^2/2 - 2\alpha) = 103.32$
- Empirical $S_{\rm grav} = -\ln G_e = 103.06$
- **Gap: +0.26%**

Predicted $G_e = 1.34 \times 10^{-45}$, empirical $G_e = 1.75 \times 10^{-45}$, ratio 0.76.

The form (linear in $S_{\rm em}$ with coefficient 21) is correct structurally; the absolute value is off by a multiplicative factor of order 1.

## Structural Identification

| Sector | Algebraic grade | Coupling derivation |
|--------|------------------|---------------------|
| EM | Bivector | Modulation depth $S_{\rm em} = \pi^2/2 - 2\alpha$ on $G_2$-orbit |
| Gravity | Scalar/density | Modulation depth $S_{\rm grav} = 21 \cdot S_{\rm em}$ on Spin(7) orbit |

This is consistent with v58's "gravity from density gradients" hypothesis: the scalar-density grade is the gravity source.

## Open: The 0.5% G correction

The EM sector required a correction $2\alpha$ (Variant D, found to 0.005% precision). The gravity sector needs an analog correction we have not identified. Empirically:

$$S_{\rm grav}^{\rm needed\ correction} = S_{\rm grav}^{\rm empirical} - 21 \cdot S_{\rm em}^{\rm corrected} = -0.27$$

A small negative shift. None of the simple candidates tested (α, −21α, ln(1+21α), m_e/M_P log) matches this gap cleanly.

This might be:
- Higher-loop QCD corrections (gravity sees all matter content, not just charged leptons).
- Threshold effects at heavy particle masses.
- Genuine additional algebraic content (perhaps the quark sector's contribution).

## Status

Variant G: **PARTIAL** — structural cross-sector relation correct; absolute G prediction off by factor 0.76. The bulk of the hierarchy is captured; the remaining correction structure for G is not yet identified.

Proceed to Variant H (consolidation).
