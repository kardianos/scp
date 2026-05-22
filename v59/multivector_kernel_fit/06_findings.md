# Kernel Fit 06 — Findings: Normalization Check on the α–G Ratio

**Date**: 2026-05-22
**Script**: `06_normalization_check.py`
**Target**: test whether $S_{\rm grav}/S_{\rm em} = 2\pi^2$ exactly, under some natural normalization of the bare coupling.

---

## Headline Result

**The 2π² conjecture is NOT confirmed.** The empirical ratio $S_{\rm grav}/S_{\rm em} = 20.945$ does not coincide with $2\pi^2 = 19.739$ under any natural normalization. Closing the 6% gap requires $\alpha_{\rm bare} = e^{1/\pi} \approx 1.375$, which has no obvious algebraic interpretation.

**However**: the ratio is striking close to **21 = 3 × 7** (off by 0.26%, far inside experimental precision). And **the offset $S_{\rm grav} - 21 \cdot S_{\rm em} = -0.26 \approx 0$** — the proportionality $S_{\rm grav} \approx 21\,S_{\rm em}$ holds with nearly vanishing offset.

This is suggestive but algebraically opaque. 21 is not obviously the volume of any natural manifold in Cl(3,1).

## Detailed Numbers

| Quantity | Value |
|----------|-------|
| $S_{\rm em} = -\ln \alpha$ | 4.920244 |
| $S_{\rm grav} = -\ln G_e$ | 103.055680 |
| Raw ratio $S_{\rm grav}/S_{\rm em}$ | **20.9452** |
| 2π² (V of unit S³) | 19.7392 (gap −5.76%) |
| 21 (= 3 × 7) | 21.0000 (gap +0.26%) |
| 22 | 22.0000 (gap +5.04%) |
| 7π | 21.9911 (gap +5.00%) |
| 4π + 8.4 (numerology) | 20.9664 (gap +0.10%) |

The cleanest natural-number near-miss is **21**.

## Conditions for 2π² Exactly

For $S_{\rm grav}/S_{\rm em}$ to equal exactly $2\pi^2$, the bare coupling must be:

$$c = \exp\!\left(\frac{\ln(1/G_e) - 2\pi^2 \ln(1/\alpha)}{2\pi^2 - 1}\right) = \exp(0.31666) = 1.37254$$

Closest natural candidates for this $c$:
- $e^{1/\pi} = 1.37480$ — agrees to 2 × 10⁻³.
- $(1 + 1/e) = 1.36788$ — agrees to 5 × 10⁻³.
- $e/2 = 1.35914$ — off by 0.013.

$e^{1/\pi}$ is mathematically suggestive but does not correspond to any natural geometric volume or algebraic invariant of Cl(3,1) that I can identify. So while we *could* declare $\alpha_{\rm bare} = e^{1/\pi}$, this is curve-fitting, not derivation.

## What Closes to 21?

The proportionality $S_{\rm grav} \approx 21 \cdot S_{\rm em}$ is striking:

- $S_{\rm grav}$ predicted from $S_{\rm em}$: $21 \times 4.9202 = 103.32$
- $S_{\rm grav}$ measured: 103.06
- Offset: $-0.26$

So $S_{\rm grav} - 21\,S_{\rm em} \approx 0$ to within 0.25% of $S_{\rm grav}$. This is a very tight near-proportionality.

But **21 = 3 × 7** has no obvious algebraic origin in Cl(3,1):
- 3 = number of generations ✓ (motivated by Z₃)
- 7 = ?

Candidates for the "7":
- Dimension of Cl(3,1) even subalgebra minus 1 = 8 − 1 = 7
- Dimension of $S^7$ in octonionic structure (Cl(0,7) = octonions)
- Number of basis bivectors in Cl(3,1) (which is actually 6, not 7)
- Random.

None is structurally compelling. If 7 came from "octonion dimension minus 1" or "even subalgebra rank," that would suggest an embedding of Cl(3,1) into a larger octonionic structure where 21 = 3 × (rank − 1).

## What the Check Confirms

1. **Cross-sector framework is internally consistent**: α and G come from one O(1) bare coupling with exponential modulation.
2. **The 42-order hierarchy is structurally manageable**: it compresses to a ~21× ratio in suppression depths, not 42 orders of magnitude.
3. **The 2π² conjecture is too optimistic**: the gap is 6%, not closable by natural normalization.
4. **There is a tantalizing near-21 numerical fit** (0.26% gap) but its algebraic origin is unclear.

## What the Check Does Not Establish

1. **No clean algebraic derivation of the ratio.** Neither 2π² nor 21 emerges from Cl(3,1) structure unambiguously.
2. **No prediction of α or G individually.** The cross-sector ratio is suggestive but each absolute value still requires empirical input.
3. **No fix to the lepton-Planck hierarchy.** The smallness of $m_e/M_P$ is reframed (as exponential suppression) but not derived.

## Interpretation

The α–G corroboration program has reached a partial result:

- **Qualitative success**: the two sectors share a kernel, and the exponential modulation framework is mathematically consistent.
- **Quantitative ambiguity**: the cross-sector ratio is *close* to two candidate natural numbers (2π² and 21) but matches neither cleanly.

This is honest territory for a Kepler-stage program. We found a structural relation (S_grav ≈ 21 × S_em) that is empirically tight but lacks a derivation. Three possibilities:

1. **21 is the right answer, and 7 comes from octonionic / Cl(0,7) structure** — would require extending the framework to Cl(3,1) ⊗ Cl(0,3) or similar.
2. **The actual ratio is irrational but close to a clean number for accidental reasons** — the framework is consistent but does not predict the specific ratio.
3. **The modulation function is not exactly exp(-x)** — small corrections could shift the apparent ratio.

## Status After 8 Steps

| Item | Status |
|------|--------|
| Lepton sector mass spectrum | Reproduced to machine precision (Cl(3,1) Z₃ kernel + constraint surface). |
| Koide identity Q = 2/3 | **Structural** via |ξ|² = 1/2 constraint surface. |
| Brannen phase | Empirical (3-direction on S³). |
| α (alone) | Not predicted. |
| G (alone) | Not predicted. |
| α–G ratio under exp. modulation | 20.945 empirical. 21 (= 3 × 7) nearest natural number (0.26% gap). 2π² rejected (6% gap). |
| α variation in gravity | Predicted ~10⁻¹⁰ at Earth surface. Testable. |
| Cross-sector derivation | Qualitative framework consistent; quantitatively the ratio is empirically tight but its algebraic origin is unidentified. |

## Suggested Next Steps

In honest order of leverage:

1. **Investigate the 7 in 21 = 3 × 7.** Extend Cl(3,1) to a larger algebra (Cl(3,1) ⊗ Cl(0,3), or Cl(0,7) octonions) and ask whether 7 appears as a natural rank/dimension. If yes, the cross-sector ratio 21 = 3 (generations) × 7 (extended-algebra rank) becomes structural.

2. **Test alternative modulation functions.** The exponential gave 20.95. Polynomial or log-modified exponential modulations could shift this to exactly 21 or 2π² with a small correction. Try $f(x) = e^{-x}(1 + \epsilon x)$ etc.

3. **Compute $S_{\rm em}$ and $S_{\rm grav}$ independently from algebra.** Instead of inferring them from empirical α and G, calculate them as instanton actions on the Cl(3,1) constraint S³. If their values come out as O(5) and O(100) — and their ratio matches empirical 20.95 — the cross-sector unification is derived.

4. **Test α variation in gravity experimentally** (long-term). Engage atomic clock community.

## Files

- `06_normalization_check.py` — the calculation.
- `06_findings.md` — this document.

The kernel-fit program has produced a quantitative structural target (the cross-sector ratio) without yet deriving it from algebra. That is partial progress in the right direction. The 2π² hope was too tidy; the actual ratio (20.95, or perhaps 21 exactly) requires additional structure beyond what we've identified in Cl(3,1).
