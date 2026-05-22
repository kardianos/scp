# Kernel Fit 03 — Findings: Cross-Sector α Test

**Date**: 2026-05-22
**Script**: `03_em_bivector_coupling.py`
**Target**: derive the fine-structure constant α = 1/137.036 from the same Cl(3,1) algebra that gave the lepton spectrum.

---

## Headline Result

**Cl(3,1) alone does not predict α.** The natural bivector-vector coupling magnitude in the algebra is O(1), the right order of magnitude for a UV gauge coupling but off by a factor of ~3 from the empirical IR coupling $g_{\rm em} = \sqrt{4\pi\alpha} \approx 0.303$. No clean dimensionless combination of the Brannen phase $\varphi = 2/9$ rad and α emerges.

The lepton sector and EM sector, in pure Cl(3,1), are **structurally decoupled**. This is information, not failure.

## What Was Computed

1. **Natural Cl(3,1) bivector coupling**: $g_{\rm nat} = 1$ in unit-normalized algebra. Compared to $g_{\rm em} \approx 0.303$, the ratio is 3.3 — same order of magnitude but not equal.

2. **Wyler's 1969 formula**: $\alpha^{-1} = (9\pi^4/16)(2\pi^5/120)^{1/4} \approx 137.036$ — agrees to ~5 digits with the empirical value. The numerical coincidence is striking; the mathematical justification (specific homogeneous space volume ratio) is widely regarded as arbitrary. Recorded for completeness, not adopted.

3. **Cross-sector ratios** between Brannen $\varphi = 2/9$ rad and α: no clean combination. $\varphi/\alpha \approx 30.45$, $\varphi^2/\alpha \approx 6.77$, $\alpha/\varphi^2 \approx 0.148$, $\exp(-\varphi/\alpha) \approx 10^{-14}$ — none of these are small rationals.

## Honest Interpretation

The "unified multivector force" claim of v58 had two distinct senses:

- **Qualitative**: one equation $\to$ multiple sectors by grade projection. *This works* — Newton-like and Maxwell-like behaviors do emerge from different grade projections of the same multivector field.

- **Quantitative**: cross-sector predictions (relating G, α, m_μ/m_e, etc.). *This does not work in Cl(3,1) alone.* The lepton sector is fixed by the vector-grade Z₃ structure with two empirical parameters (|ξ|², φ); the bivector sector has its own coupling strength that is structurally independent.

After step 5, the kernel-fit program has produced:

- **Strong result for the lepton sector**: Brannen form is structural; Koide and Brannen phase are isolated to two empirical inputs (|ξ|² = 1/2 and φ = 2/9 rad); both can be pinned by a Z₃-invariant potential.

- **Null result for cross-sector unification**: Cl(3,1) alone cannot produce a numerical relation between the lepton mass scale and α. The sectors are decoupled within this algebra.

This is exactly the right kind of result for a Kepler-stage program: positive in the place where the data lives (lepton spectrum), negative where the natural extension was speculative (cross-sector predictions).

## What Could Close the α Gap

Listed for record, in order of standard-physics plausibility:

1. **RG running from a UV-natural value**: at high energy, gauge couplings unify and are O(1). RG flow from the GUT scale to the IR gives the running value α at any specific scale. This is the standard explanation in GUT theories. It doesn't predict α uniquely without the unification scale and the full RG content, but it explains why 1/137 isn't algebraically simple.

2. **Embedding in a larger algebra**: Cl(0,7) (octonions), Cl(3,1) ⊗ Cl(0,3), or exceptional Lie algebras like E_6, E_8 might have natural structure constants that produce 1/137. The Furey/Dixon program pursues this. Not yet quantitatively successful for α.

3. **Wyler-style geometric volume ratios**: numerically near-miss, mathematically arbitrary. Possibly meaningful if the right principle for selecting the homogeneous space is found.

4. **Vacuum manifold topology**: if the vacuum manifold has specific homotopy / Chern numbers, the EM coupling can be quantized in units of α. Requires identification of the manifold.

5. **Instanton suppression**: $\alpha \sim e^{-S_{\rm inst}}$ with $S_{\rm inst} \approx 5$ would give the right order. No first-principles derivation in this framework.

## Strategic Reckoning After 5 Steps

The kernel-fit program has produced:

| Item | Result |
|------|--------|
| Lepton mass form | Cl(3,1) Z₃ kernel **derives Brannen structure**. |
| Koide identity | Pinned by Z₃-invariant potential at $\|\xi\|^2 = 1/2$. |
| Brannen phase | Pinned by cubic phase parameter $\alpha_{\rm cubic} = 2/3$ rad. (Empirical, not yet derived.) |
| Lepton masses (e, μ, τ) | Reproduced to machine precision with 2 empirical parameters. |
| α (fine structure) | **Not predicted** by Cl(3,1) alone. |
| Cross-sector unification | **Quantitatively negative** in Cl(3,1) without additional structure. |

The framework is **partially successful and honestly limited**. It works for the lepton sector. It does not work as a quantitative unification of leptons + EM without extension.

## Recommendation Going Forward

Given these results, the kernel-fit program has three viable directions:

### Direction A: Refine the lepton sector
Pursue the structural conjecture α_cubic = Q rad inside the Cl(3,1) framework. If an algebraic operation in Cl(3,1) produces the cubic-phase coefficient = 2/3 rad naturally, then Koide AND Brannen become structurally derived. This would be a major positive result confined to the lepton sector. **Cost: medium. Risk: medium. Payoff: completes the lepton-sector ellipse.**

### Direction B: Extend the algebra to address α
Move from Cl(3,1) to Cl(0,7) (octonions) or Cl(3,1) ⊗ Cl(0,3). Larger algebras have more invariants. The Furey program suggests this contains the Standard Model fermion content as an automorphism class. Worth a serious look. **Cost: high. Risk: high. Payoff: potentially predicts α, gauge group, generations.**

### Direction C: Accept partial unification
Treat the lepton sector and the EM sector as independent applications of the GA framework, each with its own parameters. The success criterion shifts from "predict α from leptons" to "match more of the TYCHO table in any sector via natural algebraic structure." This is more conservative but more honest about what we've achieved. **Cost: low. Risk: low. Payoff: solidifies the lepton result without overpromising.**

My recommendation: **Direction A first** (it tests a specific conjecture inside the existing kernel), **then C** (consolidate the lepton sector as a Kepler ellipse), **then B** (extend the algebra for cross-sector predictions). Each step is bounded and falsifiable.

## Files

- `03_em_bivector_coupling.py` — the probe script.
- `03_findings.md` — this document.

## Status

| Target | Status |
|--------|--------|
| Brannen form | Derived (step 1–3) |
| Koide $\|\xi\|^2 = 1/2$ | Pinned by Z₃ potential (step 4) |
| Brannen $\varphi = 2/9$ rad | Pinned given $\alpha_{\rm cubic} = 2/3$ rad (step 4). The relation $\alpha_{\rm cubic} = Q$ rad is conjectural. |
| Lepton masses | Reproduced to machine precision (step 3) |
| α (fine structure) | **Not predicted by Cl(3,1) alone (step 5).** Requires algebra extension or new principle. |
| G (gravity) | Untested. Likely also not predicted by Cl(3,1) alone. |
| m_p/m_e | Untested. Requires confinement-sector physics beyond Cl(3,1). |

The lepton-sector Kepler ellipse is in hand. Cross-sector unification is not.
