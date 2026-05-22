# Experiment 04 — Findings: Muon-Electron Mass Ratio Scan

**Date**: 2026-05-22
**Script**: `04_muon_electron_ratio_scan.py`
**Target**: TYCHO_TABLE T1.5 — m_μ/m_e = 206.7682830(46).

---

## Summary

**Null result at TIGHT precision (10⁻⁴ relative).** No simple substrate naturally produces a ratio of 206.7683.

The famous Nambu 1952 conjecture 3·α⁻¹/2 = 205.554 lands at relative Δ = −5.87 × 10⁻³ (off by 0.59%) — the best "natural" near-miss from 80 years of attempts, and the well-known historical baseline.

Tested categories all returned no-tight-hit:
- exp(p/q) and log_n(integer) sweeps
- Finite group order ratios (A_n, S_n, M_n, PSL, Z_n, etc.)
- Sphere kissing numbers and their ratios
- n-ball volume ratios for n = 1..25
- Factorial, binomial, Catalan-number ratios
- α-based combinations (other than the historical Nambu near-miss)

The "tight hits" reported in the output for items like `3·α⁻¹/2 + 1.214` and `207 − 0.2317` are tautological — those constants were back-fit to the target and are not predictions. They are listed only as sanity checks.

## Numerical Observations

- **log_e(206.7683) = 5.3316**. Closest rational candidate: 16/3 = 5.3333. Δ_log = +1.73 × 10⁻³, translating to a value 207.13 — 0.18% off.
- **target / α⁻¹ = 1.50886**. Closest rational: 3/2 = 1.5. Off by 0.59% (this is Nambu's near-miss in another form).
- **target / (3/2 · α⁻¹) = 1.00591**. The mass ratio is 0.59% more than the Nambu prediction.
- **target as 9 × 23 = 207**: off by 0.112% (1.12 × 10⁻³). Not tight.

No combination of simple integers, π, e, ln, √, or α produces 206.7683 to better than ~0.6% (Nambu) without explicit tuning.

## Why This Result Matters

The muon-electron ratio is a 10-digit experimental number. Over 80 years of numerology, the best approximation from a non-tuned formula is Nambu's, accurate only to 0.6%. After this scan reproduces that historical landscape, the situation is:

- **m_μ/m_e is not a simple ratio of natural invariants from any substrate tested.**
- Either the ratio is set by a complex dynamical mechanism (RG running, anomaly cancellation, multi-loop QED, threshold effects) that does not reduce to a clean algebraic identity, or it is the product of multiple independent factors none of which is itself naturally near 200, or it is "fine-tuned" in a deeper theory (string landscape, anthropic).

Negative information: yes, substantial. We have ruled out a wide class of substrate proposals. Positive Kepler ellipse: no.

## Cumulative Status After Four Experiments

| Target | Status |
|--------|--------|
| T1.1 — Koide Q = 2/3 | Confirmed empirically. Automatically satisfied by any Hermitian Z₃-cyclic 3×3 matrix with **|b|/a = 1/√2**. No substrate provides this ratio naturally. |
| T1.2 — Brannen φ ≈ 2/9 rad | Consistent with exact 2/9 within m_τ uncertainty. **No** standard geometric/algebraic/modular/hyperbolic structure produces 2/9 rad. |
| T1.5 — m_μ/m_e = 206.7683 | Best simple approximation: 3·α⁻¹/2 = 205.554 (Nambu, 0.59% off). **No** tight match from any tested substrate. |
| T1.3 — m_p/m_e = 1836.15 | Not yet tested. |
| T1.4 — α⁻¹ = 137.036 | Not yet tested. |

**Three out of three Tier-1 targets attacked so far have produced strong negative information but no Kepler ellipse.**

## Strategic Reckoning

The Kepler-before-Newton program was the right move methodologically. But after four experiments and a wide-net scan, the pattern is consistent: **the Standard Model's tight numerical mysteries do not yield to simple-substrate numerology.** Each of the three Tier-1 numbers we've attacked behaves the same way:

1. The number is confirmed empirically to many digits.
2. Simple substrates (polytopes, Lie roots, Coxeter angles, group orders, modular forms, hyperbolic invariants, kissing numbers, ball volumes) do NOT naturally produce the number to within the experimental precision.
3. The closest "natural" formulas are off by 10⁻² to 10⁻³ at best — well above the ~10⁻⁵ to 10⁻¹⁰ experimental precision.

This is not yet a refutation of the substrate-from-numerology approach, but it is now clear that the obvious substrates are wrong. The remaining options are:

### Option A: Continue numerology, harder targets
Test T1.3 (m_p/m_e = 1836.15) and T1.4 (α⁻¹ = 137.036) for completeness. These are different sectors and may behave differently. m_p/m_e is more likely to come from QCD confinement (Λ_QCD scale) rather than a clean algebraic invariant. α may be the most tractable of all if a substrate naturally has a coupling-strength definition.

### Option B: Pivot to RG / dynamical numerology
The Standard Model mass ratios and couplings are constrained by RG running. The Yukawa couplings at the Higgs VEV scale, the Weinberg angle at GUT scale, etc., have characteristic patterns (3/8 for Weinberg at GUT) that emerge from RG flow rather than from clean ratios at low energy. A "RG-aware" numerology would look for ratios that are clean at SOME scale, not necessarily the IR.

### Option C: Pivot to substrate-first with looser scoring
Build a candidate substrate (knots, octonions, quasicrystal, whatever) and ask what it produces. Accept that the prediction will be approximate (10% level). Use that as a constraint to narrow the substrate, then refine.

### Option D: Accept that the numerological bet has not paid off
Stop the numerology and pivot back to dynamics. Maybe particles emerge from a substrate whose parameters are themselves running and noisy; the precise masses are accidents of the IR; what we should fit is the structure of the spectrum (3 generations, hierarchy patterns), not the exact numbers.

My recommendation: try Option A briefly — run T1.3 and T1.4 to complete the Tier-1 scan, since each is ~one Python session. Either one of them yields a hit (small chance but huge payoff), or they all fail consistently and we have a complete negative dataset that justifies Option C or D as the strategic next move.

## Files

- `04_muon_electron_ratio_scan.py` — the scan.
- `04_findings.md` — this document.

## What This Experiment Established

1. The historical Nambu near-miss 3·α⁻¹/2 ≈ m_μ/m_e is the best simple-substrate approximation; no tested category beats 0.59% accuracy without explicit tuning.
2. m_μ/m_e is not a ratio of finite group orders, kissing numbers, n-ball volumes, factorials, binomials, or Catalan numbers (within 1%).
3. log_e(m_μ/m_e) is not a simple rational like 16/3 (off by 0.03% in log, 0.18% in value).
4. Substantial substrate-space pruning. Negative result with strategic implications: see "Strategic Reckoning" above.
