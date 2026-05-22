# Variant E — Findings: Quark Sector

**Date**: 2026-05-22
**Script**: `05_quark_sector.py`

---

## Result: Quarks do NOT satisfy a clean Koide identity

Computed Koide ratios for quark generations (PDG 2024 MS-bar masses):

- Up-type (u, c, t): **Q_up = 0.849** (gap from 2/3: +0.182)
- Down-type (d, s, b): **Q_down = 0.731** (gap from 2/3: +0.065)

Compare to charged leptons: Q = 0.6667 (within 6 × 10⁻⁶ of 2/3).

The quark Koide ratios are NOT at 2/3 within experimental precision.

## Brannen Phases Differ Across Sectors

| Sector | Brannen phase φ |
|--------|------------------|
| Charged leptons (e, μ, τ) | 0.2222 rad ≈ 2/9 |
| Up-type quarks (u, c, t) | −2.02 rad |
| Down-type quarks (d, s, b) | 0.110 rad |

The three sectors have **distinct** Brannen phases. The clean φ = 2/9 derivation in step 10 is specific to the lepton sector.

## Why This Is Expected

1. **Leptons are SU(3)_c singlets**; quarks are color triplets. Color mixing modifies the eigenvalue spectrum non-trivially.
2. **QCD running**: quark masses are quoted at different scales (MS-bar at 2 GeV for light quarks, at $m_c$ for c, $m_b$ for b, pole for t). A clean algebraic identity would need running to a common scale or use of running-invariant quantities.
3. **Higgs mixing**: quark Yukawa couplings include CKM mixing, which our minimal Furey framework hasn't yet incorporated.

## Implication

The lepton-sector Koide identity Q = 2/3 is a SPECIAL feature of color-singlet fermions. The framework's prediction of Koide structurally (via $\dim G_2/\dim \text{Spin}(7)$) **applies only to the lepton sector** in its current form. Extending to quarks requires explicit SU(3)_c content and running-mass treatment, neither of which is in the minimal Furey kernel as constructed.

This is not a failure — it correctly identifies that the lepton sector is the cleanly-derivable piece, and quark physics requires more machinery.

## Status

Variant E: **COMPLETE** with negative-for-quarks result. The structural derivation in v59 covers the lepton sector cleanly but does NOT extend to quarks without additional structure (color, running, CKM).

Proceed to Variant F (Lean formal verification).
