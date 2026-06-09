# v59→v63 Predictivity Audit — the decision number

**Date**: 2026-05-28
**Question**: is the octonionic structure **over-determined** (few inputs → many
independent, look-elsewhere-robust SM relations = predictive, worth continuing),
or is relations ≈ inputs (a sophisticated fit = stop)?
**Answer**: naive ratio **4.5×** (looks predictive); **honest ratio ≈ 0.36×** (on
the **fit** side). One genuinely striking survivor: the lepton Koide relation.

---

## Method

Codify the project's own audits (`v61/PROVEN_LEDGER.md`,
`v59/furey_construction/RIGOR_AUDIT.md`) into one number, applying three deflations
a naive match-count ignores:

- **D1 — SM-tree redundancy**: `m_W`, `m_Z` follow from `(v, g_W, θ_W)` by SM
  definition → not independent reproductions.
- **D2 — look-elsewhere**: `α` ("5 of 8 arbitrary templates hit it to 0.03%") and
  quark Koide ("12 simple rationals per RG band") carry ≈ zero evidential weight.
- **D3 — load-bearing unproven ansätze**: `v_Higgs=784a²` ("looks fitted"),
  `g_W²=5√α` ("α in disguise"), the `Z₂×Z₂` selection rule ("undriven") are
  **assumptions** — they belong on the input side.

## Result (`predictivity_audit.py`, 5/5 checks)

| | count |
|---|---|
| naive reproductions | 9 |
| naive inputs (`a_ℓ + α`) | 2 |
| **naive ratio** | **4.5×** (reads as over-determined) |
| honest independent + LE-robust relation content (Σ weights) | **1.8** |
| honest assumptions (inputs + load-bearing ansätze) | **5** |
| **honest predictivity ratio** | **0.36×** ( < 1 ) |

The naive 4.5× is an artifact of (D1) double-counting SM-tree quantities and (D2)
counting look-elsewhere-weak matches as strong. After honest deflation there are
**at least as many undriven assumptions as look-elsewhere-robust relations** — the
structure is **not** over-determined. This matches the project's own bottom line:
*"a set of striking-but-unexplained numerical coincidences anchored on two genuine
inputs."*

## The one survivor

**Lepton Koide `Q = 2/3 = dim G₂/dim Spin(7)` at `~10⁻⁵`** — a 5-digit match of a
40-year-unexplained relation to a *forced* group-dimension ratio (weight 0.80). It
does not deflate away. The honest caveat: it is conditional on the undriven
selection rule (`D_lepton=28 ⇒ t²=1/2 ⇒ Q=2/3`), so even it is "reproduces Koide
given one undriven assignment." But it is the strongest, hardest-to-dismiss piece —
the real residual value of the octonionic program.

## Verdict (against the agreed criterion)

By "ratio ≫ 1 ⇒ predictive; ≈ 1 ⇒ fit," the octonionic *parameter-derivation*
program sits on the **fit** side (0.36×). Combined with the v62 no-go (the
transcendentals `α`, `cos(2/3)` are provably unreachable algebraically), this is a
clean, defensible close: the program reproduces the lepton Koide strikingly and is
otherwise a fit anchored on `a_ℓ + α` plus a stack of unproven ansätze.

**Note on switching algebras**: a different algebraic structure (E8, h₃(𝕆),
Clifford) does **not** help — it inherits the same no-go *and* offers *more*
invariants, i.e. *more* look-elsewhere freedom (easier to fit, not harder). The
only genuinely different tracks are **dynamical/RG** (where a transcendental can be
a fixed-point output) or **reframing away from constant-derivation** (the original
soliton/simulation program, whose success criterion — does a stable particle with
the right quantum numbers emerge — is not a numerology problem).

**Artifact**: `predictivity_audit.py` (5/5 checks).
