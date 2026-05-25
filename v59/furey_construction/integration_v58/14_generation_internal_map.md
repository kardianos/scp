# The generation↔internal map: PARTIAL — the ½ reaches the amplitude, not the phase

*2026-05-24.  The make-or-break attempt: carry the L-grade complex-structure root-product `½`
(`13_octomath_half.md`, `OctoHalf.lean`) to the generation Brannen amplitude `t²`.  **Result: the
`½` IS carried to the amplitude** — `t²=½` is the reversion-norm of the L-grade complex
half-element — **but NOT to the phase**, and the amplitude grounding still rests on the
half-norm (grade-balance) assumption.  Honest partial bridge, with a clean new L/F distinction.*

## What carries (the amplitude)

Identify the Yukawa coupling `ξ` (the off-diagonal of `M=a(I+ξS+ξ̄S²)`) with the **L-grade complex
half-element** `(1+u)/2`, `u²=−1`:

> **reversion-norm** `ξξ̃ = ((1+u)/2)((1−u)/2) = (1−u²)/4 = ½` — a **scalar**, `= |ξ|² = t²`.

So `t² = ½ = ` Koide amplitude.  The `½` is exactly the OctoHalf complex-structure norm, and it is
carried to the generation amplitude.  Crucially, the **L/F distinction is clean**:

| coupling grade | `u²` | `ξξ̃` (reversion-norm) |
|---|---|---|
| **L (complex)** | `−1` | **`½` (a scalar — a genuine norm)** |
| **F (real)** | `+1` | `(1+w)/2` — **not a scalar** (idempotent; no norm) |

The L-grade (`u²=−1`) gives a genuine *scalar* reversion-norm `½`; the F-grade (`u²=+1`) doesn't
give a scalar norm at all.  Since **`mass ∈ L` is proven** (`LeptonRealityForcing`), the lepton
coupling carries reversion-norm **`½`**, hence `t²=½`.  *This grounds the amplitude `½` in the
proven L-grade complex structure* — the strongest anchor yet.

## What does NOT carry (the phase) — and a hard obstruction

- **The phase is FREE — and `45°` is the amplitude, not a phase prediction** *(corrected 2026-05-24
  after the "is 45° a missing phase shift?" check)*.  The half-element's `arg(½+½u)=45°` is the angle
  where `|Re ξ|=|Im ξ|` — the **internal balance** that *makes* `|ξ|²=¼+¼=½`.  So **`45°` encodes the
  amplitude `t²=½`, not the phase.**  Since `|ξ|²=½` holds for *any* `arg(ξ)`, the physical coupling
  is a **free rotation of the half-element on the `|ξ|²=½` circle**: `ξ_phys=(1/√2)e^{uφ}`.  The
  "missing phase shift" `φ−π/4` is exactly that free rotation — **the phase residual**, not a
  conflict.  Empirically the physical phase `φ=2/9 ≈ 12.7°` sits **near the electron-massless / chiral
  point** (`M4`; the electron is `~2.3°` short of `m_e=0`), **far from `45°`**, and the shift
  `π/4−2/9 ≈ 0.563 rad` shows no clean structure.  So: **amplitude `½` = the `45°` balance (carried);
  phase = the free shift off it (residual, near the chiral limit)** — `φ` not symmetry-fixed, not
  π-rational.
- **The generation off-diagonal `ζ=ξS+ξ̄S²` is *never* a complex structure**: `⟨ζ²⟩₀ = +2t² ≥ 0 ≠ −1`
  (verified at both `φ=π/4` and `φ=2/9`).  So the carrier is the **coupling `ξ`** (the half-element),
  *not* `ζ`.  This is a hard obstruction (the scalar part of `ζ²` is non-negative), robust to
  octonionic structure.

## The honest residual

The amplitude bridge reduces `t²=½` to `|ξ|² = ξξ̃ = ½`, i.e. **`ξ` is a half-norm L-grade complex
element**.  This norm `½` *is* the OctoHalf complex-structure signature (and the F-grade gives no
scalar norm) — so it is now anchored in the proven `mass∈L`.  But `|ξ|²=½` holds for *any* phase,
so the condition is still "`ξ` has reversion-norm `½`" = the **grade-balance**, now *carried by* a
concrete octonionic object (the L-grade complex half-element) rather than imposed by hand.  It is
**grounded, not a clean from-nothing forcing**: "why the coupling is half-normed" is the same
grade-balance, though its *value* `½` (vs the F-grade's non-norm) is now structural.

## Verdict (build it / show it can't)

> **PARTIAL BUILD.**  The `½` *does* reach the generation amplitude: `t² = ξξ̃ = ½` is the
> reversion-norm of the L-grade complex half-element, well-defined as a scalar *precisely because*
> `u²=−1` (L-grade, proven `mass∈L`), with the F-grade giving no scalar norm.  This is the most
> grounded the amplitude `½` gets — anchored in proven structure via the OctoHalf identity.
>
> **It does NOT close to a clean forcing:** (i) the amplitude still rests on "`ξ` is half-normed"
> (the grade-balance, now carried by the half-element, not derived from a deeper principle); (ii)
> the **phase `φ=2/9` is a free rotation** on the `|ξ|²=½` circle — the `45°` of the half-element is
> the *amplitude* balance, not a phase, so the phase is genuinely free (residual, near the chiral
> limit); (iii) the generation `ζ` provably cannot itself be the complex structure.

So: the amplitude side of the bridge is **built and grounded** (the Koide `½` = L-grade complex
reversion-norm); the phase is a **separate, still-open residual**; and a clean *forcing* of the
amplitude (beyond the grade-balance) is **not** achieved.  Net: `t²=½` is now "the L-grade complex
half-norm" — structural and proven-anchored — rather than a bare number, but `φ=2/9` and the
ultimate "why half-norm" remain open.

## Files
- This attempt: inline computation (octonion/Cl C-subalgebra + 3×3 generation).
- Foundations: `13_octomath_half.md`, `lean/OctoHalf.lean` (the `½` = complex-structure root-product).
