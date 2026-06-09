# FINDINGS — Gap G7: a live non-geometric origin for Q = 2/3?

**Date:** 2026-05-25 · **Folder:** `v59/gaps/lepton_phase/` · **Gap:** G7, the
origin of the lepton mass-phase magnitude `Q = 2/3` (phase `φ = Q/3 = 2/9`).

---

## Verdict

**Split the question, as the gap itself splits:**

1. **The MAGNITUDE `Q = 2/3` HAS a live non-geometric origin.** Two consistent
   non-geometric readings of `t²=1/2`:
   - **(A) Maximal-mixing / G₂-content** [thm arithmetic]: `t² = (D−dim G₂)/D`,
     a Lie-**dimension** ratio (no angle, immune to the π-rationality death). The
     lepton's `1/2` is *forced* by `dim L = 28 = 2·dim G₂`; the quark sectors give
     `3/5, 7/9` (so `1/2` is special, not generic).
   - **(B) Equipartition** [thm, new this run]: `t²=1/2` is exactly the point where
     the Brannen kernel's identity content (weight 1) balances its circulant-shift
     content (weight `2t²`): `2t²=1 ⟺ t²=1/2`. A clean extremal characterization.

   Both produce `Q=2/3` from structure with **no geometry**. This is the most
   promising avenue and it survives every test run.

2. **The phase invariant `cos(2/3)` has NO clean algebraic shortcut** [NULL]. It is
   transcendental; no Lie-dimension ratio, Casimir ratio, or `√3`-built quantity
   equals it (best near-miss `11/14`, off by `1.7×10⁻⁴ ≫` the `10⁻⁵` data precision).
   Crucially, the data does **not** prefer any algebraic `f(Q)` for `cos 3φ` — only
   the transcendental `cos(Q-as-radian)` fits at `10⁻⁵` (`second_invariant.py §3`).
   So the magnitude problem cannot be dodged by finding a "nicer" target.

3. **The HARD CORE (gap F) remains OPEN.** Even granting `Q=2/3` from A/B, the law
   `φ=Q/3` *inserts that pure dimensionless number as a radian* (`cos(2/3)`). No
   non-geometric mechanism explains why the **amplitude** number `Q` re-appears as
   the **phase** argument. Self-referential fixed points (`x=cos 3x`) and kernel
   self-consistency (`3φ=Q`) both fail — the latter provably, since `Q` is
   phase-independent (`BrannenPhase.Q_phase_independent` [thm]).

**Bottom line:** the magnitude `Q=2/3` is non-geometrically real and well-founded
(dimension ratio = equipartition point); the **`φ=Q/3` radian insert is the genuine
residual**, with no live mechanism and no algebraic shortcut. G7's open content is
sharper than before: not "explain `2/3`" but "explain why the amplitude ratio `Q`
re-enters as the cosine argument, producing the transcendental `cos(2/3)`."

---

## New results produced this run

| result | where | status |
|---|---|---|
| **Equipartition characterization** of `t²=1/2`: kernel identity vs circulant-shift content balance | `search.py §D2`, `LeptonPhaseMagnitude.equipartition_iff_half`, `equipartition_gives_koide` | [thm, new] — a non-geometric extremal reading of the magnitude |
| **Skewness reframing**: `cos 3φ = 4·⟨x³⟩/(2t)³` — the phase invariant IS the standardized third moment of `√m` | `second_invariant.py §1`, `LeptonPhaseMagnitude.cos3phi_is_third_moment` | [thm, new] — recasts `cos(2/3)` as a pure data moment, no angle |
| **No algebraic `f(Q)`** fits `cos 3φ` (only `cos Q` does) — the radian insert is real, not a disguise | `second_invariant.py §3` | [emp, new] — confirms the transcendental target is unavoidable |
| **Casimir twins**: `h∨(G₂)/h∨(Spin8) = h∨(su2)/h∨(su3) = 2/3` exactly | `search.py §C` | [obs] — numerical twins, no mechanism (MIXED) |
| **`√3` of ψ does not source the magnitude** | `search.py §E` | [NULL] — `ψ` gives the `/3`, never `Q` |
| **Sensitivity**: `2/3` is the unique low-complexity fraction at `10⁻⁵` data precision | `search.py §F` | [emp] — confirms there is something real to derive |

These complement the prior negative results (geometric route dead: π-rationality,
holonomy, naive potential, Weinberg cosine — all [thm] in `furey_construction/lean/`).

---

## Most promising avenue and what a derivation needs

**Most promising:** avenues **A+B** (the magnitude is a non-geometric dimension
ratio = equipartition point) feeding the **skewness reframing (E)** (the invariant is
the third moment of `√m`). This reduces G7 to a "second-Koide" problem:

> Koide is the structural constraint fixing the **variance** of `√m` (`⇒ Q=2/3`). G7
> needs a structural constraint fixing the **skewness** of `√m` to exactly the value
> making `cos 3φ = cos(Q-as-radian) = cos(2/3)`.

**A derivation of the hard core (F) needs** a single v59-structural object that:

1. carries the dimensionless ratio `Q=2/3` (from the G₂-content, A), **and**
2. feeds it as the **argument of a cosine** — producing the transcendental
   `cos(2/3)` as a genuine spectral/analytic quantity, not a fitted value.

The recommended **untested probe** (avoiding the π-rational trap): a **character /
partition-function** evaluation `χ(ψ)` of the order-3 sedenion automorphism `ψ`
(which carries the structural `/3`), **weighted by the G₂-content Casimir** so the
`Q` enters as a *weight*, not as `ψ`'s own (π-rational) rotation angle. If some such
`χ(ψ)` evaluates to `cos(2/3)` with a G₂-content reason for the argument, G7 closes.
If no operator/character of the v59 algebra does, G7's magnitude is — like Koide for
decades — a precise (`~10⁻⁵`) but mechanism-less empirical regularity. This is stated
as the open Lean theorem `LeptonPhaseMagnitude.magnitude_to_phase_open` (trivial
witness, mechanism flagged OPEN).

---

## Honest scope / caveats

- **Lean status:** `LeptonPhaseMagnitude.lean` is **written, not built this run**
  (the shared `furey_construction/lean` project must not be co-built; concurrent
  `lake build` conflicts). The arithmetic theorems (A, B) and the trig identities
  (the cycle/cube sums, E) mirror already-proven analogues in `BrannenKernel.lean` /
  `BrannenPhase.lean` and should check; the open core is an explicit placeholder.
- **The bridge in A does not fully close** (retraction in `MaximalMixingKoide.lean`):
  max-mixing over the *generation* structure leaves `t²` free; `(D−dimG₂)/D` is the
  so(8)-grade dimension ratio that *numerically equals* `t²=1/2`, not a derived
  dynamical amplitude. So even the magnitude is "structural number matched to data,"
  not "dynamically derived" — consistent with the v59 meta-rule (`UNIFIED_THEORY.md`).
- **The geometric route stays dead** — none of this resurrects it; A, B, E are
  dimension/extremal/moment statements, not angles.
