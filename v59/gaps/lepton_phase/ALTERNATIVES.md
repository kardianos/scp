# ALTERNATIVES — non-geometric origins for Q = 2/3 = 14/21 and cos(2/3)

**Gap G7.** The geometric route is DEAD (`README.md §3`: π-rationality, holonomy,
potential alignment, Weinberg cosine — all ruled out [thm]). The magnitude must come
from a **non-geometric** algebraic / spectral / extremal relation that produces the
transcendental `cos(3φ) = cos(2/3)`. This document opens the solution space: each
candidate states a concrete **test** and a **falsifier**. Numerical backing is in
`search.py` and `second_invariant.py`; clean identities are formalized in
`LeptonPhaseMagnitude.lean`.

Status legend: **[LIVE]** survives its test · **[NULL]** falsified here ·
**[OPEN]** untested target · **[MIXED]** partial.

---

## A. Maximal-mixing / G₂-content: `t² = (D − dim G₂)/D`  **[LIVE for the magnitude]**

**Idea.** The mass-splitting weight `t²` is the fraction of the sector outside the
`G₂ = Aut(𝕆)` automorphism core (which preserves the algebra and so carries no
generation-distinguishing information). For a `D`-dim sector with the `dim G₂ = 14`
core inert: `t² = (D − 14)/D`. This is a Lie-**dimension** ratio — **no angle**, so
it sidesteps the π-rationality obstruction entirely.

**Test** (`search.py §A`, `MaximalMixingKoide.lean` [thm]):

| sector | `D` | `t² = (D−14)/D` | `Q = (1+2t²)/3` |
|---|---|---|---|
| **lepton** | 28 = **2·14** | **1/2** (exact) | **2/3** |
| d-quark | 35 | 3/5 | 11/15 |
| u-quark | 63 | 7/9 | 23/27 |

The lepton is special **because `D_lepton = dim L = 28 = 2·dim G₂`**, making the
non-G₂ fraction *exactly* `1/2`. This gives `Q=2/3` non-geometrically.

**Falsifier** (passed): if another sector also hit exactly `1/2`, the lepton's
"specialness" would be coincidence — it does not (`3/5, 7/9 ≠ 1/2`).

**Honest limit** (`README.md §4`, retraction M): the *bridge* "maximally-mixed
L-vacuum non-G₂ weight = Brannen `t²`" does not close at the generation level
(max-mixing over the generation structure gives `t²∈{0,1}`). So this establishes
`Q=2/3` as the **so(8)-grade dimension ratio numerically equal to `t²=1/2`**, not a
derived dynamical amplitude. **Still the best magnitude candidate** because it
produces `2/3` from structure with no geometry.

**What it does NOT do:** explain the `φ=Q/3` *radian insert* — see §F.

---

## B. Equipartition / extremal characterization of `t²=1/2`  **[LIVE, new this run]**

**Idea.** A dynamics-free selection principle: the order parameter sits at the point
where the kernel's content is *balanced* between its symmetric and circulant-shift
parts — a genuine extremal / max-symmetry characterization (not an angle).

**Test** (`search.py §D2`). Write `M = a(I + ξS + ξ̄S²)`. The Frobenius² content
splits into a **diagonal/symmetric** part (the `I`, weight `1` per the normalization
`Σ|s_k|²/3a² = 1 + 2t²`) and a **circulant-shift** part (the `ξS + ξ̄S²` bands,
weight `2t²`). These are **equal** exactly when `2t² = 1`, i.e. **`t² = 1/2`**:

> `t²=1/2 ⟺ 2t²=1 ⟺` (off-diagonal/shift content) = (diagonal/identity content).

So `t²=1/2` is the **equipartition point** between the symmetric (no-mixing) and the
maximally-mixed circulant content of `M`. This is consistent with the "maximal
mixing" principle (P1) of `MaximalMixingKoide.lean`, but stated as a clean *algebraic
balance* with no reference to G₂.

**Formalized:** `LeptonPhaseMagnitude.equipartition_iff_half` [thm], and the bridge
`equipartition_gives_koide` ([thm]: `2t²=1 ⇒ Q=2/3`).

**Falsifier:** if the balanced point gave any other `t²` (e.g. the diagonal carried a
`t`-dependent weight), equipartition would not select `1/2`. It does. **[LIVE]**

**Relation to A:** A and B are two readings of the same `1/2` — A from the so(8)
dimension count (`28=2·14`), B from the kernel's symmetric/circulant balance. Both
are non-geometric. A explains *why 1/2* via dimensions; B characterizes 1/2 as an
extremum. Together they are the strongest non-geometric account of the **magnitude**.

---

## C. Casimir / dual-Coxeter ratios equal to 2/3  **[MIXED — twins, no mechanism]**

**Idea.** `Q=2/3` might be a ratio of quadratic Casimirs / dual Coxeter numbers `h∨`
of the v59 groups (a spectral, not geometric, invariant).

**Test** (`search.py §C`). Exact `2/3` ratios among `h∨` and dimensions:

- `h∨(G₂)/h∨(Spin8) = 4/6 = 2/3` (exact),
- `h∨(su2)/h∨(su3) = 2/3` (exact),
- `dim G₂/dim Spin7 = 14/21 = 2/3` (the Koide identity).

So there are **multiple** clean structural sources of the number `2/3`. But:

**Falsifier / verdict:** for a *mechanism* the Casimir ratio must be tied to the
lepton mass operator by representation theory. `h∨(su2)/h∨(su3)=2/3` involves groups
that are **not** the lepton sector's; `h∨(G₂)/h∨(Spin8)` is suggestive
(`G₂⊂Spin7⊂Spin8` is exactly the v59 chain) but no derivation links it to `t²`. So
these are **numerical twins**, evidence the value is structurally favored, **not** a
mechanism. **[MIXED]** — worth one targeted check (does the lepton mass bilinear's
Casimir genuinely carry `h∨(G₂)/h∨(Spin8)`?), else numerology.

---

## D. cos(2/3) as an algebra invariant / eigenvalue ratio  **[NULL]**

**Idea.** A spectral origin would have the *transcendental* invariant `cos(2/3)`
emerge directly as an eigenvalue or trace ratio of a real v59 operator.

**Test** (`search.py §B`, `second_invariant.py §3`). Scanned Lie-dimension ratios,
small-integer ratios, and `√3`-built quantities against `cos(2/3)=0.78589`. Best
rational near-miss: `11/14 = 0.78571` (`|Δ|=1.7×10⁻⁴`). Every algebraic candidate
(`Q`, `1−Q/3`, `7/9`, `√(2/3)`, `11/14`) **misses** at the data precision (`10⁻⁵`);
only `cos(Q-as-radian)` fits.

**Falsifier / verdict:** `cos(2/3)` is **transcendental**; no rational/Lie ratio can
equal it, and a `10⁻³`-level near-miss is meaningless at `10⁻⁵` data precision.
**[NULL]** — there is no non-geometric *algebraic shortcut* to `cos(2/3)`; a real
mechanism must produce the transcendental number itself (see §E, §F).

---

## E. The skewness / "second-Koide" reading  **[LIVE reframing — sharpest target]**

**Idea (new this run, `second_invariant.py`).** The phase invariant is a pure data
moment, no angle. Standardize the amplitudes `x_k = (√m_k − a)/a = 2t·cos(φ+2πk/3)`.
Then:

- the **2nd** standardized moment `⟨x²⟩ = 2t²` is fixed by **Koide** (`⇒ Q`);
- the **3rd** standardized moment is `⟨x³⟩ = (2t)³·(1/4)·cos 3φ` — so

> **`cos(3φ) = 4·⟨x³⟩/(2t)³`** — the phase invariant **is** the (rescaled)
> **skewness** of the `√m` distribution.

Verified numerically: from raw PDG `√m` moments, `cos(3φ) = 0.78587 = cos(2/3)` to
`1.4×10⁻⁵`. **Formalized:** `LeptonPhaseMagnitude.cos3phi_is_third_moment` [thm].

**This reframes G7 as a "second Koide":** Koide fixes the *variance* of `√m` (`⇒Q`);
G7 would be a structural relation fixing the *skewness* of `√m` to exactly the value
making `cos(3φ) = cos(Q-as-radian)`. The mystery, stated sharply: *why does the
data's skewness land precisely where the cosine of the variance-ratio (in radians)
sits?*

**Test / target [OPEN]:** find a single structural functional `F(√m_0,√m_1,√m_2)`
whose vanishing (`F=0`, a constraint like Koide's `Σm = (2/3)(Σ√m)²`) forces both
the variance to `Q=2/3` **and** the skewness to `cos(2/3)`. Koide is the
variance-constraint; the unknown is the skewness-constraint.

**Falsifier:** if the data's skewness were equally well fit by a clean *algebraic*
value (e.g. `7/9` or `11/14`), there would be no transcendental to derive and §F
would be a non-problem. It is **not** — only `cos(Q)` fits (`§D`). So the
transcendental target is real. **[LIVE reframing; OPEN derivation]**

---

## F. The hard core: why does the amplitude number `Q` re-enter as a radian?  **[OPEN]**

**The residual that even A+B+E do not touch.** Granting `Q=2/3` from the dimension
ratio (A) / equipartition (B), the law `φ=Q/3` *inserts that pure dimensionless
number as an angle in radians* (`cos(2/3)`). Why would the **amplitude**'s value
re-appear as the **phase**'s value?

**Tests run** (`search.py §G`, `second_invariant.py §3`): 

- **Self-referential fixed point** `x = cos(3x)`: gives `x=0.390 ≠ 2/9`. **[NULL]**
- **Kernel self-consistency** `3φ = Q(φ)`: impossible, since `Q` is
  phase-independent (`Q_phase_independent` [thm]) — the kernel cannot dynamically
  force `3φ=Q`. **[NULL]**
- **Algebraic disguise** `cos 3φ = f(Q)` for clean non-radian `f`: every candidate
  misses; only `cos(Q)` fits (`§D`). So the radian insert is **real**, not a
  disguised rational. **[confirmed]**

**What a derivation needs (the OPEN statement of G7):** a single structural object
that simultaneously (i) carries the dimensionless ratio `Q=2/3` (from the
G₂-content, §A) **and** (ii) feeds it as the *argument of a cosine* (i.e. produces
`cos(2/3)` as a genuine spectral/analytic quantity, not a fitted value). Candidate
shapes to look for:

- a **single eigenvalue** of a v59 operator of the form `λ ∝ cos(c)` where the
  *same* dimension count fixes both the prefactor and `c = Q`; 
- a **partition-function / character** evaluation `χ(g)` at the group element `g=ψ`
  (the order-3 sedenion automorphism) whose value involves `cos` of a
  G₂-content-weighted angle — tying the explained `/3` (ψ) to the residual `Q` in
  one expression. This is the only route that could unify §A (magnitude) with the
  structural `/3` and produce a transcendental — **untested, the recommended next
  probe** (it must avoid the π-rational trap: `ψ`'s *own* eigenphases are π-rational,
  so the `Q` must enter as a *weight/Casimir*, not as `ψ`'s rotation angle).

**Falsifier:** if no operator/character of the v59 algebra evaluates to `cos(2/3)`
with a G₂-content reason for the argument, G7's magnitude remains a precise but
mechanism-less empirical regularity (the honest default, like Koide for decades).

---

## G. The √3 of the sedenion ψ as a magnitude source  **[NULL]**

**Idea.** `ψ` (the order-3 generation automorphism) carries `√3 = 2sin60°` in its
120° block `[[−½, √3/2],[−√3/2, −½]]`. Could a `√3`-built algebra invariant equal
`2/3` or `cos(2/3)`?

**Test** (`search.py §E`): `1/√3, √3/3, √3−1, 2−√3, 1−1/√3` — none equals `2/3` or
`cos(2/3)`.

**Falsifier / verdict:** the `√3` is *π-rational structure* (it builds the 120°
rotation = the `/3`); its eigenphases are `{0,±2π/3}`. It explains the **rotation**
(the divisor 3), **not** the magnitude. **[NULL]** — consistent with the structural
vs residual split: `ψ` gives `/3`, never `Q`.

---

## Summary table

| § | avenue | target | status |
|---|---|---|---|
| A | max-mixing / G₂-content `t²=(D−14)/D` | magnitude `Q=2/3` | **LIVE** (so(8)-grade ratio; lepton `28=2·14`) |
| B | equipartition `2t²=1` (kernel balance) | magnitude `t²=1/2` | **LIVE** (extremal, non-geometric) |
| C | Casimir / `h∨` ratios = 2/3 | magnitude | MIXED (twins, no mechanism) |
| D | `cos(2/3)` as algebra invariant | invariant | **NULL** (transcendental; near-misses are numerology) |
| E | skewness / "second Koide" | invariant `cos(2/3)` | **LIVE reframing**, derivation OPEN |
| F | why `Q` re-enters as a radian | the `φ=Q/3` insert | **OPEN** (the hard core; character/eigenvalue probe recommended) |
| G | the `√3` of ψ | magnitude | **NULL** (explains `/3`, not `Q`) |

**Most promising:** A+B (magnitude is non-geometrically real) feeding E (the
skewness reframing makes `cos(2/3)` a pure data moment). The unsolved knot is F.
