# v62 — The Number-Type Map of the SCP Residuals

**Date**: 2026-05-28
**Supersedes**: `speculative/OCTO_SPECTRAL_PROJECTION_MENTAL_MODEL.md` as the
working map of the problem space. That document was an honest but spineless
enumeration (seven readings of "projection" as a flat menu, "07 with better J's"
as the top next step). This is the spine.

**One sentence**: every quantity in the SCP program is either *algebraic*
(reachable by representation theory) or *transcendental* (not), and that single
distinction tells you which "residual conjectures" are theorems-in-waiting, which
are closable by state-selection, and which are barred from the algebra forever.

---

## 1. The discriminant: algebraic vs transcendental

An algebraic construction — branching/restriction of representations, Casimirs,
group characters, dimension counts, eigenvalues of matrices with algebraic
entries — can only output **algebraic numbers** (roots of `ℚ`-polynomials; closed
under `+ − × ÷`). Therefore:

> A **transcendental** target is unreachable by any algebraic construction.

Sorting the program's quantities:

| algebraic (rep theory CAN produce) | transcendental (rep theory CANNOT) |
|---|---|
| `784 = 28²`, `dim G₂ = 14`, `dim Spin(7) = 21` | the Brannen phase invariant `cos(2/3)` |
| Koide `Q = 2/3`, `sin²θ_W = 2/9`, gauge `5/2/9` | `cos(2/9)` (cosine of the phase itself) |
| the `(1−t²)D = 14` cross-sector identity | `α` (fine-structure constant) |
| selection rule, ±2 helicity content | `f_g ~ α^{21/2}` (a function of `α`) |

The program's *successes* are exactly the algebraic column (rigid, Burnside- or
Casimir-forced). Its chronic *frustrations* are exactly the transcendental column
— it kept trying to extract them from an algebraic machine. **The frustration was
a category error, not a lack of cleverness.** (Proof for the phase: `no_go/NOGO.md`.)

---

## 2. "Projection" is four operations, not one

The speculative document's options A–G are not a menu; they are a pipeline, and
each stage handles a fixed number-type:

| op | operation | produces | number-type |
|---|---|---|---|
| **P1** | equivariant restriction / branching (G₂-covariant, Burnside) | dimensions, irrep content, selection rule | algebraic |
| **P2** | Hessian at a selected vacuum | the rank-3 Brannen mass *matrix* (2nd variation of the democratic vacuum) | algebraic structure; magnitudes set by a scale |
| **P3** | coarse-graining / integrating out (RG; the Shulga kernel's real job) | running couplings, the loop-ratio that pins the phase | **transcendental** (logs, loop integrals) |
| **P4** | soldering → spacetime (Cl(3,1) = commutant of Spin(7)) | metric as tetrad², ±2 helicity, GR | structural |

P1/P2/P4-structure is where rep theory wins. **P3 is the only home for the
transcendentals.** Forcing a P3 quantity (the phase, `α`) through P1/P2 machinery
(characters, `J`-eigenphases, Shulga modulation) is the barred move.

The rank tension (`v60/gaps/rank_tension`) dissolves here: the democratic `784`
vacuum is a P2 *vacuum* (full symmetry), the rank-3 Brannen kernel is its P2
*Hessian* (whatever rank the breaking dictates). Standard SSB, not a contradiction.

---

## 3. The flat-direction law

> Representation theory pins the symmetric skeleton of each sector and is
> constitutionally **blind to the flat (angular) directions** it leaves behind.
> Every residual conjecture sits on such a flat direction.

Demonstrated on two independent sectors (`flat_direction/flat_direction_demo.py`):

- **Brannen phase**: every `Z₃`-*symmetric* invariant (`Σ√m`, `e₂`, `Σm`, Koide
  `Q`) is `φ`-**independent**; `φ` enters only the `cos(3φ)`-class invariants.
  The symmetric structure fixes the whole skeleton (and `Q = 2/3` at `t² = 1/2`)
  and leaves `φ` as a flat coordinate it cannot see.
- **EW democracy (R1)**: the O(784) Frobenius hat fixes only the radius (1 Higgs);
  the `783` angular directions are exact Goldstones (`S^783`,
  `v61/lean/EwVevHome.lean`). "Democracy" is one point on this flat sphere.

Same shape both times: **magnitude fixed, angle flat.** The flat part is exactly
what the algebra cannot fix — which is *why* the residuals are residual.

---

## 4. The residual partition (and the input count)

v61's closeout lumps four residuals as "value/symmetry inputs, all on the same
footing as `α`." That is too coarse. By number-type and flatness they split into
**three kinds with three different fates**:

| residual | kind | fate |
|---|---|---|
| the scale `a` (`a = √v / dim L`) | **pure scale** | unavoidable units anchor — *not* a defect; stop counting it |
| EW democracy in `v = 784a²` | **degeneracy** (flat angle, rational point) | **closable** by state-selection (below) |
| Brannen phase `φ = 2/9` | **transcendental value** on a flat angle | P3 loop-ratio, or input — *barred from algebra* |
| `α`; `f_g ~ α^{21/2}` | **transcendental coupling** | input (`f_g` is downstream of `α`, see §6) |

So three of the four have real exits; only the transcendental *couplings* are
genuinely hard. The `784` itself is algebraic and Burnside-forced — a P1 win, not
a residual; only its *alignment* (democracy) was open.

**Democracy is closable** (`residual_audit/democracy_selection.py`): the
democratic vacuum is selected by *any* of three principled, non-transcendental
criteria — (1) it is the unique point on `S^783` fixed by the generation
permutation `Sₙ`; (2) it uniquely maximizes the mixing entropy; (3) the omitted
`Sₙ`-invariant quartic lifts all `n−1` Goldstones with it as the strict minimum.
Closing it removes `v = 784a²` from the residual list, leaving only the scale.

---

## 5. Why the phase is different from democracy (the crux)

Both are flat directions, but the flat *coordinate* differs in number-type:

- Democracy's flat coordinate ranges over a sphere of **rational/algebraic**
  points; selecting one is a *degeneracy* problem → state-selection solves it.
- The phase's flat coordinate has a **transcendental value** (`cos(2/3)`);
  *no* state-selection or representation move reaches a transcendental. Only a
  dynamical loop ratio `−c₃/(4c₆)` (P3) can equal it
  (`residual_audit/phase_dynamical_home.py`).

This is the whole point of the partition: same geometry (flat direction), opposite
resolution, dictated by number-type.

---

## 6. Input-count corrections to v61

- **`f_g ~ α^{21/2}` is probably not independent.** `21 = dim Spin(7)`, so the
  exponent is `dim Spin(7)/2`: algebraic structure on a transcendental base. If
  `21/2` is *derived*, `f_g` is a function of `α`, not a separate input — the
  closeout over-counts by one. (Open: confirm the exponent is forced, not fit.)
- **After democracy selection, `v = 784a²` is not an independent input** — it
  reduces to the single scale `a`.

Net genuine inputs trend toward: **one scale `a`** + **the transcendental
coupling(s) `α`** (with `φ` and `f_g` downstream of dynamics/`α` if the open P3
and exponent tasks close). This is the metric that matters — *input count* — not
tolerance-hitting on any one number.

---

## 7. Falsifiability: one source, many outputs, zero new knobs

The reframing earns its keep only if it **reduces independent inputs**. The test
is not "did we hit `2/9` within 5%" but:

- **Win**: one distinguished state + one equivariant restriction (P1/P2) yields
  `{784, 28/3, 14, 5/2/9, sin²θ_W = 2/9, Q = 2/3, Z₃ spacing}` with **zero new
  knobs**. (Most are already claimed structural; the open question is whether
  they follow from *one* state choice or several.)
- **Honest segregation**: `φ`, `α`, `f_g` are the transcendental sector — derived
  via P3 or held as inputs, never pretended to be algebra outputs.

---

## 8. What v62 proves, and the resulting effort allocation

Proven here (numerically self-checking + machine-checked Lean):

1. **The phase no-go** (`no_go/`, `lean/PhaseNoGo.lean`): `cos(2/3)`, `cos(2/9)`
   transcendental ⟹ phase unreachable by any algebraic map; cuts 01–06's fixed
   phase explained; the "07 cut" barred.
2. **The flat-direction law** (`flat_direction/`): phase and democracy are the
   same phenomenon; symmetric invariants are flat in the residual coordinate.
3. **Democracy is closable** (`residual_audit/democracy_selection.py`): three
   independent selectors pick the democratic vacuum.
4. **The phase's only home is P3** (`residual_audit/phase_dynamical_home.py`): the
   loop-ratio is the unique number-type that can equal `cos(2/3)`.

**Effort allocation (the actionable output):**

- **STOP**: making the algebra emit the phase (`J`'s, characters, Shulga kernels);
  treating the scale `a` as a defect. *Barred / not-a-defect.*
- **CLOSE (tractable)**: democracy selection (state-selection, §4); the soldering
  action principle for P4 (`v60` G9, posited not derived — a find-the-principle
  gap, not a barred one); confirm the `f_g` exponent `21/2` is forced.
- **HARD (transcendental)**: the P3 loop Lagrangian for the phase; `α` as input.

The verdict follows from number-type, not taste.
