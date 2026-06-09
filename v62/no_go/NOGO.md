# The Brannen-Phase Transcendentality No-Go

**Claim.** The lepton Brannen phase `φ = 2/9` cannot be produced as an output of
*any* algebraic construction built from the octonionic/Cl(7) structure — no
Casimir eigenvalue, no group character, no dimension count, no eigenphase of a
matrix with algebraic entries (toy `J`, Furey-pinned complex structure, or
Shulga-modulated rotation). It is either a **dynamical** (loop-ratio, P3) output
or an **input**. This closes the spectral-projection program *for the phase*.

This is the theorem-grade core of v62. v59 already *suspected* it
(`furey_construction/lean/LeptonPhaseEmpirical.lean`: "`cos 3φ = cos(2/3)` is
transcendental, non-geometric, not the Weinberg angle"; `CLOSEOUT.md` N5:
"`cos(2/3)` is not an algebra invariant … *transcendental*"). v62 turns that
remark into a proof and follows it to the consequence the mental-model document
(`speculative/OCTO_SPECTRAL_PROJECTION_MENTAL_MODEL.md`) missed: the program's
own top-ranked next steps (option D "characters/derivations", the "07 cut with
better J's") are *barred a priori*.

---

## 1. The number-type discriminant

Partition every target quantity by the kind of number it is.

- **Algebraic** — a root of a polynomial with rational coefficients. This class
  is closed under `+, −, ×, ÷` and contains every quantity an algebraic
  construction emits: dimensions (`784`, `14`, `28`), rationals (`2/3`, `2/9`,
  `5/2/9`), Casimir eigenvalues, structure constants, roots of unity, and the
  eigenvalues of any matrix whose entries are algebraic.
- **Transcendental** — not algebraic. Includes `cos q`, `sin q`, `e^q`, `log q`,
  `π` for nonzero algebraic `q`, and generic ratios of loop integrals.

An algebraic machine can only emit algebraic numbers. So **a transcendental
target is unreachable by any algebraic construction, full stop.**

---

## 2. The phase invariant is transcendental

The Brannen mass kernel uses `√mₖ = a(1 + √2·cos(2πk/3 + φ))`; the only
`φ`-dependence of any `Z₃`-invariant of the spectrum enters through `cos(3φ)`
(proven symbolically in `../flat_direction/flat_direction_demo.py`). For the
leptons `φ = Q/3 = 2/9`, so the invariant is

```
    cos(3φ) = cos(2/3) ≈ 0.7858872607769…
```

**Lemma (Lindemann–Weierstrass corollary).** *For every nonzero rational `q`,
`cos q` is transcendental.*

*Proof.* Lindemann–Weierstrass: if `α ≠ 0` is algebraic then `e^α` is
transcendental. Take `α = iq` (a nonzero algebraic, since `i` and `q` are
algebraic). Then `e^{iq}` is transcendental. Suppose `cos q` were algebraic.
Then `e^{iq}` satisfies

```
    z² − 2(cos q)·z + 1 = 0          (because e^{iq}·e^{-iq} = 1 and
                                      e^{iq} + e^{-iq} = 2 cos q),
```

a polynomial with algebraic coefficients, making `e^{iq}` algebraic — a
contradiction. ∎

Applying the lemma with `q = 2/3` and `q = 2/9`:

> **`cos(2/3)` and `cos(2/9)` are transcendental.**

(Formalized in `../lean/PhaseNoGo.lean` as `cos_two_thirds_transcendental` /
`cos_two_ninths_transcendental`. Mathlib currently ships only the *analytical
part* of L–W, so that one step is a cited `sorry`; everything downstream is
machine-checked. Numerical corroboration — `cos(2/3)` admits no integer relation
through degree 8, height 10⁶ — is in `transcendence_nogo.py` Part 2.)

This also refutes the specific algebraic identity v59 hoped for: the would-be
"deep identity" `cos 3φ = cos²θ_W = 7/9` is **false**,
`cos(2/3) − 7/9 = +8.1×10⁻³ > 0` (`AxiomCheck.lean:86`; reproduced in
`transcendence_nogo.py` Part 1).

---

## 3. Two corollaries that close the program

**Corollary A (no algebraic invariant).** No Casimir eigenvalue, dimension,
character value, or structure constant equals `cos(2/3)`, because all of these
are algebraic and `cos(2/3)` is not. *(`no_algebraic_eq_phase_invariant`.)*

**Corollary B (no algebraic eigenphase — the `J` no-go).** The v59/01–06 phase
machinery realizes the phase as the eigenphase of a rotation `R = exp(φ·J)`. A
matrix with algebraic entries has an algebraic characteristic polynomial, hence
**algebraic eigenvalues**; if `e^{iθ}` is an eigenvalue then
`cos θ = (e^{iθ} + e^{-iθ})/2` is algebraic. Since `cos(2/9)` is transcendental,
**no matrix assembled from the (algebraic) octonion structure can have eigenphase
`2/9`** — not the toy `J`'s, not "proper Furey-pinned" `J`'s, not any
Shulga-modulated combination. *(`no_algebraic_eigenphase_two_ninths`.)*

The phase rotation by `2/9` *exists* only when built by inserting the
transcendental numbers `cos(2/9), sin(2/9)` by hand; that matrix is not assembled
from algebraic generators (`transcendence_nogo.py` Part 3).

---

## 4. Why cuts 01–06 saw a kernel-independent fixed phase (and why 07 is barred)

The speculative cuts found a phase `≈ 0.742 rad` that did not move under any
compression / Shulga-kernel / mode-selection choice. This is a **corollary, not a
tuning failure**:

- The Shulga kernel and the compressions are *amplitude* maps (real positive
  weights) — a P3-style coarse-graining. A real scalar reweighting cannot rotate
  a phase (cut 01 proved this for the `Z₃` DFT).
- The phase is fixed by the *eigenphase* of the `J`-rotation, an **algebraic
  angle** pinned by the (algebraic) `J` matrices.
- No amplitude modulation can move an algebraic eigenphase onto a transcendental
  target.

So the phase is frozen at a kernel-independent algebraic value — exactly the
observed behaviour. **The recommended "07 cut with proper `J`'s" will produce a
different fixed algebraic eigenphase, still never `2/9`.** It is barred before it
is run. (`transcendence_nogo.py` Part 4.)

---

## 5. The exact-vs-approximate dilemma (robustness)

The no-go does not even depend on whether `φ = 2/9` is exact:

- **(I) `φ = 2/9` exact** ⟹ `cos(3φ) = cos(2/3)` transcendental ⟹ no algebraic
  invariant/character/eigenphase equals it. Algebraic program barred.
- **(II) `φ = 2/9` approximate** ⟹ `2/9` is a fit, not fundamental. The quark
  phases already show this texture: `φ_d ≈ 1/9` (2.3%), `φ_u ≈ −1/14` (1.5%),
  with non-uniform denominators `{9, 9, 14}` — a genuine algebraic law would be
  uniform and exact; it is neither.

Either branch closes the algebraic-phase program. (`transcendence_nogo.py` Part 5.)

---

## 6. Where the phase *can* live

The only number-type that can equal a transcendental is a **ratio of dynamical
(loop) coefficients**. The `φ`-potential is a Fourier series in `3φ`,
`V_eff(φ) = c₀ + c₃cos(3φ) + c₆cos(6φ) + …`, so its critical points obey
`cos(3φ) = −c₃/(4c₆)` — a ratio of loop coefficients that carry logarithms and is
therefore transcendental. This is the right category, and it is a *continuous*
function of the loop data (unlike a frozen algebraic eigenphase). Reproducing
`φ = 2/9` requires the ratio to equal `4cos(2/3) = 3.1436…`; doing so from a
complete Lagrangian (fixing the scale and the sector counting `N_X`) is the
genuine open task — v59's unfinished Step-6. See
`../residual_audit/phase_dynamical_home.py`.

---

## 7. Verdict

| object | number-type | algebraic construction can produce it? |
|---|---|---|
| `784`, `14`, gauge `5/2/9`, `sin²θ_W = 2/9`, Koide `Q = 2/3` | algebraic | **yes** (this is where rep theory wins) |
| Brannen phase invariant `cos(2/3)`, `cos(2/9)` | transcendental | **no** (this no-go) |

The phase is not a *missing* algebraic mechanism; it is a *non-algebraic* object.
Effort spent making the algebra emit it (better `J`'s, characters, Shulga
kernels) is barred by Corollaries A–B. The only viable path is the loop
Lagrangian (§6).

**Artifacts.** `transcendence_nogo.py` (numerical, self-verifying);
`../lean/PhaseNoGo.lean` (machine-checked, one cited `sorry` for L–W).
