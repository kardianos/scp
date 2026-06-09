# Residual Partition + Input Count

v61 `CLOSEOUT.md` lists four residual conjectures "all value/symmetry inputs, on
the same footing as `α`": `α`, `v = 784a²`, the Brannen phase `φ = 2/9`, and
`f_g ~ α^{21/2}`. **That lumping is too coarse.** By number-type and flatness
(`../THESIS.md` §§1,3) they split into three kinds with three different fates.

---

## The partition

| residual | kind | flat coordinate | fate | artifact |
|---|---|---|---|---|
| scale `a` (`= √v / dim L`) | pure scale | — | **unavoidable units anchor** (not a defect) | — |
| `v = 784a²` democracy | degeneracy | `S^783`, rational points | **closable** by state-selection | `democracy_selection.py` |
| Brannen phase `φ = 2/9` | transcendental value | `cos(3φ) = cos(2/3)` | P3 loop-ratio or input — **barred from algebra** | `phase_dynamical_home.py`, `../no_go/` |
| `α`; `f_g ~ α^{21/2}` | transcendental coupling | — | input; `f_g` downstream of `α` | — |

`784` itself is algebraic and Burnside-forced — a P1 *win*, not a residual; only
its *alignment* (democracy) was ever open.

---

## Democracy is closable (the one tractable residual)

`democracy_selection.py` shows the democratic vacuum is selected by **any** of
three principled, non-transcendental criteria (all self-verified):

1. **Equivariance** — the democratic point `(1,…,1)/√n` is the *unique* point on
   the sphere fixed by the generation permutation group `Sₙ` (the
   permutation-fixed subspace is 1-dimensional).
2. **Max-entropy** — democracy uniquely maximizes the entropy of the
   squared-component distribution (`S = log n`, attained nowhere else; 20 000
   random sphere points all fall short).
3. **Explicit breaking** — adding the leading `Sₙ`-invariant quartic the hat
   omits lifts **all `n−1` Goldstones** and pins democracy as the strict minimum
   (tangent Hessian becomes positive-definite).

Because the flat coordinate is rational/algebraic, none of this needs a
transcendental. Closing it removes `v = 784a²` from the residual list, leaving
only the scale `a`.

**Contrast with the phase**: the phase's flat coordinate has a *transcendental*
value, so no state-selection reaches it — only the dynamical loop ratio
`−c₃/(4c₆)` (P3) can (`phase_dynamical_home.py`). Same geometry, opposite
resolution, dictated by number-type.

---

## Input-count corrections

- **`f_g ~ α^{21/2}` is probably not an independent input.** `21 = dim Spin(7)`,
  so the exponent is `dim Spin(7)/2` — algebraic structure on the transcendental
  base `α`. If `21/2` is *derived* (not fit), `f_g` is a function of `α`; the
  closeout over-counts by one. **Open**: confirm the exponent is forced.
- **`v = 784a²` reduces to the scale `a`** once democracy is selected.

**Trend of genuine inputs**: one scale `a` + the transcendental coupling `α`
(with `φ`, `f_g` downstream of dynamics / `α` if the P3 and exponent tasks close).
Input count — not tolerance-hitting — is the metric that decides whether the
program is reducing assumptions or relabeling them.
