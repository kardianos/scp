# Opening the solution space for R1 / R2 / the 25 directions

**Date:** 2026-05-25 · companion to `README.md`. Each candidate mechanism for the open residuals
gets: the **claim**, a concrete **test** (mostly already run in `formalize_bridge.py`), and a
**falsifier** (what would kill it). Goal: move beyond "posit Frobenius²" and either find a real
mechanism or sharpen exactly why each stays `[conj]`.

The residuals to close (from README):
- **R1** — why the EW scale `= ‖Y‖²_F` (Frobenius², not trace / top-eig).
- **R2** — why each of the 784 components has magnitude `= a_ℓ`.
- **25-dir** — `L` is 28-dim, only 3 generations light; account for the other 25.

---

## A. so(8)-invariant Higgs potential whose minimum is `‖Y‖²_F`  (R1 + R2)

**Claim.** Promote `Y ∈ End(L) = M₂₈(ℝ)` to a dynamical field with the most general
renormalizable `SO(28)`-invariant (adjoint-action) potential
`V(Y) = −μ²·I₁ + (λ/4)·I₁² + (g/4)·I₃ + (h/4)·I₂`, where `I₁ = ‖Y‖²_F`, `I₂ = (Tr Y)²`,
`I₃ = Tr((YYᵀ)²)`. Its minimum should both (i) read `v` as a Frobenius² (R1) and (ii) make the
components democratic of common scale `a_ℓ` (R2).

**Test (run — `formalize_bridge.py` §4).**
- *Case A* (pure Frobenius hat, `g=h=0`): the minimum is the sphere `‖Y‖²_F = 2μ²/λ`, direction
  **flat**. So a quartic **can** pin `‖Y‖²_F` — but to `2μ²/λ`, an *input*, not 784 a². This
  is R1's *form* (the scale is a Frobenius²) but **relocates** `v` into `(μ,λ)`.
- *Case B* (add shape term `g≠0`): the sign of `g` controls the vacuum's *spectrum*. Numerically
  `g>0` (`+Tr((YYᵀ)²)` penalizes large singular values) flattens the spectrum toward
  **democratic** (eff. rank ~23/28); `g<0` rewards concentration (eff. rank ~13). So a
  **democratic, full-rank** vacuum is a genuine *dynamical outcome* for `g>0`.

**Verdict: partial — the closest R1 gets, but does not close it.** The potential delivers (a)
"`v` is a Frobenius²" and (b) "the vacuum is democratic" — but the **overall scale** `2μ²/λ`
remains free (R2's value undriven), and a full-rank democratic vacuum has **28** light
directions, not 3 (the 25-dir tension). So this route **trades** R1+R2 for "why is `2μ²/λ =
784 a_ℓ²`," i.e. it does not reduce the input count.

**Falsifier.** If the only `SO(28)`-invariant minima are either (i) the singlet `Y∝I`
(`‖Y‖²_F = 28 a²`, the excluded √28) or (ii) rank-1 (`784 b⁴`, wrong power), then no invariant
potential yields the bridge and this route is dead. *Result:* neither — `g>0` gives a genuine
high-rank democratic basin, so the route is *alive but insufficient* (scale + generation-count
unresolved).

---

## B. A mass sum-rule / spectral identity (R1, sidestep the potential)

**Claim.** Avoid the dynamical potential: derive `v = ‖Y‖²_F` as a **spectral identity** of the
already-fixed Brannen kernel. The gravity sector already uses `Σm = Tr(M†M) = ‖M‖²_F` over the
3-dim generation space (`gravity_charge_test.py`); the bridge is the *same* Frobenius² over the
28-dim ambient. If a single identity `(EW scale) = ‖(mass bilinear)‖²_F` holds in *both* spaces,
R1 is one statement, not two.

**Test (run — README §1, `gravity_ew_coeff`).** The cross-check `Σm_ℓ = (9Q/dim(L)²)·v =
(6/784)·v` holds at **0.07%** — i.e. the *same* `784 = dim(L)²` relates the generation-space
Frobenius² (`Σm = 9Qa²`) to the ambient-space Frobenius² (`v = 784a²`). This is a real
consistency: one `a²`, two Frobenius² readings, tied by `dim(L)²`.

**Verdict: promotes R1 to "one identity in two spaces," still `[conj]`.** It shows the Frobenius²
reading is *internally consistent* across gravity and EW (a non-trivial check it passes), but the
*physical reason* the EW vacuum energy equals a Frobenius² is still unproved. It is the most
economical framing: **one** conjecture (`scale = ‖bilinear‖²_F`) now covers both the EW bridge
and the gravity charge.

**Falsifier.** If the gravity charge that fits the data were a *different* invariant of `M` than
the EW scale is of `Y` (e.g. gravity ∝ `Tr M`, EW ∝ `‖Y‖²_F`), the "one identity" framing breaks.
*Result:* both are Frobenius² (EP-exact for gravity, README §1 for EW) — survives.

---

## C. Composite condensate / NJL route (R2 + the 25 directions)

**Claim.** `Y^{ab} = ⟨ψ̄_L^a ψ_R^b⟩` is a **composite** of two *independent* 28-legs, so it has
`28×28 = 784` independent components from the start (Option 2). A symmetry-breaking (democratic)
condensate makes all 784 components comparable; the scale `a_ℓ` is set by an NJL **gap equation**.

**Test (run — `formalize_bridge.py` §5, `obe_options_2_5.py`).** Component counting confirmed:
democratic `Y` → `‖Y‖²_F = 784 a²`; singlet `Y∝δ` → `28 a²` (excluded); rank-1 outer product →
`784 b⁴` (wrong power of scale). So 784 needs **two independent 28-legs**, which the composite
supplies naturally.

**Verdict: live as a *counting* mechanism; R2 still open + a new prediction.** Two legs give the
784 count cleanly. But (a) the **dynamical scale `a_ℓ`** is an unsolved gap equation (so R2 is
*reframed*, not closed), and (b) a composite Higgs implies a **compositeness scale** — a
*collider-testable* consequence (an honest cost, not free). The 25 directions become the broken
`so(8)` generators (Goldstones), but see the count problem below.

**Falsifier.** A composite Higgs predicts deviations in Higgs couplings / a new strong scale at
the TeV–10 TeV range; **LHC/HL-LHC Higgs-coupling precision** that finds the Higgs exactly
SM-pointlike to the relevant level would disfavor a low compositeness scale. (Does not falsify a
very-high-scale condensate, but removes its phenomenological appeal.)

---

## D. Identify the v58 vacuum field with the L-grade mass bilinear (R1's "home")

**Claim.** The R1 obstruction is "v59's actual Higgs (`XiVacuum`, on `ℍ≅ℝ⁴`) is 4-dim and
dimensionless, so it cannot *be* the 28-dim bilinear." Fix this by positing the **physical Higgs
multiplet is the L-grade bilinear `Y ∈ End(L)`** (28×28), with `XiVacuum`'s `ℍ≅ℝ⁴` the *radial +
3-Goldstone* sub-sector that actually gets gauged (the `SU(2)_L` of `SilentDirection.lean`).

**Test (partly run — `XiVacuum.lean` `[thm]`; extension needed).** `XiVacuum` already proves the
`{λ,0,0,0}` spectrum: 1 massive radial Higgs + 3 Goldstones eaten as `W±,Z`. The *extension* this
route needs: embed that 4-dim sector as the "active" block of the 28×28 `Y`, with the remaining
`28−4 = 24` directions either heavy or the missing Goldstones. **Not yet built** — the explicit
`ℝ⁴ ↪ End(L)` embedding (the deferred Witt/generation map, `SevenDAlgebra.lean`) is required.

**Verdict: the most physically-grounded route, but unbuilt.** It is the only candidate that
connects to v59's *existing, Lean-proven* Higgs mechanism (`XiVacuum`). If the `ℝ⁴ ↪ End(L)`
embedding is constructed and the radial mode's VEV comes out `= 784 a_ℓ²` as a Frobenius², R1+R2
close together. This is the recommended next build (see FINDINGS).

**Falsifier.** If every embedding `ℝ⁴ ↪ End(L)` consistent with the `SU(2)_L` gauging gives a
radial VEV that scales as `√28·a` (vector/equipartition) rather than `28·a` (bilinear), the
identification is excluded — it would mean the physical Higgs is a *vector in L*, not a
*bilinear*, reviving the dead equipartition reading.

---

## E. Account the 25 directions as eaten Goldstones of broken so(8)  (25-dir)

**Claim.** Spontaneous breaking `so(8) → H` makes `dim(so(8)/H)` directions Goldstones; choose
`H` so that exactly 3 directions stay light (or 25 are eaten), giving 3 generations.

**Test (run — `formalize_bridge.py` §6).** **NULL.** No proper subalgebra of `so(8)` has
dimension 25 (would need `dim H = 3`; the only 3-dim subalgebra is `su(2)`, but `so(8)/su(2)`
is 25-dim *as a coset*, not a symmetry that leaves 3 *generations* — and the maximal proper
subalgebras top out at `dim Spin(7) = 21`). No single-step `so(8) → H` leaves exactly 3 light
directions matching the 3 generations.

**Verdict: single-step breaking is dead.** The clean "3 light + 25 eaten" picture does not arise
from one classical subgroup. Two live alternatives remain:
1. **A breaking chain** `so(8) → … → (3 light)` with several Goldstone stages (untested,
   speculative).
2. **The vacuum is rank-3, not full-rank.** If `Y` is rank-3 (only 3 nonzero singular values =
   √m_e,√m_μ,√m_τ), there are exactly 3 light directions *by construction* — but then `‖Y‖²_F =
   m_e+m_μ+m_τ = Σm`, the **gravity charge** (`9Qa²`), **NOT** `784 a²`. The 784-count requires
   full rank.

**Falsifier / the sharp tension.** These two requirements are **mutually exclusive**: the
`784`-count needs a *full-rank* (28 light) democratic `Y`; the 3-generation spectrum needs a
*rank-3* `Y`. So **R1's "Frobenius² = 784 a²" and the 3-generation count cannot both be read off
the same `Y`.** This is the deepest open problem of the bridge (see FINDINGS) — and a clean
falsifier of the naive "one democratic 28×28 vacuum does everything" story.

---

## F. Reading `dim(L)²` as an internal volume (Option 7, R1 alternative)

**Claim.** No bilinear at all: the kernel factorizes `K = K_spatial · K_internal(M_L)` over the
28-dim L-coset; the lepton generations are the lowest harmonics on `M_L`, and the EW coupling
dilutes as `1/Vol_int ∝ 1/dim` per leg, with the *square* coming from two factors of the kernel.

**Test (run — `obe_geometric_6_7.py`).** The 1/r² law is the 4D massless zero mode
(independent of `M_L`); the zero-mode coupling scales exactly as `1/Vol_int` (1.000, 0.500,
0.250 for Vol=1,2,4). So the 28-count and the radial law live in *separate factors* of `K`.

**Verdict: live but does not produce `784` from a single principle.** It explains *one* factor
of `dim(L)` (the internal-volume dilution) cleanly, but getting `dim(L)²` needs the dilution to
act *twice* (two legs / a bilinear) — so it collapses back to C/D for the square. Useful for
keeping the 1/r² gravity law compatible with the 28-count, not for deriving R1.

**Falsifier.** If the zero-mode coupling scaled as `1/Vol^(1/2)` or `1/Vol²` rather than `1/Vol`,
the "one factor of dim per leg" reading would be wrong. *Result:* exactly `1/Vol` — survives as a
*single*-power mechanism (hence needs two legs for the square).

---

## Summary table

| route | targets | status | key cost / open piece |
|---|---|---|---|
| A. so(8)-invariant potential | R1+R2 | partial | scale `2μ²/λ` free; full-rank → 28 light, not 3 |
| B. spectral identity (gravity↔EW) | R1 | promotes to 1 conjecture | *why* scale = Frobenius² still unproved |
| C. composite/NJL (two 28-legs) | R2, count | live as counting | gap-equation `a_ℓ`; compositeness scale (testable) |
| D. v58 vacuum = L-bilinear | R1+R2 | most grounded, **unbuilt** | needs explicit `ℝ⁴ ↪ End(L)` embedding |
| E. 25 = eaten Goldstones | 25-dir | **single-step dead** | full-rank(784) vs rank-3(spectrum) tension |
| F. internal-volume (KK) | R1 alt | live, one factor only | square still needs two legs (→ C/D) |

**Most promising:** **B + D combined** — frame R1 as the single conjecture "physical scale =
Frobenius² of the mass bilinear" (B, consistent across gravity and EW), and give it a home by
identifying the v58 `XiVacuum` radial Higgs with the active block of `Y ∈ End(L)` (D). The
blocker for both is the **rank tension** (E): the count wants full rank, the spectrum wants
rank 3. Resolving that is the crux.
