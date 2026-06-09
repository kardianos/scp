# FINDINGS — Gap G2/G3 (the dimensionless couplings)

**Folder:** `v59/gaps/alpha_couplings/` · **Date:** 2026-05-25 · **Verdict document.**

All numbers below are reproduced by `verify_forms.py`, `overfit_scan.py`, `rg_fixedpoint_test.py`
in this folder; all clean identities are in `AlphaCouplingIdentities.lean` (written, not built
this run). Status tags per `README.md`.

---

## Per-gap verdict

### G2 — the `√α` form `g_W² = 5√α` → **DEAD as a law / OPEN only as a value**
- The **`5` is theorem-grade** (`h∨(Spin7) = 21−16 = 7−2`); never the problem. **[thm]**
- The **`√α` *form* is not a derivable law.** Machine-checked here
  (`AlphaCouplingIdentities.gW_sqrt_alpha_is_alphaMZ` + `…alphaMZ_satisfies_consistency`):
  given the SM-definitional `4πα = g_W²·sin²θ_W` with `sin²θ_W=2/9`, `g_W²=5√α` is **exactly
  equivalent** to the single value `α(M_Z)=25/(324π²)`. A line `18πα` and a root `5√α` meet at one
  point (`verify_forms.py [2]`). The intended `Ω²` Lagrangian mechanism is **excluded three ways**
  (`04_gW_sqrt_alpha_result.md`, avenue D) — **do not re-run.**
- **Net:** G2 contains no information beyond the G3 EW value. It is **dead as an independent gap**;
  what survives is the G3 value-match. (Note a downstream wrinkle: feeding the v59 `α(M_Z)` into
  `g_W=√(5√α)` gives `g_W=0.665` vs PDG `0.652`, a 2% slip — the relation is tightest on
  `α(M_Z)` itself, not on `g_W`.)

### G3 — α itself → **OPEN, currently a [free]/value-conjecture input**
- **EW form** `α(M_Z)=25/(324π²)`: clean identity **[thm]**, matches PDG to **0.032%** **[emp]**.
  Integers `5` and `2/9` have prior provenance ⇒ moderate evidential weight. But it is a *value*,
  not a *derivation* — there is no mechanism (it is the G2 √α-form in disguise).
- **IR form** `−ln α + 2α = π²/2`: matches CODATA to **3.6×10⁻⁵** **[emp]** — *but the tight match
  rides on the reverse-engineered `+2α`*. The mechanism-backed *pure* form `−ln α = π²/2 = 8π²/16`
  gives only **1.47%** (`verify_forms.py [3]`). The instanton template is blocked by `π₃(S⁷)=0`
  (no S⁷-based instanton) and an asserted `g²=16`. **[conj, partly fitted]**
- **Net:** α is **not derivable from v58/v59 as constituted** (consistent with `ALPHA_SCOPING.md`).
  The honest parameter count is **2 dimensionless-relevant inputs: `a_ℓ` (dimensionful) + α.**

---

## Overfitting verdict (the rigor question) — two-sided, honest

- **For v59, surprisingly favourable *within* a family:** in `α=(p/(mπ))²` the v59 `(5,18)` is the
  **only** reduced form beating 0.1% (`p≤30,m≤60`); in `−ln α + cα = aπ²/b` the v59 `(1,2,2)` is
  the **unique** small-integer hit to 0.1% (`overfit_scan.py [A],[B]`). The parametrizations are
  *tight*, not loose — earlier prose understated this.
- **Against, *across* families:** the real freedom is template choice. **5 of 8** arbitrary,
  unmotivated templates already hit `α₀⁻¹` to 0.03% with small integers (`overfit_scan.py [D]`).
  So "α has *a* clean structural form" is **cheap**; precision alone is **not** evidence.
- **Discriminator:** shared provenance of the integers (`5=h∨`, `2/9`, `16=dim Cl(3,1)` from the
  *same* algebra that fixes Koide/`sin²θ_W`). EW form scores moderate; IR form weak (free `+2α`).

---

## Most promising avenue

**Avenue C — RG / β-function (`ALTERNATIVES.md §C`, `rg_fixedpoint_test.py`).** It is the only
route that (i) sidesteps the fatal "`√α` is not a law" problem — standard running gives the correct
*linear* `g_W²=18πα` at every scale, and the v59 `α(M_Z)=1/127.91` **is consistent with the value
standard running passes through at `M_Z`** (orientation check: leptonic 1-loop alone moves `1/α`
into the `1/128` ballpark), reframing `5√α` as a benign single-scale coincidence rather than a
broken law; and (ii) has a concrete, falsifiable next step. The two live sub-tests:
1. **Boundary-value test:** impose the theorem-grade unification pattern `g_W²:g_R²:g_{B-L}²=5:5:2`
   at a unification scale, run 2-loop SM RGEs down, and check whether the **single** scale that
   reproduces `sin²θ_W(M_Z)=0.23121` *also* lands `α(M_Z)`. A simultaneous fit would be real
   evidence; a mismatch refutes the `5:5:2` origin.
2. **Fixed-point test:** the pure-SM U(1) β is positive (no zero); a genuine `β=0` requires the
   exceptional-group matter content of the `Spin(7)` embedding — **uncomputed, the real open work.**

Second choice: **Avenue B reborn** — replace the dead `π₃(S⁷)=0` instanton with the live
`π₃(Spin(7))=ℤ` / `π₃(G₂)=ℤ` sectors and check whether the minimal action forces `g²=16` (vs `14`
or `21`). If it forces `16`, the IR form's foundation is repaired; otherwise it refutes it.

**Dead/do-not-pursue:** D (`Ω²` self-coupling, excluded), and—on present evidence—E (eigenvalue:
finite-dim operators have algebraic spectra, α's forms are transcendental).

---

## What a real derivation would require

1. **A forced coefficient, not a fit.** Any α relation must come from a Lagrangian/RG computation
   that *forces* its integers (the `2` in `+2α`, the `16` in `8π²/16`, the `5:5:2` ratios) — with
   the integers traceable to the same algebra that fixes Koide and `sin²θ_W`. A formula whose
   integers are chosen to fit is numerology (the `overfit_scan.py [D]` baseline).
2. **A constructed topological object** (avenue B): an actual finite-action self-dual configuration
   with integer charge on a `π₃≠0` base inside the v59 algebra, whose action is computed (not
   borrowed) to be `8π²/g²` with an algebraically forced `g²`.
3. **Or a genuine fixed point** (avenue C): a `β(α*)=0` from the enlarged `Spin(7)`-embedded gauge
   + matter content, fixing α* with no input — or a unification boundary value that *jointly*
   predicts `sin²θ_W(M_Z)` and `α(M_Z)`.
4. **Disentangling G2 from G3 is impossible** as stated: any derivation of `g_W²=5√α` *is* a
   derivation of `α(M_Z)` and vice versa (`AlphaCouplingIdentities.gW_sqrt_alpha_is_alphaMZ`). So
   there is **one** dimensionless target (the α value/relation), not two.

---

## Bottom line

- **G2 = dead** as an independent law (it is `α(M_Z)`'s value in disguise; the `5` survives as
  `h∨`). **G3 = open**, with α a genuine value-conjecture input. The two are one gap.
- **The match precision is real but not, by itself, evidence** — the cross-template look-elsewhere
  is large. The within-family tightness and shared-provenance integers give the EW form *moderate*
  weight, the IR form *weak* weight (fitted `+2α`).
- **Highest-leverage live work:** the RG boundary-value/fixed-point test (avenue C), which can also
  *explain away* the `5√α` coincidence without needing it to be a law.
- **Honest parameter count for the lepton+EW block:** `a_ℓ` + α (two), unless avenue C or B
  succeeds.

### Lean build status
`AlphaCouplingIdentities.lean` is **written, not built this run** (no `lake build` on the shared
`furey_construction/lean` project, per task rules). It is self-contained (`import Mathlib` only) and
uses only `ring`/`norm_num`/`decide`/`field_simp`/`positivity`/`linear_combination` and standard
`Real.sqrt`/`Real.log` lemmas, in the exact idiom of the sibling modules that build cleanly. Two
`Prop`s carry `sorry`/placeholder and are flagged in-file as the genuine open steps (the IR
instanton derivation; the dead `√α`-Lagrangian origin).
