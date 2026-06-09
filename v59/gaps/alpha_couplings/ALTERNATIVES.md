# Solution avenues for G2/G3 — beyond the current conjectures

**Folder:** `v59/gaps/alpha_couplings/` · **Date:** 2026-05-25

The current v59 conjectures (IR `−ln α + 2α = π²/2`, EW `α(M_Z)=25/(324π²)`, gauge `g_W²=5√α`)
are *value-matches*, not derivations, and G2/G3 are entangled (`README.md §4`). This document
opens the solution space: avenues that would turn α from a [free]/[conj] input into something
*derived*. Each avenue states **(i)** the mechanism, **(ii)** how it would be tested, **(iii)**
what would falsify it, and **(iv)** a prior-status flag (some are already dead — recorded so we
don't re-run them).

The non-negotiable bar (from the project ethos): a real derivation must produce α (or the
relation it satisfies) from **the same algebra that already fixes Koide, `sin²θ_W`, the
generation count** — i.e. the integers must have *shared provenance*, not be chosen to fit.

---

## A. Fixed-point / self-consistency equations for α (generalizes the IR form)

**Mechanism.** The IR conjecture is itself a fixed-point equation `F(α)=0` with `F(α)=−ln α+2α−π²/2`.
Generalize: seek a fixed-point equation whose *operator* is forced by the v59 algebra, so the
root is α with no free coefficient. Candidates:
- a **gap equation** `α = G[α]` from a Schwinger–Dyson / NJL-type self-energy on the bivector
  connection `Ω`, where the kernel's normalization is the dual Coxeter `h∨` or `dim Cl(3,1)`;
- an **effective-potential stationarity** `∂V_eff/∂α = 0` where `V_eff` is the one-loop potential
  of the `Ω` self-coupling.

**Test.** Derive `F` from a Lagrangian; check that (1) the coefficient on the `α` (linear) term is
*forced* to be exactly `2` (not fitted), and (2) the constant is `π²/2` with `16` from `dim Cl(3,1)`.
`overfit_scan.py [B]` already shows that within `−ln α + c α = a π²/b`, only `(a,b,c)=(1,2,2)` hits
to 0.1% — so a derivation that *independently* forces `c=2` would be strong evidence.

**Falsify.** If the Lagrangian forces a different `c` (e.g. `c=1`, `c=0`), or forces a non-`π²`
constant, the IR conjecture is refuted. Equally, if no coefficient is forced (the `c` remains free),
the avenue is dead — α stays a fit.

**Prior status.** The `+2α` is currently reverse-engineered (`04_alpha_prediction.py`), so this
avenue is **open but unworked**. The instanton route to the `π²/2` is blocked by `π₃(S⁷)=0`
(avenue B below).

---

## B. Instanton / topological-action route (the IR template's intended mechanism)

**Mechanism.** `α = e^{−S}`, `S = 8π²/g²`, `g² = dim Cl(3,1) = 16` ⇒ `S = π²/2`. This is the BPST
action of a unit-charge instanton on the v59 constraint surface.

**Test.** Construct an actual instanton: a finite-action self-dual configuration on a 4-dimensional
base, with integer topological charge, on (or fibred over) the Furey surface `G₂\Spin(7) ≅ S⁷`,
and show its action is `8π²/16`.

**Falsify / current status — WEAK, partly DEAD.** `π₃(S⁷) = 0` ⇒ **there are no topological
instantons with S⁷ as base** (this matches the v6 gravity "Avenue C" null result, `π₃(S⁷)=0`).
So the `8π²/g²` is a *borrowed number*, not a constructed action. Reviving this requires either
(a) a *different* base with nontrivial `π₃` inside the v59 algebra (e.g. a `Spin(7)` or `G₂`
subgroup orbit with `π₃≠0` — `π₃(G₂)=ℤ`, `π₃(Spin(7))=ℤ`, so this is the live sub-avenue), or
(b) abandoning the instanton reading entirely. **Recommended concrete step:** compute the minimal
instanton action of the `π₃(Spin(7))=ℤ` or `π₃(G₂)=ℤ` sector and check whether the natural
coupling normalization yields `8π²/16` (or `8π²/14`, `8π²/21`). If a *bona fide* `π₃≠0` instanton
gives `g²=16` for an algebraic reason, the IR form's foundation is repaired. If the minimal action
is structurally `8π²/14` or `8π²/21`, the conjecture is refuted in favour of a different integer.

---

## C. RG / β-function fixed point for α

**Mechanism.** Instead of a *value*, derive a *relation* α must satisfy from running. Two sub-forms:
1. **UV/IR boundary condition.** If the v59 algebra fixes the gauge couplings at a unification
   scale via the Killing indices (`g_W²:g_R²:g_{B-L}² = 5:5:2`, all `∝√α`), then *one* boundary
   value plus the SM β-functions predicts `α(M_Z)` and `α(0)` by running. The content would be the
   *unification pattern* `5:5:2`, which is theorem-grade, not the α value.
2. **Genuine fixed point.** A nontrivial UV or IR fixed point `β(α*)=0` of an enlarged theory
   (e.g. asymptotic safety, or a Banks–Zaks-like point in the `Spin(7)` gauge sector) would *fix*
   α* with no input.

**Test.** (1) Take `5:5:2` at a chosen unification scale, run SM 2-loop RGEs down, compare
`α(M_Z)`, `sin²θ_W(M_Z)` to PDG. The *non-trivial* check is whether the **single** scale that makes
`sin²θ_W(M_Z)=0.23121` *also* lands `α(M_Z)`. (2) Compute `β` of the `Spin(7)`-embedded coupling
including the exceptional-group matter content; look for a zero.

**Falsify.** (1) If no single unification scale fits both `sin²θ_W` and `α(M_Z)`, the `5:5:2`
pattern is not the origin (it is then only the EW-scale coincidence of `README.md §4`). (2) If
`β` has no zero in the physical range, no fixed point exists.

**Prior status.** Untested in this folder. This is the **most standard-physics-respectable** route
and it sidesteps the "√α is not a law" problem entirely (running gives the *linear* `g_W²=18πα` at
every scale, which is correct; the `5√α` coincidence at `M_Z` would then be explained as the value
the running happens to pass through). **Recommended as the highest-leverage live avenue.**

---

## D. The `Ω²` bivector self-coupling (the intended G2 mechanism) — DEAD

**Mechanism (proposed).** v59's gauge field is the grade-2 connection `Ω`. A quadratic
self-interaction `λΩ² + μ⟨Ω,Ω⟩` was proposed to make a coupling that enters EM linearly appear as
`√` in the weak cross-channel, giving `g_W² ∝ √α`.

**Status — DEAD (excluded 3 ways, `04_gW_sqrt_alpha_result.md`):**
1. v58 *derives* Maxwell as `DΩ_EM = J_χ` — EM is a **primary** sector, not sourced by a weak `Ω²`
   (the balance the mechanism needs is the opposite of what was proved).
2. v58 has **no weak sector**, and the self-coupling `λ≈0.2` needed for the prefactor is **40×**
   outside v58's safe band `|λ|≤0.005`.
3. Standard gauge embedding `U(1)_em ⊂ Spin(7)` gives `e = g_W sinθ_W` ⇒ `α ∝ g_W²` (**linear**),
   never `√α`. The `A∧A` non-abelian term contributes the `g²W⁺W⁻` quartic but still through
   `e ∝ g_W`.

**Falsify (already done).** No projection of a `Spin(7)` connection's `Ω²` yields `g_W²∝√α`.
**Do not re-run.** The honest restatement is `g_W²=5√α` is the EW α-value repackaged (`README §4`).

---

## E. α as an eigenvalue of a v59 operator

**Mechanism.** Some theories obtain dimensionless constants as eigenvalues/spectral ratios of a
natural operator (e.g. the ratio of two Casimirs, or an eigenvalue of the mass kernel / the `ad`
action on `L=so(8)`). Candidate: α (or `√α`, or `−ln α`) as an eigenvalue of the bivector
self-coupling operator on `Λ²(𝒜)`, or as a ratio of structural traces.

**Test.** Enumerate the natural operators of the v59 algebra (the `ad`-action Casimirs on the 28-,
35-, 63-dim sectors; the mass-kernel eigenvalues `λ_k = a(1+2t cosθ_k)`; Frobenius norms) and check
whether any *dimensionless ratio* reproduces α to the precision of the conjectures **without** an
adjustable integer. The `KernelEigenvalues.lean` module already has the eigenvalue machinery.

**Falsify.** If every such ratio is either a small rational (no transcendental π/e content) or
requires a fitted integer to reach α, the eigenvalue route is empty. (Note: pure-rational
eigenvalues *cannot* give the transcendental `−ln α ≈ 4.92`; so this avenue must route through an
operator whose spectrum is genuinely transcendental, e.g. involving `exp`/`log` of a trace — which
is unusual and a red flag.)

**Prior status.** Untested here; **low prior** (the spectra of finite-dim algebra operators are
algebraic numbers, while α's structural forms are transcendental). Worth a quick scan, likely null.

---

## F. Dimensional transmutation

**Mechanism.** In QCD, `Λ_QCD = μ e^{−1/(b₀ g²(μ))}` turns a dimensionless coupling into a scale.
*Run it backwards:* if the v59 EW scale `v_Higgs = 28² a_ℓ²` (gap G1) and the lepton scale `a_ℓ`
are tied by a dimensional-transmutation relation, the *ratio* `v_Higgs/a_ℓ²` (a pure number, 784)
could be `∝ e^{c/α}` for a structural `c`. Then α would be fixed by the *measured* `784` instead
of conversely.

**Test.** Check whether `ln(v_Higgs/a_ℓ²) = ln 784 = 6.66` equals `c/α` or `c·(−ln α)` for a
structural `c`. (`6.66·α₀ = 0.0486`; `6.66/(−ln α₀) = 1.353` — neither is obviously structural.)

**Falsify.** If no structural `c` links `ln 784` to α, the avenue is dead. This avenue **couples
G3 to G1** (the `v_Higgs` gap) — only attractive if G1 also gets a transmutation mechanism.

**Prior status.** Speculative; flagged because G1 (`v_Higgs=28²a_ℓ²`) is the sibling load-bearing
conjecture and a *joint* transmutation story would explain both at once. Low prior, high payoff.

---

## G. Demotion (the honest null avenue): accept α as a [free] input

**Mechanism.** None — this is the conclusion if A–F all fail. α and `a_ℓ` are the two genuine
inputs; everything dimensionless else (Koide, `sin²θ_W`, `v_Higgs/a_ℓ²`, `m_W/m_Z`) is structural.

**Test / falsify.** This is the *default*; it is "falsified" only by a successful A–F. It is the
position `ALPHA_SCOPING.md` and `RIGOR_AUDIT.md` currently endorse: "α is a genuine value-conjecture,
not derivable from v58/v59 as constituted."

**Status.** This is the current honest baseline. The headline becomes **"2 inputs: `a_ℓ` + α."**

---

## Priority ranking (live → dead)

| avenue | prior | leverage | recommendation |
|---|---|---|---|
| **C. RG / β-function fixed point** | medium | high (sidesteps "√α not a law") | **work first**: run `5:5:2` from a unification scale; test the single-scale `sin²θ_W` + `α(M_Z)` joint fit |
| **B. instanton via `π₃(Spin7)=ℤ` / `π₃(G₂)=ℤ`** | low-med | high (repairs the IR form) | compute minimal action; check `g²=16` is forced (not `14`/`21`) |
| **A. forced fixed-point equation** | low-med | high | needs a Lagrangian forcing `c=2`; couple to B |
| **E. eigenvalue of v59 operator** | low | medium | quick scan, expect null (algebraic vs transcendental) |
| **F. dimensional transmutation (joint with G1)** | low | high if it works | only with a G1 transmutation story |
| **D. `Ω²` self-coupling** | — | — | **DEAD, do not re-run** |
| **G. demote α to [free]** | — | — | the honest baseline if A–F fail |
