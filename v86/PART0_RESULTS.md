# v86 Part 0 + census — results
**Date:** 2026-07-26 · **Scope:** the ε/foundation battery (N1–N6, N9, N10),
the census rungs HC-1/HC-2/HC-3, and the exchange rung EX-2. All zero-spend
(CPU/analysis only). Instruments, logs and TSVs live in `v86/n_battery/`.

Each rung was implemented against the canonical protocols in
`v86/council/grok45/N_BATTERY_REVIEW.md` §1, reviewed adversarially by a
grok-4.5 seat **before** execution (`v86/n_battery/REVIEW_n12_grok45.md` and
the N9/N4 review), and every required fix applied. Verdict tags follow
THEORY_v86: **[M]** measured · **[D]** derived · **[C]** conjecture.

---

## 0. Headline

Part 0's central result is a **virial reduction of the excess**:

> **Σ ≡ E − Σ_a ω_a Q_a = (2/3)(E_∇ − E_g)**, with residual ≤ 5×10⁻⁸ of E
> across every branch, every coupling, every flavour partition and every
> potential in the (μ,κ) family — while ε = Σ/ωQ itself is 0.45–4.34%.
> **[D, verified M]**

Read this precisely. Σ is *reduced* by the virial theorem to two thirds of the
gradient-minus-gauge imbalance; ε remains a **computed profile functional** —
you still need the shooter to get E_∇(Q,g). What is new is a sharp,
parameter-free expression with no fitted constants, plus a *derivation* of the
universally positive sign (ungauged, Σ = (2/3)E_∇ > 0 identically). Given I1
and I3, I4 is logically forced; confirming it to 5×10⁻⁸ verifies that the
energy bookkeeping, the Gauss law and the virial identity are mutually
consistent on the branch — it is not a second independent miracle.

Three further results are load-bearing for the program:

* **n(H_ω) = 1** over the accepted (box-adequate) ω range (HC-1). Both indices
  of the corrected GSS criterion are now measured *there*; outside that range,
  and in the gauged Coulomb phase, it is not established.
* **The internal-harmonic count is not a constant of the theory** (HC-1/HC-2).
  It is a function of branch position: a ladder of bound monopole modes in the
  thin-wall region, **zero bound internal modes above ω ≈ 1.35** — the region
  every production run of this program has lived in. PROPOSAL I.1's slogan is
  falsified *as a bound-mode census*.
* **The exchange channel is matter, not gauge** (EX-2): the gauge share of the
  escaping power is ~10⁻⁵, and the energy budget is consistent at the
  order-of-magnitude level under coarse sampling (114%, cadence-limited).

And one correction the documents must absorb:

* **FOUNDATIONS R2 / GROUNDING §0 / THEORY_v86 A4 state the excess with the
  wrong gauge sign.** They write "gradient + potential + **gauge** − kinetic".
  With the true integrated kinetic energy the correct statement is
  Σ = E_∇ + E_m + E_V − E_kin **− E_g**. The `+E_g` is off by 2E_g. (It becomes
  correct only if "kinetic" is read as the *throughput* ωQ/2, which by the
  Gauss identity already contains E_g.) See §1.

---

## 1. N1 — component decomposition of E and Σ  **[M]**

Instrument `n_battery/n12_decomp.py`; production grid H=0.004, RMAX=150
(matching the frozen `v69/theory/gscan.tsv`), 292 branch points over
g ∈ {0, 0.02, 0.05, 0.10}.

**Conventions (this is the part the docs get wrong).** With the shooter's
E_kin = ∫(3/2)wt²f²dV (wt = ω − χ, the *true* kinetic energy), the Gauss law
gives an exact identity,

    I1:  ω Q = 2 (E_kin + E_g)              residual ≤ 4.9×10⁻⁹ of E

so the "throughput" ωQ/2 and the true kinetic energy differ by exactly E_g.
The excess is therefore

    Σ = E_∇ + E_m + E_V − E_kin − E_g   =   E_∇ + E_m + E_V − (ωQ/2)

and never `+E_g`.

**Phase A** (cross-check vs the frozen gscan artifact) is **CONSISTENT**:
every apparent difference — including two rows that a naive 0.2% tolerance
flagged — is explained by ω-matching resolution on a near-vertical branch
(dQ/dω up to −3.7×10⁷). Re-scored in `n1_phasea.py`: 0 of 15 sampled rows need
an ω offset larger than the 1×10⁻⁴ matching tolerance.

**Phase B** (the full split) is delivered in `n1_decomp.tsv`. The raw fractions
of Σ are large and cancelling (E_m and E_kin are each ~10–100× Σ), which is
expected and is exactly why the N2 reduction below is the useful statement.

---

## 2. N2 — the virial identity, and a closed form for the excess  **[D/M]**

**Derivation.** The free-extremum functional reproducing *both* radial ODEs is
F = E_∇ + E_m + E_V − E_kin − E_g — the gauge energy enters with a **minus**
because electrostatics is a saddle point in a₀. Using I1, Σ = E − ωQ = F.
F is stationary under the Derrick family f,χ → f(r/L),χ(r/L), under which
E_∇, E_g scale as L and E_m, E_V, E_kin as L³, so

    I3 (virial):  R_vir ≡ (E_∇ − E_g) + 3(E_m + E_V − E_kin) = 0
    I4 (closed):  Σ = (2/3)(E_∇ − E_g),   ε = (E_∇ − E_g) / (3(E_kin + E_g))

**Measured residuals (production grid, relative to E, max over each branch):**

| g | n | I1 (Gauss) | R_vir (N2) | I4 (closed form) | Newton resid | ε range |
|---|---|---|---|---|---|---|
| 0.00 | 209 | 4.4e-16 | 1.45e-07 | 4.83e-08 | 6.9e-09 | 0.0045–0.0434 |
| 0.02 | 46 | 2.1e-10 | 1.44e-07 | 4.77e-08 | 8.6e-10 | 0.0094–0.0429 |
| 0.05 | 26 | 1.3e-09 | 1.39e-07 | 4.51e-08 | 9.3e-10 | 0.0159–0.0404 |
| 0.10 | 11 | 4.9e-09 | 1.22e-07 | 3.57e-08 | 1.1e-10 | 0.0213–0.0313 |

**N2 GATE: PASS** — six orders of magnitude between the identity residuals and
the physical ε signal. The pre-registered falsification conditions (residual
systematic at the ε scale; parts failing to rebuild E; Gauss residual above the
integration floor) did not fire.

**What I4 is and is not.** I1 (Gauss) and I3 (Derrick) are two independent
constraints; I4 is their algebraic elimination of the bulk terms, so given
I1+I3 it is logically forced. Verifying it to 5×10⁻⁸ therefore confirms the
*virial package* — energy bookkeeping, Gauss law and stationarity are mutually
consistent — rather than predicting ε from nothing. I4 does **not** give ε(Q)
without a profile: E_∇ must still be computed. FOUNDATIONS' own N2 target,
"ε as a derived function of measurable profile moments with no free fit
parameters", is met; "ε is closed" in the stronger sense of a standalone
formula for ε(Q) is **not** claimed.

Two consequences worth stating plainly:

* **Ungauged, Σ = (2/3)E_∇ > 0 identically** — this *derives* the universally
  positive ε that the program had only ever measured.
* Σ < 0 is not forbidden by the identity: it requires E_g > E_∇, i.e. the deep
  gauged/large-Q corner. Nothing on the scanned grid reached it.

N9 independently confirms the identity is potential-shape independent: over all
25 (μ,κ) theories, max |Σ − (2/3)E_∇|/E = 1.9×10⁻⁶.

---

## 3. N3 — scaling collapse  **[M]**

Instrument `n_battery/n3_collapse.py`, 292 production points. R frozen to
r_half; the Coulomb subtraction uses the **measured** E_g (not a point-Coulomb
estimate), which by N2 makes the subtracted excess exactly

    ε' = ε + (2/3)E_g/(ωQ) = (2/3)E_∇/(ωQ)      (verified to 5×10⁻⁸)

| collapse variable | metric (0 = perfect) |
|---|---|
| ε vs ξ/R | 0.238 |
| ε vs g²Q/R | 0.916 |
| **ε' vs ξ/R (measured E_g removed)** | **0.148** |

**Verdict: geometry dominant, Coulomb correction real.** ξ/R already organises
the raw excess; removing the measured gauge energy tightens the collapse 1.6×,
with the improvement concentrated at the largest g. g²Q/R alone explains
essentially nothing, so the gauge sector is *not* an independent scaling
variable — it enters only through the −(2/3)E_g term N2 already isolates.

**Thin-wall behaviour, with its caveats attached.** Over the 34 thinnest-wall
points, a **local** log-log fit gives ε ~ (ξ/R)^1.011, and ε falls to 0.0045 at
ξ/R = 0.048. Exponent ≈ 1 is what Σ ∝ 4πR²σ_eff with Q ∝ R³ predicts *if*
σ_eff and ω vary slowly. Three limits must travel with that number:

* it is a **local fit at the thin-wall end**, dominated by large-Q, g=0 points
  — not a multi-g proof of a surface law;
* ε → 0.0045 is "heading toward zero", not a measured asymptotic zero;
* **R_rms/r_half varies 0.76–1.36 over the branch (56% spread)**, so the
  collapse abscissa is O(1) redefinable. The *existence* of the collapse
  survives that; a quantitative σ_eff fit, or quoting 1.011 as a physical
  precision, does not until R is frozen and the sensitivity is inside the
  error budget.

Likewise the collapse metric (within-bin spread over total spread) is a
**ranking statistic, not an R²**: 0.148 still means non-negligible residual
structure (thick wall, potential shape, residual g-dependence). "Geometry
dominant, Coulomb correction real" is the supported statement; "the surface
law, measured" is not.

---

## 4. N4 + HC-3 — flavoured branch  **[M]**

Instrument `n_battery/n4_hc3_flavored.py` (ungauged, 3 components, corrected
energy functional; normalisation asserted equal to the v66 symmetric
convention to 1e-12).

**N4.** Detuning ray (ω₀ down, ω₁=ω₂ up) at fixed total charge Q* = 209.567,
held by a Newton-corrected uniform ω shift (|ΔQ|/Q* ≤ 5×10⁻⁷ on every accepted
row).

* Σ moves by **2.52%** across the scan…
* …and E_∇ moves by **2.52%**, with I4 holding to 3.0×10⁻⁷.

**Verdict: geometry-driven partition dependence.** The flavour split acts on
the excess through *profile shape*, not through a separate V(s) interference
channel — because if it acted through the potential coupling it would break
I4, and I4 does not budge. The multi-frequency Derrick derivation carries over
unchanged: Σ = (2/3)E_∇ for any partition.
*Path limitation (pre-registered):* one detuning ray at one total charge.

**HC-3.** D_ab = ∂Q_a/∂ω_b by central differences (h=10⁻³), symmetrised
(D = Dᵀ on a smooth branch; asymmetry reported as an FD diagnostic).
Every scanned partition carries eigenvalues of signature (−,+,+):
**n(D) = 1 throughout**, a contiguous region rather than an isolated VK point.
Combined with HC-1's measured n(H_ω) = 1 (§6), the corrected GSS indices
**match across the whole scanned flavoured region** — in the ungauged theorem
regime, which is the anchor GROUNDING §1 asks for.

---

## 5. N5 + N6 — action, charge and the ħ_eff triple  **[M]**

Instrument `n_battery/n56_action.py`, 46 archived diag files, averaged over the
last 50% of each run.

The **CANNOT-VERIFY flag is resolved**: the diag does expose E_phi_kin, so
T = 2E_kin is directly available at the 0.25 t.u. diag cadence — N5 is free.

A convention hazard was found and is now documented: the diag's `omega_core` is
a **point sample** of the local phase rate at the argmax-s voxel, i.e. the
gauge-shifted clock at the ball centre, not the bare branch ω. This rung
therefore inverts the N1 Gauss identity to obtain the honest bare frequency,
ω̄ = 2(E_kin + E_em)/Q, and reports the difference as a measurement.

Across the 45 single-object runs (pair/ball+cloud runs, whose Q is a *net*
charge, are flagged and excluded):

| quantity | measured |
|---|---|
| in-kernel excess ε̄ = E/(ω̄Q) − 1 | 0.0261 … 0.0489, median 0.0366 |
| gauge share of the throughput 2E_em/(ω̄Q) | 7.5×10⁻⁵ … 0.031 |
| clock shift ω_core/ω̄ − 1 | −0.0019 … −0.0407 |

**N6.** ħ_E/Q − 1 is *identically* ε̄ (since ħ_E = E/ω̄), and lands at
0.0261–0.0489. This settles a FOUNDATIONS prediction outright: **v70's "ħ_eff =
Q to 3–5%" and the shooter's "ε = 1–4%" are the same residual measured twice.**
The identity language should be dropped in favour of the measured ratio, as
GROUNDING §0 already asks.

**N6 is now COMPLETE — see §5b.** The third leg was recovered from the archived
v70/v68 boost series by measuring the momentum instead of assuming it.

The nonzero, g-dependent clock shift is the direct empirical content of
GROUNDING §0's three-way split: charge, action and the local clock are three
different quantities, measured rather than argued.

---

## 5b. N6 completed — the ħ triple with a measured momentum  **[M]**

Instruments `n_battery/sfa_momentum.c` (new) + `n6_triple.py`, on the archived
boost series `/space/scp/v68/bs0{3,5,7}.sfa` (the runs behind v70/FINDINGS §2).

**Why v70's third leg did not count.** v70 reports a column "ħ_eff = p/k", but
its p was *not* measured: for a rigidly boosted soliton p = γE₀v and k = γωv,
so **p/k ≡ E₀/ω identically**. That column is ħ_E wearing ħ_pk's name, and the
triple had only two independent legs. (Its variation across the series comes
from k departing from γωv — the lattice dispersion — not from any momentum
measurement.)

**What changed.** `sfa_momentum.c` integrates **p_i = ∫T^{0i} dV** (matter +
gauge Poynting) directly from the frames, with the kernel's own covariant link
differences. Nothing kinematic is assumed in obtaining p.

| run | Q | E_lab | p (measured) | β = p/E | v (tracked) | β/v − 1 | E₀ = E√(1−β²) |
|---|---|---|---|---|---|---|---|
| bs03 | 481.85 | 723.86 | 215.96 | 0.2983 | 0.2977 | **+0.20%** | 690.90 |
| bs05 | 481.85 | 794.45 | 393.03 | 0.4947 | 0.4864 | **+1.7%** | 690.42 |
| bs07 | 481.84 | 948.82 | 651.79 | 0.6869 | 0.6678 | **+2.9%** | 689.51 |

Two checks fall out that were not available before:

* **E₀ is constant to 0.20%** across γ = 1.05 … 1.38, recovered from measured
  (E, p) alone with no reference to Q, ω or any mass hypothesis — the boosted
  objects really are the same object.
* **E₀ = 690.28 measured vs 691.90 from the ungauged shooter** at matched
  charge (ω = 1.39, Q = 482.18) — agreement to **0.23%**.
* β = p/E tracks the independently measured COM velocity to +0.2% at v = 0.30,
  rising to +2.9% at v = 0.67 — precisely the documented 1–5% lattice
  group-velocity anomaly, growing with v as it must.

**The triple** (ω = 1.39, ε = 0.0323 from the rest branch):

| run | ħ_Q = Q | ħ_E = E₀/ω | ħ_pk = p/k | ħ_E/Q − 1 | ħ_pk/Q − 1 | ħ_E/ħ_pk − 1 |
|---|---|---|---|---|---|---|
| bs03 | 481.85 | 497.05 | 493.85 | 0.0315 | 0.0249 | 0.0065 |
| bs05 | 481.85 | 496.71 | 489.88 | 0.0308 | 0.0167 | 0.0139 |
| bs07 | 481.84 | 496.05 | 480.74 | 0.0295 | −0.0023 | 0.0319 |

**Verdict.** Every residual sits in ε's family. "ħ_eff = Q" is not an identity
and the program should stop writing it as one: ħ_E/Q − 1 *is* ε by construction
(ħ_E/Q = E₀/(ωQ)), so v70's "3–5%" and the shooter's ε were never two findings.
The one genuinely new number is ħ_pk, whose departure from ħ_E grows 0.6% →
3.3% across the series — that gap is the lattice dispersion in k, and it is the
honest precision floor on any momentum-based mass statement in this kernel.

**What this does not do:** it is kinematics, not dynamics. Recovering E₀ from
(E, p) shows the boosted object is relativistically consistent; it does **not**
decide whether the inertia that resists a *force* is E/c². That is N7.

---

## 6. HC-1 — the BdG mode catalog  **[M/D]**

Instrument `n_battery/hc1_bdg.py`. Ungauged (the theorem regime GROUNDING §1
names as the anchor). Perturbing Φ_a = e^{iωt}(f + x_a + i y_a) gives

    L₀        = −∇² + m² − ω² + P₀                    (y block, all 3 channels)
    L_x^sym   = −∇² + m² − ω² + 5P₀ + 3B              (flavour-symmetric x)
    L_x^flav  = −∇² + m² − ω² − P₀                    (2 flavour x channels)

with P₀ = μf⁴/(1+κf⁶)², B = −4κμf¹⁰/(1+κf⁶)³. Discretised as u = rφ →
symmetric tridiagonal, so the negative counts are trustworthy by construction.

**Structural checks (each a real falsification opportunity — all passed):**

1. L₀'s lowest l=0 eigenvalue is the U(1) Goldstone: **−3.5×10⁻¹²**, with
   eigenvector overlap **1.000000** against u = r·f.
2. L₀ contributes **0** negative directions (its zero mode is nodeless).
3. L_x^flav contributes **0** (A = 2P₀ < 0 ⇒ L_x^flav = L₀ − A ≥ L₀).
4. L_x^sym contributes exactly **1** — the l=0 dilational direction.

⇒ **n(H_ω) = 1 at every ω scanned.** The GSS/VK anchor is measured, and HC-3's
provisional assumption is now an input rather than a hope.

**The continuum edge — a correction the census needs.** In the co-rotating
frame the coupled (x,y) system's asymptotic dispersion is
(Ω²−G)² = 4ω²Ω² with G = m²−ω²+k², giving **Ω_c = m − ω**, not √(m²−ω²).
At ω = 1.42 that is 0.080 versus 0.483 — a factor 6. Equivalently, in the lab
frame a perturbation at co-rotating Ω radiates iff ω+Ω ≥ m, which is exactly
GROUNDING §2's band statement. The single-operator eigenvalue maps to frequency
as Ω(λ) = √(ω²+λ) − ω, which sends λ_c = m²−ω² to Ω_c = m−ω exactly.

**Box adequacy is a first-class systematic, not a footnote.** The first pass
used a fixed R_max = 40 and produced a *spurious* result at the thin-wall end:
at ω = 1.315 (Q = 3.7×10⁵) the Goldstone came out at 6×10⁻³ with eigenvector
overlap 0.78 against r·f, and the run then reported n(H_ω) = 0. The box simply
did not contain the ball. The instrument now

* sizes the box to the object, **R_max = max(40, 5 r_rms)**, from a scouting
  solve, and
* applies a **hard gate**: if |λ₀| > 100h² or the Goldstone overlap < 0.99, the
  ω is **rejected** and nothing from it — least of all n(H_ω) — is reported.

With the adaptive box, **all 8 scanned ω pass the gate and none is rejected.**

**Box-state discrimination.** Every spectrum is solved at two box sizes; a
sub-threshold frequency survives only if its shift between boxes is below
min(2% of the highest sub-threshold mode, ¼ of the larger box's local level
spacing) — the spacing term matters because near threshold the box ladder is
dense and a flat 2% tolerance admits accidental alignments.

**The branch-resolved census (box-converged):**

| ω | Q | r_rms | R_max | Ω_c = m−ω | n(H_ω) | bound internal modes |
|---|---|---|---|---|---|---|
| 1.3200 | 144006 | 22.63 | 113.2 | 0.180 | 1 | **none** |
| 1.3300 | 22207 | 12.11 | 60.6 | 0.170 | 1 | 0.0939 |
| 1.3400 | 7222 | 8.35 | 41.8 | 0.160 | 1 | 0.1352 |
| 1.3600 | 1745 | 5.31 | 40.0 | 0.140 | 1 | none |
| 1.3900 | 482 | 3.71 | 40.0 | 0.110 | 1 | none |
| 1.4200 | 210 | 3.14 | 40.0 | 0.080 | 1 | none |
| 1.4500 | 118 | 3.09 | 40.0 | 0.050 | 1 | none |
| 1.4700 | 92 | 3.43 | 40.0 | 0.030 | 1 | none |

**The box fix changed the physics answer.** At ω = 1.32 an undersized box
(R_max = 40 ≈ 1.8 r_rms) had reported *three* "box-stable" bound modes; with
R_max = 113 there are **none**. The ladder was an artifact. The corrected
census is: **at most one bound internal monopole mode, in a narrow window
ω ≈ 1.33–1.34, and zero everywhere else.**

**Scope of the n(H_ω) = 1 claim.** It holds over the accepted range
ω ∈ [1.32, 1.47], **ungauged**. It is not established outside that range, and
the production regime (g = 0.05, Coulomb phase) needs a gauged BdG — until
then, GSS statements about production objects remain heuristic, exactly as
GROUNDING §1 caveat (i) requires.

**A distinction the catalog must carry.** The per-channel catalog frequencies
are **single-operator estimates** Ω(λ) = √(ω²+λ) − ω applied to L₀, L_x^sym and
L_x^flav separately. The true BdG frequencies come from the coupled quadratic
pencil, which is solved (and box-tested) only for the l = 0 symmetric sector.
Where the two disagree, the coupled box-tested result is authoritative — the
ω = 1.32 case above is exactly such a disagreement.

---

## 7. HC-2 — multipole-first resonance audit  **[D]**

Instrument `n_battery/hc2_resonance.py`, consuming the HC-1 catalog.
Classification order as the council corrected it: multipole **first**, then gap
arithmetic only for monopole/pure-matter modes.

* Every l ≥ 1 mode is assigned a golden-rule width against the massless gauge
  channel; gap arithmetic decides nothing for it.
* Every propagating **monopole** candidate has its fundamental **below** the
  matter continuum, so its leading radiative order is **n ≥ 2** (the harmonic
  crossing table gives n = 2 … >6 per mode).
* The census answer to *"how many harmonics can a particle hold?"* is
  **at most one, in a narrow window**: two box-stable bound internal modes
  across eight ω values, both in ω ∈ [1.33, 1.34].

Two consequences the program must absorb:

1. **HC-4's pre-registration tests the wrong power.** "Width ∝ amplitude²
   = golden-rule confirmation" assumes direct (n=1) emission. Direct emission
   is blocked for every monopole mode in the catalog; the leading process is
   n-quantum and the exponent must be derived per mode from the harmonic table
   before any measured width is called a confirmation. Combined with the
   council's box-IR hazard (GROUNDING §2 fix #4), **the halving list's
   demotion of HC-4 is earned, not merely economical** — and if HC-4 runs at
   all it should run at ω ≤ 1.34, where bound lines actually exist.
2. **Monopole radiative protection now has a linear-theory root.** A monopole
   mode cannot radiate through the gauge channel at linear order *and* has its
   fundamental below the matter continuum, so the leading perturbative opening
   is multi-quantum. That is the linear-theory root of the X10-series
   measurement that radial breathing persists. It is **not** a complete
   radiation catalog: parametric decay, multipoles induced by nonlinear
   self-interaction, gauge emission once motion breaks spherical symmetry, and
   sponge/box artifacts are all real channels in the nonlinear world X10
   actually lives in, and none of them is counted here. "Protection is derived
   at linear order" is the supported statement.

**PROPOSAL I.1's slogan, disambiguated.** "The number of independently stable
harmonics = the dimension of the conserved-charge lattice" has two readings and
they fare differently:

| reading | status after HC-1/HC-2 |
|---|---|
| a count of **bound internal modes** of a fixed soliton | **FALSE on the working branch** — the count varies with ω and is zero above ω ≈ 1.35, while dim(charge lattice) is fixed at 3 (+winding) |
| a count of **stable state families / clocks** labelled by the charges | not tested here; this is what the corrected GSS statement n(H) = n(D) formalises, and it is the reading the program should keep |

So the honest claim is: *if I.1 is read as a prediction for the linear
bound-mode census, it is false on the working branch; the load-bearing
replacement is GSS plus multipole/gap classification, not dim(charge lattice).*
What has **not** been shown is that "no stable harmonics exist in the theory" —
only that there are no below-threshold internal BdG modes on the monochromatic
ball in the working window. The QRK lines, the flavoured splittings and the
observed breathing may still be weakly-above-threshold, box-protected, or
nonlinear structures; they are simply not linear bound modes.

---

## 8. N9 — the soft-window scan and the HIER bound  **[M]**

Instrument `n_battery/n9_softwindow.py` + `n9_synthesis.py`. 25 theories,
μ × {0.6…1.5}, κ × {0.5…2.0} about the standards, ungauged, VK-stable and
bound (E < mQ) points only.

The protocol's two gates were **re-assigned** on review (the 2× gate belongs to
fixed-Q comparisons across theories; the 10% gate to the along-branch window),
and the along-branch 2× test was pre-registered as pre-doomed on this grid.

| reading | measured | gate |
|---|---|---|
| (a) E/Q window, standard potential | **11.72%** | 10% → technically not blocked |
| (a) E/Q window, widest theory (μ×1.5, κ×0.5) | 42.92% | — |
| (b) E at fixed Q across all 25 theories | **≤ 1.21×** | 2× → **failed decisively** |

**Synthesis (the honest reading).** On the VK-stable branch of *any* theory in
this family, E/Q = ω(1+ε) with ω < m and ε at the percent level, so E/Q ≤
m(1+ε). The only way to buy a large single-sector mass ratio is to push ω_min
toward zero. Measured, that buys 1.15× for the frozen theory and 1.52× at the
softest corner probed. **Reading (a)'s "not blocked" is a mis-calibrated gate:
an 11.7% spread is not a hierarchy.**

**Gates and synthesis, kept separate.**

* *The literal FOUNDATIONS gates:* the 10% along-branch gate **does not fire**
  (11.72% > 10%), so on its own terms it does not declare HIER blocked. The 2×
  fixed-Q gate **fails decisively** (1.21×), so potential redesign is not an
  alternative to a second sector.
* *The structural argument (what actually decides it):* E/Q ≤ m(1+ε) on any
  VK-stable localized solution of this class, so a large single-sector mass
  ratio requires ω_min → 0 or a different gap m. Measured, retuning buys 1.15×
  (frozen theory) to 1.52× (softest corner probed).

**HIER verdict: the single sector is blocked for hierarchy purposes** — carried
by the structural bound plus the failed 2× redesign test, *not* by the 10%
gate, which is mis-calibrated for the question (an 11.7% spread is not a
hierarchy). Stating it this way avoids the goalpost-moving objection: one
adopted gate did not fire, and the verdict does not rest on it.
*Caveats: ungauged (generous — gauging shrinks the window); the bound is
airtight only inside this soliton class and this functional family — other
potentials, multi-centre/topological objects, a different m, or a second sector
are not excluded by this rung.*

---

## 9. N10 — shell-mode E = ωQ exactness  **[M, with interpretation]**

Instruments `n_battery/sfa_radial.c` (new) + `n10_shellmode.py`, on the X10c
GPU dataset (ball Q=267, g=0.05, opposite-charge cloud, 15 frames).

**The new tool.** `sfa_radial.c` mirrors `scp_sim.c`'s
`compute_energy_complex_gauge()` and `compute_charges()` term by term —
covariant link differences, compact plaquette magnetic energy, the same
periodic indexing — and bins everything radially with the charge density split
by sign. **Validated against the kernel**: on x10c frame 14 it reproduces
Q_phi = 266.6825 (diag 266.6824) and E_total = 393.6511 (diag 393.6508).

**The subtraction.** No matched ball-alone run exists at Q_N = 267, so the
protocol's option 2 (integrate outside a core radius) is used, made auditable
by the sign-resolved charge: the cut is scanned, and the ball's charge
contamination of the selected region is < 1% at every reported cut. E_elec is
**excluded** (it is the ball's Coulomb field, present with or without the
cloud) and its size is reported.

**The truncation term.** On a region r > r_cut the identity carries an exact
boundary term, Surf(r_cut) = −π r_cut² d(ρ₂)/dr, which is *measured*, not
fitted. Cuts are selected where |Surf| is smallest (the cloud's inner turning
point), and both raw and corrected deviations are reported.

**Result.** Over 15 frames, at the minimum-|Surf| ball-clean cut:

* |deviation| = 0.0023 … 0.0772, **median 0.026**
* the only nonlinear term, E_pot, is **≤ 0.010% of E_cloud**
* the deviation is **positive in 12 of 15 frames**

**Verdict: PASS-WITH-INTERPRETATION.** For *any* linear superposition with
occupations n_k at frequencies ω_k, ρ₂ = Σn_k, |Q| = Σω_k n_k, E = Σω_k²n_k and
the measured ω_cl = |Q|/ρ₂ = ⟨ω⟩, so

    E / (ω_cl |Q|) = ⟨ω²⟩/⟨ω⟩² = 1 + Var(ω)/⟨ω⟩²   ≥ 1.

For a linear cloud the residual would therefore *be* the relative frequency
variance, and would be necessarily non-negative.

**This is a hypothesis with partial support, not a completed interpretation,
and it must be reported as such.** In its favour: the cloud is linear to three
decimal places (E_pot ≤ 0.010% of E_cloud), the residual is positive in 12 of
15 frames, and X10c is known independently to breathe with a 500–600 t.u.
period and never to have reached stationarity. Against it:

* the measured **spatial** clock variance supplies only **~3%** of the
  residual — 97% of the proposed explanation is *unmeasured*;
* **3 of 15 frames have a negative deviation** (t ≈ 900, 1350, 1500). Pure
  linear spectral variance *forbids* E < ω_cl|Q|, so those frames are an open
  tension requiring one of: surface/definition error, residual core
  contamination, an omitted gauge contribution, or non-stationary bookkeeping.

The protocol asked for residual ≪ soliton ε as *exactness*; that is correctly
demoted here, and it must not then be re-promoted via an unmeasured variance.
What the data supports: **the cloud is a linear object** (measured), and its
2.6% residual is *consistent with* — not demonstrated to be — a frequency
spread. The distinguishing test is a direct spectral variance from a densely
sampled run (see §12); if the spectral Var comes out ≪ residual, this
interpretation is falsified.

---

## 10. EX-2 — sponge-flux spectroscopy  **[M]**

Instrument `n_battery/ex2_spongeflux.py` on the `sfa_radial` output. The
measurement surface is r = 49.9, just inside the sponge (R_damp = L − 3 = 52).

Frequency content is recovered from **spatial** structure, because the archived
snapshot cadence (143 t.u.) cannot Fourier-resolve ω ≈ 1.4. For a free matter
wave, ω_A = m√(E_kin/E_mass) and ω_B = m√(1 + E_grad/E_mass) are independent
estimates that must agree.

**Results — two runs, deliberately contrasting sources.** x10c is a parked
ball with a radially breathing cloud (a **monopole** source). x11pb is the
two-closure orbit (heavy flavoured ball + light soliton on a circular orbit) —
a genuine **rotating dipole**, the configuration that *should* radiate through
the gauge channel if anything does.

| quantity | x10c (monopole) | x11pb (orbiting dipole) |
|---|---|---|
| gauge share of escaping power, median | 1.4×10⁻⁵ | 2.4×10⁻⁵ |
| gauge share, max | 7.4×10⁻⁵ | 1.3×10⁻⁴ |
| ω_A | 1.434 … 1.639 | 1.464 … 3.544 |
| ω_B | 1.509 … 1.615 | 1.509 … 3.250 |
| free-wave consistency median\|ω_A−ω_B\|/ω_A | 0.013 | 0.013 |
| energy-budget closure (cadence-limited) | 114% | 78% |

The orbiting dipole raises the gauge share by only ~1.7×, leaving it four
orders of magnitude below matter carriage. (x11pb's frame 18 failed to read —
the archived file is short — and its frame 17 carries a corrupt timestamp;
neither affects the channel ratio.)

**Readings.**

1. **The exchange channel is matter.** Gauge (Poynting) carriage is four to
   five orders of magnitude below matter carriage in **both** runs. The three
   candidate channels of PROPOSAL II.2.3 are not co-equal: above-gap matter
   waves carry the exchange; the gauge energy present near the object is
   *static Coulomb near field*, not radiation.
2. **Why, structurally:** radiating through the massless channel needs a
   time-varying l ≥ 1 multipole. A parked ball with a radially breathing cloud
   is a monopole source. This is the same obstruction HC-2 derives from the
   linear spectrum — EX-2 now measures it in the escaping power itself. The
   x11pb contrast is the control: giving the source a real rotating dipole
   moves the gauge share by a factor 1.7, not by orders of magnitude, so the
   matter dominance is a property of the coupling strength (g = 0.05), not
   merely of the monopole selection rule.
3. **The escaping waves sit just above the gap** (ω ≈ 1.51 vs m = 1.5), and the
   two independent frequency estimates agree to 1.3%, so the outflow really is
   a free above-gap matter wave rather than near-field leakage.
4. **Consequence for EX-1:** the radiative losses of a boosted object should be
   dominated by above-gap matter waves. That is now a pre-registerable
   expectation rather than a guess.
5. **What 114% does and does not mean.** With 15 samples 143 t.u. apart the
   trapezoid cannot resolve bursts, and the measurement surface sits 2 units
   inside the sponge. 114% is therefore **order-of-magnitude consistency — no
   missing channel at the factor level** — not precision energy accounting.
   Mild over-closure is expected from trapezoid bias and near-field
   contamination at the surface. The alarm conditions (wrong sign, or a factor
   off) did not fire.
6. **On the 1.4×10⁻⁵ gauge share.** The conclusion "gauge ≪ matter" is robust:
   an O(1)–O(10) error in the plaquette B normalisation would not touch it,
   only a ~10⁴ error would. The number itself should not be quoted as a
   precision constant without a lattice-unit audit.

*Limit (pre-registered):* this rung delivers a mean frequency per shell and a
matter/gauge split, **not** a resolved line spectrum. A resolved spectrum needs
dense snapshots or an in-kernel flux diagnostic — a new campaign, not a
re-analysis.

---

## 10b. N7 — the inertia lock (D5′)  **[M]**

Design **N7-STRESS** (grok-4.5 seat). Instruments: `n_battery/sfa_momentum.c`,
`n7_inertia.py`, `make_profile.py`. Scout run: two same-sign gauged Q-balls
(gauged-shooter profile ω = 1.430, g = 0.05, Q = 221 continuum → 211.5 on the
lattice), anti-phase, centres at x = ±5, N = 96, L = 24, absorbing sponge,
T = 90. Charge drift −4×10⁻⁴, energy drift −8×10⁻⁴.

**Nothing in the measurement chain presupposes a relation between E, Q and M:**
X(t) is the charge-weighted centre of each half-space, P(t) = ∫T^{0i} over it,
F(t) = −∮T^{ji}n_j dA through the midplane. The analytic g²Q₁Q₂/4πD² is never
used. (The review's original circularity flag does not apply — g is a kernel
input, not fitted from energies — but the analytic force's lattice form-factor
error is the same few percent as ε, which is exactly the gap being resolved.)

**The built-in validation did its job.** The momentum balance dP/dt vs the
measured flux resolved the surface-normal sign automatically (residual 0.079 vs
2.079) and then exposed a real systematic: **∫F dt = 18.457 against ΔP = 19.299,
so the plane-flux integral runs 4.4% low** — a discretisation artifact of
standing a single voxel slab in for a surface with a staggered plaquette B.
dP/dt is therefore the calibrated force, and any F/a estimator inherits the
bias. This is precisely the failure the R_PF gate was designed to catch.

**Estimators** (the first three use only P and X — one derivative each):

| estimator | M | M/E − 1 | M/(Qω) − 1 |
|---|---|---|---|
| M_slope = dP/dv | 312.40 | **−0.0029** | +0.0330 |
| M_P = P/v | 315.03 | **+0.0055** | +0.0416 |
| M_dPa = (dP/dt)/a | 314.68 | **+0.0044** | +0.0405 |
| M_a = F/a *(excluded: known 4.4% flux bias)* | 305.75 | −0.0242 | +0.0110 |

Comparators from this run, not from the continuum shooter: Q = 211.49,
E = 313.31 (region energy minus its 0.14% share of interaction energy),
Qω = 302.43, and the hypothesis separation **ε = E/(Qω) − 1 = 0.0360**.

Derivative-window sensitivity (w = 5…21) moves the three momentum estimators by
under 1.5% and leaves the ordering untouched.

**N7 VERDICT: D5′ CONFIRMED.** |M/E − 1| = 0.29–0.55% while |M/(Qω) − 1| =
3.3–4.2% ≈ ε, on every momentum-based estimator. The inertia that resists a
*measured* force is the vacuum-subtracted energy, **M = E/c²** — not the
thin-wall skeleton Qω. H4 is dead as a mass definition, and the v85-era
demotion narrative is not revived.

*Systematic budget (pre-registered):* the lattice group-velocity anomaly (1–5%)
is the honest ceiling, so sub-percent agreement was neither expected nor
demanded; the scout has only 6.6 voxels per r_half and a 4.3% lattice charge
deficit against the continuum profile, which absolute masses inherit but the
*ratios* M/E and M/(Qω) largely do not, since E, Q and M are all measured on the
same lattice. A refinement pair (N=112 and N=160 at L=30, D=12, T=130) is
running to close the dx axis and to serve as D7-lite on this observable.

---

## 11. What this changes in the standing documents

| Document | Required change |
|---|---|
| `GROUNDING.md` §0, `THEORY_v86.md` A4, FOUNDATIONS R2 | Σ's gauge sign: **−E_g**, not +E_g (or state it against the throughput ωQ/2). Add the closed form Σ = (2/3)(E_∇ − E_g). |
| `THEORY_v86.md` A4 | ε's "three closed-form targets under active derivation" resolve to one: the virial identity delivers Σ = (2/3)(E_∇ − E_g) exactly. ε itself remains a computed profile functional; the thin-wall exponent (≈1, local fit, R-definition sensitive) is evidence for the surface-law target, not a closed σ_eff. |
| `THEORY_v86.md` A4 / A9 | ħ_eff/Q − 1 and ε are the **same** residual, measured (N6). |
| `PROPOSAL.md` §I.1 | The stable-harmonics slogan is **falsified**, not merely superseded; the m=1.5 continuum figure should be Ω_c = m − ω. (Grok Finding 6 asked for this rewrite; the measurement now forces it.) |
| `PROPOSAL.md` Part I.4 / Amendment 1 | HC-4 must be re-designed around n ≥ 2 leading order, or run at ω ≤ 1.34, or dropped per the halving list. |
| `THEORY_v86.md` Part C | N9 closes the "retune the potential instead of adding a sector" option; HIER's second-sector requirement is now measured, not argued. |
| `THEORY_v86.md` A4, Part C locks | **D5′ is decided: M = E/c² measured to 0.3–0.6%, against 3.3–4.2% for Qω.** The lock can be closed. |
| `THEORY_v86.md` A8 | n(H_ω) = 1 is measured over the box-adequate ungauged range; the GSS statement's ungauged anchor is in place *there*. The gauged (production) case stays heuristic until a gauged BdG exists. |

## 12. What is still open

* **N7 (inertia lock / D5′)** — the next rung, and still before any EX-1 boost
  interpretation. Design settled as **N7-STRESS** (grok seat, 2026-07-26): a
  two-ball Coulomb *setup* but with a **measured** force (Maxwell/matter stress
  flux through a plane between the balls, `sfa_momentum.c`) and a **measured**
  momentum (∫T^{0i} per half-space), giving two independent inertia estimators
  M = F/a and M = |P|/|v| plus the momentum-balance residual
  R_PF = |dP/dt − F|/|F| as a built-in kill criterion. The analytic
  g²Q₁Q₂/4πD² is demoted to a scale check: it is not circular via a fitted g
  (g is a kernel input), but its lattice form-factor and image errors are the
  *same few percent as ε*, which is exactly the gap between the hypotheses
  M = E and M = Qω. A force-model-free **ratio arm** (unequal Q, M_A/M_B =
  |a_B|/|a_A|) rides along. Scout running: N=96, L=24, D=10, Q=221 each,
  anti-phase, g=0.05, T=90.
* **N8 (ring mass ladder)** — deferred per the halving list; non-blocking.

* **HC-6 (converse decay)** — needs seeded GSS-unstable partitions; HC-3 has
  now identified that n(D) = 1 everywhere on the scanned ray, so HC-6's targets
  must be sought outside it.
* **HC-4** — see §7; redesign required before it is worth running.
* **EX-1/EX-3** — GPU.
* **Gauged BdG** — HC-1 is the ungauged anchor; the Coulomb-phase case needs
  the A₀ perturbation solved self-consistently. Until it exists, **every GSS
  statement about production objects (g = 0.05) remains heuristic**, and this
  is the weakest brick in the Part-0 wall.
* **HC-1 completion debt (CPU, cheap)** — the ω range must be extended with
  box-converged runs at the thin-wall end, and the bound-mode count reported
  only where the Goldstone check passes. Rejected ω values are listed in
  `hc1.log` and quoted nowhere.
* **N10 spectral confirmation** — the Var(ω) reading of the 2.6% residual is a
  hypothesis until a densely sampled run gives a direct spectral variance, and
  until the three negative-deviation frames are explained. This is the test
  that decides between "interpretation" and "excuse".
* **HC-3 ray limitation** — n(D) = 1 everywhere on the single detuning ray
  scanned, so HC-3 currently maps a stable *tube*, not the partition space.
  HC-6's decay targets must be sought **off** that ray (large detuning,
  asymmetric charge partitions, multiple Q), or HC-3 must be extended to a
  (Q₀,Q₁,Q₂) volume scan.
* **EX-2 resolution debt** — a resolved line spectrum needs dense snapshots or
  an in-kernel flux diagnostic; the present rung is order-of-magnitude.
