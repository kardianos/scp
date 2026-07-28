# v86 GROUNDING (v2, post-council) — mathematics under the census + the frozen foundation
**Date:** 2026-07-26 · **v1 reviewed by the three-seat council** (see
`council/SYNTHESIS.md`); this revision applies every converged fix and
INTEGRATES the adopted Grok FOUNDATIONS frame (R1/R2/R4/R5). Where v1 is
contradicted, this file supersedes it; v1's text survives in git history.

## 0. The frozen foundation (adopted from FOUNDATIONS R1/R4 — supersedes TOE Step 3's E=Qω-first phrasing)

- **Primary objects:** the vacuum-subtracted Hamiltonian energy **E**, the
  conserved charges **Q_a** (+ winding n), and the branch relations. Inertia
  is **M = E/c²** — E = mc² stands as the honest rest-energy relation, not
  as a "derived shadow" of E = Qω. (v85's elevation of E=Qω repeated v84's
  mistake in mirror image: it closed the D5 units gap cosmetically by
  claiming the closed formula that CREED A4 had rightly refused.)
- **ω_a = ∂E/∂Q_a are chemical potentials** — the Legendre structure is the
  exact, measured content (8×10⁻⁵). "E = Qω" is the monochromatic thin-wall
  skeleton, exact only for bound linear modes.
- **Σ ≡ E − Σ_a ω_a Q_a is the virial/surface excess** — the object formerly
  known as ε(Q)·ωQ. **DERIVED by N1/N2** (`v86/PART0_RESULTS.md` §1–2); this
  bullet is corrected here from its v1 wording:
  - **Sign correction.** With the *true* integrated kinetic energy
    E_kin = ∫(3/2)wt²f²dV, the excess is
    **Σ = E_∇ + E_m + E_V − E_kin − E_g** — a MINUS on the gauge term, not the
    "+ gauge" of the v1 text. The two are reconciled by the exact Gauss
    identity **ωQ = 2(E_kin + E_g)** (residual ≤ 5×10⁻⁹ of E): if "kinetic" is
    read as the *throughput* ωQ/2 then Σ = E_∇ + E_m + E_V − ωQ/2 and no bare
    E_g appears at all. Mixing the two conventions is an error of 2E_g.
  - **Closed form.** The Derrick/virial identity of the profile system,
    (E_∇ − E_g) + 3(E_m + E_V − E_kin) = 0, reduces the excess to
    **Σ = (2/3)(E_∇ − E_g)**, verified to ≤ 5×10⁻⁸ of E on 292 branch points
    across g ∈ {0, 0.02, 0.05, 0.10}, on the flavoured branch, and across 25
    (μ,κ) potentials. Ungauged this gives Σ = (2/3)E_∇ > 0 identically, which
    *derives* the universally positive ε.
  - **What is still not closed.** ε remains a computed profile functional —
    E_∇ must come from the shooter. The thin-wall surface-law target is
    supported (local fit ε ~ (ξ/R)^1.01, ε' collapses on ξ/R across g after the
    measured-E_g subtraction) but σ_eff is not fitted, and the collapse
    abscissa is O(1) sensitive to the definition of R.
- **Three-way split even within one sector:** gauge charge Q_G ≠ process
  action 𝒜 ≠ inertia M. "ħ_eff = Q" is demoted from identity to measured
  ratio ħ_eff/Q = 1+δ. **N6 closes the unification this bullet predicted:**
  ħ_E/Q − 1 ≡ E/(ω̄Q) − 1 is *identically* ε̄, measured at 0.026–0.049 across
  45 archived single-object runs — v70's 3–5% and the shooter's 1–4% are the
  same residual measured twice. The third leg ħ_pk = p/k still needs a boost
  series. N5 also measures the local-clock/throughput split directly: the
  ball-centre clock runs 0.2–4.1% slow against the throughput-derived bare ω,
  g-dependently.
- **Q-degeneracy consequence:** degeneracy is a *regime property of
  monochromatic solitons*, not destiny. Foundation-level attacks (ring
  mass-ladder M(n) at fixed Q; potential-window scans) run before/alongside
  the second-sector design, which stays motivated by carrier physics and
  universal ħ — not by rescuing a false fundamental law. Cloud-mass does
  not solve HIER; it restates it.

## 1. Multi-charge stability: GSS, corrected (council fix #1 — the sign error)

Setup as before: three commuting U(1)s at η=0; stationary states extremize
E − Σω_a Q_a; d(ω⃗) = E − Σω_a Q_a on the branch; ∂d/∂ω_a = −Q_a.

**v1 stated the criterion backwards** ("n(H_ω) = p(D)" with D = ∂Q_a/∂ω_b),
which fails its own single-charge reduction: the VK-stable branch has
dQ/dω < 0, i.e. p(D) = 0, while n(H_ω) = 1. **Corrected statement:** with
Hess d = −∂Q_a/∂ω_b, GSS orbital stability requires
**n(H_ω) = n(∂Q_a/∂ω_b)** — the count of NEGATIVE eigenvalues of the
charge-response matrix must equal the constrained negative index of the
linearized energy. Single-charge check: dQ/dω < 0 gives one negative
eigenvalue = n(H_ω) = 1 ✓ (VK recovered).

**What HC-3 computes:** D_ab = ∂Q_a/∂ω_b by finite differences on the
flavored branch; record its full signature; the stability claim is the
matching of negative indices — and requires n(H_ω) from the HC-1 linear
solve, which is now provisioned (council fix: v1 had no rung for n(H_ω)).

**MEASURED (v86 Part 0).** HC-1 gives **n(H_ω) = 1** at every box-adequate ω in
[1.32, 1.47], ungauged: L₀ contributes no negative directions (its zero mode is
the nodeless U(1) Goldstone, recovered at |λ₀| ≤ 10⁻¹¹ with eigenvector overlap
1.000000 against r·f), L_x^flav = L₀ − A with A = 2P₀ < 0 contributes none, and
L_x^sym contributes exactly one — the l = 0 dilational direction. HC-3 finds
**n(D) = 1** across every scanned flavoured partition. The indices therefore
match over the scanned region, in the ungauged theorem regime. Two limits stand:
the HC-3 scan is a **single detuning ray**, so it maps a stable tube rather than
the partition space (HC-6's decay targets must be sought off it), and the
production regime g = 0.05 is Coulomb-phase, where a gauged BdG does not yet
exist and GSS remains heuristic.

**Caveats (hardened per council):** (i) the gauged system is
COULOMB-phase — the massless channel reaches ω = 0, the cited vortex
results are Higgs-phase, so gauged GSS is used as a heuristic with the
ungauged theorem as anchor; verdicts are labeled accordingly. (ii) GSS is
orbital, not asymptotic — monopole-protected breathing coexists with
orbital stability (measured, v85).

## 2. Resonance decay: multipole first, then arithmetic (council fix #3)

**MEASURED (v86 Part 0), and one figure in this section is wrong.** The
matter continuum edge for the coupled (x,y) system in the CO-ROTATING frame is
**Ω_c = m − ω**, not √(m²−ω²) and not the bare m: the asymptotic dispersion
(Ω²−G)² = 4ω²Ω² with G = m²−ω²+k² gives Ω = −ω + √(m²+k²). This is the
lab-frame statement "a perturbation at co-rotating Ω radiates iff ω+Ω ≥ m",
i.e. exactly the band [m, ω_max] this section names — but the number is a
factor ~6 smaller than √(m²−ω²) at ω = 1.42, and using the wrong one promotes
box-discretised continuum states to "bound internal modes". The single-operator
eigenvalue maps as Ω(λ) = √(ω²+λ) − ω, sending λ_c = m²−ω² to Ω_c exactly.
HC-2's harmonic arithmetic runs against Ω_c.

**Classification order corrected:** for every internal mode, FIRST classify
its multipole content. Any mode with a time-varying multipole l ≥ 1
couples to the **massless gauge channel, open at all frequencies** — gap
arithmetic decides nothing for it; its width is a golden-rule computation
against the A-channel density of states (κ_box-limited in finite volume).
Only monopole (and pure-matter) excitations get the gap-arithmetic
treatment vs the matter band [m, ω_max], ω_max² = m² + 12/dx² (formula
council-verified). The lattice band-top re-protection warning stands, with
Grok's addendum: band-edge Van Hove structure makes near-edge rate
predictions unreliable in either direction — refinement tests are
mandatory for any near-edge verdict.

**Width-floor calibration (Kimi, adopted):** before any HC-4 width is
quoted, measure the numeric width floor on a mode with NO open channel
(monopole, below-gap, below-cutoff) — "zero width" is then calibrated, not
assumed.

**Box-IR hazard (council fix #4):** the QRK-1 lines (0.018–0.126) all sit
below the box IR cutoff (≈0.209 at N=64/L=15; ≈0.057 at N=384/L=55) — in
any proposed box those modes have no open radiative channel, and nonlinear
frequency-pulling can fake an ε² "confirmation." HC-4 therefore (a) runs
the width-floor calibration, (b) uses a box pair chosen to move the cutoff
across at least one line, and (c) treats ε²-scaling WITHOUT cutoff-crossing
sensitivity as NULL, not confirmation.

## 3. KAM: framing only (unchanged, council-verified).

## 4. Claims, restated at provable strength
C1 (theorem-shaped): a configuration is orbitally stable iff its
sector-constrained GSS indices match (corrected §1), and every other
localized excitation carries a width that is either golden-rule computable
(open channel) or bounded by the calibrated floor (blocked channel), with
multipole classification deciding which — lattice band-top blocks flagged
as artifacts. C2 (correspondence with real physics' lowest-state-per-
sector rule) stays interpretive, [C]-tagged. **HC-1..HC-6 as amended are
judged by the council sufficient in structure for C1**, with HC-4's
hardening load-bearing.

## 5. The ε/foundation numeric battery (FOUNDATIONS N1–N10, adopted as Part 0 of the program)
All shooter/analysis-level (zero GPU spend) except N7:
- **N1** component decomposition of E on existing gscan + branch tables →
  Σ split into E_∇, E_m, E_V, E_g vs Q (is Σ the virial imbalance as
  claimed?).
- **N2** integrated virial identity of the radial ODE (parameter-free) →
  ε(Q) as a derived function; TOE-grade if it lands.
- **N3** scaling collapse: ε vs ξ/R across g after Coulomb subtraction.
- **N4** multi-flavor Legendre excess on the flavored branch.
- **N5** throughput-vs-charge (𝒜 vs 2πQ) on flavored/breathing objects —
  the action≠charge split made empirical.
- **N6** ħ_eff triple comparison (reopen v70 calibration against ε).
- **N7** N-INV weak-acceleration inertia measurement — decides D5′
  empirically (M_inertial vs E/c² vs Qω/c²). GPU-light.
- **N8** ring spin-ladder M(n) at fixed Q (v73 ring machinery) — mass
  without charge; the first foundation-level degeneracy crack test.
- **N9** soft-window shooter scan (μ,κ) — how much can M/Q move in one
  sector?
- **N10** shell-mode E = ωQ exactness on X10 profiles (the linear-mode
  anchor, atom-facing).

## 6. EX-1 adiabaticity, corrected (council fix #5)
The v1 threshold dropped the mass factor. Boost energy per unit charge is
½·(M/Q)·v² ≈ ½·m·v² (E/Q ≈ m for the cloud), so the survival condition is
**v ≪ √(2ε₁/m) = α_f ≈ 0.053** — v = 0.05 sits at 89% of threshold with a
pre-registered sudden-kick stripping estimate of ~55% (13% at v = 0.02),
and Kimi's differential-drift mechanism strips below the bound as well.
**Protocol:** ramped boost (≥4 steps, spacing ≫ 1/ε₁) mandatory above
v = 0.02; sudden-kick at v = 0.02 retained deliberately as a stripping-
fraction measurement. EX-4 retargets clock transport to the BALL clock
(measurable at T=2000), cloud-mode transport deferred to a longer run.
