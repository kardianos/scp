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
- **Σ ≡ E − Σ_a ω_a Q_a is the virial/surface excess** — the object
  formerly known as ε(Q)·ωQ, now to be DERIVED (see §5 battery): gradient +
  potential + gauge − kinetic imbalance, with closed-form targets (virial
  identity of the profile ODE; thin-wall surface law
  Σ ≈ 4πR²σ_eff(ω) + g²Q²/8πR; scaling collapse in ξ/R).
- **Three-way split even within one sector:** gauge charge Q_G ≠ process
  action 𝒜 ≠ inertia M. "ħ_eff = Q" is demoted from identity to measured
  ratio ħ_eff/Q = 1+δ (v70's 3–5% and ε's 1–4% are the same residual
  family, to be unified by the battery).
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

**Caveats (hardened per council):** (i) the gauged system is
COULOMB-phase — the massless channel reaches ω = 0, the cited vortex
results are Higgs-phase, so gauged GSS is used as a heuristic with the
ungauged theorem as anchor; verdicts are labeled accordingly. (ii) GSS is
orbital, not asymptotic — monopole-protected breathing coexists with
orbital stability (measured, v85).

## 2. Resonance decay: multipole first, then arithmetic (council fix #3)

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
