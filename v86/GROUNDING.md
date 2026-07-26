# v86 GROUNDING — the mathematics under the harmonic census
**Date:** 2026-07-26 · **Role:** paper-math grounding for v86/PROPOSAL.md
Part I, written BEFORE the numerical rungs so the analysis has a precise
target and cannot overclaim. Reviewed by the v86 council (this file is a
review input).

## 1. Multi-charge stability: the Grillakis–Shatah–Strauss frame

**Setup.** The gauged complexified Cosserat system is Hamiltonian with (at
η=0) three commuting U(1) symmetries and conserved charges Q_a
(code-verified, v85/ANALYSIS §2E). Stationary states are constrained
extremizers: δ(E − Σ_a ω_a Q_a) = 0, so along any smooth family
**ω_a = ∂E/∂Q_a** (measured to 8×10⁻⁵ in the single-charge case).

**The stability object.** Define the dual d(ω⃗) = E − Σ_a ω_a Q_a evaluated
on the branch. Then ∂d/∂ω_a = −Q_a and the Hessian
**D_ab ≡ −∂²d/∂ω_a∂ω_b = ∂Q_a/∂ω_b** is symmetric. The
Grillakis–Shatah–Strauss (GSS, 1987/1990) theory gives the criterion:
orbital stability holds when n(H_ω) — the number of negative eigenvalues of
the constrained linearized energy — equals p(D), the number of positive
eigenvalues of D. For the single-charge Q-ball this reduces to the familiar
Vakhitov–Kolokolov dQ/dω < 0 (with the standard one-negative-direction
count). **The multi-harmonic generalization is therefore not a guess: it is
the signature condition on the 3×3 matrix ∂Q_a/∂ω_b along the flavored
branch.**

**What HC-3 must actually compute.** The flavored shooter gives
Q_a(ω_1,ω_2,ω_3) on the branch; finite-difference D_ab; record its
signature over the partition space. The claimed "stability region of
partitions" IS the region where the GSS count matches. Honest caveats to
carry: (i) GSS assumes the linearization has no unstable essential
spectrum — on the lattice this must be checked against the band (see §2);
(ii) gauge fields add constraint directions: the Gauss constraint must be
projected before counting n(H_ω) (standard for gauged GSS, e.g. abelian
Higgs vortex literature); (iii) GSS gives orbital stability, not asymptotic
stability — a stable partition may still be surrounded by long-lived
breathing (which v85 measured; monopole protection is exactly an
obstruction to asymptotic relaxation).

## 2. Resonance decay: Fermi golden rule with a lattice wrinkle

**The continuum statement (Soffer–Weinstein 1999; Sigal; Buslaev–Perelman).**
A discrete internal mode ω_int of the linearization about a soliton decays
when some small integer multiple couples into the essential spectrum:
n·ω_int > ω_edge (here the gap m). The rate is Γ ∝ |⟨F_n, ψ_res⟩|² — the
golden-rule constant built from the n-th order nonlinear coupling and the
density of states at n·ω_int — provided the coupling is nondegenerate
(the "nonlinear Fermi Golden Rule" condition, generically nonzero but not
guaranteed; it must be computed or measured, not assumed).

**Applicability conditions for OUR system (the review should attack these):**
1. Multi-component coupling: the mode couples through the product potential
   V(s), whose expansion mixes components — the golden-rule constant is a
   matrix element over component structure, computable from the HC-1
   catalog.
2. Gauge sector: the essential spectrum includes the massless gauge
   channel. A mode combination below the matter gap may still radiate
   through A if it carries a time-varying multipole — EXCEPT monopole
   breathing (no dipole; v85's measured protection). The census must
   therefore classify modes by multipole as well as frequency.
3. **The lattice band-top wrinkle (novel, honest artifact warning).** On
   the lattice the "continuum" is a BAND: ω ∈ [m, ω_max(dx)] with
   ω_max² = m² + 4d/dx² (d=3). Combinations with n·ω_int ABOVE the band top
   cannot radiate either — the lattice re-protects high harmonics that the
   physical continuum would kill. Any HC-2 "stable" verdict from
   above-band arithmetic is a lattice artifact and must be flagged; the
   dx-dependence of ω_max makes this testable (an artifact-protected mode
   destabilizes under refinement).

**Prediction for HC-4:** measured linewidths of the QRK-1 lines should
scale as the square of the excitation amplitude at lowest coupling order
(Γ ∝ ε² for a first-order-allowed channel, higher even powers when the
first channels are multipole- or arithmetic-blocked). Deviation from
even-power scaling falsifies the golden-rule reading.

## 3. KAM: framing only

Kuksin/Wayne-type KAM for PDEs proves persistence of finite-dimensional
quasi-periodic tori for parameter sets of positive measure in SPECIFIC
(mostly 1D, bounded, non-resonant) settings. Nothing we run satisfies those
hypotheses; the citation earns its place only as the structural reason
"generic extra frequencies die" — it licenses no quantitative claim in
this program and must not be cited as if it did.

## 4. The correspondence claim, stated at provable strength

Claim C1 (fabric theorem-shaped): *in this system, a configuration is
orbitally stable iff it is a GSS-signature minimizer within its conserved
sector (flavor charges × winding), and every other localized excitation
carries a golden-rule width unless blocked by multipole selection or
arithmetic (with lattice band-top blocks flagged as artifacts).*
Claim C2 (correspondence, interpretive): this reproduces the structural
rule by which real physics organizes stability (lowest state per
quantum-number sector; resonances elsewhere). C2 is an analogy claim and
stays [C]-tagged; C1 is what HC-1..5 can actually establish, and the
council review should judge whether HC-1..5 as designed suffice for C1.

## 5. EX-1 adiabaticity (pre-registered arithmetic)

Boost concern: the X10c cloud binds at ε₁ ≈ 2.1e-3 (binding time 1/ε₁ ≈
473). A boost applied as a seed-level phase tilt is impulsive; the cloud
survives only if the imparted energy per unit charge ~ ½v² ≪ ε₁ →
**v ≪ √(2ε₁) ≈ 0.065**. The proposed v = 0.02–0.05 spans gentle (0.02:
½v² = 2e-4 ≈ ε₁/10) to marginal (0.05: ½v² = 1.25e-3 ≈ 0.6·ε₁). The review
should check this estimate and whether a ramped (multi-step) boost is
needed at the upper end.
