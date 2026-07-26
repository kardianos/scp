# v86 Council Feedback — Seat: kimi_k3
**Date:** 2026-07-26 · **Scope:** v86/PROPOSAL.md + v86/GROUNDING.md, against
v85/STATE_OF_THE_THEORY.md and spot-checks in v85/QRK_RESULTS.md,
v85/X10_RESULTS.md, sfa/sim/scp_sim.c.

Verification notes up front (so findings can be weighed):
- QRK-1 line values (0.018/0.036 family; flavored 0.108 at 8%, 0.126 at 4%):
  **verified** verbatim, v85/QRK_RESULTS.md:28–31.
- ε₁ = 2.1e-3, 1/ε₁ = 473 t.u. for the X10c cloud: **verified**,
  v85/X10_RESULTS.md:104,114.
- Matter gaps: m_Φ = 1.5, m_θ = 1.6: **verified** (v85/TOE_FEEDBACK_v85.2.md:71
  and elsewhere). Note there are *two* matter gaps plus the massless gauge
  channel; GROUNDING writes "the gap m" throughout.
- Lattice Laplacian: standard 7-point stencil, `(... − 6·f)·idx2`
  (scp_sim.c:864–866, 1187–1194). The band-top formula in GROUNDING §2 is
  therefore exactly the KG lattice dispersion maximum: **verified**.

---

## A. GROUNDING correctness

**Finding 1 — The GSS signature criterion is stated with the wrong count;
internally contradictory.**
Attacked sentence: *"orbital stability holds when n(H_ω) — the number of
negative eigenvalues of the constrained linearized energy — equals p(D), the
number of positive eigenvalues of D."*
With D_ab defined one sentence earlier as D_ab = −∂²d/∂ω_a∂ω_b = ∂Q_a/∂ω_b,
the single-charge reduction of "n(H_ω) = p(D)" is: stable iff ∂Q/∂ω > 0
(p = 1 = n). That is the VK-**unstable** side. The very next sentence —
*"For the single-charge Q-ball this reduces to the familiar
Vakhitov–Kolokolov dQ/dω < 0"* — requires instead the **negative**-eigenvalue
count of D to equal n(H_ω) = 1 (since d″ = −∂Q/∂ω > 0 on the stable branch).
As written, the criterion declares the entire measured-stable Q-ball branch
unstable. Fix: either state n(H_ω) = q(D) (# negative eigenvalues of D,
equivalently # positive of d″), or redefine D with the opposite sign. This is
a one-line definitional repair, but it is load-bearing: if HC-3 implements
the printed formula it will invert the stability map.

**Finding 2 — The Gauss-constraint remark is directionally right, hand-waved
in content, and its anchor is off-target.**
Attacked sentence: *"(ii) gauge fields add constraint directions: the Gauss
constraint must be projected before counting n(H_ω) (standard for gauged GSS,
e.g. abelian Higgs vortex literature)"*
Three problems. (a) The document never states the operative simplification
available in THIS system: flavor redistributions at fixed total charge are
gauge-neutral, so HC-3's variations of (ω₁,ω₂,ω₃) should be taken at fixed
Q_tot = Q₀+Q₁+Q₂, and on that surface the A₀-mediated contribution to the
partition Hessian plausibly drops out. That argument is neither made nor
referenced — the reader cannot tell whether D_ab is the gauged or the
effectively-ungauged 3-component Hessian. (b) "abelian Higgs vortex
literature" is the wrong anchor (codimension-2 topological defects); the
relevant body is gauged Q-ball classical stability, where dQ/dω < 0 is known
to be necessary but not always sufficient once electric-energy terms enter.
CANNOT VERIFY that "standard for gauged GSS" covers the gauged 3-charge case
as claimed. (c) Nothing says how the projection is implemented in the
shooter/Hessian pipeline. Fix: state the fixed-Q_tot variation rule
explicitly and either prove the flavor-sector decoupling from the constraint
or carry it as an explicit error budget on HC-3.

**Finding 3 — The golden-rule amplitude scaling is asserted for the wrong
regime; two decay laws are conflated.**
Attacked sentence: *"measured linewidths of the QRK-1 lines should scale as
the square of the excitation amplitude at lowest coupling order (Γ ∝ ε² for
a first-order-allowed channel, higher even powers when the first channels are
multipole- or arithmetic-blocked)."*
For a channel allowed at the linear level (internal mode × background →
continuum), the Soffer–Weinstein rate is a spectral property of the
linearization: the mode decays **exponentially with amplitude-independent
Γ**. The linewidth does not scale with ε in that regime. Γ ∝ ε² is instead
the signature of a self-harmonic channel (source ∝ ε² at 2ω_int, radiated
power ∝ ε⁴, mode energy ∝ ε²), and there the decay is **power-law**
(|a|² ~ 1/t), not exponential — so a "linewidth" is fit-window-dependent and
the HC-4 protocol ("long-T for narrow lines") is not even well-posed without
specifying the windowing. The sentence is right for one regime and wrong for
the one it names. Fix: make the HC-4 prediction conditional and two-branched
— flat width + exponential decay for linearly-allowed channels; ε² effective
rate + power-law envelope for self-harmonic channels — and require a
decay-law discriminant (log-linear vs 1/t fit) before any width is quoted.
Which regime the QRK-1 lines occupy is CANNOT VERIFY until HC-2 runs; see
Finding 10 for why that threatens the rung.

**Finding 4 — Band-top formula correct; the artifact warning is
quantitatively misleading by omission and incomplete on rates.**
Attacked sentence: *"On the lattice the 'continuum' is a BAND: ω ∈ [m,
ω_max(dx)] with ω_max² = m² + 4d/dx² (d=3)."*
Formula verified against the kernel stencil. But: at the GPU centerpiece
resolution dx = 0.287, ω_max = √(2.25 + 12/0.0824) ≈ **12.2 ≈ 8.1·m**; with
ω_clock ≤ 1.48, only combinations with **n ≥ 9** land above the band top —
orders where any golden-rule coupling is presumably negligible. The
re-protection is near-vacuous at EX-1's dx and only bites at dx ≳ 1
(ω_max → 3.77 at dx = 1, n ≥ 3). As written the warning implies imminence it
does not have; it should carry a per-rung table n_cross(dx) (CPU campaigns
have no single dx on record — L is per-campaign — so this must be tabulated,
not asserted). Second, the warning is incomplete: inside the band the lattice
density of states (van Hove singularities, v_g → 0 at zone edges) distorts
golden-rule **rates**, not just the kinematic cutoff — a width measured near
a zone edge is not the continuum Γ even when the channel is open. Third, the
wrinkle is written for a single "gap m"; the system has three bands (Φ at
1.5, θ at 1.6, gauge massless with ω_max = √12/dx ≈ 12.1 at GPU dx) and
HC-2's arithmetic must use the right edge per channel.

**Finding 5 — The multipole caveat is honest but weaker than stated;
"EXCEPT monopole breathing" rests on an unstated centeredness assumption.**
Attacked sentence: *"A mode combination below the matter gap may still
radiate through A if it carries a time-varying multipole — EXCEPT monopole
breathing (no dipole; v85's measured protection)."*
The gauge channel is massless: any mode with a time-varying l ≥ 1 moment of
the **total-charge** current radiates at its own frequency with **no
threshold and no integer-combination condition**. That is a much larger
allowed class than the matter-gap arithmetic suggests, and it makes the
multipole classification (correctly demanded two sentences later) mandatory,
not optional — see Finding 9 for the proposal-side consequence. The monopole
exemption was measured for a *centered* radial oscillation; a displaced or
composite radial mode (exactly what a boosted or multi-center configuration
produces, cf. v85 §4 "monopole-protected breathing… a structural annoyance
for every future bound-state experiment") carries a dipole moment and the
exemption lapses. The sentence should carry the centeredness qualifier,
because EX-1 will violate it.

**Finding 6 — A GSS hypothesis is missing from the caveat list: the kernel
of H_ω must be exactly the symmetry directions.**
The caveat list (i)–(iii) omits the nondegeneracy assumption. With a
phase-blind product potential, exchange-near-symmetric components generate
near-zero directions beyond the three U(1)s (v85: the clock partition is
invisible to density probes). Numerically, n(H_ω) then depends on an
eigenvalue cutoff, and the "GSS count" the whole rung architecture scores
against becomes threshold-ambiguous. HC-1/HC-3 need a stated null-space audit
with a cutoff protocol. Minor but cheap to pre-register.

---

## B. PROPOSAL design

**Finding 7 — No rung computes the count that C1 is made of.**
Attacked sentence (GROUNDING §1): *"The claimed 'stability region of
partitions' IS the region where the GSS count matches."*
The GSS count has two inputs: signature of D (HC-3 delivers) and n(H_ω), the
negative-eigenvalue index of the constrained linearized energy (**no rung
delivers**). HC-1 is scoped as *"the mode catalog the whole census scores
against"* — frequencies, not indices. The BdG linear response contains the
index (negative-norm / Krein count), but the proposal never asks for it. As
designed, HC-1..5 can establish at most an empirical correlation between
D-signature and observed stability, not the signature equality C1 asserts.
Fix: add spectral-index extraction (with the Finding-6 cutoff protocol) to
HC-1's deliverables. Without this, C1 is untestable by construction.

**Finding 8 — C1's "iff" and its sector scope both exceed the rungs.**
C1 (GROUNDING §4) quantifies stability over the sector *"flavor charges ×
winding"* and asserts an **iff**. Two gaps: (a) no rung varies winding —
v85 lists ring winding on top of clocks (L_z = nQ) as established physics,
but the census never tests a winding partition's signature or decay; C1 must
either be re-scoped to flavor charges or buy a winding rung. (b) The converse
direction (signature-violating ⇒ decays) is tested only by HC-5, which seeds
*"2-frequency content in ONE component"* — excess energy inside a single
flavor sector. That tests "everything else decays," a result already
foreknown from v64 3ω decay, and says nothing about the signature map across
partitions. The sharpest falsifiable test — seed a flavored partition that
HC-3 declares signature-unstable at fixed Q_tot and watch it decay to the
signature-stable partition — is absent (see C, M2).

**Finding 9 — HC-2 as tabled contradicts GROUNDING's own multipole
requirement.**
PROPOSAL HC-2: *"enumerate n·ω_i ± m·ω_clock vs gap for the HC-1 catalog →
predicted stable/decaying classification."* No multipole column, no gauge
channel. Per Finding 5, every mode with a time-varying dipole/quadrupole of
the total current is gauge-allowed at zero threshold, so the HC-2 "stable"
class will be contaminated by gauge-radiating modes exactly when modes sit
below the matter gap — the census's most interesting region. GROUNDING §2
item 2 explicitly demands the classification (*"The census must therefore
classify modes by multipole as well as frequency"*); the proposal table drops
it. Fix: HC-1 must output multipole labels of the total-charge current
moments per mode, and HC-2 must add a gauge-channel column.

**Finding 10 — EX-1 adiabaticity: arithmetic verified; the criterion is
mis-shaped and will mislabel outcomes.**
Attacked sentence: *"the cloud survives only if the imparted energy per unit
charge ~ ½v² ≪ ε₁ → v ≪ √(2ε₁) ≈ 0.065."*
Arithmetic checked: √(2·2.1e-3) = 0.0648 ✓; ½(0.02)² = 2e-4 ≈ ε₁/10.5 ✓;
½(0.05)² = 1.25e-3 ≈ 0.60ε₁ ✓ (per-charge normalizations differ by O(1)
factors of ω, immaterial here). The physics, however: ½v² ≪ ε₁ is the
condition for *all* translational kinetic energy of the cloud to thermalize
into the relative coordinate — a worst-case bound, not an adiabaticity
estimate. A common phase tilt injects co-moving momentum into core AND
cloud; the actual failure modes are (a) tilt ≠ true boosted eigenstate
(missing Lorentz contraction, missing A_μ boost) → an impulsive transient
whose relative-mode content is a small fraction of ½v², and (b) the measured
1–5% lattice group-velocity anomaly acting **differentially** on deeply-bound
core modes vs the shallow cloud mode → slow stripping in flight at velocities
where the impulse test passes cleanly. The single-number criterion will call
a drift-stripped run "impulse damage" or a clean-impulse run "safe" while it
strips at t ≫ 473. Fix: pre-register the differential-drift observable
(core–cloud separation vs t at fixed v, compared against the anomaly
prediction), and at v = 0.05 use a ramped tilt over τ ≫ 1/ε₁ = 473 t.u. —
the document's own closing question (*"whether a ramped (multi-step) boost is
needed at the upper end"*) should be answered "yes, and pre-registered,"
not left open.

**Finding 11 — HC-4 is the rung most likely to produce a misleading result.**
Three independent reasons. (a) Floor: v85's own problem list flags
*"BC-limited decay (pathfinder)"* as load-bearing and X2 refinement as never
run; a linewidth measured in one box at one dx with an absorbing sponge
cannot be attributed to golden-rule decay without an instrumental floor.
(b) Target: the QRK-1 lines sit at ω ≈ 0.018–0.126 ≪ m (verified values);
if these are the mode frequencies (not beats of near-gap eigenmodes), the
lowest self-harmonic channel is n ≥ 12 (12 × 0.126 = 1.51 > m) and every
combination channel is ultra-high-order — HC-4 would then measure the
sponge/discretization floor and report it as a "width," and could even
"confirm" Γ ∝ ε² spuriously if the floor extraction is amplitude-sensitive.
(c) Regime: per Finding 3, the predicted scaling is wrong for a
linearly-allowed channel, so a true golden-rule line showing flat width
would be logged as a falsification. Hardening: (i) run a known-protected
mode (monopole breathing — v85 guarantees it does not radiate) through the
identical HC-4 pipeline to map Γ_floor(ω, box, sponge, fit-window);
(ii) condition targets on HC-2 identifying at least one line with an allowed
low-order channel; (iii) duplicate one point at a second dx; (iv) require
the decay-law discriminant before quoting any width.

**Finding 12 — EX-2 feasibility is assumed, not established.**
CANNOT VERIFY that existing SFA frames contain the gauge links/A_μ needed to
decompose sponge-shell outflow into (matter/gauge) × frequency bins; v85
lists per-component Q_a as "SFA-frame analysis only" (frames are
scalar-field-centric) and records async SFA tail corruption for x11pb at
t ≥ 1800 — precisely the late-time bins a transfer census wants. If links
are not in the frames, EX-2 is a re-run with an instrumented shell, not
"analysis + reuse"; the cost table should say so now.

---

## C. Missing analyses (three, in priority order)

**M1. Instrumental-floor calibration for all linewidth claims.** A dedicated
protected-mode run (monopole breathing) through the full HC-4 pipeline —
same box, sponge, diag cadence, fit windows — producing the floor curve
Γ_floor(ω, dx) that every census width must exceed, plus one HC-4 point
duplicated at a second dx. v85 calls box/sponge systematics load-bearing and
the D7 4-axis certification "arguably the most overdue chore in the
program"; v86 Part I is a width-measuring campaign with no width floor. This
is the cheapest rung in the program and gates the credibility of HC-4/HC-5.

**M2. Converse-direction partition rung (call it HC-6).** Seed a flavored
partition inside HC-3's signature-unstable region at fixed Q_tot and measure
decay to the signature-stable partition. This is the only design that can
falsify the "stable ⇐ signature" half of C1's iff; HC-5 cannot (Finding 8).
If winding stays in C1's sector, the same rung at n = 1 covers Finding 8(a).

**M3. Gauge-channel multipole audit as an explicit HC-1/HC-2 deliverable.**
For every catalog mode: the l = 1, 2 moments of the total-charge current and
a gauge-allowed/blocked classification, folded into the HC-2 stable/decaying
table. Without it the "stable" class is contaminated at zero threshold
(Findings 5, 9), and C1's clause *"unless blocked by multipole selection"* is
never actually exercised.

---

## D. Verdicts

**GROUNDING.md: SOUND-WITH-FIXES.**
Required fixes:
1. Correct the signature criterion: n(H_ω) equals the # of **negative**
   eigenvalues of D as defined (Finding 1).
2. State the Gauss-constraint handling concretely: fixed-Q_tot partition
   variations; prove or budget the flavor-sector decoupling; drop or replace
   the vortex anchor (Finding 2).
3. Re-issue the HC-4 prediction as a two-regime conditional (exponential +
   flat Γ vs power-law + Γ ∝ ε²) with a decay-law discriminant (Finding 3).
4. Quantify the band-top warning per dx (n_cross table) and add the in-band
   lattice-DOS rate distortion; write the wrinkle for all three bands
   (Finding 4).
5. Add the centeredness qualifier to the monopole exemption and the
   null-space/cutoff caveat to the GSS hypothesis list (Findings 5, 6).

**PROPOSAL.md: SOUND-WITH-FIXES.**
Required fixes:
1. HC-1 deliverables must include the spectral index n(H_ω) (with cutoff
   protocol) and per-mode multipole labels of the total-charge current
   (Findings 6, 7, 9).
2. HC-2 must add the gauge-channel/multipole column and per-dx band-top
   flags (Findings 4, 9).
3. HC-4 must be gated on M1's floor calibration, target lines with allowed
   low-order channels per HC-2, and fit decay law before width (Findings 3,
   10, 11).
4. Add the converse partition rung (M2) or re-scope C1 to what HC-1..5 test;
   resolve the winding-sector scope explicitly (Finding 8).
5. EX-1: ramped boost at v = 0.05 with τ ≫ 473 t.u. pre-registered, and the
   differential-drift observable added to the pass criteria (Finding 10).
6. EX-2: verify gauge-field content of existing SFA frames before claiming
   "analysis + reuse"; otherwise re-budget as an instrumented re-run
   (Finding 12).

Neither document is unsound: the framework choices (GSS for the census,
Soffer–Weinstein for widths, band-top flagging) are the right tools and the
honesty caveats are mostly real. But Finding 1 is a sign error in the
program's central criterion, and Findings 7–9 mean that **as designed, the
rungs cannot establish C1 as stated** — both must be repaired before the
census runs, not after.
