# v86 Council Feedback — Seat: Grok 4.5 (`grok45`)

**Date:** 2026-07-26  
**Inputs read (order required):** `v85/STATE_OF_THE_THEORY.md`, `v86/PROPOSAL.md`, `v86/GROUNDING.md`  
**Optional context used:** `v85/QRK_RESULTS.md`, `v85/X10_RESULTS.md` (ε₁/X10c numbers), `v85/ANALYSIS.md` (continuum / alignment prior), `v66/THEORY.md` (m=1.5, window)  
**Charge:** adversarial + constructive on GROUNDING correctness, PROPOSAL design, missing pieces, verdicts.

Ground rules observed: attack quoted sentences; no praise padding; CANNOT VERIFY where not checked against primary sources; findings numbered.

---

## A. GROUNDING correctness

### A1. GSS signature criterion — **sign / counting error (internal inconsistency)**

**Attacked sentence (GROUNDING §1):**
> "The Grillakis–Shatah–Strauss (GSS, 1987/1990) theory gives the criterion: orbital stability holds when n(H_ω) — the number of negative eigenvalues of the constrained linearized energy — equals p(D), the number of positive eigenvalues of D. … For the single-charge Q-ball this reduces to the familiar Vakhitov–Kolokolov dQ/dω < 0 …"

**Attacked definition (same §):**
> "D_ab ≡ −∂²d/∂ω_a∂ω_b = ∂Q_a/∂ω_b"

**Finding 1.** With *this* definition of D, the single-charge reduction does **not** work as written.

- Dual calculus is fine: d = E − Σ ω_a Q_a ⇒ ∂d/∂ω_a = −Q_a ⇒ −∂²d/∂ω_a∂ω_b = ∂Q_a/∂ω_b. So D = ∂Q/∂ω is correctly related to d.
- VK-stable single-charge branch: dQ/dω < 0 ⇒ D is a **negative** 1×1 matrix ⇒ **p(D) = 0**, **n(D) = 1**.
- Standard Q-ball action Hessian H_ω = δ²(E − ωQ) typically has **n(H_ω) = 1** (one negative direction), with kernel from U(1)/translations handled separately.
- Criterion "n(H_ω) = p(D)" then reads **1 = 0** on the VK-stable branch — it does **not** reduce to VK.

**Correct pairings (standard index form; either is fine, not both mixed):**

| Define | Stability match (under usual nondegeneracy) | Single-charge VK |
|---|---|---|
| D = ∂Q/∂ω | n(H_ω) = **n(D)** (negative counts) | n(H)=1, n(D)=1 when dQ/dω<0 |
| D′ = −∂Q/∂ω = d″ | n(H_ω) = **p(D′)** | n(H)=1, p(D′)=1 when dQ/dω<0 |

GROUNDING mixes definition D = ∂Q/∂ω with the p(D) side of the other convention.

**Finding 2.** The phrase "constrained linearized energy" for H_ω is **terminology-wrong or at least ambiguous**. In the usual GSS index formula, H_ω is the Hessian of the *augmented* action S_ω = E − ω·Q (Lagrange-multiplier form), **not** the second variation of E already reduced to the fixed-Q manifold. On the true charge-constrained manifold, orbital stability is local minimization: n(H_constrained) = 0 (modulo symmetries). Writing "n(H_ω) = p(D)" only makes sense for the *unreduced* action Hessian. As written, a reader who implements "constrained Hessian negative-count equals p(∂Q/∂ω)" will implement the wrong object.

**Finding 3.** Multi-charge statement "the multi-harmonic generalization is therefore not a guess: it is the signature condition on the 3×3 matrix ∂Q_a/∂ω_b" is **directionally right but incomplete** even after the sign fix:

- Signature of D alone is not sufficient; you need **both** n(H_ω) (Morse index of the action Hessian at that flavored stationary point) **and** the inertia of D, then match.
- Kernel dimensions: three global U(1)s + spatial translations (+ gauge orbit) must be quotiented or counted in the kernel before "n" is meaningful.
- For *ground-state* multi-charge Q-balls one often has the stronger sufficient condition that ∂Q_i/∂ω_j is **negative definite** when n(H_ω) equals the number of independent charges — but that is a special case of the index formula, not a free-standing "signature of D" test.

**Finding 4.** Winding / ring sector: PROPOSAL I.1 lists "ring winding n" in the conserved lattice; GROUNDING §1 and §4 claim C1 speak of "flavor charges × winding," but the GSS object is only written as a 3×3 in (ω_a). If n is an independent conserved charge with conjugate angular frequency Ω, the stability matrix is at least block-enlarged (4×4 or the appropriate (Q_a, L_z) Jacobian). **CANNOT VERIFY** from archive whether the program treats n as an independent GSS parameter on a joint branch; as written, GROUNDING under-specifies the sector it claims in C1.

### A2. Gauss-constraint projection — **directionally correct, hand-waved for this system**

**Attacked sentence (GROUNDING §1 caveat ii):**
> "gauge fields add constraint directions: the Gauss constraint must be projected before counting n(H_ω) (standard for gauged GSS, e.g. abelian Higgs vortex literature)"

**Finding 5.** True that physical Morse index is computed on the Gauss-law constraint surface (Coulomb gauge / Helmholtz projection of δA, δE; eliminate pure-gauge and longitudinal residual). **Not sufficient as written for this kernel:**

1. **Massless gauge continuum.** Abelian Higgs vortices are typically studied with a *massive* photon (broken U(1)) or in 2D Ginzburg–Landau energy minimizers. Here A is massless Maxwell on the lattice. Essential spectrum of the linearization reaches **down to ω = 0 in the gauge sector**. Caveat (i) ("no unstable essential spectrum — on the lattice this must be checked against the band") discusses the *matter* band starting at m; it does not say how a massless gauge continuum interacts with the GSS spectral assumption. That is the more dangerous essential-spectrum issue.
2. **What to project.** One sentence; no protocol: discrete Gauss residual already sits at 1e-13 in the nonlinear runs (STATE §1) — good for the background — but the *linearized* operator's domain, gauge fixing, and whether longitudinal E is eliminated or retained as a constrained pair are unspecified. Without that, HC-1 "BdG spectrum" is not an implementable specification for the gauged system.
3. **Literature gesture.** "e.g. abelian Higgs vortex literature" is not a derivation. Acceptable as a pointer; **not** a substitute for stating the physical subspace for *this* 3-component + U(1) lattice Hamiltonian.

**Verdict on Gauss remark:** not false; **hand-waved** relative to the load it is asked to carry for C1.

### A3. Golden rule (Soffer–Weinstein) applicability — **honest about genericity, weak on the operational continuum edge**

**Attacked block (GROUNDING §2 continuum statement):**
> "A discrete internal mode ω_int of the linearization about a soliton decays when some small integer multiple couples into the essential spectrum: n·ω_int > ω_edge (here the gap m)."

**Finding 6.** "here the gap m" is **ambiguous and likely wrong if ω_int means the QRK envelope frequencies**.

- QRK-1 (v85/QRK_RESULTS.md) reports lines **0.018 / 0.036** (symmetric) and **0.108 / 0.126** (flavored) from FFT of `s_max(t)` — these are **amplitude-modulation / internal-mode frequencies Ω**, not field oscillation frequencies near m.
- For a Q-ball in a co-rotating / BdG frame, the continuum edge for radiation of an internal mode is **not** "Ω > m = 1.5". Absolute lab frequencies of the carrier are already ~ω_clock ∈ (1.31, 1.5); the relevant gap is the distance from the discrete eigenvalue to the essential spectrum of the *linearized* operator about the soliton (often order **m − ω** or a gauge-channel threshold near 0 for multipole-allowed modes).
- Numeric check at X10c clock ω ≈ 1.41: m − ω ≈ **0.09**. Then Ω = 0.018 may sit below a matter continuum gap while Ω = 0.108 may sit above — *exactly* the interesting census — but only if HC-2 uses the **BdG continuum edge from HC-1**, not the vacuum gap m as a hard threshold on the raw QRK numbers.

PROPOSAL I.1 compounds the error:
> "every low-order integer combination n·ω_i ± m·ω_clock stays below the radiation continuum (gap m=1.5)"

Mixing (small Ω) with (± ω_clock) vs absolute m=1.5 without a specified frame is **not a well-posed arithmetic test**. This is the highest-leverage theory bug in the pair of documents.

**Finding 7.** Multipole × gauge channel (GROUNDING §2 point 2) is the right qualitative structure and matches STATE/X10c (monopole breathing radiatively protected). Remaining holes:

- Classification "by multipole as well as frequency" is demanded of "the census" but **HC-1/HC-2 never list multipole labels as deliverables** (see B).
- Selection rules for a multi-component, phase-blind product potential are **not** vacuum multipole rules; V(s) expansion couples components. Matrix elements can open channels that naive single-field multipole counting forbids, or suppress ones it allows. GROUNDING point 1 admits this; the plan never budgets a matrix-element computation beyond "catalog."
- **CANNOT VERIFY** that Soffer–Weinstein 1999 hypotheses (specific dispersive PDE classes, non-resonance, etc.) hold for gauged complex Cosserat on a lattice. The document correctly says the FGR constant is "generically nonzero but not guaranteed; it must be computed or measured." That honesty is fine; claiming C1 as "theorem-shaped" later overreaches the same caveat.

**Finding 8.** HC-4 prediction (GROUNDING §2):
> "measured linewidths of the QRK-1 lines should scale as the square of the excitation amplitude at lowest coupling order (Γ ∝ ε² …)"

- Even-power scaling is the right FGR fingerprint **if** Γ is the **energy/amplitude decay rate of a resolved linear mode**, not an FFT peak width of `s_max` under a single kick protocol.
- Monopole-protected lines (the 0.018 family is a prime suspect for radial breathing) should have **Γ ≈ 0** from radiation; any measured "width" is then dephasing, numerical dissipation, sponge coupling, or weak higher-multipole leakage — and will **not** obey clean ε² radiation scaling. Falsification clause is good; **null width is also a success mode** that the prediction text under-emphasizes.
- "higher even powers when first channels are multipole- or arithmetic-blocked" — OK as a sketch; operational distinction between Γ=0, Γ∝ε⁴, and Γ from lattice artifacts is not specified.

### A4. Lattice band-top formula — **algebra correct; warning incomplete**

**Attacked sentence (GROUNDING §2 point 3):**
> "On the lattice the \"continuum\" is a BAND: ω ∈ [m, ω_max(dx)] with ω_max² = m² + 4d/dx² (d=3)."

**Finding 9.** For the standard second-order nearest-neighbor Laplacian in d spatial dimensions,

\[
\lambda_{\max}(-\nabla^2_{\mathrm{lat}}) = \frac{4d}{\mathrm{d}x^2},
\]

so ω_max² = m² + 4d/dx² is **correct** for the free massive scalar on this lattice stencil. Numeric: at GPU dx ≈ 0.287, ω_max ≈ 12.2; at N=64 typical coarser dx, lower but still ≫ m. Harmless for low internal modes; relevant only for high harmonics / multi-linear combinations.

**Finding 10.** Artifact warning is **incomplete**:

1. **Gauge band is different.** Massless A: ω ∈ [0, 2√d / dx], not [m, √(m²+4d/dx²)]. A multipole-allowed mode can radiate into gauge **below the matter gap** and is constrained by a **different** band top. Flagging only matter-band overshoot misses half the lattice pathology (and half the physical decay channel).
2. **Background A and link-covariant derivatives** shift/deform the local dispersion relative to the free stencil; free-band ω_max is only a first estimate around a gauged soliton.
3. **Umklapp / zone-boundary** combination resonances on the Brillouin zone are not mentioned; they are the other lattice-only decay path.
4. dx-refinement test ("artifact-protected mode destabilizes under refinement") is the right experimental handle — keep it; extend it to **gauge-band** and **continuum-edge** refinement, not only ω_max(matter).

### A5. EX-1 adiabaticity arithmetic — **√(2ε₁) correct; mass-per-charge and protocol caveats**

**Attacked block (GROUNDING §5):**
> "the cloud survives only if the imparted energy per unit charge ~ ½v² ≪ ε₁ → v ≪ √(2ε₁) ≈ 0.065. The proposed v = 0.02–0.05 spans gentle (0.02: ½v² = 2e-4 ≈ ε₁/10) to marginal (0.05: ½v² = 1.25e-3 ≈ 0.6·ε₁)."

**Finding 11.** Arithmetic check (ε₁ = 2.11×10⁻³ from X10c):

| Quantity | Value | Match to GROUNDING? |
|---|---|---|
| √(2ε₁) | 0.06496 ≈ 0.065 | YES |
| ½(0.02)² | 2×10⁻⁴ ≈ ε₁/10.6 | YES (~ε₁/10) |
| ½(0.05)² | 1.25×10⁻³ ≈ 0.59 ε₁ | YES (~0.6·ε₁) |
| 1/ε₁ | ≈ 474 | YES (they cite 473) |

**Finding 12.** Physical prefactor is slightly optimistic: kinetic energy density for a slowly moving object is (γ−1)E₀ ≈ ½v² E₀, so **per unit charge** one should compare ½v² (E/Q) to ε₁ if ε₁ is a **fractional** binding (energy per cloud charge ~ ε₁ in units where the free threshold is m). With E/Q ~ O(1.4) on-branch, v ≪ √(2ε₁/(E/Q)) ≈ **0.055**, not 0.065. Order-of-magnitude unchanged; the upper edge v=0.05 is then **~ fully marginal / slightly over**, not "0.6·ε₁".

**Finding 13.** The deeper issue is **protocol**, not arithmetic:

- X10c end-state is **not stationary** — monopole breather period ~500–600 t.u. (X10_RESULTS). Boosting a breathing cloud makes "strip vs co-move" confounded with breath phase and with radiation-reaction cascade leftovers.
- Global co-boost of (ball+cloud) vs phase-tilt on the ball only are different experiments; GROUNDING's "seed-level phase tilt is impulsive" assumes the dangerous relative-kick case. PROPOSAL EX-1 does not specify which.
- Binding energy of the *retained* cloud is ε₁×|Q_cloud| with |Q_cloud|~0.43 at T=2000 — total bound energy is tiny; diagnostics must normalize carefully or "survives" will be noise-dominated.

Ramped multi-step boost at the upper end: **yes, required** if the question is adiabatic transport; impulsive v=0.05 is a stripping probe, not an adiabaticity probe.

### A6. KAM framing (GROUNDING §3)

**Finding 14.** KAM paragraph is appropriately non-quantitative. No attack. PROPOSAL I.1 still lists Kuksin–Wayne as a "theoretical body that applies" in parallel with VK and FGR — softer overclaim than citing it for a rate, but it should be demoted to "metaphor only" in the proposal to match GROUNDING §3.

### A7. Claim C1 strength

**Attacked sentence (GROUNDING §4):**
> "Claim C1 (fabric theorem-shaped): *in this system, a configuration is orbitally stable iff it is a GSS-signature minimizer within its conserved sector … and every other localized excitation carries a golden-rule width unless blocked by multipole selection or arithmetic (with lattice band-top blocks flagged as artifacts).*"

**Finding 15.** "iff" and "theorem-shaped" overclaim what HC-1..5 can deliver:

- GSS gives **sufficient** conditions under spectral hypotheses; necessity fails if other protections exist (topological, additional hidden symmetries, monopole selection already in the document).
- "every other localized excitation carries a golden-rule width unless blocked" is a **generic-FGR ideology**, not a theorem for this PDE; GROUNDING §2 already said the matrix element is not guaranteed.
- Honest C1 strength: *necessary numerical program to test whether GSS signature + FGR/multipole/arithmetic classification organizes the observed long-lived spectrum*. That is still a strong program deliverable; it is not a theorem.

---

## B. PROPOSAL design

### B1. Do HC-1..HC-5 suffice for C1 as stated?

**Short answer: No — not for C1 as written; partially for a weakened C1.**

| C1 sub-claim | Covered by? | Gap |
|---|---|---|
| Orbital stability ↔ GSS signature match | HC-3 only builds **D = ∂Q_a/∂ω_b** (or Hessian of E(Q) in PROPOSAL I.1 — dual language, inconsistent between docs) | **Never computes n(H_ω)**; cannot "match" GSS counts with one side missing |
| Sector = flavor charges × winding | HC-3: partitions (Q₀,Q₁,Q₂) only | **Winding n omitted** from the scan; fiber charge omitted |
| Mode catalog | HC-1 BdG | No gauge fixing / Gauss projection / multipole labels / continuum-edge report in the rung text |
| Golden-rule widths unless blocked | HC-2 arithmetic + HC-4 widths | HC-2 uses **wrong edge formula risk** (Finding 6); HC-4 may measure protected (Γ=0) lines |
| Lattice band-top artifacts flagged | Mentioned in GROUNDING only | **No HC rung** does dx-refinement of a "stable by arithmetic" mode |
| Dynamical orbital stability | HC-5 overload decay of *off-sector* content | Does not evolve **on-sector GSS-stable vs GSS-unstable partitions** as a controlled pair |

**Finding 16.** PROPOSAL I.1 vs GROUNDING language clash:

- PROPOSAL: "Multi-charge stability = positivity of the Hessian of E(Q₀,Q₁,Q₂)"
- GROUNDING: signature of D_ab = ∂Q_a/∂ω_b with n(H)=p(D)

These are Legendre-dual stories only when the branch is smooth and nondegenerate; "positivity of Hess E(Q)" is the **energy-minimizer** criterion (good), not automatically the same checklist as a botched n(H)=p(D). The two documents must lock one operational definition before HC-3 is coded.

**Finding 17.** PROPOSAL slogan:
> "the number of independently stable harmonics = the dimension of the conserved-charge lattice"

This does **not** follow from GSS. Dimension of the charge lattice labels **sectors**; it does not equal the number of stable *internal excitations*. Stable objects = minimizers per sector; each minimizer can still support protected bound modes (monopole breather is the archive's counterexample to "everything else decays"). STATE §4 already flags monopole-protected breathing as a structural annoyance — the slogan fights the archive.

### B2. EX-1 adiabaticity (design side)

Arithmetic: see Findings 11–13. Design hardenings:

1. Pre-register two protocols: (i) **global** velocity field on the full end-state; (ii) **relative** boost of core only. Only (ii) tests stripping against ½v² vs ε₁; (i) tests shape-mode excitation under boost.
2. Cap impulsive v at ≤0.03; use **ramped** boost over several breath periods (Δt ≳ 2π/Ω_breath ~ 500–600) for the adiabatic co-motion claim.
3. Do **not** use raw X10c t=2000 as "the atom" without a stationaryization step (or at least breath-phase tagging and multi-phase ensemble).
4. Primary observables: Q_cloud(t), cloud COM vs core COM lag, binding proxy ε_eff, multipole content of radiated flux (ties to EX-2).

### B3. Which rung is most likely to produce a misleading result?

**Finding 18. Most misleading: HC-2 (resonance arithmetic audit), with HC-4 and EX-1 as close runners-up.**

**Why HC-2:**

- Cheap, analysis-only, will produce a crisp "stable/decaying" table that the rest of the census "scores against."
- Continuum edge is underspecified / wrong if taken as m=1.5 on QRK frequencies (Finding 6).
- No multipole column in the proposed output.
- No band-top / gauge-band / umklapp columns.
- False "stable" labels will freeze into HC-4 design and into C1 rhetoric before any dynamical check.

**How to harden HC-2:**

1. **Inputs only from HC-1:** discrete eigenvalues Ω_j **and** measured essential-spectrum thresholds (matter and gauge) of the *same* linearization, same gauge fix.
2. Output table columns: mode id, Ω, multipole / component structure, lowest open channel (n, matter vs gauge), continuum edge used, predicted Γ order (0 / ε² / ε⁴ / …), **lattice-artifact flag** (channel above matter or gauge band top at this dx).
3. Mandatory **one dx-pair** on any mode classified "arithmetic-stable" before it is allowed to count toward C1.
4. Explicit frame: co-rotating BdG vs lab; write the combination rules once in GROUNDING and cite them from PROPOSAL.

**HC-4 runner-up:** QRK-1 protocol + kick amplitude on density FFT will report "linewidths" for monopole-protected lines that are not radiative Γ. Harden: (a) separate multipole-projected observables; (b) measure mode energy decay in BdG basis, not only `s_max` FFT width; (c) pre-register Γ≈0 for monopole family as the X10c-consistent null.

**EX-1 runner-up:** non-stationary breathing cloud + impulsive v=0.05 → "stripped" or "survived" will track breath phase and sponge more than transport physics. Harden: §B2 list.

**HC-3** is misleading in a quieter way: it will plot a pretty "stability region" of partitions from signature(D) alone and label it GSS without n(H_ω). Harden: either compute Morse index of H_ω at representative points (hard) or **downgrade the claim** to "VK-matrix region (necessary slope condition)" and require dynamical HC-5-style tests on both sides of the D boundary.

**HC-5** is the least misleading if it stays a pure off-sector decay test; it does not by itself prove C1.

### B4. Order of work

PROPOSAL order HC-1 → HC-2/3 → EX-2 → HC-4/5 → EX-1/4 is sensible **if** HC-1 is truly complete (gauged, multipoles, continuum edges). EX-2 before EX-1 is good (mine existing flux). EX-3 (boost-frame invariance vs dx) remains the overdue D7-class chore STATE already flags; parking it after EX-1 is acceptable but it should not be optional if EX-1 claims "movement is simulated correctly."

---

## C. Missing analyses / theory checks (up to 3)

**Finding 19. Missing item 1 — Continuum edge of the *linearized* gauged operator (matter + gauge), with gauge fix.**  
Neither GROUNDING nor PROPOSAL specifies how ω_edge is computed for this system. Without it, HC-2 is free arithmetic on the wrong number. This is a theory check + a one-time numerical spectral computation on the HC-1 background.

**Finding 20. Missing item 2 — Both sides of the GSS match: n(H_ω) (or dynamical orbital-stability pairs), not only D or Hess E(Q).**  
HC-3 as written is a slope/Hessian scan. C1 requires the index match or a dynamical substitute: evolve two nearby partitions, one with "good" D signature and one with "bad," at fixed numerical budget, and show orbital stability vs drift. Currently absent.

**Finding 21. Missing item 3 — dx-refinement / lattice-artifact battery bound to the census (and to EX-1 kinematics).**  
STATE §4 already lists unaudited resolution and undone X2/D7 as load-bearing. GROUNDING's band-top warning demands a refinement test that no rung schedules. Same battery should cover: (i) arithmetic-stable modes, (ii) QRK linewidths, (iii) process-ledger 0.25% residue (STATE/exchange prior), (iv) EX-1 co-motion at two dx. Without this, C1's "lattice band-top blocks flagged as artifacts" is a slogan with no measurement.

**Honorable mentions (not in the top 3 but real):** multipole-projected observables in the toolpath; stationaryization of X10c before boost; winding-n in the stability matrix; gauge2/fiber out of scope statement (CUDA lacks gauge2 — STATE §2 — so C1 sector language should explicitly exclude fiber until ported).

---

## D. Verdicts

### D1. `v86/GROUNDING.md` — **SOUND-WITH-FIXES**

Not unsound as a research grounding sketch: dual d(ω), FGR genericity caveat, monopole×gauge remark, KAM humility, and EX-1 order-of-magnitude are in the right intellectual neighborhood. Not sound as written because of a **hard counting error** in the GSS criterion and an **underspecified continuum edge** that poisons the resonance section.

**Required fixes (enumerate):**

1. **GSS sign/count:** With D = ∂Q/∂ω, state stability match as **n(H_ω) = n(D)** (negative eigenvalues), **or** redefine D = −∂Q/∂ω = d″ and keep n(H_ω)=p(D). Show the single-charge reduction algebraically in one line.
2. **Rename H_ω:** "Hessian of the augmented action S_ω = E − Σ ω_a Q_a," not "constrained linearized energy," unless you truly mean the reduced Hessian (then stability is n=0, and the D relation is an index formula — write it as such).
3. **State that D alone is not GSS:** HC-3 computes inertia(D); full GSS needs n(H_ω) or a dynamical proxy; list kernel (3×U(1), translations, gauge orbit).
4. **Gauss subsection:** one paragraph on physical subspace (discrete Coulomb gauge + Gauss projection) and explicit note that **massless gauge essential spectrum reaches 0**, so GSS "no unstable essential spectrum" is a nontrivial assumption here, not a lattice-band footnote.
5. **Replace "ω_edge (here the gap m)"** with: edge from HC-1 linearization (matter + gauge), with co-rotating/BdG frame specified; vacuum m=1.5 is not the HC-2 threshold for QRK Ω's.
6. **Band-top warning:** keep ω_max² = m² + 4d/dx² for free matter; add gauge band [0, 2√d/dx]; add umklapp; require dx-refinement protocol.
7. **C1:** demote "iff / theorem-shaped" to testable organizing claim; drop necessity; keep multipole/arithmetic/lattice caveats.
8. **EX-1:** keep √(2ε₁)≈0.065 as the baseline; note E/Q prefactor → ~0.055; require ramped boost + protocol split (global vs relative) + non-stationary cloud caveat.

### D2. `v86/PROPOSAL.md` — **SOUND-WITH-FIXES**

Census + exchange split matches STATE open items (mode spectroscopy, transport of bound systems, unaudited movement). Rung order is mostly right. Not sound as a plan that "establishes C1" because the rungs do not implement the GSS match, the continuum edge, multipoles, or lattice refinement that GROUNDING itself says C1 needs.

**Required fixes (enumerate):**

1. **Lock PROPOSAL I.1 to corrected GROUNDING:** one operational stability criterion (Hess E(Q) minimizer **or** corrected GSS index) — not two divergent slogans.
2. **Strike or heavily qualify** "number of independently stable harmonics = dimension of the conserved-charge lattice" (Finding 17); replace with "stable *ground states* are classified by the conserved-charge lattice; extra excitations decay unless multipole/arithmetic protected."
3. **HC-1 deliverables (expand):** gauged BdG or equivalent linearization; gauge fix; multipole / component structure; essential-spectrum edges (matter, gauge); list of Ω_j used by all later rungs.
4. **HC-2 deliverables (expand):** table per Finding 18 hardening; forbid "stable" without edge source citation and artifact flag.
5. **HC-3 deliverables (expand):** report full inertia of D(ω); sample n(H_ω) or dynamical stability pairs across the D-boundary; include winding if C1 claims winding; state fiber out of scope.
6. **HC-4:** define Γ operationally (mode energy decay vs FFT width); pre-register monopole Γ≈0; amplitude scan only on multipole-open lines for FGR ε² test.
7. **HC-5:** keep; add a sibling **HC-3b** dynamical test of partitions on both sides of the D boundary (or fold into HC-5).
8. **EX-1:** protocol split, v-ramp, stationaryization/phase control, v≤0.03 impulsive or ramped to 0.05; success metrics beyond "looks bound."
9. **Add a cheap refinement rung** (can be HC-2b / EX-0) for band-top and exchange-residue dx dependence — required by GROUNDING §2 and STATE §4.
10. **KAM:** match GROUNDING §3 — framing only, not a third "body that applies" at equal rank with VK/FGR.

### D3. Joint readiness for runs

Analysis rungs HC-1/2/3 should **not** proceed to "C1 established" language until fixes 1–5 (GROUNDING) and 1–5 (PROPOSAL) are edited into the documents. CPU HC-4/5 and GPU EX-1 should stay gated on the hardened HC-1 edge catalog. This is a documentation/design repair, not a reason to abandon the census.

---

## Finding index (quick)

| # | Topic | Severity |
|---|---|---|
| 1 | GSS n(H)=p(D) with D=∂Q/∂ω fails single-charge VK reduction | **Critical** |
| 2 | "Constrained" H_ω terminology wrong/ambiguous | High |
| 3 | D signature alone ≠ full GSS | High |
| 4 | Winding sector in C1 but not in 3×3 object | Medium |
| 5 | Gauss projection hand-waved; massless gauge ESS ignored | High |
| 6 | ω_edge = m mis-applied to QRK internal frequencies | **Critical** |
| 7 | Multipole×gauge incomplete vs HC deliverables | Medium |
| 8 | HC-4 linewidth ≠ FGR Γ for protected modes | High |
| 9 | Band-top formula algebra OK | (confirm) |
| 10 | Band-top warning incomplete (gauge, umklapp, background A) | High |
| 11–12 | EX-1 √(2ε₁) OK; E/Q prefactor tightens v_max | Medium |
| 13 | EX-1 protocol / non-stationary cloud | High |
| 14 | KAM OK in GROUNDING; over-ranked in PROPOSAL | Low |
| 15 | C1 "iff/theorem-shaped" overclaim | High |
| 16 | PROPOSAL Hess E(Q) vs GROUNDING D clash | High |
| 17 | "stable harmonics = charge-lattice dim" false slogan | High |
| 18 | HC-2 most misleading rung | High |
| 19–21 | Missing: linearized edges; n(H_ω)/dynamical GSS; dx battery | High |

---

*End of seat report — grok45.*
