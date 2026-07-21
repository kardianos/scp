# v85 ANALYSIS — deep review + analysis-only execution

**Date:** 2026-07-21
**Agent:** Claude Fable 5, dedicated deep-analysis session (read-only on kernel; no lattice runs)
**Scope:** Part 1 = adversarial + constructive critique of `v85/PROPOSAL.md` + `v85/FULLNESS.md` against the v66–v84 archive. Part 2 = execution of every analysis-only rung (X6a-branch, X6a-shedding, D12 design, X10 semi-analytic, three-U(1) code verification, T1 beat derivation).
**Tools used:** existing branch data (`v69/theory/gscan.tsv`), local diag time series (`v74/results/*_diag.tsv`), a scratchpad re-drive of the **existing** v69 gauged-shooter numerics (verbatim copy of `v69/theory/gauged_shooter_fast.py`), scratch python (KG solver), and read-only inspection of `sfa/sim/scp_sim.c`. **No kernel edits, no scp_sim executions, no remote/GPU work.** Scratch scripts live in the session scratchpad (`branch_drive2.py`, `a_analysis.py`, `shed_analysis.py`, `shed_refine.py`, `kg_coulomb2.py`, `c_scan.py`); every number quoted below traces to a command + output shown in §2.

---

## 0. Executive verdicts

| # | Claim (verbatim target) | Verdict | Load-bearing evidence |
|---|---|---|---|
| 1 | H3 "E = Qω" as the fundamental energy law | **REFUTED as exact; survives only as thin-wall asymptote.** ε ≡ E/(Qω)−1 is measured **+0.9% to +4.3%** everywhere on every computed branch (§2A step 3); the proposal's own §1 equation block carries the (1+ε). | branch tables, `a_analysis.py` |
| 1 | H4 "M = Qω/c²" closes council D5 | **Units: yes. Calibration: no.** A new gap opens — Qω and ∫T₀₀ disagree by the same measured 1–4%, and the proposal does not say which is inertial. The "E=mc² demotion" is definitionally circular (E₀≡Qω, M≡Qω/c² ⟹ E₀=Mc² is a tautology). | §1.1 |
| 2 | H5 "spectrum of persistent states is discrete" | **Contradicted by the archive as stated.** The theory's own persistent states form *continuous* families: Q ∈ (86.6, 924) gauged branch, continuous polarization spin (v73 §6.1), non-quantized ring flux (v73 magneton: 0.374 vs quantum 126). Survivorship discretizes *mode indices at fixed charges*, not charges/action. | §1.2 |
| 2 | H6 Bohr–Sommerfeld = closure | **Partially supported:** action-per-cycle = 2πQ is measured to 0.03% (v73 P1). The leap to "flux quantization / quantized exchange" is contradicted by v73's own flux measurement. X5 remains the fair test. | §1.2 |
| 3 | §2.1 kernel inventory + interference cross-term | **Mostly correct; two omissions.** Cross-term algebra verified against Vt(s) (§2F); but beat-averaging kills only the *first-order* term — a second-order rectified attraction (coefficient 9Pf₁²f₂²) survives at any Δω, and the gauge Coulomb force never averages. T2-absent verdict confirmed in code. | §2F, §2E |
| 4 | §4 CLOSED/TRANSIENT/ARTIFACT classifier | **Right idea, under-specified.** As written it mis-files every absorbing-BC run: boundary-tail loss refines with **L**, not dt/dx. Certification needs 4 axes (dt, dx, L, seed precision). Good news found in code: compute is f64 everywhere; "f32" is snapshot/seed I/O only. | §1.4, §2E |
| 5 | FULLNESS iron-peak Q*; c6 park = fullness sighting | **Interior Q\* does not exist — provably.** An interior min of E/Q requires ε=0; measured min ε>0 on all branches. E/Q min sits AT the capacity fold (Q_max = 921–924). Q=650 is curvature-featureless and the run was still shedding (dQ/dt=−0.17): the park is **kinetic**, not a branch feature. Softener: by the frozen F(Q), F(650)=0.98 — marginal binding is 98% saturated there. **Q_max as capacity: PASS. c6=Q\*: FAIL.** | §2A |
| 5 | F(Q) candidate (D11 curvature formula) | **Rejected — ill-conditioned** (normalization diverges at the bottom fold; F_curv≈1 everywhere by Q=114). **Frozen replacement delivered:** F(Q) = [ω(Q_bot)−ω(Q)]/[ω(Q_bot)−ω_fold]. | §2A step 4 |
| 6 | Three global U(1)s at η=0; only total Q gauged | **VERIFIED in code** with file:line cites (§2E). But "mismatched pairs **cannot fully annihilate**" **overclaims**: free massive Φ_a waves are legal flavor-charge carriers; the block is a *branching/rate* statement, not an absolute conservation stall. v71's "no isolated components" forces the residual to disperse or form a two-flavor Coulomb cloud — a sharper, still-testable X8 prediction. | §1.6, §2E |
| 7 | SHELL: KG-Coulomb response modes; D12 | **Spectrum is real — computed here to machine precision (X10 predictions delivered).** D12 confirmed quantitatively: at carbon Q_N=650, a₀/R = **1.04** (30% of the 1s sits inside r_half). Two viable designs delivered (§2C). New structural finding: **α_max = g²Q_max/4π ≈ 0.17–0.21, nearly g-independent** — a built-in fine-structure ceiling; supercritical shells are impossible. | §2C, §2D |
| 7 | SHELL protection interlock ("capacity" + detuning) | **Capacity protection is sign-confused.** F→1 saturation governs *same-sign accretion*; annihilation of an opposite-charge cloud quantum releases ~2m regardless of fullness. The cloud's survival rests **only** on detuning/kinetics (+flavor if engineered). v70 pos1 (slow annihilation, not never) is the measured warning. Expect F10 to fire on long T unless flavor protection is added. | §1.7 |
| 8 | D6 Bell out-of-scope declaration | **Honestly drawn.** Minor rhetorical overreach in §3.1 (delayed-choice/eraser asserted as resolved). | §1.8 |
| — | FULLNESS §3 cites "v83/e1 sign runs" as the Δω=0 phase prior | **Mis-citation.** Those runs are lock-charge tests (Q_phi=0, locks_medium_only=1, oscillon init — verified from cfg+log). The real prior is **v70 EXISTENCE** (cos(Δφ)·e^{−0.56D} force; detuned-pair beat with 4–10× Coulomb surges; lock range <0.06 at D=11). | §1.3, §2F |
| — | D8 mass defect as "Δω bookkeeping" | **Product form Qω fails by 33% relative** on the c6 defect (3.6% vs 5.3%). **Fix available and exact:** integral form ΔE = ∫ω dQ (Legendre dE/dQ=ω verified to 8×10⁻⁵ here; 0.1–0.26% in-kernel, v72). | §1.1, §2A |
| — | F6 (c12 shedding continuous?) | **AMBIGUOUS at diag resolution.** The shed is *episodic* (two smooth bursts of 185 and 330 units, 14× rate contrast), not stepped (max 2.5-t.u. drop 19.4 ≪ 88) and not clean drool. Episode sizes = 2.06 and 3.68 minimal closures — not integers. Cluster-resolved verdict **BLOCKED** on the archived 42.4 GB SFA (§5). | §2B |

---

## 1. Critique of the v85 stack

### 1.1 E = Qω fundamental, E = mc² derived; M = Qω/c² (PROPOSAL §1, H3, H4)

**Steelman.** Three measured legs support a Qω-centric ontology: (i) ℏ_eff ≡ E₀/ω = p/k = Q(1+ε) measured to 3–5% (v70); (ii) the throughput identity 2E_kin = ωQ to 0.03% (v73 P1) — the *kinetic* part is exactly ωQ; (iii) the Legendre relation dE/dQ = ω, confirmed at 0.10–0.26% in-kernel (v72 §5) and at 8×10⁻⁵ on the shooter branch here (§2A). Mass-as-closure-frequency is a coherent soliton-physics stance, and Qω/c² has genuinely correct mass units, which ENT/LQ did not. D5's *units* objection is closed.

**Attack.**
1. **H3 is false as an exact statement, by the document's own equation.** PROPOSAL line 42 writes ℏ_eff = E₀/ω = Q(1+ε) — the ε is *in the axiom's source equation* — then H3 drops it. The branch data (§2A) make ε mandatory: ε ranges **+1.6% → +4.0%** on the g=0.05 stable branch (min at the cap, max at Q≈124), +0.9% → +4.3% at g=0.02, +3.4% → +4.3% ungauged (computed range; v70's 3–5% agrees). E = Qω is the **thin-wall asymptote**, approached but never reached. The proposal must either state H3 as an asymptotic/idealized law or redefine E ≡ Qω — and the second option makes "E=mc² is derived" **circular**: with E₀ ≡ Qω and M ≡ Qω/c², E₀ = Mc² is an identity of definitions, not a derivation. The actual physical content — boosted balls satisfy E² = p² + E₀² with E₀ = ∫T₀₀ (v70) — uses the *estimator* the creed demoted, and Qω/c² fails that dispersion relation at the ε level. **A calibration lock is needed: which object (Qω or ∫T₀₀) is the inertial mass? They differ by a measured 1–4%.** (Proposed name: D5′.)
2. **The "Kepler demotion" narrative is decorative where it is not circular.** "Mass is propagation at c closed into a cycle": the ball's internal clock runs at ω = 1.41–1.48 < m = 1.5, *below* the free band — internal circulation does not "saturate the locality bound"; nothing in the kernel pins any internal speed to c. What is true and measured is lattice SR kinematics (E²=p²+E₀², γω clocks, k=γωv). That content predates v85 and is not strengthened by the demotion language.
3. **H4's composite accounting fails quantitatively in product form.** For c6: defect via E/Q = 5.3% (branch) / 4.85% (t=300, v74); via 1−ω_f/ω_i = **3.6%** — a 33% relative understatement (computation shown in §2A step 5). The error term is exactly the ε(Q) drift (ε: 0.0402→0.0211 across the fusion). **Constructive fix:** because dE/dQ = ω holds to 0.26%, the *integral* form ΔE = ∫_{Q_i}^{Q_f} ω(Q)dQ (+radiation) is exact bookkeeping in ω alone — freeze D8 on the integral form, never the product form.

**Verdict.** D5 closed as units; H3/H4 must be re-stated as "E = Qω(1+ε(Q)), ε measured, →0 thin-wall" or the axioms contradict the archive they cite. The demotion of E=mc² adds no predictive content in its current form.

### 1.2 H5 closure survivorship; H6 old quantum theory (PROPOSAL §2)

**Steelman.** "Quantization is survivorship, not granularity" is a respectable soliton stance: the *existence* windows (ω ∈ (1.3087,1.5); (1.40624,1.4975) gauged), the VK-stable segment, branch *count*, and ring winding n are all genuinely discrete structures on a continuum substrate. v72's fixed-Q flow finding discrete attractor destinations at fixed Q is real evidence for attractor discreteness. H6's core identity — one action quantum h_eff = 2πQ per re-constitution cycle — is *measured* (v73 P1, 0.03%).

**Attack.**
1. **The persistent-state spectrum of this theory is measurably continuous.** The flagship particle family is a **one-parameter continuum**: every Q in (86.6, 923.9) at g=0.05 is an exactly-closed stationary state (§2A). Nothing in H5 discretizes Q, and therefore nothing discretizes E or M.
2. **The archive already contains three anti-H6 measurements:** (i) polarization spin is *classically continuous* — "Unlike the ring's topological J = nQ… polarization spin is classically continuous" (v73 §6.1); (ii) ring flux does **not** quantize — "trapped flux is far below the London/vorton quantum 2π/g ≈ 126 — not a flux-quantizing superconducting vortex" (v73 §5.4); (iii) transfer magnitude is continuous (T2 absent — §2.1's own inventory). H6's list ("flux quantization… quantized exchange are all instances") is contradicted by (ii) and (iii) on this kernel.
3. **Well-posed retreat exists:** H5-weak — *at fixed values of all conserved charges, persistent states are isolated (countable) attractors; transients flow onto them*. That is testable (X2 flow map) and consistent with the archive. H5-strong (discrete spectrum → quantum phenomenology) requires T2, which the kernel lacks; the proposal even concedes this in §2.1 but does not propagate the concession back into H5/H6's wording.
4. **Guitar-string analogy status:** load-bearing for *mode* discreteness at fixed boundary conditions (that is exactly X5/X10's content — discrete radii/response modes at fixed Q); decorative for charge/action/energy quantization (a string's amplitude — the analog of Q here — is continuous, which is precisely the problem).
5. **H6's universality gap:** Planck's constant is universal; ℏ_eff = Q is per-object (Q=90 ball vs Q=650 nucleus). de Broglie λ = 2πQ/p is object-charge-dependent in a way real QM's λ = h/p is not. CONCEPT §7 flags this honestly ("not a derivation of quantum mechanics"); H6's prose papers over it.

**Verdict.** Keep H5-weak; strike or re-scope H5-strong and H6's second sentence. X5 remains the right decisive test — and note X5's prediction should be stated with the orbiter's own Q fixed, since ℏ_eff varies per object.

### 1.3 H7 T1/T2/T3 and the §2.1 kernel inventory (incl. cross-term correction)

**Steelman.** The inventory table is the most self-honest piece of the proposal: it correctly identifies linear superposition/evanescence as emergent T3, the mass-gap resonance kinematics as partial T1, phase-blind V(s) as the T1 violation at equal-clock dense overlap, and **T2 as absent and not expressible by config**. The 2026-07-21 correction (object-level gating via interference cross-terms) is directionally right and the per-component identity |Φ_a|² = f₁²+f₂²+2f₁f₂cos(Δωt+δ) is algebraically correct (verified against the kernel's exact potential in §2F).

**Attack.**
1. **The beat-averaging conclusion is incomplete.** FULLNESS §3: "the cross term beats at Δω and time-averages toward zero — detuned objects interact weakly at contact." Doing the 3-component product honestly (§2F): the time-average of s keeps a **rectified second-order term 9P f₁²f₂²** (P=f₁²+f₂²) with fixed attractive sign (μ<0), range e^{−2κD}, plus the **gauge Coulomb term, which never averages**. Detuned objects do *not* decouple; they lose only the O(f₁f₂) coherent force. At saturated contact (κs ≈ 2.5 at core), Vt′ is knee-suppressed by (1+κs)² ≈ 12 — the phase force flattens there for a *fullness* reason, not a detuning reason.
2. **Mis-citation:** "the v83/e1 sign_opp/sign_same runs probed exactly this sign structure" — they did not. Those configs are `locks_medium_only=1`, `init=oscillon`, `n_locks=2` with charges ±4 and **Q_phi = 0.000000e+00** in the completion log: a lock-charge Coulomb sign test with no Q-balls and no matter phase. The correct measured prior is **v70 EXISTENCE**: contact force = cos(Δφ)·e^{−0.56D} riding on Coulomb (Δφ=0 fuses through Coulomb at 2.7×; π/2 → pure Coulomb, ratio 0.94; π → enhanced escape), and the genuinely detuned pair (b2_detuned, ω=1.46/1.42, D=11) shows the force *beating* at Δω with a 4–10× Coulomb surge at anti-phase passage and **no locking (lock range < 0.06 in Δω at D≥11)**. Note the v70 quadrature ratio 0.94 is consistent with my predicted second-order residual (~3% of the coherent term at D≈8–11, §2F step 3).
3. **T3's conflation:** sub-gap evanescent tunneling (linear physics) and soliton pass-through (v73 E2, 99.75% layment — a *nonlinear conduction* phenomenon) are called "the same mechanism." They are different mechanisms with different failure modes; the rhetorical unification buys nothing and risks mis-designing X1.
4. **X1's "alignment" is undefined as written.** All propagating radiation lies *above* the gap (ω_rad > m = 1.5 for Φ, > 1.6 for θ), while every ball clock is *below* it (ω_ball ≤ 1.4975): Δω = ω_rad − ω_ball > 0 strictly — the X1-A sweep can never cross alignment. Gating, if present, will appear as resonance with the ball's *internal linear response modes* (breathing/quadrupole), not clock matching. F1a's "flat curve kills" is uninterpretable until "aligned" is operationally defined (recommend: define against the measured linear-response spectrum of the parked ball — which X10's machinery computes).

**Verdict.** Inventory: keep, with the §2F correction folded into the T1 row. X7's quantitative prediction should be the three-regime law of §2F step 3 (coherent / impulse-suppressed 2/(Δωτ) / locked below Δω_lock), with v70 b2_detuned as the anchor.

### 1.4 §4 drift-tolerance abolition and the CLOSED/TRANSIENT/ARTIFACT classifier

**Steelman.** Replacing endpoint-drift tolerances with drift-*rate asymptotics* is a real methodological upgrade, and the archive already contains its proof-of-concept: v72's exact discrete stationary state out-drifts the interpolated continuum profile (−0.038% vs −0.072%), showing the team can distinguish converging-to-closure from slowly-dying when it tries. The three-way classification with a mandatory refinement-scaling certificate (D7) is the right shape.

**Attack.**
1. **The classifier as written mis-files every real run.** On absorbing BC, the ball's evanescent tail overlaps the sponge and is eaten *forever*: v72's T=1000 run shows "slow LINEAR boundary-tail absorption −2.5e-4%/t.u. — no acceleration." That loss rate scales with **box size L** (∝ e^{−2κ(L−DW)}-class) and does **not** refine away under dt/dx refinement. Under §4.2's table this lands in ARTIFACT ("tracks BC → fix or discard") — but it cannot be "fixed" at finite L, so *every* state, including true closures, fails CLOSED certification. The classifier needs a fourth column: **BC-limited**, certified by exponential convergence of the drift rate in L at fixed (dt,dx).
2. **"There is no third category" is metaphysics smuggled into methodology.** A transient with τ ≫ any affordable T is empirically indistinguishable from CLOSED. The honest output of X2 is never "CLOSED," only "not shown TRANSIENT at refinement R" — D7 should say so explicitly.
3. **What certification would actually look like on this codebase** (concrete protocol, all config-level, CPU-feasible at N=64–128):
   - Axes: (a) dt: dt_factor × {1, ½, ¼} — expect drift-rate component ∝ dt² (velocity-Verlet, `verlet_step`, scp_sim.c:1583); (b) dx: N × {1, 1.5, 2} at fixed L — expect ∝ dx²; (c) L: L → L+4, L+8 at fixed dx — expect exponential decay with slope 2κ = 2√(m²−ω²); (d) **seed quantization**: seeds written f32 vs f64. Code fact established in §2E: the *compute state is f64 on both CPU and GPU* — `precision` (scp_config.h:123) selects only SFA output width, and the GPU f32 buffer is staging for snapshots (scp_sim.cu:244, 4092). So the "f32 floor" language in v72 is seed/snapshot I/O, not arithmetic — dt/dx refinement will not hit a float wall, but f32-written seeds inject ~10⁻⁷ relative off-shell noise that radiates; certification requires f64 seeds.
   - Fit: rate(dt,dx,L) = A·dt² + B·dx² + C·e^{−2κ(L−DW)} + D_res; **CLOSED(R)** iff D_res is statistically zero.
   - Cheap first target: rerun v72's r5_long (η=0.5, N=64, L=15) at L=19, L=23 — if −2.5e-4%/t.u. drops ×e^{−2κ·4}, the η-soliton is certified BC-limited-CLOSED and the §4.4 re-audit of v72 closes at near-zero cost.
4. **v73's 0.25% pass-through residue** already has a predicted classification: the wake is lattice dispersion (v70 measured the group-velocity anomaly) → expect ∝ dx²: **ARTIFACT(dx)**, testable with one refinement pair.

**Verdict.** Adopt the classifier with the L-axis and BC-limited category added, D7 reworded to "not shown TRANSIENT at refinement R," and f64 seed I/O mandated for X2.

### 1.5 FULLNESS: μ=dE/dQ saturation, iron-peak Q*, force map, c6 as sighting

Data verdicts in §2A. Summary of the critique:

1. **The interior iron peak does not exist — structurally.** d(E/Q)/dQ = (ω − E/Q)/Q = −ε·ω/Q < 0 wherever ε > 0. Measured min ε: 0.0157 (g=0.05), 0.0090 (g=0.02), 0.0337 (g=0 up to Q=352; gscan continues the decrease to Q=1.4×10⁵). So E/Q is strictly decreasing on every stable branch and its minimum sits **at the capacity fold** (or at infinity, ungauged). FULLNESS §1.1's "minimizer of E/Q" defines Q* = Q_max, making "fusion favorable climbing toward Q*, neutral past it" vacuous — there is no "past it."
2. **c6's park at 650 is not a branch feature.** ω, E/Q, μ, d²E/dQ² are all monotone and featureless through Q=650 (table in §2A step 5). v74 itself records dQ/dt ≈ −0.17/t.u. at T=300 ("parked Q is an upper bound") — FULLNESS §1.2's "a droplet that *stopped absorbing*" misquotes its own source; the droplet was still radiating when the run ended. The park is set by radiative-cooling kinetics. What *is* true: by the frozen F(Q) below, F(650) = 0.980 — the marginal fusion gain from 650→cap is 0.11% per unit charge — so the *qualitative* fullness reading ("nothing much left to gain") is compatible, but it is a property of the whole upper half-branch (F(500)=0.93), not of 650.
3. **Force map rows:** "Strong = closure exchange while μ falls" — supported (fusion branch measured). "Hard core = fullness" — the κ-knee mechanism is real (Vt′ suppression ×12 at core s), and my §2F adds: the knee also kills the *phase-coherent* force at contact, leaving density+gauge — a cleaner statement of "failed-exchange repulsion." "Weak = dissolve/reform" — still a metaphor; no analysis-level content yet (X9 exploratory, agreed). "Neutrino = minimal closure": the minimal closure is Q_min = 86.6–89.7 (three branches, §2A) with E/Q = 1.532 > m — VK-stable but *evaporation-metastable*; calling it a transparency-protected neutrino analog needs the X7 gating result first.
4. **X6b's "near-critical pairs bounce" prediction contradicts the measured precedent.** c12 (co-phase, 2.13×Q_max seed) merged into ONE droplet and evaporated — no bounce, no fission (v74). Fullness-as-saturation would have to overturn that precedent to survive F7; the honest prior is merge-then-evaporate for co-phase contact at any inventory, with "bounce" plausible only for *detuned* or anti-phase pairs — which is X7's variable, not X6b's.

### 1.6 Flavor charges as lepton number; annihilation block (FULLNESS §4)

**Steelman.** The symmetry claim is exactly right, and §2E verifies it in the shipped code at η=0: the potential force on Φ_a is a real, phase-blind scalar times Φ_a (`- 2.0*Vp*fu[a]*prod_rest[a]`, scp_sim.c:1196–1201 and 966–977); the gauge link phase is common to all components (1143–1145); the gauge current is the component-**sum** (1235–1238); the sponge damps (u̇,v̇) isotropically per component (1454–1456), so even boundary losses are flavor-proportional. Three independent U(1)s over one gauged diagonal U(1) is a true structural fact of this kernel, and "charge/lepton-number split" is a fair name for it.

**Attack.**
1. **"The residual flavor charge has no carrier state" is false.** Free Φ_a waves above the gap carry Q_a. A mismatched conjugate pair (partitions q⃗ and −q⃗′, Σq=Σq′) has residual q⃗−q⃗′ with zero total charge; annihilation releasing ≈ 2E_ball ≈ 3Q is energetically ample to radiate the residual as flavored waves (cost ≈ m·Σ_a|Δq_a| ≪ 2E_ball for any partial mismatch). Conservation therefore does **not** stall the channel; it constrains the *final state*.
2. **What conservation actually forbids** is annihilation into pure θ/gauge radiation. Combined with v71's measured "no isolated components" (single-component lumps don't self-bind: s = Π_a|Φ_a|² vanishes identically off the diagonal — code-verified, the potential force is ∝ prod_rest[a]), the sharp X8 prediction becomes: **annihilation completes on the matched part of the partitions; the flavored residual either (a) disperses as free flavored radiation to the sponge, or (b) — if the residual has opposite charges in different components — collapses into a two-flavor Coulomb-bound cloud** (indestructible at η=0 because Q_a are separately conserved and the components cannot contact-annihilate each other). Outcome (b) would be a genuinely new object — and note it is *itself* a SHELL-thread object (a KG-Coulomb self-bound charge pair of the same fabric).
3. **Instrumentation gap:** the kernel diagnoses only summed Q (compute_charges, scp_sim.c:2089–2100) — no per-component Q_a in diag.tsv. X8 cannot be scored without per-flavor charge tracking; that exists only in SFA-frame analysis tools (`sfa_qcomp`). Design X8 with SFA cadence sufficient for per-component ledgers, or it will be unfalsifiable in-run.

**Verdict.** Symmetry claim verified; "blocked by conservation, not by fiat" must be downgraded to "redirected by conservation"; F9's kill condition should be re-worded from "annihilation completes despite mismatch" (it can, legally) to "**the flavored residual fails to persist** (Q_a ledgers close to zero without a surviving carrier)."

### 1.7 SHELL: response harmonics, α_f = g²Q_N/4π, screening, D12/D13 (FULLNESS §5)

**Steelman.** This is the strongest new physics in v85. The ontology is legitimate: around vacuum, a single-component fluctuation is *exactly* a free KG field of mass m (the product potential contributes nothing off the diagonal — code-verified §2E), gauged with the same g; a parked ball sources A₀ = gQ_N/4πr (shooter-verified Coulomb tail, `chi ~ g^2 Q/(4 pi r)`); therefore discrete gap-protected bound modes with ε_n ≈ mα_f²/2n² exist below the gap. §2D computes the full spectrum to machine precision — X10's predictions now exist. Self-limiting screening at neutrality is Thomas–Fermi-plausible. D12/D13 are honestly registered.

**Attack.**
1. **D12 is worse than "redesign before X10/X11" — at carbon scale the current parameters put the cloud *inside* the nucleus:** a₀/r_half(650) = 5.16/4.94 = **1.04**, with 30% of the 1s charge within r_half and ⟨r⟩_1s = 7.5 ≈ the ball's own rms radius (v74 r_core 6.9). At the cap (Q_N=921): ratio 0.65, 59% inside. Confirmed and sharpened; viable redesigns in §2C.
2. **The interlock's "capacity" protection is sign-confused.** FULLNESS §5.2-1: "the nucleus is full — absorbing the cloud's quanta gains no binding (κ-knee / F→1)." Fullness/μ-saturation governs **same-sign** charge accretion. A cloud quantum is **opposite**-charge; absorbing it is *annihilation of one unit*: energy release = μ(Q_N) + ω_n ≈ 1.41 + 1.50 ≈ 2.9 per unit — maximally exothermic, *independent of F(Q)*. The archive already measured this channel open: v70 pos1 opposite-charge pairs "annihilate slowly (charge segregation −4× over 600 t.u.), not never" (CONCEPT §6). So protection #1 is void; the cloud's persistence rests entirely on protection #2 (Δω detuning + T2-absence) and, if engineered, flavor mismatch (§1.6). **Expected X11 failure mode: slow drain, pos1-style — F10 fires on long enough T.** Recommendation: design X11 *with* the §4 flavor block from the start (cloud seeded in a partition orthogonal to the core's conjugate), turning the SHELL survival question into the X8 question — which conservation actually protects.
3. **Bookkeeping error in §5.3:** "a Z=6 nucleus accumulates **−6** of response charge." Neutrality means cloud charge = −Q_N = **−650** fabric units (the gauge field sees gQ; screening completes at ∮E=0). "−6" conflates the structural label Z=6 with the conserved charge; there is no −6 anywhere in the dynamics. (Charitable reading: −6 "nucleon-equivalents" of 108 units — but with no electron quantum, per-mode occupancy is a continuous amplitude, and the sentence as written will mislead X11's success criterion. Freeze the criterion as: uptake stalls at |Q_cloud + Q_N| < tolerance.)
4. **Detuning protection needs a number, and it is not obviously favorable.** Cloud modes sit at ω_n ≈ m−ε ≈ 1.4994–1.4997; the core clock at ω_N = 1.406–1.43: Δω ≈ 0.07–0.09, beat period 70–90 t.u. — *comparable to* the v70-measured surge timescales, and the pair-force there surged 4–10× Coulomb at anti-phase passages rather than averaging quietly. A permanently co-located cloud (τ = ∞) does average — but it is also parametrically driven at the beat frequency for ~10⁴ t.u. runs. X11 should include a beat-frequency spectral monitor on the cloud charge; heating at Δω is the signature that detuning protection is failing dynamically.
5. **New structural ceiling (good news):** α_f,max = g²Q_max/4π = 0.169 / 0.183 / 0.207 at g = 0.02 / 0.05 / 0.10 (§2C step 5) — nearly g-independent because Q_max ≈ const/g². Consequences: (i) shells can never reach KG supercriticality (α = ½) from an on-branch nucleus — no fall-to-center pathologies anywhere in config space; (ii) the deepest possible 1s binding is ε ≈ mα_max²/2 ≈ 0.025 (1.7% of m) — the theory's entire "atomic physics" lives within 2% of the gap; (iii) hydrogenic scaling holds to ≤4% relativistic corrections everywhere (verified against exact KG formula, §2D).
6. **D13 stands as declared** (no Pauli, no 2/8/18); add: the *gross* level structure (n,l) is now predicted (§2D tables) — X10 can be scored on it without any statistics claims.

### 1.8 Bell / D6 out-of-scope declaration (PROPOSAL §3.5)

Honestly drawn: single-particle interference claimed, joint spacelike closure explicitly not claimed, no-EPR registered as D6, F-series does not include a Bell rung. Two cavils: (i) §3.1 asserts delayed-choice/eraser outcomes ("all in the past light cone") as if settled — that is a design-level conjecture; mark it as such or move it under D6; (ii) "a local deterministic field adds amplitudes locally (the right ingredient)" — classical field superposition adds *fields*, not probability amplitudes; whether closure statistics inherit amplitude²-addition is exactly X3-§3.3's open rate-law bet (F5 correctly registers it as survivable). The boundary is honest; keep it.

---

## 2. Analyses (Part 2, A–F)

Scratchpad root (all scripts + outputs): `/tmp/claude-1000/-home-d-code-scp/994ad21e-63ba-4f82-88cd-502be2a05cb3/scratchpad/`. `gs_fast.py` is a **verbatim copy** of `v69/theory/gauged_shooter_fast.py` (the existing validated Newton-BVP shooter); `branch_drive2.py` only re-drives it on denser ω grids and writes TSVs to the scratchpad.

### 2A. X6a-branch — E(Q) curvature → F(Q), Q*

**Step 1 — locate existing branch data.**
```
$ awk -F'\t' 'NR>1{...}' v69/theory/gscan.tsv     # block census
g=0.000: 244 rows ω∈[1.320000,1.497500]; g=0.020: 74 rows [1.360879,1.497500];
g=0.050: 41 rows [1.406245,1.497500];   g=0.100: 17 rows [1.465542,1.497500]
```
Columns: g, omega, f0, chi0, a0_0, Q, E_matter, E_field, E_total, E_over_mQ, dQdw_sign, r_half, weff0. Inference: full-window coverage exists but is coarse near the g=0.05 folds and around Q=650 (adjacent samples at Q=675→529 across one 0.0025 ω step).

**Step 2 — densify with the existing shooter.** `branch_drive2.py` (piecewise ω plans; fine steps only near folds; g-continuation 0→0.02→0.05 from the banked v66 ω=1.45 profile):
```
g=0 seed Q=118.05 res=9.6e-10
g=0.02: leg -> 1.309 reached 1.360891 [FOLD]  → branch_g002.tsv (91 rows)
g=0.05: leg -> 1.309 reached 1.406244 [FOLD]  → branch_g005.tsv (166 rows)
g=0:    mid-branch+bottom-fold only           → branch_g000.tsv (72 rows)
```
Cross-validation against the archive (`a_analysis.py` first line): at ω=1.420000, g=0.05 → Q=311.46, E=456.27, r_half=3.804 — **identical to gscan** (311.46, 456.27, 3.804). Fold locations reproduce gscan to the step size (1.406244 vs 1.406245; 1.360891 vs 1.360879). The g=0 thin-wall side (Q up to 1.44×10⁵) was *not* regenerated — gscan's 244 rows cover it; my g=0 table spans Q ∈ [86.6, 352].

**Step 3 — μ(Q), d²E/dQ², ε(Q)** (`a_analysis.py`, centered differences on the VK-stable side):
```
g=0.05: Legendre max|μ−ω|/ω = 8.0e-05 ; ε=E/(Qω)−1 ∈ [0.0157 @Q=924, 0.0404 @Q=124]
        d(E/Q)/dQ: 0 interior sign changes; E/Q min = 1.42836 AT Q=923.9 (endpoint)
        |d2E/dQ2|: 1.9e-02 @Q_bot=89.7 → 3.4e-05 @Q=467 → 7.1e-07 @Q_max
g=0.02: ε ∈ [0.0090, 0.0429]; E/Q min at Q=5199 (endpoint);  g=0: ε ∈ [0.0337, 0.0434]
```
Inferences: (i) **dE/dQ = ω** verified at 8×10⁻⁵ — the branch is the canonical family; (ii) **no interior minimum of E/Q exists on any branch**, and none *can*: d(E/Q)/dQ = (ω−E/Q)/Q = −εω/Q, and ε > 0 everywhere (min 0.9%). An interior iron peak requires ε → 0 on-branch, which this potential never does. (iii) |d²E/dQ²| falls five decades from bottom fold to cap — the branch is *flat* near capacity, diverging in curvature only at the **bottom** fold.

**Step 4 — F(Q) freeze (lock D11).** The FULLNESS §1.1 curvature candidate F_curv = 1 − |d²E/dQ²|_Q/|d²E/dQ²|_{Q_min} is **ill-conditioned**: its normalization sits at the bottom fold where |d²E/dQ²| → ∞ (value depends entirely on how close to the fold the endpoint is sampled), and it saturates uselessly (F_curv(114) = 0.97 already). **Frozen proposal:**

  **F(Q) ≡ [ω(Q_bot) − ω(Q)] / [ω(Q_bot) − ω_fold]**  — "how far μ = dE/dQ has fallen from the branch-bottom marginal cost toward its capacity asymptote." Monotone, dimensionless, first-derivative only, 0 at the minimal closure, exactly 1 at capacity.

```
Q      F_mu     F_curv   ω(Q)      E/Q       ε          (g=0.05, a_analysis.py)
 90    0.0283   0.7897   1.48253   1.53211   0.0334
114    0.3143   0.9700   1.46008   1.51876   0.0402
210    0.6746   0.9909   1.43179   1.48424   0.0366
311    0.8243   0.9956   1.42004   1.46504   0.0317
650    0.9801   0.9993   1.40781   1.43747   0.0211
921    1.0000   1.0000   1.40625   1.42843   0.0158
```
(Q=88 does not exist on the g=0.05 branch — bottom fold is Q_min=89.7. On g=0.02: F(88)=0.035, F(210)=0.497, F(650)=0.794, F(921)=0.855.)

**Step 5 — verdicts.**
- *Q\* vs c6 park:* fine scan around 650 (a_analysis.py): ω, E/Q, μ, d²E/dQ² all monotone; d²E/dQ² = −1.29e-05 at 650 vs −1.45e-05 at 630 and −1.16e-05 at 670 — **featureless**. Combined with v74's own "still slowly shedding, dQ/dt ≈ −0.17" the park is kinetic. **FAIL for "c6 park = Q\* sighting."** (Softener: F(650) = 0.980; the marginal-gain reading is qualitatively compatible but non-causal and non-specific — F > 0.9 for the entire upper half-branch.)
- *Q_max as capacity:* the fold is real, terminating the branch at ω=1.406244 with Q_max = 923.9 by my fold-refined sampling (the canonical 921.0 is gscan's sampling of the same fold; capacity = 921–924). **PASS.**
- *Overall fullness–saturation identification:* **AMBIGUOUS-leaning-FAIL** — saturation-of-μ is real (F_mu well-defined, capacity fold real) but the *interior* iron-peak narrative is structurally impossible, and the one claimed sighting dissolves under its own run log.
- *D8 check* (product vs integral bookkeeping, computed inline):
```
ω_f(650)=1.40793 (interp);  defect_E/Q(branch)=5.33%;  defect_Qω=1−ω_f/ω_i=3.57%
→ product form understates by 33% relative; ε_i=0.0402 → ε_f=0.0211 accounts for the gap.
```

### 2B. X6a-shedding — c12_light quantized-lump test (falsifier F6)

**Step 1 — locate the time series.**
```
$ wc -l v74/results/c12_light_diag.tsv  → 1181 rows (t=0→300, Δt_sample=0.25)
$ ls -la v74/results/*.sfa → c12_light.sfa is a SYMLINK → /space/scp/v74/results/ (absent)
$ rclone ls scpsfa:scpsfa/v74/ → 42,377,905,264 c12_light.sfa (cloud-archived)
$ cat /space/scp/v74/results/slices/c12_light.meta.txt → 7 frames only (render panel)
$ grep snap v74/results/c12_light.log → "snap=2.5" ... "122 frames"
```
Inference: the 0.25-t.u. diag series is local; the cluster-resolvable data (122 SFA frames at 2.5 t.u.) is **cloud-archived only** (42.4 GB).

**Step 2 — extract and characterize Q(t)** (`shed_analysis.py`, `shed_refine.py`):
```
Q_phi: 1959.88 → 1410.53 (−549.35). Rates (10-t.u. window, every 5 t.u., t>60):
  −0.11 (t=60) … −1.27 (100) … −5.9 (115) … −1.1 (145) … −1.2 (165) …
  −7.34 (210) … −1.04 (245) … −0.15 (255) … −0.42 (295)
Spike test vs 6×MAD at windows 1.0/2.5/10 t.u.: NO spikes (all windows).
Largest single-sample (0.25 t.u.) drop: −1.945 at t≈210 (= smooth rate ×0.25, not a step).
Largest 2.5-t.u. drop (100<t<250): −19.42.
Episode integrals: shed(100–175)=198, shed(175–245)=315; calm tail (245–300)=17.5
  → burst/calm rate contrast ≈ 14×.
omega_core through episodes: swings 1.27–2.28 (off-branch/hot) while r_core swells
  8.7→13.7 and recontracts to 7.4; ends ω_core=1.311, still Q=1410=1.53×Q_max.
```
Inference: the shed is **episodic, not monotone drool** — two large smooth bursts (≈185 units over t≈100–170, ≈330 units over t≈170–250; boundaries by cumulative-shed table) separated by a lull and followed by near-quiescence (−0.3/t.u.), with the core swelling during bursts (collective-mode-driven evaporation cycles, GDR-adjacent), **but nowhere stepped**: no discontinuity exceeds 19.4 units per 2.5 t.u., against a minimal-closure quantum of 88–90.

**Step 3 — candidate quantum comparison.** Minimal closure from §2A: Q_min = 89.7 (g=0.05 branch bottom; v72 η-branch bottom 88; ungauged 86.6). Episode sizes: 185/89.7 = **2.06**; 330/89.7 = **3.68**. Against the light nucleon (114.13): 1.62 and 2.89. Neither series is integer; the near-2.06 alone is not evidence (the second episode falls mid-gap).

**Step 4 — verdict on F6: AMBIGUOUS.**
- *Against LUMPY:* no steps at any window; within-episode flow is smooth; episode sizes not integer multiples of any closure candidate.
- *Against CONTINUOUS (as "featureless drool"):* the 14× episodic rate structure rejects a constant-rate evaporation null; the discharge is organized by collective dynamics.
- *Sensitivity limit (why diag can't settle it):* a detached subcritical lump (v~≪1) takes ≳10–30 t.u. to cross from core to sponge; Q_phi (whole-box) smears any step over that transit, and Q_core's adaptive radius (r_core swinging 5.5→13.7) makes it definition-noisy (dQ_core/dt MAD = 2.6 with ±70 excursions). Lumps of ~90 units emitted ≈ every 30 t.u. would be **indistinguishable from the observed ramps** in this observable.
- **The decisive test is cluster-resolved and BLOCKED** (§5): retrieve the archived SFA (122 frames @ 2.5 t.u. is sufficient cadence) and run `bin/sfa_qball_track --tsv` → per-cluster charge inventory; LUMPY ⟺ discrete departing clusters with Q ≈ 90-class each; CONTINUOUS ⟺ charge leaves only as sub-threshold radiation. No guess is substituted here.

### 2C. D12 — shell scale-separation design

**Step 1 — derivation (first principles).** Gauss on the lattice (v69 conventions, div E = +g ρ_Q; shooter tail `chi ~ g² Q/(4π r)`): a parked ball of charge Q_N sources electrostatic potential energy for a unit-charge fluctuation V(r) = −g·(gQ_N/4πr) ≡ −α_f/r, **α_f = g²Q_N/4π**. A small fluctuation of one component about vacuum is exactly gauged-KG with mass m (the product potential is flat off the diagonal — §2E): [(ω + α_f/r)² + ∇² − m²]ψ = 0. Radial form u'' + [(ω+α/r)² − m² − l(l+1)/r²]u = 0 maps to Schrödinger–Coulomb with E_eff = (ω²−m²)/2, Z_eff = ωα, l_eff(l_eff+1) = l(l+1)−α²; hence hydrogenic **a₀ = 1/(α_f m) = 4π/(g²Q_N m)**, **ε_n ≈ mα_f²/2n²**, exact KG spectrum ω_nl = m[1+α²/N²]^{−1/2}, N = n−l−½+√((l+½)²−α²). **Supercriticality** = imaginary index at l=0: α_f > ½. Current sit: α_f ≤ α_max ≈ 0.18 (step 5) — margin ≈ 2.7× below critical, everywhere, always.

**Step 2 — current parameters (g=0.05, m=1.5), R_ball = r_half from the dense branch** (`c_scan.py`):
```
Q_N   ω(Q)      r_half  α_f      a₀      a₀/r_half  ε₁(hydro)   1s-charge inside r_half
 90   1.48253   2.34    0.0179   37.2    15.92      2.40e-04    0.000
210   1.43179   3.30    0.0418   16.0     4.83      1.31e-03    0.009
650   1.40781   4.94    0.1293    5.16    1.04      1.25e-02    0.300
921   1.40625   5.58    0.1832    3.64    0.65      2.52e-02    0.592
```
Inference: D12's arithmetic confirmed and sharpened — at carbon scale the "atom" is inside its own nucleus (30% of 1s charge within r_half; ⟨r⟩_1s = 7.5 vs droplet rms r_core ≈ 6.9 from v74). Light nuclei (Q ≲ 210) are *already* separated at current parameters.

**Step 3 — scan for a₀ ≥ 5·R_ball.** Design law (exact algebra from step 1): fixing ρ ≡ a₀/R_ball, **ε₁ = 1/(2 m ρ² R_ball²)** — independent of (g, Q_N) separately. Consequences (numbers from `c_scan.py`):
```
ρ=5:  R=2.34 (Q≈90)  → ε₁ ≤ 2.44e-03, T_resolve ~ 2π/ε ≈ 2.6e3 t.u.
      R=3.30 (Q≈210) → ε₁ = 1.22e-03, T ≈ 5.1e3
      R=4.94 (Q≈650) → ε₁ = 5.5e-04,  T ≈ 1.1e4
Self-consistent ρ=5 crossing at g=0.05: Q_N* ≈ 205 (α=0.0409, ω=1.4325, a₀=16.3)
Carbon at ρ=5: α_needed = 1/(5·m·r_half(650)) = 0.0271 → g = 0.0229
   (r_half(650) is g-insensitive: 4.90@g=0.02, 4.94@g=0.05)
α_max = g²Q_max/4π: 0.1655 (g=0.02) / 0.1838 (g=0.05) / 0.2072 (g=0.10) — ceiling ≈0.17–0.21
```
Box feasibility: need L ≳ 2.5–3 a₀ + sponge for the 1s (n=2 modes have ⟨r⟩ = 4–6 a₀ — capture them only in the largest boxes or accept 1s+2p only).

**Step 4 — deliverables (both concrete and viable; no kernel change, config only):**
- **X10 pathfinder (recommended first):** g=0.05 (unchanged), **Q_N ≈ 205** (seed ω=1.4325 profile; VK-stable), α_f=0.041, a₀=16.3, ρ=5.0. Grid **N=384, L=55** (dx=0.286; ball core resolved by 11 cells over r_half). Ringdown: mode existence + profile in T ≈ 2×10³; frequency to ±10% of ε₁ needs T ≈ 1.5×10⁴ (≈2×10⁶ steps — days-class on V100; run the short version first).
- **Carbon-scale X11:** **g = 0.023, Q_N = 650** (re-run the c6_light recipe at g=0.023; Q_max(0.023) ≈ 4353 keeps 650 mid-branch), α_f=0.027, a₀=24.5, ρ=5.0, ε₁=5.5e-04. Grid **N=512, L≈90–100** (dx≈0.35–0.39); V100-32GB class. Seed the cloud in a flavor partition orthogonal to the core's conjugate (see §1.7-2) and monitor beat-frequency heating.
- **No-go direction:** raising ρ at fixed physics costs ε₁ ∝ 1/ρ² — at ρ=10 the carbon shell drops to 1.4e-04 (T~5×10⁴, unmeasurable). ρ=5 is the practical ceiling. Raising m (config-level) raises ε₁ ∝ m at fixed ρ·R — but rescales the entire v66–v74 phenomenology (window, branches, Q_max all move); flag as a program-level decision, not a default.

### 2D. X10 semi-analytic — KG-Coulomb bound spectrum (predictions)

**Step 1 — solver.** `kg_coulomb2.py`: adaptive two-sided shooting (LSODA, rtol→1e-13) on u'' = [m²−(ω+α/r)²+l(l+1)/r²]u; outward from series u ~ r^s, s = ½+√((l+½)²−α²); inward from u ~ r^ν e^{−κr}, κ=√(m²−ω²), ν=ωα/κ; log-derivative match at r_m = 0.7a₀n²; brentq on ω; residual check |miss|<1e-4 rejects pole-crossings (two levels initially mis-bracketed by pole roots were re-solved with the tightened bracket — shown in session log; both land on the closed form to 3e-15).

**Step 2 — convergence certificate** (hardest case α=0.1293, 1s):
```
rtol 1e-07 → w err −8.4e-10 ; 1e-09 → −5.1e-11 ; 1e-11 → +3.8e-13 ; 1e-13 → −7.3e-15
far-boundary 15→60 /κ: w shift < 3e-13. Error bars on all quoted ω: < 1e-12 (method),
i.e. entirely negligible vs the lattice systematics X10 will face.
```

**Step 3 — spectra.** Current parameters (validation + baseline; α = g²Q_N/4π, m=1.5):
```
Q_N=650 (α=0.1293, a₀=5.16):  ω_1s=1.487187 (ε=1.281e-02) ⟨r⟩=7.54
   ω_2s=1.496821 (3.18e-03) ⟨r⟩=30.5 ; ω_2p=1.496857 (3.14e-03) ⟨r⟩=25.7
   ω_3s=1.498593 (1.41e-03) ⟨r⟩=68.9 ; ω_3p=1.498603 (1.40e-03) ⟨r⟩=64.3
Q_N=921 (α=0.1832): ω_1s=1.473683 (2.63e-02) ⟨r⟩=5.18 — 4% relativistic shift vs hydrogenic
Q_N=210 (α=0.0418): ω_1s=1.498688 (1.31e-03) ⟨r⟩=23.9
Q_N= 90 (α=0.0175): ω_1s=1.499770 (2.30e-04) ⟨r⟩=57.1
```
**X10 ringdown predictions at the §2C design points** (these are the numbers a run must hit):
```
PATHFINDER g=0.05 Q_N=205 (α=0.0409, a₀=16.3, ρ=4.97):
  1s: ω=1.4987428  ε=1.257e-03  ⟨r⟩=24.4      2p: ω=1.4996863  ε=3.137e-04  ⟨r⟩=81.5
  2s: ω=1.4996859  ε=3.141e-04  ⟨r⟩=97.7      3p: ω=1.4998606  ε=1.394e-04  ⟨r⟩=203.7
CARBON g=0.0229 Q_N=650 (α=0.0272, a₀=24.5, ρ=4.98):
  1s: ω=1.4994446  ε=5.554e-04  ⟨r⟩=36.7      2p: ω=1.4998613  ε=1.387e-04  ⟨r⟩=122.5
  2s: ω=1.4998612  ε=1.388e-04  ⟨r⟩=147.0
```
Caveats bound to these predictions: (i) they assume the pure-Coulomb far field — inside r ≲ 1.5·r_half the mode sees the core's Vt curvature and the local clock shift g·a₀(r) (weff0 on the branch: 1.356 at the cap), so **expect the measured 1s to be shifted from these numbers by O(f_inside)** — 0.4% of charge inside for the carbon design, ~2% for the pathfinder wait — (pathfinder f_inside(r_half=3.28, a₀=16.3): 1−e^{−x}(1+x+x²/2), x=0.40 → 0.8%); (ii) l-degeneracy splitting (2s vs 2p: 4×10⁻⁷) is far below lattice resolution — report n-levels only; (iii) in-box the spectrum is additionally discretized by L (box modes at Δω ~ (π/L)²/2m-class); choose L so the 1s binding ≫ box splitting: satisfied at the design L for n=1, marginal for n≥2.

### 2E. Three-U(1) verification (read-only code inspection of `sfa/sim/scp_sim.c`)

**Claim (a) — potential couples only through s = Π_a|Φ_a|², phase-blind per component: VERIFIED.**
- Ungauged complex force (`compute_forces_complex`): s built at :951–954 (`s2_a = u_a²+v_a²`, `s = s2_0*s2_1*s2_2`); Vt′(s) = (μ/2)/(1+κs)² at :956; force on Φ_a = `−2·Vp·Φ_a·prod_rest[a]` with prod_rest = Π_{b≠a}|Φ_b|² at :958, applied identically to u_a (:966–968) and v_a (:975–977) — a real scalar multiplying the component's own (u,v): invariant under independent phase rotation of each Φ_a.
- Gauged path (`compute_forces_complex_gauge`): identical structure at :1160–1166, :1196–1201. Higgs-bag optional term (v71) `−HKAP·hh·Φ_a·prod_rest[a]` (:1198, :1201) is likewise real-scalar × Φ_a — flavor-preserving. `cfg_validate` (scp_config.h:439–467) hard-rejects every legacy real-mode coupling in complex mode — no other matter nonlinearity exists.

**Claim (b) — gauge couples to the TOTAL current: VERIFIED.**
- One link phase for all components: :1143–1145 (`cP,sP` from `link_c/link_s[d]`, no a-index) — per-component charges exist only in the *second* gauge sector (`qg2[a]`, :1069, :1256), which is off in standard configs (`gauge2_mode=0`).
- The E-kick current sums components: :1234–1238 `for (a=0;a<3;a++) Jg += fu[a]*TIp[d][a] − fv[a]*TRp[d][a]` (+Θ terms) — the gauge field cannot see flavor. Gauss residual: `divE − g·ρ_EM` at :2170–2176; E-update `E += hdt·K` and link drift `th_dot = −g·a·E` in `verlet_step` :1609–1613, :1641–1650.

**Claim (c) — at η=0 nothing else couples component phases; η>0 breakers identified: VERIFIED.**
- The **only** component-mixing terms are the η couplings: ungauged curls (∇×Θ)_a = ε_abc∂_bΘ_c entering Φ_a's force at :965/:974 (and (∇×Φ)_a into Θ at :983/:990); gauged covariant curls at :1181–1184 entering :1197/:1200/:1202–1203; and the η seagull in the gauge current (T12−T21, :1239–1249). All carry the ε_abc structure across the component index: they preserve only the joint diagonal U(1) (Φ and Θ co-rotating), breaking the three separate U(1)s. At η=0 every one of these terms is absent.
- Boundary: the absorbing sponge (`apply_damping` :1440–1470) multiplies (u̇_a, v̇_a) by the same factor — charge removal is per-component-proportional; no mixing even at the boundary.

**Discrete update: is Q_a conserved at η=0?** The semi-discrete system (space-discretized, time-continuous) has three exact per-component U(1) symmetries (all force coefficients above are real, phase-equivariant, common to (u_a,v_a)); `verlet_step` (:1583) is symplectic kick–drift–kick whose map commutes with per-component rotation, so the discrete Noether charge is conserved to the integrator's floor with no secular term (the project's measured floors: gauged total-Q exact-by-construction per v69 SPEC; ungauged Q drift at the boundary-absorption floor, v72 §3 — and the comment at :2068 "periodic-box bath runs where total Q is trivially conserved"). **Caveat for X8 instrumentation:** the diag computes only the component-**sum** (qp_v += u·v̇ − v·u̇ over a, :2089–2100 → Q_phi at :2106); per-component Q_a is *not* written to diag.tsv (thp1u/v etc. are θ point-probes, :2136–2145; Q_flux is the Gauss charge in a centered cube, :2229–2237). Per-flavor ledgers exist only via SFA-frame tools.

**Consequence for X8 design.** Protected: the *vector* (Q_0,Q_1,Q_2) at η=0, exactly, in bulk and through the sponge. Not protected: localization of the residual (see §1.6 — free flavored waves are legal carriers; single-component lumps cannot self-bind because the potential force on Φ_a carries the factor Π_{b≠a}|Φ_b|², :958). X8 must therefore (i) run at η=0 strictly, (ii) score on per-component ledgers from SFA frames (`sfa_qcomp`-class cadence), (iii) pre-register the two legal outcomes of §1.6 (dispersed flavored radiation vs two-flavor Coulomb cloud) as PASS variants of the block, with F9 firing only if the Q_a ledger itself fails.

### 2F. Beat-suppression derivation (object-level T1) — the X7 prediction

**Step 1 — setup and the exact cross-term structure.** Two balls k=1,2, standard seeds (three components co-phased within each ball), real profiles: Φ_a^(k) = f_k(x)e^{i(ω_k t+δ_k)}. In overlap, per component (all a identical):
|Φ_a|² = P + X, P ≡ f₁²+f₂², X ≡ 2f₁f₂cos β, β = Δω·t + δ.
The kernel's exact potential argument (code-verified §2E) is the **product**:
**s = (P+X)³ = P³ + 3P²X + 3PX² + X³** — not the single-field |Φ|⁴ expansion; the 3-component structure triples the first-order coefficient and generates the rectified 3PX² term.

**Step 2 — the two regimes.** Tail overlap ⇒ κs ≪ 1 ⇒ Vt(s) ≈ (μ/2)s, μ<0. Interaction density Δv ≡ Vt(s) − Vt(s₁) − Vt(s₂):
- **Δω = 0 (static):** Δv = (μ/2)[ 6P²f₁f₂cos δ + 3Pf₁²f₂² + 12Pf₁²f₂²cos²δ + 8f₁³f₂³cos³δ ]. Leading (f₂≪f₁ near ball 1): **(μ/2)·6f₁⁵f₂cos δ** — first order in overlap, sign ∝ cos δ: co-phase attracts (μ<0), anti-phase repels, quadrature ≈ zero. Matches the v70 measurement (fusion at Δφ=0 through Coulomb ×2.7; pure Coulomb at π/2, ratio 0.94; enhanced escape at π).
- **Δω ≠ 0 (beat-averaged):** ⟨X⟩=⟨X³⟩=0, ⟨X²⟩=2f₁²f₂² ⇒ **⟨s⟩ − s₁ − s₂ = 9Pf₁²f₂²** ⇒ ⟨Δv⟩ = (μ/2)·9Pf₁²f₂² < 0: a **phase-independent residual attraction, second order in overlap (range e^{−2κD}, κ_k=√(m²−ω_k²))**. It does *not* vanish. Residual/coherent ratio ≈ 1.5(f₂/f₁)/|cosδ| ~ 1.5e^{−κD}: at D≈8–11, ~2–3% — consistent with (and hidden inside) v70's quadrature ratio 0.94±.
- The **gauge Coulomb** interaction (charge-density × charge-density) carries no relative phase and is untouched by averaging at any Δω.

**Step 3 — the predicted suppression law (what X7 measures).** Three regimes for the *net impulse* of the first-order force F₀cos(Δωt+δ) over an encounter of duration τ:
1. **Coherent** (Δω·τ ≲ 1): no suppression — J ≈ F₀τ·cos δ.
2. **Averaging** (Δω·τ ≫ 1, no lock): J = (F₀/Δω)[sin(Δωτ+δ)−sin δ] ⇒ **|J/J₀| ≤ 2/(Δω·τ)** — a 1/Δω envelope, *not* exponential, with the floor set by the rectified term + Coulomb.
3. **Locked** (|Δω| < Δω_lock): Adler capture — the relative phase locks near the energy minimum (co-phase for μ<0) and the pair behaves as Δω=0 *attracting*; suppression fails entirely. Δω_lock(D) ∝ |V₀(D)|·(1/Q₁+1/Q₂) with V₀ the coherent interaction energy ⇒ **Δω_lock grows ∝ e^{−κD} toward contact**. Measured anchor: v70 b2_detuned bounds Δω_lock(D=11) < 0.06; the e^{−κ(D−11)} scaling puts the *upper bound* at ≈0.25 by D=8 — X7's primary deliverable should be the measured Δω_lock(D) curve, which no run has yet produced.

**Failure/complication regimes to pre-register for X7:** (i) fusion-speed encounters (τ ~ 20–50 t.u.) vs branch detunings (Δω ≤ 0.09 ⇒ beat period ≥ 70 t.u.) sit in regime 1–1.5 — *real merger encounters are barely averaged*, matching v70's observed 4–10× surge at anti-phase passage rather than quiet suppression; (ii) parametric response if Δω lands on an internal mode of either ball; (iii) at saturated contact (κs ≈ 2.5) the knee suppresses Vt′ by (1+κs)² ≈ 12 — the phase force flattens and *density* (gradient+mass) + Coulomb terms dominate; (iv) flavored partitions: with per-component clocks ω_a^(k), each factor gets its own β_a; first-order terms beat at each Δω_a while the 9-type rectified DC terms survive per pair — detuning protection never removes the DC floor.

**Bottom line for H7-T1 at object level:** gating exists but is *soft* — first-order-only, algebraic (1/Δωτ) not exponential, defeated by locking at small Δω and by the DC floor + Coulomb always. FULLNESS §3's "interact weakly at contact" and §5.2's detuning protection must be re-stated in these terms.

---

## 3. Updated lock & falsifier status

| Lock | Was | Now (this analysis) | Evidence |
|---|---|---|---|
| D1 (vacuum ref for δE) | open (v84) | unchanged; still needed for any E-claims | — |
| D2 (locality quantum) | open (v84) | **mooted** by v85's H-axioms (unit = one cycle); fold into D5′ | §1.1 |
| D3 (entanglement measure) | open (v84) | demoted with A3 → structure witness; keep for QFI thread only | — |
| D4 (light seed) | open (v84) | superseded by SHELL alternative + §1.6 flavor design | §1.6/1.7 |
| D5 (units of M) | **blocking** (council) | **closed as units** by M=Qω/c²; **replaced by D5′ (new): calibration** — Qω vs ∫T₀₀ differ by measured 0.9–4.3%; declare which is inertial and how ε is booked | §2A |
| D6 (Bell scope) | out of scope | unchanged — honestly held | §1.8 |
| D7 (floor certification) | "dt/dx scaling" | **re-specified**: 4 axes (dt², dx², L-exponential, seed f64); add BC-limited category; verdicts are "not shown TRANSIENT at R" | §1.4 |
| D8 (composite ω) | "Δω bookkeeping" | **frozen fix available**: integral form ΔE=∫ω dQ (exact to 0.26% in-kernel; product form fails 33% on c6) | §2A step 5 |
| D9 (slit mask) | open | untouched (no new evidence) | — |
| D10 (kernel-v4) | gated on F1b | unchanged; note X1's "alignment" must be defined vs response spectrum first or F1a/F1b are unreadable | §1.3 |
| D11 (F(Q) freeze) | curvature candidate | **frozen proposal delivered**: F(Q)=[ω(Q_bot)−ω(Q)]/[ω(Q_bot)−ω_fold]; curvature candidate rejected (ill-conditioned) | §2A step 4 |
| D12 (shell scale sep.) | "redesign needed" | **CLOSED**: carbon ratio measured 1.04; two viable configs delivered (g=0.05/Q≈205; g=0.023/Q=650); design law ε₁=1/(2mρ²R²); α ceiling 0.17–0.21 | §2C |
| D13 (Pauli/occupancy) | open | open; sharpened — gross (n,l) levels now predicted (§2D) so X10 is scorable without statistics claims; §5.3's "−6" must be restated as −Q_N | §1.7 |

| Falsifier | Status change |
|---|---|
| F1a/F1b (X1 gating/quantization) | Untested; **design defect flagged** — propagating radiation is always above-gap (Δω>0 strictly); define "aligned" against the ball's linear response spectrum (X10 machinery) before running | 
| F2 (orbit continuum, X5) | Unchanged — remains the sharpest H5 test; state with orbiter Q fixed |
| F3 (no floor convergence, X2) | Unchanged; X2 protocol now concrete (§1.4); cheap first target = v72 r5_long L-scan |
| F4 (wrong λ, X3) | **Tightened**: λ = 2πQ(1+ε)/p — the 3–4% ε systematic must be inside the kill window or F4 will false-fire |
| F5 (Born rate law) | Unchanged (survivable) |
| F6 (c12 continuous shedding) | **Part-tested: AMBIGUOUS.** Episodic (two bursts 185/330 ≈ 2.06/3.68 minimal closures — non-integer), no steps >19.4/2.5 t.u.; cluster-resolved verdict BLOCKED on archived SFA | 
| F7 (merge boundary vs F(Q), X6b) | Untested; **prior correction**: co-phase near-critical pairs should merge-then-evaporate (c12 precedent), not bounce; put the bounce prediction on detuned/anti-phase pairs (X7 variable) |
| F8 (no Δω suppression, X7) | Untested; **prediction upgraded** to the 3-regime law (coherent / 2/(Δωτ) / locked) with measured v70 anchors; correct prior is v70 b2_detuned, **not** v83/e1 sign runs (mis-citation) |
| F9 (annihilation despite mismatch, X8) | Untested; **re-worded** — legal completion channels exist (flavored radiation; two-flavor cloud); F9 should fire only on Q_a-ledger failure |
| F10 (cloud drains, X11) | Untested; **risk raised** — capacity protection is void (annihilation exothermic regardless of F); expect pos1-style slow drain unless flavor-protected; add beat-heating monitor |
| F11 (no discrete response spectrum, X10) | Untested; **predictions delivered** (§2D tables with machine-precision error bars; lattice caveats bound) |

---

## 4. Recommended next actions (ranked)

**Runs now justified (all config-level, no kernel changes):**
1. **X2-cheap (D7 certification pilot):** rerun v72 r5_long at L = 15/19/23 (N scaled to hold dx), dt_factor ×{1,½}, f64 seeds — one weekend of CPU; certifies the η-soliton BC-limited-CLOSED or TRANSIENT and validates the 4-axis fit before any other rung depends on the classifier.
2. **X10-pathfinder (SHELL go/no-go):** g=0.05, Q_N≈205, N=384, L=55, T=2×10³ (existence+profile), scored against §2D's 1s ω=1.49874, ⟨r⟩=24.4. Cheapest decisive new-physics rung in the stack: F11 kills SHELL, a pass funds everything downstream.
3. **X7 (Δω_lock(D) curve):** pairs from branch profiles (Δω up to 0.09) at D = 6/8/10, τ controlled by impact parameter; measure impulse vs Δω and the lock boundary; scores F8 with the §2F law. Reuses `gen_qball_pair`/`gen_sfa_pair` unchanged.
4. **X8 (flavor block, re-worded):** η=0, mismatched conjugate pair, per-component ledgers from SFA frames at ≤2.5 t.u. cadence; pre-registered outcomes per §1.6/§2E.

**Analysis/design still needed before running:**
5. **D5′ decision (paper only):** declare the inertial-mass convention (recommend: M ≡ ∫T₀₀/c² = Qω(1+ε)/c² with ε(Q) tabulated from §2A; H3/H4 restated as thin-wall asymptotes) and adopt integral-form D8. Without this, X2's CLOSED verdicts on fusion states are unscorable.
6. **X1 redesign (paper only):** define T1 "alignment" against the parked ball's linear response spectrum (computable with the §2D machinery generalized to the ball background — a bound-state solve in the ball's own potential well) before any X1 sweep; otherwise F1a is uninterpretable.
7. **F6 closure:** retrieve the archived c12 SFA and run the cluster tracker (§5, cheapest path to a real F6 verdict — no new simulation needed).
8. **X6b re-scope:** fold into X7 (the merge/no-merge boundary is a detuning-and-phase question, not a pure-F(Q) question — §1.5-4).

**Explicitly not recommended:** any kernel-v4 discussion (F1b has not fired — §2.1's own gate); X11/X12 before X10-pathfinder passes; any further weight on "c6 park = fullness sighting" in creed documents (it does not survive §2A).

---

## 5. BLOCKED register

| Step | What is blocked | Exact unblock condition |
|---|---|---|
| §2B step 4 (F6 cluster-resolved verdict) | Distinguishing lumpy vs continuous shedding at the ~90-unit scale — diag Q_phi smears any step over the 10–30 t.u. core→sponge transit | Retrieve `scpsfa:scpsfa/v74/c12_light.sfa` (42.4 GB; `rclone copy scpsfa:scpsfa/v74/c12_light.sfa /space/scp/v74/results/`) and run `bin/sfa_qball_track /space/scp/v74/results/c12_light.sfa 0.30 --tsv c12_clusters.tsv` (122 frames @ 2.5 t.u. suffice). Alternative (requires a sim run — outside this task's authority): re-run c12_light with a per-cluster diagnostic cadence ≤1 t.u. |
| §2B (same, c6 cross-check) | Same discrimination for the c6 tail shed (t=100–170 episode) | Same as above with `c6_light.sfa` (37.6 GB) |
| §2E X8 scoring | Per-component Q_a time series do not exist in diag.tsv (kernel emits summed Q only — scp_sim.c:2089–2106) | Score X8 from SFA frames via `sfa_qcomp`-class analysis at ≤2.5 t.u. snap cadence (config-level), or user-authorized diag extension (kernel edit — **not** requested, **not** performed) |
| §2D lattice-corrected 1s (carbon) | The pure-Coulomb prediction carries an O(0.4–0.8%) core-overlap shift not computable without the ball-background mode solve | Extend the §2D solver with the shooter's f(r), χ(r) background (analysis-only, ~1 day) — recommended before scoring X10 at better than 1% |

Nothing else in Part 2 is blocked: A (branch/F(Q)/Q*), C (D12 designs), D (X10 spectra), E (three-U(1) verification), and F (beat law) are fully executed above.

---

*Method note: every quoted number traces to a shown command/output in §2 or to a cited file:line in the repo. The scratchpad scripts (`branch_drive2.py` + verbatim `gs_fast.py`, `a_analysis.py`, `shed_analysis.py`, `shed_refine.py`, `kg_coulomb2.py`, `c_scan.py`) reproduce all tables from a clean checkout plus `v66/results/profile_omega1.4500.txt`. No repo file other than this one was created or modified; `sfa/sim/scp_sim.c`, `scp_sim.cu`, and `sfa/format/sfa.h` were read, never written.*
