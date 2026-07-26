# Bell rung B0 — working file: FABLE seat
(Owned by the fable agent. Grok may READ this but must not write it.)
**Status: COMPLETE** (all seven deliverables; final VERDICT at end).
Started 2026-07-26; worksheets under `v87/work/fable/`, all run clean end-to-end.

---

## Plan of attack (so GROK can route around me)

I am taking **Door A (relax MI) as the primary construction** — it is the door
O1∧O2 motivates on its own terms — and doing the full deliverable stack there:
explicit model, closed-form E(a,b), symbolic S-max, MI cost, Monte Carlo,
derivation-vs-postulate verdict. I will then work **Doors B and C as
closure/openness proofs against the actual v86 kernel dynamics** (which are a
local hyperbolic lattice PDE — that fact does most of the work), and the
**Tsirelson half** via a one-parameter family that interpolates through the
quantum point, to test whether 2√2 is *structural* or *imported*.

GROK: if you want to divide labour, the deepest thing I will NOT do heavily is
the retrocausal (Door C) constructive side — I intend to close it against the
kernel's initial-value formulation, not build a Wharton-style model.

Worksheet index (each referenced from findings below):
- `ws1_singlet_smax.py` — symbolic max of S for E = −a·b (Tsirelson point).
- `ws2_model_closedform.py` — the explicit Door-A model; E(a,b) in closed form
  (SymPy), including the p-family generalisation.
- `ws2b_model_closedform.mac` — Maxima cross-check of the same integral.
- `ws3_md_bound.py` — the budget laws: tight TV bound S ≤ 2 + 2D (LP-verified),
  minimal-MI curve incl. the global Blahut–Arimoto solve and closed form.
- `ws4_mutual_info.py` — S_max(p) and I(λ : a,b) along the p-family (the
  Tsirelson-is-imported scan).
- `ws5_montecarlo.py` — explicit-PRNG Monte Carlo of the model.

---

## Findings

(numbered as F1, F2, …; status tags [D]/[M]/[P]/[C] per the brief. Here
[M] = "measured" means numerically verified to stated precision.)

### F0. Which assumption the ontology relaxes — the pre-selection argument

**Claim: O1 ∧ O2 admit Door A naturally, admit Door B only by adding structure
the kernel does not have, and do not admit Door C at all for the present
dynamics.** [D from the ontology + kernel facts; argued, not proved, at this
stage — proofs/constructions below.]

- O1 (definite substrate, epistemic probability) is exactly assumption (R):
  a definite λ with ρ(λ) epistemic. O1 *keeps* (R). So the relaxation must be
  (L), (MI), or the temporal-ordering premise.
- O2 (field monism: choosers and detectors are fabric) attacks (MI) directly:
  ρ(λ | a,b) = ρ(λ) is an assumption about the *joint* configuration of
  source-fabric and chooser-fabric, and in a monist theory both are
  functionals of one field history with a shared causal past. (MI) is then a
  dynamical claim requiring proof, not an axiom. This is the standard
  observation behind superdeterminism/measurement-dependence (Hall 2010;
  Brans 1988 — C. H. Brans, Int. J. Theor. Phys. 27, 219 (1988), identifier
  UNVERIFIED, flagging per brief).
- The v86 kernel is a **local** lattice field theory: finite stencil, finite
  timestep, strictly bounded signal cone. The preferred frame it exhibits
  (1–5% group-velocity anisotropy, v86 THEORY A2) is a *dispersion* anomaly,
  not a nonlocal channel: it changes signal speed, it does not connect
  spacelike-separated regions. **Door B requires a genuinely nonlocal
  guidance structure (deBB velocity field depends on the distant
  configuration NOW on a preferred slice); a preferred frame is necessary
  for that but nowhere near sufficient.** The kernel has the frame and lacks
  the channel.
- Door C requires the fabric history to be fixed by two-time boundary data.
  The kernel is an explicit initial-value integrator (leapfrog). There is no
  formulation of the present dynamics in which future boundary data
  constrains λ at the source. Door C is closed *for this realization* — not
  for field monism in general. (Details later.)

### F1. Tsirelson point for the singlet correlator, symbolically  [D]

Worksheet `ws1_singlet_smax.py` (runs clean). For E(a,b) = −a·b:
max |S| = 2√2, attained at b ⊥ b′, a ∥ (b+b′), a′ ∥ (b−b′). Two routes:
a 200-start numeric global maximisation (agrees to 1e-12), and the exact
vector proof S = −[a·(b+b′) + a′·(b−b′)] ≤ |b+b′| + |b−b′| ≤ √(2·4) = 2√2
(Cauchy–Schwarz; the |b±b′| identity checked symbolically). This is the
target the Door-A model must hit and must NOT exceed at its quantum point.

### F2. The explicit Door-A model, closed form  [D, verified to 8e-13]

Worksheets `ws2_model_closedform.py` (SymPy) and `ws2b_model_closedform.mac`
(Maxima cross-check). **The model (p-family):**

- **λ ∈ S²** — one unit vector per run: the uncontrolled fabric degree of
  freedom (the "car phase" of O1's analogy).
- **Responses, deterministic and LOCAL:** A(a,λ) = sgn(a·λ),
  B(b,λ) = −sgn(b·λ). A never sees b; B never sees a. (L) and (R) both hold.
- **Event source (the MI relaxation):** ρ_p(λ|a,b) = c_p |λ·b|^p with
  c_p = (p+1)/4π. The source's λ-distribution is tilted toward the B-side
  analyser axis. p = 0 is the uniform (MI-respecting) model; p = 1 is the
  Degorre–Laplante–Roland density (PRA 72, 062314 (2005) — identifier
  believed right, UNVERIFIED; everything is re-derived from scratch here so
  nothing rests on it).

**Closed form, derived by machine-checked steps** (azimuthal sgn-integral →
Leibniz differentiation → Beta-function integral; each step symbolically or
30-digit verified):

    E_p(θ) = −1 + ((p+1)/π) · B((p+2)/2, 1/2) · ∫₀^θ sin^p φ dφ

with the anchors:
- **p = 0: E = −1 + 2θ/π** — the classic sawtooth; local; |S|max = 2.
- **p = 1: E = −cos θ** — the EXACT singlet correlator, hence |S|max = 2√2
  by F1. A deterministic, parameter-independent-response, no-signalling
  model sitting exactly on the quantum singlet.
- **p → ∞: E → −sgn(cos θ)** — the step correlator; |S|max → 4 (PR-box).

Verified against the definition at 20 (p,θ) points to < 8e-13 (via the
separately-verified azimuthal reduction; the naive 2D quadrature of the
sign-discontinuous integrand stalls at 1e-4 and was replaced), and E_p(π) = +1 for all p (Beta identity, symbolic). Operational
no-signalling verified: ⟨A⟩ = ⟨B⟩ = 0 independent of the far setting (ρ is
even in λ, responses odd), so the hidden-level measurement dependence never
surfaces as a signalling channel. [D]

### F3. The cost bound, and it is TIGHT  [D + M]

Worksheet `ws3_md_bound.py`. Two quantifications of "how much MI-violation
buys how much S":

**(a) Total-variation budget — analytic and tight.** For any deterministic
responses and per-setting-pair densities ρ_xy, and any reference density m:
every deterministic strategy has CHSH value ±2, so

    S ≤ 2 + Σ_xy ∫|ρ_xy − m| dλ = 2 + 2D,   D ≡ Σ_xy TV(ρ_xy, m).

An exact LP over the 16 deterministic strategies (m = the actual mixture)
achieves D = (S−2)/2 for every S ∈ [2,4]: **the bound S ≤ 2 + 2D is tight.**
S = 2√2 costs exactly D = √2 − 1 ≈ 0.4142 of summed TV; S = 4 costs D = 1.

**(b) Mutual information (Hall-style) — CORRECTED after cross-seat check.**
My first SLSQP minimisation reported I_min(2√2) = 0.0638 bits. **That was a
local minimum** (softmax parametrisation is non-convex even though the
underlying problem is convex): GROK's GOOD-tilt construction (their G5)
achieves 0.046274 bits, and when I evaluate their construction inside my own
formulation it checks exactly. I then solved the convex problem globally by
the Blahut–Arimoto fixed point ρ_xy(s) ∝ m(s)·exp(ν·c_xy A_x B_y) (ws3 §4),
which converges to the global optimum, and found the **closed form for the
entire optimal curve**:

    I_min(S) = D_KL( (1 + S/4)/2  ‖  3/4 )   bits,   S ∈ [2, 4]

verified at 7 points to 1e-14 against the fixed point:

| S target | I_min (bits) = closed form | Pinsker lower bd |
|---|---|---|
| 2.0 | 0 | 0 |
| 2.5 | 0.016006 | 0.0113 |
| **2√2** | **0.046274 = D_KL((1+1/√2)/2 ‖ 3/4)** | 0.0309 |
| 3.0 | 0.069593 | 0.0451 |
| 3.5 | 0.176808 | 0.1014 |
| 4.0 | **log₂(4/3) = 0.415037** | 0.1803 |

This CONFIRMS GROK's G5 value and slightly extends it: GROK exhibited the
tilt family and its I(t); the Blahut–Arimoto solve shows that family is the
**global optimum at every S**, hence I_min(S) = D_KL((1+S/4)/2 ‖ 3/4) is the
exact price list for the 4-setting CHSH problem, ends included. So the
Tsirelson point costs exactly **0.046274 bits of source–setting mutual
information per run**, and the PR box costs log₂(4/3) ≈ 0.415 bits — not
1 bit. (Hall PRL 105, 250404 (2010) reports ≈0.0663 bits for the FULL
singlet at all settings — a strictly harder task, and properly above our
CHSH-only 0.0463. Ordering consistent; his number UNVERIFIED against the
primary source this session.) The concrete p = 1 model of F2 is NOT
MI-optimal: at its optimal settings (b ⊥ b′) it spends I = 0.2018 bits
(ws4), ≈ 4.4× the optimum — the price of its clean geometric form and of
correlating with only one side's setting.

### F4. The Tsirelson question: in this family the bound is IMPORTED  [D/M]

Worksheet `ws4_mutual_info.py`. S_max(p) and I(p) along the family:

| p | S_max | I (bits) |
|---|---|---|
| 0 | 2.0000 | 0 |
| 0.5 | 2.4849 | 0.0731 |
| **1** | **2.8284 = 2√2** | 0.2018 |
| 1.5 | 3.0817 | 0.3296 |
| 2 | 3.2732 | 0.4427 |
| 3 | 3.5355 | 0.6191 |
| 10 | 3.9725 | 0.9734 |
| 20 | 3.9994 | 0.9993 |

**S_max(p) is smooth and strictly increasing through p = 1. Nothing in the
relaxation mechanism distinguishes the quantum point.** The measurement-
dependence door, by itself, has no Tsirelson wall: the wall is put in by
choosing p = 1, i.e. by matching quantum mechanics. This is the honest
negative result the brief anticipated in §4: the model reaches 2√2 but does
not *derive* it. (Note also the pretty coincidence S_max(3) = 5/√2 — not
pursued.) The only structural anchors in the family are p = 0 (zero cost,
local bound) and p → ∞ (1 bit, algebraic maximum 4).

### F5. Monte Carlo — deliverable 4  [M]

Worksheet `ws5_montecarlo.py`. Explicit PRNG (numpy PCG64, base seed
20260726), N = 2×10⁷ per correlator, exact inverse-transform sampling of
ρ_p (|λ·b| component: |z| = U^(1/(p+1)), sign fair, azimuth uniform, rotate
to b). Deterministic responses. Results:

| run | S (MC) | S (closed form) | deviation |
|---|---|---|---|
| p=1, quantum angles | 2.828068 ± 0.000447 | 2√2 = 2.828427 | 0.80 σ |
| p=0, its optimal angles | 1.999599 ± 0.000447 | 2 | consistent, ≤2 |
| p=8, step-adapted angles | \|S\| = 3.900662 ± 0.000447 | >2√2 | overshoot confirmed |

The symbolic, closed-form, and Monte-Carlo values agree within stated
precision at every tested point. The Maxima worksheet
(`ws2b_model_closedform.mac`) independently proves the key derivative
identity **exactly** (radcan residual 0), E₁(θ) = −cos θ, S(canonical) =
2^{3/2}, and stationarity of the canonical angles. [M]

### F6. Door B is CLOSED for the present kernel — the frame is not a channel  [D]

The claim to close: "the v86 preferred frame (1–5% group-velocity anomaly,
THEORY A2) is the structure Door B requires." It is not, and the two things
should never be conflated:

1. **What Door B needs.** In de Broglie–Bohm (and the double-solution
  picture), the guidance of subsystem 1 depends on the *instantaneous*
  configuration of subsystem 2 on a preferred simultaneity slice — an
  influence with unbounded speed, exercised at spacelike separation, every
  run. The preferred foliation is only the *bookkeeping* that makes this
  well-defined; the physics is the instantaneous dependence.
2. **What the kernel has.** `scp_sim` is a finite-stencil, finite-timestep
  lattice PDE: the state at (x, t+dt) depends only on the state within one
  stencil radius of x at times t, t−dt. Signal speed is strictly bounded
  (and measured — that is what the A2 group-velocity numbers ARE). The
  anisotropy/dispersion anomaly moves the *value* of the bound by 1–5%; it
  does not create any dependence outside the cone. Formally: for any two
  regions R₁, R₂ whose measurement events are outside each other's lattice
  cones (with the anomaly-corrected speed), the outcome in R₁ is a function
  of its own cone's data — parameter independence (L) holds *by
  construction of the update rule*.
3. **Consequence.** For any experiment realizable in the present kernel in
  which the settings are injected locally (fresh entropy after last common
  contact) and the measurements are cone-separated, Bell's theorem applies
  with (L) and (MI) both holding, so **|S| ≤ 2**. Door B in this kernel
  would require adding an instantaneous nonlocal term — a kernel change,
  policy-gated, and a *new postulate*, not a reading of the existing frame.
  And even then the 2√2 would come from importing the quantum equilibrium
  measure |ψ|², i.e. the Tsirelson half is imported in Door B exactly as in
  Door A. **The v86 "defect" frame is necessary-but-nowhere-near-sufficient
  for Door B; recording it as Door-B structure would be wrong.** [D]

### F7. Door C is CLOSED for the present kernel — no two-time structure  [D]

The kernel is an explicit initial-value integrator (leapfrog): the
configuration at the source event is a deterministic function of *past*
Cauchy data alone. A retrocausal/all-at-once model requires the source
hidden variable to be constrained by *future* boundary data (the settings) —
i.e. the history is selected by a two-time boundary-value problem (Wharton &
Argaman, Rev. Mod. Phys. 92, 021002 (2020) — identifier from the brief).
No such formulation of this kernel exists, and O1's own language ("evolving
by definite dynamics", A1 "exists fundamentally in motion") is
Cauchy-flavoured. Operationally, any retrocausal λ–settings dependence is
statistically indistinguishable from Door-A measurement dependence
(ρ(λ|a,b) ≠ ρ(λ)), so all the F3 cost accounting applies unchanged to Door C
— it buys no discount on the bits and adds a reformulation cost (two-time BC)
the programme's dynamics cannot express. **Closed for this realization;
philosophically available to a block-universe reading of field monism, but
that reading changes the theory, not just the interpretation.** [D]

### F8. The crux (deliverable 5): the structure is POSTULATED, not derived —
### and the measured fabric dynamics actively work against deriving it  [D/C]

What Door A needs, made physical: the source's uncontrolled DOF λ (natural
fabric candidate: the emitted pair's internal U(1) phase/clock at emission —
per A5 these are real, measurable DOF) must carry ≥ 0.046 bits/run of mutual
information with the *future outputs* of both setting choosers, run after
run. Three observations:

1. **Local routing only.** By F6/F7 the only route in this kernel is the
  common causal past. That is logically possible (this is superdeterminism)
  but note what it quantifies to: by the data-processing inequality, if the
  chooser applies ANY function f to its input entropy, correlating λ with
  the output f(input) requires the fabric to have correlated λ with the
  relevant input entropy at ≥ the same MI — for EVERY entropy source the
  chooser might consult (thermal noise, a PRNG, a photon from the causal
  edge). The postulate is not one correlation but a correlation with the
  entire algebra of possible setting-generating procedures. Nothing in the
  measured kernel (v66–v86) produces or maintains targeted correlations of
  this kind; the measured dynamics of this fabric are of the opposite,
  mixing type — Adler-averaging decoherence of detuned clocks, radiation of
  all unclosed flow, thermalization of the condensate. [C, anchored to
  measured behaviours]
2. **What was derived vs postulated in this rung.** Derived: the response
  functions' locality, the closed form E_p, the tight budget law
  S ≤ 2 + 2D, the exact price list I_min(S) = D_KL((1+S/4)/2 ‖ 3/4).
  Postulated: the tilt ρ_p(λ|b) itself — the map from settings to source
  distribution. The tilt is a *parametrisation of the answer*: it is chosen
  to make E come out right, precisely the thing the brief warns about.
3. **What a derivation would require** (concrete, so a future rung can
  attempt or kill it): (a) identify λ in kernel variables (emission-time
  component phases); (b) exhibit a dynamical mechanism by which
  chooser-fabric evolved from shared past data biases the emission phase
  with an effective |λ·b|-type tilt; (c) show the induced MI ≥ 0.046
  bits/run survives arbitrary fresh-entropy injection at the choosers; and
  (d) separately derive p = 1 (the Tsirelson selection — F4 shows the MD
  mechanism itself has no preference). I checked the obvious candidate for
  (d): "linear response" does not select p = 1 — the p = 1 density is an
  order-unity modulation, not a perturbative one, and the family's S(p) has
  no structure at p = 1 in any case. **My assessment: (a) is doable, (b)
  contradicts the measured mixing character of the dynamics, (c) is
  quantitatively falsified for the kernel as-is (there is no mechanism at
  all, so I = 0 and S ≤ 2), (d) has no candidate mechanism. Verdict:
  postulated.** [D for the accounting, C for the impossibility-in-practice]

### F9. Falsifiable consequences (deliverable 6)  [D → quantitative predictions]

1. **For the fabric programme (sharpest, in-kernel).** The present kernel,
  run as-is with in-fabric choosers injecting fresh entropy (e.g. chaotic
  scattering subsystems) and cone-separated measurements, MUST give
  **|S| ≤ 2 + 2D_meas**, where D_meas = Σ_xy TV(ρ̂(λ|xy), ρ̂-mixture) is
  *directly measurable* by histogramming the source variables conditioned
  on the later settings. Since no MD mechanism exists in the kernel,
  prediction: D_meas ≈ 0 and |S| ≤ 2. Observing |S| > 2 + 2D_meas in-kernel
  falsifies the locality of the implementation (bug tripwire, same class as
  a Gauss-law drift). Converges with GROK's fingerprint 3.
2. **For any Door-A completion of the theory.** The exact price list
  inverts to S_max(I) = 4·(2η*(I) − 1), η* the solution of
  D_KL(η ‖ 3/4) = I. If the choosers' settings contain F bits/run of
  entropy injected after the last common-past contact with the source,
  the available MI is reduced accordingly and **S must fall below 2√2
  along this curve**: I = 0.02 → S ≤ 2.55; I = 0.01 → S ≤ 2.40;
  I = 0.001 → S ≤ 2.13. A loophole-free CHSH experiment whose settings are
  throttled to fresh-entropy sources and still measures 2√2 at these
  timings pushes the required conspiracy MI back up to 0.046 bits/run
  routed through ever-earlier common past — the cosmic-Bell strategy
  (Handsteiner et al. 2017, PRL 118, 060401 — identifier UNVERIFIED) is
  exactly this test, and any Door-A fabric theory predicts S DEGRADES
  toward 2 when the settings' entropy genuinely postdates the last
  common-contact epoch. Quantitatively: sustaining S = 2.8 needs
  I ≥ 0.043 bits/run; S = 2.5 needs 0.016; S = 2.1 needs 0.0006.
3. **If Door B were instead added to the kernel** (policy-gated): the
  nonlocal channel respects the lattice foliation, so CHSH experiments
  with the A–B axis at different orientations/boosts relative to the
  lattice frame acquire order-(v/c_eff) asymmetries in timing-resolved
  correlations, at or above the measured 1–5% dispersion scale. (Matches
  GROK's fingerprint 2 in kind.)

### F10. Cross-seat notes (GROK engagement)

- **Convergence on the main question is genuine, not converged-for-comfort:**
  both seats independently landed Door A primary, Door B frame≠channel,
  Door C closed under the kernel IVP, Tsirelson not derivable from the
  relaxation. Different model families (GROK: finite 16-strategy GOOD-tilt +
  circle tilt; FABLE: continuous S² p-family with closed-form E_p) reach the
  same structural conclusions — that is the strongest form of the result.
- **GROK's G5 corrected my F3.** Their GOOD-tilt value 0.046274 bits beat my
  SLSQP local minimum 0.0638; I verified their construction inside my
  formulation, then closed the loop with a globally-convergent solve giving
  I_min(S) = D_KL((1+S/4)/2 ‖ 3/4) — which also proves their tilt family is
  optimal at EVERY S (they showed the family; the global-optimality
  statement is the new piece). Their closed form D_KL(η‖3/4) at η=(1+1/√2)/2
  is exact. Credit where due: this is the collaboration working as designed.
- **One quibble with G5's framing:** "0.046274 bits is the apparent minimum"
  — it is the actual global minimum (the problem is convex; Blahut–Arimoto
  certifies it; ws3 §4). GROK can upgrade "apparent" to proven.
- **One addition to G7/G10:** GROK's family table stops at t=1 giving I =
  0.415 = log₂(4/3) at S = 4; worth noting this is *below* 1 bit — the PR
  box needs less than one full bit of measurement dependence, a point Hall
  makes qualitatively (UNVERIFIED) and our exact curve makes sharp.
- **Their G6 report "discretised I(Λ:S) ≈ 0.0246 bits" is WRONG, and I can
  say exactly why it must be and what the true value is.** By the
  data-processing inequality, any continuous-λ model achieving S = 2√2
  induces a distribution on the 16 response strategies with the same
  correlators, so its I(λ:settings) is ≥ the 16-strategy global optimum
  0.046274 — 0.0246 is impossible. I re-computed the true
  I(λ:XY) = avg KL(ρ_xy ‖ mixture) for GROK's exact circle-tilt model on a
  200k grid: it is **0.046274 bits — their model exactly SATURATES the
  global optimum** (the four tilted conditionals mix back to uniform at the
  Tsirelson angles, so KL-to-mixture = KL-to-uniform = D_KL(η‖3/4)). So the
  correction runs the other way: GROK's continuous model is better than
  they claimed (optimal, not merely close), and their 0.0246 discretised
  number is an implementation error in g7's mixture handling. **GROK: your
  circle tilt is MI-optimal; fix the 0.0246 line.**

---

## Deliverable checklist (brief §5)

| # | Deliverable | Status |
|---|---|---|
| 1 | Symbolic model, ρ(λ\|a,b), responses, event source, closed-form E, symbolic S-max | **Done** — F2 (E_p closed form, machine-derived), F1 (S-max: 2√2 with proof), event source stated in ws5/F2 |
| 2 | Bound on max\|S\| | **Done** — model at p=1 reaches exactly 2√2 (F1+F2); family overshoots to 4 as p→∞ (reported as the failure mode it is, F4); tight general law S ≤ 2 + 2D (F3a) |
| 3 | Cost quantified | **Done** — exact: D = (S−2)/2 in TV; I_min(S) = D_KL((1+S/4)/2 ‖ 3/4); 0.046274 bits at Tsirelson; concrete model pays 0.2018 bits (F3) |
| 4 | Monte Carlo | **Done** — F5, seeds and code committed |
| 5 | Derived or postulated | **Done** — POSTULATED, with the derivation requirements itemised and assessed against measured fabric behaviour (F8) |
| 6 | Falsifiable consequence | **Done** — three, quantitative (F9) |
| 7 | Honest negatives | **Done** — Tsirelson not derived (F4); Doors B and C closed for the kernel (F6, F7); F3 self-correction recorded |

---

## VERDICT (FABLE seat, final unless GROK's recheck moves something)

**Which assumption does O1 ∧ O2 relax?** (MI), and only (MI). O1 keeps (R).
The kernel's locality keeps (L) — the preferred frame is a dispersion fact,
not a nonlocal channel (F6). The kernel's initial-value structure closes the
retrocausal door (F7). O2 makes (MI) a dynamical claim rather than an axiom,
which is exactly the opening Door A needs.

**Can a deterministic sub-quantum field theory obeying O1 ∧ O2 reproduce
2 < |S| ≤ 2√2?** Yes — explicitly, in closed form. λ ∈ S², A = sgn(a·λ),
B = −sgn(b·λ), ρ(λ|b) = |λ·b|/2π gives E = −a·b exactly and |S|max = 2√2
(F1, F2), no-signalling at the operational level, Monte-Carlo confirmed
(F5). The exact cost is 0.046274 bits/run of source–setting mutual
information at the optimum (F3, jointly with GROK), 0.2018 bits for this
concrete geometric model.

**Does the relaxation explain the Tsirelson ceiling?** No. The p-family
passes through the quantum point with no structure there (F4); unrestricted
MD runs to |S| = 4 at log₂(4/3) bits. The ceiling must be imported (by
matching QM, as here) or derived from an independent principle no one has
derived from fabric dynamics. Half the problem is solved; the other half is
honestly open, for us and for the field.

**Is the required structure derivable from fabric dynamics?** Not today, and
the measured character of the dynamics (mixing, decohering, radiating) points
the wrong way (F8). In this rung it is POSTULATED. The concrete derivation
requirements (a)–(d) in F8 are the checklist a future rung must meet or
formally kill.

**Bottom line for the programme:** Bell does not force the fabric ontology to
change — it prices it. The price is known exactly now: ≥ 0.046 bits/run of
setting–source conspiracy, routed through the common past, robust to fresh
entropy, plus an unexplained Tsirelson selection. The kernel as measured
pays neither and therefore predicts |S| ≤ 2 for any in-fabric Bell test —
which is itself the cleanest falsifiable statement this rung produces (F9.1).
THEORY_v86 A9's "Bell is closed by theorem" should be restated as: *closed
against derivation, open to parametrisation at a measured price* — the
theorem tells you what to buy, the fabric currently declines to pay.
