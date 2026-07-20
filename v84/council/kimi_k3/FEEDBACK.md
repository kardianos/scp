# Kimi K3 — v84 CREED feedback (Round 1)

**Date:** 2026-07-20
**Seat:** Third-seat design advisor (Kimi K3); distinct skeptical voice from Fable R3.
**Charge:** FEEDBACK_PROMPT.md — Part A forward from creed; Part B backward from QFT.
**Documents read:** v84/CREED.md, v83/CLOSEOUT.md, v73/PROCESS.md §0–2, v76/PROBLEM.md §0–3.1, CONCEPT.md §2–3, §7–8, Fable R3.

---

## 1. Executive stance: PARTIAL ACCEPT — narrower than Fable's

Accept as research basis: **A2** (energy = δ-field, as a *program* — the
vacuum-subtraction move is QFT-orthodox, Casimir-anchored, and the project
already implements it), **A5** (process closure — the one genuinely measured,
quantitative claim in the creed's ancestry; v73 P1–P3), and **A6** (charts as
re-partitions — methodologically clean).

**Reject as curreMemo written: `v84/council/kimi_k3/FEEDBACK.md` (457 lines).

**Stance: PARTIAL ACCEPT (narrower than Fable).** Accept A2 (δ-field, Casimir-anchored), A5 (process closure, already measured), A6; **reject A1 and A3 as stated.**

Distinct skeptical contributions beyond Fable R3:
- **D5 (units of M)** — elevated to the blocking definitional lock. Entanglement is dimensionless; action has units; N-M ∈ [action]⁻¹, which is not mass and cannot become mass without an imported constant. A3's algebra does not support the word "mass."
- **A1 is untestable on the fixed-lattice kernel** (v76 §3.1) — Fable accepted A1; I demote it to aspiration until a field→geometry back-reaction path is authorized (new G5).
- **QFT Part B goes deeper than Fable's "E=mc² is a theorem"**: mass and entanglement are *orthogonal* in QFT for elementary particles (|p⟩ is pure, zero entanglement, nonzero mass), and *redundant* for composites (proton mass is field energy; entanglement is unnecessary as the mechanism). A3 has no QFT home on either side.
- **Ball vs Q-ring at equal (Q,J)** proposed as the cleanest discrimination test — δE-degenerate by construction (~3% K-isomer), maximally different correlation structure.
- **QFI-on-ρ_Q candidate struck pre-audit** (provably constant N-M for coherent-state-like U(1) balls), not audited in parallel.
- **"Positronium" label rejected** under the creed's own §4.1 honesty standard (no quantized annihilation channel); **Higgs tension flagged** (A2 is closer to Higgs than to entanglement).
ct,
and v84's first job is to decide whether either is fixable on paper before any
pipeline is built. If not, A3 is retired and the creed reduces to A2 + A5 —
still a useful program (δ-field energy + process closure), just not a new mass
definition.

**Single strongest objection.** The creed proposes a new definition of mass
that, by its own algebra, is not in mass units, and then demotes the one
measurement (inertia) that *is* in mass units to "optional consistency." That
ordering is backwards regardless of whether A3 is eventually rescued.

---

## 2. Part A — forward from the creed

Assume A1–A6 as written. What follows?

### 2.1 What genuinely becomes possible

**(a) The δ-field energy program is already running and is the creed's
strongest asset.** N-δE = ⟨E⟩_box − ⟨E⟩_vac is what the heavy stack already
measures (E/Q, mass defect, fusion binding, v74 RESULTS). The creed does not
need to build this; it needs to *name* it as the primary energy observable
and stop treating ∫T₀₀ as the definition. Cheap, immediate, costs nothing.

**(b) Process closure (A5) becomes the stability criterion, replacing the
VK-only criterion.** v73 P2 measured ⟨S·r̂⟩ = 0 at every radius for the
stationary state and ≠ 0 for the drain seed, quantitatively. This is a
validated, already-implemented diagnostic (`analysis/flux_profile.c`). Under
the creed, "stable mass-form" *means* "ledger closes over a cycle." This is
the creed's one ready-to-run operational win.

**(c) The reject list (§4) is coherent and should be enforced at proposal
intake.** The v80–v83 failures share a diagnosis (testing lock kinematics,
never mass) and the creed's §4.2 codifies it. This is the most useful text in
the document.

### 2.2 What is underspecified (blocking, and Fable did not flag all of it)

**(i) The units of M — my strongest Part-A objection.** Entanglement entropy
is dimensionless (nats or bits). A "locality quantum unit" with units of
action (the only candidate with a measured value is h_eff = 2πQ, in action
units) gives

    N-M = N-ENT / N-LQ  ∈  [action]⁻¹ .

That is not mass. It is not even mass-like. To reach mass you must multiply
by a quantity with units [mass² · length² · time⁻¹]; the creed names no such
object. The only dimensionful quantities in the theory are the field mass m,
the locality speed c, and the potential parameters (μ, κ). Using any of them
to convert is importing a constant — which may be legitimate, but it must be
*exhibited* and its ontology owned, not hidden inside the word "mass." And
once exhibited, it will almost certainly reconstruct some form of E/c² or
E/(c·something), at which point A3 has smuggled dualism back in through the
constant. **Before any pipeline is built, the creed must either (α) exhibit
the conversion constant and its derivation from A1–A6, or (β) concede that
N-M is a dimensionless relational quantity and that "mass" in A3 is a
metaphor.** Fable notes the units issue in passing (§3.3 item 4); I elevate
it to a definitional lock — call it **D5** — and make it action #1.

**(ii) A1 is untestable on the current kernel.** v76/PROBLEM.md §3.1 is
explicit: a fixed-lattice PDE assumes ambient space exists independently;
fields are values *on* points; T₀₀ is bookkeeping of excitations *on* the
stage. That is the dualist architecture. The creed asserts A1 (one continuum)
but the simulator is dualist. Consequences:
- No experiment on the current kernel can *distinguish* monism from dualism,
  because the kernel cannot express the monist claim (no back-reaction of
  field on geometry).
- Therefore A1 is not a research basis in the operational sense; it is a
  target ontology requiring kernel authorization to test. The creed should
  say so.
- The "warp corollary" (§2.1 optional free-capacity scalar; A4 bullet 4) is
  *the* monist claim and it is unaddressed by any current tool. Flag it as
  post-Stage-5, not Stage 3–4. Do not let A1 drive carrier or binding design
  on the existing stack.

**(iii) A3's uniqueness is unprovable on a one-parameter family.** Along the
measured Q-ball branch (Q ∈ [88, 921], ω ∈ [1.41, 1.47]) every quantity
currently measured (E, Q, ω, E/Q, radius, Θ-fraction) is a function of
(Q, ω) alone. On a one-parameter family, *any* two monotonic functions are
proportional, so N-M will trivially correlate with δE. The creed cannot
demonstrate A3 is non-redundant on this family. It needs (α) a *two*-parameter
family and (β) a pre-registered prediction of where N-M and δE decouple.
Fable's discriminating triple (single / separated / fused) is one form. I
propose a sharper one: **ball vs Q-ring at equal (Q, J)** (v73 §5, CONCEPT
§3). They are degenerate in δE by construction (K-isomer split ~3%) and
maximally different in spatial correlation structure (ring has a hole +
real-space winding; ball is spherically symmetric, internal-only rotation).
If N-M cannot separate these two, A3 is dead — they are the minimal
configuration where A2 and A3 *can* decouple.

**(iv) The light mass-form has no creed-level source.** Fable flags that
conjugation preserves the spectrum. I add: A2 (δ-field) does not help either
— δ-field energy scales with the same charge branch that produces heavy
balls; there is no light branch in δE. Stage 4 needs m_light ≪ m_nucleus
*structurally*, and none of A1–A6 produces a scale hierarchy. This is not a
D4 detail; it is a creed-level gap. The creed should either derive a light
scale from its own axioms or declare Stage 4 blocked on a scale-generation
mechanism and **stop scheduling atom-shell experiments**.

### 2.3 What kills the creed (Part-A view, my gates)

1. D5 (units) is not closable without an imported constant → A3 is metaphor.
2. N-M ∝ δE on every measurable configuration (one-parameter family + ball/ring
   control) → A3 is redundant with A2.
3. N-M disagrees with inertial response (Fable G3) → A3 is not mass.
4. A1 stays untestable on any kernel the project will authorize → monism is
   aspiration, not science; demote to long-term.

### 2.4 How Stage 3–4 should proceed

Assuming the creed survives a D5/units audit (i.e., N-M is reformulated as a
dimensionless relational quantity, or a constant is honestly exhibited):

1. **Metrology first, on the heavy stack** — N-δE (trivial; already running),
   N-PROC (already running, v73), N-ENT (needs D3 ensemble choice), N-LQ
   (freeze h_eff = 2πQ; strike the lattice cell). No new carriers.
2. **Discrimination test on ball vs ring at equal (Q, J)** — cleanest A3-vs-A2
   test; uses only existing seeds. Pre-register the ordering.
3. **Conjugate-pair bound state** — but *drop the name "positronium"* (§3.2d).
4. **Stage 4 stays blocked** on a light-scale mechanism. Do not schedule
   atom-shell experiments.

---

## 3. Part B — backward from QFT

Start from: particles as field quanta; m² as Poincaré Casimir (pole of the
two-point function); E² = p²c² + m²c⁴; entanglement in QFT; microcausality;
renormalization; composite bound states.

### 3.1 What is a plausible structural analog

**(a) "Mass of a composite is field energy" — yes, and the creed should lean
on this, not A3.** ~99% of proton mass is gluon field energy + quark kinetic
energy + trace anomaly, not Higgs-coupled constituent mass. A2's δ-field is
structurally identical to normal-ordering / vacuum subtraction in QFT, and the
heavy-stack program (mass defect, E/Q, fusion binding) maps onto the QCD
mass-generation story cleanly. *This is the part of the creed QFT supports,*
and it is the part the creed demotes in favor of A3. The ordering should be
inverted: A2 is the QFT-anchored claim; A3 is the speculative addition.

**(b) Vacuum subtraction has a precise QFT ancestor the creed should cite.**
The Casimir energy *is* δ-field energy: E_Casimir = ⟨E⟩_plates − ⟨E⟩_vac. It
is finite, regulator-dependent in absolute value, but ratio-stable across
schemes. This is the closest real-QFT ancestor of N-δE and it is *not*
"energy of stuff on a stage" — it is a property of the vacuum arrangement.
The creed should cite Casimir as the existence proof that δ-field energy is a
well-defined monist quantity, and accept that it inherits Casimir's regulator
dependence (Fable's §3.2e point).

**(c) The process-form (v73 P1) is semiclassically respectable.** ∮2E_kin dt
= 2πQ is the adiabatic invariant of the internal rotation; in QFT the Q-ball
charge is quantized in units of ħ and the classical branch is the ℏ→0 limit
of a quantum tower (Bohr–Sommerfeld). Calling h_eff a *candidate process
quantum* is exactly right. Fable made this point; concur.

### 3.2 What conflicts with QFT-tested structure

**(a) Elementary particles have zero spatial entanglement entropy in the
relevant sense — and that is not a defect of the particle.** A single-particle
Fock state |p⟩ is pure; its entanglement entropy with any spatial region is
zero. The mass m_e is a property of the field (the pole of the propagator),
not of any entanglement inventory of "the electron state." In QFT, **mass and
entanglement are independent quantities**: a massive particle can be in a pure
state (zero entanglement) or a mixed/entangled state (nonzero entanglement)
without its mass changing. The creed's central move — *mass = entanglement
per locality quantum* — has no QFT counterpart for elementary particles. It
can only be defended for *composites* (proton = field energy, where the energy
story already suffices and entanglement is unnecessary as the mass mechanism).
So A3 proposes something QFT actively contradicts for elementary mass, and is
redundant for composite mass. This is a deeper objection than Fable's "E=mc²
is a theorem for solitons": it is "entanglement and mass are orthogonal in
QFT, and the creed has not shown they are not orthogonal in this model
either."

**(b) For classical solitons E = mc² is a theorem (Fable's point, concurred).**
If the continuum action is Lorentz-invariant (Cosserat is; the lattice breaks
it at O(Δx²)), boosting a static solution gives E = γMc² with M ≡ E_rest/c².
Derrick virial identities and the measured de Broglie tilt (v70/v73, 1–5% lag
at lattice-artifact level) confirm the model is in this regime. To lattice
accuracy, **inertia = δE/c²** no matter what N-ENT says. The creed cannot
escape this by demoting the inertia test. Fable's order-parameter restatement
is one option; mine (retire A3, keep A2 + A5) is cleaner if D5 is not
closable.

**(c) Mass is a Poincaré Casimir; entanglement entropy is frame-dependent.**
m² = P^μP_μ is a Lorentz scalar. Entanglement entropy of a spatial region is
observer-dependent (Rindler observer sees Minkowski vacuum as thermal). Any
N-M built from equal-time regional data is a rest-frame quantity. The creed
must stipulate the rest frame and derive boost behavior — silent
frame-dependence makes "mass" quietly stop being a scalar. This is not a
minor caveat; it belongs in A3's text, not in D3's appendix.

**(d) "Positronium analog" is an imported label — exactly what §4.1 forbids.**
Positronium's defining physics is the *quantized annihilation channel*
e+e- → γγ, which in QFT arises because e+ and e- are quanta of the *same*
Dirac field; particle number is not conserved (only charge). The conjugate
Q-balls are *separate classical solitons* of a U(1) theory with a conserved
charge Q; while a soliton-antisoliton pair with total Q = 0 can in principle
dissolve to field radiation classically, there is no *quantized* decay to a
discrete photon-analog final state — the gauge field A exists, but the
discrete final-state structure of QED pair annihilation does not. Calling
the conjugate-pair bound state "positronium" imports quantum field-theoretic
structure the classical model lacks. This is the same category error the
creed's §4.1 rejects elsewhere ("multi-fab multiplet L as 'light electron'").
Rename it "conjugate-pair bound state" and reserve "positronium" for a state
that actually has a quantized annihilation channel. This is not a quibble; it
is the creed's own honesty standard applied to itself.

**(e) The Higgs tension the creed does not address.** §4.1 rejects "Cosserat
multiplet Higgs as free light mass-form." But the Higgs mechanism is *the*
QFT example of mass-from-vacuum-structure: m_e = y_e v/√2, where v is the
vacuum condensate. This is structurally a δ-field move (the vacuum is
nontrivial; the particle's mass is how it couples to that arrangement). **The
creed's A2 is closer to the Higgs story than to the entanglement story.** The
creed rejects the Higgs tool without explaining why its δ-field is different
from a vacuum-VEV shift. Either (α) explain the distinction (probably: Higgs
mass is from coupling to a *spatially uniform vacuum* condensate; creed mass
is from a *localized* non-vacuum arrangement — but this is a statement about
localization, not about the mass mechanism), or (β) concede A2 is Higgs-like
and the disagreement with §4.1's rejected method is about *what condensate*,
not about the mechanism. The current text papers over this.

**(f) Renormalization / scheme dependence (Fable's point, concurred).** N-δE
and N-ENT both depend on (N, L, Δx, snapshot cadence). A serious comparison
requires a continuum-trend check (two lattice spacings; check that orderings
survive even if absolute values drift). The creed should pre-declare *which*
quantities are claimed scheme-independent (probably: orderings and ratios
within a fixed scheme; definitely not: absolute N-M values).

### 3.3 Operational definitions needed for a serious QFT comparison

Minimum set, distinct from Fable's where I can add:

1. **Resolve D5 (units) before anything else.** Exhibit the constant that
   converts [action]⁻¹ to mass from A1–A6 alone, or concede N-M is
   dimensionless-relational. Without this, no QFT comparison is meaningful —
   you cannot compare a quantity in [action]⁻¹ to a mass.
2. **Ensemble declaration** (Fable's point): Gaussian-fluctuation covariance
   is the one with a QFT counterpart (von Neumann entropy of Gaussian states
   from symplectic spectra of two-point functions). Cycle-ensemble MI is
   classical and analogical only. State which.
3. **Rest-frame declaration** + boost behavior as a measured property (Fable's
   point, but creed-text level, not appendix).
4. **Double subtraction** (Fable's point): N-ENT must itself be a δ-quantity
   (object minus vacuum, same partition geometry) — A2's contrast principle
   applied to A3's numerator.
5. **A two-parameter discrimination family**, not just a (Q, ω) scan. Ball vs
   ring at equal (Q, J) is the minimal one and it exists in the codebase
   today.
6. **Decouple "mass" claims from "annihilation/decay" claims.** Conjugate-pair
   bound state is a bound-state claim, not a positronium claim. State this in
   the creed.

---

## 4. Definitional risks (D1–D4) + D5

**D1 (vacuum reference for δ) — LOW RISK, partly vacuous.** Classically,
empty-box vacuum has exactly zero energy, so N-δE ≡ ∫T₀₀ minus escaped
radiation — identical to the bulk estimator the creed demotes. A2 gains
content over the old definition *only* via (i) the process-ledger reading
(v73 P2, real and measured) and (ii) the entanglement subtraction applied to
N-ENT (where it is load-bearing, not here). One real subtlety: with gauge on,
the charged object's Coulomb tail is part of its δ; the reference cannot be
the empty box. Include the tail inside the box; ledger the absorbing-BC loss.
Fable said this; concur.

**D2 (locality quantum) — HIGH RISK.** Fable strikes the lattice cell and
freezes h_eff = 2πQ. I go further: **both candidates are misnamed.** h_eff =
2πQ is a *process action quantum* (internal-clock adiabatic invariant), not a
locality quantum. The lattice causal cell is a *regulator*. Neither is a
"quantum of locality" in any QFT sense — QFT has no countable locality unit;
locality is microcausality, a property, not a quantum. The creed's "locality
quantum" vocabulary is borrowed from neither v76 nor QFT and does no work;
drop it. Rename N-LQ → N-AQ (action quantum); freeze N-AQ = 2πQ; accept that
M = N-ENT/N-AQ is *entanglement per unit process action*, structurally a
cousin of E/Q, and the creed must check it is not rediscovering E/Q with
extra steps.

**D3 (entanglement measure) — HIGHEST RISK.** Fable's three sub-risks (no
ensemble, degeneracy, area-law) are all correct. I add a sharpening: **the
QFI-on-ρ_Q route should be struck pre-audit, not audited in parallel.** For a
coherent-state-like U(1) object, QFI with respect to the charge generator
scales ∝ Q; with N-LQ = 2πQ, N-M = QFI/2πQ ≈ universal constant. Every U(1)
ball gets the same mass under this measure. This is not a candidate to audit;
it is a proof that the candidate is dead for the project's own object class.
Strike it now. Remaining live candidates: Gaussian-fluctuation EE (has QFT
counterpart; δ-subtract; symplectic spectrum of linearized fluctuations
about the parked soliton) and component-bipartition MI (cheap, classical,
analogical). Pick the first as primary; report the second as cross-check only.

**D4 (light mass-form seed) — MEDIUM RISK, with one gift and one creed-level
gap.** Gift: conjugate ball (Φ*, ω→−ω) is an exact equal-mass opposite-charge
partner; the bound-state experiment needs no new physics. Gap: conjugation
preserves the spectrum, so it cannot produce a *light* partner. Stage 4 needs
m_light ≪ m_nucleus and the creed has no light scale. This is not a D4 detail;
it is a creed-level absence. Flag it in the creed and block Stage 4 on it.

**D5 (units of M) — NEW, proposed by this seat.** Entanglement is
dimensionless; action has units; their ratio is [action]⁻¹, not mass. The
creed must exhibit the conversion constant (and its ontology, derived from
A1–A6) or concede N-M is dimensionless-relational and that "mass" in A3 is
metaphorical. This blocks all numerics until closed, and it is *independent*
of the D3 ensemble question — D5 would remain even if D3 were solved
perfectly.

---

## 5. Ranked next three actions (design only, no code)

**1. Close D5 on paper — the units audit.** Exhibit, from A1–A6 alone (no
imported constants), a quantity with mass units that is a function of N-ENT
and N-LQ. If this cannot be done, restate A3 as a dimensionless relational
conjecture (or Fable's order-parameter form) and proceed on A2 + A5.
Everything else is premature — there is no point building an N-ENT pipeline
for a quantity that cannot be mass by its own algebra.

**2. Strike the dead D3 candidate (QFI-on-ρ_Q) and freeze the survivor.** One
page, paper only: for each remaining candidate (Gaussian-fluctuation EE,
component MI), derive expected scaling with Q, R, ω, η along the measured
branch and the ball/ring degenerate-pair control. Strike any candidate that
is provably constant, ∝ δE·R (first-law trivial), or cut-dominated. Freeze
one primary + one cross-check with stated ensemble and units. This is Fable's
action #1 minus the dead candidate, sharpened.

**3. Pre-register the ball-vs-ring discrimination test and the conjugate-pair
bound-state design (renamed).** (a) Write the predicted N-M ordering under
creed vs E=mc² vs first-law-trivial S ∝ R·δE for the ball/ring pair *before*
measurement. This pair is cleaner than Fable's triple because δE is degenerate
by construction (~3% K-isomer split) and spatial correlation structure is
maximally different — it is the minimal configuration where A2 and A3 can
decouple. (b) Design (config-level) the conjugate-pair bound-state experiment
with closed-Gauss and closed-process-ledger criteria, and *rename* it — drop
"positronium" until a quantized annihilation channel is demonstrated. Reserve
inertia (N-INV) as the standing referee (Fable's action #3a): a creed-M that
fails to track inertia is renamed, not defended.

---

## 6. Kill gates

**G1 — Units gate (after action #1).** If D5 cannot be closed without an
imported constant, A3 is retired from "definition of mass" to "dimensionless
relational diagnostic" (or Fable's order-parameter form). v84 proceeds on A2
+ A5. No N-ENT pipeline is built for a struck A3.

**G2 — Redundancy gate (after action #3a, measured).** If N-M is proportional
to δE across the ball/ring pair AND a (Q, ω, η) scan — within two-lattice-
spacing scheme error — then A3 is empirically equivalent to A2 and should be
retired as redundant. Record the equivalence as a result (the creed's monism
is real but trivial), not a failure.

**G3 — Inertia gate (standing, reinstated as primary).** If N-M and inertial
response disagree beyond lattice error where both are measurable, N-M is not
mass. Rename it (entanglement inventory per action quantum), restore
inertia/δE as the mass definition (the QFT-backed position), and keep the
inventory only if it predicts something *else* (stability, branch ends,
binding).

**G4 — Stage-3 carrier gate.** If the conjugate-pair bound-state experiment
shows no bound state at any separation inside the validated medium-Coulomb
regime (with closed ledgers — not a protocol confound; apply the Fable R2
checklist), the blocker is the binding sector, not the mass ontology.
Escalate to a carrier/sector design discussion with explicit kernel-
authorization scope. No creed refinement addresses a missing binding sector.

**G5 — A1 testability gate (long-term).** If no kernel back-reaction path
(field → geometry) is authorized within the next two version cycles, A1 (one
continuum) is formally demoted from research basis to long-term aspiration.
Monism cannot be tested on a fixed-background simulator and the creed should
not let it drive design.

---

## 7. Stop-doing

1. **Stop using "mass" for a quantity in units of [action]⁻¹.** Either exhibit
   the conversion constant from A1–A6, or call N-M what it is (a dimensionless
   relational number, or an entanglement inventory). The word "mass" is doing
   rhetorical work the algebra does not support.
2. **Stop accepting A1 as a research basis on a fixed-lattice kernel.** A1 is
   a target ontology. The simulator is dualist by architecture (v76 §3.1).
   Until a back-reaction path exists, A1 is an aspiration, not a premise.
3. **Stop calling the conjugate-pair bound state "positronium."** It lacks the
   quantized annihilation channel that defines positronium in QFT. The creed's
   own §4.1 honesty standard forbids the label.
4. **Stop carrying the lattice-causal-cell option in D2.** It is a regulator;
   Fable said this, concur. Further: stop calling the survivor a "locality
   quantum"; rename to action quantum. The "locality" vocabulary has no QFT
   referent and does no work.
5. **Stop carrying the QFI-on-ρ_Q candidate in D3.** It is provably degenerate
   (N-M = const) for the project's own object class. Strike it now, not after
   an audit.
6. **Stop demoting N-INV (inertia) to optional.** It is the one measurement in
   mass units; it is the referee. Fable's point, strongly concurred.
7. **Stop scheduling Stage 4 atom experiments** until a light-scale mechanism
   exists. The creed does not produce m_light ≪ m_nucleus and Stage 4 silently
   inherits the gap.
8. **Stop re-opening v80–v83 failed methods** under creed vocabulary (creed
   §4.2). A renamed bag is still a bag.

---

*Kimi K3, Round 1. The creed's instinct — that this project's monism was
drifting into relabeled dualism, and that the lock-as-mass habit was an
imported label — is correct and valuable. The creed's central axiom (A3) is
not yet a definition of mass: it has a units defect (D5), a referent defect
(D3), and a QFT conflict (mass and entanglement are orthogonal in QFT for
elementary particles; redundant for composites). The honest path is: treat A2
(δ-field energy, QFT-anchored, Casimir-like) and A5 (process closure, already
measured) as the research basis; treat A3 as a conjecture to be falsified or
rescued by the D5 units audit and the ball/ring discrimination test. If A3
survives both, it earns its place. If not, v84 is still a real program — just
not a new mass definition.*
