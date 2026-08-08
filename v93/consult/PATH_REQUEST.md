# v93 — research-path consultation request (post faces A/B/C)

You are reviewing the SCP programme. This brief is self-contained. Produce a
**staged research path** (short-term concrete experiments + long-term arc)
toward the goals in §4. Be honest about what is MEASURED vs SPECULATIVE, and
respect the no-background constraint (§1). Do not write code; reason and plan.

## 1. The programme (one law, four constraints)

SCP: **energy is never destroyed; it only changes mode; space is one of the
modes.** Therefore: no background (nothing persists and is merely re-valued);
no imported field or species; everything emergent from the fabric. Matter is
converted space; motion is regular conversion, not displacement. There are no
equations of motion in the usual sense (nothing has a trajectory).

## 2. The substrate (v89+ free-cell, kernel freecell.c)

A dynamical cell roster (cells born/die/jam/repack — NO fixed lattice) with
three coupled sectors per cell:
- **field** fa1,fa2 (complex amplitude; Ee=|ψ|²) — UNITARY transport (product
  of pairwise plane rotations, pass F).
- **dense matter** Em (scalar) + clock th2 (phase).
- **space** Es (convertible; radii r∝∛Es).
Plus a **conversion door** (field↔matter): field condenses at matter on a beat
above a work-function threshold; matter evaporates above cap; the click is
quantized into action atoms ε=A₀ω/2π. Radiance (k_rad) taxes over-fixed-point
matter back to field (fixed point x̂*=0.62). Conserved: Em+Ee+Es.

## 3. v93 — the unitary dense channel (the current open programme)

Give the dense sector the field's transport algebra: within-mode dense
transport becomes UNITARY PAIRWISE PLANE ROTATIONS (pass U, a cousin of pass
F) on ψ_m=√Em·e^{i·th2}, replacing the additive magnitude "want." Each Givens
hop conserves the two-cell norm exactly ⇒ ΣEm conserved to roundoff
(conservation is a theorem of the update, not a patched ledger). The cross-
term 2·Im(ψ_i*ψ_j) is the link current J = dense momentum. Knob amp_nat
(default 0 = byte-inert). The door stays discrete/quantized. Status: L1-A
(coherent translation) RESCUED with caveats [linearize amp_logate>0 + measure
the dense-hop current fd, not the misleading tagged centroid]; L1-B
(conservation) GREEN on live matter; L1-C (anti-ignition) UNRUN on the armed
coherent bath.

### Faces A/B/C (just executed, 2026-08-08)
- **A — hop schedule:** the sequential sweep is a C6-breaking Trotter
  artifact. New knob hop_order (0=seq, 1=Strang symmetric, 2=randomized).
  LIVE vacuum ring holds winding W=+2 to t=1000 under Strang (seq degrades).
  e3b dt-spread 4%→1.3%→0.3%. Strang=best retention; random=best dt-invariance.
- **B — arg(ψ) door refill:** REFUTED as a retention rescue. Door traffic with
  a random-phase field scrambles the winding regardless of coherent handling.
- **C — condensation lane (THE LEAD):** the unitary channel + the law's load-
  detuning w2e=w2/(1+q_detune·x) is a real DNLS self-trap nonlinearity. A
  spread packet SPONTANEOUSLY CONDENSES into long-lived dense hoards (stable
  to t=300, NO lock/gate/door; the additive law cannot do this). Deep corner
  (q_detune=3.6, amp=2): Em_max~7 (3×cap), frozen envelope. 189 hoards in a
  long-tailed (power-law-like) hierarchy 0.5→7.2. Radiance (V3a) SELECTS sizes:
  truncates the spectrum at the fixed point Em*≈1.55; deep self-traps resist.
  This is creation-adjacent (echoes QUENCH: "dynamics created what gates could
  not"). Image sequence shows no visible dissipation — possibly a meta-stable
  particulate. CAVEAT: it's a stable condensed PACKET (PR≈250 cells, broad),
  not a tight soliton; it binds ENERGY not PHASE; mobility is antagonistic
  (the corner that binds does not translate).

## 4. The goal chain (the user's question)

If Face C is a meta-stable particulate, the programme wants to climb a scale
stack toward real matter:
1. **Characterize / refine the particulate** — is it a genuine stable soliton
   or a transient condensate? Can it be made tight (low PR), monodisperse,
   phase-coherent? Does it survive the bath / the door / radiance long-term?
2. **Scale correspondence to reality** — what physical scale/mass does this
   fabric particulate map to? How would one establish a correspondence (a
   gauge) between the simulation's dimensionless units and SI mass/length/
   time, WITHOUT importing real-world constants (constraint §1)? What is the
   in-fabric observable that sets the scale?
3. **Electron shells / Pauli exclusion** — can a bound particulate (nucleus-
   analog) hold a cloud of light opposite-charge field/matter in discrete
   shells that obey an exclusion-like rule? (The programme's prior stage-3
   goal; v75 three-fabric proposal exists but was never built.)
4. **Quark / proton / neutron behavior** — sub-structure of the particulate
   (component flavors), confinement-like binding, the stable charged ball +
   neutral companion. (NOTE: pre-v89 "quark" multi-component lattice fields
   on a fixed background were EXPLICITLY REJECTED — any revival must be mode
   structure that dies with the cell, not a fixed fiber with a separate
   ledger/identity. v93 §IV.2/IV.7.)
5. **Atomic binding → discernible elements** — nucleus + shells, binding
   energies, a discrete spectrum of stable configs identifiable as H/He/...

## 5. What we need from you

A staged path:
- **Short-term (next ~5–10 concrete experiments):** what to run to (a) decide
  whether Face C is a real particulate (lifetimes, bath survival, tightness,
  monodispersity, phase coherence); (b) establish the scale correspondence.
  For each: the apparatus, the bar, what outcome decides the question.
- **Long-term arc:** the credible path from a confirmed particulate through
  shells/exclusion → nuclear sub-structure → atomic binding → elements.
  Mark the load-bearing assumptions and the most likely failure points. Where
  does the no-background constraint bite hardest? Where might a kernel change
  be unavoidable (and what would it have to preserve)?
- **Recommendations:** the 2–3 highest-value next moves; what to NOT do yet.

Ground everything in what is measured (§3). Flag speculation explicitly.
