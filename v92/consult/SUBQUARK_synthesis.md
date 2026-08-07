# SUBQUARK design consultation — synthesis of three reviews (2026-08-07)

User question: can the dense cell be broken into discrete sub-components
(individually incomplete, joined into one cell) and still pass the battery?
Minimally two components; interpreted as uptake/outtake defining movement;
decay rates (translation without stable uptake/outtake modes); and three
geometries — projective line, toroidal, discrete differential geometry.

Three independent reviews solicited (the question text is in
`SUBQUARK_question.txt`): grok CLI (`SUBQUARK_grok.md`), claude CLI fable
effort max (`SUBQUARK_claude_fable.md`), and opencode's own (below, woven
into the synthesis). All three converge.

## The convergence (6 points of agreement)

1. **The state layout is not the point.** Two floats `(m1,m2)` vs `(Em,th2)`
   is a relabeling — the kernel already stores the polar form and a per-slot
   complex shadow. Complexifying the *state* changes nothing.
2. **What must change is the transport algebra:** dense within-mode transport
   becomes **unitary pairwise rotations** (pass F's cousin), replacing the
   additive magnitude want. That is the genuinely new physics.
3. **Conservation *dissolves*, not solves.** Pairwise Givens hops never
   sum-then-square, so `|ΣA|² ≠ Σ|A|²` never arises — each hop conserves the
   two-cell norm exactly (why the field is unitary to roundoff). The cross
   term `2·Re(ψᵢ*wψⱼ)` that broke additive bookkeeping IS the link current —
   it is where momentum lives.
4. **Projective line = pure notation** (RP¹ ≅ S¹ = the phase circle already
   stored). Drop it. (CP¹/Bloch would be new but needs 4 components — bigger.)
5. **Torus is real physics but at the parcel/ring level**, not the cell. Use
   it as the *acceptance language* (winding retention, phase-slip rates at
   doors, ring circulation quantization).
6. **Uptake/outtake is the interpretation of the current's moments and of the
   doors — NOT two state components.** Monopole of J = metabolism; dipole =
   momentum; intake/outtake = the two mode boundaries. Attaching physical
   jobs to m1 vs m2 samples a rotating frame (parametric pump, ignition-
   shaped) and is the first step down the rejected quark road.

## The one synthesis insight (claude, confirmed by all three)

> The three feedback forms did not fail by bad luck. They failed because the
> additive Em-ledger rejects the interference term.

Keeping additive Em-conservation and bolting phase alongside left the cross
term nowhere to live — so it surfaced as throughput amplification (xsec
7.76→7.18), got renormed away (zero-sum), and that renorm degraded the very
signal wanted (e3b |cos| 0.28→0.08). "Conservation must become a theorem of
the update (norm preservation), not a ledger patched per form." This is why
form-3 done right still failed — and why unitary hops will not.

## The decay-rate diagnosis, made exact (claude)

"Matter can't translate" = **Γ_p/Γ_E → ∞**. Γ_E (structure loss) small
(blobs at x̂*=0.62, burn down over hundreds); Γ_p (momentum loss) ≈ infinite
(B2 ceiling ~5e-3, P2 ~100× deficit). Γ_p ≈ ν_door·(1−r); QUENCH-3 measured
r ≈ 0 (door writes phase at 3.5σ, nothing retains) ⇒ momentum dies within a
few door cycles. **Pre-register: hotter objects lose momentum faster unless
r→1; P2 deficit closes as r does.**

## The resulting design direction (the next memo)

- **Build:** complex dense parcel + DEC-form currents — J a 1-cochain on
  slots, cells update by δJ, doors as source terms, **unitary pairwise
  rotations** as the conservation mechanism (replacing the additive want).
  Hodge star already in the kernel (field hop normalization); a dynamical
  star = FLOW's bed law = the implementable "space curves."
- **Acceptance:** torus/winding — slip rates at doors, winding retention,
  ring circulation quantized in 2π/N. L2 = maintenance of a live fifth
  (sustain ≠ ignite); nucleation a separate QUENCH-style experiment.
- **Interpretation:** uptake/outtake = J's moments + the doors. Open-books
  gravity (spatially separated intake/exhaust along a chain) = a later face,
  after translation exists.
- **Retention (the risky piece):** the fired atom must carry arg(ψ) —
  departure phase from the composed amplitude, NOT m·th2 (grok criterion 1;
  the field's phase-faithful clicks prove the pattern can work).
- **Anti-ignition (honest risk):** unitarity kills energy-runaway outright;
  coherence-runaway remains (bath ρ_coh 0.77). Mitigation: rotations
  EXCHANGE amplitude reversibly (no Kuramoto pull fixed-point — kappa_lock
  and the 4 ignition machines are all dissipative pulls); wall at the door;
  bath-armed staged bars. Bar shift: "nothing condenses" → "bath stats
  unchanged at inert, creation rate matches the measured quench law when live."
- **Engineering:** knob-gated parallel lane (amp_native), scalar want
  untouched at 0 (byte-inertness vs the 87 bars); tolerance-based abx
  (compounded rotations amplify the 1-ulp C/Go divergence); **comb gate
  derive-then-retire** (it may be derivable from composition — the deepest
  vindication — but earn it, don't assume it).

## Net

The "break the cell into two components" intuition was right that the dense
sector needs the field's internal structure — but the three consultations
converge that the components are modes of ONE amplitude whose transport is
unitary, and "uptake/outtake" names the current and the doors, not the two
halves. That reframing turns the failed feedback programme into a buildable
one: the conservation wall dissolves because unitary hops never sum-then-
square, and the cross term that the additive ledger rejected is exactly the
momentum current that was missing.
