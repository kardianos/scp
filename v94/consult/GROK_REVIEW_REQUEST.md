# v94 — grok review request: door-pair implementation vs design

You are reviewing the **v94 implementation** of claude's "door pair" design.
This brief is **self-contained**. You have read access to the v94 directory
(`--add-dir /home/d/code/scp/v94`); produce a **review** of the
implementation against the design and the programme constraints. Output to
`v94/consult/REVIEW_grok.md` (the orchestrator pre-creates the file path; you
just write to it).

## What v94 is

v94 is the next-theory campaign after v93. The user's three observations:
1. EMF translates at c in vacuum.
2. EMF translates slower through matter.
3. EMF alternates between two expressions (E and B), as if two alternating
   incomplete conversions — a parallel the user drew to v93's matter
   intake/outtake attempt.

The user's hypothesis: **if we can find/construct/discover a primitive
translation device for EMF, the v93 matter intake/outtake becomes adoptable.**

Claude's developmental design (`v94/consult/DESIGN_claude.md`) derives the
primitive (pass F already does it: local precession + pairwise Givens = the
two interleaved first-order conversions = the Maxwell shape on the roster),
then mirrors it at the door as a third link type — a link in the **mode**
direction. The door pair is a continuous field↔register Givens rotation
(pass P, Strang-symmetric around pass F), settling into dense matter only in
whole action atoms at the beat.

The opencode orchestrator implemented face 1 of claude's design (the forward
intake path; the stimulated outtake is deferred to face 2). The implementation
status is `v94/IMPLEMENTATION_STATUS.md`.

## The programme constraints (permanent; v94 §IV / v93 §IV)

- **No background.** Nothing persists and is merely re-valued. The roster is
  dynamical; no fixed lattice; no fixed fiber with a separate ledger/identity.
- **Energy is never destroyed — the only law.** Conservation must be a theorem
  of the update (pairwise norm-preservation), never a patched ledger.
- **No imported field or species.** No Maxwell import, no SI constants. The two
  "expressions" of EMF must be derived from the fabric.
- **The door stays special.** Quantized action atoms `ε = A₀ω/2π`, work-function
  gate, two-atom credit — these are the "atoms at mode boundaries." Do not
  dissolve the click grain into a continuous rotation.
- **No quark road.** Two state-halves with separate physical jobs is the
  rejected pre-v89 lattice.
- **No sum-then-square.** Never `|ΣA|²`; pairwise hops only.
- **No feedback / shadow / bias correction.** v92 closed three forms.
- **No phase written to the cell clock `th2` at the door.** Carry slot-borne.
- **The ratchet.** Every kernel change runs the full battery before commit;
  bars sharpen by measurement, never soften to pass.

## What was implemented (in this order, all in `v94/`)

1. **State + knob + drift meter.** New persistent per-cell register
   `dc1, dc2` (the polarization χ = dc1 + i·dc2). New knob `door_g`
   (default 0). `total_energy` now includes `|χ|²` per cell (so the drift
   meter stays honest when the register carries energy). Mirrored in C
   (`kernel/freecell.c`) and Go (`fab/sim.go`, `fab/params.go`,
   `fab/meters.go`, `fab/run.go`).
2. **Pass P (the door pair), Strang-symmetric around pass F.** Two half-steps
   bracketing pass F. Each half-step, per cell:
   - register precession by `w2e·dt/2` (matches the dense-clock direction);
   - fiber Givens by `θ = 0.5·door_g·w2·dt` (note: **bare w2, NOT load-detuned
     w2e** — claude §2.4: vacuum detunes away from its own door, bound matter
     is pulled toward resonance; this re-derives the workfn law as detuning
     physics).
   The Givens uses claude's "Maxwell-shape cross pairing":
   `fa1' = cc·fa1 + ss·dc2; fa2' = cc·fa2 − ss·dc1; dc1' = cc·dc1 + ss·fa2; dc2' = cc·dc2 − ss·fa1`
   (verified pairwise norm-preserving for `|ψ_e|² + |χ|²`).
3. **Pass-6 settlement (the click).** When `door_g > 0`, the legacy
   `if (Ee > econd_i)` condensation branch is replaced by `if (chi2 > econd_i)`
   (the workfn text is reused verbatim, applied to `|χ|²` instead of `Ee`).
   The transfer is a **solved Givens rotation** between χ and ψ_m by the angle
   θ that moves exactly d1 (the fired energy). θ_max found by 256-sample scan
   on [0, π/2] (robust to all sign cases of A≡|χ|²−Em and C≡2·Im(ψ_m*·χ));
   then 60-iteration bisection on [0, θ_max] for ΔEm(θ)=d1r. The Givens
   convention is claude's settlement convention `m1' = cc·m1 − ss·dc2;
   m2' = cc·m2 + ss·dc1; dc1' = cc·dc1 − ss·m2; dc2' = cc·dc2 + ss·m1`
   (verified `|ψ_m|² + |χ|²` invariant). The space pull `dsp = s_pull · d1`
   is debited from Es exactly as in the legacy branch; the field fa1/fa2/Ee
   are NOT touched (the field→register transfer already happened in pass P).
4. **Diag.** When `door_g > 0`, the RESULT block gains a v94 line printing
   `chi2_tot` and `door_in_total`.

## What the implementation measured (post-implementation)

| bar | result |
|---|---|
| **Invariant surface** (87 bars at `door_g=0`) | **ALL GREEN 87/87** (battery exit 0) — byte-inertness across all changes |
| **M0 C≡Go parity at strength** (`door_g=0.3, amp_nat=2`) | **byte-identical** C and Go (exceeds claude's tolerance-abx expectation; closes the v93 Go-only debt for the door pair) |
| **M2 conservation** (drift ≤ 1e-13 on live matter) | **GREEN** — T=500 with `noise_amp=1.5` (live field): `drift_rel = -3.345e-14` ≪ 1e-13 |
| **M3a anti-ignition (dark bath)** | **GREEN (structural)** — dark bath byte-identical between `door_g=0` and `0.3`; pass P has no field to couple to on a dark bath, so `chi2≡0` and settlement never fires |
| **Functional verification** | **GREEN** — register fills continuously (`chi2_tot` 91.8 at `noise=1.0`, 176 at `noise=2.0`); settlement fires in whole atoms once per-cell chi exceeds the grain (`door_in` 0 → 0.30 → 1.51 across noise 1.0 → 1.5 → 2.0) |

## What we want from you (the review)

A **physicist + reviewer's read** of the implementation, in this order:

### 1. Theory compliance — does the implementation match claude's design?

Read `v94/consult/DESIGN_claude.md` §III.1–III.5 carefully. Read the actual
kernel diff (the orchestrator's summary is in `IMPLEMENTATION_STATUS.md`).
Verify:
- The pass-P Givens stencil matches claude's §5.2 Maxwell-shape cross pairing.
- The settlement Givens convention matches claude's §2.5 (the ΔEm formula and
  the real-component stencil).
- The coupling scale is `g = door_g · w2` (bare-pitch, NOT load-detuned —
  claude §2.4's load-bearing choice that re-derives the workfn law).
- The pass schedule is Strang-symmetric around pass F (claude §5.2).
- The settlement preserves the click grain (atoms_fire / atoms_clamp verbatim;
  only fired at the beat).

### 2. Conservation proof — does `Σ(Em + Ee + Es + |χ|² + flight)` truly telescope?

Verify by direct algebra that:
- Pass P (per-cell precession + fiber Givens) preserves `|ψ_e|² + |χ|²` exactly.
- The settlement Givens preserves `|ψ_m|² + |χ|²` exactly.
- The space pull + ledger bookkeeping is conservation-neutral (as in legacy).
Identify any path by which the implementation could leak energy that the drift
meter would not catch (the drift meter includes `|χ|²`, but does it include
every chi-touching operation?).

### 3. Anti-ignition — does the structural argument (claude §4.3) hold?

Claude argues four layers (detuning, phase averaging, settlement gates,
zero-throughput coupling). The implementation has only the forward (intake)
half. Is the structural anti-ignition argument still sound for face 1? What
new risk does face 2 (stimulated outtake) introduce, and how should it be
gated?

### 4. The four open flags (claude §III.8 / `IMPLEMENTATION_STATUS.md`)

The implementation touches four programme-rule edges. Give your judgment on
each:
- **(a) Two new persistent floats per cell.** Is χ a fourth energy mode (which
  would need user signoff as a major face), or is it legitimately mode-transit
  bookkeeping (cousin of the link flight parcel `slem`, which is also a per-
  link pair of floats)?
- **(b) The settlement writes `th2 = atan2(m2n, m1n)` after the Givens.** The
  QUENCH-3 rule (v93 §IV.6) says "do not write th2 at the door; phase must be
  slot-borne." Does the settlement's write (the arg-readout of a norm-
  preserving Givens, with the phase sourced from χ not from the field) comply
  in spirit, or does it repeat the failure mode QUENCH-3 measured (R2d≈0.02)?
- **(c) Vacuum dispersion shift `O(g²/Δ²)`.** The DS optics bars live on the
  field sector. Does the implementation's `door_g=0` byte-inertness protect
  the DS bars adequately, or is there a non-byte-inert path?
- **(d) The stimulated outtake is deferred to face 2.** Is face 1 (intake
  alone) a valid test of the design's central claim (the metabolic cycle
  propagates), or does it need both directions to be meaningful?

### 5. What's missing or wrong?

Identify any place where the implementation deviates from claude's design, any
bug, any optimisation that loses physics, any constraint violation. Be
specific (file + line if possible).

### 6. Verdict

For each of the four flags: ACCEPTABLE / NEEDS-WORK / BLOCKING. For the
implementation as a whole: SOUND / SOUND-WITH-FIXES / REWORK.

## Anchor files (read these in this order)

1. `v94/consult/DESIGN_claude.md` — the design (the spec).
2. `v94/IMPLEMENTATION_STATUS.md` — the implementation summary (this brief's
   long form).
3. `v94/kernel/freecell.c` — the C kernel of record. Key sites:
   - knob: `door_g` decl + default + parser.
   - state: `dc1_, dc2_` decl + alloc.
   - drift meter: `total_energy()` (includes `|χ|²`).
   - **pass P** (the door pair, two halves): search `pass P (v94 door pair`.
   - **settlement**: search `v94 door settlement` (inside pass 6).
4. `v94/fab/step.go` — the Go campaign-kernel mirror. Same sites, same labels.
5. `v94/fab/sim.go` — Sim struct (dc1, dc2 fields), allocAll, totalEnergy.
6. `v94/fab/params.go` — DoorG field + Defaults + SetKV.
7. `v94/fab/run.go` — law header + RESULT block.
8. `v94/laws/laws_V3a.cfg` — the adopted law table (unchanged; door_g=0 is
   the default).
9. `v94/battery.log` — the post-implementation battery run (ALL GREEN 87).

## Reference (do NOT re-litigate)

- v93 closed (frozen): `v93/theory/FINAL_STATUS_ADOPTION.md` — the substrate
  v94 stands on. The v93 wall is the open problem; v94 does not relitigate
  v93's measured closures.
- The QUENCH-3 / RING_DNLS / Face C record is the measured evidence for the
  wall; v94 engages it, not re-asks it.

Be specific, be honest, flag every load-bearing assumption. The orchestrator
will fix any NEEDS-WORK / BLOCKING items before opening the M1/M4/M5/M6
campaign arms.
