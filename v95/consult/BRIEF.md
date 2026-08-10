# v95 developmental brief — the J-dependent link metric (space-direction transport partner of the door pair)

You are the **developmental session** for SCP v95. This brief is **self-
contained** but assumes the v94 record. You have read access (`--add-dir`).

## Context (from v94, 6 resume cycles)

v94 established two physics domains:

**CONVERSION (mode direction) — WON.** The door pair (pass P + settlement)
is a coherent, norm-preserving rotation between the field ψ_e and a register
χ, settling into matter ψ_m in whole action atoms at the beat. Measured:
spatial information preserved in conversion (DS-G3: r(Em,Ee) = +0.684 with
the pair vs −0.259 legacy); anti-runaway confinement (M1: door pair Em stays
at 17 vs legacy 417); metabolic park approaching (M6: chi_load=1 creates
10× slope reduction + bi-directional convergence from different initial
conditions); conservation exact (drift at FP floor); anti-ignition structural.

**TRANSPORT (space direction) — LOST (with a named confound).** The door
pair CONFINES matter (M1: no propagation). The face-2 stimulated outtake
ERODES coherent structure (M4: winding R2=0.24 vs baseline 0.62; M5: blob
anchored, no translation-by-conversion). The confound: the settlement
samples a rotating quantity (the cross-term C) at the beat, which is
stroboscopic aliasing — `door_avg` (beat-integrated cross-term) is the
named fix, partially tested (M4′: R2 improved 0.24→0.32).

**Claude's architectural synthesis (v94 face-3 design, resume 6):** the
registered architecture for v95 is **pair for the mode direction (conversion),
J-metric for the space direction (transport)**. The J-dependent link metric
is the space-direction coherent transport partner that the door pair cannot
provide alone.

## The J-dependent link metric — what it is

The v93 design (§II.5) named it: the link metric (the weight `w_ij` on each
spatial link) should respond to the link current `J = 2·Im(ψ_i* · ψ_j)`.
Today the metric is geometry-only (`w_ij = (sA/Aref)·(dref/sd)`). The
proposal: make the metric **dynamical** — `w_ij` evolves toward the flux
pattern. This is "space curves because energy is conserved and space is
convertible" — the programme's founding line, implemented as the metric
responding to the current.

This is RECKONING's recommended next face, claude's promoted v95 candidate,
and the one structural piece the v94 record says is missing for coherent
spatial transport.

## What v95 must produce

### 1. Theory — derive the J-metric from conservation + the substrate

- State the postulate: the link metric responds to the link current. Derive
  the form (a relaxation? a variational principle? a Hodge-star dynamics?).
- Show that this makes spatial transport **phase-structure-aware**: coherent
  structures (where J has a definite sign) bind differently from incoherent
  ones (where J averages to zero).
- Address the binding/mobility antagonism (RING_DNLS B.3): does the J-metric
  release it? A self-trapped blob has a definite J pattern; if the metric
  responds to J, the trap might dissolve for coherent motion while persisting
  for static structures.
- Address the phase-retention wall (M4): the J-metric makes the spatial
  channel phase-sensitive. Does this help the door pair's phase cross?

### 2. Fully derived symbolic equations

- The metric dynamics: `dw_ij/dt = ...` as a function of `J_ij`.
- The conservation proof: does the dynamical metric preserve the conserved
  sum? (The geometric metric is fixed; a dynamical metric changes the hop
  angles. Show that pairwise norm-preservation still holds.)
- The dispersion: with a J-dependent metric, what is the group velocity?
  Does the antagonism (trap kills mobility) dissolve?

### 3. Numerical simulation method (the kernel face)

- A concrete update primitive behind a new byte-inert knob (e.g. `j_metric`,
  default 0). Specify:
  - the state (a per-link weight accumulator? a low-pass on J?),
  - the operations (how w_ij evolves; how the hop angle uses the new w_ij),
  - the schedule (where in the pass order),
  - byte-inertness at 0.
- The C/Go mirror parity requirement.

### 4. Acceptance bars (the ratchet)

- **J1 (phase-sensitive binding):** a coherent blob (seed_mw=2 winding) with
  j_metric=1 retains R2 ≥ 0.5 (vs 0.24 at j_metric=0).
- **J2 (antagonism release):** a self-trapped blob with kx tilt and
  j_metric=1 translates (COM advance ≥ 5 cells) while retaining Em_max ≥ 50%.
- **J3 (conservation):** drift at FP floor.
- **J4 (anti-ignition):** bath stats unchanged at j_metric=1.

### 5. The composition with the door pair

The v95 architecture: pair (mode direction) + J-metric (space direction).
Show how they compose:
- The door pair provides coherent conversion (field ↔ matter).
- The J-metric provides coherent spatial transport (matter moves coherently).
- Together: a metabolizing, translating, phase-retaining object.

## Constraints (permanent)

Same as v94: no background, no imported species, conservation as theorem,
door stays special, no quark road, no sum-then-square, ratchet governs.

## Output

`v95/consult/DESIGN_claude.md` — theory, equations, numerical method, bars.
Print your session ID at the head. Flag every load-bearing assumption.
