# v94 — developmental design brief: the EMF translation primitive mirrored to matter intake/outtake

You are the **developmental session** for SCP v94. This brief is **self-
contained** — read top to bottom. You have read access to the v94 directory
(`--add-dir`); the kernel is `v94/kernel/freecell.c` (modifiable in v94 only;
**do not modify v93/**). Produce a **developmental design** (theory,
derivations, math, mechanics, numerical simulation method) — not code. Output
to `v94/consult/DESIGN_claude.md`.

The user has authorized: **"full authority to write and manage the numeric
kernel in v94."** You design; the opencode orchestrator will implement your
design in the kernel, then grok will review the implementation, then the runs
will be shepherded to completion.

---

# PART I — the programme

## I.1 One law, four constraints

SCP rests on one principle: **energy is never destroyed; it only changes mode;
space is one of the modes.** Therefore (v89/PRINCIPLE):

- **No background** — nothing persists and is merely re-valued. A dynamical
  cell roster (cells are born, die, jam, repack); no fixed lattice.
- **No imported field or species** — everything emergent from the fabric.
- **No equations of motion in the usual sense** — nothing has a trajectory;
  "motion" is regular conversion (mode→mode), not displacement.

Matter is converted space; space curves because energy is conserved and space
is convertible; the only law is the conservation of total energy across the
three modes.

## I.2 The substrate (v89+ free-cell, kernel `freecell.c`)

A dynamical cell roster carrying **three coupled sectors** per cell:

| sector | state | energy | role |
|---|---|---|---|
| **field** | complex amplitude `fa1, fa2` (i.e. `ψ_field`) | `Ee = fa1² + fa2² = \|ψ_field\|²` | the unitary, interfering, vacuum-speed sector |
| **dense matter** | scalar `Em` + clock phase `th2` | `Em` | the loadable, cap-defended, radiating sector |
| **space** | `Es` (convertible mode; cell radii `r ∝ ∛Es`) | `Es` | the curving footprint sector |

Plus a **conversion door** (pass 6 in the kernel):
- **intake** = field→matter condensation: at a clock beat, if `Ee` exceeds a
  work-function threshold, an action atom `ε = A₀ω/2π` clicks (A₀=`quant_A0`
  =1.15; ω = the source's load-detuned pitch). Two-atom credit; integrated-
  and-fire.
- **outtake** = matter→field radiance: graded sub-cap emission
  `O = k_rad · cap · x^{p_rad}` with `x = Em/cap` (laws_V3a: `k_rad=0.05,
  p_rad=4`); matter→field evaporation above cap.
- **space** gets / loses the footprint on each conversion (`s_k`, `s_disp`).

Conserved charge: `Em + Ee + Es`, witnessed at the FP floor (≤1e-13) by the
battery's drift bars.

## I.3 The two transport algebras (where the wall lives)

| | field sector | dense sector (v93) |
|---|---|---|
| state | `fa1, fa2` (complex amplitude) | `Em` (scalar) + `th2` (phase); v93 forms `ψ_m = √Em·e^{iθ₂}` |
| within-mode transport (the hop) | **pass F** (`kernel:2151`): **local precession** `(fa1,fa2) ← Rot(w₁e·dt)·(fa1,fa2)` **then** **pairwise Givens hops** of angle `τ = field_J·w·dt/√(s_i·s_j)` between link endpoints — exactly norm-preserving | **pass U** (`kernel:1591`, v93, knob `amp_nat`, default 0 = byte-inert): same structure — local precession of `ψ_m` by `w₂e·dt`, then pairwise Givens hops on `(m₁,m₂)` |
| conservation | unitary to roundoff (rotations) | unitary to roundoff when `amp_nat>0`; additive magnitude-want when `amp_nat=0` |
| mode-boundary (the door) | n/a (field is the far side of the door) | **pass 6 door**: quantized clicks; intake/outtake as above |

**Key fact:** pass U is *structurally identical* to pass F. The field sector
already propagates at the substrate's natural speed (its `c`) by exactly this
primitive: a local clock rotation interleaved with a pairwise Givens hop
(Strang-split). What pass U does NOT do is engage the **door** — the
field↔matter↔space mode-boundary — as a *translation* primitive. The door is a
discrete click today: a switch, not a wave.

## I.4 Where v93 ended (the open wall, frozen 2026-08-09)

Read `v93/theory/FINAL_STATUS_ADOPTION.md` for the full record. The close:

- **v93 thesis confirmed where structural:** pass U gives the dense sector the
  field's *within-mode* unitary transport. Conservation dissolves into a
  theorem (drift at FP floor). Anti-ignition is structural (random-phase bath
  self-gates).
- **But v93 hit three measured walls:**
  1. **Translation is narrow-band and seed-variant** (L1-A `cos→0.99` at one
     `amp_nat≈2` for one seed; free packet diffracts).
  2. **Binding and mobility are antagonists** (RING_DNLS B.3: the corner that
     self-traps does not translate).
  3. **The door does not carry phase** (QUENCH-3: spin retention `R2d≈0.02`
     across all arms; the door writes the *cell clock* `th2` from a decohering
     field — `arg(ψ_m)` slot-borne never delivered steady-state retention).
- **The metabolic intake/outtake never closed.** `rad_dual` (a spontaneous
  field→matter intake dual of radiance) showed real **I** on a quiet clump
  (SLOW_REFILL ~1.2–1.3× vs SLOW_DRAIN) but **no `Em*` park** — no restoring
  fixed point. The coherent cross-term `g` (Claude's qgraph §3.5) was never
  coded.
- **RECKONING:** Face C's "spontaneous condensate" is a *transient* that peaks
  ~t=300 and relaxes away. Not a particle.

**The honest one-line status:** the dense sector's *within-mode* transport is
solved; the dense sector's **mode-boundary transport (the metabolic cycle)** is
not. **This is v94's territory.**

---

# PART II — the user's three observations and the hypothesis

## II.1 The observations (verbatim from the user)

1. **EMF moves around "c" the speed of light when not moving through matter.**
2. **EMF moving through matter is slower.**
3. **EMF alternates between two expressions, as if there are naturally two
   alternating conversions going on, each incomplete.**

The user draws a parallel between (3) and the v93 matter intake/outtake
attempt. The user's stated idea:

> If we are able to find/construct/discover a **primitive translation device
> for EMF**, the v93 matter intake/outtake will become adoptable.

## II.2 Reading the observations against the substrate

In the SCP substrate, **the field sector already translates** by exactly the
primitive the user is pointing at. The two "expressions" of the field amplitude
are `fa1` and `fa2` (the real / imaginary parts of `ψ_field`). Their alternation
is the kernel's pass F:

```
local precession:   (fa1, fa2) ← Rot(w₁e·dt) · (fa1, fa2)     [in-time rotation]
pairwise hop:       (fa1,fa2)_i  ↔  (fa1,fa2)_j  via Givens τ  [in-space exchange]
```

These are **two interleaved, mutually-inducing operations**, each
*incomplete* by itself: the precession alone would just spin each cell in
place (oscillation, no translation); the hop alone (without phase evolution)
would just spread amplitudes symmetrically (diffusion, no directional wave).
The **interleaving** is the propagation — and the resulting dynamics is a
discrete Schrödinger / Maxwell-equivalent on the roster. Observation (1) is
*this* primitive in vacuum; observation (2) is the *load* `w₁e → w₁e/(1+q·x)`
slowing the precession in loaded matter (the refractive effect is in the
kernel already, `q_detune`, just for the dense clock); observation (3) is the
structure of pass F itself.

**This is the primitive.** It is **not** the local clock alone, **not** the
hop alone; it is the **mutually-induced pair**, where each conversion is
first-order and incomplete, and the pair produces propagation. Real Maxwell:

```
∂E/∂t =  +c² ∇×B       (E grows from the curl of B; first-order, incomplete)
∂B/∂t =  −   ∇×E       (B grows from the curl of E; first-order, incomplete)
```

is exactly this shape: two first-order equations, each **incomplete** (the curl
of the other field, with no source), and the pair is the wave equation at c.
Substitute E↔(fa1) and B↔(fa2) (or more carefully, the analytic-signal
decomposition) and the discrete pair above is the roster version of the same.

## II.3 The hypothesis, made precise

The user's claim is that this primitive — **two first-order, mutually-induced,
incomplete conversions interleaved** — is the *universal* translation device of
the fabric, and that the v93 wall at the door is precisely the absence of this
primitive at the **mode boundary**: today the door is a *single* one-way
switch (intake OR outtake, fired discretely), not a *mutually-inducing pair*.

**v94 asks:** can the door be rewritten so that **intake** (field→matter) and
**outtake** (matter→field) are **two first-order, mutually-inducing,
incomplete conversions** — exactly analogous to (∂E, ∂B) — such that:

- the **metabolic cycle** of matter propagates as a coherent wave (the matter
  analogue of an EM wave),
- **conservation** is a theorem of the pair (Stokes / Givens, as in v93 §II.4),
- the discrete **click grain** of the action atom survives (the pair lives
  *between* clicks; the click is the integrator-of-fire for the pair, not a
  replacement),
- the two expressions of matter that alternate are **not** m₁/m₂ (that is the
  within-mode amplitude v93 already used) but **the door's two directions** —
  or, if you can argue cleanly that the within-mode pair already *is* the
  primitive and the door should be its boundary condition, say so and derive
  the consequences.

The user explicitly invites either answer: derive the primitive for EMF, then
show how the matter intake/outtake adopts it — **whether by extending the
primitive to the door, or by showing the door should be a boundary condition
of the existing within-mode primitive, is itself part of the design question.**

---

# PART III — what to produce (the developmental design)

Produce `v94/consult/DESIGN_claude.md` with these sections. Be a physicist:
state postulates, derive, identify the math, design the numerical method,
pre-register bars. **All in fabric language** (no Maxwell import; derive
Maxwell-shape FROM the substrate, not the reverse).

## III.1 Theory (the EMF translation primitive, fabric-derived)

- Derive, from conservation + the three-sector substrate, the **two-
  expression pair** that makes the field sector propagate. Identify the two
  expressions precisely (fa1/fa2? Ee/B-analogue? a different decomposition?)
  and show how the kernel's pass F instantiates the pair.
- State the **refractive / load effect** (observation 2) in fabric terms —
  where does the slow-down live? Show it is already in the kernel (q_detune
  on the dense side; identify the field-side analogue, if any).
- Show that the pair is **mutually inducing** (each conversion sources the
  other's time-derivative) and **incomplete** (first-order; neither alone
  propagates).

## III.2 Mirror to matter (the central design)

- Identify the **matter-side analogue of the pair**. Two candidates to
  evaluate explicitly:
  - **(A) Within-mode extension:** the matter pair is already `(m₁,m₂)` under
    pass U (as in field pass F). Then the door should be a **boundary
    condition** of this primitive — state exactly how.
  - **(B) Mode-boundary extension:** the matter pair is **(intake, outtake)**
    at the door, mutually inducing. Cast the door as a coupled first-order
    pair whose mutual induction propagates the metabolic cycle across a chain
    of cells (the matter analogue of an EM wave). Identify what plays the role
    of `∂E/∂t = c²∇×B` and `∂B/∂t = −∇×E` precisely.
- Argue which is correct (or whether both are needed and at which scale). The
  measured v93 wall is: **pass U solved within-mode transport, but the door
  scrambles phase and the metabolic cycle never closes.** Your design must
  engage *this* wall, not relitigate within-mode transport.
- Address the **antagonism** (RING_DNLS B.3: binding vs mobility). Does the
  primitive release the antagonism, or repeat it?

## III.3 Physics equation derivation

- From the conservation law `d/dt(Em+Ee+Es) = 0` (closed system) and the
  substrate, derive the **coupled first-order pair** for the chosen
  extension. Show each step.
- Identify the **link current** `J` (the cross-term, as in v93 §II.4) for the
  new pair, and show what the **first moment** of `J` (the momentum) becomes
  at the mode boundary.
- State the conditions under which the pair reduces to a **second-order wave
  equation** (the matter analogue of `□ψ = 0`), and what the propagation
  speed is in vacuum and in loaded matter.

## III.4 Fundamental math & mechanics (DEC fragment)

- State the discrete algebra (0-cochains, 1-cochains, `d`, `δ`, `δJ` updates)
  for the new pair. Take the **minimal fragment** v93 §II.5 used (0-forms,
  1-forms, door source terms; 2-forms and discrete curvature are over-kit).
- **Conservation proof** (the theorem): show that under the pair,
  `Σ(Ee+Em+Es)` telescopes to boundary terms and is conserved exactly
  (Stokes), exactly as v93 §II.4 dissolved conservation into a Givens
  theorem.
- State what **anti-ignition** looks like for the new pair (v93's structural
  argument: random-phase bath self-gates; does the new pair inherit this? is
  there a new coherence-runaway risk at the door?).

## III.5 Numerical simulation method (the kernel face)

- A concrete **update primitive** for the kernel, behind a **new byte-inert
  knob** (default 0, battery-gated; name it — e.g. `amp_pair` or whatever
  fits). Specify:
  - the **state** the new pair reads / writes (the existing `fa1,fa2,Em,th2,
    Es`, or a new minimal slot-borne carrier? — justify any new state);
  - the **operations** (what rotation / hop / door rewrite, in the same
    pairwise-norm-preserving style as pass F / pass U);
  - the **schedule** (where in the pass order: before pass F? after pass U?
    replacing a pass-6 sub-step? — be precise about pass order);
  - the **interaction with the click grain** (the action atom ε stays; the
    pair is *between* clicks — derive how);
  - **byte-inertness at 0** (the exact branch that reproduces v93's V3a
    byte-for-byte when the knob is 0).
- The **C / Go mirror** parity requirement (v93 §IV.12: expect tolerance-based
  abx at strength, byte-identity at the inert default).

## III.6 Acceptance bars (pre-registered; the ratchet)

Pre-register the bars that would *decide* the design, in the v93 style. At
minimum:

- **M1 (metabolic propagation):** the intake/outtake pair propagates a
  metabolic disturbance across a chain of cells at a finite, measurable speed
  (the matter analogue of `c`). State the meter.
- **M2 (conservation):** drift at FP floor (≤1e-13) under the pair, on live
  matter — not the matterless bath.
- **M3 (anti-ignition):** bath byte-stats / glow / birth rates unchanged at
  the inert default; armed-coherent-bath test specified (v93's unrun L1-C).
- **M4 (door phase retention):** QUENCH-3-style winding retention with the
  new pair, vs the v93 baseline (`R2d≈0.02`). State the target.
- **M5 (the wall):** does the pair release the binding/mobility antagonism?
  State the discriminating experiment (RING_DNLS-style).
- **Invariant surface:** 87 V3a bars GREEN at the inert default (the ratchet).

## III.7 What NOT to do (the hard constraints, repeated)

- No background, no imported field/species, no SI constants.
- No state relabeling (two floats with no transport change = nothing).
- No per-component physical jobs (the rejected quark road).
- No unitarization of the click grain (the door stays quantized).
- No sum-then-square.
- No feedback / shadow / bias correction (v92 closed).
- No phase written to the cell clock `th2` at the door (carry it slot-borne;
  v93 §IV.6 / QUENCH-3).
- No bar softening; the ratchet governs.

## III.8 What to flag explicitly

- Where the design **completes v93's open wall** vs where it **opens a new
  risk**.
- Whether the primitive requires a **new sector** (a fourth mode) or lives in
  the existing three. (Programme rule: prefer the existing three; a new mode
  is a major face and needs explicit user sign-off — flag, don't assume.)
- The **two or three highest-value first experiments** and what outcome
  decides each.
- The **most likely failure mode**, honestly stated.

---

# PART IV — substrate reference (read-the-kernel anchors)

You have `--add-dir` on `v94/`. Key anchors:

| what | where |
|---|---|
| pass F (field transport: precession + Givens hops) | `v94/kernel/freecell.c:2151` |
| pass U (v93 unitary dense hop, the cousin of pass F) | `v94/kernel/freecell.c:1591` |
| pass 6 (the door: beat, condensation, radiance) | `v94/kernel/freecell.c:2198` |
| radiance outtake (graded sub-cap emission, `k_rad·cap·x^p_rad`) | within pass 6 |
| `q_detune` (load-flattens-pitch: `w₂e = w₂/(1+q·x)`) | local-clock update |
| action atom `ε = A₀ω/2π` (the click grain) | `atoms_w`, `atoms_fire`, `atoms_clamp` |
| law table (V2g core + radiance + workfn) | `v94/laws/laws_V3a.cfg` |
| v93 closeout (the wall) | `v93/theory/FINAL_STATUS_ADOPTION.md` (read access via v93/) |
| v93 RECKONING (Face C transient gate) | `v93/RECKONING.md` |
| v93 RING_DNLS (binding/mobility antagonism) | `v93/RING_DNLS.md` |
| v93 metabolic framing (the open cycle) | `v93/theory/METABOLIC_MATTER.md` |

Read the kernel anchors; you do not need to read the entire 4800-line file.
The relevant passages are the three passes (F, U, 6) and the door bookkeeping.

---

**Deliverable:** `v94/consult/DESIGN_claude.md` — a cohesive developmental
design covering III.1–III.8. Be a physicist: derive, prove, specify. Do not
write kernel code; specify the primitive precisely enough that the opencode
orchestrator can implement it from your design. Flag every load-bearing
assumption. Pre-register the bars.

*The user has authorized full v94 kernel authority. v93 stays frozen.*
