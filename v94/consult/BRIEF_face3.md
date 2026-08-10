# v94 face 3 — claude developmental session: the unification (matter as a self-bound state of the field via χ)

You are the **developmental session** for **v94 face 3**. This brief is **self-
contained** but builds on your face-1 design (`v94/consult/DESIGN_claude.md`)
and the face-1 implementation record. You have read access to the v94
directory (`--add-dir /home/d/code/scp/v94`); read these first:

1. `v94/consult/DESIGN_claude.md` — **your** face-1 design (the spec).
2. `v94/IMPLEMENTATION_STATUS.md` — the face-1 implementation + bar status.
3. `v94/consult/REVIEW_grok.md` — grok's review of face 1 (SOUND-WITH-FIXES;
   one blocking bug found and fixed; flags ACCEPTABLE).
4. `v94/kernel/freecell.c` + `v94/fab/step.go` — the implementation (door_g
   knob, pass P, settlement).

Produce a developmental design for **face 3** at `v94/consult/DESIGN_claude_face3.md`
— theory, fully derived symbolic equations, numerical simulation method. No
kernel code; specify the primitive precisely enough for the opencode
orchestrator to implement.

---

# PART I — the user's question (the unification)

## I.1 The symmetry the user named

In face 1, the central claim of your design was that the door is a **third
link type** — a link in the **mode** direction instead of the **space**
direction — and the same primitive (unitary pairwise rotation = Givens
bind/unbind) that propagates the field (pass F) also mediates the door pair
(pass P) and the settlement (the click). Face 1 **unified the primitive**
across spatial and modal edges of one graph.

The user's question goes one rung deeper:

> Both EMF translation and matter intake/outtake rely on a bind/unbind during
> translation in the field. The matter becomes more complete or different in
> situation, but the underlying mechanism should be the same mathematical,
> theoretical mechanism. Is there a plan to explore that symmetry?

Read literally: today the substrate carries **two distinct complex amplitudes**
per cell — `ψ_e = fa1 + i·fa2` (the field) and `ψ_m = √Em · e^{iθ₂}` (matter) —
plus the persistent register `χ = dc1 + i·dc2` you introduced in face 1.
That's **three** complex carriers per cell. The user's claim is that this
*duality of state* (field vs matter as separate sectors) is a defensive
structure, not a theoretical commitment, and that the deeper unification is
**matter as a self-bound state of the one field** — with bind/unbind as the
universal primitive, and the click grain **emerging** from binding geometry
rather than being imposed at a door.

## I.2 Where face 1 already points at this

Three places in your own design:

- **§1.3 (refractive effect):** observation 2 — light through matter is slow
  *because it alternates between two incomplete conversions with the medium*
  (the Lorentz-oscillator / polariton mechanism). That is already "matter as
  field dressing."
- **§2.3 (the register χ):** you argued χ is *mode-transit bookkeeping*
  (cousin of `slem`), not a fourth free energy mode. Grok agreed (flag (a)
  ACCEPTABLE). The cousin-of-flight framing is exactly the "incomplete
  conversion needs somewhere to stand" line.
- **§2.6 (why not registerless direct `ψ_e ↔ ψ_m`):** you rejected the direct
  rotation because it would (i) unitarize the door, (ii) leave vacuum with no
  virtual/real distinction (ignition), (iii) re-import the antagonism. **Face 3
  asks: are those rejections permanent, or are they the open work?**

## I.3 The v89 PRINCIPLE, taken literally

`v89/PRINCIPLE`: "matter is converted space; space curves because energy is
conserved and space is convertible; there are no equations of motion because
nothing has a trajectory." Taken literally, matter is a **configuration** of
the convertible mode, not a separate sector. The face-3 question is whether
the kernel can be reduced to **one amplitude ψ per cell** plus **space Es**,
with matter = a bound state of ψ (the unitary descendant of the q-ball), and
intake/outtake = the bind/unbind cycle.

The decisive test, stated sharply: **can ψ_m be retired, with matter re-
derived as a self-bound configuration of χ (which is itself sourced by ψ_e
via the door pair)?** If yes, the click grain is no longer imposed — it
emerges from χ's binding topology. The unification is complete.

---

# PART II — what face 3 must produce

## II.1 Theoretical underpinnings (the unification, derived)

- State the **postulates** of the unified substrate. Identify what plays the
  role of `ψ_e`, `ψ_m`, `χ`, `Es` in the unified picture (is it one amplitude?
  one amplitude + its conjugate? an amplitude + a binding field?).
- Derive **matter as a bound state** from conservation + the unitary
  primitive. Show the conditions under which a self-bound configuration of
  the field exists, is stable, and carries the matter-analogues of the v93
  measurements (Em inventory, cap defense, radiance outtake).
- Address the **three rejections** of §2.6 one by one. For each: either
  (a) show the unified design structurally resolves it, or (b) name the
  minimal defensive structure that must remain (e.g. why a click grain is
  still needed, and what geometry produces it).
- Address **anti-ignition** in the unified picture. Face-1 layers 1–4 (your
  §4.3) were sector-bound; the unified design must re-derive the
  random-bath-self-gates argument without a separate matter sector.

## II.2 Fully derived symbolic equations

- The **unified equation of motion** for ψ (and χ, and Es — whichever your
  postulates require), in continuous form and in the DEC fragment (0-forms,
  1-forms, `d`, `δ`, source terms — the minimal fragment you used in face 1
  §II.5 / III.4).
- The **self-binding mechanism**, derived symbolically. Candidate starting
  point: a `q_detune`-analogue on χ — `Ω_χ(χ) = w_χ / (1 + q_χ · |χ|²/cap)` —
  but you may derive a different form. Show that this admits self-trapped
  solutions; derive their existence region and stability (the analogue of
  RING_DNLS B.1's "qd≥0.6 condenses").
- The **emergent click grain**, if it can be derived. Under what conditions
  does the unified dynamics quantize conversion into whole action atoms
  `ε = A₀ω/2π` *without* an imposed integrate-and-fire? If it cannot be
  derived, name what minimal grain-imposing structure must remain and why.
- The **conservation theorem** (your §III.4.2 style), restated for the
  unified dynamics. Show the conserved sum, the link currents, and that
  drift is FP-roundoff by construction.
- The **field/matter correspondence**: a dictionary mapping observables in
  the unified picture to the v94 face-1 observables (Em, Ee, chi2, door_in,
  J_D, R2d, etc.). State which v94 measurements the unified design must
  reproduce to be accepted.

## II.3 Numerical simulation method (the kernel face)

- A concrete **update primitive** behind a new byte-inert knob (default 0;
  name it — e.g. `chi_self` or whatever fits). Specify:
  - the **state** (does the unified design retire `Em`/`th2`? add new state?
    reuse `χ` with new dynamics?);
  - the **operations** (what rotation / binding / click-emergence step, in
    the same pairwise-norm-preserving style);
  - the **schedule** (where in the pass order; interaction with pass P, pass
    F, pass 6);
  - **byte-inertness at 0** (the exact branch that reproduces face 1 byte-
    for-byte when the knob is 0).
- The C / Go mirror parity requirement (face 1 closed the v93 Go-only debt
  for the door pair; face 3 must maintain it).
- **Acceptance bars** (the ratchet), pre-registered. At minimum:
  - **U1 (existence):** the unified substrate admits a self-bound
    configuration that plays the role of matter (Em inventory, cap defense,
    radiance outtake all reproducible).
  - **U2 (conservation):** drift at FP floor under the unified dynamics.
  - **U3 (anti-ignition):** random-phase bath self-gates in the unified
    picture, structurally re-derived.
  - **U4 (click emergence or imposition):** either the click grain is
    derived, or the minimal impositional structure is named and motivated.
  - **U5 (face-1 correspondence):** the unified design reproduces the v94
    face-1 measurements (the door pair's chi-filling, atom-scale settlement,
    metabolic propagation if M1 has run by the time you read this).
- The **most likely failure mode**, honestly stated.

## II.4 What to flag explicitly

- Where face 3 **completes the unification** (matter retired as a separate
  sector) vs where it **names a permanent residue** (some structure that
  cannot be eliminated).
- The **interaction with face 2** (the stimulated outtake, in flight — see
  Part III). If face 3 needs face-2 results to settle a claim, FLAG the
  result you need as a headline; produce a tentative answer in its absence;
  the orchestrator will resume your session when face 2 lands.
- Any **programme rule** the unified design touches (no background; no
  imported species; no quark road). The state reduction to "one amplitude"
  is exactly the boundary between "unification" and "the rejected pre-v89
  multi-component lattice on a background" — argue carefully which side the
  unified design is on.

---

# PART III — in-flight work that may produce results while you run

The orchestrator is running, in parallel with this session:

- **Face 2 — stimulated outtake** (your §2.5 "Direction"; the second half of
  the metabolic cycle). Implementing the negative-`C` settlement branch
  (matter anti-phased with χ un-settles into the register, ledgered on
  outtake credits). Expected: closes the bidirectional cycle.
- **M1 — metabolic propagation** (chain/ring, arrival-time vs distance).
- **M3b — armed coherent bath** (the unrun L1-C; phase-coherent lit bath
  over a seeded clump; measurement of the coherent-pumping threshold).
- **M4 — door phase retention** (lit-bath ring with door traffic firing;
  baseline `R2d ≈ 0.02`, target `R2d ≥ 0.5`).
- **M5 — antagonism release** (RING_DNLS B.3 rerun).
- **M6 — metabolic fixed point** (the I=O park).

If you need any of these to settle a claim, **flag the result as a
headline at the top of your design output** ("NEEDS: M4 outcome to settle
§X.Y"), produce a tentative design in its absence, and the orchestrator will
**resume your session** (`claude --resume <session-id>`) with the result.
Your session ID will be in the JSON output of this run; please also print it
at the head of your design output for redundancy.

---

# PART IV — programme constraints (permanent; same as face 1)

- No background (nothing persists and is merely re-valued).
- Energy is never destroyed — the only law. Conservation as a theorem of the
  update, never a patched ledger.
- No imported field or species. No SI constants. No Maxwell import.
- The door stays special — UNLESS face 3 derives the click grain from
  binding geometry (then it would emerge, not be imposed; argue which).
- No quark road (no per-component physical jobs; no separate ledgers/
  clocks/identity that outlives the cell).
- No sum-then-square. Pairwise hops only.
- No feedback / shadow / bias (v92 closed).
- The ratchet governs; bars sharpen by measurement, never soften to pass.

**Deliverable:** `v94/consult/DESIGN_claude_face3.md` — a cohesive
developmental design covering II.1–II.4. Print your session ID at the head
of the file. Be a physicist: derive, prove, specify. Flag every load-
bearing assumption; flag every in-flight result you need; produce tentative
results in their absence.

*The user has authorized full v94 kernel authority. v93 stays frozen.*
