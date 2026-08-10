# v94 — the EMF translation primitive, mirrored to matter intake/outtake

**Opened:** 2026-08-09.
**Authority:** user-directed campaign, opened on the user's three observations
about EMF. **You have full authority to write and manage the numeric kernel in
v94** (the v93 kernel-protection rule is **lifted for v94/**; the v93 kernel
stays frozen).
**Carried substrate:** v93 (its kernel copied verbatim into `v94/kernel/`,
both `freecell.c` and the Go mirror in `v94/fab/`). v93 closed with a measured
wall; v94 is the next theory campaign.

## The user's three observations (the seed)

1. **EMF translates at "c"** in vacuum (when not moving through matter).
2. **EMF translates slower** through matter (a refractive / load effect).
3. **EMF alternates between two expressions** — as if there are two alternating
   *incomplete* conversions going on — which the user draws as a **parallel to
   the v93 matter intake/outtake** attempt.

## The user's hypothesis (the line of investigation)

> If we are able to find / construct / discover a **primitive translation
> device for EMF**, the v93 matter intake/outtake becomes adoptable.

In real Maxwell EM, propagation is precisely the **mutual induction** of two
expressions: changing **E** circulates a **B** (spatially offset), changing
**B** circulates an **E** (spatially offset); each conversion is **incomplete**
(first-order in time), and the pair, mutually induced, **is** the propagation
of the wave at c. The user reads this as a primitive "translation device" and
asks: can the v93 dense channel's **intake** (field→matter condensation at the
door) and **outtake** (matter→field radiance) be cast as **two alternating
incomplete conversions** in the same sense, so that a metabolic / propagating
cycle of matter emerges by construction — the same way an EM wave emerges from
Maxwell's two first-order equations?

This is **not** the v93 thesis (v93 was: give matter the field's *within-mode*
unitary transport algebra). v94 is one rung deeper: ask **why** the field
propagates as two alternating expressions, derive the primitive, then mirror
that primitive to the **mode-boundary** (intake/outtake) of matter — which v93
explicitly left untouched ("the door is never unitarized", §IV.3).

## What v94 must produce (the deliverable)

1. **Theory.** A cohesive, fabric-only derivation of the EMF translation
   primitive in this substrate (no Maxwell import — derive it from the
   sector-coupling the kernel already has), stated as the **two alternating
   incomplete conversions** the user named. Then the **mirror to matter**: how
   intake/outtake, written in the same primitive, becomes a propagating
   metabolic cycle (the "matter wave" / coherent temperature wave).
2. **Physics equation derivation.** From conservation + substrate to the
   coupled first-order pair, with explicit identification of the two
   "expressions" in each sector.
3. **Fundamental math & mechanics.** The discrete algebra (DEC / Givens /
   Stokes fragment v93 already used) stated for the new pair; conservation
   proof; what the link current / momentum becomes.
4. **Numerical simulation method.** A concrete kernel face — the new update
   primitive, behind a byte-inert knob (default 0), battery-gated, C and Go
   mirrored, with the same ratchet discipline v93 used.
5. **Acceptance bars** (pre-registered, sharpened by measurement, never
   softened to pass).

## What v94 must NOT do (the hard constraints — read before designing)

Inherited verbatim from v93 §IV (the project rules that were never measured
false):

- **No background.** Nothing persists and is merely re-valued. The roster is
  dynamical; no fixed lattice; no fixed fiber with separate ledger/identity.
- **Energy is never destroyed — the only law.** Conservation must be a
  theorem of the update (pairwise norm-preservation), never a patched ledger.
- **No imported field or species.** No Maxwell import, no SI constants, no
  electron species. The two "expressions" of EMF must be **derived** from the
  fabric.
- **The door stays special.** Quantized action atoms ε=A₀ω/2π, work-function
  gate, two-atom credit — these are the "atoms at mode boundaries." Do not
  dissolve the click grain into a continuous rotation.
- **No quark road.** Two state-halves with separate physical jobs is the
  rejected pre-v89 lattice.
- **No sum-then-square.** Never `|ΣA|²`; pairwise hops only.
- **No feedback / shadow / bias correction.** v92 closed three forms; the
  additive ledger rejects the interference term. v94 is a *transport-primitive
  replacement*, not a correction.
- **The ratchet.** Every kernel change runs the full battery before commit;
  bars sharpen by measurement, never soften to pass; the 87-bar V3a surface
  must stay GREEN at the inert default.

## How v94 is run (the process)

1. **(DONE)** A self-contained brief went to **claude** (developmental session):
   `consult/BRIEF.md` → `consult/DESIGN_claude.md` (12 min; theory, equations,
   math, mechanics, numerical method).
2. **(DONE)** Implemented in `v94/kernel/freecell.c` + Go mirror `v94/fab/`,
   behind a new byte-inert knob `door_g`, battery-gated. Status:
   `IMPLEMENTATION_STATUS.md`.
3. **(DONE)** Grok CLI reviewed the implementation (`consult/REVIEW_grok.md`):
   verdict **SOUND-WITH-FIXES**; one blocking bug found (ΔEm solver sign);
   fixes applied; re-verified. Status now SOUND.
4. **(partial)** Shepherded the runs: M0/M2/M3a/M3b verified GREEN; the
   primitive is observationally live (door pair settles ~20 units of energy
   into matter that the legacy door does not, at door_g=0.5). M1 directionally
   positive. M4/M5/M6 (the full campaign arms) need their own focused session
   — each requires dedicated apparatus (FCS dump of register plane, ring
   winding analyzer, RING_DNLS B.3 baseline, metabolic fixed-point sweep).

## Status as of 2026-08-09 (face 1 closeout)

**Implemented and verified:**
- `door_g` knob + `dc1, dc2` register (both kernels); `total_energy` includes
  `|χ|²`; law-header line in both.
- Pass P (the door pair, Strang-symmetric around pass F) — per-cell register
  precession by `w2e·dt/2` + fiber Givens by `0.5·door_g·w2·dt` (bare w2,
  NOT load-detuned — claude's workfn-as-detuning choice).
- Pass-6 settlement (the click): `if (door_g > 0)` replaces the legacy
  `Ee`-condensation with a solved Givens between χ and ψ_m by the angle that
  moves exactly d1 (the fired energy) — in whole atoms, workfn-gated, with
  the existing atom machinery.
- Grok's blocking ΔEm-sign fix applied + verified algebraically and on a
  coherent ring. BH χ-flush added.
- 87/87 battery GREEN at `door_g=0` (byte-inertness preserved).

**Bars verified:**
| bar | status |
|---|---|
| Invariant surface (87 bars at `door_g=0`) | **GREEN** |
| M0 C≡Go parity at strength | **GREEN** (byte-identical, exceeds tolerance-abx) |
| M2 conservation (drift ≤ 1e-13 on live matter, T≥1000) | **GREEN at door_g≤0.3**; borderline at 0.5 (drift plateaus at 1.5e-13, does not grow) |
| M3a anti-ignition (dark bath) | **GREEN (structural)** — byte-identical `door_g=0` vs `0.3` |
| M3b armed coherent bath (face 1) | **GREEN** — kx=1.1 + noise_amp=0.8 + door_g=0.5, T=500: bounded, cap defends, no runaway |
| Functional (register fills; settlement fires in atoms) | **GREEN** |
| Solver correctness (post-grok-fix) | **GREEN** — algebraic unit check ~1e-16; coherent ring end-to-end |

**Bars not yet run (face 2 / next session):**
- M1 metabolic propagation (chain arrival-time vs distance) — directionally
  positive (door pair settles ~20 units at door_g=0.5); full bar needs a 1D
  chain seed.
- M4 door phase retention (R2d ≥ 0.5 on a lit-bath ring with cond>0) — needs
  FCS dump of register plane + ring winding analyzer; claude predicts this is
  where face 1's limit shows (settlement may decohere).
- M5 antagonism release (RING_DNLS B.3 rerun).
- M6 metabolic fixed point (the closing of I=O at an isolated x*).
- The stimulated outtake (claude §2.5 "Direction") is **face 2** — required
  for the full metabolic cycle and for the M5/M6 bars.

## Document map

| path | role |
|---|---|
| `v94/README.md` | **this file** — the charter + status |
| `v94/consult/BRIEF.md` | the self-contained brief sent to claude |
| `v94/consult/DESIGN_claude.md` | claude's developmental design output (the spec) |
| `v94/consult/GROK_REVIEW_REQUEST.md` | the brief sent to grok |
| `v94/consult/REVIEW_grok.md` | grok's review of the implementation (SOUND-WITH-FIXES) |
| `v94/IMPLEMENTATION_STATUS.md` | implementation summary + bar status (post-fix) |
| `v94/kernel/freecell.c` | C kernel of record (modifiable in v94) |
| `v94/fab/step.go` | Go mirror of the kernel (campaign kernel) |
| `v94/laws/laws_V3a.cfg` | the adopted law table (V2g + radiance + workfn) |
| `v94/battery.log` | the final battery run (ALL GREEN 87) |
| `v93/theory/FINAL_STATUS_ADOPTION.md` | v93 closeout (frozen; the substrate v94 stands on) |

---

*STEP ZERO: the v94 tree carries v93's kernel verbatim; the battery must pass
87 GREEN at the inert default before any face is opened. The kernel-protection
rule of CLAUDE.md applies permanently to `v93/` and earlier; v94's kernel is
explicitly authorized for modification under this campaign.*
