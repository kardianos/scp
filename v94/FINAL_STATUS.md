# v94 final status — face 1 + face 2 + face 3 + DS-G + M-series

**Opened:** 2026-08-09. **Status:** face 1–3 implemented; DS-G confirmed;
M6 approaching park; transport physics handed to v95 (J-metric).

## Knobs implemented (all byte-inert at default 0; 87/87 battery GREEN throughout)

| knob | default | what |
|---|---|---|
| `door_g` | 0 | v94 door-pair coupling (pass P: field↔register Givens + settlement: χ→ψ_m solved Givens in atoms) |
| `door_avg` | 0 | face-3 rev Y.4: beat-integrated cross-term for the outtake (synchronous detection) |
| `chi_load` | 0 | face-3 §III: register self-detuning (includes \|χ\|² in the load x) |
| `door_out_enable` | 1 | face-2 gate: 1=stimulated outtake active; 0=face-1 intake only |
| `sw_kx`, `sw_amp` | 0 | DS-G3 apparatus: standing-wave field seed |

All in both kernels (C `kernel/freecell.c` + Go `fab/step.go`). Law-header
lines in both. Drift meter includes \|χ\|². BH flushes χ. 87/87 GREEN.

## Bars — the full scoreboard

### CONFIRMED (conversion physics won)

| bar | status | evidence |
|---|---|---|
| **Invariant surface** (87 bars) | **GREEN** | battery ALL GREEN 87/87 at all knobs=0 |
| **M0 parity** | **GREEN** | C≡Go byte-identical at door_g=0 (the amp_nat>0 debt is pre-existing v93) |
| **M2 conservation** | **GREEN** | drift ≤ 1e-13 at door_g≤0.3; plateaus at ~1.5e-13 at 0.5 (FP noise, not leak) |
| **M3a anti-ignition (dark)** | **GREEN** | structural — dark bath byte-identical door_g=0 vs 0.3 |
| **M3b anti-ignition (armed)** | **GREEN** | no runaway through noise_amp=2.0; cap+evap+face-2 defend |
| **DS-G3 (no-cancellation)** | **CONFIRMED** | r(Em,Ee)=+0.684 (door pair) vs −0.259 (legacy); standing wave, frozen field, face-1-only |
| **M1 side-finding (anti-runaway)** | **CONFIRMED** | door pair Em=17 (confined) vs legacy Em=417 (runaway box-fill) — pair cures legacy pathology |
| **M6 metabolic park** | **APPROACHING** | chi_load=1 creates 10× slope reduction + bi-directional convergence (amp 0.5→127, 1.0→192, 2.0→233) |

### NEGATIVE (transport physics lost, with named confound)

| bar | status | evidence |
|---|---|---|
| **M4 phase retention** | **RED** (confounded) | R2=0.24 (door_avg=0) vs baseline 0.62. door_avg=1 helps: R2=0.32. Confound: stroboscopic aliasing. |
| **M5 antagonism** | **RED** (confounded) | COM≈0 (anchored). door_out dominates 4:1. Same confound. |
| **M1 propagation** | **NOT OBSERVED** | pair confines, doesn't propagate. The metabolic wave needs the J-metric. |

### The confound and the fix

**Confound:** the face-2 outtake samples a rotating quantity (C = 2·Im(ψ_m*·χ))
at the beat, which is stroboscopic aliasing. This erodes coherent structure
(M4 winding, M5 translation) and prevents matter accumulation (DS-G with
outtake).

**Fix:** `door_avg` (beat-integrated D instead of instantaneous C). Partially
tested — M4′ improved R2 0.24→0.32. The full campaign at door_avg=1 for M4′/M5′
needs the door_avg accumulator performance optimised (currently 50% overhead).

**Next face:** J-dependent link metric (v95) — the space-direction transport
partner. The door pair handles conversion (mode direction); the J-metric
handles coherent spatial transport. Together: metabolizing, translating,
phase-retaining matter.

## Document map (v94)

| path | what |
|---|---|
| `consult/DESIGN_claude.md` | face-1 design (claude) |
| `consult/DESIGN_claude_face3.md` (80 KB) | face-3 unification design (claude, 6 resume cycles) |
| `consult/NOCANCELLATION_claude.md` (24 KB) | no-cancellation theory + DS-G3 confirmation |
| `consult/REVIEW_grok.md` | grok's face-1 implementation review (SOUND-WITH-FIXES, fixed) |
| `consult/FACE2_M3B_RESULTS.md` | face-2 + M3b results |
| `consult/M4_RESULTS.md` | M4 negative + door_avg fix |
| `consult/DSG_RESULTS.md` | DS-G slit + standing wave results |
| `paper/scp_theory.pdf` (323 KB) | theory paper (revised per grok, format fixed) |
| `paper/REVIEW_grok.md` (32 KB) | grok's pedagogical review of the paper |
| `paper/REVISION_PLAN.md` | claude's revision plan |
| `kernel/freecell.c` + `fab/step.go` | both kernels (door_g, door_avg, chi_load, door_out_enable, sw_kx/sw_amp) |

## Claude session record

| session | what | output |
|---|---|---|
| `aebbc7eb-4716-4e30-8c22-484e115428f2` | v94 face 1–3 + paper + no-cancellation + 6 resumes | 80 KB face-3 design, 24 KB no-cancellation, 323 KB paper |
| v95 (new) | J-dependent link metric design | in progress |
