# v96 RESUME — stage-setting (FROZEN; active programme is v97)

W-A is reinterpreted by v97 P5, not adopted. W-C remains a lawful
config-only probe if the dipole hunt wants it.

**Written:** 2026-08-09 (end of day).
**Next claude session:** FRESH START (do NOT resume `0fc5d7e7-...`). Start a
new developmental session with a fresh brief.
**Session to read:** `v96/consult/DESIGN_claude.md` (76 KB, the full v96
design including §W stage-setting from the resume).

---

## Claude's §W diagnosis (the two missing pieces)

The four measured absences (no unique size, no preferred winding, no 1/r
gravity, no decomposition spectrum) reduce to **two structural gaps**:

1. **No consumption** — the metabolic cycle has no space leg (door:
   field↔matter; radiance: matter→field; space only click-share bookkeeping
   that nets zero on a parked object). Nothing in the substrate out-eats its
   medium at accessible scale. S0 measured this: α ≈ 1.1 (volume-limited);
   the field diffusion length (~50) >> object size (~5-10). The 1/r far field
   IS the supply chain of stable matter seen from outside; its absence is why
   there's no unique N* and no gravity.

2. **No pricing** — winding slips are free. When the ring's winding changes by
   ±1, no energy is paid. So all windings are degenerate (Q0: all m equally
   stable) and decay is continuous (Q1: no fragment spectrum). Making the slip
   pay an action atom would make winding a PRICED INTEGER — the baryon-number
   analogue.

## Candidate completions (claude's recommended order)

| candidate | what | scope | decides |
|---|---|---|---|
| **W-C** (config-only, run now) | bulk topology test — vortex in filled droplets at N∈{12,30,60,100} with j_metric=1; plus the localized smash (Q1 instrument debt) | config + 2 seed generators, no kernel change | Is winding-as-number scale-blocked or structure-blocked? |
| **W-B** (small knob, kernel change) | the slip quantum — a 2π slip pays one atom A₀ to radiance, gated to J-metric-dug links | one knob, battery-gated, needs user authorization | Makes W a priced, booked integer |
| **W-A** (space-sector law face, v97) | steady Es→Em conversion as upkeep (matter consumes space). Creates: 1/r depression δEs ∝ Q̇_s/(4πD_s·r), gravitational attraction, AND accessible N* ≈ 244 cells (the space sector's slow transport puts the transport-limited regime at N~10² instead of 10⁵) | major face, needs consultation round + authorization | The missing 1/r + accessible size selection |

## User's directives for the next session

1. **Start claude fresh** — do not resume the v96 session. Write a new brief.
2. **Run W-C first** (the cheap decider). It's config-only; no kernel change.
3. **The missing 1/r gravity is the deep observation.** The user notes:
   "We are missing 1/r gravity, as such something fundamental is missing."
   W-A (the space-sector law face) is the candidate fix — it's the v97 face.
4. **Winding is a number, not a geometry.** The user notes: "The winding is
   probably more of a number and less of a geometry, which means it isn't
   fully being represented in effect." W-B (slip pricing) or W-C (bulk
   topology emergence) are the candidate fixes.
5. **Q1 results were expected negative** — confirmed. The smash test on the
   C6 ring doesn't produce a fragment spectrum.

## The v94→v96 campaign record (for the next session's context)

### Confirmed (physics that works)

| result | measurement |
|---|---|
| Coherent conversion (door pair) | DS-G3: r(Em,Ee)=+0.684 (vs legacy −0.259) |
| Phase retention (J-metric) | J1: R2=0.774 (best ever); C1-face1: R2=0.862 |
| Anti-runaway confinement | M1: door pair Em=17 vs legacy 417 |
| Metabolic park (approaching) | M6: chi_load=1 gives 10× slope reduction + bi-directional convergence |
| Translation (face-1 only) | com_dev median 1.27 (3.7× baseline) — hidden by parsing bug, corrected via scpan |
| Conservation | drift at FP floor across all configs |
| Anti-ignition | structural (dark bath byte-identical) |

### Negative (physics that doesn't work yet)

| result | measurement | diagnosis |
|---|---|---|
| Unique size (N*) | α≈1.1 at accessible scale | volume-limited; needs W-A for accessible surface limitation |
| Phase retention (full door pair) | R2=0.07 (outtake erodes) | the stimulated outtake is the single pathology |
| Translation (full composition) | COM≈0.84 (outtake drains) | face-1 only works (1.27); outtake halves it |
| Preferred winding | all m equally stable (Q0) | winding is geometry not number; needs W-B or W-C bulk emergence |
| Decomposition spectrum | ring degrades continuously (Q1) | no fragment count; the bound state doesn't shatter |
| 1/r gravity | π-flat (no monopole) | no space consumption; needs W-A |

### The structural pattern

The door pair's stimulated outtake (face 2) is the single pathology across
M4/M5/DS-G/C1. With it disabled (door_out_enable=0), the face-1 compositions
work (translation, retention, DS-G). The outtake needs either door_avg
(synchronous detection, partial fix) or a fundamental redesign.

The J-metric is an amplifier, not a generator. It needs a pre-existing
coherent current to amplify. The ring has one (winding); the blob doesn't.

### Tools built

- **`scpan`** (`v95/cmd/scpan/main.go`) — structured parsing of fabrun/freecell
  output by column name. Skill at `.opencode/skills/scp-analysis/SKILL.md`.
- **Standing wave seed** (`sw_kx`, `sw_amp`) — for DS-G spatial pattern tests.
- **`door_out_enable`** knob — disables the stimulated outtake for face-1-only tests.

### Active claude sessions

| session ID | version | topic | status |
|---|---|---|---|
| `aebbc7eb-4716-4e30-8c22-484e115428f2` | v94 | face 1-3 + paper + no-cancellation (6 resumes) | complete |
| `415d0caa-80d0-4689-b10a-ea15e29ba881` | v95 | J-metric design (1 resume) | complete |
| `0fc5d7e7-9f51-4efc-a43e-ce358e708e4c` | v96 | size selection + quark ontology (1 resume) | complete |

All three are complete. The next session starts FRESH.

### Knobs implemented (all in v95/kernel, all byte-inert at default)

| knob | default | what |
|---|---|---|
| `door_g` | 0 | door-pair coupling |
| `door_avg` | 0 | beat-integrated cross-term |
| `chi_load` | 0 | register self-detuning |
| `door_out_enable` | 1 | face-2 gate (0=face-1 only) |
| `j_metric` | 0 | J-dependent link metric |
| `j_tau` | 4 | J-metric low-pass timescale |
| `j_rel` | 0 | J-metric healing |
| `sw_kx`, `sw_amp` | 0 | standing-wave field seed |

87/87 battery GREEN at all defaults throughout.

---

*Stage set. The substrate converts, retains, and parks — but it does not yet
consume (no 1/r) or price (winding is geometry not charge). W-C is the cheap
decider; W-A is the v97 face.*
