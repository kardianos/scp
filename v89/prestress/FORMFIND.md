# FORMFIND — form-finding on the real foam (PRESTRESS stream B)

*2026-07-28/29. Solver: `prestress/formfind.py` (pure numpy/scipy; the
simulator was never run). Foam: `prestress/foam/foam_s20260727.tsv` — the
campaign standard (9741 cells, L=24, d̄=1.5053), identical to every
`seed=20260727` kernel run, so cell picks and phases transfer verbatim.
Outputs: `prestress/candidates/*.net` (+ `_ctrl.net` twins,
`*_report.txt`, `summary.tsv`). Deterministic: fixed seeds
(FSEED=20260727 + 101·candidate-index); identical reruns give identical
bytes.*

---

## 1. Verified kernel facts the solver is built on

Read from `v89/cellfab.c` (2026-07-28 state) and `battery/laws_V2g.cfg`:

* gate(ψ) = ((1+cos ψ)/2)^8; forward gate of link i→j:
  ψ_f = wrap(θ_i − ω_i·d/C − θ_j); backward mirrored with ω_j.
* ω = w2/(1+q·x), w2=2.9, q=1.2, C=1; seeded Em = x·cap (cap=2.5, the
  add + s_pull mechanics net to exactly x·cap; verified against the
  `# LUMP` voice masses in both validation logs).
* Link rule: i,j linked iff d < 1.15·(r_i+r_j). Loader reproduces the
  foam: 9741 cells, mean degree 17.63, d̄=1.5053; 85884 links vs the
  kernel's 85886 — the TSV carries 4-decimal coordinates/radii, and
  exactly 2 links sit within rounding of the rule boundary. No pick or
  gate in this document is affected (structures live ≥0.06 from any
  boundary case).
* Ring seeder (init=7, lines ~916–991) and shell seeder (init=8) pick
  loops replicated exactly: ascending id order, strict `<`, `cflag`
  all-zero at seed time (edge sinks are flagged after seeding — print
  order in both logs proves it), shell refinement = 3 sequential passes,
  search radius √1.44, core-ball exclusion r_core=0.8.
* `init=net` reader (lines ~1189–1326): `V x y z xk th2` / `E a b`,
  `#` comments ignored, xk clamped ≥0.02, θ stored mod 2π, picks =
  nearest free cell to the V position (we emit exact foam coordinates,
  so picks are exact and `worst_pick=0`), forward direction of an edge
  = a→b as listed. Parasite scan = the exact link rule over picked
  pairs. All 18 emitted `.net` files conform (V/E/# lines only).

Design algebra (LEDGER translation, confirmed against B1):

* **Strut** (both-gate link): 2ω·d = 2πm with only m=1 reachable ⇒
  d = π/ω, one shared length per component; optimal phases are
  antiphase (θ_i−θ_j = π), each gate then sits at ψ = π − ω·d. All-strut
  even cycles are automatically consonant (targets sum to 0 mod 2π):
  a cube's post-retune deficit is **pure length spread**.
* **Cable** (forward-only): any φ = ω·d; θ exact on a tree; every
  independent cycle needs its oriented targets ≡ 0 mod 2π (struts
  contribute π, cables ±ω·d) — `ring_m` generalized. Cable back gate
  = gate(−2φ).
* **Windows**: x ∈ (0.05, 0.9) usable, x ≥ ~0.25 for skirt margin
  (x_skirt = 0.0617) ⇒ ω ∈ (1.394, 2.736), safe ω ≤ 2.231; strut
  lengths π/ω ∈ (1.15, 2.25) intersected with the link window
  [1.0, ~2.07].

## 2. Method

**Picker.** Per-vertex candidate pools (cells within pool_r ≈ 1.4–1.6 of
the template target, all ≥2.8 from box faces — edge sinks live within
1.6). Search = template-rotation search (uniform random rotations, first
one identity) × multistart coordinate descent: every sweep tries every
pool candidate for every vertex and keeps improvements. This strictly
dominates the kernel shell seeder's 3 greedy passes. Objective, minimized
jointly over a 60-point ω grid:
exact strut deficits 2(1−gate(π−ωd)) plus 3× the worst edge; cable-cycle
spread proxy 2·wrap(defect)²/n; **pinned-loop terms** (forced closure
integers, UNwrapped — this couples cell choice to the intended winding;
without it the C5 tube was solved for the wrong pitch); parasites at
their **tree-predicted** gates (θ_a−θ_b ≈ π·n_struts_on_path +
ω·signed cable path length) + 0.25 structural; missing intended links
+25 each.

**Tuner.** ω band from the design (strut rung π/d̄, pinned-cycle
2πm/L, or free window); 240-point grid + golden refine, each point
scored by the exact leak of its LSQ phase solution; loads x = (w2/ω−1)/q
uniform (static lock requires uniform ω per component; Em = x·cap per
voice reported as mass).

**Phase solver.** Targets: strut π, cable ω·d (oriented). BFS-forest
init, integers per cycle from wrapping the tree solution (rounding),
then wrap-iterated weighted least squares on the lock Laplacian
(weights = live-direction lock stiffness, strut:cable = 2:1 —
−gate″(0) = 4 per direction), then L-BFGS polish of the exact leak
objective **including parasite gates**. Spectrum: weighted Laplacian
with per-link stiffness Σ_live −gate″(ψ*); eigenvalues reported, gap =
first above the gauge zero(s).

**Leak accounting** (the scoreboard; per-link table in every report):

* strut: (1−g_f) + (1−g_b) — both directions are design channels;
* cable: (1−g_f) + g_b — forward deficit plus **openness** of the
  unintended back direction;
* parasite: g_f + g_b — any open unintended flow;
* intended edges that are not foam links are excluded from solving and
  scoring and reported as MISSING (see §4 — this class exists).

Note: the task sheet's literal cable term (1 − gate(−2φ)) would score a
perfect one-way ring (comp6 back gate 1.5e−5) as maximally leaky,
contradicting B1's design intent and the measured "back gate 0.100
leaks slowly" reading; it is reported as `leak_alt` in summary.tsv, the
openness convention as `leak`. Where the two disagree (wound designs),
`leak_alt` is punishing exactly what the winding is for.

**Robustness.** Re-run pick+tune at random centers (20 for C1–C5, 12
for C6, 6 for C7; smaller descent budget), collect min/mean gate;
`robust_pctile` = fraction of jitter runs the box-center result beats.

**Controls.** Every candidate ships a `_ctrl.net` twin: same cells, same
loads, phases from a mean-tuned BFS that ignores per-edge lengths (per
edge type, drop = ω·d̄_type). The tuned-vs-ctrl pair isolates exact
retuning in the kernel.

## 3. VALIDATION (the gate for everything else)

**Ring12 (vs `v3_comp12.log`): MATCH.** Kernel pick loop replicated on
the TSV: Lring = **16.542** (ref 16.542, print precision),
closure/2π = **5.0000**, ring_m=5 ⇒ ω = 1.899152, x = 0.439164,
Em/voice = 1.0979 — equal to the log's `# LUMP` voice masses (1.0976–
1.0979) and E_dense(t=0) = 13.1747 vs our mass 13.1749.

**Cube (vs `h1_shell0_v3.log`): MATCH.** Shell seeder replicated
(3 refine passes, r_core=0.8): abar = **1.586**, ω = **1.9809**,
x = **0.3866**, kernel-BFS gates min = **0.001**, mean = **0.597** —
all at print precision.

Both mimics are asserted in `formfind.py --selftest` (PASS).

## 4. Setup findings (new, kernel-relevant)

1. **The H1 cube contained two phantom edges.** On the exact H1 v3
   cells, intended edges 0–1 (d=2.0002, ceiling 1.9261) and 2–3
   (d=1.9866, ceiling 1.9833) are **not foam links**: no channel exists;
   the kernel's 12-edge gate report scored two non-channels, and the
   "cube" was a 10-edge open lattice. The shell seeder never checks the
   link rule. (`# MISSING` lines in `c1_cube125_report.txt`.)
2. **comp12 carries one parasitic link.** v0–v10, d=1.7443, seeded gates
   g_f=0.16 / g_b=0.48 — a half-open unintended channel in the
   longest-lived object the program has produced
   (`c8_ring12_report.txt`, NETGATE P line). The old ring seeder never
   scored it.
3. The feared 12 face-diagonal parasites of the a=1.25 cube do NOT
   materialize on the H1 v3 cells: foam inflation (abar 1.586) pushes
   diagonals out of range; only 3 parasites remain and all sit in the
   gate desert (≤1e−4) at the solved phases. The P2 arithmetic applied
   to the *target* geometry, not the realized one. They remain
   1:1 lock-eligible — the kernel NETGATE P report should watch them.

## 5. C1 — the cube-deficit decomposition (H1 done right)

H1 v3 kernel seed line: 12-edge mean gate 0.597, min 0.001. Decomposed
on the same foam, same cells:

| step | mean gate | min gate | what it isolates |
|---|---|---|---|
| kernel BFS, 12 listed edges | 0.597 | 0.001 | the H1 report |
| kernel BFS, the 10 REAL edges | 0.616 | 0.001 | drop phantom edges |
| exact retune, same cells (ω*=1.9702, LSQ+polish) | **0.676** | **0.130** | phases only |
| free cell choice, same template (rot+descent) | **0.851** | **0.312** | picks |

* **Phases buy +0.06 mean but ~130× on the min gate** (0.001 → 0.130):
  BFS dumps all cycle defect on co-tree edges; least squares spreads it.
  The mean barely moves because for an all-strut cube the phase system
  is consonant — post-retune gates are gate(π−ωd_e), pure length spread.
* **Picks buy +0.175 mean** — three times the phase gain. The freepick
  cube relaxes toward the foam's natural edge (C2 realizes d̄=1.585,
  spread 1.26–1.76).
* **Irreducible remainder ≈ 0.15 mean-deficit** (min gate ~0.3–0.4):
  8-cell cubes in this foam have σ_d ≈ 0.14 at the box center — a
  40-rotation, 16-restart search does not improve it (foam-limited, not
  search-limited). Across 20 random centers: min-gate quartiles
  0.16/0.29/0.44, best 0.78. **No jittered 8-cell cube achieves a
  consonant skin** — quantifying H1 v3's verdict with best-possible
  numbers.
* Parasite leak at seed ≈ 0 (all 3 in the desert; see §4.3).

Emitted: `c1_cube125.net` = the H1 cells exactly retuned (the pure-phase
A/B against h1_shell0_v3; same cells, x 0.3933 vs kernel 0.3866) +
`_ctrl` twin.

## 6. Candidates (box-center solves; full tables in `*_report.txt`)

**C2 cube a=1.5 — best candidate by predicted leak.** ω=1.9631,
x=0.3977, mass 7.95. Gates min 0.401 / mean 0.866; leak 3.227 (best
total and best per live direction, 0.134); 0 missing; 4 parasites, all
≤1e−4 (target diagonals 2.12 clear the max ceiling 2.072; the realized
distortion lets two pairs back in range but the phases close them).
Robustness: 60th pctile; jitter median min-gate 0.29, best 0.78.

**C3 hexagonal prism (all-strut, a=h=1.5) — over-constraint confirmed.**
18 equal lengths demanded, foam supplies 1.06–1.96 (mean 1.67): min gate
0.044, mean 0.837, leak 7.839. 9 parasites and one sits EXACTLY on the
π-rung (4–6: d=1.676, ψ=0.088, gates 0.985/0.985) — a consonant
absorption candidate: declare it a 19th strut. Best mixed split of the
prism IS the C5 tube (rings → m=2 cables): leak/live-direction 0.167 vs
0.436. Robustness pctile 10 (center is a bad spot for prisms).

**C4 ring12 + 2 consonant chords — the chord algebra works.** Single-link
chords are geometrically impossible (2R=5.8; the k=2 skip ≈ 2d·cos15°
≥ 2.4 at any ring the foam supplies — beyond the 2.07 ceiling), and on a
wound m=5 ring every symmetric crossing is frustrated by exactly π/2
(quarter-ring phase 5π/2). Worked solution: **unwound strut ring (m=6,
φ=π) + two 5-hop cable chords at ω·L_chord = 4π** (the chords carry the
winding, m=2 each, hop φ=144°) **+ two cross-struts between the arched
chords** — consonant because the quarter-ring phase is 3π ≡ π.
Realized: ω=1.7652, x=0.5358, mass 26.79 (20 voices), min 0.154 / mean
0.914 (live), fwd mean 0.926; chord-cycle defects ±0.50 rad (spread
~0.045/link); **cross-struts land at gates 0.995/1.000 and 0.996/0.979 —
the weave is real**; two chord-end hops sit near φ=π and become
struts in all but name (back gates 0.84/0.60), i.e. natural consonant
absorption. 6 parasites, worst 0.13. This is the richest consonant
object found: mutual skeleton + one-way circulating chords + winding.

**C5 co-rotating wound tube (2× ring6, B1c) — the rescue holds at
moderate quality.** Pinned-loop steering (both rings forced m=2) is what
made it solvable: without it the picker optimized cells for the wrong
pitch. Realized ω=1.8979, x=0.4400, mass 13.20 — **mass-matched to
comp12 (13.17) within 0.2%**: a clean structural A/B, tube vs ring at
equal mass and near-equal pitch. Rings m=2/m=2 exact, back_sum 0.023
over 12 cables (fully one-way); axial π-struts live both directions
(0.68–1.00). Min 0.055 (one ring cable eats a face-cycle defect), mean
0.883, leak 4.018. 13 parasites (worst 0.48 — the tight d̄=1.10 rings
pull neighbors into range). `leak_alt` = 15.56: the task-sheet cable
formula penalizes exactly the one-wayness the winding is for.

**C6 quad/tri-mesh torus — INFEASIBLE on this foam (quantified).**
Minor-loop closure m=1 forces ω = 2π/L_minor to the bottom of the pitch
window: 3×8 lands x=0.828, 4×6 lands x=0.899 — at the 0.9 usable-load
limit. Row closures cannot be simultaneously integer: homologous major
rows at different latitudes differ in length by ~1.5·r_m·n₂·sin(π/n₂),
and making the row ratio a small-integer ratio (2:1 or 3:2) pushes the
outer row past the 2.07 link ceiling or the minor edge below dmin=1.0
for every (r_m, R_M) — the link window ratio 2.07 is marginally too
narrow. Measured best-effort: 1–2 missing links, 29/48 parasites, min
gates 0.005/0.022, leak 13.9/19.3. The homology species needs a bigger
box or a finer foam (window ratio ≥ 2 with margin). Best-effort .nets
emitted with warnings.

**C7 truncated octahedron n=24 — the foam does not supply it.** 36
equal struts demanded; realized spread 1.04–1.97 (mean 1.61) after an
8-restart × 14-rotation search: min 0.035, mean 0.708, leak 21.0.
Stiffness gap 0.0116 — an almost-soft mode: the structure is not held.
Robustness pctile 0 (every jitter center did as badly or better —
uniformly infeasible, not a bad center). 11 parasites, all closed.

**C8 ring12 exact re-derivation — the validation candidate.** Kernel
picks reproduced (§3); `c8_ring12.net` re-emits comp12 through init=net
(θ₀=0 gauge; per-link lock recursion; E lines oriented around the
loop). Forward gates 1.0000/1.0000 by construction; back gates
gate(−2φ_k) sum 4.339 (mean 0.36 — long links leak backward much faster
than B1's uniform-d̄ 0.100 estimate; the shortest links are fully
one-way). 1 parasite (§4.2). Ctrl twin (mean-tuned BFS): predicted fwd
gates min 0.226 / mean 0.566 vs tuned 1.000/1.000 — the sharpest
exact-retuning A/B in the fleet; on a closed ring the mean-BFS control
spreads mistuning as per-link ω(d̄−d_k), no seam dump.

## 7. Scoreboard

`candidates/summary.tsv` (rebuilt from the per-candidate reports).
Ranking by predicted leak: **C2 3.23** < C5 4.02 < C8 4.98 < C4 5.56 <
C1 6.48 < C3 7.84 < C6a 13.9 < C6b 19.3 < C7 21.0. Per live direction:
C2 0.134 < C4 0.146 < C5 0.167. All loads inside the usable window;
C6 at its edge (flagged). Spectral gaps: healthy for C2 (6.42), C3
(3.90), C4 (1.52), C5 (1.42), C8; near-zero for C7 (0.012).

**Recommended next kernel runs** (all seeds ready): (1) `c1_cube125.net`
vs `c1_cube125_ctrl.net` vs the historical h1_shell0_v3 — what exact
retuning buys in lifetime on identical cells; (2) `c2_cube150.net` —
the cleanest consonant skin this foam supplies; (3) `c5_tube6.net` vs
comp12 — structure-vs-structure at matched mass; (4) `c8_ring12.net` vs
`c8_ring12_ctrl.net` — retuning isolated on the known-best object;
(5) `c4_ring12_chords.net` — the composite. NETGATE P lines should be
tracked over time for the lock-eligible parasites (C3's on-rung one
especially, if C3 is run at all).

## 8. Run protocol and caveats

* `python3 prestress/formfind.py --selftest` (validations), `--cand cN`,
  `--all`, `--foam PATH --out DIR`, `--jitter N`.
* The single-process `--all` sweep was interrupted at the session
  checkpoint before writing its combined summary; every artifact on
  disk is from the per-candidate runs (identical seeds — per-candidate
  results are bytewise deterministic and were verified stable across
  reruns). summary.tsv was rebuilt from the on-disk reports;
  validations re-executed at rebuild time: ring MATCH, shell MATCH.
* Jitter counts: 20 centers (C1–C5), 12 (C6), 6 (C7); C8 none (it is a
  fixed rebuild, not a search).
* 4-dp TSV rounding: picks and gates unaffected (margins ≫ 1e−4);
  link census differs by 2 borderline links out of 85886.
* The load x written to each .net is %.6f-rounded; the realized pitch
  differs by O(1e−6) — negligible against every gate scale here.
