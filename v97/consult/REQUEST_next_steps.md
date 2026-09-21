# v97 consult — validate the fresh-eyes evaluation, rate the next steps, add your own (2026-08-13)

You are resumed with your prior session context. A fresh independent
review (claude, this session) just read v89→v97 end to end and produced
the findings and next-step plan below. Your three jobs:

1. **Validate or refute the findings** — especially the math in §2.
2. **Rate/critique the proposed next steps** (§4) — ordering, what to
   cut, what is missing, where the plan fools itself.
3. **Contribute your own ideas** for making H1 inhabited — especially
   for pinning J (latch vs alternatives) and for coherent intake.

Scope note from the user: **C-kernel parity is out of scope for now.
The Go kernel (`v97/fab/`) is the working substrate.** Config-only and
Go-side work is the arena; kernel *faces* still need explicit user
sign-off. Ground yourself in the actual files (repo root
`/home/d/code/scp`): `v97/{README,PRINCIPLES,METRICS,H0,H1,LEADS}.md`,
`v97/fab/promo.go`, `v97/hunt/convert.go`, `v97/runs/leads/H1.out`,
`v93/theory/qgraph/docs/UNIFIED_LAW.md`.

---

## 1. Status recap (measured, adopted)

- v97 = THE REGULAR METHOD: six principles adopted; scoreboard = the
  five-clause conversion-dipole hunt (I/O separated; Φ>0; dM/dt≈0;
  dipole of J steady; Es throughput at constant M). Battery 87 GREEN
  at defaults; carried knobs all default-0.
- **H0 PASS** (meter honesty, 5 classes; one named SKIP: snapshot J is
  not `fd` — on e3b the snapshot TJnet sloshes −2±14 while integrated
  `fd` is seed-robustly +x).
- **H1 EMPTY**: 0/203 archive streams inhabited; best 3/5. I/O
  spatial separation appeared in 2/106; never together with an Es
  cycle at constant M.
- **Promo leads** (`promo=0/1/2/3`, Go-only oriented door):
  orientation is real (polar 0.14→0.29–0.60), the space-leg is real
  (`promo=2` makes sep + Es-at-constant-M co-occur — the archive never
  had that pair). H1 still EMPTY on two clauses:
  Φ/Em = 8e-4–3.6e-3 vs bar 0.015 (~10× short), and Jdip dir-stable
  0.54–0.69 vs bar 0.80. Also measured: I/O *imbalance* both ways
  (control ⟨I⟩=0.58/⟨O⟩=0.34; hot arm I=0.27/O=0.59), p2_ring dark →
  door never fired (1/5), p2_jm melted (j_metric amplifies leaks too),
  p3_close starved cond (0.17 vs 0.53).

## 2. The fresh review's findings (validate these)

**F1 — The two failing clauses ride the two weakest meters.**
(a) Φ = min(I,O) is graded by the frame-diff mode-crossing heuristic
(`hunt/convert.go`): *net* Em/Ee change per matched cell per ~4 t.u.
frame. A cell that both condenses and evaporates within a frame nets
toward zero → the heuristic systematically under-counts **gross** door
traffic — the one error direction H0 could not catch (it only verified
known-empty streams read empty; the synthetic dipole is heuristic-
shaped). Meanwhile the **exact per-cell door ledger already exists and
is already in the FCS**: `pI`/`pO` accumulate in the Go door write
path (`fab/promo.go`: cond→pI at dest; evap+radiance→pO) and are
written every snap (`fab/fcs.go:183`) — and `hunt/` never reads them.
(b) The J clause uses snapshot bilinear J — the same signal H0
SKIPped as unreliable for direction. Caveat noted: `promoPick` also
*steers by* that snapshot compass, so part of the slosh is real door
behavior, not just meter noise.
The archive-wide EMPTY survives this critique (net-based I/O still
detects genuine polarization; sep was 2/106). The *quantitative* 10×
Φ shortfall on the leads arms does not yet deserve trust.

**F2 — The Φ/Em bar is provably unpassable at g=0 (check this math).**
At the radiance balance point of the adopted law
(O = k_rad·cap·x^p, Em = x·cap, p=4, k_rad=0.05, cap=2.5, xb=0.62):

    Φ/Em = k_rad·x*^(p−1) = 0.05·x*³.

Bar 0.015 ⇒ x* ≥ (0.3)^(1/3) = 0.669 > xb = 0.62 — the motif must
hold load **above ambient**. qgraph Theorem 3.4 (box-fill) proves
that is impossible at g=0 (ceiling 0.05·0.62³ = 0.0119 < 0.015);
Theorem 3.5 (u⁴ − g·u³ − 1 = 0) gives the threshold:
u* = 0.669/0.62 = 1.08 ⇒ **g ≳ 0.29**. So H1's Φ clause is, in
disguise, the statement "the motif has phase-matched (stimulated)
intake with g ≳ 0.3." Consequences if the math holds:
- The promo leads' Φ failure was **predetermined** — every arm ran in
  dark/random baths (g=0). No door geometry can fix Φ.
- The charter's "first face = slot door" covers only half the gap:
  slot door targets sep + J-lock; only coherent intake targets Φ.
  (rad_dual as coded is spontaneous-only, I = k·rad_dual·cap·xb^p —
  consistent with p2_dual not lifting Φ.)

**F3 — Surface/volume favors small motifs.** Φ is surface-fed, Em is
volume: the ~100-Em leads blobs dilute Φ/Em structurally. The C6 ring
(Em ≈ 5; observed cond rates ~0.5 where doors fire) is quantitatively
the right size for the bar.

**F4 — promo=2 already implements P5's space-leg as a FULL cycle**
(Es consumed at the fore cell in `applyPromoCond`, back-share returned
at the aft cell in `applyPromoOut`) — answering v96's "W-A is half a
cycle" objection at the event level.

## 3. Open question for you (Q4 below depends on it)

Can the **existing** resonance machinery (G_res pitch matching, comb
gate, workfn) deliver effective g ≳ 0.3 from a *pitch-matched
organized field* (standing wave `sw_kx/sw_amp` tuned to the ring's
w₂e), or is a kernel g-term (the stimulated cross term
k·cap·g·x^(p−1)·xb of qgraph L4/Thm 3.5, never coded) unavoidable?
Note the negative priors: v92 closed three feedback forms; P5 measured
amp_door+**noise** ≈ null (ratio ~1.00) — but noise is g=0 by
definition, so that null does not test an organized field.

## 4. Proposed next steps (rate these; propose ordering changes)

**Step 1 — Sharpen the two meters, re-grade (offline `hunt/` only).**
Teach `hunt/convert.go` to use exact `pI`/`pO` columns when the FCS
carries them (all post-leads Go runs do), heuristic fallback for the
archive. Add an H0 class "exact vs heuristic on a leads stream." Then
re-grade the leads arms. Optional (Go `fab/` edit, physics-inert,
needs user OK): an integrated directed-hop vector accumulator `pJ`
alongside pI/pO, retiring the snapshot-J SKIP for good.

**Step 2 — The lit-bath helix ring (config-only; the registered-but-
never-run arm from RING_DNLS A.5).** `exp=ring nv=6 seed_mw=2` + kx
tilt (helix: polarized cycle so Jdip isn't a pure loop) + laws_V3a +
`amp_nat=2 amp_logate=1 hop_order=1` + pitch-matched standing wave
(`sw_kx/sw_amp` at the ring's w₂e; `bath_xb` as apparatus fallback;
NOT noise_amp) + `amp_door=1` + `promo=2` (+ optionally rad_dual=1).
T ≥ 1000, 3 seeds, grade with sharpened `-h1`. Unites: the only object
that holds phase (cycle), the only condition where doors fire (lit
bath), the only intake that can beat box-fill (organized field through
the resonance gate), the only geometry that splits I/O (oriented door
+ space-leg). p2_ring failed *only* because the vacuum was dark.

**Step 3 — Two cheap theory derivations (qgraph-style, symbolic+Go).**
(a) Formalize F2 as a theorem (Φ/Em bar ⇔ g ≳ 0.3).
(b) **Edge-door selection rule**: qgraph L5 computed a vertex door has
flat support over Δm (mean 0, sd 3.49 on C12 — the measured winding
scrambler) and a wound door pumps exactly Δm=k. Compute M_ab for a
door supported on one *oriented edge* of C_N; prediction: support
concentrated at Δm=±1 (dipole selection rule). Pre-registers the
slot-door face's falsifiable bar.
(c, optional) Closure-rule graph-selection toy: iterate "promote iff a
joint cycle closes" on random contact graphs; measure survivor b₁ /
magic spectrum. Attacks qgraph's biggest named gap ("the graph is
given, not selected") = P3's core claim.

**Step 4 — Only if Step 2 fails with honest meters: the slot-door
latch face (Go kernel, user sign-off).** LEADS' own closing candidate:
the link remembers its last successful promotion ("this link is the
machine"); intake fires at the fore end and outtake at the aft end of
the SAME latched pair. Fixes Φ structurally (I/O fire as a matched
cycle, not by statistical accident) and dir-stability (the door PINS
the direction instead of re-picking max-J from a sloshing snapshot
compass every event — see `promoPick`/`cellDir` in `fab/promo.go`).
Pedigree: the ASYM bed law (conductance follows the signed net)
applied to the door's *locus* — third application of an adopted law
family; unlike j_metric it biases only where conversion fires, not
transport rates (j_metric melted p2_jm because it amplifies leaks).

## 5. Questions (answer each explicitly)

- **Q1.** Is F2's derivation sound? Check: per-cell vs per-motif
  aggregation, the surface/volume correction (does it make the bar
  *harder* than g≈0.3?), the `comp = w2/w2e` asymmetry in rad_dual
  (M3), and whether Φ = min(I,O) at balance equals O.
- **Q2.** Is F1's meter critique right — can net-vs-gross
  under-counting plausibly account for the full 10× Φ shortfall, or
  only part of it? Would you re-grade before or after running Step 2?
- **Q3.** Rate the Step 1→4 ordering. What would you cut, add, or
  reorder? Is anything here self-deception (bars softened in disguise,
  meters tuned toward passing)?
- **Q4.** §3: existing resonance machinery vs a new g-term — which do
  you predict, and what config arm would decide it fastest?
- **Q5.** Your own ideas: (a) alternatives to the one-bit latch for
  pinning J (e.g. signed vector credit from ASYM — the refinement v95
  explicitly did NOT take; phase-locked fore/aft pair; something
  else); (b) any different route to Φ ≥ bar that the plan missed;
  (c) anything in the v97 charter you'd push back on now.

Write your full review; be adversarial where the plan deserves it.
