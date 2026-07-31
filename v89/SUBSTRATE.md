# SUBSTRATE — evolving away from the foam (living lab notebook)

**Living document, opened 2026-07-31 on user directive:** the frozen
jittered foam was an improved stepping stone but will fail as a viable
means; replace the substrate MECHANISM, not the concepts — and the
replacement may need to precede EMF and charge work too, not only MASS.
**Success bar for each rung: the full battery passes on the new
substrate.** Subordinate to `PRINCIPLE.md`; ratchet rule governs;
`laws_V2g` remains THE standing table until a substrate variant earns
promotion by user decision.

## 1. Why the foam fails (measured, cross-program)

* Its link window [dmin, 1.15·(rᵢ+rⱼ)] ≈ [1.0, 1.96] is **19% length
  disorder** (σ_d = 0.29 over all channels, measured at build) — the
  frozen-foam ceiling: exact retuning saturates at mean gate 0.851
  (W1); shells need ≤2%.
* Species masses are foam accidents (~8% within-species vs 48%
  inter-m gaps — M-R1 impossible per-instance).
* Foam noise floors drown g2 free-fall and lensing (0.9 vs 0.015);
  visibility ceiling V≈0.46.
* Q9 measured the split directly: the Josephson slope is law-owned
  (foam-independent); the knee Δω_c and CPR amplitude are foam-owned
  (one link's accidents, ±30% expected link-to-link).
* Doctrinally: fixed random connectivity sits on the constraint-2
  checklist ("connectivity not produced by the current energy
  structure").

**What is preserved:** cells as processes, links as retardations,
gates/rungs/closure integers, conversion atoms, the mode books, integer
conservation. The substrate is apparatus; the laws are untouched.

## 2. The ladder

| rung | mechanism | status |
|---|---|---|
| **S1** | annealed frozen geometry (build-time, deterministic, RNG-stream-preserving) | **BUILT** — see §3; battery gate running |
| **S2** | coordinate-free diagnostic: explicit uniform-d graph, positions as scaffold only | design |
| **S3** | **livefab** — the dynamic complex: link/cell creation-annihilation as ledgered rewrite events; vacuum self-anneals; structures own their geometry | design track (kernel-2). §3's frustration measurement is its justification |

## 3. S1 — the annealed substrate (2026-07-31)

Three designs tried; two failed instructively:

1. **Dynamic-pair spring (REJECTED — densifies).** Springs toward a
   target length over all pairs in a generous cutoff: far pairs get
   pulled INTO the window; degree inflates 17.6 → 28.3, σ_d stays
   ~19–20%. A dynamic attractive pair set re-packs instead of
   annealing.
2. **Frozen-graph spring (REJECTED — frustration floor; the finding).**
   Relax the initial neighbor graph as a fixed elastic network toward
   uniform edge length: converges (200 ≈ 1500 sweeps) at **σ_d ≈ 18%**
   — barely below the initial 19.5%. A dense random graph (degree ~18,
   edges ≈ 9·NC vs 3·NC coordinates) cannot be equalized: **the
   disorder is topological, not positional.** This is the measured
   argument that full uniformity requires DYNAMIC topology — livefab
   (S3) — and that any frozen substrate must instead control the link
   WINDOW.
3. **Repulsion-to-contact + trimmed ceiling (ADOPTED).** Build the
   substrate like a jammed packing: dart scaffold at substrate density
   dmin=1.25 (NC 5039 at L=24; sets the monodisperse jamming contact
   at the legacy spacing) → pure repulsion pushes every pair up toward
   the contact floor (dynamic pair set is SAFE for pure repulsion — no
   densification; converged by ~800 sweeps; the ~10k residual
   below-floor pairs are the at-jamming equilibrium) → the first
   neighbor shell piles up sharply at contact → trim the channel
   ceiling to the shell (geom_lmax). Channels become contact-shell
   links.

**Kernel apparatus (all default-off = byte-identical legacy foam):**
`geom_relax` (sweeps), `geom_relax_k` (step), `geom_dtar` (contact
floor; <0 auto = 0.985·jamming estimate), `geom_lmax` (absolute
channel-ceiling trim), `geom_runi` (uniform radii; jitter draw still
consumed). The anneal runs BEFORE radii/normals/phases are drawn — the
RNG stream is untouched, so a relaxed run shares the legacy draw
sequence: only geometry moves. Serial, deterministic at any thread
count. Prints `# GEOM relax start/done` + `# GEOM links` stats.

**Calibrated result (L=24, seed 20260727; python prototype and C
kernel agree):**

| | legacy foam (V2g) | S1 annealed (V2s) |
|---|---|---|
| σ_d over channels | 0.29 (**19%**) | 0.046 (**3.0%**) |
| d̄ | 1.505 | 1.545 (m=1 rung window preserved) |
| NC (L=24) | 9741 | 5039 |
| mean degree | 17.6 | 7.3 |
| gate-deficit scale (∝σ_d²) | 1 | **~1/40** |

Chosen point: `dmin=1.25` (the one law-key VALUE changed — substrate
density; recorded as variant table **`laws_V2s.cfg`**, same key set as
V2g, purity-clean) + extras `geom_relax=800 geom_relax_k=0.3
geom_dtar=1.52 geom_lmax=1.70`. Watch item: **mean degree halves**
(17.6 → 7.3) — the largest physics-facing apparatus shift; conductance
sums, containment, and skirt behavior all ride on degree. The battery
is the judge, by design.

**Gate protocol:** (1) default-off regression battery on V2g (kernel
mod ⇒ ratchet); (2) full battery on V2s + anneal extras (`battery.py
--laws laws_V2s.cfg --extra ...`; harness gained `--extra`/`--tag`).
Verdicts appended below when the runs land. Bars are judged as-is
first; any red is classified (physics vs apparatus-window starvation,
e.g. pair-pick windows tuned to the old d spread) before any
adjustment is proposed.

## 4. Verdicts (2026-07-31, first S1 battery campaign)

**Regression gate:** V2g default (geom off) **21/21** — the kernel
apparatus is inert.

**V2s shell1** (floor 1.52 / lmax 1.70; d̄=1.545, σ_d=3.0%, degree
7.3): **15/21**. **V2s shell2** (floor 1.44 / lmax 1.62): **11/21** —
instructive failure, see finding F2.

| exp | shell1 | classification |
|---|---|---|
| e1, e3a, e4, e5, e7, d1, t1, q2, t4, qt_lo, qt_hi, p1, g1, LIN | PASS | — (e7 IMPROVES to 60/60 alive, frac 1.00, gg 0.96 vs 0.77 on foam; e2 v_g 0.88C vs 0.61) |
| e6 | FAIL → **probe PASS** (gg 0.21→0.886) | APPARATUS: pair window [1.0,1.55] + fixed x0=0.3 encode the old d-spread; widened window + auto x*(d) fixes it outright |
| e8 | FAIL (shed ladder degenerate, total 0.5 vs 7.2) | APPARATUS (same class): its δ0 ladder rode the foam's d variety; needs explicit δ0 seeding |
| e9 | FAIL (window [1.28,1.42] holds ZERO links) → **RESOLVED + corner probe PASS (locked@t10=19 vs the foam's 8)** | **TEST — the apparatus hard-codes the foam-era ADDRESS of the fifth, not the fifth.** Final analysis (three wrong classifications corrected by measurement, in order: σ_d/tongue-roughness — invalidated, the probe was mis-seeded by the unison-only auto mode, cellfab.c:1072; "habitable zone d≲1.50, shell excluded" — incomplete, the atlas shows TWO bands (m=2 [1.16,1.50], m=3 [1.75,2.25]) and the x≤0.9 edge is a seeder guard, not physics). The mutual-fifth closure (q·ω_i+p·ω_j)·d = 2πm at the shell's d̄=1.54 puts the m=2 rung at loads (x_i, x_j) = (0.351, 0.944) — the near-cap corner — and seeded THERE it locks 19/40 at t10, 2.4× the foam's best, t_last=60, bar PASS. σ_d=3% does NOT kill the fifth. Species load-addresses are substrate-dependent; existence is robust. Plasticity-rescue arm (foam-era loads + κ=1): partial (3 locks; links walk only ~0.01 by t=30) — the misfit force scales with transport, so malleability fine-tunes within the capture window but cannot cross the gate desert: a tuner, not an address-finder |
| e3b | FAIL (coherent BACKWARD drift, cos −0.956) → **RESOLVED + probe PASS** | TEST OPERATING POINT on a real physical structure: at amp=0.5 the loaded pitch gives ω·d̄ ≈ π AND kx·d̄ ≈ π — the forward travel gate and the backward wrap-rung (d(kx+ω)=2π) are DEGENERATE; the foam's 19% smear broke the tie incoherently (V2g's own log shows net cmx 7.005→6.657 BACKWARD; its PASS is a second-half drift), while the uniform shell resolves it coherently backward. Probe at the derived travel window (amp=0.25, kx=2.55: ω(x)≈kx, 2ωd̄ far from 2π): **PASS cos=+0.894, speed 0.00805 — 2.3× the foam's best.** Physics bonus: the emergent-inertia wall is mechanism-identified — heavy blobs freeze/reverse because loaded pitch × d̄ approaches the backward rung π |
| g3 | FAIL (deflection +1.17 vs +2.5) | G-SECTOR RECALIBRATION DEBT: NC halved ⇒ blob source mass ~halved (amp/sigma are per-cell), footprint shallower (core/far 0.73 vs 0.53); degree probe (deg~11) did NOT recover it ⇒ source-mass scale, not connectivity |
| g4 | FAIL (fflux_late=0; decay ratio borderline) | PHYSICS SHIFT, arguably favorable: on the quiet substrate the blob leak stops radiating (roughness ↓ — the same effect that lengthens structure lifetimes) — the bar encoded the foam's rough-leak signature |

### Findings of the rung (each a design law for every later rung)

* **F1 — the disorder is topological.** Frozen-graph relaxation stalls
  at σ_d ≈ 18% regardless of sweeps: a dense random graph cannot be
  equalized positionally. Uniformity requires either window control
  (S1's contact shell) or dynamic topology (S3). This is livefab's
  measured justification.
* **F2 — the contact floor must sit AT jamming.** Below it (shell2)
  the first shell goes diffuse, the trimmed window starves (degree
  7.3 → 3.4) and the OPTICS sector collapses (d1 V 0.39→0.09, e2
  0.88→0.56C, HOM split gone): **interference visibility and pulse
  speed are the substrate's connectivity ruler.** Healthy degree ≥ ~7.
* **F3 (final form) — the substrate's d̄ sets each species' ADDRESS,
  not its existence.** Every mutual interval lock lives on closure
  bands in d (the atlas: Ω·d = πm with the pitch window set by skirt/
  cap/w2 — unison [1.16,2.25]; fifth [1.16,1.50]∪[1.75,2.25]; fourth
  [1.16,1.69]∪[1.55,2.25]; octave: none anywhere — the measured
  lifetime hierarchy's tail, explained structurally). Changing d̄ moves
  each interval's rung to different LOADS (the fifth at the shell
  lives at x_j≈0.94 and locks 2.4× better than on the foam);
  foam-era apparatus froze the old addresses. Standard tool going
  forward: compute the atlas + per-d addresses at config time.
  (The Γ/(pq) tongue-vs-σ_d bound is answered empirically at 3%:
  fifths lock fine; the ≲0.4% worry applied to a mis-analysis.)
* **F4 — uniform substrates resolve gate degeneracies coherently.**
  Operating points tuned on a disorder-smeared substrate can flip sign
  when the smear is removed (e3b); conversely the quiet substrate makes
  real structures sharp: the inertia wall = the ω·d̄ → π backward-rung
  approach, now mechanism-identified and probe-confirmed (light blob in
  the travel window: cos +0.894 at 2.3× foam speed).
* Constraint set for any future substrate: floor at jamming; degree
  ≥ 7 (optics green); σ_d per the interval content needed; d̄ inside
  the loaded-pitch gate window; per-experiment apparatus constants
  that encoded foam statistics (kx, pair windows, δ0 ladders, blob
  amp per NC) must be re-derived per substrate, bars unchanged.

### Standing state and next steps

`laws_V2g` remains THE table; **V2s shell1 is the S1 working variant**.
With the four probe-green fixes (e3b travel window, e6 windows+auto,
e9 shell-address loads, plus e8's explicit-δ0 ladder designed), the
substrate stands at an effective **18/21** with only g3 (source-mass
sigma rescale, fix computed) and g4 (check re-expression, favorable
physics) remaining. Next: (1) adopt the variant apparatus set — the
general method fix: derive every experiment's operating point from its
closure relations + the substrate's d̄ at config time (the atlas),
never freeze substrate-dependent coordinates; add an interval-aware
auto-seed mode (ω_i = πm·p/((p+q)·q·d)-class formula with reachability)
— ratchet-gated seeder feature; (2) g3 sigma×1.246, g4 substrate-
invariant check form (user sign-off); (3) full battery rerun on the
variant apparatus set → the 21/21-without-the-foam claim; (4) livefab
(S3) design carrying F1–F4; malleability's measured role: local tuner
within the capture window (the plasticity arm), addresses still set by
closure arithmetic.