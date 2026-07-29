# MORPHO — the structural morphological search (catalog + backtracer + exotics)

*2026-07-29, morpho agent (stream C of the PRESTRESS campaign, LEDGER.md).
Everything here was computed on the REAL foam `foam/foam_s20260727.tsv`
(9741 cells, L=24, 85884 links, d̄=1.505, link window [1.000, 2.062],
mean foam degree 17.6) with `search.py` — a lightweight scorer built
independently of `formfind.py` (deliberate double implementation; the two
should cross-check each other's picks). Re-run: `python3 search.py`
(≈15–30 min; `--quick` for a 1–2 min smoke). Outputs: `ledger.tsv`
(every candidate, one line), `shortlist.json` (top 8 + 2 negative
controls, concrete cells/edges/phases), `nets/*.net` (kernel-ready
seeds, frozen V/E format), `report.json` (full detail).*

Physics verified against the kernel before anything was scored:
`gate_of` (cellfab.c:489) = ((1+cos ψ)/2)^8; gates ψ_f = θ_i − ω_i d/C − θ_j
(:1811); pitch w2e = w2/(1+q·(Em+flload)/cap) (:1715); link rule
d < 1.15(r_i+r_j) (:701); p:q comb tongues (:1787); the `init=net` NETGATE
+ parasite report (:1286–1326); QATOM credit machinery (:1558). Laws:
`battery/laws_V2g.cfg` (w2=2.9, q=1.2, cap=2.5, p_gate=8, γ_m=0.10).

---

## 1. The scorer (independent implementation)

Given a cell set + intended edge list embedded on the foam, at uniform
pitch ω (static lock ⇒ one ω per locked component):

* **Phases** are solved by the kernel's own lock relaxation — each edge
  pulls its receiver toward forward alignment at a rate weighted by the
  open gate (the rate-level fixed point of `kappa_lock`, cellfab.c:2034),
  from a spanning-tree init that fixes the winding basin. Dead gates pull
  nothing, so winding is protected in the scorer exactly as in the sim
  (the desert). Design basins (ring m) are pinned by oriented target
  retardations, ring_m-semantics: ω = 2πm/L_actual.
* **Leak model** (a RANKING device; the sim is the judge):
  intended strut (1−gf)+(1−gb); intended cable (1−g_fwd)+0.5·g_back
  (back-gate openness = desert-protected frustration, weight from the
  comp12 anchor: gb=0.100 and longest-lived); parasite within a locked
  component gf+gb at the solved phases; parasite across components
  2·⟨G⟩ = 0.39 (relative gauge drifts; ⟨G⟩ = C(16,8)/2¹⁶ = 0.196);
  vacuum port res₁₁(ω)·⟨G⟩·0.11 per unlisted foam neighbour, with
  res₁₁ = γ_m²/(γ_m²+(w2−ω)²) — the skirt trickle. The vacuum term
  falls with load, reproducing the measured lifetime-rises-with-load law.
* **Structural penalties**: open vertex (intended degree < 2) 0.40;
  bridge edge 0.20; gate-flux Kirchhoff imbalance (per-vertex open-gate
  inflow ≠ outflow ⇒ congestion, unlock) weight 1.0. Labels
  strut/cable/direction are DERIVED per edge (argmin of the three
  costs), so tensegrity splits of over-constrained shells are found
  automatically rather than enumerated by hand.
* **Validity domain**: anchored by measured lifetimes at per-voice
  x ∈ [0.128, 0.44] (skirt ladder, MASS.md). The ω scan is clamped to
  x ≤ 0.62; anything scored beyond is flagged `x_unmeasured` and kept
  out of the shortlist (one deliberate heavy probe is in the ledger).

**Calibration against the measured anchors** (leak/voice, foam-realized):
blob 0.30–0.36 (worst measured class ✓), open chain penalized for its 2
open ends ✓, cube_a1.25 ≈ 0.3 with ~10 parasites (died 1.2–1.8k ✓),
ring12_m6 struts > ring12_m5 wound (2221 vs alive-at-5000 ✓), paper
ring12_m5 = 0.0501/voice = exactly the designed 0.5·G(60°) back gate ✓.
The one anchor tension: pure leak ranks comp6 (one-way, tiny back gate)
above comp12, while comp12 outlived comp6 — resolved by the parasite
discovery in §5 (comp6 carried ~5 hot second-neighbour parasites), which
the scorer sees and the 2026-07-28 gate reports did not.

---

## 2. The morphological rulebook for THIS foam (derived, then verified)

Everything below follows from the kernel facts + the measured foam
statistics; each was checked numerically in `search.py`.

* **R1 — one pitch, one strut length.** Static lock ⇒ uniform ω per
  component; both-gate links need 2ωd = 2πm with only m=1 reachable
  (m=2 ⇒ d = 2π/ω ≥ 2.8 > ceiling) ⇒ ALL struts share d = π/ω.
  Seedable strut lengths: d ∈ [1.41, 1.96] (ω ∈ [1.60, 2.23],
  x ∈ [0.25, 0.68]). The foam's own mode d̄ = 1.505 sits at ω = 2.087,
  x = 0.325 — the comp12-class load. The laws engineered this seam;
  morphology should exploit it: **d* = 1.5 is the strut sweet spot.**
* **R2 — parity.** Strut cycles close automatically iff even (Σψ = 2kπ);
  odd cycles are frustrated by π. Struts-only closed structures must be
  BIPARTITE. Triangulated shells (octahedron, icosahedron), odd rings,
  antiprisms, Möbius strips (even-k) are dead as strut objects.
* **R3 — cables are the length-freedom.** Forward-only links take any
  d with φ = ωd; each independent cycle needs its oriented sum on the
  ladder: ω·L_cycle = 2πm. Back-gate leak G(−2φ) is a designed constant.
  Odd N is FINE for cable rings (parity is a strut constraint only) —
  the catalog's ring9_m3 exists because of this.
* **R4 — the back-gate ladder (wound ring family).** For a uniform ring,
  φ = 2πm/N and back gate = G(2π(N−2m)/N) = G(4πw/N), w = N/2 − m:

  | ring | φ | d at rung | x | back gate | paper leak/voice |
  |---|---|---|---|---|---|
  | ring12 m=5 (comp12) | 150° | 1.25 | 0.32 | 0.1001 | 0.0501 |
  | ring10 m=4 | 144° | 1.35 | 0.47 | 0.0263 | 0.0168 |
  | ring8 m=3 | 135° | 1.30 | 0.50 | 0.0039 | 0.0019 |
  | ring9 m=3 | 120° | 1.15 | 0.50 | 1.5e-5 | ~0.0000 |
  | ring8 m=2 | 90° | 1.126 | 0.90 | 0 (G(π)) | 0.0000 |

  Every step DOWN in φ from comp12 kills the back gate exponentially
  (gate = cos¹⁶(ψ/2)) while pushing the load UP — and the measured death
  law says heavier lives longer. φ = 90° is the perfect one-way point
  (back gate exactly 0) but needs foam-floor links (d = 1.13 at x = 0.9,
  unmeasured-heavy, headroom 0.1 ⇒ congestion risk — ledger probe only).
  **ring8_m3 / ring10_m4 / ring9_m3 are comp12's natural successors.**
* **R5 — parasite clearance is a LENGTH law.** Any unintended pair
  closer than 1.15(r_i+r_j) (up to 2.07 for the largest radii) is an
  off-rung channel. Measured on the foam radii:
  P(pair cut > 1.905) = 0.835, P(cut > 2.12) = 0.000. Hence:
  quad faces need edge d ≥ 1.46 (diagonal d√2 > 2.07 — at d=1.5 the
  diagonal 2.12 is NEVER a link); hexagon/octagon second-neighbours are
  safe at any strut length; ring second-neighbour chords need
  2R·sin(2π/N) > 2.07. Rhombic faces are dead (short diagonal 1.15·d
  always inside the ceiling ⇒ 12 built-in radiators for the rhombic
  dodecahedron). Double shells (cube-in-cube) are dead: one strut length
  forces shell separation = d, putting O(n²) inter-shell pairs inside
  the ceiling — and the outer-shell edges beyond it.
* **R6 — collimation (the anti-polyhedron law).** A cell's two planes
  intersect in ONE axis, and a channel carries only along directions in
  both planes (cellfab.c:1746-56): conductance factor
  (1−(n̂₁·û)²)(1−(n̂₂·û)²) per end. A degree-2 voice on a gentle arc
  aligns its axis with both links (ring12: 15° bend ⇒ factor 0.93).
  A cubic vertex cannot: best case (1, ¼, ¼) across its three edges;
  a quad-mesh (deg-4) vertex is ≤ ¼ on a whole ring family. **1D curves
  outconduct 2D skins per link — rings and gently-curved tubes are the
  foam's preferred morphology; compact polyhedra fight the collimation
  law.** (Consistent with measurement: every cube died young while
  rings persist; gates mean 0.6 was only part of the story.)
* **R7 — load & skirt.** x ≥ 0.25 seeded (skirt death at 0.0617);
  lifetime ordering by per-voice load is measured; validity of the
  scorer capped at x = 0.62 (heaviest measured survivor class 0.44).
  Mass = n·x·cap.
* **R8 — Kirchhoff balance.** Steady circulation needs per-vertex
  open-gate inflow = outflow. Rate-level phase solutions can otherwise
  settle into folded source/sink patterns (leak-free on paper, congested
  in the sim — the receiver caps at head→0 and the lock decays). The
  scorer penalizes imbalance; junctions (theta graphs, books) must route
  circulation consistently through their arms.

---

## 3. Sweep 1 — the catalog (results on the real foam)

**Run fidelity note (2026-07-29):** the numbers below are the completed
QUICK pass (12 embed tries/candidate — 36 for n>16, 2 restarts × 4 pitch
classes × 250 annealing steps in the free search; RNG seed 20260728).
A deeper pass (`python3 search.py --tries 60 --restarts 6 --steps 1400`)
was launched and terminated at session checkpoint before writing output;
realizability and σ_d columns are therefore coarse (±1 embed), leak
rankings are stable across both quick passes run. Paper columns are
exact (foam-independent).

Columns: realiz = fraction of random placements where every intended
edge is a real foam link; σ_d = median edge-length spread of successful
embeds; leak/voice (e/p/v) = intended-edge / parasite / vacuum-port
terms; s/c = auto-derived strut/cable split (the tensegrity demotion);
wind = Σ|m − N_c/2| (chirality content); imbal = Kirchhoff gate-flux
imbalance; paper = leak/voice on ideal geometry; `x!` = outside the
measured-validity load window.

| candidate | class | n | parity | realiz | sigma_d | omega | x | leak/voice (e/p/v) | npar | s/c | wind | imbal | paper | verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| blob12 | control | 12 | bip | 1.0 | - | 1.938 | 0.414 | 0.298 (0.00/3.53/0.05) | 9 | 0/0 | 0.0 | 0.0 | - | worst class (calibration ✓) |
| chain12 | control | 12 | bip | 0.25 | 0.238 | 1.992 | 0.380 | 0.183 (2.14/0.00/0.05) | 0 | 6/5 | 0.0 | 0.10 | 0.000 | open control (2 skirt ports) |
| ring12_m6 | anchor | 12 | bip | 0.42 | 0.215 | 1.935 | 0.416 | 0.169 (1.96/0.03/0.05) | 1 | 5/7 | 0.0 | 0.12 | 0.000 | viable (strut ring; died 2221) |
| ring12_m5 | anchor | 12 | bip | 1.0 | 0.195 | 1.829 | 0.488 | 0.064 (0.71/0.03/0.04) | 2 | 2/10 | 1.0 | 0.08 | 0.044 | **viable — the standing champion** |
| ring6_m2 | anchor | 6 | bip | 1.0 | 0.162 | 1.840 | 0.480 | 0.256 (0.11/1.40/0.02) | 3 | 0/6 | 1.0 | 0.00 | 0.376 | **parasite-frustrated** (§5: hot 2nd-nbr chords) |
| cube_a1.25 | anchor | 8 | bip | 0.92 | 0.235 | 2.133 | 0.300 | 0.241 (1.37/0.51/0.04) | 8 | 6/6 | 0.0 | 0.18 | 0.073 | parasite-dense (died ~1.2–1.8k ✓) |
| cube_a1.5 | shell | 8 | bip | 0.83 | 0.228 | 2.028 | 0.359 | 0.269 (2.11/0.00/0.04) | 5 | 3/9 | 0.0 | 0.14 | 0.000 | poor realization (paper-clean) |
| tower2_a1.5 | shell | 12 | bip | 0.33 | 0.262 | 2.011 | 0.369 | 0.501 (5.15/0.82/0.04) | 14 | 8/12 | 0.0 | 0.12 | 0.000 | parasite-dense |
| tower3_a1.5 | shell | 16 | bip | 0.58 | 0.245 | 2.009 | 0.370 | 0.441 (5.04/1.96/0.06) | 21 | 12/16 | 0.0 | 0.15 | 0.000 | parasite-dense |
| hexprism_d1.5 | shell | 12 | bip | 0.42 | 0.244 | 2.176 | 0.277 | 0.321 (2.98/0.80/0.08) | 10 | 7/11 | 0.0 | 0.15 | 0.000 | parasite-dense |
| hexprism_d1.6 | shell | 12 | bip | 0.08 | 0.252 | 2.070 | 0.334 | 0.460 (3.46/2.00/0.05) | 12 | 4/14 | 0.0 | 0.26 | 0.000 | parasite-dense, rare |
| octprism_d1.5 | shell | 16 | bip | 0.42 | 0.258 | 2.105 | 0.315 | 0.270 (3.80/0.43/0.09) | 9 | 4/20 | 0.0 | 0.18 | 0.000 | parasite-dense |
| truncoct_d1.5 | shell | 16 | bip | 0.67 | 0.222 | 1.981 | 0.387 | 0.246 (3.86/0.01/0.06) | 7 | 9/11 | 0.0 | 0.17 | 0.000 | best shell; still 4× worse than rings |
| rhombdodec_d1.5 | shell | 14 | bip | 0.08 | 0.244 | 2.096 | 0.319 | 0.647 (5.50/3.50/0.06) | 21 | 8/16 | 0.0 | 0.17 | 0.000 | **dead: rhombus diagonal 1.73 = built-in radiators (R5 confirmed)** |
| cubeincube_a1.5 | shell | 16 | bip | 0.00 | - | - | - | - | - | - | - | - | - | **NOT SUPPLIED** (outer edges 3.23 > ceiling, as derived) |
| tube6x3_d1.5 | tube | 18 | bip | 0.33 | 0.240 | 2.084 | 0.326 | 0.430 (7.01/0.65/0.08) | 15 | 8/22 | 0.0 | 0.14 | 0.000 | parasite-dense |
| tube6x3_d1.6 | tube | 18 | bip | 0.17 | 0.244 | 2.022 | 0.362 | 0.382 (4.93/1.87/0.07) | 13 | 14/16 | 0.0 | 0.12 | 0.000 | parasite-dense |
| tube8x3_d1.5 | tube | 24 | bip | 0.11 | 0.246 | 2.168 | 0.281 | 0.365 (7.74/0.88/0.13) | 20 | 19/21 | 0.0 | 0.15 | 0.000 | parasite-dense |
| tube6x4_d1.5 | tube | 24 | bip | 0.17 | 0.243 | 1.973 | 0.392 | 0.418 (9.56/0.39/0.08) | 26 | 16/26 | 0.0 | 0.21 | 0.000 | parasite-dense |
| woundtube12x2_w1 | tube | 24 | bip | 0.06 | 0.225 | 1.970 | 0.393 | 0.428 (9.40/0.79/0.08) | 16 | 13/23 | 3.0 | 0.22 | **0.057** | paper-excellent, foam-rare (see note below) |
| theta_666_d1.5 | theta | 17 | bip | 0.14 | 0.192 | 1.952 | 0.405 | 0.193 (3.21/0.00/0.07) | 3 | 7/11 | 1.0 | 0.19 | 0.000 | viable (strained junctions) |
| theta_866_d1.5 | theta | 19 | bip | 0.19 | 0.225 | 2.111 | 0.311 | 0.184 (3.36/0.01/0.11) | 1 | 13/7 | 0.0 | 0.16 | 0.000 | viable (cycle integers 7/6 realized) |
| octahedron_d1.5 | control | 6 | ODD | 0.50 | 0.257 | 2.030 | 0.357 | 0.721 (4.30/0.00/0.03) | 1 | 0/12 | 2.0 | 0.26 | 0.667 | **FRUSTRATED (parity) ✓ control** |
| icosahedron_d1.5 | control | 12 | ODD | 0.17 | 0.267 | 2.157 | 0.287 | 1.175 (13.2/0.82/0.07) | 10 | 2/28 | 6.0 | 0.33 | 0.917 | **FRUSTRATED (parity) ✓ control** |
| antiprism4_d1.5 | control | 8 | ODD | 0.42 | 0.259 | 2.094 | 0.321 | 0.854 (6.37/0.43/0.04) | 5 | 5/11 | 3.0 | 0.17 | 1.000 | **FRUSTRATED ✓** (slant rung d=0.75 unreachable) |
| mobius6_d1.5 | control | 12 | ODD | 0.17 | 0.262 | 2.044 | 0.349 | 0.373 (4.14/0.28/0.05) | 8 | 2/16 | 1.0 | 0.19 | n/a* | **FRUSTRATED ✓ control** (*ideal ribbon violates link window at n=12: doubly impossible) |
| ring9_m3 | ring | 9 | odd-OK | 0.92 | 0.191 | 1.674 | 0.610 | 0.107 (0.05/0.89/0.02) | 2 | 0/9 | 1.5 | 0.00 | **0.000** | **viable — parasite-free one-way (the comp6 discriminator, §5)** |
| ring8_m3 | ring | 8 | bip | 0.75 | 0.171 | 1.736 | 0.558 | 0.058 (0.14/0.30/0.02) | 7 | 0/8 | 1.0 | 0.02 | **0.002** | **viable — back gate 26× below comp12** |
| ring10_m4 | ring | 10 | bip | 0.92 | 0.177 | 1.704 | 0.585 | 0.067 (0.56/0.08/0.03) | 1 | 2/8 | 1.0 | 0.11 | **0.016** | **viable — back gate 4× below comp12** |
| ring16_m8 | ring | 16 | bip | 0.08 | 0.211 | 2.059 | 0.340 | 0.183 (2.84/0.01/0.08) | 1 | 7/9 | 1.0 | 0.21 | 0.000 | viable, foam-rare at n=16 |
| ring8_m2_x0.9 | ring | 8 | bip | 1.0 | 0.170 | 1.394 | 0.900! | 0.169 (0.26/1.08/0.01) | 1 | 0/8 | 2.0 | 0.00 | 0.000 | heavy probe, x unmeasured — ledger only |
| ring13_strut | control | 13 | ODD | 0.42 | 0.192 | 2.157 | 0.287 | 0.103 (1.25/0.00/0.09) | 0 | 6/7 | 0.5 | 0.16 | 0.077 | frustrated seam control (m=6.5 impossible; leak = the π defect spread) |
| hopf12x12_m5 | exotic | 24 | bip | 0.67 | 0.198 | 1.903 | 0.437 | 0.159 (1.30/2.43/0.08) | 7 | 6/18 | 2.0 | 0.18 | 0.044 | viable — zero cross-channels verified on foam |

**Class verdicts.** (1) *Every 2D skin and tube realizes parasite-dense
on this foam* — quick-pass embeds land 5–26 unintended in-ceiling links
(σ_d ≈ 0.22–0.26 forces vertices together), even for shells that are
paper-perfect (cube_a1.5, truncated octahedron, all tubes, the wound
tube: paper 0.000–0.057). The catalog's paper/foam split is the central
result: **the surface topologies are consonant in principle and
foam-starved in practice at n ≤ 32**; they need either annealed growth
or a larger foam patch, not better picks. Truncated octahedron is the
best-realizing shell (0.67, npar 7) and still scores 4× worse than
rings. (2) *The wound-ring family realizes cleanly* (realiz 0.75–1.0,
σ_d ≈ 0.17–0.20) and owns the top of the table. (3) Parity controls
frustrate exactly as predicted (paper 0.67–1.0 leak/voice ≈ dead-gate
territory). (4) ring13_strut shows the odd-seam defect π spreading into
13 slightly-lifted gates (leak 0.103, npar 0) — frustration without
parasites, the clean signature.

## 4. Sweep 2 — the free search (backtracer)

**Closure criterion (justified):** a state is CLOSED when every vertex
has intended degree ≥ 2 AND no intended edge is a bridge (2-edge-connected:
every edge lies on a cycle). Physical grounds: an end voice's unpaired
gate faces the room = a skirt-trickle port (measured: chains leak 1.7×
faster than their closed rings; full closure worth 4× vs the blob,
MASS.md). Degree ≥ 3 ("skin") is NOT required for closure — the
longest-lived measured object (comp12) is degree-2 throughout — and the
collimation law (R6) says skins pay a per-link conductance tax that
curves do not. So the backtracer requires 2-edge-connectivity, not
skins; skins compete on equal terms and lose.

Method: simulated annealing directly over foam subgraphs (states = 8–32
connected cells + intended-edge subset of their real links; moves =
add/remove cell, toggle link intended↔parasite, close-an-end; objective
= leak/voice + penalties; every internal link is always either intended
or scored as a parasite, so nothing is unaccounted). Restarts across
foam regions × pitch classes ω ∈ {2.094, 1.963, 2.231, 1.800}
(strut d* 1.50/1.60/1.41/1.75), plus restarts seeded from the best
catalog embeds. Top states re-polished with a free ω scan and a greedy
closure pass.

**How far it got:** quick fidelity — 2 restarts × 4 pitch classes × 250
annealing steps each (plus catalog-seeded restarts at ω=2.094 and a
greedy closure polish); the 6×1400-step deep pass was terminated at
session checkpoint before emitting. Findings below are lower bounds on
what the backtracer can do; re-run with `--restarts 6 --steps 1400`.

| find | n | ne | omega | x | leak/voice | npar | wind | imbal | bipartite | cycles |
|---|---|---|---|---|---|---|---|---|---|---|
| free_0 | 8 | 8 | 1.986 | 0.384 | 0.1492 | 4 | 0.0 | 0.013 | yes | 4e:m=2; 4e:m=2 |
| free_1 | 9 | 10 | 1.804 | 0.506 | 0.1569 | 2 | 0.0 | 0.018 | yes | 4e:m=2 ×3 |
| free_2 | 8 | 9 | 2.000 | 0.375 | 0.0493 | 3 | 0.0 | 0.073 | yes | 6e:m=3; 4e:m=2 |
| free_3 | 8 | 8 | 2.225 | 0.253 | 0.0501 | 3 | 0.0 | 0.070 | yes | 4e:m=2 |
| free_4 | 14 | 20 | 2.104 | 0.315 | 0.2726 | 11 | 0.0 | 0.116 | yes | quad-mesh patch (m=2/4/5) |
| free_5 | 11 | 15 | 2.225 | 0.253 | 0.2349 | 6 | 0.5 | 0.095 | NO | 5e:m=3 + quads |
| free_6 | 13 | 19 | 1.794 | 0.514 | 0.2180 | 6 | 0.0 | 0.074 | yes | 6 quad cycles |

**What the backtracer found (and what it means):** left free, the
search does not invent new shells — it converges on **small fused-cycle
1D objects**: free_2 is a hexagon+quad fused pair (a theta-graph
fragment the catalog's theta entries only had in larger form), free_3 a
folded quad-cycle octet, free_0/free_1 twin/triple quad clusters.
Everything degree-heavy it tried got eaten by parasites and collimation,
independently confirming the catalog's class verdict. Two finds
(free_2/free_3, leak ≈ 0.049–0.050) beat every shell in the catalog but
carried one open vertex / bridges at quick fidelity, so they enter the
shortlist below only via their closed siblings (free_0/free_1). A
non-bipartite mixed object (free_5: a 5-cycle m=3 fused to quads —
cable pentagon + strut quads!) appeared spontaneously and survived at
0.235: the first evidence that **odd cable cycles fused to even strut
mesh** is a workable mixed class the catalog missed — worth a designed
candidate (pentagon-prism hybrid) in the next pass.

## 5. The comp6 reinterpretation (found by the parasite scan)

The paper score of the comp6 geometry (N=6 ring, d=1.10, φ=120°)
carries **six second-neighbour chords at 1.905** — inside the link
ceiling with probability 0.835 per pair on the real radii (expected ~5
live parasites per ring). At the locked phases these channels sit at
ψ_par = ω(2d − d_chord) ≈ 0.56 rad ⇒ **gate 0.53 — half-open unintended
channels** delivering two-step-retarded energy into every voice, exactly
the H1-cube disease (P2) at smaller scale. The 2026-07-28 fleet scored
comp6's intended gates only. This offers a cleaner reading of the
endgame anchor "wound mutual comp12 (>5000) outlived heavier one-way
comp6 (3836)": comp6 may have died of parasites, not of one-way-ness.

**Discriminating experiment (pre-registerable):** ring9_m3 (d=1.15,
φ=120°, x≈0.50) — same one-way physics as comp6, but second-neighbour
chords 2.16 are NEVER links (P=0.000): parasite-free one-way ring.
If mutuality is load-bearing, ring9_m3 still underperforms comp12; if
parasites killed comp6, ring9_m3 should beat comp12 (its back-gate leak
is 3 orders smaller). Either outcome is informative; seed nets are in
`nets/`.

## 6. Sweep 3 — exotics (paper checks)

* **Interval-rung bridges (the 3:2 molecule).** Within the seedable
  pitch window (ω 1.394–2.231, ratio span 1.60) the fifth is the ONLY
  reachable non-unison interval (octave needs 2.0). The pitch pair is
  forced to the window's edges: high cluster ω_hi ∈ [2.09, 2.23]
  (x 0.25–0.36 ✓ measured territory), low cluster ω_lo = ω_hi/1.5 ∈
  [1.394, 1.487] ⇒ **x_lo 0.79–0.90 — seedable but unmeasured-heavy**.
  The bridge locks through the coincident partial Ω = 2ω_hi = 3ω_lo;
  the bridge pair-rung
  (qω_i + pω_j)d = 2πm gives d(m=1) = 0.72 (below the foam floor —
  unreachable) and **d(m=2) = 1.44 — inside the link window**: fifth
  bridges exist geometrically. But the tongue is Tenney-narrowed:
  acceptance γ_m/(pq) = 0.0167 in det = 2ω_hi−3ω_lo, which maps to an
  rms LOAD-matching tolerance of **0.0035 (0.4%) against a measured
  foam-chaos of ±30%**, with on-tooth response 1/6. Verdict: bridge
  links are ~6× weaker and ~100× more fragile than unison links —
  molecules are real in the law table but not buildable until a unison
  particle holds its load constant to <1%. Paper-only; re-check after
  certification.
* **Beat-locked pairs (Δω on the atom cadence).** NOT expressible as a
  static lock in this kernel. Pitch is an instantaneous function of load
  (w2e = w2/(1+q·(Em+flload)/cap)) — there is no independent frequency
  dynamics to entrain, so Δω ≠ 0 persists as long as the loads differ
  and ψ sweeps at rate Δω. Conversions are atomized (quant_mode=2,
  `atoms_fire`: integrate-and-fire credit, lapse at 2ε) but the credit
  ledger has NO phase variable, so nothing couples the beat phase to
  the atom cadence. A stroboscopic "lock" would be a limit cycle of
  kappa_lock kicks crossing the dead-gate desert once per beat — a
  periodic unlock/refire radiator (roughness burst at every beat), not
  a particle. The only stable Δω≠0 relations the kernel offers are the
  comb's p:q locks — i.e. exotic #1, the fifth. CLOSED unless the
  kernel ever grows a phase-bearing credit (would be a law change).
* **Wound double-tori (linked homology).** True quad-mesh tori need
  major radius > minor + parasite clearance ⇒ ≥ ~12 rings of ≥ 6 cells
  ⇒ n ≥ ~72 (two linked: ≥ ~150) — far beyond n ≤ 32. The reachable
  linked-homology object is the **Hopf pair of wound rings**: two
  ring12_m5 circles interlocked, minimum inter-ring approach 2.41 >
  worst ceiling 2.07 ⇒ ZERO channels between the components (verified
  on the embedded pair: no cross-parasites). The bond is purely
  topological: separating the rings forces a crossing ⇒ transient
  parasite storm ⇒ death — while linked, they are invisible to each
  other. Sim-worthy as the first two-body bound state with no
  interaction channel; shortlisted.

## 7. Shortlist and negative controls

Selection: rank by leak/voice among scored candidates with x in the
seedable+measured window, zero open vertices, zero bridges; ≥2 wound
entries guaranteed (satisfied naturally). Full cell ids, positions,
per-edge labels/ψ/gates, parasite lists and phases: `shortlist.json`;
kernel-ready seeds: `nets/<name>.net` (`init=net` prints the NETGATE
cross-check — the independent-solver comparison the campaign wants).

| rank | candidate | class | n | omega | x (×x_skirt) | mass | leak/voice | npar | wind | cycle m | net |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | ring8_m3 | ring | 8 | 1.7364 | 0.558 (9.1×) | 11.2 | 0.0584 | 7 | 1.0 | 5≡3* | nets/ring8_m3.net |
| 2 | ring12_m5 (comp12) | anchor | 12 | 1.8295 | 0.488 (7.9×) | 14.6 | 0.0642 | 2 | 1.0 | 7≡5* | nets/ring12_m5.net |
| 3 | ring10_m4 | ring | 10 | 1.7040 | 0.585 (9.5×) | 14.6 | 0.0668 | 1 | 1.0 | 4 | nets/ring10_m4.net |
| 4 | ring9_m3 | ring | 9 | 1.6742 | 0.610 (9.9×) | 13.7 | 0.1068 | 2 | 1.5 | 3 | nets/ring9_m3.net |
| 5 | free_0 | free-search | 8 | 1.9856 | 0.384 (6.2×) | 7.7 | 0.1492 | 4 | 0.0 | 2,2 | nets/free_0.net |
| 6 | free_1 | free-search | 9 | 1.8040 | 0.506 (8.2×) | 11.4 | 0.1569 | 2 | 0.0 | 2,2,2 | nets/free_1.net |
| 7 | hopf12x12_m5 | exotic | 24 | 1.9026 | 0.437 (7.1×) | 26.2 | 0.1590 | 7 | 2.0 | 5, 7≡5* | nets/hopf12x12_m5.net |
| 8 | ring12_m6 | anchor | 12 | 1.9348 | 0.416 (6.7×) | 12.5 | 0.1690 | 1 | 0.0 | 6 | nets/ring12_m6.net |

\* m and N−m are the two chiralities of the same winding (traversal
convention); B1's species split. The realized foam loads (x 0.42–0.61)
sit above the design table's 0.32 because foam jitter inflates loop
lengths ~8–12% and ring_m retunes ω = 2πm/L_actual downward — the same
effect MASS.md measured; all within the validity window except where
flagged.

**Negative controls (seeded frustrated, prediction: early death):**
* `mobius6_d1.5` — Möbius quad strip, k=6 (even ⇒ odd cycles): foam
  leak/voice 0.373 + 8 parasites; the ideal ribbon additionally violates
  the link window at n=12, so it is frustrated AND barely constructible
  — the strongest available "cannot work" control. net: `nets/mobius6_d1.5.net`.
* `octahedron_d1.5` — the odd shell: paper 0.667, foam 0.721 leak/voice
  (dead-gate territory), equator 4-cycles and triangles jointly
  unclosable. net: `nets/octahedron_d1.5.net`.
(The catalog also carries icosahedron, antiprism4 and ring13_strut as
additional frustrated references; a deliberately mis-phased twin is
handled by another stream.)

## 8. Ledger — every candidate evaluated, one line each

(name, class, verdict, leak/voice score; machine-readable: `ledger.tsv`)

- `blob12` [control] n=12 — worst class ✓calibration; 0.298; no intended edges; the measured −0.232%/tu class
- `chain12` [control] n=12 — open control (2 skirt ports, penalty-gated); 0.183; measured 1.7× worse than closed ring
- `ring12_m6` [anchor] n=12 — viable; 0.169; unwound strut ring m=6; died t=2221
- `ring12_m5` [anchor] n=12 — viable, standing champion; 0.064; comp12 wound cable ring m=5; alive t=5000
- `ring6_m2` [anchor] n=6 — parasite-frustrated (§5); 0.256; comp6: ~5 expected hot 2nd-neighbour chords at 1.905
- `cube_a1.25` [anchor] n=8 — parasite-dense; 0.241; H1 cube, face diagonals 1.77 in-ceiling; died ~1.2–1.8k
- `cube_a1.5` [shell] n=8 — paper-clean, foam-poor; 0.269; P2 cube, diagonals 2.12 out
- `tower2_a1.5` [shell] n=12 — parasite-dense; 0.501; capped square tube
- `tower3_a1.5` [shell] n=16 — parasite-dense; 0.441; 3-cube stack
- `hexprism_d1.5` [shell] n=12 — parasite-dense; 0.321
- `hexprism_d1.6` [shell] n=12 — parasite-dense, foam-rare; 0.460
- `octprism_d1.5` [shell] n=16 — parasite-dense; 0.270
- `truncoct_d1.5` [shell] n=16 — best shell, still 4× behind rings; 0.246; Kelvin cell
- `rhombdodec_d1.5` [shell] n=14 — DEAD by R5 (rhombus diagonal 1.73 = 12 radiators); 0.647
- `cubeincube_a1.5` [shell] n=16 — NOT SUPPLIED (outer edges 3.23 > ceiling), as derived
- `tube6x3_d1.5` [tube] n=18 — parasite-dense; 0.430
- `tube6x3_d1.6` [tube] n=18 — parasite-dense; 0.382
- `tube8x3_d1.5` [tube] n=24 — parasite-dense; 0.365
- `tube6x4_d1.5` [tube] n=24 — parasite-dense; 0.418
- `woundtube12x2_w1` [tube] n=24 — paper-excellent (0.057), foam-rare (realiz 0.06); 0.428; B1c co-rotating tube d_ring=1.417/d_ax=1.7
- `theta_666_d1.5` [theta] n=17 — viable (strained); 0.193; junction deg-3, cycles m=6
- `theta_866_d1.5` [theta] n=19 — viable (strained); 0.184; unequal arms, realized cycle integers 7/6
- `octahedron_d1.5` [control] n=6 — FRUSTRATED (parity) ✓; 0.721; NEG control
- `icosahedron_d1.5` [control] n=12 — FRUSTRATED (parity) ✓; 1.175; NEG reference
- `antiprism4_d1.5` [control] n=8 — FRUSTRATED ✓ (slant rung unreachable); 0.854; NEG reference
- `mobius6_d1.5` [control] n=12 — FRUSTRATED ✓ (odd cycles + geometric strain); 0.373; NEG control
- `ring9_m3` [ring] n=9 — viable, parasite-free one-way, the comp6 discriminator; 0.107
- `ring8_m3` [ring] n=8 — viable, back gate G(π/2)=0.0039; 0.058; top of shortlist
- `ring10_m4` [ring] n=10 — viable, back gate G(2π/5)=0.0263; 0.067
- `ring16_m8` [ring] n=16 — viable, foam-rare at n=16; 0.183
- `ring8_m2_x0.9` [ring] n=8 — heavy probe (x=0.90 unmeasured), zero back gate G(π)=0; 0.169; ledger-only
- `ring13_strut` [control] n=13 — frustrated seam (m=6.5 impossible), π defect spread over 13 gates; 0.103; NEG reference
- `hopf12x12_m5` [exotic] n=24 — viable, purely topological bond, zero cross-channels; 0.159
- `free_0` [free-search] n=8 — viable, twin quad cluster; 0.149
- `free_1` [free-search] n=9 — viable, triple quad cluster; 0.157
- `free_2` [free-search] n=8 — best free find (hexagon+quad theta fragment); 0.049; one open vertex at quick fidelity
- `free_3` [free-search] n=8 — folded quad octet; 0.050; bridges at quick fidelity
- `free_4` [free-search] n=14 — parasite-dense quad-mesh patch; 0.273
- `free_5` [free-search] n=11 — NON-BIPARTITE find: cable pentagon (m=3) fused to strut quads; 0.235; the class the catalog missed
- `free_6` [free-search] n=13 — parasite-dense quad cluster; 0.218
- exotic `fifth_bridge` [paper] — geometrically real (bridge rung d(m=2)=1.44 in-window), 6× weaker and load-tolerance 0.4% vs ±30% chaos ⇒ deferred until a particle holds load to <1%
- exotic `beat_locked_pair` [paper] — NOT expressible as a static lock (no frequency dynamics; phaseless atom credit); a beat "lock" would be a periodic unlock/refire radiator ⇒ CLOSED for this kernel
- exotic `wound_double_tori` [paper] — true tori need n ≥ ~72 (linked ≥ ~150) ⇒ out of the n≤32 mandate; reachable form = the Hopf ring pair (shortlisted)
