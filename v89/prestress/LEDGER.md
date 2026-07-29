# PRESTRESS — the possibility ledger (living document)

**Campaign directive (user, 2026-07-28):** work the math theory and geometry
of the tension-surface problem OUTSIDE the simulation first — force-density /
prestress form-finding + sparse symbolic discovery — then run EVERY good
possibility through the simulation. Every avenue tried; write each possibility
down as it comes in; write each verdict once simulated; do not stop early even
if something works. MASS.md carries the dated lab-log entries; this file is
the campaign's possibility/verdict table.

## The translation (v89-native)

The prestress program maps onto the consonance laws with no imported physics:

| tensegrity concept | v89 object |
|---|---|
| force-density matrix Q | weighted phase-lock Laplacian B·W·Bᵀ (W = gate stiffness) |
| equilibrium Q·x = 0 | all-gates-open lock: ψ_f = θ_i − ω_i·d/C − θ_j ≡ 0 |
| self-stress state | circulating one-way flux with quantized closure ω·L_cycle = 2π·m |
| strut (compression) | both-gate link — needs the π-rung: d = π/ω, one length per component |
| cable (tension) | forward-only link — free φ, back-gate leak gate(−2φ) |
| prestress magnitude | seeded load x (sets ω = w2/(1+q·x)) + winding integer m |
| mechanism / soft mode | ker of the lock Laplacian beyond gauge; load-sector drift |
| tangent-stiffness PSD | lock attractivity under load-transport feedback (S2 sign) |
| form-finding | solving (picks, ω, θ, m) on the ACTUAL jittered foam geometry |

Static-lock constraint chain (to be formalized by the theory agent): θ̇ = ω ⇒
ψ̇ = ω_i − ω_j ⇒ locked components are uniform-ω; struts then need
2ω·d = 2π·m with only m=1 in the windows ⇒ all struts share length d = π/ω ⇒
struts-only closed surfaces are bipartite and equal-edged; cables generalize
ring_m (per-cycle Diophantine closure).

## Outside-sim workstreams (agents, launched 2026-07-28)

| stream | question | deliverable | status |
|---|---|---|---|
| A theory | consonant-network mathematics: drift theorem, counting, the lock Laplacian, bipartite generalization, flux-as-prestress | prestress/THEORY.md | running |
| D stability | sign of load-transport feedback (self-repair vs runaway), modal spectrum predictions (H2 lines), death taxonomy, active-lock assessment | prestress/STABILITY.md | running |
| E regression | sparse leak/death-law discovery over the 74-log MASS corpus + pre-registered predictions for Phase-2 | prestress/regress/ | **LANDED** (see findings below) |
| B formfind | pick/tune/phase solver on the real foam (foam_s20260727.tsv); candidates C1–C8 + controls as .net seedfiles | prestress/formfind.py, candidates/ | running |
| C morpho | topology catalog + free evolutionary search + exotics; top-8 shortlist + negative controls | prestress/morpho/ | running |
| F plasticity | user proposition: cells are plastic — the frozen foam is the obstruction; harmonic misfit force realigns/morphs cells until hard; formalize (retardation plasticity / metric-from-space / node motion / work-hardening), conservation bookkeeping, convergence, experiment designs | prestress/PLASTICITY.md | running |

Toolchain: numpy/scipy/sympy + Maxima present. PySR/Julia absent (hand-rolled
STLSQ instead). Lean/lake absent — P12 formal certificate recorded as
toolchain-blocked (numeric certificates instead if wanted).

## Possibilities

Status flow: proposed → solved (outside sim) → seeded → SIMULATED → verdict.

| id | possibility | origin | status |
|---|---|---|---|
| P1 | cube a=1.25 exact-retune: per-link phases + ω from actual geometry (H1 done right) | H1 verdict + form-finding | proposed |
| P2 | parasite-free cube a≈1.5: face diagonals (1.77 at a=1.25) sit INSIDE the link ceiling (~1.96) — the H1 cube carried 12 unscored off-rung face links; a≥1.35 pushes them out and puts edges at the foam's natural d̄=1.505 | found during setup, 2026-07-28 | proposed |
| P3 | equilateral-strut mining: search the foam for near-equal-length bipartite subgraphs (the foam supplies the geometry or it doesn't — measure) | translation | proposed |
| P4 | hexagonal prism n=12 (3-regular; expected over-constrained by 2 ⇒ mixed strut/cable split) | counting | proposed |
| P5 | tube: stacked wound rings co-rotating (B1c rescue), axial struts at the π-rung | B1 design note | proposed |
| P6 | quad-mesh torus k×l with winding pair (m1,m2) — homology species labels | theory | proposed |
| P7 | ring12 + consonant chords (stiffening; chord rung feasibility TBD) | translation | proposed |
| P8 | mixed tensegrity splits: demote the over-constraint excess (generically 2 links on 3-regular shells) to cables carrying flux | counting | proposed |
| P9 | interval-rung bridges (p:q locked links joining unison clusters — molecules); acceptance narrowed by gamma_res_m | CONSONANCE comb | proposed |
| P10 | negative controls: Möbius quad strip / odd shell (predicted frustrated), mis-phased twins of every candidate | method | proposed |
| P11 | free-search champion: whatever the evolutionary backtracer finds that the catalog missed | morpho | proposed |
| P12 | formal stability certificate (Lean) for the winning quadratic form | user workflow §6 | blocked: no Lean toolchain (numeric certificate fallback available) |
| P13 | truncated octahedron n=24 (bipartite 3-regular) if the foam supplies it | catalog | proposed |
| P14 | beat-locked pairs: Δω synchronized to the action-atom cadence (stroboscopic lock) — expressibility check first | exotics | proposed |
| P15 | CELL PLASTICITY (user, 2026-07-28): the frozen foam is the obstruction — plastic cells would realign under harmonic misfit force and morph until hard. d_ij is a link retardation, not a Euclidean distance (no background to embed in) ⇒ plasticity = link-property dynamics; options: misfit gradient flow on d / metric-from-space (d from live Es — the pressure law as annealing channel) / node motion / lock-hardening hysteresis (first equilibrium?). LAW-change class: ratchet-gated, vacuum must stay frozen | user proposition | agent working |

## Agent findings (as they land)

### E regression (landed 2026-07-28) — prestress/regress/

Corpus: 74 logs → 23 unique mass runs, 9 uncensored deaths. Fitted laws:

1. **Death is the load line**: t_death = 274·(x50/0.0617)^1.066 —
   R²=0.99, leave-one-class-out CV R²=0.97; one term beat every 2-term
   rival; 10/11 censored runs consistent.
2. **Universal per-voice current**: c₀ = cap·(x50−x_skirt)/t_life =
   4.25e-4 Em/t.u./voice (MAD 3%) across rings/chains/shells, seeded
   gate min 1.5e-5→1.0, N 6–12, open/closed box. **Seeded gate quality
   does NOT set the leak in the existing corpus**; roughness ≈ 0 for
   all structured objects (roughness is the BLOB killer, not the ring/
   shell killer — sharpens the H1 mechanism statement); dM/dt ≈ −c₀·N
   is the only defensible leak law from this data.
3. **One exception in 20 rows**: wound-mutual comp12 beat the line
   ×2.5 (c_eff ≤ 0.40·c₀; t_struct 4879 vs 1696 predicted) — the only
   structural lifetime effect in the corpus.

**Pre-registered predictions (the sims score these):** the exact-retuned
min≥0.95-gate cube dies ON the load line, t≈1875 (band 1250–2810),
despite perfect gates. If it instead exceeds 4700 or plateaus, the
consonant-skin mechanism is real. The wound tube is the structural bet
at ≥4600. Corpus corrections: a1 ring6 death = 1667 by census (MASS.md's
"~1900" was loose).

**Campaign implication:** P1/P2 (exact retuning) are now a clean
discriminating experiment — the regression says gates don't matter
(load line), the H1 mechanism story says they must. One of them loses.
If the load line wins, the leak channel is gate-independent (candidate:
the mob_floor trickle — to be confirmed by the theory/stability
agents), and the design target shifts from "open gates" to "recapture
the floor trickle" — which is exactly what the comp12 exception hints.

## Verdicts

(filled as simulations complete; every simulated possibility gets a row and a
MASS.md log entry)

| id | seeded gates (min/mean) | lifetime vs control | verdict |
|---|---|---|---|

## Instruments

- `foam/foam_s20260727.tsv` — the standard foam geometry (9741 cells, L=24,
  d̄=1.505), dumped from a vacuum run; identical in every campaign sim (same
  seed) so solver picks/phases transfer verbatim.
- `init=net` seedfile mode (cellfab.c, PRESTRESS-1): V/E lines; kernel places
  voices at listed cells with listed loads/phases and prints the # NETGATE
  per-edge report including parasitic links. Ratchet-gated.
- Frozen .net format: `V x y z xk th2` (0-based order = vertex id) and
  `E a b` (intended edges, for reporting only); `#` comments.
