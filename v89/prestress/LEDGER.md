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
| A theory | consonant-network mathematics: drift theorem, counting, the lock Laplacian, bipartite generalization, flux-as-prestress | prestress/THEORY.md | **LANDED** (see findings below) |
| D stability | sign of load-transport feedback (self-repair vs runaway), modal spectrum predictions (H2 lines), death taxonomy, active-lock assessment | prestress/STABILITY.md | **LANDED** (see findings below) |
| E regression | sparse leak/death-law discovery over the 74-log MASS corpus + pre-registered predictions for Phase-2 | prestress/regress/ | **LANDED** (see findings below) |
| B formfind | pick/tune/phase solver on the real foam (foam_s20260727.tsv); candidates C1–C8 + controls as .net seedfiles | prestress/formfind.py, candidates/ | see row below (API-error resume) |
| C morpho | topology catalog + free evolutionary search + exotics; top-8 shortlist + negative controls | prestress/morpho/ | **LANDED** (quick pass complete; deep pass died at checkpoint — rerun: `python3 search.py --tries 60 --restarts 6 --steps 1400`) |
| B formfind | pick/tune/phase solver on the real foam; C1–C8 + controls as .net seedfiles | prestress/formfind.py, candidates/ | **LANDED** (see findings below) |
| F plasticity | user proposition: cells are plastic — the frozen foam is the obstruction; harmonic misfit force realigns/morphs cells until hard; formalize (retardation plasticity / metric-from-space / node motion / work-hardening), conservation bookkeeping, convergence, experiment designs | prestress/PLASTICITY.md | **LANDED**: retardation plasticity + hardening recommended (ḋ = −κ·Φ·∂V/∂d, Φ = the S2 reactive prefactor; vacuum exactly inert; anneals cube 7%→0% in ~20 t.u.; hardening = first C1-plateau candidate); metric-from-space rank-deficient (kept as piston tripwire); node motion rejected (re-imports background). Implementing as pass D, ratchet-gated |

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
| P4 | hexagonal prism n=12 (3-regular, bipartite; CORRECTED by theory: 0 demotions needed — all-strut feasible; a P3 mining target) | counting (corrected 2026-07-29) | proposed |
| P5 | tube: stacked wound rings co-rotating (B1c rescue), axial struts at the π-rung | B1 design note | proposed |
| P6 | quad-mesh torus k×l with winding pair (m1,m2) — homology species labels | theory | proposed |
| P7 | ring12 + consonant chords (stiffening; chord rung feasibility TBD) | translation | proposed |
| P8 | mixed tensegrity splits: demote the over-constraint excess (generically 2 links on 3-regular shells) to cables carrying flux | counting | proposed |
| P9 | interval-rung bridges (p:q locked links joining unison clusters — molecules); acceptance narrowed by gamma_res_m | CONSONANCE comb | proposed |
| P10 | negative controls: Möbius quad strip / odd shell (predicted frustrated), mis-phased twins of every candidate | method | proposed |
| P11 | free-search champion: whatever the evolutionary backtracer finds that the catalog missed | morpho | proposed |
| P12 | formal stability certificate (Lean) for the winning quadratic form | user workflow §6 | blocked: no Lean toolchain (numeric certificate fallback available) |
| P13 | truncated octahedron n=24 (bipartite 3-regular) if the foam supplies it | catalog | proposed |
| P14 | beat-locked pairs: Δω synchronized to the action-atom cadence (stroboscopic lock) | exotics | **CLOSED NEGATIVE on paper** (theory §6: transport carries no atoms; strobe point 25× beyond the unlock window; sustained detune IS the roughness channel) |
| P16 | **DIODE-16** (theory §2): cable ring N=16 at φ=π/2 per link — back gate EXACTLY 0 (gate(π)=0), m=4, d≈1.08, x≈0.83; parasite-safe; highest-load candidate ⇒ load line alone predicts ~4400+; theory leak ratio 0.5 vs comp12's 1 | theory discovery | proposed |
| P17 | fifth-triangle {3:2, 2:3, 1:1} — odd cycle closed by even-m interval rungs, ω₀≈2.35, d=1.337 uniform: the first NON-bipartite consonant object; Γ/6 tongue ⇒ predicted fragile (species probe, not a mass candidate) | theory §5 refutation of strict bipartiteness | proposed |
| P18 | flight-corrected seeds: up to 27% of a wound ring's rung mass is census-invisible FLIGHT inventory — comp12 should seed x=0.233 (not 0.321), ring6 0.117; P-A predicts the retune kills the early shed | theory §7 | proposed |
| P19 | MASS SHARPNESS (user M-R1): species mass must be an attractor point, not a valley segment — hardened-annealed objects across foam seeds must cluster tighter than seed scatter (hardening = the pinning mechanism; frozen foam = the null) | user mass-spectrum note, 2026-07-29 | proposed |
| P20 | M-FAMILIES (user M-R2, neutrino pattern): ≥2 coexisting stable m-classes of one topology on one foam, masses distinct beyond scatter (ring m-ladder; chirality split; gate desert = transition barrier) | user mass-spectrum note | proposed |
| P21 | INERTIA TENSOR (user M-R3, quasiparticle pattern): kick experiments must return a tensor; pre-registered: rings anisotropic (quasiparticle-class), shells isotropic (particle-class); all-anisotropic = no true particles yet | user mass-spectrum note | proposed |
| P15 | CELL PLASTICITY (user, 2026-07-28): the frozen foam is the obstruction — plastic cells would realign under harmonic misfit force and morph until hard. Formalized: ḋ_l = −κ_plast·Φ_l·∂V/∂d (retardation plasticity; d is a link property, no background to embed in) + lock-hardening. Diagnosis CONFIRMED in-model: foam supplies 16% strut spread, shells need ≤2%; flow anneals cube to gates 0.996 in ~20 t.u.; amputates frustrated seams; vacuum exactly inert (Φ=0 with no dense) | user proposition → design accepted | building pass D (ratchet-gated) |

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

### D stability (landed 2026-07-28) — prestress/STABILITY.md

Anchors reproduced first (ring12 gates 1.0000/0.1001 = measured; load
invariant under locked flow to 5 digits). Results:

1. **SIGN VERDICT: HEAVY→LIGHT — RESTORING.** On a rung the locked
   gate offsets are exact negatives (even gates ⇒ zero odd flow at any
   entrainment); mobility symmetric above the floor; the sole surviving
   asymmetry is headroom ⇒ the heavy voice feeds the light one.
   λ_detune = −0.19/t.u. (pair, exact map; analytic −0.135/−0.166 with
   flight delay). Unlock needs dE ≈ 0.95·store — detune runaway
   unreachable. **Consonant networks self-repair internally.**
2. **The common mode is neutral: passive nets die of ENVIRONMENT
   BLEED** (vacuum-side: res_vac·mob_floor·g_plain), never of
   detuning. Independent convergence with E's load line: two agents,
   two methods, one killer — the vacuum bleed, not the structure.
3. **Mode lines (pre-registered):** ring6/cube dense sector is
   OVERDAMPED — no breathing line; slowest relaxation −0.009±0.002/t.u.
   **Wound ring12 is a CHIRAL PUMP:** co-propagating load waves grow at
   n=±1: +0.035±0.01/t.u., ν=0.64±0.05 rad/t.u.; n=±2: +0.053±0.01,
   ν=1.05±0.08; ONE propagation sense (±n split = closure-integer
   signature); saturates at x_std≈0.06–0.10 then slow drain −6e-4/t.u.
   — the measured comp12 class reproduced from theory. π-rung/one-way
   rings immune.
4. **Death taxonomy:** D1 skirt slide dominant (exp decay, knee at
   x=0.0617, terminal (t_d−t)^{1–2}); D3 seam/strangulation — **the
   ā=1.586 cube's own gravitational footprint (g1 Es depression)
   shrinks local radii and closes its own edges** — a self-
   strangulation channel for long struts; D4 parasites self-tune into
   satellite channels; D2 pump overload only for wound even rings
   (saturating). NEW DESIGN CONSTRAINT for B/C candidates: strut
   length must clear the LOADED death radius (C1' geometry: ~1.48–1.6
   band) — the parasite-free a≈1.5 cube may sit near the
   strangulation edge; check per candidate.
5. **Minimal active lock: none needed for phases — load maintenance
   only** (the C1′ piston B+ arm on a π-rung ring_m-exact ring;
   feeding is life support, not stability). Bonus discriminator:
   under quant_mode=2 credit, kicked mass loss is exactly 0 (vs
   0.3–0.9 continuous) — the action atom freezes the comma tax for
   transients.

**Campaign implication:** the two landed agents agree the killer is
the vacuum bleed. Static seed phases cannot close the vacuum-facing
gates (foam sits at vacuum pitch ω=w2; Δω sweeps every foam-facing
gate through open) — which explains c₀'s universality. The routes
that remain against the bleed: (a) the chiral pump / wound-mutual
recapture (comp12's ×2.5), (b) load maintenance (piston), (c) any
plasticity/hardening channel (agent F pending).

### C morpho (landed 2026-07-29) — prestress/morpho/ (quick pass)

Deliverables: MORPHO.md (rulebook + catalog + free search + exotics +
per-candidate ledger), search.py (kernel-verified scorer), shortlist.json,
**nets/ — 10 kernel-ready init=net seeds**. Results:

1. **The back-gate ladder** (top class): each winding step down from
   comp12's φ=150° kills back-leak exponentially while RAISING load —
   and the measured death law favors load. Top-3: ring8_m3 (leak/voice
   0.058, x=0.56), ring12_m5 re-solved on this foam (0.064), ring10_m4
   (0.067). [Cross-check vs theory's diode-16: theory puts the exact
   back-gate zero at φ=π/2; morpho's ladder tops out at ring8_m3 —
   reconcile the φ convention (gate(2φ) vs gate(φ)) when the sims run
   both. Both agree on the direction: wind down, load up.]
2. **Most surprising:** comp6 carries ~5 expected second-neighbour
   parasites (chords 1.905, in-ceiling with P=0.835, gates 0.53 OPEN at
   lock) — its death may be parasites, not one-way-ness. CONVERGES with
   theory §9 (independent methods, same finding). Discriminating probe
   shortlisted: ring9_m3 — odd-N cable ring, parasite-free by geometry.
3. **Surfaces are foam-starved at n≤32**: every 2D skin/tube is
   paper-consonant but realizes parasite-dense on this foam. Free
   search independently converges to fused small cycles; found a
   non-bipartite cable-pentagon + strut-quad hybrid class the catalog
   missed.
4. Negative controls: mobius6_d1.5 (leak 0.373), octahedron_d1.5
   (0.721), + icosahedron/antiprism4/ring13_strut frustrated refs.
5. Exotics: fifth-bridge geometrically real (d=1.44) but 0.4% load
   tolerance vs ±30% foam chaos — deferred; beat-locks NOT expressible
   (phaseless atom credit) — P14 doubly closed; double-tori need n≥72;
   the **Hopf ring pair** (two linked rings, zero cross-channels — a
   purely topological bond) shortlisted instead.

### B formfind (landed 2026-07-29) — prestress/formfind.py + candidates/

18 .net seeds (incl. mis-tuned ctrl twins), 9 reports, summary.tsv.
**VALIDATION: MATCH both anchors** (ring12: Lring 16.542, closure 5.0000,
x 0.439164 = comp12's voice masses; cube: abar 1.586, ω 1.9809, gates
0.001/0.597 at print precision — the solver IS the kernel's seeder,
bit-faithful).

1. **PHANTOM EDGES (new setup finding):** 2 of the H1 cube's intended
   edges are NOT foam links (the shell seeder never checks the link
   rule d < 1.15(rᵢ+rⱼ)) — the cube was built with two nonexistent
   channels. The init=net kernel path + formfind check both score real
   links only.
2. **Cube-deficit decomposition** (H1 cells): 0.597 kernel → 0.616
   dropping phantoms → 0.676 exact phases (min gate 0.001→0.130,
   ~130×) → 0.851 better picks → **remainder ≈0.15 is irreducible
   foam length spread (σ_d≈0.14, search-saturated)**. Exact retuning
   buys mean 0.60→0.85 but cannot reach 0.95 on the frozen foam —
   P15 plasticity is the only route past the ceiling.
3. **Best by predicted leak: c2_cube150** (leak 3.23, min 0.401, mean
   0.866, parasites all dark ≤1e-4); per-live-direction: c4 ring12+
   chords (0.146; the consonant chord algebra WORKS — unwound m=6
   ring + 4π cable chords + cross-strut weave at gates 0.995–1.000),
   c5 tube (0.167; mass-matched to comp12 within 0.2% — the clean
   structural comparator; rings m=2 exact one-way).
4. **comp12 itself carries 1 half-open parasitic link (g_b=0.48)** —
   even the champion has an unscored channel.
5. **Infeasible on this foam:** C6 torus (minor-cycle m=1 forces x to
   the 0.9 ceiling; row closures can't be integer in the link window;
   29–48 parasites), C7 truncated octahedron (spread 1.04–1.97, min
   0.035, near-soft gap 0.012). C3 all-strut hex prism over-
   constrained in PRACTICE (min 0.044) — note tension with theory's
   "0 demotions" (structural feasibility ≠ foam feasibility; the
   prism's π-rung parasite is an absorption candidate).

## Program gate (2026-07-29)

**Exact mass (M-R1 / P19) is the first MASS goal** — see `MASS.md` §5c.
Wave reading: W1 frozen mechanism nulls (expect M-R1 fail); W2 pin
(plast+harden); W3 package/discriminators; W4 **P19 first**, then P20/P21.
MASS↔EMF: `EMF.md` §5. Integer ledger: `INT_LEDGER.md` (no fleet block).
Audit: `prestress/audit/AUDIT_2026-07-29.md`.

## Verdicts

(filled as simulations complete; every simulated possibility gets a row and a
MASS.md log entry)

| id | seeded gates (min/mean) | lifetime vs control | verdict |
|---|---|---|---|
| **W1 c2_cube150** (P1/P2 retune) | 0.401 / 0.866 | t_death **908** (ctrl **698**); both dead; neither in load-line band 1250–2810 nor skin ≥4700 | **Load-line / early death.** Retune buys ~30% vs mis-phased twin, not a plateau. Consonant-skin rate-level claim **fails** on frozen foam. Integer ledger sum_err=0. |
| **W1 c8_ring12** (comp12 twin) | 1.000 / 1.000 | t_death **449** (T=3000); spectra T=200 still alive | Short structural life on this foam under int kernel; needs x50 load-line score + spectra analysis. Parasite gpar_max=0.48. |
| **W1 ring8_m3** (diode ladder proxy) | 0.000 / 0.270 | t_death **1631** — **longest** Wave-1 object | Best frozen lifetime in wave; still dies; back-gate ladder direction plausible, not a particle. |
| **W1 c5_tube6** (wound tube bet) | 0.055 / 0.862 | t_death **1042** ≪ predicted ≥4600 exception | **Not** the comp12-class exception on this seed; on/near load line. |
| **W1 c4_chords** | 0.219 / 0.926 | t_death **630** | No chirality protection from chords alone. |

## Instruments

- `foam/foam_s20260727.tsv` — the standard foam geometry (9741 cells, L=24,
  d̄=1.505), dumped from a vacuum run; identical in every campaign sim (same
  seed) so solver picks/phases transfer verbatim.
- `init=net` seedfile mode (cellfab.c, PRESTRESS-1): V/E lines; kernel places
  voices at listed cells with listed loads/phases and prints the # NETGATE
  per-edge report including parasitic links. Ratchet-gated.
- Frozen .net format: `V x y z xk th2` (0-based order = vertex id) and
  `E a b` (intended edges, for reporting only); `#` comments.
