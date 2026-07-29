# RESUME — session checkpoint 2026-07-29 (PRESTRESS campaign paused)

**Paused at user request** ("agent ops, possibly restart") after the
outside-sim phase COMPLETED: all six math agents landed, both kernel
instruments built and gated, zero simulations of candidates run yet.
**The entire Phase-2 fleet is ready to launch** — 28 kernel-ready seed
files on disk with pre-registered predictions to score them against.

Read order for a cold start: `v89/PRINCIPLE.md` → `v89/README.md` →
`v89/prestress/LEDGER.md` (the campaign's possibility/verdict table,
all agent findings summarized) → this file. Lab log: `v89/MASS.md`
(dated entries through 2026-07-29 cover everything below in detail).

## What this campaign is

User directive (2026-07-28): work the math theory and geometry of the
tension-surface problem OUTSIDE the simulation first (force-density /
prestress form-finding + sparse symbolic discovery, via agents), then
run EVERY possibility through the simulation. Do not stop early even if
something works. Later additions: P15 cell plasticity (the frozen foam
is the obstruction; cells would realign by force — user proposition);
mass-spectrum requirements M-R1..M-R4 (specific mass / neutrino-style
families / quasiparticle orientation-dependence — user note, gates
"taking mass seriously"; MASS.md 2026-07-29, LEDGER P19–P21).

## State of the world

- Git: all agent deliverables + MASS.md entries staged on top of
  `6ea35b6` (init=net + scaffolding + regression, battery 21/21).
  Final commit of this checkpoint happens after the plasticity battery
  verdict (see "Pending at pause").
- Kernel `v89/cellfab.c` — two additions this session, both default-off:
  1. **init=net** (COMMITTED, battery 21/21): places externally
     form-found networks from `.net` files (`V x y z xk th2` /
     `E a b`); prints `# NETGATE` per-edge + parasite report at seed,
     `# NETG` live-gate summary each diag row, final `# NETGATE F`
     dump. Smoke-validated exact (π-rung pair gates 1.0000/1.0000).
  2. **kappa_plast / tau_harden** (P15 retardation plasticity + lock
     hardening; UNCOMMITTED pending the battery run that was in flight
     at pause): the S2 odd term on its geometric conjugate —
     `Δd = −κ·base·Sm·dV/dd` buffered in pass 2 (c==1 branch), applied
     in new pass D (between pass 5 and pass F), Jacobi, deterministic;
     `# PLAST` diag row (serial, byte-stable); floor d≥0.5; vacuum
     PROVEN exactly inert (commissioning run: links_moved=0,
     max|d−d0|=0, drift 0.000e+00). kappa_plast=0 is the exact legacy
     kernel (guard-skipped memset+print only).
- Foam standard: `v89/prestress/foam/foam_s20260727.tsv` (9741 cells,
  L=24, d̄=1.5053, seed=20260727 — same foam as all measured anchors;
  every solver pick/phase transfers verbatim to sims on this seed).
- Harness: `v89/prestress/run_net.py` (laws_V2g + standard apparatus:
  edge_sink=1.6, lump_diag=1, snapshots, renders via viz/render_slice,
  auto-scoring) and `score_net.py` (frozen verdict definitions: m0,
  t_death = largest lump < 0.3·m0, late dM/dt, rough share) →
  `runs/SCORES.tsv`.

## The six agents — landed results (full text in LEDGER.md + each doc)

1. **E regression** (`prestress/regress/`): death is the load line
   t_death = 274·(x₅₀/0.0617)^1.066 (CV R²=0.97); universal per-voice
   current c₀=4.25e-4 Em/t.u. — seeded gate quality does NOT set the
   leak in the corpus; comp12 the lone ×2.5 exception. PRE-REGISTERED:
   retuned cube dies on the line t≈1875 (band 1250–2810) despite good
   gates; consonant-skin real only if >4700 or plateau; wound tube
   ≥4600.
2. **D stability** (`prestress/STABILITY.md`): SIGN VERDICT heavy→light
   (restoring, λ=−0.19/t.u.; detune runaway unreachable) — networks
   self-repair; the killer is the VACUUM BLEED (common mode neutral;
   res_vac·mob_floor·g_plain). Ring/cube overdamped (no breathing
   line, slowest −0.009/t.u.); wound ring12 = CHIRAL PUMP (n=±1 grows
   +0.035/t.u. at ν=0.64 rad/t.u.; n=±2 +0.053 at ν=1.05; one sense;
   saturates x_std 0.06–0.10 → −6e-4 drain = the comp12 class,
   derived). Death: D1 skirt slide dominant; D3 SELF-STRANGULATION (a
   heavy structure's own g1 footprint closes its long edges — struts
   must clear the loaded death radius ~1.5-band); D4 parasites
   self-tune. Active lock: load maintenance only (piston B+ = life
   support). Kick discriminator: quant credit freezes the comma tax.
3. **A theory** (`prestress/THEORY.md`, t1–t5): drift theorem (one ω
   per component; struts m=1 only, d=π/ω ∈ (1.16,1.62), tolerance
   ±3.65% for gates≥0.9); force-density matrix Q=BWBᵀ (comma law
   ψ∝1/w); DIODE POINT φ=π/2 → back gate EXACTLY 0 → **diode-16**
   (N=16, m=4, d≈1.08, x≈0.83 — highest-load candidate, load line
   alone ⇒ ~4400+); strict bipartiteness REFUTED (fifth-triangle
   {3:2,2:3,1:1} closes odd cycles, fragile); FLIGHT-LOAD FIXED POINT
   (27% census-invisible → corrected seeds comp12 x=0.233, ring6
   0.117); gpl(N)=cos⁴(2π/N) (ring6 9× conductance-poor); parasite
   rule d_nonedge ≥ 2.08 ⇒ cube needs a≥1.47; P-A..E pre-registered
   (incl. flux moment comp12 ≈2.15 vs 0 unwound).
4. **C morpho** (`prestress/morpho/`, 10 nets): THE BACK-GATE LADDER
   (wind down from 150°, load up — ring8_m3 best 0.058 leak/voice);
   comp6 death likely parasites (5 chords at 0.53 open — CONVERGES
   with theory); ring9_m3 = parasite-free odd-N discriminator;
   surfaces n≤32 foam-starved (parasite-dense); free search found
   non-bipartite pentagon+quad hybrid; Hopf ring pair (topological
   bond, zero cross-channels) shortlisted; Möbius/octahedron negative
   controls ready. Deep pass rerun: `python3 search.py --tries 60
   --restarts 6 --steps 1400`.
5. **F plasticity** (`prestress/PLASTICITY.md`): the recommended law
   (now implemented, see kernel above); foam supplies 16% strut spread
   where shells need ≤2% (frozen-foam obstruction CONFIRMED); flow
   anneals cube 7%→0.00% in ~20 t.u.; amputates frustrated seams;
   hardening = the first true-equilibrium (C1 plateau) candidate AND
   the mass-sharpness (M-R1/P19) pinning mechanism; κ window ≥50×.
   PLAST-0 partly done (vacuum inertness exact); PLAST-1..3 designed
   with bars.
6. **B formfind** (`prestress/formfind.py`, 18 nets + 9 reports):
   validated bit-faithful against both anchors; PHANTOM-EDGE confound
   found (2 of the H1 cube's edges were not links); cube-deficit
   decomposition 0.597→0.851 with ≈0.15 irreducible (frozen-foam
   ceiling — plasticity is the route past); best: c2_cube150 (leak
   3.23, parasites dark), c4 chords (algebra works, gates 0.995+), c5
   tube (mass-matched comp12 ±0.2%); comp12 has 1 half-open parasite;
   torus + truncated octahedron INFEASIBLE on this foam.

**Convergences worth trusting** (independent methods agreeing): the
vacuum bleed as the killer (E+D); comp6's parasites (A+C); surfaces
foam-starved (B+C); parasites-dark-when-locked (B+F).
**Open cross-checks**: diode φ convention (theory's diode-16 vs
morpho's ring8_m3 G(90°) — which angle the back gate sees; sims
measure both); theory's "hex prism 0 demotions" vs B's measured
over-constraint (structural vs foam feasibility); κ_plast calibration
(probe in flight at pause); load line vs theory's λ ∝ √x·res leak
ratios (ring6:comp12:diode = 6.6:1:0.5) — the fleet decides.

## Pending at pause (two background runs)

1. **Plasticity ratchet battery** (`battery.py --laws laws_V2g.cfg`,
   task be5cn2yga): MUST be 21/21 before the kappa_plast kernel change
   commits (kappa_plast absent from laws_V2g ⇒ default 0 ⇒ expected
   green; the earlier init=net battery was 21/21).
2. **Anneal probe — LANDED, the mechanism works in the real kernel.**
   Mis-tuned pair (seed mean gate 0.913, ψ_b=0.31): annealed to BOTH
   gates 1.0000 within 10 t.u. at kappa_plast=1.0; the link walked
   d 1.4999→1.5420 to a joint (θ, ω, d) lock point (phases and pitch
   co-adapt, so the rest point is nearer than the naive d*=1.577;
   final ψ_f=−0.0006, ψ_b=+0.0004). κ*=1.0 provisional (τ_anneal <10
   t.u. — near the adiabatic floor; 0.3–0.5 would sit mid-window if
   PLAST-1 shows churn). HONEST SURPRISE to watch in A/Bs: 40185
   foam links moved (max|Δd|=0.116) — the pair's leaked crumbs make
   the SURROUNDING foam plastically active (a leaking object slowly
   rewrites vacuum geometry around it; doctrine-flavored — space
   responds to matter — but it can shift local skirt/conductance:
   PLAST-1's control comparison and the PLAST-0(c) κ*-battery arm
   are the guards). Log: scratchpad/plast_pair.log.

## Program priority (updated 2026-07-29, integer default)

**Exact mass (M-R1 / P19) is the first goal of the MASS program** —
sharp universal package, not merely long life. Full statement:
`MASS.md` §5c. PRESTRESS waves below are the mechanism → pin →
spectrum path toward that goal.

**MASS ↔ EMF:** parallel cheap EMF (EM1/EM2, C4 tools) OK anytime;
EM5/P2 only after MASS checkpoints. See `EMF.md` §5.

**Kernel default = integer ledger (mode 3):** `cellfab.c` → `cellfab`.
FP64 reference: `cellfabf.c` → `cellfabf`. Results:
`int_ledger/RESULTS.md` (20/20). Use `ledger_mode=0` only for A/B.

**Plasticity:** in kernel (κ via `--extra` only). Plan audit:
`prestress/audit/AUDIT_2026-07-29.md` — GO_WITH_FIXES.

## Next-session runbook (PRESTRESS-2, task #44 — Wave 1 COMPLETE 2026-07-29 (Go fleet.go; all 7 jobs; see LEDGER Verdicts))

Fleet = every avenue, mis-tuned controls, ≥3 foam seeds for scored
claims, visuals per run (standing practice). Launch pattern:

    python3 v89/prestress/run_net.py --name <run> --net <file> \
        [--T 3000] [--extra kappa_plast=1.0] [--threads 4]

(2 concurrent × 4 threads on the 8-core box; T=3000 ≈ 15 min each;
run in background, never foreground >10 min. **Absolute `--net` paths.**)

### Wave 1 — frozen foam mechanism nulls (κ=0) — expect M-R1 fail

Prep: rebuild/smoke cellfab; T=20 smoke_pair. Spectra first or parallel
with first long job: `c8_ring12` `--T 200 --snap_every 250`.

Runs: c2_cube150 vs c2_cube150_ctrl (log as **best frozen retune**
min~0.40/mean~0.87 — **not** PREDICTIONS min≥0.95); c8_ring12;
**ring8_m3 as diode ladder proxy** (emit true diode-16 offline, optional);
c5_tube6; c4_ring12_chords. Score: load line on **measured x50**, leak
ratios, spectra ν, flux moment if instrumented. T=3000 decides in-band
death; T≥5000 only if alive/near skin≥4700 or tube≥4600 upset.

Optional parallel: EMF EM1 or EM2 on free cores (`EMF.md` §6 Mode P).

### Wave 2 — pin hunt (plasticity + hardening) — first exact-mass chance

PLAST-1: c2_cube150 ±kappa_plast=κ* (bars: live gates ≥0.9 by t≈200,
roughness halved, leak below control; prize = C1 plateau). PLAST-3:
c1_cube125 + parasites + plast (self-seal). PLAST-2: anneal-then-kick
with tau_harden. **Any plateau → C4 packet + M0 field-channel before
"exact mass" language.** Watch foam plastic contamination vs κ=0 control.

### Wave 3 — package + discriminators

Flight-corrected comp12 (x=0.233, P-A) — mass = bound+flight; ring9_m3;
negatives (Möbius, octahedron, ctrl twins); Hopf; fifth-triangle;
free-search. Size/N as species axes only.

### Wave 4 — spectrum gates (MASS first goal lives here)

**P19 first** (hardened multi-seed mass cluster vs frozen null) — gate
for taking mass seriously. Then P20 m-families; P21 inertia tensor
(ring vs shell; quant comma-tax freeze). Prefer ≥5 seeds for C3-class
claims.

Verdicts: one row per run in LEDGER "Verdicts" + dated MASS.md entry.

## Task list at pause

- #44 PRESTRESS-2 (pending): Wave 1 fleet (+ optional EMF Mode P).
- #45 PRESTRESS-3 (pending): verdicts → LEDGER/MASS.md/memory/commits.
- #46 PRESTRESS-4 (pending): **P19 first**, then P20/P21.
- #41 MASS-C1' piston (pending): load maintenance + locus-shift tripwire.
- #39 MASS-A3/A4 (pending): unblocks on first true plateau + C4; **A4'
  P19** before "standing particle" language.
- EMF Mode P (pending, parallel): EM1 and/or EM2 apparatus-only.
- INT_LEDGER (design only): no implementation task until opened.
- #32/#33 SPEED-3/4 (deferred): only if fleet wall time hurts.
- Historical #29–#40: completed, see MASS.md.

## Standing rules that bind the next session

Ratchet (full battery per kernel/law mod; bars sharpen only); no
pre-v89 material; visuals every campaign; deterministic (byte-identical
at any thread count — the # PLAST diag is serial for this reason);
laws_V2g remains THE standing table (kappa_plast enters experiments as
an explicit extra key, never silently into the table); fleet runs in
background; absolute paths (Bash cwd resets); exact mass before g5;
EMF law changes only at MASS checkpoints.
