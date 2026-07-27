# ROADMAP — state of the v89 program, open areas, and the reality ladder

**Working document.** Subordinate to `PRINCIPLE.md`. Written 2026-07-27 at
the user's direction: record the open-areas list, the phenomena ledger, and
the program state; then execute, holding one criterion above all others —
**correspondence with reality**, with every deviation named rather than
hidden. Energy conservation remains the only law; every mechanism below is
an exactly paired ledger move or an exactly unitary rotation.

---

## 1. State of the program (as of this writing)

**Kernel** (`cellfab.c`, ~1700 lines, no lattice, no pre-v89 code):

* **Dense sector** — gated transport: tail-phase gates, retarded
  entrainment, cycle-gated flight, partial comb (Tenney-weighted p:q),
  roughness radiation with space return, mutual (sqrt) coupling.
* **Field sector** — two-component signed amplitude ψ = a₁+i·a₂ on the
  plane pair; exact unitary onsite + pairwise hop rotations with
  symmetric-normalized weights ŵᵢⱼ = wᵢⱼ/√(sᵢsⱼ); superposition,
  diffraction, interference native; c emergent (group velocity ∝ field_J).
* **Instruments** — wall (detuned mute), slits, condensing screen with
  click grains and shutter, which-path recorders (kick + tap), pair seeds,
  analytic ladder mode, event-level Bell.

**Measured and standing** (all runs conserve to ≤2e-15):

| result | where |
|---|---|
| conservation at the floating-point floor through every mechanism incl. detection | all logs |
| field packet propagation; curvature ΔA ∝ E_dense (r²≈0.99), core clock slowing; radius contraction | CELLFAB §10 |
| mass translation with zero cell motion (+2.3u cohesive under round-2 mechanisms) | CELLFAB §10 |
| pair separation ladder (qωᵢ+pωⱼ)d/C = 2πm derived from the gate law; tongue profile matched point-for-point; comma tempered to δ/2 | CONSONANCE III, IV |
| the fabric pays the comma (one-sided relaxation ratchet, shed monotone in δ₀, radiated and ledgered) | CONSONANCE VI (E8) |
| interval lifetime hierarchy: unison ≫ fifth (~20 t.u.) > octave | CONSONANCE VI (E9) |
| CHSH S = 2.826 from the joint harmonic vs 1.414 LHV control | CELLFAB §6, E5 |
| double slit round 1: exact additivity I_AB = I_A+I_B — single-phase field falsified | DOUBLESLIT §7 |
| double slit round 2 (repair): fringes at parameter-free loci, V_norm = 0.316; complementarity dial V(σ) = 0.316/0.312/0.282/0.025; 223 click grains | DOUBLESLIT §8 |

**Known deviations from reality, currently standing:**

* Fringe visibility 0.32 (real experiments ≈1) — residual foam scatter +
  source bandwidth, not mechanism.
* Field dispersion is quadratic (electron-like), not linear (photon-like).
* Vacuum is Kerr-nonlinear by default (load detuning); optics runs impose
  the linear regime (q_detune = 0).
* Which-path recorders attenuate what they dephase (kick-scattering);
  real experiments can mark paths losslessly — repaired by Tier 2 below.
* Clicks are thresholded accumulations, not indivisible quanta — repaired
  by Tier 1 below.
* Two speeds: dense flight has imposed C; field v_g is emergent and
  k-dependent. One conversion ceiling should govern both.

---

## 2. Open areas of exploration (the expanded list)

**A. Consonance dynamics (dense sector)**

1. **Frequency entrainment — the choir's correction.** Locks are
   phase-entrained only; occupancy drift kills intervals (fifth ≈ 20 t.u.).
   Missing: exchange direction correlated with detuning sign
   (∝ ∂res/∂x), a paired move. Success: interval lifetimes extend by
   orders; rungs become two-sided attractors (fixes the below-rung
   fall-through).
2. **Interval identity hysteresis.** The comb argmax flips instantly; a
   held interval should need a decisive margin to change its name.
3. **Rings — the fabric's chords.** Σωᵢdᵢ/C = 2πm for triangles; narrower
   Arnold tongues; the bridge to `construct_species.c` N_c=3 locks and
   composite species.
4. **Crystallization.** The rung ladder is a preferred spacing: a relaxing
   population of loaded cells should develop pair-correlation peaks at d*
   — a lattice constant from first principles, and the honest route to
   bound extended matter (a blob as a crystal of mutual rungs).
5. **The prime limit, measured.** Tongue survival vs ambient noise
   amplitude, against computed widths.

**B. Field sector / optics**

6. **Visibility ceiling.** Longer coherent sources, scattering mean free
   path vs normalization scheme, foam regularity as a physical parameter.
7. **Nonlinear optics.** q_detune ≠ 0 is a Kerr medium by construction:
   solitons (self-focusing vs dispersion), four-wave mixing, condensation
   as saturable absorption. Tension to keep visible: real vacuum is linear
   to extreme precision; the model must explain why its field-mode
   q_detune is tiny or be wrong.
8. **Chirality/polarization.** The ± of u± as two field species: vector
   optics, lossless which-path tagging, eraser, field-level Bell with real
   analyzers. (Executed below — Tier 2.)
9. **Dispersion decision.** Quadratic (massive/electron-analog, current)
   vs engineering a linear branch. Either commit to the electron reading
   or find the massless mode.
10. **c unification.** Tie dense C and field v_g to one ceiling.

**C. The gravity–optics bridge**

11. **Bandwidth ∝ space energy.** Derive the symmetric normalization from
    Es uniformity; then Es depletion near mass ⇒ light slows/bends there:
    lensing and Shapiro delay from the same clause that enabled
    interference. Experiment: coherent packet past an E3 blob, deflection
    and delay vs the Es profile. Also still untested: free fall of dense
    structures (steepest closure).

**D. The quantum layer**

12. **Single-quantum claim rule.** The first completed conversion cycle
    claims the whole grain; amplitude guides, energy moves once, whole.
    Mostly derived (CONSTRUCTION R3 + the cycle gate applied to a
    many-ended T-mode process); the hazard ∝ |ψ|² is the one postulate
    (same law as the absorb rate). (Executed below — Tier 1.)
13. **Where ħ comes from.** If per-cycle transferred action is universal,
    E = A₀ω/2π emerges. Falsifiable fork; audit `act` data. (Executed
    below.)
14. **No-signaling audit.** The claim rule and the Bell registry are
    nonlocal-atomic; native to a stage-free ontology, but must be
    demonstrated statistically, not assumed. (Executed below.)
15. **Two-quantum interference (HOM).** Joint amplitudes as registry
    objects — one process, two ends, ψ(x₁,x₂). Conceptually licensed;
    computationally N² (reduced geometry or GPU). The many-body wall
    beyond pairs is the program's long-term wager: processes-with-many-
    ends scale with actual entanglement, not Hilbert space.

**E. Foundations housekeeping**

16. Cell birth/death (retire the audit residue).
17. Per-cell integer ω-words fusing the comb with the species census.

---

## 3. The reality ladder for the double slit (tier ledger)

| tier | phenomenon (reality) | requirement | status |
|---|---|---|---|
| 0 | fringes, envelopes, decoherence dial, conservation | amplitude field | **done** (DOUBLESLIT §8) |
| 1 | one indivisible click per quantum; fringes as a histogram of single events (Taylor 1909, Tonomura) | claim rule: first completed cycle claims the whole grain; hazard ∝ \|ψ\|²; sampling from the wave's own absorption process is *exact* for non-interacting quanta (linearity + norm-decay = survival) | **this session** |
| 2 | lossless which-path marking; quantum eraser (fringes/anti-fringes); duality V²+D² ≤ 1 | second field component (chirality); unitary tagger at slit A; analyzer-basis records | **this session** |
| 3 | delayed choice (Wheeler) | flip analyzer basis after slit transit; statistics depend only on detection-time configuration | **this session** |
| 4 | Hong–Ou–Mandel; entangled-pair delayed-choice eraser | pairwise joint amplitudes (registry ψ(x₁,x₂)); N² cost — reduced geometry or GPU | **partial** (see §5 addendum) |
| 5 | quantitative fidelity: V→1, thousands of fringes, dispersion decision, relativistic causality of claims | scale + design decisions | future |

**Falsifiable forks logged here as they resolve:**

* **F-ħ (per-cycle action universality):** *resolved 2026-07-27 —
  NOT universal.* Measured per-cycle action of standing locks (`act`
  column, E6 round-3, 42 in-tongue locked pairs of the same species at
  the same nominal frequency): 0.084–1.584, a 19× spread, median 0.62 —
  amplitude-proportional, continuous.
  The kernel's transfers are rate-proportional, so Planck's E = ħω is
  **not emergent in the current kernel**. The route is already named in
  the inherited construction: CONSTRUCTION §1.3 R1 postulates a *deposit
  quantum δ* which the kernel approximated away with continuous rates.
  Implementing fixed-δ deposits is the (major, dense-sector) change that
  could make ħ_eff = A₀/2π emerge — future work, flagged honest.
  *Update, same day:* the action-atom route was implemented behind
  `quant_A0` (auto A₀ = e_s0·d̄/C, no new constant) and regression-proven:
  toggle-off byte-exact everywhere; toggle-on leaves field/Bell/double-
  slit/HOM bit-identical, conserves at the floor, and produces the
  **photoelectric threshold unprompted** (E1: sub-atom demand cannot
  condense at any intensity). Cost mapped honestly: few-quantum dense
  stores freeze rather than circulate (blob drift and standing-lock
  exchange throttle at 2–3 atoms) — coherent few-quantum locks need
  amplitude treatment, the dense-sector successor construction. Default
  remains off; full account in `HBAR.md` §6; baselines refreshed as
  `cellfab_runs/reg_*_{R,Q}.log`.
* **F-NS (no-signaling):** *resolved 2026-07-27 — holds exactly.* The
  total (basis-summed) screen distribution is bit-identical across
  analyzer settings — native, 45°, and delayed-choice (max bin delta 0).
  Honest scope: exact by construction at the ensemble level, because the
  analyzer decomposition never feeds back into the dynamics and claim
  sampling follows the basis-independent total absorption. A deeper test
  — per-quantum collapse with feedback — belongs to the tier-4 registry
  work.

---

## 4. Execution order and acceptance criteria

Reality-correspondence is the acceptance test at each step; deviations are
recorded, not smoothed.

1. **Tier 1 — claim rule** (`n_quanta`): accept if clicks are one per
   quantum by construction-from-axiom (documented as such), and the
   click histogram reproduces the wave fringes.
2. **ħ audit**: resolve fork F-ħ from existing data (done, above).
3. **Tier 2 — chirality + tagger + analyzer**: accept if (a) full tagging
   kills fringes with *no* energy loss at the tagger (unlike the kick
   recorder), (b) the eraser revives fringes in one analyzer ledger with
   anti-fringes in the other while the total stays fringe-free, (c) the
   measured (V, D) points respect V² + D² ≤ 1.
4. **Tier 3 — delayed choice**: accept if late-basis-flip statistics
   depend only on the basis at detection.
5. **No-signaling**: accept if the total (basis-summed) screen
   distribution is analyzer-independent within foam noise.
6. Tier 4 (HOM), the gravity–optics bridge, frequency entrainment, and
   crystallization are the next campaigns after this one.

## 5. Execution log (2026-07-27, same day; logs `cellfab_runs/t1_*, q1_*, q2_*, q3_*`)

Conservation 5.6e−16 … 1.7e−15 in every run below. Kernel additions:
`n_quanta` (claim sampling), second chirality component (fb1/fb2, active
when `tag_rate > 0`), unitary slit-A tagger, analyzer-basis screen ledgers
exA/exB with `analyzer_deg` and `t_choice`.

**Tier 1 — one quantum at a time (Tonomura): PASSED.**
500 quanta, each claimed whole by the first completed cycle (one QCLICK
per quantum — indivisibility is CONSTRUCTION R3's atomicity, documented as
axiom, not fit). Individual dots at scattered (y, z); the histogram
rebuilds the fringes: 33/41 clicks in the central-maximum bins vs 2–13 in
the minima; discrete-event visibility 0.412 vs the wave's 0.466 (shot
noise). Sampling from the wave's own absorption process is exact for
non-interacting quanta in the linear regime (identical unitary evolution;
hazard = absorb law; survival = the wave's norm depletion).

**Tier 2a — duality (Englert–Greenberger–Yasin): PASSED, with the mixed
regime included.** Unitary slit-A tagger (lossless: no energy taken, no
phase noise), D estimated from the tagged fraction (D = √(2·P_V)):

| tag_rate | V_central | P_V | D | V/V₀ measured | √(1−D²) |
|---|---|---|---|---|---|
| 0 | 0.466 | 0 | 0 | 1.000 | 1.000 |
| 0.12 | 0.432 | 0.068 | 0.369 | **0.927** | **0.929** |
| 0.25 | 0.348 | 0.184 | 0.607 | 0.747 | 0.795 |
| 0.50 | 0.173 | 0.248 | 0.705 | 0.371 | 0.710 |
| 1.00 | −0.018 | 0.164 (over-rotated) | — | ≈0 | — |

Weak tags sit on the EGY *equality* (0.927 vs 0.929); strong tags go
*mixed* (residence spread disperses the marking) and fall below the bound
(V² + D² = 0.53 < 1 at tag 0.5) — precisely the pure-vs-mixed marker
structure of the real relation. Full tagging kills the fringes entirely
with zero energy cost at the tagger.

**Tier 2b — quantum eraser (Walborn-style): PASSED.**
Tag 0.5, analyzer +45°: V₊ = **+0.403** (fringes recovered),
V₋ = **−0.292** (anti-fringes — bright at predicted dark), total
unchanged (0.173). Analyzer −45°: the ledgers swap exactly
(−0.292/+0.403, totals mirrored). Even the fully-washed mixed tag
(tag 1.0, V_total = −0.02) yields ±(0.23/−0.31) on erasure — the
coherent fraction is recoverable, the mixed fraction is not.

**Tier 3 — delayed choice (Wheeler): PASSED.**
Analyzer basis flipped at t = 40, after slit transit, while the light is
in flight: the post-choice light shows erasure (V₊ = 0.257 diluted by the
pre-choice native-basis accumulation, as integration predicts); the total
distribution is identical to the fixed-basis run. Nothing about "when the
choice was made" reaches the interference — only the configuration at
detection matters, with no retrocausal bookkeeping anywhere in the kernel.

**No-signaling (F-NS): max bin delta = 0** across native / 45° / delayed
totals — see fork resolution above.

**Standing after this campaign:** tiers 0–3 of the reality ladder are
green with signatures matching the real experiments' structure (fringe
loci parameter-free; EGY equality→inequality; fringes/anti-fringes;
choice-time invariance; conservation throughout). Known open fidelity
gaps: V₀ ≈ 0.47 ceiling (foam scatter + bandwidth), quadratic dispersion
(electron-analog), imposed linear regime, ħ not emergent (fork F-ħ:
19× per-cycle action spread — route is CONSTRUCTION R1 deposit quanta).
Next campaigns, in order: tier 4 (HOM via pair registries), the
gravity–optics bridge (bandwidth ∝ E_s ⇒ lensing), frequency entrainment
and crystallization in the dense sector.

### §5 addendum — same day: F-ħ written out; Tier 4 first result

**F-ħ brainstorm → `HBAR.md`.** The situation written in full: the cycle
gate quantizes *when*, nothing quantizes *how much*. Central reframe: the
atom of transfer is **action, not energy** — fixed action A₀ per complete
cycle gives E = (A₀/2π)·ω, so Planck's relation is the statement that the
tail-call frame has a fixed size. Candidate origins ranked: (1) restore
the construction's own integer harmonic content h(v) ∈ Z^k (the kernel is
a mean-field approximation of an integer theory — F-ħ as a fidelity bug,
strongest candidate); (2) ħ as the fabric's phase-space resolution
(tongue widths compute a minimal distinguishable action); (3) the
space-cell grain A₀ ~ e_s0·d̄/C = 1.15, with the measured median act 0.617
already at half that scale — the foam sets the scale, a mechanism must
lock the spread; (4) zero-point self-consistency [G]; (5) topological
winding exchange [G]; rejected: thermal noise floor (would make ħ
temperature-dependent). Five discriminating tests defined (linearity,
photoelectric, blackbody, shot noise, cross-sector consistency).
Recommended route: integer occupations + action-quantum R1 deposits,
dense sector first, own campaign.

**Tier 4 — HOM apparatus built and run: PARTIAL.** `init=hom`: two
waveguides in the foam, tunneling coupler (gap cells at `hom_barrier`
detune), end sinks; the two quanta are the two field components; the
pair registry computes coincidences exactly from the (anti)symmetrized
joint amplitude. Findings, in order:

1. A short strong coupler (Δ/κ ≈ 5, hard edges) is **not a beamsplitter**
   — reflections and barrier population break the ±i cross-phase and gave
   inverted bunching. Long weak coupler (window 12–13 units) fixed the
   phase behavior. Coupler split calibrated to 0.45.
2. Broadband packets (Δk/k ≈ 45%) wash the exchange term entirely —
   **as in reality: HOM needs narrowband sources.** At Δk/k ≈ 18%
   (L=45, σ=5): correct exchange-statistics ordering appears —
   **g_boson = 0.42–0.43 < 0.5 < g_fermion = 0.57–0.58**, recovering
   toward the incoherent value at large delay (0.420 at dx=6 vs 0.464 at
   dx=15). Conservation 1.2e−16 at 65k cells.
3. The dip is shallow (ideal ≈ 0.01 at this split): the floor is wire-to-
   wire mode mismatch and coupler chromaticity — the same factors that
   limit real HOM experiments. Narrower bandwidth alone (σ=7) did not
   help at this box size (the packet overlaps the coupler).

Verdict: the exchange registry works and shows the boson/fermion split
with the right signs and delay dependence; the deep dip needs a larger
apparatus (long adiabatic coupler, spatial mode filtering, reservoir
source). Queued as tier-4 fidelity work alongside the §4 list.
