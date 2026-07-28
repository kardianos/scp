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
  the linear regime (q_detune = 0). *(§6 U3 removes the imposition:
  detune sourced by dense occupation only.)*
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
   fall-through). *(May emerge natively from §6 U1 amplitude coupling —
   beat-signed energy exchange; test there before adding a mechanism.)*
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
   q_detune is tiny or be wrong. *(Resolved on paper by §6 U3: detune is
   sourced by dense occupation, so vacuum is linear with q_detune ≠ 0;
   loaded matter stays Kerr. Nonlinear-optics-in-matter remains open.)*
8. **Chirality/polarization.** The ± of u± as two field species: vector
   optics, lossless which-path tagging, eraser, field-level Bell with real
   analyzers. (Executed below — Tier 2.)
9. **Dispersion decision.** Quadratic (massive/electron-analog, current)
   vs engineering a linear branch. Either commit to the electron reading
   or find the massless mode.
10. **c unification.** Tie dense C and field v_g to one ceiling.
    *(Absorbed into §6 U6.)*

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

> **2026-07-27, standing FIRST PRIORITY: §6 — unification.** One law
> table, every experiment, no per-experiment physics switches. The items
> below were the double-slit campaign's order and are executed; §6
> governs what runs next.

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

---

## 6. UNIFICATION — one law table, every experiment (FIRST PRIORITY)

**Scoped 2026-07-27 at the user's direction.** The demand: one kernel and
**one set of law constants** under which *every* standing experiment
passes, with only **apparatus** varying between runs — geometry, sources,
seeds, instruments, durations. No per-experiment physics switches. The
`quant_A0` on/off seam posed the question; the audit below shows the seam
is much wider, and §6.3 argues that almost all of it has a single cause.

### 6.1 Law vs apparatus

A parameter is a **law constant** if it states what the fabric *is* or
how conversion *works*: frequencies of the plane pair, coupling and
conversion rates, gate/resonance structure, the action atom, saturation,
space pull. A parameter is **apparatus** if it states what we built in
the fabric or how we look at it: box, duration, seeds, source amplitude/
wavevector/width, wall/slit/screen/coupler geometry, detuned-mute
strength, analyzer angles, shutter and choice times, diagnostic cadence.
Apparatus may vary per experiment (real labs differ); law constants may
not. "Pass all simulations without changing any parameters" =
**no law constant varies by experiment.**

### 6.2 The audit (12 standing cfgs + defaults, 2026-07-27)

Law constants currently taking more than one value across the battery:

| constant | default | values in use (where) | seam class |
|---|---|---|---|
| `quant_A0` | 0 (off) | 0 everywhere; >0 for threshold/ħ physics (reg_*_Q) | the quantum switch |
| `e_cond` | 0.30 | 0.2 (e1), 0.25 (e4), **99 = disabled** (e2, e6–e9, d1, t4) | law disabled per run |
| `q_detune` | 0.35 | 0.35 (dense era), **0** (d1, t4 — forced-linear optics), **1.2** (e9 — fifth reachability) | law disabled / retuned per run |
| `w1` | 1.5 | 1.5 (e2), 1.65 (d1, t4) | fabric constant, two values |
| `w2` | 0.93 | 0.93 (e3, e4), **2.9** (e6–e9 — puts the unison rung inside the link window) | fabric constant, two values |
| `field_J` | 0.06 | 0.06 (dense-era cfgs, vestigial), **1.8** (d1, t4) | fabric constant, two values |
| `p_gate` | 4 | **2** (d1 — soft gates for diffraction), 4, **8** (e3b, e6–e9 — rim seal) | sharpness knob per regime |
| `gamma_res_m` | =gamma_res | **0.02** (e3, e6–e8), **0.12** (e9; /pq → same 0.02) | width knob per regime |
| `k_dep·k_dep_m` | 1.2 | 2.4 (e3b, e6–e8), 1.2 (e9) | rate knob per regime |
| `s_pull` | 0.5 | 0.5 (e4), **0.15** (e6–e9) | space-pull (gravity coupling), two values |
| `f_conv` | 0.25 | 0.25, 0.4 (e4) | rate knob |
| `kappa_lock` | 0.9 | 0.9, 1.0 (most dense runs) | entrainment knob |
| `kappa_align` | 0.5 | 0.5, **0.1** (d1 — under-diffraction fix) | alignment knob per regime |
| `mob_floor` | 0.01 | 0.01, **0.002** (e9 — bleed control) | floor knob per regime |
| `e_click` | 0 | 0.15 (d1 family) | instrument grain that should be physics |

Notes. (a) `k_dep` alone is dead for transport since the amplitude
repair — field wants are skipped (`c==0 continue`); only the dense
product `k_dep·k_dep_m` acts. (b) The load variable is
`x = (Em+Ee)/cap` — *field* energy detunes cells, which is the sole
reason optics must force `q_detune=0`. (c) `e_cond` is a hand threshold
on `Ee`; eight runs disable it with 99. (d) Everything not listed
(`C`, `cap`, `e_s0`, foam statistics `dmin/r0/rjit`, `comb_limit`,
`rough_k`, `gamma_rough`, `mob_sym`, `lock_floor`, `sigma_tumble`,
`f_evap`, `gamma_res`) is already single-valued — the fabric is closer
to unified than the cfg sprawl suggests.

### 6.3 The thesis — one cause under most seams

The field sector, since the amplitude repair, has **no per-run physics
switches inside itself**: the same unitary law does free packets, slits,
fringes, eraser, delayed choice, Bell, and HOM, with only apparatus
varying. Every remaining seam lives in, or protects, the **dense sector's
classical rate machinery**:

* The tail-phase gate Δ = θᵢ − ω·d/C − θⱼ **is a retarded amplitude
  overlap written as a classical gate.** Under amplitudes,
  ((1+cosΔ)/2)^p is |⟨retarded src|rcv⟩|² at p=1, and *effective*
  sharpness grows dynamically with accumulated coherence — `p_gate`
  (2 vs 4 vs 8) stops being a knob.
* Entrainment is what coupled amplitudes do natively (phase pulling).
  The hand-coded imprint machinery is a rate-era model of it — and it is
  exactly what broke under quantization (imprints riding transfer lumps:
  corr_e3b drift does not recover even at ~70 atoms). Amplitude coupling
  has no lumps to ride — `kappa_lock` and the granularity defect retire
  together.
* A resonance width is 1/coherence-time. `gamma_res_m` as a tuned number
  is the rate-era stand-in — retires under amplitudes.
* `e_cond` (and its 99-disables) is a hand threshold where reality has
  the **action atom**: conversion fires when a whole atom of demand
  exists. The E1-Q photoelectric result already demonstrated this — the
  threshold *emerged* from `quant_A0` and would make `e_cond=99` in
  optics runs unnecessary (optical demand is sub-atom there).
* `e_click` is the same confession at the screen: the click grain should
  be ε(ω) = A₀ω/2π, not a knob (tier-1 already samples MCWF; only the
  grain size is hand-set).
* `q_detune=0` in optics protects against a mis-sourced load: with
  x = Em/cap (dense occupation only), vacuum optics is linear
  **automatically** — Kerr nonlinearity becomes a property of loaded
  matter, as in reality — and `q_detune` takes one value everywhere.
* `quant_A0` itself sits on the wrong boundary: as a *transport* law it
  froze flow (E3b, E6); in reality ħ lives at **mode boundaries** —
  absorption, emission, detection — and in **closure integers** (levels).
  Intra-mode transport is amplitude-continuous.

**Architecture in one line: amplitudes within a mode; atoms at mode
boundaries; integers in closure.** (This is quantum mechanics' own
architecture: unitary evolution between events; quantized exchange at
events; quantized levels from boundary conditions.) On this thesis the
parameter seams are not fifteen independent tuning problems — they are
shadows of the dense sector still lacking its amplitude completion, plus
a small residue of genuine constants to be fixed once.

### 6.4 The construction

* **U1 — dense amplitudes.** ψ_D = (b₁, b₂) per cell on the dense plane;
  Em = |ψ_D|², th2 = arg ψ_D become derived caches (mirroring the field
  repair exactly). Dense links become pairwise norm-conserving rotations
  whose **retarded phase ω·d/C is built into the rotation axis** and
  whose strength is geometry × comb weight × mutual coupling ×
  symmetric normalization ŵ = w/√(sᵢsⱼ). Gates, entrainment, tongue
  widths, and the flight delay all become interference phenomena rather
  than rate rules; the dense flight registry retires (retardation lives
  in the hop phase); C2 roughness moves to conversion events. Expected
  knob retirements: `p_gate`, `kappa_lock`, `gamma_res_m`, `k_dep_m`
  (absorbed into one dense coupling), the imprint machinery.
* **U2 — atoms at conversions only.** All S↔D↔F conversions
  integrate-and-fire in whole atoms ε = A₀ω/2π, A₀ = e_s0·d̄/C (auto —
  no new constant), **sampled by the claim rule** (MCWF, as tier 1) not
  deterministic credit. Intra-mode transport is never quantized.
  Retires: `e_cond` (threshold = one atom of demand), `e_click`
  (grain = ε at the screen), and the `quant_A0` on/off seam itself —
  the machinery is permanently on, but it lives where ħ lives.
* **U3 — load detune from dense occupation only.** x = Em/cap. Vacuum
  optics linear automatically; `q_detune` single-valued. Reachability
  constraint: the fifth needs (1+q)/1 ≥ 3/2 ⇒ q ≥ 0.5 — propose
  **q_detune = 1.2 as the law value**; tuning curves x*(d) re-derive,
  and the reachable interval vocabulary becomes a *prediction* of q
  rather than a per-run accommodation.
* **U4 — fabric constants fixed once.** Propose **w1 = 1.65,
  w2 = 2.9, field_J = 1.8** (the optics-grade values; w2 chosen so the
  unison rung d* = π(1+q·x)/w2 sits inside the nearest-neighbour link
  window at working occupancies). e2/e3/e4 apparatus recompute (source
  kx, tilt kx = w_eff(x)) — apparatus recompute is allowed; law revalue
  is not.
* **U5 — surviving rates joint-passed.** Whatever knobs survive U1–U3
  (`s_pull`, `f_conv`, `f_evap`, `kappa_align`, `mob_floor`,
  `sigma_tumble`, `rough_k`, `gamma_rough`, dense coupling scale) get
  **one value each**, selected by running the full battery jointly, not
  per experiment.
* **U6 — one ceiling (absorbs open area B10).** The same C appears as
  the retardation in both sectors' hop phases; field_J and the dense
  coupling both derive from conductance × mode overlap. Stretch goal
  inside U1's rewrite; at minimum, document the two speeds against one
  C honestly.

### 6.5 The protocol — "no parameter changes" made operational

1. **One `LAWS` block** — a canonical text block of every law constant —
   pasted byte-identical into all cfgs (`laws_check.sh` diffs the block
   across cfgs and refuses a battery run on mismatch). Apparatus blocks
   are free.
2. **The battery** = the 14 standing experiments (e1, e2, e3a, e3b, e4,
   e5 Bell, e6, e7, e8, e9, d1, t1 Tonomura, q2 eraser, t4 HOM) **plus
   the quantum set** (photoelectric threshold, grain linearity ε ∝ ω,
   click statistics) — every one under the same LAWS block, quantum
   machinery permanently on.
3. **Acceptance is physics, not bytes** (laws change ⇒ numbers move;
   the *claims* must survive): conservation at the floor; e2 packet
   advances ≥ 0.3C; e3a sealed / e3b translates along the tilt; e4
   ΔA–E_dense linearity r² > 0.95; e5 S > 2.7 vs LHV 1.41; e6/e7 tongue
   matches the computed profile; e8 shed monotone in |δ₀| (one-sided);
   e9 fifth locks at computed d* for tens of t.u.; d1 fringes at
   parameter-free loci with V_norm ≥ 0.25 + tiers 1–3 signatures
   (clicks rebuild fringes; EGY equality; eraser ±; delayed choice;
   no-signaling); t4 g_b < 0.5 < g_f with delay recovery; threshold:
   sub-atom demand condenses nothing at any intensity; linearity:
   grain ∝ ω across modes. Deviations recorded, not smoothed.
4. **Baselines**: current set frozen as `reg_*_{R,R2,Q}` (committed).
   The unified battery writes `reg2_*` — a new baseline generation, not
   edits to the old one.

### 6.6 Execution stages

* **S1 — cheap decisive test first (no dense rewrite).** Implement U2
  (atoms → conversions only; transport unquantized) + U3 (x = Em/cap)
  on the current kernel; draft the LAWS block with §6.4's proposed
  values; run the full battery. This alone predicts: flow physics
  restored (e3b drift, e6 locks) *with* thresholds kept, optics linear
  without `q_detune=0`, `e_cond`/`e_click` switches gone. If S1 passes
  everything, U1 becomes fidelity work rather than rescue; if it fails,
  the failure localizes the thesis.
* **S2 — U1, the dense amplitude completion** (the big rewrite), walls
  pre-listed from the field repair: sign/seed pairing (κ = −k with
  exp(−iτX)), symmetric normalization against foam disorder,
  direction-aware metrics (energy centroid, not radius), linear-regime
  discipline. Re-derive the CONSONANCE laws in amplitude form (the
  ladder derivation §III carries over — the gate algebra is unchanged,
  it just becomes literal).
* **S3 — U4/U5 joint-pass sweep** → freeze THE law table; write
  `reg2_*` baselines; update CELLFAB.md §2 laws and the params table.
* **S4 — quantum battery under unified laws**: HBAR §4 discriminating
  tests (linearity, photoelectric, blackbody bath, shot noise,
  cross-sector A₀ consistency); update forks F-ħ and the tier ledger.

### 6.7 Honest risks

* **w2 = 2.9 blob re-seal is not guaranteed** (ω·d moves regime); if the
  blob cannot seal at the unified w2, that tension — containment window
  vs rung window — is a finding to record, possibly the first real
  constraint on the fabric's frequency pair.
* **E8's one-sided ratchet must re-emerge** from interference + atomic
  conversion; the current shed mechanism lives in rate machinery. If the
  comma does not get paid under U1, the thesis loses a limb.
* Claim-rule conversions add stochasticity to e4 curvature — linearity
  must survive ensemble averaging.
* The fifth's numbers move under one q_detune (they were tuned at 1.2
  with bespoke rates); e9's *claim* (locks at computed d*, finite
  lifetime, Tenney ordering) is what must survive.
* Field sector and quantum tiers 0–4 are untouched by design — any
  regression there is an implementation bug, full stop.
