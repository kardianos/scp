# v95 developmental design — the J-dependent link metric (the space-direction transport partner of the door pair)

**Session ID:** `415d0caa-80d0-4689-b10a-ea15e29ba881` (resume with
`claude --resume 415d0caa-80d0-4689-b10a-ea15e29ba881` when results land).
**Author:** claude (fable), developmental session, 2026-08-09.
**Deliverable for:** the v95 brief (J-dependent link metric: theory, equations,
numerical method, bars, composition with the door pair).
**Grounded in (read this session):** `v94/FINAL_STATUS.md` (the scoreboard),
`v94/consult/DESIGN_claude_face3.md` §Y/§Z (the frozen-phase ratchet, the
Z.4 priority order, the v95 architecture registration),
`v94/consult/M4_RESULTS.md` (R2 0.24 vs baseline 0.62; hot-bath self-heal),
`v94/consult/DSG_RESULTS.md` (door_avg suppresses only *incoherent* outtake),
`v94/RING_DNLS.md` (A.1–A.5 retention anatomy; B.1–B.4 DNLS map; B.3 the
antagonism), `v94/RECKONING.md` (R1–R3: condensation is a transient; the
raised stability bar; the J-metric named as the next face),
`v93/README.md` §II.5 (the dynamical-Hodge-star design language),
`v94/kernel/freecell.c` (pass S + FLOW bed law :1329–1409, pass 2 τ_s
:1572–1597, pass U :1624–1712, `apply_uhop` :1296–1313, slot rebirth :1148),
`v94/fab/step.go` (the Go mirror), `v94/laws/laws_V3a.cfg`.
**Status:** design + **resume update 1 (PART VII)** — J1 GREEN (exceeds
target, best in-medium phase retention ever measured), J2 RED (diagnosed:
the metric is an amplifier, not a generator), J3/J4 GREEN (partial), the
composition bars C1–C5 scoped. Where PART VII supersedes a bar or a
grading, PART VII governs.

---

## Executive summary

The design's central move is a recognition, not an invention: **the
substrate already owns a dynamical link metric — the FLOW bed-digging law**
(`bed_k`/`sbed`, adopted architecture of the v91 ASYM campaign, live in the
kernel for the *space* current at `freecell.c:1379–1409`): a per-slot
conductance with a low-pass memory of the actual moved flow, multiplicative
growth on the persistent part, a zero-sum per-cell renormalization (the
anti-ignition structure), and clamps. v95's J-dependent link metric is
**the same law applied to the dense channel's link current** — one metric
law, both currents. The v93 §II.5 language ("a dynamical Hodge star — link
weights evolving toward the flux pattern — is precisely FLOW's bed-digging
law") already stated this identification; v95 executes it.

What makes it dense version *phase-structure-aware* — the property the
space version never needed — is the same physics as face-3's `door_avg`
fix, applied one sector over: **the low-pass is synchronous detection.**
The dense link current on an incoherent contact oscillates at the pitch
detuning (period ≈ 7.7 t.u. at ring-in-bath contrast) and integrates to
≈ 0; the current carried by a coherent structure (a winding's circulation,
a mover's forward current) is DC and survives the integral. So coherent
structures dig their own channels; incoherent contacts starve — and the
zero-sum renorm turns that contrast into a competition each cell must
settle, which replicator dynamics settles exponentially in favor of the
persistent current. The measured numbers say the contrast is enormous where
it matters: a wound C6 ring's internal circulation is |J| ≈ 1.4 per link
while its bath-contact drain is ~10⁻²–10⁻³ — the ring self-insulates
(§I.5). **[RESUME 1: MEASURED — R2 0.774 vs 0.131, Em_min ×14–×670; the
self-insulation is real. §VII.2.]** A tilted packet's forward current is
O(1) *in total* but diluted over hundreds of links and transient (§I.6)
— **[RESUME 1: MEASURED INSUFFICIENT — J2 RED; the dilution+transience
estimate in §I.6 was the design's one wrong load-bearing assumption;
diagnosis and the composition road in §VII.3–§VII.4.]**

The three derived claims, each with a pre-registered bar: **(1)
phase-sensitive binding** (J1) — the metric attenuates exactly the contact
links through which RING_DNLS A.2/A.3 measured both the leak and the
dephasing entering, stretching the winding's medium lifetime **[RESUME 1:
GREEN, exceeds target]**; **(2) antagonism release** (J2) — the self-trap
criterion `δΩ ≳ z·κ` is fought on the bandwidth side (`κ → κ·w`, up to the
×5 clamp), un-trapping the law corner *on the dug track only* while static
structure keeps its binding **[RESUME 1: RED for the bare blob — no
current source; the machinery is confirmed by J1, the source moves to the
composition, §VII.4]**; **(3) conservation untouched by theorem** (J3) —
the metric enters only through Givens angles **[RESUME 1: GREEN, drift=0]**.
Anti-ignition (J4) is structural: bounded weights, zero-sum per cell,
unitarity untouched, and a derived digging floor (§II.6) that puts bath
noise ~4 orders of magnitude below engagement **[RESUME 1: dark-bath arm
GREEN (structural byte-identity); lit-bath stats + hoard-floor arms open]**.

The composition (§V, re-scoped §VII.4): **the door pair owns the mode
direction** (coherent conversion DS-G3 r = +0.684; anti-runaway 17 vs 417;
the chi_load park, 10× slope) — measured, won. **The J-metric owns the
space direction** — binding/retention half now measured won (J1);
transport half awaits a current *generator*, which is exactly what the
door pair's metabolic cycle is (C2) and what topology already is (C4).
Pair alone confines (M5/Z.2); metric alone anchors (J2); the v95 object
needs both.

Every load-bearing assumption is flagged inline and collected in §VI;
§VII.6 re-grades them against measurement.

---

# PART I — Theory

## I.1 The postulate, from the founding line

The programme's one line: *energy is never destroyed; it only changes mode;
space is one of the modes.* v89 made the static half of that line
measurable: a mass maintains a graded space depression (g1); matter is
emergently opaque (g3); and — decisively — **a footprint requires
throughput** (g4: a blob's space flux is all mass-rate bookkeeping; no
steady monopole without an internal cycle).

The J-metric is the *transport* half of the same line, stated as a
postulate:

> **P-metric.** The conductance of a spatial link is not fixed geometry;
> it responds to the energy current the link actually carries. Sustained
> current converts the channel — digs it — and the space so converted is
> conserved: a cell's total link conductance is fixed; current only
> *redistributes* it among the cell's links.

The second clause is the load-bearing one. "Space curves because energy is
conserved and space is convertible" — convertible, not creatable. A cell
has a fixed budget of channel (its share of the roster's geometry); a
persistent current through one link claims budget *from the cell's other
links*. This is simultaneously the physics (curvature as redistribution,
never net creation of conductance) and the anti-ignition structure (no
global transport amplification is expressible — the same zero-sum
architecture the v91 ASYM campaign derived as "the LEDGER internalized as
the bed-digging channel law" and the kernel already implements for the
space current).

**[Load-bearing assumption A1: the zero-sum budget is per-cell (the node's
link star), not global or per-link. This is inherited from the adopted
FLOW law, where it was measured anti-ignition-safe for the space current;
it is *assumed*, not re-proven, for the dense current — bar J4 measures
it.]**

## I.2 The form, derived three ways (they agree)

The brief asks: a relaxation? a variational principle? a Hodge-star
dynamics? All three, and they are the same law.

**(a) As Hodge-star dynamics (the DEC reading, v93 §II.5).** The dense
channel is a 0-cochain ψ_m on cells with a 1-cochain current J on oriented
slots. The hop normalization/weight is the discrete Hodge star on 1-forms
(the metric's appearance in the discrete Laplacian). A *fixed* star is a
fixed background metric — tolerated only because it was inherited from the
roster's geometry. The dynamical star `d★/dt = f(J)` is the minimal
background-free completion: the metric becomes state that responds to the
physics, exactly as the roster's positions already do (live geometry,
tumble, bond plasticity). v93 §II.5 named this "the one place 'space
curves' gets an implementable form."

**(b) As a variational principle (adaptive transport networks).** Consider
the dissipation-like functional of the dense channel at frozen amplitudes,
`D[w] = Σ_s Ĵ_s²/(w_s·c_s)` (c_s the geometric factor), minimized over the
weights subject to the per-cell budget `Σ_{s∋i} w_s = n_i`. Stationarity
gives conductance concentrated on the links carrying the largest persistent
current — the optimal-channel solution (the same variational structure as
river networks and adaptive biological transport). The *gradient flow* of
this constrained problem is multiplicative growth on |Ĵ| with per-cell
renormalization: precisely the bed law. Equivalently (§II.3): the
constrained multiplicative law **is replicator dynamics on each cell's
link simplex with fitness |Ĵ_s|** — each cell's budget flows toward its
highest-persistent-current links at a rate set by the fitness *gap*. This
is why a modest coherent/incoherent contrast suffices: the renorm
competition amplifies any persistent gap exponentially in time (§II.3).

**(c) As a relaxation with memory (the implemented form).** Operationally
the law is a first-order relaxation of a per-slot memory toward the actual
moved flow, then growth on the memory's magnitude — the exact arithmetic
shape of the adopted FLOW law, reused verbatim (§II.1). The memory's time
constant τ_J is what makes the law a *synchronous detector* (§I.4): it is
not a numerical smoothing knob but the physical statement "the channel
responds to the persistent component of the current."

**[Load-bearing assumption A2: the fitness is |Ĵ| (unsigned) — a channel
is a channel, direction-blind, mirroring the bed law. The signed
(vector-credit) refinement exists in the ASYM record; it is NOT taken
here. If J2 shows track-digging helps the *backward* slosh as much as the
forward motion, revisit.]**

## I.3 What the dense current is (kernel-exact)

Pass U's hop on slot s (canonical orientation i→j) is the unitary coupler
`ψ_i' = cosτ·ψ_i − i·sinτ·ψ_j`, `ψ_j' = cosτ·ψ_j − i·sinτ·ψ_i`
(`freecell.c:1296–1313`). The exact per-hop transfer into j:

```
ΔEm_j(τ) = sin²τ·(Em_i − Em_j) + sinτ·cosτ·J_s ,
J_s ≡ 2·Im(ψ_i·ψ_j*) = 2·(m2_i·m1_j − m1_i·m2_j)
```

— a diffusive (population-difference) term plus the coherent cross-term.
`J_s` is the v93 link current / dense momentum up to the fixed sign
convention of the −i coupling (v93's docs write `J = 2·Im(ψ_i*·ψ_j)`; the
two differ by sign; the metric law below uses only magnitudes and the
*measured* transfer, so no convention can leak into the physics — but the
implementation must pin one convention for the meter, and the Sincos-swap
lesson says: **do not compute J analytically in the kernel; measure the
actual ΔEm_j the hop moved**, which `apply_uhop` already computes for the
p1 meter).

The metric's input is therefore the **actual moved flow**, per slot per
step (summed over the Strang half-hops):

```
F_s = ΔEm_j(step) / dt        (signed, canonical i→j orientation)
```

This is the same choice FLOW made ("net memory from the ACTUAL moved
flow") and it is robust by construction: whatever the schedule, the
convention, or the gate, F_s is what the channel really carried.

## I.4 Phase-structure-awareness: the low-pass is synchronous detection

The memory is a one-pole low-pass of F_s with time constant τ_J:

```
Ĵ_s ← Ĵ_s + (dt/τ_J)·(F_s − Ĵ_s)
```

Its transfer function attenuates a component of F_s oscillating at
frequency δω by `1/√(1 + (δω·τ_J)²)` and passes DC at 1. Now read the
measured situations through it:

- **Incoherent contact (ring voice ↔ bath cell).** The relative phase
  advances at the pitch gap δω = w2e_bath − w2e_ring ≈ 2.9 − 2.086 = 0.81
  rad/t.u. (RING_DNLS common config), so the cross-term component of F_s
  oscillates with period ≈ 7.7 t.u. and integrates toward zero. What
  survives is the slow diffusive drain (small) plus a random-walk residue
  `~σ_F/√(τ_J/dt)`.
- **The winding's own circulation.** On the C6 cycle with th2 = m·φ, every
  cycle link carries `J = 2·Em·sin(m·Δφ)` = 2·0.81·sin(120°) ≈ **1.40**
  (m=2, Δφ=60°, Em=0.8125/voice) — DC, locked by the topology, maintained
  for as long as the winding lives (the winding IS a persistent current;
  this is the transport-side face of "W=+2 exact ≥200 t.u.").
- **The mover's forward current.** The tilted packet's `fd` current is
  coherent and unidirectional (+287/+297/+277 across seeds, L1-A) — DC
  while the packet moves. **[RESUME 1: in total, yes — per link, no; the
  dilution is the J2 killer. §VII.3.]**
- **A static hoard.** Net J ≈ 0 (breathing/Rabi-like internal exchange is
  zero-mean); its slow diffusive leak is ~10⁻⁴–10⁻⁵ per link per t.u.
  (RECKONING R1 decay rates spread over dozens of links) — below the
  digging floor (§II.6) by orders of magnitude.

So the metric *is* the phase-structure detector the brief asks for:
structures with definite-sign J (coherent) bind the metric; structures
whose J averages to zero (incoherent) leave it flat. This is `door_avg`'s
synchronous-detection argument (face-3 §Y.4) applied to the spatial
channel — the same fix, the other direction. One honest difference,
carried from DSG finding 2: synchronous detection distinguishes
*persistent* from *oscillating*, not "good" from "bad" — a steady drain is
also DC (§I.5 shows why the competition handles it; §VI flags the residual
risk).

**[Load-bearing assumption A3: τ_J is the design's ONE new scale. Default
τ_J = 4 t.u., chosen so that the slowest measured incoherent oscillation
(δω ≈ 0.8 → period 7.7) is attenuated ≈ 3.3× while the response is fast
enough to engage before the B.3 stall (~4 t.u.). Neither end is fully
comfortable — the tension is real and named (§VI F1) — but the replicator
amplification (§II.3) means the low-pass only needs to create a fitness
*gap*, not a clean separation. Sweep {1, 4, 16} registered in J1/J2.]**

## I.5 Phase-sensitive binding (J1): why the metric rescues the ring

RING_DNLS A.2/A.3 separated the medium's two killers, and **both enter
through the contact links**: (1) bulk leak = accretion crowding — Em
exported through new bath contacts; (2) pure contact dephasing — the bath
as a random-phase reservoir walks one cycle edge through π (kills the
winding even at keep ≥ 0.99). Both scale with the contact hop angle
`τ_contact ∝ w_contact`.

At a ring cell, the renorm competition is between its ~2 cycle links
(persistent |Ĵ| ≈ 1.4) and its bath contacts (|Ĵ| → drain ~10⁻²–10⁻³ plus
the attenuated noise residue). Replicator dynamics with that fitness gap
drives the cycle links toward the ×5 clamp and the contacts toward the
0.2× floor within a few digging e-folds (e-fold time `1/(k_J·|Ĵ|)` ≈ 0.7
t.u. at k_J=1 for the cycle links). Consequences, in order of confidence:

- **Leak cut ~5×** (τ_contact ∝ w → 0.2): direct.
- **Dephasing time stretched ≥5×**: the contact noise *amplitude* on the
  cycle phase scales with τ_contact; if the slip is a diffusive
  first-passage through π, the slip time scales as 1/τ_contact² — up to
  ~25×. Claimed range 5–25×; the bar (J1) asks only for the low end.
- **The self-insulated ring is an emergent impedance mismatch** — the
  XSEC "inert load = reflector / impedance defect" phenomenology,
  now *produced* by the metric rather than found.

The A.2 baseline dies at t ≈ 8–16; a ≥5× stretch is survival ≥ 40–80 t.u.
— i.e. R2@80 climbs off the floor. On the *lit* bath (M4 geometry) the
door_g=0 ring already self-heals via accretion + topology; there the
metric's job is to keep the healing while cutting the noise (§IV J1a).

**[RESUME 1: MEASURED GREEN — R2@76 0.774 vs 0.131; Em_min@78 0.471 vs
0.0007 (no slip sites open). §VII.2.]**

**Why this does not re-run the A.3 non-monotone failure:** q_detune
throttled contacts through a *resonance barrier* that a draining voice
destroys by retuning (the drained-cell runaway, qd3.6 worse than law).
The metric's throttle has the opposite sign of feedback: a draining
contact's flow *falls* as its weight falls (w → 0.2 cuts the drain that
feeds Ĵ), so the starvation is self-limiting on the safe side — it
converges to the floor, not to a runaway.

## I.6 Antagonism release (J2): bandwidth against the trap, on the track only

**[RESUME 1: this section's mechanism (bandwidth vs trap) is NOT refuted,
but its premise — that the bare tilted blob supplies the digging current —
is measured FALSE. The corrected reading is §VII.3; the criterion algebra
below stands and now applies to the composed movers C2/C4.]**

The B.3 measurement: a tilted packet condenses harder, advances ~1.3
cells, and stalls — the DNLS mobility/trapping tradeoff. The face-3 §II.3
criterion: self-trap when the detuning gap beats the hop bandwidth,

```
w2·q·Δx/(1+q·x̄)²  ≳  z·κ̄ ,      κ = τ_s/dt (the per-link hop rate)
```

The metric attacks the right-hand side, *selectively*: the mover's forward
current digs the along-track links (`κ → κ·w`, up to ×5), so the effective
trapping threshold in q scales up by w on the track. B.1 located the bare
threshold at `q_c ∈ (0, 0.6)`; on a fully dug track the effective `q_c`
rises to ≳ 3.0 — the law value 1.2 moves to the *mobile* side with ~2.5×
margin. Meanwhile:

- **Static structure keeps its binding.** A parked hoard has net J ≈ 0;
  its weights renorm back toward 1; the trap criterion is unchanged. The
  brief's asked-for asymmetry — "the trap dissolves for coherent motion
  while persisting for static structures" — is not a bolted-on gate; it
  falls out of the detector (§I.4).
- **The mover self-collimates.** At a track cell the renorm takes budget
  *from the transverse links* (they carry no persistent J), so transverse
  hopping weakens: the packet is confined side-to-side by its own motion —
  a self-dug waveguide. This attacks the OTHER half of the item-4
  diffraction failure (the localized winding diffracted because the linear
  channel spreads isotropically).
- **Group velocity** rises with w on the track: `v_g = (2κ̄·w)·d·sin(k·d)`
  (§II.4) — the J2 COM bar has ×5 headroom over the baseline speed that
  produced 1.3 cells before the stall.

The timing question — can the track be dug before the stall? — is the
design's sharpest quantitative risk (§VI F1): the stall is at t ≈ 4, the
digging e-fold at k_J=1, |Ĵ|≈O(1) is ~1 t.u., the memory build-up is
τ_J-limited. **[RESUME 1: answered NO for the bare blob at τ_J ∈ {1, 4} —
and the deeper reason is per-link dilution, not timing alone. §VII.3.]**

## I.7 The retention wall and the door pair's phase cross

The brief asks whether the J-metric helps the door pair's phase survive
the mode crossing (M4). Two answers, kept separate:

- **Directly: no, by design.** The metric touches only within-mode spatial
  transport. The door — timing, grain, gates, settlement — is untouched.
  The door stays special (constraint honored; nothing here rotates,
  weights, or re-times any conversion).
- **Indirectly: yes, measurably.** M4's settlement composes χ against a
  ψ_m whose phase the bath has walked; `door_avg` fixed the sampling error
  (0.24 → 0.32) but cannot fix a reference that is itself decohering. The
  J-metric's job is precisely to keep the spatial phase field ordered
  (self-insulation, §I.5). Registered prediction: at `j_metric > 0` the
  M4′ gap closes, because the beat-integrated D finally accumulates
  against a stationary matter phase. This is composed bar C1 (§V, sharpened
  §VII.4).

## I.8 What the metric is NOT (constraints check)

- **Not a background.** The weights are state on the living roster's
  slots, born at 1 on link rebirth (`freecell.c:1148` precedent), dying
  with their links, redistributed under a conservation (the per-cell
  budget). Nothing persists-and-is-merely-revalued outside the dynamics.
- **Not an imported field.** No new species, no new carrier: w_s is
  channel geometry (the same ontological class as sbed, sA, sd — the
  roster's shape), carries no energy, has no pitch, is not transported.
- **Not a feedback/bias on wants.** The v92-closed road was additive bias
  terms injected into want formation. The metric multiplies a *unitary
  angle*; it cannot create or destroy energy (theorem, §II.5), cannot
  amplify a cell's total outflow capacity (zero-sum budget), and detects
  rather than injects phase (it reads |Ĵ|; it never writes th2).
- **Not sum-then-square.** Pairwise hops only, unchanged.
- **The door stays special.** Untouched (§I.7).

---

# PART II — Fully derived symbolic equations

## II.1 The metric dynamics

Per slot s (live, canonical orientation i→j), with F_s the actual moved
dense flow (§I.3), the law — stated as the exact FLOW mirror plus one
optional healing term:

```
memory:   dĴ_s/dt   = (F_s − Ĵ_s)/τ_J                          (one-pole low-pass)
growth:   dw_s/dt   = k_J·|Ĵ_s|·w_s  −  (w_s − 1)/τ_rel        (τ_rel = ∞ default: off)
budget:   w_s ← w_s·√(f_i·f_j),   f_i = n_i / Σ_{s'∋i} w_s'    (zero-sum per-cell renorm)
clamp:    w_s ∈ [0.2, 5]
use:      τ_s = amp_nat·base_s·G_s·w_s                          (the pass-U hop angle;
                                                                G_s = gsym gated / 1 linear)
```

Discrete-time form (mirroring `freecell.c:1382–1408` clause for clause):
`Ĵ_s += (dt/τ_J)·(F_s − Ĵ_s)`; `w_s *= 1 + k_J·dt·|Ĵ_s|`; renorm; clamp.
All quantities dt-correct to first order (κ = τ_s/dt is dt-free; the
low-pass and growth carry explicit dt) — and all J-campaign runs use the
adopted Strang schedule (`hop_order=1`), which is also the dt-invariance
fix of record.

The healing term is **off by default** (exact FLOW mirror; fewest new
decisions). Its physical case — g4 says footprints need throughput; a
metric maintained only by persistent current would heal its scars when the
current dies — is recorded, and a finite-τ_rel arm is registered in J5.
Note the roster already heals scars its own way: geometry churn recycles
links, and a reborn link has no memory (w=1, Ĵ=0).

## II.2 Where it acts

Only the pass-U hop angle. Not the additive want path (amp_nat=0 ⇒
j_metric inert), not pass F (the field does not dig — flagged §VI F5), not
pass S (the space current has its own adopted bed law; the two metrics are
siblings, not shared state), not the door (§I.7). The metric therefore
engages exactly the channel whose transport/binding/retention triangle
(RING_DNLS synthesis 3) it is designed to break.

## II.3 The replicator structure (why a modest contrast wins)

Fix a cell i, let `u_s = w_s/n_i` be its budget shares. Growth +
per-cell renorm give, in the continuous limit,

```
du_s/dt = k_J·u_s·(|Ĵ_s| − ⟨|Ĵ|⟩_i) ,    ⟨|Ĵ|⟩_i = Σ_{s'∋i} u_{s'}·|Ĵ_{s'}|
```

— replicator dynamics on the link simplex with fitness |Ĵ_s|. (The
kernel's √(f_i·f_j) pairwise renorm realizes this approximately, exactly
as the bed law does; the deviation is the same as the adopted law's and is
bounded by the clamps.) Consequences:

- Share ratios evolve as `u_a/u_b ∝ exp(k_J·(|Ĵ_a|−|Ĵ_b|)·t)`: any
  *persistent* fitness gap grows exponentially. The low-pass does not need
  to null the incoherent current, only to attenuate it below the coherent
  one (§I.4's 3.3× at default τ_J is ample).
- Uniform fitness (an isotropic bath) is a fixed point at u = 1/n_i: the
  bath digs nothing *in expectation* — the renorm cancels the common mode
  of |Ĵ| exactly (this is why growth-on-|noise| does not inflate a bath:
  the bias is uniform and the budget removes it; what remains is a slow,
  clamped random walk, §II.6).
- The clamps bound the winner-take-most endpoint: no link exceeds ×5 or
  falls below 0.2× — the metric can focus, never disconnect or explode.

## II.4 Dispersion and group velocity

On a uniform chain with per-step hop angle τ per link (small-τ,
linearized channel), pass U is the discrete Schrödinger coupler with
hopping rate κ = τ/dt:

```
ω(k) = w2e − 2κ·cos(k·d)          v_g(k) = 2κ·d·sin(k·d)
```

(the L1-A reframe of record: v_g = +2J̄_d·sin(kx·d) > 0, measured by fd).
With the metric, κ → κ·w_track on dug links:

```
v_g(k) = 2κ̄·w_track·d·sin(k·d)    — up to ×5 the bare speed
```

and the transverse hopping rate falls to κ̄·w_⊥ with w_⊥ < 1 by renorm:
the dispersion surface anisotropizes around the track — a waveguide. The
self-trap criterion (face-3 §II.3) with the metric in place:

```
trap  ⇔  w2·q·Δx/(1+q·x̄)²  ≳  Σ_{s∋i} κ_s·w_s ≈ z_track·κ̄·w_track + z_⊥·κ̄·w_⊥
```

— on the track the bandwidth term is up-weighted ×w_track: the antagonism
(the corner that binds does not translate) is released *along the track*,
and only there. Static hoards (Ĵ ≈ 0 ⇒ w ≈ 1) recover the bare criterion:
binding persists exactly where the brief asks it to.

## II.5 The conservation theorem (energy), with the metric-sector corollary

Pass U moves energy only by Givens rotations; a Givens rotation of ANY
angle preserves the two-cell norm exactly (`|ψ_i'|²+|ψ_j'|² =
|ψ_i|²+|ψ_j|²` for all τ — the v93 theorem). The metric enters the update
**only through the angles** τ_s; the weights are not energy, appear in no
ledger, and are touched by no click. Therefore:

> **Corollary (metric-independence).** The conservation theorem holds
> uniformly over the entire metric sector: for ANY weight configuration
> {w_s(t)} — however it evolves, clamps, or renormalizes — ΣEm changes
> only by FP roundoff in pass U. `j_metric` cannot move the drift floor,
> by construction.

This is the face-3 §II.7 pitch-independence corollary transplanted to the
metric sector: bar J3 can only catch coding errors (an angle computed from
a stale weight is still an angle; the one implementable mistake class is
accumulating F_s inconsistently with the hops actually applied — which
touches the *meter and the metric*, never the energy). The τ ≤ 0.5
integrator guard is retained AFTER the weight multiply (a ×5 weight on a
deep-overlap link can reach it sooner; the clamp saturation count is a
standing diag, §III.6). **[RESUME 1: GREEN — drift = 0 at j_metric=1.]**

## II.6 The digging floor (why the bath cannot dig)

Engagement requires the growth rate to beat the run: digging matters when
`k_J·|Ĵ|·T ≳ 1`, i.e.

```
|Ĵ|_floor ≈ 1/(k_J·T)
```

At k_J = 1, T = 80: floor ≈ 1.2×10⁻². Measured currents against it:
ring circulation 1.4 (**110× above** — digs in ~1 t.u.); mover forward
current O(1) *in total* — **[RESUME 1: but per-link ~10⁻²–10⁻³ after
dilution over the packet's support: at or below the floor. The J2 killer;
§VII.3]**; hoard diffusive leak 10⁻⁴–10⁻⁵ (**100–1000× below** — never
engages, so the metric cannot accelerate R1's slow diffusion); bath noise
residue after low-pass + renorm common-mode cancellation: zero-mean,
variance-driven random walk with step `k_J·dt·σ_res`, clamped — J4's
dark/lit-bath stats bars measure that it stays quiet.

## II.7 The space budget (the metric's own conservation)

The renorm identity `Σ_{s∋i} w_s ≈ n_i` (exact at the fixed point of the
pairwise √(f_i f_j) iteration, approximate per step, clamp-bounded) is the
metric's conservation law: **total channel per cell is conserved; current
only redistributes it.** This is the implementable content of "space is
convertible but conserved," and it is the structural anti-ignition
guarantee: no configuration of currents can raise a cell's total
conductance, so no positive transport-feedback loop closes through the
metric. (The same argument, for the same law, on the other current, is the
ASYM campaign's zero-sum result of record.) J4 prints the identity's
residual as a standing meter.

---

# PART III — Numerical simulation method (the kernel face)

## III.1 Knobs (all in BOTH kernels; struct + default + parser + header)

| knob | type | default | meaning |
|---|---|---|---|
| `j_metric` | double | **0** | k_J, the digging gain. 0 = byte-inert (no state touched, no arithmetic executed). |
| `j_tau` | double | 4 | τ_J, the memory time constant (t.u.). Inert while j_metric=0. |
| `j_rel` | double | 0 | healing: relaxation of w toward 1 with time constant j_rel; 0 = off (FLOW mirror). |

Clamps fixed at [0.2, 5] (the adopted bed constants — shared, not new
knobs). `j_metric > 0` requires `amp_nat > 0` to do anything (the metric
weights pass-U angles only).

## III.2 State (file-scope C / Sim-struct Go; slot-indexed, NLMAX)

```
sjw[s]    double, init 1.0     — the link weight (the metric)
sjnet[s]  double, init 0.0     — the low-passed dense flow memory Ĵ_s
sjfl[s]   double, zeroed/step  — this step's actual moved flow accumulator
```

(**[RESUME 1: the implementation's names are `jbed`/`jflow` — record the
actual names in IMPLEMENTATION notes; the clauses are what is specified,
not the identifiers.]**)

Slot lifecycle: on link rebirth reset `sjw=1, sjnet=0` — same site as the
bed reset (`freecell.c:1148`, "a reborn link has no bed"; Go
`geom.go:128`). Weights die with their slots; no state survives the
roster.

## III.3 Operations, in schedule order (each guarded `P.amp_nat>0 && P.j_metric>0`)

1. **Pass 2 (where `shau_[s]` is computed, :1590–1593):** after the
   existing `taub` computation and BEFORE the 0.5 clamp:
   `taub *= sjw[s];` then clamp. Also zero `sjfl[s]` here (the loop
   already visits every slot). At `j_metric=0` neither line executes —
   byte-inert by branch, per kernel convention (never rely on ×1.0 IEEE
   identity).
2. **`apply_uhop` (:1296):** after the writeback, accumulate the actual
   transfer: `sjfl[s] += (|ψ_j'|² − |ψ_j|²);` — the identical difference
   the p1 meter already forms; under Strang both half-hops accumulate,
   giving the true per-step moved flow. (Guarded; at j_metric=0 the
   function is byte-identical to today.)
3. **Post-pass-U metric update (new block after the pass-U writeback,
   before pass 3):** for every live slot (including zero-flow slots — the
   memory must decay where flow stops, exactly as bednet does):

```c
double kb = P.j_tau > dt ? dt / P.j_tau : 1.0;
sjnet[s] += kb * (sjfl[s] / dt - sjnet[s]);
sjw[s]   *= 1.0 + P.j_metric * dt * fabs(sjnet[s]);
if (P.j_rel > 0) sjw[s] += (dt / P.j_rel) * (1.0 - sjw[s]);
/* zero-sum per-cell renorm — the bed law's clauses verbatim: */
/*   f_i = n_i / Σ sjw over live slots of i;  sjw *= sqrt(f_i*f_j);  */
/*   clamp [0.2, 5]                                                  */
```

   The renorm reuses the bed's `bedf_`-style scratch (a separate `sjwf_`
   array — do NOT share bedf_, the two metrics may run together).
   Stale-by-one discipline: the weights used in pass 2 this step were
   updated after last step's pass U — the same one-step lag as x, Em,
   flload, and the bed itself; no new race.

## III.4 Byte-inertness and battery

- `j_metric = 0`: no branch executes, no new arithmetic touches any
  existing array — byte-for-byte the v94 surface, hence byte-for-byte
  laws_V3a at all other knobs 0. **Full battery 87 GREEN, both kernels,
  before commit** (the ratchet).
- `amp_nat = 0`: pass U is skipped entirely; j_metric>0 is inert (the
  battery surface never sees it).
- Conserve-diff against the standing byte-inert reference before the
  battery, per v93 practice.

## III.5 C / Go mirror parity

Go: `P.JMetric, P.JTau, P.JRel` (`params.go`: struct + defaults + parser +
header print); `s.sjw, s.sjnet, s.sjfl` (`sim.go` alloc, `geom.go` rebirth
reset, `step.go` pass-2 multiply + `applyUhop` accumulate + post-pass-U
block). Gotchas of record: `math.Sincos` returns (sin, cos); no boolean
shortcuts (`P.JMetric != 0` style); mirror clause order exactly (the
renorm's iteration order matters for byte-parity). Parity grading per the
standing X.6 language: byte-identity on the battery surface; at
`amp_nat > 0` all campaign data are Go-only with the C side
tolerance-graded (pre-existing v93 debt, unchanged by this face).

## III.6 Meters (diag only, no physics)

- Header line: `# v95 J-metric: j_metric=%g j_tau=%g j_rel=%g`.
- Per-diag `# JMET t= wmin= wq1= wmed= wq3= wmax= ngrown= nclamp_tau=`
  (bed-meter clone :2921–2933 + the count of slots whose τ hit the 0.5
  guard).
- `# JMAP t= i j w jnet` dump behind `snap_every` (BEDMAP clone) — the
  track/insulation visualizer J1/J2 grading needs.
- The renorm-identity residual `max_i |Σ_{s∋i} w_s − n_i|` printed in the
  RESULT block (J4's structural meter).

---

# PART IV — Acceptance bars (pre-registered; the ratchet governs)

**[RESUME 1: results in — J1 GREEN (exceeds), J2 RED (diagnosed), J3
GREEN, J4 partial GREEN. Per-bar markers below; grading detail and the
re-scoped composition series in PART VII.]**

**Invariant surface first:** battery 87 GREEN at `j_metric=0`, both
kernels, before and after the edit; bars sharpen by measurement, never
soften. All strength runs: `hop_order=1` (Strang, adopted), linear channel
(`amp_logate=1`) unless the gated arm is named, Go kernel (campaign
kernel of record at strength, per X.6).

- **J0 (inertness + parity).** Battery 87 GREEN ×2 kernels; conserve-diff
  byte-identical at `j_metric=0`; battery-surface C≡Go byte-identical.
  **[RESUME 1: not explicitly reported — confirm and record (§VII.5).]**
- **J1 (phase-sensitive binding).** The brief's bar, with provenance
  pinned per arm:
  - **J1a (the brief's number):** lit-bath m=+2 C6 ring, M4 config
    (`exp=ring ring_n=6 bath=1 seed_mw=2 amp_nat=2 amp_logate=1
    hop_order=1 noise_amp=1.0 k_rad=0.05 p_rad=4 wf_on=1 q_detune=1.2
    L=16 T=80 seed=20260802`), `door_g=0.5 door_avg=1 j_metric=1`:
    **R2@80 ≥ 0.5**. **[RESUME 1: the run executed was this geometry at
    `door_g=0` — GREEN at 0.774 vs 0.131; the composed door_g=0.5 arm
    is now C1. Bar sharpened to the sustained form, §VII.2.]**
  - **J1b (cold-bath, the clean mechanism arm):** RING_DNLS A.2 config
    (`rb6_uni_m2` + `hop_order=1 j_metric=1`): winding survival ≥ 40 t.u.
    (4× the measured 8–16) AND R2@80 ≥ 0.5 (vs 0.22; noise floor 0.41 —
    only R2 with W locked counts, per the standing A-section rule).
    **[RESUME 1: not yet run — remains registered.]**
  - **J1c (stretch, not gate):** the item-4 localized blob vortex
    (`seed_mw=2` blob, no cycle support): R2d@80 ≥ 0.5 vs 0.07.
  - Sweeps inside J1: `j_tau ∈ {1, 4, 16}`, `j_metric ∈ {0.3, 1, 3}`;
    JMAP must show the insulation (contact w → floor, cycle w → ceiling)
    for the mechanism claim, not just the number.
- **J2 (antagonism release).** B.3 config (`exp=blob bath=1 L=16 dt=0.02
  T=80 sigma=2.5 amp=0.5 kx=1.1 q_detune=1.2 amp_nat=2 amp_logate=1
  hop_order=1 seed=20260802`) + `j_metric=1`: **all-cell COM advance ≥ 5
  cells** AND Em_max@end ≥ 50% of the j_metric=0 value. **[RESUME 1: RED
  at j_tau ∈ {1, 4} — COM ≈ 0, anchored. Diagnosis §VII.3: no per-link
  current source; the pre-registered failure semantics ("record and
  consult") triggered — this resume is the consult; the bar migrates to
  the composed movers C2/C4.]**
- **J3 (conservation).** `|drift_rel| ≤ 1e-13` on live matter, T ≥ 1000,
  j_metric=1. **[RESUME 1: GREEN — drift = 0 (the §II.5 corollary held;
  the run only certifies the implementation).]**
- **J4 (anti-ignition).** (i) Dark bath at `j_metric=1`: byte-identical
  where Em≡0 keeps pass U silent. **[RESUME 1: GREEN — structural
  byte-identity confirmed (no dense hops → no flow → weights stay 1).]**
  (ii) Lit bath stats within noise; `wmax` off the ceiling;
  renorm-identity residual at floor; no secular weight drift T ≥ 500.
  **[RESUME 1: open.]** (iii) The digging-floor check: hoard-only run
  (B.1 law corner, kx=0), w quartiles ≈ 1 throughout. **[RESUME 1:
  open.]**
- **J5 (endurance — RECKONING's raised bar, adopted verbatim).**
  **[RESUME 1: open — still gates ANY stability claim; upgraded to C5
  (the composed object's endurance) alongside the J1-winner static ring
  at T = 3000.]**

Run order (updated): confirm J0 → J4(ii,iii) → J1 sweeps + time-series
(§VII.2 NEEDS) → C1 → C2a → C2b/C4 → C3 → C5/J5.

---

# PART V — The composition with the door pair (the v95 architecture)

## V.1 The division of labor, measured vs designed

| direction | mechanism | status |
|---|---|---|
| **mode (conversion)** | door pair: pass P Givens + settlement in atoms; `door_avg` synchronous detection; `chi_load` self-detuning | **measured, won**: DS-G3 r=+0.684 vs −0.259; anti-runaway 17 vs 417; M6 park approach 10×; conservation FP floor; anti-ignition 5+1 layers |
| **space (transport)** | J-metric: dynamical Hodge star on pass-U angles; synchronous detection on the dense current | **[RESUME 1: binding/retention half measured, won (J1: R2 0.774, Em_min ×14–×670). Transport half needs a current generator — the composition (C2) or topology (C4). §VII.]** |

Structural orthogonality: the pair touches the mode fiber (ψ_e–χ–ψ_m) and
never a spatial angle; the metric touches spatial angles and never the
fiber. Shared state is only what physics shares (Em, th2, x). No new
coupling term exists between them — composition is by superposition of
passes, which is why the J-bars ran at `door_g=0` first and compose after.

## V.2 Why the composition closes v94's open wall

- **Z.2 (the pair confines; the metabolic wave does not propagate):** with
  the metric, propagation is *within-mode* — matter moves as matter on dug
  tracks; the pair converts in place. The grain floor (a front cell needs
  χ² ≥ ε before its first click) no longer blocks motion, because motion
  no longer requires conversion at the front.
- **M4/M4′ (phase across the door):** `door_avg` removed the sampling
  error; the metric supplies the stationary phase reference (§I.7). The
  two fixes are complementary halves of one coherence chain:
  field → (pass P, coherent) → χ → (door_avg, integrated) → ψ_m →
  (J-metric, spatially protected) → ψ_m elsewhere.
- **M5′ (translation-by-conversion):** stays registered on the pair's side
  (trailing-edge outtake), now runnable without the erosion confound
  (`door_avg=1`) *and* — if C2a is green — with the trap released (J2's
  mechanism powered by the metabolic current). If M5′/C2 still stalls,
  Z.4's verdict stands in its strongest form and C4 (topology) is the
  remaining road.
- **The v95 target object** (the brief's "metabolizing, translating,
  phase-retaining object"): a q_detune self-trapped core (mass — B.1),
  carrying a winding on closed support (identity/phase — A.1), on its
  self-dug track (mobility — C2/C4), self-insulated from contact noise
  (retention — J1, measured), metabolizing through the door pair
  (intake/outtake — faces 1–2) at the chi_load park (persistence — M6).

## V.3 Composed bars — SUPERSEDED by §VII.4 (re-scoped after J1/J2)

The original C1–C3 registrations stand in spirit; §VII.4 sharpens C1 to
the measured 0.774 line, splits C2 into the staged C2a/C2b, adds C4 (the
topological mover) and C5 (composed endurance), and pins configs.

---

# PART VI — Flags: load-bearing assumptions and failure modes

**Assumptions (inline A1–A3, collected; §VII.6 re-grades):**
- **A1** Per-cell zero-sum budget = the anti-ignition structure,
  inherited from FLOW's measured record on the space current, assumed
  transferable to the dense current (J4 measures).
- **A2** Unsigned fitness |Ĵ| (channels are direction-blind); the
  vector-credit refinement is explicitly not taken.
- **A3** τ_J = 4 default — the one new scale; chosen against measured
  timescales (contact beat 7.7 t.u., stall ~4 t.u.); swept, not trusted.
- **A4** The clamp ceiling ×5 (carried from the bed) is assumed
  sufficient for J2's release (derived margin 2.5× at the law corner —
  but the (q, Δx) threshold surface behind that derivation is
  single-point-calibrated; face-3 U1(a)'s two-parameter sweep remains
  open and would harden it).
- **A5** All strength results are Go-kernel; the C-side amp_nat>0
  tolerance debt is pre-existing and unchanged (X.6 language governs).

**Failure modes, most likely first:**
- **F1 — the timescale tension (named in §I.6).** Digging (needs fast
  response) and noise rejection (needs slow memory) pull τ_J in opposite
  directions. **[RESUME 1: resolved-fatal for the BARE blob (J2 RED at
  τ_J 1 and 4) — but the root cause is per-link dilution (§VII.3), which
  no τ_J fixes; the tension is moot for the bare blob and re-arms for the
  composed movers.]**
- **F2 — DC drain digging.** Synchronous detection cannot distinguish
  coherent circulation from a steady drain (both DC). **[RESUME 1: did
  not materialize at J1 — keep improved (0.31 → 0.49) and Em_min rose;
  the §I.5 competition analysis held. Provisional pass; watch in C1's
  door-traffic books.]**
- **F3 — accretion-channel amplification.** Condensation transients carry
  real currents into forming hoards; if those exceed the digging floor,
  the metric deepens accretion. J4(iii) is the direct probe (open).
- **F4 — scar anisotropy.** With j_rel=0 a dead structure's track
  persists until geometry churn recycles the links. Bounded by clamps and
  healed by roster churn; J5c measures the healing alternative; adoption
  judgment is the user's.
- **F5 — the field does not dig (scope cut).** Pass F keeps its fixed
  metric. If the J-metric succeeds for matter, the symmetric question
  (does light dig? — lensing, the g-lane) opens; out of v95 scope,
  flagged for the record alongside face-3 §IV.4's space-register
  speculation.

---

# PART VII — Resume update 1 (J1–J4 landed, 2026-08-09)

**Source:** orchestrator resume-1 results (J1 lit-ring at door_g=0; J2
blob; J3; J4 dark arm). Same session. Where this part supersedes a bar or
a grading above, this part governs.

## VII.1 Scoreboard

| bar | grade | headline |
|---|---|---|
| J1 (phase-sensitive binding) | **GREEN, exceeds target** | R2@76 **0.774** vs 0.131 (target 0.5); Em_min ×14 @76, ×670 @78; keep 0.49 vs 0.31 — best in-medium phase retention ever measured in the programme |
| J2 (antagonism release, bare blob) | **RED** | COM ≈ 0 at j_tau ∈ {1, 4} — anchored; diagnosis §VII.3 |
| J3 (conservation) | **GREEN** | drift = 0 at j_metric=1 (the §II.5 corollary certified) |
| J4 (anti-ignition) | **GREEN (dark arm, structural)** | byte-identical: no dense hops → no flow → weights stay 1; lit-stats + hoard-floor arms open |
| J0, J1b/c, J4(ii,iii), J5 | **open** | §VII.5 |

## VII.2 J1 graded: the self-insulation mechanism is real

**What the numbers say.** The metric starved the ring→bath contacts
(incoherent, |Ĵ| ~ 10⁻²) and grew the cycle links (coherent circulation,
|Ĵ| ≈ 1.4) — the §I.5 competition, operating as derived. The decisive
mechanism number is **Em_min**: 0.183@76 / 0.471@78 vs 0.013 / 0.0007
without the metric. Phase-slip sites are support holes (A.2/A.3: a voice
draining toward Em ~ 10⁻³ is where the winding dies); the metric keeps
every voice fed, so no slip site ever opens. R2 = 0.774 follows. This is
the best in-medium retention on record (prior: cold-bath ~0.02 dead by
t≈10; hot-bath self-heal 0.62 as a late-time spike). The vacuum reference
(A.1: 0.988, W exact) remains the ceiling.

**Grading discipline (three sharpenings, per the ratchet):**

1. **The baseline correction.** M4's door_g=0 "0.62@80" and this run's
   "0.131@76" are the *same fluctuating series* (M4's own table: 0.13@60,
   0.226@78, 0.62@80) — single-snapshot R2 on a 6-cell ring is violently
   noisy (random-phase floor ≈ 0.41). The honest baseline is "R2 at or
   below the noise floor with intermittent spikes." The j_metric=1 arm's
   0.774 must therefore be certified as *sustained*, not sampled. **J1
   bar sharpened (supersedes the @80-snapshot form): median R2 over the
   last quarter of the run ≥ 0.5, AND W locked at +2 for the certified
   span, AND Em_min > 0.1 throughout it.**
2. **NEEDS — the W(t) column.** The standing A-section rule ("only R2
   with W locked counts as retention") requires the winding number
   series for the j_metric=1 arm. Report it next resume; if W walked
   while R2 stayed high, the claim demotes from retention to re-formation
   (still a result — but a different one).
3. **NEEDS — the JMAP certificate.** The registered mechanism requirement
   ("JMAP must show the insulation, not just the number"): contact-w
   quartiles at the floor, cycle-w at the ceiling, at 2–3 snapshot times.
   The narrative asserts it; the numbers should be on file. Plus 3 seeds
   (the M4 lesson: one seed's late-time R2 is not a measurement).

**F2 (DC-drain digging) did not materialize:** keep *improved* (0.31 →
0.49) — the drain competition resolved toward insulation, as §I.5 argued.
Provisional pass; re-examine in C1 where door traffic adds a second drain
channel.

## VII.3 J2 diagnosed: the metric is an amplifier, not a generator

The failure localizes to one assumption, and it is **not** the metric's
machinery (which J1 just confirmed end-to-end — detection, competition,
renorm, clamps all behaved as derived). It is §I.6's premise that the
bare tilted blob supplies its own digging current. It does not, for three
compounding reasons:

1. **Dilution.** The measured total forward current (fd ~ O(1)·t.u.⁻¹,
   L1-A) is spread over the packet's support (PR ≈ 456 cells, hundreds of
   links): per-link |Ĵ| ~ 10⁻²–10⁻³ — at or below the digging floor
   1/(k_J·T) ≈ 10⁻². The ring wins by concentrating |J| ≈ 1.4 on six
   links; the blob loses by spreading the same order of total current a
   hundred ways.
2. **Transience.** The B.3 stall (~4 t.u.) kills the source before any
   memory can consolidate it — confirmed by the j_tau=1 arm failing too
   (faster memory cannot remember a current that is already gone).
3. **Self-defeat.** As the packet slows, it condenses harder (B.3: the
   moving packet is the strongest condenser), further cutting the current
   that would have dug its escape. The chicken-and-egg is dynamically
   *stable* on the wrong side.

**The measured law, stated for the record:** *the J-metric amplifies
persistent coherent currents; it generates none.* Mobility therefore
requires a **current generator** + the metric as **consolidator**. The
substrate's generator inventory, with status:

| generator | persistence | per-link amplitude | status |
|---|---|---|---|
| topological circulation (winding on closed support) | permanent (topology) | O(1) | **measured — powers J1** |
| linear tilt (kx seed) | ~4 t.u. (transient) | 10⁻²–10⁻³ (diluted) | **measured insufficient (J2)** |
| metabolic conversion cycle (door pair) | sustained while fed | unmeasured | **the C2 hypothesis** |
| bath accretion (hot-bath feeding) | sustained while hot | unmeasured (isotropic — needs symmetry breaking) | recorded, no arm |

**The insulation/mobility relation, one level up.** J1's mechanism
(starve contacts) and J2's goal (open forward contacts) are not in
contradiction: the metric opens exactly the contacts that carry
*persistent coherent* current. A genuine mover's forward contacts do; a
stalled blob's don't. The metric is thus a **coherence filter on
mobility**: only phase-organized motion can open its own road, while
incoherent matter is bound and insulated where it stands. That is the
transport/binding/retention triangle (RING_DNLS synthesis 3) breaking
*selectively for coherent objects* — the design intent, now half-measured:
the binding/retention half is won (J1); the mobility half awaits a
generator (C2/C4).

## VII.4 The composition bars, scoped (supersedes §V.3)

**The composition hypothesis, stated as a mechanism chain (C2):** the
door pair is the current generator the bare blob lacks. Coherent field
with spatial structure (kx) → pass P loads χ phase-locked to the drive →
settlement at `door_avg=1` composes matter whose phase carries the
field's spatial gradient (the measured basis: DS-G3 r = +0.684 —
conversion preserves spatial information; this is the one measured
coherent e→m phase channel in the programme) → freshly settled matter
carries a forward phase tilt → pass-U forward current, **replenished
every beat by continued intake** (fixing transience) and **localized at
the object** (fixing dilution) → the metric digs → bandwidth rises →
translation before re-trapping. Pair alone confines (M5/Z.2); metric
alone anchors (J2); the claim is that the composition moves.

- **C1 (composed retention — do-no-harm at the new line).** J1 geometry
  (lit-bath m=+2 ring) + `door_g=0.5 door_avg=1 chi_load=1 j_metric=1`,
  3 seeds. **Bar:** sustained R2 (last-quarter median) ≥ 0.5 with door
  traffic > 0, and not below the pure-J line (0.774 single-seed; use the
  3-seed band) by more than the seed spread. Grades §I.7's prediction
  (the metric gives the settlement a stationary phase reference) and
  watches F2's second drain channel (door books: door_out must not
  ratchet — door_avg=1 required).
- **C2a (the generator meter — run FIRST, it decides the mechanism
  cheaply).** Blob + lit bath (`noise_amp=1`) + `door_g=0.5 door_avg=1
  chi_load=1 j_metric=1 kx=1.1`, T ≈ 80. **Measure, not translate:**
  per-link Ĵ quartiles on the packet's leading/trailing shells (JMAP)
  vs the digging floor 1/(k_J·T), and their persistence (≥ 3·τ_J).
  **GREEN = the metabolic cycle sustains per-link coherent current ≥
  the floor.** If C2a is RED, C2b is unreachable by this mechanism —
  record, skip C2b, and the road narrows to C4.
- **C2b (the metabolic mover — the deadlock breaker).** Same config,
  T ≥ 200: **COM advance ≥ 5 cells AND Em ∈ [50%, 200%] of start**
  (metabolizing, neither starving nor runaway — the M1 anti-runaway
  result says the pair defends the upper bound) **AND drift ≤ 1e-13.**
  Sweeps: `j_metric ∈ {1, 3}`, `j_tau ∈ {2, 4}`, `door_g ∈ {0.3, 0.5}`.
- **C3 (attribution 2×2).** The M5′ config at `door_avg × j_metric`
  ∈ {0,1}²: attributes translation-by-conversion vs transport-by-metric
  vs synergy-only. Unchanged from the original registration; run after
  C2 so the interesting cell is grounded.
- **C4 (the topological mover — the door-free alternative).** The J1
  winner (lit-bath wound ring, `door_g=0`, `j_metric=1`) + kx tilt.
  The winding supplies the persistent current (the generator that
  measured GREEN); the tilt supplies direction. Open question, stated
  honestly: circulation digs the *cycle*, not the road ahead — whether
  the superposed tilt opens forward contacts against the ring's own
  self-insulation is exactly what the arm measures. **Bar: COM ≥ 3
  cells with W held** = the programme's first mobile coherent object,
  door-free. (If both C2 and C4 are RED, mobility needs a generator
  the substrate does not yet own — that finding would gate v96.)
- **C5 (composed endurance = J5, RECKONING's bar).** Whichever object
  moves (C2b/C4) — or, if none, the static J1 ring — at **T = 3000**:
  non-transient (Em_max/PR/rms not relaxing to seed levels; W held;
  R2 sustained). **No "stable matter" or "mobile matter" claim before
  C5 is green** — q_detune's transient (R1) is the standing warning.
  Include the `j_rel ∈ {80, 200}` healing arm (J5c) here.

Run order: J0 confirm → J1 time-series/W/JMAP/seeds (VII.2 NEEDS) →
J4(ii,iii) → **C2a** → C1 ∥ C2b → C4 → C3 → C5.

## VII.5 Updated NEEDS (next resume)

- **J1 certification pack:** R2(t) + W(t) series, JMAP contact/cycle
  quartiles, 3-seed spread for the j_metric=1 lit ring. (Sharpened bar
  form in VII.2 — sustained median, W locked, Em_min > 0.1.)
- **J0 confirmation:** battery 87 GREEN both kernels at j_metric=0 +
  byte-inert conserve-diff; record the implementation's actual
  identifiers (`jbed`/`jflow`) and knob spellings in an
  IMPLEMENTATION note.
- **J4(ii) lit-bath stats + renorm-identity residual; J4(iii)
  hoard-floor probe** (the metric must NOT engage on static
  condensation).
- **C2a before C2b** (the cheap decider); then the VII.4 order.
- **J1b (cold-bath ring)** — the clean mechanism arm, still worth one
  run: it isolates the metric from accretion healing entirely.

## VII.6 Assumptions re-graded

| flag | status after resume 1 |
|---|---|
| A1 (zero-sum anti-ignition transfers to dense) | holding (dark arm structural; lit arm open) |
| A2 (unsigned fitness) | untested against a mover (no mover yet); C2b/C4 will exercise it |
| A3 (τ_J = 4) | J1 passed at default; J2 insensitive to τ_J (root cause was dilution) — the "one new scale" carries less load than feared |
| A4 (×5 clamp sufficient for release) | untested (no current source reached the clamp on a track); C2a measures the reachable Ĵ first |
| §I.6 premise (bare tilt current suffices) | **REFUTED — the design's one measured-false assumption; replaced by the generator/amplifier law (VII.3)** |
| F2 (DC drain digs) | did not materialize at J1; provisional pass, re-watch at C1 |

*End of resume update 1. The scoreboard, honestly: the metric's
detection-competition-consolidation machinery works exactly as derived —
on the one structure that owns a persistent coherent current, it
delivered the programme's best in-medium phase retention. What it cannot
do is create the current it amplifies. The v95 architecture's division of
labor is now measured from both sides: the pair converts but confines;
the metric consolidates but cannot initiate. The composition — the door
pair as generator, the J-metric as consolidator — is the single
registered road to mobile matter, C2a is its cheap decider, and C4
(topology as generator) is the door-free hedge.*
