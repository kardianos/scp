# STABILITY — tangent stiffness and modal theory of consonant networks

**v89 prestress program, theory side. 2026-07-28.**
Derived from the exact update rules of `v89/cellfab.c` (standing law table
`v89/battery/laws_V2g.cfg`), with no kernel runs: every number below comes
from the derivation chain in S1–S2 plus the runnable reduced computations
in `theory/ring_map.py` (exact per-step map, Floquet, kick spectroscopy)
and `theory/modal_2x2.py` (closed-form per-mode matrix). Outputs quoted
inline; rerun with `python3 theory/ring_map.py all` and
`python3 theory/modal_2x2.py`.

Anchors reproduced before anything was predicted: ring12 m=5 seed gates
fwd **1.0000** / back **0.1001** (measured 1.0 / 0.100); x(ring12,m=5) =
**0.32054**; x(ring6,m=3) = 0.12823; x(cube,π-rung,ā=1.586) = 0.3867;
x_skirt = 2Γ_m/(q(w₂−2Γ_m)) = **0.06173**.

---

## S0. Scope, method, and what is frozen

The reduced state is the set of degrees of freedom a seeded dense network
actually exercises: per-cell dense clock θ (kernel `th2`), dense store
`Em`, and per-link two-direction flight slots `lem`/`lph`. The map in
`theory/ring_map.py` reproduces `step_field()` pass by pass (line refs in
its docstring). Frozen, with justification:

* **Es / pass S** — identical cells at equal load exchange no space
  (pressure differences vanish); Es enters only through radii → the
  conductance product. Quoted as a systematic (seeded vs
  footprint-relaxed Es).
* **Plane normals / pass 7** — kappa_align and tumble move normals, not
  phases; they enter as the multiplicative plane factor `gpl` on every
  channel. All rates below scale with the single conductance product
  **Φ = k_dep·k_dep_m·geo·gpl**; mode *ratios* and all *signs* are
  Φ-independent. Spectra are quoted at the seeded-geometry gpl with a
  band up to the aligned value.
* **Field sector** — Ee = 0 on a dense net; roughness output is booked
  to a radiated ledger (the real kernel recondenses ~11% of it; M0).
* **Beat conversion / pass 6 doors** — with Ee = 0 and Em < cap both
  doors are shut. This is the measured M0 fact (bulk evaporation exactly
  0) falling out of the pass structure.

Numerical floor: Floquet baselines are periodic to |dz| ≲ 1e-3 (stable
nets after 1500 t.u. relax), giving an eigenvalue floor of about
±0.003–0.005 /t.u. on Re; the two exact n=0 neutrals (uniform clock
gauge, conserved mass) land at 0 to ≤3e-3.

---

## S1. The exact discrete dynamics (what actually moves)

Per step (dt = 0.02), for free cells (cflag 0), dense sector, standing
table (kappa_freq = kappa_reac = 0 — both bias blocks are dead code at
rate level):

**Pitch (pass 0+1, lines ~1462, ~1552).**
x_i = (Em_i + flload_i)/cap, with flload_i = ½·Σ_incident(lem_fwd+lem_bwd);
ω_i = w₂/(1+q·x_i). *Bound* energy only; passing field does not re-pitch.

**Wants (pass 2, ~1573).** Per link (i,j), per direction,
F_{i→j} = k_dep·k_dep_m·dt·geo·gpl·res(det)·g_{ij}·head_j·m_i^eff, with
g_{ij} = G(θ_i − ω_i d/C − θ_j), G(ψ) = ((1+cos ψ)/2)^p_gate (`gate_of`
~482), head_j = clip(1−(Em_j+Ee_j)/cap), m_i^eff =
√(Em_i·max(Em_j, mob_floor·cap)), res = Γ_b²/(Γ_b²+det²) via the comb
(1:1 on uniform nets, asserted). The outflow limiter (pass 3) caps a
cell's total outbound at 0.98·Em — inactive at steady flow (wants·dt ≪
Em).

**The clock θ (this answers "what moves th2").** Exactly three things:

1. **The advance** (pass 6, line 2069): θ_i += ω(x_i)·dt. Nothing else
   is deterministic drift.
2. **Deposit-side entrainment** (pass 4c, 1795–1819): fires *every
   step* on each arriving resolved deposit f: err = wrap(θ_src −
   ω_src·d/C − θ_recv), mix = f/(f+Em_recv+lock_floor·cap), θ_recv +=
   kappa_lock·mix·err. Rate-level torque: kappa_lock·F·err/(Em+lockf).
3. **Delivery-completion entrainment** (pass 5, 1879–1889): fires once
   per completed flight cycle (period τ = d/C per slot) with mix built
   from the accumulated take ≈ F·τ. Average torque again ≈
   kappa_lock·F·err/(F·τ+Em+lockf).

Total entrainment gain per unit err:
**g_e = kappa_lock·F·[1/(Em′) + 1/(Em′+Fτ)]**, Em′ = Em+lock_floor·cap —
about 2× the single-channel value at light load. Tumble (σ_tumble) and
kappa_align touch plane normals only: they are conductance noise on gpl,
never direct phase noise. th1 is a diagnostic cache (line 2052).

**Transport and flight.** Deposits debit the source instantly and ride
the slot; the slot delivers its integrated content once per τ
(integrate-and-dump, not smooth delay); residual sub-cycle phase rides
again (1891–1893). A full receiver stalls the slot without loss.
**Flight is load**: the staging is exactly half-in/half-out (deposit:
source x −f/2·cap⁻¹, receiver +f/2; delivery: the other halves), so a
*uniform* locked flow leaves every x invariant — measured in the map:
x stays 0.32054/0.12823 to 5 digits through the entire flight-fill
transient. The locked ring's pitch never moves while it circulates.
This is the engineered consequence of flload and the reason rung-true
seeds stay rung-true.

**Roughness (pass 5, 1859–1877).** On every dense delivery: det =
ω_i−ω_j (link orientation), R = 2|det|Γ_r/(det²+Γ_r²), demand =
take·rough_k·R, fired through the receiver's D→F credit.
R has a |det| kink at 0: roughness contributes **no smooth linear term**
at a uniform state — it is a dry-friction (|amplitude|-linear) loss on
detune perturbations, plus the vacuum-bleed tax (S4).

**Action atoms (quant_mode 2, `atoms_fire` 1411–1424).** Conversions
fire in whole atoms ε = A₀ω/2π (ring6: ε_F = 0.460) through a per-cell
credit clipped at 2ε. Consequences for the linearized response,
derived and then measured in the map (`theory/ring_map.py tax`):

* **Persistent** flows (steady bleed): mean rate exactly preserved;
  the atom adds a dead time ε/(demand rate) and a burst ceiling 2ε per
  delivery cycle. No hysteresis in the mean.
* **Transient** perturbations: **fully forgiven below one atom per
  receiver.** A kicked ring6 (n=3 load kick, T=200) loses **0.29–0.91**
  units in the continuous limit (quant_A0=0) and **exactly 0** (1e-15 =
  roundoff) under the standing table, at every kick amplitude up to
  0.1: the decaying detune's rough demand never accumulates ε before
  the mode dies. The credit pool persists indefinitely, so *repeated*
  kicks accumulate toward the atom — a genuine memory: the response to
  kick #2 depends on kick #1. Pre-registered discriminator P-S4 below.

---

## S2. The sign theorem (campaign-deciding)

**Setup.** Locked pair (or any rung link), 1:1, standing table.
Perturb antisymmetrically: δEm_i = +a/2 (heavy → flatter, δω_i = −γa/2
with γ = w₂q/(1+qx̄)²/cap), δEm_j = −a/2 (light → sharper).

**Lemma 1 (gate-argument sum; CONSONANCE VII.1 made exact).**
ψ_{ij} + ψ_{ji} = −(ω_i+ω_j)d/C ≡ −δ_rung (mod 2π). An antisymmetric
load perturbation preserves ω_i+ω_j to first order, so the pair *stays
on the rung* and the sum of the two gate arguments stays 0.

**Lemma 2 (the locked offset is symmetric).** With entrainment gains
K_i = K_j = K, the Adler lock condition gives the locked phase
v* = −Δ(τ + 1/K) (Δ = γa/2 the pitch detune), and the two gate
arguments settle at **ψ_{ij}* = −Δ/K, ψ_{ji}* = +Δ/K — exact
negatives**. Gates are even ⇒ g_{ij} = g_{ji} = G(Δ/K) *identically*:
the gate factor carries **no odd-in-detune flow at any entrainment
strength** (the rate-level version of s2's negative lemma).

**Lemma 3 (mobility is symmetric).** Above the sympathetic floor,
m_i^eff = m_j^eff = √(Em_i·Em_j): no direction from mob_sym.

**What survives: headroom.** The only odd factor left is
head_j − head_i = (occ_i − occ_j) = a/cap. Net flow i→j:

  **net = Φ·res·G(Δ/K)·√(Em_iEm_j)·(occ_i − occ_j)/1  > 0 for i heavy**

**VERDICT: HEAVY → LIGHT. The surviving rate-level dynamics is
restoring.** The heavy (flat) voice feeds the light (sharp) voice
because the light one has headroom; the loads re-converge, and since
load is pitch, the *ratio* re-converges too — the choir's correction
survives at rate level as a headroom effect, not a gate effect. This is
the same restoring sign s2 derived at amplitude level (interference
cross-flow feeds the sharp voice); the two mechanisms are independent
and agree. Passive detune runaway does not exist for π-rung pairs.

**Numbers** (`ring_map.py pair`, Floquet + analytic):

| object | rate of detune decay λ_h | source |
|---|---|---|
| pair d=1.25, m=1, gpl=1 | **−0.191 /t.u.** | exact Floquet |
| analytic −2Φ·G·S/cap, no delay | −0.135 | modal_2x2 |
| with flight-delay factor e^{−λτ} | −0.166 | modal_2x2 |
| direct kick a=0.01 | −0.168 | ring_map demo |

Direction confirmed heavy→light at every amplitude tried (0.005–0.05).
**Unlock threshold**: the pair unlocks only at dE_asym = **0.140** —
95% of its entire store (initial detune 0.29 ≈ 3Γ_b). Detune runaway is
unreachable by physical perturbation of a healthy pair.

**Second-order drift (the common mode).** The uniform-load direction is
*neutral* at first order (conserved mass, Jordan block with the
Goldstone clock). Its slow negative drift is not part of the linear
spectrum; it is the sum of (i) the vacuum bleed through
res_vac = Γ_m²/(Γ_m²+det_vac²) channels with the mob_floor readiness
(estimated 5e-5/link/t.u. at x=0.32, reproducing the measured
−0.05%/t.u. class), and (ii) the comma tax on *sustained* internal
detunes: loss rate ≈ (2·rough_k·γ/Γ_r)·(flux)·|δx_i−δx_j| — linear in
|amplitude| (nonanalytic), and **frozen below the action atom for
transients** (S1). Consequence: **the locked network is an attractor in
every internal direction, sitting on a neutral leaking ray.** Passive
seeds do not die of detuning; they slide down the load axis. The C1
plateau is not blocked by internal instability — it is blocked by the
environment-side bleed factors (res_vac·mob_floor·⟨G⟩·gpl_vac), which
is exactly the piston's side of the ledger (S5).

**One structural exception** — the wound even ring: see S3. Winding
moves the back gate off its flat top and buys a *chiral pump* that the
sign theorem's flat-gate premise excludes.

---

## S3. Modal spectrum — pre-registered predictions

Method: exact Floquet of the reduced map over an integer delivery
window (dim 4N state: θ, Em, 2 slots/link; lph is a rigid clock),
cross-validated by time-domain kicks + matrix-pencil fits and by the
closed-form 2×2 per-mode matrix (modal_2x2.py; agrees with Floquet to
6–20% for all π-rung slow/fast branches — the derivation chain is
sound). Rates in 1/t.u., frequencies in rad/t.u., at kernel dt = 0.02.

### ring6 (m=3, φ=π, x=0.1282, seeded gpl=0.0625) — the measured death object

All internal modes **real (overdamped): a kick relaxes, it does not
ring.** Two branches per |n|:

| mode n | load branch (slow) | phase branch (fast) |
|---|---|---|
| ±1 | **−0.0077** | −0.146 |
| ±2 | −0.0302 | −0.462 |
| 3  | −0.0416 | −0.635 |

2×2 analytic: −0.0094/−0.0283/−0.0377 and −0.138/−0.415/−0.553 ✓.
Conductance band (gpl 0.0625→0.25): slow n=1 −0.0077→−0.0099 (nearly
Φ-independent — pinned by the γ-pitch feedback); n=2 −0.030→−0.087;
n=3 −0.042→−0.130; fast branch scales ≈ ∝Φ (−0.146→−0.554).
**Pre-registered: slowest internal relaxation −0.009 ± 0.002 /t.u.
(τ_relax ≈ 110 ± 30 t.u.), insensitive to conductance.** Flight-slot
modes sit at Re ≈ −6.8..−7.3 (aliased at the window Nyquist; not
observable lines).

### cube (π-rung; the task's two readings)

Given (ā≈1.25, x≈0.39, φ=π) is rung-inconsistent (π-rung at ā=1.25
gives x=0.128; x=0.39 needs ā=1.586 — the H1 v2/v3 cube). Both
computed:

* **cube ā=1.586, x=0.3867** (H1 geometry; seeded gpl=0.037,
  geo=0.047 — the shell's channels are ~150× weaker than ring12's;
  edge overlap margin 2r−ā = 0.04 only, and the footprint-relaxed Es
  ≈ 0.71 **closes the direct edges entirely** (2r = 1.52 < 1.586) —
  the H1 cube lives at the strangulation boundary; D3b below).
  All modes real. Slow branch: **−0.0096, −0.0306, −0.0437**
  (adjacency classes μ = ±1, −3); at adapted gpl=0.15: −0.0088,
  −0.025, −0.039, −0.041, −0.084, −0.127. No underdamped line.
* **cube ā=1.25, x=0.1282**: slow branch triple-degenerate
  **−0.0373** (μ = +1, −1, −3), then −0.277, −0.559. All real.

**Pre-registered (H2-style kick, either cube): monotone two-time-scale
relaxation, NO breathing line in the dense sector.** If H2 spectroscopy
shows an oscillatory line, it is *not* the (θ, Em) sector — it must
come from the space mode (Es, pass S — outside this reduced model:
the s_k equalization time ~40 t.u. gives a candidate ~0.15 rad/t.u.
scale) or the field skin. That is a sharp discriminator between
"tension" pictures: dense-sector tension is overdamped; only the
space-cushion can ring.

### ring12 (m=5, φ=150°, x=0.32054) — THE WOUND RING IS A CHIRAL PUMP

The winding displaces the back-gate argument to ψ̄_b = +π/3 (g_b =
0.1001 — the measured 0.100), putting the back channel on the gate
slope, d ln G/dψ = −p_gate·tan(ψ̄_b/2) = −4.62. The exact map shows the
locked state is **linearly unstable** to long co-propagating load
waves:

| mode | growth Re | line Im (rad/t.u.) | period | validation |
|---|---|---|---|---|
| n=±1 | **+0.035 ± 0.010** | **∓0.64 ± 0.03** | 9.8 t.u. | kick fit +0.030 |
| n=±2 | **+0.053 ± 0.010** | **∓1.05 ± 0.05** | 6.0 t.u. | kick fit +0.059 |
| n=±3 | +0.002 (marginal) | ∓1.20 | 5.2 | +0.003 |
| n=±4 | −0.117 | ∓1.18 | — | damped |
| n=0 breathing | −0.126 | ±0.14 | — | damped |

Only **one propagation sense grows per |n|** — the waves co-propagate
with the forward transport (against the phase-winding sense): **the ±n
degeneracy is lifted. This chiral line splitting is the first
spectroscopic signature of the closure integer** (the B1 "chirality
signature" made concrete).

**Attribution (surgical probes on the exact map + delayed 2×2):**
freezing g_b at 0.1001 kills the instability (all n≠0 damped ≤ −0.10);
π-rung structures (Gb′=0: ring6, cube, unwound m=6) are immune; the
one-way ring (m=2, g_b=1.5e-5) is immune; the *smooth*-delay analytic
is stable — the pump needs all three of: **winding-displaced back gate
× its slope × the impulsive integrate-and-dump delivery** (sampled
feedback, not smooth delay). Growth is threshold-like in conductance:
+0.011 (gpl 0.3), +0.059 (0.5625), +0.050 (0.75); at dt=0.01 the peak
rate is +0.037 (a finite continuum part plus a dt-enhancement — both
present in the kernel, which runs the same dt).

**Reconciliation with the measured comp12 (alive at t=5000, −3e-4).**
Three measured-in-model facts close the apparent conflict: (i) the
instability *saturates* — x_std arrests at 0.06–0.10, gates degrade to
~0.6, and the ring continues as a slowly-radiating disordered lock,
not a fragmentation; (ii) foam jitter (5% link disorder, ring_m
closure) converts the coherent pump into an early quenched-roughness
shed followed by a **slow drain at −6e-4 /t.u. fractional — the
measured comp12 class**; (iii) real foam conductances (jittered
normals) sit near the pump threshold. The prediction stands for a
*clean* wound ring: **grow-and-saturate chiral waves with the ν₁/ν₂
lines above**, and for any wound ring an early transient at those
periods. comp12's 94%-mass-loss-yet-alive endpoint is the saturated
state of exactly this pump.

### pair (d=1.25, m=1)

−0.191 (load antisymmetric — the sign-theorem rate), −2.91 (phase),
plus slot modes ≈ −6.8..−7.3. All real.

---

## S4. Death taxonomy from the linearization

Ranked by measured share and by what the tangent analysis says. M(t)
forms are for the parallel regression agent; Λ symbols are per-t.u.
fractional rates.

**D1 — the skirt slide (dominant; kills everything passive).**
Common-mode neutral ray + environment bleed. While x > x_skirt:
  M(t) ≈ M₀·exp(−Λ₁t), Λ₁ = Λ_bleed(x)·[1 + rough share],
Λ_bleed ∝ Φ_vac·res_vac(x)·√(mob_floor·cap/Em) with
res_vac = Γ_m²/(Γ_m²+det_vac(x)²), det_vac = w₂qx/(1+qx). As x falls,
res_vac grows (Lorentzian) while mobility falls — the fractional rate
stays roughly flat, then **accelerates past the skirt** (res_vac → 1
inside |det_vac| < 2Γ_m, i.e. x < 0.0617) and ends in finite-time
dissolution: terminal form **M ∝ (t_d − t)^p, p ≈ 1–2**. Satellite
demo (`suite3`): x slides 0.128→0.09→0.04→0.01 with no plateau,
near-constant fractional rate through the skirt and a cliff at
x ≈ 0.003; A1's measured deceleration-then-acceleration
(−0.00118→−0.00063→−0.00146→death) is this curve at foam conductance.
Regression form: fit log M piecewise — early slope Λ₁; knee at
x(t) = x_skirt; terminal (t_d−t)^p.

**D2 — the chiral pump (wound even rings only; saturating, not fatal
alone).** Growth e^{+λt} of co-propagating load waves at λ up to
+0.05, λ threshold-like in Φ; saturation at x_std ≈ 0.06–0.10;
signature **before** any mass anomaly: growing per-voice load
modulation at period 9.8 / 6.0 t.u. traveling forward. Post-saturation:
elevated roughness drain (adds to Λ₁). M(t): unchanged early, then a
one-time shed of a few % at saturation, then D1 with larger Λ₁.
π-rung structures: absent. One-way rings: absent.

**D3 — seam/defect concentration and strangulation (seed-quality
deaths).** (a) Quenched link disorder δ_l puts static gate offsets
∝ δ_l and a standing roughness tax: early **linear** M(t) with slope
∝ Σ|δ_l| (the 5%-jitter ring sheds 75% by t=100 at full conductance —
the open-chain −0.09%/t.u. class is this channel at foam values);
concentration: the lock recursion dumps the closure residue on the
weakest link (measured seam gate ≈ 0 in fleet v2). (b) Strangulation:
channels die at d > 2r(Es); the ā=1.586 cube sits 0.04 from the edge
and its own displacement footprint (Es → ~0.71) closes its direct
edges — predicted fragmentation into weakly-coupled vertices, matching
"gates ~0.6, dies ≈1800". M(t): step drops at link deaths, then D1 per
fragment.

**D4 — parasitic siphons (foam cells the design didn't intend).**
A vacuum neighbor within reach starts at det_vac = −0.39 (res 0.06)
but *tunes itself in* as it accumulates crumbs (its x rises, pitch
falls toward the ring's): measured in the satellite demo, det_vac
−0.28 → −0.05 → 0 — an accelerating capture that ends indistinguishable
from the skirt (the ring meets the parasites at their common pitch).
Mid-life R(det) passes through its maximum (det ≈ Γ_r = 0.5: R = 1 —
35% of everything delivered radiates). Census signature: satellite
lumps appearing at sub-threshold masses near the ring; M(t): smooth
extra Λ₁ plus census n rising before death.

Ranking for a well-seeded π-rung ring: **D1 ≫ D4 ≈ D3 ≫ D2 (absent)**.
For a wound ring: D2 fires first (early, saturating), then D1.
For jittered/small shells: D3 first.

---

## S5. Active locking assessment (one page)

The linearization's message: **phases lock themselves** (all phase
modes damped at −0.15..−0.65; Huygens acquisition measured in e7) and
**internal detunes self-repair** (S2). What passive nets cannot do is
hold the *common-mode load* against the environment bleed: the neutral
ray leaks. "Active locking" therefore means **active load maintenance**
— phase control is pointless. Interventions consistent with the laws
(no new species, no suction, exact ledgers), ranked:

1. **The C1′ piston, B+ arm** (`design/C1_piston_design.md`; EXP-2 is
   already specified). Raises ambient Es at a ledgered boundary skin.
   It cannot move the rung loci or x_skirt (bound-energy pitch law —
   ambient-invariant, a free tripwire), but it acts on **every factor
   of the bleed the network cannot reach**: deepens the
   self-insulation contrast (displacement sealing, measured e3a
   0.00174→0.00084), extends the strangulation ceiling (d_cut 1.807 at
   pin 1.2 — rescues D3b shells), and slows skirt kinetics through the
   vacuum-side conductance. Cost: apparatus energy through R_s only.
   Addresses D1 kinetics, D3b, D4 (parasites insulate too).
   **Minimal lawful intervention; first choice.**
2. **Consonant-skin seed hygiene** (zero energy cost; addresses D3,
   D2). ring_m-exact closures (standing rule, closure < 0.05), edge
   uniformity, and — new from S3 — **avoid wound even rings at high
   conductance or seed them expecting the saturation shed**; prefer
   π-rung or one-way topologies where lifetime is the goal, winding
   where the species label is the goal.
3. **Quasi-static pressure schedules** (C2′/C3′ anneal): formation-side
   version of 1; ramps t1−t0 ≥ 1000 (piston note). Addresses D1 at
   birth (heavier, better-insulated survivors).
4. **Bath/record-media environment**: measured dead end as stabilizer
   (A1c closed box died at 1600 vs 1900 open — the shed bath does not
   recondense usefully). A *tuned* bath (ambient field at the ring's
   own pitch, inside the comb window) would re-feed via condensation,
   but is periodic re-seeding in disguise — see 5.
5. **Periodic re-seeding / feeding** (last resort; addresses D1
   directly and honestly). Ledgered Em injection at the seed pitch at
   the bleed rate (~5e-4·M/t.u.). Maximal cost; defensible only as an
   instrument to *measure* x* if a true equilibrium exists under the
   piston (i.e., to locate the plateau the passive system approaches
   but cannot hold).

**The lawful path to C1**: piston B+ on a π-rung ring (immune to D2),
ring_m-exact, ā safely inside the strangulation band — with the S3
prediction that its kick spectroscopy is overdamped and its kicked mass
loss is zero under the action atom (P-S4). If B+ produces the plateau,
C1 is reachable passively-in-the-piston; if it does not, the bleed
factors (res_vac·mob_floor) are law-level and the stable particle needs
S2-full (wave-borne translation) physics, not better engineering.

---

## Pre-registered numbers (for the simulation phase to score)

* **P-S1 (sign)**: pair/ring detune kicks decay heavy→light at
  λ_h = −0.19 ± 0.03 /t.u. (gpl=1 pair; scales ∝ Φ). No passive detune
  runaway below dE_asym ≈ 0.95·store.
* **P-S2 (ring6/cube lines)**: overdamped — no dense-sector breathing
  line; slowest relaxation −0.009 ± 0.002 /t.u. (Φ-robust); fast phase
  branch −0.15..−0.65 (Φ-scaled). Any observed ring is not (θ, Em).
* **P-S3 (wound-ring chirality)**: clean comp12-class rings show
  growing co-propagating load waves, ν₁ = 0.64 ± 0.05 and ν₂ = 1.05 ±
  0.08 rad/t.u. (one sense only), growth +0.01..+0.06 /t.u.
  (conductance-thresholded), saturating at x_std ≈ 0.06–0.10 with a
  few-% shed. π-rung and one-way controls show neither.
* **P-S4 (the atom freeze)**: kick-induced extra mass loss is exactly 0
  under quant_mode=2 whenever the kick's integrated rough demand per
  receiver stays below one atom ε_F (measured true for every kick up to
  a=0.1 on ring6); the continuous-limit control loses 0.3–0.9 per kick.
  Repeated kicks accumulate credit and eventually fire — kick-history
  memory.
* **P-S5 (death forms)**: D1 knee at x = 0.0617 with terminal
  (t_d−t)^{1–2}; D3 early-linear slope ∝ measured closure defect;
  cube-1.586 edge strangulation under its own footprint.

---

## Files

* `theory/ring_map.py` — the exact reduced map; modes: fixed / floquet
  / kick / pair / unlock / tax / all. All tables above are its output
  (suite logs in the session scratchpad; rerun takes ~10 min total).
* `theory/modal_2x2.py` — closed-form fixed point, per-mode 2×2,
  delayed-arrival characteristic roots, sign-theorem rate with the
  flight-delay factor.
