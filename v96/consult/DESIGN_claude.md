# v96 developmental design — the ground state: size selection and the quark ontology

**Session ID:** `0fc5d7e7-9f51-4efc-a43e-ce358e708e4c` (resume with
`claude --resume 0fc5d7e7-9f51-4efc-a43e-ce358e708e4c` when results land).
**Author:** claude (fable), developmental session, 2026-08-09.
**Deliverable for:** `v96/consult/BRIEF.md` (the user's two questions: why is
matter a specific size; are quarks parts or emergent attributes).
**Grounded in (read this session):** `v89/PRINCIPLE.md` (the one law),
`v93/RECKONING.md` (R1–R3: condensation is a transient; the raised stability
bar), `v94/consult/DESIGN_claude_face3.md` (the fiber-minimality theorem, the
one-load DNLS, the self-trap criterion), `v94/consult/NOCANCELLATION_claude.md`
(T1–T3, the rotation-primitive uniqueness theorem, DS-G3 confirmed),
`v94/FINAL_STATUS.md` (DS-G3 r=+0.684; M1 pair confines 17 vs legacy 417;
M6 park approaching; M3b knee Ee≈1.0 at door_g=0.5),
`v95/consult/DESIGN_claude.md` incl. PART VII (J1 GREEN R2=0.774, Em_min
×14–×670; J2 RED — the metric is an amplifier, not a generator; the digging
floor; the composition roads C1–C5), `v95/laws/laws_V3a.cfg` (w1=1.65,
w2=2.9, q_detune=1.2, cap=2.5, k_rad=0.05, p_rad=4, wf_on=1, quant_A0=1.15,
f_conv=0.25, field_J=1.8, dmin=1.0, r0=0.85).
**Status:** design + **resume update 1 (§W)** — S0/Q0/Q1 landed: ℓ confirmed
to the third decimal; α ≈ 1.1 (volume-limited — the pre-registered pivot
executed); Q0 negative (no preferred winding — the pre-registered negative);
Q1 negative as instrumented. §W re-grades the design against the results,
takes up the user's two observations (missing 1/r gravity; winding is a
number, not a geometry), and derives the candidate missing structural
pieces. Where §W supersedes a bar or a grading, §W governs.

---

## Executive summary

**Question 1 (size).** The substrate's answer is a liquid-drop balance, but
with one sign inverted from the nuclear textbook, and the inversion is both
derivable and already measured. The brief's suggested curve
`E(N) = −a_vol·N + a_surf·N^(2/3)` has **no interior minimum** — its unique
stationary point is a *maximum* (the critical nucleus). With those signs
alone the fate of matter is the measured dichotomy of the whole v90–v93
record: evaporate below the barrier, box-fill runaway above it (legacy door
M1: Em=417), and *no stable size anywhere* — which is precisely what the
programme kept measuring (v90 FORGE "no stable mass," v93 RECKONING "all
transients"). Real nuclei escape this only by importing Coulomb; the
substrate has no Coulomb and needs none, because it has something the static
drop does not: **matter is a flux machine** (the user's 2026-07-31 directive,
now load-bearing). A metabolizing object's upkeep is *volumetric* — every
cell at the park load radiates at the derived rate
ℓ = k_rad·cap·x^p_rad·(1+q·x) ≈ 0.032 per cell per t.u. — while its feeding
is *transport-limited*: beyond a derived crossover size the bath cannot
deliver field to the object faster than the diffusion-limited flux
I ∝ N^(1/3). Loss extensive, income sub-extensive ⇒ the per-cell metabolic
budget b(N) is **strictly decreasing in N**, and its unique zero is the
ground state:

```
N* ≈ (a/ℓ)^{3/2} · [1 − (3/2)·Λs₀/(ℓ·N*^{1/3})]     (stable: b′(N*) ≈ −2ℓ/3 < 0)
```

with a the diffusion-limited supply coefficient (bath level × field
diffusivity), Λ the per-surface-cell dense contact leak (the J-metric's
lever), s₀ = (36π)^{1/3} ≈ 4.84. The full size flow has exactly **three
fixed points — vacuum (stable), the critical nucleus N_c (unstable), the
ground state N* (stable)** — which is the complete answer to "everything
else is meta-stable": below N_c dissolves, between N_c and N* grows to N*,
above N* starves its interior back to N*. The proton is the size at which
its flame exactly consumes what its surroundings can feed it: **r* is a
flame radius**, and it is set by law constants plus the ambient bath level.
Two parameter-free scaling predictions fall out (N* ∝ Ē_e^{3/2},
N* ∝ k_rad^{−3/2}) and one retrodiction lands hard: the N_c formula, with
the measured J-metric contact starvation (w → 0.2), moves the nucleation
barrier by a factor ~5³ = 125 — putting the 6-cell ring *below* the barrier
without the metric and *above* it with the metric, which is exactly J1's
measured rescue (R2 0.774 vs 0.131, Em_min ×670). §I.9.
**[RESUME 1: S0 measured α ≈ 1.1 — the field-fed medium is volume-limited
at every accessible scale (field diffusion length ≈ 50 ≫ object ~5–10), so
the diffusion-limited branch and its N* sit beyond N ~ 10⁵ cells. The
theory's regime condition is unmet in the field-fed configuration, exactly
per the pre-registered failure semantics; ℓ was confirmed at 0.0322. The
structural consequence — and the road to an *accessible* N* — is §W.]**

**Question 2 (quarks).** The substrate's answer is **(C)** — sharper than
(B), incompatible with (A). A bound state is ONE ledgered object; what a
break produces is fixed by three things the whole owns: (i) the **spectral
decomposition** of the shattered state into self-bound lumps + radiation
(B's content — the "bell modes," here the DNLS soliton content); (ii)
**integer conservation of winding** — the one measured topologically
protected quantity (W=+2 exact, RING_DNLS A.1) must partition into integer
fragment windings, which is why fragments look like particles with quantized
charges rather than smears; and (iii) **the nucleation barrier N_c** — a
fragment below the critical size cannot persist (surface leak dominates its
budget), so it evaporates back to field or re-dresses onto a neighbor.
Point (iii) is **confinement, derived**: free sub-critical fragments are
never asymptotic states, exactly as free quarks are never seen — they exist
transiently, during the break, which is the user's intuition stated as a
theorem of the size theory of Part I. "Hadronization" is re-dressing: any
fragment above N_c relaxes to the *same* N* regardless of the break's
details, which is why smashing a proton yields protons and pions rather
than proton-shrapnel. Where does three-ness live? The substrate does not
derive the integer 3 and this design **refuses to impose it** (that would
be the quark road). What it derives is *which invariant the valence count
is*: the ground state's winding content — the number of minimal protected
winding carriers the break must emit. The substrate's own number is a
measurement (bar Q1), and per the programme's standing convention the
parallel to the real proton is structural, not quantitative. The
falsifiable (A)-vs-(C) discriminators — multiplicity vs probe energy,
pre-break identity continuity, winding bookkeeping across the break,
universality of re-dressing — are pre-registered as Q1a–Q1d, and one of
them (winding partition, Q1c) is runnable **today** on the measured J1
ring, before the ground state itself exists. §II.6.
**[RESUME 1: Q0 negative — no preferred winding (m = 0,1,2,3 all survive
equally at j_metric=1). Q1 negative as instrumented — uniform noise
degrades the ring continuously; no fragment spectrum (localized-probe
instrument debt stands). The registered three-ness hypothesis is
unsupported at ring scale; the diagnosis — winding is currently a geometry
the law never prices, not a number the ledger counts — is the user's
observation, developed in §W.4–§W.5.]**

The prohibition boundary, stated once and used throughout: an **attribute**
is a decomposition of the one field that the analyst (or a transient
fragment) realizes — delete it and the dynamics is unchanged. A **part** is
per-component state the *update law* reads and writes — separate amplitudes,
clocks, or ledgers inside the whole. Mode structure, winding partitions,
current channels are attributes and are licensed. Anything that gives a
"quark" its own ledger is a part and is the rejected road. DIS's three
sub-structures do not require parts: structure functions are correlators of
the one internal current field, and the substrate's bound state has exactly
the concentrated internal currents (|J| ≈ 1.4 on cycle links vs 10⁻² at
contacts, measured) for a hard probe to scatter off. §II.2, §II.5.

---

# PART I — Theory: the ground state and the size selection

## I.1 The record, read as a theorem: unfed matter has no ground state

Every stability failure in the v89–v93 record is one statement measured
from different sides. List the major ones with their feeding status:

| object | fed? | fate | where measured |
|---|---|---|---|
| v90 balance-curve hoards | no door | cap-wall step; "no stable mass" | v90 FORGE |
| v93 Face C condensates | amp_nat only, no door | transient, peak t≈300, melts by t=3000 | RECKONING R1 |
| v93 tight droplet | no door | melts (rms ×2.8 by t=1000) | RECKONING R2 |
| v91 metabolic chord | radiance-fed | lives t>5000 **at fixed size** | COMBINE (lock+tuner) |
| C6 ring in medium | no door, no metric | dead t≈8–16 | RING_DNLS A.2 |
| C6 ring + J-metric | accretion-fed, insulated | R2 0.774 @76, Em_min held | v95 J1 |
| door pair on film | door-fed | grows **at the field's pattern** | DS-G3 |

The pattern: **no object without an income has ever held, and no object
with an income has ever chosen its own size** — the fed objects held the
size they were seeded at (the ring's 6, the chord's 48) or grew without
bound until a cap intervened (legacy door, M1 Em=417). The programme has
measured an *intensive* fixed point — the load x̂* ≈ 0.62 of the radiance
campaign, now approached dynamically by the M6 chi_load park — but never an
*extensive* one. v96's question is exactly the missing half: the law
selects a density; what selects a size?

The answer cannot be a static energy minimum. The v89 principle forbids it
twice over: nothing persists to *have* a static energy (motion is regular
conversion), and the measured record confirms it (every static hoard is a
transient). A particle in this substrate is a **process** — the flux
machine — and a process is sized by the balance of its *rates*. So the
correct liquid-drop object is not E(N) but the **net rate functional**

```
Ṅ ∝ b(N)·N ,     b(N) = [I(N) − L(N)] / N      (the per-cell metabolic budget)
```

and the ground state is the stable zero of b. The rest of Part I derives
I(N) and L(N) from the law, finds the fixed points, and converts the result
back into the brief's E(N) language honestly (§I.5).

## I.2 The loss side: upkeep is exactly extensive

Every cell of a parked object radiates. The adopted law (laws_V3a,
RADIANCE): per cell, per unit time, in atoms,

```
L_cell = k_rad · cap · (Em/cap)^p_rad · (w2/w2e) ,    w2/w2e = 1 + q·x
```

At the park (x = x̂* ≈ 0.62, the measured radiance fixed point at
k_rad=0.05, p_rad=4; the M6 chi_load mechanism approaches the same point
dynamically):

```
ℓ ≡ L_cell(x̂*) = 0.05 · 2.5 · 0.62⁴ · (1 + 1.2·0.62) ≈ 0.032  per cell per t.u.
```

**[RESUME 1: MEASURED — 0.0322/cell/t.u. The radiance leg of the budget
theory is now measured law, to the precision quoted.]**

Two properties matter. **(i) Extensivity is exact:** radiance is per-cell
local; a droplet of N parked cells pays ℓ·N, no shielding, no surface
discount — the field the interior radiates into is its own dressing and
leaves (measured: dying matter's radiative far field r^−0.18, ASTRO; the
glow leaves the object even when nothing else does). **(ii) ℓ is intensive
and size-independent**, because x̂* is intensive — this is what makes the
budget argument clean: the loss slope in N is a constant of the law.

To the radiance upkeep add the **dense contact leak** — the surface term
proper: boundary cells drain Em into bath contacts by pass-U diffusion and
dephasing. Per boundary cell:

```
Λ = c_b · λ_link ,   λ_link ≈ 10⁻³–10⁻²  (RING_DNLS A.2, measured range),
                     c_b ≈ 3–4 bath contacts per boundary cell (z ≈ 7.3, S1 foam)
```

and the boundary population of a compact droplet is s₀·N^{2/3} with
s₀ = (36π)^{1/3} ≈ 4.84. Total loss:

```
L(N) = ℓ·N + Λ·s₀·N^{2/3}
```

**The J-metric enters here and only here** (§I.7): the measured
self-insulation (J1) drives contact weights to the 0.2 floor, so
Λ → Λ_J ≈ Λ/5 when j_metric=1 and the object carries a coherent internal
current. Radiance is untouched (the metric weights pass-U angles; radiance
is a click channel).

## I.3 The income side: feeding is transport-limited, and that is the whole point

The object eats through the door pair: coherent conversion at
J_D = 2g·Im(ψ_e*χ), settled in whole atoms ε = A₀·Ω₂/2π at beat windings.
Per fed cell the intake rate is bounded by the beat clock:

```
j = f_click · ε_m / T_beat ,
ε_m = A₀·w2/(2π(1+q·x̂*)) ≈ 1.15·1.663/2π ≈ 0.304 ,
T_beat = 2π(1+q·x̂*)/(w2−w1) ≈ 8.8 t.u.   ⇒   j ≤ 0.035 per cell per t.u.
```

(f_click ∈ (0,1] is the fraction of beats that fire — it depends on the
local field level through the M3b knee; it is a *measured* input, bar S0.
Note j_max ≈ 0.035 against ℓ ≈ 0.032: at law constants the fed margin
j − ℓ is thin. §I.8 discusses why that is interesting rather than
embarrassing.)

But j is only the *demand* side. The supply side is the bath's ability to
deliver field to the object, and it has three regimes, in increasing N:

- **(α=1) Reaction-limited (small N):** the bath everywhere around and
  inside the object holds its ambient level Ē_e; every cell eats at j.
  I(N) = j·N. Holds while the object's total consumption is small against
  local resupply.
- **(α=1/3) Diffusion-limited (large N):** the foam field is diffusive on
  scales beyond a few cells (the roughness/churn scattering the B7 opacity
  measured — mean free path of cells order unity). A consuming droplet of
  radius r in a diffusive bath draws at most the classic screened flux

  ```
  I_diff = 4π·D_e·ρ_e·r  =  a·N^{1/3} ,
  a = 4π(3/4π)^{1/3}·D_e·Ē_e/d̄² ≈ 7.8·D_e·Ē_e   (d̄=1 units)
  ```

  independent of how hungry the interior is — the transport theorem every
  droplet-growth problem obeys. D_e is the foam's field diffusivity
  (measurable in isolation: relax a field blob in a dark bath, fit r²(t);
  bar S0 does this).
- **(α=2/3) Surface/opacity-limited (intermediate, if opacity is strong):**
  if the door absorption length in *fed matter* is shorter than r, only a
  skin of depth δ eats: I = j·s₀·δ·N^{2/3}. Whether this regime exists
  between the other two depends on measured δ (XSEC measured absorption is
  clock-rate, i.e. per-cell per-dwell — δ follows from j and the field
  transit rate; expected δ ~ tens of cells, so this regime likely sits
  *above* the diffusion crossover and never governs — but S0 measures the
  exponent rather than trusting this).

The master income curve is the minimum of the applicable regimes:

```
I(N) = min( j·N ,  a·N^{1/3} )        [with the 2/3-branch inserted if measured]
```

**[RESUME 1: MEASURED — door_in = 0.075·N^{1.11}: the accessible range is
entirely reaction-limited (α ≈ 1, in fact mildly super-linear: per-cell
intake ∝ N^{0.11}, a weak cooperative/antenna effect). The field's
diffusion length ≈ 50 box units ⇒ the α = 1/3 branch begins only at
r ≫ 50, N ≫ 1.26×10⁵ — out of numerical reach. The field-fed N* exists in
the theory and not in the box. §W.1 for what this means structurally.]**

**The theorem the whole design stands on:** the radiance loss is *exactly*
linear in N; therefore a stable finite size exists **iff the income is
sub-linear at large N** — iff feeding is transport-limited. If the bath
could always feed every cell (α=1 forever), size would be neutral (any N
parks — the M6 situation extended, no selection) or runaway (if j > ℓ with
no supply ceiling — the measured legacy box-fill). Transport limitation is
not an assumption bolted on: it is forced by locality and the finite field
diffusivity of the foam, both already measured. **Size selection is what
locality does to appetite.**

## I.4 The three fixed points

Assemble the budget. Per-cell, b(N) = I(N)/N − ℓ − Λs₀·N^{−1/3}:

```
small N (reaction-limited):   b(N) = (j − ℓ) − Λs₀·N^{−1/3}
large N (diffusion-limited):  b(N) = a·N^{−2/3} − ℓ − Λs₀·N^{−1/3}
```

- **N = 0 (vacuum): stable.** The workfn presence rule (wf_on=1: no
  conversion without Em ≥ wf_floor) plus the atom grain (the M3b knee: no
  click below threshold drive) keep j = 0 on empty foam — the bath cannot
  self-seed (measured: W-L3 empty foam dark; the anti-ignition stack).
  Creation is nucleation on seeds only (measured: W-L3b; QUENCH — creation
  under a transient beam, i.e. under transiently raised local supply).
- **N_c (the critical nucleus): unstable.** In the reaction-limited branch
  the surface leak per cell diverges as N^{−1/3} at small N, so below

  ```
  N_c = ( Λ·s₀ / (j − ℓ) )³            (exists iff j > ℓ — the fed condition)
  ```

  the budget is negative: sub-critical objects **evaporate**. Above N_c
  (and below the crossover) the budget is positive: the object **grows** —
  boundary-adjacent bath cells accrete Em past wf_floor and enlist
  (nucleation-on-seed, the measured creation channel). This N_c is exactly
  the stationary point of the brief's E(N) = −a_vol·N + a_surf·N^{2/3} —
  the brief's curve is the *barrier*, not the ground state (§I.5).
- **N* (the ground state): stable.** In the diffusion-limited branch the
  income per cell falls as N^{−2/3} while the upkeep stays ℓ; the budget
  crosses zero once, from above:

  ```
  a·N^{1/3} = ℓ·N + Λs₀·N^{2/3}
  ⇒  N* = N₀·[1 − (3/2)·Λs₀/(ℓ·N₀^{1/3})] ,   N₀ ≡ (a/ℓ)^{3/2}
  ```

  Stability: b′(N*) = −(2/3)a·N^{−5/3} + (Λs₀/3)N^{−4/3} ⇒ the total-rate
  derivative F′(N*) = (a/3)N^{−2/3} − ℓ − (2Λs₀/3)N^{−1/3} ≈ ℓ/3 − ℓ =
  **−(2/3)ℓ < 0** (using aN^{−2/3} ≈ ℓ at the fixed point). Perturbations
  relax at rate ~(2/3)ℓ ≈ 0.02/t.u. — size recovery time ~50 t.u., well
  inside the T=3000 endurance window.

Above N*, the interior runs a deficit: cells fall below the park, then
below wf_floor, and de-list — the object **shrinks back to N***. The
ordering that makes the window exist:

```
N_c  <  N_x  <  N* ,     N_x = (a/j)^{3/2}   (the supply crossover)
```

N_x < N* follows from j > ℓ. N_c < N_x is a genuine condition on law
constants + bath level — **the window of creation**. Where it fails,
nothing nucleates and the foam is inert (the anti-ignition regime); where
it holds, matter both *can form* and *must stop growing*. The substrate's
existence-of-matter condition and its anti-ignition safety are the two
sides of one inequality, which is the kind of unification the programme's
one-law posture demands.

**Meta-stability of everything else (the user's sentence, derived):**
- N < N_c: dissolves to vacuum (transient — the entire v93 Face C record).
- N_c < N < N*: grows toward N* (accreting adolescents).
- N > N*: starves toward N* (the fate the legacy door never reached
  because its intake wasn't transport-honest; the pair's confinement, M1,
  is this physics seen from the intensive side).
- N = N* with a non-ground internal configuration (wrong winding, off-park
  load, excited mode): relaxes at fixed N by radiance — redward walk
  (measured: the death spectrogram; here re-read as de-excitation). These
  are the substrate's **resonances** — the Δ-analogue: same drop, excited
  configuration, finite lifetime set by radiance.

## I.5 The E(N) curve, honestly (what the brief asked for, with the sign audit)

Define the size potential V(N) by Ṅ = −dV/dN (gradient-flow reading of the
coarse-grained size dynamics):

```
reaction branch:    V(N) = −(j−ℓ)·N + (3/5)·Λs₀·N^{5/3} · (3/2)…      [barrier piece]
diffusion branch:   V(N) = (ℓ/2)·N² + (3/5)·Λs₀·N^{5/3} − (3/4)·a·N^{4/3}   [well piece]
```

(continuously joined at N_x; constants absorbed). Shape: V rises from
vacuum to a **maximum at N_c** (the brief's −a_vol·N + a_surf·N^{2/3}
stationary point — note d²E/dN² < 0 there: it was always a barrier top,
never a ground state), descends through the growth window, and reaches its
**unique minimum at N***, rising as (ℓ/2)N² beyond. The nuclear liquid
drop needs Coulomb to terminate growth; the metabolic drop terminates
growth with transport-limited feeding — the roles of "volume" and
"surface" in the brief's template are inverted on the income side because
this matter *runs*, it doesn't *sit*:

| term | nuclear drop | metabolic drop (SCP) |
|---|---|---|
| volume ∝ N | binding (income) | **upkeep (radiance loss ℓ)** |
| surface ∝ N^{2/3} | leak (cost) | leak (cost Λ) — same |
| growth terminator | Coulomb (imported) | **transport-limited income a·N^{1/3} (derived)** |
| ground state | none without Coulomb | unique N*, stable |

## I.6 Coefficients from the law (the emergence table)

| coefficient | formula | law inputs | value at law table |
|---|---|---|---|
| ℓ (upkeep/cell) | k_rad·cap·x̂*^p_rad·(1+q·x̂*) | k_rad, p_rad, cap, q_detune; x̂* itself from the radiance balance | ≈ 0.032 /t.u. **[RESUME 1: measured 0.0322 ✓]** |
| j (demand/cell) | f_click·A₀·w2 / (2π(1+q·x̂*)·T_beat), T_beat = 2π(1+q·x̂*)/(w2−w1) | A₀, w1, w2, q_detune; f_click(door_g, Ē_e) measured (M3b knee) | ≤ 0.035 /t.u. **[RESUME 1: measured 0.075·N^{0.11} per cell — above the beat-clock bound as spec'd; the composed config's effective grain differs from the single-cell estimate; recalibrate at S0 re-run]** |
| a (supply) | 7.8·D_e·Ē_e/d̄² | field_J, roughness (D_e); ambient bath Ē_e | **[RESUME 1: L_diff ≈ 50 ⇒ the a-branch is unreachable]** |
| Λ (surface leak/cell) | c_b·λ_link·⟨w_contact⟩ | amp_nat, q_detune contrast (λ_link); **j_metric (⟨w⟩: 1 → 0.2)** | (3–40)×10⁻³ |
| s₀ | (36π)^{1/3} | geometry | 4.84 |
| x̂* (park load) | radiance⁄intake fixed point | k_rad, p_rad, door_g, chi_load | 0.62 ± 0.02 (measured) |

Everything in N_c, N_x, N* is therefore law + one ambient (Ē_e). No new
constant is introduced anywhere in Part I. The ground state's totals:

```
E* = N*·x̂*·cap ≈ 1.55·N*        r* = d̄·(3N*/4π)^{1/3}
```

**r* is a flame radius:** the size at which diffusion-limited feeding
equals volumetric burning, 4πD_e·ρ_e·r* = ℓ·(4π/3)r*³·n_c, i.e.
r* = √(3D_e·Ē_e/ℓ) in cell units — the proton-size analogue is a derived
length scale of the law plus its environment. (A candle flame, a raindrop
in vapor, and this object are sized by the same theorem; the substrate
earns the analogy rather than importing it.)

## I.7 The J-metric's role: it moves the barrier a hundredfold and the ground state a little

The metric modulates exactly one coefficient, Λ (contact starvation,
w → 0.2 on incoherent contacts — measured, J1). Effects, derived:

- **N_c ∝ Λ³** ⇒ metric on: N_c falls by ≈ 5³ = **125×**. Nucleation
  becomes radically easier for coherent objects — and *only* for coherent
  objects, since the metric starves only phase-incoherent contacts (the
  synchronous-detection selectivity, v95 §I.4). Incoherent condensates get
  no discount: the metric is a **coherence filter on existence** at small
  N, the same role it measurably plays for retention.
- **N*: small increase.** ΔN*/N₀ = (3/2)·(Λ_off − Λ_on)·s₀/(ℓ·N₀^{1/3})
  = 1.2·Λ_off·s₀/(ℓ·N₀^{1/3}) ≈ **4–40%** across the measured λ_link
  range (at N₀ ≈ 30). Sign: positive (less surface tax → slightly bigger
  drop). Bar S3 tests the sign and ordering, with magnitude read off the
  measured Λ (JMAP quartiles), not assumed.

## I.8 Two scaling laws and one thin margin (the falsifiable spine)

From N₀ = (a/ℓ)^{3/2} with a ∝ Ē_e and ℓ ∝ k_rad:

```
S1b-i :   N* ∝ Ē_e^{3/2}          (bath-level tuning; exponent 1.5)
S1b-ii:   N* ∝ k_rad^{−3/2}       (law-constant tuning; exponent −1.5)
```

Both are parameter-free exponents (the prefactors carry the uncertainty;
the exponents don't — caveat: x̂* shifts weakly with k_rad, contributing a
correction to (ii) through ℓ(x̂*); the pre-registered tolerance ±0.3 on the
fitted exponent absorbs it and S0's measured x̂*(k_rad) removes it).
Practically, (i) is also the campaign's *reachability dial*: choose Ē_e to
put N* ≈ 20–50 so the ground state fits an affordable box (§III.4).
**[RESUME 1: the dial does not reach — L_diff ≈ 50 puts the crossover
beyond any affordable box regardless of Ē_e in the field-fed medium. The
scaling spine transfers intact to the space-fed configuration, §W.3.]**

**The thin margin, flagged as physics:** at law constants j_max − ℓ ≈
0.003 — the fed margin is ~10% of either rate. N_c ∝ (j−ℓ)^{−3} is
violently sensitive to it. This is not a defect of the derivation; it is a
candidate *explanation* of the record: a substrate near j ≈ ℓ has a huge
nucleation barrier almost everywhere (matter is rare and hard to ignite —
anti-ignition for free) yet supports arbitrarily long-lived matter where
the barrier was once crossed. The QUENCH measurement — creation only under
a transient beam, then self-sustainment — is this sensitivity operating: a
beam transiently raises local supply (a↑, effectively j↑), collapses N_c,
nucleates, and withdraws. **[Load-bearing assumption A-I1: f_click at the
park is O(1) rather than ≪1; if S0 measures f_click ≪ 1 the fed condition
j > ℓ fails at this bath level and the whole fixed-point structure moves
to higher Ē_e — the theory survives, the numbers move. S0 decides.]**

## I.9 Retrodictions (the theory against the standing record)

1. **v90/v93 "no stable mass"/"all transients":** every one of those
   objects had j = 0 (no door) — with no income the budget is b = −ℓ −
   Λs₀N^{−1/3} < 0 for *all* N: no fixed point exists, everything decays
   on the radiance/leak clock. The theory does not merely permit the
   negative record; it **requires** it.
2. **M1 legacy runaway (Em=417) vs pair confinement (17):** the legacy
   magnitude door's intake was not transport-honest (it consumed field
   magnitude where it stood, DSG finding 2), realizing the α=1, j>ℓ,
   no-ceiling branch: box-fill, as measured. The pair's coherent intake +
   chi_load throttle restores the intensive brake.
3. **J1's rescue, quantitatively:** with Λ at the upper measured range,
   N_c(metric off) ~ O(10²–10³); with the metric's 5× cut, N_c(on) ~
   O(1–10). The 6-voice ring sits *between* — sub-critical bare,
   super-critical insulated. J1's measured Em_min ×14–×670 (no slip site
   opens) is the barrier crossing seen cell-by-cell. **[This is a
   retrodiction with error bars, not a fit — the honest range of λ_link
   spans it; S2's bisection turns it into a measurement.]**
4. **The v91 chord at nv=48 (lock+tuner) holding at its seeded size:** an
   asserted-identity apparatus substitutes for the door income
   (maintenance pays upkeep); with income pinned per-cell by apparatus,
   b(N) ≡ 0 at any N — size-neutral, as observed (it neither grew nor
   chose a size). Assertion hides exactly the physics v96 measures.

**[RESUME 1: retrodiction 5, new — S0 itself: α ≈ 1.1 with j·N^{0.11} vs
ℓ·N is a knife-edge between box-fill (j > ℓ) and evaporation (j < ℓ) with
no interior zero — "those signs alone reproduce the measured record
(evaporate-or-box-fill, no stable mass)," now measured as the operative
regime of the field-fed medium at accessible scales.]**

## I.10 Proton and neutron: one size, two configurations

N* depends only on {a, ℓ, Λ} — none of which reads the internal phase
configuration. Therefore **every ground-state species has the same size**
to within the surface correction: the substrate's proton/neutron
size-degeneracy is structural. What distinguishes species at fixed N* is
the internal configuration: the winding content W (and, at the fiber
level, the beat-phase relationship of the metabolic cycle — seed_mw's
role). Configuration energy differs by the winding kinetic term
(ΣJ across the cycle — order 2·Em·sin²(πm/n) per link), a small split on
top of the common bulk E* — a mass splitting hierarchically below the mass
itself, as with the real nucleon pair (structural parallel only, per the
standing convention). Species are told apart spectroscopically (B4:
species by line ratios — the standing instrument), and the de-excitation /
inter-configuration relaxation channel is radiance (the redward walk).
Whether one configuration is strictly stable and another weakly unstable
(a β-analogue) is left as an outlook question for the v96 record — it
requires the configuration-resolved spectroscopy of the S-series winners,
not new theory. **[RESUME 1: Q0 measured NO configuration preference at
ring scale — the winding kinetic term either vanishes against the noise
floor or is not, in fact, priced by the law. §W.4.]**

---

# PART II — Theory: the quark ontology

## II.1 The question, formalized: the part/attribute litmus

Let S(t) be the substrate's full state (the one fiber field + geometry +
ledgers). Consider any decomposition D that exhibits "three things" inside
a bound state — mode occupations, winding partition, current channels,
spectral components.

- D is an **attribute** iff the update law never reads or writes D's
  components separately: erase the decomposition (forget it) and S(t+dt)
  is unchanged. Attributes live in the analyst's frame — or transiently in
  the physics when a break *realizes* the decomposition as spatially
  separated fragments.
- D is a **part** iff the update law carries per-component state:
  separate amplitudes, separate clocks, separate ledgers whose identity
  survives the whole. Erase a part's ledger and the dynamics changes.

The pre-v89 quark model was parts: per-component fields Φ_a with per-
component charges Q_a and clocks ω_a on an immortal lattice — every
component separately ledgered, forever, everywhere, including in vacuum.
That is what CLAUDE.md prohibits, and the prohibition is not aesthetic:
per-source identity inside a superposition is exactly what the
no-cancellation analysis (v94 §1.2) identified as "ledgers that outlive
the cell" — it re-imports a background one level down.

## II.2 The boundary, drawn precisely (what is licensed, what is not)

**Licensed (attributes of the whole):**
- The mode fiber (ψ_e, χ, ψ_m) per cell — three *modes* of one ledgered
  energy, proven minimal (the fiber theorem), never three particles: no
  fiber component has its own conserved charge; energy moves freely among
  them under one law.
- Winding numbers of the one phase field on the one support — topological
  attributes; W=+2 was measured on ONE object.
- Spectral/soliton decompositions of a bound state — computed from S, not
  stored in S.
- Transient fragments — spatially separated lumps *after* a break are
  simply the one field in a multi-lump configuration; each fragment is a
  whole in its own right, ledgered by the same one ledger.

**Prohibited (parts):**
- Per-"quark" amplitudes co-resident on the same cells with separate
  update terms.
- Per-"quark" conserved charges enforced by the law (a flavor ledger).
- Per-"quark" clocks (separate pitch families per component).
- Any identity that the *update* preserves through the bound state's life
  (the gid apparatus is a *meter* and stays one; the moment a gid enters
  the dynamics it becomes a part).

The design below never crosses the line: every "three-ness" candidate is
an attribute, and the experiments *measure decompositions*, they do not
*install components*.

## II.3 The substrate's answer: (C) — the break-spectrum of one object, with topological bookkeeping and a confinement barrier

**(A) is rejected** twice: by the prohibition (II.2) and — more
interestingly — by *redundancy*: Part I's theory needs no parts to make a
bound state, and II.5 shows DIS-style phenomenology needs no parts either.
A hypothesis that is both forbidden and unnecessary is not a tension; it
is closed.

**(B) as stated is too weak.** Bell modes are the right intuition with the
wrong physics: linear normal modes of a struck bell are neither conserved,
nor quantized, nor particle-like — they disperse. If quarks were only
"what the vibration decomposes into," fragments would be smears of
radiation with arbitrary content. The measured substrate adds exactly the
three ingredients that turn decomposition-modes into particle-like
emergents:

1. **Nonlinearity ⇒ discrete lump content.** The dense channel at
   amp_nat>0 with q_detune is a DNLS (face-3 §II.1); a violently perturbed
   DNLS state resolves into **self-bound lumps + radiation** — a discrete
   soliton census plus a continuum, the (approximate) inverse-scattering
   picture. The number and sizes of the lumps are functionals of the
   *whole pre-break state* — determined by it, not resident in it. This is
   (B) upgraded from "modes" to "lumps."
2. **Topology ⇒ integer charges.** The one measured protected quantity is
   winding. Total winding is conserved through the break (it is a phase-
   circulation integral of the one field; T1–T3 protect the ledger, the
   topology protects the integer); fragments carry **integer** windings
   summing to the whole's W. This is why break products have quantized
   "charges" without any charge being imported: the quantization is
   topological bookkeeping of one phase field. **[Registered caveat: W
   conservation through a *shattering* event is exact only while phase
   stays defined on the relevant loops; slip events (measured at support
   holes, A.2/A.3) can change W by ±1 with a corresponding radiated
   burst. Q1c measures the bookkeeping including slips, not assuming
   exactness.]**
3. **The nucleation barrier ⇒ confinement.** Part I: a fragment below N_c
   has strictly negative budget — it **cannot persist**, no matter how it
   was made. It evaporates back to field (most debris) or re-dresses onto
   a super-critical neighbor. Sub-critical fragments are therefore
   *transient by theorem*: visible during the break, never asymptotic.
   **This is the substrate's confinement**, and it is not a new force — it
   is the same surface-vs-income bookkeeping that sizes the proton. The
   user's sentence — "we see them emerge transiently when we break them
   apart" — is this theorem, verbatim.

And the closing act, **hadronization = re-dressing:** any fragment above
N_c flows to the *same* ground state N* (Part I's attractor is universal).
Break a ground-state object hard and the asymptotic products are — ground
state objects (plus radiation), regardless of the break's violence. The
substrate predicts jets of re-dressed matter, never bare shards. That the
observable products of proton collisions are always hadrons is, in this
frame, the visible face of the size attractor.

## II.4 Where the three-ness lives — ranked, with the refusal stated

The brief's candidates, graded against the substrate:

1. **Winding content of the ground state (STRONGEST — adopted as the
   registered hypothesis).** The valence count = the number of minimal
   protected winding carriers the break must emit = the integer partition
   of the whole's W into stable fragment windings (|W_frag| = 1 carriers
   are the smallest protected objects; a W = m ground state yields m
   valence fragments + unwound debris). This is derivable *given* the
   ground state's W — but the substrate does not yet derive W(ground
   state): whether the size attractor prefers a wound configuration (the
   winding's internal current feeds the J-metric insulation, cutting Λ —
   a real energetic preference for W ≠ 0, measured in J1's Λ/5) and
   *which* W it prefers is a measurement (bar Q0). Note the honest
   direction of inference: **if the ground state measures W = 3, the
   three-ness is derived; the design does not seed W = 3 to get it.**
   **[RESUME 1: Q0 NEGATIVE — m = 0,1,2,3 all survive equally; no
   preference. The registered hypothesis is unsupported at ring scale.
   Diagnosis: the law does not price winding (§W.4); either the pricing
   is missing structure (W-B) or the preference emerges only in bulk
   objects (W-C). The refusal below stands and has now been exercised.]**
2. **The three phase relationships of the metabolic cycle (e–χ, χ–m, and
   the composed e–m beat).** Real structure — the fiber's two edges plus
   their composition — but these are *mode-direction* attributes: they do
   not partition into spatially separated fragments (every fragment
   carries the whole fiber). They are the substrate's analogue of
   *internal quantum numbers of the cycle*, not of fragment count.
   Retained as spectroscopy (the beat spectrum is measurable, B4-style),
   rejected as the origin of fragment three-ness.
3. **The three sectors (field/matter/space).** These are modes of the one
   energy — a break does distribute energy into all three (fragments,
   glow, space books), and the three-way split is measurable and
   conserved. But "three sectors" is a fixed architectural fact, not a
   variable attribute of a particular bound state; it cannot explain why
   a *baryon-analogue* would show three and something else would not.
   Rejected for three-ness; kept as the break's energy bookkeeping.
4. **m = 0, ±1 winding modes / spin.** Conflates the spin-multiplicity
   question with the constituent question; parked for the spin lane.

**The refusal:** if Q0/Q1 measure a ground state with W = 1 whose break
yields one protected fragment + debris, then the substrate's baryon
analogue is a one-valence object and the structural parallel to the
3-quark proton *fails at this rung* — and that negative would be recorded
exactly as the census/Q-ball negatives were. The programme's convention
("parallels are structural, not quantitative") governs; nothing in this
design bends the substrate toward 3.

## II.5 Deep inelastic scattering without parts: partons as current concentrations

What does a hard probe actually measure? Correlators of the target's
internal currents at short wavelength. The substrate's bound state owns
concentrated internal currents as a *measured fact*: the wound core
carries |J| ≈ 1.4 per cycle link while contacts carry 10⁻²–10⁻³ — a
two-to-three-decade contrast — and the J-metric *consolidates* this
contrast into the geometry (cycle links at the ×5 ceiling, contacts at
the 0.2 floor: the current structure is written into the roster itself).
A short-wavelength probe therefore scatters not off a smooth drop but off
the **discrete, stiff, current-carrying structure inside it** — pointlike
scattering centers at high momentum transfer, from ONE object with zero
parts. The centers' multiplicity and form factors are set by the wound
core's structure (how many links carry the circulation; their spatial
arrangement), i.e. by the same winding content that II.4 nominates for
valence. Prediction (qualitative, pre-registered): probe-energy sweeps in
Q1 should show a transition from collective response (long-wavelength:
the drop recoils/deforms as one) to center-resolved response
(short-wavelength: fragments emerge along the current concentrations) —
the substrate's scaling-onset analogue. Structure functions, in this
frame, are Fourier transforms of one field's current-current correlator —
no per-parton ledger anywhere in the story, exactly as in II.2.

## II.6 The falsifiable discriminators, (A) vs (C) — pre-registered

- **D1 — multiplicity vs probe energy.** (A): fragment count fixed by the
  part count (3), independent of E_probe once above threshold. (C):
  protected-fragment count fixed by W (constant in E_probe), while
  *total* transient multiplicity grows stepwise with E_probe (each ~N_c·
  x̂*·cap of deposited energy can fund one more transient lump over the
  barrier). Substrate predicts (C): a constant valence floor + an
  energy-dependent sea. Both counts measured separately (winding-carrying
  vs unwound fragments).
- **D2 — pre-break identity continuity.** (A): the three parts exist
  before the break ⇒ tri-partite identity structure in the *quiet* bound
  state, with continuity (the "same third" persists). (C): internal mode
  occupations exchange freely; no sub-object identity survives internal
  churn. Instrument exists: the episode-gid meter (v91 IDENTITY lane —
  built, inert, measured on exactly this question one level down; it
  found "identity preserves, does not create"). Run it *inside* the bound
  state: partition the object's cells three ways at t₀ by any rule and
  measure whether the partition is dynamically distinguished at t₀+Δ.
  (C) predicts no distinguished partition exists before the break.
- **D3 — winding bookkeeping across the break.** ΣW_frag + ΔW_slips =
  W_initial, with slips localized at measured support holes and paid in
  radiated bursts. Sharp version, runnable NOW on the J1 ring (no ground
  state needed): break W = +2 vs W = 0 vs W = +1 objects at matched N and
  E_probe; (C) predicts the protected-fragment census tracks W and the
  W = 0 break yields debris only; (A) has no mechanism tying fragment
  count to a seeded winding. **[RESUME 1: the noise-sweep arm ran — the
  ring degrades continuously at all noise levels 1.0–8.0 (Em_mean 0.59
  even at 8×: the insulation is remarkably robust), no fragment spectrum.
  Negative as instrumented; uniform noise is a bath heater, not a break —
  the localized-probe instrument (III.1) was never built and remains the
  registered debt. §W.5 re-scopes.]**
- **D4 — universality of re-dressing.** (C): fragment end-states converge
  to the same N* (and the same species set) whatever the probe energy or
  angle — hadronization is attractor dynamics. (A): end-states inherit
  the pre-existing parts' properties; varying which "part" was struck
  varies the products. Measured as the N-histogram of surviving fragments
  at T_long across probe conditions.

**Refutation conditions for (C), stated against ourselves:** a measured
distinguished tri-partition in the quiet state (D2 positive) or fragment
counts locked at a probe-independent integer with no winding correlate
(D1 anti-C + D3 null) would falsify the attribute reading and reopen the
ontology — that finding would gate any further v96 claim.

---

# PART III — Numerical simulation method

## III.1 What is and is not built (the constraint audit)

**No kernel physics change.** The entire campaign runs on the v95 kernel's
existing knob surface (amp_nat, amp_logate, hop_order, q_detune, door_g,
door_avg, chi_load, door_out_enable, j_metric, j_tau, k_rad, p_rad, wf_on,
noise_amp, bath, exp/seed family). Additions, all apparatus or meters:

1. **Seed generators (fab, sanctioned apparatus per the standing
   convention):**
   - `exp=drop`: compact droplet seed of N cells at load x_seed (N and
     x_seed as knobs; support = BFS ball on the foam from a center cell) —
     needed because S1's seed sizes {3…48} must be placed without imposing
     shape beyond compactness. Optional `drop_mw` for winding
     configuration on the drop's outer cycle (Q0/Q1 arms; default 0).
   - `probe_*`: a *field-sector* Gaussian packet (position, width, kx,
     amplitude) superposable on any experiment — the smash apparatus. The
     L1-A tilted-blob precedent is dense-sector; the probe must be
     field-sector so it travels under pass F and hits the object through
     the door, the honest DIS geometry. (Two-stage alternative if
     superposition at t=0 proves insufficient: converge, dump, re-seed
     from dump + packet — the fcs dump/read path exists in the report
     tooling; decide at implementation, no physics either way.)
2. **Diag-only meters (both kernels, print-only, no state feedback):**
   - Shell books: per radial shell from the object COM — door intake
     (clicks·ε), radiance out, dense flux through the shell boundary
     (the XSEC sector-meter/convtag precedent). This is the I(N), L(N)
     spectrometer; without it S0/S1 are number-less.
   - Object census per diag row: N(t) (largest connected component with
     Em ≥ wf_floor, union-find — R3 precedent), r_rms, Em_tot, x̄, W of
     the dominant internal cycle (analyze_hoard_phase.py logic promoted
     to a standing report script, not a kernel meter).
   - Fragment census (report script, post-hoc on fcs dumps): per-fragment
     N, Em, W, COM, fate-tracking across dumps.
3. **Battery:** all additions inert at default ⇒ battery 87 GREEN ×2
   kernels before any campaign run (the ratchet); conserve-diff
   byte-identity at defaults.

## III.2 The common configuration (the fed medium)

All strength runs: laws_V3a + `amp_nat=2 amp_logate=1 hop_order=1`
(Strang, adopted) + the composed v94/v95 machinery
`door_g=0.5 door_avg=1 chi_load=1 door_out_enable=1 j_metric=1 j_tau=4` +
lit bath `bath=1 noise_amp=<set per arm>` + `q_detune=1.2` (law), Go
kernel at strength (X.6 language; C side tolerance-graded, pre-existing
debt), 3 seeds minimum per claim (the M4 lesson). `dt=0.02`.

**The closed-box protocol (a derived requirement, new to v96):** in a
closed box the pair (N*, Ē_e) co-adjusts — the object's equilibrium size
depends on total box energy, since eating depletes the bath that feeds
back. Therefore: (i) **all S1 arms run at matched E_tot** (bath amplitude
adjusted to compensate seed-size differences — otherwise "convergence to
common N*" is confounded by different ambient levels); (ii) box large
enough that the bath reservoir ≫ E* (L=16 minimum, L=24 for the N*
confirmation arm; bath energy ~Ē_e·L³ vs E* ~ 1.55·N*); (iii) the S1b-i
scaling arm *uses* the Ē_e dependence rather than fighting it.

## III.3 The experiment ladder (cheap deciders first, v95 discipline)

- **S0 — the regime meter (run FIRST; decides everything downstream
  cheaply).** Three short measurements:
  (a) D_e: relax a field blob in dark bath, fit r²(t) — the diffusivity.
  (b) f_click and j: a single fed cell / tiny drop (N=3) in lit bath,
  shell books → intake per cell vs Ē_e (the M3b knee, re-measured in the
  composed config).
  (c) α(N): drops N ∈ {3, 6, 12, 24, 48}, short T (~200), shell books →
  I(N) log-log slope; and L(N) decomposition (radiance vs contact leak →
  ℓ, Λ directly; JMAP quartiles → Λ's j_metric modulation).
  **S0 outputs the predicted N_c, N_x, N* with error bars before any long
  run.** If α is measured ≈ 1 with no downturn across the accessible
  range at any affordable Ē_e — the reaction-limited ceiling — then size
  selection is out of numerical reach at this law table; record, and the
  campaign pivots to the S1b scalings at the largest reachable contrast
  (the theory is then tested through its exponents, not its fixed point).
  **[RESUME 1: EXECUTED — α ≈ 1.11, L_diff ≈ 50, ℓ = 0.0322. The
  reaction-limited ceiling is the measured reality; the pivot clause
  fired. §W.]**
- **S1 — convergence to the ground state (the brief's bar).** Seeds
  N ∈ {3, 6, 12, 24, 48} at matched E_tot, T=3000 (the RECKONING
  endurance bar governs — no transient may pass). Measure N(t) per arm.
  **GREEN: all arms that survive converge to a common N* (±20%, seed
  spread), where N* agrees with S0's prediction within its error bars;
  arms seeded below the measured N_c die — and that death is *part of the
  pass* (it is S2's barrier seen from below).**
  **[RESUME 1: BLOCKED by S0's regime result in the field-fed medium —
  no N* exists at accessible scale to converge to. Re-scoped to the
  space-fed configuration if W-A is authorized (§W.3), else stands as
  the recorded reachability wall.]**
- **S1b — the scaling spine.** (i) Ē_e sweep ×{0.5, 1, 2}: fitted
  N* exponent 1.5 ± 0.3. (ii) k_rad sweep ×{0.5, 1, 2} at fixed Ē_e:
  exponent −1.5 ± 0.3 (with S0's x̂*(k_rad) correction applied). These are
  the parameter-free claims — the campaign's falsifiable core even if the
  absolute N* is imprecise. **[RESUME 1: blocked with S1; transfers to
  W-A's space-fed form verbatim.]**
- **S2 — the barrier (meta-stability below).** Bisection on seed N at
  fixed bath: locate N_c between the largest dying and smallest growing
  seed. GREEN: three-regime structure observed (die / grow-to-N* /
  shrink-to-N*), N_c within a factor ~3 of the S0 formula (the (j−ℓ)^{−3}
  sensitivity forbids a tighter bar honestly). **[RESUME 1: the N_c half
  (die vs grow) remains runnable in the field-fed medium — the barrier
  exists even where N* does not; still registered.]**
- **S3 — the J-metric tuning.** Full S2 bisection + one super-critical
  arm at j_metric ∈ {0, 1}: **pre-registered ordering N_c(1) ≪ N_c(0)
  (target: ≥ 10× at the measured Λ ratio; the formula says ~125× at a
  clean 5× Λ cut) and N*(1) ≥ N*(0) (small, sign only).** The N_c arm is
  the strong test — it retrodicts J1 (§I.9.3) and it is where the
  metric's 5× becomes a 125×: the most falsifiable single number in the
  design. **[RESUME 1: N_c arm stands; N* arm blocked with S1.]**
- **Q0 — the ground state's own configuration.** Take S1's converged N*
  object: measure its W (does the attractor spontaneously carry
  winding?), its internal current map (JMAP — how many links carry the
  circulation), and its beat spectrum (spec.awk). This is the "what IS
  the proton analogue" measurement, and II.4's three-ness hypothesis
  lives or dies by it. **[RESUME 1: EXECUTED on the available objects
  (rings m=0–3): no preference — the hypothesis is unsupported at ring
  scale. §W.4.]**
- **Q1 — the smash.** Probe packets at E_probe ∈ {0.5, 1, 2, 4, 8}·E*
  into (a) the S1 ground state, (b) the J1 wound ring (runnable
  immediately, before S1 lands — D3 does not need the ground state),
  (c) **the QUENCH control: identical probe into bare bath** — the beam
  itself creates matter (measured, v91 QUENCH); the fragment census must
  be background-subtracted against this arm or every multiplicity claim
  is confounded. Measure per III.1's fragment census: multiplicity
  (protected vs debris), winding bookkeeping (ΣW + slips), fates to
  T_long = 2000 (evaporate below N_c / re-dress to N* above), and the
  D1–D4 discriminators. **[RESUME 1: run as a uniform noise sweep, not a
  localized probe — negative as far as it goes (continuous degradation,
  no spectrum). The probe-packet apparatus remains the registered
  instrument debt; re-scoped in §W.5 to the bulk-vortex target.]**

## III.4 Compute honesty

L=16 → ~4k cells; L=24 → ~14k; T=3000 at dt=0.02 = 150k steps. The
S-series is dozens of such runs; at v93/v95 measured throughput this is
days of Go-kernel CPU across the ladder — front-load S0 (short runs decide
the interesting region), run S1's five arms concurrently, and gate the
L=24 confirmation arm on S0's predicted N* actually needing it. If S0
puts N* > box at every affordable Ē_e, do not brute-force: report the
scaling-spine results and the reachability wall — that too is a measured
answer to Question 1 (the substrate selects a size; the box couldn't hold
it). **[RESUME 1: this is what happened — N* > box by three orders. The
wall is the finding.]**

---

# PART IV — Acceptance bars (pre-registered; the ratchet governs)

**Invariant surface first, always:** battery 87 GREEN ×2 kernels at all
new knobs default before and after the apparatus lands; conserve-diff
byte-identical at defaults; drift ≤ 1e-13 on every strength run;
anti-ignition control (bare lit bath, no seed, full campaign config,
T=3000: no nucleation — N(t) census empty beyond noise clusters).

| bar | statement | GREEN condition | RESUME 1 |
|---|---|---|---|
| **S0** | the regime and the coefficients | α, D_e, f_click, ℓ, Λ measured; predicted (N_c, N*) with error bars published *before* S1 runs | **EXECUTED**: α≈1.11, ℓ=0.0322 ✓, L_diff≈50 ⇒ N* unreachable field-fed |
| **S1** | size selection | all surviving seed arms converge to common N* (±20%) at matched E_tot, T=3000, 3 seeds; N* within S0's error bars | blocked (regime); transfers to W-A |
| **S1b** | the scalings | N* ∝ Ē_e^{1.5±0.3} and ∝ k_rad^{−1.5±0.3} | blocked; transfers to W-A |
| **S2** | meta-stability | three-regime structure (die/grow/shrink) observed; N_c within ×3 of formula | N_c half still runnable |
| **S3** | J-tuning | N_c(j=1) ≪ N_c(j=0) (≥10×); N*(j=1) ≥ N*(j=0); Λ modulation confirmed in JMAP | N_c arm stands |
| **Q0** | the ground state's structure | W, current map, beat spectrum of the converged object measured and on file | **EXECUTED**: no W preference (negative, pre-registered) |
| **Q1a** (=D1) | multiplicity law | protected count constant in E_probe; total multiplicity rising stepwise; both background-subtracted vs the QUENCH control | open (instrument debt) |
| **Q1b** (=D2) | no pre-existing parts | no dynamically distinguished tri-partition in the quiet object (gid-meter) | open |
| **Q1c** (=D3) | winding bookkeeping | ΣW_frag + slips = W₀ across breaks; protected-fragment census tracks seeded W ∈ {0,1,2} (ring arm — runnable now) | noise-sweep arm negative-as-instrumented |
| **Q1d** (=D4) | re-dressing universality | surviving-fragment N-histogram peaks at N* independent of probe condition; sub-critical fragments 100% non-asymptotic (the confinement bar) | open |

**Failure semantics, pre-registered:** S0 α≈1 everywhere → record, pivot
to scaling spine (§III.3). S1 no convergence at measured α<1 → the budget
theory is wrong in an interesting way — publish the measured b(N) curve
and consult. Q1b positive (a distinguished tri-partition) or Q1a/Q1c
jointly anti-(C) → the ontology section is falsified as registered
(§II.6) and the finding gates v97. No claim of "ground state," "proton
analogue," or "confinement" before the corresponding bar is green at
T=3000 — RECKONING's transient warning stands over every one of them.

---

# PART V — What NOT to do, and the flags

**Reaffirmed prohibitions (from the brief, made operational):**
- No imposed size: seeds span sizes; every size claim is about the
  *attractor*, and the S1 bar fails if arms merely keep their seeds.
- No imposed three-ness: no three-component anything; W is seeded only in
  arms whose *purpose* is to vary it (D3), and the ground state's own W
  is Q0's measurement, not an input.
- No QCD import: no color, no flavor ledger, no SU(3); the only group
  theory in this design is U(2) per edge (already the substrate's) and
  integer winding (already measured).
- No background, no separate ledgers: all additions are seeds and
  print-only meters; the gid instrument stays a meter (II.2).
- The door stays special; scp_sim/sfa.h untouched (not used in v96).

**Load-bearing assumptions, collected:**
- **A-I1** f_click = O(1) at park (else the fed window moves; S0
  decides). §I.8. **[RESUME 1: superseded by the measured per-cell intake
  0.075·N^{0.11} — recalibrate j from the S0 books at next resume.]**
- **A-I2** The foam field is diffusive at bath scales (D_e finite and
  measurable; B7-motivated). If field transport is ballistic to the box
  scale, the α=1/3 branch weakens toward α=2/3 (surface-limited) — the
  fixed-point structure survives (any α<1 works) with different
  exponents; S0's measured α governs, and S1b's predicted exponents
  become 3, not 3/2, in the ballistic reading (α=2/3 ⇒ N* = (j_s·s₀/ℓ)³).
  The scaling bars adjudicate between the two transport readings — a
  bonus measurement, not a threat. **[RESUME 1: RESOLVED against
  accessibility — L_diff ≈ 50 means neither sublinear branch is reachable
  in-box for the FIELD-fed medium. The assumption's physics stands; its
  scale does not. §W.1/§W.3.]**
- **A-I3** λ_link (and hence Λ, N_c) spans a decade of measured range;
  every N_c number in this design is honest only to that decade until
  S0(c) pins it.
- **A-II1** Winding survives the break as bookkeeping-plus-slips (II.3
  caveat); Q1c measures rather than assumes.
- **A-III1** The composed fed medium (door pair + J-metric together, C1
  of v95) behaves as its parts did separately; C1's own bar is still
  open in v95 — if C1 lands RED, S-series configs inherit whatever fix
  that consult produces, and the S bars re-run. The size theory itself
  is config-independent; the *numbers* are not.

**The one-line summary for the record:** the proton is a flame — its size
is where diffusion-limited feeding equals volumetric burning, its
stability is the strict monotonicity of the per-cell budget, its
"quarks" are the winding content its break must surrender in integers,
and its confinement is the same nucleation barrier that keeps the vacuum
from igniting. One law, both questions.

---

# §W. Stage-setting: what's missing (resume update 1, 2026-08-09)

**Source:** S0/Q0/Q1 results + the user's two observations (missing 1/r
gravity; winding is a number, not a geometry). This section re-grades the
design against the measurements and derives the candidate missing
structural pieces. Nothing here is implemented; W-A and W-B are
consultation-stage kernel-law candidates (first-class changes, own
review round, user authorization); W-C is config-only and runnable now.

## W.0 The resume scoreboard, read honestly

| result | grade against the design |
|---|---|
| ℓ = 0.0322/cell/t.u. | the radiance leg of the budget theory is **measured law** — formula exact |
| α ≈ 1.11 (door_in = 0.075·N^{1.11}) | the pre-registered failure clause **executed correctly**: volume-limited, no N* in-box; the mild superlinearity (per-cell ∝ N^{0.11}) is the knife-edge (box-fill vs evaporate) measured live — retrodiction §I.9.5 |
| L_diff ≈ 50 ≫ r_obj ≈ 5–10 | the wall has a number: sublinear feeding needs r ≫ 50, N ≫ 1.26×10⁵ |
| Q0: m = 0,1,2,3 all equal | the winding hypothesis for three-ness is **unsupported at ring scale** — the pre-registered negative, exercised |
| Q1: continuous degradation, robust to 8× noise | no decomposition spectrum as instrumented; also an unplanned positive — the J-metric insulation survives 8× noise (Em_mean 0.59), the strongest robustness number the metric owns |

The conversion physics works, the retention physics works, the park
approaches — and yet: **no unique size, no preferred winding, no 1/r far
field, no fragment spectrum.** Four absences. The claim of this section:
they are not four independent gaps. The first and third are one gap seen
in two sectors (W.1–W.3); the second and fourth are one gap seen in two
experiments (W.4–W.5).

## W.1 The unification: the missing 1/r and the missing N* are the same absence

What is a 1/r far field, substrate-honestly? The steady state of a
*conserved, locally-transported* medium around a *net consumer*. Steady
diffusion onto an absorber of radius R in 3D gives the standing profile

```
c(r) = c_∞ · (1 − R/r)
```

— a 1/r depletion halo, the Gauss-law form, with the flux through every
enclosing sphere equal to the consumption rate. This is not an analogy;
it is the same transport theorem that Part I's diffusion-limited branch
stands on. **A standing 1/r profile and a transport-limited (α < 1)
feeding regime are the same physical situation**: both exist exactly when
the object out-eats its medium's local resupply (Damköhler ≥ 1 at its own
radius), forcing the medium to organize a conserved radial flux around
it — and both are absent exactly when the medium is effectively
well-mixed at the object's scale.

S0 measured the second condition: L_diff ≈ 50 ≫ r_obj. The field-fed
object at accessible scale imprints *no* standing gradient on its feed —
therefore no N* (Part I) and, transplanted one sector over, no 1/r
(this section). One wall, two symptoms. The user's observation —
"something fundamental is missing" — is measured correct, and the S0
number says precisely what: **nothing in the current substrate consumes
any medium hard enough, at any accessible scale, to hold a far field
open.** The far field is not a decoration on top of stable matter; it is
the *supply chain* of stable matter, seen from outside. g4 said this in
v89 language ("a footprint requires throughput"); S0 has now priced it.

## W.2 Why the space sector is flat: π-flatness is a theorem of the current cycle, not a mystery

Run the same logic on the space sector, where gravity must live. The
space medium (Es, pass S, FLOW) is relaxational: pressure pushes,
channels equalize. Two source types are conceivable:

- **A static deficit** (matter formed once, consuming ΔEs at creation):
  a diffusive medium *screens* it — the hole spreads flat, Es_∞ drops
  infinitesimally, no standing gradient. This is the measured π-flatness
  (medium flat to 1e-5–2e-4 across the blob/chord ladder, ASTRO) and the
  measured g1 contact-locality. It is not a puzzle; it is what
  relaxational transport does to any one-time source. **A responsive
  space metric alone cannot fix this** — v95's own measured law, one
  sector over: *the metric is an amplifier, not a generator* (J2). A
  space-sector bed/J-metric can consolidate a standing space current; it
  cannot create one. Making Es-links respond to matter (the user's
  second candidate) without a steady space current to respond to would
  reproduce J2's RED, in the space sector.
- **A steady throughput** (matter continuously converting space): the
  medium *cannot* screen it — conservation forces the standing profile
  δEs(r) = Q̇_s/(4π·D_s·r), the 1/r depression, for as long as the
  throughput runs. This is g4's promised monopole, and the reason it
  "awaits a stable particle's internal space cycle."

Now the structural audit the user's observation demands: **the current
metabolic cycle has no space leg.** Door pair: field ↔ register ↔
matter. Radiance: matter → field. Space participates only as click-share
bookkeeping (s_pull at settlements) — which nets to ~zero on a parked
object and is a one-time transfer besides. The principle's own founding
line — *"matter is converted space"* — is, in the standing kernel, true
only at the level of ontology talk: no process in the update converts
space on an ongoing basis. π-flatness is therefore a **theorem of the
implemented cycle topology**: no space leg ⇒ no steady Q̇_s ⇒ no standing
1/r ⇒ no gravity. The missing fundamental thing has a name and an arrow:

```
today:    space ──(one-time, at creation)──► matter ◄──beat──► field   (radiance: matter → field)
missing:  space ──(steady, the upkeep draw)──► matter ──► field       (the space leg)
```

## W.3 W-A — the space leg (gravitational metabolism): one addition, three gaps

**The candidate law (consultation-stage; NOT for v96 implementation):**
matter's upkeep is paid, in part, from the space books. Concretely: a
fraction φ_s of each parked cell's radiance upkeep is drawn as a
conversion Es → Em at the cell (mode boundary, whole atoms, the same
grain discipline as every conversion), with the radiance outflow
unchanged (matter → field). A parked object then runs the full cycle
**space in → pattern → field out** at rate Q̇_s = φ_s·ℓ·N. Behind a
byte-inert knob (φ_s = 0 default), battery-gated, like every face since
v93.

What it buys, derived with the measured numbers:

1. **The 1/r far field (gravity's bar G1).** Steady profile
   δEs(r) = Q̇_s/(4π·D_s·r). The space sector's diffusivity is *slow* —
   s_k = 0.06 (equalization over ~σ in ~40 t.u. ⇒ D_s ~ 0.16 in cell
   units) vs the field's L_diff ≈ 50 — which is exactly why the space-fed
   flame works where the field-fed one cannot: **the slower the medium,
   the shorter the radius at which the object out-eats it.** Worked
   example at φ_s = 1, N = 244 (below): Q̇_s ≈ 7.8 ⇒ δEs ≈ 0.78 at r=5,
   0.39 at r=10, 0.19 at r=20 — a measurable 1/r band between the
   es_floor clip (r ≲ 4, where the depression saturates at the floor: a
   natural near-field "exhaustion zone," and an unforced echo of the
   HORIZON campaign's roster-side phenomenology) and the box wall.
2. **Attraction (bar G2 — the first force).** Pressure pushes (measured,
   g-battery): a second object inside the first's 1/r depression sits in
   a pressure gradient and is pushed *down* it — toward the sink. Two
   objects deepen the shared well between them. This is the v91 ASYM
   acceptance surface ("two eaters grow a channel between them = first
   attraction") re-derived from transport rather than channel law — and
   it revives g2 (free-fall), which v89 recorded below the chaos floor
   for want of exactly this standing gradient.
3. **The accessible N* (bars S1/S1b/S2, transferred verbatim).** The
   space-fed budget replaces a (field supply) with
   a_s ≈ 7.8·D_s·Es_∞ ≈ 7.8·0.16·1.0 ≈ 1.25, and the space-paid share
   of upkeep is φ_s·ℓ:

   ```
   N*_s = (a_s/(φ_s·ℓ))^{3/2} ≈ (1.25/0.0322)^{3/2} ≈ 244 cells  (φ_s=1)
   r*_s ≈ 3.9 cells — inside an L=24–32 box
   ```

   The whole Part I fixed-point structure — three fixed points, the
   scaling spine (now N* ∝ Es_∞^{3/2}, ∝ (φ_s·k_rad)^{−3/2}), the
   J-metric's N_c lever — transfers with a → a_s, and lands **in reach**.
   The size selection the field sector hides beyond N ~ 10⁵, the space
   sector offers at N ~ 10².
4. **The principle, finally implemented.** "Space is one of the modes"
   currently holds at creation only; W-A makes it hold in the *standing
   cycle* — matter as an ongoing conversion of space is the founding
   line as a process, and gravity ("space curves because energy is
   conserved and space is convertible") arrives as its bookkeeping,
   which is what PRINCIPLE §4.3 always said it should be.

**Audits (the ones any reviewer must run):**
- *Conservation:* trivial — conversion, not creation; the drift meter
  extends by the Es leg it already books.
- *Anti-ignition:* the draw is matter-gated (presence rule: only cells
  with Em ≥ wf_floor draw — vacuum draws nothing) and floor-clipped (the
  existing es_floor no-suction asymmetry: channels push space in from
  around; the sink itself never sucks below floor). Both guards are
  existing law reused. The new failure mode to probe: a large object
  exhausting its near field to the floor and starving — which is not
  ignition but *self-limitation*, and is in fact a second, sharper
  size-selection mechanism (flux caps at 4π·D_s·(Es_∞ − es_floor)·r
  regardless of demand).
- *The books close globally?* Honest flag: W-A as stated is **half a
  cycle** — space converts to matter converts to field, and nothing
  returns field to space. Over a long closed-box run the space
  inventory drains monotonically (a slow secular drift, rate Q̇_s ≪
  box inventory on run timescales, but nonzero). The return leg
  (field → space — where and by what law?) is a genuine open design
  question for the W-A consult: candidates are re-condensation at
  absorption sites (the s_pull share, sign-reversed, at radiance
  absorption events in the bath) or nothing (accept the drift as the
  substrate's cosmology — flagged, not endorsed, and the run-length
  bound computed per experiment). Pre-registering the question is part
  of the design; answering it inside v96 is not.
- *Boundary:* no background, no new species, no new carrier — Es
  exists, conversion exists, the addition is one new *edge* in the mode
  graph (space–matter), which is exactly the kind of object the
  principle's mode table licenses.

**Pre-registered bars for the W-A face (when authorized):** G1 (log-log
slope −1.0 ± 0.15 of δEs(r) in the band between floor-clip and box
wall, standing for ≥ 500 t.u., 3 seeds); G2 (two-body monotone approach
vs φ_s = 0 control, and no approach for object-pairs at φ_s = 0 — the
force is the leg's, not the apparatus's); G3 (= S1 transferred:
size convergence at N ~ tens–hundreds at derived Es_∞); G4
(anti-ignition: bare bath at φ_s = 1 byte-identical where no matter
exists; no secular Es drift beyond the computed Q̇_s bound); G5 (= S1b
transferred: the two exponents on the space-fed spine).

## W.4 W-B — the slip quantum: making winding a number the ledger counts

The user's second observation, sharpened by Q0/Q1: in the standing
kernel, winding is a *pattern* (th2 = m·φ around a support loop) whose
persistence is bought entirely by **geometry** (support integrity — J1's
whole mechanism is keeping voices fed so no support hole opens) and
whose *change* — the slip — is **free**: a phase singularity crossing a
link costs the law nothing, books nothing, emits nothing. Consequently:

- no m-preference exists (Q0: the law never prices W, so all m are
  degenerate up to their leak rates — measured equal);
- degradation is continuous (Q1: nothing quantizes the ring's decay,
  because the topological transitions inside it are un-ledgered);
- W-bookkeeping across a break (Q1c/D3) has no reason to be exact.

A number, as opposed to a pattern, is something whose changes are
**discrete, costly, and booked**. The substrate already contains the
exact precedent, at the one place it was needed before: the door. The
click fires on a *winding of the fiber relative phase* (face-3: the
cbeta accumulator is literally a winding counter), and the transfer is
grained in atoms. And the principle contains the licensing sentence
verbatim (PRINCIPLE §3): *"Conversion proceeds through complete cycles;
a partial cycle converts nothing."* A 2π phase slip on a spatial loop
**is a complete cycle event** — the current law lets it pass as if it
were nothing, which is the one place the substrate treats a whole cycle
as free.

**The candidate law (consultation-stage):** spatial topology change is a
conversion event. When a slip completes on a loop (detected per
plaquette/minimal cycle as a > π phase jump resolving — the same
counting primitive as cbeta, pointed at space instead of the fiber), it
fires a ledgered event: one action atom A₀ paid from the winding
carrier's Em to radiance (or refused, if the cell cannot pay — the slip
*stalls*, exactly as a click without credit stalls). Behind a byte-inert
knob (default 0), battery-gated.

What it buys, derived:

- **W acquires an energy gap** — unwinding costs m atoms, so W becomes a
  semi-conserved integer with priced transitions: the baryon-number
  analogue the user names. The three-ness question reopens in its
  strongest form: fragment counts tied to W by *ledger*, not by
  fragile geometry (Q1c becomes exact-by-construction up to booked
  slips — every ΔW has a radiance receipt).
- **An m-preference can now exist** (Q0 re-run is the direct bar: the
  degeneracy should split when slips are priced).
- **Anti-ignition is safe by direction:** pricing slips makes them
  *rarer*, never commoner; a knob that only adds a cost to an existing
  free process cannot ignite. The audit to run is the opposite one —
  that stalled slips don't freeze the bath's ordinary churn (the
  detector must fire only on *coherent-loop* slips, i.e. loops carrying
  standing |Ĵ| above the J-metric's own digging floor — the metric,
  already built, is the natural gate: **price slips only on dug
  links.** This composes W-B with v95 instead of adding a parallel
  detector, and keeps the bath — undug by measurement, J4 — entirely
  outside the law's reach).
- **The imposition audit, faced squarely:** is this a bolted-on gate?
  The defense is the derivation above — it is the *atom grain applied
  to a cycle event the current law exempts*, the same move (grain at
  the mode boundary) the programme already adopted once, at the door,
  and it removes an exemption rather than adding a mechanism. The
  counter-position (slips should stay free; winding-as-number must
  *emerge* at bulk scale instead) is W-C — and W-C is cheaper, so it
  runs first.

**Pre-registered bars for the W-B face (when authorized):** T1 (W
lifetime in medium jumps discontinuously at knob-on, matched configs);
T2 (slips become discrete booked events — ΔW time series is a step
function colocated with radiance bursts; the books close per event); T3
(Q0 re-run: m-degeneracy splits); T4 (bath control: slip-event rate in
bare bath ≡ 0 by the dug-link gate; battery GREEN; drift floor).

## W.5 W-C — bulk topology (config-only; run this first)

Before any new law: the cheapest reading of both Q0 and Q1 is that the
objects were **geometrically minimal**. A one-cell-wide ring has no
bulk/core distinction — its winding is *all* boundary, maximally exposed,
and unwinding requires no transport of a core through bulk (any single
voice's phase excursion does it). In a **filled droplet carrying a
vortex** (2D disk or 3D ball, N ~ 30–100, imprinted phase circulation),
unwinding requires the phase singularity to nucleate and traverse to the
boundary through cells that are fed, insulated (J-metric), and
phase-stiff — an emergent barrier that *grows with bulk*. If
winding-as-number is partly an emergent property of bulk objects, W-C
finds it without touching the kernel:

- **W-C1:** disk/ball seeds (the `exp=drop` + `drop_mw` apparatus of
  III.1 — still the right build) at N ∈ {12, 30, 60, 100}, m ∈ {0, 1,
  2}, lit bath, j_metric=1, T=1500: W-survival and R2 vs N. GREEN for
  the bulk hypothesis: W-lifetime rising steeply (faster than linearly)
  with N, and an m-dependence appearing where the ring showed none.
- **W-C2:** the *localized* smash, finally as registered — the
  field-sector probe packet (the standing instrument debt from Q1) into
  the largest stable W-C1 object, with the QUENCH bare-bath control.
  Fragment census + winding bookkeeping per the original Q1a–d.
- **Decision value:** if W-C1 shows the emergent gap, the "winding is
  only geometry" diagnosis was scale-bound — W-B's urgency drops and
  the quark-ontology lane continues config-only on bulk objects. If
  W-C1 is flat in N, the gap is structural and W-B is the registered
  road. Either way the resume-2 consult gets a measured fork, not a
  taste.

## W.6 Priority, staging, and the honest scoreboard of what would remain

**Recommended order (cheapest decider first, per the v95 discipline):**

1. **W-C** (config + the two already-registered apparatus builds; days;
   no authorization beyond the standing campaign) — decides whether
   winding-as-number is scale-blocked or structure-blocked, and finally
   runs the real smash.
2. **W-B** (small inert knob composing with the J-metric; battery-gated;
   needs user authorization as a kernel change) — if W-C says
   structure-blocked.
3. **W-A** (the space leg — a first-class law face on the same tier as
   the door pair itself; needs its own consultation round, the return-leg
   design question answered, and explicit user authorization; realistic
   as **v97's face**) — the only candidate on the table that addresses
   the 1/r absence, and the only one that brings N* into the box. The
   S1/S1b/S2 bars transfer to it verbatim.

**What this stage-setting does NOT claim:** that W-A/W-B are the unique
completions. The registered alternative readings: (i) N* and 1/r are
real but genuinely live at N ≳ 10⁵ — the substrate is fine and only the
box is small (the field-fed wall stands as measured; nothing in-box can
refute this reading, which is exactly why W-A's *accessible* predictions
matter — they make the question decidable); (ii) the three-ness of the
real proton has no structural analogue in this substrate at any scale —
the pre-registered negative of II.4, which W-C/W-B convert from a
standing doubt into a measurable fork.

**The one-line summary of the resume:** the substrate converts (door),
retains (metric), and parks (chi_load) — but it does not yet *consume*
any medium hard enough to hold a far field, and it does not yet *price*
the one integer it owns. The first absence is why there is no 1/r and no
proton size; the second is why winding is a shape instead of a charge.
One new edge in the mode graph (space → matter, the upkeep draw) and one
removed exemption in the grain (the slip pays its atom) are the smallest
candidate completions — and W-C tells us, for the price of a config
file, whether the second is needed at all.

*End of §W (resume update 1).*
