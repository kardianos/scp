# v94 developmental design — the translation primitive, mirrored to the door

**Author:** claude (fable), developmental session, 2026-08-09.
**Deliverable for:** PART III of `v94/consult/BRIEF.md` (the EMF translation
primitive mirrored to matter intake/outtake).
**Grounded in:** `v94/kernel/freecell.c` (pass F :2151, pass U :1591, pass 6
:2198, local clock :1386, `apply_uhop` :1265), `v93/theory/
FINAL_STATUS_ADOPTION.md`, `v93/RECKONING.md`, `v93/RING_DNLS.md`,
`v93/theory/METABOLIC_MATTER.md`, `v94/laws/laws_V3a.cfg`.
**Status:** design only. No kernel code here. Everything specified behind a
new byte-inert knob, default 0, battery-gated.

---

## Executive summary (the design in six sentences)

The field sector already translates by one primitive: **two first-order,
mutually-inducing, incomplete conversions, interleaved** — the quadrature
pair `(fa1, fa2)` under local precession + pairwise Givens hops (pass F).
The dense sector got the *within-mode* copy of that primitive in v93 (pass
U) and it worked exactly as far as the primitive reaches: transport yes,
mode-boundary no. **The v94 claim: the door is a third link type — a link
in the *mode* direction instead of the *space* direction — and the same
primitive applies to it.** Today's door is the stroboscopic, magnitude-only
shadow of that link: the beat accumulator `cbeta += (w1e − w2e)·dt` is
already the carrier frequency of a two-level field↔matter rotation, but the
kernel samples it once per cycle, transfers a magnitude, and discards the
phase cross-term — *structurally the same amputation* that made the additive
dense hop reject the interference term before pass U restored it. The design
below restores the cross-term at the door: a **door register** `χ` (the
fabric's polarization — energy in mode-transit, the exact cousin of the link
flight parcel `slem`, which is energy in space-transit), coherently and
unitarily coupled to the field spinor between clicks, settling into dense
matter **only in quantized action atoms at the beat** — the click grain is
never unitarized; the pair lives *between* clicks, exactly as the brief
hypothesizes. The payoff, derived below: refraction (observation 2) becomes
a theorem instead of an absent mechanism, the workfn law becomes detuning
physics instead of a threshold rule, the metabolic cycle becomes a
propagating wave (the matter analogue of the EM wave — a polariton on the
roster), momentum crosses the door (the standing ~100× radiation-pressure
deficit is engaged), and the binding/mobility antagonism has a structural
release: **the object translates by cycling, not by hopping** — motion as
regular conversion, the programme's founding line, finally implemented
literally.

---

# III.1 — Theory: the EMF translation primitive, fabric-derived

## 1.1 What the kernel's pass F actually is

Write the field spinor as one complex number per cell,

```
ψ_e(i) = fa1(i) + i·fa2(i),        Ee(i) = |ψ_e(i)|².
```

Pass F (freecell.c:2151) is exactly two operations per step:

**(a) Local precession** (freecell.c:2152–2159):

```
fa1 ← cc·fa1 + ss·fa2
fa2 ← −ss·fa1 + cc·fa2         with (cc,ss) = (cos, sin)(w1e·dt)
```

which is `ψ_e ← e^{−i·Ω₁·dt}·ψ_e` with `Ω₁ ≡ w1e = w1/(1 + q_detune·x)`
(the field clock is **already load-detuned** — freecell.c:1387; hold that
fact, it is load-bearing for observation 2 and for the door gate in §III.2).

**(b) Pairwise Givens hop** (freecell.c:2169–2194): for each live link
`s = (i,j)` in canonical order, angle `τ_s = field_J·w_s·dt/√(s_i·s_j)`:

```
fa1_i ← cc·fa1_i + ss·fa2_j        fa1_j ← cc·fa1_j + ss·fa2_i
fa2_i ← cc·fa2_i − ss·fa1_j        fa2_j ← cc·fa2_j − ss·fa1_i
```

In spinor form this is precisely

```
ψ_i ← cos(τ)·ψ_i − i·sin(τ)·ψ_j
ψ_j ← cos(τ)·ψ_j − i·sin(τ)·ψ_i
```

i.e. `exp(−i·τ·σ_x)` on the two-cell subspace. The generator of the whole
pass is Hermitian:

```
i·(dψ/dt) = H·ψ,     H = diag(Ω₁(i)) + K,     K_ij = field_J·w_ij/√(s_i·s_j)
```

— a discrete Schrödinger evolution on the roster graph. Norm (= energy) is
conserved per rotation, pairwise, exactly: **conservation is a theorem of
the update.** This is the fact v93 transplanted within-mode; v94 transplants
it across-mode.

## 1.2 The two expressions and their mutual induction

The two expressions of the field (observation 3) are the **quadratures**
`fa1` and `fa2`. In quadrature form the generator reads, per cell i:

```
d(fa1_i)/dt  =  +Ω₁·fa2_i  +  Σ_j K_ij·fa2_j        (fa1 sourced ONLY by fa2)
d(fa2_i)/dt  =  −Ω₁·fa1_i  −  Σ_j K_ij·fa1_j        (fa2 sourced ONLY by fa1)
```

Compare the Maxwell shape:

```
∂E/∂t = +c²·∇×B          ∂B/∂t = −∇×E
```

The correspondence is structural, and exact in the sense that matters:

- **first-order**: no second time derivatives anywhere;
- **incomplete**: neither quadrature has a self-term — `fa1`'s entire time
  derivative is sourced by `fa2` (on-site and neighbors), and vice versa
  with the opposite sign. Freeze one quadrature and the other is inert;
- **mutually inducing**: the pair, eliminated to second order, gives
  `d²ψ/dt² = −H²·ψ` — the wave equation on the roster graph.

Caveat stated honestly: the continuum curl is the 3D antisymmetric
derivative acting on vector fields; the roster field is one complex scalar
per cell, so the neighbor coupling is the symmetric graph operator `K`, not
a literal curl. What carries over — and what the user's observation names —
is the **shape**: two first-order cross-sourced halves whose interleaving is
propagation. (The analytic-signal decomposition of a 1D Maxwell wave gives
exactly this quadrature pair.)

**Why each half alone does not translate.** Precession alone (`K = 0`):
every cell spins in place, flat band, group velocity zero — oscillation, no
transport. Hops alone (`Ω₁ = 0`): a generic (single-quadrature, zero-phase-
gradient) seed spreads **symmetrically** in ±k — dispersion, no direction.
A *directional* wave is a phase gradient (momentum) being advanced by the
precession while the hop converts phase gradient into displacement of the
envelope: `v_g = dω/dk` with `ω(k) = Ω₁ + ε(k)`, `ε(k)` the graph band of
`K` (on a regular chain `ε(k) = 2·K̄·cos(k·d)`, `|v_g| ≤ 2·K̄·d`). The
substrate's `c` — observation 1 — is the maximum band slope

```
c_sub ≈ 2·K̄·d̄ = 2·field_J·w̄·d̄/√(s̄²)
```

fixed by law constants and roster geometry, no import.

## 1.3 The refractive effect (observation 2), in fabric terms

Two mechanisms, one present, one missing:

**(i) Load-detuned on-site pitch (in the kernel today).** `w1e =
w1/(1 + q_detune·x)` (freecell.c:1386–1388): loaded matter lowers the field
precession. A spatial region of load is a step in the on-site term of `H` —
a refractive **potential**. A packet at conserved total ω entering the
region must re-solve `ω = Ω₁(x) + ε(k)` with a shifted `Ω₁`, changing `k`
and hence `v_g` — reflection, refraction, impedance mismatch. This is real
but is the *smaller* half of refraction.

**(ii) Conversion mixing (absent today — supplied by this design).** In
continuum physics the refractive index of matter is *not* mainly a
potential effect; it is the field continuously and incompletely converting
into matter polarization and back (the Lorentz-oscillator / polariton
mechanism). Light through glass is slow **because it alternates between two
incomplete conversions with the medium** — which is the user's observation
2 *and* observation 3 read together, and it is precisely the mechanism the
kernel lacks: today's door is a discrete click, so the field in loaded
matter never virtually dresses itself with matter. §III.2's pair supplies
exactly this, and the slow-down becomes derivable (§III.3.4): the group
velocity of the mixed wave is the bare velocity times the field weight of
the mixed mode,

```
v_mixed = v_bare · |⟨field|polariton⟩|²  <  v_bare   inside loaded matter.
```

So: observation 1 = the pair in vacuum at the band slope; observation 2 =
the same pair dressed by the door coupling inside matter; observation 3 =
the pair's two-quadrature structure itself, now to be mirrored at the door.

---

# III.2 — Mirror to matter: the central design

## 2.1 Candidate (A): within-mode extension, door as boundary condition

Under (A), the matter pair is `(m₁, m₂)` under pass U — already built, and
this is the correct within-mode statement: pass U (freecell.c:1591) is
byte-for-byte pass F's algebra on `ψ_m = √Em·e^{iθ₂}`. The door would then
be a "boundary condition" of that primitive — an absorbing/emitting edge in
the mode fiber.

**Evaluation: (A) alone cannot engage the wall.** The measured v93 wall is
*at the door*: (i) the door scrambles phase (QUENCH-3 `R2d ≈ 0.02`; the
`amp_door`/`qp_phase` heuristics — trig-formula phase mixes written into
`th2` — never delivered retention); (ii) the metabolic cycle never closes
(`rad_dual` gives SLOW_REFILL, no `Em*` park; the coherent cross-term `g` —
FINAL_STATUS item 10 — was never coded); (iii) treated as a *pure* boundary
condition of a unitary within-mode pair, the door would have to be
unitarized — continuous `Ee ↔ Em` transfer — which is **forbidden** (III.7:
the click grain stays; and rightly so: ħ-linearity, the two-atom credit,
and the entire battery's atom bars live on the quantized door). (A) is
necessary and is *done*; it is not the open problem.

## 2.2 Candidate (B): mode-boundary extension — the door as the third link

The substrate has, after v93, **two kinds of links carrying the same
primitive**:

| link type | endpoints | primitive | pass |
|---|---|---|---|
| spatial, field | (cell i, cell j) within field mode | precession + Givens | F |
| spatial, dense | (cell i, cell j) within dense mode | precession + Givens | U |
| **modal (this design)** | **(field, matter) within cell i** | **precession + Givens, settling in atoms** | **D (new)** |

The claim: **conversion is transport in the mode direction**, and the
translation device for it is the same device. The door pair is:

```
i·dψ_e/dt  =  Ω₁·ψ_e  −  g·χ           (field sourced by the register)
i·dχ/dt    =  Ω₂·χ    −  g·ψ_e         (register sourced by the field)
```

where `χ = dc1 + i·dc2` is the **door register** — a per-cell two-component
carrier holding energy **in mode-transit** (defined fully in §2.3), `Ω₂ =
w2e` (the register is matter-pitched: it is matter-in-the-making), and `g`
is the door coupling (its scale fixed in §2.4 with no imported constant).
Each equation is first-order; each is incomplete (the coupling term has no
self-part); each sources the other — the Maxwell shape, verbatim, in the
mode direction. The role-mapping requested by the brief:

```
∂E/∂t = +c²∇×B    ↦    dψ_e/dt ∝ −i·g·χ     (outtake induction: register drives field)
∂B/∂t = −∇×E      ↦    dχ/dt   ∝ −i·g·ψ_e   (intake induction: field drives register)
```

**Intake and outtake are the two half-cycles of one rotation.** The pair
`(ψ_e, χ)` is a driven two-level system with detuning equal to the **beat**:

```
Δ(x) = Ω₁ − Ω₂ = (w1 − w2)/(1 + q_detune·x)
```

— the *identical* quantity the kernel already integrates in `cbeta`
(freecell.c:2241) to gate the door. The occupancy oscillates at the
generalized Rabi frequency `Ω_R = √(Δ² + 4g²)` with maximum transfer
fraction

```
F(x) = 4g² / (Δ(x)² + 4g²)   <  1        (conversion is INCOMPLETE),
```

and the direction of net transfer alternates each half-period — **two
naturally alternating conversions, each incomplete**, which is observation
3 restated at the door. Today's pass-6 door is what remains of this pair
after two amputations: sample only once per beat period (stroboscopic), and
keep only the magnitudes (discard `Im(ψ_e*·χ)`). The v93 diagnosis of the
additive dense hop ("the additive Em-ledger rejects the interference term")
applies to the door verbatim; this design is the same repair, applied at
the same joint, one level up.

## 2.3 The register χ: what it is and why it is not a fourth mode

`χ` is **energy in transit across the mode boundary** — the exact cousin of
the link flight parcel `slem` (energy in transit across space, belonging to
neither endpoint cell's inventory, advanced by `slph += dt·C/d` at the
substrate speed, settled on arrival). The symmetry:

| | space-transit (exists) | mode-transit (this design) |
|---|---|---|
| carrier | `slem[2s+dir]` (+ Phase-M shadow `sre_,sim_`) | `χ(i) = (dc1, dc2)` |
| clock while in transit | flight phase advance, transit rotation | precession at `Ω₂` |
| crossing time | `d/C` | beat/Rabi period `2π/Ω_R` |
| settlement | delivery on `slph ≥ 1`, capacity-checked | **click in atoms `ε` at the beat** (§2.5) |
| conservation | flight ledger in the conserved sum | `|χ|²` in the conserved sum |

It is **not** a fourth energy mode: it introduces no new conversion
channel-ends (everything entering `χ` comes from `Ee` unitarily and leaves
to `Em` in atoms or back to `Ee` unitarily), no new pitch law, no transport
of its own (`χ` does not hop — it is strictly cell-local, the flat band of
§III.3.4), and its energy is booked explicitly in the conserved total
`Σ(Em + Ee + Es + |χ|² + flight)`. It is door-flight. **Flag per III.8:**
it *is* two new persistent floats per cell, and per programme rule the user
should sign off on that explicitly; the alternative (a registerless direct
`ψ_e ↔ ψ_m` rotation) is evaluated and rejected in §2.6.

QUENCH-3 compliance: the phase of the conversion is carried **in `χ`**,
slot-borne in the sense that matters — a dedicated carrier, never a
heuristic write into the cell clock `th2`. (Settlement composes into the
dense spinor by a norm-preserving rotation, the same operation class as
pass U's own write-back; see §2.5 and the explicit flag in §III.8.)

## 2.4 The coupling scale g: the workfn law derived

No SI constant and no new frequency scale may be imported. The candidates:

- `g ∝ Δ(x)` (beat-scaled): then `g/Δ` is load-independent and the door is
  equally open in vacuum and matter — **rejected** (ignition-unsafe,
  physics-wrong).
- `g ∝ w2e` (detuned-pitch-scaled): both `g` and `Δ` scale as
  `1/(1+q·x)`, ratio again load-independent — **rejected**, same reason.
- **`g = door_g · w2` (bare-pitch-scaled) — adopted.** The door's stiffness
  is a property of the *fabric* (the bare law constant `w2`), not of the
  local load; only within-mode clocks are load-detuned. Then

```
g/|Δ(x)| = door_g · w2 · (1 + q_detune·x) / (w2 − w1)
```

**grows with load**: vacuum is detuned away from its own door; bound matter
is pulled toward door resonance. Numbers at the law table (`w1 = 1.65`,
`w2 = 2.9`, `q_detune = 1.2`): `|Δ(0)| = 1.25`; at the V3a radiance fixed
point `x* ≈ 0.62`: `|Δ| = 0.717` — the transfer fraction `F(x)` rises ≈3×
from vacuum to law-loaded matter (small-g limit), monotonically in `x`.
**This is the workfn law — "empty foam demands a far threshold, bound
matter converts freely" — re-derived as detuning physics** rather than
imposed as an `econd` rule. It is a soft gate, so the adopted hard workfn
stays layered at settlement (§2.5); but the design predicts the workfn
*shape* from the law's own constants, which is a falsifiable claim (bar
M3's dark-bath arm tests its safety; a dedicated `F(x)` sweep can test its
shape).

## 2.5 Settlement: the click grain survives (the pair lives between clicks)

The `(ψ_e, χ)` rotation is continuous and reversible — **virtual**
conversion. Real conversion — the ledger move `|χ|² → Em` — happens **only
at the beat, only in whole action atoms**, exactly on today's machinery:

- **When:** `beat_fire` (existing `cbeta` accumulator, unchanged).
- **Gate:** the adopted workfn text, unchanged (`econd_i` from `wf_on/
  wf_floor/wf_far`), applied to the register inventory `|χ|²` in place of
  `Ee`.
- **Grain:** demand `d1 = f_conv·(|χ|² − econd_i)`, then `atoms_fire` /
  `atoms_clamp` with the same `eF = A₀·w1e/2π`, `eD = A₀·w2e/2π` credits —
  ħ-linearity and the two-atom credit preserved verbatim.
- **How (the coherent part):** the transfer executes as a **Givens rotation
  between `χ` and `ψ_m`** (the dense spinor reconstructed from `(Em, th2)`
  as pass U does) by the angle `θ` that moves exactly the fired energy:

```
ΔEm(θ) = sin²θ·(|χ|² − Em) + sinθ·cosθ·C,     C = 2·Im(ψ_m*·χ)
```

  solved for `ΔEm(θ) = d1` on `θ ∈ (0, π/2]`; integrate-and-fire honestly —
  if the geometry cannot move a whole `d1`, clamp to the largest reachable
  whole-atom multiple, else no fire. Norm `Em + |χ|²` is exactly conserved
  by the rotation; the fired energy is exactly quantized; **the phase
  crosses the door as algebra, not as a formula patch.**
- **Direction (stimulated outtake for free):** the cross-term `C` is
  signed. When the phase geometry drives `ΔEm < 0` — matter in anti-phase
  with its own register — the same solved rotation moves atoms **out** of
  `Em` into `χ` (fired on the outtake atom credits, ledgered separately),
  and the register unloads to the field through the pair. This **is** the
  coherent cross-term `g·x^{p−1}·xb` of qgraph Theorem 3.5 (FINAL_STATUS
  item 10, "highest-upside unfinished, not coded") — arriving not as a rate
  formula but as the phase geometry of one primitive. Spontaneous radiance
  (`k_rad·cap·x⁴`) is untouched: it remains the adopted spontaneous-
  emission channel; the pair adds the stimulated channel beside it, exactly
  the spontaneous/stimulated split of real emission.
- **Space book:** `s_pull` backsplash bookkeeping at settlement identical
  to the legacy branches, per direction. `Es` never enters the continuous
  rotation — the space footprint moves only at clicks, so the space ledger
  keeps its click grain too.

Between clicks: unitary, reversible, phase-carrying. At clicks: quantized,
gated, ledgered. **The pair lives between the clicks; the click is the
integrator-of-fire for the pair** — the brief's requirement, met by
construction.

## 2.6 Why not a registerless direct pair (ψ_e ↔ ψ_m)?

Rotating the field spinor directly against the dense spinor every step
would (i) transfer `Ee ↔ Em` continuously — unitarizing the door, breaking
the atom bars and III.7; (ii) leave vacuum with no virtual/real distinction
— every fluctuation converts a little, ignition-unsafe; (iii) couple the
fast field directly to the `q_detune`-trapped dense hoard, re-importing the
antagonism. The register is what makes "virtual conversion" a state the
substrate can hold: incomplete conversion **needs somewhere to stand**.
That is arguably the deepest reading of observation 3 — the two alternating
expressions exist *because* conversion has a mid-point, and real EMF in
matter is field + polarization, not field + settled matter.

## 2.7 The verdict the brief asked for, and the antagonism

**Both candidates are correct at their own scale, and they compose.** (A)
within-mode pairs (F, U) are the primitive on spatial links — done. (B) the
door pair is the *same* primitive on the modal link — this design. The door
is *neither* a boundary condition of (A) *nor* a replacement of the click:
it is the third edge-type of one graph (§III.4's prism complex), with the
unique feature that its settlement is quantized.

**The antagonism (RING_DNLS B.3: the corner that binds does not
translate).** Measured: `q_detune` self-trapping binds `Em` by resonant
decoupling — a loaded cell detunes from its neighbors and its *dense hops*
die; a moving packet condenses and stalls (COM advances ~1.3 cells, then
anchors). Under the door pair the object has a second way to move that the
trap cannot veto:

> **Translation by conversion.** At the leading edge, ambient/own field is
> resonantly settled (intake clicks — the front cells' load grows); at the
> trailing edge, stimulated outtake un-settles matter into the field, which
> crosses the object at field speed and re-settles at the front. The `Em`
> distribution advances **without any dense hop carrying it** — the trap
> can keep its within-mode veto and the object still moves, because motion
> is the cycle, not displacement.

This is not a metaphor bolted on: it is the group velocity of the mixed
(lower-polariton) branch, §III.3.4, which is nonzero even when the bare
dense band is flat (trapped). The antagonism is *released*, not repeated,
**if and only if** the door traffic preserves phase — which is what the
coherent settlement is for, and what bar M4/M5 decides. Honest converse: if
M4 fails, the antagonism survives (§III.8, most-likely-failure).

---

# III.3 — Physics derivation

## 3.1 The pair from conservation

Postulates (fabric only):

- **P-cons:** `d/dt Σ_i (Em + Ee + Es + |χ|² + flight)(i) = 0` for the
  closed roster, at machine precision, with no compensating ledger.
- **P-grain:** every `Em`-changing conversion fires in whole atoms
  `ε = A₀ω/2π` at the source pitch, beat-gated (the adopted law).
- **P-prim:** all within-step energy motion is by norm-preserving rotation:
  per-cell `U(1)` (precession) or per-edge `SO(2)` Givens (hop), edges
  being spatial (F, U) or modal (D).
- **P-noscale:** every frequency in the update is built from `{w1, w2,
  q_detune, field_J, A₀, dt}` and roster geometry.

P-prim + P-cons force the door's continuous part to be a rotation in some
two-plane whose squared coordinates are energies. The only pair that (a)
does not touch `Em` continuously (P-grain) and (b) couples the field to the
matter-bound sector is `(ψ_e, χ)` with `|χ|²` a booked inventory —
uniqueness up to the choice of `χ`'s precession pitch, fixed at `Ω₂ = w2e`
by the requirement that the pair's detuning be the beat (so that the
existing beat gate is the door's own carrier, not a second clock).
Linearized per cell:

```
dψ_e/dt = −i·Ω₁·ψ_e − i·Σ_j K_ij·ψ_e(j) + i·g·χ
dχ/dt   = −i·Ω₂·χ + i·g·ψ_e                          g = door_g·w2
d|ψ_e|²/dt + d|χ|²/dt = 0 under the g-terms  (cross-terms cancel exactly)
```

with click settlement `|χ|² → Em` in atoms at beats (§2.5). Each g-term is
first-order and incomplete; the two are mutually inducing; conservation of
`Ee + |χ|²` under the coupling is an identity (the two cross-terms are
`+2g·Im(ψ_e*χ)` and `−2g·Im(ψ_e*χ)`).

## 3.2 The link current and its moments at the mode boundary

Within-mode (v93, unchanged): spatial link current
`J_s = 2·K_s·Im(ψ_i*·ψ_j)`; first moment `p = Σ_s J_s·d_s·û_s` = momentum
(center-of-energy theorem).

**Modal link current (new):**

```
J_D(i) = 2·g·Im(ψ_e*(i)·χ(i))     [energy per time crossing the mode boundary at i]
```

- **Zeroth moment** `Σ_i J_D(i)` = net instantaneous conversion rate = the
  metabolic flux `Φ = I − O` of METABOLIC_MATTER §1 — now a *measured
  algebraic quantity*, not a bookkeeping difference of two rate formulas.
- **First moment in space.** The modal link has zero spatial extent
  (`d = 0` within a cell), so the door itself carries no momentum — but the
  register field `χ(i)` across neighboring cells inherits the field's
  spatial phase gradient (`arg χ` tracks `arg ψ_e` through the coupling),
  and coherent settlement hands that gradient to `ψ_m`, whose spatial link
  current `J_s^{(m)} = 2·Im(ψ_m,i*·ψ_m,j)` is then immediately nonzero.
  **Momentum crosses the door because phase does.** Today's magnitude door
  erases `k` at conversion — which is the mechanistic shape of the standing
  v89 anomaly, the **~100× radiation-pressure deficit at the conversion
  door** (recorded, not gated). Pre-registered prediction **P-p2**: the
  door pair reduces that deficit; the p2 meter run with `door_g > 0` shows
  momentum delivery at the door significantly above the v89 floor. (Not a
  pass/fail bar for v94 — a scored prediction.)

## 3.3 Reduction to second order (the matter analogue of □ψ = 0)

Eliminate `χ` in the linearized pair (uniform background, load `x`,
Fourier on a regular patch of the roster):

```
[ω − Ω₁ − ε(k)]·[ω − Ω₂] = g²
```

the avoided-crossing (polariton) dispersion:

```
ω_±(k) = ½[Ω₁ + ε(k) + Ω₂] ± √( ¼[Ω₁ + ε(k) − Ω₂]² + g² )
```

Away from resonance the branches are field-like (`ω ≈ Ω₁ + ε(k)`, wave
equation at `c_sub`) and register-like (flat). Near the crossing
`Ω₁ + ε(k₀) = Ω₂` both quadrature pairs participate and the second-order
form is the **gapped wave equation of the mixed mode** — the matter
analogue of `□ψ = 0` with an index:

```
v_±(k) = dω_±/dk = v_bare(k) · W_e^{±}(k),      W_e = |⟨field|±,k⟩|² ∈ (0,1)
```

**Propagation speed, vacuum:** `x = 0`, detuning maximal, `W_e ≈ 1` on the
field branch — speed `c_sub`, unchanged (observation 1 preserved; the
vacuum correction is `O(g²/Δ²)`, measurable, small).
**Propagation speed, loaded matter:** detuning shrinks, `W_e < 1` — the
wave is part-field, part-register(-matter), slower by exactly its matter
fraction (observation 2, now derived from the door instead of absent). The
**metabolic wave** of bar M1 is the matter-weighted branch: a disturbance
in the cycle (extra load at one cell) propagates at

```
v_M ≈ c_sub · g² / [ Δ(x)²/4 + g² ]        (heavy-detuning limit)
```

— finite, `door_g`-tunable, load-dependent: a **pre-registered quantitative
scaling** (M1 measures `v_M(door_g)` and this curve is the claim). With
clicks included, `v_M` is duty-cycled by the settlement fraction per beat —
the correction is measurable and belongs to M1's meter, not to a new knob.

## 3.4 Vacuum and matter speeds, summarized

| regime | carrier | speed |
|---|---|---|
| vacuum field | `(fa1, fa2)` | `c_sub` (band slope; observation 1) |
| field in loaded matter | `(ψ_e, χ)` dressed | `c_sub·W_e < c_sub` (observation 2) |
| dense within-mode | `(m₁, m₂)` pass U | slow band; killed by trap (v93) |
| **metabolic cycle** | **`(ψ_e, χ)` + clicks** | **`v_M` — new, finite, M1's meter** |

---

# III.4 — Fundamental math and mechanics (DEC fragment)

## 4.1 The complex

Minimal fragment, as v93 §II.5 (0-forms, 1-forms, source terms; no 2-forms,
no discrete curvature — over-kit). The v94 addition is one product
structure:

- **Base:** the roster graph `G = (cells, links)` — dynamical, as always.
- **Fiber:** the three-point mode path `e — χ — m` per cell.
- **Total complex:** the prism `G × {e, χ, m}` with three edge families:
  spatial-e (pass F edges), spatial-m (pass U edges), and **fiber edges**
  `(i,e)—(i,χ)` (continuous, coupling `g`) and `(i,χ)—(i,m)` (quantized,
  click).
- **0-cochains:** `ψ_e, χ, ψ_m` — one complex value per prism node (with
  `|ψ_m|² = Em`, `arg ψ_m = th2` the pass-U identification). `Es` remains a
  scalar 0-cochain on cells (space does not join the spinor algebra; it
  moves only by click-ledger, as adopted).
- **1-cochains:** the currents — `J^{(e)}_s`, `J^{(m)}_s` on spatial edges,
  `J_D(i)` on e–χ fiber edges, and the click ledger (atoms) on χ–m fiber
  edges.

## 4.2 The update as edge rotations; conservation as a theorem

Every continuous operation in the step is one of:

1. per-node `U(1)`: `ψ ← e^{−iΩdt}ψ` — norm exact;
2. per-edge `SO(2)` Givens on a spatial edge (F hop, U hop) — two-node norm
   `|ψ_i|² + |ψ_j|²` exact;
3. per-edge `SO(2)` Givens on the e–χ fiber edge — node-local sum
   `Ee(i) + |χ(i)|²` exact;

and every discrete (click) operation is a ledgered move that was already
conservation-exact in the adopted law (condensation/evap/radiance
bookkeeping, `s_pull` splits), now including the χ–m solved-Givens whose
norm exactness is `Em + |χ|²` invariance (§2.5).

**Theorem (conservation).** Over one full step, `Σ_prism |ψ|² + Σ Es +
Σ flight` telescopes: each rotation contributes zero by (1)–(3); each click
contributes zero by its ledger identity; boundary terms vanish on the
closed roster (cell death flushes `|χ|²` by `field_inject` before roster
removal — §III.5). Hence drift = FP roundoff. This is v93 §II.4's
dissolution of conservation into the Givens theorem, extended to the fiber
— **the same proof, one more edge family.** (Stokes reading: the discrete
divergence of the three current families equals the node inventory rates;
summing, all internal edge contributions cancel in pairs.)

## 4.3 Anti-ignition for the new pair

v93's structural argument (random-phase bath self-gates) must be
re-established, not assumed. Four layers, three structural, one inherited:

1. **Detuning (structural, §2.4):** vacuum sits at maximal `|Δ| = w2 − w1 =
   1.25` with `g = door_g·w2`; the virtual admixture is bounded by
   `F(0) = 4g²/(Δ² + 4g²)` and is **reversible** — an unsettled register
   oscillates back at `Ω_R`. Nothing accumulates from coherent vacuum
   dressing alone.
2. **Phase averaging (structural):** a random-phase bath drives `χ`
   incoherently; the vector sum executes a bounded random walk with
   `⟨|χ|²⟩ ≤ F(0)·⟨Ee⟩` — no secular pumping (this is the v93 argument
   transplanted: the cross-term averages to zero over random phase).
3. **Settlement gates (inherited, unchanged):** workfn `econd` at the
   click, atom threshold `ε`, `f_conv` fraction, cap + evaporation above
   cap. A register that never settles is optics, not matter.
4. **Zero-throughput coupling:** the g-rotation is norm-preserving — it
   cannot amplify total throughput, only redistribute it (the L-1 lesson,
   here by construction rather than by renorm).

**The new risk, named honestly (III.8): coherent runaway.** A *phase-
coherent* lit bath over existing matter defeats layer 2: resonant pumping
with the settlement gate open is positive feedback (more load → smaller
`|Δ|` → larger `F` → more intake), bounded only by local `Ee` supply and
cap/evap. This is not obviously a bug — a bounded version of it is
*particle formation* — but it is exactly the configuration v93 left as the
unrun L1-C. It is pre-registered as bar M3b and must be run **early**.

---

# III.5 — Numerical simulation method (the kernel face)

## 5.1 Knob, state, defaults

- **Knob:** `door_g` (double, **default 0**; parse alongside the v93 knob
  family). `door_g = 0` ⇒ every v94 addition is skipped: no allocation
  touched in the hot loop, pass 6 takes the original branch — **byte-for-
  byte laws_V3a** (the 87-bar surface). One knob; no second scale.
- **New state:** `dc1[NC], dc2[NC]` (doubles, zero-initialized; born-zero
  on cell birth). Include `Σ(dc1² + dc2²)` in the conservation/drift meter
  and in the diag energy columns. FCS dump: add the register plane behind
  the same knob (orchestrator's format call; flag in the dump header).
- **Cell death / roster events:** before a cell leaves the roster, flush
  its register to the field: `field_inject(i, dc1² + dc2²)`, zero it. (The
  register is unsettled mode-transit; return-to-source is the flight-parcel
  rule α analogue.) Repack copies `dc1/dc2` with the cell.
- **Recommended campaign companions (not defaults):** `amp_nat = 2,
  amp_logate = 1, hop_order = 1` (the v93 working transport recipe) — the
  door pair is designed for the dense spinor to exist; see §5.4 for the
  `amp_nat = 0` degraded mode.

## 5.2 Pass D: placement and operations

Pass order becomes: … pass U → passes 3–5 → **pass D(τ/2)** → pass F →
**pass D(τ/2)** → pass 6 → … i.e. the door rotation is **Strang-
symmetrized around pass F** from day one. (RING_DNLS A.4 measured the cost
of getting Trotter order wrong — the vacuum retention ceiling was the
kernel's own sweep order. Do not repeat that discovery; symmetrize now.)

Each half-application, per cell `i`, `door_g > 0`:

```
/* (a) register precession, half-step: χ ← e^{−i·w2e[i]·dt/2}·χ  */
/* (b) fiber Givens, half-angle: θ = 0.5 · door_g · w2 · dt      */
fa1' = cc·fa1 + ss·dc2        dc1' = cc·dc1 + ss·fa2
fa2' = cc·fa2 − ss·dc1        dc2' = cc·dc2 − ss·fa1
       with (cc, ss) = (cos θ, sin θ)   — the pass-F/U Givens stencil,
       applied on the fiber edge (same-cell, cross-sector)
```

then recompute `Ee[i] = fa1² + fa2²` after the second half (the existing
recompute at freecell.c:2195 moves after D's second half). Notes for the
implementer: the stencil is `apply_uhop`'s (freecell.c:1265) with `(i, j)`
replaced by the two sectors of one cell — reuse the LUT trig discipline of
pass U; the rotation is unconditionally applied (no gate — linearize
lesson, FINAL_STATUS item 18: gates on rotation angles melt structure; the
*detuning* is the gate).

## 5.3 Pass 6 rewrite (the settlement), door_g > 0 branch only

Inside the existing `beat_fire` block, the **condensation branch is
replaced** (not duplicated — no double-intake) by register settlement:

```
if (door_g > 0) {
    chi2 = dc1² + dc2²
    econd_i = (existing workfn/e_cond text, verbatim)
    if (chi2 > econd_i) {
        d1 = f_conv · (chi2 − econd_i)
        d1 = atoms_fire(d1, eF, eD, &qcnvD[i])          /* unchanged credits */
        θ  = solve ΔEm(θ) = d1                           /* §2.5; ψ_m from (Em, th2) */
        if no θ reaches d1: clamp d1 to the largest whole-atom reachable; else no-fire
        apply Givens(ψ_m, χ; θ); Em, th2 ← |ψ_m|², arg ψ_m; dc ← rotated
        space book: dsp = s_pull·d1 from Es, exactly as the legacy branch
        ledger: door_in_total += d1 (+ existing cond_total, click diag, tags)
    }
    /* stimulated outtake: if C = 2·Im(ψ_m*·χ) < 0 and the solved branch
       gives ΔEm < 0 with |ΔEm| ≥ one atom on the OUTTAKE credits
       (eF2/eD2 swapped, atoms via qcnvF[i]): apply the negative-θ Givens;
       ledger door_out_total; backsplash bookkeeping as radiance/evap. */
}
/* evap-above-cap: unchanged. radiance: unchanged. rad_dual: leave as an
   apparatus knob, default 0 — the pair supersedes its role. */
```

The legacy `Ee`-condensation branch, `amp_door`, and `qp_phase` are
untouched at `door_g = 0` and bypassed at `door_g > 0` (one door, one
algebra — no stacking of phase heuristics on top of the coherent door).

## 5.4 Click-grain interaction, and the amp_nat = 0 degraded mode

The atom `ε = A₀ω/2π` machinery is used verbatim (`atoms_fire`,
`atoms_clamp`, both credit counters, `qatom_diag`). Between beats the
register holds and precesses; the beat samples it; the atoms quantize it —
**derivation of "the pair lives between clicks" is the schedule itself.**
At `amp_nat = 0` there is no dense spinor: settlement degrades to
magnitude-only (`Em += d1`, `th2` untouched, phase remains in `χ` and is
eventually returned fieldward) — legal for A/B tooling, **not** a campaign
configuration; document it as such.

## 5.5 C / Go mirror parity

- `door_g = 0`: **byte-identity** C ≡ Go (the standing bar).
- `door_g > 0`: tolerance-based abx (the v93 §IV.12 expectation at
  strength).
- **Standing debt folded in:** FINAL_STATUS item 5 flags that all
  `amp_nat > 0` campaign data are Go-only — C pass-U parity at strength is
  unverified. v94's abx surface must cover `(amp_nat > 0, door_g > 0)`
  jointly before any campaign conclusion is graded (bar M0 below).

---

# III.6 — Acceptance bars (pre-registered; the ratchet)

**Invariant surface:** the full battery, 87 V3a bars GREEN, at
`door_g = 0` (byte-inert), both kernels, before and after every kernel
edit. Bars below never soften; they sharpen by measurement.

- **M0 (parity):** C ≡ Go byte-identical at `door_g = 0`; tolerance-abx
  agreement on the M1 chain observables at `(amp_nat = 2, door_g > 0)` —
  closing the v93 Go-only debt. *Meter:* existing abx harness + `v_M`.
- **M1 (metabolic propagation):** on a quiet 1-D loaded chain (ring/chain
  seed at uniform `x₀ ≈ 0.3–0.6`), a localized field pulse at one end
  produces an `Em` disturbance arriving at distance `n` with finite speed
  `v_M`: threshold-crossing arrival times linear in `n` over ≥ 8 cells,
  ≥ 3 seeds. *Sharpener:* `v_M(door_g)` follows §3.3's scaling shape
  (monotone rising, saturating toward `c_sub·W`). *Decides:* the primitive
  propagates the cycle at all.
- **M2 (conservation):** `|drift_rel| ≤ 1e-13` with `door_g > 0` on live
  matter (M1 chain and a 3-D blob), `T = 1000+`, register energy included
  in the sum. *Decides:* the theorem survived implementation.
- **M3a (anti-ignition, dark):** dark V3a bath, `door_g > 0`, `T = 5000`:
  birth/glow/cond byte-stats within noise of `door_g = 0`; `Σ|χ|²` bounded
  (no secular growth). **M3b (armed coherent bath — the unrun L1-C, run
  EARLY):** phase-coherent lit bath over a seeded clump, `door_g > 0`:
  measure the growth curve; bar = no unbounded `Em` growth past cap defense
  (evap engages; `Em_tot` bounded); **record** the coherent-pumping
  threshold as a measured number, whatever it is.
- **M4 (door phase retention):** QUENCH-3-style winding on a lit-bath ring
  **with door traffic actually firing** (`cond > 0` verified — the A.5
  lesson: v93's refill arm never fired, so the baseline is honest only if
  clicks flow). Baseline `R2d ≈ 0.02`. **Target: `R2d ≥ 0.5` at `t = 80`
  with active traffic.** *Decides:* phase crosses the door.
- **M5 (the wall — antagonism release):** RING_DNLS B.3 rerun, condensing
  corner (`amp = 0.5`, law `q_detune`, `kx = 1.1`), `door_g > 0` vs `0`.
  Baseline: COM advances ~1.3 cells then anchors. **Bar: sustained COM
  advance ≥ 5 cells over `T = 80` with `Em_max` retention ≥ 50%** —
  translation *with* binding. *Decides:* translation-by-conversion is real.
- **M6 (metabolic fixed point):** quiet clump + lit bath at ambient `xb`,
  full V3a + `door_g > 0`, `T = 5000`: `Em` parks — late slope `|dEm/dt|`
  below the SLOW_REFILL slope (~2e-3/t.u.) by 10×, AND restoring: a ±20%
  `Em` perturbation at `t = 3000` decays back. Baselines: `rad_dual`
  SLOW_REFILL (climbing, no park), radiance-only SLOW_DRAIN. *Decides:* the
  cycle closes — `I(x*) = O(x*)`, `I′ < O′` — the METABOLIC_MATTER
  attractor, measured. *Score `x` on core cells, PR, `Em_max`, phase —
  never `Em_tot` alone (FINAL_STATUS language discipline).*

---

# III.7 — Hard-constraint compliance (point by point)

| constraint | how the design complies |
|---|---|
| No background / imported field / SI constants | `g = door_g·w2`; every frequency from `{w1, w2, q_detune}`; χ is transit bookkeeping on the existing roster |
| No state relabeling | `χ` changes transport: it closes a rotation pair, carries current `J_D`, and produces new dynamics (M1's `v_M` exists or the design dies) |
| No per-component jobs | one register, no flavor |
| **No unitarization of the click** | `Em` changes **only** at beats, in whole atoms, via `atoms_fire`/`atoms_clamp` verbatim; the continuous pair touches `Ee ↔ |χ|²` only |
| No sum-then-square | every energy move is a rotation (norms) or an atom ledger; coherent *addition* of amplitudes into an energy never occurs |
| No feedback / shadow / bias | no bias term anywhere; the g-rotation is norm-preserving (cannot amplify throughput — §4.3 layer 4); distinguish from v92: those were *biases on wants*; this is unitary algebra with booked energy |
| No phase written to `th2` at the door | phase is carried in `χ` (dedicated carrier); settlement touches `th2` only as the arg-readout of a norm-preserving Givens — the same operation class as pass U's own write-back, not a mixing formula. **Flagged for user judgment in III.8** since the letter of the rule could be read more strictly |
| No bar softening | §III.6 targets stated against measured baselines; the 87-bar surface is the invariant |

---

# III.8 — Flags, first experiments, failure modes

## What this completes vs what it opens

**Completes (v93's open wall):** the metabolic cycle as a first-class
transport channel (the open wall by name); the coherent cross-term `g`
(FINAL_STATUS item 10 — never coded — here as algebra, not a rate); the
door-phase problem (items 7/19's amp_door reworked from formula to
rotation); the untested lit-bath refill arm (RING_DNLS A.5); the standing
`rad_dual` no-park failure (M6 is its redesign with restoring balance); the
refraction gap (observation 2 had no mechanism); and it engages the v89
radiation-pressure deficit (P-p2, scored prediction).

**Opens (new risks):** (1) two new floats per cell — **explicit user
sign-off requested**: I argue §2.3 that `χ` is mode-transit bookkeeping
(cousin of flight), not a fourth mode — no new pitch, no transport of its
own, energy fully booked — but the judgment is the user's, per programme
rule. (2) Coherent-runaway (§4.3) — M3b runs early. (3) The `th2`-at-
settlement reading of the QUENCH-3 rule (III.7 table) — flagged, user
judges. (4) Vacuum dispersion shift `O(g²/Δ²)` — small, but the DS-bar
optics live on the field sector; any `door_g > 0` optics claim must re-run
the DS loci off-battery before being graded.

## The three first experiments (value order)

1. **E1 = M1 + M2 + M0** (one campaign): the chain. Cheapest, decides
   existence of the metabolic wave and the conservation theorem, and pays
   the C≡Go debt. *Outcome that decides:* arrival-time linearity → the
   primitive is real; no propagation → the design fails at the root (see
   failure 1).
2. **E2 = M5**: the antagonism rerun. One arm against a fully-measured v93
   baseline. *Outcome:* COM ≥ 5 cells with retention → the wall is
   released; stall → the trap vetoes conversion too (failure 2).
3. **E3 = M6 (+M3b first, for safety)**: the park. The programme's oldest
   open ledger — `I = O` at an isolated `x*`. *Outcome:* restoring park →
   METABOLIC_MATTER's attractor exists and v94 has its adoption case;
   climbing or draining → the coherent channel alone does not close the
   cycle and the failure mode localizes the missing structure.

## Most likely failure mode, honestly

**Settlement decoherence.** The pair will propagate in the `(ψ_e, χ)`
subspace (M1 likely passes — it is linear physics), but the *click* may
still scramble phase: the beat sampling is not commensurate with `Ω_R`, and
each solved-Givens composes `χ` against a `ψ_m` whose phase the bath has
meanwhile walked (RING_DNLS: contact noise dephases small supports in ~10
t.u.). Then M4 fails at the baseline, M5 stalls, and v94 has built a
correct optics of virtual conversion (refraction, dispersion — real
results) without metabolic matter — the wall moved from the door to the
settlement, one joint deeper. Mitigations inside this design: settle in the
register's rotating frame (compose phases at the beat instant, which the
solved-Givens does by construction), Strang placement, and the M4 lit-ring
geometry where `ψ_m` support never vanishes. If M4 still fails, the honest
conclusion is that phase retention needs the J-dependent link metric
(RECKONING's named next face) *in addition to* the door pair — the two
designs are compatible, and that composition is the registered fallback.

---

*End of DESIGN_claude.md. The orchestrator implements §III.5 exactly;
grok reviews against §III.7's table and the battery invariant; the first
runs are §III.8's E1–E3 in order, M3b before any lit-bath campaign.*
