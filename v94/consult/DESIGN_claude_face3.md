# v94 developmental design — face 3: the unification (matter as a self-bound state of the field via χ)

**Session ID:** `aebbc7eb-4716-4e30-8c22-484e115428f2` (resume with
`claude --resume aebbc7eb-4716-4e30-8c22-484e115428f2` when in-flight results land).
**Author:** claude (fable), developmental session, 2026-08-09.
**Deliverable for:** `v94/consult/BRIEF_face3.md` PART II (II.1–II.4).
**Grounded in:** `v94/consult/DESIGN_claude.md` (face 1, mine),
`v94/IMPLEMENTATION_STATUS.md`, `v94/consult/REVIEW_grok.md` (SOUND after the
ΔEm sign fix), `v94/kernel/freecell.c` (pass 1 :1401–1409, pass P :2170/:2249,
beat gate :2326–2329, workfn :2331–2333, radiance :2582–2600, settlement
:2333–2420), `v93/RING_DNLS.md` (B.1 condensation map), `v89/PRINCIPLE.md`,
`v94/laws/laws_V3a.cfg` (w1=1.65, w2=2.9, q_detune=1.2, cap=2.5, wf_on=1,
k_rad=0.05, p_rad=4).
**Status:** design only. No kernel code here. The one physics knob specified
(§III) is byte-inert at default 0, battery-gated.

---

## NEEDS (in-flight results, flagged as headlines)

- **NEEDS: face-2 stimulated outtake + M6 outcome** to grade bar **U1(b)**
  (the bound state proper is the metabolic park; intake alone cannot park —
  tentative answer given in §II.5, marked tentative).
  **[RESUME 1: face 2 LANDED — dependency half-resolved; U1(b) now awaits
  M6 only. See §X.]**
  **[RESUME 3: M6 first evidence — `chi_load = 1` slows refill 10×,
  approaching the park; mechanism = the beat throttle. Confirmation
  pending (T ≥ 2000 + perturbation). See §Z.1.]**
- **NEEDS: M4 outcome (door phase retention, `R2d ≥ 0.5` target)** to scope
  the "one algebra on every edge" claim (§I.6). Tentative: design assumes M4
  meets target; if M4 fails, the claim demotes to magnitude-sector-only and
  the registered fallback is the J-dependent link metric composition
  (face-1 §III.8).
  **[RESUME 2: MEASURED — NOT MET, and worse than baseline (the face-2
  outtake erodes the winding). (L1)'s phase half demoted as pre-registered;
  diagnosis, the frozen-phase ratchet, and the design-faithful fix
  (beat-integrated cross-term, `door_avg`) in §Y. `chi_load` does NOT
  rescue this — wrong sign; §Y.3.]**
- **NEEDS: M3b measured coherent-pumping threshold** to calibrate the fold
  prediction `Ee_c ≈ Δ₀²·cap/(4g²q(1+qx₀))` (§II.4; prefactor trusted only
  to O(1) — the *shape* in `door_g` and `x₀` is the claim).
  **[RESUME 1: MEASURED — knee `Ee ≈ 1.0` at `door_g=0.5` = 2.6× the
  vacuum prediction (O(1) as pre-registered); bounded, not runaway.
  Calibration + the surviving `1/g²` slope test in §X.2.]**
- **NEEDS: M1 `v_M(door_g)` curve** to fix the polariton weight `W_e` in the
  correspondence dictionary (§II.9). Tentative: face-1 §3.3 scaling assumed.
  **[RESUME 3: MEASURED — the metabolic wave does NOT propagate; the pair
  CONFINES (and cures a legacy runaway, 417 vs 17). `v_M` and the `W_e`
  dictionary entry are moot as specified. See §Z.2.]**

None of these blocks the face-3 kernel work (§III): the knob is inert at 0
and its first bars (U0, U2, U3, U4) are independent of the in-flight arms.

---

## Executive summary

The user's question, taken at its sharpest — *can `ψ_m` be retired, with
matter re-derived as a self-bound configuration of `χ`?* — has a definite
answer, and it is **no, by a theorem, not by habit** (§I.3): any substrate
that (H1) exhibits refraction (continuous, reversible, phase-coherent
partial conversion — the polariton mechanism face 1 restored), (H2) obeys
the atom grain (conversion in whole `ε = A₀ω/2π`, ħ-linearity measured to
~1e-8), and (H3) updates locally and Markovianly (no background, no
persistent reference) **must** carry a mode fiber with at least two edges —
one continuous (the virtual half) and one stroboscopic-quantized (the real
half) — and a persistent complex carrier at the joint between them, distinct
from the settled inventory. That carrier is `χ`; the settled inventory is
`ψ_m`. The three-node fiber `e — χ — m` of face 1 is not defensive
scaffolding; it is the minimal structure compatible with the measured
physics. Retiring `ψ_m` is also, independently, the wrong side of the
programme boundary: matter-as-a-profile-in-a-field is the rejected pre-v89
road (§I.5); the v89 principle makes the mode *boundary* — conversion — the
ontological content of matter, and erasing it erases the principle's
referent.

What the unification **does** deliver — and it is substantial — is unity of
**law**, not erasure of **structure** (§I.6): one fiber field `Ψ = (ψ_e, χ,
ψ_m)` per cell under one algebra (norm-preserving rotations on every edge),
one pitch functional `Ω_a = w_a/(1+q_detune·x)`, and — the face-3 physics
knob — **one load**: the v89 law "flight is load" extended to door-flight,
`x = (Em + flload + |χ|²)/cap` behind byte-inert `chi_load` (default 0).
That single booking line closes the last asymmetry between space-transit and
mode-transit, makes the register self-detuning with **zero new constants**
(the brief's candidate `Ω_χ = w_χ/(1+q_χ|χ|²/cap)` realized with `w_χ = w2`,
`q_χ = q_detune`, derived rather than posited), and yields derived,
falsifiable physics: the Kerr effect (§II.4), a closed-form coherent-pumping
fold `Ee_c` (the substrate's pair-production threshold, §II.4), a derived
beat-slowing throttle that joins the anti-ignition stack as layer 5
(§II.8), and a park-stability inequality showing the metabolic fixed point
survives the Kerr feedback by ~30× margin at law constants (§II.5).

On the click grain (§II.6): face 3 shows the door is **more derived than
face 1 assumed**. The click's *timing* is a topological event — the kernel's
`cbeta` accumulator is literally a winding counter of the mode-fiber
relative phase, and the door fires on winding — and the grain's
*integer-ness* is winding topology (integers for free, as `W = +2` exact on
the C6 ring already demonstrated for spatial winding). What remains imposed
is exactly two things: the action scale `A₀` per winding (the substrate's
ħ — its one quantum postulate) and the two-edge fiber structure (proven
necessary). The door stays special, and its specialness now has three named
irreducibles instead of being one opaque mechanism.

---

# PART I — Theory: the decisive test, run honestly (II.1)

## I.1 The inventory and the reduction candidates

Per cell today (face 1, `door_g > 0`): `Es` (space, scalar, click-ledgered),
`ψ_e = fa1 + i·fa2` (field), `χ = dc1 + i·dc2` (register / polarization),
`ψ_m = √Em·e^{i·th2}` (settled matter). Three complex carriers + one scalar.
The candidate reductions, in decreasing ambition:

- **(R1) One amplitude:** `ψ` per cell; matter = a bound configuration of
  `ψ`. Field and matter are envelope/spectral components of one carrier.
- **(R2) Retire `ψ_m` only:** matter = the self-bound sector of `χ`; the
  fiber becomes `e — χ` with the grain at... somewhere (that is the
  problem).
- **(R3) Retire nothing; unify the laws:** one algebra, one pitch
  functional, one load booking; the carriers remain three.

The brief's sharp form is (R2) (with (R1) behind it). §I.2 steelmans them;
§I.3 proves what kills them; §I.6 develops (R3), which survives and has
content.

## I.2 The steelman: matter as self-bound χ

Give `χ` the brief's self-detuning, realized with no new constants by
booking `|χ|²` into the one load (§III.1): `Ω_χ(x) = w2/(1+q_detune·x)`,
`x ⊇ |χ|²/cap`. Then:

- A locally intense `χ` detunes its own cell. The within-mode analogue is
  *measured*: RING_DNLS B.1 — `q_detune = 0` melts (PR 456→853); any
  `qd ≥ 0.6` condenses (peak ×3.5, envelope-frozen at deep `qd`). The same
  DNLS mechanism transplants to the register: self-trapping exists as
  physics in this substrate. So "bound χ" is *dynamically* constructible.
- The bound component and the driven (virtual) component are spectrally
  distinct: the virtual part is slaved to the field (`χ_v ≈ g·ψ_e/Δ`,
  oscillating at the field's pitch); the bound part is free (oscillating at
  its own pitch `Ω_χ`). One complex number can carry the superposition.
- "Real matter" would then be the self-trapped spectral component; "virtual
  dressing" the slaved component; the virtual/real distinction becomes
  dynamical (a bifurcation) instead of structural (a second carrier).

This is the strongest available form of the user's hypothesis, and it gets
remarkably far. It fails in exactly three places, and each failure is
informative — they become the three irreducibles of §II.6.

**Failure 1 — a bifurcation is a threshold, not a ladder.** The
self-trapping transition is sharp (that is what B.1 measured) but it
quantizes *existence*, not *transfer*: above threshold a continuous family
of bound norms exists. The measured door physics is a **ladder** — whole
atoms `ε = A₀ω/2π`, ħ-linearity to ~1e-8, the two-atom credit. No classical
bifurcation of a smooth flow produces an exact integer ladder in the
transferred quantity; ladders come from counting structures (topology or
imposed integrate-and-fire). Binding geometry alone therefore cannot
reproduce the atom bars.

**Failure 2 — the spectral split is not locally decidable.** To *act* on
the bound component only (settle it, cap-test it, radiate it) the update
must project it out. A frequency projection requires temporal depth — a
filter, a rotating-frame reference, a history register. Any such structure
is persistent state that is merely re-valued: background-adjacent, and in
practice at least one more real pair per cell — i.e. `ψ_m` rebuilt under
another name. Under (H3) locality/Markov, the split must be *structural* if
anything is to read it.

**Failure 3 — the gates lose their home.** Face-1 anti-ignition layer 3 is
the settlement gate stack (workfn `econd`, atoms, `f_conv`, cap+evap) — all
predicated on the register/inventory distinction ("a register that never
settles is optics, not matter"). In (R2) there is no settlement event to
gate: energy in `χ` *is* matter the moment it self-traps. The adopted
workfn (presence rule `Em ≥ wf_floor`, freecell.c:2332) is not expressible:
there is no `Em`. Every gate would have to be rebuilt as a condition on a
spectral decomposition that (Failure 2) is not locally available.

(R1) inherits all three failures plus one more: with field and matter one
carrier, the load-detuning `q_detune` acts on the *field's own* pitch
through its *own* intensity even in vacuum flight — vacuum optics become
generically nonlinear, contradicting the measured linear vacuum optics
(v89: "vacuum optics linear, Kerr = loaded matter").

## I.3 The fiber-minimality theorem

The three failures generalize into a structure theorem. Stated carefully,
with its hypotheses exposed for attack:

**Hypotheses.**
- **(H1) Refraction realism.** Loaded matter slows the field by continuous,
  reversible, phase-coherent partial conversion (the Lorentz-oscillator /
  polariton mechanism; face-1 §1.3 observation 2 — the mechanism face 1
  supplied and the DS optics of the field sector demand in matter).
- **(H2) Grain realism.** The settled inventory changes only in whole atoms
  `ε = A₀ω/2π`, beat-timed (measured: ħ-linearity ~1e-8, two-atom credit;
  the battery's atom bars).
- **(H3) Locality/Markov.** The update is one-step, per-cell/per-edge; no
  persistent auxiliary reference (no background; no temporal filters).

**Claim.** Any dynamics satisfying (H1)–(H3) contains: (a) at least two
mode-fiber edges — one carrying continuous reversible transfer and **no**
inventory change, one carrying **all** inventory change, stroboscopically
and in atoms; and (b) a persistent two-real-component carrier at the node
between them, distinct from the settled inventory.

**Proof sketch.** (a): a single edge would have to carry transfer that is
both continuous-in-time (H1: the dressing exists between clicks — Rabi
period `2π/Ω_R ≫ dt`) and exactly integer-summed into the inventory (H2).
Continuous transfer with integer inventory increments means the inventory
is constant between clicks while the *other* end changes continuously — so
the edge's two ends are not a conserved pair, contradicting the unitary
primitive (P-prim); alternatively transfer ≡ 0 between clicks, contradicting
(H1). Hence ≥ 2 edges with the roles as stated. (b): the mid node holds the
in-transit energy *and its phase* across many steps and beats (H1 dressing
persists; M4's claim is that phase survives to settlement), hence persistent
≥ 2 real components; it cannot be identified with the inventory, else (H2)
applies to it and kills (H1) again — the argument of (a) one node over. ∎

**Corollary (the decisive test answered).** `ψ_m` cannot be retired. The
escapes are exactly: drop (H2) — an a-quantal substrate, contradicting the
battery; drop (H3) — temporal filters, i.e. background-adjacent memory;
or quantize the substrate's symplectic structure (second quantization) — a
different programme. **Refutation condition, pre-registered:** exhibit a
local Markovian one-carrier-per-mode update that passes both a refraction
bar (index physics in loaded matter) and the ħ-linearity bar. If such a
design exists, this theorem's (a)-step is wrong and face 3 must be
reopened.

## I.4 The three rejections of face-1 §2.6, revisited one by one

The brief asks whether they are permanent or open work. Verdict:

1. **"It would unitarize the door."** **Permanent, now derived.** Face 1
   asserted it defensively; §I.3(a) proves the continuous edge and the
   grained edge cannot be one edge. The door's specialness is reduced to
   three named irreducibles (§II.6) — smaller than face 1 assumed, but not
   zero.
2. **"Vacuum would have no virtual/real distinction (ignition)."**
   **Permanent, now derived.** The distinction must be structural (Failure
   2): virtual = `χ` (reversible, unsettled, gate-invisible), real = `ψ_m`
   (settled, gated, capped). The adopted workfn presence rule makes vacuum
   settlement *impossible* (`econd = wf_far = 99 ≫ cap` for `Em <
   wf_floor`) — a gate that has no expressible form without the
   register/inventory split (Failure 3).
3. **"It would re-import the binding/mobility antagonism."** **Resolved
   structurally — this one was open work, and face 1 already did it.** The
   field couples to the trapped hoard only through `χ` with bare-pitch `g`;
   the trap keeps its within-mode veto; translation-by-conversion (face-1
   §2.7) bypasses it. M5 measures whether the release is real. Face 3 adds
   the quantitative dressing physics (§II.4) but no new structure here.

## I.5 The programme boundary: which side is "one amplitude" on?

The brief flags it exactly: state reduction to one amplitude is the boundary
between "unification" and "the rejected pre-v89 multi-component lattice on a
background." The argument that (R1)/(R2) fall on the **rejected** side:

- Pre-v89's constitutive failure was matter-as-profile: shapes *in* a field
  that persists and is re-valued — the q-ball in its background field.
  "Matter = a bound configuration of the one field ψ" is that architecture
  verbatim, with the roster substituted for the lattice. The roster's
  dynamism does not save it: the *field* becomes the thing that persists
  everywhere and is merely re-valued, i.e. the background, and matter
  becomes bookkeeping on its profile.
- The v89 principle — "matter is converted space; energy only changes
  mode" — locates matter's identity in the **mode boundary**, i.e. in
  conversion events. Erase the boundary (make everything field
  rearrangement) and "energy changes mode" has no referent; the one law
  loses its content. The mode boundary *is* the principle, not an
  implementation detail.
- What the principle licenses — demands, even — is that the *same law*
  govern both sides of the boundary and the crossing itself. That is (R3),
  and it is what faces 1–3 jointly build.

So the unified picture's answer to "what plays the role of `ψ_e, ψ_m, χ,
Es`" (brief II.1): **not** one amplitude, **not** an amplitude plus its
conjugate (the conjugate carries no independent freedom; smuggling freedom
into it doubles the symplectic structure — a second carrier covertly), but
**one section `Ψ(i) = (ψ_e, χ, ψ_m)(i)` of a rank-3 mode fiber over the
roster**, plus the base scalar `Es`. One object, one algebra, three fiber
components — the mathematically honest reading of "matter and field are one
thing" that does not erase conversion.

## I.6 The unification face 3 delivers: one law on one fiber field

**(L1) One algebra.** Every continuous energy motion is a norm-preserving
rotation: per-node `U(1)` (precession), per-edge `SO(2)` Givens — spatial-e
(pass F), spatial-m (pass U), fiber e–χ (pass P). Every discrete motion is
a ledgered atom click (fiber χ–m settlement; radiance; evap). Faces 1–2
built this; face 3 states it as law and adds nothing to it. *Scope caveat:*
"one algebra including phase" is graded by M4 (NEEDS flag); magnitude-sector
unity is already measured (M0/M2/M3 green).

**(L2) One pitch functional.** `Ω_a(i) = w_a/(1 + q_detune·x_i)` for
`a ∈ {1,2}`, applied to `ψ_e` (w1), `ψ_m` (w2), and `χ` (w2 — the register
is matter-pitched, already implemented, grok §1.7). No carrier has a private
clock law. (The field/dense `e^{∓i}` sign split is a fixed convention of the
substrate, not a second law; the idealized equations below use one sign.)

**(L3) One load — the face-3 knob.** v89: "load flattens pitch and flight
is load: `x = (Em + flload)/cap`." Space-transit energy (`flload`) loads the
pitch; mode-transit energy (`|χ|²`) today does not — the last asymmetry
between the two transit ledgers. Face 3 completes the law:

```
x_i = (Em_i + flload_i + |χ_i|²) / cap        (chi_load = 1)
```

This is not a new law; it is the existing law applied to the new transit
ledger. Consequences derived in Part II: register self-detuning with zero
new constants, the Kerr effect, the pumping fold, the beat throttle.

**(L4) One grain, minimally imposed.** Conversion into or out of the
settled inventory fires only on windings of the mode-fiber relative phase
(topological timing), in whole source-cycles each carrying action `A₀`
(imposed scale). §II.6.

## I.7 Postulates of the unified substrate (restated)

- **P-cons:** `d/dt Σ_i (Em + Ee + Es + |χ|²)(i) + flight + Eh = 0` on the
  closed roster at machine precision, no compensating ledger.
- **P-prim:** all continuous energy motion is norm-preserving rotation on a
  node or an edge of the prism complex (§II.2); all discrete motion is an
  atom-ledgered click.
- **P-grain:** inventory-changing conversion fires only at windings of the
  fiber relative phase, in whole atoms `ε = A₀ω/2π` at the source pitch,
  two-atom credit.
- **P-pitch:** one pitch functional `w_a/(1+q_detune·x)`; **one load `x`
  booking holdings + all transit** (space and mode) against `cap`.
- **P-noscale:** every rate in the update is built from `{w1, w2, q_detune,
  field_J, door_g, A₀, cap, dt}` + roster geometry. Face 3 adds **no
  constant** (`chi_load` is a 0/1 booking switch, not a scale).

---

# PART II — Fully derived symbolic equations (II.2)

## II.1 The unified equations of motion

Continuous part (idealized single-sign form; kernel conventions per grok
§1.7 — field `e^{−i}`, dense/register `e^{+i}` — are a fixed isomorphism,
not different physics). Per cell `i`, with `x_i` per (L3), `Ω₁ = w1/(1+q·x)`,
`Ω₂ = w2/(1+q·x)`, `g = door_g·w2` (bare — the door stiffness is the
fabric's, not the load's; face-1 §2.4):

```
i·dψ_e/dt = Ω₁(x_i)·ψ_e  −  Σ_j K_ij·ψ_e(j)  −  g·χ_i
i·dχ_i/dt = Ω₂(x_i)·χ_i  −  g·ψ_e(i)                       (no spatial hops: χ is transit)
i·dψ_m/dt = Ω₂(x_i)·ψ_m  −  Σ_j K^m_ij·ψ_m(j)              (pass U; between clicks)
Es_i      : click-ledgered only (s_pull splits at settlements; unchanged)
```

Discrete part, at beat windings (`|∫(Ω₁−Ω₂)dt| = 2π`): the solved χ–m
Givens moving `d1 = f_conv·(|χ|² − econd_i)` re-quantized to whole atoms
(face-1 §2.5 with grok's sign fix: `ΔEm(θ) = sin²θ·(|χ|²−Em) − sinθcosθ·C`,
`C = 2·Im(ψ_m*χ)`); the face-2 negative branch (stimulated outtake) on
outtake credits; radiance `k_rad·cap·(Em/cap)^{p_rad}·(w2/w2e)` in atoms;
evap above cap. All unchanged by face 3 except through `x`.

The nonlinearity structure: with `chi_load = 1` the system is a **DNLS on
the prism** — the pitches depend on the local norms through the one load.
The three derived regimes: vacuum (linear, `x ≈ 0` — vacuum optics stay
linear, preserving the v89 measurement, because `χ ≈ 0` there); loaded
matter (self-trapped `ψ_m`, measured in RING_DNLS; dressed field, §II.4);
and the boundary (the fold, §II.4).

## II.2 The DEC fragment (unchanged complex, one new 0-form definition)

The face-1 prism complex stands verbatim (base roster graph × fiber path
`e — χ — m`; 0-cochains `ψ_e, χ, ψ_m` on prism nodes, `Es` on cells;
1-cochains = the three current families: spatial `J^e_s, J^m_s`, fiber
`J_D(i) = 2g·Im(ψ_e*χ)`, and the atom ledger on χ–m edges). Face 3 adds
**no cell, no edge, no cochain**. It redefines one derived 0-form: the load

```
x = ⟨inventory + incident transit⟩/cap :  cells → ℝ⁺
    inventory = Em ;  transit = flload (space-flight 1-forms incident on i,
    pulled back to the cell) + |χ|² (mode-flight 0-form at i)
```

and pipes it through the one pitch functional. The aesthetic point is also
the technical point: `flload` is the cell-restriction of edge-transit;
`|χ|²` is the cell-restriction of fiber-transit; the load is the total
transit-plus-holdings density at the node. `chi_load = 1` makes the
definition uniform over the two transit families.

## II.3 Self-binding I: the within-mode bound state (existence + stability)

The settled sector's binding is the dense DNLS, already measured; here is
the derivation the brief asks for, with the existence region.

Linearize pass U about a lump: cell at load `x`, neighbor at `x′`. The
pitch gap is

```
δΩ(x, x′) = Ω₂(x) − Ω₂(x′) ≈ −w2·q·(x − x′)/(1 + q·x̄)²
```

Self-trapping (resonant decoupling) when the gap exceeds the hop bandwidth
`B = z·J̄_m` (z = coordination, `J̄_m` the mean dense hop rate): the lump's
frequency leaves the band of its surroundings and linear transport off the
lump is suppressed as `(J̄_m/δΩ)²`. Existence region for a lump of contrast
`Δx = x − x′`:

```
w2·q·Δx / (1 + q·x̄)²  ≥  z·J̄_m            (self-trap criterion)
⇒  Δx_c(q) ≈ z·J̄_m·(1+q·x̄)²/(w2·q)        — the 1/q hyperbola at small q
```

**Calibration against measurement:** RING_DNLS B.1 (amp = 0.5, `amp_nat=2`)
measured the threshold in `q` between 0 and 0.6 with condensation at all
`qd ≥ 0.6` and depth monotone in `q` — the hyperbola crossed from the `q`
side at fixed `Δx`. This pins `z·J̄_m/w2` empirically; the **U1(a) sharpener**
(§III.6) is the two-parameter check: the threshold locus in the `(q, Δx)`
plane follows `q·Δx/(1+q·x̄)² = const`.

Stability of the trapped lump: the trap deepens with contrast
(`d(δΩ)/d(Δx) > 0` at fixed `x̄`) — the rich-get-richer mechanism RING_DNLS
A.4/B.1 measured — so the bound state is an attractor of the within-mode
dynamics up to cap (where evap clips it) and down to the melt threshold.
What the trap does *not* provide is persistence against the door (radiance
drains it: v91, `t_half` 80–140) — persistence is the metabolic fixed
point, §II.5.

## II.4 Self-binding II: the register under `chi_load` — Kerr, and the fold

The register's own nonlinearity, with zero new constants. Adiabatic
(slaved) response of `χ` to a quasi-static local field (`g ≲ |Δ|`; at
`door_g = 0.3`, `g = 0.87` vs `Δ₀ = 1.25` — marginal, so prefactors below
are O(1)-trusted only):

```
χ_v ≈ g·ψ_e/Δ(x)   ⇒   |χ|² ≈ g²·Ee/Δ(x)²,    Δ(x) = Δ₀/(1+q·x),  Δ₀ = w2−w1
```

**Kerr, derived.** With `chi_load = 1`, the dressing loads the pitch:
`x = x₀ + |χ|²/cap` ⇒ the local index (both `Ω₁` and the polariton weight)
depends on `Ee` — intensity-dependent optics **in loaded matter only**
(in vacuum `wf_far` blocks settlement and `x₀ ≈ 0` keeps the dressing both
small and virtual — and a small `|χ|²` in vacuum contributes `O(g²Ee/Δ₀²
cap)` to `x`: the measured-linear vacuum optics acquire a correction that
must be *measured* small, bar U5's DS arm, not assumed small). This is the
v89 line "Kerr = loaded matter" produced by the door instead of postulated.

**The fold (pair-production threshold), derived.** Self-consistency of the
slaved solution with `u = |χ|²/cap`, `s = 1 + q·x₀`, `β = g²·Ee/(Δ₀²·cap)`:

```
u = β·(s + q·u)²
```

Two branches merge at the fold `∂u(RHS) = 1`: `2βq(s+qu) = 1`, giving

```
β_c = 1/(4·q·s)      ⇒      Ee_c = Δ₀²·cap / (4·g²·q·(1+q·x₀))
```

Above `Ee_c` the slaved (perturbative, reversible) dressing has no bounded
solution — the register runs to the nonperturbative regime, bounded only by
the Rabi ceiling (`|χ|² ≤` local field supply; norm-preservation) and, on
the settled side, by cap+evap. **Numbers at law constants** (`Δ₀ = 1.25`,
`cap = 2.5`, `q = 1.2`):

| door_g | g | Ee_c (vacuum, x₀=0) | Ee_c (x₀=0.325 chain) |
|---|---|---|---|
| 0.15 | 0.435 | **4.30** | 3.09 |
| 0.3 | 0.87 | **1.08** | 0.77 |
| 0.5 | 1.45 | 0.39 | 0.28 |

Pre-registered claims (U3 sharpener): the fold exists; its location scales
as `1/g²` in `door_g` and falls with ambient load `1/(1+qx₀)`; the measured
M3b threshold calibrates the O(1) prefactor. **Interpretation, stated
plainly:** the fold is not an ignition bug — it is the substrate's
**stimulated pair-production threshold**, the intensity at which coherent
field over matter dresses nonperturbatively and (where the workfn presence
rule permits, i.e. over existing matter only) converts. At the adopted law
the vacuum cannot settle at any intensity (`econd = wf_far = 99 > cap ≥
|χ|²_max`) — genesis requires a matter seed, which is the correct physics
(stimulated conversion on existing matter) and the structural anti-ignition
backstop.

## II.5 The bound state proper: the metabolic soliton (tentative — NEEDS face 2)

**[RESUME 3 (§Z): M6 first evidence — `chi_load = 1` slows the refill
slope 10× (SLOW_REFILL → near-park). The measured restoring mechanism is
the §II.8 layer-5 beat throttle on the dressing-booked load, not the
radiance wall this section leans on — see §Z.1 for the calibration and
the corrected mechanism attribution.]**

Matter that *persists* is the composition of the two bindings with the
cycle closed: a self-trapped `Em` core (§II.3), dressed by its own register
(§II.4), parked at the metabolic fixed point `I(x*) = O(x*)` with
`I′ < O′` (restoring). Face 2's outtake is the missing half; intake-only
face 1 provably cannot park (its `I > 0` always at `chi2 > econd`, and its
only outtake is radiance — the measured v93 SLOW_REFILL/SLOW_DRAIN pair,
no restoring crossing). **Tentative** pending M6, here is the stability
inequality with the face-3 Kerr term included, so the bar is pre-derived:

Rate-level idealization: `I(x) ≈ ν_beat(x)·E_click(x)` with `ν_beat =
|Δ(x)|/2π = Δ₀/(2π(1+qx))` and per-click settlement rising with the
transfer fraction `F(x) = 4g²/(Δ(x)²+4g²)`; `O(x) = k_rad·cap·x_h^4·(1+qx)`
(radiance on holdings `x_h = Em/cap` with the clock-compensation factor,
freecell.c:2584) + evap above cap. Using the derived identity

```
dF/dx = 2q·F·(1−F)/(1+qx)
```

the park is restoring iff

```
O′/O − I′/I > 0 :   4/x*  +  q/(1+qx*)   >   2q(1−F(x*))/(1+qx*)  −  q/(1+qx*)
```

(LHS: radiance steepness + comp factor; RHS: Kerr door-opening − beat
throttle). At `x* = 0.62`, `q = 1.2`, `door_g = 0.3` (`F(x*) = 0.85`):
LHS ≈ 6.45 + 0.69 = 7.14; RHS ≈ 0.20 − 0.69 < 0. **The park survives the
Kerr feedback with ~30× margin, and the margin *grows* with `door_g`**
(larger `F` ⇒ smaller `1−F`). The Kerr term cannot destabilize the
metabolic fixed point at law constants for any `door_g`; the genuine risk
is upstream — the fold (§II.4) sitting inside campaign field intensities —
not the park slope. This inequality is the **U1(b) sharpener**.

## II.6 The click grain: what emerges, what remains (U4 answered)

**Derived — the timing is topology.** The kernel's beat gate is, read
exactly, a winding counter: `cbeta_i = ∫(w1e − w2e)dt` is the mode-fiber
relative phase; the door may fire only when it completes a full winding
(`|cbeta| = 2π`, freecell.c:2326–2329). Face 3 elevates this from mechanism
to law (P-grain): **conversion fires on windings of the fiber phase.**
Winding numbers are integers by topology — the same topology that held
`W = +2` exact on the C6 ring for 200 t.u. So the click's *when* and the
grain's *integer-ness* are derived, not imposed: a partial winding is not
an event, exactly as a partial loop is not a winding.

**Derived — the throttle.** The winding rate `ν_beat = Δ₀/(2π(1+qx))`
falls with load: heavily loaded (or, under `chi_load`, heavily dressed)
matter clicks *less often* per unit time. Anti-ignition layer 5, free of
charge (§II.8) — and independently measurable (click-rate vs `x`, a U4 arm).

**Imposed — the scale.** What one winding *carries*: whole source-cycles,
`A₀` of action each, `ε = A₀ω/2π` at the source pitch, two-atom credit.
No classical continuous flow generates that exact ladder (§I.2 Failure 1);
in continuum physics it comes from canonical quantization, which is not
this substrate. **`A₀` is the substrate's ħ — its one quantum postulate —
and it is a permanent residue.** Deriving it would mean second-quantizing
the fabric: a different programme, explicitly not proposed.

**Imposed-but-proven-necessary — the two-edge structure.** §I.3. The
integrate-and-fire at the χ–m edge is the minimal counting structure
compatible with (H1)+(H2)+(H3).

Scorecard for U4: timing **derived** (topology); integer-ness **derived**
(topology); throttle **derived** (pitch law); scale `A₀` **imposed**
(the ħ postulate); two-edge fiber **imposed and proven minimal**.

## II.7 The conservation theorem, restated (with a stronger corollary)

Face-1 statement (grok-verified §2.5): over one step,

```
Σ_i (Em + Ee + Es + |χ|²)(i) + Σ_slots flight + Eh
```

changes only by FP roundoff — each continuous operation is a node-`U(1)` or
edge-`SO(2)` rotation (norm-exact per (1)–(3) of face-1 §4.2), each click a
ledger identity, boundary terms vanish (χ death-flush + BH χ-flush; the
former is a §III.2 prerequisite).

**Face-3 corollary (pitch-independence).** Every rotation is norm-preserving
for *any* rotation angle; the pitch functional enters only through angles.
Therefore the theorem holds **uniformly over all pitch laws** — in
particular `chi_load` (which touches only `x`, which touches only pitches
and gates) **cannot move the drift floor, by construction**. U2 verifies
the implementation, but a drift regression under `chi_load = 1` would
indicate a coding error, never a physics choice. This is the strongest form
of "conservation as a theorem of the update, never a patched ledger": the
theorem quantifies over the entire pitch sector of the design space.

## II.8 Anti-ignition, re-derived under the unified law (U3)

The brief demands the re-derivation "without a separate matter sector";
§I.3 concluded the sector stays, so the obligation becomes: re-derive the
stack under the unified law (L1–L4), where the register self-detunes and
every gate reads the one load. The stack, now five layers — four inherited
and re-checked, one new and derived:

1. **Vacuum detuning (structural).** `|Δ(0)| = Δ₀ = 1.25` maximal; `g`
   bare; virtual admixture `F(0)·Ee` reversible. Unchanged by `chi_load`
   in vacuum (`χ ≈ 0` ⇒ `x` unchanged ⇒ byte-identical dark behavior —
   and this is now *structurally* checkable: U0's dark-bath arm at
   `chi_load = 1` must be byte-identical, since every new term multiplies
   `|χ|² = 0`).
2. **Phase averaging (structural).** Incoherent drive ⇒ `χ` executes a
   bounded random walk, `⟨|χ|²⟩ ≈ F·⟨Ee⟩`, no secular term. `chi_load`
   feeds `⟨|χ|²⟩` into `x` — a small *mean* shift, not a pump; bounded
   because layer 4 bounds `|χ|²` by local supply.
3. **Settlement gates (inherited, sharpened).** Workfn presence rule:
   vacuum (`Em < wf_floor`) demands `wf_far = 99 > cap` — **vacuum
   settlement is impossible at any field intensity**; genesis requires a
   matter seed. Atoms, `f_conv`, cap+evap unchanged. Note `econd` reads
   `Em`, not `x` — the gate is `chi_load`-blind by construction (checked
   against freecell.c:2331–2333).
4. **Norm preservation (structural).** The g-rotation redistributes, never
   amplifies; `|χ|²` ≤ local field supply always (the L-1 lesson as
   algebra).
5. **Beat throttle (NEW, derived).** `ν_beat ∝ 1/(1+qx)`: load — including
   dressing, under `chi_load` — *slows the click clock*. The two coherent
   feedbacks at the door (`F` rises with `x`; `ν_beat` falls with `x`)
   compete; per-time intake `∝ ν_beat·F` is **non-monotone** in `x` when
   `2g > Δ(x)` (already true at `door_g = 0.3` in vacuum: `2g = 1.74 >
   1.25`), i.e. deep dressing throttles its own click rate. The net sign
   is regime-dependent; U3's sweep measures it. (Radiance is exempt from
   the throttle by its comp factor — an *asymmetry favoring outtake* at
   high load, which is anti-ignition-correct and pre-dates face 3.)

The one genuine coherent hazard is the fold (§II.4): above `Ee_c`, dressing
over *existing* matter goes nonperturbative and intake feedback engages
until cap defense. That is bounded creation on a seed — candidate physics,
not runaway — but it must be *measured* before any lit-bath campaign at
`chi_load = 1` (bar U3, extending face-1 M3b exactly as face 1 extended
L1-C).

## II.9 The correspondence dictionary (U5's contract)

| unified object (theory) | kernel observable (v94 face 1) |
|---|---|
| fiber field, e-component `ψ_e` | `fa1, fa2`; `Ee` |
| fiber field, χ-component (polarization / door-flight) | `dc1, dc2`; `chi2_tot` diag |
| fiber field, m-component (settled inventory) | `Em, th2` |
| the one load `x` (holdings + all transit) | pass-1 `x` (chi_load=1 adds `|χ|²`) |
| fiber current `J_D = 2g·Im(ψ_e*χ)` | metabolic flux meter; rate of `door_in` (− `door_out`, face 2) |
| beat winding `∮dφ_D` | `cbeta` wraps; click count |
| atom ladder | `qcnvD/qcnvF` counters; ħ-linearity bars |
| polariton weight `W_e(k)` | `v_M/c_sub` from M1; loaded-arm DS index shift |
| self-trap threshold (within-mode) | RING_DNLS B.1 locus; U1(a) `(q, Δx)` sweep |
| pair-production fold `Ee_c` | M3b/U3 measured dressing-jump threshold |
| Kerr coefficient | loaded-arm index shift vs `Ee` (U5 DS arm) |
| metabolic bound state | M6 park (`x*`, `Em*`, restoring slope) + `chi2` profile locked to `Em` profile |
| park stability margin (§II.5 inequality) | M6 perturbation-decay rate |

**Must-reproduce set (U5 pass condition):** at `chi_load = 0` — everything,
byte-for-byte (87 bars + face-1 M0/M2/M3 evidence). At `chi_load = 1` —
dark bath byte-identical (structural); M1 chain `v_M` within the face-1
scaling curve with the predicted Kerr correction (sign: `v_M` shifts *up* —
dressing raises `x`, lowers `Δ`, raises the matter weight; magnitude from
§II.4 slaved formula); DS vacuum loci unchanged within the measured vacuum
dressing correction; atom-scale settlement and chi-filling (face-1
functional numbers) reproduced within tolerance stated per-run.

---

# PART III — Numerical simulation method (II.3)

## III.1 The knob: `chi_load` (byte-inert, default 0)

- **Name/type/default:** `chi_load`, double, **default 0**. `0` = face-1
  behavior byte-for-byte. `1` = the lawful booking (L3). Values in `(0,1)`
  legal for off-law diagnostics only (threshold-hunting sweeps), never a
  campaign law.
- **The one physics change** (C kernel, pass 1, freecell.c:1401–1409; Go
  mirror in `fab/step.go` pass 1):

```c
double x = (Em[i] + flload[i]) / P.cap;
if (P.chi_load > 0 && P.door_g > 0)
    x += P.chi_load * (dc1_[i]*dc1_[i] + dc2_[i]*dc2_[i]) / P.cap;
double det = 1.0 + P.q_detune * x;
w1e[i] = P.w1 / det;
w2e[i] = P.w2 / det;
```

- **Byte-inertness at 0:** the added term is behind an explicit branch
  (kernel convention — do not rely on `+0.0` IEEE identity). At
  `chi_load = 0` the pass-1 arithmetic is untouched: byte-for-byte face 1,
  hence byte-for-byte laws_V3a at `door_g = 0`. **Structural dark-bath
  identity:** at `chi_load = 1` with `χ ≡ 0` (dark bath) the branch adds
  exactly `0.0/cap` — assert byte-identity as a *measured* bar (U0), not
  an assumption.
- **Downstream consumers, enumerated (the spec is: ONE x, consistently):**
  `w1e/w2e` (all precessions, all beat rates — pass F, pass U, pass P, the
  `cbeta` integrand), the atom sizes `eF/eD ∝ w1e/w2e`, the radiance comp
  factor `w2/w2e = 1+qx`, and any diag that prints `x`. **Explicitly NOT
  changed:** the workfn `econd` (reads `Em` — presence rule, untouched);
  the cap/evap test (reads `Em + Ee` — transit is not holdings; the same
  rule that exempts `flload` exempts `|χ|²`; grok item 6 resolved by
  documentation, here); radiance's holdings base `x_h = Em/cap`
  (freecell.c:2583 — decision R-D1 stands: flight is not radiatable).
- **State: none added. Operations: none added. Passes: none added.** Face 3
  deliberately has the smallest kernel footprint of any v94 face: the
  theorem (§I.3) says the structure is already right; what was missing is
  the law statement, one booking line, the derived thresholds, and the
  bars.

## III.2 Prerequisites folded in (carried debt, from grok NEEDS-WORK)

Land with (or before) the `chi_load` edit, both kernels:

1. **Cell-death/repack χ flush** (`field_inject(i, |χ|²)`, zero `dc1/dc2`,
   before roster removal) — required by the conservation theorem's boundary
   clause (§II.7); currently unexercised (NC fixed mid-run) but the theorem
   must be complete in code.
2. **FCS register plane** (`dc1/dc2` in the dump behind `door_g > 0`, header
   flagged) — needed for U1(b)/M4/M6 diagnosis (the `chi2`-profile-locks-to-
   `Em`-profile bar needs per-cell χ).
3. Confirm the ΔEm sign fix + BH χ flush are in (they are, per
   IMPLEMENTATION_STATUS — verify at battery time).

## III.3 Schedule

No new pass. Pass order unchanged: … pass U → 3–5 → P(τ/2) → F → P(τ/2) →
6 → … `chi_load` acts entirely through pass-1 coefficients computed once
per step (the standard one-step lag of `x` behind the state it summarizes —
same lag discipline as `Em`/`flload` today; no new race).

## III.4 The winding-grain meter (U4 tooling; diag only, no physics)

Per-run counters (behind `door_g > 0`, print in the RESULT block):
`beat_wraps_total` (increment where `cbeta` wraps), `clicks_total`
(settlement fires), `click_atoms_total` (atoms moved), plus the click-rate
× load histogram (`ν_click` binned by `x` deciles) for the layer-5 throttle
curve. Bar U4 asserts: clicks ⊆ wraps (every click at a wrap — should be
exact by construction; the meter certifies no other path fires atoms), and
`ν_beat(x) ∝ 1/(1+qx)` measured on the histogram.

## III.5 C / Go mirror parity

- `chi_load = 0`: byte-identity C ≡ Go (standing bar, both kernels).
- `chi_load = 1`: byte-identity is *expected to hold* (the change is one
  deterministic arithmetic line in pass 1, no reordering, no trig-LUT
  divergence) — assert byte-identity first, fall back to tolerance-abx only
  if a pre-existing 1-ulp C↔Go divergence (the P-ID1 precedent) surfaces;
  record which.
- The abx surface must cover the triple `(amp_nat = 2, door_g > 0,
  chi_load = 1)` before any campaign conclusion is graded (extends face-1
  M0's closure of the Go-only debt).

## III.6 Acceptance bars (pre-registered; the ratchet)

**Invariant surface:** full battery, 87 V3a bars GREEN at `door_g = 0`
(byte-inert), both kernels, before and after every edit. Bars sharpen by
measurement, never soften.

- **U0 (inertness + parity):** battery 87 GREEN at `chi_load = 0` both
  kernels; C ≡ Go per §III.5; **dark bath at `chi_load = 1, door_g = 0.3`
  byte-identical to `chi_load = 0`** (the structural χ≡0 check).
- **U1 (existence of the self-bound state):**
  **(a)** within-mode threshold locus: off-law `(q, Δx)` sweep (battery
  precedent for off-law q arms: the xs bars) matches the §II.3 hyperbola
  `q·Δx/(1+qx̄)² = const` in shape, anchored at the measured B.1 point.
  **(b) [NEEDS face 2 + M6]** metabolic park at `chi_load = 1`: M6
  criteria (late slope 10× under SLOW_REFILL; ±20% perturbation decays)
  AND the register profile locked to the matter profile (spatial
  correlation `corr(chi2_i, Em_i) ≥ 0.9` on core cells) AND the measured
  restoring rate consistent in sign with the §II.5 inequality.
- **U2 (conservation):** `|drift_rel| ≤ 1e-13`, `chi_load = 1`, lit matter
  (M1 chain + 3-D blob), `T ≥ 1000`, both kernels. Per §II.7 this can only
  catch coding errors — which is exactly why it runs first.
- **U3 (anti-ignition + the fold):** (i) incoherent lit bath at
  `chi_load = 1`: birth/glow/cond stats within noise of `chi_load = 0`,
  `Σ|χ|²` bounded, no secular growth; (ii) **the fold sweep**: coherent
  field amplitude ramp over a seeded clump at `door_g ∈ {0.15, 0.3, 0.5}`;
  measure the dressing-jump threshold `Ee_c(door_g)`; **bar:** `Em`
  bounded by cap defense on both sides of the fold; **sharpener:** the
  measured `Ee_c` scales as `1/g²` (two decades of `door_g²` → slope −1 on
  the log-log). The threshold *number* is recorded physics (the
  pair-production intensity), whatever it is.
- **U4 (grain accounting):** clicks ⊆ beat-wraps exact; atoms integer on
  the ledger; `ν_beat(x)` follows `1/(1+qx)` (the layer-5 throttle,
  measured); ħ-linearity spot-check re-run off-battery at
  `chi_load = 1, door_g > 0` (the atom sizes moved with `w_e` — linearity
  in `A₀` must survive; this is the U4 arm that can *fail*, not just
  certify).
- **U5 (face-1 correspondence):** the §II.9 must-reproduce set. In
  particular the M1 chain re-run at `chi_load = 1`: `v_M` within the
  face-1 curve + predicted upward Kerr shift; and the DS loci at strength
  (`door_g > 0, chi_load = 1`): vacuum loci unchanged within the measured
  vacuum-dressing correction (record it; it is the O(g²/Δ₀²) number face 1
  flagged), loaded-arm shift recorded as the Kerr measurement.

Run order: U0 → U2 → U4 → U3 → U5 → U1(a); U1(b) queues behind face 2 +
M6. U3(ii) **before** any lit-bath campaign at `chi_load = 1` (the M3b
discipline, inherited).

## III.7 Most likely failure mode, honestly stated

**The fold sits inside campaign intensities.** At `door_g = 0.3` the
derived `Ee_c ≈ 1.1` per cell (vacuum-adjacent; ≈ 0.8 over chain-loaded
matter) — and face-1 functional runs already push per-cell field energies
of that order in lit baths. If the measured fold is at or below campaign
amplitudes, every `chi_load = 1` lit run near matter is bistable: dressing
(and hence index, beat rate, atom size) jumps hysteretically between
branches, and optics/metabolic measurements acquire a branch label that
must be tracked, or the campaign drops to `door_g ≤ 0.15` (fold at ≈ 4.3,
comfortably above traffic) and pays 4× in every `g²` rate. Second failure:
the adiabatic estimate is O(1)-wrong *and* the shape is wrong (the fold is
not observed as a jump at all because beat-sampling + non-slaved dynamics
smear it) — then §II.4's fixed-point analysis is the wrong reduction, the
Kerr claim survives only qualitatively, and U3(ii)'s sharpener fails while
its safety bar still passes. Third (inherited, scoped): M4 fails — phase
does not survive the grained edge — in which case the "one algebra" claim
demotes to the magnitude sector, U1(b)'s profile-lock bar likely fails with
it, and the registered fallback is the J-dependent link metric composition
(face-1 §III.8). None of these failures touches the theorem of §I.3, which
is falsified only by the refutation condition stated there.

---

# PART IV — Flags (II.4)

## IV.1 What face 3 completes vs the permanent residue

| | status |
|---|---|
| One primitive on every edge (spatial-e, spatial-m, fiber e–χ, fiber χ–m) | **complete** (faces 1–2 build, face 3 states as law; phase scope graded by M4) |
| One pitch functional, all carriers | **complete** (was already true; now law) |
| One load booking holdings + all transit (space + mode) | **complete at `chi_load = 1`** — the face-3 knob; adoption is the user's |
| Click timing + integer-ness | **derived** (winding topology; the kernel already implements it — face 3 re-reads, adds the U4 meter) |
| Click scale `A₀` | **permanent residue** — the substrate's ħ; deriving it = second quantization = a different programme |
| Two-edge fiber `e — χ — m` (register ≠ inventory) | **permanent residue, proven minimal** (§I.3); `ψ_m` is NOT retired |
| Matter as bound state | **re-derived as configuration *of the dense sector*** (self-trap §II.3 + dressing §II.4 + park §II.5), not as configuration of the field — the distinction is the programme boundary (§I.5) |
| `Es` in the spinor algebra | **not attempted** — see IV.4 |

## IV.2 Interaction with face 2 (and the other in-flight arms)

- **Face 2 (stimulated outtake):** load-bearing for U1(b) (no park without
  outtake — §II.5) and for M5/M6. Face 3's knob composes with it trivially
  (outtake reads the same `x` through the same pitches). **The fold sweep
  U3(ii) must be re-run once outtake lands** — outtake phase-matched to a
  coherent bath is the amplifying half the face-1 M3b could not exercise
  (grok §3.2); the §II.4 threshold prediction applies to the *dressing*
  fold either way, but the post-fold fate (bounded conversion vs sustained
  pump-cycle) is a face-2+3 joint measurement.
- **M1:** supplies `v_M(door_g)` → fixes `W_e` in the dictionary; U5
  re-runs the chain at `chi_load = 1`.
- **M3b:** its measured threshold is the calibration point for `Ee_c`.
- **M4:** scopes the phase half of (L1) — see NEEDS.
- **M5/M6:** unchanged by face 3 at `chi_load = 0`; re-run at
  `chi_load = 1` only after U0–U4 are green.

## IV.3 Programme rules touched (argued, for the user's judgment)

1. **The one-amplitude boundary.** Face 3's principal finding: the
   reduction to matter-as-field-configuration is the *rejected* side
   (§I.5) — matter-as-profile on a persistent field is the pre-v89
   architecture, and the mode boundary is the v89 principle's content. The
   lawful unification is bundle-level (one fiber field, one algebra), not
   carrier-erasure. **No rule is broken; the rule is *derived*** (§I.3
   turns "the door stays special" from a constraint into a theorem with
   named irreducibles).
2. **`chi_load` touches the law's `x`.** The load definition sits inside
   every pitch in the substrate — this is a law-adjacent change, behind a
   byte-inert knob, default 0. The argument that it is the *lawful
   completion* rather than a new law: v89's own text ("flight is load")
   booked space-transit into `x`; `|χ|²` is mode-transit (the `slem`
   cousin, grok-ACCEPTED); the current `chi_load = 0` state is the
   asymmetry. Adoption nonetheless requires the user + full battery + a
   reckoning, per the ratchet.
3. **No new constants** (P-noscale): `chi_load` is a booking switch;
   `q_detune`, `cap`, `w2`, `door_g` are reused. The brief's candidate
   `Ω_χ = w_χ/(1+q_χ|χ|²/cap)` is realized with `w_χ = w2`, `q_χ =
   q_detune` — derived from (L2)+(L3), not posited with free parameters.
4. **No feedback/shadow/bias:** face 3 adds no term to any want or rate —
   only the load booking; the pitch functional is the adopted one. The
   Kerr feedback is the *law's own* nonlinearity applied to the new
   ledger, in the open, with its fold derived and its bars pre-registered.
5. **Cap semantics documented (grok item 6):** transit (`flload`, `|χ|²`)
   loads *pitch*, never *holdings*; cap/evap remain on `Em + Ee`. One
   rule, both transit families, stated here as the standing text.

## IV.4 The next rung of the symmetry (flagged, NOT for implementation)

The user's symmetry, applied once more: space-transit has a register
(`slem`) and mode-transit now has one (`χ`); but the m–s conversion (the
`s_pull` backsplash, condensation's space book) is today a magnitude-only
click ledger — the exact structural position the e–m door occupied before
face 1. The pattern predicts a **space register** `σ` (an m–σ–s fiber
extension: virtual space-conversion, with `Es` curvature as its settled
sector) — the polariton mechanism for *gravity*, where the v89 g-bars
(footprint g1, opacity g3, the missing 1/r far field g4) live. This is
speculation beyond the current evidence; the far-field lane's measured
prerequisite (medium-carried books, v91 ASTRO) is consistent with it. **Do
not implement**; recorded so the symmetry's direction is on file when the
g-lane reopens.

---

*End of the original DESIGN_claude_face3.md. Orchestrator: implement
§III.1–III.5 (one booking line + meters + the two carried prerequisites),
run §III.6's U-bars in the stated order, U3(ii) before any `chi_load = 1`
lit-bath campaign; resume this session (ID at head) when face 2 / M1 / M3b
/ M4 land, per the NEEDS block. §X below is resume update 1.*

---

# §X. Resume update 1 (face-2 + M3b landed, 2026-08-09)

**Source:** `v94/consult/FACE2_M3B_RESULTS.md` (read in full). Same session
(`aebbc7eb-4716-4e30-8c22-484e115428f2`), first resume. This section
calibrates §II.4, updates §II.5, and re-grades the affected bars. Nothing
above is edited except the NEEDS block (resolution markers); where §X
supersedes a number or a grading, §X governs.

## X.1 What landed

Face 2 (stimulated outtake) is implemented in both kernels behind the same
`door_g` knob, byte-inert at 0, battery 87/87 GREEN — **the bidirectional
cycle is closed.** The outtake fires when the corrected-convention
cross-term `C = 2·Im(ψ_m*χ) > 0` drives `ΔEm < 0` (note: face-1 §5.3's
draft text said `C < 0`; that predates grok's ΔEm sign fix — the results
file states the corrected condition, and this doc's §II.1 language "the
negative branch" should be read with `C > 0`). M3b ran as a
`noise_amp` ramp (0.3 → 2.0) over a seeded blob at `door_g = 0.5`,
`amp_nat = 2`, `wf_on = 1`, Go kernel: onset of `door_in` at `na ≈ 0.6`
(per-cell `Ee ≈ 0.36`), channel-dominance knee at `na ≈ 1.0` (`Ee ≈ 1.0`),
outtake onset at `na ≈ 1.2`, **no runaway through `na = 2.0`**
(`Ee ≈ 4.0`), drift at the FP floor throughout (post-fix, X.7).

## X.2 §II.4 calibrated: the fold is real, bounded, and O(1)-located

The §II.4 prediction at the tested coupling: `Ee_c(door_g=0.5, x₀=0) =
0.39` (0.28 over chain-loaded matter). Measured knee: `Ee ≈ 1.0`.

- **Prefactor: measured/predicted ≈ 2.6.** This is within the
  pre-registered O(1) trust band, and the direction is expected: at
  `door_g = 0.5`, `g = 1.45 > Δ₀ = 1.25` — the slaved (adiabatic)
  approximation behind the closed form is outside its validity range
  (`g ≲ |Δ|`), so an O(1) underestimate of the threshold is the predicted
  failure of the estimate, not of the mechanism. **Calibrated fold:**
  `Ee_c ≈ 2.6·Δ₀²·cap/(4g²q(1+qx₀))` (prefactor pinned at one point).
- **Boundedness: confirmed as designed.** §II.4 claimed the post-fold fate
  is bounded conversion on a seed, defended by cap+evap — and now, with
  face 2, by the stimulated outtake as a *derived* return channel
  (`door_out` engages at `na ≈ 1.2` and scales with drive: 0.86 → 5.3 →
  22.4 across `na` 1.2 → 1.5 → 2.0). The defense is structural, not an
  imposed gate — exactly the §II.8 stack doing its job, with face 2 as an
  emergent sixth layer: **stimulated return** (above the outtake onset,
  the same phase geometry that pumps intake also opens the drain).
- **Supralinearity is the fold's signature.** `door_in ∝ na^3.5 ≈ Ee^1.75`
  above onset — steeper than the linear-response `∝ Ee`. Below the fold
  the self-consistent dressing amplifies as `du/dβ > 1` (diverging at the
  fold before saturation cuts it); an exponent between 1 and the fold's
  formal divergence is qualitatively the §II.4 branch structure. Recorded
  as consistent, not claimed as confirmed — the exponent is one number
  from one geometry.
- **The surviving sharpener (the shape test), sharpened with calibrated
  numbers.** The `1/g²` slope is still unmeasured (one `door_g` point).
  Discriminating pre-registration, now with the calibrated prefactor: the
  **onset** (`Ee ≈ 0.36`) is grain-set (X.3) and should be roughly
  `door_g`-independent, while the **knee** should move as `1/g²`:
  predicted knees `Ee_c ≈ 2.8` (`na ≈ 1.7`) at `door_g = 0.3` and
  `Ee_c ≈ 0.51` (`na ≈ 0.7`) at `door_g = 0.7`. A 3-point ramp at
  `door_g ∈ {0.3, 0.5, 0.7}` spans reachable drive (`na ≤ 2`) and tests
  slope −1 on the log-log plus the onset/knee decoupling. This replaces
  U3(ii)'s original sweep spec.

## X.3 The onset is the atom grain — a U4 confirmation for free

The measured `door_in` onset (`Ee ≈ 0.36` per cell) coincides with the
atom grain at the loaded dense pitch: `ε = A₀·w2e/2π = 1.15·2.086/2π ≈
0.38` (blob cells at `x ≈ 0.325`; results file: "just below atom grain
0.37"). So the *first* gate met in practice is exactly (L4)/P-grain — no
click until one whole source-cycle's action is available in the register —
while the *knee* is the §II.4 dressing fold. Two distinct thresholds, both
now measured, matching the two distinct pieces of the design (imposed
scale `A₀`; derived fold). This is a quantitative down-payment on bar U4
before the winding meter even exists.

## X.4 §II.5 updated: the cycle is closed; the park inequality with the stimulated channel

The §II.5 tentative claim "intake-only cannot park" is now moot in the
right way: the outtake half exists and fires. Updates to the park
analysis:

- The stimulated channel enters **both** sides of the balance with the
  same Kerr factor (`I_stim` and `O_stim` both carry `F(x)`'s
  door-opening term `2q(1−F)/(1+qx)` and the beat throttle), so their
  contributions to the restoring criterion partially cancel; the radiance
  `x⁴` wall (`4/x*` ≈ 6.45) still dominates the inequality, and the ~30×
  margin conclusion stands unchanged.
- The measured **7:1 intake:outtake ratio** at `na = 2.0` says this blob
  geometry is deeply intake-dominant under strong coherent drive — the
  balance at a quiet park will be struck mostly by radiance + evap against
  intake unless the park's phase geometry favors outtake. Pre-registered
  for M6: measure the park's outtake **composition** (spontaneous radiance
  vs stimulated `door_out` vs evap) — the design predicts radiance-led
  balance at law constants (`k_rad = 0.05`, `p_rad = 4`), with `door_out`
  a minority channel whose share grows with local phase order. If M6
  parks with `door_out` as the *majority* channel, §II.5's rate-level
  idealization needs rebuilding.

## X.5 Bars re-graded

| bar | previous | now | basis |
|---|---|---|---|
| U3 safety (bounded coherent pumping) | pre-registered | **GREEN at `door_g = 0.5` through `na = 2.0`** (Em bounded, drift FP floor, outtake + cap/evap defend) | M3b ramp |
| U3 sharpener (`1/g²` knee slope) | pre-registered | **OPEN — spec replaced by X.2's calibrated 3-point test** | one coupling measured |
| U1(b) (metabolic park) | NEEDS face 2 + M6; tentative "impossible without outtake" | **NEEDS M6 only; upgraded to "possible, awaiting measurement"** — both half-cycles now exist and fire (in/out both ledgered) | face 2 + functional table |
| U4 (grain accounting) | pre-registered, meter unbuilt | pre-registered, **plus the X.3 onset≈ε quantitative anticipation** | M3b onset |
| U2 (conservation at `chi_load = 1`) | pre-registered | unchanged — `chi_load` not yet implemented; note the face-2 drift record (FP floor at all `na`, post-fix) covers the `chi_load = 0` face-2 surface | results table |
| U0 / U5 | pre-registered | unchanged (`chi_load` not yet implemented) | — |

The face-3 kernel work (§III) remains exactly as specified; nothing in
these results changes the `chi_load` spec. Run order unchanged, with
X.2's 3-point fold test folded into U3(ii).

## X.6 Parity language corrected (M0 / §III.5)

The measured C↔Go divergence at `amp_nat > 0` (e.g. `cond` 1168 vs 1116,
5%, at `door_g = 0`) confirms the **pre-existing v93 debt** — it predates
and is independent of the door pair (face-1 M0's byte-identity result was
config-specific, not a general `amp_nat > 0` guarantee). Consequences for
this design's §III.5: (a) the battery surface (`door_g = 0, amp_nat = 0`)
stays C ≡ Go byte-identical — the standing bar; (b) **all v94 campaign
data are Go-only** (`fabrun`), consistent with v93; (c) §III.5's
"byte-identity expected at `chi_load = 1`" claim is scoped to *within the
Go kernel* (`chi_load = 1` vs `0` comparisons) and to the battery surface —
cross-kernel abx at `amp_nat > 0` is tolerance-based until the C-side
`amp_nat` debt is paid, which is upstream of and orthogonal to face 3.

## X.7 The conservation theorem earned its keep (the bs-cap bug)

The Go outtake bug (backsplash `bs` uncapped at post-Givens `emNew`; the
`Em = 0` clamp then created energy, drift `2.4e-5`) is exactly the failure
class §II.7 names: the theorem's continuous part is exact by construction,
so the only places implementation can break it are the **click-ledger
clamp clauses and boundary flushes** — and the drift meter, which the
theorem makes a tripwire rather than a tuning aid, caught it at 8 orders
above the floor. Post-fix `-2.5e-13` at the same config. Standing
implementation rule, restated for the face-3 edit and all future clicks:
*every clamp inside a click ledger must be two-sided (never clamp a state
variable without re-deriving the ledgered amount from the clamped move)* —
the C kernel already obeyed it; mirrors must be checked clause by clause.

## X.8 Updated NEEDS (for the next resume)

- **M6** → grades U1(b) (park + `chi2`/`Em` profile lock + X.4's
  composition split).
- **M4** → scopes the phase half of (L1); needs the FCS χ-plane (§III.2
  prerequisite — now also blocking M4 tooling, raise its priority).
- **M1** → `v_M(door_g)` fixes `W_e` in the dictionary; U5 chain re-run
  at `chi_load = 1` once the knob lands.
- **M5** → antagonism release; face-2-dependent, now unblocked.
- **U3(ii) 3-point fold test** (X.2) → the `1/g²` slope; runnable now,
  before the `chi_load` edit, and worth running first — it hardens the
  §II.4 calibration that `chi_load`'s own bars will lean on.

*End of resume update 1.*

---

# §Y. Resume update 2 (M4 negative, 2026-08-09)

**Source:** `v94/consult/M4_RESULTS.md` (read in full). Same session,
second resume. M4 is RED and the failure is *instructive*: it localizes a
sampling error in face 2's outtake that this design's own §II.6 language
already names. This section records the result, diagnoses the mechanism,
answers the orchestrator's direct question about `chi_load` (answer: **no
— and the suggested mechanism has the wrong sign**), and registers the
design-faithful fix with pre-registered bars. Where §Y supersedes §X or
the original text, §Y governs.

## Y.1 The result, taken at face value

- **M4 NOT MET.** Lit-bath m=+2 ring, `door_g = 0.5`: R2 dips to 0.15 by
  t=8 (identical to baseline — the early dephasing is bath contact, not
  the door) and then *decays* to 0.24 at t=80 while `Em_mean` collapses
  0.81 → 0.07. Books: `door_out = 8.56` vs `door_in = 1.83` — the outtake
  drains the ring 5:1.
- **The baseline is a discovery in its own right.** At `door_g = 0` the
  same ring on the same hot bath dips to 0.14 and **self-heals to R2 =
  0.62 at t=80 with `Em_mean` 0.81 → 1.16** — the bath *feeds* the ring
  through the legacy condensation door, and topology + the unitary
  channel re-establish W=+2. This supersedes the v93 cold-bath pessimism
  (`R2d ≈ 0.02`) on ring geometry: **hot-bath accretion + closed topology
  is a phase-retention mechanism the substrate already owns.** Record it;
  it resets the M4 bar's honest baseline from 0.02 to 0.62.
- So the door pair did not fail to *help* — it actively broke a healing
  process that works without it. The "do no harm" form of the bar (Y.4)
  is the right re-registration.

## Y.2 Diagnosis: the frozen-phase ratchet (face 1's amputation, recommitted)

The results file attributes the drain to the winding's C>0 geometry. The
mechanism underneath is sharper, and it is one this design already named:

**The cross-term is not a DC quantity.** `ψ_m` precesses at `w2e`; `χ` is
driven by the field at `w1e`; their relative phase advances by exactly 2π
per beat period *by construction* (that is what the beat is). So
`C(t) = 2·Im(ψ_m*χ)` **oscillates through a full cycle every beat**, and
its beat-instant value is a stroboscopic sample at a per-cell phase that
is *persistent* (`cbeta` is initialized once, random per cell, and wraps
coherently thereafter). A cell whose frozen sampling phase lands in the
`C > 0` half-cycle fires outtake **every beat**; the winding correlates
these phases around the ring and phase-slip sites reseed them. Face 2's
outtake thus samples a rotating quantity once per cycle and treats the
sample as geometry — which is precisely the amputation face 1 diagnosed
in the magnitude door ("sample once per beat, discard the cross-term's
dynamics"), recommitted one edge over.

**The gate asymmetry amplifies it.** Intake is grain-gated on the register
(`chi2 > ε ≈ 0.37`) which is `g²`-fed and rare at `na = 1.0`; outtake
draws on `Em ≈ 0.81` (≥ 2 atoms always available) and is gated only by
the sampled `C` sign — the click itself is *not* `g`-scaled (the solved
Givens moves atoms regardless of how slowly `χ` refills). Order-of-check:
`ν_beat = Δ₀/(2π(1+qx)) ≈ 0.143/t.u.` at ring load → ~11.5 beats/cell
over T=80; 6 cells × ~O(½) firing × ε ≈ 0.35 gives ~12 vs measured 8.56 —
the ratchet model is quantitatively consistent. **Testable from data
already in hand:** the M4 sweep ran `door_g ∈ {0.3, 0.5, 1.0}`; the
ratchet predicts `door_out` roughly *flat* across that sweep while
`door_in` rises ~`g²`. Orchestrator: please report those two columns —
it decides this diagnosis at zero cost.

## Y.3 The direct question answered: `chi_load` does NOT rescue M4

The results file (and the resume prompt) suggest the face-3 self-detuning
"would detune its own door, reducing the outtake rate exactly where it's
destructive." Honestly: **no — the sign is wrong.** Under (L2)+(L3), load
*shrinks* `|Δ| = Δ₀/(1+qx)` — that is the design's own derivation of the
workfn ("bound matter converts freely, vacuum demands a far threshold",
face-1 §2.4). A bound, dressed ring has a *more open* door, not a
detuned-shut one; `chi_load` cannot be simultaneously the workfn
derivation and an outtake suppressor. What `chi_load` does offer here:

- the **beat throttle** (`ν_beat ∝ 1/(1+qx)`): dressing adds
  `|χ|²/cap ≈ 0.1–0.2` to `x` at this bath strength → ~10–15% fewer
  clicks. Mitigation, not rescue.
- and the throttle **weakens as the drain proceeds** (`Em` ↓ ⇒ `x` ↓ ⇒
  beats *faster*) — the gradient is destabilizing, so the ratchet
  accelerates into the kill.

Verdict, stated for the record: `chi_load` changes M4's numbers by tens
of percent at most and does not change its sign. The §III spec is
unmodified, but no face-3 claim may lean on `chi_load` for phase
retention, and I withdraw the results file's suggestion (3) on my
design's behalf.

## Y.4 The design-faithful fix: integrate the cross-term (synchronous detection)

The fix is not a geometry gate bolted on (an `Em_min` floor or a
winding-gradient veto would be a hand-tuned bias — the v92 closed road).
It is the design's own §2.5/§II.6 language implemented literally: *"the
click is the integrator-of-fire for the pair"* — face 2 integrated the
magnitude but point-sampled the phase. Registered revision (face-2 rev,
one knob):

- **State:** one persistent per-cell accumulator `dcc_i` (double),
  the cross-term sibling of `cbeta`:
  `dcc_i += 2·(m1·dc2 − m2·dc1)·dt` each step (values read post-pass-P,
  pre-settlement; zeroed at each beat wrap after use). Same state class
  as `cbeta` (a per-cell phase-cycle accumulator — the precedent), one
  new float: flagged for user signoff like `dc1/dc2` were.
- **Knob:** `door_avg` (default 0 = face-2 instant-C behavior, preserved
  for reproducibility; 1 = the beat-integrated law). Byte-inert at 0;
  battery untouched either way (`door_g = 0` skips all of it).
- **Law at `door_avg = 1`:** at the beat, the outtake branch keys its
  direction and its demand on `D = dcc_i` (the beat-integrated
  coherence), not on the instantaneous `C`: fire outtake only if `D > 0`
  and the solved Givens can move ≥ 1 atom; scale the demand by the
  coherent fraction (`|D|` against its full-lock bound
  `2|ψ_m||χ|·T_beat`), clamped by atoms as always. Intake keeps its
  magnitude gate (`chi2 > econd`) unchanged.

**Why this is the physics and not a patch.** Over one full beat the
relative phase completes one cycle, so the reversible (oscillatory)
component of `C` integrates to ≈ 0 — `D` performs synchronous detection
at the beat frequency. Only *pitch-locked* coherence survives the
integral: a `χ` dragged by its own matter (a mover's trailing edge, a
coherent drive resonant with the dense pitch) accumulates `D` linearly;
a `χ` fed by a beating or incoherent field averages out. Consequences,
derived and pre-registered:

1. **Bound static structures stop draining** (`D ≈ 0` on the ring —
   the frozen-phase ratchet is exactly what the integral kills).
   **Bar M4′ (replaces M4, sharpened to the measured baseline):** lit-bath
   m=+2 ring, `door_g = 0.5, door_avg = 1`: `R2@80 ≥ 0.62` (do no harm —
   the pair must not underperform the legacy door) with door traffic
   verified `> 0`.
2. **Movers keep their outtake** (translation-by-conversion needs the
   trailing edge): the packet's own emitted field is phase-locked to its
   matter, `D` accumulates. M5 is the measurement; if M5 stalls at
   `door_avg = 1` but moved at `door_avg = 0`, the lock claim is wrong
   and the fallback (below) engages.
3. **Reciprocity restored:** `D ∝ |χ| ∝ g` — outtake becomes `g`-fed
   like intake, closing the Y.2 gate asymmetry (and the Y.2 sweep check
   verifies the asymmetry existed).
4. **Anti-ignition:** strictly safer on the drain side (a coherence
   requirement was added to a channel that returns energy fieldward).
   One honest cost: the M3b "stimulated return" defense (§X.2's sixth
   layer) is weakened where the bath's drive is incoherent — but M3b's
   boundedness never depended on it (`door_out = 0` for all `na ≤ 1.0`
   and `Em` stayed bounded by cap+evap alone). **Re-run the M3b ramp at
   `door_avg = 1`** as the closing bar; expect unchanged boundedness with
   smaller `door_out`.
5. **M6 risk, pre-registered before it runs:** under `door_avg = 0` the
   same frozen-phase ratchet will slow-drain any quiet clump whose `χ` is
   bath-fed (predict: if M6 fails to park, its books will show
   `door_out`-dominant outtake composition — X.4's meter). If so, the
   cause is Y.2, and M6 re-runs at `door_avg = 1` before any design
   conclusion is drawn.

**Fallback if M4′ also fails:** the J-dependent link metric (RECKONING's
named next face, face-1 §III.8's registered fallback) — which the results
file independently re-names as still necessary. The two are compatible;
`door_avg` is the cheaper and more principled first move because it
removes a demonstrated *error* (sampling a rotating quantity) rather than
adding structure.

## Y.5 Claims and bars re-graded

| item | previous | now |
|---|---|---|
| M4 | target `R2d ≥ 0.5` with traffic | **RED** (0.24 vs baseline 0.62; pair worse than no pair) |
| M4 baseline | v93 cold-bath `R2d ≈ 0.02` | **superseded: hot-bath ring self-heals to 0.62 at `door_g = 0`** — recorded as standing physics (accretion + topology heal phase) |
| M4′ (new) | — | pre-registered: `R2@80 ≥ 0.62` at `door_g = 0.5, door_avg = 1`, traffic > 0 |
| (L1) "one algebra", phase half | assumed pending M4 | **demoted as pre-registered** (NEEDS block): currently magnitude-sector only; restoration rests on M4′ (+ M5 at `door_avg = 1`), else the J-metric composition |
| §I.3 theorem | proven-sketch | **untouched** — it claims structural minimality, not retention; M4 measures the dynamics, not the hypotheses (H1)–(H3) |
| U1(b)/M6 | awaiting M6 | awaiting M6 **with the Y.4(5) drain risk pre-registered**; M6 books must be read against it |
| §X.2 sixth layer (stimulated return) | emergent defense | **scoped:** operative at `door_avg = 0`; partially traded away at `door_avg = 1` (M3b re-run is the closing bar) |
| `chi_load` spec (§III) | as written | **unchanged** — but explicitly *not* a phase-retention mechanism (Y.3); its bars (U0–U5) unaffected |

## Y.6 Updated NEEDS (third resume)

- **M4-sweep columns already on disk:** `door_in`/`door_out` at
  `door_g ∈ {0.3, 1.0}` — decides the Y.2 ratchet diagnosis for free.
  Please include in the next resume regardless of what else lands.
- **M4′** after `door_avg` is implemented (with the Y.4 spec; one float +
  one knob, user signoff on the accumulator flagged).
- **M5** — now doubles as the mover-side test of the `D`-lock claim
  (run both `door_avg` values if cheap; the pair decides Y.4(2)).
- **M6** — read its outtake composition against Y.4(5) before grading
  U1(b).
- **M1** — unchanged (propagation is pass-P physics, upstream of the
  settlement sampling question; `v_M` still fixes `W_e`).
- **M3b re-run at `door_avg = 1`** — the anti-ignition closing bar for
  the revision.

*End of resume update 2.*

---

# §Z. Resume update 3 (M6 / M1 / M5, 2026-08-09)

**Source:** orchestrator resume-6 results (M6 park approach at
`chi_load = 1`; M1 confinement-not-propagation; M5 no antagonism
release). Same session, third resume. This is the biggest single update
to the design's grading: the face-3 knob scores its first physics win
(§Z.1), two of the design's central dynamical predictions fail as run
(§Z.2, §Z.3), and the cross-cutting diagnosis (§Z.4) says both failures
and M4 share one cause — the still-unimplemented Y.4 sampling fix —
while the door pair's *conservative* character emerges as a measured
virtue nobody pre-registered.

## Z.1 M6: the park approached — and the face-3 knob is the restoring mechanism

The measured arms (blob + lit noise, `door_g = 0.5`, outtake enabled,
T ≈ 520):

| t | `chi_load = 0` | `chi_load = 1` |
|---:|---:|---:|
| 0 | 186 | 186 |
| 210 | 220 | 186 (recovered to initial) |
| 520 | 282, slope ≈ +0.5/t.u. (SLOW_REFILL) | **192, slope ≈ +0.05/t.u. — 10× slower** |

The `chi_load = 0` arm reproduces v93's `rad_dual` SLOW_REFILL (climb,
no park). The `chi_load = 1` arm recovers to its initial inventory and
flattens. **The only difference between the arms is the face-3 booking
line** — `|χ|²` in the load — so the restoring physics is the
dressing-aware pitch law, exactly the knob §III specified.

**Mechanism attribution, corrected against my own §II.5.** The §II.5
rate analysis treated the radiance `x⁴` wall as the dominant restoring
slope and the Kerr/dressing terms as small corrections. The measured
differential says otherwise: both arms share the radiance law; the
*difference* — the thing that turns SLOW_REFILL into near-park — is the
dressing entering `x`. The §II.8 derivation contains the mechanism and
even flagged its sign as regime-dependent: per-time intake
`∝ ν_beat(x)·F(x)` with `ν_beat ∝ 1/(1+qx)`, and for `2g > Δ` (at
`door_g = 0.5`, `2g = 2.9 ≫ Δ₀ = 1.25` — deep in the over-coupled
regime) the beat-throttle factor *wins*: intake per unit time **falls**
as the object loads and dresses. At `chi_load = 1` a feeding blob's own
dressing throttles its click clock — a negative feedback with no
counterpart at `chi_load = 0`, where the dressing is invisible to every
pitch. §II.8's "layer 5" is not just an anti-ignition layer; it is the
park's restoring force. **[measured at one coupling; the sign flip
predicted for `2g < Δ` (door_g ≲ 0.2) is now a registered sharpener —
a `door_g` sweep of the late slope should show the park weaken as the
coupling drops below the over-coupled regime]**

**Grading discipline:** this is *evidence of approach*, not a confirmed
park — U1(b) stays open pending T ≥ 2000, the ±20% restoring
perturbation test, and the outtake-composition split (X.4). But the
10× slope reduction on the programme's oldest open ledger, delivered by
the one line of kernel physics face 3 added, is the design's first
measured payoff. U1(b) upgrades: "possible, awaiting measurement" →
**"approaching; mechanism identified; confirmation pending."**

## Z.2 M1: the door pair confines — the metabolic wave is NOT observed

| t | `door_g = 0` (legacy) | `door_g = 0.5` (pair) |
|---:|---:|---:|
| 0 | Em = 20 | Em = 20 |
| 100 | **Em = 417** (runaway box-fill) | **Em = 17** (confined, stable) |

Two findings in one table, and the unregistered one first:

- **The legacy door under lit noise is a runaway accreter** (20 → 417:
  the magnitude-only conversion ratchets bath energy into matter with
  no reversibility). The door pair, same drive, holds ~17. The pair is
  not merely no-worse than legacy on ignition — it is the *cure* for a
  legacy runaway mode. This is face-1 §4.3's layer structure (norm
  preservation, detuning, reversible dressing) doing measured work, and
  it belongs in any adoption case for `door_g > 0`.
- **The metabolic wave (`v_M`, face-1 §3.3) does not propagate.** No
  arrival front; the pair *confines*. The face-1 derivation was linear
  physics of the continuous pair; the real channel is atom-grained and
  beat-gated, and two graining effects the linear theory ignores are
  each sufficient to kill a travelling front: (i) **the grain floor** —
  a front cell must accumulate `chi2 ≥ ε ≈ 0.38` before its first
  click, so the wave's leading tail (which is exponentially dim) can
  never advance it — the M3b onset result applied spatially; (ii) **the
  outtake erosion** — measured `door_out` dominance (Z.3) drains any
  structure the intake cascade builds. The `v_M(door_g)` curve, the
  dictionary's `W_e` entry, and U5's Kerr-shift prediction on `v_M`
  are all **moot as specified** — amended to "not measurable under
  face-1+2 as implemented."

The honest summary sentence for the record: **the door pair is a
conservative, confining mechanism — the opposite of the legacy door's
dissipative spreading** — which is precisely what the no-cancellation
analysis (NOCANCELLATION §3–§4) says the rotation primitive should be,
and which is a *virtue* for stability and a *wall* for transport.
Propagation needs a mechanism the pair alone does not supply.

## Z.3 M5: the antagonism stands — and the test is confounded

Self-trapped blob, `door_g ∈ {0, 0.5}`: COM deviation 0.002 vs 0.001 at
t = 80 — anchored both ways; translation-by-conversion (face-1 §2.7)
did not engage. Books: `door_in = 20.3` vs `door_out = 87.8` — the
outtake drains the blob 4:1 before any front/rear asymmetry could move
it. This is the M4 frozen-phase ratchet operating on a blob: the
instant-C outtake fires on the oscillating cross-term's frozen sampling
phase, isotropically, and eats the object. **M5 as run is therefore RED
but confounded** — it does not grade translation-by-conversion, because
the trailing-edge mechanism it was designed to test requires outtake
that fires on *persistent* phase geometry, which is exactly what the
unimplemented Y.4 revision (`door_avg`, beat-integrated D) provides.
Re-registered as **M5′: rerun after `door_avg` lands**, with the honest
caveat carried from §Y: if M5′ still stalls, translation-by-conversion
is wrong as a mechanism and the J-metric face is the only registered
road to mobile matter.

## Z.4 The structural pattern, and the corrected priority order

Across M4 / M5 / M1 / DS-G-with-outtake, one cause recurs: **the
face-2 instant-C outtake erodes every phase-structured object it
touches**, while everything face-1-only is clean (DS-G3 `r = +0.684`;
M1 confinement; M3b boundedness). Two fixes are on the books and they
are not competitors — they are sequential:

1. **`door_avg` (Y.4, synchronous detection) — implement FIRST, still
   pending.** It removes a demonstrated *sampling error* (firing on a
   frozen sample of a rotating quantity), costs one accumulator float
   and one knob, and un-confounds every phase-sensitive bar: M4′, M5′,
   the M6 outtake-composition split, and DS-G with outtake enabled all
   wait on it. No conclusion about the metabolic cycle's phase half is
   gradeable until it lands.
2. **The J-dependent link metric (RECKONING's face) — the propagation
   candidate.** `door_avg` stops erosion but does not create transport;
   Z.2 says the pair alone confines. If mobile, metabolizing matter
   exists in this substrate, its transport plausibly needs the link
   metric to make the *spatial* channel phase-structure-aware, with the
   door pair supplying conversion and confinement. That composition —
   pair for the mode direction, J-metric for the space direction — is
   now the registered architecture hypothesis for v95.

## Z.5 Bars and claims re-graded (cumulative through resume 6)

| item | grade after resume 6 |
|---|---|
| U1(b) metabolic park | **approaching** — 10× slope reduction at `chi_load = 1`; mechanism = layer-5 beat throttle (measured differential); confirmation needs T ≥ 2000 + perturbation + composition split |
| `chi_load` (the face-3 knob) | **first measured physics win** (the park differential); its `door_g`-sweep sign-flip sharpener registered (Z.1) |
| §II.5 rate analysis | **corrected** — radiance wall is not the operative restoring slope; the dressing throttle is (§Z.1) |
| M1 metabolic wave / `v_M` / dictionary `W_e` | **not observed; moot as specified** — pair confines; grain floor + outtake erosion each sufficient (Z.2) |
| M1 side-finding | **pair cures a legacy runaway** (417 vs 17) — adoption-case evidence, unregistered discovery |
| M5 antagonism release / translation-by-conversion | **RED as run, confounded by the ratchet** → M5′ after `door_avg`; if M5′ fails too, J-metric is the only road (Z.3) |
| Face-1 §2.7 (translation-by-conversion) | demoted to **conditional** on `door_avg` (and possibly J-metric) |
| `door_avg` (Y.4) | **elevated to the single highest-priority kernel item** — four bars blocked behind it (Z.4) |
| J-metric face | **promoted** from fallback to registered v95 architecture hypothesis (space-direction transport partner of the pair) |

*End of resume update 3. The design's scoreboard, honestly: the
conversion physics (coherent door, no-cancellation, park approach,
anti-runaway) is winning; the transport physics (wave, translation,
phase retention) is losing with one common confound; the confound has a
named, cheap, still-unimplemented fix; and the propagation question is
now explicitly handed to the next face.*
