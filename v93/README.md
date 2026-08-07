# v93 — THE UNITARY DENSE CHANNEL (design document; docs-only, NOT opened)

**STATUS: a design document, pre-registered before any kernel edit. v92
remains the active programme; v93 is the pending next-version design that
the v92 consultation (`v92/consult/SUBQUARK_synthesis.md`) converged on.
Nothing here is code, a decision, or authorized for implementation — it is
the design the ratchet will test, awaiting explicit user sign-off.**

## The one-sentence thesis

Give the dense (matter) sector the field sector's transport algebra —
**within-mode dense transport becomes a product of unitary pairwise plane
rotations** (a cousin of the field's pass F), replacing the additive
magnitude want — so coherent dense transport (translation, the live fifth,
the far field) emerges with **exact conservation by construction**, and the
first moment that was missing (momentum) lives in the rotation cross-term
that the additive Em-ledger kept rejecting.

---

# PART I — BACKGROUND (everything needed to read this document standalone)

## I.1 The programme's one law and four constraints

The SCP substrate rests on one principle: **energy is never destroyed; it
only changes mode; space is one of the modes.** Therefore (v89/PRINCIPLE):
no background (nothing may persist and be merely re-valued); no imported
field or species; everything emergent from the fabric. Matter is converted
space; motion is regular conversion, not displacement; there are no
equations of motion in the usual sense because nothing has a trajectory.

## I.2 The substrate (v89+ free-cell, kernel `freecell.c`)

A dynamical cell roster (no fixed lattice — cells are born, die, jam, repack)
carrying three coupled sectors per cell:
- **field** `fa1, fa2` (a two-plane complex amplitude ψ_field; energy
  `Ee = fa1² + fa2² = |ψ_field|²`);
- **dense matter** `Em` (a scalar magnitude) + clock `th2` (a phase);
- **space** `Es` (the convertible mode; radii `r ∝ ∛Es`).

Plus the **conversion door** (field↔matter, pass 6): field condenses at
matter on a clock beat above a work-function threshold; matter evaporates
above cap; the click is quantized into whole **action atoms** ε = A₀ω/2π
(A₀ = `quant_A0` = 1.15, a fixed universal grain; ω = the source's local
load-detuned pitch; an integrate-and-fire with a two-atom credit).

The conserved charge is `Em + Ee + Es`; the battery's drift bars witness it
at the FP floor (≤1e-13).

## I.3 The asymmetry that is the open frontier

The two sectors transport energy by **different algebras**:

| | field sector | dense sector |
|---|---|---|
| state | complex amplitude `fa1,fa2` | scalar `Em` + phase `th2` |
| transport (pass F) | **product of plane rotations** (local precession + pairwise Givens hops, `freecell.c:1933`) — each exactly norm-preserving | **additive magnitude want** `swant`, debited/credited (`Em[src]−=d1; Em[rcv]+=d1`), phase used only to GATE |
| conservation | unitary to roundoff (rotations) | exact by ledger (debit=credit) |
| phase | carried by the amplitude, composes | cargo/gate only; not transported |

The field is unitary, carries phase, interferes, detects (Born ±2%). The
dense sector is **scalar**: it moves magnitudes, gates on phase, and
**does not coherently transport**.

## I.4 The measured absences (the motivation — every one downstream of §I.3)

All confirmed by experiment, all backtraced to the scalar dense transport:
- **B2 / nothing translates** (v_COE ≤ 1e-3 of closing; e3b drift a marginal,
  seed-variant residual ~0.002).
- **the fifth dies** (nv=6 boundary gate gg→0; the kr=1 negative template).
- **no far field** (the π-flatness theorem: closed-book matter loads the
  medium with identically zero; ASTRO π-flat ±1e-5).
- **P2 ~100× radiation-pressure deficit at the door.**
- **spin dies at the door** (QUENCH-3: τ_φ ≈ 3 t.u. flying; door writes
  phase at 3.5σ but retention zero).

`BLACKHOLE.md §0` pins all of these to one line of code: the conversion door
reads `Ee = fa1²+fa2²` and writes Em as a pure scalar — the first moment
dies at every mode boundary. The "amplitude completion" (make dense
transport coherent like the field) is the programme's central open problem.

## I.5 The v92 trail (what was tried, what failed)

v92 (Phase L) pursued the completion via the Phase-M **shadow amplitude**
(a per-slot complex meter, `sre_/sim_`, that rides the dense deposits and
measures ρ_coh). Three "feedback" driver forms were implemented and tested
(`v92/L1_FINDINGS.md`):
1. magnitude-bias `(c0−c1)/(c0+c1)` — blind to the tilt's phase gradient;
2. phase-current `Re(A·e^{−iθ_rcv})` — engages but perturbative;
3. zero-sum phase-current (FLOW architecture) — clean (fixed xsec
   over-absorption) but **still did not deliver L1-A** (e3b |cos|
   0.28→0.08, 0.94→0.83, 0.93→0.93).

The chord probe (Path B) was negative (the dead fifth has no deposits ⇒ no
shadow ⇒ no feedback ⇒ chicken-and-egg). A "Path A" design (native complex
transport via a per-cell received-amplitude R_j) was written and reviewed.

## I.6 The v92 consultation convergence (the basis of THIS design)

Three independent reviews (grok CLI, claude CLI fable/max, opencode) of
"break the dense cell into 2 sub-components (uptake/outtake), decay rates,
projective/toroidal/DEC geometry" **converged** (`v92/consult/`):
- the state layout (two floats vs `Em,th2`) is a **relabeling** — the kernel
  already stores the polar form + per-slot shadow;
- **what must change is the transport algebra** (unitary pairwise rotations
  replace the additive want);
- **conservation dissolves** (Givens hops never sum-then-square; the cross
  term `2·Re(ψᵢ*wψⱼ)` that broke additive bookkeeping IS the link current —
  it is where momentum lives);
- projective line = notation (drop); torus = parcel-level (acceptance
  language); **DEC = the implementation language**;
- **uptake/outtake = the current's moments + the doors, NOT two state
  components** (component-asymmetric coupling is a parametric pump —
  ignition-shaped — and the first step down the rejected quark road).

**The synthesis insight (claude, confirmed by all three):** the three
feedback forms did not fail by bad luck — they failed because **the additive
Em-ledger rejects the interference term.** The cross term surfaced as
throughput amplification (xsec 7.76), was renormed away (zero-sum), and that
renorm degraded the very signal wanted (e3b 0.28→0.08). *"Conservation must
become a theorem of the update (norm preservation), not a ledger patched
per form."*

---

# PART II — THE DESIGN

## II.1 The thesis

Within-mode dense transport is promoted from an additive magnitude want to a
**unitary channel**: a product of pairwise plane rotations on a complex
dense amplitude, exactly analogous to the field's pass F. The dense sector
finally obeys the v89 synthesis — *"amplitudes within a mode; atoms at mode
boundaries"* — which the gated magnitude want violates in spirit.

This is **not** "complexify Em" (a state relabel) and **not** "feedback the
shadow" (the closed hypothesis). It is a change to the **transport
primitive** itself.

## II.2 The state (enabling, not the point)

A complex dense amplitude per cell, `ψ_m = m1 + i·m2`, with `Em = |ψ_m|² =
m1²+m2²` and `th2 = arg(ψ_m)`. **Note carefully:** this is isomorphic to the
`(Em, th2)` the kernel already carries — the state layout is not the
innovation. The two components are a **U(1) frame** (they sweep into each
other under the local clock, exactly as fa1/fa2 do). They carry **no
separate physical jobs** (see safeguards §IV.1–§IV.2). The state change
exists only to *enable* the rotation algebra; the physics is in the hops.

## II.3 The transport algebra (the core change — pairwise unitary rotations)

Mirror pass F (`freecell.c:1933`). For each live slot `s` (link i–j), the
dense amplitude hops by a **pairwise Givens rotation** of angle `τ_s`:

```
m1[i]' = cc·m1[i] + ss·m2[j]
m2[i]' = cc·m2[i] - ss·m1[j]
m1[j]' = cc·m1[j] + ss·m2[i]
m2[j]' = cc·m2[j] - ss·m1[i]
```
with `cc=cos τ_s, ss=sin τ_s`. **This is exactly norm-preserving for the
pair** (`|ψ_i|²+|ψ_j|²` unchanged — the cross terms cancel), so `Σ|ψ_m|² =
Σ Em` is conserved to FP roundoff by the same mechanism that makes the field
unitary. No sum-then-square ever occurs; no renorm is needed.

The hop angle `τ_s` carries the dense physics that currently lives in the
want (chart order, closure, mobility, headroom), folded into a phase:
`τ_s ∝ base · g_ij · head · mobility · (chart factor)` — i.e. the existing
want-computation (`freecell.c:1347`) maps onto the hop angle rather than a
magnitude to debit. The **gate** (cos^p closure) survives as the hop-angle
envelope (it was already a rotation-ready real factor).

## II.4 Conservation (dissolves — the key result)

Under II.3, conservation is a **theorem of the update** (pairwise norm
preservation), not a patched ledger. The cross term `2·Re(ψ_i*·w·ψ_j)` that
broke additive bookkeeping — that surfaced as xsec over-absorption and had
to be renormed away — is no longer a bug: **it is the link current J_s, the
dense momentum.** The conservation wall (`|ΣA|² ≠ Σ|A|²`) that blocked Path
A is simply never encountered, because Givens hops never form `ΣA`.

## II.5 The DEC framing (the implementation language)

Discrete exterior calculus is the natural language and is closest to what
the kernel already does (claude's strongest point):
- cells = 0-cells; slots = 1-cells; `ψ_m` a 0-cochain; the in-flight
  amplitude / current `J` a 1-cochain on oriented slots;
- `d` (along-link difference) = want-formation (reading imbalance — the
  intake operator); `δ` (signed incident sum) = deposit-collection (the
  outtake operator); **J is the bilinear between them**;
- **cells update only by `δJ`** — then `Σ` over the complex telescopes to
  boundary terms and conservation is exact (Stokes). *The zero-sum renorm
  of v92 was Stokes' theorem derived by hand, once, for one form; DEC makes
  it structural.*
- the field's hop normalization `ŵ = w/√(s_i·s_j)` is already a discrete
  Hodge-star (mass matrix). A **dynamical** star — link weights evolving
  toward the flux pattern — is precisely FLOW's bed-digging law: *the bed
  law is the metric responding to the current*, the one place "space curves
  because energy is conserved and space is convertible" gets an
  implementable form.

DEC is combinatorial (incidence relations, no coordinates), so it is
background-free in exactly the programme's sense; a dynamical roster is its
home turf (cell birth/death = the complex changing; conservation degrades
gracefully to boundary flux at doors and roster events). **Take only the
minimal fragment** (0-forms, 1-forms, `d`, `δ`, door source terms); 2-forms
and discrete curvature are over-kit for L-1.

## II.6 The current & its moments — uptake/outtake done right

The link current `J_s = Im(ψ_i* · w · ψ_j)` (the rotation cross-term). Its
moments around a cell are the dense dynamical quantities:
- **monopole** (net `δJ`) = **metabolism** (eat vs shine — the scalar book
  balance the programme already measures at the radiance fixed point
  x̂*=0.62);
- **dipole** (first moment of J) = **momentum** — *"translation IS the
  current"* as a mechanism, not a slogan. A mover is a cell/parcel with
  balanced monopole and nonzero dipole: uptake fore, outtake aft, **a
  configuration of one current field**, not two state components.

**Uptake/outtake is the interpretation of J's moments and of the doors**, not
of the Re/Im pair: intake = the conversion/space-sector door (the g1
footprint side); outtake = the radiance channel; the in-mode amplitude
current runs between them. (Open-books gravity — spatially separated
intake/exhaust along a chain breaking π-flatness — is a *later face*, after
translation exists; see §II.9.)

## II.7 The door stays special — and the retention fix

The conversion door (atoms, cap, work-function, radiance) is **not
unitarized** — it remains the discrete, quantized mode boundary (the v89
"atoms at mode boundaries"). Two consequences:
- **conservation** holds across the door exactly as today (field loses d1,
  matter gains d1, space gains the footprint; `Em+Ee+Es` invariant).
- **retention (the risky piece):** the fired atom must carry `arg(ψ_m)` —
  the departure phase from the composed slot amplitude — **not `m·th2`** (the
  grok criterion-1 / D-3 lesson; QUENCH-3 showed the naked cell clock churns
  at ~0.2 rad/t.u. and qp_phase's retention was zero). The phase is carried
  slot-borne (protected from delivery churn). The field's own phase-faithful
  clicks (single-quantum clicks rebuild fringes) prove the pattern can work.

## II.8 Geometry — torus as acceptance language; projective dropped

- **Projective line `[m1:m2]`** = pure notation (RP¹ ≅ S¹ = the phase circle
  already stored). **Drop it.** (CP¹/Bloch would be new physics but needs
  four real components — a separate, larger campaign, not this one.)
- **Torus** = real content at the **parcel** level, not the cell. The
  programme's persistent matter *is* rings (nv=48/6/comp12); a ring carries
  winding around the ring-cycle × the clock-cycle, and a *moving* ring is a
  helix on that torus. Key theorem: **winding is topologically protected iff
  `|ψ|>0` along the cycle, and every conversion event is a phase-slip site**
  — which retrodicts both halves of the winding record (exact closure lives
  longest; winding dies in transit at the doors). Use the torus as the
  **acceptance language**: winding retention, phase-slip rates at doors,
  ring circulation quantized in steps of `2π/N` (pre-register: measurable).

## II.9 Anti-ignition — the honest risk, and the structural mitigation

**Unitarity kills the entire energy-runaway class outright** (transport
cannot amplify `Σ|ψ|²`). The remaining risk is **coherence-runaway** (grok's
named risk), armed by the bath's measured ρ_coh ≈ 0.77 (Phase M).

The structural mitigation: **rotations exchange amplitude reversibly** —
they have no fixed point at alignment. The programme's four measured
ignition machines (cantus v1/v1.1, registry, identity-maturity) are all
**dissipative Kuramoto-style phase pulls** (cf. `kappa_lock`,
`freecell.c:1533`); replacing pull with composition arguably *shrinks* the
ignition surface. But this is an argument, not a measurement. **The wall
sits at the conversion door** (two-atom credit, chart gate if needed); the
bars are the measured ones (bath byte-stats, glow/birth rates), staged per
grok. **Bar shift:** "nothing ever condenses" is no longer the right bar
(QUENCH measured lawful creation); the bar is *bath statistics unchanged at
inert, creation rate matching the measured quench law when live*.

## II.10 The decay-rate framing — Γ_p/Γ_E → ∞

"Matter can't translate" made exact (claude): the sector has **two decay
rates, maximally split**. Γ_E (structure loss) is small (blobs sit at
x̂*=0.62, burn down over hundreds of t.u.); Γ_p (momentum loss) is
effectively infinite (B2 ceiling ~5e-3, P2 ~100× deficit).

`Γ_p ≈ ν_door·(1−r)`, where ν_door is the conversion-event rate (∝
metabolism) and r is the per-event phase retention. QUENCH-3 measured
r ≈ 0 ⇒ today momentum dies within a few door cycles. **Pre-register:** under
the unitary channel + the arg(ψ) door, hotter objects lose momentum faster
unless r→1; the P2 deficit closes as r does.

---

# PART III — IMPLEMENTATION PLAN (pre-registered; nothing authorized)

## III.1 Knobs, defaults, byte-inertness

- `amp_nat` (default **0** = the additive want path, byte-identical to v92's
  V3a; >0 = the unitary dense channel engaged).
- `amp_nat` gates the **whole transport face swap** (rotation hops replace
  the additive want when >0). At 0, `swant` and the additive door run exactly
  as v92.
- `amp_tau` (existing) — the slot-borne amplitude window (retention).

The lane is **parallel and knob-gated**: the scalar want is *untouched* at
0, so byte-inertness against the 87-bar V3a surface is achievable.

## III.2 Acceptance bars (sharpened, staged)

- **L1-A (transport):** e3b |cos_to_kdir| → 1, speed ≥ 2.6e-3, seed-robust.
  **Staged** (per grok): median |cos| improves and seed variance *falls* at
  small `amp_nat` — not a cliff (two forms already "engaged but didn't
  deliver").
- **L1-B (conservation):** drift at floor (worst ≤ 1e-13) — the theorem.
- **L1-C (anti-ignition):** bath byte-stats + glow/birth rates **unchanged**
  at `amp_nat>0` vs `=0`; the armed ρ_coh≈0.77 test. If it ignites, the
  design records a strike and returns to the design round.
- **L2 (chord, maintenance):** nv=6 boundary fifth gg ≥ 0.9 at T=100 (a
  **live** fifth sustained; sustain ≠ ignite — nucleation is separate).
- **Spin/winding retention:** the m=+2 texture outlives its driving field
  (QUENCH-3 rerun on the unitary substrate); phase-slip rate at doors
  measured.
- **Invariant surface:** the 87 V3a bars; C≡Go.

## III.3 Engineering guards

- **Byte-inertness:** `amp_nat=0` reproduces v92's V3a byte-for-byte (header
  line only); battery ALL GREEN 87.
- **C≡Go at strength:** compounded rotations **amplify** the 1-ulp C/Go FP
  divergence the v92 feedback already exposed. Expect to move the abx bars
  to a **tolerance-based** comparison (relative, not byte) for production
  sweeps; the physics-state columns stay the identity carrier.
- **The comb gate — derive-then-retire:** whether the consonance comb
  (p·q≤6) can migrate from links to doors without losing the interval
  hierarchy should follow the κ_freq precedent (imposed → derived as the
  interference cross-flow → retired). Keep the comb gate at the door
  initially; add a bar that phase-matched selectivity reproduces its measured
  selectivity; retire it only by measurement. It *may* be derivable from
  composition (the deepest vindication) — earn it, don't assume it.

## III.4 Order of work (after sign-off)

1. Implement the unitary dense hop (pass F's cousin) behind `amp_nat` in C;
   verify byte-inertness + the conservation theorem (drift floor) + the e3b
   (now working).
2. Mirror in Go; C≡Go (tolerance abx).
3. The `arg(ψ)` door (retention); QUENCH-3 winding rerun.
4. Full acceptance (L1-A, L1-B, L1-C, L2 maintenance, invariant surface).
5. DEC refactor (J as 1-cochain, `δJ` updates) — the language cleanup, once
   the physics lands.

---

# PART IV — SAFEGUARDS: WHAT NOT TO DO, AND HOW NOT TO IMPLEMENT

These are hard prohibitions, each tied to a measured failure or a permanent
project rule. **Read before any code.**

## IV.1 Do NOT make it a state relabeling
Two floats `(m1,m2)` in place of `(Em,th2)` is notation — the kernel already
stores the polar form (`freecell.c:903`) and the per-slot shadow. If the
only change is the state layout, you have built nothing. **The change must be
to the transport primitive (the hops), not the state container.**

## IV.2 Do NOT assign physical jobs to m1 vs m2
The two components are a rotating U(1) frame (they sweep into each other
under the local clock). A law that reads m1 *as* intake samples a rotating
frame — a **parametric pump** (lock-in) that rewards whatever is phase-
aligned, i.e. it hands the m=1 bath (ρ_coh≈0.77) a rectifier → **ignition-
shaped**. And per-component dynamics with separate identities is exactly the
**rejected quark road** (the pre-v89 multi-component lattice on a
background). **Uptake/outtake belongs to the current's moments and the doors,
never to the two state halves.**

## IV.3 Do NOT unitarize the door
The conversion door (atoms ε = A₀ω/2π, cap, work-function, radiance) stays
discrete and quantized — the "atoms at mode boundaries." Folding the door
into a continuous rotation destroys the click grain, PAULI-0, and the XSEC
structure. **The door is a source term in the DEC picture, not a hop.**

## IV.4 Do NOT sum-then-square (the |ΣA|² ≠ Σ|A|² trap)
Never compute a delivered energy as `|Σ A_k|²` from a coherent sum of
arrivals — that breaks conservation (cross terms). The whole point of
pairwise Givens hops is that **you never form ΣA**; each hop conserves the
two-cell norm. If a design step requires summing amplitudes before squaring,
it is wrong.

## IV.5 Do NOT use a feedback / bias / shadow-driven correction
The feedback hypothesis is **closed** (three forms failed; `L1_FINDINGS`).
Any "bias the magnitude want by the shadow" is form-3 and will fail the same
way (the additive ledger rejects the interference term). The unitary channel
is a transport-primitive replacement, not a correction on the existing want.

## IV.6 Do NOT write the cell clock th2 at the door
QUENCH-3 measured the naked cell clock churns (~0.2 rad/t.u. differential
rotation); `qp_phase`'s door-write had retention zero. **Phase is carried
slot-borne** (protected from delivery churn), and the fired atom carries
`arg(ψ)` — not `m·th2`.

## IV.7 Do NOT reintroduce a background
The pre-v89 multi-component "quark" fields were rejected because they lived
on a fixed lattice. A fixed two-component **fiber** per cell is acceptable
**only as mode structure that dies with the cell and indexes nothing
immortal** (like the field's fa1/fa2). The quark line is crossed the moment
the components acquire **separate ledgers, separate clocks, or an identity
that outlives the cell**. Guards: one norm, one clock, U(1)-covariant
dynamics, no per-component observable in the battery.

## IV.8 Do NOT try to revive a dead fifth
Sustain ≠ ignite (the chord probe showed chicken-and-egg: dead channel ⇒ no
deposits ⇒ no current). The unitary channel fixes the **carrier** half
(rotations move amplitude regardless of a live-gate prerequisite), but a
standing fifth still needs a **seed**. L2 acceptance is **maintenance of a
live fifth**; nucleation is a separate QUENCH-style experiment.

## IV.9 Do NOT soften bars to pass (the ratchet)
Every moved bar gets an explicit sharpen/re-gauge/reject decision with the
user; bars sharpen by measurement, never soften to pass; a green test leaves
the suite only by explicit user decision.

## IV.10 Do NOT modify the kernel without explicit user sign-off
CLAUDE.md's kernel-protection rule is permanent. This design is a **major**
kernel change (the dense transport primitive, in both C and Go, full
conservation re-derived). No code until the user explicitly authorizes v93
implementation.

## IV.11 Do NOT pre-claim the comb is derivable
It *may* be that the consonance comb emerges from composition (the deepest
vindication) — but that is a result to earn by measurement, not an
assumption to build on. Gate derive-then-retire; do not assert derivation
before the bar lands.

## IV.12 Do NOT expect C≡Go byte-identity at strength
Compounded rotations amplify the 1-ulp C/Go FP divergence. Build the
tolerance-based abx from the start; do not chase byte-identity by distorting
the rotation algebra.

---

# PART V — SIGN-OFF GATE

This document is the design. The path to code is:
1. User reviews this design (and `v92/consult/`).
2. User **explicitly authorizes v93 implementation** (a major kernel change).
3. STEP ZERO (v93 tree carry from v92 + baseline battery) → then the unitary
   dense hop behind `amp_nat`, face by face, each battery-gated.

**Nothing in v93 is opened, adopted, or implemented. v92 (laws_V3a + the
safe L-1 apparatus) remains the active, committed substrate.**
