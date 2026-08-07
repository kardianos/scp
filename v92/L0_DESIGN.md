# L-0 — the Phase L design round (v92; registered prerequisite before any L-1 code)

**Phase M closed with: "Phase L proceeds only from a design round over this
map."** This is that round. It specifies the L-1/L-2/L-3 law faces — knob
names, inert defaults, the update rule at the code-structure level, the
anti-ignition gate, the protected-carriage mechanism, the sharpened bars,
and the invariant surface — PRE-REGISTERED before any kernel edit. L-1 is
implemented first (smallest measurable); L-2 and L-3 land face-by-face, each
battery-gated. Nothing here is code or a decision; it is the design the
ratchet will test. **L-1 implementation is a kernel modification and awaits
explicit user sign-off.**

## 0. The design inputs (the five registered constraints)

1. **The coherence map (AMPLITUDE Phase M, `v91/AMPLITUDE.md` §4.1).** The
   shadow amplitude's delivery coherence ρ_coh = |ΣdA|/Σ|dA| per class:
   warm bath 0.641/0.766/0.853/0.910 (q25/50/75/90); frozen-chord bonds
   0.436/0.495/0.495/0.495; tilted blob 0.651/0.756/0.840/0.903. Two
   inversions of the registered predictions: the bath is NOT incoherent
   (unison 1:1 composes at 0.77), and the deficit is on the CHORD'S OWN
   BONDS (the m=2 fifth composes at half strength). Coherence deficit
   scales with CHART ORDER — Phase L would act hardest where the flux
   machine lives.
2. **The ignition asymmetry (registered BEFORE any law code).** An
   |A|²-driven transport law rewards the bath's coherent unison flows more
   than the chord's half-coherent fifth — ignition-shaped. CANTUS measured
   it twice: any closure-gated self-growing coherence ignites a bath-wide
   Kuramoto transition (glow −10× to −50×). **Phase L must gate on
   something the bath's unison cannot buy.** The candidate: chart order.
3. **The QUENCH-3 requirements card (`v91/QUENCH.md` §7.7).** The door CAN
   carry phase (3.5σ imprint at qp_phase=1) but the naked cell clock loses
   it in ~10–20 t.u. (differential rotation ~0.2 rad/t.u.). What is missing
   is phase PROTECTED from delivery churn: carried on the parcels/slots
   (delivery coherence measured 0.77–0.95), not exposed in the cell clock.
   Requirements: survive ~0.2 rad/t.u. differential rotation; do not
   inherit speckle deliveries; drive on m≥2 charts.
4. **The kr=1 negative template (`carried/S2_CHANNEL.md` §1).** Any design
   that kills the nv=6 boundary fifth (gg→0) or bloats the D-ring past
   x̄ 0.4 fails. The raw rate-level reactive term (−κ·g·sin ψ) is the
   wrong-level object: a unison instrument that suppresses the fifth
   wherever ψ oscillates.
5. **The ignition ledger.** Four statistical self-growth laws ignited
   (cantus v1 cell-borne, cantus v1.1 link-borne, registry, identity-
   maturity). ONE structural law did not: FLOW's bed, whose zero-sum
   per-cell budget starves symmetric churn by construction. **The
   architecture that survives is structural (zero-sum), not statistical.**

## 1. The apparatus L-1 inherits (already in the kernel, as the Phase-M meter)

The shadow amplitude is ALREADY built (`kernel/freecell.c`, guarded by
`amp_tau`, default 0 = byte-inert, physics-silent):

- **per slot-direction** (the want's own granularity): at deposit
  `sre_[2s+dir] += √f·cos(m·θ_src)`, `sim_ += √f·sin(m·θ_src)`, where
  `m = slq[s]` or `slp[s]` — the link's **chart order** (freecell.c:1428).
  The phase is folded by the chart. A unison link folds by m=1; a 3:2
  fifth by m=2/3.
- **in transit:** advanced by `w·d/C` (phase rotation), composed into the
  dst chart frame at arrival (freecell.c:1490).
- **per cell:** a shadow account decays on the par-style window
  (`amp_tau`); the meter prints ρ_coh per class.

So the dense transport (`swant[2s+dir]`, the magnitude want) ALREADY moves a
complex amplitude (the shadow) alongside the magnitude — chart-folded,
path-composed, cell-accumulated. **L-1 promotes this shadow from meter to
driver.** The chart-order structure (`slq`/`slp`), the path composition, and
the slot-borne carriage are not new design — they are the Phase-M apparatus
switched from read to write.

## 2. The three design decisions (resolved here from the evidence; flagged for review)

**D-1: deep completion, not shallow modulation.** A magnitude boost
(`swant *= 1 + k·|A|²`) is shallow — it cannot deliver the fifth-as-amplitude
(L-2) or the exchange sign (L-3), and the kr=1 sweep shows rate-level
magnitude pulls kill the fifth. The design memo's "translation IS the
current" demands the dense transfer CARRY phase that composes. **Decision:
the shadow amplitude drives a directional transport current (redistribution
toward the coherent direction), not a scalar boost.** This is the registered
S2 §2 item 1.

**D-2: the anti-ignition gate is chart order m≥2, AND the feedback is
zero-sum per cell.** Two independent ignition failure modes must be blocked:
(a) rewarding bath unison — blocked by gating on `slq[s]≥2 ∨ slp[s]≥2`
(the link is a chord, not unison); the bath's churn is overwhelmingly
m=1, so it is structurally untouched. (b) runaway growth — blocked by
making the feedback a REDISTRIBUTION of each cell's existing transport
budget (zero-sum per cell, renormalized), exactly FLOW's structural
anti-ignition (ignition-ledger item 5). **Both gates are structural, not
statistical** — the lesson of the four ignitions. Decision: chart-order
gate (D-2a) AND zero-sum renormalization (D-2b), together.

**D-3: protected carriage is slot-borne (inherit the shadow).** The QUENCH-3
requirements card says phase must survive ~0.2 rad/t.u. differential cell
rotation. The shadow IS slot-borne (it lives on `sre_/sim_[2s+dir]`, not
`th2[i]`), so it is already protected from the cell's own delivery churn by
construction — the cell clock can rotate freely under the slot's amplitude.
Decision: L-1 reads and feeds back the slot-borne shadow; it does NOT write
the cell clock (qp_phase's mistake — a driven texture, retention zero).

## 3. L-1 — the transport face (momentum first; smallest measurable)

**Thesis.** Where the shadow amplitude composes coherently on a chord
(m≥2) link, it biases the dense transport along the coherent direction —
"translation IS the current." A tilted blob's forward window acquires a
coherent amplitude current and translates (the e3b door). This is the
missing first moment (BLACKHOLE §0): momentum as the first moment of
conversion, carried by the amplitude current.

**Knobs (default-inert, purity-pinned, battery-gated, both kernels).**
- `amp_drv` (default **0** = byte-inert; the Phase-M shadow stays a meter).
  Strength of the amplitude-driven transport redistribution.
- `amp_tau` (default **0** stays; >0 = the shadow window, already wired).
  L-1 requires `amp_tau > 0` to have a shadow to read.
- `amp_mmin` (default **2**; the chart-order gate floor). Amplitude feedback
  fires only on links with `max(slq,slp) ≥ amp_mmin`. 2 = fifths and above;
  unison (m=1) untouched.

**Update rule (pass 2, where wants `w_ij`/`w_ji` are set; freecell.c:1347).**
For each live slot `s` with chart order `m = max(slq[s],slp[s])`:
```
if amp_drv > 0 and amp_tau > 0 and m >= amp_mmin:
    read the slot's coherent shadow projection  c_s = sre·û + sim·v̂
        (the shadow amplitude projected on the link unit vector)
    bias  b_s = amp_drv · c_s / (1 + |c_s|)        (saturating, signed)
    w_ij *= (1 + b_s);  w_ji *= (1 - b_s)          (redistribute, not create)
then per cell i renormalize ZERO-SUM:
    W_i = Σ_incident w;  f_i = (budget_i) / W_i;  w_s *= √(f_i·f_j)
```
The bias is SIGNED (redirects transport toward the coherent direction),
SATURATING (bounded), and the per-cell renormalization keeps `Σ_incident w`
constant (FLOW's architecture — no global quantity to grow). Energy is
redistributed, never created: `Σ w_ij` per cell is invariant ⇒ `drift`
untouched by construction (mirrors the FLOW P-F0 proof).

**Acceptance bars (sharpened; pre-registered).**
- **L1-A (transport):** e3b tilted blob forward speed ≥ 2.6e-3, cos ≥ 0.95,
  seed-robust (≥3 seeds) — structure TRANSLATES. (V2g/V3a: marginal,
  ±1e-3 residual; the S2 acceptance.)
- **L1-B (conservation):** full battery drift at floor (worst ≤ 1e-13); the
  zero-sum budget holds; `Σ w` per cell invariant (printed witness).
- **L1-C (anti-ignition, the critical bar):** warm bath ρ_coh distribution
  AND the V3a glow/birth-rate UNCHANGED at `amp_drv>0` vs `=0` (the
  chart-order gate + zero-sum starve the bath). If this fails the design
  records a strike and returns to L-0 — the fifth ignition lesson.
- **L1-D (P2, the S2-full criterion, travels with L1):** with translation
  as the current, the ~100× radiation-pressure deficit at the conversion
  door closes (P2). Measured after L1-A lands.
- **L1-E (byte-inertness):** `amp_drv=0` reproduces V3a byte-for-byte
  (header line only); battery ALL GREEN 87; C≡Go byte-equal firing.

## 4. L-2 — the chord-chart face (after L-1; the flux-machine interior)

**Thesis.** The p:q chart transports amplitude with the interval's phase
map (θ → (q/p)θ). The shadow ALREADY folds phase by `m` (freecell.c:1428)
and composes into the dst frame — L-2 makes this composition LOAD-BEARING
for the chord: the boundary fifth survives coherence, the UUD triad holds
its circulation, a flavored composite RADIATES.

**Knobs.** `amp_chart` (default 0; the chart-composition feedback strength).
L-2 uses the same shadow but reads its CHART-FRAMED composition (the
arrival-side rotation at freecell.c:1490, currently meter-only).

**Acceptance bars.**
- nv=6 boundary fifth `gg ≥ 0.9` at T=100 with D-ring `x̄ ≤ 0.35` (the kr=1
  negative template inverted — the fifth survives).
- UUD `ggm ≥ 0.9` (the flux-machine circulation holds per link).
- a flavored composite RADIATES (rough > 0 sourced at the boundary fifth).
- the comb still admits only chords at x̂*=0.62 (the geometry wall holds).

## 5. L-3 — the exchange face (after L-2; the exclusion precursor)

**Thesis.** Two-quantum amplitudes with a sign under exchange — the HOM
template (v89 measured g_b 0.42 < 0.5 < g_f 0.58). The shadow, carried per
parcel (IDENTITY's slot-borne gid), gains an exchange phase; identical
arrivals refuse one state in occupancy (the B5 precursor), distinct from
cap saturation (measured pitch-blind).

**Knobs.** `amp_exch` (default 0; the two-quantum exchange-sign strength).

**Acceptance bars.** Y-junction occupancy statistics diverge: identical vs
fifth-consonant arrivals separate beyond the cap-saturation control.

## 6. The invariant surface (must NOT change — the S2 §2 item 5)

The atoms machinery (grains ε = A₀ω/2π, two-atom credit, ħ-linearity
~1e-8), the conversion door (cap, headroom, evap — PAULI-0 and XSEC), the
space sector (pressure pushes, footprint g1, the skirt), the contact rule,
and the adopted V3a law (radiance + workfn). **The 87 V3a bars are the
invariant surface** — L-1 lands as a new mode of the existing want channel,
gated by the same ratchet; any bar that moves gets an explicit decision.
The kernel's energy conservation is untouched by construction (zero-sum
redistribution ⇒ drift floor).

## 7. Pre-registered predictions (L-1, before any run)

- e3b forward speed rises above 2.6e-3 with cos ≥ 0.95 (the current is
  carried by the coherent amplitude; the ±1e-3 residual is the incoherent
  part L-1 coheres).
- the bath is UNCHANGED at amp_drv>0 (chart-order gate + zero-sum) — the
  critical anti-ignition prediction. registered BEFORE the run.
- P2 (radiation pressure) does NOT close at L-1 alone — it is the L1-D
  criterion that travels with translation; the ~100× deficit closes only
  when the current reaches the door. Registered as a deferred measurement.
- drift stays at the floor (zero-sum ⇒ no energy touched).

## 8. Open design questions (for the user, before L-1 sign-off)

1. **D-1 confirmation:** deep (directional current) over shallow (scalar
   boost). Recommended: deep — the evidence (kr=1 kills the fifth; L-2/L-3
   need phase) demands it. The shallow form is a fallback if L-1 shows the
   current is unstable.
2. **The zero-sum budget class** for L-1's renormalization: FLOW uses a
   per-cell geometric-mean renorm of conductance. L-1's want-redistribution
   can use the same form (proven ignition-immune) or a simpler per-cell
   sum-normalization. Recommended: inherit FLOW's form (defense-in-depth
   via the one survivor).
3. **L-1 scope:** transport-face only at first (e3b + anti-ignition +
   byte-inertness), OR bundle the P2 re-measure (L1-D) into the first run.
   Recommended: transport + anti-ignition + inertness first (the smallest
   measurable unit); P2 re-measure as the immediate sequel once L1-A lands.
4. **The shadow's fidelity:** Phase M measured ρ_coh at the existing
   shadow resolution. L-1 reads the same shadow. If L1-A lands marginal,
   a higher-fidelity shadow (more phase bits / tighter window) is the
   registered escalation — but start with what is built.

## 9. Order of work (after sign-off)

L-1 (transport) → L1-A acceptance → L1-D (P2) → winding-retention rerun
(QUENCH-3's standing vortex on the L-1 substrate: spin becomes carried or
back to L-0) → L-2 (chord charts) → open-books + beds (gravity) → L-3
(exchange). Each face battery-gated, inert-default, C≡Go re-verified, the
ratchet governing.
