# PLASTICITY — the foam anneals (P15)

**User proposition (2026-07-28):** the FIXED jittered foam is itself what
prevents consonant structures. Physically the cells are plastic — under
tension and harmonic misalignment they would realign their harmonic planes
and morph their shells until hard. On a frozen foam the slightest
misalignment destroys a structure; real cells would constantly realign BY
FORCE.

**Verdict up front:** the proposition survives quantification and comes out
sharper. The frozen foam denies exactly the one thing a consonant shell
needs (one shared strut length); the misalignment force is computable from
quantities the kernel already evaluates every step (∂ψ/∂d = −ω/C — the
misfit gradient is the reactive term S2 already derived, acting on its
OTHER conjugate variable); the flow is conservation-exact because geometry
is not a ledger; and in the reduced models it converges to exact consonance
on bipartite shells in ~20 t.u. (80× faster than death), pays the comma in
geometry instead of energy, silences parasites by parity, and amputates
frustrated links instead of tempering them. Recommended: option (i)
retardation plasticity with option (iv) hardening as its second phase; one
new constant (κ_plast, a soft mobility with a ≥50× stability window).

Tags as in CONSONANCE.md: [D] derived, [P] postulated, [G] guess,
[M] measured (here: measured in the reduced numpy models or in the foam
geometry, never in the simulator — no kernel run was made for this
document). Models: `prestress/theory/plast_foam.py`, `plast_pair.py`,
`plast_net.py`, `plast_net2.py`. Constants throughout: `battery/laws_V2g.cfg`.

---

## 0. The kernel audit — what geometry IS in cellfab.c [M, code]

Read before designing; every claim checked against `v89/cellfab.c`.

### 0.a Positions are scaffold. Dynamics reads only the link table.

`cx/cy/cz` are declared "scaffold positions (init + diagnostics only)" and
the code honors it: no dynamics pass reads a position. Build uses them once
to dart-throw cells and construct the link table; seeders use them to pick
cells; every diagnostic (centroids, shells, census, RAD, CLICK coordinates)
uses them as reconstruction. The physics of every step reads only:

| array | meaning | mutability today |
|---|---|---|
| `li/lj` | link endpoints (topology) | frozen at build |
| `ld` | link length = retardation ω·d/C, flight time d/C, conductance 1/d | **frozen at build** |
| `lux/luy/luz` | link direction cosines (axis collimation, plane alignment) | frozen at build |
| `lA` | live lens-overlap area A(d, cr_i, cr_j) | **recomputed every step** (pass 2) |
| `cr` | live radius = cr0·(Es/e_s0)^{1/3} | **recomputed every step** (pass 1) |

**Consequence [D]:** "plasticity" is precisely *make d_ij respond to
physics* — a link property, not node motion. There is no background to
move nodes in, and none is needed: d IS the physical content (retardation,
latency, conductance); the embedding is bookkeeping. v89 already lives
this way — it just froze d at the dart-throw accident.

### 0.b Radii: the geometry is already HALF-alive

`cr` already evolves: pass 1 sets `cr = cr0·cbrt(Es/e_s0)` from the live
space store, and `lA` follows every step through the lens formula. So the
space law already morphs channel AREAS (conductance, liveness/pinch-off,
strangulation — a measured effect, CONSONANCE III.3) while channel LENGTHS
stay frozen. **The current kernel is asymmetric with no principled
defense: Es re-sizes the contact but not the retardation.** Option (ii)
below is the repair of that asymmetry; option (i) is the direct force.

The candidate-link set is built once from `cr0` with ceiling
1.15·(r_i+r_j). Candidates can never be added; but *liveness is already
dynamic*: lA hits 0 when d ≥ cr_i+cr_j and can return. On the standard
foam [M]: 85,884 candidates, 57,688 live at vacuum (67%). Dormant
candidates at d = 1.10·(r0_i+r0_j) wake if Es grows 1.33×; the ceiling
link needs 1.52×. So option (i) inherits a bounded, rebuild-free topology
plasticity for free: d moving across the pinch is a link dying or being
born *within the closed candidate universe*.

### 0.c Every runtime use of `ld` (the response spec for any d-dynamics)

1. **Dense gate retardation** (pass 2): ψ_f = qθ_i − qω_i d/C − pθ_j and
   ψ_b likewise; also the κ_reac arguments (inactive at kr=0).
2. **Geometric conductance** (pass 2): geo = (A/Aref)·(dref/d).
3. **Lens area** lA(d, cr): pass 2 (and init) — liveness itself.
4. **Space transport** (pass S): w = (lA/Aref)·(dref/ld).
5. **Deposit entrainment target** (pass 4): err = wrap(m_s(θ_src −
   ω_src d/C) − m_r θ_r).
6. **Flight advance** (pass 5): adv = dt·C/ld — the cycle takes d/C; plus
   the delivery entrainment err.
7. **Field sector** (pass F): fsum weights and hop τ, both (lA/Aref)(dref/ld).
8. **Diagnostics**: pair δ and act = fl·2d/C, DBG gates, RAD field current,
   ladder mode; lump connectivity via lA.
9. **Seeders**: all phase seeding via ω·d/C (build-time only).
10. **The action atom**: A0 auto mode uses build-time d̄ once. In
    laws_V2g `quant_A0 = 1.15` is a pinned constant. **Spec: A0 must never
    track live d̄** — the atom is law, not geometry; annealing must not
    drift ħ_eff. (With the pinned value this is already true; the `-1`
    auto seam must stay a build-time evaluation.)

Every one of 1–8 reads `ld[l]` fresh each step. **Therefore a single
dynamical `ld` propagates coherently through retardation, conductance,
liveness, flight time, space transport and field hops with no further
code paths** — the kernel is already shaped for it. Directions `lux` stay
frozen (they encode relative channel orientation; nothing in the
no-background reading requires them to move, and the alignment pass
already adapts the *planes* to them). In-flight energy is unaffected by
d changing under it: `lem` is a ledger, `lph` a cycle fraction; a length
change alters remaining flight time, never energy. Conservation-neutral.

### 0.d The adaptation channels that already exist

| channel | what it adapts | driver | rate |
|---|---|---|---|
| `kappa_align` (pass 7) | plane NORMALS n1/n2 — rotates the planes' intersection line onto the transfer direction, sign-matched to partners | flux-weighted (lflux) | 0.5, ~few t.u. under traffic |
| `kappa_lock` (passes 4/5) | clock phases toward the retarded tail | per-deposit mix f/(f+mob+lockf) | locks in ~10 t.u. [M, e6/e7] |
| `lp/lq` pair temper | the link's interval identity p:q | comb best-score each step | instant (no hysteresis — the open "interval identity" item) |
| load detune | pitch ω = w/(1+q·x) | energy exchange (F1: one-sided walk) | 20–80 t.u. [M, e8] |
| `cr` ← Es | contact areas, liveness | space transport s_k + s_pull | s_k equalization ~40 t.u. |
| `sigma_tumble` | normal diffusion (noise) | — | 0.01 |
| `rough_k` ledger | (not adaptive — the misfit ENERGY leaving) | off-comb deliveries | 52% of a blob's leak [M, M0] |

**The user's proposition, mapped [D]:** "realign their harmonic planes" —
the orientation half EXISTS (kappa_align; planes follow traffic); the
phase half EXISTS (kappa_lock); the pitch half EXISTS one-sidedly (load
walk). "Morph their shells" — the radius half EXISTS (cr ← Es). **The
missing degree of freedom is exactly one: the link length.** Every
adaptive channel the kernel owns steers TOWARD consonance and one frozen
number (d, set by a dart-throw at build) vetoes it. That is the
obstruction, named.

---

## 1. Diagnosis — the frozen foam is the obstruction, quantified [M]

From the real campaign foam (`foam/foam_s20260727.tsv`, 9741 cells,
d̄ = 1.505; `plast_foam.py`):

**What consonance requires.** A bipartite shell of unison voices needs
every strut on the π-rung at ONE shared length d* = πC/ω (LEDGER,
translation table). The gate cost of missing it, at p_gate = 8, is
G²(δ) = [(1+cos(δ/2))/2]^16 with δ = 2ω·(d−d*)/C; small-δ:
**1 − G² ≈ δ² ≈ (2π·spread)²** — churn is quadratic in the length spread:

| edge spread | δ (rad) | G² |
|---|---|---|
| 1% | 0.063 | 0.996 |
| 2% | 0.126 | 0.984 |
| 5% | 0.314 | 0.906 |
| 10% | 0.628 | 0.673 |
| 15% | 0.942 | 0.408 |

**What the foam supplies.** Live links: fractional spread σ_d/d̄ = 16%.
Cube picks with the kernel's own seeder algorithm (greedy + 3 refinement
passes), 40 placements at a = 1.5: median edge spread **16.6%** (best
placement 10.3%), median per-edge G² ≈ 0.44, median MIN edge gate 0.02.
The center-of-box picks that H1 actually used: a=1.25 → spread 20.3%,
a=1.5 → 18.9%, a=1.5 core-excluded → 15.5% with max |δ| = 2.3 rad. BFS
phase seeding concentrates cycle defects on co-tree edges: modeled seed
gates mean 0.63–0.76, min 0.002–0.016 — **reproducing the measured H1
seed report (mean ≈ 0.6, min ≈ 0) from pure geometry.** The foam also
inserts parasites (median 6 in-ceiling face diagonals per placement) and,
at a = 1.5, throws individual edges beyond the pinch (up to 1.95 > 1.7):
some "edges" of a seeded cube are not live channels at all.

**Refinement cannot fix it.** The 3-pass edge-uniformity swap search —
already in the kernel — bottoms out at 10–15% spread. The foam's blue
noise (dart throw at dmin = 1.0) simply does not contain 12 equal lengths
meeting at 8 specific cells. Form-finding SELECTS from what exists; at
cube scale what exists is a 0.4–0.6-gate skin. This is the measured H-series
wall: "no jittered mini-cube achieves it (gates mean 0.6)" — and the
H1-v3 verdict says the skin must be CONSONANT because off-rung skins
radiate (roughness = 52% of the M0 leak; the leak is the killer via the
skirt road, death at t ≈ 1600–1900 for light rings, 2200–3800 heavy).

**What plasticity must deliver.** Per-edge corrections |Δd| = |d − d̄_set|:
mean ≈ 0.20, max ≈ 0.45–0.66 (13–30% of the strut) — an order of
magnitude beyond what pick refinement achieves — on a timescale long
against phase lock (~10 t.u., so it does not fight entrainment) and short
against death (≪ 1600): a target window τ_plast ≈ 20–300 t.u. To reach
gate 0.98 the residual spread must anneal to ≤ 2%: an 8× tightening.
Halving the spread quarters the churn; annealing 16% → 2% cuts the
structural roughness channel ~64×.

---

## 2. The four formalizations

### The rank trilemma [D]

The misfit lives on links: E gate residuals. The knobs decide the option:

| option | knobs | rank vs E misfits | background? |
|---|---|---|---|
| (i) d per link | E | **full rank** — every residual individually tunable | none (d is relational) |
| (ii) d from cell Es | V (cells) | rank-deficient (cube: 8 knobs / 12 misfits; foam: V ≪ E = 8.8V) | none |
| (iii) node motion | 3V | cube: enough (18 ≥ 12); foam bulk: deficient (3V < E); and d constrained to Euclidean-embeddable vectors | **yes — promotes the scaffold to physics** |
| (iv) hardening | rate modifier on (i) | — | none |

Only (i) has full rank without re-importing the embedding. (iii) is
doubly wrong: fewer knobs than misfits on the bulk foam AND the reachable
d-vectors are exactly the "triangle-inequality" cone of a background v89
says does not exist. **v89 has no background to embed in — d_ij is a
retardation, and unconstrained-by-embedding is the honest geometry.**
This is the conceptual point of the whole campaign: frustration that is
GEOMETRIC in a Euclidean foam (can't have 12 equal edges there) is only
COMBINATORIAL here (cycle parity, §3.3). Plasticity removes the false
obstruction and leaves the true one — the integers.

### 2.1 Option (i) — RETARDATION PLASTICITY [recommended]

**Law [P].** Each link's length is a dynamical property driven by the
harmonic misfit of the channel:

```
ḋ_l = − κ_plast · Φ_l · ∂V_l/∂d_l

V_l   = 1 − G(ψ_f) G(ψ_b)                    (the both-gate misfit, ∈[0,1])
ψ_f   = wrap(q θ_i − q ω_i d/C − p θ_j)      (kernel pass-2 quantities)
ψ_b   = wrap(p θ_j − p ω_j d/C − q θ_i)

∂ψ_f/∂d = − q ω_i / C ,   ∂ψ_b/∂d = − p ω_j / C          [D]
⇒ ∂V/∂d = (1/C) [ q ω_i G'(ψ_f) G(ψ_b) + p ω_j G(ψ_f) G'(ψ_b) ]
G'(ψ) = − (p_gate/2) h^{p_gate−1} sin ψ ,  h = (1+cos ψ)/2

Φ_l   = k_dep·k_dep_m · geo(d) · gpl · res · sqrt(m_i^eff · m_j^eff)
```

Φ_l is the **gate-free sympathetic urge** — literally `base·Sm`, the
prefactor the kernel already computes for the S2 reactive term (pass 2,
kappa_reac block). At entrained equilibrium (ψ_f = ψ_b = −δ/2, the comma
tempered evenly) the flow reduces to the rung-offset form [D]:

```
ḋ = − κ_plast · Φ · (p_gate/2) · h(δ/2)^{2p−1} · sin(δ/2) · ω_sum/C
```

Restoring across the whole inter-rung interval: δ > 0 (link too long)
⇒ ḋ < 0; δ < 0 ⇒ ḋ > 0; unstable only at the anti-rung δ = ±π (basin
boundary; §3.5). Numerical guard: floor d ≥ 0.5 (adv = dt·C/d and geo ∝
1/d stay bounded); no ceiling needed — the pinch d ≥ cr_i+cr_j closes lA,
which zeroes geo hence Φ hence the force: **a link that strangles itself
freezes sealed** (existing law does the sealing; plasticity only walks it
there).

**The S2 connection [D].** CONSONANCE VII.1: at rate level the gates are
even, g_ij ≡ g_ji on a rung — no odd-in-detune flow exists; the restoring
mechanism is interference energy, and its raw rate-level form (kappa_reac
= 1) breaks tilted transport (16/17, the standing S2-full criterion).
The plastic force is that same odd, sin(ψ)-class term acting on its
*other* conjugate variable: the phase ψ is conjugate both to the flows
(where rate compression cannot host it) and to the length d (where
nothing else lives). The homeless derived term finds a home that touches
no ledger. κ_plast is then not a new coupling but a new *mobility*: how
fast geometry yields to a stress the theory already derived.

**Conservation — where does the work go? [P]** Nowhere, and that is the
principled answer, not an evasion. d is not an energy ledger; the kernel
already moves geometry for free every step (cr from Es, normals from
kappa_align — neither books work). The misfit energy that drives the flow
is *already* leaving through exactly paired moves (roughness D→F with its
space share returned; churn during convergence — the choir's shed
disagreement, I.4). Plasticity redirects future conversion; it never
touches Es/Em/Ee/lem, so drift stays at the 1e−15 floor **by
construction**, and the E8-class prediction flips: an annealing pair pays
the comma in geometry, not energy (measured in the model, §3.1). A
"plastic toll" variant (charge each |Δd| a fraction of link traffic to
field) was considered and rejected: adds a constant and a channel with no
measured need; revisit only if PLAST-1 shows unphysical free-lunch
behavior (it should not: the flow is dissipative in V, bounded by V ≤ 1).

**Vacuum exactness [D].** Φ ∝ sqrt(m_i·max(m_j, floor)) with m = Em: in
vacuum Em = 0 on both ends ⇒ Φ = 0 *exactly* (the mob_floor arm is
multiplied by the zero live side). Field amplitude never enters Φ (dense
sector only), and dense transfer at instrument cells is already skipped
(`c==1 && cflag` guard) — so ld is frozen, bit for bit, in: vacuum, every
optics experiment (e2, d1, t1, q2, t4, p1), and everything Em-free.
Implementation must use an exact-zero skip (`if (Phi == 0) continue;`) so
no −0.0 write breaks byte identity. The battery-safety argument is then
structural, and the ratchet run proves it (PLAST-0).

**Determinism [D].** One new pass (pass D, between 5 and F, or fused into
pass 2): per-link Jacobi update from the pre-step snapshot (th2s, w2e,
Em, ld) — each `ld[l]` written only by its own link, `schedule(static)`
irrelevant to the result, byte-identical at any thread count. The
existing snapshot discipline (th1s/th2s) already provides the read set.

**The constant [P].** One new law-table entry, `kappa_plast`. Seam
justification: the force normalization is inherited from the S2 coupling
(the unitarity point fixed kappa_reac = 1; kappa_plast is the yield rate
of geometry under that stress — a mobility, dimensionally new, so a seam
is honest). Two candidate derivations to test against measurement rather
than assume: κ_plast ~ s_k (geometry flows at the space-transport
conductance — the (i)≡(ii) unification hook; predicts anneal time ≈ the
measured ~40 t.u. Es-equalization scale, mid-window), or κ_plast from
kappa_align's timescale (geometry class). The stability window is ≥ 50×
wide (§3.5), so the constant is soft — the benign seam class.

### 2.2 Option (ii) — METRIC-FROM-SPACE (d from live Es)

**Law [P].** No new DOF: d_l(t) = d_l(0) · [(Es_i+Es_j)/(Es_i(0)+Es_j(0))]^{1/3}
— the retardation is the space between; the same cbrt law that already
drives cr. Zero new constants. Space transport (s_k, already law) becomes
metric flow; the existing kernel's area/length asymmetry (§0.b) is
repaired; "matter is converted space" acquires its metric face: the g1
footprint (core Es ≈ 0.5× vacuum) becomes a **21% metric contraction
inside masses** (0.5^{1/3} = 0.794) — curvature as shortened rulers, not
just shrunk contacts.

**Does the pressure law already anneal? [D + M-model] Mostly no, as a
tuner.** Three findings:

1. **Rank [D]:** Es is a cell quantity — all links of a cell rescale
   together. V knobs cannot zero E residuals (trilemma above). (ii) can
   set the *scale* of a region, never the relative spread that kills
   shells.
2. **Sign [D]:** the two couplings are misfit-blind. Loading (s_pull,
   s_disp) lowers local Es ⇒ shortens ALL a mass's links, long and short
   alike; roughness's space return (back_s) raises Es ⇒ lengthens all,
   restoring only for the too-short half. No systematic restoring
   direction per link.
3. **Magnitude [M-model]:** back_s per rough unit = s_pull/(1+s_pull) =
   0.13; O(0.1) Es shifts give Δd/d = ΔEs/3Es ≈ 3–5% — below the required
   13–30% corrections.

**Is there an already-present annealing channel we never looked for?**
The channel exists (Es moves, cr moves, lA moves — strangulation is
measured) but it acts on areas only; d frozen means today's kernel has NO
metric response at all, so no signature can be in existing logs — nothing
to mine. **The discriminating instrument exists already**: the C1' piston
design proved rung LOCI ambient-invariant under the current pitch law
(bound-energy-only detune). Under (ii), ambient pressure changes Es
changes d changes every rung locus — **a measured locus shift under the
piston falsifies the current law and is (ii)'s fingerprint; a static
locus falsifies (ii)-at-strength.** Free tripwire, no new apparatus.

**Verdict:** keep (ii) as the principled long-run geometry unification
(it is what "space is a mode" ultimately wants) and as the piston
tripwire prediction; reject it as the annealer — wrong rank, blind sign,
weak magnitude.

### 2.3 Option (iii) — NODE MOTION [rejected]

Positions dynamical with re-linking. Costs, honestly:

* **Ontology:** promotes cx/cy/cz from scaffold to physics — the exact
  backslide (a permanent reference geometry) v89 exists to prevent. The
  no-background reading of "cells morph" is d-plasticity; node motion is
  its background-contaminated shadow.
* **Rank:** 3V < E on the bulk foam; and reachable d-vectors constrained
  to the Euclidean-embeddable cone — reintroducing the very frustration
  plasticity should dissolve (12 equal cube edges are impossible in the
  embedding, trivially possible relationally).
* **Machinery:** candidate table, CSR incidence, edge coloring, and the
  4NL flight slots all rebuild on re-link; in-flight lem on a deleted
  link needs a conserved relocation rule (ugly); determinism survivable
  (rebuild at fixed cadence, serial) but the cost is a new kernel, not a
  new law.
* **Redundancy:** everything node motion buys (lengths adapting), (i)
  buys with E clean knobs; the one thing (i) cannot do (change lux
  directions, acquire genuinely new neighbors beyond the candidate
  ceiling) has no identified physics need — collimation follows planes,
  which already adapt.

Minimal scoped fallback if ever needed: seeded cells only, moved every
N steps by the (i) force projected through lux, links recomputed among
the seeded set only — noted for completeness, not proposed.

### 2.4 Option (iv) — SHELL HARDENING ("morph until hard")

**Law [P].** (i) plus a per-link consolidation state: locked time T_l
accumulates while gg > g_thr (decays otherwise), and the plastic rate
falls with it:

```
κ_l(t) = κ_plast / (1 + T_l/τ_h)        (or exp(−T_l/τ_h))
T_l   += dt   if  G(ψ_f)G(ψ_b) > g_thr   (g_thr ≈ 0.9)
T_l    = max(0, T_l − dt)  otherwise      — a kick re-softens (work-hardening
                                            with annealing on unlock)
```

One more constant (τ_h ~ 50–100) and one per-link double. Deterministic
(same Jacobi pass). Conservation untouched.

**What hardening is FOR [D].** Plain (i) has a flat direction: the lock
manifold is the whole tuning-curve valley ω(x)·d = πmC — a pair can slide
along it (d tracking x) with gates open, so load drift is benign but the
mass within a species is a continuum segment, and nothing pins the
annealed geometry against slow wander. Hardening freezes d once the lock
has held, which (a) re-pins x*(d) — restoring force against load drift
returns, two-sided now because the geometry was annealed first; (b)
restores the discrete mass spectrum (species = closure integers × the
hardened lengths, not the foam's accident); (c) realizes the user's
"strong glass": stiff (frozen d + entrainment) and brittle-but-healing (a
kick that unlocks long enough re-softens the very links that were hit —
plastic exactly where wounded, hard elsewhere). It is the C1-plateau
candidate: anneal removes the roughness channel (52% of the leak [M,M0]),
hardening removes the slow-wander channel, and what remains is the
skirt-road load leak — which the annealed skin, no longer radiating,
feeds ~an order of magnitude more slowly. Whether that crosses the
sharpened C1 bar (|dM/dt| ≤ 1e−4 late, or t_half ≥ 1e4) is exactly
PLAST-1/2's question — the first configuration with a mechanism for a
true equilibrium rather than a slow death.

In the reduced models hardening changed nothing during clean anneals
(edges lock in ~20 t.u. and stop moving anyway — κ falls after the work
is done) [M-model]: it is a late-time/identity law, invisible in the
anneal, load-bearing under drift and perturbation. That is the correct
shape for a consolidation rule.

---

## 3. Stability and convergence [M-model]

Reduced models with laws_V2g constants (w2 = 2.9, q = 1.2, p_gate = 8,
Γ_m = 0.10, rough_k = 0.35, cap = 2.5, mob_floor, lens geo(d) with
r0 = 0.85, kernel dt = 0.02). Reference voice x = 0.32 ⇒ ω = 2.0954,
d* = 1.4993 (the foam's own d̄ = 1.505 — P2's observation), pinch at
2r = 1.70 vacuum, 1.64 after s_pull seeding.

### 3.1 The pair: two-sided attractor, comma paid in geometry (`plast_pair.py`)

Coupled (θ₁, θ₂, x₁, x₂, d), slightly unequal loads (the generic case),
Euler at kernel dt.

**Frozen (κ_p = 0) reproduces the measured one-sided ratchet** [cross-check
against E6/E8]: δ₀ < 0 (overloaded) pairs shed energy and walk to the rung
(ret 0.36–0.83); δ₀ > 0 pairs freeze off-rung holding energy (ret 0.78–1.0,
gg ≈ 0.00–0.04); on-rung pairs bleed by det-churn and slide off positive.

**Plastic (κ_p = 0.5): the attractor is two-sided and total.** Every
capture-range δ₀ ∈ [−1.13, +0.75] row locks to gg = 1.000, δ_f = 0.000,
in t_lock ≈ 1–2 t.u. (model units; the kernel-calibrated statement is
"within the entrainment decade"), retention 0.88–0.98 on the frozen side
where the frozen pair holds gg ≈ 0. The δ₀ = +1.13 row (d₀ = 1.77 >
pinch) stays frozen — **the capture range is the live range**; plasticity
cannot reach through a pinched channel (no urge), which is the correct
no-suction physics and a seeding rule: never intend an edge beyond
d ≈ 1.64.

**The comma, re-paid [M-model].** Overloaded protocol (E8 analog,
δ₀ = −0.85): frozen sheds 0.32 of the pair energy to reach the rung
(energy pays); plastic at slow κ_p sheds 0.14–0.22 with the rest paid by
d stretching 6–9% (geometry pays); at fast κ_p the geometry path races to
the pinch and the pair seals dissonant (retention 1.0, gates 0, silent) —
a failure of *rate*, not of law, fixed by the timescale window below.
Pre-registered kernel prediction: an E8 rerun with plasticity shows
**shed decreasing with κ_plast while final gg stays ~1** — the defect
ledger moves from the energy books to the geometry.

**Valley sliding [M-model].** A plastic pair that bleeds does not detune:
d tracks the valley (observed: x 0.34→0.107 with δ pinned at 0.000, gates
1.0 throughout). Detuning death (E9's fifth-killer) is abolished at the
configuration level; what remains is pure energy leak, which stops when
det → 0 (equalized pair, R = 0). This is the mechanism behind the §2.4
spectrum remark and the reason PLAST-1 predicts the leak's roughness
component vanishing rather than merely shrinking.

### 3.2 Networks: exact consonance on bipartite graphs (`plast_net.py`, `plast_net2.py`)

Phases + plastic lengths, uniform pitch, ±12–15% foam-measured jitter,
ambient phase noise 0.02.

* **ring4** (all live): gg = 1.0000 every link, max|ψ| = 1e−4, **cycle
  closure Σωd/2πC = 2.0000 exactly** — the flow lands on the integer.
* **ring6**: same when all links live; a seed with one link thrown beyond
  the pinch anneals as a locked open CHAIN (5/5 live gates 1.0,
  non-integer loop sum — correct: a chain has no cycle constraint). The
  B1-fleet-v2 accident (open chains) is thus a *generic outcome class*,
  now understood: the pinch, not the flow, decides ring vs chain.
* **cube, seeded-lock phases (kernel practice), 3 seeds** — the headline:

  | | gg min | gg mean | d spread | t(all > 0.99) | delivery-churn/link |
  |---|---|---|---|---|---|
  | frozen | 0.00–0.45 | 0.49–0.70 | 6.4–7.4% | — | ~0.08 (half-open) |
  | plastic | **0.987** | **0.996** | **0.00%** | **13–24 t.u.** | ~0.0036 |

  **The flow manufactures the one shared strut length the foam denies**
  (spread 7% → 0.00%), 80× inside the death time, and the residual is the
  noise floor, not the geometry (gg saturates at 0.99 against phase noise
  — CONSONANCE's tongue-vs-noise limit, now for the annealed state).
  Delivery-weighted churn drops ~22× per link.

### 3.3 Frustration: the flow amputates, it does not temper [M-model]

With d free, full lock requires per link ωd/C = πm and around every cycle
Σ±m even [D]. Bipartite + all-live ⇒ solvable exactly (measured, §3.2).
Odd cycles (all-m = 1 unreachable): the pre-registered equal-split
prediction (each link carries |δ| = π/N) was **WRONG, informatively**:
the measured flow concentrates the whole defect on ONE link and locks the
rest exactly —

* ring3: gg = [0.000, 0.978, 1.000], ωd/π = [0.998, 1.000, 1.000]
* ring5: gg = [0.000, 0.984…1.000], seam link stretched to 1.104π

The high-power gate potential is flat at large |ψ|, so one dark link
costs less total V than N tempered ones: **frustrated topology
self-truncates to its largest consonant subgraph** — a seam forms and
goes silent, exactly the defect-concentration the kernel showed
spontaneously (B1 v2 seam, H1 co-tree gates 0). Temperament (even
splitting) is what *phases alone* do (E6's measured δ/2 sharing);
amputation is what *geometry* does. Odd shells don't die under
plasticity — they shed an edge and become open consonant structures.
(Negative-control prediction for P10's Möbius/odd candidates: one dark
seam, remainder locked — sharp and checkable.)

### 3.4 Parasites: darkness by parity, takeover from cold [D + M-model]

Under a LOCKED bipartition, same-parity vertices are in phase, so a
parasite's gate argument is ψ = −ωd/C mod 2π regardless of jitter.
Gate > 0.01 needs |ψ| < 1.45 ⇒ d < 0.69 (below dmin = 1.0 — unreachable)
or d > 2.31 (beyond the ceiling) [D]. **Every candidate-range parasite of
a locked shell is dark**: gates 1e−5 … 1e−21 across d = 1.0–1.95
(computed). The model confirms: seeded-lock cube + 4 live diagonals ⇒
edges anneal to gg 0.99 while every parasite stays dark (gg_live =
0.000), unmoved, delivering nothing. **The annealed skin silences its
parasites by parity — self-sealing without expulsion.** (H1's parasites
radiated because the H1 skin never locked — gates 0.6 leave everything
half-open. Sealing the intended edges IS the parasite fix.)

The honest converse [M-model]: from RANDOM phases (natural formation, no
seeded bipartition), short live parasites can win the entrainment contest
— one seed annealed the diagonals to gg 0.995 and collapsed the intended
edges (mean 0.58): the flow found a *different* consonant structure than
designed. So plasticity is not a license for sloppy topology (P2
parasite-free geometry stands for engineered seeds), while for
condensate-annealing routes (C-tree) the same behavior is the feature:
the fabric form-finds for itself. PLAST-3 measures which face shows in
the kernel.

### 3.5 The κ_plast window and failure modes [M-model]

Pair model, frozen-side start (δ₀ = +0.75):

| κ_p | outcome |
|---|---|
| 0.05 | anneals too slowly — pair bleeds out first (retention 0.01) |
| 0.1–5 | clean: lock in ≤ 10 t.u., retention 0.61–0.98, plateau 0.2–1.0 |
| ≥ 10 | overshoot instability: d chases phase transients faster than entrainment follows; overshoots the rung (d_f 1.13–1.20 d*), gates 0 |

A ≥ 50× stable window bounded below by the bleed race (τ_plast < τ_leak)
and above by adiabaticity (τ_plast > τ_lock ≈ 10 t.u.). Wide-window ⇒
soft constant [supports the §2.1 seam argument]. Additional structural
failure modes, all understood: the **pinch stall** (overloaded pairs whose
rung lies beyond strangulation seal themselves dissonant-but-silent —
§3.1; mitigated by slow κ_p letting the energy path share the comma);
the **dead-gate desert** (links caught near ψ = ±π feel ~zero force AND
~zero phase pull — G and G′ both ~h^{p−1}; noise-free gradient flow can
trap there; the kernel's ambient churn is the physical unfreezer — e7-P5
measured a 60-t.u. saddle escape; the reduced model at noise 0.02
under-escapes, so cold random-phase ring acquisition is *not established*
— seeded-lock practice sidesteps it entirely).

### 3.6 What annealing does to species [D, consequence]

Plain (i) trades the foam-pinned point spectrum for **integer topology ×
continuous valley segments** (any (d, x) on ω(x)d = πmC locks). Discrete
labels survive (N, interval content, closure integers m — ring6's m = 3,
comp12's m = 5); the mass within a label becomes a bounded segment (cap
and skirt cut the valley). Hardening re-pins the flat direction at the
annealed point, restoring sharp masses whose values now come from the
*structure's own history* rather than the dart-throw. If the theory is
right, hardened-annealed masses across foam seeds should cluster far
tighter than seeded masses ever did — an ensemble prediction (PLAST-1,
C3-style) distinguishing plasticity from the frozen program.

---

## 4. Experiment designs (ratchet-gated; run after implementation)

All under laws_V2g + the two new keys (`kappa_plast`, `tau_harden`;
`kappa_plast=0` MUST be the exact legacy kernel). Every experiment
includes the standard conservation row (drift ≤ 2e−15) and the full
battery runs FIRST (ratchet rule 1). New cheap instruments wanted:
runtime `# NETGATE` rows (per-edge gates each diag window, not just at
seed — the D3 debt) and a `# PLAST` row (per-window: mean|Δd|, max|Δd|,
n_links moved, Σ|ḋ|, hardened fraction).

### PLAST-0 — vacuum and optics inertness (the gate)

* **Arms:** (a) kappa_plast = 0: FULL battery + one vacuum run
  byte-compared against the pre-change kernel binary's output. Bar:
  **byte-identical logs, 20/20 green.** (b) kappa_plast = κ* (the chosen
  production value), init=vacuum T = 100: bar: **every ld bit-identical
  to seed** (the exact-zero skip proves Φ ≡ 0), drift ≤ 2e−15. (c)
  kappa_plast = κ*, full battery: bars: all optics/vacuum experiments
  (e1, e2, d1, t1, q2, t4, p1, qt_lo/hi, e5) **log-identical to arm (a)**
  (their free cells never carry Em); dense experiments (e3a/b, e4,
  e6–e9, g1, g3, g4) within their existing bars — green, values may
  move. **Kill:** any optics log differs, or any battery bar fails, or
  vacuum ld moves by one ulp.

### PLAST-1 — the A/B: the cube anneals (the claim)

* **Seed:** P2 parasite-free cube via init=net (a ≈ 1.45–1.5, all 12
  intended edges live: d < 1.60 at seeded load — reject picks otherwise;
  the NETGATE parasite scan must report its in-ceiling diagonals, which
  the parity argument predicts dark). Load at the π-rung tuning curve,
  seedlock phases. T = 5000, closed box (the A1c protocol), lump_diag on.
* **Arms:** kappa_plast = 0 (frozen control — the measured baseline
  class: gates ~0.6, death by skirt road) vs kappa_plast = κ*.
* **Pre-registered predictions:** (P-a) runtime NETGATE gates anneal from
  ~0.6 to mean ≥ 0.9 within t ≤ 200 (model: 0.996 at ~20 t.u.; kernel
  slower is fine, direction is the claim); (P-b) `# PLAST` shows Δd
  converging to a plateau with mean|Δd| ≈ 0.1–0.3 (the §1 requirement)
  and near-zero flow after anneal; (P-c) roughness ledger rate during
  [200, 1000] at most HALF the control's (the 52% channel closing);
  (P-d) leak: late |dM/dt| below the control's at matched epoch;
  **success** = P-a AND P-c AND P-d; **the prize** (not required for
  success): a C1-sharpened plateau — late |dM/dt| ≤ 1e−4 or extrapolated
  t_half ≥ 1e4 — the first equilibrium candidate. **Kill:** gates fail
  to improve (P-a false) — the force is mis-derived; or leak WORSENS
  while gates improve — annealed skins radiate through a channel the
  theory missed (that result would be as informative as success).
* **Ensemble rider (C3 discipline):** ≥ 5 foam seeds; the §3.6
  prediction: annealed-mass dispersion across seeds < seeded-mass
  dispersion.

### PLAST-2 — anneal-then-kick (hardening and identity)

* **Seed:** the PLAST-1 winner, tau_harden ∈ {0 (off), 50}. Run to
  t = 1000 (annealed + hardened), then the H2 protocol: weak aux field
  packet strike (sub-fracture), continue to t = 3000.
* **Pre-registered:** (P-e) hardened arm's post-kick gate dip recovers to
  ≥ 0.9·pre-kick within ~100 t.u. (re-softening exactly where wounded —
  the work-hardening signature: `# PLAST` shows κ_l reawakening only on
  the struck links); (P-f) hardened arm out-lives the never-hardened arm
  under repeated kicks (the "strong glass": frozen-d restoring force in
  x); (P-g) breathing-mode spectroscopy (STABILITY agent's H2 lines) is
  SHARPER on the annealed shell than the frozen control (uniform struts
  ⇒ degenerate modes ⇒ narrow lines) — the cleanest falsifiable
  spectral signature of annealing. **Success** = P-e + C4-class identity
  (returns to mass plateau). **Kill:** hardening prevents healing
  (brittle shatter with no re-softening) ⇒ the T_lock decay rule is
  wrong, iterate the consolidation law.

### PLAST-3 — the naive a=1.25 cube WITH parasites (self-sealing)

* **Seed:** the original H1 geometry (a = 1.25 class, its 6–12
  in-ceiling diagonals as measured), seedlock phases, both arms as
  PLAST-1.
* **Pre-registered (the parity mechanism, §3.4):** (P-h) with plasticity,
  intended-edge gates anneal high AND every parasite NETGATE stays
  < 0.01 the whole run (dark by parity — they never wake because the
  skin locks); (P-i) total roughness rate below the frozen control
  DESPITE the parasites (control: half-open everything radiates);
  (P-j) honest takeover watch: if instead edges collapse while diagonal
  gates rise (the §3.4 cold-start mode), plasticity amplified the wrong
  structure — records that engineered seeds MUST be parity-clean and
  P2 remains load-bearing. **Success** = P-h + P-i. Either outcome
  settles the "does plasticity feed or seal parasites" question with a
  mechanism attached.

Sequencing: PLAST-0 (gate) → PLAST-1 (claim) → PLAST-3 (cheap, same
tooling) → PLAST-2 (needs hardening implemented; can trail).

---

## 5. Verdict and recommendation

**Diagnosis confirmed and sharpened.** The frozen foam denies consonant
shells the single shared strut length they need (16% spread where ≤ 2% is
required; gates 0.6/min 0 reproduced from pure geometry); no pick
refinement reaches it; and every adaptive channel the kernel owns steers
toward consonance except the one frozen number that vetoes them. The
user's mechanism is the missing DOF, and its no-background form is exact:
**cells do not move — the retardation between them yields.**

**Recommended mechanism: (i) retardation plasticity, hardened by (iv).**
It is the unique full-rank, background-free option (§2 trilemma); its
force is computable from pass-2 quantities with ∂ψ/∂d = −ω/C and is the
S2-derived reactive term acting on its geometric conjugate (the term
rate-transport provably cannot host finds the variable nothing else
owns); it costs one soft constant (κ_plast; window ≥ 50×) plus τ_h for
hardening; conservation is exact by construction (geometry is not a
ledger — precedented by cr and the normals); vacuum and the whole optics
battery are frozen bit-for-bit by the Φ = 0 structure. (ii) is retained
as the long-run geometric unification and its piston tripwire; (iii) is
rejected as background reintroduction with insufficient rank.

**What it predicts that the frozen program cannot:** two-sided rung
capture (the frozen ratchet is one-sided [M]); the comma paid in geometry
(E8's shed falling with κ_plast at constant final gg); strut spreads
annealing 7% → ~0 with gates at the noise floor in ~20 model-t.u.;
frustrated topology amputating a seam instead of tempering; parasites
dark by parity under an annealed skin; degenerate (sharpened) breathing
lines; annealed-mass clustering across foam seeds; and — the prize — the
first true equilibrium candidate (roughness channel closed by annealing,
wander channel closed by hardening).

**The consistency check (form-finding vs plasticity).** The form-finding
solver can only SELECT geometry; plasticity CORRECTS it. If the theory is
right they converge to the same consonant states: same bipartite
topology, same closure integers, same tuning-curve relation — the
annealed (d, x, θ) of a sloppy seed should match the solver's solution on
the same foam link-by-link up to valley degeneracy. Divergence between
them is a falsifier for one side or the other; agreement is the theory's
first non-trivial self-consistency loop. (Concretely: run formfind on the
foam; seed its cube DELIBERATELY mis-tuned; anneal; compare.)

**Smallest implementable step:** pass D in cellfab.c — ~25 lines: per-link
Jacobi update `ld[l] += -kappa_plast * base*Sm * dVdd(psi_f,psi_b) * dt`
with the exact-zero skip, floor d ≥ 0.5, one `# PLAST` diagnostic row,
`kappa_plast=0` default → PLAST-0 gate (byte-identity + 20/20) →
PLAST-1. Hardening (τ_h, one array) lands as a second ratchet-gated
change only after PLAST-1 reads out.

**LAW-change class notice:** everything here is design; any
implementation is a kernel/law change and runs the FULL battery before
commit (ratchet rule 1). The vacuum-inertness argument in §2.1 is the
reason to expect green; the battery is the proof.

---

*Models and numbers reproducible via:*
`python3 prestress/theory/plast_foam.py` (diagnosis),
`plast_pair.py` (pair: two-sidedness, comma-in-geometry, κ window),
`plast_net.py` + `plast_net2.py` (networks: exact bipartite convergence,
seams, parasites, noise floor). No simulator runs were made.
