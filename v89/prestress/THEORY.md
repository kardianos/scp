# CONSONANT NETWORK THEORY — discrete prestress / form-finding for v89

*Stream A deliverable (campaign LEDGER.md). 2026-07-28. Derived from
`v89/PRINCIPLE.md`, `v89/CONSONANCE.md`, `v89/MASS.md`,
`v89/design/B1_winding_compensation.md`, `v89/battery/laws_V2g.cfg`, and the
kernel `v89/cellfab.c` — no pre-v89 material. Every kernel fact cited below
was re-verified against the source at the current line numbers; where the
tasking brief or the ledger was wrong, the correction is stated in place and
collected in §8. Runnable scripts with quoted outputs: `theory/t1…t5*.py`
(+ `.out` files).*

Notation: `gate(ψ) = ((1+cos ψ)/2)^p = cos^{2p}(ψ/2)`, `p = p_gate = 8`;
pitch `ω(x) = w2/(1+q·x)`, `x = (Em+flload)/cap`; laws_V2g constants
throughout (`w2=2.9, q=1.2, cap=2.5, Γ_m=0.10, p=8, C=1, dmin=1.0, r0=0.85,
rjit=0.06, s_pull=0.15, k_dep·k_dep_m=2.4, mob_floor=0.004,
lock_floor=0.005, rough_k=0.35, Γ_r=0.5, A0=1.15`).

---

## 0. Kernel verification (what the theory stands on)

Verified in `cellfab.c` (current numbering; the `init=9` net seeder added by
stream B shifted post-1188 lines by +148, physics passes untouched — diff
audited):

| fact | where | status |
|---|---|---|
| `gate_of(ψ) = ((1+cosψ)/2)^p`, even, p=8 | 489 | ✓ |
| comb: coprime p:q, p·q ≤ 6, simplest-first tie-break; width Γ_m/(pq), amplitude 1/(pq) | 443, 1786–1795 | ✓ |
| foam: dart-throw `dmin=1.0`; radii `r0·(1±0.06·U)` (uniform, not gaussian); candidate links iff `d < 1.15·(cr0_i+cr0_j)` (NOMINAL radii); LIVE transport needs lens overlap of LOADED radii `cr = cr0·(Es/e_s0)^{1/3}` — `A≤0` links carry nothing, flight on pinched links holds | 593–719, 1698, 1724 | ✓ (brief said "radii r0·(1±0.06·U)" — confirmed) |
| dense gates: `ψ_ij = q·θ_i − q·ω_i d/C − p·θ_j`, `ψ_ji = p·θ_j − p·ω_j d/C − q·θ_i` (1:1 ⇒ `θ_i − ω_i d/C − θ_j`) | 1796–1797 | ✓ |
| transport chain: `want_ij = k_dep·k_dep_m·dt·geo·gpl·res·gate_ij·head_j·mob_i^eff`; `geo=(A/A_ref)(d_ref/d)`, `gpl=axi·axj·(n̂2_i·n̂2_j)²`, `res` = comb Lorentzian, `head_j = 1−(Em_j+Ee_j)/cap`, `mob^eff = sqrt(mob_i·max(mob_j, mob_floor·cap))` | 1719–1877 | ✓ |
| **transport is continuous** — atoms fire ONLY at mode boundaries: pass-5 roughness, pass-6 condensation/evaporation (`atoms_fire` call sites) | 1917–1920, 2010–2016, 2226–2266 | ✓ |
| flight: deposits enter `lem`; phase advances `dt·C/d`; delivery of the whole slot on cycle completion; receiver-full ⇒ wait a cycle; roughness `R = 2|det|Γ_r/(det²+Γ_r²)·rough_k` of the take, receiver-side, D→F with space return `s_pull/(1+s_pull)` | 1969–2042 | ✓ |
| clock: `th2 += w2e·dt` (pass 6) — exactly ω(x) per step, x from step start | 2216 | ✓ |
| **everything else that moves th2**: receiver-only entrainment at deposit (`mix = f/(f+mob+lock_floor·cap)`) and at delivery (`mix = take/(take+mob_prev+lockf)`), both `θ += κ_lock·mix·wrap(ms(θ_src−ω_src d/C) − mr·θ)/mr`; seeding. Alignment (pass 7) moves plane normals only. th1 is a diagnostic cache of the field amplitude phase | 1902–1966, 2023–2036, 2199–2201 | ✓ (brief's "exactly ω(x) per step" is TRUE of the free advance but INCOMPLETE without entrainment) |
| pitch: `x = (Em + flload)/cap`, `flload_i = ½·Σ_incident(lem_fwd+lem_back)` | 1603–1617, 1707 | ✓ |
| seeders: pair/ring/shell/net all use lock recursion `θ_next = θ − ω·d/C`; ring `init=7` is SEQUENTIAL around the loop (seam = last link); `ring_m ≥ 0` sets uniform `ω = 2π·m·C/L_ring` from the ACTUAL loop; shell `init=8` = BFS tree over cube edges, `ω = π·C/ā`; net `init=9` (stream B) reads `V x y z x_load θ` / `E a b`, scores intended edges AND parasitic candidate links (`# NETGATE`, `# NETGATE P`) | 841–1327 | ✓ |

Independent cross-check: stream D's `theory/ring_map.py` (exact reduced
per-step map) lists the identical pass structure.

---

## 1. THE DRIFT THEOREM

**T1 (drift).** For a static link (d, ω_i, ω_j constant) with free-running
clocks, every directed gate argument drifts at the detune:

    ψ̇_ij = q·ω_i − p·ω_j =: det ,   ψ̇_ji = −det          (1:1: ψ̇ = ω_i − ω_j)

*Proof.* `th2` advances at exactly `ω·dt` per step (pass 6, line 2216); the
retardation term `q·ω_i d/C` is constant; nothing else moves θ in a locked
steady state — the entrainment kicks are proportional to `wrap(lock error)`
and vanish identically at the lock (§0). ∎

**Corollary (static lock ⇒ rational pitch).** A statically open gate needs
det = 0. On a 1:1 link that is `ω_i = ω_j`; comb links allow `ω_i/ω_j = p/q`
(p·q ≤ 6). A connected component whose links are all 1:1-locked therefore has
**one ω** — and hence, through the pitch law, **one load x(ω)** on every
voice. (Non-uniform ω is either a comb lock — §4 — or drift = churn.)
Entrainment does not evade this: it is a finite-bandwidth phase-locked loop
(capture ≤ min(0.30·R₀, 2Γ_m), t5 §4); within capture it holds a lock *with
residual error* (the tempered comma), outside it the link slips at det with
mean gate `⟨gate⟩ = C(16,8)/2¹⁶ = 0.19638` (t1-A).

**T2 (the strut ladder).** Both directions of a link statically open ⇔

    ψ_ij + ψ_ji = −(q·ω_i + p·ω_j)·d/C ≡ 0  (mod 2π)
    ⇒ (q·ω_i + p·ω_j)·d/C = 2π·m            (1:1, uniform ω:  d = π·m/ω)

and the forward lock pins `q·θ_i − p·θ_j = q·ω_i d/C (mod 2π)`. On-rung the
per-link phase drop is **always the half-period**: `φ = q·ω_i·d/C = π·m` —
for every interval, not just 1:1 (lock ⇒ q·ω_i = p·ω_j ⇒ φ = 2πm·q·ω_i/(2q·ω_i)
= πm). Odd m = antiphase-grade, even m = inphase-grade (this drives §4).

**T3 (reachable m, 1:1).** With the pitch window ω ∈ [1.3942, 2.7000]
(x ∈ [x_skirt, 0.9]; x_skirt = 2Γ_m/(q(w2−2Γ_m)) = **0.061728**, derived: the
voice enters the room's ±2Γ 1:1 acceptance and dissolves — the measured e7/A1
law) and the live-contact ceiling `d < r_i(x)+r_j(x)`, `r(x) = r0·(1 −
s_pull·x·cap/((1+s_pull)e_s0))^{1/3}` (seeder bookkeeping, verified):

* m=1: d = π/ω ∈ [1.1636, 2.2533] — reachable, but the contact ceiling cuts
  it: solving `π/ω(x) < 2r(x)` gives the **strut window**

      x ∈ (0.0617, 0.413),  d = π/ω ∈ (1.164, 1.620)     [nominal radii]
      x ∈ (0.0617, 0.348),  d ∈ (1.164, 1.535)           [worst −6% jitter]

* m≥2: d ≥ 2.327 > worst-case candidate cut 2.0723 = 1.15·2·r0·1.06 —
  **no channel exists**. (t1-B. Brief's "only m=1, d ∈ ~(1.14, 2.3)"
  confirmed and sharpened: the loaded-contact ceiling, not the candidate cut,
  is binding.)

So: **all 1:1 struts of a component share one length d = π/ω and are
antiphase** — with the caveat of T2 that comb links can carry even-m rungs at
in-window lengths (§4).

**T4 (jitter tolerance).** Length error δd on a strut (target d* = π/ω):

* seed/tree lock (kernel seeders): forward gate 1 exactly, back gate carries
  the whole round trip: `g_b = gate(2ω·δd) = cos¹⁶(ω·δd)`;
* entrained equilibrium (E6-P2's even tempering, = the LS solution §3):
  both gates `cos¹⁶(ω·δd/2)` — **tempering doubles the tolerance** (product
  e^{−4(ωδd)²} vs e^{−8(ωδd)²}).

Thresholds (t1-C), ω-independent as a fraction of strut length:

| criterion | tree lock | tempered |
|---|---|---|
| gate ≥ 0.9 | δd ≤ 0.1146/ω = **3.65 % of d** | δd ≤ 0.2293/ω = **7.30 % of d** |
| gate ≥ 0.99 | δd ≤ 0.0354/ω = **1.13 %** | δd ≤ 0.0709/ω = **2.26 %** |

Foam nearest-pick placement error is ~10–15 % of d ⇒ raw seeded shells sit at
mean gate 0.6–0.7 with a dead co-tree edge — exactly the measured H1 cube
(0.60–0.64, min ≈ 0); the annealed/LS assignment at the same jitter reaches
0.73–0.84 mean, min 0.23–0.48 (t3-D table). The E6 tongue formula
`G²(δ) = cos³²(δ/4)` (δ = rung offset = 2ω·δd) reproduces the measured E6
points to ≤0.01 in-tongue (t1-D). *(Correction: B1's "gate ≥ 0.9 for |ψ| ≤
0.16" is actually gate(0.16) = 0.95; the 0.9 point is |ψ| = 0.229.)*

---

## 2. CABLES (forward-only links)

**T5 (tree lock).** ψ_f = 0 on every cable of a spanning tree is exactly
solvable by the lock recursion `θ_j = θ_i − ω d_ij/C` (what every kernel
seeder does). Per-link drop φ = ω d/C is free within the **band**

    φ ∈ B(ω) = [ω·dmin, ω·(r_i(x)+r_j(x))]    (x = x(ω) — one load per component)

with B = [1.394, 2.111] at ω=1.394 … [2.700, 4.559] at ω=2.700 (t1-B table).

**T6 (cycles).** Each independent cycle adds one closure condition on the
*signed* drop sum: traversing a cable along its lock direction contributes
+φ, against it −φ:

    Σ_c ±φ_l ≡ 0 (mod 2π)   ⇔   ω·Λ_c/C ∈ 2π·Z,  Λ_c = signed length.

* **Balanced cycles (Λ_c = 0) impose nothing** — free tunability. This is a
  designable resource: parallel same-sense paths of equal length close at
  any ω.
* Unbalanced cycles quantize ω = 2πm_c C/Λ_c. F cycles ⇒ F−1 commensurability
  conditions on geometry: theta-graph L=15 vs 9 tunes at ω = 2.0944 (m=5,3),
  15 vs 10 at 2.5133 (m=6,4), 15 vs 11 — never; a 0.5 % foam jitter (15.08 vs
  9) already breaks exact common ω (t2-E) ⇒ **multi-cycle nets must either be
  built from commensurable loops or absorb per-cycle commas ≤ N·0.229 rad
  (all gates ≥ 0.9, tempered; §3)**. The standing seed-acceptance
  "closure integer to < 0.05" (0.314 rad) is compatible for N ≥ 2.
* Torus nets N₁×N₂, uniform φ: φ = 2πm₁/N₁ = 2πm₂/N₂ ⇔ m₁N₂ ≡ 0 (mod N₁) —
  homology pair (m₁, m₂) as the species label (ledger P6).

**Back-gate leak and the diode point.** At forward lock the back gate is
exactly `gate(−2φ)` (verified: ψ_ji = −(ω_i+ω_j)d/C = −2φ at lock — the
brief's claim is confirmed):

    φ = π (m even rungs): back gate 1 (strut-grade, mutual)
    φ = 5π/6 (comp12):    0.1001      φ = 2π/3 (comp6, wound-3): 1.53e−5
    φ = π/2:              **exactly 0** — the DIODE POINT (cos(π)= −1 zero)

Diode window: φ = π/2 needs d = π/2ω ≥ dmin ⇒ **ω ≤ π/2 ⇒ x ≥ 0.705, d ∈
[1.0, 1.127]** — heavy voices only; the second zero φ = 3π/2 needs d ≥ 1.745 >
contact — unreachable (t1-B). The leak-minimizing distribution of φ around a
closed cycle subject to Σφ = 2πm: uniform φ = 2πm/N is the stationary point,
and a local minimum of Σ gate(2φ_l) wherever cos¹⁶ is convex (|φ−kπ| >
0.253); the global optimum saturates at back gate ≡ 0 iff φ = π/2 is
reachable and N = 4m — mixed {π/2, 3π/2} solutions require both diode lengths
and are geometrically excluded. **What actually flows through a partial
gate**: the gate multiplies the continuous want (a rate) — nothing is
probabilistic; the flow enters flight and delivers per cycle; roughness fires
only at DELIVERY and only from comb detune (`R(det)`), NOT from gate offset —
a uniform-ω cycle with imperfect gates churns and slows but radiates nothing
internally (B1's "winding cost is purely in the gates", verified). The
radiative losses of a uniform-ω structure are entirely at its **rim** (vacuum
detune Δω = w2·qx/(1+qx)), §5.

---

## 3. THE FORCE-DENSITY MATRIX (central deliverable)

**Setup.** Per component: unknowns θ ∈ T^n (mod one gauge) and ω (1). Write
every directed gate constraint as a row: cable → one row `θ_i − θ_j = φ_l`;
strut → **two** rows (`±(θ_i −θ_j) = φ_l`): a strut is a 2-cycle of the
directed constraint system whose defect is the pair-rung residual. With
integer comb multipliers the row is `q·θ_i − p·θ_j = q ω_i d/C` (matrix
M_pq ∈ Z^{rows×n}).

**T7 (W from gate″).** Near an operating point the delivered flow through a
row is `K_l·gate(ψ_l)` with K_l = k_dep k_dep_m·geo·gpl·res·head·mob (the
pass-2 chain sans gate). Since `gate(ψ) ≈ 1 − (p/4)ψ²` (gate″(0) = −p/2;
Gaussian top e^{−pψ²/4}, valid |ψ| < ψ_infl = 2 atan(1/√(2p−1)) = 0.505),
maximizing total flow = minimizing

    Σ_l w_l ψ_l²,   w_l = K_l·p/4      ["phase stiffness" = force density]

**T8 (the Laplacian).** Stationarity of the weighted LS on the torus is

    Q θ = B W φ,   Q = B W Bᵀ   (B = signed incidence over the rows)

— the v89 force-density matrix. Solved by wrapped Gauss–Newton in the seed's
winding sector (t3), which is also what kernel entrainment does dynamically
(E6-P5 acquisition ~10 t.u.). Verified structure (t3-A):

* **ker Q = the gauge**: for 1:1 the constant vector (global clock shift);
  for interval networks the **ω-vector** (θ → θ + ω·τ = time translation;
  `M_pq·ω = det = 0` on lock). Extra kernel beyond gauge = *mechanisms*:
  free retunings (disconnected parts, balanced-cycle ω-freedom).
* **Residual ψ\* lies in W⁻¹·ker(Mᵀ)** — the weighted cycle space
  (|MᵀWψ*| ≈ 1e−15). Cycle defects are invariant: `Σ_cycle ψ* = −δ_c`
  exactly. coker Q ≅ ker Q (symmetric); the defect data lives in the cycle
  space, of dimension rows − (n−1) = (graph cycles) + (# struts).

**T9 (the 1/w comma law).** For a single cycle with defect δ:

    ψ*_l = −δ·(1/w_l)/Σ_k(1/w_k),    min Σ w ψ² = δ²/Σ(1/w)

— commas distribute in **inverse proportion to stiffness** (series
compliance). Verified numerically to 4 digits; uniform w ⇒ even tempering =
the measured E6-P2 comma split; w→0 on one link ⇒ the seam dump measured in
the fleet-v2 fixed-ring_x seeds (seed: seam ψ = −δ gate 0.48; LS: ψ = −δ/12
each, gates 0.995 — t3-B). A strut with length error is the 2-cycle special
case: LS splits the round trip evenly — tempered gates beat the tree lock's
product (t3-C). Consequence for design: **soft links attract comma** — a
deliberately soft "sacrificial" link (small A or bad alignment) can shield
the working links of an imperfect cycle; conversely a hard net distributes
evenly.

**T10 (counting).** Per component: unknowns n (θ mod gauge + ω). Forward
locks: E rows, rank n−1; B = E−n+1 cycle conditions with integer freedoms.
Struts contribute (i) their θ-row pair, (ii) the geometric rung d = πm/ω —
with m=1 only (1:1), E_s equal-length conditions absorbed by form-finding
(3n−6 placement dof, generically for n ≥ 6), (iii) a π-residue in every
cycle through them. On a fixed foam the generic deficit = (#unbalanced
cycles − 1) quantizations + E_s length equalities + the payment conditions
(§4) — the discrete parts are **not** removable by continuous tuning.

**Dictionary** (ledger table, corrected and completed):

| FDM tensegrity | v89 consonant network |
|---|---|
| node coordinates x | clock phases θ (+ pitch ω) |
| force density q = t/ℓ | phase stiffness w_l = K_l·p/4 (pass-2 chain × gate curvature) |
| equilibrium Dx = f | Q θ = B W φ |
| external load f | retardation targets φ = ω d/C |
| **self-stress** (member forces, equilibrium, no load) | **divergence-free circulating flux σ ∈ ker B** — winding (§5; steady state + conservation make it automatic) |
| **mechanism** | ker Q beyond the ω-gauge: free retunings; load-sector drift |
| prestress stabilizes mechanisms | flux adds entrainment damping + gyroscopic **precession**, no static stiffness; protection = the gate desert (§5) |
| form-finding: choose q, solve x | choose (topology, m, intervals); solve (picks, ω, θ, x*) on the actual foam |

**The circle-placement reduction** (for the formfind/morpho agents): since a
component has one ω, the whole static design is — place n phases on the
circle such that every strut pair is **antipodal** and every cable pair's
separation lies in **B(ω) ∪ (2π−B(ω))** (orientation = sign choice); then
d_l = |drop|/ω realizes the geometry, x = x(ω) the loads. All of §§1–4 are
constraints on this placement; Q is its quadratic relaxation; `init=9`
verifies it edge by edge (`# NETGATE`), including **parasites** (candidate
links you did not intend — see §6 rule 9).

---

## 4. BIPARTITENESS, GENERALIZED — and one refutation of the brief

**T11 (1:1 case).** All-strut components: every drop is π (m=1 forced by
T3), so every cycle needs an even number of struts ⇒ **1:1 strut networks
are bipartite, equal-length, antiphase** — the H-series selection rule,
derived. Odd cycles are frustrated exactly as the ledger says.

**T12 (the payment lemma, mixed cycles).** A cycle with n_s struts and n_c
cables requires `Σ_cables ±φ ≡ π·n_s (mod 2π)` with each φ ∈ B(ω) — and ω
pinned to the strut window when any strut exists. Feasibility scan (t2-A):

    triangle:  n_s = 0 only  (all-cable wound m=1, ω ≤ 2.094)
    pentagon:  n_s ≤ 2
    heptagon:  n_s ≤ 4

⇒ *"demote 2 links to cables" is wrong in both number and remedy*: one or
two cables can never pay an odd triangle (their φ-sum cannot reach 0 mod 2π
in-window); the triangle must be all-cable. **Frustration index (E − maxcut)
is only a lower bound and a bad one**: K4/tri-prism/pentag-prism have index
2, but the true minimum demotions are — K4: no consonant realization at all
(the two circle-diagonals always land in the band's forbidden hole; anneal
confirms at all ω, t2-C); tri-prism: 6 (both triangles all-cable — the
co-rotating triangular tube, verified design t2-D); pentagonal prism: 10
(the m=2 wound tube). **Hexagonal prism (ledger P4 "expected
over-constrained by 2"): REFUTED — it is bipartite (all faces even); 0
demotions, all-strut feasible** (explicit phase assignment verified).
Icosahedron: infeasible even all-cable (anneal, t2-F). Dodecahedron:
frustration index 6; pentagons cap struts at 2 each — feasible mixed forms
exist on paper but no verified assignment yet.

**T13 (interval rungs CAN beat bipartiteness — corrected claim).** The
exhaustive interval scan (t4, all comb assignments, both windows, rung
integers, and the integer left-null parity `Σ y_k m_k ≡ 0 mod 2`) **found
in-window odd strut-grade cycles**. The loophole is T2's half-period drop
with **even m**: an interval link's rung sum `q ω_i + p ω_j = 2 q ω_i` is
twice the naive scale, so m=2 lands at in-window length while contributing
drop 2π ≡ 0 — a ferromagnetic bond in the antiphase Ising picture.
Smallest example, **the fifth-triangle** (design card, t4):

    v0 —(3:2, m=2)— v1 —(2:3, m=2)— v2 —(1:1, m=1)— v0
    all three links both-gate on-rung; equal lengths d = π/ω0
    ω0 ∈ [2.09, 2.5] ⇒ x0 = x2 ∈ [0.13, 0.32], x1 = x(2ω0/3) ∈ [0.62, 0.90]
    e.g. ω0 = 2.35: d = 1.337, contact margin +0.27

Cost (why it is not a mass candidate): the 3:2 legs run at comb amplitude
1/6 and acceptance Γ_m/6 = **0.0167 rad/t.u.** (6× weaker, 6× narrower than
1:1), and E9 measured the fifth dying of detuning drift in ~20 t.u.
Pre-registrable prediction: forms, holds while the legs stay in the Γ/6
tongue, dies by interval drift — lifetime ≫ naked fifth is NOT expected.
Worth one fleet probe as the first non-bipartite species; the **1:1-only
bipartite rule stands** wherever interval content is absent.

**T14 (beat-consonance / stroboscopic lock: REFUTED).** Non-uniform ω
without a rational lock drifts at det and cannot be rescued by the action
atoms (t4-b):

1. transport within the dense mode is **continuous** (atoms only at D↔F
   boundaries — pass-4 comment + `atoms_fire` call sites): there is no
   transfer quantum for the gate window to synchronize with;
2. the delivery-time entrainment *is* beat-sampled, but stroboscopic
   freezing needs det·d/C = 2πk ⇒ det ≥ 2π/d = 5.03 — **25× beyond the
   measured unlock boundary 2Γ_m = 0.2**, where comb res = 4e−4 and mean
   gate 0.196: nothing left to lock;
3. the atoms quantize amounts, not timing; the credit lapses at 2 atoms and
   stores no phase;
4. what a sustained det actually does is **radiate**: rough fraction
   `R(det)·rough_k` up to 0.34 of every delivered parcel (peak at det =
   Γ_r = 0.5) — the measured particle-killer (M0: 52 % of blob loss), not a
   lock. The only discrete beat-locks are the comb's own p:q rungs, already
   in the law table, already measured fragile (E9 hierarchy).

---

## 5. FLUX AS PRESTRESS

**What plays tension.** The circulating dense flux σ_l (energy/time through
link l, = the resolved want rate at lock). At steady state with no
conversions (uniform-ω interior: det = 0 ⇒ no roughness; Em < cap; Ee ≈ 0)
conservation forces per-voice in = out: **σ is divergence-free — σ ∈ ker B,
exactly the FDM self-stress** (automatic from the paired-ledger bookkeeping
at steady state, not instantaneously; a leaking structure's divergence = its
bleed). C5's "internal cycling at constant mass" = a self-stress with zero
divergence including the rim. (t5-2.)

**T15 (flight-load fixed point — a seeding correction).** Flight is load
(pass 0): a slot's sawtooth inventory averages σ·d/2C, so a circulating ring
carries `flload = (σ_f+σ_b)·d/2C` per voice, and the seeder's `Em = x*·cap`
wakes up sharp of the rung by flload/cap and sheds to it (F1's one-sided
ratchet, measured as the early fast M drop). Fixed points (t5-1):

| class | gpl(N) | σ_f | flload | flight % of rung mass | corrected seed x |
|---|---|---|---|---|---|
| ring6 (m=3) | 0.0625 | 0.023 | 0.028 | 8.9 % | 0.117 (was 0.128) |
| comp12 (m=5) | 0.5625 | 0.319 | 0.219 | **27.3 %** | 0.233 (was 0.321) |
| comp6 (m=2) | 0.0625 | 0.068 | 0.038 | 3.5 % | 0.421 |
| diode8 | 0.25 | 0.229 | 0.120 | 6.4 % | 0.702 |

Up to 27 % of a wound ring's rung mass is **census-invisible** (in `lem`) —
part of the measured "early leak" is inventory transfer, not loss (the M0
−0.039 residual's mechanism; the owed tET check will see it). The **ring
plane-alignment factor is a design variable**: seeded gpl = cos⁴(2π/N)
(planar ring geometry, derived t1-E): N=6 → 0.0625, N=12 → 0.5625, N=16 →
0.73; κ_align drifts it toward cos⁴(π/N). N ≤ 5 rings are conductance-dead
at seed (N=4: gpl = 0 exactly) unless axes are pre-set to the bisector.

**T16 (gyroscopic, derived).** Linearized clock perturbations on a ring obey
`u̇_k = R_f(u_{k−1}−u_k) + R_b(u_{k+1}−u_k)` with R = κ_lock·σ·[1/(Em+lockf)
+ 1/(Em+lockf+σd/C)] (both entrainment sites, receiver-only — kernel-exact).
Modes: `λ_j = −(R_f+R_b)(1−cos 2πj/N) ± i(R_f−R_b) sin 2πj/N`. **Decay is
set by the total exchange; net circulation contributes only precession.**
Circulation does NOT stiffen the lock — winding protection is not
gyroscopic. Numbers: comp12 R_f = 0.86/t.u., slowest decay 0.129/t.u.,
precession 0.376 rad/t.u.; ring6 an order weaker (R_f = 0.14). (t5-3.)

**T17 (the desert — B1 quantified).** One link against a detune:
`ψ̇ = δω − R₀·gate(ψ)·ψ` (restoring rides the flux, which rides the gate).
Max restoring at ψ* = 0.495 (gate 0.61): capture `|δω| ≤ min(0.30·R₀, 2Γ_m)`
— comp12: 0.20 (comb-clipped), ring6: 0.047. Beyond ψ ≈ 1.2 the restoring is
< 1 % (gate(1.5) = 6.7e−3, gate(2.618) = 4.05e−10 — B1's "~4e−10 at ψ≈−2.6"
verified to three digits, and gate(π) = 0 exactly): **the desert has no
restoring force AND no motor** — a slip needs a sustained one-sided δω for
T ≈ (2π−2.4)/δω (~19 t.u. even at the unlock scale 0.2), but a uniform ring
leaks uniformly, det stays 0 pairwise, and no motor appears. Winding is
protected by desert + detune symmetry; slip risk concentrates where symmetry
breaks (unequal voices, rim contacts, strikes). (t5-4.)

**Leak floor (F1).** A rung voice bleeds through every candidate link to the
vacuum choir:

    λ_voice ≈ n_rim · k_dep·k_dep_m · geo · gpl_rim · res(Δω) · ⟨g⟩ · sqrt(Em·mob_floor·cap)
    Δω = w2·q·x/(1+q·x),  res = Γ_m²/(Γ_m²+Δω²);  rough share R(Δω)·rough_k radiates

An upper-bound family (measured ring6 sits ×3.7 below its low band —
recruitment ⟨g⟩, crumb headroom and recondensation all overestimated;
comp12's late rate matches within ~2× at its late Em). The load-bearing part
is the **scaling law** `λ ∝ sqrt(x)·res(Δω(x))`: ring6 : comp12 : diode8 =
6.6 : 1 : 0.5 per unit mass — heavy rungs are quieter against the vacuum in
the resonance factor alone, on top of being farther from the skirt. The
floor is nonzero while any candidate rim link exists (mob_floor·cap = 0.01
guarantees a sqrt(Em·0.01) trickle). Flux-moment instrument prediction
(MASS open item): loop-tangential lem imbalance = N(σ_f−σ_b)d/2C — comp12:
**2.15**, comp6: 0.23, unwound rings: 0 (the chirality signature). (t5-5.)

---

## 6. DESIGNER'S CHECKLIST (the spec for formfind + morpho)

Component-level invariants:

1. **One ω per component** (T1); every voice at x = (w2/ω−1)/q; design floor
   x ≥ 2·x_skirt = **0.123** (skirt death at 0.0617); ceiling x ≤ 0.9
   (headroom 1−x·cap/(cap) vanishes; evaporation door at cap).
2. **Struts** (mutual, both-gate): only in x ∈ (0.123, 0.41), all lengths
   equal d = π/ω ∈ (1.24, 1.62); antiphase; **1:1 strut subgraphs must be
   bipartite** (T11). Length tolerance per strut: ±3.65 % of d (tree) /
   ±7.3 % (annealed) for gates ≥ 0.9; ±1.13 %/±2.26 % for ≥ 0.99 (T4).
3. **Cables** (one-way): φ = ωd ∈ B(ω) (t1-B table); back gate = gate(2φ)
   exactly; **diode point φ = π/2 (back ≡ 0)** at x ≥ 0.705, d = π/2ω ∈
   [1.0, 1.127]. Odd cable cycles are fine (wound); every cycle needs
   Σ±φ = 2πm (T6); balanced cycles (Λ=0) are free.
4. **ω/m recipe** (rings and belts): pick the loop, then ω = 2πm·C/L_loop
   from the ACTUAL picked loop (`ring_m` practice); accept seeds only with
   closure integer to < 0.05 (standing rule). Ring (N, m) feasibility table
   in t2-D; recommended: N=12 m=5 (back 0.10), N=12 m=4 (back 1.5e−5),
   **diode-16 m=4** (back 0, x ≈ 0.83, d ≈ 1.08, parasite margin +0.05;
   diode-8 and diode-12 have parasite exposure — rule 9).
5. **Mixed cycles**: a cycle with n_s struts needs its cables to pay π·n_s
   (mod 2π): triangles pay only n_s = 0 (all-cable, ω ≤ 2.094, φ = 2π/3
   each); pentagons n_s ≤ 2; heptagons n_s ≤ 4 (T12). Tube designs (verified
   phases, t2-D): tri-tube 3 struts + 6 cables (ω ∈ [1.94, 2.09]); pent-tube
   5 + 10 (m=2 belts); hex prism and cube all-strut.
6. **Interval links**: only 2:1 and 3:2 fit the pitch window (ratio cap
   ω_max/ω_min = 1.94–2.05); they carry 1/(pq) exchange and Γ_m/(pq)
   acceptance (fifth: ×6 penalties, tongue 0.0167). Even-m interval rungs
   close odd cycles (fifth-triangle card, T13) — species axis, not a
   stability axis.
7. **Seeding correction (T15)**: seed Em = x*·cap − flload_eq (t5-1 table;
   comp12: seed x 0.233 for rung 0.321) or accept the F1 shed of up to 27 %.
   Prefer N ≥ 8 rings for gpl (N=6 pays ×9 conductance; N ≤ 5 dead at seed
   unless axes pre-set to the chord bisector).
8. **Commas**: per cycle, the LS/annealed state distributes δ_c ∝ 1/w_l
   (T9); keep |δ_c| ≤ N·0.229 for all gates ≥ 0.9; soft links attract comma
   — usable as sacrificial absorbers, dangerous as accidents (the seam
   effect).
9. **Parasites** (unintended candidate links): any seeded pair within
   d < 1.15·(cr0_i+cr0_j) (worst case 2.073) exchanges. Rule: non-edge
   pair distance ≥ 2.08, or design the parasite ONTO a rung/dead phase.
   Ring second-neighbour chord = 2d·cos(π/N): safe for N=6 at d ≥ 1.20,
   N=12 at d ≥ 1.08; cube face diagonal √2·a safe for a ≥ 1.47. comp6 at
   d = 1.10 carries six half-open off-rung parasites (gate 0.53) — a
   candidate cause for its death at 3836 vs comp12's survival. `init=9`
   scores all of this (`# NETGATE P` rows must show g < 0.05).
10. **Expected leak floor**: λ_voice from F1 with the class table (t5-5);
    design target: minimize sqrt(x)·res(Δω(x)) (go heavy), phase-decouple
    the rim, keep every internal link on-rung (det ≡ 0 ⇒ zero internal
    roughness — the M0 52 % channel closed by construction).
11. **Perturbation budget**: relock rate of the slowest ring mode =
    (R_f+R_b)(1−cos 2π/N) (t5-3 numbers); capture |δω| ≤ min(0.3·R₀, 0.2).
    Keep strike transients below capture or the link walks the desert.
12. **Interchange**: nets go to the kernel as `V x y z x_load th2` +
    `E a b` files (`init=9`, `net_file=`); the kernel is the referee via
    `# NETGATE` (intended) and `# NETGATE P` (parasites); acceptance:
    intended g_f ≥ 0.99 (tree exact), parasites < 0.05, closure < 0.05.

---

## 7. What this predicts (pre-registrable, before any run)

P-A. **comp12 retune**: seeding x = 0.233 (not 0.321) lands the ring on-rung
   with no F1 shed; early census M drop shrinks by ≈ the 27 % flight share.
P-B. **diode-16** (m=4, x ≈ 0.83, d ≈ 1.08): back gates exactly 0, highest
   sqrt(x)·res quietness — predicted to beat comp12's lifetime per unit
   mass; its flux moment ≈ N σ_f d/2 with σ_f from t5-1.
P-C. **fifth-triangle**: forms at ω0 ≈ 2.35, holds ~fifth-lifetime (~20–100
   t.u.), dies by interval drift — the first deliberately non-bipartite
   object.
P-D. **flux-moment**: comp12 ≈ 2.15 loop-tangential lem sum; unwound ≈ 0 —
   the chirality instrument's calibration numbers.
P-E. **hex prism all-strut** (P4 corrected): feasible iff the foam yields 18
   near-equal picks within ±7 % (annealed) — a direct P3-mining target.

---

## 8. Corrections to the brief / ledger (explicit)

1. "th2 advances at exactly ω(x) per step" — true of pass 6, but receiver
   clocks are ALSO moved by entrainment at deposit and delivery (κ_lock·mix·
   err, both sites verified); alignment never touches clocks.
2. Practical pitch window: exact [1.394, 2.700] for x ∈ [0.0617, 0.9]
   (brief's ~(1.36, 2.75) was loose); x = 1.0 hard floor gives 1.318.
3. Strut existence is contact-limited well below the candidate cut: x ≤
   0.41, d ≤ 1.62 (not the full d ≤ 2.25 of the bare ladder).
4. "3-regular shell over-constrained by exactly 2 ⇒ demote 2" (ledger P4) —
   refuted: hex prism needs 0; triangles cannot be paid by 1–2 cables;
   tri-prism needs 6, pentagonal 10; K4 has no realization at all.
5. "Interval rungs cannot beat bipartiteness" (expected) — refuted: even-m
   interval rungs close odd cycles in-window (fifth-triangle); 1:1-only
   networks remain bipartite.
6. B1's "gate ≥ 0.9 for |ψ| ≤ 0.16" → gate(0.16) = 0.95; the 0.9 point is
   0.229. B1's desert number 4e−10 at ψ ≈ 2.6 confirmed (4.05e−10 at 2.618).
7. The naive leak-floor prefactor over-predicts ring6 by ~4× (upper-bound
   family); the scaling law λ ∝ sqrt(x)·res(Δω(x)) is the robust content.
8. Ledger translation table: "prestress magnitude ↔ seeded load x + winding
   m" — refine: the true prestress analog is the circulating flux σ (the
   self-stress); x and m are its potentials.

*Scripts: `theory/t1_gates_windows.py` (gates, windows, tolerances),
`t2_cycles.py` (payment lemma, shells, rings, tunability),
`t3_laplacian.py` (Q = BWBᵀ, 1/w law, jittered cube, Z2),
`t4_oddcycles.py` (interval odd cycles, strobe refutation),
`t5_flux.py` (flload fixed point, gyroscopic, desert, leak floor).
Each has a `.out` transcript; all pure math/numpy — no simulator runs.*
