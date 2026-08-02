# CHARGE — the charge construction (effects-first design)

**Design document, opened 2026-07-31 on user directive:** quantify what
charge is in QFT **in effect** (not ontology), map the effects into this
theory, then provide options for the ontology and the numerical means.
Names do not matter ("charge", "spin" are labels); the EFFECTS are the
requirements. Motivation: protons and neutrons contain charged
substructure — **MASS cannot be completed for nucleon-class objects
without the charge construction** (the p–n mass split is partly an
electromagnetic self-energy effect in reality; see E10).

Subordinate to `PRINCIPLE.md`. Standing constraint 3 binds hard here:
**no imported field** — charge must be a property of the existing modes
(space / dense / field) or the construction fails. Ratchet rule governs
any kernel/law change. Correspondence targets below are STRUCTURAL
(CONCEPT-doc discipline), not numerical fits.

---

## 1. What charge DOES in QFT — the effect catalog

Operational content only; each row is something measured about reality.

| id | effect | quantitative content |
|---|---|---|
| **E1** | **Exact additive conservation** | never violated in any observed process; additive over subsystems; both signs; net ≈ 0 universe-wide |
| **E2** | **Quantization + universality** | free particles carry integer multiples of e; \|q_p + q_e\| < 10⁻²¹ e despite utterly different internals (composite vs point); constituents carry thirds (±e/3, ±2e/3) and are never free |
| **E3** | **Gauss law** | flux of the sourced field through ANY closed surface counts enclosed charge exactly — charge is measurable from the boundary without local inspection; field lines end only on charge |
| **E4** | **Coulomb response, both signs** | 1/r² force, attractive AND repulsive (unlike gravity); superposition; strength α = e²/4πε₀ħc = 1/137.036 |
| **E5** | **Holonomy (Aharonov–Bohm)** | interference shifts by qΦ/ħ around enclosed flux with zero local field; flux quantization h/2e in superconductors; Dirac: monopole existence ⇒ qg = nħ/2 ties quantization to topology |
| **E6** | **Magnetic companion under motion** | F = q(E + v×B); magnetism is the frame-partner of Coulomb — comes with relativity, not separately |
| **E7** | **Radiation at the vertex** | accelerated charge radiates (Larmor P = q²a²/6πε₀c³); q sets emission/absorption amplitude; pair production/annihilation strictly ± paired |
| **E8** | **Spin–charge interplay** | moment μ = g(q/2m)S, g ≈ 2 for point fermions (+α/2π = 0.00116 vertex correction — the g−2 program we just archived measures this to 127 ppb); spin-statistics: half-integer spin ⇒ exclusion ⇒ shell structure of matter |
| **E9** | **Neutral-with-substructure** | the neutron: q_n < 10⁻²⁰ e yet μ_n = −1.913 μ_N ≠ 0 and ⟨r²⟩_ch = −0.116 fm² (positive core, negative skin) — internal charge distribution with exact global cancellation; proton r_ch = 0.841 fm; DIS sees point-like thirds with weight ratio 4/9 : 1/9 |
| **E10** | **Charge contributes to mass** | M_n − M_p = 1.293 MeV decomposes (lattice QCD+QED, BMW-class) ≈ +2.5 MeV quark-flavor − 1.0 MeV electromagnetic — the charge texture's self-energy is a load-bearing part of the nucleon mass budget |
| **E11** | **Screening / running** | vacuum polarization: α(0)=1/137 → ≈1/128 at M_Z; IR charge universal and exact (Ward identity) — the far-measured charge is protected against internal complexity |

**Load-bearing subset for the MASS program** (nucleon-class): E1–E4
(conserved, quantized, Gauss-countable, ± force), E9 (neutral objects
with internal charge structure), E10 (charge self-energy in the mass
books). E5 is the cheapest decisive TEST once anything exists. E6
arrives free with the Maxwell/relativity completion (EMF.md 2026-07-31
entry). E8 spin-statistics is noted honestly as the hardest distal item
(the HOM exchange-registry machinery is the only statistics foothold).

## 2. Mapping the effects into v89 — what each requires

| effect | v89 requirement | exists today? |
|---|---|---|
| E1 conservation | a label that survives all conversions — i.e. TOPOLOGICAL, not an amount of stuff | winding machinery exists (two-plane phase, closure integers, gate desert barrier ~4e-10) |
| E2 quantization/universality | integers by construction ⇒ winding/wrapping numbers; composite charge = sum of windings EXACTLY (homotopy is additive) — this solves \|q_p+q_e\| exactness structurally, which no "amount" ontology can | closure-integer machinery exists |
| E3 Gauss | the structure must imprint a FAR TEXTURE whose closed-surface index equals the label — an index theorem, not a dynamical accident | **missing** — needs polarization-structured field (EM5) |
| E4 Coulomb ± | texture overlap energy must depend on relative orientation: like-windings cost, unlike-windings cancel | missing; falls out of texture energetics if E3 lands |
| E5 AB | transported packets acquire geometric phase around the texture | interference instruments EXIST (slit/eraser) — cheap test |
| E6 magnetism | Lorentz-covariant field sector | EM5 + the lock-Lorentz result |
| E7 vertex | conversion door already fires in atoms; charge modulates door amplitude | partial (door exists; coupling absent) |
| E9 neutral substructure | internal ± distribution with zero total — a closed machine with balanced handedness books but nonzero handedness MOMENTS | flux-machine picture ready to host it |
| E10 mass contribution | texture self-energy enters the balance books | follows once E3/E4 exist |

## 3. Ontology options (ranked; O1/O2 may be one thing seen from two sides)

### O1 — Poincaré-texture charge (field-side: the hedgehog)

The field chiral pair (u₊, u₋) has, at each place, a normalized state =
a point on the **Poincaré sphere S²** (polarization Bloch vector). A
point defect whose surrounding sphere is WRAPPED once by that vector
(hedgehog) carries π₂(S²) = ℤ:

* **Gauss = index theorem**: the wrapping number through any enclosing
  surface is the same integer — E3 exact by topology, surface-
  independent by construction.
* **± = orientation** (inward/outward hedgehog); C-conjugation = swap
  the planes (handedness flip).
* **Annihilation/pair production**: opposite wrappings sum to zero and
  can decay to pure radiation; creation is strictly paired — E7's ±
  pairing automatic.
* **AB analog**: a packet transported around/past the texture acquires
  a Pancharatnam–Berry phase — E5 as a fringe-shift experiment on the
  existing slit machinery.
* **Coulomb sign structure**: two like hedgehogs = total wrapping 2
  (higher texture energy → repel); unlike = 0 (relax → attract). The
  1/r² falloff must be MEASURED (it is the far-texture gradient energy
  under the Maxwell cone; on a band it will be wrong — EM5 first).
* Needs: EM5 (link-borne chiral pair — the polarization frame is the
  link's own transverse plane pair, PRINCIPLE §3) + a texture seeder.

### O2 — Helicity-throughput charge (dense-side: the chiral machine)

The flux-machine principle (MASS.md 2026-07-31): a particle is
structurally forced intake = outtake. Those books have a SCALAR part
(energy — mass) and a HANDEDNESS part (which chirality flows in vs
out). **Charge = the conserved chiral imbalance of the throughput
books**: a machine that cycles both handednesses symmetrically is
neutral; one with net handedness throughput is charged.

* Sign = circulation sense (the chiral pump's one-way ν already
  measured in the comp12 class; the B1 chirality split m = N/2 ∓ w was
  pre-registered as "the first charge-like species axis").
* Neutral-with-substructure (E9) is NATIVE: a globally balanced
  machine with internal circulating handedness has zero charge but a
  nonzero handedness MOMENT — the neutron's μ ≠ 0, ⟨r²⟩_ch < 0 pattern.
* Gauss: a steady handedness source must show as a conserved helicity
  flux through enclosing shells (the field-side face of the same
  books — likely O1's texture viewed from the source). DERIVATION
  OWED: helicity continuity in the two-plane field sector.
* Spin: the machine's circulation integer is its spin candidate; a
  charged spinning machine is automatically a helicity dipole — the
  moment μ ∝ q·(circulation) shape of E8, with g computable once real.

### O3 — Fractional shares on odd closed cycles (the substructure pattern)

Inside either ontology: a closed 3-voice cycle can carry total winding
2π distributed as per-voice shares (2π/3 each, or 2·2π/3 + ... mixed by
interval content). **Thirds exist only inside closed triples; only
integer totals close freely** — the quark pattern (±1/3, ±2/3 confined;
integers asymptotic) as a closure statement, using machinery that
already exists (odd cycles closed by even-m interval rungs — the
fifth-triangle). Nucleon = 3-voice machine; p vs n = same topology,
different interval/handedness content; their mass split = different
texture self-energies (E10).

**→ Deepened into the working ontology O3′ (§7, user direction
2026-07-31): charge = the un-closed tail call.** O3′ subsumes O1 and O2
(the deficit forces the handed throughput and dresses itself with the
far texture) and moves a large part of the program off the EM5
critical path.

## 4. Numerical means

1. **Prerequisite: EM5** (link-borne chiral pair + curl coupling —
   EMF.md 2026-07-31). Charge lives in polarization structure; today's
   per-species scalar amplitudes cannot hold a hedgehog. This is the
   concrete content of "MASS is blocked on the EMF charge
   construction": the dependency is real in-theory, not just
   phenomenological.
2. **Substrate**: charge is topology; topology needs quiet geometry.
   Run charge work on the annealed substrate (S1, MASS.md 2026-07-31)
   — foam jitter at σ_d ≈ 9% would shred textures with defect noise.
3. **Seeders**: (a) hedgehog seeder — wind the Poincaré vector over
   shells around a cell; (b) helicity-pump seeder — chiral ring with
   handedness-asymmetric gates (diode machinery + the two planes).
4. **Instruments**: (a) **Gauss meter** — per-shell wrapping-number /
   helicity-flux counter (the g4 radial machinery pattern applied to
   Stokes parameters of the field pair); (b) Berry interferometer —
   slit apparatus with a texture between the paths, score fringe
   SHIFT (E5); (c) two-charge force rig — centroid drift vs
   separation, sign and power law (E4); (d) annihilation cell —
   opposite textures released, books → 0 + radiation burst (E7).
5. **Integer charge ledger**: winding numbers join the integer books —
   per-structure charge account, conservation exact by construction
   (the max_sum_err=0 standard extended to charge).

## 5. Experiment ladder (all pre-registerable)

| id | question | bar |
|---|---|---|
| Q1 | seeded texture survives and its shell index is INTEGER | wrapping stable, same integer every shell |
| Q2 | Gauss surface-independence | index identical across nested/deformed surfaces, exact |
| Q3 | AB analog | fringe shift ∝ enclosed winding with zero local field energy on the paths |
| Q4 | two-charge force | like repels / unlike attracts; power law measured (1/r² only expected under EM5 cone) |
| Q5 | annihilation | ± pair → zero books + radiated energy conserved to ledger floor |
| Q6 | composite additivity | assembled structure's far index = sum of parts, exact |
| Q7 | neutral substructure | balanced machine: far index 0, near-field handedness moment ≠ 0 (neutron pattern) |
| Q8 | charge self-energy in the mass books | two same-topology machines differing only in handedness content have systematically split balance points (p–n split analog, structural) |

## 6. Sequencing

Substrate S1 (cheap, unblocks everything) → EM5 design + variant table
(doubly motivated: Maxwell cone AND polarization substrate — law event,
scheduled per EMF.md §5 doctrine) → seeders + Gauss meter (apparatus,
buildable against the variant table) → Q1–Q5 → then the nucleon program
(O3 substructure, Q6–Q8) rejoins MASS. The flux-machine balance work
(MASS) and this construction share one bookkeeping: mass = the scalar
balance, charge = the conserved chiral asymmetry, spin = the
circulation integer of the same books.

**Amended by §7:** the slip sector (Q9–Q11) runs on the STANDING table
with a load-pinning apparatus only — before EM5. Only the far dressing
(Gauss meter, Coulomb power law) waits for the cone.

---

## 7. O3′ — the deficit ontology: charge = the un-closed tail call
### (user direction, 2026-07-31; symbolic numerics verified this date)

**The proposition (user):** perfection of tail calling is not
universal. Some structures are intrinsically OFF closure — and that
deficit IS charge. The electron in particular is not a closed machine
at all: it is the dispersed harmonic field generated by a tail call
that is never fully adhered to — "electrons aren't really electrons,
they are the field itself of incomplete harmonics." Electric current
is the flow of the un-tail-called harmonics themselves, and it is
chiral because the residue of an incomplete call lives in the two
handed planes.

### 7.1 The holonomy algebra (charge quantization derived)

An interval lock p:q is a **multivalued circle map**: the locked
relation q·θ_j − p·θ_i = α_l defines θ_j only up to a branch,
θ_j = (p·θ_i + α_l + 2πn_l)/q_l, n_l ∈ ℤ_{q_l}. A closed cycle needs
frequency closure Π(p_l/q_l) = 1 (exact rational arithmetic — the comb).
Composing the maps around the cycle gives the monodromy. Verified
symbolically for the third-harmonic triangle {3:1, 1:3, 1:1}
(comb_limit=6 admits it):

    M(θ) − θ = (α₁+α₂)/3 + α₃  +  2πn/3 ,   n ∈ ℤ₃

Two parts, utterly different in kind: a **continuous comma** (the α
sum — annealable by plasticity/tempering, dumped on seams, payable in
atoms: all measured machinery) and a **discrete branch 2πk/Q** that no
continuous deformation can remove. **CHARGE := k/Q — the ℤ_Q holonomy
class of the cycle.** Conserved topologically (changing k means
re-deriving the cycle through the gate desert), quantized in 1/Q,
additive over composites, signed by traversal orientation (= the
handedness of the planes: C-conjugation is the mirror swap).

The comb's small integers generate the small quantum numbers: 3:1
content ⇒ ℤ₃ ⇒ **thirds** (the quark pattern); 2:1 (octave) content ⇒
ℤ₂ ⇒ **double-valued closure — a 4π-return structure, the spinor
pattern** (nomenclature per user: we care about the effect; "spin-½"
is an octave holonomy).

**Confinement theorem (structural):** the vacuum is unison — Q = 1.
A fractional branch is defined only relative to interval content; only
integer holonomy is transportable through a unison medium. Fractions
are confined to interval-rich structures (the hadron pattern); integers
roam (the lepton pattern). E2's free-charge integrality is derived, not
imposed.

### 7.2 The kink realization (the slip soliton in the standing gate law)

The locked chain's energy is E[{θ}] = Σ_l w_l·V(ψ_l) with
V(ψ) = 1 − ((1+cos ψ)/2)^{p_gate} — **the standing gate law is a
Frenkel–Kontorova / sine-Gordon system**, and a cycle with boundary
twist 2πk/Q carries its charge as k discrete **phase-slip solitons
(kinks)**: localized on links, mobile, chiral (kink vs antikink = the
plane handedness), pair-created and annihilated in ± pairs. Core
energy E = √(2A)·I(p_gate)/Q with I(1) = 4.0000 (sine-Gordon exact),
**I(8) = 5.405** (the standing p_gate) — fractional cores cost 1/Q of
an integer core; the far dressing adds ~q². Max gate torque at p=8 is
**1.233·ε at ψ = 0.505** (ε = coupling from the measured tongue
widths, E6 — parameter-free).

> **Scope note (2026-08-02).** Of this paragraph, three things are
> established and one is not. Established: I(1)=4 and I(2) exactly, I(8)
> by quadrature; the **integer** core energy |k|·√(2A)·I(p), measured to
> a known lattice correction (§8.1); and the max-torque constant, now a
> closed form — cos ψ\* = (p−1)/p, T_max = 1.232578445 at p=8 (§8.1).
> **Not established: the 1/Q reduction of the core cost for a fractional
> charge.** It is a structural conjecture. The scan that appeared to
> confirm it was run on a unison ring, where a 2πk/Q twist is not a
> soliton, and it is withdrawn (§8.4). Testing it needs a ℤ_Q interval
> chain. Likewise "the far dressing adds ~q²" is asserted by continuum
> analogy; no artifact measures any far dressing.

### 7.3 The electron

The dressed **delocalized** kink: no dense core, no closed loop of its
own. Its substance is (i) the slip and (ii) the far phase dressing
∮∇χ·dl = 2πk/Q the slip forces onto the surrounding field — which is
O1's texture (Gauss = the dressing's index) and O2's handed throughput
(the residue is chiral), so O3′ generates both. Its mass is dressing +
core energy — field-gradient scale, naturally light against a
nucleon's dense closure energy (the mass-hierarchy direction). Its
position is a collective coordinate dispersing through the field band;
it localizes only at conversion events — and v89 has ALREADY measured
that machinery end to end (claim rule R3, MCWF-exact clicks, Tonomura
fringe statistics): **the measured double-slit sector IS electron
phenomenology under this ontology** — "statistical fields where they
may travel" is the demonstrated click statistics. (Convergence note:
the v85-era conclusion "no electron soliton; shells = bound response
harmonics" recurs here derived, v89-native: shells = bound harmonics
of the deficit field around a charged core.)

Coulomb sign in three lines: dressing energy ∝ ∫|∇χ₁ + ∇χ₂|² =
self₁ + self₂ + 2∫∇χ₁·∇χ₂; the cross term carries sign q₁q₂ — like
charges repel, unlike attract. The power law (1/r²) needs the EM5 cone
(on the measured band it screens wrong) — that is the ONLY part of
this sector gated on EM5.

### 7.4 The electric sector (derived correspondences)

* **Voltage = pitch differential** (load flattens pitch: the load
  landscape IS the potential landscape).
* **Josephson relation, exact:** across a weak link, ψ̇ = Δω ⇒
  asymptotic slip rate ν = Δω/2π — one slip per 2π of accumulated
  phase, torque-independent (verified numerically). Pitch difference ↔
  slip current, with A₀ in ħ's role.
* **Superconducting branch = the lock:** no slips below the critical
  detune Δω_c = 1.233·ε (the tongue edge — the Arnold tongue IS the
  critical current, already measured as E6/E7 apparatus).
* **The I–V characteristic (computed, parameter-free once ε is read
  off the tongue):** ν = 0 below Δω_c; √(Δω²−Δω_c²)/2π near onset
  (saddle-node universality — verified to 4 digits); → Δω/2π far.
  Gate-torque (p=8) curve tabulated 2026-07-31.
* **Ohm/Joule = the roughness toll:** each slip pays the R-window
  through the conversion door in whole atoms ⇒ dissipation is
  quantized per slip, P = ν_slip·ε_rough. Resistance is the per-slip
  roughness payment; dissipationless flow below threshold.
* **Persistent current = ring holonomy** (the ring_m circulation,
  already standing).
* **Magnetism = the helicity current of handed slips:** the chiral
  part of the existing momentum-current instrument,
  J^χ_ij = −2J_hop·Im[u₊ᵢ*u₊ⱼ − u₋ᵢ*u₋ⱼ]; the right-hand rule is the
  plane-handedness convention. Correspondence bar: the screw-locked
  slip's μ/L ratio (g ≈ 2 for the point pattern).
* HBAR tie (flagged, not claimed): the holonomy ledger is integer
  harmonic content — HBAR.md's strongest F-ħ candidate; per-slip
  action vs A₀/Q is measurable in the QATOM books.

### 7.5 Ladder additions (Q9–Q12)

| id | question | bar | needs |
|---|---|---|---|
| **Q9** | weak-link I–V: two load-pinned reservoirs, one bridging link; ν_slip(Δω) | locked below Δω_c; √ onset; Δω/2π asymptote; all parameter-free from the measured tongue | load-pinning apparatus only — STANDING TABLE |
| **Q10** | slip chirality: single-slip handedness in the two planes; kink↔antikink under plane swap | slips are handed; mirror swap flips sign | standing table |
| **Q11** | fractional confinement: export a single ⅓ slip from a ℤ₃ cycle into unison vacuum | fails or pair-creates; integer composite exports cleanly | standing table + ℤ₃ cycle seeder |
| **Q12** | dressed-kink interference: slip texture through the two-slit machinery | click statistics with conserved holonomy — the in-sim "electron double-slit" | after Q9–Q11 |

Q9–Q11 precede and de-risk EM5. The Gauss meter and Coulomb power law
(Q1–Q5) remain the EM5-gated half.

### 7.6 Q9 EXECUTED — the v89 weak link measured (2026-07-31)

**Apparatus (kernel, ratchet-gated, battery 21/21 with everything
default-off):** `pin_net` — holds init=net voices' Em at their seeded
targets via ONE reservoir account inside the conservation sum (exact
int64 moves; `# PIN`/`# RESULT pin`); `slip_diag` — per-step unwrapped
forward-gate phase over net edges + the integral of the bare pitch
difference (`# SLIP`/`# RESULT slip`). DICHOTOMY GUARD honored: no
charge state anywhere — energy books + derived phase counts only.
Harness: `charge/q9_gen.py` (bridge pick + symmetric detune at exact
rung SUM: ω_a+ω_b = 2πC/d — Δω is a pure difference-mode drive),
`q9_run.py`, `q9_score.py` → `runs/q9_iv.png`.

**Bridge:** foam 20260727 cells 6395↔6928, d = 1.4080 (the m=1 rung
window), seed gates 1.0000/1.0000, zero parasites. x̄ = 0.2498; scan
Δω_seed 0→0.55 (15 points, T=600). Conservation: `max_sum_err=0`,
drift ~1e-16 on every run with active pinning; the reservoir drains at
~1.6e-3 Em/t.u. at Δω=0.21 — the driven junction's dissipation, fully
booked (the pin reservoir IS the voltage source's energy account).

**Results — all three pre-registered bars PASS:**

| bar | verdict |
|---|---|
| B1 locked branch | **PASS** — 6 points Δω_meas ≤ 0.1210 all locked: ν ≤ 4.1e-5 (noise floor), slips = 0. The dissipationless branch exists. |
| B2 knee | **PASS** — **Δω_c = 0.127 ± 0.006** (locked 0.1210 / running 0.1329). Onset STEEPER than sine-Adler (48% of the line at +5% past the edge, vs 22% for sine) — the measured signature of the p_gate=8 torque shape (§7.2: max 1.233·ε at ψ=0.505); sine-form extrapolation underestimates the edge (0.085), as the gate torque predicts. |
| B3 Josephson line | **PASS** — far branch ν → Δω_meas/2π: point ratio 0.998 at the top; far-branch slope fit 1.014 (n=7, R²=0.9997; 9-point fit incl. knee 1.028). Parameter-free, torque-independent — the law-owned prediction. |

Bonus (unplanned but measured): the **current–phase relation** — on
the locked branch the antisymmetric deflection ψ_a = (ψ_f−ψ_b)/2 grows
0 → 0.036 → 0.082 → 0.154 → 0.289 as Δω → the edge (gates easing
0.95 → 0.79/0.90): the static torque curve traced point-by-point, the
junction's CPR. Also measured: the flight-load correction (Δω_meas
1–8% below seeded, converging at large detune — flight-loads-pitch
visible in the instrument), and a common-mode ψ offset (+0.152 at
Δω=0) from the vacuum environment's asymmetry.

**Foam-owned vs law-owned (substrate directive, user 2026-07-31 —
"we may need to evolve away from the foam cells prior to EMF and
charge experiments as well, or full success"):** Q9's own decomposition
makes the split measurable. LAW-OWNED (foam-independent): the
locked/running dichotomy, the saddle-node onset class, and the
Josephson slope = 1 (measured 1.4%). FOAM-OWNED (this link's
accidents): Δω_c, the CPR amplitude (the coupling ε through contact
area/degree/local Es — expect ±30% link-to-link), and the ψ baseline
offset. TIMELINE CONSEQUENCE: single-link probes (Q9, Q10) are
meaningful on the frozen foam because their law-owned content
separates cleanly; **Q11 (ℤ₃ interval cycles) and the EM5 texture work
are NOT foam-safe** — fractional holonomy needs interval rungs whose
tolerance windows the foam's ±30% chaos eats (the fifth-bridge
deferral generalizes). So the substrate evolution (S1 at minimum)
precedes Q11 and EM5-texture experiments; livefab (S3) precedes any
"full success" claim. Universal-Δω_c claims would need multi-link
ensembles even for Q9-class runs.

Next on the ladder: Q10 (slip chirality under plane swap — same
apparatus, add the two-plane handedness readout) on the current foam;
Q11 waits for S1.

### 7.7 Resumption results (2026-07-31, post-S1)

**Q10-lite — C-symmetry of the slip sector: EXACT.** The weak-link I–V
is antisymmetric under voice swap (kink ↔ antikink): ν(−Δω) = −ν(Δω)
measured at three points — running branch −0.032043 vs +0.032046
(1:10⁴), −0.063162 vs +0.063157 (8:10⁶), locked branch symmetric
(|ν| ≤ 3e-5 both signs at |Δω|=0.056). No built-in rectification;
charge conjugation is an exact symmetry of the transport law. (Full
Q10 — the two-plane handedness registration of single slips — still
open; needs a per-slip plane-resolved counter.)

**Q11 status:** substrate-unblocked (S1 green) but seeder-blocked —
the .net format carries no interval column (E a b defaults 1:1), so a
ℤ₃ cycle (3:1/1:3/1:1 content) cannot yet be placed. Next apparatus
item: per-edge `E a b p q` extension (ratchet-gated), then Q11 runs on
V2s.

**M-B1 cross-link:** the balance-curve instrument (MASS.md this date)
measured the pair drain floor B_min = 8.85e-4 ≈ 2·c₀·cap — the same
mob_floor bleed that sets Q9's dissipation books; one microscopic
channel now measured from both programs' sides.

### 7.8 Q11 first light (2026-07-31): the circulator, and the instrument gap

The fifth-triangle (P17 address: cells 3928/4407/5902, sides
1.33–1.38, phase-chain closure residual 0.11 rad) placed via the new
interval-net edges (`E a b p q`, battery-gated). Unpinned: loads drift
off 3:2, the triangle runs as a phase-slip circulator whose books
close like Kirchhoff (ν₀ = +0.15822 vs ν₁+ν₂ = 0.15849 — 0.2%).
Pinned: a clean two-edge circulator (ν₀ = +0.121 ≈ −ν₁; the 1:1 edge
SILENT at ν₂ = 0.005 — A and C hold coherence across the cycle).
**Locked-vs-free on the 3:2 edges is undecidable with the current
instrument**: slip_track uses the unison phase, and a locked p:q edge
advances that variable exactly like a free one. Required for Q11 v2
(small, ratchet-gated): comb-aware slip books — Φ = q·θᵢ − p·θⱼ from
the link's lp/lq — plus the same in the # NETG gate report. Then the
ℤ₆ holonomy count is readable directly.

### 7.9 Q11 v2 (2026-07-31): the dual-book instrument decides — partial
### lock, comb-chart circulation, and the tongue-chatter discovery

The §7.8 gap is closed, in two instrument steps that each taught
physics (battery-gated; q9 unison regression EXACT at every step —
all archived digits reproduced, fixed book ≡ live book on 1:1 nets).

**Step 1 (comb-aware live books) exposed the architecture:** lp/lq is
not a register but the comb DETECTOR's live state, re-derived each
step from actual pitches (cellfab.c C1 comb scan). Even the PINNED
fifth-triangle ends at pq=1:1 — and the final instrument counts why:
**flips=867/865 per interval edge over T=600 (1.4/t.u.)** — the
circulator's own flight load detunes the effective pitches and the
3:2 tongue (width γ/6, the Tenney-height narrowing) is marginal
against it, so the detector chatters. A slip sum over a chattering
chart is not a holonomy count (live books read ±41 "slips" — chart
mixing, not physics). FINDING, stated once: **interval locks are
fragile under their own transport — complexity narrows tongues, so
higher intervals are hard to HOLD, not just to reach.** (Unison, with
the full-width tongue, is what the corpus locks on — same reason.)

**Step 2 (dual books): slip_diag now tracks BOTH** the live chart
(kernel-faithful, flip-counted) and the net's DECLARED p:q chart
(frozen at seed, orientation-matched to the link) — clean ℤ-countable
holonomy in one chart plus visibility of detector state. Measurement
only; battery 21/21 (gate 2) and rerun for gate 3 with the dual books.

**The fifth-triangle verdict (fixed 3:2/2:3 chart, window [300,600]):**

| edge | declared | fix_nu | fix_dwb/2π | ratio | fix_slips | flips |
|---|---|---|---|---|---|---|
| A–B | 3:2 | −0.0144 | −0.0326 | **0.44** | −4.33 | 867 |
| B–C | 2:3 | +0.0145 | +0.0327 | **0.44** | +4.36 | 865 |
| C–A | 1:1 | +0.00005 | +0.00002 | — | +0.01 (LOCKED) | 0 |

(1) **Locked-vs-free is now decidable, and the answer is PARTIAL
LOCK**: pinned interval edges slip at 0.44× their bare comb detune —
intermittent phase-locking at the tongue edge — vs the free arm's
0.84–1.00 (free-running; its detector settles to 1:1 with 7/3/0
flips). The old unison instrument read ±0.117 on these edges: a
chart artifact, 8× the true comb-chart rate. (2) **First measured
ℤ-holonomy circulation in the interval chart**: −4.33/+4.36/+0.01
around the cycle — antisymmetric to 0.7%, Kirchhoff-closed to 0.04
slips, the unison edge holding coherence exactly (both books, 0
flips). The fifth-triangle is a comb-phase circulator. (3) The
residual comb detune ON the pinned address, 0.205 rad/t.u., is the
flight-load detuning measured directly (the pin holds Em; flight
adds to x; the coincidence 2ω_A = 3ω_B shifts). NEXT (the Q11
payload proper, §7.5): the fractional-export experiment — kick a
single ⅓-class slip out of an interval cycle into unison vacuum;
confinement or pair-creation. Needs a transient-kick apparatus
(pin-target step or edge Δω pulse) — small, ratchet-gated, now
designable against a measured baseline. Files:
charge/runs/q11_*_v2.log, scratchpad q9_regress2 (exactness).

### 7.10 Q11h (2026-07-31): the harmonic scan — trapped winding IS the
### charge state; bistability, loop-rate conservation, creep floor

The scan the charge question wanted ("holonomy as a function of the
interval address, and of size"): four pinned cycles, dual-book
instrument, all Kirchhoff-closed ≤0.5 slip (full table:
charge/q11h_score.py; runs q11h_*.log; nets with closure arithmetic
in headers).

| cycle | harmonic (Tenney) | state found | per-edge fix_nu | flips |
|---|---|---|---|---|
| fifth_tri P17 | 3:2 (6) | CIRCULATOR (winding at vertex B) | ∓0.0145 | 867 |
| fifth_tri2 | 3:2 (6) | LOCKED (all edges, live=3:2 held) | ~2e-4 | 0 |
| oct_tri | 2:1 (2) | CIRCULATOR (winding at vertex A, through the unison edge) | ∓0.0145 | 0 |
| fifth_hex (ℤ₆) | 3:2 (6) | CIRCULATOR, two counter-arcs + two quiet edges | ∓0.0063 | 0 |

FINDINGS. (1) **BISTABILITY = the charge state.** The same address
supports a quiet neutral state and a persistent circulator; the
control nailed it: re-seeding the P17 cells with kink-free chained
phases collapses the circulation — live book Σ|slips| 76.4 → 16.0
(4.8×), declared chart ±4.33/4.36 → ≤1.53 (2.9×), e0 locked, chatter
867→107 — with pin/loads/geometry identical. Trapped integer
winding — set by the initial phase texture, conserved thereafter —
is the state variable. In O3′ terms: the un-closed tail-call count
the cycle carries. Neutral = 0, charged = n≠0, sign = circulation
sense (C-symmetric). (2) **The trapped quantum's current divides
with loop size**: triangle Σ|ν| 0.0290 units/t.u. vs hexagon 0.0252
(87% conserved) at 0.43× per running edge — near-conserved total
current, fluxoid-in-a-ring phenomenology; on the hexagon the flow organizes as two counter-arcs
with two quiet edges (texture, not uniform current). (3) **The
loop rate is harmonic-independent** (fifth and octave circulators
both |ν|≈0.0145/edge on triangles) — the winding precession scale is
a property of the cycle+pin machinery, not of the interval; origin
unhunted (candidate scales: γ_res_m=0.10 ≈ 2π·0.0145=0.091 rad/t.u.).
(4) **Frustration sets a creep floor, not the circulation**: the
closure arithmetic (fifth-triangle defect = −ω_h·perimeter mod 2π
with ℤ₃ branch freedom ⇒ no zero-defect fifth triangle exists
in-band; octave triangles and hexagons can close) predicts which
cycles can fully quiesce: tri2 (closable-ish, 0.77) locks dead;
rephased P17 (0.94) creeps at ≤1.5 slips/300 with tongue chatter —
but a zero-defect octave triangle still runs as a full circulator
when seeded wound. (5) **Tongue fragility is Tenney-driven,
cleanly**: under identical transport the octave (γ/2) never flips
its live detector once; the running fifth chatters at 1.4/t.u.;
LOCKED fifths (tri2) also hold 3:2 with zero flips — the chatter is
a consequence of transport (flight-load detune), not of the interval
per se. Combined with M-B4's vacuum-tongue drain spikes (fifth 114×,
octave 11× over the desert floor at their x-coincidences) the
charge picture has: trapped
ℤ winding (state), quantized persistent current (its observable),
tongue-fragility + drain channels (its interactions with the
medium). NEXT: the kick apparatus (pin_kick_*, built, gate 4
running) injects/extracts single quanta — the export/confinement
experiment (Q11 payload) and the winding-number instrument
(per-vertex θ-winding counter) decide fractional confinement.

### 7.10b Kick experiments (2026-07-31): winding conservation
### demonstrated dynamically; the mixed-interval cycle is PUMPABLE

Three single-vertex kick runs (pin_kick apparatus, gate 4; kick =
vertex-B pitch excursion calibrated to one comb-quantum scale,
dv=0.5 Em × 10.33 t.u.): **inject** (locked tri2 + kick) → RE-LOCKS,
±0.06 slips, nothing trapped; **deinject** (P17 circulator −
counter-kick) → circulation UNTOUCHED (fix −4.33/+4.36/+0.02,
identical to the un-kicked run to 3 digits); **frac3** (⅓-scale) →
re-locks clean. WHY, exactly: the loop winding is Σ_cycle Φ̇_e, and
for this net's reciprocal interval pattern (3:2 then 2:3 then 1:1) a
single-vertex excursion telescopes to zero — vertex B enters edge AB
with weight −p_AB = −3 and edge BC with weight +q_BC = +3. **A local
excursion can neither create nor destroy trapped winding: charge
conservation is dynamically protected, the loop-EMF statement of
this medium** (a local potential fluctuation is not a Faraday loop).

THE PUMP (derived, then built): the telescope BREAKS for coordinated
two-vertex excursions on mixed-interval cycles. Kick B and C
simultaneously with Δω_C = (q_BC/p_BC)·Δω_B = 1.5·Δω_B: edge BC's
comb detune (3Δω_B − 2Δω_C) nulls, AB advances at −3Δω_B, CA at
+1.5Δω_B ⇒ net loop injection −1.5·Δω_B/2π per unit time. A UNISON
cycle can never be pumped this way (weights telescope for any vertex
set); the ℤ₃ covering structure (LIT.md novelty #2) is therefore
DYNAMICAL: interval cycles are chargeable by coordinated pumping,
sense = charge sign, and the quantum arrives in τ = 2π/(1.5·|Δω_B|)
= 20.6 t.u. at dv_B=0.5/dv_C=0.333. Apparatus: pin_kick2_* second
window (gate 5 running). Experiments staged: pump1 (+1 quantum
sense), pumpR (reverse sense), pump3 (⅓ pulse — the fractional
confinement question, now properly posed: does a partial pump
relax, or trap a fractional state?).

### 7.10c Pump results (2026-07-31): the THIRD state — pinned
### fractional strain; mobile charge is integer, bound charge is not

The RESULT rate-window read "all locked" — the SLIP series showed
why that was the wrong instrument: the pump WORKED, and what it
creates is static. Final loop windings (Σ unw around the cycle,
t=595, all edges lock-quiet from pulse end to T):

| run | pulse | loop winding | predicted (with ∂ω/∂Em nonlinearity) |
|---|---|---|---|
| pump1 | +1 design, 20.6 t.u. | **+0.883·2π, STATIC** | +0.89 (coefficient softens as x_B rises) |
| pumpR | reverse | **−1.184·2π, STATIC** | −1.15 (coefficient stiffens as x_B falls) |
| pump3 | ⅓ pulse | **+0.236·2π, STATIC** | +0.30 |

The transient is exactly the §7.10b arithmetic (pump1: AB +1.83·2π,
CA −0.92·2π, BC nulled throughout; freeze at pulse end; zero further
slips for 380 t.u.). READINGS: (1) **A third state class exists**:
not neutral, not the running circulator — PINNED STRAIN, a static
multi-branch lock texture holding CONTINUOUS (non-quantized) loop
winding; each gate absorbs displacement within its capture range and
whole branch shifts distribute the rest (the large-inductance-SQUID/
pinned-fractional-vortex analog; LIT.md #5). It is stable because it
is transport-free: no flight load, tongues intact, nothing to slip.
(2) **Fractional charge = bound polarization**: the ⅓ pulse traps
⅓-scale strain indefinitely — it does NOT relax; but nothing
fractional MOVES: mobile circulation (P17/oct circulators) has only
been observed in whole-quantum units with the universal loop rate.
Current picture: neutral / polarized (pinned, continuous) /
conducting (running, integer) — the dielectric–conductor distinction
native to interval cycles. (3) The forward/reverse asymmetry
(0.88 vs 1.18) is pure pitch-law nonlinearity and matches to 3% —
the pump is now a CALIBRATED charge-injection instrument once τ is
corrected by the measured ∂ω/∂Em(x) along the excursion.
INSTRUMENT NOTE (next increment): a final WINDING report row (unw
per edge + loop sum) — the data already lives in # SLIP; the
RESULT-row rate window is blind to static states by construction.
OPEN (sharp): does pinned strain DEPIN into a running circulator at
a threshold (LAMH barrier crossing — drive with transport or piston
noise and watch for ignition); can a pump pulse timed to end
OUTSIDE capture leave a running integer circulator (charge
injection into the conducting state); ℤ₃ branch bookkeeping of
which gates hold which branch (the covering-space census).

### 7.10d Ignition null + the octave reappraisal (2026-07-31, honesty
### pass)

**Ignition null:** ending the pump mid-slip (τ=28.1, loop 1.22·2π,
AB stranded outside capture) does NOT start a circulator — the
system re-pins with more strain (+1.178·2π static, zero slips in
[300,600]). On tri2, every pump endpoint — adiabatic or stranded —
relaxes to PINNED strain; the running branch is not pump-adjacent
(LAMH-consistent: its basin is not reachable by quasi-static
deposition; a genuinely dynamical preparation — moving texture, not
displaced texture — appears required, which is exactly what the P17
historical seed was). **Octave reappraisal:** oct_tri's circulator
SELF-STARTED from a kink-free chained seed — and the band arithmetic
FORCES the octave's light vertices under the skirt (ω-ratio 2.0 in a
2.2-wide band ⇒ x_light ≤ 0.03 < x_skirt = 0.0617; no in-band octave
triangle exists with all vertices above the skirt). Its circulation
is therefore probably SKIRT-DRIVEN (vacuum-dissolution physics at
the empty vertices, pin-sustained), not trapped-winding — the
in-band octave is intrinsically ionization-adjacent, a finding in
itself. §7.10's oct row is reinterpreted accordingly; the winding
claims now rest on the fifth-class pair (P17 runs / tri2 locks /
rephased-P17 quiets), the deinject invariance, and the pump strain
ledger — all of which stand. The identical ±0.0145 loop rate across
P17 and oct now reads as apparatus/geometry-owned (pin protocol),
strengthening the "rate is not interval physics" verdict.

VERIFICATION (2026-08-01): every number in §7.10–§7.10d re-derived
from the run logs in a clean pass. Pump loop windings recomputed
from the # SLIP series: pump1 +0.8825, pumpR −1.1843, pump3
+0.2354, ignite +1.1779 (·2π; tri2 un-pumped baseline −0.043);
P17 baseline = q11_pinned_v2.log (fix ∓0.01444/+0.01453, flips
867/865), deinject matches it to 3–4 digits. Provenance: the kick
and pump apparatus each passed the full battery before use (gates
4 and 5, 21/21 each, all observables byte-identical to the
pre-apparatus kernel with the knobs off), and the q9 weak-link
regression is exact at every kernel state of the arc.

### 7.11 Q-rate (2026-08-01): the loop-rate origin — γ eliminated;
### the fifth circulator rides a γ-free RSJ working point

The γ_res_m scan on the P17 circulator (0.05/0.10/0.20, runs
charge/runs/qrate_g*.log vs q11_pinned_v2): fix_nu = 0.01438 /
0.01444 / 0.01466 — a 4× stiffness lever moves the rate 1.9%.
**γ_res_m is NOT the origin** (the §7.10 candidate refuted). The
data hands over the real structure: the bare comb-detune rate
fix_dwb2pi = 0.0326 is IDENTICAL across all three (0.1%), so
ν/dw_bare = 0.443 — the RSJ depinning factor (LIT.md transfer #1,
measured at 0.44 before) is the γ-free invariant: **the fifth
circulator's rate = (its own bare detune, set by loads+geometry
through the pitch law) × the universal RSJ working point 0.443.**
Corollary test (future): reshape the load pattern → ν must track
the computed bare detune. And the octave check lands the §7.10d
reinterpretation independently: oct_tri has ν/dw_bare = −3.7
(wrong magnitude AND sign) — its circulation does not obey the
winding law at all, as expected of skirt-driven flow. The
"universal ±0.0145" across harmonics was two mechanisms
coinciding, and the ratio ν/dw_bare is the diagnostic that
separates them.

### 7.12 Q13 (2026-08-01): three-site assembly of one quantum —
### integer deposits PIN regardless of texture; the pump generalizes

The third kick window (pin_kick3_*, gated 21/21 byte-identical off)
enables the three-vertex coordinated pump. Loop algebra from the
tri2 net (E rows: 3:2 / 2:3 / 1:1): Λ = 2(δ_A − δ_C) — the heavy
vertex ALWAYS telescopes out; equal-|dv| three-site drive (A −0.3,
C +0.3, B +0.3) deposits across all three edges. MEASURED: the
transfer factor is 25% of the naive Λ design (constant in τ —
pin-ramp loss + the light vertices shunting drive into the
vacuum-unison tongue as x_A → skirt; an apparatus constant, so one
calibration run fixes it). Calibrated τ=23: **loop deposit
+0.998·2π — one quantum to 0.2% — and the endpoint is PINNED**
(all edges re-locked, zero flips; final texture strongly
non-uniform e0 +23.6 / e1 −3.7 / e2 −13.7 rad, loop exactly +1 —
arithmetic exact, texture free). VERDICT: with pump1 (+0.883,
two-site), ignite (+1.178, stranded) and now tri3 (+0.998,
three-site symmetric), EVERY quasi-static deposition — any
texture, any site count, integer or fractional — lands in the
pinned-strain state. Integer assembly from distributed thirds does
NOT ignite the running branch: the basin separation is
deposit-texture-independent and genuinely DYNAMICAL. Fractional
deposits at three sites (+0.263, +0.027 arms) are bound, as
before. Files: charge/runs/q13_*.{cfg,log}.

### 7.13 Q12 launcher (2026-08-01): the ℤ₃ wall is MOBILE — runs on
### clean track, pins at defects, evaporates at open ends

The kink launcher (q12_gen.py): a 14-voice fifth chain (alternating
0.195/0.709, ω ratio 1.5000 exact) with a +2π/3 ℤ₃ branch wall
seeded mid-chain (real topology on 3:2 edges; any 2π kink on a
unison chain is gauge). TWO TRACKS, one lesson each:
(1) v1 (unfiltered walk, edges to 1.91): the wall CHURNS at its
seeded site for 600 t.u. (+0.017/t.u. extra drift, same flip class
as background) — pinned at a defective edge; the track also had an
intrinsic defect edge (1.91 length, ~1000 flips in kink AND ctrl).
(2) v2 (track bounded to [1.30, 1.70], span 14.3): the kinked
chain's endpoint is IDENTICAL to control at the 1e-5 level, and
the wall's exit is written in the books: e12 (the open end) shows
kink−ctrl = **−2.094 rad = −2π/3 EXACTLY**, first crossing at
t≈80 — the wall transited ~6 edges and evaporated at the boundary
(open paths have no topological protection; e11's difference is
contaminated by intrinsic end chatter present in both runs).
VERDICT: **mobile dressed textures exist** — walls run on clean
track (transit ≥0.12 length/t.u. lower bound), pin at defects
(LAMH consistent, matching the cycle-side pinned strain), and die
at open boundaries. Combined with §7.12: deposits pin, launched
textures move — the charge sector's static/dynamic dichotomy is
now measured from both sides. Q12-proper needs: a longer
instrumented track (wall velocity + dressing profile), then the
two-path geometry; boundary management (the wall must reach the
slit, not an end). Files: charge/q12_gen.py, charge/nets/q12_*,
charge/runs/q12_*.{cfg,log}.

---


## 8. The mass–size–charge law (symbolic + numerical)

**User direction:** charge formation conditions quantize the size, and the
size quantizes the mass — "charges will only be able to form in mass under
specific size conditions"; electrons are unbound field harmonics. This
section makes that chain explicit (symbolic in `charge_sym.mac`, numerical
in `charge_fk.c`, `charge_fk2.c`).

> **CORRECTION NOTICE — 2026-08-02.** The first version of this section
> (2026-08-01) is retracted in part. Independent re-derivation found four
> defects in its artifacts; all four are repaired in the current files and
> each repair is documented at the head of the file it belongs to (D1–D11).
> Two of them changed conclusions, not just digits:
>
> * **§8.4's "the quark pattern, measured" was an artifact and is
>   WITHDRAWN** — see §8.4 below. Nothing in this document now claims a
>   measured fractional-charge energy.
> * **§8.1's agreement with the Bogomolny bound was measured on the
>   Peierls–Nabarro saddle, not the minimum.** The corrected ratio is
>   0.972, not 0.994.
>
> The monodromy of §7.1 was also being computed with a wrong composition
> law (the accumulated offset was not multiplied by `p_l`). The headline
> — charge = the ℤ_Q branch class — survives, but the comma and the branch
> coefficients printed by the old artifact were wrong on the *physical*
> fifth-triangle. Repaired in `charge_sym.mac` PART 2.

### 8.1 The kink energy, verified

Charge = k/Q lives as k phase-slip solitons (kinks) on the locked cycle
(§7.2). The continuum static kink energy is the Frenkel–Kontorova /
sine-Gordon Bogomolny bound for the **standing gate potential**
V(ψ)=1−((1+cosψ)/2)^p:

```
E_core = sqrt(2A) · I(p) · |k| ,   I(p) = ∫_0^{2π} sqrt(V(ψ,p)) dψ
```

| p | I(p) (quadrature) | measured E_core (A=1) | ratio to the bound |
|---|---|---|---|
| 1 | 4.000000 (exact, √V = \|sin(ψ/2)\|) | 5.615890 | 0.9928 |
| 2 | 4.591174 (exact, 2√2+2·asinh 1) | 6.436195 | 0.9913 |
| 4 | 5.057061 | 7.066302 | 0.9881 |
| **8** | **5.405336** (the standing p_gate) | **7.432761** | **0.9723** |
| 16 | 5.658564 | 7.657656 | 0.9569 |

**The discrete minimum sits BELOW the continuum bound**, by a lattice
correction that grows with p (the core narrows) and shrinks with A (the
core widens): at p=8 the ratio runs 0.8702 at A=0.25 → 0.9967 at A=4.
This is the expected direction and magnitude for a chain whose continuum
core width √(2A/V″(0)) = 0.71 site is *sub-lattice* at the standing
parameters.

The 2026-08-01 table reported 0.994 at p=8 because its relaxation started
from a single uniform-ramp seed, which on an even-length chain is
symmetric about a *bond* and therefore converges to the bond-centred
stationary state — the PN saddle (E = 7.598210, max|grad| ~ 1e-14: a
genuine stationary point, just not the ground state). `charge_fk.c` now
relaxes from several seeds and keeps the lowest (D4).

**Max gate torque is now a closed form, not a scan** (`charge_sym.mac`
PART 1): T′=0 ⇒ **cos ψ\* = (p−1)/p** for general p (the residual
simplifies to 0 symbolically), giving

```
T_max(p) = (sqrt(2p−1)/2) · ((2p−1)/(2p))^(p−1)
```

At the standing p=8: ψ\* = acos(7/8) = 0.50536051, T_max = **1.232578445**.
So §7.2's critical-detune constant 1.233·ε and the Q9 knee it predicts are
parameter-free *and* exact.

### 8.2 The mass–size–charge relation

A charge-carrying cycle of N voices, interval cover Q, trapped winding k:

```
M(N,Q,k) = N·e_balance(L)  +  [charge core]  +  (k/Q)²·E_dress
           \_____________/                       \___________/
            dense balance                         far dressing
            (MASS program)                        (EM5 cone)
```

The charge-core term, honestly:

* **integer k on a unison cycle: MEASURED.** E_core = |k|·7.432761 at
  p=8, A=1 — additive to all printed digits at k=1,2,3, because the cores
  are contact-range (κ = √(p/2A) = 2/site) and do not interact.
* **fractional k/Q: NOT MEASURED.** §7.2's 1/Q reduction of the core cost
  is a structural conjecture. It cannot be tested on a unison ring, where
  a 2πk/Q twist is not a soliton at all (§8.4). Testing it requires a
  **ℤ_Q interval chain** whose gate variable is the comb phase
  Φ = q·θ_i − p·θ_j — the same instrument Q11 v2 built in-kernel, applied
  to a straight chain rather than a cycle.
* **E_dress: not measured by anything.** Gated on the EM5 cone.

So one of the three terms is closed, and this section no longer claims
otherwise.

### 8.3 The size condition (the user's thesis, measured — `charge_fk.c`)

A charge-|k| state on a periodic cycle of N voices is a clean soliton only
when each core has room. Measured (p=8, A=1), and the result **collapses
onto N/k**:

| k | N | N/k | E_total | E/(k·E₁) | max step |
|---|---|---|---|---|---|
| 1 | 10 | 10.0 | 7.432778 | 1.0000 | 1.4006 |
| 1 | 6 | 6.0 | 7.452555 | 1.0027 | 1.4102 |
| 1 | 4 | 4.0 | 7.926439 | 1.0664 | 1.5885 |
| 1 | 3 | 3.0 | 8.579706 | 1.1543 | 2.0945 ≈ 2π/3 |
| 2 | 12 | 6.0 | 14.905110 | 1.0027 | 1.4102 |
| 2 | 8 | 4.0 | 15.852878 | 1.0664 | 1.5885 |
| 2 | 6 | 3.0 | 17.159411 | 1.1543 | 2.0945 |
| 3 | 15 | 5.0 | 22.613786 | 1.0141 | 1.4468 |
| 3 | 9 | 3.0 | 25.739117 | 1.1543 | 2.0945 |

The k=1,2,3 rows at equal N/k agree **to every printed digit** — the cores
are independent and each needs its own room. Measured floor: N/k ≥ 6 costs
<0.3%, N/k = 4 costs 6.6%, N/k = 3 costs 15.4% (and there the "core" is a
single hard 2π/3 jump — the whole cycle is uniform strain).

**Reading.** The winding is forced by the cycle BC and stays integer
regardless of N — so "charge cannot form" is not "winding→0"; it is "the
smooth core cannot fit, and the charge becomes a high-energy squeezed
defect." **Only cycles with N ≳ 6|k| host a clean charged state.**

**A caution about ξ.** Three different numbers have been called the core
width and they are not interchangeable: the continuum √(2A/V″(0)) = **0.71
site**; the FWHM of the step profile = **4 sites**; the capacity crossover
of `charge_pack.c` = **4 sites per charge**; the energy floor above =
**6 sites per charge**. The kink is deep in the discrete regime, so the
continuum formula is simply not the size scale here. Do not quote a single
ξ, and in particular do not carry a unison-chain ξ onto an interval cycle
of 3–8 voices — that was the unjustified step in the first species table.

The sharp criterion is not the energy excess but the **Peierls–Nabarro
barrier** — whether the core can move (`charge_fk2.c` N2, §8.4b).

### 8.4 Fractional charge — the 2026-08-01 measurement is WITHDRAWN

The first version reported a scan of E_core against w = k/Q on a unison
FK ring and concluded that integer charges scale linearly while fractional
charges are sub-linear (⅓ at 0.60 of the linear share, ⅔ at 0.90),
reading this as "the quark pattern, measured."

**That was an artifact of the apparatus, and the conclusion is withdrawn.**
2π/Q is not a period of V, so a twist of 2πk/Q imposed on a *unison* ring
is not a topological state and there is no soliton to measure. The relaxed
configuration is an elastic strain sitting on the ring's artificial
closing bond. `charge_fk.c` section 6 now prints the diagnostic as an
explicit negative control:

| w = k/Q | E_relax | load on the seam bond | max interior step | non-vacuum sites (of 400) |
|---|---|---|---|---|
| 1 | 7.432761 | 0.0000 | 1.4006 | 5 |
| **1/3** | **1.510164** | **1.3948** | 0.2895 | **2** |
| **2/3** | **4.541719** | **1.2362** | 1.2362 | 4 |

At w=1 the twist is carried by a localised core and the seam is empty. At
w=1/3 the entire twist sits on the seam and 398 of 400 sites are at the
vacuum — there is no core anywhere. The energies scale roughly as w²
(strain), not as the "0.60 / 0.90 of linear" that was reported.

What survives: §7.1's **confinement argument is untouched** (it is an
algebraic statement about which holonomies a unison medium can transport,
and it does not depend on this scan) — and it in fact *predicts* this
null. The pump experiments of §7.10c remain the real evidence that
fractional deposits are bound: those ran in-kernel on genuine interval
cycles, not on a unison FK ring.

What is now owed: the fractional core energy, measured on a ℤ_Q interval
chain. Until then no number for the cost of a third of a charge exists.

### 8.4b Mobility: the Peierls–Nabarro barrier (`charge_fk2.c` N2)

The PN barrier is the gap between the two symmetric stationary states of
the translation path: the **minimum** (symmetric about a site, one site
parked at π) and the **adjacent saddle** (symmetric about a bond, obtained
by shifting the relaxed minimum half a lattice spacing). At p=8, A=1:

```
E_site = 7.432761 (the minimum)   E_bond = 7.598210 (the saddle)
E_barrier = 1.654491e-01 = 2.2259e-02 of E_min = 2.23% of E_core
```

both states stationary to max|grad| ~ 2e-14.

**Two things the corrected instrument shows that the old one could not.**

*(a) Which symmetric state is the minimum flips with (p, A).* At p=8, A=1
the site-symmetric state wins; at p=8, A=0.5 and at every p=1 row the
bond-symmetric one does. That is structure of the PN landscape, and it is
why the barrier has to be taken as a |gap| rather than by assuming one
class always wins. It also gives a free cross-check: the winning class
agrees with `charge_fk.c`'s independently relaxed E_core at every (p,A)
the two artifacts share. The old apparatus, pinned at π by construction,
could not see the flip at all.

*(b) The mobility size floor is sharp, and it is not where the energy
criterion is.* Scanning N at p=8, A=1:

| N | 400 | 200 | 40 | 20 | 12 | 10 | 8 | 6 | 4 | 3 |
|---|---|---|---|---|---|---|---|---|---|---|
| barrier/E_min | 2.2259e-02 | 2.2259e-02 | 2.2259e-02 | 2.2259e-02 | 2.2260e-02 | 2.2264e-02 | 2.2405e-02 | **unresolved** | 6.3825e-03 | **unresolved** |

The barrier holds its large-N value to 4 digits all the way down to N=8,
then **collapses below the FP floor at N=6** — the two symmetric states
have become the same state, i.e. there is no longer a distinct kink
position to pin. At N=4 a small barrier reappears with the minimum
flipped to the bond class; at N=3 it is gone again (the cycle is uniform
strain, max step 2π/3).

So in mobility form the floor is **N ≥ 8**, against **N/k ≥ 6** from the
energy excess of §8.3. The two criteria answer different questions — "is
the core still a pinnable object" versus "does it still cost the additive
core energy" — and they are reported separately rather than averaged into
one ξ. (The 2026-08-01 version reported this scan as "constant 1.27% down
to N=8, drops to 0.59% at N=4, jumps to 5.5% at N=3". The level was wrong
by 1.75× throughout, and the N=3 "jump" does not exist: there is no
resolvable barrier there at all.)

The 2026-08-01 figure (0.094482, 1.27%) was low by a factor 1.75: it swept
a pinned site's phase over π ± 0.5 rad, but π *is* the minimum and the
saddle sits at ≈0.79 rad away — the sweep never left the bottom of the
well (D6). The correction matters beyond the digit: the earlier "the naive
PN threshold is 2× below the measured depinning, the gap is core
distortion" needs no such excuse once the barrier is right.

Note the mobility criterion and the §8.3 energy criterion are different
questions and do not have to agree on a single N.

### 8.5 The electron (§7.3 restated — ontology, not a measurement)

The electron is the **dressed delocalized kink with no dense core**: no
closed cycle, hence no N·e_balance term, so

```
M_electron = E_dress(field gradient)   — the field-gradient scale, naturally light
```

Its substance is the slip plus the far phase dressing ∮∇χ·dl = 2πk/Q; it
localizes only at conversion events, and the measured double-slit sector
(claim rule R3, MCWF-exact click statistics, Tonomura fringe statistics)
is electron phenomenology under this ontology. The mass hierarchy
(nucleon ≫ electron) is immediate.

**This is a reading of the ontology, not a measured result.** No artifact
in this document measures E_dress, or any far field at all.

### 8.6 Path to quark / proton / neutron masses

The charge piece of the mass budget is **not** closed. What is closed:

* the integer core quantum √(2A)·I(p) (measured, 2.8% lattice correction
  at the standing parameters);
* the holonomy arithmetic — which cycles exist and what quantum each
  carries (`charge/species_enum.py`; see the caveat below);
* the mobility criterion (§8.4b).

The remaining inputs, unchanged in kind but now correctly scoped:

1. **e_balance(L)** — the dense balance energy per voice. The MASS
   program's quantity; blocked on a stable balance point (livefab).
2. **The fractional core cost** — the 1/Q scaling, on a ℤ_Q interval
   chain. Newly owed as of the §8.4 withdrawal.
3. **The allowed (N,Q,k)** — closure arithmetic intersected with a size
   floor. The closure half is enumerated; the size half needs (2) first,
   because the floor for a fractional core is not the floor for an
   integer one.
4. **E_dress** — gated on the EM5 cone.

**Species table caveat.** The 2026-08-01 enumeration reported 175 classes
and a proton/neutron slot count. It enumerated sorted edge *multisets*
only, so it saw one cyclic ordering of each — but Q depends on the
ordering (its own docstring says so, and the example it cites was among
the cases it silently dropped). The corrected ordered enumeration gives
**17922 classes** over N ≤ 8; at N=4 alone the census goes from
{Q=1:4, 2:2, 3:5, 4:1} to {Q=1:15, 2:11, 3:9, 4:1, **6:1**}. The
structural conclusion survives intact and is the one thing that file now
claims: **every charge quantum is 1/(2^a·3^b)** — the comb carries only
the primes 2 and 3, so no fifth-of-charge exists at comb_limit=6, and the
small quantum numbers are the small harmonic numbers. No mass table and
no slot assignment is emitted any more (D9).

### 8.7 Artifacts

| file | role | status |
|---|---|---|
| `charge_sym.mac` | I(p) with exact I(1), I(2); corrected N-cycle monodromy → charge = k/Q; exact torque locus cos ψ\*=(p−1)/p; V″(0)=p/2; tail rate κ=√(p/2A) | verified 2026-08-02 |
| `charge_fk.c` | FK charge soliton: E_core vs p and A (multi-seed, minimum kept); size condition vs N/k; integer additivity; the fractional-twist negative control | verified 2026-08-02 |
| `charge_fk2.c` | two-kink force law (sign = sign q₁q₂, contact range); PN barrier via the symmetric stationary states | verified 2026-08-02 |
| `charge_pack.c` | packing capacity of one cycle; merged-jump vs relaxed-lattice cost | verified 2026-08-02 |
| `charge/species_enum.py` | closure + holonomy census (ordered cycles) | verified 2026-08-02 |

Build: `gcc -O2 -o charge_fk charge_fk.c -lm` (same for `charge_fk2`,
`charge_pack`). Run: `maxima -b charge_sym.mac`, `./charge_fk`,
`./charge_fk2`, `./charge_pack`, `python3 charge/species_enum.py`.
Logs in `charge/runs/charge_{fk,fk2,pack}.log`, `charge/runs/species_enum.log`.

**Standing discipline for this section.** Every one of the defects above
was invisible in prose and visible in the run log — a seed that converged
to the wrong branch, a sweep that never reached the saddle, an enumerator
that skipped 99% of its domain. Artifacts here now print their convergence
witness (max|grad|) and their negative controls. A number in this document
without a log line behind it is not a result.
