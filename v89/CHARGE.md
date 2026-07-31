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
