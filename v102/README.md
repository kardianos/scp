# v102 — SPIN-ONE FIBER

**Opened:** 2026-10-07.
**Authority:** user, this session. After `v99/reopen/OPTIONS.md`, the
user accepted three of the five A+B premises as rulings ("I'm okay
with those changes. Go a bit further."). They are R1, R2, and R5
below. The other two are gates. **This document is the plan.** Gates:
`PREREG.md`. A pass or a stop is claimed only in `DECIDE.md`.

**No kernel is carried.** `v97` is not edited. Opening v102 is not a
law-table change and does not authorize a mesh. The rulings change
what the theory may assume. They do not change the code.

**Status:** the gates were run on 2026-10-07, ending C.

- G-wedge **passes**: the antiferromagnetic, bilinear link exchange is
  derived from R1 and R2.
- G-two **fails**: a cell's occupation equals the lock's occupancy
  ceiling, whose value nothing fixes.

One sentence would fix it: "the occupancy ceiling of a lock is one
quantum per link end." Phase 1 is not opened. See `DECIDE.md`.
Reproduce:
`python3 -I v102/gates/cas/wedge.py`, then
`python3 -I v102/gates/cas/two.py`.

---

## Thesis under test

A lock is the gapped quantum of a fiber chain. Each cell of the chain
carries quanta of a two-plane (Jones-vector) field. Each link converts
field into space through the link's handedness. If that is the whole
story, the fiber is a spin-1 antiferromagnetic ring (family 268 in
`github.com/openai/math`), and its correlation length (ξ ≈ 6.03
grains) is forced, with no continuous coefficient.

Two things stand between the rulings and that claim:

1. **G-wedge.** Is the link exchange antiferromagnetic and bilinear
   because of the premises, or by choice?
2. **G-two.** Is there one spin per cell, `N = 2` quanta, because of
   the premises, or by choice?

---

## Rulings (user, 2026-10-07)

| | ruling | what it reverses or retires |
|---|---|---|
| R1 | Each cell's field is two complex amplitudes, one per transverse plane (a Jones vector), with a U(2)-invariant coupling (the wording accepted in `v99/reopen/FEASIBILITY.md` §7 A). Every coupling of the field is SU(2)-invariant, and transport is U(2)-invariant. A conversion changes the number of field quanta, so it cannot also keep the U(1). | The kernel's single complex amplitude per cell. |
| R2 | Transport is quantized. The quanta are those of the classical wave, so they are bosons. | v89 law 1, "transport is never quantized." |
| R5 | Integer `W` is retired as the lock. The lock candidate is the gapped quantum of the fiber chain. | v98 §4's "winding is an integer" as the lock, and the v99 wrap thesis. |

Not ruled at opening: the occupation per cell and the form and sign
of the link exchange. Those were the gates.

| | ruling, after the gates | closes |
|---|---|---|
| R6 | **The occupancy ceiling of a lock is one quantum per link end** (user, 2026-10-07: "That's okay, can you continue with that modification?"). | `DECIDE.md` §2's open sentence. With it, W4's ceiling exists, `N = 2` on the ring, and an end cell of an open chain (degree 1) holds one quantum. |
| R7 | **S1–S4 of `law/DECIDE.md`** (user, 2026-10-07: "Sure, follow those rabbit trail S1-S4 and try to make a lock work."). A link's space mode is a discrete oscillator at `Ω = w2`. Conversion is a coherent pair amplitude `g`. Transport inside a lock is zero. | `law/DECIDE.md` ending C. The lock is tested under this "S-law" in `lock/` (`lock/PREREG.md`). The rulings change the theory, not the code: no kernel implements S1–S3, so no battery has run. |
| R8 | **The 3+4 widening (W)** (user, 2026-10-08: "Allow the law change for 1."). Each field quantum carries spin ½ × isospin ½, four modes per cell. Conversion pairs a quantum across a link into a singlet in both (P = ε ⊗ ε). Charge is read as `I^z + ½` on an open slot. | `perm/DECIDE.md` §7.1. R1's single Jones doublet becomes two doublets. The theory changes, not the code: no kernel step is authorized. Decisions 2–3 (transport and space at the lock) are settled by test in `attract/`. |
| R9 | **T-floor** (user, 2026-10-08: "accept decisions 2+3"). Transport inside a lock is allowed with amplitude t < t_max, where t_max lies between 0.1 and 1 J0 (`attract/DECIDE.md` §2). Transport may not empty a lock cell. **Working value t = 0.1 J0** (9.6e-4 in rate units), the largest t at which every gate passed. | `options/` T0 (t = 0 by label). The hole walks; hole pairs cannot collapse into fiber cuts. |
| R10 | **S-boundary** (user, same). Space may leave a lock only through a link with one end outside the lock, at the table's `s_k`. | A limited exception to reading V. It is the formation route: light makes an open slot. |
| D11 | **The working regime is small t; the ladder is deferred** (delegated to me: "Make the best decision yourself here"). At t = 0.1 J0 the nucleon is firm, walks, forms, and binds into a nucleus. The only ladder found needs t ≈ 7 J0, where pairs unbind and the nucleon deforms. v98 puts unique size first. A ladder at small t stays an open search. | `attract/DECIDE.md` §7.2. |
| K | **Kernel step authorized** (user, same: "(3) approved"): put R8 + R9 + R10 into a real simulation. Plan: `kernel/PLAN.md`. | v97's "do not edit the kernel" no longer applies to this step. The ratchet still applies: full battery, new behaviour off by default. |
| R12 | **Singlet core** (user, 2026-10-09: "Go with singlet core, keep going"). Two N = 1 cells that share a link (two adjacent holes, or a hole beside an end cell) may not hold their pair in the conversion singlet: that component is excluded. Their other channels carry no term. It is reading V taken strictly: both exact branches of a singlet (1,1) pair hold about half a space quantum on the link (`bind/PREREG.md` §2). No new constant. | The leading-order lock's "no term" on (1,1) links (K2b §5 wall 1). In the kernel: `lockq_p11=2`, with the excluded component projected out of each jump (`kernel/K2e_PREREG.md`). |

Phase 1 is opened under R1, R2, R5, and R6: `PHASE1.md`, with gates
in `phase1/PREREG.md`.

**Phase 1 stopped** (2026-10-07; `phase1/DECIDE.md`).

- **The lock is a phase.** Transport must stay below roughly
  t/Δ ≈ 0.16–0.26; above that, the field converts into a ferromagnet.
- **Inside the phase, the fiber is a Haldane chain** with one line and
  no bound pair.
- **ξ ≈ 6 grains only at leading order.** It drifts continuously with
  g/Δ and t/Δ: 6.06, 5.79, and 5.07 at g/Δ = 0.05, 0.1, and 0.2.
- **Two registered bars were defective** (P1 at n = 4, and P6(c)).
- **The two physics misses are real:** the transport threshold, and
  ξ's drift with g/Δ.

**Law step** (2026-10-07; `law/DECIDE.md`): **ending C.**

- The current table does not host the lock.
  - Space has no pitch, so a field pair resonates with it and decays.
  - g has no counterpart: the table's conversion is a rate law.
  - Even borrowing the retired dense rung, 3-D transport is about 10×
    too fast.
- The lock needs four new law statements (S1–S4). Each is the user's
  ruling, and a law-table change.
- The two defective Phase 1 bars were re-registered on new cases, and
  both passed.

**Lock step under R7** (2026-10-07; `lock/DECIDE.md`): **ending A,
in a model.** Under S1–S4 (the "S-law"), a fiber ring at its ceiling
is a working lock:

- **It holds its field.** Over 100/J, from every unlocked start, space
  stays at or below 0.103 per link. With space refused a pitch, the same
  ring converts at the golden-rule rate.
- **Its quantum walks only by conversion.** It moves at rate
  J ∝ g², and the field never travels.
- **Its quantum is massive and single.** Δ_∞ = 0.436 J, there is one
  line, and space-carrying states sit 32× higher.

Costs: transport must stay under t/Δ ≈ 0.16–0.18, and the space
band's half-width under Δ. The lock cannot cool itself. `g` is
unset. Size sits on the smallest even cycle. No kernel or battery
carries any of it.

**Options step** (2026-10-08; `options/DECIDE.md`): **a measurement,
no ruling.** Each open question was given 3–5 options. Each option was
tested alone and in combination, and the bundles were scored by a rule
fixed in advance.

- **The lock needs one label: which cells are in it.**
  - The geometric ceiling leaks the lock, or makes spin z/2 cells.
  - The dynamic ceiling freezes adjacent free light.
- **Given the label, t = 0 on lock links is the robust transport
  rule.** The local "frozen at ceiling" rule is identical alone, but
  collapses under any space flow on lock links.
- **A lock forms only through boundary conversion, and only from light
  colder than its gap:** about 10⁸ photons at g/Δ = 0.1. Resonant light
  heats it or destroys it.
- **The lock is metastable.** Per quantum it sits far above free light.
- **Touching locks merge.**
- **Space may flow on a lock only coherently,** with 2t_s < 0.2275.
  Rate-law flow creates energy.
- **The door spares a lock that carries no space deficit.** Removing
  the door costs 10 bars.
- **The g window is g/Δ ≤ 0.1,** with the table's `s_k` reading on its
  edge.
- **No size is selected.**

Best bundle: **B-label-cool** (T0, C-lock, X1, no space flow on lock
links, door exempt). It needs one label.

The user then noted that "cold" is a layer mismatch: the fiber and its
links are what energy is made of, and they cannot be hot. Read that
way, the cooling requirement was about emptying the fiber of its lock
quanta, not about making a lock. `options/DECIDE.md` §12 sets out
which results that touches.

**Fabric step** (2026-10-08; `fabric/DECIDE.md`): **ending B.** The
user's two ideas were pursued toward "a firm proton or neutron from the
underlying fabric, not cellular clumping".

- **Both ideas work.**
  - Reading V (the link holds no energy) reproduces the S-law's low
    spectrum to 0.03 J.
  - Under Q (count lock quanta), one resonant photon makes exactly one
    lock quantum 12–28% of the time.
- **The fabric hosts one firm spin-½ object: a hole (one open slot).**
  Its size is set by the fabric (decay length 5.8 grains, the same at
  two lengths), and it has a gap.
- **It falls down in four places:**
  1. **There is no neutron.** At its ceiling a fiber holds an even
     number of quanta, so spin ½ needs an open slot. If charge is the
     deficit, no neutral spin-½ object exists.
  2. **It is pinned at t = 0.** It walks only if transport is allowed
     inside the lock (0 < t ≲ 9 J).
  3. **It forms only if space leaves the boundary link,** which V
     forbids.
  4. **It has no excited ladder.**

**Permutation step** (2026-10-08; `perm/DECIDE.md`): **a measurement, no
ruling.** The re-reading (charge as isospin, I, or v98's 3+4 widening,
W) and the decisions (hole transport, T; space leaving the boundary
link, S) were tested in all 12 combinations. The long chains used
TeNPy, a dependency the user approved.

- **One cell reaches the Full tier: W + T + S.**
  - The widened fiber is a lock (gap 1.32 J0, three times 3+2's).
  - Its open slot is a firm (½, ½) quartet: p and n, both spin ½, one
    object of one size (2.3 grains, stable 24 → 32).
  - It walks under T, and light makes it under S (≈ 5% per resonant
    photon).
  - The moving W hole's size is carried from the static one.
- **I alone trades spin for p/n.** Every I cell fails.
- **No cell has an excited ladder,** and no cell has a nucleus: two
  holes convert, or under T collapse into a fiber cut.
- **T at t ≳ J0** keeps a length-independent cloud that fails the
  registered size fit; at 3 J0 it is a Nagaoka precursor.
- **T + S together** fills a lit fiber with holes, and pairs of them cut
  it.

**Attraction step** (2026-10-08; `attract/DECIDE.md`): R8 (the 3+4
widening) is ruled by the user. Decisions 2–3 were settled by the
registered tests:

- **T-floor:** transport inside a lock is allowed with t < t_max,
  where t_max lies between 0.1 and 1 J0, and it may not empty a lock
  cell. Without the floor, holes collapse into fiber cuts.
- **S-boundary.**

Item 4 reaches **ending LN**, with **the fiber's own correlations as
the attraction**:

- **A nucleus.** Two open slots bind into a pair of fixed size
  (about 3 grains) in the all-singlet channel: W (S, I) = (0, 0) with
  B = −0.50 J0, 3+2 S = 0 with B = −0.13 J0, at t = 0.1 J0. Every other
  channel is free.
- **A ladder.** A bound (3/2, 3/2) or S = 3/2 level appears only at
  t ≈ 7 J0 (spin polaron). There the pairs unbind.
- The two never coexist. The fourth-order terms at larger g give no
  ladder.

**Kernel step, K1** (2026-10-08; `kernel/PLAN.md`, `kernel/K1_DECIDE.md`):
**Go.**

- **Formation is one quantum jump.** It matches the exact absorption
  with fidelity 1.0000 in both fabrics, because the space band is far
  wider than the lock's internal spread. The kernel can implement it as
  an event.
- **A closed lock never binds two holes:** the binding energy has
  nowhere to go.
- **Cold light captures them** (radiative capture): light colder than
  about 0.2 J0 brings the pair to 100%, but at about 10¹⁰ photons per
  pair.
- **Next: K2a,** a C module `lockq` in `v97/kernel`, off by default.

**Kernel step, K2a** (2026-10-08; `kernel/K2a_PREREG.md`,
`kernel/K2a_DECIDE.md`): **Go (K2b).**

- **Off, it is byte-inert:** the battery on both kernels is ALL GREEN
  (87 bars) and identical, and kernel stdout is byte-identical.
- **On, it matches the Python lock model** (`lo.py`) to every printed
  digit, for the ground energies and for 20 000 steps of hole and spin
  dynamics.
- **Light makes holes**, one atom per event, at the golden-rule rate
  onto the space band. The ledger closes at the kernel's own floor
  (drift ≤ 3.5e-15).
- **One gate failed as first written:** a hole formed beside another
  in W shifts the space quantum out of the band. The window now uses
  the exact energy change of each jump.
- **The kernel's light is one mode,** so every hole it makes has the
  same charges, and the all-singlet nucleus cannot form. That is K2b's
  job: four field modes per cell.

**Kernel step, K2b** (2026-10-08; `kernel/K2b_PREREG.md`,
`kernel/K2b_DECIDE.md`): **Go, and stop before K3.**

- **The kernel field now carries 1, 2 or 4 complex modes** (`fmodes`),
  mode-blind.
  - Off, it is byte-inert: battery 87 ALL GREEN and identical.
  - Four modes with an equal split are byte-identical to one mode,
    step for step, and permuting the modes changes nothing.
- **Light in mode ν makes holes carrying ν's charges,** checked in the
  kernel on every event, with the ledger at 3.6e-15. Under mixed light
  the lock reaches the nucleus's zero-charge two-hole sector.
- **Two bugs were found and fixed,** one a latent overflow in K2a's
  code. K2a's recorded results stand.
- **Two walls before K3:**
  - bright light makes adjacent holes, which the rulings leave unpriced
    (perm §5: they bind by converting, or merge);
  - pairing needs a sink, and the only one known (cold-light capture)
    takes about 10¹⁰ photons per pair.

**Bind step** (2026-10-09; `bind/PREREG.md`, `bind/DECIDE.md`):
**ending S.** A measurement in models; it rules nothing.

- **Energetically, the lock forces a size: pairs (N\* = 2).** Under
  zero, core and upper, in both fabrics and at two lengths, the second
  hole binds (−0.13 J0 in 3+2, −0.50 J0 in W, ⟨d⟩ ≈ 3), the third is
  free and the fourth makes a second pair. Score (now R12) gives the
  same, its second length run as a labelled supplement.
- **The binding parks in the formation space quantum (P1).** A
  complementary hole formed at d = 3 is born into the pair with weight
  0.66 (3+2) or 0.52 (W), inside the window. P2–P4 are closed or cost
  about 10¹⁰ photons.
- No field-reading ground state holds an adjacent pair, so the (1,1)
  pricing matters only for what formation can make.
- Lower (the space branch) clumps at d = 1 and violates V.

**Kernel step, K2d** (2026-10-09; `kernel/K2d_PREREG.md`,
`kernel/K2d_DECIDE.md`): **Go.**

- The five (1,1) pricings of `bind/` (zero, core, score, upper, lower)
  are kernel options (`lockq_p11`). They are byte-inert at the default
  (battery 87 ALL GREEN) and match `lo.py` to 1e-11.
- In bright light, zero, upper and lower fill the lock with adjacent
  holes, far from bound.
- Core and score (with K2a's mean-energy window) keep holes apart. The
  3+2 lock jams at three holes, and a complementary pair at d = 3 is
  born bound (weight 0.51 under score).

**Kernel step, K2e** (2026-10-09, under R12; `kernel/K2e_PREREG.md`,
`kernel/K2e_DECIDE.md`): **Go.**

- `lockq_proj` projects each jump onto R12's allowed subspace, so only
  the singlet part of an adjacent jump is refused. It is byte-inert at
  the default (battery 87 ALL GREEN) and matches an independent
  projector at every printed digit.
- **Under R12 exact, light builds blocks.** Holes form beside holes in
  the allowed channels. Bright 3+2 ends with 6 of 8 sites holed, mean
  B = +15 to +17 J0; W ends at +4.5 to +12 J0. Nothing is bound.
- K2d's born-bound pairs came from the mean-energy window, not from
  R12.
- Formation sees adjacency only through the lock energy it costs, about
  3.6 J0 per N = 1 neighbour, inside a ±12.5 J0 window. It cannot see
  binding (0.13–0.50 J0).

**Kernel step, K3** (2026-10-09, under R12; `kernel/K3_PREREG.md`,
`kernel/K3_DECIDE.md`): **Go.** 26 registered runs, T = 2000; every bar
passes.

- **Light fills the lock with unbound blocks.**
  - Bright 3+2 always ends at seven holes in eight bulk sites, because
    the eighth is refused by R12 and by the window.
  - Bright W ends at three or four, by the order the holes formed.
  - Dim light is light-limited: W makes one hole, beside the end.
  - No run has any weight on bound levels.
- **A single hole walks in the 3+2 bulk** (D = 0.63) and is pinned
  beside an end. In W n = 7, which has no bulk, it is pinned.
- **A seeded bound 3+2 pair holds** for 2000 t.u. (bound weight 0.64,
  ⟨d⟩ = 3.0). W n = 7 has no bound two-hole level.
- **The size the energetics force (pairs) never forms from light.**
  Formation sees adjacency but not binding, and a closed lock cannot
  relax. The candidates (a reverse boundary event, core in place of
  R12) each need a ruling.
- **A seeding defect:** the lock's xorshift draw, started from seeds 1–5,
  makes each run's first event non-random. Part A was re-run with
  well-mixed seeds as a labelled supplement; the picture stands.

## What is carried

- `v89/PRINCIPLE.md`, in full. Where the kernel and the principle
  disagree, the principle governs. That matters here, because the
  kernel converts on a cell while §3 converts through a link.
- `v89/README.md` law 2, "a channel is a joint process of its two
  ends."
- v98 §1, "locked field cannot travel as itself," and the unfold
  filter.
- The v99 representation split.
- v100 §6: no hedgehog and no radial profile.
- v101 rules 4 and 5: no continuous ratio to break a tie, and no cell
  count read off a result.

## What is not carried

The v97 kernel as the law. The two-atom credit as an occupation. A
scanned coupling ratio. Integer winding as the lock.

---

## Order

| step | what | opens when |
|---|---|---|
| G-wedge | link exchange: operator, form, sign, net sign | `PREREG.md` written |
| G-two | occupation per cell | `gates/WEDGE.md` decided |
| Phase 1 | derived chain: line, ξ, no second line, slit regression | both gates pass; not authorized by this page |

## What not to do

1. Do not consult any version before v89.
2. Do not edit the v97 kernel or `v97/fab`.
3. Do not read ξ as a cell count, and do not take `n = 6` from it.
4. Do not scan `g/Δ`, `t/U`, or a biquadratic angle to place a line.
5. Do not use the two-atom credit (a credit that lapses at two atoms)
   as an occupation.
6. Do not map anything onto 938 MeV, 0.84 fm, or 1836.
7. Do not seed or score a hedgehog or a radial profile.
