# v99 reopen — options C, A, B, and A+B

**Written:** 2026-10-07, after `FEASIBILITY.md`.
**Follow-on:** A+B was opened as `v102/` the same day, under rulings
R1, R2, and R5.
**Authority:** user, this session: look at C briefly, then at A, B,
and A+B, each in isolation.
**Status:** analysis only. **Nothing is opened.** No kernel or law
change. Numbers: `cas/out/options.out`
(`python3 -I v99/reopen/cas/options.py`); `cas/out/gap.out` for X
and W.

The options, from `FEASIBILITY.md` §7:

- **C:** the size gate becomes "one isolated line, plus ratios that
  do not depend on the law's one free number."
- **A:** a link field with two complex amplitudes and a U(2)-invariant
  coupling; integer `W` retired as the lock.
- **B:** spin ½ on each link end, which means reversing v100 G2.

---

## 0. Verdict

| option | forced scale? | selects ring size `n`? | new premises | verdict |
|---|---|---|---|---|
| C | not a mechanism | — | gate change only | no recorded result passes it; it only matters for A |
| A | **no**: it is a ferromagnet (E1) | no | 2–3 | not promising alone |
| B | only with a free coefficient (E3) | no | 3–4 | not promising alone |
| A+B | **yes**: ξ ≈ 6.03 grains, Δ ≈ 0.41 J, no continuous coefficient | **no** (E2: `E₀` concave in `n`) | 4, two of which reverse the record | the only coherent one; expensive; its two decisive gates are of the kind v100 G2 already failed |

---

## 1. C, briefly

C re-reads results; it creates none. Against C:

| result | one isolated line? | a ratio independent of the free number? |
|---|---|---|
| v98 1B trough, `E = αE_w + βn` | one `n*` for each `ρ`, but `ρ = 1` ties {4, 6} | no. It works at the grain scale, so every observable moves with `ρ` |
| v99 `S_rel` | no line | — |
| v100 hedgehog | one minimizer | no (the rotor carries the Skyrme parameter); refused anyway |
| v101 succession | zoo | — |

Nothing on record passes C.

**What C requires in practice:** a separation of scales, `ξ ≫ grain`,
so that grain details drop out. Only a transmuting mechanism provides
that, which means A with a continuous coupling. Under C, A's
exponential sensitivity stops mattering: β only sets the unit. C has
no effect on A+B, because A+B has no continuous coupling to vary.

**C and the record:** C is a literal reading of v98 §I.4: "the only
lawful comparison is a pure number the theory is forced to produce."
Demanding `n*` in grain units is stricter than §I.4. The price of C
is that the theory never predicts a lock's size in grain units, so
the grain stays unobservable.

---

## 2. A in isolation

**Field.** A is arguably a more literal reading of `v89/PRINCIPLE.md`
§3 than the kernel. A link actuating in two transverse planes has a
Jones vector `z ∈ C²` (one complex amplitude per plane). The kernel
keeps a single complex amplitude, so it keeps only one circular
component. With a coupling that overlaps both components equally
(the natural lift of the kernel's hop, `z_i' = cos τ z_i − i sin τ z_j`
applied to each component), the field sector is U(2)-invariant.

**How each regime fails:**

| regime | result | why it fails |
|---|---|---|
| classical (transport not quantized, as v89 states) | the 2-D Dirichlet energy has no scale | v98 1B: no trough |
| thermal (a bath) | an O(4) model with `β = J/T` | size depends on the environment, which makes it a star (v98 §8) |
| quantized with bosonic quanta (atoms) | **ferromagnet** | **E1** |

**E1.** Two-mode boson rings with a spin-blind hop of either sign and
a spin-blind on-site load term (the analog of the law's load
nonlinearity). In every interacting case, the ground state has
**maximal total spin**. That includes the fabric's own hop sign (`+t`)
on odd rings. The free `+t` odd ring has a degenerate, frustrated
level with mixed spin, and any `U > 0` restores the ferromagnet. This
is the Eisenberg–Lieb theorem (2002) for spin-blind bosons. A
ferromagnet has gapless magnons whose energy grows as `k²`. It
generates no scale.

**Other consequences of A:**

- Two-mode bosons at fixed occupation `N` form one spin `S = N/2` per
  cell (Schwinger). The on-site singlet that an O(4) rotor needs does
  not exist. 215's O(4) result therefore does not apply to a
  quantized A. The reachable gapped model is the O(3) one (268),
  and only if the exchange is antiferromagnetic.
- Antiferromagnetism needs a link energy that rewards the **wedge**
  `z_i ∧ z_j = z_i¹z_j² − z_i²z_j¹`: the pair's handedness, which
  vanishes when the two ends are parallel. By the C² Lagrange
  identity, `|z_i∧z_j|² = |z_i|²|z_j|² − |z_i†z_j|²`. Rewarding the
  wedge penalizes overlap (antiferromagnetic). Rewarding overlap,
  which is what the hop does, is ferromagnetic. No term in the law
  rewards the wedge.
- π₁(S³) = 0, so integer `W` is gone. Handedness becomes a direction
  on the Poincaré sphere, not a sign.

**Verdict:** A alone yields a ferromagnet, a hot blob, or no scale.
Not promising.

---

## 3. B in isolation

The field stays U(1). B adds spin ½ as a separate non-abelian label.
Until it is derived, that is an imported species (v100 G2: no
derivation from the two-atom credit).

**Where the spins sit decides everything:**

- **One per link.** That is spin ½ per unit cell, so
  Lieb–Schultz–Mattis applies: the ring is gapless or degenerate. E2
  shows the S = ½ ring's gap closing as `1/n` (`n·gap` ≈ 4.0 → 4.34
  for n = 4 → 18). No scale.
- **One per link end, so two per cell.** A rule is then needed for
  how a cell's two halves combine. **E3:** on an alternating ring
  with `J_link = +1`, the gap is a continuous function of
  `J_cell/J_link`. It is 1.0 at `J_cell = 0` (isolated link singlets),
  0.78 → 0.20 as `J_cell` goes from −0.5 to −8, tends to 0.103 as
  `J_cell → −∞` (the spin-1 chain), and closes at `J_cell = +1` (the
  uniform S = ½ chain). This is `ρ` again (Hida 1992: the whole
  ferromagnetic side is one Haldane phase, continuously connected).

**Other consequences of B:**

- The spins do not interact with the U(1) wrap. `W` and the spin
  chain are separate sectors unless a coupling is added.
- Ring size is not selected (§4, E2).

**Verdict:** B alone gives a scale only through a free coefficient.
It does not select `n`, and it stays an import until G2 is reversed.
Not promising.

---

## 4. A+B

**The synergy.** Each option supplies what the other lacks.

- A's bosonic two-mode cell removes B's coefficient. Two quanta in
  two modes can only form `S = 1` (states |20⟩, |11⟩, |02⟩). The
  singlet does not exist, which is `J_cell = −∞` exactly, with
  nothing chosen.
- B's half-per-link-end becomes **one quantum of A's field**. It is a
  spinor of A's internal (Poincaré-sphere) SU(2), not a spatial
  spinor. Under a physical rotation about the link axis, one quantum
  carries helicity ±1. So A+B does not reverse G2. Spatial spin ½
  stays underived, and A+B does not need it.
- A single bilinear antiferromagnetic exchange across each link,
  projected onto `S = 1` cells, is `(J/4) S_i·S_j`: the pure spin-1
  Heisenberg ring, which is family 268's object.

**What it gives:**

- A gap Δ ≈ 0.4105 J that stays positive as `n` grows. E2 on even
  rings: gap 1.00, 0.72, 0.59, 0.52, 0.48 for n = 4 … 12, while
  `n·gap` rises (an S = ½ ring's would stay flat).
- ξ ≈ 6.03 grains (White–Huse 1993). Neither has a continuous
  coefficient; J is the unit.
- One isolated line: the triplet magnon below the two-magnon
  continuum.
- Odd rings are frustrated. In E2, odd rings sit +0.08 to +1.2 above
  the extensive line `n·e∞`, while even rings sit below it. That
  matches v90's measured F-odd ("odd rings π-frustrated"). The match
  is a consistency check, not a derivation.

**What it does not give: a ring size.**

- On even rings, `E₀(n)` is **concave** in `n` (second differences
  −0.102, −0.038, −0.018). `E₀(n) + c·n` has no interior minimum for
  any `c`. The ring sits on the smallest even cycle or hoards, the
  same endings as v99 Phase 1.
- 268 itself proves the gap is uniform in `n`.
- The forced scale is ξ, a property of the excitation. It is not `n`.
  ξ = 6.03 is a decay length and is not an integer. v101 rule 5
  forbids reading a cell count off it, and that applies here.

**The premises it costs (each is the user's decision):**

1. Two complex amplitudes per link with U(2)-invariant transport.
   This is a law and kernel change.
2. Quantized transport. This reverses the v89 synthesis "amplitudes
   within a mode; transport is never quantized."
3. Exactly two quanta per cell. Reading the two-atom credit as an
   occupation is the move `v100/phase0/SPIN.md` §4 refuses (a ledger
   count turned into a label). If `N` is free, the result is a ladder:
   odd `N` is gapless, `N = 2` gives ξ ≈ 6, `N = 4` gives ξ ≈ 50
   (literature), and so on. That is a zoo.
4. A bilinear, antiferromagnetic wedge exchange on each link. This is
   a new law term, and its form is a choice. A biquadratic admixture
   moves ξ continuously across the Haldane phase, from ∞ (θ = −π/4)
   through 6.03 (Heisenberg) and 0.91 (AKLT) to ∞ (θ = π/4).
5. Integer `W` is retired as the lock. The lock becomes the gapped
   triplet quantum of the fiber chain, against v98 §4 and the v99
   thesis.

---

## 5. If A+B were opened: the gates that would come first

Paper gates, pre-registered, in this order. Either failure ends it.

- **G-wedge.** Derive from `v89/PRINCIPLE.md` §3 that the link energy
  of two Jones vectors is `−J|z_i∧z_j|²` (antiferromagnetic,
  bilinear), with no biquadratic term and no choice of sign. "The
  relative phase is the handedness" is a reading. It is not that
  derivation.
- **G-two.** Derive two quanta per cell from the conversion law
  without equating the two-atom credit (atoms booked per completed
  cycle) with cell occupation.

Only then would exact diagonalization or DMRG on the derived chain be
worth running. The bars would be: one line; ξ measured and not read
as a cell count; no second line; the slit recovered from the
linearized field. The record suggests both gates are the kind v100 G2
failed. That is the real feasibility estimate for A+B.

---

## 6. What is now true, and what is not

**True:**

1. With spin-blind transport, a two-mode bosonic fiber is a
   ferromagnet and generates no scale (E1).
2. Spin ½ per link end, without A, leaves a continuous coefficient
   that sets ξ (E3).
3. A+B with a bilinear antiferromagnetic exchange is the spin-1
   Heisenberg ring: a forced ξ ≈ 6.03, Δ ≈ 0.41 J, and one line.
4. No option selects ring size. `E₀(n)` is concave on even rings (E2).
5. No option on record passes C.

**Not true:**

- `n = 6` or any cell count was derived. ξ is a decay length.
- A+B's premises follow from v89. Two of them reverse the record.
- Spatial spin ½ was derived. A+B uses an internal SU(2).
- Anything was opened or authorized.
