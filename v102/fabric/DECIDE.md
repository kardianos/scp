# v102 fabric — decision (a firm nucleon from the fabric?)

**Written:** 2026-10-08, after the outputs in `cas/out/`. The
criteria, the algebra and the endings are in `PREREG.md`.
**Status:** a measurement in models. It rules nothing. No kernel and
no law table was edited.

---

## 1. The answer: **ending B**

**Both of the user's ideas work.**

- **V (the link holds no energy)** reproduces the S-law's whole low
  spectrum to 0.01–0.03 J. It removes exactly the states behind the
  layer-mismatched failures, which are 90–95% of the microscopic
  states.
- **Q (count lock quanta)** turns resonant light from "heating" into a
  clean source. One photon leaves exactly one lock quantum 12–28% of
  the time.

**The fabric hosts one firm spin-½ object: the hole,** one open slot in
the fiber.

- It is a defect of the fabric, not a cluster built for it.
- Its spin-½ cloud has **one size**, a decay length of 5.8 grains that
  does not change with the fiber's length.
- It is the lowest state of its sector, with a gap.

**It falls down in four places:**

1. **There is no neutron.** At its ceiling a fiber always holds an
   even number of quanta (the handshake lemma), so a "complete" object
   always has integer spin. Spin ½ requires an open slot. With charge
   read as the deficit, a neutral spin-½ object cannot exist. This
   holds for any fiber graph.
2. **It does not walk by conversion.** At t = 0 its position is
   exactly conserved. It moves only if field transport is allowed
   inside the lock. With transport it forms a spin-½ band (0.14 J
   wide at t = 0.1 J, 0.45 J at 0.3 J), and the ground state stays
   spin ½ up to t ≈ 9 J. Beyond that it polarizes the fiber (S = 9/2
   on the six-ring).
3. **It cannot form under V.** A hole is made when a photon and a lock
   quantum convert at the boundary and the space leaves (12–25% per
   resonant photon). That needs the boundary link to hold space as it
   leaves, which V forbids.
4. **It has no ladder.** Every excitation above it is the hole plus a
   free lock quantum, and the hole–quantum channel is repulsive. No
   bound excited level appears at any length tried.

## 2. The two ideas

### F1. V against the S-law (`f1_v.out`)

| ring | S-law states | V states | lowest 30 levels, max \|diff\| | relative to the first gap |
|---|---:|---:|---:|---:|
| n = 6 | 1831 | 141 | 0.0095 J | 1.3% |
| n = 8 (S-law Σc ≤ 3) | 21827 | 1107 | 0.029 J | 4.9% |
| n = 6, one hole | 928 | 96 | 0.019 J | 1.6% |
| n = 8, one hole (Σc ≤ 3) | 12273 | 750 | 0.032 J | 3.5% |

- **V is the S-law's low spectrum without the space-carrying states.**
  The n = 8 differences include the S-law's own truncation at
  Σc ≤ 3, which is −0.9% in J (`../lock/` amendment).
- The failures that ran through real link occupation have no states
  to run through under V:
  - resonant heating into converted levels;
  - space flow on lock links;
  - T4's collapse under space flow;
  - L6's 10–20% space weight.

  That is true by construction, not measured. What is measured is that
  removing those states leaves the low spectrum intact.

### F2. Counting lock quanta (`f2a_*.out`, `f2b.out`)

**(a) One photon on the n = 4 lock in its ground state.** Final states
per photon:

| g/Δ | boundary | source | ground | **one lock quantum** | more | space-carrying (V's mismatch) | hole made (space left) |
|---:|---|---|---:|---:|---:|---:|---:|
| 0.1 | no flow | resonant | 0.836 | **0.120** | 0.041 | 0.002 | 0 |
| 0.1 | flow | resonant | 0.869 | 0.008 | 0.003 | 0.0002 | **0.121** |
| 0.2 | no flow | resonant | 0.599 | **0.283** | 0.106 | 0.012 | 0 |
| 0.2 | flow | resonant | 0.663 | 0.062 | 0.023 | 0.003 | **0.249** |
| 0.1 / 0.2 | either | broadband | 0.98–0.998 | 0.001–0.011 | | | 0–0.014 |
| any | either | cold | 1.0000 | 0 | 0 | 0 | 0 |

- **About 70% of every non-ground outcome is exactly one lock
  quantum.** Under Q, the resonant photon that `../options/` E5 called
  "heating" is a clean single-quantum source.
- Cold light does nothing to a lock in its ground state, as it should.
- The last column is hole creation. It returns in §3 (F4e).

**(b) Sudden starts** (effective rings, reading V), in lock quanta per
cell (energy above the ground state over n Δ_∞):

| start | n = 8 | n = 10 |
|---|---:|---:|
| Néel | 0.98 | 0.96 |
| Jones 1–3 | 2.4–4.0 | 2.5–3.1 |

A lock quantum's own extent is about ξ ≈ 6 grains. One to four quanta
per cell is not a gas of separate quanta; it is an overlapping soup.
**Under Q, separate quanta come only from gentle, one-at-a-time
processes,** such as (a). A sudden start never makes them.

## 3. The candidates

### F3. The lock quantum (magnon)

From `../lock/` and F2:

- it is an excitation of the fabric (N1);
- it is massive and single: Δ_∞ = 0.436 J (N4);
- it walks by conversion (N6, L3);
- it is made one at a time by resonant light (N8, F2a).

**It is a triplet: N3 fails, by symmetry.**

### F4. The hole (one open slot)

**(a) Rings** (effective, n = 6 … 14; `f4a_rings.out`):

| n | ground S | first excitation (S = 3/2) | one-quantum threshold | second S = ½ |
|---:|---:|---:|---:|---:|
| 6 | ½ | 1.175 J | 0.738 J | 2.07 J |
| 10 | ½ | 0.775 J | 0.539 J | 1.32 J |
| 14 | ½ | 0.625 J | 0.474 J | 0.95 J |

- The ground state is a doublet at every n. Lieb–Mattis predicts this
  on a bipartite ring.
- The S = 3/2 level sits above the threshold by 0.44 → 0.15 J,
  shrinking roughly as 1/n. On the chains the excess is 0.097, 0.05 and
  0.031 J at L_bulk = 24, 40 and 56. That is a repulsive scattering state of the
  hole and one lock quantum, not a bound level.

**(b) Long chains** (effective, DMRG χ = 64, screened ends, hole in the
middle; `f4b_analyze.out`):

| L_bulk | decay length of the hole's cloud | moment within d ≤ 3 | residual at the ends | E(3/2) − E(½) |
|---:|---:|---:|---:|---:|
| 24 | (fit range too short) | +0.663 | −0.074 | 0.592 J (threshold 0.495) |
| 40 | **5.91** | +0.661 | −0.019 | 0.510 J (threshold ≈ 0.46) |
| 56 | **5.80** | +0.661 | −0.005 | 0.478 J (threshold ≈ 0.447) |

- **N2 passes: 1.9% change.** At L_bulk = 40 and 56, the profile agrees
  to three decimals out to d = 20.
- The hole cell itself carries −0.147, against the doublet.
- The neighbours carry +0.794, and a staggered cloud decays over about
  ξ.
- This is the integer-spin fabric's own correlation length dressing a
  missing quantum. It is a size set by the fabric, not by cells placed
  for it. It is not read as a cell count.

**(c, d) Microscopic one-hole ring, n = 6, with transport t**
(`f4cd_micro.out`):

| t | t/J | ground S | band width W of the 6 positional states | lowest S ≠ ½ |
|---|---:|---:|---:|---|
| 0 | 0 | ½ | **0 (exactly degenerate)** | 3/2 at 1.171 J |
| 0.1 J | 0.1 | ½ | 0.141 J | 3/2 at 1.107 J |
| 0.3 J | 0.3 | ½ | 0.451 J | 3/2 at 1.007 J |
| J | 1 | ½ | ≥ 0.60 J (the sixth member lies above the 3/2 level) | 3/2 at 1.035 J |
| 3 J | 3 | ½ | (mixes with 3/2) | 3/2 at 0.489 J |
| 0.1 Δ | 5.2 | ½ | | 3/2 at 0.628 J |
| 0.17 Δ | 8.8 | ½ | | 3/2 at 0.091 J; 7/2 at 0.174 J |
| 10 J | 10 | **9/2** | | — |

- **At t = 0 the hole cannot move (N6 fails),** as `PREREG.md` §3.3
  said.
- With small t it is a mobile spin-½ band.
- Near t ≈ 9 J, high-spin states crowd in. By 10 J the hole has
  polarized the whole ring: the Nagaoka mechanism, seeded by a real
  hole.
- The microscopic t = 0 ground doublet and its 3/2 level agree with
  V's (F1, F4a).

**(e) Formation.** A resonant photon makes a hole with probability
0.121 (g/Δ = 0.1) or 0.249 (0.2), when the boundary link's space can
flow away (F2a). Under V it cannot, so the rate is **0. N8 fails under
V.**

## 4. The scorecard

| | criterion | lock quantum | hole |
|---|---|---|---|
| N1 | from the fabric | ✔ an excitation | ✔ a defect (one missing quantum) |
| N2 | one size | ✔ | **✔ 5.91 → 5.80 grains** |
| N3 | spin ½ | ✘ (triplet) | **✔** |
| N4 | firm (lowest of its sector, gap > 0.1 J) | ✔ | ✔ the first excitation is hole + one quantum: 0.51 J at L_bulk = 40, tending to Δ_∞ = 0.44 J |
| N5 | ladder | — | ✘ no bound level; repulsive channel |
| N6 | translation by conversion at t = 0 | ✔ | ✘ pinned exactly; moves for 0 < t ≲ 9 J |
| N7 | p and n, one lock with two occupancies | ✘ | ✘ **parity: no neutral spin-½ object** |
| N8 | forms under V | ✔ (12–28% per resonant photon) | ✘ under V; 12–25% if space may leave |

**Ending B.** The hole passes N1–N4. It needs a cost for N6 (transport
inside the lock) and for N8 (space leaving the boundary link). N7
fails by algebra.

## 5. Where it falls down, as walls

1. **Parity.** At its ceiling the fabric holds an even number of
   quanta, so spin ½ is always "one open slot". If charge is the
   deficit, the proton has a place and the neutron does not. The wall
   is algebra; no rule tested here moves it.
2. **Pinning.** R6's local law conserves a hole's position exactly. A
   hole that walks needs field transport inside the lock: the margin
   0 < t ≲ 9 J, or about 0.17 Δ. That reopens S4, which `../options/`
   closed with t = 0 for robustness. The two pull against each other.
3. **Formation against V.** The only hole-making process found needs
   space to leave through a link, which is the occupation V removes.
4. **No ladder.** The hole has no excited level of its own; above it
   is only hole plus free quantum.

## 6. One re-reading that would move wall 1 (not tested, the user's call)

The fiber's SU(2) is the Jones (polarization) symmetry. It was read as
spin. **Read as isospin instead,** several things line up:

- the hole's ground doublet is one object with two occupancies: which
  of the two planes holds the open slot;
- that is `v98/UNFOLD.md` §4's isospin row, "two occupancies of one
  pinch";
- the lock quantum is an isospin triplet;
- the S = 3/2 channel above the hole is hole plus quantum, and it
  decays.

This is an observation of shape only, at the level of UNFOLD's
classes. It is not a claim.

Its costs:

- spin ½ would then have to come from somewhere else (v98 says "the
  handedness of the collapse fiber"). A 1-D fiber model has nothing to
  supply it;
- it conflicts with v98's charge row ("open slot vs complete"), since
  both members would have an open slot;
- the p/n split would need a small symmetry breaking between the two
  planes.

## 7. What is now true, and what is not

**True (measured in models):**

1. Reading V is the S-law's low spectrum without its space-carrying
   states (≤ 0.03 J).
2. Under Q, one resonant photon makes exactly one lock quantum 12–28%
   of the time. Sudden starts make an overlapping soup, not separate
   quanta.
3. One missing quantum in the fiber is a firm spin-½ object. Its size
   is set by the fabric (decay length 5.8 grains, the same at two
   lengths), and it has a gap.
4. At t = 0 it cannot move. With 0 < t ≲ 9 J it is a mobile spin-½
   band. Beyond that it polarizes the fiber.
5. It has no bound excited level.
6. By parity, no complete (neutral) object in the fabric has spin ½.

**Not true:**

- A neutron, as a neutral spin-½ object of this fabric.
- A hole that walks by conversion alone, or forms under V.
- A discrete excited ladder.
- Anything in a kernel.
- Any mapping onto measured nucleon numbers.

## 8. What would move past this page (each is the user's call)

1. **Charge.** Keep "charge = deficit" and accept that the fabric's
   spin-½ object is proton-like only. Or re-read the fiber SU(2) as
   isospin (§6), and find spin ½ elsewhere.
2. **Hole transport.** Allow transport inside a lock (0 < t ≲ 9 J). The
   hole then walks as a spin-½ band, and S4 changes from "t = 0" to a
   bounded hop.
3. **Hole formation.** Allow space to leave through a lock's boundary
   link, a limited exception to V. Or find another hole-making process.
4. **The ladder.** A Δ-like level would need an attractive
   hole–quantum channel. Nothing in the S-law at g/Δ = 0.1 supplies one.

## 9. Deviations from `PREREG.md`

1. **F4(b)'s L_bulk = 24 decay-length fit** has only three points
   (d = 2–4) clear of the ends, so it is not used for N2. N2 uses 40
   and 56, as registered ("the two longest").
2. **The F2(a) final-state classes** (one quantum = S = 1 below the
   two-quantum threshold) were set in the code before it ran, and
   follow the registered wording. The S labels come from recurrence of
   each energy across S^z sectors.
3. **F4(d)'s band width W** is the spread of the six lowest S^z = ½
   states. From t ≈ 3 J up, S = 3/2 levels enter that window, so W is
   reported only while all six are S = ½.
4. **F2(b)** uses the effective rings (reading V), as registered. Its
   unit Δ_∞ = 0.4357 J is the infinite-chain gap from `../lock/` L4.

## Reproduce

From `v102/fabric/cas/`:

| step | command | output |
|---|---|---|
| F1 | `python3 -I f1_v.py` | `out/f1_v.out` |
| F2 | `python3 -I f2_quanta.py b`; `python3 -I f2_quanta.py a {0.1,0.2} {F0,F1}` | `out/f2b.out`, `out/f2a_*.out` |
| F4(a) | `python3 -I f4a_rings.py 14` | `out/f4a_rings.out` |
| F4(b) | `./run_f4b.sh`, then `python3 -I f4b_analyze.py` | `out/f4b_*.out`, `out/f4b_analyze.out` |
| F4(c, d) | `python3 -I f4cd_micro.py` | `out/f4cd_micro.out` |

Instruments:

- `eff.py`: reading V. The effective chain in an S^z sector, with the
  hole-centred (2,1,2) term.
- `../../options/cas/graph.py`, `e5_cool.py`: carried.
- `../../lock/cas/clusters.py`, `mpo3.py`: carried.
- `../../phase1/cas/dmrg.py`: carried.
