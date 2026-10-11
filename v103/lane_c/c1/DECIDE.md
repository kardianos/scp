# v103 lane c, step c1 — decision (voids in a chain medium)

**Written:** 2026-10-11, after the registered runs. Gates: `PREREG.md`.
A measurement in models; it rules nothing.

- Tool: `cas/chain.c` (C). Runs: `cas/run_t0.sh`, `cas/run_grid.sh`.
- Fits: `cas/fit.go` (Go), output in `cas/out/fit.out`.
- 40 runs: eight rule sets × {2-D, 3-D} × {two voids, four voids},
  plus the 1-D ring for the four bipartite sets. Each run is 8 seeds ×
  10⁶ sweeps.

## Ending: **E** (charge-like)

**Chains re-forming over geometric distance do give the long-range
law.** In the exact ground state of a medium that re-forms chains,
two voids interact as Coulomb charges on that very lattice:
- logarithmic in 2-D;
- 1/r in 3-D;
- (in 1-D, the strict medium confines them).

**But the law is signed, not universal.** A void carries the sign of
the sublattice its missing quantum came from. Unlike voids attract and
like voids repel. So in this model the user's mechanism produces
**electric-like charge, not gravity.**

| | needed? | evidence |
|---|---|---|
| **b**, chain conservation (the lack must recurse) | **yes** | every b-off set is short range, in both dimensions |
| **c**, bipartite fabric | **yes** | b on a non-bipartite fabric (kagome, pyrochlore) is short range |
| **a**, exclusion per link end | no, but it strengthens the law | coupling β: 1.43 → 2.10 (2-D), 2.44 → 2.90 (3-D) |

## 0. Tool gate T0: **PASS**, with one disclosed statistical fail

- **Lattices and Green's functions:** all pass. Links N·z/2, the
  2-colouring as registered, residuals ≤ 1.4e-13.
- **Gate 1, Monte Carlo against exact enumeration:** 17 of 18 cases
  pass. Cubic 2×2×2 under bc gave χ²/dof = 1.84, against a bar of 1.5.
- **Gate 1b, added before production:** the sampler's full transition
  matrix, on every enumerable state space. Every column sums to 1 to
  ≤ 1.2e-15, in all 15 cases including that one, so the uniform measure
  is exactly stationary.
- **The failing case, rerun** with independent seeds and 4× the
  statistics: χ²/dof = 1.08.
- **Ledgers:** zero violations in every production run.

Recorded in `PREREG.md` §6 before production.

## 1. The pair potential (two voids)

V(r) = −ln[P(r)/P_ideal(r)], fitted on shells with r ≥ 2:
- F0: a constant;
- F1: Coulomb on the same lattice, V = α + β·Ĝ(r);
- F2: screened Coulomb.

| rules | 2-D | 3-D |
|---|---|---|
| none, a, b, c, ab, ac | short: F1 never improves on F0 by Δχ² ≥ 25 | short |
| **bc** | **long**: Δχ² = 89 680, χ²/dof 1.46, β = −1.429 ± 0.005 | **long**: Δχ² = 1129, χ²/dof 0.87, β = −2.44 ± 0.07, best screening m = 0.2 not significant (Δχ² 1.4) |
| **abc** | **long**: Δχ² = 248 419, β = −2.097 ± 0.004 | **long**: Δχ² = 1644, χ²/dof 0.96, β = −2.90 ± 0.07, m = 0 |

- β < 0 means the voids (always on opposite sublattices in a
  two-void state) attract.
- The 3-D profile under abc: V = −0.195 at r = 1, −0.071 at √3, −0.025
  at 3, 0.000 at 4.6, +0.015 at 8. It rises monotonically, as the
  lattice's 1/r does.

## 2. The sign (four voids)

Pair correlations g at the two nearest allowed shells:

| rules | d | g_same | g_opp | verdict |
|---|---|---|---|---|
| bc | 2-D | 0.364 (r = 0), 0.682 (√2) | 1.623 (1), 1.333 (√5) | **signed** |
| bc | 3-D | 0.466 (0), 0.914 (√2) | 1.176 (1), 1.059 (√3) | **signed** |
| abc | 2-D | 0.238, 0.574 | 1.954, 1.510 | **signed** |
| abc | 3-D | 0.382, 0.894 | 1.220, 1.077 | **signed** |
| every other set | both | contact suppressed (g ≈ 0.4–0.56 at r = 0), slight excess at r = 1 | (as same) | none |

All errors are ≤ 0.004. Like-sublattice voids repel and
unlike-sublattice voids attract: a two-sign Coulomb gas. No set shows
universal attraction.

## 3. Anemia (b off; reported)

In the sets where a lack can stay local, the number of unbonded quanta
per cell rises next to a void and then heals within one or two cells:

| set (3-D) | at r = 1 | at r ≥ √3 |
|---|---|---|
| none | 0.606 | 0.574–0.576 |
| ac | 0.645 | 0.590–0.600 |

The anemia is local. It does not carry the void's influence outward.

## 4. 1-D control

| rules | result |
|---|---|
| abc | the two voids are always adjacent: confined, which is the 1-D Coulomb law at infinite slope |
| bc | separations are flat: allowing two bonds per link frees the chain |
| c, ac | flat |

## 5. What I expected (written before running), against what happened

| expectation (PREREG §0, §5) | outcome |
|---|---|
| ending E, carried by chain conservation on a bipartite fabric | **held** |
| a not needed | **held**, though a strengthens the coupling |
| non-bipartite and b-off sets short range | **held** |

## 6. What lane c now says

1. **The user's mechanism is real in this model, and it gives a 1/r
   law.**
   - Make the lack recurse: chains may end only at voids.
   - Lay the medium on a fabric with two sublattices.
   - Then a void's disturbance spreads geometrically, and two voids
     feel Coulomb's law in exactly the dimensions predicted: log in
     2-D, 1/r in 3-D.
2. **The law it gives is charge, not gravity.**
   - Every void is the same kind of lack, but a void inherits a sign
     from where it sits.
   - Neutrality is enforced: a lone like pair cannot exist without
     partners.
   - In this reading, **the two sublattices are the two signs of
     charge**. An antiparticle is the same void on the other
     sublattice, or an excess on the same one.
3. **Anemia does not reach.** Where the lack can be absorbed, it
   heals within about one cell. Nothing unsigned propagates.
4. **For gravity, the unsigned part must come from somewhere else.**
   The candidates already named are:
   - the steady draw: a void that converts substrate into an outgoing
     mode, with net intake;
   - or a coupling to an unsigned quantity that is itself conserved
     and long range.

   Here, the only long-range field is the signed flux.

## 7. Limits

- **Sizes.** 2-D 32² (square) and 3-D 12³ (cubic), with kagome and
  pyrochlore of similar size. The 3-D fit uses 28 shells.
- **Equilibrium only** (the RK ground state). No dynamics, no draw.
- **q = 2 (chains) only.** q = 1 (dimers) is known to behave
  similarly, but was not run.

## Then

- **c2, the steady draw.** A void that converts chain into an outgoing
  mode, with a conserved substrate flowing in. Does an unsigned 1/r
  deficit appear, and does it act on other voids universally?
- **The sign as charge.** Use lane a's tools to ask whether a
  sublattice-charged void also carries spin ½ and a wave, which would
  make it a charged quantum particle.
