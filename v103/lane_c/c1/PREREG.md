# v103 lane c, step c1 — pre-registration (voids in a chain medium)

**Written:** 2026-10-11, before any run of this step.
**Authority:** the user, 2026-10-11:
- "Field void interaction.... the more void you have the field becomes
  animic and bonds lack of chains recurse, exclusion principle in field
  over space, The gravity 1/r gradient is a product of chains
  re-forming over geometric distance";
- then: "write this new approach down, then work on it in three
  different lanes".

The plan is `../../README.md`.

**A measurement in models; it rules nothing.**

## 0. The question

Take a medium of field quanta that bond in chains, under exclusion,
and remove a quantum somewhere: a void. Then:
1. Does the lack the void leaves re-form over geometric distance, so
   that two voids feel a long-range law: 1/r in 3-D and log r in 2-D?
2. Is that law **universal**, with every void attracting every other
   as gravity needs? Or is it **signed**, with like voids repelling and
   unlike voids attracting, as charge behaves?
3. Which ingredients does it need?

**What known physics says** (written before running). In constrained
bond models on bipartite lattices, a local conservation rule makes
defects interact as Coulomb charges, which carry the sign of their
sublattice. Without the rule, or on a non-bipartite lattice,
correlations are short-range. So the expected ending is **E**
(charge-like; §5). The test is whether the user's three ingredients
reproduce that, and whether any of them adds a universal part.

## 1. The medium

- **Cells** sit on a periodic lattice: no edge.
- Each cell i holds Q_i field quanta.
  - A normal cell holds q = 2: a cell in a chain has two bonds.
  - A void is a cell with a quantum missing (Q_i = 1). Two missing
    make Q_i = 0.
  - Voidness is read from Q_i, a field count, and nothing else
    (R13).
- **Bonds.** Each link l = (i, j) carries n_l bonds. A bond pairs one
  quantum at each end of the link.
- **The degree** deg_i = Σ_{l∋i} n_l is the number of cell i's quanta
  that are bonded. Always deg_i ≤ Q_i: a quantum bonds at most once.
- **The quantum medium.** Its kinetic term re-forms chains: it moves
  one bond unit along a cell, and moves a quantum between neighbours,
  each with equal amplitude. At its Rokhsar–Kivelson point, the ground
  state is the equal-amplitude superposition of all valid
  configurations. So every diagonal ground-state correlation is the
  uniform-ensemble average computed here.
- **No constant enters.** Every weight is 1.

## 2. The candidates

Each is a switch. All eight combinations are run.

| | candidate | on | off |
|---|---|---|---|
| **a** | **exclusion per link end** ("exclusion principle in field over space") | a link carries at most one bond: n_l ∈ {0, 1} | n_l ∈ {0, …, q} |
| **b** | **chain conservation** ("bonds lack, chains recurse") | every quantum is bonded, deg_i = Q_i. A lack cannot be absorbed where it is, so chains can end only at voids | deg_i ≤ Q_i. A cell may leave quanta unbonded, so the lack can stay local |
| **c** | **bipartite fabric** | square (2-D), simple cubic (3-D), ring (1-D) | kagome (2-D), pyrochlore (3-D). Same coordination as their bipartite partners (z = 4, 6), but with triangles. 1-D has no non-bipartite partner |

## 3. Systems and runs

| d | c on | c off |
|---|---|---|
| 1 | ring, N = 128 | (none) |
| 2 | square 32 × 32, N = 1024 | kagome 18 × 18 cells, N = 972 |
| 3 | simple cubic 12³, N = 1728 | pyrochlore 5³ conventional cells, N = 2000 |

**Ensembles.**
- **E2:** two missing quanta, the void pair.
- **E4:** four missing quanta, for the sign test.

Every (rule set, system, ensemble) runs 8 independent seeds
(splitmix64). The sweep counts are fixed after T0's timing, before any
production run, and recorded in §6.

**Moves.** Every move keeps the uniform measure: a symmetric proposal,
accepted if the result is valid.
- **b on:** the void hop. A missing quantum at cell v takes one more
  bond onto a random link (v, k). Cell k then gives up one of its other
  bond units, (k, l), chosen at random. The missing quantum is now at
  l, two steps away on the same sublattice. Bonds re-form along the
  way. The counts in each direction are equal, so detailed balance
  holds (§4 T0-1 checks it).
- **b off:** two kinds of move, half the time each:
  - a link flip, n_l → n_l ± 1;
  - a quantum hop, a quantum moving from a neighbour k into the void
    cell v.

## 4. Gates (T0, before any registered run)

1. **Exact enumeration against Monte Carlo.** On small periodic graphs,
   every valid (configuration, void-position) state is enumerated, and
   the exact distribution of void-pair positions is compared with the
   sampler's, for every rule set whose enumeration fits:
   - square 2 × 4 and 4 × 4;
   - kagome 1 × 2 and 2 × 2 cells;
   - cubic 2 × 2 × 2 (a multigraph);
   - one pyrochlore conventional cell.

   Pass: χ²/dof ≤ 1.5 and max |z| ≤ 4.5 in every case run. Cases too
   large to enumerate are listed, not run.
2. **Lattices.**
   - Every cell has degree z, there are no self-links, and there are
     N·z/2 links.
   - A 2-colouring succeeds for every c-on lattice and fails for every
     c-off lattice.
3. **Green's functions.** For each lattice, Ĝ solves the graph
   Laplacian L·Ĝ = δ − 1/N to a residual ≤ 1e-10. The screened Ĝ_m
   solves (L + m²)·Ĝ_m = δ to the same residual.
4. **Ledger.** Along every run, Σ_i Q_i and the deg ≤ Q (or deg = Q)
   condition hold at every check (every sweep), with zero
   violations.

## 5. Readout, bars and endings

**E2: the pair potential.**
- Pairs are binned in shells of equal minimum-image distance.
- V(shell) = −ln[P_obs/P_ideal], where:
  - P_ideal counts the position pairs that the rule set allows (only
    opposite sublattices under b·c);
  - V is fixed up to a constant;
  - errors come from the scatter over the 8 seeds.
- Fits on shells with r ≥ 2 nearest-neighbour distances:
  - F0: V = α;
  - F1: V = α + β·Ĝ(r), Coulomb on this very lattice. This gives the
    dimension fingerprint for free: log in 2-D, 1/r in 3-D;
  - F2: V = α + β·Ĝ_m(r), screened, with m fitted on a fixed grid
    (m = 0 is F1).
- **Long range** (per run): F1 beats F0 by Δχ² ≥ 25, F1 has
  χ²/dof ≤ 3, and F2 does not beat F1 by Δχ² ≥ 9 with a screening
  length 1/m < L/4.
- **Short range:** otherwise.

**E4: the sign.**
- Pair correlations g_same(r) and g_opp(r), normalized to ideal:
  - "same" and "opposite" mean sublattice on bipartite lattices;
  - all pairs count as "same" on non-bipartite lattices.
- Read at the two nearest allowed shells:
  - **signed:** g_opp > 1 and g_same < 1, each by ≥ 5σ;
  - **universal attraction:** every g > 1 by ≥ 5σ;
  - **universal repulsion:** every g < 1 by ≥ 5σ;
  - **none:** otherwise.

**Anemia (b off only):** the mean number of unbonded quanta per cell,
against distance to the nearest void. Reported, not scored.

**Endings:**

| ending | condition |
|---|---|
| **G**, gravity-like | some rule set is long range in both 2-D and 3-D, and its E4 sign is universal attraction |
| **E**, charge-like | long range appears, but only with signed E4 |
| **S**, short only | no rule set is long range in both 2-D and 3-D |

Each rule set's necessary ingredients are read from the eight
combinations. 1-D is reported as a control: the pair's spread of
separations, and whether the voids are confined.

## 6. Disclosed before running

- **Small systems.** L ≤ 32 in 2-D and 12 in 3-D. On these, 1/r and a
  long exponential are told apart by the fit, not by eye.
- **The RK point is a choice of point, not a coupling.** It is the one
  point where the quantum ground state is known exactly, and it has no
  free constant.
- **q = 2 (chains) is the user's picture.** q = 1 (dimers) is not run
  here.
- **E4's sign is read only at short range.** Its long-range sign
  follows from E2 if the system is a Coulomb gas.
- **Sweep counts:** written here after T0 and before production:
  (pending).

## 7. Reproduce

From `cas/`:
- `gcc -O3 -march=native -Wall -Wextra -o chain chain.c -lm`;
- `./run_t0.sh`, then `./run_grid.sh`;
- tables: `go run fit.go`.
