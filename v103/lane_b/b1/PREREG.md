# v103 lane b, step b1 — pre-registration (where classical begins)

**Written:** 2026-10-11, before any run of this step.
**Authority:** the user, 2026-10-11, thought experiment (b): "Can we
use this inversion to help unify quantum/classical mechanics?" Then:
"work on it in three different lanes". The plan is
`../../README.md`.

**A measurement in models; it rules nothing.**

## 0. The question

The inversion suggests two regimes from one field:
- **classical:** the bulk medium, at large occupation;
- **quantum:** a unit deficit, one void.

For that to unify anything, both must be true at once. As the medium's
occupation n̄ grows:
1. the bulk becomes a classical field, coherent in phase across the
   system;
2. one missing quantum stays a discrete thing: here or there, not
   smeared into the bulk.

Known physics gives a reason to expect a trade-off. Number and phase
are conjugate, so a medium that locks its number (and keeps voids
sharp) gives up its phase, and a phase-coherent medium dissolves a unit
deficit. This step measures whether that trade-off holds, and whether
any ingredient escapes it.

## 1. The medium

- **Bosons on a periodic ring** of L cells, with filling n̄:
  H = −t Σ_{i,μ} (b†_{i+1,μ} b_{i,μ} + h.c.) + (U/2) Σ_i N_i(N_i − 1).
- **Units.** t = 1, and there is no twist.

## 2. The candidates

| | candidate | on | off |
|---|---|---|---|
| **a** | **number lock (saturation)** | cap = n̄: the medium is full | no cap |
| **b** | **stiffness (the hat)** | U = 4t, the free kinetic bandwidth: neither term dominates by construction | U = 0 |
| **c** | **two-plane field** | M = 2, SU(2)-invariant | M = 1 |

**Ladders** (a measurement of trend, not a scan):

| ring | rule sets | n̄ |
|---|---|---|
| L = 4 | every rule set | 1, 2, 3, 4, 6, 8 |
| L = 6 | M = 1 sets (c off) | 1, 2, 3, 4, 6, 8 |
| L = 6 | M = 2 sets | 1, 2, 3 (a 4 × 10⁶-state sector at n̄ = 3) |

## 3. Measurements (`../../tools/bx.c`, mode b1)

- **Medium phase coherence:** g1 = Σ_μ ⟨b†_{i,μ} b_{i+L/2,μ}⟩ / n̄,
  averaged over translations, in the ground state of N = n̄L. A
  classical field has g1 → 1.
- **Medium number fluctuation:** σ_N/n̄ per cell.
- **Void discreteness**, from the half-ring number N_half (averaged
  over every origin):
  β = [Var₁(N_half) − Var₀(N_half)] / [p(1 − p)].
  - Var₀ is taken in the medium's ground state.
  - Var₁ is taken in the lowest state with one quantum removed.
  - p = ⟨N_half⟩₀ − ⟨N_half⟩₁ is the deficit's share in the half
    (= ½).
  - β = 1: the missing quantum is a discrete thing, here or there.
  - β = −1: it is spread through a coherent bulk.
- **Void quasiparticle weight:** Z at the void band's minimum, as in
  a1.

## 4. Gates

These are the same tool gates as a1 (`../../tools/out/t0.out`): free
fermions, free bosons, Hermiticity, SU(2), Lanczos and spectral routines
against dense diagonalization, at 1, 4, 8 and 16 threads.

## 5. Bars and endings

A point (rule set, L, n̄) is:
- **classical-bulk** if g1 ≥ 0.9;
- **discrete-void** if β ≥ 0.5;
- **two-regime** if both.

| ending | condition |
|---|---|
| **U**, unified | some rule set is two-regime at the top of its ladder, at L = 4 and at L = 6 where run. The classical limit keeps a discrete void |
| **P**, partial | some point is two-regime, but none at the top of a ladder |
| **C**, complementary | no point is two-regime: wherever the bulk is classical the void is smeared, and wherever the void is discrete the bulk is not classical |

The trends of g1 and β with n̄ are reported for every rule set.

**What I expect** (written before running):
- **a on:** g1 = 0 and β = 1 (number-locked).
- **a off, b off:** g1 = 1 and β = −1, a free condensate.
- **b** moves between the two.
- **c** changes little.

The expected ending is **C**.

## 6. Disclosed before running

- **Small rings.** L = 4 has a half-ring of two cells. g1 at L/2 is
  the longest distance there is.
- **U = 4t is a registered natural unit,** not tuned to place a
  crossover.
- **This step uses a quantum medium.** Whether a medium with no ħ put
  in can produce quantum voids is a later step. There, Bell is the
  wall.

- **Tool faults the gates caught, all fixed before any registered
  run:**
  1. **An off-by-one in the Lanczos step count.** Every restart ended
     after one pass, and a one-step pass handed LAPACK a size of 0.
  2. **Fragile restarts without reorthogonalization.** With threaded
     reductions, a restart could lock onto an excited state. Lanczos
     now reorthogonalizes fully, small sectors (≤ 2500 states) are
     diagonalized densely, and only the residual ends a restart.
  3. **An absolute breakdown threshold.** Once the Krylov space inside
     a probe's momentum sector was exhausted, roundoff carried the
     recurrence into other sectors, and a spurious, weightless lower
     pole was reported. The threshold is now relative, and the lowest
     pole is the lowest one that carries weight.

  Gates 5 and 6 found these faults. The final binary passes every gate
  at 1, 4, 8 and 16 threads (`../../tools/out/t0.out`).

## 7. Reproduce

From `cas/`: `./run_b1.sh`; tables from `go run tab.go`.
