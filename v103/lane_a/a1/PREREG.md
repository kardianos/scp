# v103 lane a, step a1 — pre-registration (the void as a quantum particle)

**Written:** 2026-10-11, before any run of this step.
**Authority:** the user, 2026-10-11, thought experiment (a): "Can
existing quantum physics survive if we invert what 'stuff' is?" Then:
"work on it in three different lanes". The plan is
`../../README.md`.

**A measurement in models; it rules nothing.**

## 0. The question

In known quantum mechanics, particles are excitations of an empty
vacuum. Invert that: the vacuum is a full medium, and "stuff" is a
void, one quantum missing. Does the void carry the properties quantum
mechanics needs a particle to carry? The ledger asks five things:

| | property | what it means for a void |
|---|---|---|
| **A1** | a wave | a band with sharp quasiparticle weight, and Aharonov–Bohm interference around the ring |
| **A2** | exclusion | two voids do not share a cell |
| **A4** | spin ½ | in a singlet medium, the void carries spin ½ |
| **A5** | an antiparticle | an excess (one extra quantum) exists as a sharp particle |
| **A6** | a breathing line (the Higgs reading) | the medium's own structureless breathing has a sharp line below the void-plus-excess continuum |

A3 (exchange statistics) needs 2-D and is deferred to a later step.

## 1. The medium (known quantum mechanics)

- **Bosons on a periodic ring** of L cells (no edge), with filling
  n̄ = 2 per cell.
  - n̄ = 2 is R6 on a ring: one quantum per link end, and every cell
    has two links.
  - It is a standard lattice quantum model:
    H = −t Σ (e^{iθ} b†_{i+1} b_i + h.c.) + (candidate terms).
- **Units and twist.** t = 1 sets the unit. The twist θ = φ/L threads
  a flux φ through the ring.
- **A void** is the state left when one quantum is removed:
  b_k|G⟩, where |G⟩ is the medium's ground state. Nothing labels it
  (R13).

## 2. The candidates

| | candidate | on | off |
|---|---|---|---|
| **a** | **exclusion (saturation)** | cap = n̄ = 2: the medium is full, with no headroom | cap = 3: one quantum of headroom per cell |
| **b** | **two-plane field** | two modes per quantum (M = 2), SU(2)-invariant (R1's field) | one mode (M = 1) |
| **c** | **conversion exchange** | −J A†A on every link. A = Σ P b_i b_j: P = ε (the singlet) for M = 2, P = 1 for M = 1. J = J0/2 = 5t, from R9's ratio t = 0.1 J0 and J0 = 2g²/Δ | none |

There is no on-site interaction (U = 0). Every combination is run, at
L = 6 and L = 8. The largest sector is 2.2 × 10⁶ states (L = 8,
M = 2, cap 3).

## 3. Measurements (`../../tools/bx.c`, mode a1)

| | measurement | how |
|---|---|---|
| 1 | **the medium** \|G⟩ | lowest state of the N = n̄L sector (2S^z = 0). Its gap to the next state, ⟨S²⟩, and its number fluctuation |
| 2 | **the void band** | for every k, the lowest pole E_v(k) of the spectral function of b_{k↑}\|G⟩, and its weight Z(k) |
| 3 | **bandwidth and minimum** | W = max E_v − min E_v; Z at the band minimum |
| 4 | **void spin** | ⟨S²⟩ of the lowest-pole state at the band minimum |
| 5 | **Aharonov–Bohm** | dAB = [E_{N−1}(π) − E_N(π)] − [E_{N−1}(0) − E_N(0)], the lowest energies of each sector |
| 6 | **the excess** (a off only) | the lowest poles of b†_{k↑}\|G⟩, their weights, and Z at the minimum |
| 7 | **breathing** | the lowest pole of O\|G⟩, with O = Σ_i (N_i − n̄)² − ⟨O⟩ and the ground state projected out (k = 0, spin 0, an amplitude probe), plus its weight. The threshold is min_k [ε_v(k) + ε_ex(k)] (a void and an excess of total momentum 0), and R = E_breath / threshold |
| 8 | **two voids** | the lowest state of the N − 2 sector: its binding against two free voids, and P_double = ⟨Σ_i [N_i ≤ n̄ − 2]⟩ |

## 4. Gates (T0, before any registered run)

All are in `../../tools/out/t0.out`, run at 1, 4, 8 and 16 threads.

1. Hard-core bosons (M = 1, cap 1) on rings L = 4 … 10, at every
   filling and two fluxes, equal free fermions, with the boundary twist
   set by parity, to 1e-8.
2. Free bosons: E0 = −2tN cos(φ/L), to 1e-8.
3. Hermiticity: ≤ 1e-12, with random vectors.
4. SU(2): ‖[H, S⁺]x‖/‖S⁺Hx‖ ≤ 1e-12 (M = 2, with U and J).
5. Iterative (Lanczos) lowest energies against dense diagonalization,
   to 1e-9. This covers every shape of rule set at L = 4, with
   sectors up to 1500 states.
6. The spectral routine's lowest pole and weight against dense
   diagonalization (b_k|G⟩, every rule-set shape at L = 4): energy to
   1e-8, weight to 1e-6.

## 5. Bars and the ledger

A property **survives** for a rule set only if its bar holds at both
L = 6 and L = 8.

| | survives when |
|---|---|
| **A1** | the medium is non-degenerate (gap ≥ 1e-6), W ≥ 1e-3, Z(band minimum) ≥ 0.5, and \|dAB\| ≥ 1e-8 |
| **A2** | P_double ≤ 0.01 |
| **A4** | M = 2, the medium is a non-degenerate singlet (⟨S²⟩ ≤ 1e-6), and the void's ⟨S²⟩ = 0.75 ± 1e-6 |
| **A5** | an excess exists with Z(band minimum) ≥ 0.5 |
| **A6** | R < 1, with breathing weight ≥ 0.1 |

**Endings:**

| ending | condition |
|---|---|
| **F**, full survival | some rule set passes A1, A2, A4, A5 and A6 |
| **P(k)**, partial | the best rule set passes k of the five. The ledger names which properties break, under which ingredients |

**What I expect** (written before running):
- **Saturation (a on)** gives sharp voids (A1), but no excess (A5
  breaks), and the full medium cannot breathe (A6 breaks:
  O|G⟩ = 0).
- **Headroom (a off)** allows an excess, but softens the void.
- **Spin ½ (A4)** needs b and c: a two-plane field in a singlet medium.
- **Exclusion (A2)** is not given by saturation at n̄ = 2, because a
  cell can hold zero quanta.

So I expect P(3) or less. The tension would be that the saturation
that makes voids sharp is the same thing that forbids antiparticles
and breathing.

## 6. Disclosed before running

- **Small rings** (L = 6, 8). Small Z or W differences near a bar are
  finite-size effects.
- **J = 5t is not a scan.** It is R9's ratio, the only standing
  ruling that relates the two scales.
- **Removing an up quantum** puts the void in the 2S^z = −1 sector. Its
  spin is read from ⟨S²⟩.
- **A5 and A6 cannot pass under a on**, by construction: no excess
  state exists, and O|G⟩ = 0. This is part of what is measured, and it
  is reported, not hidden.

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

From `../../tools/`:
- `gcc -O3 -march=native -fopenmp -Wall -Wextra -o bx bx.c -l:liblapack.so.3 -l:libopenblas.so.0 -lm`;
- `./bx mode=t0`.

Then, from `cas/`: `./run_a1.sh`, and tables from `go run tab.go`.
