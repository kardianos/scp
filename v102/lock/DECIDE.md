# v102 lock — decision (S1–S4 ruled)

**Written:** 2026-10-07, after the outputs in `cas/out/`. Gates and the
amendment: `PREREG.md`. Ruling R7: `../README.md`.
**Status:** this page closes the lock step. It does not edit a kernel or
the law table in the tree, and it does not open a kernel step.

---

## 1. The decision: **ending A**, under the S-law

L1, L2, L3 and L4 all pass. Under S1–S4, a fiber ring at its
ceiling is a working lock:

1. **It keeps its field.** Conversion stays virtual. The pitchless
   control converts for real, so S1 is what makes the hold.
2. **Its quantum walks only by conversion.** It moves at a rate set
   by `J ≈ 2g²/Δ`, and the field quanta never move.
3. **Its quantum is massive and single.** There is a gap, one line,
   and space-carrying states far above.

This holds in a model, not in the kernel. No kernel implements S1–S3.
The battery has not run. The costs owed before the lock can live in
the world are in §3, and two of them are measured there.

## 2. Results by gate

### L1. Spectrum: **PASS**

The first 12 excitations, compared with multiplicity against the
spin-1 Heisenberg ring at g/Δ = 0.1, with one fitted `J`
(`l1_spectrum.out`, `l1_exact8.out`):

| ring | dim (S^z = 0) | J_fit | J_fit/J_h2 | max miss | ground S |
|---|---:|---:|---:|---:|---:|
| n = 6, no truncation | 1831 | 0.009151 | 0.991 | 3.4% | 0 |
| n = 8, Σc ≤ 3 (gate) | 21827 | 0.009086 | 0.984 | 2.9% | 0 |
| n = 8, no truncation (check) | 30654 | 0.009167 | 0.993 | 3.3% | 0 |
| n = 8, Σc ≤ 2 (not converged; unused) | 13283 | 0.008507 | 0.922 | 1.8% | 0 |

- `J_h2 = 0.009231` comes from the two-cell cluster; `2g²/Δ` is
  0.009586.
- **The S-law differs from Phase 1's linear model at fourth order.**
  Its conversion strings are cheap (0.2275 per further adjacent
  conversion). That gives a biquadratic `K/J = −0.019` (Phase 1:
  −0.00016) and a three-site remainder of 0.057 J. The Σc ≤ 2
  truncation misses these strings, which is why it is unconverged.

### L2. The lock holds its field: **PASS**

**(a)** Eight unlocked starts evolved to T = 100/J (about 11,000 rate
time units), with energy conserved to at most 1.4e-14 of the total
(`l2_hold.out`):

| ring | max ⟨c_l†c_l⟩ (bar 0.30) | lowest field fraction | second-half mean (bar ≥ 0.85) |
|---|---:|---:|---:|
| n = 6, exact | 0.079–0.103 | 0.897–0.954 | 0.922–0.962 |
| n = 8, Σc ≤ 3 | 0.081–0.092 | 0.914–0.949 | 0.926–0.957 |

The starts were the Néel product and three random Jones products.
None had space or pairing at t = 0.

**(b) The control.** Refuse S1, so space has no pitch and is a flat
continuum around the pair. The same n = 4 ring then converts for real
at the golden-rule rate.

| t | P(converted) | golden rule |
|---|---:|---:|
| 1/Γ | 0.601 | 0.632 |
| 2.5/Γ | 0.893 | 0.918 |
| 5/Γ | **0.973** | 0.993 |

With S1's oscillator, the same ring's ⟨c_l⟩ never exceeds 0.094 over
[0, 11058]. **S1 is the hold.**

### L3. The quantum walks by conversion: **PASS**

On the n = 10 ring (Σc ≤ 3), a kick `S⁺₀|0⟩` spreads around the ring.
Its spin profile matches the Heisenberg ring's, with time in units of
1/J (`l3_walk.out`):

| g/Δ | J_fit | sup over r, τ ∈ [0, 4] of the deviation (bar 0.05) |
|---:|---:|---:|
| 0.05 | 0.002373 | 0.0046 |
| 0.10 | 0.009024 | 0.0167 |

- The two couplings collapse onto each other to within 0.012.
- At g = 0, the same kicked state does not move (1.6e-12).
- The field never travels. Each cell's occupation leaves 2 only by
  exactly twice its links' space occupation. That is the local law
  `N_i + c_{i−1} + c_i = 2`, which is R6 at t = 0.

This is v98 §2's sentence made literal: "the lock walking by
conversion, not a store sliding". The lock's quantum moves because
virtual pair conversion swaps polarization across a link. The quanta
that carry it stay in their cells.

### L4. Massive and single: **PASS**

**(a) A mass exists.** Effective chain from the S-law clusters
(h2 + r3), with R6's spin-½ end cells, χ = 64
(`l4_collect.out`):

| L_bulk | S-law gap / J_h2 | Heisenberg gap / J | ratio |
|---:|---:|---:|---:|
| 24 | 0.4951 | 0.4674 | 1.059 |
| 32 | 0.4713 | 0.4486 | 1.051 |
| 48 | 0.4506 | 0.4310 | 1.046 |
| **∞**, fit `D + a/(L+1)²` | **0.4357** | 0.4190 | 1.040 |

- **Δ_∞/J_h2 = 0.436 > 0.2: PASS.** In rate units the lock quantum's
  rest energy is Δ_∞ = 0.0040 at g/Δ = 0.1.
- **(c) Reported.** The ratio is 6.1% above 0.4105, and 4.0% above the
  Heisenberg chain run by the same procedure. That procedure itself
  sits 2.1% above the literature at these lengths; P2's longer chains
  gave 0.6%.
- The excess is the S-law's fourth-order terms (K/J = −0.019,
  r3 = 0.057 J). It is the gap's counterpart of P3's ξ drift. It
  shrinks as g/Δ → 0, and it is not a choice.
- Convergence: χ = 96 at L = 48 moves E₀ by 1e-7 J and the gap by
  2e-6 J (0.4505782 against 0.4505800).
- The DMRG MPO was rebuilt as an exact compact automaton (`mpo3.py`,
  bond dimension 21 instead of up to about 156). It was checked
  against the term-by-term builder and against the dense Hamiltonian
  to 6e-15. The physics is unchanged; only the speed is.

**(b) One line** (`l4b_line.out`, n = 8, Σc ≤ 3), with the n = 8 gap
Δ₈ = 0.005549:

- the lowest S^z = 2 state sits at 1.47 × 2Δ₈;
- the lowest non-triplet singlet sits at 2.16 × 2Δ₈;
- two levels lie below 2Δ₈ (the ground state and the magnon), and
  neither carries space (total ⟨Σc⟩ is 0.34 over 8 links).

Reported for the n = 6 ring without truncation: the first level
carrying a real space quantum is index 141, right after all 141
spin-manifold states. It sits at 0.91 Δ = 32 × 2Δ₆.

## 3. The costs, two of them measured (reports L5–L7)

### L5(a) Space transport

The table moves space (`s_k`). S1's oscillator survives it only while
the space band stays above the pair. One bond, singlet start,
g/Δ = 0.1 (`l5_costs.out`):

| t_s/Δ (1-D lead) | band edge − pair | bound field weight | late field survival |
|---:|---:|---:|---:|
| 0 | +0.479 | 0.949 | 0.903 |
| 0.3 | +0.192 | 0.934 | 0.872 |
| 0.45 | +0.048 | 0.890 | 0.791 |
| 0.50 | 0.000 | 0.836 | 0.698 |
| 0.55 | −0.048 | 0.410 | 0.157 |
| 0.60 | −0.096 | 0 | 0 |

The condition is exact: **the lock keeps a bound field state while the
space half-bandwidth stays below Δ.** Reading `t_s = s_k = 0.06` (a
diffusive conductance read as an amplitude, so a reading):

| lead | outcome |
|---|---|
| 1-D fiber | holds (edge +0.36 above the pair) |
| tree with the foam's link coordination z_l = 12.6 | holds marginally (edge +0.07) |
| lattice-edge estimate, z_l t_s = 0.75 > Δ | decays |

So S1 and the table's space flow are compatible only marginally, and
only under some readings.

### L5(b) The transport margin

On the n = 6 ring with the S-law energies, the lock gives way at
**t/Δ ∈ (0.16, 0.18]**, under both truncations. Beyond that point:

- the ground state is the converted ferromagnet, S = 3–4;
- Phase 1's linear model gave (0.18, 0.20];
- 3-D field transport is t/Δ = 1.88, about 10× over.

So S4 is not a convenience. Without it there is no lock.

### L5(c) Owed under the ratchet, not computed

1. **Detection.** A lone quantum cannot convert, so a click needs a
   partner (P6). The slit's click claims (R3 atomicity, Born clicks)
   must be re-derived with pair detection.
2. **Two conversion laws.** The kernel's cell door converts field
   above `e_cond` at a rate. Acting on a lock cell, it would convert
   the lock regardless of S1. Locked field must be exempt from it,
   or the door must be replaced by S3.
3. **Space flow.** L5(a)'s condition has to be met by whatever `s_k`
   becomes in the S-law's representation. The g1/g3/g4 claims ride
   on that flow.
4. **Context.** S4 makes transport depend on whether field is locked.
   R6's ceiling counts the lock's link ends, and it is undefined for a
   lock cell's links to cells outside the lock.
5. **`g`.** No value is set. The lock's mass scale is ∝ g²; its
   ratios are set at leading order and drift at order (g/Δ)².
6. **Ends.** A two-cell fiber converts downhill (Δ₁₁ = −0.024), so it
   is not a lock. An end link of a longer open chain has
   Δ_end = 0.2275 and `J_end ≈ 1.96 J`.

### L6. Formation: **there is no channel**

A sudden start leaves the fiber hot (n = 6; `l2_hold.out`):

| start | weight in the spin manifold | internal energy | canonical T |
|---|---:|---:|---:|
| Néel | 0.81 | 3.8 gaps | 1.2 J |
| Jones 1 | 0.86 | 8.4 gaps | 3.5 J |
| Jones 2 | 0.89 | 10.8 gaps | 15.5 J |
| Jones 3 | 0.90 | 11.8 gaps | above the spin-manifold mean (∞) |

The S-law has nothing that could carry this excess away:

- the first space-carrying level lies 0.30 (0.63 Δ) above the top of
  the spin manifold;
- transport is zero;
- space is discrete.

**A lock, once made, keeps its internal energy.** The clean line
exists only in a fiber that was made cold. No process in the S-law
makes one cold.

The obvious route: unlocked light scattering inelastically off the
lock, through virtual conversion across a link from the lock to the
outside, and carrying energy away at the J scale. That route needs
item 4 above (the ceiling and S4 at the lock's boundary) to be ruled
first.

### L7. Size: **the smallest even cycle, as E2 found**

E₀(n) on S-law rings (`l7_size.out`), per cell:

| even n | E₀/n | odd n | E₀/n |
|---:|---:|---:|---:|
| 4 | −0.022872 | 3 | −0.018467 |
| 6 | −0.022291 | 5 | −0.021149 |
| 8 | −0.022118 | 7 | −0.021686 |
| 10 | −0.021941 | 9 | −0.021778 |

- Odd rings are frustrated.
- Even rings are lower per cell the smaller they are, so two 4-rings
  beat one 8-ring (by 0.0060 = 0.66 J), and a 4-ring plus a 6-ring
  beats a 10-ring (by 0.0058).
- That is `v99/reopen/OPTIONS.md` E2's ending, "the ring sits on the
  smallest even cycle or hoards." It is a boundary set by the
  shortest even cycle the geometry allows, not an interior minimum
  set by the law.

The S-law does not select a ring size in the registered sense.

## 4. What is now true, and what is not

**True:**

1. Under the S-law, a fiber ring at its ceiling holds its field
   against conversion at every internal energy. The hold comes from
   S1, the discrete, pitched space oscillator. Remove it and the same
   ring converts at the golden-rule rate.
2. The lock's quantum walks only by conversion, at a rate J ∝ g².
   Field never travels.
3. The quantum is massive and single (L4).
4. The lock exists only inside two margins:
   - transport t/Δ < 0.16–0.18;
   - space half-bandwidth < Δ.

   The table's 3-D transport is 10× outside the first. Its space
   flow sits on the edge of the second, depending on how it is read.
5. The lock cannot cool itself.

**Not true:**

- The lock lives in the kernel, or the battery passes with it.
- `g` is set by anything.
- A ring size is selected (the smallest even cycle is a boundary).
- A cold lock forms.
- ξ, or any length, is a cell count.

## 5. What would move past this page

Each item is the user's call.

1. **Formation.** Rule the lock's boundary: R6's ceiling, and S4, for
   a lock cell's links to the outside. Then test whether unlocked
   light can cool a hot fiber by inelastic scattering through those
   links. This is the next physics question. Without it, the line
   never appears from a hot start.
2. **A kernel step.** Implementing the S-law in a kernel needs
   explicit authorization, and it is a law-table change under the
   ratchet. It must settle L5(c) items 1–3 (pair detection, the
   cell door, space flow) before the battery can be run against it.
3. **`g`.** Either rule its value, or accept that the lock's mass
   scale is a free unit, with its ratios set at leading order.
4. **Size.** Decide whether "the smallest even cycle" counts as a
   size law. `v99/reopen/OPTIONS.md` E2 and this page say it does
   not.

## Reproduce

From `v102/lock/cas/`:

| order | command | writes |
|---|---|---|
| 1 | `python3 -I l1_spectrum.py` | `out/l1_spectrum.out` |
| 2 | `python3 -I l1_exact8.py` | `out/l1_exact8.out` |
| 3 | `python3 -I l2_hold.py` | `out/l2_hold.out` (L2 and L6), about 40 min |
| 4 | `python3 -I l3_walk.py` | `out/l3_walk.out`, about 10 min |
| 5 | `python3 -I l4_mass.py {slaw,heis} {24,32,48} {0,1} 64`, then `python3 -I l4_collect.py` | `out/l4_*.out`, `out/l4_collect.out` |
| 6 | `python3 -I l4b_line.py` | `out/l4b_line.out` |
| 7 | `python3 -I l5_costs.py` | `out/l5_costs.out` |
| 8 | `python3 -I l7_size.py` | `out/l7_size.out` |

Instruments:

- `slaw.py`: the S-law ring Hamiltonian, with the local-constraint
  basis.
- `clusters.py`: the S-law des Cloizeaux clusters.
- `mpo3.py`: the compact exact MPO.
- `../../phase1/cas/dmrg.py`: the DMRG validated in P2.
