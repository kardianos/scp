# v102 lock — pre-registration (S1–S4 ruled)

**Written:** 2026-10-07, after `../law/DECIDE.md` and before any
result in this directory. The only numbers on this page are table
constants and the algebra that follows from them.
**Authority:** user, this session: "Sure, follow those rabbit trail
S1-S4 and try to make a lock work."

---

## 1. Ruling R7

The four sentences of `../law/DECIDE.md` are ruled for v102, as
written there.

| | sentence | form taken here |
|---|---|---|
| S1 | A link's space mode is a discrete oscillator with its own pitch `Ω`. | one harmonic mode `c_l` per link; no self-load |
| S2 | `Ω` exceeds the pair it converts. | `Ω = w2 = 2.9`, the bare dense rung |
| S3 | Conversion on a link is a coherent pair amplitude `g`, with `g/Δ` small. | `g Σ_l (c_l† A_l + h.c.)`, the G-wedge operator |
| S4 | Transport inside a locked fiber is slow, or zero. | **zero** (v98 §1 read strongly) |

**Choice inside S4.** `../law/DECIDE.md` offered t = 0 or a fiber
transport constant ≤ 0.08–0.12. t = 0 is taken: it adds no constant,
and it is v98 §1's own sentence ("locked field cannot travel as
itself"). The margin for t > 0 is reported (L5b), not used.

**Scope.** Like R1–R6, R7 changes what the theory may assume. It does
not change the code.

- `laws_V3a` in the tree is not edited.
- The v97 kernel is not edited.
- No battery is run, because no kernel implements S1–S3.

The ratchet debt this creates is listed in L5c. It is owed before any
kernel carries the S-law.

## 2. The S-law (what is simulated)

A fiber is a ring of `n` cells, with links `l = (i, i+1)`.

| part | statement | source |
|---|---|---|
| cell field | two boson modes `b_i1, b_i2`; occupation `N_i ≤ d_i` (2 on a ring) | R1, R2, R6 |
| field energy | `E_f(N)`, the table's load law: `E_f = N w1/(1 + q ħ E_f/cap)`, self-consistent, `ħ = A0/2π` | reading R-b, `../law/cas/law.py` |
| link space | one harmonic mode `c_l` at pitch `Ω = 2.9` | S1, S2 |
| conversion | `g Σ_l (c_l† A_l + h.c.)`, `A_l = b_i1 b_j2 − b_i2 b_j1` | S3, G-wedge |
| transport | none | S4 |

    H = Σ_i E_f(N_i) + Ω Σ_l c_l† c_l + g Σ_l (c_l† A_l + h.c.)

**Table constants.** `w1 = 1.65`, `q = 1.2`, `cap = 2.5`,
`A0 = 1.15`. So `E_f(1) = 1.46217` and `E_f(2) = 2.67252` (rate
units).

### Algebra fixed before any run

1. **R6 becomes a local conservation law.** With t = 0, every cell
   conserves `N_i + c_{i−1,i} + c_{i,i+1} = 2`. Each link end holds
   one quantum, either as field in its cell or as half of its link's
   space quantum.
2. **Every converted configuration is uphill.** A cell with deficit
   `d_i = 2 − N_i` contributes `Ω d_i/2 − [E_f(2) − E_f(N_i)]`. That is
   0.2397 for d = 1 and 0.2275 for d = 2, both positive.
   - The first conversion on a ring costs **Δ = 0.4793**.
   - Each further adjacent conversion costs 0.2275. Phase 1's linear
     model charged 0.479 for each.

   So the S-law is not Phase 1's model. Multi-conversion is cheaper,
   and that enters at fourth order in `g`.
3. **Ends.** A degree-1 end cell holds one quantum (R6).
   - An end link converts from cells holding (1, 2):
     Δ_end = 0.2275.
   - A two-cell fiber converts from (1, 1): **Δ₁₁ = −0.0243 < 0**. The
     converted state is lower, so a two-cell fiber is not a lock under
     the S-law.

   Every test below uses rings. Open chains of at least 3 cells appear
   only inside the effective-chain instrument (L4).
4. **g.** S3 gives no value. The registered point is g/Δ = 0.1
   (g = 0.04793), Phase 1's. g/Δ = 0.05 and 0.2 are reported
   companions. No result is placed by choosing `g`.

## 3. Gates

All energies are in rate units. `J` means a ring's fitted exchange
(L1 procedure) unless stated otherwise. `λ_max` is the largest
eigenvalue of `A†A` on two full cells; the instrument computes and
prints it.

### L1. The S-law ring's spectrum is the Haldane ring's

- **Setup.** Microscopic rings at g/Δ = 0.1, t = 0, sector S^z = 0:
  - n = 6, exact (no truncation of space quanta);
  - n = 8, with Σc ≤ 3. The Σc ≤ 2 result is also reported, to show
    convergence.
- **Measure.** The first 12 excitation energies, with multiplicity,
  against the spin-1 Heisenberg ring of the same `n`, with one fitted
  `J`. This is the R-P1′ procedure.
- **Pass:** every miss below 5%, and a ground state with total spin 0.
- **Reported:** `J_fit` against the two-cell cluster's `J`, and against
  `2g²/Δ`.

### L2. The lock holds its field (real time)

**(a) S-law.** The initial states are unlocked: field at the ceiling,
no space, and no pairing.

- the Néel product `|+1, −1, …⟩`;
- three seeded random Jones-vector product states: each cell holds
  `(α b†₁ + β b†₂)² |0⟩/√2` with random normalized `(α, β)`.

Each is evolved on n = 6 (exact, all S^z) and n = 8 (Σc ≤ 3, all S^z)
to T = 100/J. Pass if all three hold:

- relative energy drift below 1e-8;
- `⟨c_l† c_l⟩ ≤ 5 (g/Δ)² λ_max` at every link and every sampled time
  (four times the sudden-start Rabi amplitude, plus margin);
- the mean field fraction over the second half of the run is at least
  `1 − 2.5 (g/Δ)² λ_max`.

**(b) Control R-c: space without pitch (S1 refused).**

- **Setup.** An n = 4 ring, Σc ≤ 1. Each link's single mode is
  replaced by M = 400 modes, spread uniformly over ±1.0 around the
  resonant pair energy, with couplings `g/√M`. That is a flat
  continuum carrying the same total weight `g²`.
- **Prediction.** The golden-rule rate
  `Γ = Σ_l (π g²/W) ⟨A_l† A_l⟩`, evaluated on the initial Néel state.
- **Pass:** the probability of one real conversion exceeds 0.9 by
  t = 5/Γ.
- **Same-ring check:** the S1 version of this n = 4 ring stays under
  (a)'s bound over the same time.

**L2 passes if (a) and (b) both pass.** That shows the hold, and that
S1 is what makes it.

### L3. The lock's quantum walks by conversion

- **Setup.** Microscopic ring n = 10, Σc ≤ 2 (amended below to
  Σc ≤ 3), at g/Δ = 0.05 and 0.1.
  - Ground state |0⟩ in S^z = 0.
  - Kick: |φ⟩ = S⁺₀|0⟩/‖·‖.
  - Spin profile: `m_r(t) = ⟨φ(t)| S^z_r |φ(t)⟩`. It sums to 1
    exactly.
- **Comparison.** The same kick on the spin-1 Heisenberg ring n = 10,
  at τ = J t, with `J` fitted on this n = 10 ring by the L1 procedure.
- **Pass, if all three hold:**
  - `sup_{r, τ ∈ [0, 4]} |m_r^S(τ/J) − m_r^H(τ)| < 0.05`, at both g/Δ.
    This is collapse in units of 1/J, so the walking rate is set by
    the conversion exchange `J ∝ g²`.
  - The field never travels: no cell's `⟨N_i(t)⟩` leaves `2 ± 2·max⟨c⟩`.
    This is exact by algebra item 1; it is printed as a check.
  - At g = 0 nothing moves. This is exact; it is printed.

### L4. The lock's quantum is massive and single

**(a) A mass exists.**

- **Setup.** The effective chain from S-law clusters: the
  des Cloizeaux h2 and r3 from 2-cell and 3-cell clusters (exact to
  fourth order), with the loaded `E_f` and the local constraint. Ends
  are degree-1 cells (cap 1), as in P2.
- **Measure.** DMRG gap `E(S^z = 1) − E(S^z = 0)` at L_bulk = 24, 32,
  48, extrapolated as `D + a/(L+1)²`. The Heisenberg chain at the same
  L is run as a reference, and its extrapolation is reported beside it.
- **Pass:** `Δ_∞/J_h2 > 0.2`.

**(b) One line.** On the n = 8 micro ring (Σc ≤ 3):

- the lowest S^z = 2 state, and the lowest S = 0 excitation that is
  not a triplet member, both lie at ≥ 0.98 × 2Δ₈;
- the first level with ⟨Σc⟩ > 0.5 lies above 2Δ₈.

**(c) Reported, not gated.** `Δ_∞/J` against 0.4105. Phase 1 P3
already found drift at order (g/Δ)², so this is not re-gated.

**L4 passes if (a) and (b) both pass.**

## 4. Reports (not gates)

### L5. How the lock could still fail in the world (the costs of S1–S4)

**(a) Space transport.** The table moves space (`s_k = 0.06`,
pressure-driven). That would broaden S1's discrete oscillator into a
band.

- **Setup.** One bond (two full cells) whose link oscillator couples
  to a semi-infinite lead of space modes at pitch Ω, with hop `t_s`.
  Single-conversion sector, singlet start.
- **Reported:**
  - the long-time field survival against `t_s/Δ ∈ {0.1, 0.2, 0.3, 0.4, 0.6}`;
  - the predicted condition: the lock holds iff the band edge stays
    above the pair, `2t_s < Δ` for a 1-D lead and
    `2√(z_l − 1) t_s < Δ` for a z_l-coordinated lead;
  - what the reading `t_s = s_k` implies. The foam degree is 7.28
    (`v89/MASS.md`), so the link graph has z_l ≈ 2(7.28 − 1) = 12.6.

  `s_k` is a diffusive conductance, not an amplitude, so this is a
  reading.
- **Expectation, not a gate:** marginal.

**(b) The t > 0 margin under the S-law.** On the n = 6 ring (Σc ≤ 2,
general basis, loaded `E_f`), find the transport at which the ground
state stops being the S = 0 lock. Compare it with Phase 1's linear-Δ
value.

**(c) Regressions owed under the ratchet.** These are listed, not
computed.

- **Detection.** A lone quantum cannot convert (P6), so a click needs
  a partner.
- **Two conversion laws.** The cell door's rate conversion
  (`Ee > e_cond`) must not act on locked field, or it converts the
  lock regardless of S1.
- **Space flow.** S1's discreteness has to be reconciled with the
  `s_k` space flow that carries g1, g3 and g4.
- **Context.** S4 makes transport depend on whether field is locked.

### L6. Formation

Can an unlocked fiber become a cold lock without an imported cooler?

- **Measured.** From L2's runs, the internal energy above the ground
  state in units of the gap, and the canonical temperature that
  matches it on the n = 6 ring.
- **Argued.** Which channels in the S-law could carry the excess away.
- **Expectation:** none.

### L7. Size

E₀(n) for n = 3 … 10 (n ≥ 9 with Σc ≤ 2; amended below to Σc ≤ 3).
Reported:

- E₀(n)/n;
- the second difference;
- the fission energy `E₀(a) + E₀(b) − E₀(a+b)`.

**Expectation:** no ring size is selected.

## 5. Endings (exactly one is claimed in `DECIDE.md`)

| ending | condition | consequence |
|---|---|---|
| A | L1–L4 pass | the lock works under the S-law: it holds its field, it walks by conversion, and its quantum is massive and single. `g` remains the one unset constant, and the L5 costs are owed |
| B | L2 passes; L1, L3 or L4 misses | the field is held, but its quantum is not the derived line |
| C | L2 misses | the S-law does not hold a lock |

## Amendment (2026-10-07, after L1 and before L2–L7 ran)

L1's convergence check found that Σc ≤ 2 is not converged.

| n = 8 | J_fit | vs no truncation |
|---|---:|---:|
| no truncation | 0.009167 | — |
| Σc ≤ 3 | 0.009086 | −0.9% |
| Σc ≤ 2 | 0.008507 | −7.2% |

The cause is the S-law's cheap conversion strings: each further
adjacent conversion costs only 0.2275, against `g√λ_max ≈ 0.12`.

So L3 and L7 run with **Σc ≤ 3**, not Σc ≤ 2. The Σc ≤ 2 numbers are
not used. No bar changes. The n = 8 run without truncation also
passes L1: 3.3% miss, singlet (`cas/out/l1_exact8.out`).

## 6. What does not count

- Choosing `g`, or a reading, because a gate passes.
- Tuning `w1`, `w2`, `q`, `cap` or `A0`.
- Reading ξ, or any length, as a cell count.
- Softening a bar after a result.
- Using L5–L7 to rescue a gate.
