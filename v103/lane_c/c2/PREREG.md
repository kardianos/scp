# v103 lane c, step c2 — pre-registration (two things: a sign, and something without one)

**Written:** 2026-10-11, before any run of this step.
**Authority:** the user, 2026-10-11:
- "I think there are two things going on, a signed electric charge, and
  a void-void interaction which can only be expressed through 1/r^2,
  and is thus without a sign. The math would work if the links are
  tweaked, still shapes the field without the sign."
- On the candidates: "your proposals sound reasonable for initial
  attempts".
- "Please run through the proposed candidates combinations (all 16), be
  sure to write down the results as always."

**A measurement in models; it rules nothing.**

## 0. Frame check (R14, `../../FRAME.md`)

| question | answer for c2 |
|---|---|
| **What are the places?** | Cells of a ring. Each cell is a piece of the medium holding n_i quanta (full: n̄ = 2). |
| **Can a place exist with no medium in it?** | **No.** Every cell always holds at least one quantum (n_i ≥ 1). A cell can be depleted, never emptied. R14 forces this rule, and it applies in all 16 sets. One consequence: at most one gap fits in a cell, which is exclusion over space. |
| **How is distance counted?** | **In quanta.** The distance between two cells is the number of quanta strictly between them, along whichever way round is shorter. A gap shortens the medium: less where. |
| **What is a gap, and how is it read?** | A missing quantum, read from the cell's own count (n_i = 1). No label (R13). An **excess**, n_i = 3, is the opposite sign (the antiparticle). |
| **Container caveat** | In one dimension, the moves of the cell picture and of a sequence of quanta with gaps are the same. What differs is the distance, and here it is counted in quanta. c1, a1 and b1 counted it in fixed sites. |

**What one dimension can and cannot decide.**
- **It can decide:**
  - whether like gaps attract;
  - whether like excesses attract too, so the attraction does not care
    about sign;
  - which of the four candidates carries the attraction;
  - whether an isolated gap can move at all (Mach).
- **It cannot decide the 1/r law.** That needs three dimensions and a
  medium that holds its own geometry, which is a later step.

## 1. The medium

- **A periodic ring** of L cells (no edge), n̄ = 2, t = 1.
- **Sectors:** N = 2L (the medium), 2L ∓ 1 (one gap or one excess), and
  2L ∓ 2 (two gaps or two excesses).

## 2. The candidates (four, all 16 combinations)

| | candidate | on | off |
|---|---|---|---|
| **a** | **the lack recurses** (chain constraint) | cap 2 = n̄: no cell can hold an excess, so a gap can never be absorbed, only passed along | cap 3: a gap can meet an excess and both vanish |
| **b** | **tweaked links** (rank-2: the re-forming conserves the medium's first moment P = Σ_i i·n_i mod L) | moves T3 (n_{i−1}, n_i, n_{i+1}) ↔ (+1, −2, +1) and T4 (n_i … n_{i+3}) ↔ (+1, −1, −1, +1), each with amplitude t times the bosonic factors | ordinary re-forming: one quantum between neighbouring cells, amplitude t |
| **c** | **inertia** (coherent re-forming) | the state is the quantum ground state of H = K_b + E_d | the state is the stationary measure of the random process with the same moves, averaged over starts: the Gibbs weight e^{−E_d/t} over the whole sector |
| **d** | **flat level** | E_d = t Σ_i (n_i − n_{i+1})²: only differences between neighbours cost | stiff: E_d = 2t Σ_i (n_i − 2)². Leaving the level costs |

**Choices registered here:**
- **Matched scales for d.** U = 4t, the natural unit from lane b. A
  single deviation δ in one cell costs 2tδ² in both forms, so they
  differ only in what they charge for long, smooth changes.
- **The moves for b** are every first-moment-conserving re-forming of
  range up to three bonds. Under a on, T3 can never act (it would empty
  a cell or overfill one), so T4 carries b there.
- **c off, structurally:** averaged over starts, the stationary measure
  of any symmetric random process is the Gibbs measure on the whole
  sector, whatever the moves. **So b cannot change any c-off static
  measurement.** The four c-off pairs that differ only in b must agree
  exactly, which is a built-in check. b still changes which states can
  reach which (§3, item 4).

## 3. Measurements

The tool is `cas/gap1d.c` (C), with sizes L = 12 and L = 16. The tables
come from `cas/tab.go` (Go).

1. **States.**
   - **c on:** the lowest-energy states. If the lowest level is
     degenerate (as with fragmentation under b), the state is the equal
     mixture over the whole degenerate lowest manifold.
     - Under b on, sectors are split by P. Translation by one cell shifts
       P by N, so P = 0 and P = 1 cover every sector up to translation
       for N = 2L ∓ 1 and 2L ∓ 2. N = 2L needs every P.
     - The lowest energy is the minimum over P. If P sectors tie, the
       state is the mixture of their lowest manifolds.
   - **c off:** the Gibbs weight, enumerated exactly.
2. **Binding.**
   - Like gaps: B_gg = E(2L−2) − 2E(2L−1) + E(2L).
   - Like excesses: B_ee = E(2L+2) − 2E(2L+1) + E(2L), under a off only.
   - E is the lowest energy (c on) or the free energy −t ln Z (c off).
     B < 0 means the pair is together.
3. **Pair shape**, against counting.
   - Each pair is conditioned on its minimal composition:
     - gg: exactly two gaps and no excess, in the 2L−2 state;
     - ee: exactly two excesses and no gap, in the 2L+2 state;
     - ge: exactly one gap and one excess, in the 2L state. This one
       is **reported, not scored**: in a quantum state a gap and an
       excess are created together, so their nearness measures creation
       as well as force.
   - **Reference:** the uniform measure over the same composition
     (counting alone).
   - **Readings:**
     - ρ = ⟨D⟩/⟨D⟩_ref, with D in quanta;
     - the profile V(D) = −ln[p(D)/p_ref(D)];
     - the weight of the minimal composition in the state.
4. **Mach (structure only).** These depend only on a and b. Each is
   read from the move graph of a placed probe configuration:
   - can an isolated gap reach the next cell?
   - can two adjacent gaps translate as a pair?
   - can two gaps L/2 apart come together?

## 4. Gates (T0, before any registered run), in `cas/out/t0.out`

1. **Hermiticity:** ⟨x|Hy⟩ = ⟨Hx|y⟩ to 1e-12, with random vectors, in
   all 16 sets.
2. **First-moment conservation:** under b on, every generated move
   conserves P mod L, counted exhaustively at L = 8.
3. **Exact limit.** Set a on, b off, d off, c on. The gaps are then
   hard-core particles hopping with amplitude 2t.
   - Two-gap sector: E = 4t − 4t·[cos(q1) + cos(q2)], with the
     antiperiodic two-fermion momenta q = ±π/L.
   - One-gap sector: E = 2t − 4t.
   - Both to 1e-10, at L = 8, 12, 16.
4. **Iterative against dense:** lowest energies and mixture
   observables, in every set at L = 8 and every sector, to 1e-9.
5. **Sector sizes** against the counting formula (trinomial
   coefficients for cap 3, binomial for cap 2).
6. **The c-off check:** the b-on and b-off pairs agree exactly (§2).

## 5. Bars and the ending

Every bar must hold at both L = 12 and L = 16.

- **Attraction of a like pair** (gg or ee): B ≤ −1e-6 and ρ ≤ 0.95.
- **Repulsion:** ρ ≥ 1.05.

| ending | condition |
|---|---|
| **U**, unsigned | some set has both like pairs, gaps with gaps and excesses with excesses, attracting. The attraction does not depend on which sign the pair carries (a off only) |
| **S**, sign-dependent | some set has one like pair attracting while the other does not, and no set reaches U |
| **A**, gaps only | gap pairs attract in some set, but only where excesses cannot exist (a on), so the sign cannot be tested |
| **N**, none | no set has an attracting like pair |

**The ledger** reports, for every set:
- B_gg, B_ee, ρ_gg, ρ_ee, ρ_ge (reported only), and the Mach row;
- the V(D) profiles, which are reported and not scored, because range
  in one dimension is not the 1/r test;
- for each candidate, whether it is present in every attracting set
  (necessary) or only strengthens the attraction.

**What I expect** (written before running):
- **d (flat)** makes like deviations cluster, because a boundary between
  levels costs and a cluster has fewer boundaries. So under c off, gap
  pairs and excess pairs attract at contact, and gap–excess pairs
  repel.
- **b with c** binds gaps by mobility, which is Pretko's mechanism. Under
  a on, an isolated gap is frozen. Two gaps can only breathe between
  contact and three cells apart (the pair cannot translate). They bind
  by B = −4t.
- **b off, d off** gives no attraction. With a on and c on, this is free
  fermions: a Pauli cost, ρ > 1.
- **Expected ending: U,** carried by d alone in the random process and
  by b with c in the coherent one. All of it is at contact; there is
  nothing at range in 1-D.

## 6. Disclosed before running

- **One dimension, two small rings.** Range is reported and not scored.
- **The constants are not scanned.** t = 1 and U = 4t (lane b's natural
  unit). The matched d forms follow from U. The temperature of the
  random process is t, the only energy unit.
- **Unlike pairs** (gap with excess) are reported, not scored (§3).

## 7. Reproduce

From `cas/`: build per the header of `gap1d.c`, then run `./run_t0.sh`
and `./run_c2.sh`, then `go run tab.go`.
