# v103 gate g1 — pre-registration (can a medium with no container hold three dimensions?)

**Written:** 2026-10-11, before any run of this step. **Nothing has
been run:** no gate, no shakedown, no registered run.

**Authority:** the user, 2026-10-11:
- "Was this analysis done in 3d? Some of these properties will only
  occure in full dimensional analysis."
- After the 3-D gate was proposed: "Set it up, register it, note the
  current status of the run so any additional agent will find the
  status, commit (no coauthor) and push. Do not yet run."

The gate was named before: c2 DECIDE ("The 3-D geometry gate (R14)
before any 1/r test"), and in the README's Next list ("It gates every
3-D step, including a 3-D rerun of lane a").

**A measurement in models; it rules nothing.**

## 0. Frame check (R14, `../../FRAME.md`)

| question | answer for g1 |
|---|---|
| **What are the places?** | **The quanta themselves.** The medium's whole state is which quanta touch which. There is no lattice of places: the tool stores contact lists and nothing else. |
| **Can a place exist with no medium in it?** | **No.** A place is a quantum. When a quantum is removed, no slot is left behind: its six contacts are re-paired among its neighbours at once. An unpaired contact would be an edge inside the medium, and there is no edge. So the medium closes over the gap. |
| **How is distance counted?** | **Through the medium:** the number of contacts on the shortest path between two quanta. Volume is the number of quanta. |
| **What is a gap, and how is it read?** | **A missing quantum.** It is read from the medium's own shape. Where the medium closed over it, its contacts fall out of the flat pattern (squares, triangles, links; §5). No label (R13): the tool never records where the quantum was removed. An **excess** is one quantum more, the opposite sign. |
| **Container caveat** | **None in g1's model.** The flat start is built from bookkeeping coordinates that are dropped at once; no move and no measurement uses a coordinate. Earlier results: c1, a1 and b1 used containers (caveated). c2 was frame-faithful only in one dimension. |

## 1. Why a gate, and what three dimensions add

**a1, b1 and c2 ran on one-dimensional rings** (a1 DECIDE, dimension
note). Spin as rotation, exchange, the sharpness of a quasiparticle, the
threshold for binding and any 1/r law are dimension-dependent. **A 3-D
step cannot use a fixed grid:** that is the container R14 rules out.
The medium must hold its own three dimensions.

**One thing only more than one dimension can show: a gap in a medium
of bare quanta.**
- In one dimension, a ring of N − 1 quanta is just a ring. Nothing
  marks where the missing one was. c2 had to give each piece content
  (n = 2) for a gap to be readable.
- In three dimensions, the way the medium closes over a missing
  quantum can leave a mark in the medium's own shape.

**What known physics says** (structure and method only, R8):
- **Free local re-forming makes expanders.** Random local flips turn a
  connected regular graph into an expander: a medium with no dimension
  at all. (Allen-Zhu, Bhaskara, Lattanzi, Mirrokni and Orecchia,
  "Expanders via local edge flips", SODA 2016, arXiv:1510.07768.)
- **Euclidean dynamical triangulations in three dimensions** show only
  a crumpled phase and a branched-polymer phase, with a first-order
  transition between them. Neither is an extended 3-D geometry
  (hep-lat/9712011, and later work).
- **A crystal holds its geometry, but its forces confine.** Kleinert's
  "world crystal" (1987) gets curvature from defects. Kleinert and
  Zaanen (2004, gr-qc/0307033): in a crystal with ordinary
  (first-gradient) elasticity, the force between curvature sources
  grows with distance, which confines. A 1/R force needs dislocations
  to proliferate, giving a nematic phase that keeps orientation but not
  position. **So holding geometry and carrying a 1/R law are different
  properties, possibly of different phases.** g1 measures the phase of
  any medium that holds.
- **Fracton–elasticity duality** (Pretko and Radzihovsky, PRL 120,
  195301, 2018; arXiv:1711.11044). The elasticity of a crystal is dual
  to the fracton tensor gauge theory. Its glide constraint comes from
  conservation of vacancies and interstitials. c2's tweaked links have
  exactly this structure, so a medium that holds its geometry is where
  they live.

## 2. The medium

- **N quanta and their contacts, nothing else.**
  - Every quantum has six contacts.
  - Every contact is paired: there is no unpaired contact, because
    that would be an edge inside the medium.
  - Two quanta touch at most once, and no quantum touches itself.
- **Why six.** Six is the contact count at which a 3-D packing of
  frictionless grains is just rigid (Maxwell's count, 2d). The medium
  is told six. It is not told how to arrange them.
  - Six alone does not build in three dimensions. With no candidate
    rule the medium becomes an expander (§1).
  - The 2-D triangular arrangement also has six contacts (§3).
- **Distance and volume.**
  - Distance is contacts along the shortest path.
  - Volume is the number of quanta.
- **Re-forming: local flips.**
  - Along a path i–j–k–l, the contacts ij and kl become ik and jl. A
    flip is refused if it would double a contact or close i onto l.
  - A flip is local: every quantum it touches is within one contact of
    the path.
  - It keeps six contacts per quantum.
  - **It keeps the medium in one piece.** After the flip, i reaches j
    through k, and k reaches l through j. So every path that used ij or
    kl can be rerouted.
- **Proposal.** Pick a directed contact j→k uniformly, then i uniformly
  among j's other five contacts, and l uniformly among k's other five.
  - Each flip arises in two ways, each with probability
    (1/2B)·(1/5)², where B = 3N contacts. Its reverse arises the same
    way, so the proposal is symmetric.
  - The acceptance is Metropolis at temperature t, so the stationary
    measure is the Gibbs weight e^{−E/t}.
- **The random process only.** The geometry is not made coherent
  (quantum) here, because the space of contact patterns is out of
  reach. Lanes a, b and c will put coherent gaps on whatever medium
  passes this gate.

## 3. The candidates (three, all eight combinations)

Energy is in units of t, at 2t per unit of deviation squared. That is
c2's matched scale: a single unit deviation in any one count costs 2t.

| | candidate | on | off |
|---|---|---|---|
| **a** | **flat squares** (position) | 2t (s_e − 4)² for every contact e, where s_e is the number of squares (4-cycles) through e. Four squares per contact is the discrete zero curvature of a medium of cube-like pieces: four pieces around every edge | no cost |
| **b** | **octahedral link** (orientation) | at every quantum i, for each contact j: 2t (q_ij − 4)², where q_ij is the number of i's other contacts that share a square with i and j. With q = 4 everywhere, each contact has exactly one opposite, and the six contacts form the octahedron | no cost |
| **c** | **no triangles** | 2t for each triangle (three mutually touching quanta) | no cost |

**Exact zeros, checked by gate T0 before any run:**
- **The cubic pattern** on a closed 3-torus (L ≥ 5) has E = 0 in every
  set. Every contact has four squares, every link is octahedral, and
  there are no triangles.
- **The 2-D triangular arrangement** also has E_a = E_b = 0: each
  contact lies in four rhombi, and each link is the octahedron. But it
  has 2 triangles per quantum. **So among the zeros of a and b, only c
  separates two dimensions from three.**

**A structural fact about one gap, also checked by T0.**
- A perfect flat medium of N − 1 quanta exists: a sheared closure, the
  circulant C_n(1, a, b) for suitable a and b, with zero energy in
  every set.
- So in a closed medium, one gap could in principle be absorbed into
  the medium's global shape: the whole medium simply has one quantum
  less.
- **A lone gap is a localized thing only while the medium holds its
  shape against global re-forming.** That is the frame's own statement:
  there is no elsewhere, and the medium closes over its gaps.

## 4. Sectors and schedule

**Sizes** L = 6, 8, 10, 12, 16, so N = L³ = 216, 512, 1000, 1728,
4096. The bars use L = 12 and L = 16. The ladder is reported for its
trend; it is not a scan.

**Units.**
- A sweep is B = 3N proposals.
- A measured window is 2000 sweeps at t, with a sample every 100
  sweeps (20 samples).

**Sectors, in this order, in one run per (set, L):**

| sector | start | schedule |
|---|---|---|
| **M**, the medium | the flat start: the cubic pattern of contacts, closed | 2000 sweeps at t, then the measured window |
| **G**, one gap | M's state after its 2000 sweeps. Remove one quantum chosen uniformly. Pair its six contacts among its six neighbours by the Gibbs weight over every valid pairing (at most 15; a pairing that would double a contact is invalid) | 2000 sweeps at t, then the measured window |
| **X**, one excess | the same state. Choose a quantum p uniformly. The new quantum takes three pairwise-disjoint contacts u–v, both ends within distance 2 of p, splitting each into u–q–v. It thereby has six contacts, and nobody else's count changes. The triple is chosen by Gibbs weight over every valid triple | 2000 sweeps at t, then the measured window |
| **A**, assembly | the flat start, melted by 2000 sweeps with every rule off (infinite temperature) and read once (A0) | cooled through T = 64, 32, 16, 8, 4, 2, 1 × t, 300 sweeps each, then the measured window at t |

- Removal and insertion are exact reverses: six contacts become three,
  and three become six.
- **Seed.** The registered seed is 1000003·L + 7·(4a + 2b + c) + 1,
  expanded by splitmix64 into xoshiro256**.
- **Rigidity.** Before M, every flip of the flat medium is enumerated
  and its exact energy cost recorded: the cheapest flip, how many reach
  it, and how many cost ≤ 0.

## 5. Measurements (every sample)

**Order.**
- ψ_T: the share of contacts with s_e = 4 (position).
- ψ_O: the share of quanta whose six links all have q = 4
  (orientation).
- Triangles per quantum.
- E/N.

**Geometry.**
- ⟨d⟩: the mean distance over pairs of distinct quanta.
- The diameter.
- The number of pieces.
- The mean sphere count S(r).
- λ2: the lowest nonzero eigenvalue of the contact Laplacian, by
  Lanczos with a convergence flag.

**Mismatch.**
- A contact is mismatched if s_e ≠ 4, or it is in a triangle, or q ≠ 4
  at either end.
- D is the set of quanta touching a mismatched contact.
- Clusters: two mismatched contacts are joined when their ends are
  within distance 2.
- Spread: the largest distance between two members of D.

**Footprint (G and X).** The centre is the set of quanta whose largest
distance to D is smallest (read from the state, no label). By distance
r from the centre:
- the number of quanta;
- the square deficit Σ(4 − s_e);
- the mismatched contacts;
- the mean of m(x), the mean distance from x to every other quantum;
- the medium's Green's function from a unit source at the centre
  (zero mean, by CG).

**Reference (M).** The same profiles from one quantum.

**Derived by `cas/tab.go`.**
- **Two-size readings, 12 → 16:**
  - d_H = ln(N16/N12) / ln(⟨d⟩16/⟨d⟩12);
  - d_s = −2 ln(N16/N12) / ln(λ2_16/λ2_12).
  - The flat torus reads d_H = 3.00 and d_s = 3.05; T0 prints both.
- **Ladder fits** over all five sizes (reported).
- **The gap's signed "where":**
  - δ⟨d⟩ = ⟨d⟩_G − ⟨d⟩_M. Negative is less where.
  - Δm(r) = m̄(r) − m_ref.
  - Its decay exponent p, from Δm ∝ r^−p over 2 ≤ r ≤ L/2.
- **Field footprint:** ΔG(r) = G_gap(r) − G_ref(r).
- **Defect radius r_def:** the largest r with at least half a
  mismatched contact on average.
- **The same readings for the excess.**

## 6. Gates (T0, before any registered run), in `cas/out/t0.out`

1. **Invariants.** Every set, at L = 6, flat and melted starts, 10⁵
   proposals:
   - six distinct contacts per quantum;
   - no self contact;
   - the adjacency bits match the contact lists;
   - one piece.
2. **Local equals full.** For every proposed flip in every set (4000
   per start), the local energy change equals the full recomputation
   exactly (integers).
3. **Detailed balance.**
   - **N = 8:** the 105 states, K8 minus a perfect matching, are all
     alike. 10⁷ proposals must give the uniform measure, with total
     variation (TV) < 0.01.
   - **N = 9:** four kinds of state, whose complements are C9, C6+C3,
     C5+C4 and C3+C3+C3, with exact counts 20160, 5040, 4536 and 280.
     In every set, at a test temperature of a quarter of the energy
     span, the kind frequencies must match count × e^{−E/T}, with
     TV < 0.01.
4. **The flat medium**, at L = 6, 8, 10, 12, 16:
   - E = 0;
   - ⟨d⟩ = (3L/4)·N/(N − 1) exactly;
   - diameter 3L/2;
   - λ2 = 2 − 2cos(2π/L) to 1e-8;
   - and its two-size d_H and d_s are printed.
5. **The 2-D competitor and the instruments.**
   - The triangular arrangement at L = 6, 8 has E_a = E_b = 0 and
     2L² triangles.
   - On triangular media of N = 1764 and 4096, the instruments must
     read d_H and d_s in [1.8, 2.2].
   - On the random six-contact medium (the flat start after 2000 free
     sweeps, N = 1728 and 4096), they must read d_H > 4.
6. **Making a gap and an excess**, on the flat medium at L = 6, in
   every set:
   - 15 valid pairings;
   - every candidate's local energy change equals the full one;
   - the medium is intact afterwards, with N − 1 or N + 1 quanta and
     six contacts each.
7. **Lanczos against dense** λ2 at N = 216 (flat and melted) and 512
   (melted), to 1e-8.
8. **The structural fact (§3):** a zero-energy circulant of n = L³ − 1
   quanta is found for L = 6 and 8.

## 7. Bars and the ending

Every bar must hold at both L = 12 and L = 16.

**H, the medium holds three dimensions.** On sector M:
- H1: d_H in [2.7, 3.3];
- H2: d_s in [2.7, 3.3];
- H3: one piece and a converged λ2 in every sample;
- H4: steady: ψ_T and ψ_O each move by at most 0.01 between the first
  and second half of the window.

**P, a readable, persistent gap.** On sector G, for a set that holds:
- P1: the medium M has on average fewer than one mismatched contact
  per sample. The medium is quiet enough to read a gap in.
- P2: in every G sample, the mismatched contacts form exactly one
  cluster, with spread ≤ 6.

**Reported with the ending, not scored:**
- **Phase** (M, L = 16):
  - crystal if ψ_T ≥ 0.99 and ψ_O ≥ 0.99;
  - nematic if ψ_T < 0.99 and ψ_O ≥ 0.99;
  - disordered otherwise.
- **Footprint class:** contact if r_def ≤ 3 at both sizes and does not
  grow by more than one from 12 to 16; extended otherwise.
- **Assembly:** H1–H4 on sector A.
- **The excess:** P1–P2 on sector X. The signs of δ⟨d⟩ and of the
  square deficit, for G against X.

| ending | condition |
|---|---|
| **G**, the gate passes | some set passes H and P. That set, with its phase and footprint, is the 3-D substrate for later lanes |
| **H**, holds without a gap | some set passes H, but no set passing H passes P |
| **F**, fails | no set passes H |

**What I expect** (written before running):
- **none and c melt.** Flips with no cost make an expander: ⟨d⟩ grows
  like ln N, and λ2 stays finite. So d_H and d_s read far above 3, and
  H fails. That is "free re-forming loses three dimensions", measured.
- **Every set with a or b is rigid.** Every flip of the flat medium
  costs tens of t (my estimate is about 110t for a straight flip under
  a), so at temperature t nothing is accepted. The medium holds by
  rigidity: H passes, and the phase is crystal.
- **Nothing assembles** from the melt within the schedule (a glass).
  Without c, I expect patches of triangles there (the 2-D zero).
- **The gap is readable and persistent, and frozen** where it was
  made. Its defect footprint is contact (r_def ≤ 2 or 3). It is less
  where: δ⟨d⟩ < 0.
- **The excess is more where** (δ⟨d⟩ > 0). Both break squares, so the
  square deficit has the same sign for both. The shape's response is
  signed in distance and unsigned in squares.
- **The distance footprint Δm(r) may not decay.** Counted in contacts,
  a crystal's metric is taxicab-like, and its geodesics are massively
  degenerate. One shortcut where the medium closed over a gap then
  shortens a finite share of all paths, so I expect p ≈ 0 out to the
  torus scale. The field footprint ΔG(r) has no such degeneracy.
- **Expected ending: G (crystal, contact) for every set with a or b.**
  That reproduces the open problem as measured: a rigid medium holds
  three dimensions and its gaps, and closes over them locally. By
  Kleinert–Zaanen, such a crystal is expected to confine rather than
  carry 1/R. That points to g2 (Then).

## 8. Disclosed before running

- **Nothing has been run.** The tool and the table code compile (C with
  `-Wall -Wextra`, and `go vet`). No gate, no shakedown and no
  registered run has executed.
- **Constants are not scanned.**
  - t = 1, and 2t per unit deviation squared, are c2's scale.
  - Six contacts is a standing choice of this step (§2).
  - The sweeps, the cooling schedule, the sizes, the test temperature
    of T0.3 and the seeds are all registered here.
- **The flat start is the pattern being tested,** so "holds" means
  stability, not self-assembly. Assembly is reported separately (A).
- **Frozen dynamics are expected.** For the rigid sets, "holds" will
  mean that no re-forming is ever accepted. The rigidity number (the
  cheapest flip) is reported so that this is visible.
- **"Six contacts" is enforced exactly.** A gap or an excess changes
  N, never anyone's contact count.
- **The flips** are local and keep the medium in one piece, so the
  "one piece" bar (H3) can fail only through a tool defect.
- **Footprints are read around a centre computed from the state.** The
  removed quantum's identity is never used.

## 9. Reproduce

From `cas/`:
1. `./run_t0.sh`. It must end with `T0 PASS`.
2. `./run_g1.sh` runs 40 processes, 12 at a time, and touches
   `out/run.done`.
3. `go run tab.go > out/tab.out`.
4. Write `../DECIDE.md`.

## Then (after DECIDE)

- **If G:** the passing medium is the 3-D substrate for:
  - a 3-D rerun of lane a (spin, exchange, sharpness, binding);
  - c3: two gaps in it, the 1/r test;
  - b2.
- **If the passing medium is a crystal with a contact footprint:**
  - g2, a medium that keeps orientation but lets dislocations move
    (the Kleinert–Zaanen nematic);
  - the closest candidate here is b without a.
- **If H or F:** the 3-D lanes wait. The design moves to a medium with
  a soft metric on a held topology (quanta in the contacts).
