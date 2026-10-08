# v102 options — pre-registration (where the lock fails, option by option)

**Written:** 2026-10-08, after `../lock/DECIDE.md` and before any
result in this directory. The only numbers on this page are table
constants and the algebra that follows from them.
**Authority:** user, this session: "Continue to test the where the
locks may fail, as for the questions, on each one take a range of
options, both in isolation and in combination and test them; see
which ones make the most sense at the end of the day, don't just rely
on your intuistion or mine, see what the effect of each are."

**Scope.** This is a measurement campaign. It does not rule anything.

- No kernel is edited.
- No law table is edited.
- Battery runs with changed constants are run from scratch builds,
  through the battery's own `-extra` flag. They measure what a global
  form of an option would cost the table. They are not adoptions.

Which options "make sense" is decided by the scoring rule in §5,
which is fixed before any run.

---

## 1. The reference point (R7, as run in `../lock/`)

| axis | value at R7 |
|---|---|
| transport inside the lock | t = 0, by label (S4 read strongly) |
| boundary | none: the ring is alone |
| space flow on lock links | none |
| cell door | absent from the model |
| g | g/Δ = 0.1 |
| geometry | free rings |

The S-law, the table constants, and the loaded energies are those of
`../lock/PREREG.md` §2: `E_f(1) = 1.46217`, `E_f(2) = 2.67252`,
`Ω = 2.9`, `Δ = 0.47931`.

Further loaded energies that the options below need, computed with the
same load law:

| N | E_f(N) | E_f(N) − E_f(N−1) | Δ_N = Ω − 2[E_f(N) − E_f(N−1)] |
|---:|---:|---:|---:|
| 1 | 1.46217 | 1.46217 | −0.02435 |
| 2 | 2.67252 | 1.21035 | 0.47931 |
| 3 | 3.72861 | 1.05609 | 0.78781 |
| 4 | 4.67769 | 0.94908 | 1.00184 |
| 6 | 6.35356 | 0.80665 | 1.28670 |

Δ_N is the first conversion cost on a ring whose cells hold N quanta.
The boundary pair (a ring quantum with a lone outside quantum at
kinetic energy ε) is resonant with Ω at `ε_res = Ω − E_f(2) = 0.22748`.
That lies inside the outside band `[−2t, 2t]` for any `t > 0.114`.

---

## 2. The six questions and their options

Each option is tested **in isolation**: R7 with only that axis
changed. Selected **combinations** are then tested in §4.

### Q1. The boundary: R6's ceiling and the field's moves at a lock-to-outside link

A lock cell embedded in the world has links to cells outside the
ring. Three things must be said about such a link.

**Ceiling rule** (what R6's "one quantum per link end" counts):

| code | rule | needs a label? |
|---|---|---|
| C-lock | a ring cell's cap counts its ring links only (cap 2); cells outside the ring have no ceiling | yes: ring membership |
| C-geo | every cell's cap is its geometric degree | no |
| C-held | a cell's cap is the number of its links whose far cell holds field; a cell with no such link has no ceiling | no |

**Transport rule** (one rule for every link, inside and at the
boundary; `t` is the hop amplitude where a hop is allowed):

| code | rule | t | needs a label? |
|---|---|---:|---|
| T0 | no hop on any link touching a lock cell | 0.9 elsewhere | yes |
| T1 | uniform small hop, refused only into a cell at cap | 0.05Δ, 0.10Δ | no; a new constant |
| T2 | the table's hop, refused only into a cell at cap | 0.9 | no |
| T3 | refused into or out of a cell with `N = cap` | 0.9 | no |
| T4 | refused into or out of a cell with `N + Σ_l c_l = cap` (every link end held, by field or by half a space quantum) | 0.9 | no |

**Algebra fixed here.** On an isolated ring at cap 2, T4 refuses
every hop, because R6's local law `N_i + c_{i−1} + c_i = 2` holds at
all times. So T4 inside the ring is identical to T0, without a label.
It is run anyway, as a check.

**Boundary conversion:**

| code | rule |
|---|---|
| X0 | no conversion on a lock-to-outside link |
| X1 | S3 on that link too, with its own oscillator at Ω |

### Q2. Space flow on lock links (S1 against the table's `s_k`)

| code | rule |
|---|---|
| F0 | a lock link's space mode does not flow |
| F1 | coherent flow: the link oscillator hops to neighbouring space modes with amplitude `t_s`, read as `s_k = 0.06`, in three geometries: along the fiber (link to adjacent link), a tree lead with `z_l = 12.6`, and the lattice-edge estimate `z_l t_s` |
| F2 | rate-law flow: each lock link's space quantum leaves at rate `γ` with no pitch check (a Lindblad jump `c_l`), read as `γ = s_k` |
| F3 | pitch only where both ends of a link hold field; every other link's space is the pitchless continuum (R-c) |

### Q3. The cell door on lock cells (two conversion laws)

| code | rule |
|---|---|
| D0 | the table's door acts on lock cells: per step it converts `f_conv = 0.25` of the field above the work function (`wf_floor = 0.01` in matter) |
| D1 | lock cells are exempt (by label, or by T4's predicate `N + Σc = cap`) |
| D2 | the door is replaced everywhere by S3. A lone quantum then converts only with a partner quantum across a link |

### Q4. g

| code | value |
|---|---|
| G-grid | g/Δ ∈ {0.025, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5}, to measure the admissible window |
| G-read | table constants that have the dimensions of an amplitude, read as `g`: `s_k = 0.06`, `mob_floor = 0.004`, `lock_floor = 0.005`, `gamma_res = 0.25`, `ħΩ = A0 Ω/2π = 0.531`, `t = field_J/2 = 0.9` |

The window is measured first. The readings are then placed against it.
No result is placed by choosing g.

### Q5. Size

| code | geometry or rule |
|---|---|
| Z0 | free rings (R7) |
| Z1 | coalescence: two rings adjacent, with conversion across the links between them (X1 between full cells), against the same rings apart |
| Z2 | frame holonomy: a non-planar ring carries a twist `θ` of the Jones frame on one link (A_l → b_iᵀ ε U(θ) b_j) |
| Z3 | ends and free light: energy per quantum of rings, of open chains with R6's cap-1 ends, and of free light at its band bottom |

### Q6. Formation (what can cool a hot lock)

Not a separate rule. The cooling test of §3 E5 is run under each Q1 and
Q2 option where a channel exists.

---

## 3. The experiments

All at g/Δ = 0.1 unless the experiment varies g. Times are in rate
units; `J` is the ring's fitted exchange.

**E1. Transport inside the ring (Q1 transport rule).** n = 6 ring,
general basis, Σc ≤ 3, S^z = 0, rules T0, T1 (0.05Δ, 0.10Δ), T2, T3, T4.
Report the ground state's total S, its ⟨Σc⟩, the gap over J, and the
L1-style miss of the first 12 levels against the spin-1 ring.

**E2. Leak into empty outside (Q1).** n = 4 ring at its ground state.
One outside lead of L sites hangs from cell 0, empty. Every
(transport rule × ceiling rule) pair. The truncation keeps at most one
escaped quantum, so this measures the first-escape rate. Report the
ring's retained field and the escape probability at T = 100/J.

**E3. A photon arrives (Q1).** As E2, with one photon (one quantum,
random polarization) launched down the lead toward the ring. Two
kinetic energies: cold (near the band bottom) and resonant (`ε_res`).
X0 and X1. Report capture (the photon ends inside the lock or stuck
beside it), escape of ring field, and the ring's disturbance.

**E4. Static ceiling (Q1 C-geo).** Under C-geo, a ring cell at its
ceiling holds `z` quanta, with spin `z/2`, where `z` is its degree.
Report the S-law exchange for cells at cap 3 and 4 from the two-cell
cluster, and the gap of the resulting spin-S chain for S = 1, 3/2, 2
(effective chain DMRG, with the half-integer case expected gapless). Then
a chain of mixed degrees drawn from the foam's degree distribution.

**E5. Cooling by light scattering (Q1, Q4, Q6).** An n = 4 ring
(exact) with one boundary link from cell 0 to the end site of a
semi-infinite photon lead (`t_out = 0.9`).

- Exact multichannel S-matrix: the ring's eigenstates are the
  channels, and the lead is integrated out as a self-energy.
  Unitarity is printed.
- A transition matrix per photon for a cold photon source (kinetic
  energies uniform within the lowest `δ = 0.05` of the band) and for
  a resonant source.
- A master equation over photons, from L6's hot Néel start.
- Report the number of photons needed to reach ground-state
  population 0.9, and the cumulative destruction (absorption plus
  capture) by then.
- Run under X0 and X1, with the boundary link's space mode not
  flowing (F0) and flowing into a 1-D space lead at `t_s = s_k` (F1),
  at g/Δ ∈ {0.05, 0.1, 0.2}.

**E6. Space flow (Q2).**

- (a) F1 along the fiber: inter-link hopping `t_s ∈ {0.03, 0.06, 0.12, 0.24}`
  on the n = 6 ring (general basis). Report the ground S, ⟨Σc⟩, gap/J
  and L1 miss.
- (b) F1 leads: carried from `../lock/` L5(a).
- (c) F2: the first-jump rate `γ⟨Σc⟩` on rings n = 4, 6, 8. Report the
  lifetime in units of 1/J, and the energy that each jump fails to
  account for.
- (d) The table side: the battery with `s_k = 0`.

**E7. The door (Q3).**

- (a) D0: the door's own arithmetic on a lock cell. Report the number
  of steps to halve its field.
- (b) The table side of D2: the battery with `f_conv = 0`.
- (c) Detection under D2: a photon on a screen cell holding m = 1 or
  2 quanta, where the screen link's space is F0 (oscillator), F1 (band
  `t_s = s_k`) or R-c (continuum). Report the absorption probability
  against the photon's kinetic energy across the band.

**E8. g (Q4).** On the G-grid:

- cluster K/J, r3/J, J/(2g²/Δ), and the des Cloizeaux separation;
- the n = 6 exact ring: ground S, gap/J, L1 miss, ⟨Σc⟩;
- the transport margin t_c(g) on the n = 6 ring (Σc ≤ 2);
- the hold: max ⟨c_l⟩ over a sudden-start run of 100/J on n = 6;
- the cooling rate (from E5's g points).

The **admissible window** is the set of g that pass K1 in every one of
these. Report which G-read values fall inside it.

**E9. Size (Q5).**

- (a) Z1: a cube (two 4-rings joined by four rungs, 12 links, every
  cell keeping 2 quanta) against two separate 4-rings, both exact.
  Also the hexagonal prism (two 6-rings) against two separate 6-rings,
  with Σc ≤ 3. Report the energy per cell and the gap.
- (b) Z2: E₀(n, θ) for n = 4, 6, 8, with θ ∈ {0, π/4, π/2, π}.
- (c) Z3: energy per quantum of rings n = 4..10, open chains n = 3..10
  (cap-1 ends), and a free quantum at the band bottom (`E_f(1) − 2t`
  in 1-D, `E_f(1) − 6t` on the cubic lattice).

**E10. The table side of T1.** The battery with
`field_J ∈ {0.18, 0.09}`, which is `t = 0.09` and `0.045`. Those are
about 0.19Δ and 0.09Δ, either side of the S-law margin.

---

## 4. Combinations

The isolated runs cover each axis alone. These combinations are run
because their axes interact.

| code | combination | why |
|---|---|---|
| K-a | every (transport × ceiling) pair, in E2 and E3 | a transport rule's effect at the boundary depends on what the ceiling counts |
| K-b | X1 × {F0, F1 at the boundary link} × g, in E5 | the same link that cools may also destroy |
| K-c | X1 × Z1, in E9(a) | boundary conversion between full cells decides whether locks stay apart or coalesce |
| K-d | D2 × {F0, F1, R-c} at screens, in E7(c) | pair detection needs somewhere for space to go |
| K-e | the candidate bundles below, scored on every test they touch | |

**Candidate bundles** (each is a full set of answers to Q1–Q5):

| bundle | Q1 | Q2 | Q3 | Q4 | Q5 |
|---|---|---|---|---|---|
| B-label | T0, C-lock, X0 | F0 | D1 (label) | window | Z0 |
| B-label-cool | T0, C-lock, X1 | F0 | D1 (label) | window | Z0 |
| B-local | T4, C-held, X1 | F0 inside, F1 elsewhere | D1 (T4 predicate) | window | Z1 as measured |
| B-local-lock | T4, C-lock, X1 | F0 inside, F1 elsewhere | D1 (T4 predicate) | window | Z1 as measured |
| B-geo | T4, C-geo, X1 | F1 | D1 (T4 predicate) | window | Z0 |
| B-table | T2, C-geo, X1 | F2 | D0 | G-read | Z0 |

---

## 5. The scoring rule (fixed now)

Every option and bundle is scored on five criteria:

| | criterion | measured by |
|---|---|---|
| K1 | **The lock lives.** The ground state has S = 0, ⟨Σc⟩/n ≤ 0.1, and an L1 miss below 5%. In every dynamic test the ring keeps at least 0.85 of its field (mean over the second half of the run). | E1–E4, E6, E8 |
| K2 | **The lock can form.** From L6's hot start, light scattering brings the ground-state population to ≥ 0.9 with cumulative destruction ≤ 0.1. | E5 |
| K3 | **The table keeps its claims.** No battery bar moves under the option's global form (law-purity header excluded). A lone photon is detected at a screen. | E6(d), E7(b, c), E10 |
| K4 | **Economy.** The number of new labels, then of new constants. | the tables in §2 |
| K5 | **Margin.** The ratio of the operating point to the nearest measured failure boundary. | E1, E6, E8 |

**Ranking.**

1. A bundle is **admissible** if it passes K1 in every test it
   touches, and K3 wherever K3 is measurable for it.
2. Admissible bundles are ranked by K2 (forms, then does not), then
   by K4 (fewer labels, then fewer constants), then by K5.
3. Where K3 cannot be measured without a kernel step, the bundle is
   ranked and flagged "K3 unmeasured".

---

## 6. What does not count

- Choosing g, a reading or a rule because a test passes.
- Tuning `w1`, `w2`, `q`, `cap` or `A0`.
- Softening a criterion after a result.
- Reading ξ, or any length, as a cell count.
- Mapping anything onto 938 MeV, 0.84 fm or 1836.
- Counting a battery run with changed constants as an adoption.
