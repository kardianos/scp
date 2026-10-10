# v102 vacancy step — pre-registration (the particle as a vacancy)

**Written:** 2026-10-10, before any run of this step.
**Authority:** user, 2026-10-10: "That's a fascinating take. The
particle as a vacancy." Then: "Use the permutation modes to try to find
answers to your questions. Please proceed with this 'particle as
vacancy' idea." The standing rule from the same day: every decision is
made from experiments over all combinations of its candidates (a, b, c,
ab, ac, bc, abc). New code is Go or C; Python only where a C library
(LAPACK, BLAS) does the work.

**A measurement in models. It rules nothing.** No kernel change, no law
change, no `v97/fab` change. Constants as in every v102 step: the
S-law, g/Δ = 0.1, J0 = 2g²/Δ = 0.0095862, the T-floor t = 0.1 J0 (R9),
R12 on (1,1) links.

## 0. The working reading

**A nucleon is a vacancy:** one quantum missing from a background of
locked field. The background (the lock, a fiber at its R6 ceiling) is
not the particle. It is what the particle is missing from. This is the
reading the fabric and permutation steps already measured as "the
hole". What is new is taking it as the identity, and asking what
follows.

Three questions were left open when the reading was named:

1. Is the particle the vacancy, the whole lock, or something else?
2. If it is the vacancy, what sets its size?
3. In this reading the lock's own ripple (the lock quantum) plays the
   role of a Higgs boson: a vibration of the filled background, not an
   object in it. Is there a ripple with spin 0 and isospin 0, and is it
   a line of its own, or does it fall apart into two ripples?

A fourth is reported, not decided: what the background is made of.

## 1. Models

| name | what | where it runs |
|---|---|---|
| **LO** | the leading-order effective lock of `../perm/cas/lo.py` (reading V): exchange −(g²/Δ_type)A†A on each link, T-floor hop, R12 as the singlet core on (1,1) links | exact (C) on rings and open chains; DMRG on long chains |
| **AO** | all orders in g. On small rings and chains: the microscopic S-law of `../perm/cas/graphm.py` (field energies E_f(N), link oscillators at Ω, conversion g, hop t), exact, with space quanta on links. On long 3+2 chains: the des Cloizeaux terms h2 + r3 (`../fabric/cas/eff.py`, exact to fourth order), as in Phase 1 P3 and fabric F4b | exact (C); DMRG (3+2) |

Open chains carry R6 end cells (N = 1). The R6 ceiling holds in AO as
N_i + Σ_{links at i} c_l = cap_i at t = 0, and the T-floor hop may not
empty a cell or fill one past its cap.

**Tools.**
- `cas/vx.c` (new, C, OpenMP, LAPACK `dstev`): exact diagonalization in
  a charge sector, LO or AO, ring or open chain, static or mobile
  vacancy. Matrix-free Lanczos with full reorthogonalization. Spin-0
  (and isospin-0) states are found by adding λ(S² + I²) with
  λ = 5 J0, which lifts every other multiplet by ≥ 10 J0, and are
  checked by ⟨S²⟩ and ⟨I²⟩. It reports profiles, correlations, link
  space content and overlaps.
- DMRG on long chains reuses the drivers already used for these chains
  (TeNPy for LO, `../phase1/cas/dmrg.py` for 3+2 AO). Both are numpy
  over BLAS/LAPACK. They gain a correlation readout only.
- Fits and tables: `cas/fit.go` (new, Go).

**Tool gate T0.** Before any new number is read, `vx` must reproduce
these recorded values to their printed digits:

| case | recorded | source |
|---|---|---|
| LO W ring n = 6, one mobile vacancy, t = 0 | E0 = −20.51293 J0; lowest (3/2, ½) at +1.9917 J0 | `../perm/cas/out/x3_rings_M4_n6.out` |
| same, t = 0.1 J0 | E0 = −20.73466 J0; band width 0.38794 J0 | same |
| LO 3+2 ring n = 8, one mobile vacancy, t = 0 / 0.1 J0 | E0 = −19.66294 / −19.72715 J0; band 0.12258 J0 at 0.1 | `x3_rings_M2_n8.out` |
| AO 3+2 ring n = 4, sector 0 | lowest −0.091490, −0.081985, −0.063837 | `w_lock.out` |
| AO W ring n = 4, sector (0, 0) | lowest −0.127331, −0.108256 | `w_lock.out` |
| AO 3+2 ring n = 8, Σc ≤ 3 | Δ₈ = 0.005549; lowest S = 0 excitation 0.024014 = 2.1639 × 2Δ₈ | `../lock/cas/out/l4b_line.out` |

A miss stops the step until it is found.

## 2. Q1 — which object is the particle?

**Candidates, all seven combinations:**

| | candidate |
|---|---|
| a | the vacancy (one missing quantum and its cloud) |
| b | the background (the whole lock) |
| c | the background's ripple (one lock quantum) |
| ab | the background holding one vacancy, taken as one object (v98's "lock with an open slot") |
| ac | a vacancy carrying a ripple (an excited vacancy) |
| bc | the background carrying one ripple |
| abc | the background holding a vacancy and a ripple |

**Scoring:** the registered nucleon criteria N1–N8 of
`../fabric/PREREG.md`, read under R8–R10 as the permutation step read
them: N3 is the (½, ½) quartet in W; N7 is p and n both present as
(½, ½), charge I^z + ½; N6 and N8 are scored at t = 0.1 J0 (R9) and
with the S-boundary (R10).

**Sources.** Every cell is scored from the records (`../fabric/`,
`../perm/`, `../attract/`, `../lock/`) or from Q2 and Q3 of this step.
No run is registered for Q1 alone. Each cell names its source.

**Endings.** **V**: a scores strictly highest. **L**: a candidate that
contains b scores at least as high as a. **X**: no candidate passes
N1–N4.

## 3. Q2 — does the particle take its size from the background?

The vacancy's size ℓ_vac is the registered N2 quantity: the decay length
of the vacancy-centred profile, a pure exponential fit of |⟨S^z⟩| (a
static vacancy) or of the conditional profile C(d) (a mobile one). The
background's correlation length ξ_bg is P2's quantity: −1/slope of
log(|⟨S_i·S_{i+r}⟩|·√r) in the same model without the vacancy.

**Factors, all eight cells:** a = W (against 3+2); b = the T-floor
t = 0.1 J0 (against t = 0); c = AO (against LO).

| cell | ℓ_vac | ξ_bg |
|---|---|---|
| 3+2 LO t=0 | record: 6.074, static, own DMRG L_bulk = 40 (`../perm/cas/out/analyze_chains_ownDMRG.out`) | record: LO is the Heisenberg chain; Phase 1 P2, 6.137 (L = 80, r ∈ [6, 24]) |
| 3+2 LO t=0.1 | record: the conditional profile, 5.380 at L_bulk = 40 (`analyze_tp.out`) | as above |
| 3+2 AO t=0 | record: fabric F4b, 5.80 (L_bulk = 56) | record: Phase 1 P3, 5.79 (L = 48, r ∈ [6, 20]) |
| W LO t=0 | record: 2.299, own range d = 2..14, TeNPy L_bulk = 32 (`analyze_tp.out`) | **new:** TeNPy L_bulk = 32, r ∈ [3, 12] |
| W LO t=0.1 | **new:** TeNPy with the 14-state mobile site, L_bulk = 24, the conditional profile | as W LO t=0 |
| 3+2 AO t=0.1, W AO t=0, W AO t=0.1 | no long-chain AO model with these. **Ring scale:** on exact rings, the ratio of the vacancy profile's decay to the background correlation's decay, fitted over the same window, in AO and in LO at the same n | same rings |

**Disclosure.** The three 3+2 long-chain cells are records. Their ratios
could be read off when this page was written: 6.074/6.137 = 0.990,
5.380/6.137 = 0.877, 5.80/5.79 = 1.002. The conditional profile is
chosen for the mobile cell because the raw C(d) carries the vacancy's
particle-in-a-box weight (`../perm/DECIDE.md` §7.4). The raw value
(5.318, ratio 0.867) gives the same verdict.

The ring-scale cells test only that AO keeps the LO relation at the
largest ring each can run (3+2 n = 8, W n = 6). They do not measure
ℓ_vac or ξ_bg.

**Bar.** A long-chain cell **inherits** iff |ℓ_vac/ξ_bg − 1| ≤ 10%. A
ring-scale cell **keeps** iff its AO ratio is within 10% of its LO
ratio.

**Covariance (supplement, 3+2 AO, static).** ℓ_vac at g/Δ = 0.05 and
0.2 (DMRG L_bulk = 40), against P3's ξ_bg = 6.06 and 5.07 at the same
g/Δ. It passes iff ℓ_vac/ξ_bg stays within 10% of 1 at all three g/Δ.
These are the g/Δ values Phase 1 already ran. No value is selected and
no line is placed: the test asks whether ℓ_vac moves with ξ_bg.

**Endings.** **I**: every long-chain cell inherits, every ring-scale
cell keeps, and the covariance passes. **P**: some do. **N**: no
long-chain cell inherits.

## 4. Q3 — is there a Higgs-like line?

At each size: E0; the ripple Δ_n (3+2: the lowest S = 1 state; W: the
lowest state of sector (2S^z, 2I^z) = (2, 2)); ε_s, the lowest
excitation with S = 0 (3+2) or (S, I) = (0, 0) (W); and **R = ε_s/(2Δ_n)**.

**Factors, all eight cells:** a = W (against 3+2); b = AO (against LO);
c = open chain with R6 ends (against ring).

| | rings (n) | open chains (L_bulk) |
|---|---|---|
| LO 3+2 | 6, 8, 10, 12, 14 | 6, 8, 10, 12 |
| LO W | 4, 6, 8 | 3, 4, 5, 6 |
| AO 3+2 | 4, 6, 8 (Σc ≤ 3 at 8) | 4, 6 |
| AO W | 4, 6 | 3, 4 |

AO runs use the full closure where it fits, and otherwise Σc ≤ 3. The
truncation is reported, with the shift of E0 and Δ_n from Σc ≤ 2.

**Verdict per cell.** **Line** iff R < 1 at both of the two largest
sizes and 1 − R at the largest is at least half of 1 − R at the next.
**No line** iff R ≥ 1 at the largest size. Otherwise **unclear**.

**Endings.** **H**: every W cell is a line. **h**: some W cell is a
line. **0**: no W cell is a line; the spin-0, isospin-0 ripple falls
apart into two ripples. The 3+2 cells are reported beside them.

Also reported: the ripple's own multiplet; whether the lowest (0, 0)
excitation shares its energy with other (S, I) members (W's larger
symmetry); and ε_s and R on the largest W ring with one vacancy present
(the ripple beside a particle).

## 5. Q4 — what the background is made of (reported)

The S-law converts one space quantum on a link into one conversion
singlet pair, one quantum in each end cell. The state

  |CS⟩ = Π_links A†_l |0⟩,  A†_l = Σ P_μν b†_iμ b†_jν,

is "one converted space quantum per link". It meets R6 by construction:
every cell holds one quantum per link end. Reported, on the largest
rings and chains of Q3, LO and AO:

- the per-cell fidelity |⟨CS|ψ0⟩|^(2/n) of the ground state;
- with one static vacancy, the same for b_hμ|CS⟩ against the vacancy's
  ground state in the matching sector;
- in AO, the space left on each link, ⟨c_l⟩, and the check
  ⟨N_i⟩ + Σ⟨c_l⟩ = cap_i.

## 6. Order

1. T0. 2. Q3 (exact, C). 3. Q2 (DMRG in the background, ring scale in
C). 4. Q4 (from the Q3 states). 5. Q1 (scored last, from the records and
Q2–Q3). 6. `DECIDE.md`.

## 7. What this step does not do

- It does not make a background from field. The lock is still placed.
  How a background gathers is the next question; this step only
  reports what the background is made of.
- It does not map anything onto 938 MeV, 0.84 fm, 1836, quark species,
  the pion, or the Higgs boson's 125 GeV. "Higgs-like" names a role (a
  ripple of the filled background), not a particle of the world.
- It does not take n from ξ, and reads no length as a cell count.
