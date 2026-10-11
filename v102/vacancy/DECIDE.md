# v102 vacancy step — decision (the particle as a vacancy)

**Written:** 2026-10-10, after the registered runs. Gates: `PREREG.md`.
A measurement in models; it rules nothing. No kernel, law or `v97/fab`
change.

- Exact diagonalization: `cas/vx.c` (C, OpenMP, LAPACK/BLAS). Runs:
  `cas/run_t0.sh`, `cas/run_q3.sh` (the Q3 grid and the Q2 ring-scale
  runs), `cas/run_q3s.sh` (a labelled supplement).
- Long chains: `cas/run_q2.sh` (`cas/vq2_tp.py`, a thin TeNPy driver;
  `../perm/cas/tp_chain.py`; `../attract/cas/a2_allorder.py`).
- Fits and tables: `cas/fit.go` (Go) → `cas/out/fit_{t0,q2,q3,q4}.out`.
- Process notes, including two faults found and fixed during the runs:
  §6.

## Endings

| question | ending | in one line |
|---|---|---|
| **Q1** which object is the particle | **V** | the vacancy scores 7 of 8 registered criteria; the next best scores 4 |
| **Q2** what sets its size | **P** | the vacancy's size is the background's correlation length, and follows it when g/Δ moves it; three cells miss the bar (§2) |
| **Q3** a Higgs-like line | **h** | 3+2: no line in any cell. W: the spin-0, isospin-0 state sits at the two-ripple threshold to within 1% at leading order; only the all-orders rings put it about 10% below |
| **Q4** what the background is made of | reported | "one converted space quantum per link" is the background to 0.96–0.997 per cell at leading order |

**What the reading delivers.** The vacancy is the particle, and its
size is not a choice: it is the correlation length of the locked field
around it. That length is a property of the fabric (6.0 grains for 3+2,
2.1 for W at leading order). No test of the lock's own length ever
selected a size (`../options/` §6, `../kernel/K3_DECIDE.md` §6). The
vacancy reading moves "size" from the lock, where nothing fixes it, to
the vacancy, where the fabric fixes it.

## 0. Tool gate T0: **PASS**

All 15 checks reproduce the recorded values to their printed digits
(`cas/out/fit_t0.out`). Among them: LO W ring n = 6 with one mobile
vacancy, E0 = −20.512926 J0 (recorded −20.51293), band 0.387936 J0
(0.38794); AO 3+2 ring n = 8 (Σc ≤ 3), Δ₈ = 0.005549 and the lowest
S = 0 excitation at 2.16389 × 2Δ₈ (2.1639). `vx` runs these 100–300×
faster than the Python tools that made the records.

## 1. Q3 — is there a Higgs-like line? **Ending h**

R = ε_s/(2Δ_n): ε_s is the lowest excitation with S = 0 (W: S = 0 and
I = 0); Δ_n is the lowest ripple. R < 1 means the spin-0 state lies
below two free ripples. Every ε_s state was checked: ⟨S²⟩ = ⟨I²⟩ = 0 to
1e-6. Every run converged, with residual ≤ 1e-8 rate units.

| cell | sizes | R | verdict |
|---|---|---|---|
| 3+2 LO ring | n = 6 … 14 | 2.20, 2.26, 2.14, 1.94, 1.78 | no line |
| 3+2 LO open | L_bulk = 6 … 12 | 1.72, 1.62, 1.53, 1.46 | no line |
| 3+2 AO ring | n = 4, 6, 8 | 1.89, 2.08, 2.15 | no line |
| 3+2 AO open | L_bulk = 4, 6 | 1.76, 1.69 | no line |
| W LO ring | n = 4, 6, 8 | **1.000, 0.993, 1.009** | no line |
| W LO open | L_bulk = 3 … 6 | 1.26, 1.11, 1.135, 1.10 | no line |
| W AO ring | n = 4, 6 | **0.899, 0.890** | **line** |
| W AO open | L_bulk = 3, 4 | 1.235, 1.091 | no line |

**Supplement (labelled, not registered; `cas/run_q3s.sh`):** the next
size in two W cells (`cas/out/fit_q3s.out`):

| cell | size | Δ | ε_s | R | the cell's sequence |
|---|---|---|---|---|---|
| W AO ring | n = 8 (D = 1.26e7, full closure) | 1.3526 | 2.4735 | **0.914** | 0.899, 0.890, 0.914 |
| W LO open | L_bulk = 7 (D = 5.2e6) | 1.6523 | 3.5809 | **1.084** | 1.26, 1.11, 1.135, 1.10, 1.084 |

- **The AO ring stays below the threshold.** It still reads as a line
  by the registered rule. But its margin 1 − R went 0.101, 0.110, 0.086
  for n = 4, 6, 8. From 6 to 8 that shrinks about as 1/n (0.110 × 6/8 =
  0.083), the way a finite-size attraction shift would. It does not
  settle the question.
- **The LO open chain approaches 1 from above,** as two ripples in a box
  do.

**Reading.**
- **3+2: the spin-0 ripple falls apart.** In every cell it sits 1.46–2.26
  times above twice the ripple. At leading order R falls with size (rings
  from n = 8, open chains throughout), as a two-ripple scattering state's
  does; the small all-orders rings have not reached that regime. This
  agrees with the lock step's 2.16 at n = 8.
- **W sits on the threshold.** At leading order the lowest (0, 0) state
  is at 2Δ_n to within 1% on every ring: exactly at n = 4, where Δ = 2
  and ε_s = 4 J0. It is above 2Δ_n in open chains, approaching from
  1.26 to 1.10. That is two ripples with almost no interaction in the
  spin-0, isospin-0 channel.
- **Only the all-orders rings go below,** by about 10% at n = 4 and 6.
  The fourth-order terms make the channel weakly attractive on small
  rings. The open chains at all orders do not show it. In one
  dimension any attraction binds two particles, so a weak line is
  possible. These sizes cannot separate a weak line from a finite-size
  shift; the n = 8 supplement (R = 0.914, margin shrinking from n = 6)
  does not separate them either. The registered rule reads the AO ring
  as a line, and the ending is h.
- **The ripple itself is not spin 0.** It is a triplet (3+2) or the
  15-fold (1, 1) ⊕ (1, 0) ⊕ (0, 1) (W): ⟨S²⟩ = 2 in every cell, and
  ⟨I²⟩ = 2 in W. A Higgs-like scalar could only be a pair of ripples.

**The ripple beside a particle** (LO W ring n = 8, one static vacancy):
nothing lies below the vacancy plus one ripple. The vacancy's lowest
excitations are at +1.66 J0 (spin 3/2 members) and +1.68 J0 (½, ½),
against a threshold of +1.37 J0. There is no excited vacancy at t = 0.
This agrees with the attract step's "no ladder at small t".

## 2. Q2 — does the particle take its size from the background? **Ending P**

ℓ_vac is the registered N2 decay length (pure exponential, the
vacancy-centred or conditional profile). ξ_bg is P2's
Ornstein–Zernike fit of the background's ⟨S_i·S_{i+r}⟩.

| cell | ℓ_vac | ξ_bg | ratio | verdict |
|---|---|---|---|---|
| 3+2 LO, static | 6.074 (record) | 6.137 (record) | 0.990 | inherits |
| 3+2 LO, mobile t = 0.1 J0 | 5.380 (record, pre-R12) | 6.137 | **0.877** | **misses** |
| 3+2 AO, static | 5.80 (record) | 5.79 (record) | 1.002 | inherits |
| W LO, static | 2.299 (record) | **2.116** (new, L_bulk = 32; I^z: 2.127) | 1.087 | inherits |
| W LO, mobile t = 0.1 J0 | **2.090** (new, R12, L_bulk = 24) | 2.116 | 0.988 | inherits |
| 3+2 AO, mobile (ring scale, Σc ≤ 3) | AO/LO ratio of ratios | | 1.085 | keeps |
| W AO, static (ring scale) | | | **1.329** | **does not keep** |
| W AO, mobile (ring scale) | | | **1.263** | **does not keep** |

**Covariance (3+2, all orders, static), passes.** As g/Δ moves the
background's ξ_bg, the vacancy's size follows it:

| g/Δ | ξ_bg (Phase 1 P3) | ℓ_vac (L_bulk = 40, new except 0.1) | ratio |
|---|---|---|---|
| 0.05 | 6.06 | 6.133 | 1.012 |
| 0.1 | 5.79 | 5.91 (fabric F4b; 5.80 at L_bulk = 56) | 1.021 (1.002) |
| 0.2 | 5.07 | 5.172 | 1.020 |

ξ_bg moves by 16%, and the ratio stays within 2% of 1.

**The three misses.**
1. **The moving 3+2 vacancy is 12% tighter** than ξ_bg. The moving W
   vacancy is also tighter than the static one (2.09 against 2.30),
   but stays within the bar. Motion tightens the cloud by about 10% in
   both fabrics. The 3+2 record was run without R12; in W, R12 moved
   the size by 0.6% (2.077 → 2.090).
2. **The two W ring-scale cells miss on an instrument fault.** On a
   6-ring the W vacancy's profile is a standing pattern. Its magnitude
   at d = 3 exceeds that at d = 2 (0.82 against 0.53 of the d = 1
   value), as in the long chain where d = 2 and 3 are equal. A
   three-point "decay" fit is then ill-conditioned: it reads 9.9 (LO)
   against 12.6 (AO), from profiles that agree point by point to
   0.037. The 3+2 cells (an 8-ring) keep: 1.06 and 1.085. The verdict
   is the registered one; the profiles are in `cas/out/fit_q2.out`,
   labelled as a supplement.

The ring-scale cells test only that AO keeps LO's relation (`PREREG.md`
§3). The pointwise supplement says it does, to within 0.04 of the d = 1
value in every cell.

## 3. Q4 — what the background is made of (reported)

|Ψ_CS⟩ = Π_links A†_l|0⟩ is "one converted space quantum per link".
Per-cell fidelity |⟨CS|ψ0⟩|^(2/n):

| | 3+2 ring | 3+2 open | W ring | W open |
|---|---|---|---|---|
| LO | 0.971–0.975 (n = 6 … 14) | 0.984–0.990 | 0.961–0.969 (n = 4 … 8) | 0.990–0.997 |
| AO | 0.930–0.944 | 0.939–0.940 | 0.906–0.907 | 0.915–0.917 |

- **The background is converted space,** one quantum per link, to 96–99%
  per cell at leading order. The R6 ceiling (one quantum per link end)
  is what that conversion makes. It is the same statement read from the
  other side.
- **In AO about 4–6% of each link is still space:** ⟨c_l⟩ = 0.044 (3+2
  rings) and 0.060 (W rings), up to 0.14 next to an open end. The R6
  balance holds exactly: ⟨N_i⟩ + Σ⟨c_l⟩ = cap_i (3+2 n = 8: 1.9127 +
  2 × 0.0436 = 2.0000).
- **A vacancy is that state with one member removed.** b_hμ|CS⟩ against
  the static vacancy's ground state: 0.962 per cell (LO 3+2 n = 8) and
  0.944 (LO W n = 6); 0.916 and 0.883 in AO.

## 4. Q1 — which object is the particle? **Ending V**

The fabric step's criteria, read for W (R8) at t = 0.1 J0 (R9) with the
S-boundary (R10). Sources: the records and §1–§3.

| | a vacancy | b background | c ripple | ab background + vacancy | ac vacancy + ripple | bc background + ripple | abc |
|---|---|---|---|---|---|---|---|
| N1 from the fabric | ✔ a defect | ✘ the fiber itself, its extent a label | ✔ | ✘ contains the label | ✔ | ✘ | ✘ |
| N2 one size | ✔ 2.37 → 2.30 | ✘ its extent is the label | — not measured in W (3+2 ✔) | ✘ the label | ✘ no bound level (§1; `../attract/` §4) | ✘ | ✘ |
| N3 (½, ½) | ✔ | ✘ (0, 0) | ✘ 15-fold | ✔ | ✘ nothing bound (§1) | ✘ | ✘ |
| N4 firm | ✔ gap ≈ Δ | ✔ gap 1.32 J0 | ✔ | ✔ | ✘ | ✘ a band | ✘ |
| N5 ladder | ✘ | ✘ | — | ✘ | — | ✘ | ✘ |
| N6 moves by conversion | ✔ band 0.388 J0 (n = 6) | ✘ | ✔ | ✘ the fiber stays put | — | ✘ | ✘ |
| N7 p and n | ✔ | ✘ | ✘ | ✔ | ✘ | ✘ | ✘ |
| N8 forms from field | ✔ ≈ 5% per resonant photon | ✘ placed (≈ 10⁸ cold photons in models) | ✔ 12–28% (3+2) | ✘ the background is placed | ✘ | ✘ | ✘ |
| **score** | **7** | 1 | 4 (5 with 3+2's N2) | 3 | 1 | 0 | 0 |

- **The vacancy wins on size.** Taking the background together with its
  vacancy (v98's "lock with an open slot", ab) gives the same quantum
  numbers and the same p/n pair. But its size is the lock's label,
  which nothing fixes. The vacancy alone has a size the fabric fixes.
- **The ripple loses on spin.** It is a triplet or the 15-fold
  multiplet, never (½, ½), and it carries no p/n pair.
- **N5 fails for every candidate.** There is no excited level at small
  transport (§1).

## 5. What the vacancy reading now says

1. **The particle is a vacancy.** One quantum missing from a background
   of converted space (§3), carrying (½, ½): p and n by the isospin of
   the gap.
2. **Its size is the background's correlation length.** It is fixed by
   the fabric and follows the fabric when the coupling moves it (§2).
   A moving vacancy is about 10% tighter.
3. **The background's ripple is not a Higgs boson.** It is never spin 0.
   The spin-0, isospin-0 state of two ripples sits at their threshold:
   above it in 3+2, on it in W, slightly below only on small rings at
   all orders. If W has a scalar line, it is weak; these sizes do not
   decide it.
4. **The background is still placed.** Nothing here makes it from field.
   That is now the one link of field → particle the models do not
   carry.

## 6. Process notes

1. **Solver changes during the step, each checked before use.**
   - Restarted Lanczos (the λ(S² + I²) penalty widens the spectrum).
   - S² applied as S⁻(S⁺ψ) through the neighbouring sector, instead of
     n² cell pairs.
   - Reorthogonalization by BLAS `dgemv`.
   The new binary reproduced the finished runs to every printed digit
   (3+2 LO rings n = 8 and 10, AO n = 8, W LO n = 4) before it
   replaced the old one.
2. **A relaunch fault, caught and undone.** In the first parallel
   relaunch, the job runner dropped each job's first argument
   (`model=`), so 28 AO jobs ran the LO model. The output headers showed
   it. Those outputs were deleted and rerun; every AO output now heads
   `model=ao`. The LO jobs were unaffected, because LO is the default.
3. **Thread contention.** Two processes outside this step held two
   cores. Small jobs therefore ran single-threaded, 12 at once.
4. **The AO 3+2 mobile ring did not fit as a full closure.** With
   transport the closure kept growing past 38 minutes. It was stopped
   and run at Σc ≤ 3 (D = 2.9e6), as `PREREG.md` §4 allows. Against
   Σc ≤ 2: E0 is 0.49 J0 lower, and the profile moves by ≤ 0.005.
5. **R12 in the long chains.** The W mobile cell was first run with the
   perm step's default (no term on (1,1) links). It was rerun under
   R12 as registered: 2.077 → 2.090. The 3+2 mobile record (5.380)
   predates R12.

## 7. Limits

- Exact sizes are small: W rings up to n = 8 (n = 6 registered at all
  orders; n = 8 as a supplement). The W Higgs question sits at R ≈ 1, where finite-size
  shifts of a few percent decide the verdict.
- ξ_bg and ℓ_vac use different fit forms (Ornstein–Zernike for the bulk
  correlation, a pure exponential for the vacancy's cloud). That is the
  registered pairing; a pure exponential on both would read W at 1.26.
  In 3+2 the same pairing gives 0.99–1.00.
- Nothing here maps onto 938 MeV, 0.84 fm, 1836, quarks, the pion or
  125 GeV. "Higgs-like" names a role.

## Then

The vacancy reading leaves one open link: **the background is placed**.
Under the principle ("matter is converted space") the background is
converted space (§3), so the question is how space and field gather
into it. A next step would test, in all combinations, candidates such
as:

| | candidate | what it would do |
|---|---|---|
| a | the reverse boundary event | space arriving at a background's edge converts into one locked quantum plus one traveling quantum (light): the principle's direction, refilling vacancies |
| b | conversion in an unlabelled fiber | a run of links each converting one space quantum into a pair, with the released energy (about Ω − E_f(2) = 0.23 per link) leaving as light |
| c | membership by saturation | a cell belongs to the background when every link end holds a converted pair; no label |

Each would need a ruling before it entered the law; the test itself
would rule nothing.
