# v102 kernel step, K2b — decision (the widened field)

**Written:** 2026-10-08, after the runs in `cas/out/k2b_*`. Gates:
`K2b_PREREG.md`. Code:
- `v97/kernel/freecell.c`: `fmodes`, `fseed_mode`, the mode-blind
  field;
- `v97/kernel/lockq.c`: per-mode events, `lockq_fmode`.

## Ending: **Go, and stop before K3**

The kernel field now carries 1, 2 or 4 complex modes.

- **At one mode it is byte-identical to before.**
- **At four modes with an equal split it is byte-identical to one mode,
  step for step.**
- **Light in mode ν makes holes with mode ν's charges, every time.**

Two bugs were found during K2b and fixed before the commit (§B4). One
of them was in K2a's committed code; K2a's recorded results are
unaffected.

K3, as `PLAN.md` describes it ("nucleons form, walk, and pair"), meets
two physics walls (§5). Both need a ruling before K3 is designed.

## B1. Inert at default: **PASS**

- **Kernel stdout,** K2b binary against K2a binary, V3a law, T = 20:
  `exp=bath`, `blob`, `slit`, `pulse` are **byte-identical**, and so is
  an every-step `pulse` run (1001 rows).
- **The full battery, both kernels, on the final binary**
  (`cas/out/k2b_battery.out`): **ALL GREEN (87 bars)**, bar table
  identical to the pre-K2a baseline.
  - 70 of 71 run logs are byte-identical.
  - `p2lc` (the Go kernel) differs only in wall-clock timings.
- There are no new compiler warnings (9 before, 9 after).

## B2. Exact reduction (F = 4 with an equal split, against F = 1): **PASS**

- `exp=pulse` (amp 1, σ 2.2) and `exp=blob` (L = 16), T = 20: every
  `t=` row and every `# RESULT` line are byte-identical.
- The every-step check (`diag_every=1`, 1001 rows plus the results) is
  also byte-identical.
- **Why it is exact.**
  - Each mode amplitude is exactly half the one-mode amplitude, because
    scaling by ½ is exact in floating point.
  - The pairwise energy sum (m0 + m1) + (m2 + m3) of four equal
    quarters is exact.
  - Every field operation is the same linear map on each mode.
- **F = 2 with the split (reported).** 1/√2 rounds, so 120 of 1001
  rows differ in print, and only in the drift column, at the 1e-16
  level. Every other printed value is identical (largest relative
  deviation 0 for values above 1e-8).

## B3. Mode permutation is a symmetry: **PASS**

- F = 4, `exp=pulse`, light seeded in mode 0, 1, 2 or 3: every `t=`
  row and every `# RESULT` line are byte-identical across the four,
  every step.
- **Energy outside the seeded mode: exactly 0.** Over T = 20 no cell
  went dark and was re-lit by emission, so the empty-cell rule never
  fired (`cas/out/k2b_b3s_m*_result.out`: all 110.5709616526 in the
  seeded mode).

## B4. Per-mode events: **PASS, after two bug fixes**

**The registered run found two bugs.**

1. **A heap overflow in `lockq.c`** (K2a code).
   - The fence array `lq_rawA` was sized by `NSLOT`. That is the slot
     table's high-water mark, not its capacity (`NLMAX`).
   - A new slot touching the lock above the old high-water mark was
     written past the end. One W run (mode 0, seed 2) died of it at
     t ≈ 7, and AddressSanitizer located it.
   - Fixed: sized by `NLMAX`.
   - **K2a's recorded results are unaffected.** All 12 of K2a's G4
     configurations, rerun with the fixed binary, give event lines and
     results identical to `k2a_g4_*.out`
     (`cas/out/k2b_fixcheck.out`).
2. **A ledger leak in the lock's installation** (new K2b path).
   - It zeroed only mode 0 of the lock cells. Light seeded into modes
     1–3 stayed in the lock cells and came back into `Ee` at the next
     pass F.
   - This gave drift of 1.6–3.3e-3 in every run seeded outside mode 0.
   - Fixed: all modes are zeroed.

AddressSanitizer then ran clean on five configurations (F = 1, 2, 4;
both fabrics; seeded and split light).

**The re-run** (24 runs: 3+2 n = 8 with F = 2, W n = 6 with F = 4;
light in each mode and split; three seeds; `cas/out/k2b_M*`):

| lock | light | events (seeds 1/2/3) | holes made, by mode | final charges |
|---|---|---|---|---|
| 3+2 | mode 0 | 6 / 5 / 6 | all mode 0 | (+6), (+5), (+6) |
| 3+2 | mode 1 | 6 / 5 / 6 | all mode 1 | (−6), (−5), (−6) |
| 3+2 | split | 6 / 6 / 6 | 3+3, 4+2, 3+3 | (0), (+2), (0) |
| W | mode ν (each of 4) | 4 / 2 / 4 | all mode ν | ν's charges × N_h |
| W | split | 3 / 3 / 3 | (1,1,0,1), (0,2,1,0), (0,2,1,0) | (1,−1) in all three |

- **The in-kernel checks fail 0 times in all 24 runs.** Each event took
  exactly one atom from mode ν, and moved the lock by mode ν's
  charges.
- No event fell outside the band.
- |drift_rel| ≤ 3.6e-15.

**The nucleus's sector is reached.** Under split light every W run
passes through two holes at zero charge, (S^z, I^z) = (0, 0), with the
holes on neighbouring sites 1–2:
- seed 1: modes (+,+) then (−,−);
- seeds 2 and 3: (+,−) then (−,+).

- **These W visits are adjacent-hole configurations,** the case of
  Wall 1 (§5).
- 3+2 seed 2 reaches (0) with its holes two sites apart (sites 1 and
  3).
- The bright light then adds a third hole within 0.5–4 time units.

## B5. Per-mode A/B: **PASS**

- Seeded jumps with photon mode ν (W ν = 1, 2, 3 at site 2; 3+2 ν = 1
  at site 3), then 2000 steps.
- The hole, S^z and I^z profiles agree with `cas/k2a_ab.py`
  (extended with ν) **to every printed digit** (max |C − Python| = 0
  at 10 digits).
- The sectors match: W ν = 1 → (1, −1, 1), ν = 2 → (−1, 1, 1),
  ν = 3 → (−1, −1, 1); 3+2 ν = 1 → (−1, 1).

## 5. Why K3 should not start yet: two walls

**Wall 1: adjacent holes are outside the ruled physics.**

- Under bright light the lock fills with neighbouring holes: patterns
  like 011110 occur in most bright runs, W and 3+2.
- The leading-order lock gives a link between two holes, a (1,1) link,
  no exchange. Its conversion cost, Δ₁₁ = Ω − 2E_f(1) = −0.024, is
  below zero.
- `../perm/DECIDE.md` §5 measured this at all orders:
  - adjacent holes **bind by converting**, using real link occupancy,
    which reading V forbids;
  - with transport, two holes **merge into an empty cell**, which R9
    forbids.
- So whenever two holes touch, the kernel's lock is in a configuration
  the rulings leave unpriced. Dim light (amp 1) makes one or two
  separated holes and stays inside the model.

**Wall 2: pairing needs a sink, and the kernel has none fast enough.**

- K1-2: a closed lock never binds two holes. The binding (W: 0.50 J0)
  is smaller than the lock's own gap, so it has nowhere to go.
- K1-3: cold light captures a pair, at about 10¹⁰ photons per pair.
- The kernel lock's only door is the formation event. So "nucleons
  pair" cannot happen in any reachable kernel run under the present
  rules.

**What K3 can do now** is formation and walking with dim light:
- separated holes made by light;
- carried by the T-floor hop (one hop per about 10³ time units);
- with the space quanta delivered to the foam.

That is a demonstration of what K2a/K2b already validate. It is not
the nucleus.

## 6. Limits stated

- **Mode-0 only.** The FCS snapshot and the loop-phase meters
  (`fa_phase`) read mode 0 only.
- **Refused at F > 1:** `amp_door`, `qp_phase`, `door_g`, `bath_xb`,
  `bh_r`.
- **The Go kernel stays single-mode.**
- **No field charge ledger.** Matter carries no spin or isospin, so
  the classical kernel does not conserve the field's S^z or I^z. Only
  the lock's door exchanges them exactly.
