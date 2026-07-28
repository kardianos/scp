# Root docs + new git changes — accuracy review

**Date**: 2026-07-08  
**Scope**: Unstaged root docs and related new work; spot-check of older CONCEPT claims.

## Scope reviewed

- **Unstaged**: `CONCEPT.md`, `DISCOVERIES.md`, `FUTURE.md`, `sfa/Makefile`, `sfa/seed/eta_relaxer.c`, untracked `v72/`, `sfa/seed/eta_qflow.c`
- **Root docs**: `README.md`, `CLAUDE.md`, `EM_THEORY.md` (currency)
- **Cross-check**: v72 diag TSVs/logs vs claims; sample older CONCEPT numbers vs v69–v71 sources

**Bottom line:** The v72 physics story is real and the main drift/VK numbers match the data. There are a few **documentation overclaims**, one **false “fixed” claim** about tools, and **stale era pointers** in README/CLAUDE.

---

## Verified against data (v72)

Recomputed drifts from `v72/results/*_diag.tsv` with the same 0→60 definition as `drift.py`:

| run | claimed | measured |
|---|---|---|
| r_qflow2 (η=0.25 stationary) | −0.038% | **−0.0377%** |
| r_floor (η=0) | −0.072% | **−0.0719%** |
| r_theta0 | −1.361% | **−1.361%** |
| r_bvp | −0.111% | **−0.111%** |
| rg_qflow2 | −0.072% | **−0.072%** |
| rg_theta0 | −1.427% | **−1.427%** |
| r5_qflow (η=0.5) | −0.014% | **−0.0143%** |
| r5_theta0 | −5.63% | **−5.634%** |

Also correct:

- **ω = 1.444635**, kernel `omega_core` at t=30 is **1.44464** (5 decimals)
- **VK branch** Q=88→210, ω=1.4679→1.4167 matches `vk_branch.txt`
- **Dressing** Q_θ/Q ≈ 1.57% (η=0.25), 6.16% (η=0.5) → rounded 1.6% / 6.2%; ratio ≈ η²
- **BVP removes ~92%** of drain (docs say ~95% — mild)
- **SFA semantics fix** `{1,3,2}→{0,1,2}` is real (`SFA_POSITION=0`, `ANGLE=1`, `VELOCITY=2`, `ACCEL=3`)
- **√3 off-shell factor**: Φ_a = f/√3 ⇒ s off by **27×**; Q ~ factor **3** — correct diagnosis

Older CONCEPT headline numbers also check out (spot-checked): ω window **1.3087**, A\* **0.464**, Q_max **921**, He **~1.56%** / Li **~3.5%**, flavor drift **−0.13%/1000 t.u.**, GDR **ω≈0.14**, force ≤6%, m_A bound.

---

## Issues (accuracy / self-consistency)

### 1. High — “maxF → 5e-15 in 14k iters” not in banked artifacts

`v72/FINDINGS.md` / `DISCOVERIES.md` claim:

> maxF 2.5e-2 → **5e-15** (machine floor) in **14k** iterations

Banked log `qflow_w145_eta025.log` only runs **8000** iters and ends at:

```
maxF = 6.4e-10   (iter 8000)
```

No other log supports 5.4e-15. Extrapolating the late-time ~22× per 2k-iter ratio gives ~**6e-14 at 14k**, not 5e-15. Convergence is clearly excellent; the **floor number and iter count are not evidenced** by what’s in `v72/results/`.

**Fix:** Quote the banked log (`6.4e-10 @ 8k`, or re-run to default 20k and bank that log), or drop the precise floor figure.

### 2. High — “(b),(c) fixed” overstates what was fixed

`DISCOVERIES.md` V72 item 1 ends with **“(b),(c) fixed.”**

| defect | status in code |
|---|---|
| **(c)** wrong SFA semantics in `eta_relaxer` | **Fixed** (unstaged diff) |
| **(b)** Φ_a = f/√3 seed | **Still present** in `eta_relaxer.c` and `radial_eta_soliton.c` |
| Correct per-component f | **Only in new** `eta_qflow.c` (`u[a]=f`) |

So the June-26 tooling bugs are **diagnosed** correctly, and the **new tool** is right, but claiming the old tools had (b) fixed is false. `FINDINGS.md` (“both tool bugs fixed (`eta_relaxer.c`)”) has the same problem.

**Fix:** “(c) fixed in `eta_relaxer`; (b) fixed in `eta_qflow` only; legacy tools still seed f/√3.”

### 3. Medium — CONCEPT η bullet slightly overclaims “gauged or ungauged”

New CONCEPT §3 text:

> drift … (η=0.25: **−0.04%**; η=0.5: **−0.014%** … **below** η=0’s −0.07%), **gauged or ungauged**

Data:

- Ungauged η=0.25: **−0.038%** (below floor) ✓
- Gauged η=0.25: **−0.072%** (= floor, **not** −0.04% and not “below”)
- η=0.5 gauged: **not tested**

“No intrinsic drain once stationary” is fair; attaching the **−0.04%** numbers to “gauged or ungauged” is not.

### 4. Medium — stale root era maps (miss v72)

| file | problem |
|---|---|
| `README.md` | History ends at **v70–v71**; no v72 |
| `CLAUDE.md` | “Current era … v71”; version dirs “v28/ … v71/” |
| `CONCEPT.md` §10–11 | Nuclear/substructure through v71; **§3 has η-soliton** but seed list omits `eta_qflow`; no v72 in infrastructure |

`EM_THEORY.md` banner is correctly marked historical — fine.

### 5. Low — rounding / soft wording

- BVP “~95%” vs **~92%**
- Dressing “1.6%” vs **1.57%** (acceptable)
- CONCEPT “sheds ~1.4% energy” for Θ=0 seed: **−1.36%** energy drift — OK as ≈, but it’s total E drift, not a pure radiation budget
- He “1.6% of rest mass **radiated**”: mass defect on E/Q is ~1.56%; merger ΔE is larger because charge is also shed — same looseness as `v71/NUCLEI.md`

### 6. Low — `Makefile` still incomplete for U(1) seeds

Unstaged change correctly adds `eta_qflow eta_relaxer radial_eta_soliton`, but `SEED_TOOLS` still omits all production U(1) generators (`gen_qball_*`, `gen_blob_field`, …) that live under `sfa/seed/` and in `bin/`. Pre-existing; the new patch doesn’t fix the broader install gap.

### 7. Low — FUTURE narrative vs RESOLVED

FUTURE’s June-26 “intrinsic η-drain floor” / η\*≈0.25 QFI section is still written in present tense; a later **RESOLVED (v72)** block correctly supersedes it. Chronological FUTURE style allows that, but a one-line banner at the top of the η-drain subsection (“superseded by v72”) would prevent skimming errors. Line ~849 (v67 θ-boundedness) is historical drain-with-massless-Θ — different claim; OK.

---

## What is solid

- Physics conclusion of v72 is supported: stationary η-seed drifts at or below the η=0 floor; Θ=0 is the drain; VK branch is monotone; ω matches the Lagrange multiplier.
- Root THEORY style of CONCEPT (cohesive, null results as definition) is intact; new bullet belongs in §3.
- DISCOVERIES/FUTURE chronological + RESOLVED structure is the right split for lab vs theory.
- Code change to SFA semantics is correct and necessary.
- Pre-v72 CONCEPT numbers that were spot-checked against v69–v71 match the measurement docs.

---

## Suggested doc fixes (priority)

1. Soften/replace **maxF 5e-15 / 14k** with banked log numbers (or re-run and bank).
2. Correct **“(b),(c) fixed”** → only semantics fixed in `eta_relaxer`; norm fixed in `eta_qflow`.
3. Split CONCEPT drift numbers: ungauged −0.04% / −0.014%; gauged η=0.25 matches floor −0.072%.
4. Bump **README / CLAUDE / CONCEPT §10–11** to include v72 and `eta_qflow`.
5. Optionally fix **f/√3** in `eta_relaxer` + `radial_eta_soliton` so the tools match the docs (or mark them deprecated).
