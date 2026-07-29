# INT_LEDGER — tentative integer-based simulator design

**Status:** implemented as `cellfabi.c` (2026-07-29). Production remains
`cellfab.c` (FP64). Full battery **20/20** on modes 0, 1, 3 after poly
trig + field-snap fix — see `int_ledger/RESULTS.md`. MASS waves were
blocked during numerics by user request; unblocking is a separate call.

Subordinate to `PRINCIPLE.md`. Motivation: energy is never destroyed;
atoms at mode boundaries; exact-mass (M-R1) wants a package that cannot
numerically evaporate. See also MASS §5c, EMF §5.

---

## 1. Goal

Replace **energy accounts** with **exact integer ledgers** so that:

1. Global conservation is **identically zero drift** (not ~1e−15).
2. Mode-boundary conversions move **integer multiples of the action atom**
   (in energy units derived from ε(ω)).
3. FP is retained only where the model is inherently non-linear
   (gates, pitch, geometry factors, unitary field kinematics).

Non-goals (v1): full fixed-point phases; FP16/BF16 state; dual-int16
"fused scale+value" packs; consumer-GPU FP64 parity.

---

## 2. Why this fits v89

| Law / practice | Integer reading |
|---|---|
| Paired ledger moves | `+Δ` / `−Δ` on two `int64` accounts — exact |
| Atoms at boundaries | conversion quantum `n_ε = round(ε/u)` or exact `ε/u` if commensurate |
| Transport continuous | sub-atom moves accumulate in residual or use unit `u ≪ ε` |
| Quant credit (mode 2) | credit counters already discrete-friendly |
| Exact mass package | total package = sum of integer accounts (+ defined flight) |

---

## 3. Unit system

Pick a single energy quantum **u** (joule analog) for the whole run:

| choice | u | E_tot~3×10⁴ in units | notes |
|---|---|---|---|
| coarse | ~1e−9 | ~3×10¹³ | fits int64 easily; atoms ε~0.5 → n_ε~5×10⁸ |
| preferred band | **1e−12 … 1e−10** | ~3×10¹⁶ … 3×10¹⁴ | headroom for local accumulators |
| fine | 1e−15 | ~3×10¹⁹ | may need int128 for **global** sums |

**v1 recommendation:** `u = 1e-12` (config key `ledger_u`), with
`int64` per-cell/per-link accounts and **hierarchical int64 or int128**
only for diagnostic global sums if needed.

**Commensurability with atoms:** prefer `u` such that
`n_ε(ω) = floor(ε(ω)/u + 1/2)` is stable across the pitch window, or
store atom sizes as `int64` per voice from `ε(ω_i)/u` at fire time
(source-sized V2 rule — already standing).

---

## 4. State partition

### 4.1 Integer (ledger) — must be exact

| account | meaning | type |
|---|---|---|
| `Es[i]` | space energy | int64 |
| `Em[i]` | dense energy | int64 |
| `Ee[i]` | field energy (mode total, or per-component split) | int64 |
| `lem[l,chan,dir]` | in-flight energy on link | int64 |
| quant credit pools | two-atom credit | int64 |
| conversion fire counters | diagnostics | int64 |

**Invariant (checked every diag):**  
`sum(Es+Em+Ee) + sum(lem) + sinks = E0` **exactly** (integer equality).

### 4.2 Floating (kinematic) — FP64 v1, optional FP32 later

| quantity | role |
|---|---|
| phases θ, field plane angles | gates, hops |
| ω_eff, x_load | pitch law (x from integer Em via `Em*u/cap`) |
| gate values, res(Δω), roughness R | rate factors ∈ [0,1] |
| geometry: d, A, normals | conductance, retardation |
| field amplitudes (fa,fb,…) | may stay FP with energy scaling bridged to Ee ledger |

**Bridge rule:** any rate that moves energy computes a real desire
`Δ_real`, then  
`Δ_n = floor(Δ_real/u + residual)` with **residual carry** per edge
(or stochastic rounding — prefer deterministic residual for ratchet
byte-stability).

### 4.3 Not in ledger

Scaffold positions `cx,cy,cz` (build + viz only). Snapshots may quantize
to float32 for I/O.

---

## 5. Update algorithm (one step sketch)

Keep the existing pass structure of `cellfab.c` (1…F, D plast, …).
For every energy-affecting event:

1. **Compute** kinematic factors in FP (gates, geo, ω, want).
2. **Propose** real energy move `Δ_real` (same formulas as today).
3. **Quantize** to `Δ_n` with residual:
   - `raw = Δ_real/u + resid[edge]`
   - `Δ_n = (int64)floor(raw)` (or trunc toward 0 with documented policy)
   - `resid[edge] = raw - Δ_n`
4. **Apply paired:** `A += Δ_n; B -= Δ_n` (or fire `k * n_ε` for conversions).
5. **Conversion door:** if quant_mode requires whole atoms, only fire when
   credit/pool allows `k≥1` atoms; never fractional atom across modes.
6. **Floors:** `es_floor`, `mob_floor` expressed as integer thresholds
   `floor(floor_real/u)`.

Plasticity (`ld` real lengths) stays FP; it does not hold energy.

---

## 6. Field sector options (v1 vs v2)

| option | idea | pros | cons |
|---|---|---|---|
| **F1 (v1)** | Ee ledger integer; amplitudes FP; scale \|ψ\| so E matches ledger after hops | smaller change | hop unitarity approx at u-scale |
| **F2** | fixed-point amplitudes with shared energy norm | closer to exact | hard unitary |
| **F3** | field fully FP64 until deposit/extract at conversion door | simplest | flight continuous FP |

**v1 pick: F1 + F3 hybrid** — continuous field hops in FP64 with
periodic or end-of-hop **pull to integer Ee/lem** via residual, so the
conversion door and global sum stay exact.

---

## 7. Performance sketch (CPU / GPU)

### CPU (OpenMP graph, current host)

* int64 add is cheap; **not** the bottleneck vs gather/scatter + trig.
* Expected win: modest, unless it enables coarser FP elsewhere.
* Dual int32 software big-int: **worse** than native int64 — reject.
* int16 / dual-int16: **reject** for energy (range + atom resolution).

### GPU (future port only)

* Consumer NVIDIA: weak FP64 → integer ledgers **more** attractive.
* Datacenter: FP64 OK; integer still best for exact conservation.
* Tensor cores irrelevant to pair ledgers on foam graphs.
* atomics on int64 accounts need careful colorings (already edge-color
  hops on CPU) — reuse determinism story.

---

## 8. Implementation plan (when opened)

| step | deliverable | gate |
|---|---|---|
| 0 | This design + config `ledger_u`, `ledger=0|1` | doc only |
| 1 | Shadow mode: FP64 dynamics, **integer shadow accounts** updated in parallel; report max \|E_fp/u - E_int\| | no physics change |
| 2 | Integer apply path for dense+space only; field still FP bridge | e1 exact int drift 0; battery |
| 3 | Flight lem integer | same |
| 4 | Conversion atoms exact n_ε | ħ-linearity; LIN; quant tests |
| 5 | Optional FP32 kinematics | battery + bit-stability policy |
| 6 | GPU port uses ledger path first on consumer parts | later |

**Equivalence policy:** ledger=1 need not be bit-identical to FP64; it
must pass the **same physical bars** (battery + MASS C1/P19 when claimed).
Shadow mode (step 1) quantifies FP64's hidden numerical leak.

---

## 9. Risks

* Residual carry can bias slow secular flows — monitor vs FP64 shadow.
* ε(ω)/u not integer → rounding policy must be frozen and tested.
* Thread reductions: global int sum order must be deterministic
  (mirror current serial diag / fixed gather order).
* Plastic foam rewrite + integer energy: orthogonal but both touch
  determinism; enable separately.

---

## 10. Relation to MASS / EMF

* **MASS exact mass:** integer package makes C1/M-R1 failures physical,
  not numerical; P19 clustering is cleaner.
* **EMF:** clicks and atoms already discrete; integer conversion door
  clarifies P2 bookkeeping when joint exams run.
* **Does not replace** hardening pin, multi-seed, or PRESTRESS waves.

---

## 11. Log

* **2026-07-29** — tentative design opened from numerics discussion
  (FP64 necessity; reject dual-int16 speed path; prefer int64 energy +
  FP kinematics). No code path yet.
