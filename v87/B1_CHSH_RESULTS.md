# v87 B1 — the CHSH test: protocol, analytic layer, and the in-kernel run
**Date:** 2026-07-26 · Follows `CRANK2_RESULTS.md` (as corrected by
`crank2_review.md`) and the three-seat rung. Instruments:
`work/geom/bell_grid.c`, `work/kernel/bell_kernel.c`.

Tags: **[D]** derived · **[M]** measured · **[P]** postulated · **[C]** conjecture

---

## 0. What was asked for

CRANK2 §6 closed with one falsifiable and cheap recommendation, unchanged from
the seats':

> Run the in-kernel CHSH with two objects from a common origin and a
> phase-overlap readout. Predictions: S ≤ 2 for an unbiased readout (the fabric
> is then a Cauchy-slice local model and Bell applies); S → 2√2 *only* if the
> source ensemble itself is setting-correlated; and if S > 2 appears with a
> detector-side bias, that is the **closed detection loophole** and must be
> rejected, not celebrated.

This rung builds that test in two layers and runs both.

* **Layer 1 — the analytic construction** (`bell_grid`). λ from an RNG. This is
  where the p-family ladder of CRANK2 §1–§3 lives, and where the *measurement
  protocol* is calibrated before it is pointed at the fabric.
* **Layer 2 — the fabric** (`bell_kernel`). λ from the PDE: the measured
  internal clock phases of two gauged Q-balls built from one common
  construction and evolved by an unmodified `sfa/sim/scp_sim.c`.

Layer 1 is not a warm-up. It produces the null band without which Layer 2's
number cannot be read at all (§2.2).

---

## 1. The three protocol requirements

A CHSH number is `max |S|` over a grid of settings. That is a **maximum of a
noisy statistic**, and three things can manufacture a violation that is not
there. Each gets a named control, and both instruments implement all three.

| # | requirement | what it defends against |
|---|---|---|
| 1 | **Analytic control** | The search is run on an exact `E(θ)` whose `max|S|` is known in closed form. Bounds what the search machinery alone can produce from a noiseless table — i.e. separates grid discretisation from statistics. |
| 2 | **Search-bias null** | `max` over settings is biased **up** even when the true `S` equals the bound exactly. Measured directly by replicating ensembles whose true `S` is known, then quoting `c = excess·√N`. Nothing counts as exceeding the bound unless it exceeds the null band. |
| 3 | **Discard fraction** | Report η = mean readout weight and the fraction thrown away, against the CHSH threshold 2(√2−1) = 0.828427 and the Eberhard/CH threshold 2/3. A violation bought by discarding is the detection loophole, and must be labelled as such rather than celebrated. |

Requirement 2 also needs a **pre-registered** setting quadruple — settings
fixed before looking at the data, evaluated with no maximisation, hence
unbiased. The obvious choice, "take the analytic argmax", turns out to be a
trap and is documented as one in §2.2.

---

## 2. Layer 1 — the analytic construction  **[D,M]**

`bell_grid`, exact-integer streaming, N = 10¹⁰ samples, grid = 65536 settings
(0.00549°), 262144 histogram bins and correlation-table angles. Arm 1 is
uniform λ (respects (MI)); arm 2 carries the p = 1 tilt ρ₁(λ|b) ∝ |λ·b|.
Full log: `work/geom/big_run2.log`.

### 2.1 Analytic control — passed exactly

| arm | exact target | search on the exact table | excess |
|---|---|---|---|
| 1 (triangle) | 2 | 2.000000000000 | **+4.34×10⁻¹⁹** |
| 2 (cosine) | 2√2 = 2.828427124746 | 2.828427124746 | **0** |

The search adds nothing to a noiseless table at this grid resolution. Any
excess in the sampled runs is therefore statistics, not discretisation. **[M]**

### 2.2 Search-bias null — the finding, and it is not a formality

12 independent replicas per arm at N = 10⁹, quoting `c = (max|S| − truth)·√N`:

| arm | searched `c` | range | fixed-settings `c` |
|---|---|---|---|
| 1 (triangle) | **+3.954 ± 1.076** | [+1.929, +5.919] — **12/12 positive** | +0.440 ± 1.650 |
| 2 (cosine) | +0.604 ± 1.869 | [−2.718, +3.656] | +0.481 ± 1.844 |

At production N = 10¹⁰ the same thing shows up directly as
`max|S| − |S|_pre-registered`:

| arm | searched max\|S\| | pre-registered \|S\| | **search bias** |
|---|---|---|---|
| 1 | 2.000047949200 (+4.79 σ) | 1.999975954400 (−2.40 σ) | **+7.199×10⁻⁵** |
| 2 | 2.828412372683 (−1.48 σ) | 2.828412054160 (−1.51 σ) | **+3.185×10⁻⁷** |

**The search bias is a property of the waveform, not of the search.** A factor
of 226 separates the two arms at production N, and a factor of ~29 separates
them in the replica means. The reason is degeneracy:

* The **cosine** attains 2√2 at an *isolated* setting quadruple. The search
  cannot wander off it, so `searched ≈ fixed` and the bias is statistically
  zero. **[D,M]**
* The **triangle** attains 2 on a *plateau*. The run's analytic argmax landed
  at (a′, b, b′) = (0.0165°, 180°, 45.011°) with value exactly
  2.000000000000, while the canonical (90°, 225°, 135°) gives exactly
  2.000000000000 too. A maximum taken over a plateau of independently noisy
  entries is biased up by construction. **[D,M]**

This is load-bearing for Layer 2, because **the waveform the fabric is expected
to produce is the triangle** — precisely the arm where a naive `max|S|` drifts
above 2 for purely statistical reasons. Without this null, a spurious
"violation" in the kernel run would have been the default outcome.

**Verdict at production N:** arm 1's measured excess +4.795×10⁻⁵ sits *inside*
its own null band (< 7.182×10⁻⁵). Arm 1 is consistent with S = 2 exactly, as
it must be. **[M]**

**The pre-registration trap.** Reading the pre-registered settings off the
analytic argmax is wrong for the triangle: the search returns the degenerate
point a′ = b′ = 0, where S = 2·E(π) with the two noisy terms cancelling
identically. That estimator has *zero variance* and tests nothing. Both
instruments therefore pin the canonical CHSH quadruple
(a, a′, b, b′) = (0°, 90°, 225°, 135°), which attains the maximum for **both**
arms with four distinct, independently noisy entries. **[D]**

### 2.3 Discard fraction — Reading 2 stays refuted, with the number attached

| arm | η = ⟨w⟩ | discarded | vs 2(√2−1) = 0.828427 | vs Eberhard/CH 2/3 |
|---|---|---|---|---|
| 1 | 1.000000 | 0 | n/a — unweighted | n/a |
| 2 | **0.636618** | 0.363382 | **BELOW** | **BELOW** |

η measured at 0.636618 against the analytic 2/π = 0.636620. **[M]** Arm 2
reaches Tsirelson only at an efficiency below *both* thresholds, so read
detector-side it is inside the detection loophole — CRANK2 §4 Reading 2,
unchanged and now with the instrument reporting the verdict line itself rather
than leaving it to the reader.

### 2.4 Why arm 1 can exceed 2 at all — a caveat on this instrument  **[D]**

CHSH is supposed to be unbreakable by a local model, so arm 1 sitting
+4.79 σ above 2 deserves an explanation rather than a shrug. It is an artifact
of `bell_grid`'s **estimator**, not of the physics.

`bell_grid` builds a 1D table `E(θ)` from a single λ histogram and reads all
four CHSH terms out of it by shifting. That is only valid if the λ histogram is
exactly translation invariant — the four terms need λ reparametrised by *b* in
two of them and by *b′* in the other two, and one histogram cannot be both. At
finite N the histogram is non-uniform at O(1/√N), the per-bin bracket is no
longer confined to ±2, and the maximum over settings picks up the resulting
slack. Hence +3.95/√N, always positive.

This is a deliberate speed trade — the 1D table is what makes N = 10¹⁰ at
65536 settings affordable — but it means **`bell_grid`'s arm 1 is not a strict
CHSH estimator**, and its null is not optional bookkeeping: it is the only
thing that makes the number readable. `bell_kernel`'s 2D table has no such
defect (§3.2).

---

## 3. Layer 2 — the in-kernel test

Two runs, both `scp_sim.c` unmodified, absorbing BC, η = 0, g = 0.05, D = 13:

| run | N | L | dx | T | snap_dt | frames | role |
|---|---|---|---|---|---|---|---|
| `b1_val_N64` | 64 | 19 | 0.60317 | 40 | 0.25 | 167 | §3.2 result, coarse point of the dx pair |
| `b1_prod_N128` | 128 | 19 | **0.29921** | 150 | 0.5 | ~300 | 2× refinement, 3.6× baseline (§3.4) |

### 3.1 Protocol, and what is honestly free

Source: two gauged Q-balls from the v86 N7 profile (ω = 1.430, g = 0.05,
Q = 221 each), stamped by `gen_qball_pair` from one common construction at
internal phases δ = 0 and δ = π, centres at x = ∓6.5 (D = 13), evolved by an
**unmodified** `scp_sim.c` at η = 0, absorbing BC.

Readout, mirroring the phase construction but with *measured* phases:

```
    phi_A = arg  sum_{|r-cA|<R} sum_a Phi_a        A(a) = sgn cos(phi_A - a)
    phi_B = arg  sum_{|r-cB|<R} sum_a Phi_a        B(b) = sgn cos(phi_B - b)
```

**Gauge.** `arg Φ` is not gauge invariant and the wings are far apart, so the
naive difference is not either. The invariant object is
Δ_gi = φ_A − φ_B − Σ_links θ_x along the straight x-path joining the centres.
Both are computed; the invariant one governs.

**What is free.** The settings a, b are chosen by the analyst *after* the run
is on disk. That is the strongest available form of fresh entropy at the
choosers (GEOM Geo10.1): the simulation cannot have correlated its output with
settings that did not exist while it ran.

**So this is a tripwire, not a discovery.** (MI) holds by construction, the
model is local and deterministic, Bell's theorem applies, and |S| ≤ 2 is
guaranteed. It is worth running for three things it *can* deliver:

1. **A code check, and only that.** |S| > 2 would mean one outcome function had
   acquired a dependence on the *other* wing's setting — a bug in this analysis
   tool, not a statement about the kernel. §3.2 shows this is the whole of what
   the S-value can test, and why.
2. **The waveform.** Geo10.2's fingerprint: does the fabric's E(θ) track the
   classical triangle, or −cos θ? The two differ by +0.207 at θ = π/4.
3. **Transport.** The phase-lock visibility |R| of the pair clock, block by
   block in time. This is a *new* number, and it is the one CRANK2 §5.2 gap 3 /
   GROK G14 asks for: does a common-past correlation survive free evolution?

It does **not** address CRANK2 §5.4 gap 1 — there is still no setting that is a
fabric degree of freedom. The settings here are analyst-side labels applied to
recorded data.

### 3.2 The result that changes the experiment's status  **[D]**

Before any number: **the S-value of this experiment was decidable in advance,
and not because "Bell's theorem applies" in the loose sense.** It is an
algebraic identity of the readout.

A(a) and B(b) are ±1-valued and **each depends on its own setting only** —
A never sees b, B never sees a. (That, not spatial disjointness, is what the
argument needs; the gauge-invariant variant folds a Wilson line spanning both
wings into B, and the bound is untouched because the line is
setting-independent.) For **each individual frame**, writing α = A(a),
α′ = A(a′), β = B(b), β′ = B(b′), the CHSH bracket

```
    A(a)B(b) + A(a)B(b') + A(a')B(b) - A(a')B(b')
```

is a sum of four ±1 products algebraically confined to {−2, +2}. A sample mean
of numbers in [−2, 2] cannot leave [−2, 2]. So |S| ≤ 2 holds **exactly, at
every M, with no statistical slack whatsoever** — and it holds no matter what
the PDE did in between. Running the simulation cannot affect it.

The consequence is sharper than "the prediction was confirmed":

> **An in-kernel CHSH with an offline per-object dichotomic readout cannot be
> evidence about the fabric's locality, in either direction.** Its S-value is
> fixed by the shape of the readout before the kernel is invoked. CRANK2 §6
> recommended this test as "falsifiable and cheap"; it is cheap, and the S ≤ 2
> half is *not falsifiable* — no dynamics could have produced anything else.
> **[D]**

What the run therefore does measure is three other things, all real, and the
tool is scoped to them.

**Measured — validation run** (N=64, L=19, dx=0.60317, T=40, snap_dt=0.25,
M=167 frames, D=13, g=0.05, η=0; `work/kernel/results/val_N64.log`):

| quantity | value | reading |
|---|---|---|
| frames used / discarded | 167 / **0** | η = 1, no post-selection |
| Δ circular mean (gauge-inv.) | **179.959°** | the seeded π survives |
| Δ range across the run | 3.141403 – 3.141948 rad | spread **±0.016°** |
| phase-lock visibility \|R\| | **1.000000** | clocks perfectly locked |
| wing coherence cohA / cohB | 0.994315 – 1.000000 | "one clock per object" holds to 0.6% |
| Wilson line \|W\| | mean 7.6×10⁻⁴, max 2.4×10⁻³ rad | gauge correction ≈ 0.14°, negligible at g = 0.05 |
| fast phase φ_A vs uniform | D_KS = **0.0170** (crit₀.₀₅ = 0.1051) | λ equidistributes; semi-analytic arm licensed |
| analytic control | 2.000000000, excess +2.7×10⁻¹⁵ | search adds nothing |
| **empirical max\|S\|** | **2.000000000000** | at the bound, exactly |
| semi-analytic max\|S\| | 2.000000000 | agrees |
| search-bias null (12 reps) | **+4.44×10⁻¹⁶ ± 0** | zero, as the identity demands |

**Transport, block-resolved** — the number CRANK2 §5.2 gap 3 / GROK G14 asked
for. Over t = 0 → 41.5 at D = 13:

| block | t range | ⟨Δ⟩ deg | \|R\| | S_semi |
|---|---|---|---|---|
| 0 | 0.00 – 6.50 | 180.000 | 1.00000 | 2.000000 |
| 1 | 6.75 – 13.25 | 180.000 | 1.00000 | 2.000000 |
| 2 | 13.50 – 20.00 | 180.006 | 1.00000 | 2.000000 |
| 3 | 20.25 – 26.75 | 179.962 | 1.00000 | 2.000000 |
| 4 | 27.00 – 33.50 | 179.918 | 1.00000 | 2.000000 |
| 5 | 33.75 – 41.50 | 179.884 | 1.00000 | 2.000000 |

A common-past phase relation between two Q-balls at D = 13 **survives free
evolution with visibility 1.000000**, drifting −0.116° over 41.5 time units
(≈ −0.0028°/unit, i.e. ~9.4 clock periods with no measurable loss). **[M]**

This cuts *against* the pessimistic reading of G14 for this configuration: the
kernel's mixing does not scramble a relative clock phase on these timescales.
It does **not** rescue Door A, because a locked Δ is not a setting-correlated
ρ(λ|b) — it is one fixed number, carrying zero bits about anything chooseable.

### 3.3 Arm 2 — the p = 1 tilt on the fabric's own phases  **[M]**

Applying w = |cos(φ_B − b)| — the circle form of ρ₁(λ|b) ∝ |λ·b| — to the
measured fabric phases:

| quantity | value |
|---|---|
| max\|S\| | **2.836696** |
| Tsirelson | 2.828427 |
| excess | +0.008269 |
| **arm 2 null spread** (12 reps) | +0.051755 ± 0.025537, band **+0.128365** |
| verdict | **consistent with exactly 2√2** — the excess is noise |
| η = ⟨w⟩ | **0.636620** (analytic 2/π = 0.636620) |
| discarded | **0.363380** |

Two things this pins down, both of which needed the null to state honestly:

* The fabric's own phases **do** reach Tsirelson under the p = 1 tilt — the
  construction is not an artifact of synthetic λ. **[M]**
* They reach it **only** by discarding 36.3% of pairs, at η = 0.6366, below
  *both* 2(√2−1) = 0.8284 and Eberhard/CH 2/3. That is CRANK2 §4 Reading 2,
  the closed detection loophole, now demonstrated on fabric data rather than
  argued from a mean. **[M]**

Without the arm-2 null the raw 2.836696 > 2.828427 would have read as
"exceeds Tsirelson". It does not: unlike arm 1, arm 2 has genuine statistical
slack because the b-dependent denominators break the pointwise argument, and
the measured excess is 6.4% of the null band.

### 3.4 Resolution and duration — production run  **[M]**

N=128, L=19, **dx = 0.29921** (2× finer), T=150 (3.6× longer), M=305 frames,
0 discarded. `work/kernel/results/prod_N128.log`.

| quantity | validation (dx 0.603, T 40) | **production (dx 0.299, T 150)** |
|---|---|---|
| empirical max\|S\| | 2.000000000000 | **2.000000000000** |
| semi-analytic max\|S\| | 2.000000000 | 2.000000000 |
| arm-1 null bias | +4.44×10⁻¹⁶ | +4.44×10⁻¹⁶ |
| discarded frames | 0 | **0** |
| φ_A vs uniform D_KS | 0.0170 (crit 0.1051) | 0.0200 (crit 0.0778) |
| Δ circular mean | 179.959° | 179.700° |
| \|R\| overall | 1.000000 | **0.999849** |
| arm 2 max\|S\| | 2.836696 (+0.008269) | 2.829767 (**+0.001340**) |
| arm 2 null band | +0.128365 | +0.072723 |
| arm 2 η | 0.636620 | 0.636620 |

The tripwire is clear at both resolutions and the identity of §3.2 holds
exactly, as it must. Arm 2's excess over 2√2 shrinks from +0.0083 to +0.0013 as
M grows — consistent with pure sampling noise, which is what the null said.

**The transport number sharpens, and it is not linear.** Block-resolved over
t = 0 → 152:

| block | t range | ⟨Δ⟩ deg | \|R\| |
|---|---|---|---|
| 0 | 0.0 – 24.5 | 179.996 | 1.00000 |
| 1 | 25.0 – 49.5 | 179.933 | 1.00000 |
| 2 | 50.0 – 74.5 | 179.895 | 1.00000 |
| 3 | 75.0 – 99.5 | 179.797 | 0.99998 |
| 4 | 100.0 – 124.5 | 179.475 | 0.99980 |
| 5 | 125.0 – 152.0 | 179.156 | **0.99944** |

The pair clock **does** stay locked — |R| = 0.99944 after 35 clock periods —
but the drift **accelerates**: −0.0013°/t.u. in block 1 rising to −0.012°/t.u.
in block 5, roughly an order of magnitude. So the honest statement about
CRANK2 gap 3 is not "the correlation survives" full stop; it is *the phase
relation survives with visibility 0.9994 over the interval probed, degrading
super-linearly.* Whether that is physical decoherence or accumulated
integrator/sponge error is **not settled here** — separating them needs a
dt-refinement pair, which this rung did not run. Quoted as measured, not
interpreted.

Note the coarse run's drift rate (−0.0028°/t.u. over t≤41.5) is *larger* than
the fine run's over the same window (−0.0013°/t.u.), i.e. the two resolutions
bracket rather than contradict — the coarse grid overstates early drift, and
the late acceleration is only visible on the longer baseline.

---

## 4. Status after B1

| question | status |
|---|---|
| Is the CHSH protocol calibrated (control, null, discards)? | **Yes** [D,M] — both instruments, all three requirements |
| Does the search inflate max\|S\|? | **Depends on the estimator, and on the waveform** [D,M] — 1D translation-invariant table: +3.95/√N, 12/12 positive on the triangle, ~0 on the cosine. 2D consistent-sample table: exactly 0 |
| Does the fabric violate Bell in-kernel? | **The question is malformed as posed** [D] — with an offline per-object readout, \|S\| ≤ 2 is an algebraic identity of the readout, not a result about the fabric |
| Does the fabric's pair clock stay locked under transport? | **Yes, with a caveat** [M] — \|R\| = 0.99944 after 35 clock periods at D = 13, but the drift accelerates ~10× from block 1 to block 5. Physical decoherence vs integrator/sponge error is NOT separated (needs a dt pair) |
| Does the fabric reach 2√2 under the p = 1 tilt? | **Yes, and only by discarding 36.3%** [M] — η = 0.6366, inside the detection loophole |
| Is Tsirelson derived? | **No** — unchanged |

### 4.1 What B1 changes in the CRANK2 gap table

CRANK2 §6 listed eight gaps and recommended the in-kernel CHSH as the cheap
falsifiable next step. B1's finding re-orders them:

* **Gap 1 (what is a setting in the kernel?) is not an independent gap — it is
  a prerequisite for the recommended experiment.** An in-kernel CHSH is a
  genuine locality test only if the settings are *dynamical*: fields present in
  the initial data, at each wing, that the outcome depends on through the PDE.
  Then A *could* in principle depend on b, and S ≤ 2 becomes a real statement
  about finite propagation speed instead of an identity. With offline settings
  it is neither falsifiable nor informative about the fabric. **[D]**
* **Gap 3 (survival under transport) is measured and, for this configuration,
  open in the favourable direction.** |R| = 0.99944 over 35 clock periods, with
  a super-linearly accelerating drift whose origin is not yet separated from
  numerics. This does not close it —
  what must survive for Door A is a setting-*correlated* tilt, not a fixed
  relative phase — but the kernel's mixing is not the obstruction here that
  G14's pessimism suggested. **[M]**
* **Gap 4 (continuous → dichotomic without a detection loophole) is now
  quantified on fabric data.** The unbiased map costs nothing (η = 1) and gives
  2; the p = 1 map gives 2√2 and costs 36.3%. **[M]**

### 4.2 The experiment that would actually be a locality test  (B2, designed, not run)

Requirements, all satisfiable without kernel modification:

1. A **dynamical analyzer** at each wing — an additional fabric object whose
   internal phase encodes the setting, present in the initial data.
2. The outcome read from the *interaction* at that wing (A6 contact:
   co-phase fuses, anti-phase repels ⇒ a genuine dichotomy), not from a
   phase the analyst post-processes.
3. Wing separation kept spacelike for the readout window: with the measured
   contact scale κ ≈ 0.5117 and D = 13, direct coupling is e^{−κD} =
   1.3×10⁻³ (CRANK2 §5.2 quotes 2.2×10⁻³ at D = 12), so the wings are
   dynamically isolated over the relevant interval.
4. One run **per setting pair** — the settings are initial data, so they
   cannot be swept in post-processing. A 4×4 grid is 16 runs.

Cost at the resolution used here: ~2.5 h per run local CPU, so ~40 h for a
4×4 grid, or a few hours on a GPU. Prediction: |S| ≤ 2 still, but now as a
*measurement* of the kernel's finite propagation speed — and a violation would
be a genuine tripwire of the same class as a Gauss-law drift, which the present
test cannot be.

---

## 5. Reproduction

```
# Layer 1 — analytic construction, production + null
cd v87/work/geom
gcc -O3 -march=native -fopenmp -o bell_grid bell_grid.c -lm
./bell_grid 10000000000 65536 262144 262144 --null 12 1000000000

# Layer 2 — seed, run, analyse
./bin/gen_qball_pair 128 19.0 \
    v86/n_battery/n7_profile_w1430_g005.txt 1.430 0.0  -6.5 0 0 \
    v86/n_battery/n7_profile_w1430_g005.txt 1.430 3.14159265 6.5 0 0 \
    v87/work/kernel/seeds/pair_D13_prod.sfa
# sim: N=128 L=19 T=150 snap_dt=0.5, complex_gauge=1 g=0.05 eta=0, bc_type=0
cd v87/work/kernel
gcc -O3 -march=native -fopenmp -o bell_kernel bell_kernel.c -lzstd -lm
./bell_kernel /space/scp/v87/b1_prod_N128.sfa --radius 5.0 --grid 360 \
              --snapdt 0.5 --null 20 --out b1_prod
```

Note the geometry convention, which cost two discarded runs on this rung:
in both the kernel and the seed generators **`L` is the HALF-domain** and
`dx = 2L/(N−1)`. `N=128, L=19` is a box [−19, 19] at dx = 0.29921, not a box
of side 19.
