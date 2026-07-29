# v89 LEAK-LAW REGRESSION — sparse discovery over the existing mass corpus

*2026-07-28. Zero-spend: no simulator runs; every number below is parsed from
the standing MASS/H-campaign logs. Re-run end to end with
`python3 v89/prestress/regress/fit.py` (numpy-only; STLSQ and best-subset
hand-rolled). Raw printout: `fit_output.txt`. Table: `corpus.tsv`.*

---

## 1. Corpus

74 logs in the scratchpad + `v89/cellfab_runs/` (skimmed: double-slit era, no
census — excluded). Census-instrumented (`# LUMP`/`# CONV`) mass runs: 29
logs. Deduplication (verified by LUMP-series diff — the kernel is
deterministic, shorter runs are exact physics prefixes):

* `viz_ring6` ⊂ `m1_ring6` ⊂ `m2_ring6` ⊂ **`a1_ring6`** (one run, T=5000)
* `m1_ring12` ⊂ **`m2_ring12`**; `viz_blob` ⊂ **`m0_blob`**
* `v3_comp12|comp6|unwound12` ⊂ **`v3_*_L`** (T=5000)
* `h1_bubble` (v1) LUMP-identical to **`h1_shell0`** (the empty-hole accident)
* `h1_pocket`, `h1_pocket_v2`: no dense object ever (census empty) — excluded.

**Unique fit rows: 23** — 3 blobs (m0_blob, m1_ctrl, h1_solid), 8 ring_auto
(a1, a1c closed-box, a2_s1..s5, m2_ring12), 2 ring_naive (m1_ring6w,
b1_naive12), 3 chains (b1_comp12/comp6/unwound_match — the fleet-v2
self-rejected seeds: seam gate ~1e-22, open chains), 3 ring_exact
(v3_comp12_L/comp6_L/unwound12_L), 4 shells (h1_shell0, _v2, _v3,
h1_bubble_v2). **Uncensored deaths: 9** (+1 blob death, h1_solid, excluded
from parametric fits for lack of a voice count N). Everything else is
right-censored at its horizon. The honesty rule applies throughout: with 9
deaths and 20 structured rows, only 1–2-term laws are defensible.

## 2. Definitions (all in `fit.py` header; the load-bearing ones)

* **M_sum(t)** — sum of listed census lumps (the object; arcs of one ring are
  threshold artifacts, so the sum, not the top lump).
* **t_death** — last t at which the 21-sample rolling median of M_sum ≥ 0.15
  (just above census thr 0.1; the median kills the bath-flicker tail).
  Censored if within 30 t.u. of the horizon.
* **t_struct** — same with threshold 1.0 (one voice-mass): "structured-phase
  end"; robustness variable against the parked-crumb objection (§5.3).
* **x_voice(t) = M_sum/(N·cap)**, cap 2.5; **x50** = x_voice(50) (settled
  load; the t=0 seed surplus radiates in the first ~50 t.u.).
  **t_skirt** — first crossing of x_skirt = 0.0617.
* **leak_early** — 100·ΔM_sum/(Δt·M_sum(50)) over [50, min(300, t_death)], %/t.u.
* **c_eff = cap·(x50 − x_skirt)/(t_death − 50)** — whole-life absolute leak
  current per voice (upper bound for censored rows).
* **rough_share** — Δrough/ΔM_sum from the cumulative CONV ledger.
* Seed gates: shells print measured min/mean; ring gates computed from the B1
  algebra — g(ψ)=((1+cosψ)/2)^8 with φ̄=2π·closure/N; lock-recursion seeds:
  fwd ψ=0, back ψ=−2φ̄, seam carries the whole defect 2π·(closure−m);
  ring_m seeds: exact (fwd 1.0, back g(2·2πm/N)); naive wind=w: fwd ψ=−2πw/N,
  back ψ=−(2φ̄−2πw/N). Reproduces the design-note numbers: comp12 back
  0.1001, comp6 back 1.53e-5, naive6 all ≈0.100, chains' seam ≈0.
  *Deviation from the task brief:* naive12's predicted gates are ~0.57 fwd /
  0.44 back, not 0.100 — the "all gates 0.100" figure is the **N=6** naive
  value (ψ=60°); at N=12 the kick is 30°. The brief's number does not match
  the B1 algebra and the algebra is used.

## 3. Parser sanity vs MASS.md (all reproduced or explained)

Deaths reproduced within ~1–2%: a1c 1592/[1600], v3_unwound12_L 2201/[2221],
v3_comp6_L 3815/[3836], comp12_L censored@5000, h1_shell0 1217/[1222],
h1_shell0_v3 1811/[1814], h1_bubble_v2 1736/[1749]. Leaks reproduced on the
top-lump basis: ring6 −0.0577/[−0.058], ring12 −0.0459/[−0.051]; ctrl on the
Emfree basis −0.2275/[−0.232] (MASS mixed bases). Corrections to the record:

* **a1_ring6 "died ≈1900" is loose** — the census shows the last lump at
  t=1667 and sustained n=0 from 1672. Corpus value 1667.
* **h1_solid 1158 / h1_shell0_v2 237 were early calls** — census lumps
  persist to 1268 / ~300 (flicker tail). Corpus 1268 / 300.
* **m1_ring6w "−0.119 %/t.u." is NOT recoverable** from the census under any
  basis (top −0.063, sum −0.084, Emfree −0.071). The "naive winding leaks 2×"
  claim is **unverified** in this corpus; on consistent bases the naive-6
  penalty in [50,300] is 1.1–1.5×, and b1_naive12's whole-life current bound
  (≤4.1e-4) sits exactly at the universal value (§5.2).
* MASS.md's quoted per-voice loads for fleet v3 (0.21/0.32/0.44) do not
  reproduce; measured x50: unwound12 0.424, comp12 0.395, comp6 0.707. The
  corpus uses measured loads throughout.

## 4. Methods

Library: x50, margin=x50−x_skirt, N, 1/N, gmin, gmean, Dgate=Σ(1−g),
Bback=Σg_back, gmean·N, Dgate/N, margin·gmean (+ m, defect, wind on the
ring-only matrix). Fits: (i) exhaustive best-subset k≤2 (with n≈20,
exhaustive search over ≤2-term models is stricter and more honest than
tuned STLSQ), (ii) STLSQ (ridge λ=1e-2 + hard threshold sweep) as
specified. **CV = leave-one-topology-class-out** (blob/ring_auto/
ring_naive/chain/ring_exact/shell), pooled out-of-fold R².

## 5. Results

### 5.1 FIT 1 — early %-leak: NO transferable law (honest failure)

Best 2-term model CV_R² = **−0.03** (fit R² 0.27); every STLSQ configuration
CV-negative; ring-only matrix worse. Cause is visible in the a2 ensemble:
at *identical* topology/gates, foam seed alone moves leak_early from −0.039
to −0.152 %/t.u. (±60%) — the [50,300] %-rate is foam-chaos dominated.
**Verdict: under-determined; no %-rate law is claimed.** The transferable
quantity is the whole-life current (FIT 3).

### 5.2 FIT 2 + FIT 3 — the laws that survive

**Death law (the load line).** On the 9 uncensored structured deaths:

    t_death = 274 · (x50 / 0.0617)^1.066        [Law A]
    fit R² = 0.991, LOTO-CV R² = 0.969, σ_ln = 0.077 (×/1.08)

One term, three topology classes; nothing else survives (adding N, gmean or
gmin changes CV by <0.03 and their coefficients are noise). All ten shorter-horizon
censored rows are consistent. Physics reading: the intercept 274 t.u. is the
bare sub-skirt dissolution time (shell0_v2, which *starts* at the skirt,
dies at exactly ~300); the exponent ≈1 means the absolute current is
load-independent, which is FIT 3:

**The universal per-voice current.** c_eff = cap·(x50−x_skirt)/(t_death−50):

    c0 = 4.25e-4 Em/t.u. per voice   (median; spread 3.55–4.63e-4, MAD 3%)

across ring6 open AND closed box, a 5-seed ensemble member, one-way wound
comp6 (gmin 1.5e-5), mutual unwound12 (gmin 1.0), and three cube shells
(gmean 0.6) — i.e. **constant over 5 decades of seeded gate quality, N=6–12,
x50 0.14–0.71, and environment**. No descriptor explains the residual
(best one-term CV_R² ≈ 0). Equivalent linear form:

    t_death = 50 + cap·(x50 − 0.0617)/c0        [Law B]

Laws A and B agree within 5% over the observed range (see the scoring curve
in PREDICTIONS.md).

**The single exception — structure exists.** v3_comp12_L (exact m=5 wound
mutual ring) **beat Law A ×2.5** (alive at 5000 vs 1981 predicted;
c_eff ≤ 1.68e-4 ≤ 0.4·c0). It is the only violation among 20 rows.

### 5.3 The parked-crumb caveat, and why the exception survives it

Sub-voice remnant crumbs parked on favorable foam cells are common (a2_s1/s3/
s4, b1, m2_ring12 all end as flat 0.3–1.1 crumbs) and inflate t_death.
Robustness check on **t_struct** (last t with M_sum ≥ 1 voice-mass):
t_struct = 122·(x50/0.0617)^1.42, R²=0.844, and comp12_L's t_struct = 4879
still beats that law ×2.9. The comp12 excess is not a crumb artifact —
from t≈2500 to 5000 it holds a parked *arc* (rg 0.71, ~1.1 Em) with **no
secular decay over the final 2500 t.u.** — the nearest thing to an
equilibrium the program has produced. Caveat kept: single seed, and the
same foam cells host long-lived arcs in other N=12 runs at shorter horizons.

### 5.4 What the corpus says about the mechanisms (revisions to the record)

* **Roughness is the blob's killer, not the structured objects'.** CONV:
  m0_blob rough_share 45–52% ✓ (M0's number), but **every ring/chain/shell
  has rough ≈ 0** (worst: comp6_L 8.9%). The structured leak exits as dense
  trickle into bath/sink (E_transfer accounting closes: a1's lost 4.33 of
  5.16 is all sink-recorded dense), NOT as roughness radiation. The H1-v3
  verdict sentence "the leak leaves as field" does not hold for the
  shells' own ledger.
* **Seeded gate quality does not set the current** (c0 flat from gmin
  1.5e-5 to 1.0). The consonant-skin program (gates → no leak) is targeting
  a mechanism the rate-level data does not show; what gates/closure DID do
  in the record (chains vs rings, naive penalties) shows up in *early
  transient* rates, which FIT 1 shows are foam-dominated and not lawful.
* **Winding+mutuality is the only lever found beyond load** (comp12; the
  one-way wound comp6 sits exactly ON the load line, so winding alone or
  closure alone is not it).
* **Closed box = open box** for the current (a1c 4.33e-4 vs a1 4.13e-4) —
  confirms MASS's starvation ruling at rate level.
* t_skirt ≈ 0.6–0.8·t_death everywhere: the N-normalized load crossing the
  skirt precedes census death by ~a third of the life (consolidation keeps
  surviving cells above the local skirt while the object shrinks).

## 6. Verdict

* dM/dt = f(descriptors): **partial law only** — dM/dt ≈ −c0·N over a whole
  life (c0 = 4.25e-4 ±13%), descriptor-free; the %-rate refinement is
  under-determined (FIT 1 CV < 0) and is not claimed.
* t_death = g(descriptors): **t_death = 274·(x50/0.0617)^1.066** (CV R²
  0.97) with the pre-registered structural exception class "wound mutual
  exact ring" (≥2.4× above the line).
* One-term honesty beat five-term fits exactly as briefed: every multi-term
  model lost to the load line under class-held-out CV.

Files: `corpus.tsv` (34 rows × 43 cols, all definitions above),
`fit.py` (parser + fits, re-runnable), `fit_output.txt` (full printout),
`PREDICTIONS.md` (Phase-2 pre-registration).
