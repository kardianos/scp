# v64 FINDINGS — H1, the two-carrier hypothesis

**Date**: 2026-06-05
**Hypothesis**: [`PLAN.md`](PLAN.md) §1c (H1). The 6-field theory carries two forces
split by **carrier mass**: the **strong/nuclear force** = time-dilation lapse of the
*massive* φ-density sector (`m²=2.25`, short range); **gravity** = the *massless* sector
(θ-torsion `m_θ²=0`, and/or the v6 1/r density wave, long range).

This file records results as they land — positive and null — per the project's
lab-notebook convention.

---

## T-H1a — Is "strong force = time dilation of the massive φ sector" structural?

**Artifact**: `mass_range_check.py` (all self-checks pass).
**Question**: if the strong/nuclear force is the lapse well `α(P)` of the massive φ
sector, its spatial *range* must be the carrier Compton length `1/m`. Scan `m²`, measure
the range two independent ways, and check the scaling.

### POSITIVE result — the range is structurally carrier-set

| measure | log-log slope vs `m²` | expected (`R ∝ 1/m`) |
|---|---|---|
| Yukawa length `1/m` (definitional) | **−0.500** | −0.5 |
| convolved-response tail e-fold | **−0.485** | −0.5 |
| tail e-fold, source widths w = 0.1 / 0.5 / 1.0 | −0.500 / −0.485 / −0.438 | −0.5 (source-indep.) |

- The force **range is the carrier Compton length** — short range is *forced* by `m_φ`,
  not fitted — and it is **source-independent** (the soliton's own size does not set the
  tail decay; the carrier mass does).
- At standard `m²=2.25`: range ≈ **0.31–0.37 fm** — nuclear, not gravitational
  (cf. nuclear force ~1 fm, prior BLV ~0.8 fm).
- **Corollary (free)**: the massless sector (`m_θ²=0`) has *infinite* range — exactly
  the long-range carrier H1 assigns to gravity. The mass-split of H1 is **internally
  consistent**.

### NULL / methodological result — the "well half-width" is a bad range measure

- A first attempt measured the lapse-well *half-width* (radius where the response halves
  from its peak). It came back essentially **flat** (slope −0.05), i.e. mass-independent.
- This is **not** evidence against H1 — it is a broken metric. The screened response is
  `∝ 1/r` at the origin (divergent peak), so "half of the peak" is set by the
  mass-independent 1/r core, not the carrier. The physically meaningful range is the
  **tail decay length**, which is cleanly carrier-set (above). Recorded so the dead
  measure is not re-tried.

### What T-H1a does NOT establish (carried open)

- **Magnitude**: the strong-side coupling is ~O(1) `Φ/c²` (correct for nuclear); the
  *gravity*-side (massless) candidate is still ~10⁴⁰× too strong (v6). Hierarchy
  relocated, not solved.
- **Tensor**: gravity still needs h=±2; the massless scalar/density tail is spin-0.
- T-H1a is kinematics (range), not dynamics — it shows the *split is consistent*, not
  that the real 6-field runs produce it. That is T-H1b.

---

## T-H1b — does the long-range tail live ONLY in the massless sector?

**Stage 1 — analytic precheck DONE** (`th1b_sector_range.py`, all self-checks pass).
Per the guardrail (prechecks before GPU), this establishes the falsifiable law and the
measurement protocol the simulation must follow. Linear screened response of
`(-∇² + m²)φ = ρ` on a radial grid; far-field window r ∈ [4, 8] code lengths.

### POSITIVE result — the tail requires masslessness AND a monopole source

**(A) Carrier-mass scan** (monopole source), far-field falloff power `p = d ln|φ|/d ln r`:

| `m²` | `m` | power `p` | regime |
|---|---|---|---|
| 0.00 | 0.000 | **−1.000** | 1/r long-range (gravity carrier) |
| 0.04 | 0.200 | −2.148 | near-massless |
| 0.25 | 0.500 | −3.870 | screened |
| 1.00 | 1.000 | −6.739 | screened |
| **2.25** | 1.500 | **−9.609** | Yukawa-screened (strong/nuclear) |

The massless carrier gives a clean 1/r tail (`p=−1.000`); the standard massive φ sector
(`m²=2.25`) is steeply screened (`p=−9.6`). **Massless ⇒ long-range, massive ⇒
short-range** — the H1 mass-split, confirmed at the linear-response level.

**(B) Source-monopole structure** (massless carrier, the v6 Gauss bypass):

| source | monopole `Q` | far-field `\|φ(6)\|` |
|---|---|---|
| Gaussian (v6 `\|ω\|²` channel) | `+0.096` | `1.6×10⁻²` (real 1/r tail) |
| `∇²`(Gaussian) (total divergence) | `≈0` (`−2×10⁻¹⁶`) | `5.6×10⁻¹⁷` (no tail; ratio `3.5×10⁻¹⁵`) |

A 1/r tail requires **both** masslessness **and** a nonzero-monopole source. A
total-divergence source (zero monopole) produces **no** 1/r tail even when massless —
this is exactly the v6 Gauss obstruction (the algebraic channel falls as 1/r⁶, not
gravity). The tail's *coefficient is the monopole*.

### NULL / methodological result — "power" is the wrong tail diagnostic for Q=0

- A zero-monopole source's true far-field response is exponentially small; a naive
  power-law refit on it returns a spurious `p=−1.000` (numerical roundoff `~10⁻¹⁶`,
  divided by r, looks like 1/r). Existence of a 1/r tail is an **amplitude** question
  (the monopole coefficient), not a power-law one. Recorded so it isn't re-tried; the
  test now compares far-field amplitude (A/B `< 10⁻³`).

### What the precheck does NOT establish (the real T-H1b, needs the sim)

- It shows the *kinematics* permit the split. It does **not** show the *real nonlinear
  θ-dynamics* generate a nonzero-monopole source. v6 claims they do (`Q_eff = E₂/ρ₀`),
  but that is the load-bearing dynamical fact and must be measured in a run.
- **Magnitude unchanged**: even with the 1/r tail confirmed, the v6 coupling is
  `G_eff/G_N ≈ 3.6×10⁴⁰` — still ~10⁴⁰× too strong. **Tensor unchanged**: the tail is a
  scalar monopole, not h=±2.

### Stage 2 — the GPU run DONE (V100, 2026-06-05)

**Setup**: 192³ neutron template stamped into a ±10 (code-length) box (kernel resampled
to N=128, dx=0.157), full 6-field kernel, standard params, absorbing BC (`damp_width=3`),
`T=60`. Two runs: **massless** `m_θ=0` and **massive-θ control** `m_θ=1.5`. Radial
profiles via `shell_analysis` (φ-rms, θ-rms, θ/φ in 50 shells); fit window r∈[4,7]
(outside the extended soliton body, inside the absorbing shell). Data in
`/space/scp/v64_h1b/` (`neutron_massless.sfa`, `neutron_massive.sfa`, `diag_*.tsv`).
Provisioning required fixing the runner first (see `RUNNER_FIX` note below).

#### POSITIVE — shape segregation is in the H1 direction
The θ/φ ratio vs r separates the two runs by **shape**:

| r | θ/φ (massless) | θ/φ (massive) |
|---|---|---|
| 2.1 | 0.067 | 0.643 |
| 4.1 | 0.042 | 0.924 |
| 5.1 | 0.078 | 1.224 (peak) |
| 5.9 | 0.115 | 1.279 (peak) |
| 6.9 | 0.176 | 1.007 |
| 8.5 | 0.368 | 0.774 |
| 9.1 | 0.550 | 0.707 (falling) |

- **Massless θ/φ rises monotonically to the edge** (0.04 → 0.55): the massless sector
  out-ranges the massive φ at large r — the longer-range *shape* H1 predicts.
- **Massive-θ θ/φ peaks at r≈5 then falls**: a localized (bound) cloud, not a tail.
  Giving θ a mass converts the rising (long-range) shape into a turn-over (localized) —
  the H1 control signal, on shape.

#### NULL / SURPRISE — the dynamical amplitude is OPPOSITE the static prediction
The naive precheck picture (massless → large long-range field; massive → small screened
field) is **wrong on amplitude**, because of radiation vs binding:
- **Massless θ radiates away** (propagates at c → absorbed at the boundary): small
  residual interior amplitude (`θ_rms≈0.002` in diag), with an outgoing-wave bump still
  in transit (θ/φ non-monotonic mid-range; θ pow fit `+0.85`, i.e. *rising*, not a
  static tail).
- **Massive θ cannot radiate → binds** into a ~5× larger quasi-static cloud
  (`θ_rms≈0.008`), localized around the soliton.

So the massive run has *more* θ in the far field in absolute terms (θ_massless/θ_massive
≈ 0.09–0.29 over r∈[4,7]) — the reverse of "mass screens the field." The static
linear-response precheck (Coulomb vs Yukawa) does **not** capture this; the real
dynamical distinction is **radiate (massless) vs bind (massive)**, which is a *different*
and genuine physical signature. Recorded as the key surprise of the run.

#### What this run does NOT establish (confounds — a clean T-H1b needs a redesign)
1. **φ is not a clean Yukawa baseline.** φ is massive (`m=1.5`) in both runs yet falls as
   a power law (`φ pow ≈ −1.7…−2.0`), not exponentially — because `φ_rms` is dominated by
   the *extended* neutron body + background, not a linear Yukawa tail. The "massive ⇒
   exponential" leg is therefore **not** cleanly demonstrated.
2. **Not steady state.** Both configs lose energy to the absorbing boundary (E_total
   20.8→15.9 massless; 16.4→14.5 massive over t56–60); the stamped neutron is relaxing.
3. **θ still in transit** (outgoing-wave bump), so the massless profile is not asymptotic.
4. The **monopole question** (does the real θ source carry the v6 `Q_eff` monopole?) was
   not isolated — the extended source mixes multipoles.

**Verdict**: SUGGESTIVE, in the H1 direction *on shape* (massless out-ranges massive;
mass localizes), but **not a clean quantitative confirmation**. The run's most valuable
output is the surprise: dynamically the sectors differ by **radiation vs binding**, not by
the static Coulomb-vs-Yukawa range the precheck assumed. A clean T-H1b requires a
**redesign**: a *localized controlled source* (Gaussian density bump via a new
`sfa/seed/` generator, not the extended neutron), evolved to steady state (or with the
static field solved directly), giving a clean massive-Yukawa φ baseline and an isolated
monopole. H1's two unsolved problems (10⁴⁰ magnitude, scalar-not-tensor) remain untouched.

### Stage 2b — clean rerun with a localized impulse source (V100, 2026-06-05)

Stage 2's confounds (extended decaying neutron; φ not a clean Yukawa baseline) motivated a
**compact source**: a Gaussian φ bump (`gen_oscillon -sigma 1.5 -A 0.8`, 192³, ±15 box,
`A_bg=0` clean vacuum). Finding en route: **the bump is not a stable oscillon — it
disperses** (φ_max 0.8→0.07 by t≈24, E drift −73%). That *reframes the experiment as an
impulse*: watch which sector carries the released energy outward. Two runs, `m_θ=0` vs
`m_θ=1.5`, `T=20`, `snap_dt=2` (fine time resolution). Data: `impulse_massless.sfa`,
`impulse_massive.sfa` in `/space/scp/v64_h1b/`.

#### POSITIVE — clean, quantitative H1 segregation via energy-weighted radius
The robust observable is the volume-weighted mean radius of each sector's field energy,
`⟨r⟩ = Σ r·f_rms²·n_vox / Σ f_rms²·n_vox` (no stable source or exponent fit needed):

| t | ⟨r⟩_θ massless | ⟨r⟩_θ massive | ratio | ⟨r⟩_φ massless | ⟨r⟩_φ massive |
|---|---|---|---|---|---|
| 4 | 3.45 | 2.57 | 1.34 | — | — |
| 6 | 4.10 | 3.28 | 1.25 | 2.94 | 2.69 |
| 8 | 6.17 | 4.11 | 1.50 | — | — |
| 10 | 6.38 | 4.85 | 1.32 | 4.08 | 3.65 |

Three things, all in the H1 direction and **mutually consistent**:
1. **Massless θ out-ranges massive θ** at every time (⟨r⟩ ratio 1.25–1.50), and the gap
   **grows with time** — the massless sector migrates energy outward (long-range), the
   massive sector holds it nearer the core (short-range/bound).
2. **φ control validates the method**: φ is massive in *both* runs, and ⟨r⟩_φ is nearly
   identical between them (2.94/2.69, 4.08/3.65). The θ difference is therefore genuinely
   the θ **mass**, not a source/back-reaction artifact.
3. **The split matches H1's force assignment**: massless ⟨r⟩_θ (6.38) **exceeds** the
   massive φ (4.08) — massless out-ranges the "strong" sector — while massive ⟨r⟩_θ (4.85)
   ≈ φ (3.65) — a massive θ is as confined as φ. Profile shape confirms it: at t=10 the
   massless θ is *depleted at the core and peaks in an outgoing shell* (shell/core ≈ 3.2),
   the massive θ is a *core-peaked bound cloud* (shell/core ≈ 0.48).

#### What Stage 2b still does NOT establish
- It confirms **massless reaches/migrates further than massive** (the propagation/range
  half of H1), via a transient impulse. It does **not** show a **static** long-range tail
  (the gravitational 1/r) — the dispersing impulse cannot form one. The static-tail claim
  remains at the analytic-precheck level (Stage 1).
- **Magnitude (10⁴⁰) and tensor (h=±2) untouched** — H1's two real blockers (T-H1c).
- No clean falloff *exponent* (transient, not steady state); ⟨r⟩ is an integral measure.

**Net T-H1b verdict (Stages 1+2+2b)**: the carrier-mass split of H1 is **confirmed at the
range/propagation level** — massless sector longer-range than massive, robustly and
quantitatively, across an extended-soliton source (Stage 2, shape) and a compact-impulse
source (Stage 2b, ⟨r⟩). The *gravitational* reading still needs (i) a static 1/r tail from
a stable source, (ii) the 10⁴⁰ magnitude, (iii) the h=±2 tensor content — none addressed.

### RUNNER_FIX (incidental, 2026-06-05)
Provisioning was blocked by two bugs in `sfa/runner` (now fixed, rebuilt):
1. **Missing `rentable:{eq:true}` in the offer query** (`vast.go` `resolveGPUSpec`) — Vast
   `/bundles/` was returning non-rentable catalog machines, so every create hit
   `404/3603 no_such_ask … not available`. **This was the real root cause** (not field
   type, not pure contention).
2. **Stale-offer reuse** (`remote.go` Provision loop) — searched offers once then looped a
   frozen snapshot; now re-queries fresh offers before each create attempt.
Also added vast-python create-payload parity (`target_state:"running"`, `runtype:"ssh"`,
`client_id:"me"`, float `disk`). Requires an MCP server reconnect to load the new binary.
