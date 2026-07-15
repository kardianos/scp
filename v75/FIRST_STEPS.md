# v75 FIRST STEPS — Concrete actions (fidelity-aware)

**Updated**: 2026-07-13  
**Read first**: `PROPOSAL.md` (why), `FINDINGS.md` (what we know)  
**Kernel policy**: multi-fabric kernel **authorized**. Path locked:
`MULTIFABRIC_SPEC.md` — **B(1) Shape β** (C/Q/L arrays, Φ_Q≡Φ_C), then **B(2)**.
E-lite Step 1–2 done; next is B1 implementation + G0–G3.

---

## Goal of the first sprint

Answer, in **full 3D (L0)**, with the **existing kernel**:

> Can opposite-charge Q-balls attract under Coulomb without immediately
> annihilating, and what is the force?

That is the **positronium force gate** (FUTURE #5 force half) before multi-fabric
kernel work. Sprint uses **Option E-lite** (\(\pm\omega\)).

**Step 1 result: PASS** — see FINDINGS F8–F9. Opposite attract, same repel at
D=20, two clusters over T=100.

**Step 2 result: PARTIAL PASS** — see FINDINGS F10. Coulomb sub/circ/super orbit
scan behaves classically at D=20 (2 clusters, T=400, ~0.12 rev on circular);
head-on merges at contact. Multi-fabric still deferred for wide orbits; required
before overlapping N+e.

---

## Step 0 — Doc freeze (done when this file lands)

| Deliverable | Path | Status |
|-------------|------|--------|
| Architecture + options + principles | `v75/PROPOSAL.md` | done |
| Fidelity ladder L0–L3 | `v75/PROPOSAL.md` §6 | done |
| Setup + findings log | `v75/FINDINGS.md` | done |
| This checklist | `v75/FIRST_STEPS.md` | done |
| CLAUDE default → v75 | `CLAUDE.md` | done |

**Still to freeze (short text, no sim):**

- [ ] Charge table for Option B: \((q_C,q_Q,q_L)\), private \(G_C\) vs U(1)\(_\mathrm{EM}\)
- [ ] “Engagement off-list” (forbidden operators) as a one-page table

---

## Step 1 — L0: opposite-charge force law (3D, existing kernel)

### 1.1 Seed

```bash
# Two light balls, opposite charge, large separation, co-centered box
# Profile: v74 light nucleon
PROF=v74/profiles/f_w146_g005.txt
N=192
L=48          # larger than v74 nuclear L=36 for cleaner force
D=16          # separation; keep D < L/2 for boundary cleanliness (v70)
W=1.46

bin/gen_qball_multi $N $L /space/scp/v75/seeds/pm_force_D${D}.sfa \
  $PROF  $W   0.0   $((D/2))  0  0 \
  $PROF -$W   0.0  -$((D/2))  0  0
```

(Use python for centers if shell arithmetic is awkward; pattern matches
`gen_qball_multi` docs: **negative omega = opposite charge**.)

### 1.2 Config template

```
N = 192
L = 48
T = 80
dt_factor = 0.025
m = 1.5
m_theta = 1.6
eta = 0
complex_phi = 1
complex_gauge = 1
g_gauge = 0.05
bc_type = 0
damp_width = 4.0
damp_rate = 0.01
init = sfa
init_sfa = /space/scp/v75/seeds/pm_force_D16.sfa
output = /space/scp/v75/results/pm_force_D16.sfa
diag_file = /space/scp/v75/results/pm_force_D16_diag.tsv
precision = 0
snap_dt = 2.0
diag_dt = 0.25
```

**Must set N,L in cfg** to match seed (v74 lesson).

### 1.3 Hardware

- Prefer **≥24 GB GPU** if available (nested later); V100-16GB OK for N=192.  
- Runner: `sim_setup` remote + `sim_build` cuda sources as in v74.  
- Auto-download diags immediately; SFA can wait.

### 1.4 Analysis (success criteria)

| Metric | Pass | Fail |
|--------|------|------|
| Force sign | Attraction (D decreases or \(F\cdot\hat r < 0\)) | Repulsion |
| Rough exponent | \(F\sim 1/r^n\), \(n\approx 2\) at mid D | Yukawa/wrong |
| Gauss | `gauss_max ~ 1e-13` | Drift ⇒ bug |
| Cluster count | 2 for most of T | Instant merge to 1 |
| Charge | \(Q_\mathrm{tot}\approx 0\), \(\|Q_\pm\|\) stable early | Instant cancel |

Tools: diag.tsv; `sfa_qball_track` / cluster tools if available; slice renders
at t=0, mid, end (`sfa_slice` + `render_slices.py`).

### 1.5 Scan (optional follow-up)

D ∈ {12, 16, 20, 24} at fixed ω=1.46; one D with ω=1.42 (heavier).

**Record all results in FINDINGS.md §Phase1.**

---

## Step 2 — L0: positronium dynamics (orbit attempt)

### 2.1 Seed

Same ± pair, D≈16–20, add **equal opposite boosts** tangential (need seed
tool support: `gen_qball_boost` per ball or extend multi seed with velocity —
if no velocity API, use kernel only after checking pair generators).

Check first:

```bash
# If gen_qball_pair or multi has no velocity: document gap and either
# (a) small seed-only patch (user-authorized if in seed/, not kernel), or
# (b) start from rest and use Coulomb fall as “radial orbit” null.
```

### 2.2 Success criteria

| Outcome | Interpretation |
|---------|----------------|
| Multiple revolutions, slow radiation | **Strong** go for atom program |
| One pass + capture/annihilation | Partial; need multi-fabric isolation |
| Immediate merge | Same-fabric bag wins; multi-fabric required |
| Immediate disperse | Window/charge issue |

T ≥ 200–400; large L; absorb boundaries.

### 2.3 Results (2026-07-13) — PARTIAL PASS

Seed: `bin/gen_qball_pair_boost`, D=20, ω=±1.46, N=192 L=48 T=400.

| Run | vt/vr | Result |
|-----|-------|--------|
| vc | vt=0.0193 (naive circ.) | D 20.02→20.46, φ=43°, 2 clusters — near-circular arc |
| v05 | vt=0.010 | D→15.95 — inspiral |
| v15 | vt=0.029 | D→26.55 — outspiral |
| headon | vr=0.01 | Merge by T=400 (s_max→0.13, E_pot→−40) |

True circular from Step 1 a_C: **vt≈0.0154**. Period ~3300 → multi-rev needs
T≳3500 (Step 2b). Details: FINDINGS F10; `v75/analysis/*_Dt.tsv`.

---

## Step 3 — L2 radial scout (parallel, optional)

**Only after Step 1 is queued or if 3D queue is full.**

- Spherical BO: fixed nuclear ρ_Q background from gauged profile + light radial
  f_L(r) + A_0(r).  
- Stamp every plot: **not topology-safe**.  
- Purpose: choose \((m_L, g, Q_L)\) for Step 2–4 scale co-design.

Do **not** claim shells, spin, or multi-center from L2.

---

## Step 4 — L1 nested grid design note (no code yet)

Write `v75/NESTED_ATOM.md` sketch:

- Root: L large, coarse dx — EM between N and L  
- Child: fine on nuclear core  
- Optional second child on dense electronic peak  
- Subcycling as in `MULTI_RESOLUTION.md`

Implement only when Step 1–2 demand \(a_0 \gg r_N\) beyond single-grid memory.

---

## Step 5 — Decision gate (before multi-fabric kernel)

| If Step 1–2 show… | Then… |
|-------------------|--------|
| Long-lived ± binding | Prefer Option **A/E** implementation path; delay full C/Q/L |
| Attraction but rapid annihilation | Multi-fabric (private bag) justified — draft kernel spec for user auth |
| No Coulomb attraction | Gauge/seed bug hunt first |
| Numerics unclean | Bigger box / GPU / nested before theory change |

**Kernel multi-fabric: explicit user “modify scp_sim” authorization required.**

---

## Step 6 — L3 n-bar (explicitly later)

Do **not** start in this sprint. Prerequisites:

1. Measured \(V_{\mathrm{Coulomb}}(r)\) from L0  
2. Measured merge/annihilation cross-section or rate  
3. Engagement graph frozen  

Then: n-body with only allowed bars (EM on, fusion off or probabilistic).

---

## Directory layout for v75 runs

```
v75/
  PROPOSAL.md
  FINDINGS.md
  FIRST_STEPS.md
  cfg/           # pm_force_*.cfg, pm_orbit_*.cfg
  seeds/         # → prefer /space/scp/v75/seeds/
  results/       # → /space/scp/v75/results/ for SFAs
  analysis/      # force fit, track scripts
```

---

## Immediate next command block (when starting sims)

```text
1. mkdir -p /space/scp/v75/{seeds,results}
2. Build seed pm_force_D16 (Step 1.1)
3. Write v75/cfg/pm_force_D16.cfg (Step 1.2)
4. sim_setup remote (GPU filter: V100 or larger)
5. sim_build cuda sources (same set as v74)
6. sim_run id=pm_force_D16 + auto_download diags
7. Status loop until complete → analyze → FINDINGS.md
8. Render slices → v75/results/pm_force_D16.png
9. Teardown when downloads done
```

No sim launched in the doc-recording session unless user says “run Step 1.”

---

## Definition of done for “first sprint”

- [x] Fidelity ladder + first steps + findings recorded  
- [x] Step 1 force run complete with force-sign recorded (F8)  
- [ ] Slice render of ± pair (optional; 2/4 SFAs good: pp_D16, pm_D20)  
- [x] FINDINGS updated with pass/fail + multi-fabric impact (F9)  
- [x] Step 2 orbit/capture complete (F10 PARTIAL PASS)  
- [x] Go/no-go: E-lite continues for wide orbits; multi-fabric for contact-scale N+e  

## Step 2b / next

1. **Primary:** implement B(1) Shape β per `MULTIFABRIC_SPEC.md`; gates G0–G3  
2. After G3 PASS: lighter L and/or B(2) unlock — outcome-driven (§6 of SPEC)  
3. Optional non-blocking: E-lite multi-rev vt≈0.0154, T≥3500  

G3 success = head-on contact **without** merge (inverts F10 same-fabric merge).