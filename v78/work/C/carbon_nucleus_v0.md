# Carbon Nucleus Package v0 — Phase C (CP-C)

**Agent:** C  
**Campaign:** v78  
**Date:** 2026-07-19  
**Status:** **PASS** (Stage 2A Z-carbon); A=12 documented super-critical at g=0.05  
**Stamp:** **CP-C ADOPT** (see `logs/C_carbon_nucleus.log`)  
**Sources:** `v74/CARBON_MAP.md`, `v74/RESULTS.md`, `v74/cfg/*`, `v74/results/*`

---

## 1. Verdict (one screen)

| Target | Result | Status |
|--------|--------|--------|
| **Z-carbon** (`c6_light`) | 6 light free nucleons → single on-branch droplet, **Q→649.83** (0.71×Q_max) | **PASS** — Stage 2A fabric carbon nucleus |
| **A=12** (`c12_light`) | 12 light free nucleons → one droplet, **Q→1410.5** still **> Q_max=921** | **Documented** — super-critical control; not a parked isotope at g=0.05 |
| **Z=6 mapping** | Formation history = 6 nucleon units; identity = final Q + branch + history (not retained A) | Spec closed |
| **True A=12 park** | Requires lower g (e.g. g=0.02) or droplet+N⁰ channel | Path only — not run in v74 |

**Honest product name:** fabric **Z-carbon nucleus** (bare liquid-drop). Not an atom (Stage 3–4 / Phase L+A).

---

## 2. Z-carbon `c6_light` — PASS numbers

### 2.1 Setup

| item | value |
|------|-------|
| Recipe | 6× light nucleon, ω=1.46, co-phase δ=0 |
| Geometry | regular octahedron, edge D=10 (axial a=D/√2) |
| Profile | `v74/profiles/f_w146_g005.txt` |
| Physics | complex_phi=1, complex_gauge=1, **g=0.05**, m_θ=1.6, η=0, absorbing BC |
| Grid | N=192, L=36, dt_factor=0.025, T=300 |
| Hardware | GPU V100; wall 38.1 min (~71.8 ms/step) |
| Config | [`v74/cfg/c6_light.cfg`](../../v74/cfg/c6_light.cfg) |
| Seed | `/space/scp/v74/seeds/c6_light_N192.sfa` (repo copy: `v74/seeds/`) |
| Data | `/space/scp/v74/results/c6_light.{sfa,diag.tsv,log}`; render `v74/results/c6_light.png` |

Free-nucleon baseline (gscan light, ω=1.46): **Q_N=114.13**, **E/Q=1.5187**, r_{1/2}=2.63.

### 2.2 Endpoint table

| quantity | t=0 (seed) | t=300 (end) | notes |
|----------|------------|-------------|--------|
| Q | 907.84 | **649.83** | −28.4% radiated/absorbed |
| E | 1360.7 | 939.0 | radiation + sponge |
| E/Q | 1.4988 | **1.4450** | |
| Branch E/Q(Q) pred | — | **1.4377** | gscan g=0.05 interp |
| On-branch Δ | — | **+0.51%** | still-hot; he2 was +0.1% |
| gauss_max | 9.8e-14 | 1.0e-13 | floor held entire run |
| ω_core | 1.460 | 1.365 | drops during merger |
| r_core | 7.61 | 6.94 | |
| Q_core | 591 | 592 | core inventory ~stable late |

### 2.3 Inventory bookkeeping

| quantity | value |
|----------|-------|
| 6× free inventory | 6 × 114.13 = **684.8** |
| Seed Q | 907.8 → **+32.6%** co-phase interference excess |
| Final Q vs Q_max | 650 / 921 ≈ **0.71** (mid-high branch, subcritical) |
| Mass defect vs free light N | \(1 - 1.4450/1.5187\) = **4.85%** of rest mass per unit charge |
| Late dQ/dt (t>250) | ≈ **−0.17 / t.u.** — still slowly shedding; parked Q is upper bound |

### 2.4 Morphology (render OK)

- t=0: octahedron (z-slice shows 4 equatorial balls; poles off-plane)
- t=15–50: co-phase infall → single central droplet + multipole radiation
- t=100–300: single compact charge blob + Coulomb |E| ring; **no fission**
- Components co-located; θ inert

**Kill-gate pass:** single connected droplet; Gauss floor; on-branch within ~0.5%; no fatal radiation.

---

## 3. A=12 super-critical theorem (g=0.05)

### 3.1 Charge-budget theorem

At **g=0.05**, static gauged branch has:

- **Q_max = 921** (ω_min ≈ 1.40624)
- Lightest VK-stable free nucleon: **Q_N^min ≈ 90** (ω ≈ 1.4825)

Therefore:

\[
A_{\mathrm{free}} \cdot Q_N^{\min} \;\xrightarrow{A=12}\; 12 \times 90 \approx 1080 \;>\; Q_{\max}=921.
\]

**Theorem:** No choice of free gauged nucleons on the stable branch can fuse twelve of them into a **subcritical** static droplet at g=0.05. Always \(12\,Q_N^{\min} > Q_{\max}\).

### 3.2 `c12_light` control numbers

| item | value |
|------|-------|
| Recipe | 12× light (ω=1.46), co-phase, regular icosahedron D=10 |
| Config | [`v74/cfg/c12_light.cfg`](../../v74/cfg/c12_light.cfg) |
| Grid / T | N=192, L=36, T=300 (same as c6) |
| Seed Q | **1959.9** (2.13× Q_max); free 12×114.13=1370 → **+43%** interference |
| Final Q (t=300) | **1410.5** (1.53× Q_max; still **+489** above cap) |
| Final E/Q | 1.4265 |
| Late dQ/dt (t>250) | **−0.30 / t.u.** |
| Linear time-to-cap (rough) | ~1600 t.u. more from t=300 |
| gauss_max | ~2e-13 floor held |
| Fate | merge → **one** droplet → **evaporate toward Q_max** (li3 channel, larger excess) |
| Mass defect vs free light N | 6.07% (still hot / off-branch above Q_max) |

**Confirmed:** free A=12 cannot park as a static branch isotope at g=0.05. Honest A=12 product = **hot evaporating droplet**; stable parked carbon = **Z=6 (`c6_light`)** path.

### 3.3 Stage 2A summary table

| run | A | seed Q | final Q | vs Q_max | fate |
|-----|---|--------|---------|----------|------|
| **c6_light** | 6 | 908 | **650** | 0.71× | **parks mid-branch** (Z-carbon) |
| **c12_light** | 12 | 1960 | **1411** | 1.53× | **still evaporating** (hot A=12) |

---

## 4. How Z=6 maps (conserved bookkeeping)

Parallels to real \(^{12}\)C are **structural**, not quantitative (CONCEPT.md).

| Real carbon | Fabric analog | Bookkeeping |
|-------------|---------------|-------------|
| Atomic number **Z=6** | Six units of **nucleon inventory in formation history** (six light free nucleons fused) | History + total Q |
| Mass number A=12 | Twelve nucleon fusions when charge budget allows | A is history; fused droplet **forgets** A (v71 doctrine) |
| Nucleus | Liquid-drop **fused charge droplet** on gauged Q-ball branch | Final Q, E/Q, on-branch match |
| Electrons (6) | Light opposite-charge solitons | **Phase L / Stage 3** — out of scope for C |
| Atom | Nucleus + 6 opposite charges via A | **Phase A / Stage 4** |

**Why Z-carbon ≠ real \(^{12}\)C mass number:** at g=0.05 the charge budget forces the **stable** primary target to be **A_free=6** inventory (Z map), not free A=12. The droplet’s identity is:

1. **Q** parked mid-branch (~650 at T=300)
2. **E/Q** on-branch (~1.445 vs pred 1.4377)
3. **Formation history** (six co-phase light nucleons, octahedron)
4. **Single cluster** morphology

Not: retained count of six distinct nucleon centers (Stage 2B multi-center is optional research).

**Carbon-specific labels (from CARBON_MAP):**

| if… | we call it… |
|-----|-------------|
| 6 light nucleons → stable on-branch droplet | **Z-carbon nucleus** ← **this package** |
| 12 light nucleons → single droplet evaporating to Q_max | **A=12 hot carbon** (budget-limited) |
| multi-center static bound state with A retained | multi-center carbon (Stage 2B) |
| nucleus + 6 opposite-charge light partners | carbon atom (Phase A; blocked on L) |

---

## 5. Path to true A=12 (parked)

### Path A — Lower g (raise Q_max)

| g | Q_max (approx) | 12× standard Q_N(~311) | 12× light Q_N(~114) | status |
|---|----------------|------------------------|---------------------|--------|
| 0.05 | 921 | always super-critical | always super-critical free | **done (control)** |
| **0.02** | **~5300** | 12×311 ≈ 3732 < 5300 ✓ | easy | **deferred** — need g=0.02 gauged profile shoot |

Recommended later run: **c12_g02** — 12× standard (or light) at g=0.02; park if Q_max holds.

**Blockers today:** no g=0.02 radial profiles in tree; shooter + seed + GPU production not yet run.

### Path B — Droplet + N⁰ (isotope knob)

At fixed **Z-scale** parked charge (c6-like Q inventory):

1. Park **Z-carbon** droplet (done at g=0.05).
2. Attach or co-form **neutron-analogs (N⁰)** — EM-neutral nuclear species (Phase P / v75 C-only bag) — to raise **A** without raising EM charge budget in the same way as free charged nucleons.
3. Requires: CP-P proton/neutron package + multi-fabric or flavor-neutral construction; **not available as a pure free-nucleon fusion at g=0.05**.

### Path C — Cooler / longer park of super-critical blob (not true A=12)

c12_light will continue evaporating toward Q_max≈921. That endpoint is a **capped hot droplet**, not an A=12 isotope (inventory history is erased; final Q is Q_max-scale, not 12×Q_N). **Do not** call that parked A=12.

### Recommended order

1. **Adopt Z-carbon as Stage 2A / CP-C nuclear core** for atom recipes (A phase).
2. Optional: longer T on c6 for dQ/dt→0 (polish, not gate).
3. True A=12: shoot g=0.02 profiles **or** droplet+N⁰ after CP-P.
4. Stage 2B multi-center only if research priority shifts.

---

## 6. Config pointers (`v74/cfg`)

| path | role |
|------|------|
| [`/home/d/code/scp/v74/cfg/c6_light.cfg`](/home/d/code/scp/v74/cfg/c6_light.cfg) | **Production Z-carbon** — N=192 L=36 T=300 g=0.05; init=`/space/scp/v74/seeds/c6_light_N192.sfa` |
| [`/home/d/code/scp/v74/cfg/c12_light.cfg`](/home/d/code/scp/v74/cfg/c12_light.cfg) | **A=12 super-critical control** — same grid/physics; 12-ball seed |
| [`/home/d/code/scp/v74/cfg/c6_light_smoke.cfg`](/home/d/code/scp/v74/cfg/c6_light_smoke.cfg) | Local CPU smoke N=128 L=28 T=20; seed `v74/seeds/c6_light_N128.sfa` |

### Shared physics block (both production cfgs)

```
N=192  L=36  T=300  dt_factor=0.025
m=1.5  m_theta=1.6  eta=0  mu=-41.345  kappa=50
complex_phi=1  complex_gauge=1  g_gauge=0.05
bc_type=0  damp_width=3.0  damp_rate=0.01
init=sfa  precision=0  snap_dt=2.5  diag_dt=0.25
```

**Note:** N,L **must match the seed**; cfg load path does not auto-override grid from SFA header.

### Seed generation

```bash
# from repo root
./v74/analysis/make_carbon_seeds.sh 192 36
# uses bin/gen_qball_multi + profiles/f_w146_g005.txt
# → v74/seeds/c6_light_N192.sfa, c12_light_N192.sfa
```

### Profiles

| path | use |
|------|-----|
| `v74/profiles/f_w146_g005.txt` | light nucleon (carbon campaign default) |
| `v74/profiles/f_w142_g005.txt` | standard nucleon (He/Li history; too heavy for free A≥4 at g=0.05) |
| `v74/profiles/f_w1485_g005.txt` | edge profile (present; free nucleons: stay on VK-stable ω≤1.4825) |

### Results / renders

| path | role |
|------|------|
| `v74/RESULTS.md` | full Stage 2A writeup |
| `v74/CARBON_MAP.md` | Stage 0 mapping + theorem |
| `v74/STATUS.md` | campaign status (c6+c12 done) |
| `v74/results/c6_light.png` | Z-carbon slice panel |
| `v74/results/c12_light.png` | A=12 control panel |
| `/space/scp/v74/results/c6_light.*` | full SFA (~35 GB) + diag |
| `/space/scp/v74/results/c12_light.*` | full SFA (~42 GB) + diag |

---

## 7. Observables checklist (met for both runs)

| # | observable | c6 | c12 |
|---|------------|----|-----|
| 1 | Q(t), E(t), E/Q(t) from diag | ✓ | ✓ |
| 2 | On-branch E/Q vs gscan | +0.51% | off-branch (above Q_max) |
| 3 | Mass defect vs free N | 4.85% | 6.07% (hot) |
| 4 | Single droplet | ✓ | ✓ |
| 5 | Gauss floor ~1e-13 | ✓ | ✓ |
| 6 | Fragmentation / fission | none | none |

---

## 8. Handoff for A / L / U / P

**FOR_A (atom core):** Use **Z-carbon** as nuclear core: Z≈6 formation inventory, parked Q~650 at g=0.05. Atom = this core + 6 light opposite charges (needs CP-L). Do **not** require free A=12 for Stage 4 atom structural analog.

**FOR_L:** Electron/positronium sector still blocks atoms; C package does not change that.

**FOR_P:** N⁰ (EM-neutral) is the clean **isotope / A-knob** on a Z-carbon core without violating charge budget (Path B).

**FOR_U:** Phase C complete for board: Z-carbon PASS; A=12 super-critical theorem + control PASS-as-documented. Unlock **A core** per CAMPAIGN_MAP. Stamp co-agree: C **ADOPT**; await U co-stamp for full CP-C board green.

**Not claimed:** chemistry, quantitative 12 u, multi-center A retention, spontaneous abundance peak (Stage 5), or kernel changes.

---

## 9. File index (this package)

| path | role |
|------|------|
| `v78/work/C/carbon_nucleus_v0.md` | this document |
| `v78/work/C/README.md` | work-dir index |
| `v78/logs/C_carbon_nucleus.log` | append-only agent log + CP-C stamp |
