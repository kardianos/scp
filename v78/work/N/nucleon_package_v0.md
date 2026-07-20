# Nucleon Package v0 — Phase N (CP-N)

**Agent:** N  
**Date:** 2026-07-19  
**Status:** **ADOPT** — freeze for P/C downstream  
**Sources:** `v78/GOALS.md` §3.1, `v74/CARBON_MAP.md` §2, `CONCEPT.md` §3, `v69/theory/gscan.tsv`, `v71/NUCLEI.md`, `v74/RESULTS.md`

---

## 1. Definition

A **fabric nucleon** (Stage-1 particle / free “baryon”) is a **gauged diagonal-U(1) Q-ball**:

| property | value |
|----------|--------|
| Ansatz | \(\Phi_a = f(r)\,e^{i\omega t}\) for **all three** components equal; \(\Theta=0\) |
| Potential | \(V_t(s)=(\mu/2)\,s/(1+\kappa s)\), \(s=\Pi_a\|\Phi_a\|^2\) |
| Standard params | \(m=1.5\) (\(m^2=2.25\)), \(\mu=-41.345\), \(\kappa=50\) |
| Gauge | `complex_gauge=1`, \(g=0.05\) (default campaign) |
| \(\theta\)-drain closure | \(m_\theta=1.6>\omega\) |
| Curl coupling | \(\eta=0\) for particle experiments (no \(\Theta\) partner required) |
| Charge | \(Q=\int 3\,\omega_{\mathrm{eff}}\,f^2\,\mathrm{d}V\) (Noether); \(Q_a=Q/3\) |
| Identity | Branch point \((Q,E/Q,\omega)\), not a chemistry nucleon number |

**Not this package:** flavored multiplet \((\omega_a,Q_a)\) (v71 research), \(\eta>0\) stationary \(\eta\)-soliton (`eta_qflow`), multi-fabric P/N⁰ (Phase P / v75 B2). Those build **on** this free-nucleon template.

**Profile convention (v72):** \(f(r)\) is the **per-component** amplitude. Seed writers stamp \(\Phi_a=f\) for each \(a\), not \(f/\sqrt{3}\).

---

## 2. VK stability window (g = 0.05)

Vakhitov–Kolokolov: classical stability on the branch where \(\mathrm{d}Q/\mathrm{d}\omega < 0\).

| | |
|--|--|
| **Existence (gauged)** | \(\omega \in (\omega_{\min},\, m) = (1.40624,\, 1.5)\) |
| **VK-stable** | \(\omega \in [1.40624,\, 1.4825]\) with \(\mathrm{d}Q/\mathrm{d}\omega < 0\) (`dQdw_sign = -1` in gscan) |
| **Unstable edge** | \(\omega \gtrsim 1.485\): sign flips (`dQdw_sign = +1`); **do not seed free nucleons there** |
| **\(Q_{\max}\)** | **921** at \(\omega_{\min}=1.406245\) — Coulomb self-repulsion cap; no larger static ball |
| **\(Q_N^{\min}\)** (stable) | **≈ 90** at \(\omega=1.4825\) |

**Ungauged reference** (theory only; not default seed path): \(\omega \in (1.3087, 1.5)\), no \(Q_{\max}\) at finite g=0.

**Branch data:** `/home/d/code/scp/v69/theory/gscan.tsv` (columns: `g omega … Q E_total … dQdw_sign r_half`).  
**E/Q for on-branch tests:** use \(E_{\mathrm{total}}/Q\) from gscan (matter + field).

---

## 3. Light / standard nucleon table (g = 0.05)

From gscan rows + `v74/CARBON_MAP.md` (E/Q = E_total/Q):

| label | \(\omega\) | \(Q_N\) | \(E/Q\) | \(r_{1/2}\) | \(f(0)\) | profile | use |
|-------|------------|---------|---------|-------------|---------|---------|-----|
| **standard** | 1.42 | **311.46** | **1.4649** | 3.80 | 0.6354 | `v74/profiles/f_w142_g005.txt` | v71 He (he2); heavy; too large for A≥4 free at g=0.05 |
| **light (default)** | 1.46 | **114.13** | **1.5187** | 2.63 | 0.6061 | `v74/profiles/f_w146_g005.txt` | **default free nucleon** for carbon / multi-ball (v71 li3s, v74 c6/c12) |
| **lightest VK-stable** | 1.4825 | **90.01** | **1.5321** | 2.34 | 0.5670 | shoot on demand | minimize \(A\cdot Q_N\) |
| edge (unstable) | 1.485 | 89.70 | 1.5316 | 2.32 | 0.5608 | `v74/profiles/f_w1485_g005.txt` | **research only** (v75 L hierarchy); not free-nucleon default |

**Campaign defaults (GOALS.md):** light nucleon \(\omega=1.46\), \(Q_N\approx 114\), profile `f_w146_g005`.

### Bookkeeping quantities (every free-nucleon run)

1. Global \(Q(t)\), \(E(t)\), \(E/Q(t)\)  
2. On-branch: final \(E/Q\) vs gscan at measured \(Q\)  
3. \(r_{1/2}\) or core radius  
4. Per-component \(Q_a \approx Q/3\)  
5. `gauss_max` ~ \(10^{-13}\) floor  

---

## 4. Seed recipe

### 4.1 Radial profile

**Preferred (gauged, matches gscan):** use frozen files under `v74/profiles/`, or re-shoot with `v69/theory/gauged_shooter_fast.py` and extract 2-column `r f`.

**Ungauged only** (not for gauged particle seeds): `bin/radial_qball -omega W -profile out.txt`  
— window \(\omega\in(1.3087,1.5)\); **wrong** Coulomb envelope at \(g=0.05\).

Profile format for seed tools: comment lines OK; data rows `r f` strictly increasing \(r\).

### 4.2 Single static nucleon SFA

`gen_qball_boost` with \(v_x=0\) (or one-ball `gen_qball_multi` at origin):

```bash
# Light default — N,L match later cfg
bin/gen_qball_boost 128 24 \
  v74/profiles/f_w146_g005.txt 1.46 0.0 \
  /space/scp/v78/seeds/nucleon_light_N128.sfa

# Standard
bin/gen_qball_boost 128 24 \
  v74/profiles/f_w142_g005.txt 1.42 0.0 \
  /space/scp/v78/seeds/nucleon_std_N128.sfa
```

Equivalent multi:

```bash
bin/gen_qball_multi 128 24 out.sfa \
  v74/profiles/f_w146_g005.txt 1.46 0.0  0 0 0
```

**Output:** 24-column complex matter SFA. Gauge \(\mathbf{E}\) is built by the kernel **init Gauss projection** (`complex_gauge=1`).

### 4.3 Kernel config skeleton (`init=sfa`)

```
N = 128          # MUST match seed
L = 24           # MUST match seed
T = 100
dt_factor = 0.025
m = 1.5
m_theta = 1.6
eta = 0
mu = -41.345
kappa = 50
complex_phi = 1
complex_gauge = 1
g_gauge = 0.05
bc_type = 0
damp_width = 3.0
damp_rate = 0.01
init = sfa
init_sfa = /path/to/nucleon_light_N128.sfa
init_frame = 0
```

GPU production often uses N=192, L=36 (v70/v71 nuclear grid). Cfg path does **not** auto-override N,L from the SFA header — match by hand.

### 4.4 Multi-nucleon (downstream C; not Phase N product)

`bin/gen_qball_multi` — groups of `profile omega delta x y z` (max 16 balls).  
Co-phase fusion: \(\delta=0\), nearest-neighbour \(D\approx 10\). See `v74/analysis/make_carbon_seeds.sh`.

Negative \(\omega\) = opposite charge (v75 force/orbit prep).

---

## 5. Stability evidence

### 5.1 Theory (`CONCEPT.md` §3)

- Complex Q-ball is the unique minimal cure of Derrick + harmonic radiation.  
- **VK:** \(\mathrm{d}Q/\mathrm{d}\omega<0\) branch classically stable.  
- **Charge retention:** bare ball retains **99.99985%** of charge over **1000 t.u.** under absorbing BC.  
- **Package:** gauged diagonal U(1) + \(m_\theta>\omega\) closes \(\theta\)-drain; \(\eta=0\) keeps \(\Theta\) inert for free nucleons.  
- Static evaporation bound \(E<mQ\) does **not** govern dynamics (balls with \(E>mQ\) are classically immortal).

### 5.2 Static branch (`v69` gscan)

- Full g=0.05 scan: \(Q(\omega)\) from 921 → ~90 on the stable side; \(\mathrm{d}Q/\mathrm{d}\omega\) sign table as §2.  
- Coulomb field energy fraction small for light balls (\(E_{\mathrm{field}}/E_{\mathrm{total}}\sim 0.25\%\) at \(\omega=1.46\)).

### 5.3 Dynamics as free nucleon building blocks (`v71`)

| run | free nucleon | fate | implication |
|-----|--------------|------|-------------|
| **he2** | 2× standard \(\omega=1.42\), \(Q_N=311.5\) | Fuse → on-branch He droplet Q=604.6, Δ(E/Q)≈0.1% | Standard nucleon is a legitimate free unit |
| **li3s** | 3× **light** \(\omega=1.46\), \(Q_N=114.1\) | Clean mid-branch Li, mass defect ~3.5%/charge | **Light default validated** |
| **li3** | 3× standard | Super-critical evaporation toward \(Q_{\max}\) | Budget law, not nucleon failure |

Quark bookkeeping: \(Q_a=Q/3\) through fusion.

### 5.4 Carbon campaign (`v74`)

| run | free nucleons | seed / final Q | fate |
|-----|---------------|----------------|------|
| **c6_light** | 6× light | 908 → **650** | Parks mid-branch (0.71 \(Q_{\max}\)); on-branch Δ≈+0.5%; **Z-carbon** |
| **c12_light** | 12× light | 1960 → 1411 | Super-critical hot droplet (budget theorem) |

Both use the same light free-nucleon seed class; Gauss floor held entire T=300.

### 5.5 Phase N pass bar (GOALS.md)

| criterion | status |
|-----------|--------|
| Spec of object + bookkeeping | **this package** |
| Profile + seed recipe | **§3–4 + `recipes.sh`** |
| Stability evidence | **§5 (CONCEPT / v69 / v71 / v74)** |
| Kernel edits | **none required** |

---

## 6. Explicit non-goals (this package)

- Proton vs neutron EM distinction (Phase **P**)  
- Carbon droplet identity (Phase **C** — consumes this free nucleon)  
- Light opposite sector / atom (Phases **L**, **A**)  
- \(\eta>0\) or multi-fabric kernel changes  

---

## 7. Files

| path | role |
|------|------|
| `v78/work/N/nucleon_package_v0.md` | this freeze |
| `v78/work/N/recipes.sh` | exact commands |
| `v74/profiles/f_w146_g005.txt` | light nucleon \(f(r)\) |
| `v74/profiles/f_w142_g005.txt` | standard nucleon \(f(r)\) |
| `v74/profiles/f_w1485_g005.txt` | edge (use with care) |
| `v69/theory/gscan.tsv` | gauged branch table |
| `bin/gen_qball_boost` | single (boosted or static) seed |
| `bin/gen_qball_multi` | N-ball / nuclear seeds |
| `bin/radial_qball` | ungauged radial shooter only |

---

## 8. Handoff

**STAMP CP-N: ADOPT**

Unlocks: Phase **P** (P/N⁰ package), Phase **C** (carbon uses light free nucleon as inventory unit).

**FOR_U:** freeze board row CP-N = ADOPT; copy light default into `PARTICLE_LADDER.md` / `RECIPES.md` when written.  
**FOR_P:** free nucleon = C-bag charge inventory template; EM proton/neutron is multi-fabric, not a different \(\omega\) in this package.  
**FOR_C:** primary free unit remains light \(\omega=1.46\); budget theorem unchanged.
