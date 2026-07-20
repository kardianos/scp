# v78 RECIPES — Runnable path Field → Nucleon → P/N → C → L → Atom

**Campaign:** v78  
**Date:** 2026-07-19  
**Integrator:** U (from W/N/P/C/L/A packages)  
**Policy:** no `scp_sim` / `sfa.h` edits without human OK; prefer existing profiles/seeds.

---

## 0. Global physics block

```text
m = 1.5
m_theta = 1.6
eta = 0
mu = -41.345
kappa = 50
complex_phi = 1
complex_gauge = 1
g_gauge = 0.05
bc_type = 0          # absorbing
damp_rate = 0.01
dt_factor = 0.025
```

| Use | N | L | damp_width |
|-----|---|---|------------|
| Single nucleon smoke | 128 | 24 | 3 |
| Nuclear liquid-drop (v74) | 192 | 36 | 3 |
| Multi-fab Z2 | 128 | 32 | 4 |
| Multi-fab Z6 atom | **192** | **48** | **5** |

**N,L must match the seed** — cfg path does not auto-override from SFA header.

---

## N — Free nucleon

**Package:** `work/N/nucleon_package_v0.md` · script `work/N/recipes.sh`

### Profiles (gauged g=0.05)

| Label | ω | Q | E/Q | Profile |
|-------|---|---|-----|---------|
| **light default** | 1.46 | 114.13 | 1.5187 | `v74/profiles/f_w146_g005.txt` |
| standard | 1.42 | 311.46 | 1.4649 | `v74/profiles/f_w142_g005.txt` |
| edge (not free default) | 1.485 | ~90 | — | `v74/profiles/f_w1485_g005.txt` |

### Seed (static light)

```bash
# print + gscan fingerprint
./v78/work/N/recipes.sh

# write seeds (needs /space or set SEED_OUT)
SEED_OUT=/space/scp/v78/seeds N=128 L=24 ./v78/work/N/recipes.sh run-light

# manual
bin/gen_qball_boost 128 24 v74/profiles/f_w146_g005.txt 1.46 0.0 \
  /space/scp/v78/seeds/nucleon_light_w146_N128.sfa
```

### Config skeleton

```text
N=128 L=24 T=100
# + global physics block
init=sfa
init_sfa=/space/scp/v78/seeds/nucleon_light_w146_N128.sfa
```

### Checks

- gscan: `./v78/work/N/recipes.sh gscan-light`
- Branch: `v69/theory/gscan.tsv` g=0.05
- Do **not** use ungauged `radial_qball` for g=0.05 production seeds

---

## P — Proton / neutron analogs (B2)

**Package:** `work/P/pn_package_v0.md`

### Freeze

```text
n_fabrics=3
mf_lock_CQ=0
mf_stage=2
q_C=0  q_Q=1  q_L=-1
```

| Species | Seed |
|---------|------|
| **p** | C+Q co-located |
| **n** | C only (Q empty at site) |

### Seed tool

```text
bin/gen_pn_core N L profNuc omega out_C out_Q out_L \
  nZ  xz yz zz ... \
  nN  xn yn zn ... \
  nL  xl yl zl ... \
  [profL omegaL]
```

### Configs (`v75/cfg/pn/`)

| Config | Role |
|--------|------|
| `f17d_b2n.cfg` | unit neutron |
| `f17e_b2p.cfg` | unit proton |
| `f17f_z1n1.cfg` | Z1N1 |
| `f17g_z2n0.cfg` / `f17h_z2n2.cfg` | Z2 isotope pair |
| `p2_z2n0_nuc.cfg` / `p2_z2n2_nuc.cfg` | Z2 park |
| `p2_a_z2n0.cfg` / `p2_a_z2n2.cfg` | Z2 + L=2 |
| `z6_n0_nuc.cfg` / `z6_n6_nuc.cfg` | Z6 nuclear |
| `z6_a_n0.cfg` / `z6_a_n6.cfg` | Z6 + L=6 atom |

### Score

```bash
python3 v75/analysis/score_pn_park.py <diag.tsv> <track.tsv> --Z 6 --N 6
```

### Gates

| ID | Pass |
|----|------|
| PN-S1 | n: \|Q_flux\|/\|Q_φ\| < 0.05 |
| PN-S2 | p: Q_flux ≈ Q_φ |
| PN-N | isotope Q_flux fixed under +N |
| PN-L | n_L = Z not A |

---

## C — Carbon nucleus (single-fabric Stage 2A)

**Package:** `work/C/carbon_nucleus_v0.md`

### Z-carbon production

```bash
# seeds
./v74/analysis/make_carbon_seeds.sh 192 36
# → c6_light_N192.sfa, c12_light_N192.sfa

# run (via scp-runner preferred)
# config: v74/cfg/c6_light.cfg
#   init_sfa=/space/scp/v74/seeds/c6_light_N192.sfa
#   N=192 L=36 T=300 g=0.05
```

| Run | Config | Expected |
|-----|--------|----------|
| **c6_light** | `v74/cfg/c6_light.cfg` | Q→~650 park; on-branch ~0.5% |
| **c12_light** | `v74/cfg/c12_light.cfg` | Q→~1411 super-crit evaporate |
| smoke | `v74/cfg/c6_light_smoke.cfg` | N=128 T=20 local |

### Geometry

- **c6:** octahedron, D=10, 6× light co-phase  
- **c12:** icosahedron, D=10, 12× light co-phase  

### Multi-fab nuclear (isotope cores)

Prefer **Z6N6** (`z6_n6_nuc.cfg`) for park; Z6N0 soft. See P package.

---

## L — Light opposite sector

**Package:** `work/L/light_sector_v0.md`

### Rules

```text
q_L = -1
# same-sign Noether ω on C and L (fabric q flips EM)
n_L = Z   # never A
private bag s_L — no merge with s_C
```

### Profiles

| Use | Profile | ω |
|-----|---------|---|
| Default cloud | `f_w146_g005.txt` | +1.46 |
| Soft hierarchy scout | `f_w1485_g005.txt` | +1.485 (contact-fragile) |
| Heavy H center | `f_w142_g005.txt` | +1.42 |

### Tools

| Tool | Role |
|------|------|
| `gen_qball_multi` ±ω | E-lite same-fabric ± force/orbit |
| `gen_qball_pair_boost` | ± with boosts |
| `gen_mf_pair_boost` | multi-fab C–L pair |
| `gen_mf_shell_orbit` | multi-L shell + optional vt |
| `gen_pn_core` | production Z + L=Z |

### Proven packages

| Package | Evidence |
|---------|----------|
| E-lite force D=16/20 | F8 PASS |
| E-lite orbit D=20 | F10 PARTIAL |
| MF force + isolation | F11–F13 PASS |
| B4 single-C + L6 | F15 PASS |
| Z2 + L2 | F18 PASS (−0.6% massL) |
| Z6 + L6 R=22 | F19 sector partial (−12.5%) |

### Grid safety (Z6)

L shell R=22 → **L=48**, damp_width=5. Do **not** use L=32 box.

---

## A — C₁₂-scale atom

**Package:** `work/A/c12_atom_package_v0.md`  
**Checklist:** `work/A/recipe_checklist_multifab.md`

### Assembly

```text
Core: Z≈6 p-analogs (+ optional N n-analogs) — prefer Z6N6 parked droplet
Cloud: n_L = Z (=6), ω_L=1.46, shell R=22
Physics: B2 freeze (above) + global block
Grid: N=192 L=48 damp_width=5
```

### Primary long-T candidate

```text
# seeds via gen_pn_core (geometry F19)
# config: v75/cfg/pn/z6_a_n6.cfg
#   init_sfa / init_sfa_Q / init_sfa_L
# T=400 smoke; T≥2000 long-T claim
```

### Pass bar

| Gate | Threshold |
|------|-----------|
| PASS_nuc | c_Q_park ≤ 0.15; c_Qem_park ≤ 0.20; Gauss floor |
| L hold | c_L ≤ 0.15; massL_end > 0.5 massL_0 |
| Isolation | no C–L bag merge |
| Isotope | Q_flux nuclear identical under +N; L track identical |

**Do not** treat net Q_flux drop alone as L death.

### Score + visual

```bash
python3 v75/analysis/score_pn_park.py diag.tsv track.tsv --Z 6 --N 6
bin/volview -snapshot N -view 9 -out charge.webp run.sfa
```

### Run matrix (when multi-fab GPU available)

| Priority | ID | T | Gate |
|----------|-----|---|------|
| 1 | z6_n6_nuc | 400 | PASS_nuc |
| 2 | z6_a_n6 | 400 | sector L + park |
| 3 | z6_n0_nuc | 400 | retune if soft |
| 4 | z6_a_n0 | 400 | isotope L |
| 5 | long_a_n6 | ≥2000 | ideal bar |
| 6 | h1_orbit | ≥2000 | P1 multi-rev |

Ops: scp-runner; **no sleep inside sim_exec**; poll `v75/analysis/remote_status_snippet.sh`.  
Wall: ~48 min / T=400 N=192 multi-fab V100-16GB.

### Product status

| Layer | Status |
|-------|--------|
| Recipe package | **READY** |
| Short-T F19 smoke | **PARTIAL** (park yes; PASS_atom false) |
| Ideal long-T C₁₂ | **NOT READY** — BLOCKED on GPU re-run / retune |

---

## W — World / relations (no sim)

| Artifact | Path |
|----------|------|
| Relation matrix | `PHYSICS_RELATIONS.md` |
| World freeze | `work/W/world_freeze_v0.md` |

---

## Runner sketch

```text
sim_setup(executor="remote")   # or local
sim_build(sources=[scp_sim.cu or .c, sfa.h])
sim_run(config=<cfg content>, id=<run_id>)
sim_run_status(id=...)
sim_download → /space/scp/v78/...
# score locally
```

---

## Data roots

| Path | Content |
|------|---------|
| `/space/scp/v74/` | c6/c12 seeds + full SFA |
| `/space/scp/v75/pn/` | F17–F19 multi-fab |
| `v75/results/pn_z6/` | score JSON |
| `v74/results/*.png` | nuclear slice panels |
| `images/C6_atom_sheet.png` | F19 atom sheet |

---

## Explicit non-recipes

- Kernel multi-fabric redesign without auth  
- Fake neutron via flavor multiplet alone  
- Calling Q→Q_max evaporation “parked A=12”  
- Z6 L@R=22 in L=32 box  
