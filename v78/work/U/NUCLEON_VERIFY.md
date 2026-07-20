# Nucleon local verify (U R4)

**Date:** 2026-07-19  
**Purpose:** Confirm light nucleon template numbers used by ladder/recipes without re-running full 3D GPU.

## 1. Gauged branch (authoritative)

From `v69/theory/gscan.tsv` at **g=0.05**, **ω=1.460000**:

| Quantity | Value |
|----------|-------|
| f0 | 0.606126 |
| Q | **114.13** |
| E_matter | 172.90 |
| E_field | 0.432 |
| E_total | 173.33 |
| E_over_mQ | 1.01246 |
| dQ/dω sign | **−1** (VK stable) |
| r_half | **2.630** |
| weff0 | 1.449003 |

E/Q on matter+field total: \(E_{\mathrm{tot}}/Q = 173.33/114.13 \approx 1.5187\) — matches CARBON_MAP light nucleon row.

## 2. Profile file

| Path | Check |
|------|-------|
| `v74/profiles/f_w146_g005.txt` | Present; 2-col r,f |
| f(0) | 6.061256640e-01 ≈ gscan f0 |
| r max | 60.0; f→0 at large r |
| Convention | f(r) = **per-component** amplitude (v72+) |

## 3. Window (g=0.05)

| Edge | ω | Note |
|------|---|------|
| Thick-wall / Q_max | 1.40624 | Q_max = **921** |
| Light default | 1.46 | Q≈114 |
| Lightest VK-stable | ~1.4825 | Q_N^min≈90 |
| Unstable thin wall | ≳1.485 | do not free-seed |

**Charge-budget theorem:** \(12 \times Q_N^{\min} > Q_{\max}\) ⇒ free A=12 always super-critical at g=0.05.

## 4. Standard package (frozen)

```
complex_phi=1 complex_gauge=1 g_gauge=0.05
m=1.5 m_theta=1.6 eta=0
mu=-41.345 kappa=50
Vt(s)=(mu/2) s/(1+kappa s), s=∏|Φ_a|²
```

## 5. Local tool note

- Ungauged `v66/radial_qball` exists for ω-window smoke; **production light template is gauged gscan + `f_w146_g005.txt`**.
- Seed writers: `gen_qball_multi`, `gen_qball_boost`, `gen_pn_core` — read profile, do not re-shoot unless g changes.

## 6. Verdict

**PASS** — light nucleon template bookkeeping consistent across gscan, profile, and v74/v75 usage.
