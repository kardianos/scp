# PRESTRESS Wave 3 results — discriminators / package / free-search (2026-07-30)

**Kernel:** integer default `cellfab` (`ledger_mode=3`).  
**Foam:** `foam_s20260727.tsv`, seed 20260727, L=24.  
**Laws:** `battery/laws_V2g.cfg` (+ `kappa_plast=1` on two arms).  
**Fleet:** `prestress/fleet.go` wave 3 (completed with waves 2–4).  
**Scores:** `runs/SCORES.tsv` + `runs/WAVE3_SCORED.tsv`.

---

## 1. What Wave 3 was for

| arm | pre-reg |
|-----|---------|
| flight233 | less early shed than c8 at higher x (load down) |
| ring9/10/12 m-ladder | odd-N / m-family ranking (frozen ladder) |
| hopf | linked-ring package |
| mobius / octa | negatives die fast / frustrated |
| free0 / free1 | free formfind vs engineered |
| mobius_k1 / hopf_k1 | negatives under plast |
| torus / truncocta | infeasible-class still run as negatives (T=2000) |

---

## 2. Scored table

| run | role | nv | x50 | Law A | t_death | t_d/LawA | gmin | par | verdict |
|-----|------|---:|----:|------:|--------:|---------:|-----:|----:|---------|
| **free1** | free-search | 9 | 0.446 | 2255 | **2572** | **1.14** | 0.765 | 2 | **ON_BAND** longest in program |
| **ring12_m5** | m-ladder | 12 | 0.479 | 2437 | **1940** | **0.80** | 0.001 | 2 | **ON_BAND** best engineered ring |
| torus3x8 | infeas (T=2000) | 24 | 0.792 | 4161 | 1672 | 0.40 | 0.003 | 29 | EARLY (heavy, parasite-rich) |
| truncocta | infeas (T=2000) | 24 | 0.492 | 2504 | 1551 | 0.62 | 0.035 | 11 | MID |
| ring12_m6 | m-ladder | 12 | 0.400 | 2009 | 1402 | 0.70 | 0.200 | 1 | ON_BAND |
| hopf | package | 24 | 0.420 | 2117 | 1093 | 0.52 | 0.000 | 7 | MID |
| **mobius_k1** | neg+κ | 12 | 0.331 | 1642 | **1074** | 0.65 | 0.003 | 8 | MID (plast helps) |
| free0 | free-search | 8 | 0.374 | 1868 | 855 | 0.46 | 0.809 | 4 | EARLY |
| ring10_m4 | m-ladder | 10 | 0.577 | 2972 | 846 | 0.28 | 0.000 | 1 | EARLY |
| ring9_m3 | m-ladder | 9 | 0.565 | 2903 | 761 | 0.26 | 0.000 | 2 | EARLY |
| octahedron | negative | 6 | 0.301 | 1485 | 718 | 0.48 | 0.000 | 1 | EARLY |
| mobius | negative | 12 | 0.335 | 1665 | 521 | 0.31 | 0.003 | 8 | EARLY |
| flight233 | light load c8 | 12 | 0.203 | 976 | **347** | 0.36 | 0.400 | 1 | EARLY **worse than W1 c8** |
| hopf_k1 | hopf+κ | 24 | 0.426 | 2151 | **76** | 0.04 | 0.000 | 7 | EARLY **plast kill** |

Integer ledger: `max_sum_err=0` on scored full runs.

---

## 3. Decisions vs pre-registration

### flight233 (P18-lite)

| Claim | Result |
|-------|--------|
| Lighter flight load sheds less / lives longer than heavy c8 | **Fail** — t_death **347** < W1 c8 **449** and ≪ W2 c8_k1 **1269** |
| Lower x50 helps | x50=0.20 is light; still dies early (0.36×LawA) |

**Verdict:** Load-down alone is not a rescue. The flight-corrected seed is not an exception path on this foam.

### Ring m-ladder (frozen)

Ranking by t_death: **m5 (1940) > m6 (1402) > m4 (846) > m9_m3 (761)**.  
Not monotone in m or N. ring12_m5 is the **best pure engineered ring** in the program (0.80×LawA). W1 ring8_m3 (1631) still sits near the top of the frozen ladder but was killed under κ.

### free-search

**free1** is the **longest-lived object in the entire W1–W4 campaign** (2572, **1.14× Law A** — ON_BAND, only run above the Law A point estimate). Seed gates already good (min 0.76 / mean 0.95). free0 is ordinary (855). Free formfind can beat hand-engineered seeds; one lucky draw is not a species.

### Negatives

mobius / octa die mid-early (521 / 718) — frustrated but not instant death. Under κ: **mobius helped** (521→1074), **hopf destroyed** (1093→76). Same topology-dependent plast rule as Wave 2.

### Infeasible class (torus, truncocta)

Still live hundreds–thousands of t.u. with many parasites (torus par=29). “Infeasible” formfind is not “dies at once” under the top-lump death definition. They do **not** look like particles (high x50, rough_share, no plateau).

---

## 4. Cross-cutting

1. **Longest lives are not C1.** free1 and ring12_m5 still die; no parked mass.  
2. **Gate quality still ≠ lifetime** (flight233 mid gates, short life; ring12_m5 gmin≈0, long life).  
3. **Plasticity remains double-edged** on package objects (mobius↑, hopf↓).  
4. **Best frozen topology to carry into harden/P20:** ring12_m5 / ring12_m6 / free1 — but W4 shows harden can **kill** m5 (see Wave 4).

---

## 5. Artifacts

| path | content |
|------|---------|
| `runs/WAVE3_SCORED.tsv` | x50 / Law A / verdicts |
| `runs/w3_*.log` | full diags |
| `WAVE2_RESULTS.md` | plasticity context |
