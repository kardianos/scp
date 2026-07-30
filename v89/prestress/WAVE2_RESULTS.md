# PRESTRESS Wave 2 results — plasticity + topology under κ (2026-07-30)

**Kernel:** integer default `cellfab` (`ledger_mode=3`).  
**Foam:** `foam_s20260727.tsv`, seed 20260727, L=24.  
**Laws:** `battery/laws_V2g.cfg` + apparatus extras (`kappa_plast`, optional `tau_harden`).  
**Fleet:** `prestress/fleet.go` wave 2 (max 2 concurrent × 4 threads).  
**Scores:** `runs/SCORES.tsv` + `runs/WAVE2_SCORED.tsv`.

Plasticity is **live**: `# PLAST` shows thousands of links moved per diag
tick, `max_cum` ~0.3–0.5 by late times; `kappa_plast=0` remains the exact
legacy path (W1). Integer ledger: `max_sum_err=0` on all scored runs.

---

## 1. What Wave 2 was for

Pre-registration (`LEDGER.md` Waves 2–4 table):

| arm | pre-reg |
|-----|---------|
| PLAST-1 cube κ=1 vs W1 κ=0 | gates mean ≥0.9 by t≈200; t_death > W1 908 or lower late leak; prize C1 |
| PLAST-1 κ=0.5 | mid-window; slower anneal OK if gates lock |
| PLAST-3 c1+κ | parasites darken / self-seal vs frozen c1 |
| Topo+κ (ring8, c8, tube, chords, hex) | plast helps best topologies most |
| Harden τ=50 | less late gate churn; mass tighter (P19 prep) |

---

## 2. Scored table (measured x50)

Definitions (same as Wave 1):

- **x50** = Emfree(t≈50) / (nv · cap), cap=2.5  
- **Law A** = 274 · (x50/0.0617)^1.066  
- **t_death** = first t with largest lump mass < 0.3·m0  
- **c_eff** = cap · (x50 − x_skirt) / t_death  
- **g50 / g200** = live `# NETG` mean gate at t≈50 / t≈200  

| run | role | nv | x50 | Law A | t_death | vs W1 | t_d/LawA | g50 | g200 | verdict |
|-----|------|---:|----:|------:|--------:|------:|---------:|----:|-----:|---------|
| **c2_plast1_k05** | PLAST-1 cube κ=0.5 | 8 | 0.389 | 1952 | **1287** | 908→**1.42×** | 0.66 | 0.926 | 0.281 | **MID** best cube |
| c2_plast1_k1 | PLAST-1 cube κ=1 | 8 | 0.389 | 1952 | 988 | ~1.09× | 0.51 | 0.934 | 0.249 | MID |
| c2_plast1_ctrl_k1 | ctrl + κ=1 | 8 | 0.389 | 1951 | **1071** | 698→**1.53×** | 0.55 | 0.843 | 0.135 | MID |
| c2_harden | cube κ=1 τ=50 | 8 | 0.389 | 1952 | 1017 | ~1.12× | 0.52 | 0.934 | 0.206 | MID |
| **c1_frozen** | PLAST-3 ctrl | 8 | 0.376 | 1884 | **1013** | — | 0.54 | 0.311 | 0.250 | MID |
| **c1_plast3_k1** | PLAST-3 +κ | 8 | 0.387 | 1942 | **30** | 1013→**0.03×** | 0.02 | 0.797 | 0.347 | **EARLY fail** |
| **c8_k1** | topo +κ | 12 | 0.405 | 2038 | **1269** | 449→**2.83×** | 0.62 | 0.271 | 0.331 | MID **big lift** |
| **ring8_m3_k1** | topo +κ | 8 | 0.545 | 2796 | **14** | 1631→**0.01×** | 0.005 | 0.100 | 0.074 | **EARLY kill** |
| ring8_harden | +κ +τ | 8 | 0.545 | 2794 | **14** | same | 0.005 | 0.099 | 0.187 | EARLY kill |
| **tube_k1** | topo +κ | 12 | 0.431 | 2178 | **1613** | 1042→**1.55×** | **0.74** | 0.169 | 0.213 | **ON_BAND** longest W2 |
| chords_k1 | topo +κ | 20 | 0.528 | 2705 | **127** | 630→**0.20×** | 0.05 | 0.437 | 0.253 | EARLY hurt |
| hex_frozen | topo κ=0 | 12 | 0.480 | 2439 | 986 | — | 0.40 | 0.578 | 0.305 | EARLY |
| hex_k1 | topo +κ | 12 | 0.486 | 2474 | 1160 | 986→1.18× | 0.47 | 0.879 | 0.168 | EARLY modest |

No run reaches the exception class (≥2.3× Law A) or C1 plateau. Absolute
Law A ratios remain early-to-mid (0.005–0.74), same calibration flag as W1.

---

## 3. Decisions vs pre-registration

### PLAST-1 — cube anneals (primary claim)

| Claim | Result |
|-------|--------|
| Live mean gate ≥0.9 by t≈200 | **Fail.** κ=0.5/1.0 both peak early (g50≈0.93) then collapse (g200≈0.25–0.28) |
| t_death > W1 908 | **κ=0.5 yes** (1287); **κ=1 marginal** (988) |
| Prize C1 (parked / skin) | **No** |

**Verdict:** Plasticity **moves geometry** (moved~10k links by t=200,
max_cum~0.3) and briefly **opens gates**, but does **not lock** a
consonant cube. Mild life extension at mid κ; full κ≈1 is not better
than 0.5 on this object. Retune still matters less under κ than the
control gap under freeze (ctrl+κ 1071 vs retune+κ 988 — phases scramble
with d-motion).

### PLAST-3 — parasite self-seal

| Claim | Result |
|-------|--------|
| Parasites darken / self-seal | **Fail catastrophically** — t_death **30** vs frozen **1013** |
| Early gate boost | Real (g50 0.80 vs 0.31) then structure dies immediately |

**Verdict:** On the naive c1 cube, κ=1 **accelerates death**. The reduced
numpy self-seal picture does **not** survive the full kernel on this seed.
Do not treat “plast seals parasites” as measured fact.

### Topo under plasticity

Pre-reg: “plast helps best topologies most.”

| object | W1 | W2+κ | ratio |
|--------|---:|-----:|------:|
| tube | 1042 | **1613** | **+55%** best absolute life |
| c8 ring12 | 449 | **1269** | **+183%** largest relative lift |
| hex | 986 | 1160 | +18% |
| cube (κ=0.5) | 908 | 1287 | +42% |
| chords | 630 | 127 | **−80%** |
| ring8_m3 | **1631** | **14** | **−99%** W1 champion killed |

**Verdict:** Help is **topology-dependent**, not monotone in “best W1.”
Chirality/winding objects (tube, c8) gain; the high-load back-gate ring and
chord-heavy graph are **destabilized**. Plasticity is not a universal
upgrade — it is a mechanism that can open a longer corridor **or** pull
struts through the death radius faster.

Tube under κ reaches **0.74× Law A** (ON_BAND under ×/1.5) — closest to
the load line in either wave — but still far from the ≥4600 exception
class.

### Harden (τ=50)

| Claim | Result |
|-------|--------|
| Hardening engages | **Yes** — cube `hard_n` ~7600 by t=3000; first hard links ~t=71 |
| Less late gate churn / longer life | **No clear win** — cube harden 1017 ≈ κ=1 988 |
| Saves ring8 | **No** — dies at 14 before harden can pin |

**Verdict:** Instrument works; on these short-lived objects it does not
change the death story. Remains relevant for W4 P19-lite **if** something
lives past the anneal window.

---

## 4. Cross-cutting conclusions

1. **No C1.** Exact mass (M-R1) still blocked.  
2. **Plasticity is real and non-null** — geometry moves; lifetimes shift
   by factors of 0.01–2.8 depending on net.  
3. **The reduced-model promise (gates lock ≥0.9, self-seal, pin to mass)
   does not hold at rate level** on foam 20260727 under provisional κ*∈{0.5,1}.  
4. **Best W2 object:** `tube_k1` (t_death 1613, 0.74×LawA). Best relative
   rescue: `c8_k1` (2.8× W1).  
5. **Worst anti-result:** `ring8_m3` and `c1` under κ=1 — plasticity can
   destroy the W1 longest-lived and the parasite-test cube.  
6. **κ pin hunt:** mid κ (0.5) beat κ=1 on the cube; provisional κ*=1 is
   **not** justified as default for all nets. Topology × κ grid is still
   open if anyone reopens pin search.  
7. **Next mechanism tests** are Wave 3 discriminators (flight load,
   m-ladder, negatives) and Wave 4 harden/m-family on survivors — not
   more frozen topology hunting, and not “turn κ on for everything.”

---

## 5. Artifacts

| path | content |
|------|---------|
| `runs/SCORES.tsv` | score_net rows (all waves) |
| `runs/WAVE2_SCORED.tsv` | x50, Law A, g50/g200, plast, verdicts |
| `runs/w2_*.log` | diags + `# PLAST` + `# LEDGER` |
| `runs/w2_*.{slab,3d}.png` | last-frame renders |
| `fleet.go` | orchestrator (waves 2–4 continuing) |

---

## 6. Next

- **Wave 3:** flight233, ring m-ladder, hopf/mobius/octa/free, torus/truncocta negatives ±κ.  
- **Wave 4:** P19-lite harden T=4000; P20 m-family under harden; tube harden re-test.  
- Do **not** promote κ_plast into `laws_V2g` (apparatus only).  
- Optional: κ scan on tube/c8 only (winners), not on ring8/c1.
