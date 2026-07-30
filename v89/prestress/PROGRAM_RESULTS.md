# PRESTRESS full program results — Waves 1–4 (2026-07-30)

Single-foam campaign on seed **20260727**, integer ledger `cellfab` mode 3,
laws_V2g + apparatus (edge_sink, lump_diag; plast/harden only where noted).
Orchestrator: `fleet.go`. Wave writeups: `WAVE1_RESULTS.md` …
`WAVE4_RESULTS.md`. Scores: `runs/SCORES.tsv` + `WAVE*_SCORED.tsv`.

---

## Campaign questions → answers

| Question | Answer |
|----------|--------|
| Consonant skin on frozen foam (cube ≥4700 / load line)? | **No.** Cube 908 (0.47×LawA). |
| Tube exception (≥4600)? | **No.** Best 1613 under κ. |
| Plasticity unlocks particles? | **No.** Real geometry motion; lifetime ×0.01–2.8 by topology; **no C1**. |
| Gates lock under κ by t≈200? | **No.** Early peak then collapse. |
| Parasite self-seal (PLAST-3)? | **No** — catastrophic early death. |
| Exact mass M-R1 / P19-lite? | **No** on this foam+seeds. |
| m-family coexistence P20? | **No** under harden. |
| Anything beat Law A? | **free1 only** (1.14×); still dies at 2572. |

---

## Lifetime leaderboard (t_death)

| rank | run | wave | t_death | note |
|-----:|-----|------|--------:|------|
| 1 | **w3_free1** | 3 | **2572** | free formfind; 1.14×LawA |
| 2 | w3_ring12_m5 | 3 | 1940 | best engineered ring (frozen) |
| 3 | w3_torus3x8 | 3 | 1672 | infeas class, T=2000, 29 parasites |
| 4 | w1_ring8_m3 | 1 | 1631 | W1 longest; **killed by κ** |
| 5 | **w2_tube_k1** | 2 | **1613** | best plast help among W1 nets |
| 6 | w3_truncocta | 3 | 1551 | infeas, T=2000 |
| 7 | w4_tube_harden | 4 | 1437 | harden < κ-only |
| 8 | w3_ring12_m6 / w4 m6 | 3/4 | 1402/1401 | m6 stable under harden |
| … | w2_c8_k1 | 2 | 1269 | 2.8× W1 c8 |
| … | w4_p20_ring12_m5 | 4 | **26** | best frozen ring **killed by harden** |
| … | ring8 under κ/harden | 2/4 | **14** | W1 champion destroyed |

---

## Mechanism summary

1. **Load line / vacuum bleed** still frames death directionally; absolute Law A ratios often early (death definition / integer / apparatus calibration open).  
2. **Seeded gate quality does not set lifetime** (repeated across waves).  
3. **Plasticity is topology-conditioned:** helps tube/c8/mobius; destroys ring8/c1/hopf/m5 under harden.  
4. **Harden engages** (hard_n large) but does not create mass attractors here.  
5. **Conservation:** integer `max_sum_err=0` throughout scored runs.

---

## What is closed vs open

**Closed on foam 20260727 (this apparatus):**
- Frozen-foam “gates alone → skin”
- Universal κ upgrade
- PLAST-3 self-seal as stated
- Single-foam P19-lite / P20 under κ=1, τ=50

**Open / next (if mass program continues):**
- Multi-foam formfind ensemble (true P19)
- κ pin scan **only** on free1 / ring12_m5 / tube (not blanket κ=1)
- Death-definition A/B vs corpus M_sum (Law A calibration)
- P21 kicks deferred until something lives past ~5k
- EMF Mode P continues separately (`emf/`)

---

## Job counts

| wave | jobs | status |
|------|-----:|--------|
| 1 | 7 | complete + scored |
| 2 | 13 | complete + scored |
| 3 | 14 | complete + scored |
| 4 | 9 | complete + scored |
| **total** | **43** | fleet log: `[fleet] complete wave=2-4` |
