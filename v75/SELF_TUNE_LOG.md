# Self-tune Stage 1 — living log

**Status: CLOSED — DEFINITIVE_SUCCESS**  
**Winning θ:** B4_full — single heavy C (ω=1.42) + 6×L shell (ω=1.46, R=18)

| Item | Value |
|------|--------|
| Design | `SELF_TUNE_C.md` |
| Controller | `analysis/self_tune_controller.py` |
| Ledger | `/space/scp/v75/self_tune/ledger.jsonl` |
| VERDICT | `/space/scp/v75/self_tune/VERDICT.md` |
| GPU | `v75st` V100-32GB |

**Scorecard note:** FAIL_Q often means seed→park \(\Delta Q\) vs seed \(Q_0\), not
“L died.” PASS requires cost≤0.15, c_Q≤0.15, c_L≤0.15, gauss floor.

---

## Full trials (T=400) — closed

| ID | θ | cost | c_Q | c_L | massL | Qc | Physics read |
|----|---|------|-----|-----|-------|-----|--------------|
| B1_full | Z6 R18 | 0.259 | 0.296 | 0.163 | 472→395 **−16%** | 908→639 | L shed + park |
| B1a_full | Z6 R22 | 0.251 | 0.299 | **0** | 472→474 **stable** | 908→636 | R fixes L |
| B1b_full | Z6 softL R22 | 0.203 | 0.300 | 0.024 | 366→358 −2% | 908→636 | best Z6 |
| B2_full | Z4 R18 | 0.161 | 0.250 | 0.012 | 314→311 −1% | 562→422 | best multi-ball until B3 |
| B2a_full | Z4 R22 | 0.187 | 0.252 | **0** | 314→315 **stable** | 562→421 | L perfect |
| B3_full | Z2 R18 | 0.109 | 0.191 | **0** | 157→158 **stable** | 247→200 | PARTIAL |
| B3a_full | Z2 R22 | 0.081 | 0.192 | **0** | 157→158 **stable** | 247→199 | best multi-ball PARTIAL |
| **B4_full** | **Z1+L6 R18** | **0.112** | **~0** | **0** | 472→474 **stable** | **315→315** | **PASS** — stop_success |

Tree stopped on first PASS (B4a, B5 not run).

## Winning θ (freeze)

```
n_C=1, omega_C=1.42, n_L=6, omega_L=1.46, R_shell=18
g_gauge=0.05, N=192, L=48, B1 lock, q_Q=+1, q_L=-1
```

| Metric | t=0 | t=400 |
|--------|-----|-------|
| massC | 222.1 | 225.5 |
| massL | 472.1 | 473.5 |
| Qc | 315.4 | 315.4 |
| Ql | 689.3 | 689.3 |
| E | 1505 | 1673 (+11% → c_E, still PASS) |
| gauss_max | — | 9.4e-14 |

## Claims

1. **DEFINITIVE SUCCESS** under Stage-1 scorecard at soft θ B4.
2. **L packaging works** with multi-L shell around a single heavy center.
3. **Multi-ball nuclei (Z2–Z6)** never reached PASS (seed→park c_Q); single-C does.
4. **Larger R** fixes L erosion on multi-ball setups.
5. **Option C not required** for this scale claim (single heavy + L cloud).
6. **Caveat:** concentric rest state (COM D~0), not a demonstrated bound *orbit*; energy drifts +11%; not multi-Z “carbon” nucleus.

---

## Append-only chronology

| UTC | Event |
|-----|--------|
| 14:02 | Campaign start v75st |
| 15:07 | B1_full FAIL_Q 0.259 massL −16% |
| 16:12 | B1a_full FAIL_Q 0.251 **massL stable** |
| 17:12 | B1b_full FAIL_Q 0.203 |
| 18:22 | B2_full FAIL_Q 0.161 |
| 19:22 | B2a_full FAIL_Q 0.187 massL stable |
| 20:30 | B3_full PARTIAL 0.109 |
| 21:36 | B3a_full PARTIAL **0.081** (best multi-ball) |
| 22:00 | B4_screen PARTIAL 0.032 c_Q=0 |
| **22:43** | **B4_full PASS 0.112 → DEFINITIVE_SUCCESS** |
