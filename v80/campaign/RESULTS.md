# v80 campaign RESULTS — multi-fab mechanics graph (overnight V100)

**Queue:** 2026-07-19 02:05 → 06:28 UTC (~4.4 h wall)  
**Instance:** Vast 45259142, V100-32GB  
**Jobs:** 16/16 OK · 0 FAIL · `FORCE_OK=1` · `PAIR_OK=1` · `QUEUE_DONE`  
**Artifacts:** `campaign/results/*_{diag.tsv,runlog}` (output SFAs pruned on remote after each job)

---

## Goal function score

| Component | Weight | Score | Evidence |
|-----------|--------|------:|----------|
| **S_force** | 0.25 | **0.55** | All force jobs COMPLETE; elite opposite net \(Q\approx0\), same \(Q\sim230\); \(E_{\mathrm{em}}\)(same)> \(E_{\mathrm{em}}\)(opp) at D=12; **no pair-track \(a_{\mathrm{rel}}\)** (SFAs pruned) |
| **S_pair** | 0.25 | **0.65** | C–L rest T=400: Q held, \(E_{\mathrm{em}}\) held; head-on: \(E_{\mathrm{em}}\) 0.74→0.31 (interaction); rest E drift +14% soft |
| **S_orbit** | 0.20 | **0.50** | All 3 vt COMPLETE; low-vt (0.03/0.05) mild; **vt=0.08 E drift −27%** / \(E_{\mathrm{em}}\) drop; multi-rev not measured without tracks |
| **S_Lhold** | 0.15 | **0.70** | S4 T=800: \(E_{\mathrm{em}}\) **held** 0.71→0.75 (≠ v79 atom collapse); \(Q_\phi\) held; **\(Q_{\mathrm{flux}}\) 114→4.7**, \(Q_{\mathrm{core}}\) 114→0.46 (window/COM issue) |
| **S_ledger** | 0.10 | **1.00** | All runs `gauss_max` ~1e-13 floor |
| **S_morph** | 0.05 | **0.25** | No volview package this night |
| **G** | 1.00 | **0.62** | |
| hard_fail | | **false** | |

\[
G = 0.25(0.55)+0.25(0.65)+0.20(0.50)+0.15(0.70)+0.10(1.00)+0.05(0.25) = \mathbf{0.62}
\]

### Threshold decision

| Band | This night |
|------|------------|
| \(G \ge 0.55\) continue product hydrogenoid | **YES** |
| \(0.30\le G < 0.55\) retune only | — |
| \(G < 0.30\) soft-kill multi-fab atoms | — |

**Morning verdict:** Multi-fab **product path still alive** at the force/pair/ledger level. Orbit is **partial** (high-\(v_t\) fragile). Not a cold multi-rev atom claim. Representation (v80 thesis) remains separate daytime work.

---

## Run inventory

| Step | Job | Wall | Gauss end | Notes |
|------|-----|-----:|-----------|-------|
| S0 | smoke_mf | 1.7 min | 9.0e-14 | PASS |
| S1 | mf_force D12/16/20 | ~7.7 min ea | ~1e-13 | B1 C–L; Q~114 |
| S1 | elite_opp D12/16/20 | ~4 min ea | ~1e-13 | net Q≈0 |
| S1 | elite_same D12/16/20 | ~4 min ea | ~1e-13 | Q~230 |
| S2 | mf_pair_rest_D20 | 15.3 min | 9.7e-14 | E drift +14% |
| S2 | mf_headon_D20 | 10.4 min | 9.7e-14 | \(E_{\mathrm{em}}\) halves |
| S3 | orbit vt0.03 | 20.4 min | 9.6e-14 | mild E drift +2.8% |
| S3 | orbit vt0.05 | 20.4 min | 9.0e-14 | mild |
| S3 | orbit vt0.08 | 20.5 min | 9.3e-14 | **E −27%** |
| S4 | mf_h1_orbit_N192 | **87.7 min** | 9.5e-14 | Q hold; flux/core window collapse |

---

## Physics takeaways

1. **Ledger integrity** — multi-fab B1 on V100 is numerically healthy for this matrix (Gauss floor).  
2. **Charge bookkeeping** — E-lite opposite vs same is correct (neutral vs \(2Q\)).  
3. **E_em survival** — long-T S4 does **not** reproduce v79 atom \(E_{\mathrm{em}}\to0\); hydrogenoid-class seed holds EM energy.  
4. **Diagnostics caution** — \(Q_{\mathrm{core}}/Q_{\mathrm{flux}}\) can collapse while \(Q_\phi\) holds (COM/window); do not equate flux drop with charge annihilation without tracks.  
5. **Orbit** — low tangential speed gentler; \(v_t=0.08\) radiates hard. Need `mf_pair_track` / sfa tracks on **kept** SFAs next time.  
6. **v79 Z6+L6 null** remains: this night did **not** re-test hand-placed Z=6 shells.

---

## Recommended next (product)

1. Re-run one D=20 force pair with **SFA kept** + `mf_pair_track` / `sfa_qball_track` → real \(a_{\mathrm{rel}}\).  
2. Orbit scan with track TSVs; publish multi-rev count.  
3. Hydrogenoid T-extend only if tracks show capture class.  
4. Parallel: v80 representation toy (CPU) — not blocked by this G.

## Teardown

Instance **45259142** destroyed after download verification (see STATUS / poll log).
