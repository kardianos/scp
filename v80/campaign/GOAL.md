# Campaign gates — final scoreboard (AND layers, no composite \(G\))

**Updated:** 2026-07-19 (tracks campaign complete)  
**Rulebook:** `SCORECARD_AND.md`  
**S4 window gate:** `S4_COM_WINDOW_ANALYSIS.md` — **PASS**  
**Tracks detail:** `tracks/TRACK_RESULTS.md`

Weighted \(G = 0.62\) is **historical only**. Not used to advance stages.

---

## Layer 1 — Integrity

| Gate | Status | Note |
|------|--------|------|
| L0 Gauss | **PASS** | Floor all runs |
| L1a Diags | **PASS** | All force/elite/orbit diags |
| L1b Tracks | **PASS** | keep-SFA + pair tracks for F/R/O claims |

---

## Layer 2 — Stage-3 mechanics

| Gate | Status | Note |
|------|--------|------|
| F Force | **PASS** | Opposite C↔L attract all \(D\in\{8,12,16,20,24\}\); \(|a_{\mathrm{rel}}|\) monotone ↓ with \(D\) |
| R Same-sign | **SOFT** | \(D=20\) repel vs opposite attract; \(D=12,16\) bag attract |
| P Pair dwell | **SOFT** | Rest holds; E drift (overnight) |
| O Orbit | **FAIL** | Short \(T=400\) arcs held; long \(T\sim4000\) soft: \(Q\)/\(E_{\mathrm{em}}\) death after \(t\sim2\times10^3\); hard \(v_t=0.5\): L shreds (~0.4 rev). See `tracks/orbit2/ORBIT2_FAIL.md` |
| H Long-T hold | **PASS** | S4 COM window analysis |

---

## Layer 3 — Stage-4 shell

| Gate | Status | Note |
|------|--------|------|
| N Multi-L | **N/A** | Do not re-park Z6+L6 |

---

## Verdict

```text
VERDICT = STOP_AND_RECONSIDER   # product multi-fab atom ladder
  F PASS — rest attraction real
  O FAIL — no durable C–L bound system on long T
  do not: Z6+L6 re-park; more soft/hard orbit grids
  next: leave multi-fab product path; Stage-3 redesign (kernel auth given)
        capture v80 free/bound + locks thesis in numeric form
```

Historical composite (do not use for decisions): \(G = 0.62\).
