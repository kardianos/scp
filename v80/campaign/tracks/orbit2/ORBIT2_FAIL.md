# Orbit2 — Gate O FAIL (closed early)

**Date:** 2026-07-19  
**Instance:** 45320919 (destroyed mid soft job)  
**Decision:** Assume **O FAIL**. Soft \(T=4000\) not taken to completion after ledger death.

## Jobs

| Job | Status | Outcome |
|-----|--------|---------|
| hard `vt=0.50 T=1200` | COMPLETE | L mass **192→0 by t~170**; ~0.44 rev while L alive; **E drift −53%** |
| soft `vt=0.03 T=4000` | **KILLED ~78%** | Healthy to \(t\sim500\text{–}1000\); by \(t\sim2000\): \(Q\sim82\), \(E_{\mathrm{em}}\sim0.035\); by \(t\sim3000\): \(Q\sim9\), \(E_{\mathrm{em}}\sim10^{-4}\) |

## Scorecard

| Gate | Status |
|------|--------|
| F Force (short rest) | **PASS** (prior tracks campaign) |
| O Orbit / long hold | **FAIL** |
| Multi-L / Stage 4 | **blocked** — do not re-park |

**Verdict:** product multi-fab C–L bound system fails on long \(T\). Force without durable bound state.  
**Next:** leave multi-fab product ladder; Stage-3 redesign (kernel authorized) capturing v80 free/bound thesis — design feedback first.

## Artifacts

- Hard track (f16-fixed): `results/mf_orbit_R16_vt0p50_T1200_track.tsv`
- Soft partial diag: `results/mf_orbit_R16_vt0p03_T4000_diag.tsv`
- Tracker fix: `sfa/analysis/mf_pair_track.c` (f16 support)
