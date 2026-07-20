# work/P — Proton / neutron package

**Owner:** Agent P  
**Log:** [`../../logs/P_proton_neutron.log`](../../logs/P_proton_neutron.log)  
**Checkpoint:** **CP-P ADOPT** (see package)

| Deliverable | Path |
|-------------|------|
| **P/N package freeze** | [`pn_package_v0.md`](pn_package_v0.md) |
| Evidence (v75) | `v75/FINDINGS.md` F17–F19 · `v75/STATE.md` · `v75/PN_EXPERIMENT.md` |
| Configs | `v75/cfg/pn/` |
| Scorer | `v75/analysis/score_pn_park.py` |
| Seed | `bin/gen_pn_core` ← `sfa/seed/gen_pn_core.c` |

**One-liner:** proton = C+Q, neutron = C-only (B2); \(q_C=0,q_Q=+1,q_L=-1\); \(Q_{\mathrm{em}}\propto Z\).
