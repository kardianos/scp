# v80 overnight V100 campaign — multi-fab mechanics graph

**Started:** 2026-07-19  
**Hardware:** Vast V100-32 (`ssh7.vast.ai:19142`, contract 45259142)  
**Kernel:** existing multi-fab `scp_sim_cuda` (no kernel edits)  
**Data:** seeds+results on `/space/scp/v80/`; diags also under `v80/campaign/*/results/`

## Scoring (AND layers — no composite \(G\))

**Standing rule:** `SCORECARD_AND.md` (integrity ∧ Stage-3 mechanics ∧ Stage-4 shell).

Historical overnight used weighted \(G=0.62\); **retired as decision device**.

**Pre-force gate (done):** `S4_COM_WINDOW_ANALYSIS.md` — \(Q_{\mathrm{core}}\)/\(Q_{\mathrm{flux}}\) collapse is fixed-center windowing; global \(Q_\phi\) and \(E_{\mathrm{em}}\) hold.

**Next:** force grid with kept SFAs + `mf_pair_track` (gates F, R, L1b).

## Steps (each folder)

| Dir | Step | Grid | T | Est wall |
|-----|------|------|---|----------|
| `S0_smoke/` | Multi-fab load+Gauss+T=40 | N=128 L=32 | 40 | ~5–10 min |
| `S1_force/` | B1 C–L force D=12,16,20 + E-lite same D | N=128 | 100 | ~1–1.5 h |
| `S2_ql_pair/` | C–L rest longer + head-on scout | N=128 | 400 | ~1–1.5 h |
| `S3_orbit/` | shell orbit n_L=1 vt scan | N=128 | 400–600 | ~1.5–2 h |
| `S4_hydrogenoid/` | best seed long-T N=192 **or** L-retune | N=192 | 800 | ~2 h |

**Side paths:** encoded in `scripts/run_queue.sh` (force dead → skip orbit, do nucleus ledger fill).

## Ops

- Remote queue: `scripts/run_queue.sh` (nohup on GPU)  
- Poll: every **15 min** → status + rsync diags  
- Teardown Vast when campaign complete and downloads verified  
- Large SFA → `/space/scp/v80/results/` only if needed; prefer diags+tracks first  

## Relation to v80 thesis

This is the **product line** (force/orbit graph on multi-fab). It does **not** implement free-substrate representation. Morning score feeds whether product multi-fab still has legs vs soft-kill + representation toy.
