# v75 Status

**Updated**: 2026-07-16 — **State written, F19 scored; push checkpoint**

## Canonical snapshot

**Read first:** [`STATE.md`](STATE.md) — equations, charge table, P/N freeze, F17–F19 numbers, open problems.

## Active goal

**`C12_ATOM_GOAL.md`** — multi-fabric Phases **1–3** → time-stable C₁₂ analog.

| End state | Definition |
|-----------|------------|
| **Ideal** | Time-stable **C₁₂** core + L cloud; firm Z/N; no merge/fission/disperse |
| **Stretch** | Isotope **+N at fixed Z** + calculated/simulated decay |

## Closed baseline

| Item | Result |
|------|--------|
| F11–F16 | B1 isolation, Coulomb, B4 packaging, pair kinematics |
| **F17** | **P/N firm:** B2; n \(Q_\mathrm{flux}=0\); p charged |
| **F18** | Z2 park + L=Z + isotope \(Q_\mathrm{flux}\) identical **PASS** |
| **F19** | Z6: isotope EM **PASS**; Z6N6 park PASS; Z6N0 soft; L −12.5%; **nuclear multi-ball → single droplet** by t∼400 |

**Data:** `/space/scp/v75/pn/{,p2/,z6/}` · scores `v75/results/pn_z6/` · sheet `images/C6_atom_sheet.png`  
**Instance:** `v75f16` idle (campaign complete).

## Phase checklist

| Phase | Focus | Status |
|-------|--------|--------|
| **P1** | Multi-rev hydrogenoid, shell-radius, binding | OPEN |
| **P2.0** | Firm p vs n | **DONE** |
| **P2.1–2.5** | Park + L-from-Z + isotope (Z2) | **DONE** (F18) |
| **P2 Z6** | Carbon-class scale | **F19 scored — partial park/L** |
| **P3** | A≈12 long-T visual C₁₂ | NOT STARTED |

## Freeze (do not re-litigate without new data)

```
n_fabrics=3  mf_lock_CQ=0  q_C=0 q_Q=1 q_L=-1
p = C bag + Q charge   n = C only   L count = Z (not A)
gen_pn_core … [profL omegaL]
Grid Z6: N=192 L=48 damp_width=5  L shell R=22
```

## Next (when work resumes)

1. Retune Z6N0 park (\(c_Q\le 0.15\))  
2. Improve L hold at carbon-class shell  
3. Shell-radius diagnostic (COM D≈0 insufficient)  
4. Long-T visual C₁₂ package  
5. Larger GPU only if N≳256 needed (`PN_GRID.md`)

## Ops

- **No `sleep` inside `sim_exec`** — see `analysis/remote_status_snippet.sh`  
- Park-aware scoring for multi-ball; isotope metric = \(Q_\mathrm{flux}\) (nuclear) + L track
