# v80 — Representation gap: free substrate + locks + charts

**Date:** 2026-07-19  
**Status:** design / brainstorm — **no kernel work**, no multi-fab reruns  
**Inherits:** v76 monist free-capacity; v77 Maxwell + dual-channel RC1; v78 R1–R10 freeze; v79 multi-fab long-T atom null

## One-sentence goal

Name and develop a **simulation representation** adequate to v78’s free/bound continuum story — something that is **not** “another multiplet of continuum fields on a fixed Eulerian grid.”

## Why this version exists

v78 freezes **numeric relations** (locality \(c\), \(m\leftrightarrow E/c^2\), free/bound budget, dual sources, Coulomb, …).  
v79 runs the **multi-fabric atom recipe** (C/Q/L on fixed N³) and finds: Gauss OK, core OK, **atom product fails** (L nulls \(E_{\mathrm{em}}\)).

The working hypothesis for v80 is that multi-fab (and fixed-grid Cosserat generally) is an **emulation** of what the theory is getting at. The theory has an expression; it still lacks a **primary state type** for a simulator. Relations→field also fails. Multi-fab is getting at dual channels / isolation, but as **N copies of the same representation**.

## Documents

| File | Role |
|------|------|
| `CONTEXT.md` | Lineage v74–v79; what is frozen; what failed; two stacks |
| `REPRESENTATION.md` | Core thesis: real substrate + locks + POV charts (C/M/E) |
| `SHAPES.md` | Candidate substrates (focus 3, 4, 6, 10; full list) |
| `DUAL_CHARTS.md` | Reality as real + point of view; M-const / E-const / c-const |
| `FEEDBACK_REQUEST.md` | Brief for external review |
| `KIMI_FEEDBACK.md` | External agent feedback (kimi CLI, **model k3** / `kimi-code/k3`) |
| `KIMI_FEEDBACK.raw.txt` | Raw CLI transcript (includes thinking) |

**External review:** `kimi -m kimi-code/k3 -p …` (2026-07-19). Session resume: `kimi -r session_a833d053-3c4d-44ab-a756-6da1cf7bbf78`.

## Non-goals (this version)

- Do **not** modify `scp_sim` / `sfa.h` without explicit later authorization  
- Do **not** re-run multi-fab atom matrix as the primary probe of monism  
- Do **not** “derive a field from R1–R10” as the main path  
- Do **not** require pure novelty unused in all of physics — filter is **primary representation**, not fashion

## Success for v80 (design bar)

1. Written atlas: real state / step / lock / charts / multi-fab-as-projection  
2. External critique recorded  
3. Clear kill gates for candidates that collapse back to φ(x) on fixed N³  
4. Optional next: one thin sandbox of a **new state type** (not a carbon atom run)

## Next-step council (four independent advisors)

**See `council/`** — same brief to internal subagent, kimi k3, claude fable, agy Gemini 3.1 Pro.  
Consensus primary: hydrogenoid C–L **force/orbit with kept SFAs + pair tracks**; no Z6 re-park; representation toy parallel CPU.

## Product overnight campaign (parallel)

**Done:** `campaign/` — multi-fab force/orbit graph on V100 (not representation toys). **G=0.62.**

| Item | Value |
|------|--------|
| Plan | `campaign/CAMPAIGN.md` |
| Goal function | `campaign/GOAL.md` |
| Steps | `S0_smoke` … `S4_hydrogenoid` |
| Remote | `/root/v80/work` queue PID; binary `/root/v80/bin/scp_sim_cuda` |
| Poll | every 15 min → `campaign/scripts/poll_and_sync.sh` + scheduler |
| Seeds | `/space/scp/v80/seeds/` |
| Results | `campaign/results/` + per-step `*/results/` |

Teardown Vast **45259142** when `campaign/logs/COMPLETE` and downloads verified.
