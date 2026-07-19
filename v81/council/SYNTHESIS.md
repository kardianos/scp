# v81 design council synthesis

**Date:** 2026-07-19  
**Advisors:** claude fable · crush (zai/glm-5.2) · kimi k3  
**Full reports:** `claude_fable/`, `crush_glm52/`, `kimi_k3/`  
**Brief:** `../DESIGN_BRIEF.md`

---

## Consensus (unanimous)

| Claim | Support |
|-------|---------|
| Multi-fab product path is **dead for Stage 3** | All three; O FAIL + v79 \(E_{\mathrm{em}}\) null |
| Force medium (shared U(1)/Maxwell) **works** | All three — keep medium, change light **carrier** |
| Light charge as **field hump / Q-ball** cannot be durable | All three — soft evaporate / hard shred are same defect |
| **Primary path: free medium + typed locks** (v80 shapes 4+10, 6 literal) | All three (different names: lock-carrier / PIC-monist / link-gauge locks) |
| **Sandbox first** (standalone C, ~2 weeks), **not** kernel rewrite day 1 | All three — auth ≠ obligation |
| First metric: **durability + ledger honesty**, then binding | All three |
| **No pairwise Coulomb** — force only through live medium | All three (kill if bypassed) |
| **No multi-fab re-runs**, no Z6+L6, no soft \(v_t\) grids | All three |
| Positronium first, not carbon | All three |

---

## Three redesign options (aligned)

| Rank | Path | What it is | Kernel? | Risk |
|------|------|------------|---------|------|
| **1** | **Locks on gauge medium** | Discrete locks \(q=\pm1\), \(E^\star\), continuous \(x,p\); Boris push; Esirkepov deposit; shared Yee/U(1) free medium; optional capacity/sequestration | Port later as `scp_locks` after sandbox green | Looks like PIC — discipline is: no hump-as-definition, no pair potential, C/M/E charts real |
| **2** | **Path-measure \(\ell\) + locks** | Option 1 + hyperbolic \(\ell\) thickens around locks (shape 3) | After Option 1 green | Double unproven change; elliptic collapse / GRIN kill |
| **3** | **Token/CA pure substrate** | Pure update-budget currency; vortices = charge | Never in kernel first | Research probe only; week-scale positronium unlikely |

Claude’s “integer-flux endpoints” and crush’s “defect-primary gauge” sit as **later upgrades** of Option 1 charge honesty, not as first steps.

---

## Recommended 2-week plan (merged)

```text
Week 1 — monist_pic / lock_medium sandbox (CPU)
  D1–3  Yee (or Maxwell) + lock list + charge-conserving deposit + Boris
        Ledger close: field + lock + shell flux <1–2%
  D4–7  Rest-pair F(D) scan D=8..24 (reproduce F PASS without multiplet L)
        Soft long-T: locks cannot evaporate by construction; E_em floors

Week 2 — durability + decision
  D8–10 Hard vt / multi-rev: ≥5 revs if possible; RR inspiral rate
  D11–12 C/M/E chart readouts as first-class output
  D13–14 FINDINGS + kill-gate scorecard
         GREEN → port scp_locks into scp_sim (authorized)
         RED   → write which gate; pivot Option 2 only if Option 1 collapses to humps
```

**Day-14 product:** two-lock positronium analog with intact locks and honest energy ledger — or a named kill.

---

## Mapping to v80 thesis

| Shape | In Option 1 |
|-------|-------------|
| 4 free+locks | Architecture |
| 10 gauge/defect | Integer \(q\) + Gauss-protected deposit |
| 6 update budget | CFL as engine law; lock \(v\) from ledger |
| 3 path measure | Slot later (Option 2); optional capacity field in Claude A |
| C/M/E charts | Readouts of one state, day-one |

---

## What not to do (merged)

- Multi-fab Z6+L6 / orbit grids / fourth fabric  
- Electron = small Q-ball or radial profile seed  
- Pairwise Coulomb / springs / damping  
- Kernel rewrite before sandbox green  
- Variable \(c\) as free-frame gravity (v76 GRIN kill)  
- Carbon before positronium  
- Charts without a real state change  

---

## Decision for implementer

**Yes — step away from multi-fab product.**  
**Yes — three independent Stage-3 redesigns exist; all three advisors converge on the same primary.**  
**Next code step (when you green-light):** scaffold `v81/sandbox/` lock+medium CPU sim per Week 1 plan — **not** immediate Cosserat multi-fab edits.
