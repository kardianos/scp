# v75 SELF_TUNE_C — Option-C-ready self-tuning with predictive backtracking

**Status:** Stage 1 implemented (B1 soft θ, frozen action per trial).  
**Date:** 2026-07-14  
**Does not modify** `scp_sim.c` / `scp_sim.cu` / `sfa.h`.

---

## 1. Goal

Find whether **outer-loop soft parameters** on the current B1 multi-fabric package
can produce a **persistent nucleus + light cloud** (atom packaging), or whether
the failure surface is **hard** (same wall for all soft θ → need Option C / B2 /
lower-g nuclear redesign).

Success is **definitive**, not the first lucky partial:

| Verdict | Meaning |
|---------|---------|
| **DEFINITIVE SUCCESS** | ≥1 soft-θ trial with T≥T_full, massC and massL loss ≤15%, both masses >0.5× seed, Gauss floor held, no L-into-C bag (massL stays in L track) |
| **DEFINITIVE FAIL (soft θ)** | After the planned tree is exhausted, **every** branch ends FAIL_Q / FAIL_L / FAIL_E with the **same wall class**; best cost still above PASS threshold |
| **INCONCLUSIVE** | Budget cut short (GPU death, etc.) — not used as a physics claim |

---

## 2. Architecture layers

```
┌─────────────────────────────────────────────────────────────┐
│  Layer 3 — Auto-tuning ledger (JSONL + summary TSV)         │
│  trial_id, parent_id, θ, cost, pred, outcome, action        │
└──────────────────────────▲──────────────────────────────────┘
                           │
┌──────────────────────────┴──────────────────────────────────┐
│  Layer 2 — Predictive backtracking (parameter tree)         │
│  Early abort from diag predictors; restore last good parent │
│  Propose θ' from branch rules — NEVER free mid-Verlet ε     │
└──────────────────────────▲──────────────────────────────────┘
                           │
┌──────────────────────────┴──────────────────────────────────┐
│  Layer 1 — Self-tuning soft weights θ (outer loop only)     │
│  θ = {n_C, ω_C, n_L, ω_L, R_shell, D_nuc, g_gauge, …}       │
│  Frozen for each full trial; optional adiabatic ramp only   │
│  if predeclared in cfg (not used in Stage 1)                │
└──────────────────────────▲──────────────────────────────────┘
                           │
┌──────────────────────────┴──────────────────────────────────┐
│  Layer 0 — Frozen L0 action (B1 Shape β)                    │
│  n_fabrics=3, mf_lock_CQ=1, private bags, shared A          │
│  Hard zeros: L ∉ nuclear bag, K_13-class portal off         │
│  q_C=0, q_Q=+1, q_L=−1, U^q transport, same-sign ω          │
└─────────────────────────────────────────────────────────────┘
```

**Hard rule:** Soft θ updates happen **between trials**, not inside Verlet.
Mid-run free \(K_{st}(t)\) or free \(\varepsilon(t)\) is **forbidden** — that
would invalidate energy/ledger as theory tests.

---

## 3. Stage 1 vs Stage 2

### Stage 1 — B1 soft θ (this implementation)

No Option C kernel. Knobs that already exist via seeds + config:

| Symbol | Meaning | Stage-1 domain |
|--------|---------|----------------|
| `n_C` | Nuclear ball count | 1, 2, 4, 6 |
| `omega_C` | Nuclear frequency / profile | 1.42, 1.46 |
| `n_L` | Light ball count | 1, 2, 4, 6 |
| `omega_L` | Light frequency / profile | 1.46, 1.485 |
| `R_shell` | L shell radius | 16, 18, 22, 24 |
| `D_nuc` | Nuclear nn edge (octa/tetra) | 10 (fixed Stage 1) |
| `g_gauge` | Coupling | **0.05 fixed** (profiles are g-matched) |

Controller: `v75/analysis/self_tune_controller.py`  
Ledger: `/space/scp/v75/self_tune/ledger.jsonl`

### Stage 2 — Option C weights (only if Stage 1 is definitive FAIL *and*
the ledger attributes the wall to **engagement graph**, not scale)

Requires kernel: three Cosserat blocks + sparse \(K_{st}\), hard \(K_{13}=0\).
Then extend θ with \(K_{12},K_{23}\). **Not implemented until Stage 1 closes.**

---

## 4. Cost and predictors

### 4.1 End-of-trial cost (from `mf_pair_track` + diag)

\[
\begin{aligned}
c_Q &= \mathrm{clip}\big((Q_{C0}-Q_{C})/Q_{C0},\,0,1\big) \\
c_L &= \mathrm{clip}\big((m_{L0}-m_L)/m_{L0},\,0,1\big) \\
c_E &= \mathrm{clip}\big(|E-E_0|/E_0 / 0.20,\,0,1\big) \\
c_G &= 1 \text{ if } \max|g_{\mathrm{gauss}}| > 10^{-10} \text{ else } 0 \\
c &= 0.35\,c_Q + 0.35\,c_L + 0.20\,c_E + 0.10\,c_G
\end{aligned}
\]

**PASS** if \(c \le 0.15\) and \(c_Q\le 0.15\) and \(c_L\le 0.15\) and \(c_G=0\)
at \(T \ge T_{\mathrm{full}}\).

### 4.2 Predictive abort (during / after screen)

On diag rows (cheap, no SFA required mid-screen):

| Signal | Threshold | Code |
|--------|-----------|------|
| \(Q_\phi(t)/Q_\phi(0) < 0.80\) before \(T_{\mathrm{screen}}\) | nuclear evaporating | `PRED_Q` |
| \(E(t)/E(0) < 0.88\) | energy dump | `PRED_E` |
| \(s_{\max} > 1.5\) | fusion/rearrange spike | `PRED_S` |
| `gauss_max > 1e-10` | implementation tripwire | `PRED_GAUSS` |

On end track (after short SFA):

| Signal | Threshold | Code |
|--------|-----------|------|
| massL loss > 25% | light shed | `PRED_L` |
| massL / massL0 < 0.05 | L gone | `FAIL_L_ABSORB` |

On predict-fail: **abort trial**, log, **backtrack** to parent, try next sibling.
Do **not** mutate θ mid-trajectory.

### 4.3 Screen → full

1. Every candidate runs **screen** \(T=T_{\mathrm{screen}}\) (default 150).  
2. If predictors fire or cost_screen > 0.35 → FAIL branch, no full.  
3. If screen looks viable → **full** \(T=T_{\mathrm{full}}\) (default 400).  
4. Full outcome decides PASS / FAIL for that leaf.

---

## 5. Parameter tree (predictive backtracking)

Root starts from F14 failure class and walks **fixes**, not random walk:

```
ROOT  (virtual: F14 Z6@1.42 + L6@1.46 R~16 FAIL known)
 │
 ├─ B1  n_C=6 ω_C=1.46  n_L=6 ω_L=1.46  R=18   # parked nuclear hypothesis
 │   ├─ B1a  R=22
 │   └─ B1b  ω_L=1.485 R=22
 ├─ B2  n_C=4 ω_C=1.46  n_L=4 ω_L=1.46  R=18   # smaller nucleus
 │   └─ B2a  R=22
 ├─ B3  n_C=2 ω_C=1.46  n_L=2 ω_L=1.46  R=18
 │   └─ B3a  R=22
 ├─ B4  n_C=1 ω_C=1.42  n_L=6 ω_L=1.46  R=18   # single heavy + multi-L
 │   └─ B4a  R=22 ω_L=1.485
 └─ B5  n_C=6 ω_C=1.46  n_L=2 ω_L=1.46  R=20   # heavy multi, few light
```

**Backtrack rule:** if a node fails with `PRED_Q` / nuclear class, prefer siblings
that **reduce n_C or raise ω_C** (lower seed Q). If fails with `PRED_L` only,
prefer siblings that **raise R or soften ω_L**. If fails both → escalate to
smaller n_C.

**Exhaustion:** after all listed leaves evaluated (or pruned by parent FAIL class
with no remaining relevant sibling), declare **DEFINITIVE FAIL (soft θ)** if no
PASS; else **DEFINITIVE SUCCESS** with best θ.

---

## 6. Ledger schema (JSONL one object per line)

```json
{
  "trial_id": "B1a_screen",
  "parent_id": "B1_screen",
  "phase": "screen|full",
  "t_sim": 150.0,
  "theta": {
    "n_C": 6, "omega_C": 1.46, "n_L": 6, "omega_L": 1.46,
    "R_shell": 22.0, "D_nuc": 10.0, "g_gauge": 0.05
  },
  "metrics": {
    "Qc0": 0, "Qc": 0, "massC0": 0, "massC": 0,
    "massL0": 0, "massL": 0, "E0": 0, "E": 0,
    "gauss_max": 0, "s_max_peak": 0, "cost": 0
  },
  "pred": "PRED_Q|PRED_L|none",
  "outcome": "PASS|FAIL_Q|FAIL_L|FAIL_E|FAIL_GAUSS|FAIL_PRED|PARTIAL",
  "action": "accept|reject|backtrack|promote_full|stop_success|stop_fail",
  "cfg_hash": "...",
  "seed_note": "octa C + shell L",
  "wall_s": 0
}
```

Summary TSV: `ledger_summary.tsv` for quick human scan.

---

## 7. Implementation map

| Path | Role |
|------|------|
| `v75/SELF_TUNE_C.md` | this design |
| `v75/analysis/self_tune_controller.py` | outer loop, tree, ledger, seed build |
| `/space/scp/v75/self_tune/` | work dir: seeds, cfg, diag, track, ledger |
| GPU binary | existing `scp_sim_mf_cuda` / rebuilt `scp_sim.cu` |
| Track | `bin/mf_pair_track` |

**Kernel policy:** Stage 1 = **no** kernel edits. Stage 2 only with explicit
user authorization for Option C.

---

## 8. Relation to Option C

Option C (triple Cosserat + \(K_{st}\)) is the **implementation substrate** for
Stage 2 weights. Stage 1 answers first:

> Is atom packaging blocked by **soft scale choices** on B1, or by the
> **engagement structure** itself?

If Stage 1 finds a PASS θ → document and freeze; Option C not required for the
atom claim at that scale.  
If Stage 1 definitive FAIL with wall = nuclear supercriticality only → next is
lower-\(g\) / sub-\(Q_{\max}\) redesign, not necessarily C.  
If Stage 1 FAIL with wall = “L always erodes when any multi-C rearranges” even
for parked c6_light → engagement / BO / Option C becomes the next design move.

---

## 9. How to run

```bash
# Local smoke (small N, short T) — controller self-contained
python3 v75/analysis/self_tune_controller.py --mode local --N 64 --L 20 \
  --T_screen 20 --T_full 40 --work /tmp/scp_self_tune_smoke

# Production GPU (on instance or with --sim pointing at CUDA binary)
python3 v75/analysis/self_tune_controller.py --mode gpu \
  --sim ./scp_sim_mf_cuda --N 192 --L 48 \
  --T_screen 150 --T_full 400 \
  --work /space/scp/v75/self_tune
```

Campaign stops on definitive PASS or when the tree is exhausted (definitive FAIL).
