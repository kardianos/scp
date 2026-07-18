# Phase 2 Score — FINAL

**Agent:** TU · **Date:** 2026-07-19 · **Round:** P2-R5  
**Status:** **FINAL** — Phase-2 terminal recommendation  
**Map:** `v77/CAMPAIGN_MAP.md` · **Board:** `CHECKPOINT_BOARD.md`  
**History (not reopened):** Phase-1 CONVERGENCE STOP (A) composition-tier  

---

## 0. Verdict (one line)

**STOP Phase 2 (A-depth): M1 ∧ RC1 green; R-compose closed at RC1 co-field; not DISPROVE.**

---

## 1. Terminal bar (CAMPAIGN_MAP §0)

| Goal | Done when | Result |
|------|-----------|--------|
| **P2-M** | Comprehensive **dynamic** full Maxwell (M1 mandatory) | **MET** — `m1_claim=true` |
| **P2-RC** | Close **R-compose** (RC1 mandatory) | **MET** — `rc1_claim=true`; one grid ψ + dynamical Maxwell |
| **P2-U** | Cohesive agreement on joint model | **MET** — CP-U-FINAL ADOPT |

**STOP Phase 2 when M1 ∧ RC1 green with co-agreement** → **YES**.

---

## 2. Checkpoint scorecard (mandatory path)

| Checkpoint | Status | Evidence |
|------------|--------|----------|
| CP0 | **ADOPTED** | O-009 |
| CP-M1-SPEC | **ADOPTED** | TE+NE dual |
| CP-M1-NUM | **ADOPTED** | `work/NE/outputs/m1_result.json` m1_claim=true; O-011 |
| CP-RC1-SPEC | **ADOPTED** | TE, NE, TM, TD(ψ), NM |
| CP-RC1-NUM | **ADOPTED** | `work/NM/outputs/rc1_result.json` rc1_claim=true; RG0–6 |
| **CP-U-FINAL** | **ADOPTED** | TU-036 + this document |

Optional CP-M2 / CP-RC2 / CP-M3: **not required** for Phase-2 STOP.

---

## 3. M1 evidence (dynamic full Maxwell 2D)

| Gate | PASS | Note |
|------|------|------|
| M1-R0 | ✓ | M0 regression |
| G1 vacuum | ✓ | |
| G2 true-2D beam | ✓ | not 1D TEM; v/c≈1.044 |
| G3 off-unit | ✓ | \(c_{\mathrm{th}}=0.5\) |
| G4 energy | ✓ | drift ~0.5% |
| G5 Gauss dyn | ✓ | Cont-safe; mild cleanse documented |
| G6 div B | ✓ | ~1e-15 |
| G7 Faraday | ✓ | |
| G8 Ampère adversary | ✓ | adv kills free wave |
| G9 BC honesty | ✓ | multi-BC documented |

**Claims:** `m1_claim=true`; `full_maxwell_2d_dynamic=true`.  
**Sandbox:** `work/NE/sandbox_m1_2d.py`.

---

## 4. RC1 evidence (co-field — R-compose closed)

### 4.1 What RC1 proves

On **one continuum grid** (2D):

1. Fixed multi-locks carry \((\rho_b,\rho_Q)\).  
2. \(\rho_b \to \psi\) via free-capacity F1 (2D Laplace).  
3. \(\rho_Q \to (\mathbf{E},\mathbf{B})\) via **dynamical** NE `Maxwell2D` Yee (**not** Φ-only claim path).  
4. Shared \(c\); budget deficit; sibling independence (TE-IA1).  
5. Force diagnostics \(F^\psi\) + \(Q\mathbf{E}\) on fixed locks.  
6. Configs: vacuum / neutral / same-sign / opposite all gate-clean.

### 4.2 Gate table

| Gate | PASS |
|------|------|
| RG0 joint state / not lite-only | ✓ |
| RG1 dual source / sibling indep / ≥2 locks | ✓ |
| RG2 shared \(c\) | ✓ |
| RG3 force taxonomy N/S/O/V + dynamical qE | ✓ |
| RG4 TE-IA1 (neutral no E; no ψ≡E id) | ✓ |
| RG5 dynamical Maxwell nstep / B / M1 inherit | ✓ |
| RG6 vacuum / energy split / hierarchy soft | ✓ |

**rc1_claim = true**  
**Maxwell API:** `import_NE_sandbox_m1_2d.Maxwell2D` · `used_NE_Maxwell2D=true` · `m1_claim_inherited=true`  
**Sandbox:** `work/NM/sandbox_rc1_cofield.py`

### 4.3 Key numbers (parent run)

| Config | \(\psi_{\mathrm{mid}}\) | \(\max\|E\|\) | Note |
|--------|-------------------------|---------------|------|
| vacuum | 0 | 0 | clean |
| neutral | ≈0.937 | 0 | mass → ψ only |
| same_sign | ≈0.937 | ≈0.230 | both channels on |
| opposite | ≈0.937 | ≈0.278 | charge structure |

Shared \(c\): \(C_{\mathrm{LOCAL}}=1=1/\sqrt{\varepsilon\mu}\). Free deficit ≈0.169 on mass configs.

---

## 5. R-compose status

| Tier | Status |
|------|--------|
| RC0 composition (Phase-1) | LIVE history — separate sandboxes |
| **RC1 co-field** | **CLOSED** — one grid, dynamical Maxwell + ψ, fixed locks |
| RC2 co-evolution (moving locks) | **Not required** — optional stretch |
| RC3 3D co-evolution | stretch after M2 |

**Phase-1 residual “not one binary”** is **closed at RC1** for the mandatory depth bar.  
Remaining honesty: 2D (not 3D) co-field; locks fixed (not RC2 motion).

---

## 6. Residuals (honest, non-blocking for A-depth STOP)

| ID | Residual | Blocks Phase-2 STOP? |
|----|----------|----------------------|
| **R-2D-ψ** | Co-field ψ is **2D Laplace** (log exterior class); 3D \(1/r\) path-cost remains seed/v76 + diagnostic Green ref | **No** — labeled composition |
| **R-fixed** | Locks **fixed** (RC1); no co-evolution under \(F^\psi+F^{\mathrm{EM}}\) | **No** — RC2 optional |
| **R-force-sign** | Force diagnostics / soft hierarchy; full dynamics-derived \(\alpha_\psi\) open | **No** — virtual-work class documented |
| **R-dim-Max** | Maxwell dynamics 2D (M1); not 3D Yee multipoles | **No** — M2 optional |
| **R-iso** | Poisson-form isomorphism residual | Soft ontology (Occam tags) |
| **R-J5α** | Raw inertia triad deferred | V77-3 β already closed |
| **R-scope** | No real \(G,\alpha\); no scp_sim monism proof; no full GR | Out of Phase-2 bar |

**Kills K-M1 / K-RC1:** **not fired**.

---

## 7. G1–G3 / program relation

| Goal | Phase-2 reading |
|------|-----------------|
| G1 widely applicable idea | **Held** — free/bound + \(c\) + sibling channels organize Maxwell + locks |
| G2 multi-focus unity | **Strengthened** — joint co-field, not only composition |
| G3 build or kill | **BUILD → STOP (A-depth)** |

Phase-1 composition STOP (A) remains valid history. Phase-2 is **depth exam pass**, not a reopen of v76.

---

## 8. Recommendation to O

```text
RECOMMEND: STOP Phase 2 (A-depth)
  M1 ∧ RC1 green; R-compose CLOSED at RC1; CP-U-FINAL ADOPT
NOT: DISPROVE / NULL / V77-K
OPTIONAL LATER (not required): M2 3D Maxwell; RC2 moving locks; M3 charged matter
```

### Reproduce

```bash
cd /home/d/code/scp/v77/work/NE && python3 sandbox_m1_2d.py --quick
cd /home/d/code/scp/v77/work/NM && python3 sandbox_rc1_cofield.py
```

### Read order for human / O

1. This file (`PHASE2_SCORE.md`)  
2. `CHECKPOINT_BOARD.md`  
3. `work/NE/outputs/m1_result.json`  
4. `work/NM/outputs/rc1_result.json`  
5. `work/TM/rc1_joint_state_v0.md` · `work/TE/full_maxwell_monist_v0.md`  

---

## 9. Document control

| Version | Status |
|---------|--------|
| PHASE2_SCORE draft (P2-R3) | superseded |
| **PHASE2_SCORE FINAL (P2-R5)** | **this file** — terminal recommendation |
