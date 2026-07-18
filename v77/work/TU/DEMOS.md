# DEMOS — Live / Dead / Planned Registry

**Owner:** TU · **Round:** P2-R5 FINAL · **Source of truth for demo IDs**

**Status vocabulary:** LIVE / LIVE_PASS / LIVE_PARTIAL / DRAFT / PLANNED / DEAD

---

## 1. Registry

### Seed (v76)

| Demo ID | Status | Evidence |
|---------|--------|----------|
| **D-ψ-3D** | **LIVE** | v76 goal2_PC3D_workable |
| D-ψ-local-GRIN | **DEAD** | long-range fail |
| D-ψ-hand-Poisson | **DEAD** | not monist proof |
| D-ψ-2D-log | **DEAD as GR exterior** | log multipole |

### EM — M0 / M1

| Demo ID | Status | Evidence |
|---------|--------|----------|
| **D-EM-full-maxwell** | **LIVE_PASS** | NE R2 M0 FM1–FM7 |
| **D-EM-M1-suite** | **LIVE_PASS** | NE `m1_result.json` **m1_claim=true**; CP-M1-NUM |
| D-EM-M1-beam2d | **LIVE_PASS** | true_2d_packet |
| D-EM-M1-energy | **LIVE_PASS** | G4 |
| D-EM-M1-gauss-dyn | **LIVE_PASS** | G5 |
| D-EM-M1-ampere-adv | **LIVE_PASS** | G8 |
| D-EM-maxwell-lite / coulomb / wave / divB / faraday / continuity | **LIVE_PASS** | R1–R2 subsets |
| D-EM-sibling-psi | **LIVE_PASS** | Phase-1 dual lite; superseded depth by RC1 |

### RC1 co-field (Phase-2 depth)

| Demo ID | Status | Evidence |
|---------|--------|----------|
| **D-RC1-cofield** | **LIVE_PASS** | NM `rc1_result.json` **rc1_claim=true**; RG0–6; NE Maxwell2D + 2D F1 ψ same grid; CP-RC1-NUM |
| D-JOINT-yee-psi | **LIVE_PASS** (alias) | same as D-RC1-cofield at RC1 bar |

### Dual-channel / Matter (Phase-1 + RC1)

| Demo ID | Status | Evidence |
|---------|--------|----------|
| **D-DUAL-channel** | **LIVE_PASS** | NM r2 Maxwell-lite joint (regression) |
| D-MAT-lock-S0 / dual0 / force-tax / hier | **LIVE_PASS** | R1 G0–G3 |
| D-MAT-force-maxwell | **LIVE_PASS** | RC1 dynamical \(Q\mathbf{E}\) on co-field |
| D-SCP-vocab-bridge | **DRAFT** theory-only | no kernel monism claim |

### Dynamics

| Demo ID | Status | Evidence |
|---------|--------|----------|
| **D-DYN-j5-formfactor** | **LIVE_PARTIAL** | J5-β; V77-3 MET |
| D-DYN-free-wave-c | **LIVE_PASS** | \(v=c\) |
| **D-DYN-dual-channel-c** | **LIVE_PASS** | shared \(c\) waves |
| D-DYN-timedep / rods | DRAFT / PLANNED | soft residual |

### Unification

| Demo ID | Status | Evidence |
|---------|--------|----------|
| **D-UNIFIED-package** | **LIVE** | Phase-1 composition |
| **D-PHASE2-M1-RC1** | **LIVE** | Phase-2 depth: M1 ∧ RC1; `PHASE2_SCORE.md` FINAL |

---

## 2. Tier board

| Tier | Status |
|------|--------|
| Phase-1 V77-1…V77-3 / composition (A) | MET (history) |
| Phase-2 M1 | **MET** |
| Phase-2 RC1 | **MET** — R-compose closed |
| Phase-2 STOP (A-depth) | **RECOMMENDED** |
| DISPROVE / V77-K | **NO** |
| RC2 / M2 | optional stretch |

---

## 3. Honesty

```text
RC1 LIVE:  one 2D grid · fixed locks · F1 ψ + dynamical Yee E,B (NE Maxwell2D)
Not claimed: 3D co-field multipoles; moving locks (RC2); real G/α; scp_sim monism
```

---

## 4. Changelog

| Stamp | Change |
|-------|--------|
| Phase-1 | composition LIVE |
| P2-R3 | M1 LIVE; RC1 planned |
| **P2-R5** | **D-RC1-cofield LIVE_PASS**; Phase-2 STOP recommend |
