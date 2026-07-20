# v78 Particle Ladder — Field → Nucleon → P/N → C nucleus → Atom

**Campaign:** v78  
**Status:** frozen integration (U)  
**Date:** 2026-07-19  
**Peer packages:** W `PHYSICS_RELATIONS.md` + `world_freeze_v0`; N `nucleon_package_v0`; P `pn_package_v0`; C `carbon_nucleus_v0`; L `light_sector_v0`; A `c12_atom_package_v0`

Structural ladder only — parallels to real particles are **not** quantitative (CONCEPT.md).

---

## 0. World foundation (rung 0)

| Item | Content |
|------|---------|
| **Matter** | Complex Cosserat \(\Phi_a,\Theta_a\) (× fabrics C/Q/L) |
| **EM** | Diagonal U(1) \(A,E\); \(g=0.05\) standard |
| **Monist parallel** | Free capacity \(\psi\), free Maxwell (v76/v77) — language only for R1,R2,R4,R9 |
| **Relations** | R1–R10 — root `PHYSICS_RELATIONS.md` (W) · U freeze `work/U/RELATION_MATRIX.md` |
| **Claim** | Particles = stable localized field configurations |

**Primitive stack (W freeze L0–L6):** continuum \(c\) → Cosserat → private bags \(V_t\) → shared gauge + \(q_f\) → ledgers → optional \(\psi\) → VK/park gates.

**Gate R0:** Relation matrix assigned · **PASS** (CP-W: W ADOPT + U ADOPT)

---

## 1. Nucleon (rung N)

| | |
|--|--|
| **Object** | Gauged Q-ball; three equal \(\Phi_a\); \(\Theta=0\); \(\eta=0\); \(m_\theta=1.6\); \(g=0.05\) |
| **Default light** | \(\omega=1.46\), \(Q_N=114.13\), \(E/Q=1.5187\), \(r_{1/2}=2.63\) |
| **Profile** | `v74/profiles/f_w146_g005.txt` |
| **Stable** | VK \(\mathrm{d}Q/\mathrm{d}\omega<0\) on \(\omega\in[1.406,1.4825]\); \(Q_{\max}=921\) |
| **Evidence** | gscan; CONCEPT charge retention; v71 he2/li3s; v74 free unit |
| **Package** | `work/N/nucleon_package_v0.md` · `work/N/recipes.sh` |

**Gate N:** Spec + profile + seed + stability · **PASS** (CP-N: N ADOPT + U ADOPT)

---

## 2. Proton- and neutron-analogs (rung P / N⁰)

**Package:** `work/P/pn_package_v0.md`

B2 unlock (`mf_lock_CQ=0`), charges \(q_C=0,\,q_Q=+1,\,q_L=-1\):

| Species | Content | EM |
|---------|---------|-----|
| **p** | C bag + Q co-located | \(Q_{\mathrm{flux}}\neq 0\) |
| **n** | C-only | \(Q_{\mathrm{flux}}=0\) |
| **L** | fabric L | \(q_L=-1\) |

\[
Z=\#\{\mathrm{p}\},\quad N=\#\{\mathrm{n}\},\quad A=Z+N;\quad n_L=Z.
\]

**Not P/N:** flavored \(\Delta\omega\) alone; cancel multiplet; B1 lock.

**Evidence:** F17 S1–S5 PASS; F18/F19 isotope \(Q_{\mathrm{flux}}\) identical at Z=2 and Z=6.

**Gate P:** Species + isotope EM · **PASS** (CP-P: P ADOPT + U ADOPT)

---

## 3. Carbon nucleus (rung C)

**Package:** `work/C/carbon_nucleus_v0.md`

### Single-fabric liquid-drop (v74)

| Run | Final Q | Fate |
|-----|---------|------|
| **c6_light** (Z-carbon) | **650** (0.71×Q_max) | **PASS** park; on-branch +0.51%; mass defect 4.85% |
| **c12_light** | **1411** (1.53×Q_max) | Super-critical evaporating droplet |

**Theorem (g=0.05):** \(12\times Q_N^{\min}>Q_{\max}\) ⇒ free A=12 never parks.

### Multi-fabric cores (v75 F19)

| Setup | Park | Note |
|-------|------|------|
| Z6N6 | \(c_Q=0.046\) PASS | Multi-ball → **droplet** by t∼400 |
| Z6N0 | \(c_Q=0.184\) soft | Same \(Q_{\mathrm{flux}}=990\) |

**Gate C:** Z-carbon PASS + A=12 characterised · **PASS** (CP-C: C ADOPT + U ADOPT)

---

## 4. Light opposite sector (rung L)

**Package:** `work/L/light_sector_v0.md`

| Path | Status |
|------|--------|
| E-lite ± force | F8 **PASS** |
| E-lite orbit | F10 **PARTIAL** (arcs; head-on merge) |
| MF isolation / force | F11–F13 **PASS** |
| Soft orbit + multi-L shell | A1/B4/F16 **PASS** |
| Z-matched cloud Z=2 | F18 **PASS** (−0.6% massL) |
| Z-matched cloud Z=6 | F19 **PARTIAL** (−12.5% massL) |
| Multi-rev immortality | **OPEN** |
| True \(m_e\)-scale single ball | **BLOCKED** at g=0.05 (Q_min~90) |

**Gate L:** Prep + hydrogenoid path · **PASS** (CP-L: L ADOPT + U ADOPT)

---

## 5. C₁₂-scale atom (rung A)

**Package:** `work/A/c12_atom_package_v0.md` · checklist `recipe_checklist_multifab.md`

| | |
|--|--|
| **Core** | Z≈6 (+ optional N); prefer **Z6N6** parked droplet |
| **Cloud** | \(n_L=Z\); private bag; Coulomb via \(A\) |
| **Pass bar** | PASS_nuc + \(c_L\le0.15\) + Gauss + no merge |
| **F19 smoke** | Nuclear park PASS (N6); PASS_atom **False** |
| **Ideal long-T visual C₁₂** | **Not claimed** |

**Gate A:** Package complete; product ideal bar unmet · **BLOCKED** (product) with **CP-A ADOPT** (package ready) — A ADOPT + U ADOPT

Authorization: multi-fab GPU long-T re-run / retune (no new kernel required to start).

---

## 6. Ladder diagram

```text
R0 World (PHYSICS_RELATIONS R1–R10) ── CP-W ──┐
         │                                    │
         ▼                                    │
N  Light gauged Q-ball (ω=1.46, Q≈114)  CP-N  │
         │                                    │
         ▼                                    │
P  B2: p=C+Q , n=C-only               CP-P  │
         │                                    │
         ├──────────► C  Z-carbon c6  CP-C    │
         │               A=12 free: super-crit│
         ▼                                    │
L  Opposite L (q_L=-1, n_L=Z)         CP-L  │
         │                                    │
         ▼                                    │
A  Core(Z)+L=Z ── recipe READY ── product BLOCKED ── CP-A
         │
         ▼
U  CONVERGENCE / STOP                   CP-U
```

---

## 7. Success criteria (structural)

| Check | Required |
|-------|----------|
| Bookkeeping | \(Q\), \(Q_a\), \(\omega\), Gauss flux |
| On-branch / park | Q-ball or fusion branch; c_Q gates |
| Isolation | Private bags for atoms (R8) |
| No fatal radiation | Gauss floor; no uncontrolled death |
| Chemistry | **Not** required |

---

## 8. Open residuals

1. Free A=12 park at g=0.05 impossible (charge budget).  
2. Multi-ball nuclear → droplet (not multi-center crystal).  
3. Multi-rev immortal ± / H orbit open.  
4. F19 atom soft (L −12.5%); ideal C₁₂ not green.  
5. Monist \(\psi\) not inside `scp_sim` particle runs.

---

## 9. Document map

| Doc | Role |
|-----|------|
| This file | Ladder |
| `RECIPES.md` | Commands / configs |
| `CONVERGENCE.md` | Phase scores |
| `PHYSICS_RELATIONS.md` | R1–R10 (W) |
| `work/U/CHECKPOINT_BOARD.md` | CP stamps |
| Peer `work/{W,N,P,C,L,A}/*` | Phase freezes |
