# Unified Monist Package v1 — Cohesive Theory + Model Stack

**Agent:** TU · **Round:** 3 · **Date:** 2026-07-18  
**Status:** **cohesive agreement package** (composition tier)  
**Seed:** v76 F1-3D `goal2_PC3D_workable`  
**Stop intent:** Recommend campaign terminal **(A)** at this tier — see `UNIFICATION_SCORE_r3.md`  
**Not claimed:** Nature is monist; full GR/PPN/cosmology; real \(G,\alpha_{\mathrm{EM}}\); one production binary co-evolving full Maxwell + ψ on locks

---

## 0. One-paragraph claim

There is **one continuum**. Bound locks hold mass-form ledger \(\rho_b\) and optional gauge ledger \(\rho_Q\). Free sibling channels on that continuum are: **free capacity** \(\psi\) (F1-3D path cost), **free gauge** \((\mathbf{E},\mathbf{B})\) (full classical Maxwell), and free radiation of either. Locality \(c\) is shared: \(c_{\mathrm{path}}=1/\sqrt{\varepsilon\mu}\). Mass for path cost is \(M=E_\star^{\mathrm{path}}/c^2\) with \(E_\star^{\mathrm{path}}=\int\rho_b\); inertial mass defaults to free stress \(m_{\mathrm{inertial}}=\xi U/c^2\) (**J5-β**). Numerics demonstrate each law with kill-gates and compose under one ontology and constitutive table. **Composition residual:** dual-channel multi-lock forces use **static Maxwell-lite** \(\Phi\) on the joint lock medium; **dynamical full Maxwell** is a **separate** Yee sandbox (locks/ψ not co-evolved there). Congruence is **interface + shared \(c\) + dual-source independence**, not a single joint binary.

---

## 1. Shared primitives

| Primitive | Definition | Evidence / owner |
|-----------|------------|------------------|
| **Continuum** | Only medium; free + bound states | WORLD_v1; all foci |
| **Budget** | \(\rho_f+\rho_b=\rho_{\mathrm{tot}}\) (strong or integral) | v76; NM deficit PASS |
| **Lock** | Compact free-depleting mass-form; may carry \(Q\) | TM lock_ontology; NM G0 |
| **Free capacity \(\psi\)** | Path-cost free DOF | v76 F1-3D LIVE |
| **Free gauge \(\mathbf{E},\mathbf{B}\)** | Maxwell free DOF | TE full_maxwell; NE Yee LIVE |
| **Charge \(\rho_Q,\mathbf{J}_Q\)** | Lock–gauge ledger (signed) | TE Cont; NE FM6; TM Supp |
| **\(c\)** | Free locality = free-signal bound | JC1; ND dual-c PASS; NE off-unit |
| **Energy** | Continuum ledger (bound + free stresses \(U[\psi]\), \(u_{\mathrm{EM}}\)) | TD/ND J5-β |
| **Mass (path)** | \(M=\int\rho_b/c^2\) → multipole \(\psi\) | F1; ray=ledger |
| **Mass (inertial)** | \(m=\xi U/c^2\) (J5-β); naïve \(\int\rho_b\) killed | TD lock; ND PASS |
| **Warp / path cost** | \(\ell=\ell_0+\gamma\psi\) at fixed local \(c\) | v76 rays |

**Forbidden collapses:** \(\psi\equiv\Phi_{\mathrm{EM}}\); \(\rho_Q\equiv\rho_b\) always; local GRIN as monist gravity; hand Poisson as monist *proof*; foreign dualist stage; idle free DOF labeled monist.

---

## 2. Free capacity — F1-3D (v76 seed)

\[
-\nabla\cdot(\sigma\nabla\psi)=s\rho_b,\qquad
\ell=\ell_0+\gamma\psi,\qquad
\psi\sim\frac{s E_\star}{4\pi\sigma_0 r}\ \text{(3D)}.
\]

| Item | Status |
|------|--------|
| Theory | `v76/work/A/THEORETICAL_PACKAGE_v1.md` |
| Numeric | D-ψ-3D LIVE — multipole \(1/r\), free deficit, rays, Occam |
| Dead | local GRIN; hand Poisson as monist proof |

---

## 3. Full Maxwell free gauge (TE + NE)

### 3.1 Theory (complete set)

From TE `full_maxwell_monist_v0.md`:

| Law | Content |
|-----|---------|
| M1 | \(\nabla\cdot\mathbf{B}=0\) |
| M2 | \(\nabla\times\mathbf{E}+\partial_t\mathbf{B}=\mathbf{0}\) (Faraday) |
| M3 | \(\nabla\cdot(\varepsilon\mathbf{E})=\rho_Q\) |
| M4 | \(\nabla\times(\mathbf{B}/\mu)-\partial_t(\varepsilon\mathbf{E})=\mathbf{J}_Q\) |
| Cont | \(\partial_t\rho_Q+\nabla\cdot\mathbf{J}_Q=0\) |
| Wave | \(c=1/\sqrt{\varepsilon\mu}\); coupled \(\mathbf{E},\mathbf{B}\) |

Ontology: \(\mathbf{E},\mathbf{B}\) = free-gauge continuum stress; charge = lock ledger; sibling to \(\psi\).

### 3.2 Numeric (parent re-verified)

Sandbox: `work/NE/sandbox_full_maxwell_r2.py`  
Export: `work/NE/outputs/r2_result.json`

| Gate | Result | Note |
|------|--------|------|
| FM1 vacuum | **PASS** | max E,H = 0 |
| FM2 wave unit | **PASS** | \(v/c=1.0\); coupled E,B |
| FM2 wave off-unit | **PASS** | \(c_{\mathrm{th}}=0.5\), \(v=0.5\) |
| FM3 div B | **PASS** | \(\sim10^{-15}\) |
| FM4 Faraday | **PASS** | Yee discrete residual \(\sim10^{-16}\) |
| FM5 Gauss static | **PASS** | TE projection |
| FM6 Cont | **PASS** | \(dQ_{\mathrm{rel}}\sim10^{-7}\) |
| FM7 Coulomb 3D | **PASS** | R1 M3 recovery |

**`full_maxwell_claim = true`**

Honesty: dynamics embedding **2D** Yee TE+TM (full 2D Maxwell); Coulomb **3D** recovered from quasistatic sector; not 3D Yee radiation multipoles.

---

## 4. Dual-channel locks (TM + NM)

### 4.1 Theory

STATE-DC / DUAL-1 (TM `dual_channel_v0`, TE `dual_channel_joint_v0`):

\[
\rho_f+\rho_b=\rho_{\mathrm{tot}},\quad
\mathrm{supp}(|\rho_Q|)\subseteq\mathrm{supp}(\rho_b),
\]
\[
-\nabla\cdot(\sigma\nabla\psi)=s\rho_b,\quad
\nabla\cdot(\varepsilon\mathbf{E})=\rho_Q,\quad
c=\frac{1}{\sqrt{\varepsilon\mu}}=c_{\mathrm{path}}.
\]

Forces: \(F=F^\psi+F^{\mathrm{EM}}\) (quasistatic Coulomb Tier A; Lorentz Tier B when \(\mathbf{B},\mathbf{v}\) live).  
Hierarchy: constitutive \(\Xi\lambda_q^2\), not channel deletion.

### 4.2 Numeric joint medium (V77-2)

Sandbox: `work/NM/sandbox_r2_dual_channel.py`  
Export: `work/NM/outputs/r2_dual_result.json`

| Config | \(\psi\) | \(\Phi\) / \(E\) | Forces |
|--------|----------|-----------------|--------|
| vacuum | 0 | 0 | 0 |
| neutral | monopole on | \(\Phi=0\), \(E=0\) | \(F^\psi\) attract only |
| same-sign | monopole | monopole \(E\sim1/r^2\) | \(F^\psi\) attract, \(F^C\) repel |
| opposite | monopole | dipole | both attract |

| Gate | Result |
|------|--------|
| KG7 sibling independence | **PASS** |
| KG8 / JC1 shared \(c\) | **PASS** (\(c=1\)) |
| Budget deficit | **PASS** |
| V77-2 joint numeric | **PASS** |
| `full_maxwell_claim` on this sandbox | **false** (static Maxwell-**lite** \(\Phi\)) |

R1 multi-lock G0–G3 LIVE remain.

---

## 5. Dynamics — J5-β (TD + ND)

| Claim | Status |
|-------|--------|
| J5-β **locked** default | TD `j5_beta_default_v0.md` |
| \(m_{\mathrm{inertial}}=\xi U/c^2\) | DEFAULT |
| \(M_{\mathrm{ray}}=\int\rho_b/c^2\) | HELD (F1) |
| Naïve \(m=\int\rho_b/c^2\) | **KILLED** (ND-G1) |
| Universal \(s^*\) all shapes | **REJECTED** |
| Form factor (baseline) | \(ff\approx 0.141421\) |
| Free wave \(\partial_t^2\psi=c^2\nabla^2\psi\) | D-DYN-free-wave-c **PASS** |
| Shared-\(c\) dual waves | D-DYN-dual-channel-c **PASS** (DC1–DC5) |

**V77-3:** MET under β reading (TU-014).

---

## 6. Shared \(c\) (cross-focus)

| Demo | What | Result |
|------|------|--------|
| NE FM2 | Maxwell wave tracks \(1/\sqrt{\varepsilon\mu}\) unit+off | PASS |
| ND free-wave | free capacity \(v=c\) | PASS |
| ND dual-channel-c | \(c_\psi=c_{\mathrm{em}}\) phase agreement; off-unit both track | PASS |
| NM joint | \(C_{\mathrm{LOCAL}}=1/\sqrt{\varepsilon\mu}\) | PASS |

**JC1** satisfied at composition level.

---

## 7. Demo stack LIVE (authoritative)

| Demo | Focus | Status |
|------|-------|--------|
| D-ψ-3D | seed | **LIVE** |
| D-EM-full-maxwell | EM | **LIVE_PASS** (`full_maxwell_claim=true`) |
| D-EM-gauss-coulomb / vacuum / wave-c | EM | **LIVE_PASS** |
| D-EM-divB / faraday / continuity | EM | **LIVE_PASS** (FM gates) |
| D-DUAL-channel / D-EM-sibling-psi | joint | **LIVE_PASS** (Maxwell-lite joint) |
| D-MAT-lock-S0 / dual0 / force-tax / hier | Matter | **LIVE_PASS** |
| D-DYN-j5-formfactor | Dynamics | **LIVE_PARTIAL** (β closed) |
| D-DYN-free-wave-c | Dynamics | **LIVE_PASS** |
| D-DYN-dual-channel-c | Dynamics+EM | **LIVE_PASS** |

DEAD (stable): local GRIN; hand Poisson monist proof; 2D log as GR exterior.

---

## 8. Model architecture (honest composition)

```text
                    THEORY: one WORLD (this package)
                              │
        ┌─────────────────────┼─────────────────────┐
        ▼                     ▼                     ▼
  [Sandbox A]           [Sandbox B]           [Sandbox C]
  F1-3D ψ / path        Full Maxwell Yee      Joint dual locks
  cost (v76, NM ψ)      2D TE+TM + Cont       ψ + Maxwell-lite Φ
  D-ψ-3D LIVE           D-EM-full-maxwell     D-DUAL-channel LIVE
                        LIVE_PASS             full_maxwell=false
        │                     │                     │
        └──────── shared c, tags, dual-source ──────┘
                  congruence (not one .py binary)
```

| Layer | Content | Single binary? |
|-------|---------|----------------|
| **Theory** | One continuum, JC1–5, M1–M4, F1, J5-β, DUAL-1 | Yes (this doc) |
| **Static dual-channel** | Locks + \(\psi\) + \(\Phi\) forces | Yes (NM r2) |
| **Dynamical Maxwell** | Coupled \(\mathbf{E},\mathbf{B}\) evolution | Yes (NE Yee) — **no locks/ψ** |
| **Shared-c wave probe** | \(\psi\) wave + EM wave same \(c\) | Yes (ND) — linear vacuum |
| **Full joint production** | Locks + \(\psi\) + dynamical \(\mathbf{E},\mathbf{B}\) | **NO** — residual R-compose |

---

## 9. Residuals (do not hide)

| ID | Residual | Blocks terminal (A)? |
|----|----------|----------------------|
| **R-compose** | Dual-channel forces = Maxwell-lite static \(\Phi\); full Maxwell = separate Yee; not one co-evolving joint binary | **No for composition-tier (A)**; **Yes for “one binary” stretch** |
| **R-iso** | Poisson-form isomorphism (ψ and Φ) | Soft — Occam + free-origin tags |
| **R-dim** | Full Maxwell dynamics 2D; Coulomb 3D recovered | Soft — labeled |
| **R-J5α** | Raw triad coefficient 1 deferred | No — V77-3 β closed |
| **R-C1** | Rods/clocks operational incomplete | Soft |
| **R-force** | \(\alpha_\psi\) virtual-work not free-dynamics-derived | Soft |
| **R-scope** | No full GR, BH, real \(G/\alpha\), scp_sim monism proof | Out of campaign bar |
| **R-Lorentz** | Multi-lock Lorentz \(Q(E+v\times B)\) not joint-numeric | Soft Tier B |

---

## 10. Tier summary

| Tier | Verdict |
|------|---------|
| V77-1 multi-demo draft | **MET** |
| V77-2 dual-channel free medium | **MET** (joint ψ + Maxwell-lite + shared \(c\)) |
| V77-3 J5 | **MET** (J5-β) |
| Full Maxwell (human) | **MET** at 2D Yee + 3D Coulomb recovery |
| V77-4 / terminal (A) | **MET at composition tier** (this package) — see score r3 |
| V77-K / (B) | **NO** |

---

## 11. Reproduction (key numerics)

```bash
# Full Maxwell
cd /home/d/code/scp/v77/work/NE && python3 sandbox_full_maxwell_r2.py
# Dual-channel joint
cd /home/d/code/scp/v77/work/NM && python3 sandbox_r2_dual_channel.py
# Shared-c dual waves + J5
cd /home/d/code/scp/v77/work/ND && python3 sandbox_dual_channel_c.py
# Seed path-cost
cd /home/d/code/scp/v76/work/B && python3 offline_round3_3d.py
```

---

## 12. Document control

| Version | Role |
|---------|------|
| WORLD_v0/v1 | Maps / free-channel tables |
| **UNIFIED_PACKAGE_v1** | Terminal cohesive theory+model claim (composition) |
| UNIFICATION_SCORE_r3 | (A)/(B) recommendation |

**Interfaces frozen:** TE M1–M4; TM DUAL-1; TD J5-β; no new ontology forks without reopening (A).
