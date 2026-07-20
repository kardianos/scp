# World Freeze v0 — Field Primitive Stack (N → P → C → atom)

**Agent:** W  
**Date:** 2026-07-19  
**Status:** **FROZEN for CP-W** (ADOPT)  
**Parent:** `/home/d/code/scp/v78/PHYSICS_RELATIONS.md`  
**Campaign:** `/home/d/code/scp/v78/GOALS.md`

---

## 0. One-paragraph freeze

The v78 world is a **single continuum** with free and bound occupancy. **Particles** are stable localized bound configurations (Q-balls and multi-fabric composites) engineered in the SCP gauged complex Cosserat kernel. **Monist free channels** (capacity \(\psi\), Maxwell) supply locality, mass-energy language, Gauss/Coulomb ontology, and optional path-cost geometry. The particle ladder **N → P/N⁰ → C nucleus → L cloud → atom** is built only from field primitives listed here — no imported chemistry species.

---

## 1. Primitive stack (layers)

```text
Layer 0  Continuum / locality c
Layer 1  Matter fields  Φ_a, Θ_a  (complex Cosserat; × fabrics)
Layer 2  Short-range bag  V_t(s_f)  private per fabric
Layer 3  Shared gauge  A, E  (diagonal U(1); multi-fabric q_f)
Layer 4  Conserved ledgers  Q, Q_a, ω;  Q_em = Σ q_f ρ_f
Layer 5  Free-capacity ψ  (monist; parallel bookkeeping — not in scp_sim)
Layer 6  Stability gates  VK / park / Gauss / long-T
```

| Layer | Symbols | Source | Nuclear-T required? |
|-------|---------|--------|---------------------|
| 0 | \(c\), hyperbolic update | monist R1 + kernel CFL | **yes** |
| 1 | \(\Phi_a=u_a+iv_a\), \(\Theta_a\) | CONCEPT / kernel-v3 | **yes** |
| 2 | \(s_f=\prod|\Phi_f^a|^2\), \(V_t\) | CONCEPT; multi-fab private | **yes** |
| 3 | \(A_i,E_i\), \(g\), \(q_f\) | v69+; v75 charges | **yes** (EM species) |
| 4 | \(Q,Q_a,\omega,Q_{\mathrm{em}}\) | Noether + Gauss | **yes** |
| 5 | \(\psi,\sigma,\rho_b,\rho_f\) | v76 F1-3D; v77 RC1 | **no** (R4 optional) |
| 6 | VK, park, `gauss_max` | measured / diagnostics | **yes** |

---

## 2. Standard parameter freeze (particle package)

| Param | Value | Role |
|-------|-------|------|
| \(m^2\) | 2.25 | matter mass gap \(m=1.5\) |
| \(\mu,\kappa\) | −41.345, 50 | saturating \(V_t\) |
| \(m_\theta\) | 1.6 | close θ-drain (\(> \omega\)) |
| \(\eta\) | 0 | standard particle path |
| \(g_{\mathrm{gauge}}\) | 0.05 | Coulomb + \(Q_{\max}=921\) |
| complex_phi / complex_gauge | 1 / 1 | U(1) era |
| BC | absorbing | compact objects |

**Q-ball windows:** ungauged \(\omega\in(1.3087,1.5)\); gauged \(g=0.05\): \(\omega\in(1.406,1.5)\), \(Q_{\max}=921\).

**Light nucleon default:** \(\omega=1.46\), \(Q_N\approx 114\) (`f_w146_g005`).

---

## 3. Multi-fabric charge freeze (atom path)

When multi-fabric (v75 Option B) is available:

| Fabric | Role | \(q_f\) | Bag \(s_f\) |
|--------|------|---------|------------|
| **C** | nuclear bag / binding | 0 | yes |
| **Q** | nuclear EM charge | +1 | optional self-bag |
| **L** | light opposite sector | −1 | private; never \(s_C\) |

\[
\rho_{\mathrm{em}}=q_C\rho_C+q_Q\rho_Q+q_L\rho_L.
\]

**B2 unlock** (`mf_lock_CQ=0`) required for true P/N (independent \(\Phi_Q\)).

---

## 4. Ladder map (primitives → targets)

### 4.1 Nucleon (N)

| Need | Primitive |
|------|-----------|
| Existence | L1 equal-\(\Phi_a\) Q-ball ansatz \(f(r)e^{i\omega t}\) |
| Binding | L2 \(V_t(s)\) + charge pressure |
| EM self-energy | L3 gauged diagonal U(1) |
| Inventory | L4 \(Q_N\), \(E/Q\), \(r_{1/2}\) |
| Survival | L6 VK \(\mathrm{d}Q/\mathrm{d}\omega<0\), long-T |

**Monist gloss:** N = compact lock (bound ledger); free channels carry radiation only.

### 4.2 Proton / neutron (P / N⁰)

| Species | Field content | R3 EM |
|---------|---------------|-------|
| **P** | C bag + Q co-located | \(Q_{\mathrm{em}}\neq0\) |
| **N⁰** | C bag only | \(Q_{\mathrm{em}}\approx0\) |

| Need | Primitive |
|------|-----------|
| Isolation of bag vs charge | L2 private + L8 multi-fab B2 |
| EM species split | L3–L4 fabric \(q_f\) |
| Not sufficient alone | R7 flavor multiplet |

### 4.3 Carbon nucleus (C)

| Map | Primitive path |
|-----|----------------|
| **Z-carbon** (primary) | 6× light N co-phase fusion → parked droplet (c6_light class, \(Q\sim650\)) |
| **A=12 free** | 12× N → super-critical vs \(Q_{\max}\); evaporates (documented fate) |
| Binding | R5 \(V_t\) + phase coherence |
| Charge inventory | R3 \(Z\) via P-count or \(Q_{\mathrm{em}}\) proxy |

### 4.4 Light sector (L)

| Need | Primitive |
|------|-----------|
| Opposite charge | L fabric \(q_L=-1\) |
| No nuclear merge | private \(s_L\); no \(s_C\) engagement |
| Scale | lighter \(\omega\) / profile hierarchy vs nuclear |

### 4.5 Atom (C₁₂-scale)

| Piece | Composition |
|-------|-------------|
| Core | Z≈6 nuclear (P count) + optional N⁰ (isotope) |
| Cloud | 6× L opposite charges |
| Bind | **R6 Coulomb via shared \(A\)** only |
| Forbid | bag merge (R8 fail), fatal L radiation (R10 fail) |
| Block today | multi-fabric kernel path + L long-T (Stage 3); ask before kernel edit |

---

## 5. Relation coverage checklist

| Relation | Carried by freeze layers |
|----------|--------------------------|
| R1 Locality \(c\) | L0 |
| R2 Mass ↔ energy | L1–L2 energy \(T_{00}\); monist L5 mass language |
| R3 Charge ↔ Gauss | L3–L4 |
| R4 Path-cost | L5 optional |
| R5 Nuclear binding | L2 (+ multi-ball L1) |
| R6 Coulomb | L3 |
| R7 Flavor | L1 + L4 \((\omega_a,Q_a)\) |
| R8 Multi-fabric | L1×fabrics + L2 private + L3 shared |
| R9 Free/bound | L0 continuum + bound lumps; monist \(\rho_f+\rho_b\) |
| R10 Stability | L6 |

**All R1–R10 assigned.** Evidence updates belong to N/P/C/L/A/U; identity freeze is this document + `PHYSICS_RELATIONS.md`.

---

## 6. Explicit non-primitives (out of freeze)

- Real-field braid/oscillon “particles” (historical ≤v53; unstable)
- Chemistry import (orbitals by hand, external e⁻ species)
- \(\psi\) as required dynamical field inside nuclear recipes
- Hand-placed point masses without field content
- Kernel edits without human authorization

---

## 7. Downstream contracts

| Phase | Must respect freeze |
|-------|---------------------|
| **N** | Standard package params; equal-component gauged Q-ball; VK bar |
| **P** | B2 P/N definition; R3 not R7 for EM neutrality |
| **C** | Z-carbon primary; A=12 honest super-critical note |
| **L** | \(q_L=-1\); isolation from \(s_C\) |
| **A** | Core + 6L; Coulomb-only bind; document BLOCKED-AUTH if kernel missing |
| **U** | Score phases against R1–R10 matrix; co-stamp CP-W |

---

## 8. Stamp

**WORLD_FREEZE_v0 = ADOPT** (W, 2026-07-19)  
Revision requires W re-stamp and U co-agree on CP-W.
