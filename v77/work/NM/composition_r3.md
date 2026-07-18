# Round 3 — Composition: dual-channel joint + full Maxwell dynamical

**Agent:** NM  
**Date:** 2026-07-18  
**O theme:** UNIFY OR FAIL — compose joint ψ+gauge with full Maxwell LIVE  
**Export:** `outputs/r3_unified_export.json` via `export_r3_unified.py`

---

## 1. What is composed (sibling free-gauge stack)

| Layer | Owner | Content | Claim |
|-------|-------|---------|-------|
| **L1 Dual-channel joint** | NM R2 | Same locks: \(\rho_b\to\psi\) (F1) + \(\rho_Q\to\Phi,\mathbf{E}\) (Maxwell-lite quasistatic) | `V77_2_joint_numeric=PASS`; KG7/KG8; TE-IA1 \(\psi\neq\Phi\) |
| **L2 Full Maxwell dynamical** | NE R2 | Yee 2D TE+TM: time-dep \(\mathbf{E},\mathbf{B}\); Faraday; Ampère–Maxwell; div B; Cont; wave \(v=c\) unit+off-unit; 3D Coulomb recovery | `full_maxwell_claim=true` (FM1–FM7) |
| **Shared \(c\)** | TE JC1 + NE + NM + ND | \(c=1/\sqrt{\varepsilon\mu}=C_{\mathrm{LOCAL}}\) | Unit \(c=1\); off-unit \(\varepsilon=4\Rightarrow c=0.5\) both track |

**Composition rule (honest):**

> Free-capacity \(\psi\) and free-gauge Maxwell are **sibling free channels of one continuum**.  
> NM proves **joint dual-source locks** (static multipoles + forces).  
> NE proves **full dynamical Maxwell** (E,B evolution) as the free-gauge channel.  
> They compose under WORLD / TE-IA1 without identifying \(\psi\equiv\Phi\).  
> A single grid with simultaneous 3D Yee radiation + F1 \(\psi\) is **not** required for congruence: constitutive shared \(c\), independent sources \((\rho_b,\rho_Q)\), and non-contradictory gates suffice for O-004 “compose under WORLD.”

---

## 2. Congruence checks (no new ontology)

| ID | Check | Status |
|----|-------|--------|
| C1 | Same \(C_{\mathrm{LOCAL}}=1\) in NM dual + NE unit wave | PASS |
| C2 | NE off-unit \(c_{\mathrm{th}}=0.5\) tracks \(1/\sqrt{\varepsilon\mu}\) (free-gauge constitutive) | PASS |
| C3 | NM dual \(\Phi/E\) uses same \(\varepsilon_0,\mu_0\) language as NE Maxwell-lite/full | PASS |
| C4 | TE-IA1: neutral mass → \(\psi\neq0\), \(E=0\); opposite Q → \(\psi\) monopole, \(\Phi\) dipole | PASS (NM KG7) |
| C5 | NE full Maxwell does not force \(\psi\) identification or kill path-cost channel | PASS (no \(\psi\) in Yee; sibling attach only) |
| C6 | Budget identity \(\rho_f+\rho_b=\rho_0\) retained on dual-channel joint | PASS |
| C7 | Coulomb regression: NE FM7 recovers R1 3D \(E\sim1/r^2\); NM same-sign \(E\sim1/r^2\) | PASS |

**No contradiction found** between L1 and L2.

---

## 3. What is *not* claimed

1. **Single production sim** with 3D Yee \(\mathbf{E},\mathbf{B}\) + F1 \(\psi\) on identical voxels with radiation + multi-lock forces (deferred residual; optional future).  
2. Real \(G\), \(\alpha_{\mathrm{EM}}\), particle masses.  
3. scp_sim monism proof.

---

## 4. Unified model sketch (for TU)

```text
continuum
  ├── free capacity  ψ, σ, s, γ     ← F1-3D; path cost; mass-form locks ρ_b
  ├── free gauge     E,B (Yee) / Φ  ← Maxwell M1–M4; charge locks ρ_Q
  ├── budget         ρ_f + ρ_b = ρ0
  └── c              = 1/√(εμ) = free locality (shared)
locks
  ├── E★ = ∫ ρ_b  → sources ψ; feels −∇ψ (virtual work)
  └── Q  = ∫ ρ_Q  → sources E; feels qE; Supp|ρ_Q|⊆Supp ρ_b
```

---

## 5. Demo stack (NM view)

| Demo | Status | Source |
|------|--------|--------|
| D-DUAL-channel / D-EM-sibling-psi | LIVE joint numeric | NM R2 |
| D-MAT-dual0 / force-tax | LIVE | NM R1–R2 |
| D-EM-maxwell-full / wave-EB / divB / faraday | LIVE_PASS | NE R2 |
| D-UNIFIED-compose-r3 | LIVE composition package | this note + `r3_unified_export.json` |

---

## 6. Recommendation to O / TU

Composition supports **unified monist model (A-path)** when TU WORLD package cites:

- NM dual-channel joint PASS  
- NE full Maxwell PASS  
- shared \(c\), TE-IA1, no channel collapse  

Residual list for honesty: no single 3D Yee+ψ production grid; J5-β renorm residual (Dynamics); Poisson isomorphism Occam residual (not fatal).
