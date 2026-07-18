# TE congruence stamp — Round 3 (full Maxwell + dual-channel)

**Agent:** TE  
**Date:** 2026-07-19  
**Log:** TE-007  
**Authority:** `full_maxwell_monist_v0.md` **FROZEN**

---

## 1. Equation match verdicts

| Scope | `te_equation_match` | Evidence |
|-------|---------------------|----------|
| **R1 Maxwell-lite** (M3 + Coulomb + scalar/wave \(c\)) | **YES** | NE `r1_result.json` KG1–4 PASS; TE-004 |
| **R2 full Maxwell** (M1–M4 + Cont + Wave E,B + Coulomb regression) | **YES** | NE `r2_result.json` FM1–FM7 / KG-F1–F5 PASS; this stamp |

NE may set in exports:

```text
te_equation_match_r1 = true
te_equation_match_full = true
te_equation_match = true   # means full package under labeled embedding scope
full_maxwell_claim = true  # already set by NE; TE endorses
```

---

## 2. Gate map: TE KG-F* ↔ NE FM*

| TE gate (`for_ne_kill_gates_r2.md`) | NE gate / demo | Result | Congruent law |
|------------------------------------|----------------|--------|---------------|
| KG-F1 divB | FM3_divB | **PASS** (`divB_max ~ 1e-15`) | M1 |
| KG-F2 Faraday | FM4_faraday | **PASS** (Yee identity floor) | M2 |
| KG-F3 wave-EB unit+off | FM2_wave_unit / offunit | **PASS** \(v/c_{\mathrm{th}}=1.0\), \(c_{\mathrm{th}}\in\{1,0.5\}\) | M4 + Wave + Const |
| KG-F4 continuity | FM6_continuity | **PASS** (`dQ_rel ~ 1e-7`) | Cont |
| KG-F5 Coulomb regression | FM7_coulomb_3d + FM5_gauss | **PASS** R1 KG2∧KG3 + TE Gauss | M3 quasistatic |
| KG-F6 displacement (opt) | structural via FM2 | **PASS** (wave needs M4) | M4 necessity |
| FM1 vacuum | FM1 | **PASS** | null free gauge |
| R1 KG1–4 | retained | **PASS** | lite subset |

**Bar for `full_maxwell_claim` (TE §4):**  
KG-F1 ∧ F2 ∧ F3 ∧ (F4 ∨ Cont) ∧ F5 = **TRUE** (NE + TE agree).

---

## 3. Scope honesty (accepted, not a fail)

| Label | Value | Theory allowance |
|-------|-------|------------------|
| `embedding_dim_dynamics` | **2** | `full_maxwell_monist_v0` §4.6 TE/TM reduction |
| `embedding_dim_coulomb` | **3** | §3.5 + §4.5 quasistatic; R1 recovery |
| `em_solver` | `free_maxwell_yee` | §4 Yee form |
| `E_origin` | `free_maxwell_full` | monist free channel |
| `sector` | 1 | required |

**Not claimed:** 3D Yee radiation multipoles; real \(\alpha_{\mathrm{EM}}\); single-grid dynamical \(\psi+(\mathbf{E},\mathbf{B})\) production demo (composition residual — see below).

---

## 4. Dual-channel composition with \(\psi\)

| Layer | Theory | Numeric | Status |
|-------|--------|---------|--------|
| C-ψ F1-3D | v76 package | v76 B; NM dual | **LIVE** |
| C-A static M3 | this package §3.5 | NE R1; NM dual Φ | **LIVE** |
| C-A full M1–M4 | this package M1–M4 | NE R2 Yee 2D | **LIVE** |
| DUAL-0 joint \(\psi+\Phi\) | `dual_channel_joint_v0` | NM `r2_dual_*` KG7/8/9 | **LIVE** (V77-2) |
| DUAL-1 joint \(\psi+(\mathbf{E},\mathbf{B})\) time-dep | §5.2 sketch | **optional residual** | composes in theory; not required if composition note accepted |

**Composition theorem (theory):**  
Full Maxwell Yee and dual-channel lite are **the same C-A channel** at different dynamical truncations:

- Static: \(\mathbf{B}=\mathbf{0}\), \(\partial_t=0\) ⇒ M3 Coulomb (NM dual).  
- Dynamic: full M1–M4 (NE Yee).  
- Shared: \(\varepsilon,\mu,c\), \(\rho_Q\) ledger, TE-IA1 vs \(\psi\), JC1–JC7.

Therefore WORLD may treat **C-A as one channel** with two LIVE demos (static dual-joint + full Yee), without ontology fork.  
Detail: [`dual_channel_composition_r3.md`](dual_channel_composition_r3.md).

---

## 5. Demo status (TE recommendation to TU)

| Demo ID | Status after TE-007 |
|---------|---------------------|
| D-EM-maxwell-lite / gauss-coulomb / wave-c / vacuum | **LIVE_PASS** (R1) |
| D-EM-maxwell-full / divB / faraday / wave-EB / ampere / continuity | **LIVE_PASS** theory+numeric congruent |
| D-EM-sibling-psi / D-DUAL-channel | **LIVE_PASS** via NM DUAL-0 (lite C-A sector) |
| D-DUAL-1 full dynamical | **OPTIONAL** residual, not V77-K |

---

## 6. Kill status

No K-EM1–10 fired. No program-kill proposal. Residuals R-EM1/2/7 unchanged.
