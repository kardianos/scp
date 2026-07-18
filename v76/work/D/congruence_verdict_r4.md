# Congruence Verdict Round 4 — Approach D (final)

**Date:** 2026-07-18  
**Theme (O-005/O-006):** Close J5 inertia + final goal-(2) view  
**Code:** `congruence_score_r4.py`  
**Tables:** `results/round4_*`  
**Prior:** R3 3D free-capacity rays retained; B `round4_*` **ABSENT** at D close  

---

## 1. What was scored

| Package | Source | Status |
|---------|--------|--------|
| 3D free-capacity rays + multipole | B `round3_*` | **Kept / re-scored** |
| Dualist 3D Poisson twin | B R3 dualist rows | **Kept** |
| Non-iso 2D log adversary | D synthetic multipole stress | **Kept** |
| softE | D | **Kept** |
| Independent inertia triad | B round4 | **ABSENT** |
| Ray↔ledger mass | B R3 \(\alpha_M/\alpha_{\mathrm{eff}}\) | **Present** |

---

## 2. Occam + fit — monist still preferred?

### Ray scores (B monist 3D free map)

| Model | \(L_{\mathrm{fit}}\) | \(S\) | Prefer? |
|-------|----------------------|-------|---------|
| **monist_3d_free_capacity** | ≈0 | **≈0** | **YES** |
| dualist_3d_poisson (iso) | ≈0 | ≈1.5 | no (Occam) |
| dualist_2d_log (non-iso) | large | large | no (fit+Occam) |
| softE Poisson-as-monist | ≈0 | **100** | killed |

**Answer:** **Yes.** Monist package remains preferred under combined score.  
**Honest residual:** pure \(L_{\mathrm{fit}}\) **ties** monist free-capacity vs dualist 3D Poisson (same Green). Preference is **Occam + free-origin tags + free deficit + null gravity_solver**, not ray residual alone.

### Multipole / free deficit / solver

| Check | Result |
|-------|--------|
| multipole \(\|\alpha\|_4/\|\alpha\|_1\) | **0.250** inv_r_3d |
| free_deficit_core | **0.174 > 0** |
| gravity_solver (monist) | **null** |
| B monist_3d_1r_pass | **True** (O-005 parent-verified) |

---

## 3. Inertia triad (J5)

| Quantity | Value |
|----------|--------|
| \(m_{\mathrm{ledger}}\) | 6.299787 |
| \(m_{\mathrm{ray}}=\alpha_M/\alpha_{\mathrm{eff}}\) | 6.299787 (**rel err 0**) |
| \(m_{\mathrm{inertial}}\) independent | **None** (R4 export absent; R2 deferred tautology guard) |

| J5 status | **PARTIAL** |
|-----------|-------------|
| Independent inertial | **FAIL / ABSENT** |
| Ray↔ledger | **PASS** |

---

## 4. Final D view on goal (2)

Aligned with C-018/C-019 tiers and O-005 residual list:

| Tier | D verdict | Notes |
|------|-----------|--------|
| **goal2_minimal** | **PASS / MET** | Dynamical free-response monism (M2 log + R3 continues) |
| **goal2_PC3D_partial** | **PASS** | 3D free Green 1/r + rays + deficit + D Occam |
| **goal2_PC3D_workable** | **PASS / MET_QUALIFIED** | Theory F1-3D + numerics congruent on free-response path-cost channel; J5 deferred optional (O-005) |
| **goal2_full / complete** | **NOT MET** | Independent inertia triad open |
| **PROBLEM.md full monist bar** | **NOT claimed** | Weak-field skeleton only; residuals public |

### Does monist package remain preferred with residuals honest?

**Yes.**

Preferred because:

1. Fit to B 3D free-capacity rays at machine-level weak Born consistency  
2. Non-isomorphic dualists (2D log, compact GRIN) lose on \(L_{\mathrm{fit}}\)  
3. Isomorphic dualist loses on **ontology score** (sectors, free–bound link)  
4. softE disqualifies Poisson labeled monist  
5. Free deficit + free-origin tags present on monist channel  

Not preferred *only* because rays look good — that would fail X4 (ray-fit-only monism).

---

## 5. Honest residuals (public)

1. **Poisson-form isomorphism:** free-capacity PDE has Poisson math; monism is free/bound identity + single sector + free-origin, not a magically different equation.  
2. **Fit does not separate monist vs dualist 3D twin** — tags/Occam required.  
3. **J5 independent inertia not demonstrated** — ray↔ledger only.  
4. **Analytic Green + Born dominate B R3 export**; mini SOR parent-checked (~0.93).  
5. **PPN / Einstein factor-4 not claimed.**  
6. **Scope:** weak-field free-response skeleton, not full gravity phenomenology.

---

## 6. Alignment with A/B/C/O

| Party | Alignment |
|-------|-----------|
| O-005 | PARTIAL goal (2); PC3D package works; J5 residual — **agree**; workable tier **YES** |
| C-018/019 | goal2_PC3D_workable MET_QUALIFIED; full NOT — **agree** |
| A F1-3D | Numeric congruence on exterior 1/r + \(G_{\mathrm{eff}}\) sketch — **agree** |
| B R3 | monist_3d_1r_pass — **re-scored PASS**; R4 inertia not shipped — **J5 PARTIAL** |

---

## 7. Bottom line (quotable)

> **D final:** The monist F1-3D / free-capacity package **remains preferred** under Occam+fit with residuals stated.  
> **goal2_minimal = MET.**  
> **goal2_PC3D_workable = MET_QUALIFIED** (ray↔ledger OK; independent inertia deferred).  
> **goal2_full = NOT MET.**  
> Do **not** claim PROBLEM.md monist fabric fully proved.

---

## 8. Files

- `/home/d/code/scp/v76/work/D/congruence_score_r4.py`
- `/home/d/code/scp/v76/work/D/congruence_verdict_r4.md`
- `/home/d/code/scp/v76/work/D/results/round4_*`
