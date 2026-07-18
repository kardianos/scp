# Congruence Verdict Round 3 — Approach D (3D free-response)

**Date:** 2026-07-18  
**Orchestrator:** O-003 / O-004  
**B data:** `outputs/round3_{rays,path_cost,free_deficit,result}.*` **INGESTED**  
**Code:** `congruence_score_r3.py`  
**Tables:** `results/round3_*`

---

## 1. Question

> Does monist **3D free response** pass **fit + Occam + C checklist** enough to claim **partial goal (2)**?

C v1 tiers (C-015):

| Tier | Bar |
|------|-----|
| **goal2_minimal** | Dynamical monist free-response (2D log OK) |
| **goal2_PC3D** | + 3D free Green 1/r + PC-3D exterior (O-003) |
| **goal2 full** | PC3D + inertia triad + closed theory congruence J* |

---

## 2. B maps scored

| Map | sector_tag | phi_origin | multipole \(\|\alpha\|_4/\|\alpha\|_1\) | \(M\) | free deficit |
|-----|------------|------------|----------------------------------------|-------|--------------|
| **B monist 3D free** | monist_1sector | free_capacity_3d_green | **0.250 = inv_r_3d** | 6.300 | 0.174 |
| B dualist twin | dualist_2sector | dualist_poisson_label | 0.250 (iso) | 6.300 | 0.174 |

B reports: `monist_3d_1r_pass=True`, `gravity_solver=None` on monist channel,  
\(\psi\propto 1/r\) with \(R^2(1/r)=1\), Born rays \(\mathrm{defl}\sim -2\alpha_M/b\).

---

## 3. Winner models

### 3.1 Monist 3D free-capacity map

| Model | \(L_{\mathrm{fit}}\) | \(S\) |
|-------|----------------------|-------|
| **monist_3d_dynamical_free** | ≈0.005 | **≈0.005** |
| dualist_3d_poisson | ≈0.005 | ≈1.505 |
| dualist_2d_log | ≈0.58 | ≈2.08 |
| monist_local_GRIN | ≈0.95 | ≈0.95 |
| softE Poisson-as-monist | ≈0.005 | **≈100** |

**Fit:** monist 3D recovers B rays (weak \(1/r\) Born).  
**Occam:** beats dualist 3D Poisson despite **L_fit tie** (iso Green).  
**Non-iso:** 2D log **loses pure L_fit** (dimension diagnostic).  
**softE:** dead.

### 3.2 Dualist twin (same rays)

Rays + Occam alone still **false-award monist** (3D Green isomorphism).  
**Recovery:** `sector_tag=dualist_2sector` + `gravity_solver=poisson_3d_tagged_dualist` (B export).  
softE if dualist field labeled monist without dualist tag: \(S\sim 100\).

---

## 4. C checklist / goal tiers

| Item | Status | Evidence |
|------|--------|----------|
| N2 free deficit | **PASS** | 0.174 |
| N3 no gravity_solver (monist) | **PASS** | null |
| N6b inv_r_3d | **PASS** | multipole 0.250; B \(R^2=1\) |
| N8 free dynamics / free-origin K | **PASS (partial)** | phi_origin=free_capacity_3d_green; analytic Green class (SOR code present) |
| N9 dualist adversarial | **PASS** | Occam + non-iso log + softE |
| N13 3D free Green measured | **PASS** | B path_cost + multipole_fit |
| J3b PC-3D exterior | **PASS (weak Born)** | defl∝1/b, delay falloff |
| J4 K ontology | **PARTIAL→PASS** | free-capacity state; dualist twin tagged |
| J5 inertia triad | **FAIL** | still deferred |
| J6 iso twins | **PARTIAL** | tags required |
| J7 softE | **PASS** | |
| **goal2_minimal** | **PASS** | R2 M2 + continues |
| **goal2_PC3D (partial)** | **PASS (qualified)** | see §5 |
| **goal2 full** | **NOT MET** | J5 + strong Einstein coefficients optional |

---

## 5. Partial goal (2) claim — **YES (qualified)**

### Claim

**Partial goal condition (2) / goal2_PC3D is claimable as a project milestone:**

There exists a **congruent pair**:

1. **Theory (A+C):** F1-3D free-capacity / free Laplace → exterior Green \(1/r\);  
   \(G_{\mathrm{eff}}\) from free constants; M2 mechanism endorsed; dimension note (2D log vs 3D 1/r).  
2. **Numerics (B+D):** 3D free-capacity \(\psi\) with budget identity; exterior \(\psi\propto 1/r\);  
   rays without second gravity solver; monist_1sector tags; D score_winner monist with  
   Occam over dualist 3D Poisson and non-iso 2D log loss on \(L_{\mathrm{fit}}\).

### Qualifications (honest)

1. **Green isomorphism remains:** dualist 3D Poisson is ray-identical; monism is **not** proved by fit alone — sector tags + free-origin story are part of the claim.  
2. **Inertia triad (J5) open** → full goal2 / goal2_strong **not** met.  
3. **Analytic exterior + Born** primary in B export; full SOR grid code exists for parent confirmation — not a dualist bolt-on, but not a long dynamical formation history either.  
4. **Einstein numerical factor** (e.g. exact \(4GM/c^2b\)) is SHOULD, not MUST for partial PC3D.

### Not claimed

- Full PROBLEM.md monist gravity of the cosmos  
- Soft-Poisson as monism without free state  
- 2D log as GR-class exterior  

---

## 6. Bottom line

| Question | Answer |
|----------|--------|
| Monist 3D free response passes fit+Occam+C 1/r multipole? | **YES** (B maps + D scores) |
| Partial goal (2) / goal2_PC3D? | **YES — qualified milestone** |
| Full goal (2)? | **NO** (inertia + stronger theory formalization) |
| Dualist 3D Poisson beaten on pure fit? | **NO (iso)** — beaten on **Occam + tags** |
| 2D log dualist beaten on pure fit? | **YES** |

**O-003 path closed for multipole dimension:** free Laplace in **3D** delivers **1/r**; Round 2 “failure” was embedding dimension, not monist free-response death.

---

## 7. Cross-approach

**FOR_A:** Numeric congruence with F1-3D: \(\alpha_{\mathrm{eff}}=\gamma\kappa/(4\pi)=0.03979\), exterior \(1/r\). Publish \(G_{\mathrm{eff}}\) match note.

**FOR_B:** Ingest complete. Keep dualist twin exports. Next: independent inertia (J5); optional full SOR vs analytic residual table.

**FOR_C:** goal2_PC3D partial PASS on N6b/N13/J3b; J5 FAIL. Multipole + iso-twin lessons confirmed on real B data.

**FOR_O:** Recommend marking **goal2_PC3D partial** achieved; full goal (2) waits on inertia triad + any remaining T* formal items.
