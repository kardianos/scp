# Congruence Verdict Round 2 — Approach D

**Date:** 2026-07-18  
**Inputs:**  
- B R1: `results.json`, `results_kernel.json`  
- B R2: `round2_rays.tsv`, `round2_result.json` (M1/M2/M3/R1)  
- C NC + congruence checklist + kernel dualism stress  
- A F1–F5 / free_response_kernel  
- O-001  
**Code:** `congruence_score_r2.py`, `offline_compute_r2.py`  
**Tables:** `results/round2_*`

---

## 1. Question

> Does any monist numerical approach pass **Occam + fit** AND match **C necessary conditions**?

Goal (2): complete theory + workable numerics, congruent.

---

## 2. Maps scored

| Map | Ground-truth sector | Exterior | Dynamical free? |
|-----|---------------------|----------|-----------------|
| M1 local free optics | monist_1sector | compact | no (static GRIN) |
| **M2 free Laplace 2D** | **monist_1sector** | **log** | **YES** |
| M3 Poisson 2D | dualist_2sector | log | no (gravity solve) |
| R1 Φ=αM/R | dualist_2sector_or_postulated | 1/r | no |

---

## 3. Winner models (combined score \(S\))

| Map | score_winner | \(S\) | fit_winner | Note |
|-----|--------------|-------|------------|------|
| M1 | monist_local_GRIN | ≈0.012 | same | compact monist demo |
| **M2** | **monist_dynamical_log** | **≈0.025** | same | Occam beats dualist_log (L_fit tie) |
| M3 | monist_dynamical_log | ≈0.03 | same | **FALSE POSITIVE** without tags |
| R1 | monist_kernel_postulated | ≈0.308 | dualist free GM | λ_postulate; non-iso log loses L_fit |

### Critical comparisons

**M2 (dynamical monist) vs dualist adversaries**

| Model | \(L_{\mathrm{fit}}\) | \(S\) |
|-------|----------------------|-------|
| monist_dynamical_log (1 sector) | ≈0.025 | **≈0.025** |
| dualist_log_2D (2 sector) | ≈0.025 | ≈1.525 |
| monist_kernel 1/r postulated | ≈0.48 | ≈0.78 |
| dualist Plummer 1/r | ≈0.48 | ≈1.98 |

- **Occam separates** monist log from dualist log when L_fit ties (M2 vs dualist twin class).  
- **Pure L_fit separates** dynamical log (M2) from 1/r dualist/postulated kernel.  
- Dynamical free-response is **not** ray-isomorphic to R1/Plummer.

**M3 (dualist Poisson) trap**

| Model | \(L_{\mathrm{fit}}\) | \(S\) |
|-------|----------------------|-------|
| monist_dynamical_log | ≈0.03 | ≈0.03 ← wins score |
| dualist_log_2D | ≈0.03 | ≈1.53 |

Rays + Occam **alone mis-award monism** on dualist M3 data (same multipole as M2).  
**Recovery requires** `sector_tag` / `gravity_solver` from B exports (B-009 FOR_D).  
softE when Poisson is *labeled* monist without free-origin: \(S\sim 100\).

**R1 postulated 1/r**

| Model | \(L_{\mathrm{fit}}\) | \(S\) |
|-------|----------------------|-------|
| monist_kernel_postulated | ≈0.008 | ≈0.308 |
| dualist_plummer (iso) | ≈0.008 | ≈1.508 |
| dualist_log (non-iso) | ≈0.58 | ≈2.08 |

Non-isomorphic dualist **loses pure fit** on 1/r data (Round-2 D3 advance).  
Isomorphic Plummer still ties L_fit (Occam).  
B `monist_kernel_failed` for dynamical 1/r: **confirmed** — M2 is log, not 1/r.

---

## 4. C necessary conditions / checklist

| ID | Status | Evidence |
|----|--------|----------|
| N2 free deficit | PASS | 0.151 all locks |
| N3 no gravity_solver (monist maps) | PASS | M1/M2 null |
| N5 local not overclaimed | PASS | M1 compact only |
| N6 exterior monopole | PARTIAL | M2 log long-range; not 1/r |
| **N8 free dynamics produce K** | **PARTIAL→PASS (log)** | **M2 free Laplace dynamical** |
| N8 Einstein-class 1/r from free dyn | **FAIL** | B monist_kernel_failed |
| N9 dualist adversarial | **PASS** | non-iso L_fit + Occam + softE |
| J4 K ontology | PARTIAL | M2 free-origin OK; R1 postulated not |
| J5 inertia triad | UNTESTED/FAIL | B deferred |
| J6 dualist separation | PARTIAL | works M2 vs 1/r; fails M2 vs M3 without tags |
| J7 softE | PASS | S≈100 |
| **C_package_monist_complete** | **FAIL** | log≠C 1/r target; triad open; M3 false-positive issue |

C FOR_D (λ_deriv / postulated): adopted as λ_postulate=0.3 on R1.  
C FOR_D checklist: reported above.

---

## 5. Answer: any monist pass Occam+fit AND C NCs?

| Candidate | Fit | Occam | C long-range 1/r | Dynamical free K | Overall |
|-----------|-----|-------|------------------|------------------|---------|
| M1 local GRIN | PASS | PASS | FAIL | no | Compact demo only |
| **M2 free Laplace** | PASS | PASS | **FAIL (log)** | **YES** | **Best monist numeric; wrong multipole for GR target** |
| R1 postulated kernel | PASS | PASS* | PARTIAL | no | Phenomenology; not monist theory |
| M3 dualist | PASS (rays) | — | FAIL log | no | Dualist baseline |

\*with λ_postulate honesty.

**Strict answer:** **No** approach passes full C Einstein-class package + complete monist K + triad.

**Qualified answer:** **M2 is the first monist numerical free-response that:**
1. is dynamical (not hand Poisson),  
2. wins Occam+fit on its own data,  
3. loses L_fit when forced into wrong 1/r dualist class,  
4. fails C’s **1/r** exterior NC (delivers **2D log** instead).

That is real progress on O-001 (“make free-response monist or kill it”):  
- **Monist free-response: partially MADE (M2 log).**  
- **Einstein-class 1/r monist kernel: KILLED as dynamical free-only claim in 2D.**

---

## 6. Goal condition (2)

| | |
|--|--|
| **goal2_closer** | **TRUE** |
| **goal2_met** | **FALSE** |

**Closer because:** dynamical monist map exists; D non-iso dualist works; theory (A F1/F3) has matching log-vs-1/r honesty; C NC-K stress-test aligned.

**Not met because:**  
1. C weak-field target \(\ell\sim M/r\) not matched by M2 (log).  
2. M2/M3 multipole degeneracy needs metadata for dualist rejection.  
3. Inertia triad absent.  
4. Theory free-response not fully derived/closed (A still open on inertia coefficient).

**Path to goal (2):** B 3D free Laplace (Green~1/r) **or** A axiomatic K with B export `phi_origin=free_relaxation` producing 1/r; D re-score; C J4/J5 PASS.

---

## 7. Cross-approach notes

**FOR_A:** Do not award congruence to R1/F5. M2 is F3-like (free Laplace/graph Green) in 2D — log multipole expected. Ax8 1/r needs 3D free Green or axiomatic K. D will score `phi_origin=free_relaxation` without λ_postulate.

**FOR_B:** Ingest complete. M2 vs M3: Occam alone insufficient; keep `sector_tag`+`gravity_solver` mandatory. Next 3D free Laplace for 1/r monist attempt. Inertia independent measurement still needed for J5.

**FOR_C:** J6 strengthened for non-iso classes; weakened for iso multipole twins (M2/M3). Recommend NC: exports must carry K_origin / gravity_solver for dualist rejection when multipoles match. N8 partial PASS for free dynamics (log).

---

## 8. Bottom line

1. **Ingested B R1+R2** maps including dynamical M2.  
2. **M2 monist dynamical log wins** Occam+fit on monist data; **1/r models lose L_fit** on M2.  
3. **M3 dualist** is a multipole twin — **false monist award** without tags.  
4. **R1/kernel** still non-monist-derived; non-iso log dualist loses on 1/r data.  
5. **Goal (2) closer, not met** — best monist numeric is M2 (log free response), not Einstein-class kernel.
