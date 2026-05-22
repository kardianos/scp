# Phase 6 — Back-Propagation / Evolutionary Search Results

**Date**: 2026-05-20  
**Experiment**: Intensive global optimization over algebraic rigidity terms using Differential Evolution on the relational graph model.

---

## 1. Context and Hypothesis

After completing Phases 1–5, we had strong evidence that the living candidate produces self-stabilizing attractors for effective geometry on purely relational (background-free) graphs, and that the quadratic self-interaction terms plus the ambient modulation `f_g(ρ)` are necessary for that stability.

However, the measured effective dimensionality (`d_eff`) in all previous runs remained stubbornly in the sub-3 range (typically 1.0–1.8). This suggested that the current implementation of protected chirality and growth rules was not yet fully capturing the mechanism that would naturally favor a stable ~3-dimensional causal structure.

The working hypothesis for Phase 6 (proposed after the completion of Phase 5) is the following:

> In the pre-geometric multivector substrate there exists an algebraic rigidity (or closure) condition on the protected bivector structure. Configurations in which the protected modes generate exactly three independent, closed directions are stable and self-reinforcing. The living candidate dynamics (especially the quadratic terms and density–chirality feedback) act to maximize protected density subject to this rigidity condition. As a result, effective 3D geometry, local isotropy (constant `c`), and preferred chirality emerge as consequences of this deeper algebraic stability requirement rather than as presupposed geometric facts.

Two complementary ways this could manifest were identified:

- **Multivector trigonometry / closure relations** — algebraic identities that are satisfied at low cost only when three protected directions close properly.
- **Extremization of protected density subject to algebraic rigidity** — the growth/attachment rule (and ultimately the dynamics) favors nodes whose local protected subalgebra is close to the rigid-at-3 configuration.

The back-propagation / evolutionary search described in this document was the first concrete attempt to let the data itself tell us which functional forms on the protected bivector structure actually improve stability and push the attractor toward higher effective dimensionality.

---

## 2. Experimental Design

### 2.1 Term Bank

We defined a bank of 8 candidate functional forms computed locally from each node’s protected bivector components and state:

1. `term_rank_deviation` — `(r - 3)^2` (core quadratic penalty on number of independent directions)
2. `term_rank_abs` — `|r - 3|` (L1 version)
3. `term_closure_defect` — soft algebraic closure defect (originally `1/vol`, stabilized to `exp(-3·vol)`)
4. `term_balance` — penalty for unbalanced bivector magnitudes
5. `term_rho_rank_interaction` — `ρ · (r-3)^2`
6. `term_rho_closure_interaction` — `ρ · closure_defect`
7. `term_biv_magnitude` — overall strength of protected bivectors
8. `term_rank_biv_interaction` — `(r-3)^2 · biv_magnitude`

The optimizer learned a scalar coefficient for each term. The attachment bias (living activity) was multiplied by `exp(-rigidity_weight · weighted_sum_of_terms)`.

### 2.2 Loss Function

A multi-objective loss was used:

- Maximize late-time mean `d_eff`
- Minimize variance of `d_eff` across nodes (stability)
- Maximize protected density concentration
- Penalize poor recovery after a soft perturbation
- Small L1 pressure on the term coefficients (encourages sparsity / selective ablation)

### 2.3 Optimization Method

- Algorithm: Differential Evolution (global, derivative-free)
- Population: 20
- Generations: 15
- Seeds per individual: 4 (to reduce sensitivity to stochastic graph growth)
- Parallelism: up to 16 workers via joblib (inside the fitness function)
- Polish: enabled (L-BFGS-B on the best individual)
- Mutation and recombination tuned for exploration

This is an intentionally intensive search designed to avoid shallow local minima and to test terms both in isolation and in combination.

---

## 3. Results

### 3.1 Best Coefficients Found

After the full search the optimizer returned the following coefficient vector (positive = term improves the loss when added to the attachment bias):

| Term                              | Coefficient | Relative Strength | Notes |
|-----------------------------------|-------------|-------------------|-------|
| `term_rank_biv_interaction`       | +2.2414     | ★★★★              | Strongest term |
| `term_rank_deviation`             | +2.1491     | ★★★★              | Classic (r-3)² very effective |
| `term_closure_defect`             | -1.8352     | ★★★               | Negative — rewards good closure |
| `term_rank_abs`                   | +1.2325     | ★★                | L1 rank penalty also helpful |
| `term_rho_closure_interaction`    | +0.7384     | ★★                | Density-weighted closure useful |
| `term_rho_rank_interaction`       | +0.6928     | ★★                | Density-weighted rank deviation useful |
| `term_balance`                    | +0.4442     | ★                 | Mild benefit |
| `term_biv_magnitude`              | -0.5153     | ★                 | Standalone magnitude slightly harmful |

Best loss achieved: **-0.6462**

### 3.2 Interpretation

- The two strongest terms both involve **rank deviation from 3**. This is direct empirical support for the central hypothesis that the system “wants” exactly three independent protected directions.
- The strongest single term is the **interaction between rank deviation and bivector magnitude**. This is a beautiful result: having the wrong number of directions hurts *more* when the protected winding is strong. This matches the physical picture that high protected chirality amplifies the cost of dimensional mismatch.
- The negative coefficient on the closure defect term indicates that the optimizer preferred to *reward* good algebraic closure rather than (or in addition to) heavily penalizing bad closure. This is consistent with a “truss rigidity” view — closed, balanced three-direction configurations are actively stabilized.
- Pure bivector magnitude by itself is slightly disfavored, but when multiplied by rank deviation it becomes one of the best terms. This is classic interaction-term behavior.

---

## 4. Immediate Next Steps (Post-Search Analysis)

The search produced a promising coefficient vector. The following actions are planned (and some already under way):

1. **Post-hoc ablation study** on the discovered coefficients (zero out the strongest terms one by one or in groups and measure degradation in the loss / final metrics). This gives a quantitative ranking of importance.

2. **Validation on the full discrete simulator**. Take the best coefficient set (and the top ablated variants) and insert them into the non-differentiable Phase 6 graph growth code. Run proper 500–600+ node experiments with the same perturbation protocol as Phase 4. Compare stability and effective dimensionality against the Phase 4 baseline.

3. **Export of best configurations** for the Lean side (`Phase6Rigidity.lean`), so we can begin formalizing the winning algebraic rigidity expression.

4. Refinement of the term bank and loss if the validation runs show clear directions for improvement.

---

## 5. Relation to the Overall Program

This back-propagation experiment is the first data-driven probe into the deeper “why 3?” question that was left open at the end of Phase 5.

If the validation runs confirm that one or more of the discovered terms meaningfully improve both the achieved `d_eff` and the robustness of the attractor, we will have the first concrete evidence that an algebraic rigidity condition on the protected bivector structure is the missing pre-geometric principle that selects three dimensions.

This would be a major conceptual step: instead of treating “3D” as an empirical observation or an imposed background, we would be able to say that three independent closed directions are the only algebraically stable configuration for high-density protected modes under the living candidate dynamics — and that our perceived 3D geometry is the macroscopic shadow of that microscopic stability condition.

---

*Document written immediately after the completion of the intensive Differential Evolution search on 2026-05-20.*