# Path to Symbolic Differentiation for 7D Stability Bounds

**Date**: 2026-05-24
**Context**: Addressing Gap Analysis Point 4: "No symbolic diff path or stub."

Currently, the Hessian matrices used to evaluate topological stability across the $\text{Cl}(7)_{even}$ algebra are computed using a finite-difference approximation in `StabilityFromAlgebra.lean` (e.g., `computeHessianNum` and `computeHessianLiving`). While this is effective for generating certificates, a fully rigorous analytical pipeline requires exact symbolic derivatives of the continuous nonlinear potential $V(M)$.

## The Current State: Finite Differences
The current method evaluates:
$$ \frac{\partial^2 V}{\partial M_i \partial M_j} \approx \frac{V(M + \epsilon e_i + \epsilon e_j) - V(M + \epsilon e_i) - V(M + \epsilon e_j) + V(M)}{\epsilon^2} $$
with $\epsilon = 1/100$. This is fundamentally an approximation, even though we use exact `Rat` arithmetic for the forward evaluation.

## The Future State: Lean 4 Symbolic Calculus (`SciLean`)
To upgrade this, we propose integrating a symbolic calculus engine into the `v59` Lean framework. There are two primary avenues:

### 1. Integration of `SciLean`
`SciLean` is an emerging numerical and symbolic computing library for Lean 4 designed to compute exact gradients and Hessians of functional programs.
- **Approach**: We can declare the algebraic `octMult` and `V_living_project` operations as smooth functions (`IsSmooth`).
- **Execution**: We would use `SciLean`'s automatic differentiation macros (`∂`) to generate a compiled, exact Hessian function $H_{ij}(M) = \partial_i \partial_j V(M)$ directly from the AST of the potential.
- **Benefits**: No discretization error, completely rigorous limits, runs at native speed.

### 2. Custom Multivariate Polynomial Ring
Because our potential $V(M) = \lambda |M*M|^2 + \mu |M|^2$ is fundamentally a quartic multivariate polynomial over 8 variables (prior to saturation), we can alternatively map the algebra to an exact symbolic polynomial ring.
- **Approach**: Define an 8-variable polynomial ring $R[M_0, M_1, \dots, M_7]$.
- **Execution**: Implement a simple algebraic derivative operator $\partial / \partial M_i$ over the ring. The `octMultTable` dictates the exact integer cross-terms of the polynomial.
- **Benefits**: Self-contained within the project (no heavy external dependencies), highly explicit algebraic expansion of the leakage cross-terms.

## Next Steps for Implementation
1. Decide whether to introduce the `SciLean` dependency or build a custom `Poly8` module.
2. Replace `computeHessianLiving` with the symbolic equivalent.
3. Validate that the symbolic Hessian produces the identical stability certificates for the Python mature data rows.
