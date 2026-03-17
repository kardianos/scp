# V23-A: Hessian Splitting Analysis

## Thesis

The triple-product potential V(P) = (μ/2)P²/(1+κP²) with P = φ₁φ₂φ₃ has a
Hessian (second derivative matrix in field space) whose eigenvalues DIFFER for
symmetric vs antisymmetric perturbations around the oscillon background. If the
antisymmetric ("shear-like") eigenvalue is more negative than the symmetric
("compression-like") eigenvalue, then the shear sector has a smaller effective
mass gap inside the oscillon core. This is the prerequisite for all downstream
gravity proposals.

**This is the foundational calculation. If no splitting exists, proposals 2+7
cannot work. If splitting exists, it determines the parameters for crossgrad
and critical investigations.**

## Mathematical Setup

### Background

The symmetric oscillon has φ₁ = φ₂ = φ₃ = f(r)cos(ωt), where f(r) is the
radial profile from v21 (data in `/home/d/code/scp/v21/data/`).

### The Hessian

The potential is V = (μ/2)P²/(1+κP²), P = φ₁φ₂φ₃.

Second derivatives:

    ∂²V/∂φ_a∂φ_b = μ [(∂P/∂φ_a)(∂P/∂φ_b)(1 - 3κP²) + P(∂²P/∂φ_a∂φ_b)(1+κP²)] / (1+κP²)³

where:
    ∂P/∂φ₁ = φ₂φ₃,  ∂P/∂φ₂ = φ₁φ₃,  ∂P/∂φ₃ = φ₁φ₂
    ∂²P/∂φ₁∂φ₂ = φ₃,  ∂²P/∂φ₁² = 0,  etc.

On the symmetric background φ_a = f:

    P = f³
    ∂P/∂φ_a = f²  (all equal)
    ∂²P/∂φ_a∂φ_b = f  (for a≠b),  0  (for a=b)

So the Hessian matrix at the symmetric point is:

    H_{ab} = μf⁴(1-3κf⁶)/(1+κf⁶)³  for a=b    (diagonal)
    H_{ab} = μf⁴(1-3κf⁶)/(1+κf⁶)³ + μf⁴(1+κf⁶)/(1+κf⁶)³  for a≠b

Wait — let me redo this carefully. With φ_a = f for all a:

    (∂P/∂φ_a)(∂P/∂φ_b) = f⁴  for all a,b
    P(∂²P/∂φ_a∂φ_b) = f³·f = f⁴  for a≠b
    P(∂²P/∂φ_a∂φ_a) = f³·0 = 0   for a=b

So:
    H_{aa} = μf⁴(1-3κf⁶) / (1+κf⁶)³
    H_{ab} = μ[f⁴(1-3κf⁶) + f⁴(1+κf⁶)] / (1+κf⁶)³   for a≠b
           = μf⁴(2-2κf⁶) / (1+κf⁶)³
           = 2μf⁴(1-κf⁶) / (1+κf⁶)³

### Eigenmodes

The Hessian is a 3×3 matrix with diagonal d and off-diagonal o:

    H = | d  o  o |
        | o  d  o |
        | o  o  d |

Eigenvalues:
    λ_sym  = d + 2o   (eigenvector: (1,1,1)/√3 — symmetric/compression)
    λ_anti = d - o     (degenerate, multiplicity 2 — antisymmetric/shear)

Substituting:
    d = μf⁴(1-3κf⁶) / (1+κf⁶)³
    o = 2μf⁴(1-κf⁶) / (1+κf⁶)³

    λ_sym  = μf⁴[(1-3κf⁶) + 4(1-κf⁶)] / (1+κf⁶)³
           = μf⁴(5 - 7κf⁶) / (1+κf⁶)³

    λ_anti = μf⁴[(1-3κf⁶) - 2(1-κf⁶)] / (1+κf⁶)³
           = μf⁴(-1 - κf⁶) / (1+κf⁶)³
           = -μf⁴(1 + κf⁶) / (1+κf⁶)³
           = -μf⁴ / (1+κf⁶)²

### The Splitting

The effective mass² for each sector:

    m²_sym(r)  = m² + λ_sym(r)  = m² + μf⁴(5-7κf⁶)/(1+κf⁶)³
    m²_anti(r) = m² + λ_anti(r) = m² - μf⁴/(1+κf⁶)²

For μ < 0 (attractive):
    λ_sym  = |μ|f⁴(7κf⁶-5)/(1+κf⁶)³   — sign depends on whether 7κf⁶ > 5
    λ_anti = |μ|f⁴/(1+κf⁶)²            — ALWAYS POSITIVE (reduces m²_anti)

**Key result**: With μ < 0, the antisymmetric effective mass is ALWAYS reduced:

    m²_anti(r) = m² - |μ|f⁴/(1+κf⁶)²

This is LESS than m² wherever f(r) > 0, i.e., inside the oscillon core.

The maximum reduction occurs at the oscillon center where f = f(0) = A (peak
amplitude). At the v21 production parameters (μ=-20, κ=20, m=1.0, A≈0.89):

    m²_anti(0) = 1.0 - 20·(0.89)⁴/(1+20·(0.89)⁶)²

This must be computed numerically.

**The splitting Δm² = m²_sym - m²_anti is NONZERO even at η=0** (no cross-
gradient needed for the splitting to exist). The cross-gradient amplifies it.

## What to Compute

1. **Analytic Hessian eigenvalues** λ_sym(r) and λ_anti(r) as functions of
   the oscillon profile f(r). Use the v21 1D profile data.

2. **Effective mass profiles** m²_sym(r) and m²_anti(r). Plot vs r.

3. **Check for sign change**: Does m²_anti(r) go negative anywhere? If so,
   the antisymmetric mode is TACHYONIC (unstable) — important to know.

4. **Time-averaged Hessian**: The oscillon oscillates as f(r)cos(ωt). The
   time-averaged Hessian involves ⟨f⁴cos⁴⟩ = (3/8)f⁴. Compute the
   time-averaged effective masses.

5. **Perturbation eigenvalue problem**: Solve the Sturm-Liouville problem
   for the antisymmetric mode in the potential well created by the oscillon:
     -u'' + [m²_anti(r) - ω²]u = 0
   Find bound states (if any) and the scattering length.

6. **Correlation length** ξ of the antisymmetric mode: how far does the
   shear-softened region extend beyond the oscillon core?

## Reference Code

- v21 1D oscillon solver: `/home/d/code/scp/v21/src/triad1d.c`
- v21 3D oscillon solver: `/home/d/code/scp/v21/src/triad3d.c`
- v21 production profile data: `/home/d/code/scp/v21/data/`
- v12 normal modes solver: `/home/d/code/scp/v2/src/normal_modes.c`
  (Sturm-Liouville eigenvalue problem — similar structure needed here)
- v12 bound state solver: `/home/d/code/scp/v2/src/bound.c`
  (eigenvalue search in effective potential — directly adaptable)

## Output Structure

- `src/hessian.c` — main computation code
- `data/` — output data files
- `RESULTS.md` — results and analysis

## Parameters

Use v21 production values: μ=-20, κ=20, m=1.0
The oscillon profile f(r) can be loaded from v21 data or recomputed from
the 1D solver.

Compile: `gcc -O3 -Wall -o hessian src/hessian.c -lm`
