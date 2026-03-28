# η₁ First-Principles Derivation

## Goal
Determine whether η₁ in η(P) = η₀ + η₁|P| is uniquely fixed by the
existing parameters (m², μ, κ, η₀) through a self-consistency condition,
or whether it is a genuinely new free parameter.

## Approach

### Step 1: Symbolic virial self-consistency (Maxima)
Derive the virial identity for the modified Lagrangian with η(P).
Express η₁ in terms of profile integrals. Show the functional form
of the self-consistency equation.

File: step1_virial.mac → step1_virial_output.txt

### Step 2: Dimensional/naturalness analysis (analytical)
Enumerate all possible combinations of (m², μ, κ, η₀) that have the
correct dimensions for η₁. Check which ones give the right magnitude.

File: step2_naturalness.md

### Step 3: Numerical profile integrals (Python + V45 SFA)
Compute the required integrals from actual simulation data:
- ∫|P|² d³x (binding density)
- ∫|∇×φ|² d³x (curl energy)
- ∫|P|×|φ·(∇×θ)| d³x (topology-weighted coupling)
- ∫|∇φ|² d³x (gradient energy)

Evaluate the self-consistency equation numerically.

File: step3_integrals.py → step3_integrals_output.txt

### Step 4: Cross-check the relationship η₁ = |μ|/η₀
Verify whether the numerical result matches the conjectured relationship.

File: step4_crosscheck.md

## Deliverables
- eta_derivation/RESULTS.md — consolidated findings
- All intermediate files preserved for reproducibility
