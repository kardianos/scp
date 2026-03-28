# η₁ Derivation Results

## Summary

η₁ in η(P) = η₀ + η₁|P| is NOT a free parameter — it is constrained
by the force balance condition at the braid core. The derivation gives
η₁ = (1/η₀ - η₀)/A_core³ where A_core is the braid core amplitude,
itself determined by (m², μ, κ).

## Methods and Results

### Step 1: Virial self-consistency (Maxima)
- Derived the virial identity: E₁ + 2E₂ + 3E₃ = 0
- Expressed η₁ = -R/(2η₀J₁) where R is the virial residual
- Derived alternative: force balance at core gives η₁ = (1/η₀-η₀)/A³

**Output**: step1_virial.mac → step1_virial_output.txt

### Step 2: Naturalness analysis
- Surveyed all parameter combinations
- Leading candidate: η₁ = |μ|/η₀ = 82.7 gives η_eff(core) = 8.8
- Physical interpretation: amplification proportional to binding
  strength divided by EM coupling

**Output**: step2_naturalness.md

### Step 3: Numerical integrals (Python + V45 data)
- Computed virial residual from global diagnostics
- **Problem**: background subtraction error dominates (136k ± 200)
- Virial approach inconclusive from global data — needs per-voxel SFA
- Force balance: η₁ = 35-200 depending on A_core (0.2-0.35)

**Output**: step3_integrals.py → step3_integrals_output.txt

### Step 4: Cross-check
- Conjecture η₁ = |μ|/η₀ = 83 is in the correct range (35-200)
- NOT rigorously derived from virial (background subtraction problem)
- Force balance IS a first-principles result
- A_core is determined by (m², μ, κ) → η₁ is determined
- η₁ is constrained, not free, but the exact value needs the
  self-consistent radial ODE solution

**Output**: step4_crosscheck.md

## Key Findings

| Finding | Status |
|---------|--------|
| η₁ is constrained (not free) | **ESTABLISHED** — force balance determines it |
| η₁ range | **35 - 200** (A_core = 0.35 - 0.20) |
| η₁ = \|μ\|/η₀ = 83 conjecture | **PLAUSIBLE** but not rigorously derived |
| Virial derivation | **INCONCLUSIVE** — needs per-voxel SFA computation |
| Force balance derivation | **ROBUST** — independent of background |

## Per-Voxel Virial Computation (Step 5)

The per-voxel approach was implemented (virial_soliton.c) and run on
two independent datasets:
- V43 proton template (64³, pre-converged)
- V41 UUD proton (192³, T=200)

### Method
1. Compute P, ∇φ, ∇×φ, ∇×θ at each voxel from SFA data
2. Integrate energy densities only over soliton region (|P| > threshold)
3. Sweep threshold to isolate core contribution
4. Differential shell analysis to cancel background contamination

### Result
**η₁ ≈ 90-120** from two independent methods:
- Core shell virial (template, |P|=0.05-0.10): **η₁ = 118**
- V41 extrapolation to P_opt: **η₁ = 119**
- Force balance at A_core=0.25: **η₁ = 96**

### Approximate relationship
**η₁ ≈ √2 × |μ| / η₀ ≈ 117**

This connects η₁ to existing parameters with no new free constants.

## Recommended Value for Simulation

**η₁ = 115** (central value from virial analysis)

Gives:
- Far field: η = 0.5 (EM unchanged)
- Braid core (|P|=0.1): η = 12 (strong nuclear)
- Sufficient for binding (V45 threshold: η ≈ 8)

## Files

- step1_virial.mac — Maxima symbolic derivation
- step1_virial_output.txt — Maxima output
- step2_naturalness.md — Parameter combination analysis
- step3_integrals.py — Numerical computation
- step3_integrals_output.txt — Numerical results
- step4_crosscheck.md — Analysis and conclusions
- RESULTS.md — this file
