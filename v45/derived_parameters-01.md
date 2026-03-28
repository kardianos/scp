# Derived Parameters for SCP Field Theory
## First Principles
The SCP medium is a single 3-component vector field φ with dynamics governed by a relativistic wave equation plus a nonlinear potential V(P) where P = φ₀φ₁φ₂ is the local volume form. All observed structures (braids, protons, binding) are attractors of this equation. Any "free" parameters must ultimately be fixed by internal consistency conditions of the field itself rather than external tuning.
Empirically tuned values (m²=2.25, μ=-41.345, κ=50, η=0.5, A_bg=0.1, δ={0,3.0005,4.4325}) produce stable braids and gravity but failed to produce nuclear binding energy in V45. This indicates the parameters are not yet in the correct relational web.
## Why Parameters Must Be Derived
Free parameters break predictive power and hide missing physics. A mature theory should have:
- Only dimensionless ratios or quantities fixed by symmetries/conservation laws.
- All dimensional parameters emerging from a single scale (e.g. background oscillation frequency).
- Self-tuning mechanisms so that stable structures automatically select optimal values.
V45 showed energy is minimized at infinite separation. Binding requires a deeper minimum at finite D, which can only appear if the potential and coupling terms are related in a specific way that current independent tuning misses.
## Proposed Relations
1. **Mass term m²**  
   Fixed by background stability: the uniform oscillating background φ_a = A_bg cos(kz + 2πa/3) must be linearly stable. This sets m² ≈ k_bg² where k_bg = π/L (simulation scale). Current m²=2.25 is close but should be derived as m² = ω_bg² - k_bg² with ω_bg chosen for minimal radiation.
2. **Potential parameters μ and κ**  
   The potential minimum occurs at P_opt = 1/√(3κ). Set κ so P_opt ≈ 0.082 matches measured core triple-product (0.11-0.15). Then derive μ from equilibrium binding density:  
   μ = - (m² P_core) / (dV/dP at P_core)  
   This ties μ directly to κ and m², eliminating one free parameter. Binding energy scale then follows automatically.
3. **Coupling η**  
   η controls θ/φ energy ratio (measured 28% at far field for η=0.5). Set η so far-field θ/φ = α_phys ≈ 1/137 scaled to code units, or from force balance condition that strong:EM → 1:1 exactly at nuclear equilibrium distance (as partially seen in V42). Make η locally dependent on |∇×φ| for self-tuning.
4. **Background amplitude A_bg and phases δ**  
   A_bg fixed by requiring total energy density ρ_bg ≈ 0.03 and braid depletion depth matching gravitational strength. Phases δ derived from minimizing destructive interference in 3-braid composites (exact 2π/3 multiples preferred over empirical 3.0005/4.4325).
## Dynamic Parameters (Stronger Proposal)
Instead of constants, promote key parameters to fields or functions:
- κ(P) or μ(P) so potential stiffens exactly where needed for binding.
- η = f(|curl φ|, local breathing phase) so EM coupling activates only at correct separation.
- A_bg obeys a slow relaxation equation enforcing global average P=const, creating self-consistent depletion wells that naturally bind composites.
This turns the "compression then release while spinning" idea into a global constraint: initialize with high uniform compression + net angular momentum; the relaxation dynamics automatically select all scales and produce bound states as conserved quantities force the system into lower-energy configurations.
## Next Experiments (v46)
- Implement derived relations in a new config generator.
- Test compressed-spin seed (high A_bg sphere with rotation) on N=384 grid.
- Measure resulting binding wells vs V45 baseline.
- Update CONCEPT.md to replace tuned parameters with these relations.
These changes restore predictivity: the theory should then contain no free parameters beyond one overall scale.