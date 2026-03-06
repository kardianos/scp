# HFKT v18: Robust Integrator and Dynamic Skyrme Tensor
**Addressing the Implementation Gaps of V17**

The rigorous review of the V17 execution revealed three critical technical flaws that prevented the codebase from genuinely validating the Harmonic Field Knot Theory (HFKT). While the structural boundary hacks and Coulomb insertions were removed, the numerics themselves were insufficient.

### The Failures of V17:
1.  **Approximate Skyrme Term:** The `_skyrme_contribution` implemented a standard static-energy discretization rather than the full dynamical Maurer-Cartan commutator bracket ($[L_\mu, L_\nu]$) required by the complete Skyrme-Faddeev equation of motion. Consequently, the potential barrier was too weak to counteract radiation in 3+1D, leading to the core radiating away by step 200 (which the raw data honestly exposed).
2.  **Unstable Integrator:** The explicit forward Euler integration method introduced artificial numerical damping and phase errors. For highly non-linear topological wave equations, a symplectic integrator (like Leapfrog or Velocity Verlet) is strictly necessary to preserve stability. Furthermore, the total energy diagnostic omitted the quartic Skyrme energy, and the grid lacked necessary absorbing boundaries.
3.  **Static Superposition Hack:** Test 3 computed interaction energy by taking a static spatial snapshot of a quaternion product ($R_1 * R_2$). Because it failed to dynamically evolve the interacting solitons through the PDE, it proved nothing about true emergent Newtonian gravity over time, falling back into methodological circularity.

### The V18 Pivot
V18 will close these three final implementation gaps by deploying a mathematically exact, dynamically stable execution suite:
1.  **Full Rotor-Commutator Skyrme Term:** Implementing the precise dynamical Maurer-Cartan commutator expansion.
2.  **Symplectic Leapfrog Integrator:** Advancing the field with a high-order, energy-preserving scheme, equipped with exact topological Hopf-charge ($Q$) conservation diagnostics and proper absorbing boundary layers.
3.  **Dynamic Soliton Interaction:** Testing gravitational and geometrical forces by freely evolving two distinct Hopfions through the rigorous PDE across time, permanently retiring static snapshot evaluations.
