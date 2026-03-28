# Derived Parameters for SCP Field Theory — Revision 03

## 1. The Critique: Why Previous Derivations Failed

The `v45` null result for nuclear binding is the most important "failure" in the project's history. It exposes a fundamental disconnect in the previous parameter strategies (`-01` and `-02`):

### The Virial Mismatch (The "Exploding Proton" Problem)
In `-02`, the virial theorem was proposed as a consistency check. Applying the virial scaling to our potential $V(P)$ reveals the crisis. For a stable 3D soliton in this theory, the energy components must satisfy:
$$E_{grad} + E_{mass} + 6 \int \frac{V(P)}{1 + \kappa P^2} dV \approx 0$$
Using the `v45` diagnostic data ($E_{grad} \approx 65,000$, $E_{mass} \approx 42,000$, $E_{pot} \approx -3,200$):
$$107,000 - \text{(approx 10,000)} = 97,000 \neq 0$$
**Conclusion:** The current particles are not in equilibrium; they are "exploding solitons" held together only by the carrier-wave phase-lock and the limited simulation time ($T=500$). They have so much internal pressure that they cannot possibly bind to another particle; they are effectively "hard-sphere" repulsive because their internal gradient energy dwarfs their binding potential by 20:1.

### The Phase Incoherence Problem
Previous derivations treated parameters as static constants. But in a 3-field system, the interaction energy $V(P_1 + P_2)$ depends entirely on the **relative phase** of the carrier waves. If two baryons are placed with random or "best-guess" phases, they will almost always interfere destructively in $P = \phi_0\phi_1\phi_2$, leading to repulsion ($P \to 0$).

---

## 2. The New Path: The Topological Quench Model

Instead of "placing" particles, we must **precipitate** them from a global high-energy state. This ensures they are phase-locked to the same vacuum and satisfy the virial theorem from birth.

### First Principle 1: The Machian Energy Constraint
Total energy in a closed volume must be conserved. As matter forms, it "borrows" energy from the background.
**Proposal:** Promote $A_{bg}$ from a constant to a global Lagrange multiplier:
$$A_{bg}(t) = \sqrt{\frac{E_{total} - E_{particles}(t)}{Volume \cdot \omega_{bg}^2}}$$
This creates an **automatic attractive force**: particles that share a depletion well (bind) "release" background energy back to the system, which is a lower-energy state for the matter-field composite.

### First Principle 2: The Coriolis Stabilization (Spin)
The UDD neutron fails because its color phases converge. We need a force to keep the "colors" (field components) rotating.
**Proposal:** Introduce a global angular momentum density $\Omega$. In the field equations, the time derivative transforms:
$$\partial_t \to \partial_t + \Omega \times$$
This adds a Coriolis-like term $(\Omega \cdot \nabla)\phi$ that prevents phase convergence. The "Neutron" becomes a stabilized **vortex-cluster** where the rotation prevents the topological collapse.

### First Principle 3: The Adler Phase-Lock
Binding is **Synchronization**. Two baryons bind when their internal breathing modes lock into phase.
The effective potential between two baryons at separation $D$ is:
$$E_{int}(D) = \mu \int P_1(\mathbf{x}, \psi_1) P_2(\mathbf{x}, \psi_2) dV$$
where $\psi$ is the breathing phase. We derive $\mu$ by requiring that the synchronization energy (Adler locking) is greater than the thermal noise:
$$\mu \geq \frac{k_B T_{kin}}{\Delta P_{breathing}}$$
This ties $\mu$ directly to the particle's internal dynamics.

---

## 3. The v46 Implementation Plan: "Vortex Nucleosynthesis"

We will move away from "seed files" and toward **"Quench Scripts."**

### Step 1: The Compressed Spin Seed
Initialize a large $N=384$ or $512$ grid with:
1.  **High Compression:** $A_{bg} = 0.4$ (4x the standard density).
2.  **Global Spin:** A velocity field $v = \Omega \times r$ imprinted on the $\phi$ carrier waves.
3.  **Thermal Noise:** $A_{noise} = 0.05$ to break symmetry and seed fragmentation.

### Step 2: The Annealing Relaxation
Run the simulation with a temporary "damping" term $\gamma \partial_t \phi$ to bleed off excess kinetic energy. 
*   **$T=0 \to 100$:** High damping. The field "curdles" into high-$|P|$ regions.
*   **$T=100 \to 200$:** Reduce damping. The "curds" organize into vortices.
*   **$T=200 \to 1000$:** Zero damping. Observe whether the resulting baryons are naturally bound.

### Step 3: Parameter Correction
Based on the virial analysis, we must use:
*   **$m^2 = 1.0$** (Lowering confinement to allow longer-range interaction)
*   **$\mu = -250.0$** (Increasing potential strength 6x to balance the gradient energy)
*   **$\kappa = 20.0$** (Shifting the saturation point to allow larger $|P|$ cores)

## 4. Expected Breakthrough
By letting the system "find its own parameters" through the quench, we expect to see:
1.  **Phase-Locked Deuterium:** Baryons that emerge already meshed.
2.  **The Meson Bridge:** A persistent high-$|P|$ filament between bound baryons.
3.  **Mass Defect:** A genuine drop in $E_{total}$ compared to the "placed" baseline, because the quench-formed bound state is a true local minimum of the field's action.

This is the transition from **Kinematics** (moving objects) to **Dynamics** (the field becoming the objects).

---

## 5. Brainstorm: Speculative & Divergent Paths

In formulating Revision 03, several more radical or unlikely ideas were considered. They are recorded here to ensure no "lost paths" are forgotten.

### A. The "Liquid Space" Model (Entropic Binding)
*   **The Idea:** Space is a dense "plasma" of microscopic braids. Binding is not a potential well, but an entropic force: clusters of baryons minimize the "surface area" of the depletion zone, maximizing the volume of the undisturbed background.
*   **Status:** *Parked.* Requires high-density simulations that are currently beyond our grid resolution ($N=1024+$ needed to see the statistical "pressure" effects).

### B. Inverse Braid Theory (Theta-Matter)
*   **The Idea:** Matter is formed in the $\theta$ (angle) fields, and the $\phi$ fields are the "stiffness" of space. This would make matter fundamentally "rotational" rather than "displacive."
*   **Status:** *Rejected.* Contradicts the V33-V42 evidence that $\phi$-depletion is the primary carrier of mass and gravity. Reversing this would require a total rebuild of the theory.

### C. Non-Local Integral Potentials
*   **The Idea:** Replace $V(P)$ with a global constraint $V(\int P dV)$ to enforce absolute volume conservation across the entire universe.
*   **Status:** *Rejected.* Violates relativistic causality (non-local). While it would solve the "exploding soliton" problem instantly, it is physically unmotivated in a local field theory.

### D. The Lord Kelvin Knot (Closed Topology)
*   **The Idea:** A baryon is not a cluster of 3 open braids, but a single, complex **closed knot** (e.g., a Trefoil or Hopf link).
*   **Status:** *Investigating.* We have not yet found a stable, self-knotting solution that doesn't "un-knot" or collapse within $T=100$. Open-ended braids are the current numerical attractors.

### E. Discrete Phase-Space (Digitized Vacuum)
*   **The Idea:** Space itself is a discrete lattice of field values, and the "equation of motion" is actually a Cellular Automaton (CA) rule.
*   **Status:** *Unlikely.* Our PDE solver already approximates this, and moving to a pure CA loses the Lorentz invariance that is a cornerstone of the project.

### F. Multi-Scalar Determinant Potential
*   **The Idea:** Instead of $P = \phi_0\phi_1\phi_2$, use $P = \det(\partial_i \phi_j)$. This makes the potential a direct measure of the Jacobian of the field mapping.
*   **Status:** *High Interest but Complex.* This is mathematically much "cleaner" (Skyrme-like) but requires computing $3 \times 3$ determinants at every grid point, increasing compute cost by 4x. We select against it for `v46` to keep iteration speed high.

---

## Summary of Selection Logic

The **Topological Quench (Revision 03)** was selected as the primary path because it is the only one that addresses the **Phase-Incoherence Problem** (the root of the `v45` failure) without introducing new non-local physics or discarding the validated 6-field structure. It relies on **Dynamics** rather than **Seed Tuning**, which is the hallmark of a mature field theory.
