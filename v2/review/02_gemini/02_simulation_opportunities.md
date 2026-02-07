# CHPT Simulation Plan: From Math to Code

To validate the "loose" conceptual areas and demonstrate the theory, we propose a tiered simulation strategy.

## Tier 1: The Scalar Soliton (Proof of Stability)
**Goal**: Demonstrate that the nonlinear potential $V(\Psi) = (|\Psi|^2 - \rho_0^2)^2$ creates stable localized structures (solitons) in lower dimensions.

### Simulation 1A: 1D Kink
*   **Equation**: $\partial_t^2 \phi - \partial_x^2 \phi + \lambda \phi (\phi^2 - 1) = 0$
*   **Domain**: 1D Line.
*   **Expected Result**: Stable "Kinks" (transitions from -1 to +1) that behave like particles (bounce off each other).
*   **Tech**: Python (`numpy`, `matplotlib` animation).

### Simulation 1B: 2D Vortex
*   **Equation**: Complex scalar field $\psi$ in 2D.
*   **Domain**: 2D Grid.
*   **Topological Charge**: Winding number around a point (Vortex).
*   **Expected Result**: Vortices that are stable and interact (repel/attract).
*   **Rellevance**: This is the 2D analog of the 3D Knot.

## Tier 2: The Vector Field (Maxwell Integration)
**Goal**: Verify that the Bivector component $\mathbf{F}$ propagates as a wave at speed $c$.

### Simulation 2: 3D Wave Propagation
*   **Equation**: Linearized multivector wave equation $\square \Psi = 0$.
*   **Init**: Gaussian pulse in $\mathbf{F}$ component.
*   **Measurable**: Speed of wavefront, polarization preservation.
*   **Tech**: C++ / Python (`scipy.sparse` or FDTD mesh).

## Tier 3: The Hopfion (The Golden Test)
**Goal**: Simulate a true 3D Knot ($Q=1$) in the full nonlinear field. Current math predicts it *should* be stable.

### Simulation 3: 3D Nonlinear Evolution
*   **Equation**: Full $Cl(3,0,1)$ equation of motion.
*   **Init**: Construct a "Hopfion" ansatz (a field configuration with Hopf index 1).
    *   Map $\mathbb{R}^3 \to S^3$.
    *   Initialize $\Psi(x, 0)$ with this map.
*   **Evolution**: Run the PDE solver.
*   **Success Criteria**:
    1.  Does the energy density remain localized? (Particle stability)
    2.  Does it disperse? (Theory Failure)
    3.  Does it collapse? (Derrick's Theorem failure)

## Recommended Stack
1.  **Python Prototype**: Rapid testing of 1D/2D cases.
2.  **C++ / Compute Shader**: For 3D $256^3$ grids, we need GPU acceleration (WebGPU or CUDA).

## Immediate Action Item
Write a Python script to simulate **Simulation 1B (2D Vortex)**.
*   It visualizes "Phase" (Topology).
*   It demonstrates the mass/energy localization mechanism.
*   It uses the *exact* potential term from `spec/math/03_dynamics.md`.
