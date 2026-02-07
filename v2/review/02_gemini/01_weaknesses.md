# CHPT Phase 2 Review: Weaknesses and Gaps

This document identifies the "loose" areas in the CHPT specification where the mathematical formalization is either incomplete, hand-wavy, or in tension with known physics.

## 1. The Fractional Charge Problem (Quarks)

**The Issue**:
Current CHPT defines Electric Charge ($Q$) as the **Hopf Index** or **Winding Number** of the field map $\Psi: S^3 \to S^3$.
By definition, homotopy groups are integers: $\pi_3(S^3) = \mathbb{Z}$.
Therefore, $Q$ must be an integer.
**Observation**: Quarks have $Q = \pm 1/3, \pm 2/3$.

**Why it's loose**:
We have no topological mechanism for "1/3 winding". A map wraps around the sphere a whole number of times, or it doesn't.
*   **Proposed Fix**: Redefine the "fundamental" charge unit. If the Electron ($Q_{obs} = -1$) corresponds to a winding of $Q_{top} = 3$, then a Quark could have $Q_{top} = 1$.
*   **Risk**: This contradicts the intuitive identification of the simplest knot (Hopfion, $Q=1$) with the simplest charged particle (Electron). If the Electron is $Q=3$, what are the $Q=1$ and $Q=2$ particles? Are they hidden? Unstable?

## 2. The Mass Spectrum Prediction

**The Issue**:
We defined Mass as the eigenvalue of the Hamiltonian: $M = \int \mathcal{H} d^3x$.
We *claimed* this produces a discrete spectrum (Electron, Muon, Proton...).
**Observation**: We have not calculated a single spectral ratio.

**Why it's loose**:
We assume that the "fundamental knot" has the mass of an electron/proton, but we have no proof that the equation supports stable solutions with valid mass ratios.
*   Why is $M_{proton} / M_{electron} \approx 1836$?
*   Does the theory predict a stable particle with mass 1500x electron? Or 2.5x?
*   Without finding the actual soliton solutions, this is just a wish, not a prediction.

## 3. The Weak Force Mechanism

**The Issue**:
We attribute the Weak Force to "Chirality" and "Parity Violation" in `spec/05_chirality.md`.
**Observation**: We haven't defined the term in the Lagrangian that breaks Parity.

**Why it's loose**:
The current Lagrangian $\mathcal{L} = \langle (\nabla \Psi)^2 \rangle - V(\Psi)$ is likely Parity-Invariant (depending on exactly how $\Psi$ transforms).
To get the Weak Force (which ignores Right-handed matter), we need:
1.  **Explicit Symmetry Breaking**: A term like $\Psi \langle \dots \rangle_4$ in the action?
2.  **Spontaneous Symmetry Breaking**: A definition of $\rho_0$ that selects a specific chiral frame.
Currently, the mechanism is "hand-waved" as "Chirality $\approx$ Weak". We lack the $SU(2)_L$ gauge symmetry structure entirely.

## 4. Spin-Statistics Theorem

**The Issue**:
We asserted that "Knots (Fermions) repel" and "Waves (Bosons) pass through".
**Observation**: In QFT, the Spin-Statistics theorem is a deep result of Lorentz Invariance and path integrals.

**Why it's loose**:
In CHPT, we are imposing this behavior manually via "geometric incompatibility".
We need to prove that a $Q=1$ knot *actually* requires an antisymmetric wavefunction when two are exchanged. Currently, there is no multi-particle wavefunction in the classical field description, so "statistics" is ill-defined until we verify the semi-classical limit.

## 5. Gravity Quantization

**The Issue**:
Gravity is defined as an "effective metric" arising from density gradients.
**Observation**: This is a classical description.

**Why it's loose**:
Does this theory have Gravitons?
*   If we linearize the density perturbation, do we get a spin-2 wave?
*   `spec/math/04_electromagnetism.md` showed we get spin-1 (Vector/Bivector) waves.
*   Where is the spin-2 component? If gravity is purely scalar (density gradient), it conflicts with GR (which requires spin-2 for light bending and quadrupole radiation).
*   **Critical Risk**: Scalar Gravity typically gives the wrong light bending prediction (0 or half of GR).

## Implementation Strategy
To address these, we cannot just "write more spec". We need to **Simulate**.
Numerical solutions of the field equation will tell us if $Q=1$ is stable, what its mass is, and if it repels another $Q=1$ knot.
