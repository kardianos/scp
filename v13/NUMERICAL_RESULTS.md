# HFKT v13: Numerical Simulation of 3-Knot Condensation

## Overview
This document logs the results of the `v13/simulate_3knot.py` prototype, which models the high-energy "Big Bang" injection and subsequent frictionless evolution of topological structures in $Cl^+(3,0,1)$.

## Initialization
*   Grid: 16x16x16 (periodic boundary)
*   Integrator: Leapfrog finite-difference with explicit non-linear $S^3$ projection
*   Injection: A high-energy randomized "hot soup" core with $E_{scale} = 50.0$ injected into the vacuum.

## Evolution Profile (500 steps, $dt = 0.05$)
The evolution demonstrates a classic rapid condensation from a chaotic high-energy state into a long-lived, low-energy steady resonance.

| Step | Peak Energy Density | Shear Amplitude Proxy | Phase |
| :--- | :--- | :--- | :--- |
| **0000** | 22265.52 | 44.061 | Injection (Hot Soup) |
| **0100** | 1.47 | 0.598 | Rapid Cooling / Expansion |
| **0200** | 0.17 | 0.582 | Condensation |
| **0300** | 0.10 | 0.612 | Stable Resonance |
| **0400** | 0.13 | 0.592 | Stable Resonance |

## Analysis & Physical Implications

### 1. Spontaneous Condensation vs. Frictionless Dispersal
In a linear wave equation without friction, an injected spike of energy simply radiates outward at $c$, dropping the core energy to zero. 

However, because the $Cl^+(3,0,1)$ field is topologically constrained to the $S^3$ manifold (the rotor magnitude $|R|^2 = 1$ is invariant), the acceleration paths are highly non-linear. The results show that the field does *not* disperse to zero. Instead, the non-linear self-interference of the topology traps a residual fraction of the energy.

### 2. The Pseudo-Stable Resonance Phase (Diagnostic Test)
By step 200, the intense initial explosion stabilized. The Peak Energy plateaued around $\sim 0.10 - 0.17$, while the Shear Proxy locked into a steady oscillation near $\sim 0.60$.

*Initially*, this implied that the field "froze out" into a stable topological configuration. 

However, to address the *Avenue 4 Discrepancy* (the failure of explicit 3-quark combinations), we ran a dedicated topological diagnostic (`v13/verify_soup_topology.py`) which counted the number of explicit geometric twists ($R_w < -0.5$).
> **Raw Data Artifact:** The continuous count of topological cores vs. excited volume is streamed to `v13/results/soup_topology_diagnostic.tsv`.

**Diagnostic Result:**
*   Step 0: 12 distinct localized topological cores.
*   Step 200: 1 topological core remaining.
*   Step 300+: **0 topological cores remaining.**

The "steady state" observed in the Hot Soup was **not a composite particle**. The topology completely dissolved ($B=0$). The remaining "vibration" is just a residual, topologically trivial (scalar) wave sloshing around the periodic boundary conditions like water in a bathtub. 

This brings the Avenue 4 ("Hot Soup") result into perfect alignment with the Avenue 4 ("Explicit UUD") result: **Under pure $Cl^+(3,0,1)$ kinematics, topological knots always decay/dissolve.**

---

## Part II: Explicit UUD (Proton) Entanglement Dynamics
To track the **mutual knot wiggling** predicted algebraically, we bypassed the "Big Bang" and explicitly initialized three distinct fractional knots in an Isosceles triangle ($U-U-D$). 

*   $U_1, U_2$: Isospin $+1/2$, Phase offset separated by $120^\circ$
*   $D_1$: Isospin $-1/2$, Phase offset separated by $240^\circ$

**Evolution Profile (`sim_3knot_vibration.py`):**
> **Raw Data Artifact:** The continuous pair-wise distances ($D_{long}$, $D_{short1}$, $D_{short2}$) between the three knot centers are streamed to `v13/results/sim_3knot_vibration_raw.tsv`

| Step | Pair-Wise Distances ($d_1, d_2, d_3$) | Structural Phase |
| :--- | :--- | :--- |
| **0000** | 5.00, 6.32, 6.71 | Proper Topological Init (Isosceles) |
| **0200** | 4.00, 6.32, 6.32 | Geometric Condensation (Tension) |
| **0400** | 3.61, 3.61, 7.21 | Prolonged Transverse Wiggle |
| **0700** | 6.71, 12.17, 18.03| Asymmetric Stretching |
| **0900** | 13.93, 13.93, 26.00 | **Decay / Expulsion of Knots** |

### 4. Avenue 4 Implication: The Necessity of a Confinement Term
The generic "Hot Soup" expansion reliably froze out a generic knot, but the *explicit structural isolation* of the UUD proton archetype **decayed**. 

To rule out linear-superposition artifacts during grid initialization, the simulation was rigorously patched to perform exact point-wise $S^3$ geometric quaternion multiplication when composing the initial UUD twist. The result remained identical: while the proper geometry allowed the configuration to hold together and "wiggle" together much longer (up to step 400), the underlying pure continuous topology could not permanently bind them. The knot distances stretched from ~5 to ~26 and the particle dissolved.

This definitively implies that pure $Cl^+(3,0,1)$ topology alone is **insufficient** for 3-quark confinement. To keep the Down quark bound to the Up-Up pair forever, the V13 Lagrangian *must* logically incorporate an explicit spatial rigidity term (such as the Skyrme $L_4$ term: $Tr([L_\mu, L_\nu]^2)$) to prevent the topological medium from infinitely stretching and breaking the color-singlet.
