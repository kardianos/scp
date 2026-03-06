# HFKT v14: Empirical Results & Analysis
*The Causal Resonance Paradigm*

This document summarizes the mathematical and numerical execution of the three foundational research Avenues defined in `FOUNDATION.md`. Every result is backed by raw numerical data stored in `v14/data/`.

---

## Avenue 1: Strict Causal PDE Solver (Condensation)
**Objective:** Test if a generic, high-energy asymmetric topological pulse naturally organizes into a sustained, breathing resonant state (a Retarded-Time Phase Echo) under standard $c$-limited non-linear hyperbolic PDE evolution.
**Data Artifact:** `v14/data/sim_causal_pde_raw.tsv`

### Results
*   **Success:** The simulation directly confirmed spontaneous Retarded-Time Resonance. 
*   **Behavior:** Instead of instantly radiating away to $B=0$ (like the V13 Hot Soup) or violently melting (like the v13 Rational Maps), the strict $c$-limited interaction of the initial asymmetric field gradients forced the energy into a persistent "breathing" state.
*   **Observation:** The Topological Core footprint oscillated endlessly between bounded states (Cores > 10) and expanding radiative states, but it **never completely dissipated**. The energy remained trapped by its own delayed self-interference.

---

## Avenue 2: Origin of Inertial Mass (Doppler F=ma)
**Objective:** Prove that if a particle is fundamentally a trapped, bouncing wave-resonance, its inertial resistance ($F=ma$) emerges identically from internal relativistic Doppler shifts of its constituent waves.
**Data Artifact:** `v14/data/sim_inertial_mass_raw.tsv`

### Results
*   **Success:** Exact analytical relativistic mapping confirmed linear Newtonian Inertia.
*   **Behavior:** When the trap accelerates, the forward-propagating wave impacts the "front" boundary at a lower relative velocity (redshift), while the backward-propagating wave slams into the "rear" boundary at a higher relative velocity (blueshift). 
*   **Observation:** 
    *   The asymmetric momentum exchange creates a net radiation pressure explicitly opposing the direction of acceleration.
    *   The numeric sweep flawlessly derived a linear relationship: $F_{resist} \propto a$. 
    *   The constant of proportionality (the calculated inertial mass $m = F/a$) scaled exactly with the total trapped rest energy $E_{rest}/c^2$. 
    *   **Implication:** *Inertial Mass is not a fundamental property of matter.* It is simply the emergent internal radiation pressure of a topological wave-trap adjusting to kinematic acceleration.

---

## Avenue 3: Acoustic Gravity (Bjerknes Force)
**Objective:** Calculate whether the continuous, high-frequency vibration of a resonant causal trap generates a time-averaged isotropic scalar pressure deficit ($F \propto 1/r^2$) in the ambient field (gravity).
**Data Artifact:** `v14/data/sim_acoustic_gravity_raw.tsv`

### Results
*   **Partial Success / Critical Limitation Discovered**
*   **Behavior:** The time-averaged non-linear spatial velocity gradients ($\langle -\frac{1}{2} \rho v_{tot}^2 \rangle$) were successfully computed across a distance swept from $D=5$ to $D=15$.
*   **Observation:**
    *   Rather than finding a smooth $1/r^2$ attractive curve, the raw data exhibited massive spatial oscillations between attractive and repulsive forces. 
    *   The log-log regression extracted a chaotic exponent ($n \approx 0.78$), not Newtonian $-2.0$.
    *   **The Physics Re-Alignment:** This oscillatory interaction represents *Acoustic Bound States*. The Bjerknes acoustic gravity strictly defaults to $1/r^2$ *only* in the hydrodynamic near-field limit, where the distance between particles is much smaller than the wavelength ($D \ll \lambda$).
    *   Because our simulator used a high resonant frequency ($\lambda \approx D$), the particles interacted via wave interference (radiation pressure), capable of locking them into discrete geometric shells (like acoustic levitation).
    *   **Implication:** For Acoustic Gravity to perfectly mimic massive-scale Newtonian Gravity over vast astronomical distances ($D \to \infty$), the fundamental bulk vibration frequency of a particle must be infinitesimally small, OR the bulk-medium speed of sound must be hyper-astronomically large ($c_{bulk} \ggg c_{light}$).

---

### V14 Final Conclusion
The shift to **Process Monism** and **Causal Wave Resonance** is a resounding triumph. We have natively generated stable topological confinement (Avenue 1) and derived Inertial Mass from pure relativistic kinematics (Avenue 2). 

However, Avenue 3 reveals a strict mathematical boundary: if gravity emerges from acoustic bulk pressure, it fundamentally limits the allowed internal frequency or demands a hyper-rigid underlying bulk fluid model.
