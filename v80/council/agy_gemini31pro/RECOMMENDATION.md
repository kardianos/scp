### A. Timeline reading

*   **v74–v75:** Successfully established the Stage 2A liquid-drop Z-carbon nucleus and designed the multi-fabric architecture. E-lite forces were measured (opposites attract, sames repel), but orbit/capture remained partial or classical.
*   **v76–v77:** Shifted toward a monist free-capacity path-cost representation; implemented full Maxwell + dual-channel RC1 co-field with a shared speed of light ($c$).
*   **v78:** Froze R1–R10 physics relations and particle recipes, but the ideal $C_{12}$ atom product was blocked.
*   **v79:** Attempted a multi-fab $Z=6$ atom (long-T). Resulted in a failure: the L shell nullified $E_{em}$ and flux, creating a neutralizing radiation bath rather than a stable, parked multipole.
*   **v80:** Overnight V100 campaign (Goal = 0.62, passing threshold). Crucially established that a minimal $C+L$ hydrogenoid-class system *holds* $E_{em}$ and $Q_{\phi}$ over long T, contrasting sharply with the v79 $Z=6$ collapse. 
*   **Open Issues (v80):** Measured Coulomb relative acceleration ($a_{rel}$) pair tracks are missing. Multi-rev orbits have not been measured (only mild tangential velocities tested with some energy drift). Monist free-substrate representation is theorized but not implemented as a primary state.

### B. Long-term goal status

*   **Stages 0–2A (Carbon map, Nucleon template, Liquid-drop nucleus):** Complete. 
*   **Stage 3 (Light opposite-charge sector):** Currently the primary blocker. v80 proved that a minimal $C+L$ hydrogenoid pair can hold its electromagnetic energy without collapsing (unlike heavier configurations), but stable orbital mechanics and force-law confirmation are still unverified.
*   **Stage 4 (Atom = C-nucleus + 6 opposite light charges):** Blocked by Stage 3. v79 demonstrated that naive placement of 6 L-charges around a $Z=6$ core results in charge death/radiation bath.
*   **Real progress criteria:** To unblock the path to the carbon atom, we must demonstrate a stable, multi-revolution orbit of a single L-charge around a single core charge, supported by explicitly measured Coulomb relative acceleration tracks that confirm the force mechanics.

### C. Recommended next step (primary)

**Action:** Execute a dedicated Hydrogenoid ($C+L$) Orbital & Force Tracking Campaign.

*   **What to run/build:** 
    1. Run `mf_pair_track` to explicitly extract and measure $a_{rel}$ (Coulomb relative acceleration) between a single C and L charge to confirm force laws without SFA pruning.
    2. Run `gen_mf_shell_orbit` sweeping a wider range of tangential velocities ($v_t$) to achieve and measure a sustained, multi-rev orbit for the $C+L$ system.
*   **Why this before alternatives:** v80 proved the $C+L$ system does not immediately collapse (holding $E_{em}$), making it the most viable candidate. We cannot scale to $C_6+L_6$ (Stage 4) without first proving that a single L-charge can sustainably orbit a core without drifting or decaying into a radiation bath. 
*   **Success / kill criteria:** 
    *   *Success:* `mf_pair_track` yields a clean inverse-square acceleration profile; orbit sweep produces at least one configuration with $>2$ revolutions and stable $E_{em}$ / radius.
    *   *Kill:* Orbits consistently decay into a radiation bath or drift to infinity regardless of $v_t$ tuning, indicating the multi-fab mechanics cannot support stable bound states.
*   **Estimated effort:** Short GPU runs (sweeping $v_t$) and CPU analysis for the pair tracks.
*   **Kernel auth needed?** No. Rely entirely on existing configs and seed generators.

### D. Explicitly defer or kill

*   **KILL: Re-running $Z_6+L_6$ (or $C_6+L_6$) parked atoms.** v79 proved this acts as a neutralizing bath. Until the hydrogenoid ($C+L$) orbital mechanics are completely solved, trying to park 6 electrons is a waste of compute.
*   **DEFER: Multi-center carbon (Stage 2B).** Focus must remain on solving the opposite-charge light sector (Stage 3) which is currently blocking the core definition of an atom.
*   **KILL: Modifying `scp_sim` or `sfa.h`.** Do not attempt to fix orbit issues by hacking the kernel. The v78 physics relations are frozen; solutions must be found via initial state configs and seed generators.

### E. Product vs representation

Over the next 1–2 weeks, split effort approximately **75% Product / 25% Representation**:

*   **Product (75%):** Focus heavily on the multi-fab product ladder using the frozen v78/v80 kernel. The immediate goal is extracting the missing pair tracks and achieving multi-rev orbits via `gen_mf_shell_orbit`. This is the critical path to unblocking Stage 3 and moving toward Stage 4.
*   **Representation (25%):** Address the kimi k3 review thesis by building the recommended "2D free-medium + typed locks" toy model. This should be a completely separate, lightweight CPU script—isolated from the main kernel—used solely to validate the density-only monist math for a potential future (v81+) architecture rewrite.

### F. Confidence

**Score: 4.5 / 5**

The data strongly points to the $C+L$ hydrogenoid system as the vital stepping stone. The v79 failure vs v80 success in holding $E_{em}$ dictates that we must stay at the $N=1$ scale to solve the orbital mechanics. The deduction to use existing seed generators for force tracking directly addresses the explicitly noted missing evidence ("no $a_{rel}$ pair tracks", "multi-rev not measured") without violating the strict "no kernel modification" constraint.
