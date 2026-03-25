# V40 Postulates: Composite Particle Formation from Chiral Sub-Structures

## Core Postulate

**Any stable composite particle in the Cosserat field theory that is not single-axis
oriented (i.e., not a simple braid) and is composed of a UUD or UDD chiral
configuration will require grid sizes of 350-512^3 to resolve.**

This follows from:
1. Each chiral sub-unit has a core width of ~6 code units (from braid measurements)
2. A 3-body composite with non-collinear orientation needs ~3× the single-braid extent
3. Plus buffer for depletion halos, interaction regions, and absorbing BC
4. Minimum: 3×6 (cores) + 3×10 (halos) + 2×15 (buffer) ≈ 78 code units per axis
5. At dx=0.3 (fine resolution): N = 78/0.3 ≈ 260. With safety: 350-512.

## Requirements for Composite Particle Candidates

### R1: Dynamic System (not static)

Every component of the particle, and the particle as a whole, must be a
time-dependent dynamical system. No static equilibria.

- Each sub-unit must oscillate (traveling wave, breathing, or rotation)
- The composite must have collective dynamics (not just independent parts)
- Stability comes from dynamic balance, not energy minimization to a fixed point
- **Metric**: amplitude oscillation persists for T > 500 without decay

### R2: Energy Binding

The composite must bind energy — its total energy must be LESS than the sum of
its isolated components.

- E_composite < E_unit1 + E_unit2 + E_unit3 (binding energy < 0)
- The binding must be robust: removing one unit raises the total energy
- **Metric**: binding energy ΔE/E > 1% (measurable above numerical noise)

### R3: Depletion Response (Gravity Analog)

The composite must respond to ρ gradients in a manner analogous to the braid's
gravity response:

- The composite depletes ρ in its neighborhood (creates a gravitational footprint)
- Other structures are attracted toward the depletion zone
- The depletion profile should be power-law (not Yukawa), consistent with V34 results
- **Metric**: δρ(r) measurable at r > 3× core radius

### R4: Non-Collinear Structure

The particle must have genuine 3D structure, not a 1D braid:

- At least two of three sub-units must have non-parallel orientation
- The composite must have a non-trivial moment of inertia tensor (I_1 ≠ I_2 ≠ I_3)
- **Metric**: aspect ratio eigenvalues differ by > 10%

### R5: Chiral Assignment

Each sub-unit must carry a definable chirality (U or D):

- Chirality = sign of the helical winding (k_z direction relative to envelope rotation)
- UUD or UDD configurations have net chirality ≠ 0 (like proton/neutron)
- The chiral asymmetry must persist dynamically (not average to zero)
- **Metric**: time-averaged net chirality > 0.5 (on a [-1, +1] scale)

### R6: Stability Under Perturbation

The composite must recover from small perturbations:

- Add 5% random noise to all field values → composite re-forms within T=50
- Remove one sub-unit → remaining structure either re-captures or ejects cleanly
- The binding potential well must have a positive second derivative (restoring force)
- **Metric**: survival probability > 80% under 5% perturbation over T=200

---

## What Makes the Braid Stable

### Analysis of Known Stability Mechanisms

The V34 single braid survives for T > 500 with the following properties:

1. **Traveling wave self-reconstruction**: The helical phase structure (k_z oscillation
   with phase offsets δ={0, 3.0, 4.43}) creates a traveling wave that continuously
   reconstructs itself. Energy that radiates outward is replenished by the wave's
   forward propagation.

2. **V(P) triple-product binding**: The three fields' product P = φ_0φ_1φ_2 enters
   V(P) = (μ/2)P²/(1+κP²). With μ < 0, regions where all three fields are nonzero
   have lower potential energy. This creates an energy well that binds the three
   field components together at the braid core.

3. **Phase offset stabilization**: The specific offsets δ={0, 3.0, 4.43} (from V28
   CMA-ES optimization) maximize the time-averaged |P| at the core. Random offsets
   produce weaker binding.

4. **Elliptical cross-section**: The ellipticity ε=0.3325 breaks rotational symmetry,
   preventing the braid from collapsing to a rotationally symmetric (and unstable)
   configuration.

5. **Background field coupling**: The A_bg=0.1 background provides a "medium" that
   supports the braid's traveling wave. Without background, the braid disperses faster.

6. **Theta coupling**: The η curl(θ) term transfers energy between φ and θ channels,
   providing additional stabilization through mode coupling. But it also provides a
   radiation channel — too much η causes dissolution (V39 result).

### Generalization: Stability Hierarchy

For a composite particle, stability requires **nested stabilization**:

```
Level 0: Individual field oscillations (mass term provides restoring force)
Level 1: Triple-product binding (V(P) couples the three fields)
Level 2: Traveling wave coherence (phase structure self-reconstructs)
Level 3: Sub-unit binding (depletion interaction between sub-units)
Level 4: Composite coherence (collective dynamics of bound sub-units)
```

Each level must be satisfied for the next level to be possible. A composite particle
needs all 5 levels, whereas a single braid needs only levels 0-2.

---

## Stability Quantification

### Proposed Stability Metric: S(t)

```
S(t) = w_E × (E_bind/E_total) + w_P × (P_int/P_int(0)) + w_θ × (θ_rms/θ_rms_eq) + w_A × (1/aspect_spread)
```

Where:
- E_bind/E_total = fractional binding energy (more negative = more stable)
- P_int/P_int(0) = triple product retention (1.0 = fully coherent)
- θ_rms/θ_rms_eq = theta equilibration (approaches 1.0 when theta reaches steady state)
- aspect_spread = max(I_i/I_j) - 1 (0 = spherical, higher = more anisotropic)
- Weights: w_E=0.4, w_P=0.3, w_θ=0.2, w_A=0.1

A single braid has S ≈ 0.6. The target for a composite is S > 0.7.

### Death Criteria

A structure is declared dead when:
- P_int drops below 10% of initial (field coherence lost)
- E_pot becomes positive (no binding left)
- phi_max drops below 2× background (dispersed to background level)
- Any of these sustained for T > 20 (rolling window, not instantaneous)

### Birth Criteria (Spontaneous Formation)

A structure is declared to have formed when:
- P_int exceeds 5× background level (locally coherent)
- E_pot is negative and deepening (binding increasing)
- phi_max > 3× background (locally concentrated)
- Sustained for T > 30

---

## Thermodynamic Argument: Why Particles Lower Energy

In a high-density ρ field, the system has excess energy (kinetic + gradient + mass).
Forming a localized structure (particle) can LOWER the total energy because:

1. **V(P) is negative**: Where three fields concentrate (P ≠ 0), the potential
   energy is negative (μ < 0). Concentrating field energy into a structure with
   large |P| extracts energy from the system.

2. **Depletion reduces background energy**: When a particle forms, it depletes ρ
   in its neighborhood. The depleted region has lower gradient energy, lower mass
   energy, and lower kinetic energy than the uniform high-density state.

3. **The net effect**: E_structured < E_uniform for the same total field integral.
   The difference is the binding energy, which is radiated away as waves.

This means: **in a sufficiently dense ρ field, particle formation is thermodynamically
favored**. The system will spontaneously nucleate structures that lower its energy.
The RATE of nucleation depends on the density and the energy barrier to forming
coherent phase structure.

### Prediction: Cooling Drives Particle Formation

If we start with a hot (random) high-density field and slowly cool it (absorbing BC
or explicit damping), particles should spontaneously form when the cooling rate is
slow enough for coherent structures to nucleate before the energy is radiated away.

The optimal cooling rate is: τ_cool ≈ 10 × τ_braid (where τ_braid ≈ 2π/ω ≈ 4
code time units). Faster cooling radiates everything before nucleation. Slower
cooling allows multiple nucleation events and potential composite formation.

---

## Evolutionary Search Strategy

### Phase 1: Braid Analysis (Baseline)

Analyze existing braid SFA to extract:
- Radial profiles of all 12 field arrays
- Phase structure (extract k_z, δ_a, envelope shape)
- Energy decomposition as function of radius
- Stability metric S(t) time series
- Identify which features correlate with stability

### Phase 2: Sub-Structure Search

For each candidate modification:
1. Start from a known stable braid (warm start)
2. Add a perturbation (second braid, rotated component, amplitude modulation)
3. Run for T=50 (short probe)
4. Measure S(t=50) and death/birth criteria
5. Select top-K candidates for longer runs

### Phase 3: Evolutionary Refinement

Tournament selection over modifications:
1. Take the top-K candidates from Phase 2
2. For each, generate M variants (parameter perturbation ±10%)
3. Run each for T=100
4. Select the top-K' by S(t=100)
5. Repeat until S plateaus or a death-free candidate emerges

### Phase 4: Composite Assembly

Once a stable sub-structure is found:
1. Place 2-3 instances in a larger box with specific orientations
2. Run at N=256-512 for T=200
3. Measure binding energy between units
4. Test UUD and UDD chiral configurations
5. Measure depletion response (R3)

---

## Multi-Resolution Strategy

Use the nested grid framework (scp_multi) for Phase 4:
- Root domain (N=64, L=50): captures long-range depletion and inter-unit coupling
- Child domains (N=128, L=10): one per sub-unit, resolves internal structure
- This reduces memory from 512^3 × 18 × 8 = 18 GB to ~2 GB
- Enables running on a single V100 (16 GB) instead of requiring A100

---

## Files

| File | Purpose |
|------|---------|
| `v40/POSTULATES.md` | This document |
| `v40/analyze_braid.py` | Braid SFA stability analysis (Phase 1) |
| `v40/evolve.sh` | Evolutionary search orchestrator |
| `v40/agent_search.md` | Auto-research agent instructions |
| `v40/agent_monitor.md` | Monitor agent instructions |
