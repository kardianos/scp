# V52 Plan: UUD vs UDD Characterization

## Background

V51 showed two UUD protons forming a bound oscillating state, which
contradicts real physics (the diproton is unbound). Algebraic analysis
(THEORY.md) proved that UUD+UDD has 19-69% less hardening repulsion
than UUD+UUD while feeling identical binding attraction. The mechanism
is chirality cancellation in curl_z: braid 1's opposite chirality in
UDD causes destructive interference in the ∂φ_x/∂y contribution.

V52 tests whether this algebraic prediction manifests in simulation
and whether it produces the correct physical hierarchy:
  |E_bind(UUD+UDD)| > |E_bind(UUD+UUD)|

## Theoretical Predictions (from predict_collision.mac)

| Quantity | UUD+UUD | UUD+UDD | Ratio UD/UU |
|----------|---------|---------|-------------|
| P at midpoint | P₀ | P₀ | 1.000 (exact) |
| \|curl(φ)\|² at overlap | 0.100 | 0.056 | 0.555 |
| E_hardening density | 1.19e-3 | 3.89e-4 | 0.325 |
| Free-flight t_contact | 125.0 | 125.0 | 1.000 |
| Breathing period | 2.2 | 2.2 | 1.000 |

Predictions to verify:
1. t_contact(UU) ≈ t_contact(UD) ≈ 125 (same D, same v)
2. D_min(UD) < D_min(UU) — neutron penetrates deeper
3. D_eq(UD) < D_eq(UU) — deuterium binds tighter
4. |E_bind(UD)| > |E_bind(UU)| — deuterium binds stronger
5. M²_peak(UU) > M²_peak(UD) — proton-proton has more strain
6. Measured |curl|² ratio ≈ 0.31-0.81 (D-dependent)
7. Breathing period unchanged across all experiments

## Prerequisites

### UDD (Neutron) Template

No pre-converged UDD template exists for C4 equations. Generate one:
1. Run gen_phase_confined with UDD chirality at N=192, L=25
2. Evolve under C4 equations for T=500 to equilibrate
3. Extract last frame as neutron_template.sfa (or use V43 neutron_template)

Alternative: use V43's neutron_template.sfa (evolved under V44 equations,
not C4). This is acceptable for initial tests since the template only
sets initial conditions — the C4 equations evolve from there.

### Collision Seed Generator

Modify gen_collision.c (or write gen_collision_pn.c) to stamp:
- Position 1: UUD proton template
- Position 2: UDD neutron template (with opposite chirality)
- Same D, same v_kick as the UU version

## Experiments

### Exp A: UUD+UUD at D=25, v=0.1 (control — repeat V51)

**Purpose**: Establish UU baseline with full post-collision tracking.
V51 ran T=300 but with 5-unit snap spacing. V52 reruns with denser
diagnostics and longer aftermath.

**Config**:
- N=384, L=50, T=500, absorbing BC
- snap_dt=5, diag_dt=0.25
- Seeds: 2× proton_template.sfa at ±12.5, v_kick=±0.1
- auto_download to /space/scp/v52/expA_uu/

**Predicted timeline**:
- t=0-85: free flight approach
- t≈85: interaction zone (D < 2R = 8)
- t≈125: first contact
- t=125-500: aftermath — bound orbit or separation

**Measure**:
- t_contact, D_min, D_eq (time-averaged D for t>250)
- E_bind = E_total(t>250) - E_total(t=0)
- T_orbit (FFT of centroid separation vs time)
- M² profile at D_min (mismatch_profile tool)
- Breathing period (FFT of P_int)

### Exp B: UUD+UDD at D=25, v=0.1 (the key test)

**Purpose**: Compare proton-neutron collision to proton-proton.
This is the experiment that tests the chirality-dependent hardening.

**Config**:
- Same as Exp A except right particle is UDD neutron
- Seeds: proton_template.sfa at -12.5, neutron_template.sfa at +12.5
- v_kick=±0.1

**Predicted differences from Exp A**:
- D_min should be SMALLER (19-69% less hardening barrier)
- D_eq should be SMALLER (tighter binding)
- |E_bind| should be LARGER
- M²_peak should be SMALLER (less strain at overlap)
- t_contact should be approximately the same (same D, v)

### Exp C: UUD+UUD at D=25, v=0.3 (high energy)

**Purpose**: Test whether higher energy disrupts proton structure.
At v=0.3, E_kin is 9× larger than v=0.1.

**Config**:
- Same as Exp A except v_kick=±0.3
- T=300 (shorter — faster dynamics)

**Predicted**:
- t_contact ≈ 42 (125/3)
- Higher probability of disruption or merger
- If protons survive: wider orbit (more energy to dissipate)
- If protons disrupt: characterize fragments

### Exp D: UUD+UDD at D=25, v=0 (gravitational infall)

**Purpose**: Test pure gravitational binding without kinetic energy.
The depletion attraction should pull them together slowly.

**Config**:
- Same as Exp B except v_kick=0
- T=2000 (very long — infall is slow)
- snap_dt=20, diag_dt=1.0

**Predicted**:
- Infall time >> 125 (no initial velocity)
- Gentle approach — minimal disruption
- Should form the most stable bound state (least excess energy)
- D_eq from this run is the "natural" equilibrium separation

### Exp E: UUD+UUD at D=25, v=0 (gravitational infall, control)

**Purpose**: Control for Exp D. Same setup but proton-proton.

**Config**: Same as Exp D but both UUD.

**Predicted**:
- Same infall time as D (same depletion attraction — P is identical)
- Higher D_eq than Exp D (harder shell prevents close approach)
- If UU does NOT bind at v=0 but UD does: this is the smoking gun

## Analysis Pipeline

For each experiment:

1. **Cluster tracking** (bin/cluster_profile):
   - Every frame: centroid positions, P_int, separation
   - Track particle identity across frames

2. **Mismatch profiling** (v52/mismatch_profile):
   - At D_min: extract M², |curl|², θ², P along x-axis
   - Compare UU vs UD profiles quantitatively

3. **Energy decomposition**:
   - E_total, E_pot, E_strain from diagnostics
   - Binding energy: ΔE = E(bound) - E(separated)

4. **Spectral analysis**:
   - FFT of P_int(t) → breathing frequency
   - FFT of D(t) → orbital frequency
   - Compare across experiments

5. **Maxima verification**:
   - Extract |curl|² from simulation data at D_min
   - Compare to Maxima prediction at same D
   - Compute residual: how much of the interaction is NOT
     captured by the midpoint algebra

## Run Schedule

| Exp | GPU | T | Est. wall time | Priority |
|-----|-----|---|---------------|----------|
| A | V100-16GB | 500 | ~35 min | HIGH (control) |
| B | V100-16GB | 500 | ~35 min | HIGH (key test) |
| C | V100-16GB | 300 | ~20 min | MEDIUM |
| D | V100-16GB | 2000 | ~2.5 hr | MEDIUM |
| E | V100-16GB | 2000 | ~2.5 hr | LOW (only if D binds) |

A and B can run sequentially on one GPU session (~70 min + overhead).
C can follow on the same session.
D and E are long runs — separate session or overnight.

Total GPU time: ~4 hours. Cost: ~$1.50 on V100.

## Success Criteria

**Strong success**: D_min(UD) < D_min(UU) AND |E_bind(UD)| > |E_bind(UU)|
with measured |curl|² ratio matching Maxima prediction within 2×.

**Weak success**: Any measurable difference between UU and UD collision
dynamics (D_min, D_eq, E_bind, M²_peak) with correct sign.

**Failure**: UU and UD collisions are indistinguishable. This would mean
the midpoint algebra doesn't capture the real interaction, and the
distributed shell structure dominates over the chirality effect.

**Unexpected success**: UU does NOT bind at v=0 while UD does. This
would be the Pauli exclusion analog emerging from hardening alone.

## Source Files

- `v52/THEORY.md` — Theoretical results and proofs
- `v52/predict_collision.mac` — Numerical predictions
- `v52/predict_collision_output.txt` — Raw Maxima output
- `v52/core_interaction.mac` — Two-proton overlap algebra
- `v52/uud_udd_interaction.mac` — UU vs UD comparison
- `v52/CoreInteraction.lean` — Breathing invariance proofs
- `v52/UUDvsUDD.lean` — Hardening reduction proofs
- `v52/mismatch_profile.c` — Mismatch extraction tool
