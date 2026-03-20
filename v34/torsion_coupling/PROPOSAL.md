# Torsion Coupling: Adding Angular Stress to the Field Equation

**Goal**: Add the simplest possible cross-gradient coupling between the
three fields and observe how it affects a braid. Not trying to get EMF
yet — just understanding how the field responds when it can "feel" twist.

---

## The Modification

### Current equation (three independent scalars):

    ∂²φ_a/∂t² = ∇²φ_a - m²φ_a - ∂V/∂φ_a                    (1)

Each φ_a evolves based on its OWN Laplacian. The fields only talk to
each other through V(P) at each point. No cross-field gradient coupling.

### Modified equation (elastic solid with shear):

    ∂²φ_a/∂t² = ∇²φ_a - m²φ_a - ∂V/∂φ_a + γ ∂(div φ)/∂x_a  (10)

    where div φ = ∂φ₀/∂x + ∂φ₁/∂y + ∂φ₂/∂z

This identifies field index a with spatial index i:
    φ₀ ↔ x-displacement,  φ₁ ↔ y-displacement,  φ₂ ↔ z-displacement

The new term γ ∂(div φ)/∂x_a is the gradient of the divergence —
the standard elastic cross-coupling that separates longitudinal
(compression) waves from transverse (shear) waves.

### What this gives:

- Longitudinal waves (compression): speed c_L = √(1 + γ)
- Transverse waves (shear): speed c_T = 1 (unchanged from current)
- At γ = 0: recovers the current equation (all waves at c = 1)
- At γ > 0: compression faster than shear (normal solid behavior)
- At γ < 0: shear faster than compression (exotic but worth testing)

### Why this is natural:

CONCEPT.md says: "The triple product P = φ₀φ₁φ₂ is the VOLUME FORM
of the three-field displacement gradient." If the fields ARE a 3D
displacement, the elastic cross-coupling should be there. We dropped
it for simplicity, but it's part of the original physical picture.

The identification φ_a ↔ x_a is not arbitrary — the braid has a
helical structure that TWISTS the field-space displacement relative to
physical space. This twist is torsion. Without cross-gradient coupling,
the equation can't respond to it.

---

## Implementation

### Force computation change:

Two-pass approach:

**Pass 1**: Compute div(φ) at each interior point:
```c
double div = (phi[0][ip_x] - phi[0][im_x]) / (2*dx)    // ∂φ₀/∂x
           + (phi[1][jp_y] - phi[1][jm_y]) / (2*dx)    // ∂φ₁/∂y
           + (phi[2][kp_z] - phi[2][km_z]) / (2*dx);   // ∂φ₂/∂z
div_array[idx] = div;
```

**Pass 2**: Compute forces including ∂(div)/∂x_a:
```c
// For field a=0: ∂(div)/∂x
double ddiv_dx = (div_array[ip_x_idx] - div_array[im_x_idx]) / (2*dx);
acc[0][idx] = lap - MASS2*phi[0][idx] - mPd2*dPda + GAMMA * ddiv_dx;

// For field a=1: ∂(div)/∂y
double ddiv_dy = (div_array[jp_y_idx] - div_array[jm_y_idx]) / (2*dx);
acc[1][idx] = lap - MASS2*phi[1][idx] - mPd2*dPda + GAMMA * ddiv_dy;

// For field a=2: ∂(div)/∂z
double ddiv_dz = (div_array[kp_z_idx] - div_array[km_z_idx]) / (2*dx);
acc[2][idx] = lap - MASS2*phi[2][idx] - mPd2*dPda + GAMMA * ddiv_dz;
```

### Memory: need one extra array for div(φ)

Extend allocation from 9 to 10 arrays (phi[3] + vel[3] + acc[3] + div[1]).

### Energy: the torsion coupling adds elastic energy

    E_elastic = (γ/2) × (div φ)² × dV   (integrated over volume)

This is the compression energy. Total energy = E_kin + E_grad + E_mass + E_pot + E_elastic.

---

## Experiment

### Phase 1: γ scan on a single braid

For γ in {0, 0.1, 0.5, 1.0, 2.0, 5.0, -0.1, -0.5}:
- Single braid, N=128, L=20, T=300
- Record: E_total (including E_elastic), E_pot, braid position, winding
- Save snapshots at t=0, 100, 200, 300

**Key questions**:
1. Does the braid survive? (E_pot retention)
2. Does the braid change shape? (compare snapshots)
3. Does the energy stay conserved? (check E_total including E_elastic)
4. Does the braid's internal oscillation frequency change?
5. At what γ does the braid break?

### Phase 2: Compare modes around the braid

For γ that gives a surviving braid:
- Compute the longitudinal and transverse components of the field:
  φ_L = ∇(div φ) / |∇(div φ)|  (compression component)
  φ_T = φ - φ_L  (shear component)
- Map their radial profiles separately
- Does the shear component extend further than the compression component?

### Phase 3: Moving braid torsion wake (if Phase 1 succeeds)

If the braid survives with γ > 0:
- Give the braid a velocity kick
- Measure the shear/torsion pattern around the MOVING braid
- Is it circular (perpendicular to motion)?
- Does it depend on winding number W?

---

## Expected Outcomes

**γ = 0.1**: Minor perturbation. Braid should survive with small
changes. The elastic energy is a small correction.

**γ = 1.0**: Significant coupling. Longitudinal waves at c_L = √2.
The braid's internal dynamics will be modified — compression modes
propagate faster than shear modes, which could change the braid's shape.

**γ = 5.0**: Strong coupling. c_L = √6 ≈ 2.45. May destabilize the
braid if the compression speed exceeds the braid's internal coherence.

**γ < 0**: Exotic. Shear waves faster than compression. The braid
might become MORE stable (shear binding strengthens) or LESS stable
(compression weakens). Unknown territory.

## What we learn regardless of outcome

If the braid survives: we know the elastic cross-coupling is compatible
with braid existence. This opens the door to shear-mediated forces (EMF).

If the braid breaks: we learn what the coupling threshold is, and
whether negative γ (enhanced shear) helps.

If the braid changes shape: we learn how the internal structure
responds to the shear/compression asymmetry, which tells us about
the braid's elastic properties.
