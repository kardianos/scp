# T12 Deep Characterization: M2, M4, M7

## What Step 0 Missed
- Only 6 time snapshots (coarse)
- No full field dumps (can't do post-hoc analysis)
- No radial profile of depletion rate dρ/dt(r)
- No comparison of depletion SHAPE (1/r? exponential? ring?)
- No energy flow analysis (where does consumed energy GO?)
- No conceptual analysis of what the equations MEAN

## What We Want

### Per mechanism (M2, M4, M7): Fine-grid 3D run
- **N=128, L=30, T=500** (double previous resolution and time)
- **50 diagnostic snapshots** (every T=10)
- **At each snapshot**:
  - Radial profile ρ(r) with 100 bins (dr=0.3)
  - Differential profile δρ(r) = ρ_braid(r) - ρ_control(r)
  - Energy decomposition: E_kinetic(r), E_gradient(r), E_mass(r), E_pot(r)
  - Core metrics: fc, |P|, winding
  - Total energy, energy in shells (r<5, 5-10, 10-20, 20-30)

- **Full field dumps** at T=0, 100, 200, 300, 400, 500:
  - Save φ_a(x,y,z_mid) for z=z_mid (2D xy slice through center)
  - Save φ_a(x_mid, y_mid, z) for all z (1D z-axis profile)
  - Save ρ(x,y,z_mid) energy density slice
  - Format: binary float arrays with header (N, L, t)
  - These enable post-hoc visualization and analysis

- **Depletion rate measurement**:
  - At each snapshot, compute dρ/dt(r) = [ρ(r,t) - ρ(r,t-dt)] / dt
  - Track whether |dρ/dt| → 0 (stabilization) or persists
  - Fit depletion profile at late time: δρ ~ A/r^α → extract α

### Also run control (no mechanism modification, just braid + background)
This is the "M0" baseline for comparison.

## Conceptual Analysis of the Equations

### M2: Back-Pressure (Gradient Force)

Standard EOM:
    ∂²φ_a/∂t² = ∇²φ_a - m²φ_a - ∂V/∂φ_a

M2 adds:
    ∂²φ_a/∂t² = ∇²φ_a - m²φ_a - ∂V/∂φ_a + β × ∂ρ/∂x_a

where ρ(x) = local energy density (smoothed).

**What this means physically**: Each field feels a force ALONG the
gradient of the total energy density. Where ρ drops (depletion zone),
the gradient points inward → force pushes field INTO the depleted zone
→ opposes further depletion.

**This is analogous to**: Pressure in a fluid. Low-density regions
have lower pressure, and the pressure gradient drives flow toward them.
The coefficient β is the "compressibility" of the field medium.

**Limitation**: The force ∂ρ/∂x_a uses field index a as spatial index.
This means φ₀ feels the x-gradient, φ₁ the y-gradient, φ₂ the z-gradient.
This is the ELASTIC interpretation (field a = displacement in direction a).
It's natural for the strain/torsion framework but may introduce anisotropy.

**Lagrangian?**: This force does NOT derive from a potential in the usual
sense. It's a dissipative-like force (depends on ρ which depends on
velocities). NOT Hamiltonian.

### M4: Saturating Potential (Depletion Energy Cost)

Standard EOM:
    ∂²φ_a/∂t² = ∇²φ_a - m²φ_a - ∂V/∂φ_a

M4 adds:
    ∂²φ_a/∂t² = ∇²φ_a - m²φ_a - ∂V/∂φ_a - ∂V_depl/∂φ_a

where V_depl = λ(ρ₀ - ρ)² and ρ = ½Σφ_a².

    ∂V_depl/∂φ_a = -2λ(ρ₀ - ρ) × ∂ρ/∂φ_a = -2λ(ρ₀ - ρ) × φ_a

So the additional force is: +2λ(ρ₀ - ρ) × φ_a

**What this means physically**: The field has a preferred local amplitude
|φ|² = ρ₀. Depletion (|φ|² < ρ₀) creates a restoring force that pushes
|φ| back toward ρ₀. This is a RADIAL force in field space, not a
spatial force.

**This IS Lagrangian**: Derives from V_depl = λ(ρ₀ - ½Σφ²)².
Energy is well-defined and conserved (up to boundary effects).

**Connection to Higgs**: V_depl is a "sigma-model" potential that prefers
a specific field magnitude. Unlike the Mexican hat (which has V(ρ-ρ₀)²
with ρ = |φ|), this uses ρ = ½Σφ², which includes ALL three fields.

**Limitation**: The force is proportional to φ_a → stronger in the braid
core where φ is large. This could interfere with the braid structure.
The Step 0 results suggest it doesn't (fc=0.68), but finer resolution
may reveal issues.

### M7: Two-Component (Structural + Background)

Two separate field sets: S_a (braid) and B_a (background).

    ∂²S_a/∂t² = ∇²S_a - m²S_a - ∂V_S/∂S_a - g × B²_total × S_a
    ∂²B_a/∂t² = ∇²B_a - m²B_a                - g × S²_total × B_a

where B²_total = ΣB_a², S²_total = ΣS_a².

**What this means physically**: S and B are two populations of the
same field that interact weakly (coupling g). The S fields form the
braid; the B fields form the background. The coupling g×S²×B converts
B energy into S energy at the braid core (consumption). The reverse
coupling g×B²×S converts S back into B (radiation).

At equilibrium: consumption rate = radiation rate → self-limiting.

**This is Lagrangian**: L = L_S + L_B - g × S²_total × B²_total.
The interaction is a quartic coupling between the two sectors.

**Advantage**: The braid structure is protected because S and B are
separate. Modifications to B (depletion) don't directly alter S.
The coupling is weak (g=0.01) → the braid barely notices the
background being consumed.

**Disadvantage**: This doubles the field content (6 → 12 real fields).
And the S-B separation may be artificial — in reality, the braid IS
the field, not a separate component.

## Combining Terms: M2+M4, M4+M7

### M2+M4 (Gradient + Saturating):
Back-pressure provides the spatial gradient force (field flows toward
depletion), while the saturating potential provides the amplitude
restoring force (field magnitude returns to ρ₀). These are complementary:
- M2 controls WHERE the field goes (spatial redistribution)
- M4 controls HOW MUCH field there is (amplitude regulation)

Combined EOM:
    ∂²φ_a/∂t² = ∇²φ_a - m²φ_a - ∂V/∂φ_a
                  + β × ∂ρ/∂x_a           [M2: spatial flow]
                  + 2λ(ρ₀ - ½Σφ²) × φ_a   [M4: amplitude regulation]

### M4+M7 (Saturating + Two-Component):
Two-component provides the braid protection, while M4's saturating
potential acts on the B fields to limit depletion.

    ∂²S_a/∂t² = ∇²S_a - m²S_a - ∂V_S/∂S_a - g×B²×S_a
    ∂²B_a/∂t² = ∇²B_a - m²B_a - g×S²×B_a + 2λ(ρ₀-½ΣB²)×B_a

The B fields have both the coupling to S AND the restoring force toward
ρ₀. The braid (S) is unmodified except through the weak coupling.

### M2+M7 (Gradient + Two-Component):
Back-pressure acts on B fields, driving them toward depleted regions.

    ∂²B_a/∂t² = ∇²B_a - m²B_a - g×S²×B_a + β×∂ρ_B/∂x_a

## Output Format for Full Field Dumps

Binary format for each dump:
```
Header: int32 N, float64 L, float64 t, int32 n_fields
Data: float64[N*N] for each field's xy-slice at z=N/2
```

This allows loading in Python/C for post-hoc analysis (contour plots,
FFT, correlation functions, etc.).
