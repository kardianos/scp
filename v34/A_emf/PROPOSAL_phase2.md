# Phase 2: EMF from Braid Helical Charge and Torsion Waves

## What Phase 1a Taught Us

The linear modes around the uniform background are all massive and nearly
degenerate. The phase mode IS distinct from the amplitude mode, but both
have m_eff ≈ 1.5 in the far field. A simple "different linear mode = EM"
picture doesn't work.

But the mode splitting near the braid is enormous. The braid creates a
region where amplitude perturbations are nearly massless and phase
perturbations are strongly confined. This suggests the braid acts as
a fundamentally different medium for different perturbation types.

## The New Idea: Torsion Waves

The braid has a HELICAL structure — the three fields wind around each
other with winding number W = ±1. This winding is a topological property
that can't be changed by smooth perturbations.

**Torsion** ω_ij = ½(∂_i φ_j - ∂_j φ_i) measures the antisymmetric
part of the field gradient — how the field TWISTS in space. The braid
has nonzero torsion (V28: torsion flux Φ_T = 1.02).

A TORSION WAVE is a propagating twist perturbation: a region where the
helical structure temporarily tightens or loosens, traveling outward.

Key difference from amplitude waves:
- Amplitude wave: changes |φ| → changes ρ → gravity
- Torsion wave: changes the TWIST of φ → changes ω_ij → no ρ change

Torsion waves could be the EM analog:
- They don't create depletion (Σφ² unchanged to leading order)
- They interact with braids (which have specific twist structure)
- They can carry angular momentum (spin-1, like photons)
- Left and right-handed braids (W=±1) would interact oppositely
  with torsion waves → positive and negative "charge"

## Experimental Design

### Experiment 2a: Torsion Profile of a Settled Braid

**Goal**: Map the torsion field ω_ij(r) around a settled braid.

**Method**: Read a settled braid snapshot (from phonon test, T=200).
At each grid point compute:

    ω_ij = ½(∂_i φ_j - ∂_j φ_i)   for all i,j pairs

    |ω|² = Σ_{i<j} ω_ij²   (torsion magnitude squared)

    Torsion flux through xy-plane: Φ_T = ∫ ω_12 dx dy

Average |ω|² in spherical shells → torsion profile |ω|²(r).

**Questions**:
1. How does |ω|²(r) decay with distance? (Yukawa or power-law?)
2. Does the torsion profile differ from the depletion profile?
3. Is the torsion localized to the core or extended?

If |ω|²(r) has a DIFFERENT decay profile than δρ(r), then torsion
and amplitude perturbations propagate differently → they're distinct
force carriers.

### Experiment 2b: Torsion Wave Packet Propagation

**Goal**: Test whether a torsion perturbation propagates coherently.

**Method**: Initialize a traveling torsion pulse in the background
(no braid). The pulse locally twists the field without changing Σφ²:

    δφ_a(x) = ε × A_bg × sin(kz + 2πa/3) × envelope(x) × cos(k_x × x)

This is a 90° phase shift (cos→sin) of the background, modulated by
a Gaussian envelope in x, propagating in x. It changes the local
twist without changing the amplitude.

To be precise: at each point, rotate the field vector (φ₀,φ₁,φ₂) by
a small angle δθ(x) around an axis in field space. The rotation
changes the relative phases but not Σφ².

**Measure**:
1. Does the torsion pulse maintain shape as it propagates?
2. What is its speed? (c if massless, <c if massive)
3. Does it radiate amplitude waves (convert to gravity)?
4. Does Σφ² remain constant along the pulse? (no depletion = no gravity)

### Experiment 2c: Torsion Wave Hitting a Braid

**Goal**: Does a torsion wave interact with a braid differently than
an amplitude wave?

**Setup**: Settled braid at origin. Two separate runs:
1. Amplitude pulse from x=-20 → hits braid
2. Torsion pulse from x=-20 → hits braid

**Measure** for each:
1. Braid displacement after collision (momentum transfer)
2. Reflected/transmitted wave amplitude
3. Braid internal state change (E_pot, winding, asymmetry)
4. Is the interaction cross-section different?

If torsion waves scatter DIFFERENTLY than amplitude waves off braids,
the braid has distinct "gravitational" and "electromagnetic" responses.

### Experiment 2d: Left vs Right Braid Response

**Goal**: Do W=+1 and W=-1 braids respond differently to torsion waves?

**Setup**: Initialize two braids with opposite winding (δ → -δ gives
opposite helicity). Fire the SAME torsion pulse at each.

**Measure**: Is the scattering cross-section, momentum transfer, or
internal excitation DIFFERENT for opposite windings?

If YES: the braid's winding number acts as "charge" under torsion
interactions. W=+1 is "positive charge," W=-1 is "negative charge."
Torsion waves mediate the force between charges = electromagnetism.

---

## Implementation

### Torsion computation (for 2a)

The antisymmetric gradient tensor at each grid point:

```c
// ω_ij = ½(∂_i φ_j - ∂_j φ_i) for field components j and spatial dirs i
// With 3 fields and 3 spatial dimensions: 9 components,
// antisymmetric part has 9 independent components (3 fields × 3 antisym pairs)
// But the key quantity is ω² = Σ over all (i,a) < (j,b) combinations

for (int a = 0; a < 3; a++) {
    double dphi_dx = (phi[a][ip] - phi[a][im]) / (2*dx);
    double dphi_dy = (phi[a][jp] - phi[a][jm]) / (2*dx);
    double dphi_dz = (phi[a][kp] - phi[a][km]) / (2*dx);
    grad[a][0] = dphi_dx;
    grad[a][1] = dphi_dy;
    grad[a][2] = dphi_dz;
}
// Torsion: ω_{ij} = ½(grad[j][i] - grad[i][j]) for field-space indices
// This is a 3×3 antisymmetric matrix in field space
double omega_01 = 0.5*(grad[0][1] - grad[1][0]);  // actually this mixes
// field and space indices — need to be careful about what "torsion" means here
```

Actually, for three scalar fields, the "torsion" is the antisymmetric
part of the Jacobian J_ai = ∂φ_a/∂x_i:

    ω_ai,bj = ½(J_ai - J_bi)  — mixes field and space

The simplest torsion scalar is:

    T² = Σ_{a<b, i} (∂_i φ_a - ∂_i φ_b)²  — how differently the fields vary

Or the determinant form: det(J) = ε_{ijk} ∂_i φ₀ ∂_j φ₁ ∂_k φ₂
This is the 3D "volume form" — the Jacobian determinant.

For the braid, det(J) is nonzero and has a specific sign related to
the winding number W. It's the natural "charge density" of the braid.

### Torsion wave initialization (for 2b)

The simplest torsion perturbation: rotate the field vector at each
point by a small angle δθ(x) around the z-axis in field space:

    φ₀' = φ₀ cos(δθ) - φ₁ sin(δθ)
    φ₁' = φ₀ sin(δθ) + φ₁ cos(δθ)
    φ₂' = φ₂  (unchanged)

with δθ(x) = ε × exp(-(x-x₀)²/2σ²) × cos(k_wave × x)

This preserves Σφ² = φ₀² + φ₁² + φ₂² exactly (rotation preserves norm).
But it changes P = φ₀φ₁φ₂ (rotation in the 01-plane changes the product).

The velocity perturbation for a traveling wave:
    δ(∂_t φ_a) follows from the time derivative of the rotation.

### Resource estimate

- 2a: pure analysis on existing snapshot, ~5 minutes
- 2b: N=256, L=40, T=100, single run, ~15 minutes
- 2c: N=256, L=40, T=200, two runs, ~30 minutes each
- 2d: N=256, L=40, T=200, two runs, ~30 minutes each

Total: ~2 hours if sequential, ~1 hour parallel.

---

## What Success Looks Like

**Minimum**: Torsion profile |ω|²(r) decays differently from δρ(r).
This confirms amplitude and torsion are distinct at the nonlinear level.

**Strong**: Torsion wave propagates coherently, at c, without creating
depletion. This is a photon candidate.

**Very strong**: Torsion wave scatters differently off braids than
amplitude wave. The braid has distinct gravitational and EM responses.

**Grand slam**: W=+1 and W=-1 braids respond OPPOSITELY to the same
torsion wave. Winding number IS electric charge. EM emerges.
