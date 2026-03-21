# Electron Shell Plan: Multi-Scale Approach

## The Scale Problem

The electron orbit is ~50,000× the nuclear radius. A uniform 3D grid
can't bridge this. Solution: extract the braid's effective potential
from 3D simulation, then solve for bound states in 1D.

## Step 1: Extract Effective Potentials (from existing data)

We already have the braid's radial profiles from V34:
- δρ(r): depletion profile (gravity potential, power-law ~1/r^1.2)
- θ_φ(r): azimuthal θ field (EM potential)
- θ²(r): total θ energy (coupling strength)

From these, construct the effective radial potentials:
- V_grav(r) = -G_eff × M / r^1.2  (from the depletion profile)
- V_em(r) = -α_eff × Q / r^n      (from the θ force measurement)

The force measurements give us the constants:
- G_eff: from F = -C × ∇ρ (C ≈ 186)
- α_eff: from the 27% force enhancement at η=0.5

## Step 2: Radial Schrödinger Equation (1D, no 3D grid needed)

For a test perturbation in the θ potential well:

    -ℏ²/(2m_eff) × d²ψ/dr² + [V_eff(r) + ℏ²l(l+1)/(2m_eff r²)] ψ = E ψ

    V_eff(r) = V_em(r)  (the θ-mediated Coulomb-like potential)

This is the standard radial Schrödinger equation. Solve on a 1D grid
with N_r = 100,000 points from r=1 to r=500,000 (logarithmic spacing).
Trivial computation — runs in seconds.

The "electron" is a bound eigenstate of this equation. Its radius is
the Bohr radius analog. Different l values give different orbital shapes.

## Step 3: Where Does ℏ Come From?

Three options to investigate (not mutually exclusive):

### 3a: Field Quantum (ε)

The substrate has a minimum field step ε. This creates:
- Minimum momentum: p_min ~ ε per grid cell
- Uncertainty: Δx × Δp ~ ε × dx
- Effective ℏ_sim = ε × dx

Test: add field rounding (Track F from original plan) to the Cosserat
simulation. Does it create a minimum orbital radius?

### 3b: Grid Discreteness (dx)

The grid itself provides a minimum wavelength λ_min = 2dx (Nyquist).
A perturbation can't be localized to less than dx → Δx ≥ dx.
This gives ℏ_grid = p_min × dx where p_min is the minimum momentum
supportable on the grid.

Test: run the 1D radial equation with different dx and check if the
bound state radius converges or depends on dx.

### 3c: Intrinsic from the Dynamics

The braid's oscillation frequency ω sets a natural action scale:
ℏ_intrinsic = E_braid / ω. If the braid has E ≈ 5000 and ω ≈ 1.5,
then ℏ_intrinsic ≈ 3300. This is enormous compared to the field
values, which would make the Bohr radius very large — possibly
matching the physical ratio.

## Option A: θ-Sector Solitons

Can the θ field form its own localized structures (θ-braids)?

Currently θ has NO self-interaction — it's a free massless field
sourced by curl(φ). Without nonlinearity, θ can't form solitons.

Test: add a small θ self-interaction, e.g.:
    V_θ = (μ_θ/2)(θ₀θ₁θ₂)² / (1 + κ_θ(θ₀θ₁θ₂)²)

Same triple-product form as V(P) but for the angle fields.
If θ can form its own braids:
- θ-braids would be lighter (massless field → lighter binding)
- They could orbit in the φ-braid's potential well = electron
- Their winding would be the electron's charge

Experiment: add V_θ to the Cosserat equation, initialize a θ
perturbation near a φ-braid, see if it binds.

## Option B: Frequency Modes

The electron isn't at a specific radius — it's at a specific FREQUENCY.

The θ field around the braid oscillates at multiple frequencies.
If there are DISCRETE resonant frequencies (like a drum), each
resonance is an "energy level" and the frequency determines the
"orbital radius" through E = ℏω → r ~ ℏ²/(mE).

Test: Fourier-analyze the θ field at various radii from the braid.
Look for discrete peaks in the frequency spectrum. If the spectrum
is continuous → no quantization. If discrete → natural energy levels.

## Computational Efficiency

### For Steps 1-3 (radial equation): Already efficient
- 1D computation: N_r = 100,000 points → seconds on a laptop
- No 3D grid needed for the electron problem
- The 3D simulation only needs to characterize the braid (done)

### For Option A (θ-braids): Same 3D code, moderate grid
- N=128 is sufficient — the θ-braid would be at braid scale, not Bohr scale
- The ORBIT would be at Bohr scale, but the θ-braid itself is compact

### For larger 3D simulations if needed:
- GPU: the force computation is embarrassingly parallel. A CUDA port
  would give 50-100× speedup. N=512 would run in minutes, not hours.
- AMR: adaptive mesh refinement — fine grid at braid, coarse elsewhere.
  Complex to implement but handles multi-scale naturally.
- Spectral methods: expand in spherical harmonics for the electron problem.
  Naturally separates scales.

## Priority Order

1. **Extract V_eff(r) from existing θ data** (pure analysis, 1 hour)
2. **Solve radial eigenvalue problem** (1D code, 1 hour)
3. **Frequency analysis of θ field** (Option B, from existing SFA data)
4. **Add θ self-interaction** (Option A, code modification + simulation)
5. **Add field quantization ε** (Step 3a, code modification + simulation)
