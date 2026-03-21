# braid3(xyz): Three-Axis Crossed Braid with Chirality

## The Concept

Three braid3 structures along x, y, z axes, crossing at the origin.
Each braid carries all three fields with specific phase offsets.
The chirality (twist direction) of each braid determines UUD/UDD variants.

## Step-by-Step Construction

### Step 1: Rotate braid3 to each axis

The proven braid3(z) has three strands wound around z:
```
strand_a along z: center at (R_h cos(kz + 2πa/3), R_h sin(kz + 2πa/3), z)
φ_a = A × env(d_to_strand) × cos(kz + δ_a) × exp(-z²/2σ²)
```

Rotate this to x and y axes:
```
strand_a along x: center at (x, R_h cos(kx + 2πa/3), R_h sin(kx + 2πa/3))
φ_a += A × env(d) × cos(kx + δ_a) × exp(-x²/2σ²)

strand_a along y: center at (R_h cos(ky + 2πa/3), y, R_h sin(ky + 2πa/3))
φ_a += A × env(d) × cos(ky + δ_a) × exp(-y²/2σ²)
```

The three braids are superimposed: φ_a(total) = φ_a(x-braid) + φ_a(y-braid) + φ_a(z-braid) + background.

### Step 2: Chirality via phase sign

Right-handed (U): δ = {0, +3.0005, +4.4325}
Left-handed (D): δ = {0, -3.0005, -4.4325}

Each axis independently chooses U or D:

| Config | x-braid | y-braid | z-braid | Expected analog |
|--------|---------|---------|---------|----------------|
| UUU | U | U | U | Maximally symmetric |
| UUD | U | U | D | Proton-like |
| UDD | U | D | D | Neutron-like |
| DDD | D | D | D | = UUU by P-symmetry |

### Step 3: What happens at the crossing

At the origin (0,0,0), all three braids contribute to each field.
The total field is:
```
φ_a(0) ≈ 3 × A × cos(δ_a)  (if all same chirality)
       ≈ A × [cos(δ_a) + cos(δ_a) + cos(-δ_a)]  (if UUD)
       ≈ A × [cos(δ_a) + cos(-δ_a) + cos(-δ_a)]  (if UDD)
```

For UUU: φ₀(0) ≈ 3A, φ₁(0) ≈ 3A cos(3.0), φ₂(0) ≈ 3A cos(4.4)
For UUD: φ₁(0) ≈ 2A cos(3.0) + A cos(-3.0) = 3A cos(3.0) (same!)
For UDD: φ₁(0) ≈ A cos(3.0) + 2A cos(-3.0) = 3A cos(3.0) (same!)

Wait — at the origin, the AMPLITUDES are the same regardless of chirality
because cos(+δ) = cos(-δ). The chirality difference appears in the
VELOCITY field (through sin(±δ)) and in the CURL structure (which
determines the θ coupling and effective charge).

### Step 4: Where chirality matters — the curl

curl(φ) involves spatial derivatives. For a right-handed helix:
```
∂φ_a/∂z = -Ak sin(kz + δ_a) × env → positive curl component
```
For a left-handed helix:
```
∂φ_a/∂z = -Ak sin(kz - δ_a) × env → negative curl component
```

The curl DOES distinguish U from D. This means:
- The θ field (sourced by curl(φ)) differs between UUD and UDD
- The winding number differs
- The charge-dependent force differs

### Step 5: Mass difference

UUD and UDD have the same |φ|² (energy from amplitudes) but different
curl(φ) (energy from gradients). The gradient energy depends on
chirality because the phase derivatives differ:

E_grad(U) ∝ ∫ |∂φ/∂z|² with sin(kz + δ)
E_grad(D) ∝ ∫ |∂φ/∂z|² with sin(kz - δ)

At the CROSSING POINT, the gradient energies from the three axes
interfere differently for UUD vs UDD → MASS DIFFERENCE.

In real physics: m_neutron / m_proton = 1.00138.
The mass difference (1.3 MeV out of 938 MeV) comes from the
electromagnetic self-energy of the quark configuration.

In our theory: the mass difference should come from the θ-field
energy difference between UUD and UDD crossings.

## Variants to Test

### A: Simple superposition (easiest)
Three braid3 structures added linearly. Each braid independent.

### B: Amplitude-reduced superposition
Reduce each braid's amplitude to A/√3 so the total amplitude at
the crossing ≈ A (same as single braid). Prevents overdriving V(P).

### C: Weighted by axis
The x-braid carries 40% of φ₀, 30% of φ₁, 30% of φ₂.
The y-braid carries 30%, 40%, 30%.
The z-braid carries 30%, 30%, 40%.
This creates asymmetry between fields → richer spectrum.

### D: Different k per axis
Each braid has a slightly different wavenumber:
k_x = k, k_y = k×1.1, k_z = k×1.2.
This breaks the triple symmetry → richer dynamics.

## Testing Plan

1. **UUU** (all same): baseline, should be most symmetric/stable
2. **UUD** (proton): test mass and charge vs UUU
3. **UDD** (neutron): test mass and charge vs UUU and UUD
4. **UD-** (two-axis only, x+y): simpler test of crossing

For each: N=128, L=15, T=200, measure E_total, E_pot, aspect, θ_rms,
winding number.

## Expected Outcomes

**If crossing is stable**: A truly 3D, spherical, compact braid with
no preferred axis. The depletion profile should be spherical → the
gravity exponent may improve toward 1/r or 1/r².

**If UUD ≠ UDD mass**: The theory has a mass spectrum from chirality
alone — no additional parameters needed. The proton/neutron mass
difference emerges from geometry.

**If crossing is unstable**: The three braids may annihilate or repel.
Would need to test at different separations (braids offset from
origin rather than all crossing at one point).
