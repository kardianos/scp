# Computation 5: Trapped θ Wave Packet as Electron

## The Idea

The electron is NOT a θ-soliton (ruled out — θ lacks mass gap for
localization). Instead, it's a θ WAVE PACKET trapped in the φ-braid's
potential well, like a particle in a box.

The φ-braid creates a potential well through two mechanisms:
1. The depletion zone: lower ρ → different effective dynamics for θ
2. The curl coupling: η×curl(φ) directly drives θ near the braid

A θ wave packet placed in this well would bounce back and forth
(standing wave). Its effective mass comes from the confinement:

    m_eff = E_kinetic / c² = ℏ²k² / (2c²)

where k ~ π/R_well is set by the well size R_well.

## Key difference from Computation 4

Comp 4 (θ-soliton): needed θ to self-bind through its OWN potential.
Failed because θ has no mass gap → uniform condensate.

Comp 5 (trapped packet): θ is CONFINED by the φ-braid's potential well.
No θ self-interaction needed. The φ-braid provides the "box."

This is EXACTLY how a real electron works: it doesn't self-bind.
It's trapped by the Coulomb potential of the nucleus.

## Method

### Step 1: Measure the θ trapping potential

In the existing 6-field Cosserat simulation (no θ self-interaction),
the φ-braid modifies the local θ dynamics through:
- η×curl(φ): direct source term (strongest near braid)
- The background field: θ propagates differently where φ is structured

Measure: place a small θ wave packet at various distances from the
braid. Does it stay, drift inward, drift outward, or disperse?

### Step 2: Initialize a localized θ wave packet

NOT a helical braid pattern (that's a soliton ansatz). Instead:
- A Gaussian envelope in θ at some distance from the φ-braid
- θ_a(r) = A_θ × exp(-(r-r_0)²/2σ²) × cos(k·z + 2πa/3)
- Choose r_0 = 8-15 (outside the braid core, in the depletion zone)
- Choose σ = 2-3 (localized but not point-like)
- Choose k = π/L (one half-wave, low frequency)

### Step 3: Track the wave packet

Run the Cosserat simulation (η=0.5, m_θ=0, NO θ self-interaction):
- Does the θ packet maintain coherence?
- Does it orbit the φ-braid?
- Does it disperse outward (escapes) or inward (captured)?
- What is its oscillation frequency?

### Step 4: Measure m_eff from confinement

If the θ packet is trapped:
- The oscillation frequency ω_electron gives the energy: E = ℏω
- The confinement momentum: p ~ ℏ/σ_final
- The effective mass: m_eff = p/v = ℏk/v where v is the packet velocity

Or more directly: m_eff = E_θ_packet / c² where E_θ_packet is the
total energy of the trapped θ perturbation.

## Implementation

Use the EXISTING v33_cosserat code (no θ self-interaction needed).
Modify init to place a θ wave packet at (r_0, 0, 0) relative to braid.

The key modification in init_braid or a new init_theta_packet function:

```c
// After standard φ-braid init, add θ wave packet
for (ix...) {
    double r = sqrt((x-theta_x)² + (y-theta_y)²);
    double env = A_theta * exp(-r²/(2*sigma²));
    for (int a = 0; a < 3; a++) {
        g->theta[a][idx] += env * cos(k*z + 2*PI*a/3);
        g->theta_vel[a][idx] += omega * env * sin(k*z + 2*PI*a/3);
    }
}
```

## Experiments

### 5a: Packet at various distances (does it trap?)

Place θ packet at r_0 = {5, 8, 12, 16, 20} from the braid center.
Track the θ centroid over time. Does it stay at r_0, spiral in, or fly away?

N=80, L=25, T=200, η=0.5, m_θ=0, A_θ=0.05, σ=2

### 5b: Packet amplitude scan

At the best r_0 from 5a, vary A_θ = {0.01, 0.02, 0.05, 0.1, 0.2}.
Does amplitude affect trapping? (It shouldn't for a linear field,
but the curl coupling introduces nonlinearity.)

### 5c: Long-term stability

At the best (r_0, A_θ), run T=1000. Does the packet:
- Maintain coherence? (localization measure)
- Settle into a standing wave? (frequency spectrum)
- Radiate away gradually?

### 5d: Frequency measurement

FFT the trapped θ packet's position and amplitude over time.
The oscillation frequency → energy level.
If multiple frequencies appear → multiple orbital modes.

## Expected Outcome

If the φ-braid's potential well traps the θ packet:
- m_eff = E_packet / c²
- Combined with ℏ_eff = E_braid/ω ≈ 23,000
- Bohr ratio = ℏ²/(m_eff × α_eff × r_braid)
- If this matches ~53,000: electron from trapped θ wave packet

If the θ packet disperses: the potential well is too shallow or
the massless θ field propagates too fast to be trapped. Would need
a small m_θ² > 0 to slow it down.
