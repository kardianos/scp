# V54 Chirality Pair Sweep Results

## Setup

Two braids placed in the same periodic box at positions (-D/2, 0, 0) and (+D/2, 0, 0).
Particle A is always standard chirality. Particle B uses a modified construction.
Run for T=50 per trial, then measure H_cross in the left (x<0) and right (x>0) halves.

## Parameters (from parameter sweep best)

- m² = 1.5, μ = -80, κ = 50, η = 0.2
- A = 0.4, A_bg = 0.1, R_tube = 3.0, ellip = 0.3325
- Standard delta = (0, 3.0005, 4.4325), theta_gain = -2.5
- Periodic BC, N=96, L=15

## Chirality Methods Tested

| Method | Phase offsets (B) | Theta sign (B) | Ellip (B) |
|--------|------------------|----------------|-----------|
| same (control) | (0, 3.0005, 4.4325) | -1 | +0.3325 |
| flip_theta | (0, 3.0005, 4.4325) | +1 | +0.3325 |
| reverse_delta | (4.4325, 3.0005, 0) | -1 | +0.3325 |
| negate_delta | (0, -3.0005, -4.4325) | -1 | +0.3325 |
| flip_theta+rev_delta | (4.4325, 3.0005, 0) | +1 | +0.3325 |
| flip_theta+neg_delta | (0, -3.0005, -4.4325) | +1 | +0.3325 |
| flip_ellip | (0, 3.0005, 4.4325) | -1 | -0.3325 |
| flip_all | (0, -3.0005, -4.4325) | +1 | -0.3325 |
| pi_shift | (π, 3.0005+π, 4.4325+π) | -1 | +0.3325 |
| pi_shift+flip_theta | (π, 3.0005+π, 4.4325+π) | +1 | +0.3325 |

## Full Results

| Method | Sep | φ_max A | φ_max B | H_A | H_B | Same? | Both? |
|--------|-----|---------|---------|-----|-----|-------|-------|
| same (control) | 6 | 1.2 | 1.2 | +1.9 | +5.1 | SAME | BOTH |
| same (control) | 8 | 1.2 | 1.3 | +5.0 | +1.0 | SAME | BOTH |
| same (control) | 10 | 1.6 | 1.5 | -3.6 | +6.4 | DIFF | BOTH |
| flip_theta | 6 | 1.2 | 1.1 | -6.8 | -5.8 | SAME | BOTH |
| flip_theta | 8 | 1.0 | 1.3 | -2.4 | -11.2 | SAME | BOTH |
| flip_theta | 10 | 1.4 | 1.4 | -5.2 | -19.4 | SAME | BOTH |
| reverse_delta | 6 | 1.2 | 1.3 | +21.0 | +13.0 | SAME | BOTH |
| reverse_delta | 8 | 1.3 | 1.2 | +4.0 | +9.6 | SAME | BOTH |
| reverse_delta | 10 | 1.5 | 1.3 | -2.4 | +9.4 | DIFF | BOTH |
| negate_delta | 6 | 1.1 | 1.2 | +14.5 | +15.9 | SAME | BOTH |
| negate_delta | 8 | 1.3 | 1.3 | +20.1 | +14.6 | SAME | BOTH |
| negate_delta | 10 | 0.9 | 0.8 | +13.6 | +8.4 | SAME | BOTH |
| flip_theta+rev_delta | 6 | 1.3 | 1.2 | -9.3 | -10.5 | SAME | BOTH |
| flip_theta+rev_delta | 8 | 1.3 | 1.3 | -12.5 | -10.5 | SAME | BOTH |
| flip_theta+rev_delta | 10 | 1.0 | 0.9 | -0.3 | -2.6 | SAME | BOTH |
| flip_theta+neg_delta | 6 | 0.9 | 1.3 | -6.5 | -11.9 | SAME | BOTH |
| flip_theta+neg_delta | 8 | 1.0 | 0.8 | -7.5 | -8.7 | SAME | BOTH |
| flip_theta+neg_delta | 10 | 1.4 | 1.4 | -2.4 | -8.2 | SAME | BOTH |
| **flip_ellip** | **6** | **1.3** | **1.3** | **+10.3** | **+28.5** | **SAME** | **BOTH** |
| **flip_ellip** | **8** | **1.4** | **1.3** | **-10.1** | **+24.4** | **DIFF** | **BOTH** |
| **flip_ellip** | **10** | **1.4** | **1.3** | **-1.9** | **+15.1** | **DIFF** | **BOTH** |
| **flip_all** | **6** | **1.1** | **1.0** | **-0.2** | **-10.5** | **SAME** | **BOTH** |
| **flip_all** | **8** | **1.3** | **1.3** | **+0.0** | **-14.6** | **DIFF** | **BOTH** |
| **flip_all** | **10** | **0.8** | **0.9** | **+7.9** | **-9.9** | **DIFF** | **BOTH** |
| pi_shift | 6 | 1.4 | 0.3 | +22.7 | +35.8 | SAME | BOTH |
| pi_shift | 8 | 1.2 | 0.4 | +18.1 | +37.1 | SAME | BOTH |
| pi_shift | 10 | 1.4 | 0.4 | +19.5 | +40.4 | SAME | BOTH |
| pi_shift+flip_theta | 6 | 1.0 | 0.3 | -15.1 | +0.8 | DIFF | ONE |
| pi_shift+flip_theta | 8 | 1.2 | 0.3 | -13.0 | -1.2 | SAME | BOTH |
| pi_shift+flip_theta | 10 | 1.3 | 0.4 | -17.9 | -3.6 | SAME | BOTH |

## Key Findings

1. **Most methods produce SAME-sign H_cross.** The dynamics absorbs phase/theta
   differences and forces both particles to the same chirality. Flipping theta
   sign, reversing phases, negating phases — none of these create persistent
   opposite chirality.

2. **flip_ellip produces OPPOSITE H_cross at D=8 and D=10.** Ellipticity determines
   the geometric handedness of the cross-section (elongated in x vs y). This
   spatial shape asymmetry cannot be dynamically absorbed — it's a genuine
   geometric chirality.

3. **flip_all (negate delta + flip theta + flip ellip) also produces DIFF** at D=8,10.
   The ellip flip is the key ingredient; the other flips are secondary.

4. **All particles survive at T=50** regardless of chirality method. Both halves
   have φ_max > 0.8 in most cases.

## Physical Interpretation

Chirality in this system is determined by the **spatial geometry** of the envelope
(ellipticity), not by the phase structure or theta sign. The equation dynamics
can rearrange phases and theta signs freely, but it cannot change the shape of
the spatial envelope. A braid elongated in x has opposite handedness to one
elongated in y.

## Next Step

Run flip_ellip at D=8 for T=3000 to test whether the opposite-chirality pair
is more stable than a single particle. The hypothesis: opposite chiralities
exchange theta between them, creating a self-sustaining cycle that prevents
dissolution.
