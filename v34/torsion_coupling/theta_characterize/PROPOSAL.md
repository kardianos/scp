# θ Field Characterization: Numerical Confirmation of EM Structure

Three experiments to determine if the θ field behaves as electromagnetism.

## Experiment 1: Radial Decay of θ (Biot-Savart test)

Does θ_circular fall as 1/r from the braid axis?

**Method**: Read a settled frame from sfa_hires.sfa. Decompose θ into
cylindrical components around the braid's z-axis:
- θ_r (radial), θ_φ (azimuthal/circular), θ_z (axial)
Average θ_φ² in cylindrical shells at distance r from the braid axis.
Fit to power law: θ_φ ∝ 1/r^n. If n≈1: Biot-Savart confirmed.

## Experiment 2: Winding Reversal (charge conjugation)

Does W=-1 produce opposite θ circulation?

**Method**: Run two simulations with opposite winding:
- W=+1: δ = (0, +3.0005, +4.4325)
- W=-1: δ = (0, -3.0005, -4.4325)
Same η=0.5, m_θ=0, N=80, L=25, T=50. Compare the θ_φ component.
If θ_φ(W=+1) = -θ_φ(W=-1): winding IS charge.

## Experiment 3: Two Parallel Braids (current attraction)

Do same-winding braids attract through θ?

**Method**: Two braids at D=15, both W=+1, in the 6-field equation
(η=0.5, m_θ=0). Compare ΔD to the 3-field baseline (η=0, same setup).
If ΔD is MORE negative with η>0: the θ field adds an attractive force
between same-winding braids = parallel currents attract.

Also run: two braids, opposite winding (W=+1 and W=-1). If ΔD is
LESS negative (or positive): opposite "currents" repel through θ.
