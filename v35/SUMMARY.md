# V35 Summary: Complete

## What was accomplished

### Parameter Space Mapping
- 12,000-point sweep found exact hydrogen match (Bohr ratio 52,957)
- One-parameter family: ℏ²/m_eff = const
- Best ℏ candidate: E_braid/ω ≈ 22,727

### Negative Results (all informative)
- Field quantum ε: substrate too fine for atomic ℏ (6000× gap)
- θ self-interaction: uniform condensate, no solitons (no mass gap)
- θ frequency spectrum: continuous/dispersive (no natural quantization)
- Trapped θ packet: classical collapse without ℏ (Rutherford problem)

### Multi-Scale Architecture
- Built: 3D spherical core + 1D Schrödinger wedge (760 lines C)
- Spherical absorbing BC: braid survives, θ reaches R_match
- Standalone wedge: stable eigenstate over 44 orbital periods
- CN solver: norm conserved to 2×10⁻⁷

### Infrastructure
- SFA format: streaming mode with pre-allocated JMPF slots
- SFA compression study: BSS+zstd-3 optimal, 7 codecs benchmarked
- GPU rental: Vast.ai at $0.20/hr, OVH at $0.88/hr
- Multi-scale architecture document for V100 deployment

## Key Insight

The theory has TWO levels:
1. **Classical field** (3D Cosserat): provides matter, gravity, EM, and V_eff(r)
2. **Quantum mechanics** (from braid action ℏ=E/ω): provides orbital confinement

The electron is a stationary θ shell at the Bohr radius, stabilized by ℏ.
The classical simulation provides the potential; quantum mechanics prevents
collapse. Dynamics require perturbations (other atoms, external fields).

## Open for V36

1. CUDA port of force computation
2. Multi-atom simulation (H₂ molecule)
3. Orbital transitions (photon emission)
4. Full 3D electron visualization
