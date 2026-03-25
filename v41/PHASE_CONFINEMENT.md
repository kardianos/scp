# V41 Phase Confinement Experiment

## Hypothesis

Three braids with carrier phase offsets of 0, 2π/3, 4π/3 will resist merging
because their overlap creates destructive V(P) interference at short range,
while the depletion interaction (gravity) still attracts at medium range.
This produces CONFINEMENT: attracted but unable to merge.

The phase offset plays the role of COLOR CHARGE in QCD. Three equally-spaced
phases are the "color neutral" configuration — the only stable composite state.

## Physics of Phase Confinement

### Why braids merge without phase offset

Currently, all three braids use the same phase relationship δ = {0, 3.0, 4.4}.
When two braids overlap spatially, their contributions to P = φ₀φ₁φ₂ add
CONSTRUCTIVELY — the merged state has higher |P| than the separated state.
Since V(P) < 0, higher |P| means lower energy. Merging is energetically
favored. Nothing prevents it.

### How phase offset prevents merging

With per-braid carrier phase offsets Δ = {0, 2π/3, 4π/3}:

Each braid's fields become:
```
Braid k: φ_a(x,t) = A(x) × cos(k_z × z_rotated + δ_a + Δ_k)
```

When braids k and k' overlap at the same spatial point:
```
P = (φ₀_k + φ₀_k') × (φ₁_k + φ₁_k') × (φ₂_k + φ₂_k')
```

The cross-terms involve cos(... + Δ_k) × cos(... + Δ_k'). With Δ_k - Δ_k' = ±2π/3,
these cross-terms average to cos(2π/3) = -1/2. The overlap REDUCES |P| compared
to same-phase braids.

At the center where all three overlap:
```
P_total ∝ Σ_k cos(Δ_k) = cos(0) + cos(2π/3) + cos(4π/3) = 1 - 1/2 - 1/2 = 0
```

The triple overlap has P = 0! No binding at the center point. This is the
CONFINEMENT mechanism: three phase-offset braids cannot merge because their
combined P vanishes at the overlap center.

But at MEDIUM range (separation ~3-5 code units), each braid's depletion profile
still creates an attractive gradient for the others. The braids are attracted
but confined to their equilibrium separation.

### Analogy to QCD color charge

| QCD | Cosserat with phase offset |
|-----|---------------------------|
| Color charge (R, G, B) | Carrier phase (0, 2π/3, 4π/3) |
| Color neutral = R+G+B = white | Phase neutral = sum of phases cancels P |
| Gluon exchange (strong force) | V(P) coupling (binding potential) |
| Confinement (can't isolate quarks) | Phase repulsion at short range |
| Asymptotic freedom (weak at short r) | P → 0 when all three overlap |
| 1/3 fractional charge | 1/3 of the total phase circle (2π/3 each) |

### Fractional charge (1/3 θ per braid)

Each braid sources θ through η × curl(φ). The curl of φ depends on the braid's
spatial gradients and chirality. With phase offset Δ_k, each braid's θ contribution
is phase-shifted by Δ_k:

```
θ sourced by braid k ∝ curl(φ_k) ∝ A × k_z × sin(k_z × z + Δ_k)
```

The NET θ from all three braids:
```
θ_total ∝ sin(Δ_0) + sin(Δ_1) + sin(Δ_2) = sin(0) + sin(2π/3) + sin(4π/3) = 0
```

**Three phase-offset braids produce ZERO net theta.** Each braid carries 1/3 of
the theta circle, and the vector sum cancels. This is exactly the 1/3 fractional
charge of quarks — each quark carries 1/3 of the color circle, and the baryon
(R+G+B) is color-neutral.

## Experiment Design

### UDD Phase-Confined (Neutron Analog)

Three braids along x, y, z axes with:
- Braid 1 (z-axis): chirality U (+1), carrier phase Δ₁ = 0
- Braid 2 (x-axis): chirality D (-1), carrier phase Δ₂ = 2π/3
- Braid 3 (y-axis): chirality D (-1), carrier phase Δ₃ = 4π/3

Net charge: U + D + D = +1 -1 -1 = -1
Net θ (from phase cancellation): **0** (color neutral)

**Prediction**: The UDD with zero net θ should be LESS stable than UUD because
it has no net electromagnetic coupling to maintain coherence. The three braids
may oscillate around equilibrium positions initially, but without net θ coupling,
the binding relies entirely on V(P) depletion. Expect:
- Initial confinement (braids maintain separation due to phase repulsion)
- Slow decay via θ emission (the phase-cancelled θ field can still produce
  radiation at the sum/difference frequencies)
- The decay products would be a theta blob (photon analog) + remaining structure
- This IS beta-decay behavior: neutron → proton + electron + neutrino

### UUD Phase-Confined (Proton Analog)

Three braids along x, y, z axes with:
- Braid 1 (z-axis): chirality U (+1), carrier phase Δ₁ = 0
- Braid 2 (x-axis): chirality U (+1), carrier phase Δ₂ = 2π/3
- Braid 3 (y-axis): chirality D (-1), carrier phase Δ₃ = 4π/3

Net charge: U + U + D = +1 +1 -1 = +1
Net θ: NOT zero (two same-chirality braids don't cancel)

**Prediction**: The UUD with net θ ≠ 0 should be MORE stable because the net
electromagnetic coupling provides an additional binding channel. The proton
analog should:
- Maintain confinement (phase repulsion prevents merging)
- Retain net θ coupling (stabilizing)
- NOT decay via beta process (net charge is preserved)
- Be the stable baryon

### Parameters

- N = 192, L = 25, T = 200
- A = 0.3 (near P_opt), R_tube = 4.5 (wide for overlap)
- Pre-loaded θ (theta_init = 0.5)
- Contracting velocity profile (from V41 stability signatures)
- Carrier phase offsets: Δ = {0, 2π/3, 4π/3}

## Success Criteria

1. **Confinement**: The three braids maintain distinct centroids (>2 code units
   apart) for T > 100. They do NOT merge into a single blob.

2. **Stability ordering**: UUD (proton) survives longer than UDD (neutron).
   Ideally UUD survives T=200, UDD decays before T=200.

3. **Phase preservation**: The carrier phase offsets remain measurable at T=100
   (the phases don't randomize).

4. **Theta emission from UDD**: The UDD neutron analog should emit θ radiation
   that the UUD proton analog does not. This would be the beta decay signature.

## If This Works → V42: Nuclear Binding

If the phase confinement experiment confirms:
- Three-body confinement via phase offset
- UUD stable, UDD unstable (or less stable)
- Distinct braid identities maintained

Then V42 targets nuclear binding — multiple confined baryons bound by
the residual depletion interaction:

| Nucleus | Structure | N_grid | Notes |
|---------|-----------|--------|-------|
| ²H (Deuterium) | UUD + UDD | 384-512 | Simplest nucleus, 1p + 1n |
| ³He (Helium-3) | UUD + UUD + UDD | 512 | 2p + 1n |
| ⁴He (Helium-4) | UUD + UUD + UDD + UDD | 512+ | 2p + 2n, very stable (alpha) |

The inter-baryon binding comes from the RESIDUAL depletion field — the
equivalent of nuclear force (pion exchange in QCD). Each confined baryon
depletes ρ in its neighborhood, and another baryon is attracted to the
depletion zone. The multi-resolution simulation framework (scp_multi)
would be essential here: coarse grid for inter-baryon spacing, fine grids
for internal braid structure.

These simulations require N=512 (each baryon needs ~100 code units diameter
plus separation), which means Tesla_V100 32 GB or multi-GPU. Estimated cost:
~$1-5 per nucleus at T=200 on a 4×V100 instance.
