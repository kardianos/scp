# V24-A Results: 2D Vortex Pair

## Parameters

| Parameter | Value |
|-----------|-------|
| mu | -20 |
| kappa | 20 |
| mass | 1.0 |
| Gamma | 1.0 |
| a_core | 2.0 |
| d_sep | 10.0 |
| phi3_noise | 0.01 |
| Grid | 512x512, L=40 (dx=0.156) |
| dt | 0.044 |
| t_final | 2000 |

## Key Result: VORTEX DOES NOT SURVIVE

The vortex-antivortex pair **disperses** in the massive three-field system. There is no
stable propagating vortex solution.

## Detailed Evolution

### Phase 1: Rapid Dispersal (t=0-100)

The initial vortex structure breaks apart almost immediately:

| t | peak |phi1| | peak |phi2| | E_total | peak vorticity |
|---|---------|---------|---------|----------------|
| 0 | 0.449 | 0.298 | 15.68 | 0.497 |
| 25 | 0.158 | 0.143 | 15.72 | 0.076 |
| 50 | 0.112 | 0.098 | 13.83 | 0.041 |
| 100 | 0.104 | 0.066 | 9.30 | 0.032 |

Peak amplitudes drop by 4x and peak vorticity by 15x within t=100.
Energy drops from 15.7 to 9.3 (41% absorbed by boundaries).

### Phase 2: Slow Radiation Decay (t=100-2000)

After the initial burst, residual waves slowly leak out:

| t | peak |phi1| | peak |phi2| | E_total | peak vorticity |
|---|---------|---------|---------|----------------|
| 200 | 0.062 | 0.043 | 5.27 | 0.019 |
| 500 | 0.064 | 0.050 | 2.84 | 0.013 |
| 1000 | 0.043 | 0.018 | 1.57 | 0.006 |
| 1500 | 0.023 | 0.022 | 1.03 | 0.004 |
| 2000 | 0.026 | 0.016 | 0.72 | 0.002 |

Energy at t=2000 is 0.72 (4.6% of initial), mostly kinetic radiation bouncing in the box.
Peak vorticity at t=2000 is 0.002 (0.4% of initial) — vortex structure is gone.

### Phase 3: phi3 Behavior

phi3 (initialized as noise at amplitude 0.01) stays small throughout:
- peak |phi3|: 0.010 (t=0) -> 0.019 (t=25, brief growth) -> 0.001 (t=2000)
- E_pot remains essentially zero: the triple product P = phi1*phi2*phi3 never
  becomes significant because phi3 never grows coherently.

The noise seed does NOT get amplified. The triple-product coupling has no effect
because the vortex disperses too quickly for P to build up.

## Why the Vortex Fails

1. **Mass term kills the tail**: A vortex requires long-range 1/r flow, but the mass
   term m^2*phi forces exponential decay. The initial Lamb-Oseen profile has 1/r tails
   which are immediately expelled as radiation (visible in the energy spike at early times).

2. **No topological protection**: In 2D with two scalar fields, there is no topological
   charge. The winding number of a vortex in (phi1, phi2) space is not conserved because
   |phi| is not constrained to a circle (unlike the XY model / Ginzburg-Landau).

3. **Triple product is impotent**: P = phi1*phi2*phi3 requires ALL THREE fields nonzero.
   With phi3 = noise, P ~ 0 and the coupling provides no binding force. Even if phi3
   were initialized large in the core, the mass term would disperse it independently.

4. **Divergence constraint not enforced**: The initial condition is divergence-free, but
   the equations of motion do not preserve this. Each field evolves as an independent
   massive Klein-Gordon equation (plus the tiny P-coupling), so divergence-free character
   is lost immediately.

## Comparison with Oscillon (v21)

| Property | Oscillon (v21) | Vortex (v24-A) |
|----------|----------------|----------------|
| Topology | None (breathing) | Winding (no protection) |
| Init profile | Gaussian (localized) | Lamb-Oseen 1/r (delocalized) |
| Mass term | Traps below gap | Kills 1/r tail |
| P coupling | phi1*phi2*phi3 binds | Requires phi3 coherence |
| Outcome | **Survives >2000** | **Disperses by t~50** |

The oscillon works because it is a localized, co-located Gaussian in all three fields
with large P from the start. The vortex fails because its spatial structure (circulating,
divergence-free) is incompatible with massive scalar fields.

## Vortex Centroid Tracking

The vorticity centroid wanders between y=10-25 (bouncing within the absorbing boundary
region), confirming there is no coherent propagation. A real vortex pair would translate
monotonically in x at v = Gamma/(2*pi*d).

Measured vortex_x stays at 0 throughout — NO x-propagation observed. This is consistent
with the pair dissolving rather than translating.

## Conclusion

**The vortex pair does not survive in the massive three-field system.** The mass gap
is fundamentally incompatible with the 1/r structure of a vortex. This is a definitive
negative result: the triple-product coupling cannot stabilize vortex-like configurations.

For vortex physics to work in this framework, one would need:
- A constraint |phi|^2 = const (sigma-model), removing the mass term
- Or a Mexican-hat potential V(|phi|^2) providing spontaneous symmetry breaking
- Or a gauge field (as in Abrikosov vortices) to supply the long-range interaction

None of these are present in the current Lagrangian.

## Output Files

- `data/vortex_ts.tsv` — 11361 time-series records (13 columns)
- `data/vortex_profile_t{T}.tsv` — 2D field snapshots at t=0,10,50,100,200,500,1000,2000
- `src/vortex2d.c` — full 2D solver (512x512 grid, Velocity Verlet, absorbing BC)
