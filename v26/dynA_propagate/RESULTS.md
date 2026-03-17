# V26-DynA: Propagating Helical Wave Braid — RESULTS

## Parameters

| Parameter | Value |
|-----------|-------|
| mu | -20.0 |
| kappa | 20.0 |
| mass | 1.0 |
| A0 | 0.8 |
| R_tube | 3.0 |
| N | 128 |
| L | 20.0 |
| dx | 0.3150 |
| dt | 0.06299 |
| tfinal | 500 |
| k | 0.1571 |
| omega | 1.0123 |
| v_g = k/omega | 0.1552 |
| BC | periodic z, absorbing x,y |
| Couplings | triple product ONLY (lambda_pw=0, eta=0, lambda_L=0) |

## Summary Comparison

| Quantity | Propagating | Static Control |
|----------|------------|----------------|
| E(0) | 1153.5 | 592.9 |
| E(500) | 281.3 | 208.1 |
| E retained | 24.4% | 35.1% |
| fc(0) | 0.363 | 0.360 |
| fc(500) | 0.276 | 0.268 |
| Survived? | YES | YES |
| |P|_max final | 0.574 | 0.006 |
| Pz(0) | 172.6 | 0.0 |
| Pz(500) | 35.8 | 0.0 |
| Pz retained | 20.8% | N/A |
| Breathing? | NO | NO |
| l=0 fraction | 0.414 | 0.910 |
| l=1 fraction | 0.320 | 0.022 |
| l=2 fraction | 0.267 | 0.068 |

## Key Findings

### 1. Core fraction (fc): Propagating does NOT increase fc

Both runs end at fc ~ 0.27. The propagating braid has higher fc at intermediate
times (peaks at 0.446 near t=349) vs static (peaks at 0.475 near t=116), but
both converge to similar final values. **No clear advantage from propagation.**

### 2. Triple product |P|: Propagating maintains it FAR better

- Propagating: |P|_max = 0.574 at t=500 (strong oscillations throughout)
- Static: |P|_max = 0.006 at t=500 (collapses to ~0 by t~300)

The initial velocity sustains the nonlinear braid coupling. The static braid
loses its triple-product structure entirely, becoming three decoupled fields.
The propagating braid maintains large |P| oscillations (0.01-0.93) throughout.

### 3. Momentum Pz: Slow decay, not conserved

Pz decays from 172.6 to 35.8 (79% loss over t=500). The decay is continuous
and roughly exponential — energy/momentum radiate outward through the absorbing
x,y boundaries. This is expected: the tube radiates laterally.

### 4. l=2 content: MUCH higher for propagating

- Propagating: l=2 fraction = 26.7% (with l=1 at 32.0%)
- Static: l=2 fraction = 6.8% (with l=0 dominating at 91.1%)

The propagating braid breaks spherical symmetry much more strongly. The large
l=1 component reflects the z-propagation breaking parity. The static braid
is nearly isotropic (l=0 dominant) because the twist averages out.

However, the l=2 fraction of 26.7% is LOWER than V26 Phase 4's 41.5%.
The propagation redistributes power into l=1 (dipolar) rather than concentrating
it in l=2 (quadrupolar).

### 5. Breathing: Neither breathes

Both configurations show NO breathing mode (relative variance of rho_center
is sub-threshold). The propagating braid has omega_phi = 0.94, close to the
mass frequency omega = 1.01, consistent with the field oscillating at the
dispersion frequency rather than a breathing resonance.

### 6. Energy budget

The propagating braid starts with nearly 2x the energy (1153 vs 593) due to
the initial kinetic energy from vel. Both lose energy to radiation through
the absorbing boundaries, but the propagating case loses a larger fraction
(75.6% vs 64.9%) — the momentum drives more lateral radiation.

## Conclusion

Propagation along z via initial velocity produces:
- **Sustained triple-product coupling** (|P| remains O(1) vs collapsing to ~0)
- **Significant l=1 dipolar strain** (32%) from the broken z-symmetry
- **Moderate l=2 quadrupolar strain** (26.7%), but lower than the static V26 result (41.5%)
- **No improvement in core fraction** (fc ~ 0.27 for both)

The propagating braid is a qualitatively different object — it maintains
its nonlinear braid structure indefinitely, while the static braid decays
to three uncoupled fields. However, the propagation does not enhance the
quadrupolar strain that was the key result of V26.
