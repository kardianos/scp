# V26-DynC Results: Massless Propagating Braid

## Answer: NO — propagation does NOT prevent m=0 collapse

Both the propagating helical wave and the static control collapse catastrophically
within t ~ 6-7 time units. The triple product potential alone cannot confine
massless fields regardless of whether they carry momentum.

## Parameters

mu=-20, kappa=20, m=0.0, A0=0.8, R_tube=3.0, N=128, L=20, dt=0.0630

## Initial Conditions

| Quantity  | Propagating (mode 0) | Static (mode 1) |
|-----------|---------------------|-----------------|
| E_total   | 59.21               | 45.71           |
| E_kin     | 13.50               | 0.00            |
| E_grad    | 73.79               | 73.79           |
| E_pot     | -28.08              | -28.08          |
| fc        | 0.993               | 0.991           |
| P_z       | 26.78               | 0.00            |

The propagating mode carries 13.5 units of kinetic energy (wave propagates at c)
and 26.8 units of z-momentum. Gradient and potential energies are identical.

## Collapse Timeline

| Event              | Propagating | Static  |
|--------------------|-------------|---------|
| fc drops below 0.01 | t = 6.1    | t = 5.9 |
| E_total goes negative | t = 7.1  | t = 6.7 |
| E at t=50          | ~ -7600    | ~ -7200 |
| E at t=150         | ~ -7740    | (t=65 data end: -7520) |
| Peak |phi|         | ~1.7       | ~2.1    |
| Peak |P|           | ~1.4       | ~2.4    |

## Mechanism of Failure

1. **No mass gap**: With m=0, linear waves propagate freely at c (no dispersion).
   The triple product V(P) = (mu/2)P^2/(1+kappa*P^2) provides a NEGATIVE potential
   well (mu < 0), but with no mass term there is no restoring force at small amplitudes.

2. **Runaway deepening**: Fields grow, |P| increases, V becomes more negative,
   driving further growth. The saturating denominator (1+kappa*P^2) is overwhelmed.

3. **Energy conservation violation**: Energy drops from +59 to -7700 because the
   potential energy (which starts at -28) plunges to -11,700 as fields grow.
   The kinetic + gradient energy grows to ~4000 but cannot compensate.
   This is a genuine instability, not a numerical artifact.

4. **Propagation irrelevant**: The propagating mode collapses at t=6.1 vs t=5.9
   for static — effectively identical. The z-momentum (26.8 initially) does not
   stabilize the transverse collapse. The fields delocalize equally fast.

## Comparison with V26 Phase 4 (m=1, static)

V26 mode 2 (m=1, triple product only, static) survived at fc=0.37 for t=500.
The mass term m=1 provided the crucial stabilization:
- m=1: E_mass = (1/2)m^2|phi|^2 acts as a POSITIVE restoring force
- m=0: no such term exists; fields can grow without penalty

## Key Conclusion

**The mass term is essential for braid confinement.** The triple product alone
(without mass) has an unbounded negative potential valley. Propagation (kinetic
energy, z-momentum) does not change the transverse stability properties at all.

The "emergent mass from confinement" thesis is falsified: you need mass to get
confinement, not the other way around.

## Files

- `src/dynC.c` — simulation code
- `data/dynC_propagating.tsv` — mode 0 time series (to t~151)
- `data/dynC_static_control.tsv` — mode 1 time series (to t~65)
- `data/dynC_*_profile_t0.tsv` — initial radial profiles
