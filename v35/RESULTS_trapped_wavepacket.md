# Computation 5 Results: Trapped Theta Wave Packet

## Summary

**The theta wave packet is NOT trapped by the phi-braid's potential well.**

A massless theta wave packet disperses at speed c and fills the periodic box
within ~25 time units (L/c). Even with theta mass (m_t^2 up to 4.0), the
packet disperses completely by T=200. The phi-braid's curl coupling acts as
an antenna/source for theta waves, not a confining potential. There is no
electron analog in this classical setup.

## Code

- Source: `src/v33_trapped.c` (adapted from `v33_theta_self.c`)
- Build: `gcc -O3 -march=native -fopenmp -o v33_trapped src/v33_trapped.c -lm`
- Key modifications: removed theta self-interaction, added Gaussian wave
  packet init with command-line flags (-A_theta, -theta_x, -theta_y,
  -sigma_theta), added theta^2-weighted centroid tracking diagnostic

## Experiment 5a: Distance Scan (Massless Theta)

N=80, L=25, T=200, eta=0.5, A_theta=0.05, sigma=2.0, m_theta^2=0

| r_0 | pkt_r(T=200) | spread(T=200) | pkt_E(T=200) | peak(T=200) |
|-----|-------------|--------------|-------------|------------|
|   5 |    0.26     |    21.21     |    49.5     |   0.232    |
|   8 |    0.24     |    21.24     |    49.3     |   0.235    |
|  12 |    0.20     |    21.27     |    48.7     |   0.239    |
|  16 |    0.17     |    21.26     |    48.9     |   0.238    |
|  20 |    0.18     |    21.23     |    48.9     |   0.236    |

**Result**: All five initial positions converge to the SAME final state by
t~2. The packet disperses at speed c and fills the box (spread ~21.2, close
to the theoretical L*sqrt(2/3)=20.4 for uniform distribution in [-25,25]^3).

Key observation: at t=0.95, the centroid has already moved from r_0 to
r~0.9-4.4 (with initial position becoming irrelevant by t~2). The theta^2
centroid then tracks the braid's own curl-induced theta halo, not the
injected packet. The pkt_E jumps from 0.32 to ~30-50, reflecting energy
from the braid coupling, not the original packet.

## Experiment 5b: Mass Scan

N=80, L=25, T=200, eta=0.5, A_theta=0.05, sigma=2.0, r_0=8

| m_t^2 | spread(t=10) | spread(t=50) | spread(T=200) | pkt_E(T=200) | peak(T=200) |
|-------|-------------|-------------|--------------|-------------|------------|
|  0    |   16.0      |   20.9      |    21.2      |    49.3     |   0.235    |
|  0.1  |    7.4      |   18.9      |    19.7      |    66.6     |   0.300    |
|  0.5  |    9.8      |   18.5      |    19.6      |    85.8     |   0.356    |
|  1.0  |    5.6      |   15.8      |    20.2      |   114.6     |   0.327    |
|  4.0  |    8.2      |   13.0      |    17.8      |    20.8     |   0.095    |

Also tested r_0=3 with m_t^2=1.0: final spread=20.2, same as r_0=8.

**Result**: Mass slows the initial dispersion but does NOT prevent it.

The m_t^2=1.0 case has a transient "compression" phase where spread=5.6 at
t=10 (vs 16.0 for massless), but by t=50 the spread reaches 15.8, and by
T=200 it reaches 20.2.

Heavy mass (m_t^2=4.0) actually SUPPRESSES the theta field: the peak drops
to 0.095 (vs 0.235 massless) because the heavy theta responds less to the
curl coupling source. The tradeoff is clear: lighter theta disperses faster
but couples more strongly; heavier theta couples less but disperses slower.
Neither regime produces trapping.

## Experiment 5c: Long-Duration Run

N=80, L=25, T=500, eta=0.5, r_0=8, m_t^2=1.0

| t    | spread | pkt_E | peak  |
|------|--------|-------|-------|
|    0 |   2.0  |   2.7 | 0.060 |
|   48 |  17.0  | 124.9 | 0.317 |
|   95 |  17.8  | 138.0 | 0.430 |
|  142 |  19.1  | 187.2 | 0.540 |
|  190 |  20.2  | 105.0 | 0.276 |
|  237 |  20.1  | 113.5 | 0.536 |
|  285 |  20.7  | 147.3 | 0.489 |
|  332 |  20.3  | 115.4 | 0.288 |
|  380 |  20.2  | 126.9 | 0.463 |
|  427 |  20.0  | 147.5 | 0.448 |
|  475 |  21.4  | 112.7 | 0.404 |

The spread reaches ~18 by t~50 and then saturates at ~20 for the remainder
of T=500. No late-time trapping emerges. The theta peak oscillates
(0.28-0.54) due to energy exchange with the braid coupling, but this is a
volume-filling standing wave pattern, not a localized trapped packet.

## Energy Conservation

Energy drift over T=200: approximately -0.5% to -0.7% across all runs.
Symplectic Verlet integrator, dt = 0.10 * dx.

## Physical Interpretation

### Why the packet disperses

The theta wave equation is:
```
d^2 theta_a/dt^2 = lap(theta_a) - m_t^2 theta_a + eta * curl(phi)_a
```

This is a LINEAR wave equation with a SOURCE term (eta * curl(phi)). The
source drives theta near the braid core, but the theta field propagates
freely as a wave. There is no potential well in the Schrodinger sense.

The curl coupling acts as an ANTENNA: it excites theta waves near the braid,
which then radiate outward at the group velocity v_g = k/omega. For massless
theta, v_g = c everywhere. The "potential well" posited in the plan does
not exist -- curl(phi) is a source, not a confining potential.

### Why mass does not help enough

A massive theta field has dispersion omega = sqrt(k^2 + m^2). The group
velocity v_g = k/sqrt(k^2 + m^2) < c is reduced but still nonzero for any
k > 0. The packet spreads at the group velocity, which only vanishes in the
k -> 0 limit (uniform field).

### What would be needed for trapping

True trapping would require one of:
1. A position-dependent effective mass: m_eff^2(r) that creates a genuine
   bound state in the wave equation (like a delta function or square well
   potential), not just a source term
2. Nonlinear self-interaction in theta (ruled out in Computation 4 --
   triple product coupling does not create localized solitons)
3. A coupling of the form V(r) * theta that acts as a multiplicative
   potential rather than an additive source, e.g., from theta propagating
   on a curved effective metric created by the phi-braid

The eta * curl(phi) term generates theta waves. It does NOT create a
restoring force proportional to theta's displacement from the braid.

## Conclusion

The theta wave packet experiment definitively rules out the "particle in a
box" model for the electron in the current Cosserat framework. The phi-braid's
curl coupling radiates theta waves but does not trap them.

**This closes Computation 5 with a negative result.**

## Files

- `src/v33_trapped.c` -- simulation code
- `data/trapped/r{05,08,12,16,20}/` -- Experiment 5a distance scan
- `data/trapped/r08_mt{01,05,10,40}/` -- Experiment 5b mass scan
- `data/trapped/r03_mt10/` -- r=3 with mass control
- `data/trapped/r08_mt10_long/` -- Experiment 5c long-duration run
