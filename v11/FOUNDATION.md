# V11: Born-Infeld Magnetic Hopfion

## Question

Does Born-Infeld nonlinearity stabilize electromagnetic configurations
with Hopf topology in flat space?

## Answer: NO

BI nonlinearity provides transient slowing of dispersal but does NOT
prevent it. The hopfion disperses at speed c regardless of b.

## Method

- Pure magnetic initial data from Hopf map pullback (area form of S3 -> S2)
- Born-Infeld electrodynamics with field strength parameter b
- 3D Cartesian grid, 4th-order spatial derivatives, RK4 time integration
- State variables: (D, B), with E from Newton inversion of D(E,B) = D_target
- Sponge absorbing boundary layer

## Key Physics

The standard Ranada hopfion is a NULL field (E^2 = B^2, E.B = 0), so BI
corrections vanish identically on it. We use a PURE MAGNETIC hopfion (E=0)
with B from the Hopf map area form. This has S = B^2/2 > 0, so BI
nonlinearity is active from t=0.

Initial B normalized so peak |B| = b_param (default 1.0, independent of BI
parameter b). This means C0 = b_param * a^2 / 16 with peak at origin.

## Numerical Scheme

Staggered leapfrog (Yee time staggering) is unconditionally unstable on a
collocated spatial grid (proved by von Neumann analysis: growth factor
|g| > 1 for all k, any dt > 0). Fixed by switching to RK4 method-of-lines:

State: u = (D, B), 6 components per point
RHS: f(u) = (curl(H(E(D,B),B)), -curl(E(D,B)))

Each RHS evaluation requires Newton inversion (D,B) -> E at every grid point,
then forward constitutive (E,B) -> (D,H), then two curl operations. RK4
requires 4 evaluations per step.

CFL: dt = 0.4 h / sqrt(3)

## Energy Formula

The naive u = (b^2 + B^2)/Gamma - b^2 suffers catastrophic cancellation for
large b (b^2 + O(1) - b^2 = 0 in floating point). Cancellation-free form:

u = (B^2(1+B^2/b^2) + E^2 + (E.B)^2/b^2) / (Gamma * (1+B^2/b^2 + Gamma))

Reduces to (B^2+E^2)/2 for b -> infinity, and to b^2(sqrt(1+B^2/b^2)-1) for E=0.

## Results

### B-scan (N=64, L=6, T=10, a=1, initial B_peak ~ 1.0)

| b | R_eff(t=1) | B_peak(t=1) | R_eff(t=5) | B_peak(t=5) |
|---|-----------|-------------|-----------|-------------|
| 1.0 | 1.365 | 0.454 | 4.901 | 0.046 |
| 2.0 | 1.391 | 0.301 | 4.975 | 0.050 |
| 5.0 | 1.399 | 0.284 | 4.994 | 0.051 |
| 10.0 | 1.400 | 0.282 | 4.997 | 0.051 |
| 1e10 | 1.401 | 0.281 | 4.998 | 0.051 |

At t=1 (early phase): b=1 retains 61% more B_peak than Maxwell (0.454 vs 0.281).
At t=5 (late phase): difference shrinks to <2%. All values of b give same dispersal.

### N=128 production (b=1.0, L=6, T=10)

Energy conservation: dE/E0 < 1e-4 before sponge absorption (t < 3)
div B: ~5e-3 (4th-order truncation of initial data, preserved by scheme)
B_peak: 0.987 -> 0.009 (full dispersal by t=10)
R_eff: 1.015 -> 6.92 (grows at ~c)

### Comparison at N=128

| Quantity | b=1 (BI) | b=1e10 (Maxwell) |
|----------|---------|-----------------|
| B_peak(t=1) | 0.43 | 0.27 |
| B_peak(t=5) | 0.045 | 0.051 |
| R_eff(t=5) | 5.01 | 5.04 |
| E(t=10)/E(0) | 1.1% | 0.9% |

BI nonlinearity modifies signal speed (slower for strong fields), causing
transient slowdown of the initial dispersal. Once |B| << b, corrections
become negligible and field disperses at c.

## Conclusion

Born-Infeld nonlinearity does NOT stabilize electromagnetic hopfions.
The topology alone cannot confine EM energy in flat space. This closes
Option B (Phase 2) from V9/NEXT.md.

The result is not surprising: BI is still a relativistic theory with
maximum signal speed c. Any initial configuration in flat space must
eventually disperse at c unless confined by external forces or coupling
to a massive field.

## Files

- `src/bi_hopfion.c` — full simulation code (~600 lines)
- `src/Makefile`
- `data/evolve_b1.0000.dat` — b=1 evolution timeseries
- `data/bscan.dat` — b-scan summary
- `data/maxwell_test.dat` — Maxwell validation

## What Would Be Needed for Stabilization

1. **Massive mediator**: coupling to a massive scalar/vector field
   (as in V4 DBI Q-balls, which DO self-trap)
2. **Gravity**: curved spacetime can confine (V8-V10 geons), but
   the coupling constant gap is ~10^40
3. **Confining potential**: explicit topology-dependent potential
   (not present in BI)
