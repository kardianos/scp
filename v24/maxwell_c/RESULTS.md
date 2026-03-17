# V24-MC Results: Topological Current Coupling

## Setup

2D three-field oscillon (phi_1, phi_2, phi_3) with triple-product coupling:

    L = sum_a [ (1/2)(dt phi_a)^2 - (1/2)|grad phi_a|^2 - (m^2/2)phi_a^2 ]
      - (mu/2) P^2 / (1 + kappa P^2),   P = phi_1 phi_2 phi_3

Parameters: mu=-20, kappa=20, m=1.0, A=0.8, sigma=3.0
Grid: Nx=Ny=256, L=20 (dx=dy=0.157), dt=0.044, Nt=45355

Topological charge density measured:

    rho_top = eps_{ij} eps_{abc} phi_c (d_i phi_a)(d_j phi_b)

## Key Result: Q_top = 0 Identically

The breathing oscillon carries **zero topological charge** at all times, to machine
precision (max |Q_top| = 0.0, max |rho_top| = 0.0).

This is not merely a numerical coincidence -- it is an exact algebraic identity.

## Why rho_top = 0 for the Symmetric Oscillon

The initial condition is phi_1 = phi_2 = phi_3 = A exp(-r^2/2sigma^2). All three
fields are identical. The equations of motion preserve this symmetry exactly: since
the Lagrangian is symmetric under permutation of (phi_1, phi_2, phi_3), if the
initial data has phi_1 = phi_2 = phi_3, this relation holds for all time.

With phi_1 = phi_2 = phi_3 = f(x,y,t):

    d_i phi_a = d_i f   for all a

The cross-product entering rho_top is:

    J_{ab} = (d_x phi_a)(d_y phi_b) - (d_y phi_a)(d_x phi_b)
           = (d_x f)(d_y f) - (d_y f)(d_x f)
           = 0

Since J_{ab} = 0 for all (a,b), rho_top = eps_{abc} phi_c J_{ab} = 0 identically.

**The gradients of all three fields are parallel vectors at every point, so their
cross products vanish.** This is independent of the spatial profile shape.

## What Would Give Nonzero rho_top

Topological charge requires the three field components to be *distinct* functions
with linearly independent gradients. The canonical example is a baby Skyrmion:

    phi = rho_0 (sin f(r) cos(B*theta), sin f(r) sin(B*theta), cos f(r))

where f(0)=pi, f(infinity)=0, and B is the winding number. Here phi maps R^2 -> S^2
with degree B, and rho_top integrates to B. This is a fundamentally different object
from the oscillon -- it requires broken symmetry between the field components.

## Oscillon Properties

The 2D oscillon is confirmed:
- Peak frequency: omega = 0.864 < m = 1.0 (below mass gap --> oscillon)
- Long-lived: survives to t=2000 with slow radiation loss
- Energy: E = 21.5 initially, decays to ~8.2 by t=2000 (absorbing boundary)
- Core fraction: >99% of energy remains localized (f_core > 0.99 until ~t=1800)
- Breathing amplitude: phi_center oscillates between approximately -0.85 and +0.85

## Implications for Maxwell Coupling

Since rho_top = 0 identically for the symmetric oscillon:
- Phase 2 (coupling to A_mu) is moot: there is no source current
- The breathing oscillon cannot generate electromagnetic fields through the
  topological current mechanism
- To get EM from topology in this model, one must use a configuration with
  genuine topological charge (baby Skyrmion), which is a different initial condition

## Files

- `src/maxwell_c.c` -- 2D solver with topological charge measurement
- `data/maxwell_c_ts.tsv` -- time series (phi, energies, Q_top, rho_top)
- `data/maxwell_c_spectrum.tsv` -- frequency spectrum of phi_1(center,t)
- `data/maxwell_c_profile_t*.tsv` -- 2D field snapshots at t=0,10,50,100,200,500,1000,2000
