# Normal Mode Analysis: B=1 Hedgehog Skyrmion

## Summary

Computed breathing mode (K=0) eigenfrequencies and moment of inertia for the
B=1 hedgehog soliton. **No breathing resonance exists in the sigma model.**

## Method

Linearize the Skyrme EOM around the hedgehog profile f_0(r):

  f(r,t) = f_0(r) + epsilon * g(r) * cos(omega * t)

This gives the Sturm-Liouville eigenvalue problem:

  -(P(r) g')' + W(r) g = (omega^2/c^2) m(r) g

where:
- P(r) = 2r^2 + 4c_4 sin^2(f_0)   (stiffness)
- m(r) = r^2 + 2c_4 sin^2(f_0)     (inertia, = P/2)
- W(r) = second variation of E_2 + E_4 + E_V

Key: near r=0, the regular solution behaves as g ~ r^nu where
nu = (-1 + sqrt(1 + 2W_0))/2 = 3.79 (universal for all e, rho_0).

Code: `src/normal_modes.c` — RK4 shooting + bisection.

## Results: Sigma Model (massless pions, e=1, rho_0=1)

### Box eigenvalues (g(R_max=15) = 0)

| n  | lambda  | omega  | E (MeV) | nodes | lambda_free | ratio |
|----|---------|--------|---------|-------|-------------|-------|
|  0 |  0.163  | 0.404  |   3.7   |   0   |   0.088     | 1.86  |
|  1 |  0.421  | 0.649  |   5.9   |   2   |   0.351     | 1.20  |
|  2 |  0.810  | 0.900  |   8.2   |   2   |   0.790     | 1.03  |
|  3 |  1.408  | 1.187  |  10.8   |   4   |   1.404     | 1.00  |
|  4 |  2.210  | 1.487  |  13.5   |   5   |   2.193     | 1.01  |
|  5 |  3.205  | 1.790  |  16.3   |   6   |   3.158     | 1.01  |
|**6**|**4.386**|**2.094**|**19.1**|   6   |   4.299     | 1.02  |
|  7 |  5.750  | 2.398  |  21.8   |   7   |   5.615     | 1.02  |
|  8 |  7.296  | 2.701  |  24.6   |   8   |   7.106     | 1.03  |
|  9 |  9.021  | 3.004  |  27.3   |   9   |   8.773     | 1.03  |

### Comparison with free-particle box modes

Free particle: lambda_n = 2(n*pi/R_max)^2 for the spherical standing wave j_0(kr).

- **n=0**: 86% above free (strong soliton interaction)
- **n=3 and above**: within 1% of free (soliton perturbation negligible)
- **n=6 (19.1 MeV)**: only 2% above free

### Key result: the 19 MeV spectral peak

The FFT of the B+Bbar scattering data (scatter_anti_v0.50.dat) showed a dominant
peak at omega = 2.09 code = 19.1 MeV. This matches box mode n=6 **exactly**.

Mode n=6 has period T = 2*pi/2.09 = 3.0 code time units, which equals the
soliton light-crossing time 2R_rms = 2*1.495 = 2.99 code units.

However, n=6 is only 2% above the free-particle eigenvalue — it is **not** a
soliton-specific resonance. It is simply the standing wave whose wavelength
matches the soliton diameter. The soliton collision excites this wavelength
preferentially because it matches the spatial scale of the perturbation.

**Conclusion: no breathing mode resonance exists in the massless sigma model.**
The Skyrmion's effect on the K=0 continuum is a weak perturbation (except for
the lowest 2-3 modes).

## Effective potential shape

| Region      | r range  | W/m ratio | Character         |
|-------------|----------|-----------|-------------------|
| Core        | 0-0.3    | +36 (huge)| Strong repulsion  |
| Well        | 0.3-0.8  | -8 to 0  | Attraction        |
| Transition  | 0.8-2.0  | 0 to +3  | Rising barrier    |
| Asymptotic  | >2.0     | 4/r^2->0 | Free propagation  |

The well depth (W/m ~ -8) is too shallow to support bound states: the
continuum threshold is at lambda = 0 (massless pions), and no eigenvalues
with lambda < 0 were found.

## With pion mass (m_pi = 15.3 code = 139.6 MeV)

**Warning**: used sigma model profile (solved WITHOUT pion mass). The pion mass
changes the equilibrium profile, so these results are approximate.

Adding V = m_pi^2 * integral(1 - cos(f)) * d^3x to the energy:
- W gains a term m_pi^2 * r^2 * cos(f_0), which is NEGATIVE near the soliton
  center (cos(pi) = -1), deepening the potential well to W/m ~ -18.
- Continuum threshold rises to lambda = m_pi^2 = 234.
- 7 bound states found below threshold.
- **Ground state has lambda < 0** (omega^2 < 0 = UNSTABLE).
  This instability is an artifact of using the wrong equilibrium profile.
  The sigma model profile is NOT a minimum of the energy with pion mass.

To properly study massive pion modes, the radial solver needs to include the
pion mass term in the equilibrium ODE.

## Moment of Inertia and Delta-N Splitting

Lambda = (8*pi*rho_0^2)/(3c^2) * integral r^2 sin^2(f_0) [1 + c_4(f'^2 + sin^2(f)/r^2)] dr

| Parameter | e=4, rho_0=1 | e=1, rho_0=1 |
|-----------|-------------|-------------|
| c_4       | 0.125       | 2.0         |
| Lambda    | 2.323       | 141.6       |
| Delta_E code | 0.646    | 0.0106      |
| 1 code E  | 36.4 MeV   | 9.1 MeV     |
| Delta_E phys | 23.5 MeV | 0.1 MeV     |
| Experiment | 293.7 MeV  | 293.7 MeV   |

The Delta-N splitting is highly parameter-dependent. The classical Skyrme model
significantly underpredicts it at our chosen parameters.

Standard Skyrme result (ANW 1983): e_std = 5.45, F_pi = 129 MeV gives Delta_E = 294 MeV.
Our constraint (M_p, r_p matching) gives e_std ~ 3.69, which underpredicts Delta_E.

## Implications for Phase 12 (Boson Identification)

1. **No breathing resonances** in the massless sigma model. Meson-like features
   will NOT appear in the K=0 (radial perturbation) channel.

2. **Meson candidates require angular modes**: The rho meson (775 MeV) should
   appear in the K=1 (grand spin) channel, and f_2(1270) in K=2. These are
   substantially more complex to compute (coupled radial equations).

3. **Pion mass essential**: Without a pion mass, all radiation is massless and
   there are no genuine resonances. Adding m_pi creates a mass gap and may
   support discrete modes. Requires re-solving the profile ODE first.

4. **The 19 MeV spectral peak is physical but mundane**: It's the soliton's
   natural ringing frequency (light-crossing time), not a meson signature.
   It corresponds to the soliton "vibrating at its natural size."

## Files

- `src/normal_modes.c` — Breathing mode eigenvalue solver
- Profile used: `data/profiles/profile_sigma_e1.dat` (e=1, sigma model)
- Profile used: `data/profiles/profile_B1.dat` (e=4, sigma model)
