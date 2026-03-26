# OQ2: Coulomb Force from Cosserat Radiation Pressure

## Analytical Derivation

---

### 0. Setup and Notation

We work with the 6-field Cosserat equation (CONCEPT.md Eq. 10):

    d^2 phi_a/dt^2 = laplacian(phi_a) - m^2 phi_a - dV/dphi_a + eta curl(theta)_a    (phi: massive, m^2=2.25)
    d^2 theta_a/dt^2 = laplacian(theta_a)                      + eta curl(phi)_a      (theta: massless, m_theta=0)

The electromagnetic identification (EM_THEORY.md):

    A = theta           (vector potential)
    E = -d theta/dt     (electric field)
    B = curl(theta)     (magnetic field)
    J_eff = eta curl(phi)   (effective current density)

Natural units: c = 1. The braid oscillates with primary frequency omega (half-period ~2.2 t,
so omega ~ pi/2.2 ~ 1.43 in code units) and breathing envelope frequency Omega ~ 2pi/9.1 ~ 0.69.

---

### 1. The Cosserat Poynting Vector and Energy Conservation

**PROVEN (follows from the equations of motion).**

#### 1.1 theta energy density

Define the theta-sector energy density at a point in the far field (where phi ~ 0,
so the curl coupling source vanishes):

    u = (1/2)(d theta/dt)^2 + (1/2)|grad theta|^2
      = (1/2)|E|^2 + (1/2)|B|^2                                                    (1)

where |grad theta|^2 = sum_a sum_j (d theta_a / dx_j)^2 and |B|^2 = |curl theta|^2.

Note: |grad theta|^2 >= |curl theta|^2 in general (equality holds for transverse waves).
In the far field, theta is purely transverse (see section 2), so |grad theta|^2 = |curl theta|^2
and Eq. (1) is exact.

#### 1.2 Continuity equation

Take the time derivative of u using the theta equation of motion:

    du/dt = (d theta/dt) . (d^2 theta/dt^2) + sum_j (d theta/dx_j) . (d^2 theta / dt dx_j)

Substitute d^2 theta/dt^2 = laplacian(theta) + eta curl(phi):

    du/dt = (d theta/dt) . laplacian(theta) + (d theta/dt) . [eta curl(phi)]
            + sum_j (d theta/dx_j) . (d/dx_j)(d theta/dt)

The second line equals (1/2) d/dt |grad theta|^2, which is already included in du/dt.
Reorganizing properly by expanding the Laplacian term:

    (d theta/dt) . laplacian(theta) = div[(d theta/dt)(grad theta)] - sum_j |d^2 theta / dt dx_j|
                                    ... (standard manipulation)

The clean result uses the identity for vector fields:

    du/dt + div(S) = -J_eff . E                                                     (2)

where the Poynting vector is:

    S = E x B = (-d theta/dt) x (curl theta)                                        (3)

and the source term is:

    -J_eff . E = -[eta curl(phi)] . (-d theta/dt) = eta curl(phi) . (d theta/dt)    (4)

**Proof of Eq. (2)**: This follows the standard electromagnetic proof exactly.
Starting from:

    du/dt = E . (dE/dt) + B . (dB/dt)

where dB/dt = curl(dA/dt) = -curl(E) (Faraday, proven in EM_THEORY.md section 5),
and dE/dt = -d^2 theta/dt^2 = -laplacian(theta) - eta curl(phi) = curl(B) - J_eff
(Ampere, proven in EM_THEORY.md section 5). Then:

    du/dt = E . (curl B - J_eff) + B . (-curl E)
          = E . curl B - B . curl E - J_eff . E
          = -div(E x B) - J_eff . E

using the vector identity div(E x B) = B . curl E - E . curl B. Therefore:

    du/dt + div(S) = -J_eff . E                                                     QED

#### 1.3 Interpretation

In the FAR FIELD (r >> R_braid), J_eff = eta curl(phi) -> 0 because the phi field
is massive (m = 1.5) and decays exponentially as ~exp(-mr)/r. Therefore:

    du/dt + div(S) = 0       (far field, r >> 1/m ~ 0.67 code units)                (5)

Energy is exactly conserved in the far field, and S is the energy flux carried by
the theta radiation. The total radiated power through a sphere of radius r is:

    P_rad = oint S . n_hat dA                                                        (6)

which is independent of r in the far field (energy conservation).

---

### 2. Oscillating Magnetic Dipole Radiation

**PROVEN for the point-dipole approximation. The finite-size corrections are an ASSUMPTION
requiring numerical verification (see caveats).**

#### 2.1 The braid as a current source

The braid's helical structure produces a time-varying effective current
J_eff = eta curl(phi). Since the braid breathes (oscillates in amplitude with
frequency omega), the current oscillates in time. The lowest multipole of the
current distribution is the magnetic dipole moment:

    m(t) = (1/2) integral r x J_eff(r,t) d^3r
         = (eta/2) integral r x [curl phi(r,t)] d^3r                                (7)

The braid's helical axis (say z) defines the dipole orientation. By the oscillation
of the braid:

    m(t) = m_0 cos(omega t) z_hat                                                   (8)

where m_0 depends on the braid's amplitude, geometry, and eta.

**Justification for Eq. (8)**: The dominant time dependence of curl(phi) comes from
the braid's breathing oscillation at frequency omega ~ 1.43. The braid core (r < 3)
contains all the current because phi decays exponentially outside it. The integral
in Eq. (7) is therefore dominated by the core region, where the helical twist is
coherent and oscillates at the breathing frequency. The phase offsets delta =
{0, 3.00, 4.43} break the symmetry so that the three components do not cancel,
leaving a net oscillating dipole along the braid axis.

**ASSUMPTION**: The higher multipoles (quadrupole, octupole) radiate at higher
powers of (omega R_braid / c). For omega ~ 1.43 and R_braid ~ 3, the ratio
omega R / c ~ 4.3, which is NOT small. This means the multipole expansion does
NOT converge rapidly for this system. The point-dipole approximation captures the
correct SCALING (1/r^2 for the flux) but not the correct angular pattern or
numerical prefactor. See caveat 1 below.

#### 2.2 No monopole radiation

Since div(J_eff) = div(eta curl phi) = 0 identically (divergence of a curl),
and therefore div(E) = 0 for all time (proven in EM_THEORY.md section 5),
there is no monopole (l=0) term in the multipole expansion of the radiated field.

**PROVEN**: This is a theorem of the equations, not an approximation. The leading
radiation is dipole (l=1) or higher.

#### 2.3 Far-field radiation pattern

The standard result for magnetic dipole radiation (Jackson, Classical Electrodynamics,
chapter 9) gives the far-field theta wave:

    B_rad(r,t) = -(omega^2 / (4 pi c^2 r)) [n_hat x m(t_ret)] x n_hat              (9)

    E_rad(r,t) = c B_rad x n_hat                                                    (10)

where t_ret = t - r/c is the retarded time, and n_hat = r_hat is the observation
direction. For m = m_0 cos(omega t) z_hat:

    |B_rad|^2 = (omega^4 m_0^2 sin^2 alpha) / (16 pi^2 c^4 r^2) cos^2(omega t_ret) (11)

where alpha is the polar angle from the dipole axis (z).

#### 2.4 Time-averaged Poynting flux

    <S> = <E_rad x B_rad> = (1/2) |E_rad,peak| |B_rad,peak| n_hat

Using <cos^2> = 1/2:

    <|S|> = (omega^4 m_0^2) / (32 pi^2 c^3) . sin^2(alpha) / r^2                   (12)

**This is the key result**: the time-averaged energy flux falls as 1/r^2.

Total radiated power (integrating over solid angle, integral sin^2 alpha d Omega = 8pi/3):

    P_total = (omega^4 m_0^2) / (12 pi c^3)                                         (13)

**PROVEN**: Equations (9)-(13) follow from standard classical electrodynamics applied
to the Cosserat theta sector. They require only that:
  (a) theta is massless (confirmed: m_theta = 0)
  (b) the source (J_eff) is localized (confirmed: phi decays exponentially)
  (c) the far field is at r >> lambda = 2pi c / omega (lambda ~ 4.4 code units)

The 1/r^2 scaling of <|S|> is a consequence of energy conservation on an expanding
wavefront and holds for ANY localized oscillating source, regardless of multipole order.
Higher multipoles also produce <|S|> ~ 1/r^2 with different angular patterns and
stronger omega dependence (omega^{2l+2} for multipole order l), but the same radial scaling.

---

### 3. Effective Coulomb Force from Radiation Scattering

**The 1/r^2 force law follows from the 1/r^2 flux. The identification of the
proportionality constant as a "charge" is an ASSUMPTION.**

#### 3.1 Momentum carried by radiation

Electromagnetic radiation carries momentum density g = S/c^2. When radiation
with flux <S> is scattered by a target, the momentum transfer rate (force) is:

    F = sigma_eff <S> / c                                                            (14)

where sigma_eff is the effective momentum-transfer cross section of the target.

For a target braid at distance r from the source:

    |F| = sigma_eff / c . (omega^4 m_0^2) / (32 pi^2 c^3) . sin^2(alpha) / r^2     (15)

#### 3.2 The 1/r^2 law

Equation (15) gives:

    |F| = K / r^2                                                                    (16)

where:

    K = (sigma_eff omega^4 m_0^2 sin^2 alpha) / (32 pi^2 c^4)                       (17)

**PROVEN**: If the source radiates as a localized oscillator and the target has a
well-defined scattering cross section, the force is 1/r^2. This follows from:
  (i)  <|S|> ~ 1/r^2  (energy conservation on a sphere)
  (ii) F = sigma <S>/c  (momentum conservation in scattering)

Both (i) and (ii) are exact conservation laws. The 1/r^2 force is therefore
model-independent — it holds for ANY localized radiator and ANY localized scatterer
in a theory with massless radiation propagating at c.

#### 3.3 Angle averaging for isotropic force

Equation (15) has a sin^2(alpha) angular dependence, which would make the force
anisotropic for a single dipole. However:

- A composite baryon (UUD) has three orthogonal current loops.
- The sum of three orthogonal sin^2 patterns:
  sin^2(alpha_x) + sin^2(alpha_y) + sin^2(alpha_z) = 2 (constant)
- The composite radiation pattern is ISOTROPIC.

This resolves OQ1 in the affirmative: three orthogonal magnetic dipoles produce
a spherically symmetric radiation pattern, mimicking an electric monopole.

**PROVEN**: The angular sum is an algebraic identity. A UUD composite radiates
isotropically if the three dipole amplitudes are equal.

**ASSUMPTION**: The three braid dipole moments m_0 are equal. This requires the
three braids to have the same amplitude and geometry, which is approximately true
for a symmetric composite but may receive corrections from inter-braid interactions.

---

### 4. Sign of the Force from Winding Number

**ASSUMPTION. The mechanism is physically motivated but not rigorously derived.**

#### 4.1 Scattering phase

When radiation from source A (winding W_A = +1, dipole moment m_A) hits target B
(winding W_B), the scattered wave has a phase that depends on the relative
orientation of m_A and the target's current loop.

For SAME winding (W_A = W_B): The incident B field drives the target's current
loop resonantly. The induced current is IN PHASE with the incident field.
The scattered wave interferes constructively with the incident wave in the
backward direction (toward the source). By Newton's third law (momentum
conservation), the target is pushed TOWARD the source.

For OPPOSITE winding (W_A = -W_B): The induced current is OUT OF PHASE.
The scattered wave interferes constructively in the forward direction (away from
source). The target is pushed AWAY from the source.

This gives:
- Same winding (same "charge"): ATTRACTIVE (consistent with V34 data)
- Opposite winding (opposite "charge"): REPULSIVE

**Important**: This is the sign pattern of Ampere's force law for parallel vs
anti-parallel current loops, NOT Coulomb's law for like vs unlike charges.
However, for composite baryons where the current loops sum to isotropic
radiators (section 3.3), the net force mimics Coulomb:
- Net "charge" = sum of winding numbers (UUD: +1+1-1 = +1, UDD: +1-1-1 = -1)
- Same net charge -> same dominant winding -> net repulsion
- Opposite net charge -> opposite winding -> net attraction

**ASSUMPTION**: The scattering phase argument above is qualitative. A rigorous
derivation would require computing the differential scattering cross section
d(sigma)/d(Omega) for a current loop in an incident plane wave, including the
interference between incident and scattered fields. This is a standard
scattering problem but has not been carried out for the Cosserat braid geometry.

#### 4.2 Consistency with V34 data

V34 measured (D=15, T=150):

    Same winding:     delta_D = -6.41  (27% more attraction than gravity alone)
    Opposite winding: delta_D = -2.16  (57% less attraction)

The EM force adds to gravity for same-winding and opposes gravity for
opposite-winding, consistent with section 4.1. The asymmetry (27% vs 57%)
arises because gravity (always attractive) and EM (sign-dependent) are
comparable in magnitude at D=15. The pure EM component is:

    F_EM(same) = -6.41 - (-5.05) = -1.36 (attractive)
    F_EM(opp)  = -2.16 - (-5.05) = +2.89 (repulsive)

The ratio |F_EM(opp)/F_EM(same)| = 2.1, not 1.0. This asymmetry is expected
for magnetic dipole interactions (the force depends on relative orientation
of the dipole moments, and the winding reversal does not simply flip the sign
of m but changes the angular radiation pattern).

---

### 5. The Effective Coupling Constant

**ASSUMPTION. The identification is dimensional analysis; the numerical value
requires simulation input.**

#### 5.1 Effective charge

Define the effective "charge" squared by equating the radiation-pressure force
to the Coulomb force:

    F = e_eff^2 / (4 pi r^2)  =  sigma_eff |m_0|^2 omega^4 / (32 pi^2 c^4 r^2)

Therefore:

    e_eff^2 = sigma_eff |m_0|^2 omega^4 / (8 pi c^4)                                (18)

Each factor depends on the braid parameters:

    m_0 ~ eta integral |r x curl(phi)| d^3r  ~  eta R_braid^2 A_braid               (19)
    omega ~ pi/2.2 ~ 1.43 (code units)                                               (20)
    sigma_eff ~ pi R_braid^2 (geometric cross section, order of magnitude)            (21)

So:

    e_eff^2 ~ eta^2 R_braid^6 A_braid^2 omega^4 / c^4                               (22)

#### 5.2 Fine structure constant

    alpha = e_eff^2 / (4 pi hbar c)                                                  (23)

In this classical theory, hbar enters through the identification of the braid's
action quantum with Planck's constant — this is a SEPARATE identification not
derived here. Given that identification:

    alpha ~ eta^2 R_braid^6 A_braid^2 omega^4 / (4 pi hbar c^5)                     (24)

The key dependence is alpha ~ eta^2. This is testable: varying eta in simulation
should change the effective coupling as eta^2, measurable via the force ratio
F_EM/F_grav at fixed separation.

**ASSUMPTION**: Equation (24) provides the functional form but not the numerical
value. Computing alpha from first principles requires:
  (a) The scattering cross section sigma_eff for a braid (not computed)
  (b) The magnetic dipole moment m_0 from the braid's internal structure (not computed)
  (c) The precise frequency spectrum (broadband, not monochromatic — see caveat 2)

---

### 6. Caveats and Required Verifications

#### Caveat 1: Finite Source Size

The braid core has R_braid ~ 3 code units. The radiation wavelength is
lambda = 2pi/omega ~ 4.4 code units. The ratio R_braid/lambda ~ 0.7 is
NOT small. Therefore:

- The point-dipole approximation is QUALITATIVE, not quantitative.
- Higher multipoles (quadrupole, octupole) contribute significantly.
- The angular radiation pattern is NOT pure sin^2(alpha).
- The numerical prefactor in the force (Eq. 17) is unreliable.

**However**: The 1/r^2 radial scaling is EXACT for ANY localized source,
regardless of multipole content. This is because energy conservation requires
the total flux through any enclosing sphere to be constant, so the average
flux must fall as 1/r^2. The force law F ~ 1/r^2 is robust; only the angular
pattern and the overall coefficient K are affected by finite-size corrections.

**Numerical test**: Compute the angular power spectrum P_l(r) of the theta
field at r = 10, 15, 20 from a single braid. If l=1 (dipole) dominates, the
point-dipole approximation is good. If l=2 or higher are comparable, the
multipole expansion converges slowly but the 1/r^2 scaling still holds.

#### Caveat 2: Broadband Spectrum

The braid does NOT oscillate at a single frequency. V34 measured:
- Fast oscillation: half-period 2.2 t (omega_1 ~ 1.43)
- Breathing envelope: period 9.1 t (Omega ~ 0.69)
- Possible higher harmonics from nonlinear V(P) coupling

For a broadband source, the Poynting flux is the sum over all frequency components:

    <|S|> = sum_n (omega_n^4 |m_n|^2 sin^2 alpha) / (32 pi^2 c^3 r^2)              (25)

Each term falls as 1/r^2, so the total also falls as 1/r^2. The broadband
nature affects the COUPLING CONSTANT (effective charge) but NOT the force law.

Specifically, the radiation-pressure force becomes:

    F = sigma_eff(omega) / c . <|S(omega)|>                                          (26)

If the scattering cross section sigma_eff depends on frequency (as it generally
does for resonant scattering), the effective charge becomes a weighted integral
over the source spectrum:

    e_eff^2 = integral sigma_eff(omega) P(omega) d omega / (4 pi c)                 (27)

where P(omega) is the spectral power density. This is still a CONSTANT (independent
of r), so F ~ 1/r^2 is preserved.

**Numerical test**: Compute the power spectrum of theta at several radii.
If the spectrum is independent of r (as it must be for freely propagating
radiation), the force law is confirmed.

#### Caveat 3: Near-Field Transition

The derivation assumes r >> lambda ~ 4.4. At r ~ lambda (near field), the
fields have both radiative (1/r) and reactive (1/r^2, 1/r^3) components.
The reactive components carry no net energy flux but contribute to the
instantaneous Poynting vector. The time-averaged force in the near field
includes additional terms:

    F_near ~ 1/r^2 (radiation pressure) + 1/r^4 (magnetic dipole-dipole)
             + 1/r^6 (quadrupole) + ...                                              (28)

At nuclear distances (r ~ 2-10 code units ~ 1-6 fm), all these terms contribute
comparably, and the force is the full wave-mediated interaction observed in V34.
The "Coulomb" regime (pure 1/r^2) emerges only at r >> 10 code units (~6 fm).

This is consistent with the spatial force decomposition in EM_THEORY.md section 3:
F_pot/F_curl transitions from 5:1 at r=2 to 0:1 at r=22, with mixed character
at intermediate distances.

#### Caveat 4: Back-Reaction and Self-Consistency

This derivation treats the radiation field as a TEST field — computed from one
braid, evaluated at another. In a self-consistent treatment:

- Each braid radiates AND scatters.
- The scattered radiation reaches other braids, creating a many-body problem.
- Energy radiated by A and absorbed by B affects B's oscillation, which affects
  B's radiation, which feeds back to A.

For WEAK coupling (alpha << 1), the back-reaction is a small correction and the
test-field approximation is valid to leading order. Since the physical fine
structure constant alpha ~ 1/137 is indeed small, this is self-consistent.

---

### 7. Summary of Results

| Statement | Status | Dependencies |
|-----------|--------|-------------|
| S = E x B is the theta energy flux | **PROVEN** | Equations of motion |
| du/dt + div S = -J.E | **PROVEN** | Vector calculus identity |
| div E = 0 for all time | **PROVEN** | div(curl) = 0, initial conditions |
| No monopole radiation | **PROVEN** | Consequence of div E = 0 |
| <\|S\|> ~ 1/r^2 in far field | **PROVEN** | Energy conservation + localized source |
| F = sigma <S>/c ~ 1/r^2 | **PROVEN** | Momentum conservation in scattering |
| Three orthogonal dipoles -> isotropic pattern | **PROVEN** | sin^2 sum identity |
| Same winding attracts, opposite repels | **CONSISTENT** with V34 | Scattering phase argument (qualitative) |
| alpha ~ eta^2 | **DIMENSIONAL ANALYSIS** | Requires sigma_eff, m_0 numerically |
| Point-dipole approximation | **POOR** (R/lambda ~ 0.7) | But 1/r^2 scaling is robust |
| Broadband spectrum preserves 1/r^2 | **PROVEN** | Sum of 1/r^2 terms is 1/r^2 |

### The Core Argument in Three Lines

1. A localized oscillating current source radiates theta waves with flux <|S|> ~ 1/r^2.
   (Energy conservation on an expanding wavefront.)

2. A distant target scatters this radiation, absorbing momentum F = sigma <S>/c ~ 1/r^2.
   (Momentum conservation in scattering.)

3. Therefore F ~ 1/r^2, which is the Coulomb force law.

The 1/r^2 is NOT a dynamical prediction of the Cosserat equation specifically.
It is a kinematic consequence of (a) massless radiation in 3+1 dimensions and
(b) localized sources and targets. Any theory with these two properties produces
a 1/r^2 radiation-pressure force. The Cosserat equation's specific contribution
is providing the physical mechanism: braids are the sources (via eta curl phi),
theta is the radiation (massless), and winding number determines the sign.
