# EDMS Soliton Particle Framework

## 1. Thesis

Particles are Soler solitons: self-confined solutions of the nonlinear Dirac
equation, coupled to the Maxwell and Einstein fields. The complete classical
theory is the **Einstein-Dirac-Maxwell-Soler (EDMS) system**:

    S = integral d^4x sqrt(-g) [
        R / (16 pi G)                          (gravity)
      + psi-bar (i gamma^mu D_mu - m) psi      (Dirac matter)
      + (lambda/2) (psi-bar psi)^2              (self-interaction / confinement)
      - (1/4) F_{mu nu} F^{mu nu}              (electromagnetism)
    ]

where D_mu = nabla_mu + i e A_mu is the gauge-covariant derivative.

This action contains four sectors:

| Sector | Field | Spin | Role |
|--------|-------|------|------|
| Matter | psi (Dirac spinor) | 1/2 | Particle content |
| Confinement | (psi-bar psi)^2 | 0 | Self-interaction that localizes psi |
| Electromagnetism | A_mu (gauge field) | 1 | Long-range 1/r Coulomb force |
| Gravity | g_{mu nu} (metric) | 2 | Long-range 1/r gravitational force |

Four dimensionful parameters: m (fermion mass), lambda (self-coupling),
e (electric charge), G (Newton's constant). Three dimensionless combinations
control the physics:

    lambda m^2    = self-coupling strength (determines soliton shape)
    alpha_EM      = e^2 / (4 pi) ~ 1/137  (electromagnetic fine structure)
    alpha_G       = G m^2 ~ 10^{-38}       (gravitational fine structure)

The extreme smallness of alpha_G (the hierarchy problem) is not explained by
this framework. It is an input parameter, as in the Standard Model.

### What this framework provides

1. **Spin-1/2 particles** from Dirac structure (not imposed, follows from the equation)
2. **Self-confinement** from the (psi-bar psi)^2 nonlinearity (no external potential)
3. **Finite size** with calculable internal structure (form factors, radii)
4. **Electromagnetic interaction** via standard U(1) gauge coupling (Coulomb 1/r)
5. **Gravitational interaction** via Einstein's equations (Newton 1/r)
6. **Nuclear-range forces** from scalar density exchange (Yukawa e^{-kappa r}/r)

### What it does not provide

- Mass spectrum (why m_p/m_e = 1836): requires multi-flavor extension
- Charge quantization (why e is discrete): requires topological or anomaly argument
- Weak/strong forces: requires non-abelian gauge extension (SU(2) x SU(3))
- Gravity hierarchy (why G is small): fundamental open problem


## 2. Lessons from Previous Approaches

### v2: Skyrme Model (2023-2025)

The Skyrmion (topological soliton of the SU(2) sigma model) was investigated
as a particle model. Code: `/home/d/code/scp/v2/`.

**What worked:**
- Nuclear physics predictions within 3% (M_p = 938 MeV, r_p = 0.84 fm)
- Multi-baryon solutions (B=1-4) via rational map ansatz
- Soliton-soliton scattering: repulsive B+B, inelastic B+Bbar
- Parameter fitting: 1 code E = 9.098 MeV, 1 code L = 0.5624 fm
- 3D verification to <0.1% accuracy

**What failed (all 6 gravity paths):**
- Path 1 (finite-lambda): P/m=2 is algebraic identity, cannot vary
- Path 2 (L6 sextic): nuclear-range only (0.55 fm), 10^37x too strong
- Path 3 (B^0 p coupling): 1/r mechanism works but coupling g_top is FREE
- Path 4 (hopfion metric): anisotropic at core scale, not 1/r
- Path 5 (constraint gravity): three WZW obstructions, Lie algebra not semisimple
- Path 6 (variable signature): Yukawa eigenvalue (0.46 fm), not 1/r

Four quantum avenues (one-loop, WZW rotation, anomaly, spectral) all null.
The topological obstruction pi_3(S^7) = 0 forces epsilon_0^2 = 0 exactly.

**Structural conclusion:** The Skyrme model produces the MECHANISM for 1/r
gravity (B^0 p scalar coupling) but the coupling strength is a free parameter
at all orders. Gravity cannot emerge from the nuclear sector alone.

**Lesson for v5:** Accept gravity as a fundamental field (Einstein). Focus on
the matter sector where the Skyrme model's success (nuclear spectroscopy,
soliton dynamics) can be transferred to the Soler framework.

### v4: Self-Trapped Waves (2025-2026)

Multiple models tested for self-trapping and gravity. Code: `/home/d/code/scp/v4/`.

**Models tested:**

| Model | Self-trapping | Spin | Gravity |
|-------|:---:|:---:|:---:|
| Scalar Q-ball | Via potential | 0 | NULL (massive) |
| DBI Q-ball | BI refractive index | 0 | NULL (decouples) |
| YM breather | None (disperses) | 1 | NULL |
| Boson star | Gravity itself | 0 | External (GR) |
| Soler soliton | (psi-bar psi)^2 | **1/2** | External (GR) |
| EBId soliton | BI + dilaton | 0 | **1/r dilaton** |

**Key finding:** The Soler model is the ONLY model that produces:
- Natural spin-1/2 (from Dirac structure)
- Self-confinement (from nonlinear self-interaction)
- Working EM coupling (charged Soler solitons)
- Working gravitational coupling (Einstein-Dirac-Soler)
- Full EDMS solution (all three combined)

The EBId model produces 1/r gravity via a massless dilaton, but has spin-0
solitons and the dilaton coupling is free (hierarchy unsolved).

**Confinement-range dilemma (resolved):** The mass that confines the soliton
also screens inter-soliton forces to Yukawa range. Resolution: gravity and EM
are mediated by SEPARATE massless fields (graviton, photon) that do NOT
participate in confinement. The Soler self-interaction confines; the gauge and
metric fields mediate long-range forces.


## 3. The EDMS Action

### Full Lagrangian

    L = L_grav + L_Dirac + L_Soler + L_Maxwell

    L_grav   = R / (16 pi G)
    L_Dirac  = psi-bar (i gamma^mu nabla_mu - m) psi
    L_Soler  = (lambda/2) (psi-bar psi)^2
    L_Maxwell = -(1/4) F_{mu nu} F^{mu nu}

Gauge coupling enters through nabla_mu -> D_mu = nabla_mu + i e A_mu in L_Dirac.

### Field equations

**Dirac** (matter):

    (i gamma^mu D_mu - m + lambda (psi-bar psi)) psi = 0

The effective mass m_eff = m - lambda (psi-bar psi) is REDUCED inside the
soliton, creating a self-consistent potential well.

**Maxwell** (EM):

    nabla_mu F^{mu nu} = e psi-bar gamma^nu psi

The Dirac current j^nu = e psi-bar gamma^nu psi sources the EM field.

**Einstein** (gravity):

    G_{mu nu} = 8 pi G (T^Dirac_{mu nu} + T^Maxwell_{mu nu})

Both matter and EM stress-energy source the gravitational field.

### Symmetries

- Poincare invariance (Lorentz + translations)
- U(1) gauge invariance (local phase rotation of psi)
- Diffeomorphism invariance (general covariance)
- Discrete: C (charge conjugation), P (parity), T (time reversal)
- Global U(1)_V: psi -> e^{i alpha} psi (conserved fermion number Q)

### Conserved quantities

- Energy-momentum: P^mu (from Poincare)
- Angular momentum: J^{mu nu} (from Lorentz, intrinsic spin-1/2 from Dirac)
- Electric charge: Q_em = e integral j^0 d^3x (from U(1) gauge)
- Fermion number: Q = integral psi-dagger psi d^3x (from global U(1)_V)


## 4. Soler Matter Sector

### The flat-space soliton

In flat space (G=0, e=0), the Soler equation reduces to:

    (i gamma^mu d_mu - m + lambda (psi-bar psi)) psi = 0

Spherically symmetric ansatz for the ground state (kappa = -1, j = 1/2):

    psi = e^{-i omega t} [ f(r) chi_{1/2}, i g(r) sigma.r-hat chi_{1/2} ]^T

Radial ODEs:

    f' = (m + omega - lambda S) g
    g' + 2g/r = (m - omega - lambda S) f

where S = f^2 - g^2 (scalar density, psi-bar psi in flat space).

Boundary conditions:
- f(0) = f_0 > 0 (shooting parameter), f'(0) = 0
- g(0) = 0, g'(0) = (m - omega - lambda f_0^2) f_0 / 3
- f, g -> 0 exponentially as r -> infinity

### Established numerical results (v4, m=1, lambda=1)

| omega | f(0) | Q | M = omega Q | R_rms | E_self |
|-------|------|---|-------------|-------|--------|
| 0.50 | 1.381 | 384 | 192 | 2.61 | 46.0 |
| 0.70 | 1.361 | 111 | 78 | 2.26 | 23.5 |
| 0.80 | 1.275 | 69 | 55 | 2.33 | 23.1 |
| 0.90 | 1.065 | 48 | 43 | 2.78 | 11.7 |
| 0.95 | 0.838 | 47 | 44 | 3.60 | 7.0 |
| 0.99 | 0.419 | 74 | 73 | 7.48 | — |

Properties:
- Frequency range: omega in [omega_min, m), omega_min ~ 0.13 m
- Minimum fermion number: Q_min ~ 47 (at omega ~ 0.95 m)
- Exponential tail: f, g ~ e^{-kappa r}/r, kappa = sqrt(m^2 - omega^2)
- Lower component: g < 0 throughout, |g| << f
- Self-interaction energy: E_self = 2 pi lambda integral r^2 S^2 dr = 8-18% of total

### Physical structure

The soliton has definite spin J = 1/2 from the Dirac spinor structure. The
upper component f(r) gives the "charge density" profile; the lower component
g(r) gives the "magnetization density." Their ratio |g/f| ~ kappa/(m+omega)
reflects the relativistic nature of the bound state.

The scalar density S = f^2 - g^2 > 0 everywhere (upper component dominates).
This means m_eff = m - lambda S < m inside the soliton: the effective mass is
REDUCED, creating a self-consistent attractive potential well.


## 5. Electromagnetic Sector

### Coupling mechanism

The U(1) gauge field A_mu couples minimally to the Dirac current:

    D_mu psi = (nabla_mu + i e A_mu) psi

For static spherically symmetric solutions, only A_0 = -phi(r) is nonzero
(electrostatic potential). The Dirac frequency shifts: omega -> omega_eff(r) =
omega - e phi(r) = omega + e V(r), where V = -phi >= 0 for positive charge.

### Established results (v4, flat space, omega=0.9, lambda=1)

| e | f(0) | Q | M | R_rms | E_C/M |
|---|------|---|---|-------|-------|
| 0 | 1.065 | 48.3 | 43.5 | 2.78 | 0% |
| 0.05 | 0.955 | 42.2 | 38.0 | 2.98 | 2.1% |
| 0.10 | 0.644 | 33.8 | 30.4 | 3.71 | 5.3% |
| 0.15 | 0.396 | 23.4 | 21.1 | 4.50 | 6.8% |
| 0.20 | 0.242 | 16.3 | 14.7 | 5.19 | 7.1% |
| 0.22 | 0.121 | 4.8 | 4.3 | 5.42 | ~7% |

Key properties:
- Critical charge: e_crit ~ 0.23 (omega=0.9). No solution above this.
- Coulomb self-energy caps at ~7% of mass before unbinding.
- Soliton "puffs up" with charge: R grows from 2.78 to 5.42.
- Self-limiting: larger e reduces f_0 and Q, reducing the source for V.
- Non-perturbative: even e=0.05 reduces Q by 13%.

### Inter-soliton EM force

Two well-separated charged Soler solitons interact via Coulomb:

    V_EM(R) = e_1 e_2 / (4 pi R)

This is the standard 1/r electromagnetic interaction, mediated by the massless
photon. It is repulsive for same-sign charges, attractive for opposite signs.

The soliton's charge is Q_em = e * Q_fermion, where Q_fermion is the
conserved fermion number. At e=0.22, Q_fermion drops to 4.8, giving very
small total charge. The connection between fermion number and electric charge
is through the coupling constant e, which is a free parameter.


## 6. Gravitational Sector

### Coupling mechanism

The Soler soliton couples to gravity through its stress-energy tensor, which
sources Einstein's equations:

    G_{mu nu} = 8 pi G T_{mu nu}

For static spherically symmetric solutions, use the FSY metric:

    ds^2 = -T(r)^{-2} dt^2 + A(r)^{-1} dr^2 + r^2 d Omega^2

The coupled Einstein-Dirac-Soler system is 4 first-order ODEs for
(alpha, beta, A, T), where alpha = r f, beta = -r g.

### Established results (v4, omega=0.9, lambda=0.5, N=2)

| G | alpha_1 | Q | M_total | M_ADM | E_grav/M | GM/R |
|---|---------|---|---------|-------|----------|------|
| 0 | 1.065 | 48.3 | 98.6 | 98.6 | 0% | 0 |
| 10^{-6} | 1.065 | 48.3 | 98.6 | 98.6 | ~0 | ~0 |
| 10^{-3} | 1.141 | 48.6 | 100.7 | 97.2 | 3.5% | 0.039 |

Mass decomposition: M_ADM = N omega Q + E_self - E_grav

Key properties:
- Gravity contracts the soliton: R_rms drops 7-11% at G=0.001.
- Gravitational binding: E_grav/M = 2.4-4.8% at GM/R = 0.024-0.064.
- Maximum compactness: solutions fail above GM/R ~ 0.1-0.2.
- For real nucleons: G_code ~ 10^{-38} -> corrections COMPLETELY negligible.

### Inter-soliton gravitational force

Two solitons interact via Newton:

    V_grav(R) = -G M_1 M_2 / R

This is the standard 1/r gravitational attraction. It is 10^{38} times weaker
than the Coulomb interaction for unit-charged solitons. The weakness is an
INPUT (the value of G), not derived.

### Full EDMS (gravity + charge combined)

The complete 5-ODE system (alpha, beta, A, T, phi) with simultaneous gravity
and EM has been solved numerically via SCF iteration (v4, src/edm_soler.c):

- SCF method: fix phi -> Dirac bisection -> Green's function Maxwell -> under-relax
- Simultaneous integration FAILS (Coulomb accumulates, omega_eff runaway)
- Combined G=0.001, e=0.1: M_ADM=55.53, E_EM=3.13, GM/R=0.015
- Gravity contracts, charge inflates: competing effects verified


## 7. Particle Properties and Observables

### Mass

The soliton mass depends on omega, m, lambda:

    M_sol = omega * Q + E_self

where Q = 4 pi integral (f^2 + g^2) r^2 dr is the fermion number and
E_self = 2 pi lambda integral r^2 (f^2-g^2)^2 dr is the self-interaction
energy.

The rest mass of a free fermion quantum is m. The soliton mass M_sol < m * Q
(the soliton is bound). The binding fraction is:

    B = 1 - M_sol / (m * Q)

For omega=0.90: B = 1 - 43.5/48.3 = 10%. Substantial binding.

### Charge radius and mass radius

From v4 (standard Q-ball form factors, adapted for Soler):

    R_ch^2 = integral rho_ch(r) r^2 d^3x / Q
    R_M^2  = integral T_{00}(r) r^2 d^3x / M

The Soler soliton has R_ch ~ 2.3-3.6 (depending on omega) in units of 1/m.
The mass radius R_M > R_ch by ~13% because gradient energy peaks at the
soliton edge, not center.

### Spin and magnetic moment

The Dirac structure gives J = 1/2 automatically. The magnetic moment:

    mu = (e / 2m) * g_factor

For a point Dirac particle, g = 2 exactly. For the Soler soliton, the
extended structure modifies this. Computing the anomalous magnetic moment
(g - 2) for the Soler soliton is an OPEN problem and a key v5 target.

### Excited states

The ground state has no radial nodes (f > 0 everywhere). Excited states with
radial nodes exist but have higher energy and are expected to be unstable.
Angular excitations (higher kappa) produce solitons with J = 3/2, 5/2, ...
These could model Delta, N* resonances.


## 8. Inter-Soliton Forces

### Force channels

Between two Soler solitons at separation R >> R_sol:

| Channel | Mediator | Range | Sign | Strength |
|---------|----------|-------|------|----------|
| Coulomb | Photon (massless, spin-1) | 1/R | repulsive (same Q) | alpha_EM ~ 1/137 |
| Newton | Graviton (massless, spin-2) | 1/R | attractive | alpha_G ~ 10^{-38} |
| Scalar | psi-bar psi fluctuation (massive, spin-0) | e^{-kappa R}/R | attractive | lambda^2 |
| Pseudoscalar | psi-bar gamma_5 psi (massive, spin-0) | e^{-kappa R}/R | depends on spin | lambda^2 |

The scalar exchange (psi-bar psi channel) is the Soler analog of pion exchange
in nuclear physics. The range is 1/kappa = 1/sqrt(m^2 - omega^2). At omega=0.9:
kappa = 0.436 m, range = 2.29/m ~ 2.3 code lengths.

### Nuclear binding

Two solitons can bind if the scalar attractive force exceeds the Coulomb
repulsion at some separation. The binding condition is roughly:

    lambda^2 * e^{-kappa R} / R > alpha_EM / R

This gives a critical separation R_bind ~ (1/kappa) * ln(lambda^2 / alpha_EM).
For lambda = 1, alpha_EM = 1/137: R_bind ~ 5/kappa ~ 11 code lengths.
Nuclear binding is expected at separations of a few soliton radii.

### Atom formation

An "atom" = soliton + anti-soliton (opposite charge) bound by Coulomb:

    V_atom(R) = -alpha_EM / R + V_scalar(R)

The Bohr radius a_B = 1/(alpha_EM * mu_reduced) where mu_reduced = M_sol/2.
For alpha_EM = 1/137, M_sol ~ 43 code: a_B ~ 6300 code lengths >> R_sol.
The atom is much larger than the soliton, justifying the point-particle
approximation at atomic scales.


## 9. Multi-Flavor Extension

### Motivation

A single Dirac field produces one type of particle. Real physics has multiple
particle species (electron, proton, neutron, ...) with different masses.

### Two-flavor EDMS

Introduce two Dirac fields psi_1, psi_2 with different masses m_1, m_2:

    L_matter = sum_i [ psi_i-bar (i gamma^mu D_i_mu - m_i) psi_i
                       + (lambda_i/2) (psi_i-bar psi_i)^2 ]
             + lambda_12 (psi_1-bar psi_1)(psi_2-bar psi_2)

The cross-coupling lambda_12 mediates inter-species nuclear forces.
Each species has its own soliton solutions with characteristic M_i, R_i, Q_i.

### Particle identification (schematic)

| Soliton | Physical particle | m_i | Properties |
|---------|-------------------|-----|------------|
| psi_1 (heavy, large lambda) | Proton/neutron | m_N | R ~ 0.84 fm, M ~ 938 MeV |
| psi_2 (light, small lambda) | Electron | m_e | R ~ 0 (point-like limit), M ~ 0.511 MeV |

The electron limit: for m_e very small or lambda_e very small, the Soler
soliton becomes very extended (R -> infinity as omega -> m). In the limit
lambda -> 0, psi_2 is a free Dirac field (point particle). The electron
may be better modeled as a fundamental Dirac quantum rather than a soliton.

### Non-abelian gauge extension

Replace U(1) with SU(N):

    D_mu = nabla_mu + i g_s T^a A^a_mu

SU(2): weak force (W, Z bosons). Requires symmetry breaking for masses.
SU(3): strong force (gluons). Asymptotic freedom, confinement.

The Soler self-interaction lambda (psi-bar psi)^2 could be an EFFECTIVE
description of the non-abelian gauge dynamics at low energy. In QCD, the
quark condensate <psi-bar psi> != 0 generates a similar self-interaction
through the NJL (Nambu-Jona-Lasinio) mechanism. This connection is
speculative but worth exploring.


## 10. Parameter Space and Physical Matching

### Dimensional analysis

The theory has 4 parameters: m, lambda, e, G. In natural units (hbar=c=1):

    [m] = mass = 1/length
    [lambda] = 1/mass^2 = length^2
    [e] = dimensionless
    [G] = 1/mass^2 = length^2

The dimensionless coupling lambda m^2 controls soliton shape. All soliton
properties scale with m (mass sets the overall scale):

    M_sol = m * F_M(lambda m^2, omega/m)
    R_sol = (1/m) * F_R(lambda m^2, omega/m)

where F_M, F_R are dimensionless functions.

### Matching to nucleon

Physical nucleon: M_p = 938.3 MeV, r_p = 0.841 fm.
Dimensionless product: M_p * r_p / (hbar c) = 938.3 * 0.841 / 197.3 = 4.00.

The Soler soliton at m=1, lambda=1, omega=0.9 has M_sol * R_sol = 43.5 * 2.78 / (1/1) = 121.
This does NOT match 4.00, so m=1, lambda=1 is NOT the physical point.

The matching condition M_sol * R_sol = 4.00 constrains lambda m^2 (one equation,
one unknown after fixing omega/m). A systematic scan of lambda m^2 is needed
to find the physical parameter point.

### Parameter count

| Parameter | Determines | Physical value |
|-----------|-----------|----------------|
| m | Overall mass scale | ~ 10 MeV (if nucleon is Q~3 soliton) |
| lambda m^2 | Soliton shape (M*R) | To be determined |
| omega/m | Specific soliton in family | To be determined |
| e | Electromagnetic coupling | e^2/(4pi) = 1/137 |
| G | Gravitational coupling | G m^2 ~ 10^{-38} |

Given M_p and r_p, we can determine m and lambda m^2. Then e and G are
known from experiment. This leaves omega/m as the one free parameter
that selects which soliton in the family represents the nucleon.


## 11. Connection to Standard Model

### What matches

| Feature | EDMS Soler | Standard Model |
|---------|-----------|----------------|
| Spin-1/2 | From Dirac equation | From Dirac equation |
| EM coupling | U(1) gauge, Coulomb 1/r | U(1) gauge, Coulomb 1/r |
| Gravity | Einstein, Newton 1/r | Einstein, Newton 1/r |
| Nuclear force | Scalar exchange, Yukawa | Pion exchange, Yukawa |
| Confinement | Soler self-interaction | QCD (SU(3) gauge) |
| Finite size | Soliton radius R_sol | Proton radius r_p |

### What differs

| Feature | EDMS Soler | Standard Model |
|---------|-----------|----------------|
| Confinement mechanism | Scalar (psi-bar psi)^2 | Non-abelian gauge SU(3) |
| UV behavior | Finite (soliton has no UV divergence) | Renormalizable (but with UV cutoff) |
| Particle number | Q >= 47 (m=1,lambda=1) | Q = 3 (quarks in proton) |
| Mass spectrum | Continuous (omega free) | Discrete (QCD eigenvalues) |
| Chiral symmetry | Broken by m only | Spontaneously broken (pions as Goldstone) |
| Asymptotic freedom | No (lambda is fixed) | Yes (g_s runs with scale) |

### The NJL connection

The Nambu-Jona-Lasinio (NJL) model IS a Soler-type theory:

    L_NJL = psi-bar (i slash-d - m_0) psi + G_NJL [(psi-bar psi)^2 + (psi-bar i gamma_5 tau psi)^2]

It is the low-energy effective theory of QCD after integrating out gluons.
The NJL coupling G_NJL plays the same role as lambda in the Soler model.
NJL solitons have been studied (Kahana & Ripka, 1984; Alkofer et al., 1995)
and produce nucleon-like solutions with Q = 3 (three quarks).

The key difference: NJL includes the CHIRAL channel (psi-bar i gamma_5 psi)^2,
which gives rise to pions as (pseudo-)Goldstone bosons. This is essential for
nuclear physics. The pure Soler model (scalar channel only) does NOT have
chiral symmetry breaking and therefore does NOT produce pions.

**Implication for v5:** Extending Soler to include the pseudoscalar channel
(chiral NJL) is essential for realistic nuclear physics. This is a natural
next step.


## 12. Open Questions

### Immediate (addressable numerically)

1. **Parameter matching**: Find lambda m^2 such that M_sol * R_sol = 4.00.
   What is Q at the physical point? Is Q ~ 3 achievable?

2. **Anomalous magnetic moment**: Compute g-2 for the Soler soliton.
   The extended structure should give g != 2.

3. **Chiral extension**: Add (psi-bar i gamma_5 psi)^2 to the Soler model.
   Does this change the soliton properties? Does it produce pion-like modes?

4. **Two-soliton scattering**: Numerical time evolution of two Soler solitons
   in 1+1D or 3+1D. Elastic? Inelastic? Bound state formation?

5. **Stability**: Are Soler solitons stable under ALL perturbations (not just
   radial)? The VK criterion (dQ/d omega < 0) gives a necessary condition.

### Structural (requires analytical work)

6. **Charge quantization**: Is there a topological argument that forces
   Q_em = integer * e? In the current model, charge is continuous.

7. **Mass hierarchy**: Can two-flavor EDMS naturally produce m_p/m_e = 1836?
   Or is this just two free parameters?

8. **Pion mass**: In the chiral extension, the pion mass m_pi comes from
   explicit chiral symmetry breaking (m != 0). Can we reproduce
   m_pi = 140 MeV from the soliton parameters?

9. **Asymptotic freedom**: The Soler coupling lambda does not run. Can the
   non-abelian gauge extension produce running coupling?

10. **Graviton from soliton**: Is there any limit of the EDMS system where
    the graviton (spin-2 massless mode) emerges from the matter sector?
    All v2 evidence says NO, but the question remains open.


## 13. Research Program

### Phase 1: Parameter Space Exploration

Systematically scan the (lambda m^2, omega/m) plane for the Soler soliton.
For each point, compute M, R, Q, E_self. Find the curve M*R = 4.00
(physical nucleon constraint). Determine whether Q ~ 3 is achievable.

**Code**: Port soler.c from v4 with lambda as a scan parameter.

### Phase 2: Chiral Extension (NJL Soliton)

Add the pseudoscalar channel (psi-bar i gamma_5 psi)^2 to the Soler model.
Derive the modified radial ODEs. Solve numerically. Key questions:
- Does the soliton still exist?
- What is the pion mass (lowest-lying fluctuation mode)?
- Does Q decrease (can we get Q ~ 3)?

**Code**: New solver njl.c extending soler.c with the chiral channel.

### Phase 3: Electromagnetic Properties

Compute form factors F_ch(q), F_M(q) for the Soler soliton at the physical
parameter point. Extract:
- Charge radius r_ch
- Magnetic radius r_M
- Anomalous magnetic moment g-2

**Code**: Extend with form factor computation (as in v4 qball.c).

### Phase 4: Multi-Soliton Dynamics

Two-soliton scattering in the Soler model. Options:
- 1+1D Gross-Neveu model (dimensionally reduced Soler): faster, proof of concept
- 3+1D full Soler: expensive but realistic

Key observables: elastic cross-section, inelastic channels, bound state formation.

**Code**: New PDE solver for time evolution.

### Phase 5: Nuclear Physics

With the chiral extension and multi-flavor setup:
- Two-nucleon potential from soliton-soliton interaction
- Deuteron binding energy
- Nuclear matter (many-soliton) at finite density
- Compare with Skyrme model results from v2

### Phase 6: Gravitational Predictions

At the physical parameter point:
- Post-Newtonian corrections to soliton structure
- Gravitational form factor
- Comparison with neutron star equation of state
- Test whether EDMS gives different predictions from point-particle GR


## 14. Code Architecture

### Directory structure

    /home/d/code/scp/v5/
        FOUNDATION.md       (this document)
        RESULTS.md          (numerical results registry)
        src/
            soler.c         (port from v4, with lambda scan)
            njl.c           (chiral extension, Phase 2)
            charged_soler.c (port from v4)
            ed_soler.c      (port from v4)
            edm_soler.c     (port from v4)
            scatter.c       (soliton scattering, Phase 4)
            Makefile
        data/               (output profiles and scan data)
        tasks/              (task tracking)

### Dependencies

All code is standalone C with only -lm. No external libraries.
Build: gcc -O2 -Wall -Wextra -std=c11 -lm.

### Porting from v4

The following solvers transfer directly:
- soler.c: core Soler solver (Phase 1 starting point)
- charged_soler.c: EM coupling (Phase 3)
- ed_soler.c: gravitational coupling (Phase 6)
- edm_soler.c: full EDMS (Phase 6)

New code needed:
- njl.c: chiral NJL extension (Phase 2, major new work)
- scatter.c: time-dependent PDE solver (Phase 4, major new work)
