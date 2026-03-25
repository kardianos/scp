# Stability Characterization Research: Methods & Metrics

## 0. Context and Data Summary

We have three surviving structures from Gen 4 long-run (T=200, N=192, L=25):
- **S20** (braid 2.0x): S_final=2.34, E_pot_final=-431, P_int retention=40%
- **UDD_R4** (UDD wide tubes): S_final=2.30, E_pot_final=-392, P_int retention=75%
- **CB15** (counter-braid 1.5x): S_final=0.82, E_pot_final=-134, P_int retention=17%

All remain alive at T=200 but with strongly oscillating metrics. The key observation
from the JSON data: these structures undergo quasi-periodic oscillations in E_pot
and P_int with a period of roughly 20-30 time units, with occasional deep dips (e.g.,
CB15 at t=135 has S=0.64, while at t=10 it was S=5.27). This oscillation is
fundamental to the braid's traveling wave nature and is NOT instability per se.

Available per-voxel data in SFA frames (12 columns):
- phi_0, phi_1, phi_2 (SFA_POSITION, components 0-2)
- theta_0, theta_1, theta_2 (SFA_ANGLE, components 0-2)
- phi_vel_0, phi_vel_1, phi_vel_2 (SFA_VELOCITY, components 0-2)
- theta_vel_0, theta_vel_1, theta_vel_2 (SFA_VELOCITY, components 3-5)

All computable forces can be reconstructed from the field values using the EOM:
```
  phi_acc_a = lap(phi_a) - m^2 phi_a - dV/dP * dP/dphi_a + eta * curl(theta)_a
  theta_acc_a = lap(theta_a) - m_theta^2 theta_a + eta * curl(phi)_a
```
with V(P) = (mu/2) P^2 / (1 + kappa P^2), P = phi_0 * phi_1 * phi_2.

---

## 1. Literature-Informed Approaches

### 1.1 Derrick/Virial Analysis (Soliton Stability)

**Background**: Derrick's theorem (1964) states that in d >= 2 spatial dimensions,
static finite-energy solutions of scalar field theories with only quadratic kinetic
terms do not exist. The standard workaround is topological charge (Skyrme model) or
time-dependent structures (oscillons, Q-balls). Our braids are explicitly
time-dependent traveling waves, so Derrick's static theorem does not apply directly.
However, the virial method generalizes: define the scaling functional

    E[lambda] = lambda^(3-2) E_grad + lambda^3 E_mass + lambda^(3-4) E_pot_4

where the scaling exponents come from phi(x) -> phi(lambda*x). At equilibrium
dE/dlambda|_1 = 0, giving the virial constraint:

    E_grad - 3*E_mass - E_pot = 0  (approximate, for our system)

**Mapping to our system**: The triple-product potential V(P) scales as lambda^(-3)
under spatial dilation (three fields contribute P = phi_0*phi_1*phi_2, each
contributing lambda^(-3/2) from the compression). The kinetic energy (time
derivatives) provides the dynamic stabilization absent in Derrick's static case.

**Concrete metric**: Per-cluster virial ratio

    R_virial = (E_grad + E_coupling) / |E_pot|

where E_grad, E_pot are integrated over the cluster volume. For a balanced structure,
R_virial should be near a fixed value (estimated ~1.0-1.5). Deviation indicates
expansion (R_virial >> 1, gradient energy dominates) or collapse (R_virial << 1,
potential energy dominates).

**Stability indicator**: A structure that maintains R_virial within +/- 30% of its
time-averaged value is vibrationally stable. Secular drift in R_virial signals
impending death.

### 1.2 Q-criterion (Turbulence/Coherent Structure Detection)

**Background**: In CFD, the Q-criterion (Hunt et al., 1988) identifies coherent
vortical structures by examining the velocity gradient tensor:

    Q = (1/2)(|Omega|^2 - |S|^2)

where Omega_ij = (dv_i/dx_j - dv_j/dx_i)/2 is the rotation rate tensor and
S_ij = (dv_i/dx_j + dv_j/dx_i)/2 is the strain rate tensor. Regions with Q > 0
are vortex-dominated (rotation > strain), indicating coherent structures.

**Mapping to our system**: Define the velocity field from the phi sector:

    v_a(x) = phi_vel_a(x)  (the field time derivative IS the velocity)

The velocity gradient tensor Gamma_{ai} = d(phi_vel_a)/dx_i is a 3x3 tensor
(field component a, spatial direction i). This is NOT the same as a fluid velocity
gradient (which maps spatial directions to spatial directions), but the decomposition
into symmetric (strain-like) and antisymmetric (rotation-like) parts is still
meaningful:

    S_{ai} = (Gamma_{ai} + Gamma_{ia})/2   (symmetric: compressive deformation)
    W_{ai} = (Gamma_{ai} - Gamma_{ia})/2   (antisymmetric: internal rotation)

The analog Q-criterion:

    Q_phi = (1/2)(|W|^2 - |S|^2)

Q_phi > 0 indicates regions where the field velocities have a rotational character
(coherent vortex-like structure) rather than compressive/extensional (dispersive).

**Caveat**: Since a=0,1,2 are field components not spatial components, this tensor
mixes internal and spatial indices. A cleaner version uses the energy-momentum tensor
approach (Section 1.4). However, Q_phi is still useful as an empirical discriminant.

**Stability indicator**: High Q_phi concentrated at braid cores = coherent.
Low or negative Q_phi at cores = strain-dominated = dispersing.

### 1.3 Lambda-2 Criterion (Vortex Core Detection)

**Background**: The lambda-2 criterion (Jeong & Hussain, 1995) is more robust than
Q for identifying vortex cores. It uses the eigenvalues of S^2 + W^2 (where S and
W are the symmetric and antisymmetric parts of the velocity gradient tensor). A
vortex core exists where the second eigenvalue lambda_2 < 0.

**Mapping**: Define the 3x3 tensor per voxel:

    M_{ab} = sum_i (S_{ai} S_{bi} + W_{ai} W_{bi})

where the sum is over spatial directions i=x,y,z. Compute eigenvalues lambda_1 >=
lambda_2 >= lambda_3 of M. Regions where lambda_2 < 0 have a pressure minimum due
to vortical motion — in our context, they are regions where the field's internal
rotation creates a binding "well."

**Stability indicator**: Surviving structures should have persistent lambda_2 < 0
regions in their cores. Dispersing structures will see lambda_2 become positive as
the rotational coherence is lost.

**Implementation note**: This requires computing 18 first derivatives (3 field
components x 3 spatial directions x 2 for symmetric/antisymmetric), then a 3x3
eigenvalue decomposition per voxel. The eigenvalue solve is O(1) per voxel (closed
form for 3x3 symmetric matrices). Total: ~150 lines of C.

### 1.4 Energy-Momentum Tensor Analysis (Field Theory Standard)

**Background**: The canonical approach in field theory is the energy-momentum tensor.
For our system with Lagrangian density L = (1/2)|dphi|^2 - (1/2)m^2|phi|^2 - V(P)
+ theta terms:

    T^{mu nu} = sum_a (d_mu phi_a)(d_nu phi_a) + (theta terms) - eta^{mu nu} L

The spatial components T^{ij} give the stress tensor, whose eigenvalues determine
whether the field is under tension (positive eigenvalue, wants to contract) or
pressure (negative eigenvalue, wants to expand).

**Concrete decomposition per voxel**:

    T^{ij} = sum_a (d_i phi_a)(d_j phi_a) + sum_a (d_i theta_a)(d_j theta_a) - delta^{ij} L

The trace gives the isotropic pressure:

    P_iso = (1/3) Tr(T^{ij}) = (1/3)(sum_a |grad phi_a|^2 + |grad theta_a|^2) - L

The traceless part sigma^{ij} = T^{ij} - P_iso delta^{ij} gives the deviatoric
(shear) stress.

**Stability indicator**: In a stable soliton, the stress tensor should be dominated
by negative pressure (tension) in the core (binding the structure) and positive
pressure at the surface (the outward radiation pressure balanced by the binding).
The ratio |sigma|/|P_iso| measures anisotropy: near-spherical stable structures
should have low shear.

### 1.5 Correlation Function Methods (Condensed Matter)

**Background**: In condensed matter, order parameter correlation functions
G(r) = <phi(x)phi(x+r)> characterize the spatial coherence of a phase. For our
three-field system, the natural order parameter is the triple product P(x).

**Concrete metrics**:

(a) Spatial autocorrelation of P:

    C_P(r) = <P(x) P(x+r)> / <P^2>

where the average is over x within a cluster. The correlation length xi defined by
C_P(xi) = 1/e gives the coherent size of the triple-product structure.

(b) Phase-phase correlation: Each phi_a can be decomposed into amplitude and phase
locally. Near the braid core, phi_a ~ A_a(x,y) cos(k_z z + delta_a + phase_a(x,y)).
The phase coherence:

    C_phase(r) = <cos(phase(x) - phase(x+r))>

measures whether the traveling wave maintains its phase structure across the cluster.

**Stability indicator**: Long correlation length (xi > R_core) means the structure
is a coherent soliton. Short xi (xi < R_core) means the triple-product arrangement
is disordered — the structure is fragmenting.

### 1.6 Topological Analysis (Helicity/Linking)

**Background**: The Cosserat system has an analog of magnetic helicity:

    H_phi = integral phi . curl(phi) dV
    H_theta = integral theta . curl(theta) dV

In the Skyrme model, topological charge (baryon number) is the integral of
the topological current B^0 = epsilon^{ijk} epsilon_{abc} (d_i phi_a)(d_j phi_b)(d_k phi_c).
Our braids are not topologically stabilized in the Skyrme sense (they lack the
sigma-model constraint |q|=const), but the helicity integral measures the degree
of twist/linking in the field.

**Concrete metric**: The phi helicity density per voxel:

    h(x) = phi_0(curl_phi)_0 + phi_1(curl_phi)_1 + phi_2(curl_phi)_2

where (curl_phi)_a is computed from the OTHER two phi components in the usual
curl formula. This is already computed in the coupling energy term of analyze_sfa.c
(the eta coupling term E_coupling = -eta * integral phi . curl(theta)).

**Cross helicity** (phi-theta linking):

    H_cross = integral phi . curl(theta) dV = -E_coupling / eta

This measures the degree of coupling between the phi and theta sectors.

**Stability indicator**: Structures with high |H_phi| have twisted field
configurations that resist unwinding. The time derivative dH_phi/dt should be
near zero for stable structures and rapidly decreasing for dispersing ones.

---

## 2. Per-Voxel Stability Metrics

For each voxel at position (i,j,k), the following quantities can be computed from a
single SFA frame containing the 12 columns (phi_a, theta_a, phi_vel_a, theta_vel_a).

### 2.1 Field Density and Triple Product

**rho(x) = phi_0^2 + phi_1^2 + phi_2^2**

- What it measures: Local field energy content (proportional to mass term energy)
- Compute: Sum of squares of three phi fields
- Stability: High rho (> 0.5) at cluster cores, decaying radially.
  Stable if core rho is maintained across frames.
- Implementation: 3 multiplies + 2 adds per voxel. ~10 lines.

**P(x) = phi_0 * phi_1 * phi_2**

- What it measures: Triple product entering the binding potential V(P)
- Compute: Product of three phi values
- Stability: |P| > 0.01 defines the "bound" region. The spatial volume where
  |P| > threshold is the cluster's binding volume. Shrinking binding volume = death.
- Implementation: 1 multiply per voxel. ~5 lines.

**V(x) = (mu/2) P^2 / (1 + kappa P^2)**

- What it measures: Local potential energy density (mu < 0 means V < 0 = binding)
- Compute: From P, 2 multiplies + 1 divide
- Stability: V < 0 at cores means binding. Deeper V = stronger binding.
  Typical stable core: V ~ -0.1 to -1.0 per voxel.
- Implementation: ~5 lines.

**sigma_theta(x) = theta_0^2 + theta_1^2 + theta_2^2**

- What it measures: Local theta field amplitude (Cosserat torsion)
- Stability: Theta grows from zero via the curl coupling. Steady theta_rms
  (~0.01-0.03) indicates equilibrium. Growing theta_rms may indicate energy
  draining from phi to theta (dissolution channel noted in V39).
- Implementation: Same as rho. ~10 lines.

### 2.2 Gradient-Based Quantities

**|grad(rho)| = sqrt((d rho/dx)^2 + (d rho/dy)^2 + (d rho/dz)^2)**

- What it measures: How sharply the field density varies. Large gradients at
  cluster surfaces; small gradients inside cores and in the background.
- Compute: Central differences on rho field. 6 neighbor lookups per voxel.
- Stability: Sharp, well-defined boundaries (high |grad rho| at surface, low
  inside) indicate a coherent structure. Diffuse boundaries (low gradient
  everywhere) indicate dispersal.
- Implementation: ~20 lines.

**div(v) = d(phi_vel_0)/dx + d(phi_vel_1)/dy + d(phi_vel_2)/dz**

- What it measures: The divergence of the field velocity field, interpreted as a
  "compressibility" indicator. Positive div(v) = local expansion, negative =
  contraction. Note: this treats the three field components as a spatial vector
  field, which is physically meaningful because the braid has helical structure
  that maps field components to spatial directions.
- Compute: Central differences on the three velocity fields. 6 lookups.
- Stability: Time-averaged div(v) ~ 0 in stable cores (oscillating around zero).
  Persistent positive div(v) = expanding = dispersing.
- Implementation: ~20 lines.

**curl(phi)_a and curl(theta)_a**

- What it measures: The rotational content of the field. Curl(phi) is the source
  for theta acceleration; curl(theta) sources phi acceleration (Cosserat coupling).
- Compute: Standard curl formula with central differences. Already implemented
  in the simulation code (curl_component function).
- Stability: Coherent curl patterns (organized, smooth) in stable cores. Random,
  rapidly varying curl in unstable regions. The magnitude |curl(phi)| at the core
  indicates the strength of the driving force for theta growth.
- Implementation: ~30 lines (3 components x 4 difference terms each).

### 2.3 Energy Density Decomposition

Per-voxel energy densities (these are ALL computable from a single frame):

**e_kin(x) = (1/2) sum_a phi_vel_a^2**

- Kinetic energy density. High in oscillating cores and in radiation zones.

**e_grad(x) = (1/2) sum_a |grad phi_a|^2**

- Gradient energy density. High at boundaries, low in cores and background.

**e_mass(x) = (1/2) m^2 sum_a phi_a^2 = (1/2) m^2 rho**

- Mass energy density. Proportional to field amplitude squared.

**e_pot(x) = V(P) = (mu/2) P^2 / (1 + kappa P^2)**

- Potential energy density. Negative where binding occurs.

**e_theta_kin(x) = (1/2) sum_a theta_vel_a^2**

**e_theta_grad(x) = (1/2) sum_a |grad theta_a|^2**

**e_coupling(x) = -eta * sum_a phi_a * curl(theta)_a**

- Phi-theta coupling energy density.

**Total: e_total(x) = e_kin + e_grad + e_mass + e_pot + e_theta_kin + e_theta_grad + e_coupling**

**Stability metric from energy**: The binding fraction at each voxel:

    f_bind(x) = e_pot(x) / (e_grad(x) + e_mass(x))

This is negative in bound regions (where V(P) < 0 compensates gradient + mass cost).
The "binding radius" R_bind = max distance from centroid where f_bind < -0.01
defines the cluster's extent.

**Implementation**: ~80 lines total for all 7 energy terms.

### 2.4 Acceleration and Force Decomposition

**Reconstructed acceleration**: Given phi_a and theta_a at a single frame, the
acceleration is exactly computable:

    a_phi_a(x) = lap(phi_a) - m^2 phi_a - (dV/dP)(dP/dphi_a) + eta curl(theta)_a

Breaking this into four force terms per field component:

    F_lap_a(x)     = lap(phi_a)                   [Laplacian / wave propagation]
    F_mass_a(x)    = -m^2 phi_a                   [Mass / restoring force]
    F_pot_a(x)     = -(dV/dP)(dP/dphi_a)          [Triple-product potential force]
    F_curl_a(x)    = eta curl(theta)_a             [Cosserat coupling force]

where dV/dP = mu P / (1 + kappa P^2)^2 and dP/dphi_a is the product of the other
two fields.

**What this reveals**:
- In stable cores: F_pot should be comparable to F_mass + F_lap (balanced forces)
- In dispersing regions: F_lap dominates (pure wave propagation, no binding)
- The coupling force F_curl is typically 5-10% of F_pot (from the eta=0.5 and
  theta_rms ~ 0.01-0.03 data)

**Net acceleration magnitude**:

    |a(x)|^2 = sum_a (F_lap_a + F_mass_a + F_pot_a + F_curl_a)^2

Low |a| in cores = quasi-equilibrium (forces balanced). High |a| at surfaces =
active dynamics (radiation, reshaping).

**Force balance ratio**:

    R_force(x) = |F_pot| / (|F_lap| + |F_mass|)

R_force ~ 1 means the triple-product coupling is actively binding at this voxel.
R_force << 1 means the voxel is in the free-wave regime (no binding).

**Implementation**: ~100 lines (Laplacian, potential derivative, curl, force balance).

### 2.5 Time-Derivative Metrics (Multi-Frame)

These require comparing TWO consecutive frames (dt = snap_dt = 5.0):

**drho/dt = (rho(t+dt) - rho(t)) / dt**

- Positive = field accumulating. Negative = field dispersing.
- Time-averaged drho/dt ~ 0 in stable cores. Persistently negative = dying.

**dP/dt = (P(t+dt) - P(t)) / dt**

- Rate of change of triple product. Oscillates in braids (traveling wave).
- Variance of dP/dt: high variance = strong oscillation = dynamically active.
  Low variance = either dead (P~0) or frozen (not a braid).

**Implementation**: ~30 lines (requires loading two frames).

---

## 3. Cluster-Level Metrics

These are computed by first identifying clusters (connected components where
|P| > threshold, already implemented in spatial_analysis.c), then integrating
over each cluster's voxels.

### 3.1 Cluster Virial Ratio

**Definition**:

    R_virial = (E_grad_cluster + E_theta_grad_cluster) / |E_pot_cluster|

where each energy is integrated over voxels belonging to the cluster (|P| > threshold).

**Computation**:
1. Run BFS to identify cluster voxels (existing code)
2. Sum e_grad(x) and e_pot(x) over cluster voxels
3. Divide

**Physical meaning**: A stable bound structure should have R_virial near a
characteristic value. For the Skyrme model, E_2 = E_4 at equilibrium (Derrick).
For our time-dependent braids, the virial relation includes the kinetic energy:

    E_kin + E_grad - E_mass - E_pot = 0  (time-averaged)

so the relevant ratio is (E_kin + E_grad) / (E_mass + |E_pot|).

**Stability prediction**: R_virial consistently in [0.8, 1.5] = stable.
R_virial drifting outside this range = impending restructuring.

**Implementation**: ~40 lines (on top of existing cluster detection).

### 3.2 Angular Momentum Content

**Definition**: The angular momentum of the field within a cluster:

    L_i = epsilon_{ijk} integral (x_j - c_j) * p_k dV

where p_k = sum_a phi_a * phi_vel_a * (x_k component) ... but more precisely,
the momentum density of the phi field is:

    pi_a(x) = phi_vel_a(x)  (conjugate momentum to phi_a)

For the "spatial" angular momentum (treating the field energy flow as a fluid), the
field momentum density vector is:

    g_i(x) = -sum_a (d phi_a / dx_i) * phi_vel_a

(This is T^{0i} from the energy-momentum tensor.) Then:

    L_i = epsilon_{ijk} integral (x_j - c_j) * g_k(x) dV

**Physical meaning**: Braids are traveling waves with helical structure. Their angular
momentum content reflects the degree of rotational coherence. A perfectly helical
braid along z has L_z proportional to its winding number and amplitude.

**Stability prediction**: Conserved |L| (or slowly varying) = stable coherent
rotation. Rapidly changing L = structure breaking apart or reforming.

**Implementation**: ~50 lines (momentum density computation + cross product
integration).

### 3.3 Phase Coherence Index

**Definition**: For a braid along axis n, each phi_a should have a traveling-wave
form phi_a ~ A_a(r_perp) cos(k n.x + delta_a - omega t). The triple product
P = phi_0 phi_1 phi_2 then has a specific spatial pattern. The phase coherence
index measures how well the actual field matches this traveling-wave template.

**Practical computation** (template-free version): Compute the normalized triple
product:

    P_norm(x) = P(x) / (rho(x))^{3/2}

This is bounded by |P_norm| <= 1/(3 sqrt(3)) (AM-GM inequality, equality when
all three |phi_a| are equal). The ratio:

    C_phase = <|P_norm|> / (1/(3 sqrt(3)))

measures how efficiently the three fields combine into a triple product. C_phase = 1
means the three fields have equal amplitudes everywhere (perfect phase relationship).
C_phase << 1 means one field dominates at most locations.

**Alternative (spectral)**: Fourier transform P(x) along the braid axis. A coherent
braid will have a sharp peak at k = k_braid. The spectral concentration:

    C_spec = |P_hat(k_braid)|^2 / sum_k |P_hat(k)|^2

measures what fraction of the triple-product power is in the fundamental mode.

**Stability prediction**: C_phase > 0.3 and/or C_spec > 0.5 = coherent braid.
Below these thresholds = phase disorder = unstable.

**Implementation**: Template-free version ~30 lines. Spectral version ~60 lines
(needs 1D FFT along the braid axis, or approximate with sum over z-slices).

### 3.4 Surface-to-Volume Ratio

**Definition**: Using the |P| > threshold cluster mask:

    V_cluster = N_voxels * dV  (volume: count of above-threshold voxels)
    S_cluster = N_surface * dA (surface: voxels with at least one below-threshold neighbor)

    SV_ratio = S_cluster / V_cluster^{2/3}

Normalize by the sphere value 4.836 to get the "sphericity" index:

    Sigma = SV_ratio / 4.836

**Physical meaning**: Sigma = 1 for a sphere. Sigma > 1 for elongated or irregular
shapes. Braids are inherently Sigma >> 1 (elongated along the propagation axis).
Sigma increasing over time means the structure is fragmenting (developing more surface
area per unit volume).

**Stability prediction**: Sigma stable or decreasing = compact, stable.
Sigma increasing = fractal fragmentation = unstable.

**Implementation**: ~25 lines (count surface voxels by checking 6-neighbors).

### 3.5 Triple-Product Concentration (Q-criterion analog)

**Definition**: The fraction of the cluster's total |P| that is concentrated in
the inner core:

    Q_P = integral_{r < R_half} |P| dV / integral_{cluster} |P| dV

where R_half is the half-mass radius (the radius enclosing 50% of the cluster's
total |P|).

**Physical meaning**: Q_P measures whether the triple-product binding is concentrated
in a tight core (high Q_P = compact, efficient binding) or spread diffusely
(low Q_P = the binding is weak and distributed).

**Computation**:
1. Compute centroid of |P| within cluster (already done in spatial_analysis)
2. Sort cluster voxels by distance from centroid
3. Find R_half where cumulative |P| reaches 50% of total
4. Integrate |P| within R_half (which is trivially 50% by definition, so instead
   use a FIXED radius R_core = 3.0 code units, which is the expected braid core width)

Better formulation with fixed core radius:

    Q_P = integral_{r < R_core} |P| dV / integral_{cluster} |P| dV

with R_core = 3.0 (from the known braid core width of ~3 code units at R_tube=3).

**Stability prediction**: Q_P > 0.7 = tight core binding. Q_P < 0.3 = diffuse,
weakly bound. Time series of Q_P should be stable for survivors.

**Implementation**: ~35 lines.

### 3.6 Binding Energy Decomposition

**Definition**: For each cluster, decompose the total energy into:

    E_self = integral_cluster (e_kin + e_grad + e_mass + e_theta terms) dV
    E_bind = integral_cluster e_pot dV  (should be negative for bound structures)
    E_total = E_self + E_bind

The binding fraction:

    f_b = E_bind / E_self

**Physical meaning**: f_b < 0 means the cluster is bound (potential energy reduces
total energy). For nuclear physics analogy: f_b ~ -0.01 to -0.1 is typical for
nuclei (binding energy is 1-10% of rest mass).

From the JSON data:
- S20 at t=200: E_pot = -431, E_total = 3449. Rough f_b ~ -431/3880 ~ -0.11
- UDD_R4: E_pot = -392, E_total = 3459. f_b ~ -392/3851 ~ -0.10
- CB15: E_pot = -134, E_total = 5273. f_b ~ -134/5407 ~ -0.025

CB15 has the weakest binding fraction, consistent with its lower S_final.

**Stability prediction**: |f_b| > 0.05 = adequately bound. |f_b| < 0.02 = marginal.
|f_b| > 0.15 = strongly bound.

**Implementation**: ~30 lines (sums over cluster voxels, already partially done).

### 3.7 Theta Equilibrium Ratio

**Definition**:

    R_theta = E_theta / E_phi

where E_theta = E_theta_kin + E_theta_grad and E_phi = E_kin + E_grad. Integrated
over the cluster.

**Physical meaning**: The phi-theta coupling transfers energy from phi to theta
(and back) via the curl terms. At equipartition (thermal equilibrium), we expect
R_theta to approach a specific value determined by eta and the mode structure. From
the data, theta_rms stabilizes around 0.01-0.03, while phi amplitudes are 0.6-1.0,
suggesting R_theta ~ (0.02/0.8)^2 ~ 0.0006, i.e., theta carries << 1% of the energy.

**Stability prediction**: R_theta should converge to a steady value. If R_theta is
still growing at T=200 (theta still absorbing energy from phi), the structure has
not reached equilibrium and may still be evolving. If R_theta has plateaued, the
structure has equilibrated its internal mode coupling.

**Implementation**: ~20 lines.

### 3.8 Oscillation Frequency and Quality Factor

**Definition**: From the time series of E_pot(t) for each cluster, extract the
dominant oscillation frequency omega_osc via peak-to-peak timing, and the quality
factor:

    Q_osc = omega_osc * E_mean / (dE/dt)_dissipation

where the dissipation rate is estimated from the secular decay of E_total.

**Physical meaning**: The braid oscillates because it is a traveling wave — the
triple product P varies as the phase structure rotates. The oscillation frequency
is related to the braid's internal clock (omega ~ k*c ~ pi/L * c for the fundamental
mode). The quality factor measures how many oscillation cycles the structure survives.

From the data, the oscillation period is roughly 20-30 time units and E_total
decays by roughly 50% over T=200, giving Q ~ (2*pi/25) * 200 / ln(2) ~ 36 cycles.

**Stability prediction**: Q_osc > 20 = long-lived. Q_osc < 5 = will die soon.

**Implementation**: ~40 lines (peak detection in time series + exponential decay fit).

---

## 4. Composite Stability Index

### 4.1 Proposed Per-Cluster Score: S_stability

Combine the above metrics into a single stability score per cluster per frame:

    S_stab = w1 * norm(R_virial) + w2 * norm(f_b) + w3 * norm(C_phase)
           + w4 * norm(Q_P) + w5 * norm(R_theta_eq) + w6 * norm(1/Sigma)

where norm() maps each quantity to [0, 1] via:
- norm(R_virial) = exp(-(R_virial - R_target)^2 / (2 * 0.3^2)),  R_target ~ 1.2
- norm(f_b) = min(1, |f_b| / 0.10)
- norm(C_phase) = min(1, C_phase / 0.4)
- norm(Q_P) = Q_P  (already in [0,1])
- norm(R_theta_eq) = 1 - |R_theta - R_theta_target| / R_theta_target
- norm(1/Sigma) = min(1, 1/Sigma)

**Suggested weights**: w1=0.25, w2=0.25, w3=0.20, w4=0.15, w5=0.10, w6=0.05

The weights emphasize energy balance (virial, binding fraction) and phase coherence
as the primary stability indicators, with geometric measures as secondary.

### 4.2 Time-Integrated Stability

The per-frame S_stab is noisy (braids oscillate). The time-integrated version:

    S_integrated(T) = (1/T) integral_0^T S_stab(t) dt

and the minimum over any window of width W:

    S_min_W = min over [t, t+W] of { (1/W) integral_t^{t+W} S_stab(s) ds }

with W = 30 (roughly one oscillation period). S_min_W captures the worst-case
stability — a structure that has S_min_W > 0.3 never goes through a phase where
it is at risk of dissolution.

---

## 5. ML/LLM Prediction Approaches

### 5.1 Feature Engineering for t=0 Prediction

**Goal**: From a single snapshot at t=0, predict whether the structure survives
to t=200 (binary classification) or predict S_final (regression).

**Feature vector per initial condition** (computed from the t=0 frame):

| # | Feature | Formula | Rationale |
|---|---------|---------|-----------|
| 1 | rho_max | max(phi_0^2 + phi_1^2 + phi_2^2) | Peak field density |
| 2 | P_int_0 | integral |P| dV | Total triple-product content |
| 3 | E_pot_0 | integral V(P) dV | Initial binding energy |
| 4 | f_b_0 | E_pot_0 / E_self_0 | Initial binding fraction |
| 5 | R_virial_0 | E_grad_0 / |E_pot_0| | Initial virial ratio |
| 6 | R_core_0 | R_rms of |P| distribution | Spatial compactness |
| 7 | V_bind_0 | Volume where V(P) < -0.01 | Binding volume |
| 8 | C_phase_0 | <|P_norm|> / max(P_norm) | Phase coherence |
| 9 | Q_P_0 | Core-fraction of |P| | Core concentration |
| 10 | sigma_0 | SV_ratio / 4.836 | Shape regularity |
| 11 | E_kin_0 | integral e_kin dV | Initial kinetic energy |
| 12 | A_max_0 | max amplitude of any phi_a | Peak single-field amplitude |
| 13 | aspect_0 | I_max / I_min | Anisotropy |
| 14 | L_z_0 | Angular momentum along braid axis | Rotational content |
| 15 | K_balance | E_kin_0 / E_grad_0 | Kinetic-gradient balance |

**Target variable**: S_final at T=200, or binary alive/dead at T=200.

### 5.2 Training Data Requirements

**Current data**: We have ~36 distinct initial conditions across Gen 1-4 (12 per
generation, 3 generations), each run for T=30 (Gen 1-3) or T=200 (Gen 4 top 3).
This gives us:
- 36 samples with t=0 features and T=30 outcomes
- 3 samples with t=0 features and T=200 outcomes

This is FAR too little for any ML approach. We need at minimum:
- **Simple logistic regression / random forest**: 50-100 samples
- **Neural network**: 500+ samples
- **3D CNN on raw voxel data**: 5000+ samples (infeasible at N=192)

### 5.3 Feasible ML Strategy

**Option A: Feature-based predictor with expanded training set**

Generate 100-200 random initial conditions (varying A, R_tube, ellip, delta,
chirality, N_braids, orientations) and run each for T=100 at N=128 (cheaper than
N=192). Extract the 15 features at t=0 and label each as alive/dead at t=100.

Training a random forest or gradient-boosted tree on 15 features x 100 samples
is standard and should identify the 3-4 most predictive features.

**Estimated compute cost**: 200 candidates x 1 min/candidate = 200 min ~ 3.3 hours
on V100 at $0.13/hr = $0.43. Very feasible.

**Expected result**: The most predictive features will likely be f_b_0 (initial
binding fraction) and R_virial_0 (virial ratio), based on the physics. If so, we
can bypass ML entirely and just use these two numbers as the stability predictor.

**Option B: Autoencoder anomaly detection on voxel snapshots**

Train a 3D convolutional autoencoder to reconstruct snapshots of STABLE structures.
At test time, high reconstruction error indicates the structure is unlike any known
stable configuration. This is a one-class classification approach.

Architecture: 3D conv encoder (N=192 -> 96 -> 48 -> 24 -> 12 -> latent 64-dim)
followed by a symmetric decoder. Train on frames from the T=50-200 period of
surviving structures.

**Problem**: We only have 3 surviving structures x ~30 useful frames = ~90 training
samples of 192^3 x 6 voxels each. This is marginal for a 3D CNN. Aggressive data
augmentation (rotations, reflections, small noise) could push it to ~1000 effective
samples, which might be enough for a simple autoencoder.

**Implementation effort**: Requires Python + PyTorch, not C. ~200 lines of Python.
NOT recommended as a first approach — try Option A first.

**Option C: LLM-based pattern recognition (not recommended)**

Feed textual summaries of the 15 features to an LLM and ask it to predict stability.
This is essentially using the LLM as a very expensive lookup table. Not recommended
because: (a) the feature space is numeric and low-dimensional, perfect for
traditional ML; (b) LLMs have no physics intuition for this specific system;
(c) the cost per prediction is orders of magnitude higher than a random forest.

### 5.4 Recommendation

**Start with Option A** (feature-based random forest). The required steps are:
1. Implement the 15 features as a C tool (~200 lines)
2. Generate 200 random initial conditions with the seed generator
3. Run each for T=100 on GPU (~3 hours, $0.43)
4. Extract features at t=0, label at t=100
5. Train random forest in Python (~50 lines with scikit-learn)
6. Identify top predictive features
7. Use those features as the analytical stability predictor (no ML needed at runtime)

---

## 6. Prioritized Implementation Plan

Ordered by expected information value per implementation effort:

### Priority 1: Per-Voxel Energy Decomposition + Force Balance
- **Files needed**: phi_a, theta_a, phi_vel_a, theta_vel_a
- **Output**: 7 energy density fields + 4 force fields per phi component
- **Key metric**: f_bind(x) = e_pot(x) / (e_grad(x) + e_mass(x))
- **Key metric**: R_force(x) = |F_pot| / (|F_lap| + |F_mass|)
- **Implementation**: ~180 lines of C
- **Expected insight**: Map of where binding is strong vs weak within each cluster.
  This is the most direct answer to "what makes a region stable?"

### Priority 2: Cluster Binding Fraction + Virial Ratio
- **Input**: Energy decomposition from Priority 1, cluster masks from spatial_analysis
- **Output**: Per-cluster f_b, R_virial, E_bind time series
- **Implementation**: ~60 lines (on top of Priority 1)
- **Expected insight**: Quantitative answer to "how bound is each cluster?" and
  "is it in virial equilibrium?"

### Priority 3: Phase Coherence Index
- **Input**: phi_a fields
- **Output**: Per-voxel P_norm, per-cluster C_phase
- **Implementation**: ~40 lines
- **Expected insight**: Whether the three fields maintain their phase relationship
  (essential for V(P) binding to work)

### Priority 4: Triple-Product Concentration Q_P
- **Input**: phi_a fields, cluster centroids
- **Output**: Per-cluster Q_P
- **Implementation**: ~35 lines
- **Expected insight**: Whether binding is concentrated in a tight core (stable)
  or diffusely spread (unstable)

### Priority 5: Multi-Frame Time Derivatives
- **Input**: Two consecutive frames
- **Output**: drho/dt, dP/dt per voxel; oscillation statistics per cluster
- **Implementation**: ~50 lines
- **Expected insight**: Whether structures are expanding, contracting, or
  stably oscillating

### Priority 6: Theta Equilibrium Ratio
- **Input**: theta_a fields, phi_a fields
- **Output**: Per-cluster R_theta
- **Implementation**: ~20 lines
- **Expected insight**: Whether the phi-theta energy transfer has equilibrated

### Priority 7: Velocity Gradient Tensor (Q-criterion, lambda-2)
- **Input**: phi_vel_a fields
- **Output**: Per-voxel Q_phi, lambda_2
- **Implementation**: ~150 lines (including 3x3 eigenvalue solver)
- **Expected insight**: Vortex identification in the field — but the physical
  interpretation is less clear than for fluid dynamics because the velocity
  gradient tensor mixes field and spatial indices

---

## 7. Specific Stability/Instability Signatures (Predictions)

Based on the physics and the available data, here are concrete predictions for
what the analysis should reveal:

### Stable Regions Will Have:
- f_bind < -0.05 (significant binding fraction)
- R_force > 0.5 (triple-product force is substantial)
- |P| > 0.01 (active triple-product coupling)
- V(P) < -0.01 per voxel (energy lowered by binding)
- |curl(theta)| < 0.1 (moderate, not runaway theta)
- Balanced oscillation: drho/dt has mean ~ 0, std ~ 0.01-0.1

### Unstable Regions Will Have:
- f_bind > -0.01 (weak or no binding)
- R_force < 0.1 (dominated by wave propagation, not binding)
- |P| < 0.001 (triple product too small to bind)
- V(P) ~ 0 (in the V(P) ~ mu/2 * P^2 linear regime, not the saturated regime)
- High |div(v)| (persistent expansion or contraction)
- drho/dt persistently negative (energy draining away)

### The Critical Transition:
The boundary between stable and unstable is expected at:
- |P| ~ 0.005 (where V(P) transitions from quadratic to saturated regime,
  kappa P^2 ~ 1 gives P ~ 1/sqrt(kappa) = 0.14, but the binding force dV/dP
  peaks at P = 1/sqrt(3*kappa) = 0.082)
- rho ~ 0.3 (roughly twice the background level rho_bg ~ 0.03 * 3 fields ~ 0.1)
- R_force ~ 0.3 (transition from free-wave to bound regime)

The V(P) force is strongest (most binding) when dV/dP is maximized, which occurs at:

    P_peak = 1 / sqrt(3 * kappa) = 1 / sqrt(150) = 0.0816

At this P value: V(P_peak) = mu/(2*3*kappa) * (1/(1+1/3)^2) = mu/6kappa * 9/16

For mu=-41.345, kappa=50: V(P_peak) = -41.345/(6*50) * 0.5625 = -0.0775

So the binding energy density at peak efficiency is about -0.08 per voxel volume.
Regions with |V(P)| > 0.05 per dV are actively and significantly bound.

---

## 8. Summary of Recommended Tool Output

The analysis tool should produce, for each frame of each SFA file:

**Per-voxel maps** (saved as SFA or analyzed in-memory):
1. rho(x), P(x), V(x) — field structure
2. e_kin(x), e_grad(x), e_mass(x), e_pot(x) — energy decomposition
3. f_bind(x), R_force(x) — stability indicators
4. |a(x)| — net acceleration magnitude

**Per-cluster statistics** (printed as TSV or JSON):
1. Cluster ID, N_voxels, centroid, R_rms
2. E_self, E_bind, f_b, R_virial
3. C_phase, Q_P, Sigma (surface/volume)
4. theta_rms_cluster, R_theta
5. |L| (angular momentum magnitude)
6. phi_max_cluster, P_max_cluster

**Global time series** (one row per frame):
1. N_clusters, N_alive_clusters
2. Total E_bind, total E_self
3. Mean and std of per-cluster f_b, R_virial, C_phase
4. Theta equilibrium status (R_theta trend)

This output enables both per-frame snapshot analysis and time-series trend analysis
for predicting long-term stability.
