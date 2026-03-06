# V19: Three-Mode Resonance on S³ — Numerical Method

## 1. Objective

Test whether a **non-topological, dynamic three-mode resonance** on the S³ sigma model
with Skyrme stabilization can:

1. Form a long-lived localized oscillation (B=0 oscillon) whose lifetime exceeds
   the naive dispersal time by at least 10×
2. Exhibit a far-field time-averaged field deviation that falls as a power law r^{-n}
   with n ≤ 2 (interesting) or n = 1 (gravity-relevant)
3. Phase-lock three orthogonal quaternion modes through the intrinsic S³ nonlinearity

If Phase 1 succeeds, Phase 2 adds a conserved density field ρ (v6 architecture) and
tests whether the resonance sources a 1/r density depletion.

---

## 2. Mathematical Framework

### 2.1 Field and Target Space

The dynamical field is a unit quaternion at each spacetime point:

$$
q(x,t) = (q_0, q_1, q_2, q_3) \in S^3 \subset \mathbb{R}^4
$$

subject to the constraint:

$$
|q|^2 \equiv q_0^2 + q_1^2 + q_2^2 + q_3^2 = 1
$$

The vacuum is $q_\mathrm{vac} = (1, 0, 0, 0)$.

The three independent degrees of freedom per point are the tangent-space components
$(q_1, q_2, q_3)$, with $q_0 = \sqrt{1 - q_1^2 - q_2^2 - q_3^2}$ determined by the
constraint.

### 2.2 Lagrangian Density

The Skyrme sigma-model Lagrangian on S³:

$$
\mathcal{L} = \frac{1}{2} |\dot{q}|^2 - \frac{1}{2} |\nabla q|^2
  - \frac{c_4}{4} \left[ (\mathrm{Tr}\,D)^2 - \mathrm{Tr}(D^2) \right]
$$

where the dot denotes $\partial_t$ and the spatial strain tensor is:

$$
D_{ij} = \partial_i q \cdot \partial_j q
  = \sum_{a=0}^{3} (\partial_i q_a)(\partial_j q_a)
$$

Expanding the terms explicitly:

$$
\mathrm{Tr}\,D = D_{11} + D_{22} + D_{33}
  = \sum_{i=1}^{3} \sum_{a=0}^{3} (\partial_i q_a)^2
  = |\nabla q|^2
$$

$$
\mathrm{Tr}(D^2) = \sum_{i,j=1}^{3} D_{ij}^2
  = \sum_{i,j=1}^{3} \left( \sum_{a=0}^{3} (\partial_i q_a)(\partial_j q_a) \right)^2
$$

The Skyrme coupling constant:

$$
c_4 = \frac{2\rho_0^2}{e^2}
$$

In code units ($e = 1$, $\rho_0 = 1$): $c_4 = 2$.

### 2.3 Kinetic, Dirichlet, and Skyrme Energies

The total energy is $E = E_\mathrm{kin} + E_2 + E_4$ with:

$$
E_\mathrm{kin} = \frac{1}{2} \int |\dot{q}|^2 \, d^3x
  = \frac{1}{2} \int \sum_{a=0}^{3} \dot{q}_a^2 \, d^3x
$$

$$
E_2 = \frac{1}{2} \int |\nabla q|^2 \, d^3x
  = \frac{1}{2} \int \mathrm{Tr}\,D \, d^3x
$$

$$
E_4 = \frac{c_4}{4} \int \left[ (\mathrm{Tr}\,D)^2 - \mathrm{Tr}(D^2) \right] d^3x
$$

The Skyrme energy density can be rewritten using the Maurer-Cartan currents
$L_i^a = \mathrm{Im}(\bar{q} \cdot \partial_i q)_a$:

$$
\mathcal{E}_4 = \frac{c_4}{4} \sum_{i<j} |[L_i, L_j]|^2
$$

where $[L_i, L_j]^a = \epsilon^{abc} L_i^b L_j^c$ is the su(2) commutator.

### 2.4 Equations of Motion

Varying the Lagrangian subject to $|q| = 1$ (Lagrange multiplier λ):

$$
\ddot{q}_a = \nabla^2 q_a + c_4 (\nabla \cdot G)_a - \lambda \, q_a
$$

where the Skyrme stress vector $G$ is:

$$
G_i^a = \mathrm{Tr}(D) \, \partial_i q_a - \sum_{j=1}^{3} D_{ij} \, \partial_j q_a
$$

The constraint multiplier absorbs radial acceleration:

$$
\lambda = q \cdot F + |\dot{q}|^2
$$

where $F_a = \nabla^2 q_a + c_4 (\nabla \cdot G)_a$ is the unconstrained force. This gives
the projected EOM:

$$
\ddot{q}_a = F_a - (q \cdot F + |\dot{q}|^2) \, q_a
$$

Written out for each quaternion component $a \in \{0,1,2,3\}$:

$$
\ddot{q}_a = \underbrace{\nabla^2 q_a}_\text{wave}
  + \underbrace{c_4 \sum_{i=1}^{3} \partial_i G_i^a}_\text{Skyrme}
  - \underbrace{\left(\sum_{b=0}^{3} q_b F_b + \sum_{b=0}^{3} \dot{q}_b^2 \right) q_a}_\text{constraint}
$$

### 2.5 Divergence of Skyrme Stress (Explicit)

The divergence of $G$ required for the force:

$$
(\nabla \cdot G)_a = \sum_{i=1}^{3} \partial_i G_i^a
$$

Each $G_i^a$ involves first derivatives of $q$, so $\nabla \cdot G$ involves second
derivatives. On the lattice, we compute $G_i^a$ at all grid points first, then take
the divergence by finite differences (two-pass method, as in v18).

### 2.6 Topological Charge (Baryon Number)

$$
B = -\frac{1}{2\pi^2} \int \det(L_x, L_y, L_z) \, d^3x
$$

where

$$
\det(L_x, L_y, L_z) = L_x \cdot (L_y \times L_z)
  = \sum_{a,b,c} \epsilon_{abc} L_x^a L_y^b L_z^c
$$

For our three-mode resonance, $B = 0$ at initialization and must remain zero
throughout. Any drift signals numerical error.

### 2.7 Linearized Modes Around Vacuum

For small perturbations $q = (1 - \frac{1}{2}|\varepsilon|^2, \varepsilon_1, \varepsilon_2, \varepsilon_3)$
around vacuum, the EOM reduces to:

$$
\ddot{\varepsilon}_a = \nabla^2 \varepsilon_a + O(\varepsilon^3)
$$

This is three copies of the massless wave equation. Solutions:

$$
\varepsilon_a(\mathbf{x}, t) = A_a \, f(r) \, Y_l^m(\hat{x}) \, \cos(\omega t)
$$

where $f(r)$ satisfies the radial wave equation:

$$
f''(r) + \frac{2}{r} f'(r) + \left(\omega^2 - \frac{l(l+1)}{r^2}\right) f(r) = 0
$$

with solutions $f(r) = j_l(\omega r)$ (spherical Bessel functions). For $l = 0$:

$$
f(r) = j_0(\omega r) = \frac{\sin(\omega r)}{\omega r}
$$

For $l = 1$:

$$
f(r) = j_1(\omega r) = \frac{\sin(\omega r)}{(\omega r)^2} - \frac{\cos(\omega r)}{\omega r}
$$

These are the building blocks for our initial conditions, modulated by a Gaussian
envelope to ensure localization.

### 2.8 Nonlinear Coupling Mechanism

At cubic order in $\varepsilon$, the S³ constraint generates coupling:

$$
q_0 = \sqrt{1 - \varepsilon_1^2 - \varepsilon_2^2 - \varepsilon_3^2}
  \approx 1 - \frac{1}{2}(\varepsilon_1^2 + \varepsilon_2^2 + \varepsilon_3^2)
  - \frac{1}{8}(\varepsilon_1^2 + \varepsilon_2^2 + \varepsilon_3^2)^2 - \cdots
$$

Substituting into the EOM for $\varepsilon_a$:

$$
\ddot{\varepsilon}_a = \nabla^2 \varepsilon_a
  - |\varepsilon|^2 \nabla^2 \varepsilon_a
  + (\varepsilon \cdot \nabla^2\varepsilon) \varepsilon_a
  - |\dot{\varepsilon}|^2 \varepsilon_a
  + O(\varepsilon^5)
$$

The cubic terms couple the three modes. The key cross-terms are:

- $\varepsilon_b^2 \nabla^2 \varepsilon_a$ (mode $b$ modulates the dispersion of mode $a$)
- $(\varepsilon_b \nabla^2 \varepsilon_b) \varepsilon_a$ (mode $b$'s Laplacian drives mode $a$)
- $\dot{\varepsilon}_b^2 \varepsilon_a$ (mode $b$'s kinetic energy acts as a potential for mode $a$)

For three-wave resonance with $\omega_1 + \omega_2 = \omega_3$, the product
$\varepsilon_1 \varepsilon_2$ oscillates at $\omega_1 + \omega_2 = \omega_3$, which
resonantly drives $\varepsilon_3$. Simultaneously, $\varepsilon_3 \varepsilon_1^*$
drives $\varepsilon_2$, and $\varepsilon_3 \varepsilon_2^*$ drives $\varepsilon_1$.
This closed triad is the standard three-wave parametric resonance.

The Skyrme term adds additional quartic coupling (at order $\varepsilon^3$ in the force)
that further entangles the modes but does not change the resonance structure.

### 2.9 Time-Averaged Static Deformation (DC Rectification)

When three modes oscillate, the constraint generates a **time-independent** (DC)
correction to $q_0$:

$$
\langle q_0(x) \rangle_t = 1 - \frac{1}{2} \sum_{a=1}^{3} \langle \varepsilon_a^2(x,t) \rangle_t
$$

For $\varepsilon_a = A_a(r) \cos(\omega_a t + \phi_a)$:

$$
\langle \varepsilon_a^2 \rangle_t = \frac{1}{2} A_a(r)^2
$$

So:

$$
\langle q_0 \rangle_t - 1 = -\frac{1}{4} \sum_{a=1}^{3} A_a(r)^2
$$

This DC deformation is the "shape distortion" visible in the far field. If the mode
amplitudes $A_a(r)$ fall as $r^{-n}$, then the DC deformation falls as $r^{-2n}$.

**The critical question**: does the nonlinear phase locking, the S³ constraint, or the
Skyrme term modify the mode amplitude falloff from the linear expectation? If $A_a(r)$
falls slower than $r^{-1}$ (the free-wave 1/r radiation falloff), the DC deformation
could fall slower than $r^{-2}$ — potentially as slow as $r^{-1}$.

This is what the numerical experiment must determine.

---

## 3. Phase 1: B=0 Oscillon on S³ + Skyrme

### 3.1 Experiment 1A — Spherically Symmetric Bump (Control)

Initialize a non-topological Gaussian bump around vacuum:

$$
q_a(\mathbf{x}, 0) = A \, e^{-r^2 / 2\sigma^2} \cdot \frac{x_a}{r} \quad (a = 1,2,3)
$$

$$
q_0(\mathbf{x}, 0) = \sqrt{1 - q_1^2 - q_2^2 - q_3^2}
$$

$$
\dot{q}_a(\mathbf{x}, 0) = 0 \quad (a = 0,1,2,3)
$$

where:
- $A = 0.8$ (strong nonlinearity; $q_0(0) = \sqrt{1 - 3 \cdot 0.64} \approx 0.33$)
- $\sigma = 1.5$ (comparable to Skyrmion core size)
- $r = |\mathbf{x}| = \sqrt{x^2 + y^2 + z^2}$
- At $r = 0$: use L'Hopital, $x_a/r$ is replaced by the unit vector $\hat{x}_a$;
  in practice set $q_a(0,0) = A/\sqrt{3}$ for $a=1,2,3$

This is a **hedgehog-like** bump but with B = 0 because the profile $A e^{-r^2/2\sigma^2}$
does NOT wind from $\pi$ to $0$. The field starts at $q_0 \approx 0.33$ at the origin and
returns to $q_0 = 1$ (vacuum) at large $r$, never wrapping S³.

**Purpose**: Establish the baseline dispersal timescale $T_\mathrm{disp}$ for a generic
B=0 perturbation. Expect rapid radiation to the boundary.

### 3.2 Experiment 1B — Three-Mode Kicked Resonance

Initialize at vacuum with three orthogonal velocity kicks in different quaternion
components, each with a distinct radial structure:

$$
q(\mathbf{x}, 0) = (1, 0, 0, 0) \quad \text{(vacuum everywhere)}
$$

$$
\dot{q}_1(\mathbf{x}, 0) = \omega_1 A \, g_1(r) \, \frac{x}{r}
$$

$$
\dot{q}_2(\mathbf{x}, 0) = \omega_2 A \, g_2(r) \, \frac{y}{r}
$$

$$
\dot{q}_3(\mathbf{x}, 0) = \omega_3 A \, g_3(r) \, \frac{z}{r}
$$

$$
\dot{q}_0(\mathbf{x}, 0) = -q_0^{-1} \sum_{a=1}^{3} q_a \dot{q}_a = 0
  \quad \text{(tangent to S³; trivially zero since $q_a = 0$)}
$$

The radial envelopes have different numbers of nodes:

$$
g_n(r) = j_0(n\pi r / R_0) \cdot e^{-r^2 / 2\sigma^2}
  = \frac{\sin(n\pi r / R_0)}{n\pi r / R_0} \cdot e^{-r^2 / 2\sigma^2}
$$

where $R_0$ is the effective confinement radius (set by the Gaussian). Specifically:

- Mode 1 ($q_1$, along $\hat{x}$): $g_1(r)$ with $n = 1$, frequency $\omega_1 = \pi c / R_0$
- Mode 2 ($q_2$, along $\hat{y}$): $g_2(r)$ with $n = 2$, frequency $\omega_2 = 2\pi c / R_0$
- Mode 3 ($q_3$, along $\hat{z}$): $g_3(r)$ with $n = 3$, frequency $\omega_3 = 3\pi c / R_0$

This satisfies the three-wave resonance condition:

$$
\omega_1 + \omega_2 = \omega_3
$$

Parameters:
- $A = 0.5$ (moderate kick amplitude)
- $\sigma = 2.0$ (Gaussian cutoff)
- $R_0 = 3.0$ (sets the frequency scale; $\omega_1 = 1.047$, $\omega_2 = 2.094$, $\omega_3 = 3.142$)
- $c = 1$ (code units)

**Purpose**: The three modes start with no displacement but with kinetic energy. The
S³ nonlinearity and Skyrme coupling must transfer energy between modes. If the triad
locks, the localized oscillation persists; if not, it radiates away.

### 3.3 Experiment 1C — Large-Amplitude Breathing Triad

Same as 1B but with displacement initialization (zero velocity):

$$
q_1(\mathbf{x}, 0) = A_1 \, g_1(r) \, \frac{x}{r}
$$

$$
q_2(\mathbf{x}, 0) = A_2 \, g_2(r) \, \frac{y}{r}
$$

$$
q_3(\mathbf{x}, 0) = A_3 \, g_3(r) \, \frac{z}{r}
$$

$$
q_0(\mathbf{x}, 0) = \sqrt{1 - q_1^2 - q_2^2 - q_3^2}
$$

$$
\dot{q}(\mathbf{x}, 0) = 0
$$

where:
- $A_1 = 0.4$, $A_2 = 0.4$, $A_3 = 0.56$ (mode 3 gets the combined energy of 1+2)
- $g_n(r)$ same as Experiment 1B
- At $r = 0$: replace $x_a/r$ with $(1,0,0)$, $(0,1,0)$, $(0,0,1)$ respectively

The asymmetric amplitudes ($A_3 > A_1, A_2$) seed the triad resonance by pre-loading
mode 3 with the energy it would receive from the $\omega_1 + \omega_2 = \omega_3$
coupling.

### 3.4 Experiment 1D — Single-Mode Control

Same setup as 1C but with **only one mode** active:

$$
q_1(\mathbf{x}, 0) = A \, g_1(r) \, \frac{x}{r}, \quad q_2 = q_3 = 0
$$

$$
q_0(\mathbf{x}, 0) = \sqrt{1 - q_1^2}
$$

**Purpose**: Direct comparison with 1C. If the three-mode version persists significantly
longer than this single-mode version, the phase locking is real, not an artifact.

---

## 4. Numerical Method

### 4.1 Grid

3D uniform Cartesian grid on $[-L/2, L/2]^3$:

$$
x_i = -\frac{L}{2} + i \cdot h, \quad h = \frac{L}{N-1}, \quad i = 0, \ldots, N-1
$$

Volume element: $dV = h^3$.

Parameters for Phase 1:
- $N = 200$ (200³ grid, ~64M points)
- $L = 16$ (box large enough for far-field measurement)
- $h = L/(N-1) = 0.0804$
- Points across the core ($\sigma = 2.0$): $2\sigma/h \approx 50$ (excellent resolution)

### 4.2 Finite Differences

Spatial first derivatives (central, 2nd order):

$$
\partial_i q_a \Big|_p = \frac{q_a(p + \hat{e}_i) - q_a(p - \hat{e}_i)}{2h}
$$

7-point Laplacian (2nd order):

$$
\nabla^2 q_a \Big|_p = \frac{1}{h^2} \left[
  \sum_{d=1}^{3} \left( q_a(p + \hat{e}_d) + q_a(p - \hat{e}_d) \right) - 6 \, q_a(p)
\right]
$$

Divergence of Skyrme stress (two-pass):
1. Compute $G_i^a(p)$ at every grid point from $\partial_j q$ and $D_{ij}$
2. Compute $(\nabla \cdot G)_a = \sum_i (G_i^a(p+\hat{e}_i) - G_i^a(p-\hat{e}_i)) / (2h)$

### 4.3 Time Integration — Velocity Verlet

State: position $q_a(t)$, velocity $v_a(t) = \dot{q}_a(t)$.

At each timestep $\Delta t$:

**Step 1 — Half-kick:**

$$
v_a^{n+1/2} = v_a^n + \frac{\Delta t}{2} \, a_a^n
$$

where $a_a^n = F_a^n - (q^n \cdot F^n + |v^n|^2) q_a^n$ is the constrained acceleration.

**Step 2 — Drift:**

$$
q_a^{n+1} = q_a^n + \Delta t \, v_a^{n+1/2}
$$

**Step 3 — Project back to S³:**

$$
q^{n+1} \leftarrow \frac{q^{n+1}}{|q^{n+1}|}
$$

$$
v^{n+1/2} \leftarrow v^{n+1/2} - (q^{n+1} \cdot v^{n+1/2}) \, q^{n+1}
$$

**Step 4 — Second half-kick at new position:**

Recompute $a^{n+1}$ from $q^{n+1}$, then:

$$
v_a^{n+1} = v_a^{n+1/2} + \frac{\Delta t}{2} \, a_a^{n+1}
$$

**Step 5 — Absorbing boundary damping:**

$$
v_a^{n+1} \leftarrow \alpha(\mathbf{x}) \, v_a^{n+1}
$$

where $\alpha(\mathbf{x}) = 1$ in the interior and ramps to 0.05 in the outer 12% of
the box (see Section 4.5).

### 4.4 CFL Condition

The maximum wave speed is $c = 1$ (code units). With the Skyrme term, local wave speeds
can be higher. The CFL stability condition for the 7-point Laplacian:

$$
\Delta t \leq \frac{h}{\sqrt{3} \cdot c_\mathrm{max}}
$$

For safety with the nonlinear Skyrme term, use:

$$
\Delta t = 0.25 \times \frac{h}{\sqrt{3}}
$$

With $h = 0.0804$: $\Delta t = 0.0116$.

### 4.5 Absorbing Boundary Layer

To prevent artificial reflections, apply a velocity damping mask $\alpha(\mathbf{x})$:

$$
\alpha(\mathbf{x}) = 1 - 0.95 \, (1 - m(\mathbf{x}))
$$

where $m(\mathbf{x}) = \prod_{d=1}^{3} m_d$ and:

$$
m_d = \begin{cases}
  x_d / b & \text{if } x_d < b \\
  (N - 1 - x_d) / b & \text{if } x_d > N - 1 - b \\
  1 & \text{otherwise}
\end{cases}
$$

with $b = \lfloor 0.12 N \rfloor$ grid points in the absorbing layer.

### 4.6 Periodic vs Absorbing Boundaries

Use **absorbing boundaries** (not periodic). Periodic boundaries create image solitons
that contaminate the far-field measurement. The absorbing layer lets radiation escape.

---

## 5. Diagnostics

### 5.1 Energy Conservation

Compute at each diagnostic step:

$$
E_\mathrm{kin} = \frac{h^3}{2} \sum_p \sum_{a=0}^{3} v_a(p)^2
$$

$$
E_2 = \frac{h^3}{2} \sum_p \mathrm{Tr}\,D(p)
$$

$$
E_4 = \frac{c_4 h^3}{4} \sum_p \left[ (\mathrm{Tr}\,D(p))^2 - \mathrm{Tr}(D(p)^2) \right]
$$

$$
E_\mathrm{total} = E_\mathrm{kin} + E_2 + E_4
$$

Energy should decrease monotonically (radiation absorbed at boundary). The rate of
decrease is the radiation power $P_\mathrm{rad} = -dE/dt$.

### 5.2 Topological Charge

$$
B = -\frac{h^3}{2\pi^2} \sum_p \det(L_x, L_y, L_z)\Big|_p
$$

Must remain $|B| < 0.01$ throughout. Any larger value indicates the resonance has
spontaneously nucleated topology (interesting but different from what we're testing).

### 5.3 Localization — Core Energy Fraction

Define the core as the sphere $r < R_c$ (with $R_c = 2\sigma$):

$$
f_\mathrm{core}(t) = \frac{\sum_{|x_p| < R_c} \mathcal{E}(p) \, h^3}{E_\mathrm{total}(t)}
$$

where $\mathcal{E}(p) = \frac{1}{2}|v(p)|^2 + \frac{1}{2}\mathrm{Tr}\,D(p) + \frac{c_4}{4}[(\mathrm{Tr}\,D)^2 - \mathrm{Tr}(D^2)]$.

If $f_\mathrm{core}(t)$ stays above 0.5 for time $T$, and $T > 10 T_\mathrm{disp}$
(from the single-mode control), the resonance is "persistent."

### 5.4 Mode Decomposition

Project the field onto the three initial mode shapes to track individual amplitudes:

$$
c_a(t) = h^3 \sum_p q_a(p,t) \cdot g_a(r_p) \cdot \frac{(x_a)_p}{r_p}
$$

where $g_a$ is the initial radial profile of mode $a$, and $(x_1, x_2, x_3) = (x, y, z)$.

Then extract amplitude and phase:

$$
A_a(t) = \sqrt{c_a(t)^2 + (c_a'(t)/\omega_a)^2}
$$

where $c_a'(t)$ is computed from finite differences of $c_a(t)$.

**Phase locking diagnostic**: compute:

$$
\Phi(t) = \phi_3(t) - \phi_1(t) - \phi_2(t)
$$

where $\phi_a(t) = \mathrm{atan2}(-c_a'(t)/\omega_a, \, c_a(t))$. If $\Phi(t)$ remains
approximately constant (within $\pm \pi/4$) for time $T > 10/\omega_1$, the triad is
phase-locked.

### 5.5 Far-Field Structure — Radial Shell Averages

Compute shell-averaged quantities as functions of radius:

$$
\langle Q \rangle_\Omega(r) = \frac{1}{N_\mathrm{shell}(r)} \sum_{|r_p - r| < h/2} Q(p)
$$

for the following quantities $Q$:

**Field deviation (instantaneous):**

$$
\delta(p) = 1 - q_0(p) = 1 - q_0(p)
$$

Note: $\delta \geq 0$ since $q_0 \leq 1$, and $\delta = 0$ at vacuum.

For small perturbations: $\delta \approx \frac{1}{2}(q_1^2 + q_2^2 + q_3^2)$.

**Energy density:**

$$
\rho_E(p) = \frac{1}{2}|v|^2 + \frac{1}{2}\mathrm{Tr}\,D + \frac{c_4}{4}[(\mathrm{Tr}\,D)^2 - \mathrm{Tr}(D^2)]
$$

**Maurer-Cartan current magnitude:**

$$
|L|(p) = \sqrt{\sum_{i,a} (L_i^a)^2}
$$

### 5.6 Far-Field Structure — Time-Averaged

The key measurement. After the resonance has settled (say $t > T_\mathrm{settle}$),
compute time-averaged shell profiles:

$$
\overline{\delta}(r) = \frac{1}{T_\mathrm{avg}} \int_{T_\mathrm{settle}}^{T_\mathrm{settle}+T_\mathrm{avg}}
  \langle 1 - q_0 \rangle_\Omega(r, t) \, dt
$$

$$
\overline{\rho_E}(r) = \frac{1}{T_\mathrm{avg}} \int_{T_\mathrm{settle}}^{T_\mathrm{settle}+T_\mathrm{avg}}
  \langle \rho_E \rangle_\Omega(r, t) \, dt
$$

In practice, accumulate running sums every $M$ timesteps.

**Power-law fit**: in the region $r \in [r_\mathrm{min}, r_\mathrm{max}]$ where
$r_\mathrm{min} = 2R_c$ (outside core) and $r_\mathrm{max} = L/2 - 2b \cdot h$
(inside absorbing layer), fit:

$$
\log \overline{\delta}(r) = n_\delta \log r + C_\delta
$$

$$
\log \overline{\rho_E}(r) = n_E \log r + C_E
$$

using linear regression. Report the exponents $n_\delta$ and $n_E$ with $R^2$ values.

**Expected results**:
- Free radiation (no binding): $n_\delta = -2$ (energy spreads as $1/r^2$)
- Static Skyrmion far field: $n_E \sim -8$ (from $f \sim 1/r^2$ profile)
- **Interesting**: $n_\delta = -1$ would mean 1/r DC deformation
- **Very interesting**: $n_E = -2$ would mean Coulomb-like energy density

### 5.7 Per-Mode Energy Tracking

Decompose the total energy into contributions from each quaternion component to detect
energy sloshing between the three modes. Define per-mode energies:

**Mode kinetic energy:**

$$
E_{\mathrm{kin},a}(t) = \frac{h^3}{2} \sum_p \dot{q}_a(p,t)^2 \quad (a = 1,2,3)
$$

**Mode gradient energy (approximate — ignores cross-terms from constraint):**

$$
E_{2,a}(t) = \frac{h^3}{2} \sum_p |\nabla q_a(p,t)|^2 \quad (a = 1,2,3)
$$

**Total per-mode energy:**

$$
E_a(t) = E_{\mathrm{kin},a}(t) + E_{2,a}(t)
$$

Note: these do not sum exactly to $E_\mathrm{kin} + E_2$ because the S³ constraint
couples $q_0$ to $(q_1, q_2, q_3)$, and the Skyrme energy $E_4$ involves cross-terms
between all components. But for small-to-moderate amplitudes, $E_1 + E_2 + E_3$ is a
good proxy for the total non-vacuum energy.

**Diagnostic**: plot $E_1(t)$, $E_2(t)$, $E_3(t)$ on the same axis. Look for:

- **Periodic energy exchange**: $E_1$ and $E_2$ decrease while $E_3$ increases (and
  vice versa) with period $T_\mathrm{slosh} \sim 2\pi/|\omega_3 - \omega_1 - \omega_2|$.
  This is the hallmark of strong triad coupling (analogous to Manley-Rowe relations
  in three-wave mixing).
- **Monotonic decay**: all three decrease together → no coupling, just radiation.
- **Equipartition**: all three fluctuate around $E_\mathrm{total}/3$ with no coherent
  exchange → thermalized, not locked.

The sloshing period and amplitude directly measure the nonlinear coupling strength.

### 5.8 Radiation Power Spectrum

To check for frequency locking, compute the power spectrum of $q_1(0,t)$ (field at origin):

$$
\hat{q}_1(\omega) = \int_0^{T} q_1(0,t) \, e^{-i\omega t} \, dt
$$

If the triad is locked, we should see sharp peaks at $\omega_1$, $\omega_2$, $\omega_3$
with $\omega_1 + \omega_2 = \omega_3$ (within numerical broadening). If not locked,
the spectrum will be broad or have incommensurate peaks.

### 5.9 Gravity Proxy — Passive Poisson Solution

Given the time-averaged energy density $\overline{\rho_E}(r)$, solve the Poisson equation:

$$
\nabla^2 \Phi = 4\pi G \, \overline{\rho_E}
$$

For a spherically-averaged source, this reduces to the 1D ODE:

$$
\frac{1}{r^2} \frac{d}{dr}\left(r^2 \frac{d\Phi}{dr}\right) = 4\pi G \, \overline{\rho_E}(r)
$$

Integrating once:

$$
\frac{d\Phi}{dr}(r) = \frac{G \, M(r)}{r^2}, \quad M(r) = 4\pi \int_0^r \overline{\rho_E}(r') \, r'^2 \, dr'
$$

This will **always** give $1/r^2$ force at large $r$ (by Gauss's theorem) as long as
$M(r) \to M_\mathrm{total}$ converges. This test is a **sanity check**, not the main
result. The interesting physics is in the power-law exponents of Section 5.6.

---

## 6. Phase 2: Conserved Density Coupling (Conditional)

**Only proceed to Phase 2 if Phase 1 produces a persistent resonance** (Experiments 1B
or 1C with lifetime > 10× control).

### 6.1 Extended Field Content

Add an independent density field $\rho(\mathbf{x}, t)$ satisfying:

- **State**: $(q, \rho)$ where $q \in S^3$ and $\rho \geq 0$
- **Conservation**: $\partial_t \rho + \nabla \cdot (\rho \, \mathbf{v}) = 0$

### 6.2 Extended Lagrangian (v6 Architecture)

$$
\mathcal{L} = \frac{1}{2c_s^2} \dot{\rho}^2
  + \frac{\rho}{2} |\dot{q}|^2
  - \frac{\rho}{2} |\nabla q|^2
  - \frac{c_4}{4} \left[(\mathrm{Tr}\,D)^2 - \mathrm{Tr}(D^2)\right]
  - \frac{\alpha}{2}(\rho - \rho_0)^2
$$

Wait — the potential $\frac{\alpha}{2}(\rho - \rho_0)^2$ makes $\rho$ **massive**
(Yukawa, not 1/r). For a massless density mode we need **no potential**:

$$
\mathcal{L}_\rho = \frac{1}{2c_s^2} \dot{\rho}^2 - \frac{1}{2}|\nabla\rho|^2
$$

$$
\mathcal{L}_q = \frac{\rho}{2} |\dot{q}|^2 - \frac{\rho}{2}|\nabla q|^2
  - \frac{c_4}{4}\left[(\mathrm{Tr}\,D)^2 - \mathrm{Tr}(D^2)\right]
$$

$$
\mathcal{L} = \mathcal{L}_\rho + \mathcal{L}_q
$$

The coupling is through the $\rho \cdot |\nabla q|^2$ term: the twist energy depends
linearly on density.

### 6.3 Equations of Motion (Phase 2)

**Density EOM** (varying $\mathcal{L}$ with respect to $\rho$):

$$
\frac{1}{c_s^2}\ddot{\rho} = \nabla^2\rho - \frac{1}{2}|\dot{q}|^2 + \frac{1}{2}|\nabla q|^2
$$

Wait — this gives $\frac{1}{c_s^2}\ddot{\rho} = \nabla^2\rho + S(x,t)$ where the
source is:

$$
S(x,t) = -\frac{1}{2}|\dot{q}|^2 + \frac{1}{2}|\nabla q|^2
$$

For a stationary configuration ($\dot{q} = 0$), $S = +\frac{1}{2}|\nabla q|^2 > 0$,
so the source is positive. But this means $\rho$ is *repelled* from twists, not
attracted. That's the depletion mechanism: density flows AWAY from the soliton core,
creating a deficit.

Actually, to get the signs right we need to be more careful. The full variation:

$$
\frac{\partial \mathcal{L}}{\partial \rho} = \frac{1}{2}|\dot{q}|^2 - \frac{1}{2}|\nabla q|^2
$$

$$
\frac{\partial \mathcal{L}}{\partial \dot{\rho}} = \frac{\dot{\rho}}{c_s^2}
$$

$$
\frac{\partial \mathcal{L}}{\partial (\partial_i \rho)} = -\partial_i \rho
$$

Euler-Lagrange:

$$
\frac{1}{c_s^2}\ddot{\rho} - \nabla^2\rho + \frac{1}{2}|\nabla q|^2 - \frac{1}{2}|\dot{q}|^2 = 0
$$

Rearranging:

$$
\frac{1}{c_s^2}\ddot{\rho} = \nabla^2\rho - \frac{1}{2}|\nabla q|^2 + \frac{1}{2}|\dot{q}|^2
$$

For a stationary soliton ($\dot{q} = 0$):

$$
\frac{1}{c_s^2}\ddot{\rho} = \nabla^2\rho - \frac{1}{2}|\nabla q|^2
$$

Source term is $-\frac{1}{2}|\nabla q|^2 < 0$. In the static limit ($\ddot{\rho} = 0$):

$$
\nabla^2 \delta\rho = \frac{1}{2}|\nabla q|^2 = \frac{1}{2}\mathrm{Tr}\,D
$$

This is a Poisson equation with positive source ∝ twist energy density. By the divergence
theorem:

$$
\delta\rho(r) = \frac{1}{4\pi r} \int \frac{1}{2}|\nabla q|^2 \, d^3x' + \cdots
  = \frac{E_2}{4\pi r} + \cdots
$$

Wait — this gives $\delta\rho > 0$, which means density INCREASES near the soliton.
Let me re-examine the sign.

The background is $\rho_0$ (constant). Set $\rho = \rho_0 + \delta\rho$. The density EOM:

$$
\frac{1}{c_s^2}\ddot{\delta\rho} = \nabla^2\delta\rho - \frac{1}{2}|\nabla q|^2 + \frac{1}{2}|\dot{q}|^2
$$

Static limit: $\nabla^2 \delta\rho = +\frac{1}{2}|\nabla q|^2 > 0$.

This means $\delta\rho$ is POSITIVE everywhere (density piles up near the soliton, not
depletes). That's because the twist energy $\rho|\nabla q|^2$ means the system GAINS
energy by increasing $\rho$ near twists — but that costs gradient energy in $\rho$.

For **depletion** (the desired gravitational mechanism), we need the coupling to have
the opposite sign. The v6 Lagrangian used $E_2 = \frac{1}{2}\rho|\nabla q|^2$, which
means more density = more energy = system wants to deplete. The static condition is:

$$
\nabla^2\delta\rho = -\frac{1}{2}|\nabla q|^2
$$

giving $\delta\rho < 0$ (depletion). The sign depends on the field equation derivation.

Let me re-derive carefully. The v6 energy functional:

$$
E = \int \left[ \frac{1}{2}|\nabla\rho|^2 + \frac{\rho}{2}|\nabla q|^2 + \frac{c_4}{4}[\cdots] \right] d^3x
$$

For a static configuration, minimize E with respect to $\rho$:

$$
\frac{\delta E}{\delta \rho} = -\nabla^2\rho + \frac{1}{2}|\nabla q|^2 = 0
$$

$$
\nabla^2\rho = +\frac{1}{2}|\nabla q|^2
$$

So $\nabla^2\rho > 0$ everywhere, meaning $\rho$ has a MINIMUM at the soliton core
(where $|\nabla q|^2$ is largest). With boundary condition $\rho \to \rho_0$ at infinity:

$$
\rho(r) = \rho_0 - \frac{1}{4\pi r} \int \frac{1}{2}|\nabla q(\mathbf{x}')|^2 \, d^3x' + \cdots
  = \rho_0 - \frac{E_2}{4\pi r}
$$

Yes — depletion ($\rho < \rho_0$) falling as $1/r$. The sign is correct when derived
from the energy minimum, not from the wave equation directly (the wave equation's
static limit gives $\nabla^2\delta\rho = -\text{source}$ because the kinetic term
enters with a sign flip).

Corrected density wave equation (from the proper Lagrangian with $\rho$ kinetic term
$+\frac{1}{2c_s^2}\dot{\rho}^2$):

$$
\frac{1}{c_s^2}\ddot{\rho} + \nabla^2\rho = \frac{1}{2}|\nabla q|^2 - \frac{1}{2}|\dot{q}|^2
$$

**Note the sign of $\nabla^2\rho$**: the spatial kinetic energy of $\rho$ is
$-\frac{1}{2}|\nabla\rho|^2$ in the Lagrangian (standard convention), giving
$+\nabla^2\rho$ in the EOM.

Static limit: $\nabla^2\rho = +\frac{1}{2}|\nabla q|^2$, consistent with the energy
minimization above.

### 6.4 Phase 2 Corrected Equations

Lagrangian density:

$$
\mathcal{L} = \frac{1}{2c_s^2}\dot{\rho}^2 - \frac{1}{2}|\nabla\rho|^2
  + \frac{\rho}{2}|\dot{q}|^2 - \frac{\rho}{2}|\nabla q|^2
  - \frac{c_4}{4}\left[(\mathrm{Tr}\,D)^2 - \mathrm{Tr}(D^2)\right]
$$

Density EOM:

$$
\boxed{\frac{1}{c_s^2}\ddot{\rho} = -\nabla^2\rho + \frac{1}{2}|\dot{q}|^2 - \frac{1}{2}|\nabla q|^2}
$$

Actually — I need to be careful with the sign convention. Let me redo this from scratch.

The Lagrangian density for the density sector is:

$$
\mathcal{L}_\rho = \frac{1}{2c_s^2}\dot{\rho}^2 - \frac{1}{2}(\nabla\rho)^2
$$

The coupling term is:

$$
\mathcal{L}_\mathrm{coupling} = \frac{\rho}{2}|\dot{q}|^2 - \frac{\rho}{2}|\nabla q|^2
$$

Euler-Lagrange for $\rho$:

$$
\frac{d}{dt}\frac{\partial\mathcal{L}}{\partial\dot{\rho}}
  - \frac{\partial\mathcal{L}}{\partial\rho}
  + \sum_i \frac{\partial}{\partial x_i}\frac{\partial\mathcal{L}}{\partial(\partial_i\rho)}
  = 0
$$

$$
\frac{\ddot{\rho}}{c_s^2}
  - \frac{1}{2}|\dot{q}|^2 + \frac{1}{2}|\nabla q|^2
  + \nabla^2\rho = 0
$$

Therefore:

$$
\boxed{\frac{1}{c_s^2}\ddot{\rho} = -\nabla^2\rho + \frac{1}{2}|\dot{q}|^2 - \frac{1}{2}|\nabla q|^2}
$$

Static limit ($\ddot{\rho} = 0$, $\dot{q} = 0$):

$$
\nabla^2\rho = \frac{1}{2}|\nabla q|^2 > 0
$$

With $\rho = \rho_0 + \delta\rho$ and $\delta\rho \to 0$ at infinity, the Green's
function solution is:

$$
\delta\rho(\mathbf{x}) = -\frac{1}{4\pi}\int \frac{|\nabla q(\mathbf{x}')|^2 / 2}{|\mathbf{x} - \mathbf{x}'|} \, d^3x'
$$

At large $r$:

$$
\delta\rho(r) \approx -\frac{E_2}{4\pi r}
$$

This is **negative** (depletion) and falls as **1/r**. The total depletion deficit:

$$
\Delta Q = \int \delta\rho \, d^3x = -\frac{E_2}{4\pi} \cdot 4\pi \cdot \int_0^\infty \frac{1}{r} \cdot r^2 \, dr
$$

which diverges — reflecting the long-range nature of the depletion. In practice, the
density field equilibrates through propagation at speed $c_s$, and the depletion
wavefront expands at $c_s$.

### 6.5 Quaternion EOM (Phase 2)

The coupling also modifies the quaternion dynamics. The $\rho$-dependent terms:

$$
\frac{\partial\mathcal{L}}{\partial q_a} = 0 \quad \text{(no explicit $q$ dependence in coupling)}
$$

$$
\frac{\partial\mathcal{L}}{\partial\dot{q}_a} = \rho \dot{q}_a
$$

$$
\frac{\partial\mathcal{L}}{\partial(\partial_i q_a)} = -\rho \partial_i q_a + c_4 G_i^a
$$

Euler-Lagrange:

$$
\frac{d}{dt}(\rho\dot{q}_a) + \nabla\cdot(\rho\nabla q_a) + c_4(\nabla\cdot G)_a = \lambda q_a
$$

Expanding the time derivative:

$$
\rho\ddot{q}_a + \dot{\rho}\dot{q}_a = \nabla\cdot(\rho\nabla q_a) + c_4(\nabla\cdot G)_a - \lambda q_a
$$

$$
\rho\ddot{q}_a = \nabla\cdot(\rho\nabla q_a) + c_4(\nabla\cdot G)_a - \dot{\rho}\dot{q}_a - \lambda q_a
$$

where $\nabla\cdot(\rho\nabla q_a) = \rho\nabla^2 q_a + \nabla\rho \cdot \nabla q_a$.

So:

$$
\boxed{\ddot{q}_a = \nabla^2 q_a + \frac{\nabla\rho}{\rho}\cdot\nabla q_a
  + \frac{c_4}{\rho}(\nabla\cdot G)_a - \frac{\dot{\rho}}{\rho}\dot{q}_a
  - \frac{\lambda}{\rho} q_a}
$$

The new terms compared to Phase 1 are:
- $(\nabla\rho/\rho)\cdot\nabla q_a$: density gradient bends the wave (refraction)
- $-(\dot{\rho}/\rho)\dot{q}_a$: time-varying density modulates the amplitude (parametric)

### 6.6 Phase 2 Parameters

- $c_s = 1$ (density wave speed = twist wave speed, simplest case)
- $\rho_0 = 1$ (background density)
- Initialize: $\rho(\mathbf{x}, 0) = \rho_0 = 1$ everywhere, $\dot{\rho}(\mathbf{x}, 0) = 0$
- Initialize $q$ as in Phase 1 Experiments 1B or 1C

### 6.7 Phase 2 Diagnostics

All Phase 1 diagnostics, plus:

**Density depletion profile:**

$$
\overline{\delta\rho}(r) = \frac{1}{T_\mathrm{avg}} \int \langle \rho - \rho_0 \rangle_\Omega(r,t) \, dt
$$

**Power-law fit of density depletion:**

$$
\log |\overline{\delta\rho}(r)| = n_\rho \log r + C_\rho
$$

**Success criterion**: $n_\rho = -1$ (1/r depletion) with $R^2 > 0.95$.

**Gravitational acceleration from depletion:**

A second particle (test mass) in this density field experiences a force proportional
to $\nabla\rho$. For $\delta\rho \sim -Q_\mathrm{eff}/(4\pi r)$:

$$
|\nabla\delta\rho| = \frac{Q_\mathrm{eff}}{4\pi r^2}
$$

which is the $1/r^2$ force law, with effective charge $Q_\mathrm{eff} = E_2 / (2\rho_0)$
(from v6 results).

---

## 7. Test Matrix

### Phase 0 (1D Toy Model — Frequency Calibration)

Before running expensive 3D simulations, calibrate the triad frequencies using a
1D radial reduction. The spherically-symmetric sector of the S³ sigma model + Skyrme
reduces to an ODE for the profile function $f(r,t)$:

$$
\ddot{f} = f'' + \frac{2}{r}f' - \frac{\sin 2f}{r^2}
  + c_4 \left[ \frac{2f'\sin^2\!f}{r^2}\left(f'' - \frac{\sin 2f}{2r^2}\right)
  + \frac{(f')^2\sin 2f}{r^2} - \frac{2f'\sin^2\!f}{r^3} \right]
$$

(the full radial Skyrme EOM; see v2 `src/radial.c` for the static version).

**Procedure:**
1. Initialize $f(r,0) = A \exp(-r^2/2\sigma^2)$ (non-topological bump, $B=0$)
2. Evolve the 1D radial PDE on a fine grid ($N_r = 2000$, $r_\max = 30$)
3. Record $f(0,t)$ and compute its power spectrum $\hat{f}(\omega)$
4. Identify the dominant resonant frequency $\omega_*$ and its harmonics
5. Use $\omega_* / n$ for $n = 1,2,3$ as the exact frequency ratios for the 3D triad

This ensures the 3D initialization uses frequencies that are **natural resonances** of
the nonlinear system, not arbitrary multiples of $\pi c/R_0$. Runtime: seconds.

**1D mode energy diagnostic**: track $E_\mathrm{kin}(t)$ and $E_\mathrm{pot}(t)$ of
the 1D system. If the bump oscillates quasi-periodically (oscillon), the energy sloshes
between kinetic and potential — this is the 1D precursor of the triad energy exchange.

### Phase 1 (Pure sigma model + Skyrme)

| Test | Init | Modes | A | σ | R₀ | c₄ | N | L | Steps | Purpose |
|------|------|-------|---|---|----|----|---|---|-------|---------|
| 1A   | Bump | 3 (hedgehog) | 0.8 | 1.5 | — | 2 | 200 | 16 | 6000 | B=0 dispersal baseline |
| 1B   | Kick | 3 (triad) | 0.5 | 2.0 | from Phase 0 | 2 | 200 | 16 | 6000 | Triad locking from velocity |
| 1B'  | Kick | 3 (triad, rational ω) | 0.5 | 2.0 | tuned | 2 | 200 | 16 | 6000 | Exact commensurate ω from 1D toy |
| 1C   | Disp | 3 (triad) | 0.4/0.4/0.56 | 2.0 | from Phase 0 | 2 | 200 | 16 | 6000 | Triad locking from displacement |
| 1D   | Disp | 1 | 0.56 | 2.0 | from Phase 0 | 2 | 200 | 16 | 6000 | Single-mode control |
| 1E   | Bump | 3 (hedgehog) | 0.8 | 1.5 | — | 0 | 200 | 16 | 6000 | No Skyrme control |
| 1F   | Kick | 3 (triad) | 0.5 | 2.0 | from Phase 0 | 0 | 200 | 16 | 6000 | Triad without Skyrme |

**High-resolution follow-up** (run if any Phase 1 test shows persistence > 10× control):

| Test | Base | N | L | Steps | Purpose |
|------|------|---|---|-------|---------|
| 1X   | Best of 1B/1B'/1C | 256 | 20 | 10000 | 2+ decades in far-field fit |
| 1Y   | Best of 1B/1B'/1C | 300 | 24 | 12000 | Confirm power-law exponent |

The larger box ($L = 20$–$24$) with absorbing layer pushes the usable far-field range
to $r \in [4, 10]$, giving ~2.5 decades in log-log — sufficient to distinguish $r^{-2}$
from $r^{-1}$ reliably.

### Phase 2 (With conserved density)

| Test | Init q | c_s | ρ₀ | N | L | Steps | Purpose |
|------|--------|-----|-----|---|---|-------|---------|
| 2A | Best from Phase 1 | 1.0 | 1.0 | 200 | 20 | 10000 | Density depletion profile |
| 2B | Best from Phase 1 | 0.5 | 1.0 | 200 | 20 | 10000 | Slower c_s (stronger depletion) |

---

## 8. Success / Failure Criteria

### Phase 1

| Outcome | Criterion | Interpretation |
|---------|-----------|----------------|
| **Strong success** | Triad lifetime > 100× control; phase lock holds; far-field $n_\delta \leq 2$ | Dynamic resonance with slow far-field falloff |
| **Moderate success** | Triad lifetime > 10× control; partial phase lock | Phase locking works but eventually breaks |
| **Marginal** | Triad lifetime 2-10× control | Skyrme term extends lifetime but no locking |
| **Null** | Triad disperses at same rate as control | No triad locking effect; abandon approach |

### Phase 2

| Outcome | Criterion | Interpretation |
|---------|-----------|----------------|
| **Strong success** | $\overline{\delta\rho}(r) \sim -C/r$ with $R^2 > 0.95$; $C = E_2/(4\pi)$ | V6 mechanism confirmed for dynamic resonance |
| **Moderate success** | 1/r depletion forms but $C \neq E_2/(4\pi)$ | Mechanism works but coefficient modified by dynamics |
| **Null** | Depletion falls faster than 1/r | Massless density mode not sourced coherently |

### Known Issues and Risks

1. **B=0 oscillons may not exist** on S³ + Skyrme. The Skyrme term stabilizes topology,
   not energy. Without topological charge, there may be no barrier to dispersal. This is
   the primary risk.

2. **Radiation loss**: even if phase-locked, the oscillon radiates at higher harmonics
   ($2\omega_1$, $2\omega_2$, etc.) not covered by the triad lock. This slowly drains
   energy. The question is whether the drain rate is small enough for meaningful far-field
   measurement.

3. **Grid artifacts**: at N=200, the absorbing layer eats 24 points per side, leaving an
   effective box of ~152 points (L_eff ≈ 12.2). Far-field measurement range is limited
   to $r \in [4, 6]$ — only 1.5 decades. Power-law fits over this range have large
   uncertainties. Consider N=256 or N=300 if Phase 1 shows promise.

4. **Phase 2 sign**: the density equation has source $-\frac{1}{2}|\nabla q|^2 + \frac{1}{2}|\dot{q}|^2$.
   For a dynamic resonance, the kinetic and gradient energies are comparable
   ($\langle|\dot{q}|^2\rangle \sim \langle|\nabla q|^2\rangle$), and partial cancellation
   could weaken the effective source. Need to check whether $\langle S \rangle_t \neq 0$
   or if it averages to zero.

---

## 9. Implementation Notes

### Building on v18

The v18 `skyrme3d.c` provides the complete PDE infrastructure:
- Quaternion field storage and S³ projection
- Skyrme stress computation (two-pass G → div G)
- Velocity-Verlet integrator with absorbing boundaries
- Energy and baryon number diagnostics
- Radial shell-averaging and power-law fitting

**New code needed for v19:**
1. New initialization functions (Experiments 1A-1F)
2. Mode decomposition diagnostic (Section 5.4)
3. Time-averaged shell accumulator (Section 5.6)
4. FFT of central field value (Section 5.7)
5. Phase 2: density field arrays, density EOM stepper, coupled integrator

### Compilation

```
gcc -O3 -fopenmp -o triad3d v19/src/triad3d.c -lm -lfftw3
```

### Runtime Estimate

At N=200 with v18 performance (~500 steps/s on 8 cores): 6000 steps ≈ 12 seconds.
Phase 2 adds ~30% overhead for the density field. Total Phase 1 test suite: ~2 minutes.
