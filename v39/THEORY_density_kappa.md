# Theoretical Analysis: Density-Dependent Saturation (Variant F)

## 1. Setup and Notation

The standard SCP equation of motion derives from the Lagrangian density:

    L = (1/2) sum_a [(d_t phi_a)^2 - |grad phi_a|^2] - (m^2/2) sum_a phi_a^2 - V(P)

with V(P) = (mu/2) P^2 / (1 + kappa P^2), P = phi_0 phi_1 phi_2.

Parameters: m^2 = 2.25, mu = -41.345, kappa = 50.

The saturation ceiling is V_max = |mu/(2 kappa)| = 0.4135 per volume element.
This bounds the binding energy density and prevents unbounded collapse.

**Variant F proposes**: kappa_eff = kappa_0 / (1 + gamma Sigma), where
Sigma = sum_a phi_a^2 is the local field intensity.

---

## 2. Lagrangian Formulation and Energy Conservation

**Question 5 answered first, since everything else depends on it.**

### 2.1 The modified potential

With density-dependent kappa, the potential becomes:

    V(phi) = (mu/2) P^2 / (1 + kappa_eff(Sigma) P^2)

where kappa_eff = kappa_0 / (1 + gamma Sigma) and Sigma = sum_a phi_a^2.

This is a well-defined function V(phi_0, phi_1, phi_2) of the three field
values at each point. Explicitly:

    V = (mu/2) P^2 / (1 + [kappa_0/(1+gamma Sigma)] P^2)
      = (mu/2) P^2 (1 + gamma Sigma) / (1 + gamma Sigma + kappa_0 P^2)

This expression is smooth, bounded from below (mu < 0 makes V <= 0), and
differentiable for all field values.

### 2.2 Valid Lagrangian: YES

The Lagrangian density:

    L = T - U
    T = (1/2) sum_a (d_t phi_a)^2
    U = (1/2) sum_a |grad phi_a|^2 + (m^2/2) sum_a phi_a^2 + V(phi)

is a standard Lagrangian of the form L(phi, d_mu phi). The potential V
depends only on the field values (not derivatives), and Sigma depends only
on the field values at the same point (no nonlocality). Therefore:

- The Euler-Lagrange equations are well-defined.
- Noether's theorem applies: time-translation invariance gives a
  conserved energy E = integral of [T + U] d^3x.
- **Energy IS conserved.** The modification is physically consistent.

### 2.3 Modified equation of motion

The force on field component a is:

    d^2 phi_a / dt^2 = laplacian(phi_a) - m^2 phi_a - dV/d(phi_a)

The derivative dV/d(phi_a) has TWO contributions because phi_a appears
in both P and Sigma:

    dV/d(phi_a) = (dV/dP)(dP/d phi_a) + (dV/dSigma)(dSigma/d phi_a)

where dP/d phi_a = product of the other two fields, and
dSigma/d phi_a = 2 phi_a.

**Term 1** (standard V'(P) coupling, modified denominator):

    dV/dP = mu P / [(1 + gamma Sigma + kappa_0 P^2)^2] * (1 + gamma Sigma)^2

More precisely, defining D = 1 + gamma Sigma + kappa_0 P^2:

    dV/dP = mu P (1 + gamma Sigma) (1 + gamma Sigma + 2 kappa_0 P^2) / (D^2 ... )

Let me compute this carefully. Write V = (mu/2) P^2 (1+gamma Sigma) / D
where D = 1 + gamma Sigma + kappa_0 P^2.

    dV/dP = mu P (1+gamma Sigma) / D - (mu/2) P^2 (1+gamma Sigma) (2 kappa_0 P) / D^2
           = mu P (1+gamma Sigma) / D^2 * [D - kappa_0 P^2]
           = mu P (1+gamma Sigma) / D^2 * (1 + gamma Sigma)
           = mu P (1+gamma Sigma)^2 / D^2

So: dV/dP = mu P (1 + gamma Sigma)^2 / (1 + gamma Sigma + kappa_0 P^2)^2

At gamma=0 this reduces to mu P / (1 + kappa_0 P^2)^2 — correct.

**Term 2** (NEW density-coupling force):

    dV/dSigma = (mu/2) P^2 * d/dSigma [(1+gamma Sigma)/D]
              = (mu/2) P^2 * [gamma D - (1+gamma Sigma) gamma] / D^2
              = (mu/2) P^2 * gamma * [D - 1 - gamma Sigma] / D^2
              = (mu/2) P^2 * gamma * kappa_0 P^2 / D^2
              = (mu gamma kappa_0 / 2) P^4 / D^2

and dSigma/d phi_a = 2 phi_a, so:

    (dV/dSigma)(dSigma/d phi_a) = mu gamma kappa_0 P^4 phi_a / D^2

This is a **self-interaction** proportional to phi_a P^4. It acts as
an additional attractive (mu < 0) force that grows as P^4 and is
proportional to the local field amplitude phi_a.

### 2.4 Summary of the force

    F_a = lap(phi_a) - m^2 phi_a
          - [mu P (1+gamma Sigma)^2 / D^2] * (dP/d phi_a)     ... (I: triple-product)
          - [mu gamma kappa_0 P^4 / D^2] * phi_a               ... (II: density feedback)

where D = 1 + gamma Sigma + kappa_0 P^2.

Term (II) is the NEW contribution from Variant F. Note:
- It is proportional to P^4 (fourth power of the triple product).
- It is proportional to gamma (vanishes at gamma=0).
- It has the SAME sign as the standard V'(P) term (both proportional to mu < 0).
- It acts as an additional **attractive self-interaction**: where P is large,
  the force pulls phi_a further from zero, amplifying P.

**This is the self-reinforcing mechanism.** At high density, kappa_eff drops,
the saturation ceiling rises, and Term (II) provides extra attraction.

---

## 3. Stability Analysis

### 3.1 Small-amplitude stability (vacuum)

In the vacuum, phi_a ~ A_bg cos(...) with A_bg ~ 0.1, so:
- Sigma ~ 3 * 0.01 = 0.03
- P ~ 0.001
- gamma Sigma ~ 0.003 (at gamma=0.1) to 0.3 (at gamma=10)
- kappa_0 P^2 ~ 5e-5

At these values, D ~ 1 + gamma Sigma >> kappa_0 P^2, and Term (II) is
negligible (proportional to P^4 ~ 10^-12). The vacuum equation is
indistinguishable from the standard equation.

**Vacuum is stable for all gamma.** The density-dependent kappa only
activates where P is large (braid cores).

### 3.2 Large-amplitude behavior (braid core)

At the braid core: phi_a ~ 0.8, Sigma ~ 3 * 0.64 = 1.92, P ~ 0.512.
kappa_0 P^2 ~ 13.1.

| gamma | gamma*Sigma | kappa_eff | D | V_max_eff | Term II / Term I |
|------:|------------:|----------:|---:|----------:|-----------------:|
|   0   |       0     |    50     | 14.1 |   0.41  |       0          |
|   0.1 |       0.19  |    42     | 14.3 |   0.49  |    0.010         |
|   1.0 |       1.92  |    17     | 16.0 |   1.22  |    0.063         |
|   5.0 |       9.60  |    4.7    | 23.7 |   4.40  |    0.15          |
|  10   |      19.2   |    2.5    | 33.3 |   8.27  |    0.19          |

V_max_eff = |mu / (2 kappa_eff)| is the effective saturation ceiling.

At gamma=10, the binding ceiling rises from 0.41 to 8.27 — a factor
of **20x**. This is enormous. It means the V(P) well can now capture
20 times more energy per unit volume than the standard equation allows.

### 3.3 Does collapse to a point occur?

Consider a spherical blob of radius R with uniform field amplitude A.
Energy components (3D):

    E_grad ~ A^2 R         (gradient: surface/volume ~ 1/R, times R^3)
    E_mass ~ m^2 A^2 R^3   (mass: volume integral)
    E_pot  ~ V_eff A^2 R^3 (binding: volume integral, V_eff < 0)

Under Derrick scaling x -> x/lambda (shrink by lambda):

    E_grad -> lambda E_grad
    E_mass -> lambda^-3 E_mass     [increases as blob shrinks]
    E_pot  -> lambda^-3 E_pot(lambda)

The critical question: how does E_pot scale?

**Standard kappa (constant)**: V saturates at V_max = 0.41, so E_pot is
bounded. As the blob shrinks (lambda increases), E_pot/volume is capped.
E_pot -> lambda^-3 * V_max_vol, while E_grad -> lambda * E_grad.
At some lambda, E_grad wins. Collapse is arrested.

**Density-dependent kappa**: As the blob shrinks, Sigma = A^2 * (R_0/R)^n
increases (field amplitude grows to conserve energy, or density grows
due to compression). With higher Sigma, kappa_eff drops, V_max rises.

The question is whether V_max rises FAST ENOUGH to overcome the gradient
energy growth.

Under uniform compression (fixed particle number, decreasing R):
- Sigma ~ A_0^2 (R_0/R)^3  (energy conservation for massive field)
- kappa_eff = kappa_0 / (1 + gamma Sigma) ~ kappa_0 / (gamma A_0^2 (R_0/R)^3)
  for strong compression
- V_max ~ |mu| gamma A_0^2 (R_0/R)^3 / (2 kappa_0)

So V_max grows as R^{-3} under compression. The binding energy:

    |E_pot| ~ V_max * R^3 ~ |mu| gamma A_0^2 R_0^3 / (2 kappa_0) = CONSTANT

The binding energy is **independent of R** under compression! Meanwhile:

    E_grad ~ A^2 / R ~ A_0^2 R_0^3 / R^4  (grows as R^{-4})

**Result**: Gradient energy grows as R^{-4}, binding stays constant.
**Gradient ALWAYS wins at small R.** Collapse to a point does NOT occur.

However, this analysis assumes uniform compression. In practice, the
central density grows faster than the periphery, and the gradient energy
is concentrated at the surface. The actual dynamics will produce a
core-halo structure: a dense core where V(P) is near its (raised)
ceiling, surrounded by a gradient-supported shell.

### 3.4 The arrest mechanism

The collapse is arrested when the gradient pressure at the surface of the
dense core balances the enhanced V(P) binding. This occurs at:

    dE_grad/dR = dE_pot/dR

For the 3D scaling above, this gives a stable equilibrium at finite R.
The enhanced V_max permits a MORE COMPACT equilibrium than the standard
equation (smaller R, denser core), but NOT a singularity.

**Conclusion on stability**: Variant F does NOT produce singularities.
The gradient energy provides an absolute floor. But it DOES permit
much denser, more compact structures than the standard equation.

---

## 4. Derrick's Theorem Revisited

### 4.1 Standard Derrick argument

For a static localized solution in 3+1D with energy E = E_grad + E_mass + E_pot,
under the scaling phi(x) -> phi(lambda x):

    E(lambda) = lambda E_grad + lambda^{-3} (E_mass + E_pot)

Setting dE/d lambda = 0 at lambda=1:

    E_grad = -3(E_mass + E_pot)

Since E_grad > 0 and E_mass > 0, this requires E_pot < -E_mass/3 < 0.
Taking d^2E/d lambda^2:

    d^2E/d lambda^2 = 12(E_mass + E_pot) = -4 E_grad < 0

This is a MAXIMUM, not a minimum. Static solutions are unstable to
Derrick scaling — they are saddle points.

### 4.2 Does density-dependent kappa change this?

With V = V(phi; Sigma(phi)), the potential is still a function of phi
alone (Sigma = sum phi_a^2 is not an independent variable). Under
Derrick scaling phi(x) -> phi(lambda x), all local quantities transform
together:

    P -> P(lambda x), Sigma -> Sigma(lambda x)

The potential energy transforms as:

    E_pot(lambda) = integral V(phi(lambda x)) d^3x = lambda^{-3} E_pot

This is IDENTICAL to the constant-kappa case. The density-dependent kappa
changes the VALUE of E_pot (it can be more negative) but not its SCALING.

**Derrick's theorem applies equally to Variant F.** Static localized
solutions remain unstable saddle points under Derrick rescaling.

### 4.3 Evasion mechanisms

Derrick's theorem is evaded the same way as in the standard equation:

1. **Time dependence**: Braids are NOT static. They oscillate. The
   virial theorem for oscillating solutions differs from Derrick's
   static argument. The braid's oscillation frequency provides an
   additional energy component that breaks the scaling argument.

2. **Topological obstruction**: If the braid has nontrivial winding,
   Derrick scaling changes the winding density. The winding cannot be
   continuously deformed away, so the scaled configuration is not
   topologically accessible.

3. **Lattice discretization**: On a grid, continuous Derrick scaling is
   not possible. The lattice spacing provides a minimum length scale.

The density-dependent kappa does NOT introduce new Derrick problems
beyond what already exists. Nor does it solve them. It simply changes
the depth of the potential well available to dynamic solutions.

---

## 5. Critical Density Estimate

### 5.1 Definition

Define rho_crit as the field density Sigma at which the V(P) binding
energy density EQUALS the gradient + mass energy density for the most
compact possible configuration at that density.

### 5.2 Braid core energy balance

At the braid core (radius R_core), the energy densities are:

    e_grad ~ A^2 / R_core^2          (gradient energy density)
    e_mass ~ m^2 A^2 / 2             (mass energy density)
    e_pot  ~ V(P_max; Sigma)          (binding energy density)

where P_max = A^3 / (3 sqrt(3)) (maximum P for given amplitude A, from
the constraint that phi_a^2 <= A^2/3 each when Sigma = A^2).

Actually, for the braid with A ~ 0.8 per component:
P ~ 0.8^3 = 0.512, Sigma = 3 * 0.64 = 1.92.

The binding exceeds dispersive forces when |e_pot| > e_grad + e_mass:

    |V(P; Sigma)| > A^2/R_core^2 + m^2 A^2/2

### 5.3 V as a function of gamma

    V(P; Sigma) = (mu/2) P^2 (1 + gamma Sigma) / (1 + gamma Sigma + kappa_0 P^2)

At the braid core (P=0.512, Sigma=1.92):

    V = (-41.345/2)(0.262)(1 + 1.92 gamma) / (1 + 1.92 gamma + 13.1)
      = -5.42 (1 + 1.92 gamma) / (14.1 + 1.92 gamma)

| gamma |    V    | e_grad (R=3) | e_mass | |V| > e_grad + e_mass? |
|------:|--------:|-------------:|-------:|:----------------------|
|    0  |  -0.384 |     0.071    |  0.72  | NO (0.384 < 0.791)    |
|    0.1|  -0.395 |     0.071    |  0.72  | NO (0.395 < 0.791)    |
|    1  |  -0.492 |     0.071    |  0.72  | NO (0.492 < 0.791)    |
|    5  |  -0.913 |     0.071    |  0.72  | YES (0.913 > 0.791)   |
|   10  |  -1.361 |     0.071    |  0.72  | YES (1.361 > 0.791)   |
|   50  |  -4.36  |     0.071    |  0.72  | YES (4.36 > 0.791)    |

### 5.4 Solving for gamma_crit

At the braid core values, set |V| = e_grad + e_mass:

    5.42 (1 + 1.92 gamma) / (14.1 + 1.92 gamma) = 0.791

Let u = 1.92 gamma:

    5.42 (1 + u) / (14.1 + u) = 0.791
    5.42 + 5.42 u = 0.791 * 14.1 + 0.791 u
    5.42 + 5.42 u = 11.15 + 0.791 u
    4.629 u = 5.73
    u = 1.238
    gamma = 0.645

**gamma_crit ~ 0.6 at the standard braid core density (Sigma=1.92).**

For gamma > 0.6, the binding exceeds dispersive forces at the core.
This is the onset of the self-reinforcing collapse regime, analogous
to the V34 inv_3-5 results but achieved through kappa reduction
rather than mass reduction.

### 5.5 Critical density as a function of gamma

More generally, for a blob of amplitude A (so Sigma = 3A^2, P ~ A^3)
and core radius R:

    |mu|/2 * A^6 (1 + 3 gamma A^2) / (1 + 3 gamma A^2 + kappa_0 A^6) = A^2/R^2 + m^2 A^2/2

Defining Sigma_crit as the Sigma value where this is marginally satisfied
at the optimal R (R -> infinity eliminates gradient, so consider mass only):

    |mu|/2 * A^4 (1 + gamma Sigma) / (1 + gamma Sigma + kappa_0 A^6) = m^2/2

For large gamma Sigma >> 1 and kappa_0 A^6:

    |mu| A^4 gamma Sigma / (2 gamma Sigma) = m^2/2
    |mu| A^4 / 2 = m^2/2

This gives A^4 = m^2/|mu| = 2.25/41.345 = 0.0544, so A = 0.483.
Sigma_crit = 3 * 0.233 = 0.70.

But this is the asymptotic limit where gradient is negligible. Including
gradient energy at a finite core size R:

    Sigma_crit(gamma, R) = 3 A_crit^2

where A_crit satisfies:

    |mu|/2 * A^4 * (1 + gamma * 3A^2) / (1 + 3 gamma A^2 + kappa_0 A^6) = m^2/2 + 1/R^2

For the braid (R~3):

| gamma | A_crit | Sigma_crit | Comment |
|------:|-------:|-----------:|---------|
|   0.1 |  0.88  |    2.32    | Above current braid (1.92) |
|   0.6 |  0.80  |    1.92    | Matches current braid |
|   1.0 |  0.74  |    1.64    | Below current braid |
|   5.0 |  0.58  |    1.01    | Easily exceeded |
|  10   |  0.53  |    0.84    | Low threshold |
|  50   |  0.49  |    0.72    | Near asymptotic limit |

**Physical meaning**: For gamma=5, any region with Sigma > 1.01 (three
fields with RMS amplitude > 0.58) enters the self-reinforcing collapse
regime. The current braid core at Sigma=1.92 is nearly 2x above
this threshold.

---

## 6. Comparison with Variant D (Hybrid Mass Coupling)

### 6.1 Variant D mechanism

    m_eff^2 = m_const^2 + alpha / (1 + beta Sigma)

This reduces the effective mass at high density. The core becomes
"softer" — lower mass means less restoring force, so the field can
displace further and bind more tightly through V(P).

V34 results (inv_3-5): the positive feedback loop produces
kappa-saturated collapse. The blob hits V_max = 0.413 everywhere
in the core and then slowly drains because V_max is fixed.

### 6.2 Variant F mechanism

    kappa_eff = kappa_0 / (1 + gamma Sigma)

This raises the saturation ceiling at high density. The core can
bind MORE DEEPLY — the V(P) well deepens without limit as density
grows.

### 6.3 Structural comparison

| Property | Variant D (inv mass) | Variant F (density kappa) |
|----------|---------------------|--------------------------|
| Vacuum stability | m_const^2 provides floor | Standard m^2 unchanged |
| Core softening | YES (lower m_eff) | NO (m^2 unchanged) |
| Ceiling removal | NO (kappa constant) | YES (kappa drops) |
| Collapse mode | m_eff -> 0, V(P) saturates | V(P) -> unbounded |
| Arrest mechanism | kappa saturation (hard wall) | Gradient energy (soft) |
| New force term | No (same EOM structure) | YES (Term II: P^4 phi_a) |
| Lagrangian | Standard (just m(Sigma)) | Modified (new dV/dSigma) |
| Draining | YES (blob radiates) | Uncertain (deeper well) |

### 6.4 Which is more promising?

**Variant D** has a hard ceiling (kappa is constant). Once the blob
hits V_max everywhere, it has nowhere deeper to go. The excess energy
must radiate. The blob drains. This was observed in V34 inv_3-5.

**Variant F** has a RISING ceiling. As the blob compresses, kappa_eff
drops, V_max rises, and the potential well deepens. This provides a
mechanism to ABSORB the energy released during compression into deeper
binding, rather than radiating it.

The key advantage of Variant F: the energy released by compression is
captured by the deepening well, not radiated. This is analogous to
gravitational binding in GR: as matter falls in, the well deepens,
capturing more matter. In Variant D, the well has a fixed depth
(V_max = 0.413), so excess energy overflows as radiation.

**Variant F is structurally superior for producing stable collapsed
structures**, provided the gradient arrest is gentle enough to form
an equilibrium rather than a hard bounce.

### 6.5 Risk assessment

The danger of Variant F: if the deepening well captures energy TOO
efficiently, the structure may undergo a runaway to the grid scale
(the lattice spacing becomes the arrest mechanism, not the physics).
This would be a numerical artifact, not a physical prediction.

Mitigation: at very small R (approaching the lattice spacing), the
discrete Laplacian provides O(h^2) corrections that break the
continuum scaling analysis. If the equilibrium R is only a few grid
cells, the result is grid-dependent and unphysical.

**Prescription**: Run at multiple resolutions (N=128, 256, 512).
If the collapsed core radius R_core scales with the grid spacing,
the result is a lattice artifact. If R_core converges to a fixed
physical size, the collapse is genuine.

---

## 7. The "Black Hole Threshold"

### 7.1 GR analog

In GR, the Schwarzschild radius r_s = 2GM/c^2 marks the boundary
where the gravitational binding equals the rest mass energy. Inside
r_s, no physical process (including light pressure) can prevent
further collapse.

### 7.2 SCP analog

In the SCP theory with Variant F, the analog condition is:

    |E_pot(R)| > E_grad(R) + E_mass(R) + E_kin(R)

where all energies are evaluated for a ball of radius R containing
the collapsing structure. This is the condition that the well is deep
enough to capture ALL the energy.

From Section 5.3, this is satisfied when gamma > gamma_crit ~ 0.6 at
the current braid core density. But this is a NECESSARY condition,
not sufficient. The dynamics must also allow the energy to be captured
(not radiated away during compression).

### 7.3 Trapping condition

A more precise analog of the event horizon: the condition that
WAVE PERTURBATIONS cannot escape from the collapsed core.

A small perturbation delta phi at radius R propagates outward with
group velocity v_g. The perturbation escapes if v_g > 0 (outward).
In the modified theory, the effective wave speed is:

    v_eff^2 = 1 - (2/c^2) * d^2 V / d phi_a^2

(from the linearized equation around the collapsed background).

The trapping condition is v_eff^2 < 0, which means the effective
"mass squared" for perturbations is so large and negative that
the wave equation becomes elliptic (no propagation).

From the modified potential:

    d^2V/d phi_a^2 ~ mu (1+gamma Sigma)^2 / D^2 * [complicated]

At large Sigma and P (deep in the collapsed core):

    d^2V/d phi_a^2 ~ mu gamma^2 Sigma^2 / (gamma Sigma)^2 = mu

So the effective mass correction saturates at mu = -41.345. The
perturbation equation becomes:

    d^2(delta phi)/dt^2 = laplacian(delta phi) - (m^2 + mu)(delta phi)
                        = laplacian(delta phi) - (2.25 - 41.345)(delta phi)
                        = laplacian(delta phi) + 39.1 (delta phi)

This is a TACHYONIC equation! With m_eff^2 = -39.1, perturbations
grow exponentially rather than propagating. The growth rate is:

    omega = sqrt(39.1 - k^2)  for k < sqrt(39.1) = 6.25

Perturbations with wavelength > 2 pi / 6.25 = 1.0 code units are
trapped (they grow rather than propagate out).

### 7.4 Critical radius

The "event horizon" analog: the radius at which the local effective
mass squared transitions from positive (propagating) to negative
(trapped). This occurs where:

    m^2 + d^2V/d phi_a^2 |_bg = 0

This depends on the background profile phi_bg(r) of the collapsed
structure. For a step-function model (uniform core of density Sigma_0
out to R, vacuum outside):

    Inside (r < R): m_eff^2 = m^2 + mu f(Sigma_0, P_0) < 0  (trapped)
    Outside (r > R): m_eff^2 = m^2 > 0                        (propagating)

The transition at r = R is the horizon analog. But unlike GR, this
transition is SMOOTH (continuous profile) and WAVELENGTH-DEPENDENT
(only long-wavelength modes are trapped). Short-wavelength modes
(k > 6.25) always propagate. There is no absolute horizon.

**The SCP "black hole" is frequency-selective**: it traps long-wavelength
modes while leaking short-wavelength modes. This is analogous to a
Hawking radiation-like process, where the temperature is set by the
tachyonic mass scale: T ~ sqrt(|m_eff^2|) ~ 6.25 (code units).

### 7.5 Estimates of rho_crit for trapping

The trapping onset (m_eff^2 = 0) requires the V'' contribution to
equal m^2 = 2.25. In the deep core limit:

    V'' ~ |mu| (1 + gamma Sigma)^2 / (1 + gamma Sigma + kappa_0 P^2)^2 * [O(1)]

Setting V'' = m^2 = 2.25 and solving is parameter-dependent.
For the simpler "is the core energy self-trapping" condition from
Section 5.3, the critical densities are:

| gamma | Sigma_crit | A_crit (per component) | V_max_eff |
|------:|-----------:|----------------------:|----------:|
|   0.6 |    1.92    |         0.80          |   0.85    |
|   1.0 |    1.64    |         0.74          |   1.22    |
|   5.0 |    1.01    |         0.58          |   4.40    |
|  10   |    0.84    |         0.53          |   8.27    |

---

## 8. Recommended Parameter Ranges

### 8.1 Weak regime (gamma = 0.1-0.5)

- kappa_eff at core: 42-34 (small reduction)
- V_max_eff: 0.49-0.61 (20-50% increase over standard)
- Expected behavior: mild enhancement of braid binding, no collapse
- Useful for: fine-tuning braid stability, small corrections

### 8.2 Marginal regime (gamma = 0.5-2.0)

- kappa_eff at core: 34-17 (significant reduction)
- V_max_eff: 0.61-1.22 (ceiling nearly tripled)
- Expected behavior: enhanced binding, possible slow contraction
- **gamma ~ 0.6 is the onset of self-reinforcing collapse** (Section 5.4)
- Useful for: testing the phase transition, mapping ρ_crit

### 8.3 Strong regime (gamma = 2-10)

- kappa_eff at core: 17-2.5 (near removal of saturation)
- V_max_eff: 1.22-8.27 (ceiling raised 3-20x)
- Expected behavior: rapid compression, possible core formation
- **This is where the "black hole" analog should appear**
- Risk: collapse to grid scale (check resolution convergence)
- Useful for: observing the collapse dynamics, measuring R_core

### 8.4 Extreme regime (gamma > 10)

- kappa_eff at core: < 2.5
- V_max_eff: > 8.27 (essentially unsaturated)
- Expected behavior: rapid collapse arrested only by gradient
- High risk of grid artifacts
- Useful for: establishing the limiting behavior

### 8.5 Recommended first experiments

1. **gamma = 0.6** (marginal): start a standard braid and observe.
   Does it slowly contract? Does it find a new equilibrium?
   Compare with gamma=0 control.

2. **gamma = 2.0** (moderate): start a standard braid. Time the
   collapse. Measure R_core at different N to check convergence.

3. **gamma = 5.0** (strong): observe the collapse dynamics. Does the
   core oscillate (pulsating "black hole") or freeze? What fraction
   of the initial mass/energy is captured vs radiated?

---

## 9. Summary of Key Results

1. **Energy IS conserved.** The density-dependent kappa derives from a
   valid Lagrangian. The modification adds a new force term (Term II)
   proportional to P^4 phi_a, but the total energy is a conserved
   Noether charge.

2. **Singularities do NOT form.** Under compression, the gradient energy
   grows as R^{-4} while the binding energy (even with rising ceiling)
   stays constant or grows as R^{-3}. Gradient always wins at small R.
   The collapse arrests at finite radius.

3. **Derrick's theorem is unaffected.** The density-dependent kappa changes
   the depth of the potential well but not its scaling under Derrick
   transformations. Static solutions remain saddle points. Dynamic
   solutions (oscillating braids) evade Derrick as before.

4. **gamma_crit ~ 0.6** at the standard braid core density (Sigma=1.92).
   Above this, the V(P) binding exceeds all dispersive forces at the
   core, and the self-reinforcing collapse mechanism activates.

5. **Variant F is structurally superior to Variant D** for producing
   stable collapsed structures. The rising ceiling absorbs compression
   energy into deeper binding, whereas Variant D's fixed ceiling causes
   energy overflow and draining.

6. **The "black hole" is frequency-selective.** The tachyonic effective
   mass in the deep core (m_eff^2 ~ -39.1) traps long-wavelength modes
   but leaks short-wavelength perturbations. There is no absolute horizon.

7. **Start at gamma=0.6, sweep to gamma=5.0** to map the phase transition
   from stable braid to collapsed structure.
