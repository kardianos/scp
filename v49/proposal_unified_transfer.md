# V49 Proposal: Unified Energy Transfer with Integer Field Values

## The Problem with Current Equations

The current system has TWO separate mechanisms bolted together:

    V(P) = (mu/2) P² / (1 + kappa P²)    ← saturates (dies at high P)
    L_gap = -(1/2) lambda_theta P² |θ|²   ← added separately (V48)

When P is large, V(P) goes flat — the binding force turns off and
energy pools in a dead zone. Lambda_theta partially catches this, but
it's an additive band-aid, not a unified transfer mechanism. Energy
isn't REDIRECTED from phi to theta — V(P) just stops working while
lambda_theta independently adds a theta mass.

A correct theory should have:
- ONE unified coupling where phi-binding energy smoothly BECOMES
  theta energy as the local geometry demands
- No dead zones — energy always flows, always transfers
- Exact conservation per interaction, not just globally
- The ratio of phi-energy to theta-energy determined by local geometry

## Design Principles for Integer Fields

Integer-valued fields change the computational model fundamentally:

**Continuous (current):** Fields are doubles. Forces are accelerations.
The integrator nudges values smoothly. Saturation means the force
function goes flat — energy hides in the dead zone.

**Integer (proposed):** Fields are discrete states. Changes are TRANSITIONS
between states. Each transition has exact energy cost/release. There is
no "almost" — either the transition happens (energy transfers completely)
or it doesn't (energy stays). No dead zones because there's no smooth
function to go flat.

Key properties:
- Conservation is EXACT per transition (integer arithmetic, no float drift)
- Transfer is ALL-OR-NOTHING per quantum (but many quanta per timestep)
- Geometry determines WHICH transitions are allowed and their energy cost
- The field can't "hide" energy in an unreachable continuous dead zone

---

## Variant A: Discrete Cosserat Lattice

### The fields

Each voxel stores two integer vectors:

    phi[3]:   int16    displacement state (-32768 to +32767)
    theta[3]: int16    rotation state (-32768 to +32767)

A fixed-point representation where 1 integer unit = A_unit (a conversion
constant, e.g., A_unit = 1/256 gives phi range ±128.0 in code units).

The total energy at a voxel is computed from the integer state:

    E_mass   = (m²/2) × sum(phi[a]²) × A_unit²
    E_grad   = (1/2) × sum(phi[a] - phi_neighbor[a])² × A_unit² / dx²
    E_pot    = V(P_int) where P_int = phi[0]*phi[1]*phi[2] (integer triple product)
    E_theta  = similar for theta
    E_curl   = eta × sum(phi[a] × curl_int(theta, a)) × A_unit²

### The update rule

Instead of force → acceleration → velocity → position (continuous):

1. For each voxel, compute the energy cost of each possible transition:
   - phi[a] → phi[a] ± 1
   - theta[a] → theta[a] ± 1
2. The transition happens if it LOWERS the total energy of the voxel
   + its neighbors (microcanonical evolution)
3. The released energy goes into the kinetic degrees of freedom
   (separate integer counters for velocity)

This is a lattice Boltzmann / cellular automaton hybrid. Each step
is exactly energy-conserving because the transition energies are
computed from integer arithmetic.

### The unified transfer

There is no separate V(P) and lambda_theta. There is ONE energy
function:

    E_voxel = E_mass(phi) + E_grad(phi) + E_grad(theta)
            + E_interact(phi, theta)

    E_interact = -J × |P_int| × (1 - |theta_sum|² / theta_max²)

Where J is a coupling constant and theta_max is the maximum theta
magnitude per voxel. The key: as |theta| grows toward theta_max,
E_interact weakens — theta "absorbs" the interaction energy. As
|theta| → 0, E_interact is strongest — all coupling is in the phi
sector. Energy flows between phi-binding and theta-excitation through
this single term.

### Pros/Cons

+ Exact conservation (integer arithmetic)
+ No dead zones (discrete transitions)
+ Unified phi-theta coupling
+ Natural UV cutoff (integer resolution)
- Completely new simulation architecture (not a PDE integrator)
- Isotropy on cubic lattice is approximate (lattice artifacts)
- Massive rewrite of all tools

---

## Variant B: Fixed-Point PDE with Transfer Potential

### Keep the PDE structure, but redesign the potential

This variant keeps the continuous (float) field values and the Verlet
integrator, but replaces V(P) + lambda_theta with a single unified
potential that transfers rather than saturates.

### The unified potential

    W(P, Θ) = (mu/2) × P² × Θ² / (Θ² + Θ_c²)

Where Θ = |theta|² (local theta energy density) and Θ_c is a scale
constant.

Behavior:
- Low theta (Θ << Θ_c): W ≈ (mu/2) P² × Θ/Θ_c² — weak, theta-dependent
  binding. The binding force on phi is proportional to how much theta is
  present. No theta → no binding.

- High theta (Θ >> Θ_c): W → (mu/2) P² — full binding, theta-independent.
  The binding saturates but BECAUSE theta has absorbed its share, not
  because of an artificial cap.

The force on phi_a:

    F_phi_a = -dW/dphi_a = -mu × P × dP/dphi_a × Θ² / (Θ² + Θ_c²)

The force on theta_a (back-reaction):

    F_theta_a = -dW/dtheta_a = -mu × P² × Θ_c² × theta_a / (Θ² + Θ_c²)²

The second equation IS the energy transfer: when P is large (strong
binding), theta_a feels a force that PUSHES it toward higher amplitude.
The binding energy in P literally drives theta excitation. And as theta
grows, the binding force on phi strengthens (up to the full mu P²).

### The S-curve

The transfer follows a sigmoid:

    binding_fraction = Θ² / (Θ² + Θ_c²)

    Θ << Θ_c: binding ≈ 0      (theta hasn't absorbed yet)
    Θ =  Θ_c: binding = 0.5    (half transferred)
    Θ >> Θ_c: binding → 1.0    (fully transferred)

This is the S-curve the user described — not linear, not step function,
but a smooth sigmoid transfer from "energy in theta" to "energy in
binding." The Θ_c parameter controls the midpoint.

### What happens physically

1. A braid forms (P grows from background fluctuation)
2. The growing P creates a force on theta (F_theta ∝ P²)
3. Theta is excited — energy flows from phi kinetic into theta
4. As theta grows, the binding strengthens (the sigmoid climbs)
5. Equilibrium: theta at the value where the sigmoid is ~0.5,
   binding force balanced against theta radiation

The braid CAN'T exist without theta — the binding requires theta to
be present. And theta CAN'T grow without a braid — the force on theta
requires P ≠ 0. They're mutually dependent through one potential.

### Energy conservation

    E_total = (1/2)|dphi/dt|² + (1/2)|nabla phi|² + (1/2)m²|phi|²
            + (1/2)|dtheta/dt|² + (1/2)|nabla theta|²
            + W(P, Θ)
            + eta × phi · curl(theta)

W(P, Θ) is a single potential energy that spans both sectors. There's no
"V(P) + lambda_theta P² theta²" — just W. The Verlet integrator conserves
E_total to O(dt²) as before.

### Parameters

    mu:    binding strength (same as current, free)
    Θ_c:  transfer midpoint (new, free — replaces both kappa and lambda_theta)

This REDUCES the parameter count: current theory needs mu, kappa,
lambda_theta (3 params). This variant needs mu, Θ_c (2 params). The
saturation and the theta coupling are unified into one mechanism.

### Pros/Cons

+ Minimal code change (same PDE integrator, just replace the potential)
+ Fewer parameters (2 vs 3)
+ Unified transfer (no dead zones)
+ S-curve behavior naturally
+ All existing analysis tools still work
- Still continuous (no exact integer conservation)
- Theta must be present for binding (changes initialization — can't
  start with theta=0)
- The Θ_c parameter needs calibration

---

## Variant C: Quantized Transfer with Integer Quanta

### Hybrid: continuous PDE evolution + discrete energy transfer events

This is the middle ground — keep the PDE for wave propagation but
quantize the TRANSFER of energy between phi and theta.

### The fields

    phi[3]:   double    (continuous, standard PDE)
    theta[3]: double    (continuous, standard PDE)
    Q[3]:     int32     (integer binding quanta per voxel per axis)

Q is a NEW integer field that counts "how many quanta of binding energy
have been transferred from phi to theta at this voxel." It's the ledger
of the transfer.

### The evolution

Standard PDE for phi and theta (Laplacian, mass, curl coupling). But
the potential is:

    V(P, Q) = mu_eff(Q) × P²

Where mu_eff depends on the integer quantum count:

    mu_eff(Q) = mu_base + Q × delta_mu

Each quantum Q increases the effective binding strength by delta_mu.
More quanta → stronger binding. But each quantum also creates a theta
mass contribution:

    m_theta_eff² = m_theta² + Q × delta_m²

### The transfer rule

At each timestep, for each voxel, check:

    E_available = |P| × |dV/dP|        (available binding energy)
    E_quantum = delta_E                 (fixed energy per quantum)

    If E_available > E_quantum AND |theta| > theta_min:
        Q += 1                          (absorb one quantum)
        phi_kinetic -= delta_E          (debit from phi sector)
        theta_kinetic += delta_E        (credit to theta sector)

    If E_available < E_quantum/2 AND Q > 0:
        Q -= 1                          (release one quantum)
        phi_kinetic += delta_E
        theta_kinetic -= delta_E

Conservation is EXACT per transfer event — delta_E is moved between
sectors as an integer multiple. The phi and theta kinetic energies are
adjusted by adding/subtracting velocity increments:

    |delta_v| = sqrt(2 × delta_E / m_eff)

### The geometry dependence

The transfer direction (which velocity component gets the kick) is
determined by the local gradient of P:

    direction = nabla P / |nabla P|     (toward increasing binding)

This means the energy transfer has geometric content — it pushes phi
toward higher P (compression) and theta along the curl (rotation).

### Parameters

    mu_base:   base binding strength (replaces mu)
    delta_mu:  binding strength per quantum
    delta_E:   energy per quantum
    delta_m²:  theta mass per quantum

### Pros/Cons

+ Exact conservation per transfer event (integer accounting)
+ Continuous PDE for wave propagation (smooth, no lattice artifacts)
+ Natural quantum structure (binding comes in discrete packets)
+ Ledger field Q is directly observable (counts binding events)
- More complex than Variant B (3 fields + Q ledger)
- The transfer rule is ad-hoc (not derived from a Lagrangian)
- Requires careful tuning of delta_E to avoid too-frequent or too-rare
  transitions

---

## Recommendation

**Variant B first** (unified transfer potential) — minimal code change,
testable immediately, reduces parameters. If the S-curve transfer gives
qualitatively different behavior (binding without dead zones), it validates
the core idea.

**Then Variant C** (quantized transfer) as the next step — adds the
integer conservation that makes the energy accounting exact per event.
The Q ledger is observable and interpretable.

**Variant A** (full integer lattice) is the long-term vision but requires
a complete rewrite. It should only be attempted after B and C demonstrate
that unified transfer + discrete quanta produce better physics than the
current additive approach.

## Quick Implementation Test for Variant B

Replace the force computation in compute_forces:

```c
// OLD (current):
double P2 = P*P;
double den = 1.0 + keff*P2;
double dVdP = MU * P / (den*den);
double mtheta2_eff = MTHETA2 + LAMBDA_THETA * P2;
...
phi_acc[a] = lap - me2*phi[a] - dVdP*dPda + eta_eff*ct
           - LAMBDA_THETA*P*dPda*theta2;
theta_acc[a] = lapt - mtheta2_eff*theta[a] + eta_eff*cp;

// NEW (Variant B):
double P2 = P*P;
double T2 = theta2;  // |theta|² already computed
double Tc2 = THETA_C * THETA_C;
double sigmoid = T2 / (T2 + Tc2);

// Unified potential: W = (mu/2) P² × sigmoid
double dWdP = MU * P * sigmoid;
double dWdtheta_a = MU * P2 * Tc2 * theta[a] / ((T2 + Tc2)*(T2 + Tc2));
...
phi_acc[a] = lap - me2*phi[a] - dWdP*dPda + eta_eff*ct;
theta_acc[a] = lapt - dWdtheta_a + eta_eff*cp;
```

One new parameter: THETA_C (replaces both kappa and lambda_theta).
