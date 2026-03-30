# V49: Integer Field Proposal

## Core Idea

Replace continuous field values + force integration with:
- Integer field values (fixed-point)
- Energy-cost evaluation of proposed moves
- Accept/reject based on conservation law

No forces. No signs to get wrong. No dead zones. No integration error.

## The State

Each voxel stores:

```c
struct Voxel {
    int16_t phi[3];      // displacement field (fixed-point, unit = A_UNIT)
    int16_t theta[3];    // rotation field (fixed-point)
    int16_t phi_v[3];    // velocity (momentum / mass, fixed-point)
    int16_t theta_v[3];  // theta velocity
};
```

Total: 24 bytes per voxel (vs 96 bytes for 12 doubles in current sim).
At N=384: 1.3 GB (vs 5.2 GB). 4x memory reduction.

A_UNIT = 1.0/256 means phi=256 represents 1.0 in code units.
Range ±127.99 in code units (int16 range ±32767 / 256).

## Energy Computation (Integer Arithmetic)

All energies are computed in integer units. Define E_UNIT as the
energy quantum (smallest representable energy change).

```c
// Energy at voxel (i,j,k) — all integer arithmetic
int64_t voxel_energy(int i, int j, int k) {
    int64_t E = 0;

    // Mass: (m²/2) × sum(phi[a]²)
    for (int a = 0; a < 3; a++)
        E += M2_FIXED * (int64_t)phi[a] * phi[a];  // M2_FIXED = m² × scale

    // Gradient: (1/2) × sum((phi[a] - phi_neighbor[a])²) / dx²
    for (int a = 0; a < 3; a++) {
        for (int d = 0; d < 3; d++) {  // 6 neighbors
            int diff = phi[a] - neighbor_phi(d, a);
            E += GRAD_FIXED * (int64_t)diff * diff;
        }
    }

    // Kinetic: (1/2) × sum(phi_v[a]²)
    for (int a = 0; a < 3; a++)
        E += KIN_FIXED * (int64_t)phi_v[a] * phi_v[a];

    // Triple product binding: V(P) where P = phi[0]*phi[1]*phi[2]
    int64_t P = (int64_t)phi[0] * phi[1] * phi[2];  // 48-bit, fits int64
    int64_t P2 = (P >> SCALE_SHIFT) * (P >> SCALE_SHIFT);
    E += MU_FIXED * P2 / (DENOM_FIXED + KAPPA_FIXED * P2);

    // Same for theta (mass, gradient, kinetic)
    // ... (identical structure)

    // Curl coupling: eta × phi · curl(theta)
    // curl is computed from neighbor differences (integer)
    for (int a = 0; a < 3; a++) {
        int curl_t = curl_int(theta, a, neighbors);
        E += ETA_FIXED * (int64_t)phi[a] * curl_t;
    }

    // Transfer: |V(P)| × f(Theta) — unified binding-theta coupling
    int64_t V_abs = abs(MU_FIXED * P2 / (DENOM_FIXED + KAPPA_FIXED * P2));
    int64_t Theta = 0;
    for (int a = 0; a < 3; a++)
        Theta += (int64_t)theta[a] * theta[a];
    E += V_abs * Theta / (Theta + THETA_C_FIXED);  // linear proportion

    return E;
}
```

Note: int64 arithmetic throughout. No floating point anywhere.
The scale factors (M2_FIXED, GRAD_FIXED, etc.) are pre-computed
integer constants that absorb the dx², dt², A_UNIT conversions.

## The Update Rule

Instead of: compute force → integrate velocity → integrate position
Do: propose move → compute energy change → accept if valid

```c
void update_voxel(int i, int j, int k) {
    // Phase 1: PROPAGATION (wave dynamics via velocity)
    //
    // The velocity field carries momentum. Each timestep, velocity
    // shifts the field values. This is the integer analog of the
    // drift step in Verlet.

    for (int a = 0; a < 3; a++) {
        // Proposed new phi value from current velocity
        int16_t phi_new = phi[a] + (phi_v[a] >> VEL_SHIFT);

        // Compute energy of current state vs proposed state
        int64_t E_old = local_energy(phi[a], ...);
        int64_t E_new = local_energy(phi_new, ...);
        int64_t dE = E_new - E_old;

        // The velocity "pays" for the energy change
        // dE > 0: move costs energy → deduct from kinetic
        // dE < 0: move releases energy → add to kinetic
        int64_t KE_available = KIN_FIXED * (int64_t)phi_v[a] * phi_v[a];

        if (dE <= KE_available) {
            // Move accepted: apply position change
            phi[a] = phi_new;
            // Adjust velocity to conserve energy exactly:
            // new_KE = old_KE - dE
            // |new_v| = sqrt(2 * new_KE / m)
            // In integer: use lookup table or Newton's method
            phi_v[a] = isqrt_fixed(KE_available - dE);
            // Preserve velocity sign
            if (phi_v[a] > 0 && phi_new < phi[a]) phi_v[a] = -phi_v[a];
        }
        // else: move rejected, field stays, velocity reverses (bounce)
    }

    // Phase 2: ACCELERATION (forces from neighbors)
    //
    // Compute the energy gradient with respect to velocity.
    // This is the integer analog of the force computation.
    // Instead of F = -dE/dx, we compute:
    //   "If I increase v by 1 unit, what happens to neighbor energies?"

    for (int a = 0; a < 3; a++) {
        int64_t E_here = local_energy_with_neighbors(...);
        // Try v + 1
        phi_v[a] += 1;
        int64_t E_plus = local_energy_with_neighbors(...);
        phi_v[a] -= 2;
        int64_t E_minus = local_energy_with_neighbors(...);
        phi_v[a] += 1;  // restore

        // Discrete gradient: which direction lowers total energy?
        if (E_plus < E_minus)
            phi_v[a] += 1;  // accelerate in + direction
        else if (E_minus < E_plus)
            phi_v[a] -= 1;  // accelerate in - direction
        // else: no change (at energy minimum for this component)
    }

    // Phase 3: TRANSFER (phi ↔ theta energy exchange)
    //
    // This is where the unified transfer happens. No sigmoid needed.
    // Just: can a quantum of energy move from phi to theta?

    int64_t E_phi_sector = phi_kinetic(...) + phi_potential(...);
    int64_t E_theta_sector = theta_kinetic(...) + theta_potential(...);
    int64_t E_binding = binding_energy(P, Theta);

    // Compute transfer energy quantum
    int64_t quantum = TRANSFER_QUANTUM;  // fixed size, e.g. 1 energy unit

    // Try transferring one quantum: phi → theta
    // Temporarily: reduce phi_v magnitude, increase theta_v magnitude
    int64_t E_after_transfer = compute_total_if_transferred(...);

    if (E_after_transfer == E_before_transfer) {
        // Energy conserved exactly → apply transfer
        apply_transfer(phi_v, theta_v, quantum);
    }
    // If not exactly conserved, don't transfer. No approximation.
}
```

## Why This Works Where Continuous Fails

### No sign errors
Every computation is: compute E_old, compute E_new, compare.
If E_new < E_old, the move is energetically favorable. Period.
No sign of mu, no sign of V_base, no chain rule, no partial
derivatives. Just total energy comparison.

### No dead zones
At high P, the integer V(P) is a specific number. The move
"decrease P by 1 unit" has a specific energy cost. Either that
cost is payable from kinetic energy, or it isn't. There's no
"force going flat" — the energy landscape is discrete but the
steps are always well-defined.

### Exact conservation
E_total = sum over all voxels of voxel_energy(). This is an
integer. It changes by exactly dE when a move is accepted, and
the kinetic energy is adjusted by exactly -dE. The books balance
to the last bit. No float drift, no Verlet error, no energy
creep.

### No parametric resonance
The strobing problem (P² oscillating at 6ω) doesn't exist in
the integer model because there's no "effective mass" that
oscillates. The energy cost of a theta move is computed from the
CURRENT integer state, not from a continuous function of an
oscillating variable. If P is currently large (high |phi[0]*phi[1]*phi[2]|),
the transfer cost is one value. If P is currently small, it's
another. No averaging, no time-scale separation needed.

## The Transfer Mechanism (Integer)

The key to the unified transfer. No separate V(P) and lambda_theta.
One computation:

```c
// "Should energy move from phi to theta at this voxel?"
//
// Compute: if I decrease |phi_v[a]| by 1 and increase |theta_v[b]| by 1,
// does total energy stay the same?
//
// The binding energy V(P) depends on phi values.
// The curl coupling depends on both phi and theta.
// If reducing phi velocity lets P relax toward zero (less binding),
// AND that energy appears as theta kinetic energy that drives curl
// coupling, the transfer is self-consistent.

for (int a = 0; a < 3; a++) {
    for (int b = 0; b < 3; b++) {
        // Proposed: phi_v[a] shrinks, theta_v[b] grows
        int16_t new_phi_v = phi_v[a] - sign(phi_v[a]);
        int16_t new_theta_v = theta_v[b] + sign_from_curl(a, b);

        int64_t E_old = total_energy();
        int64_t E_new = total_energy_with(new_phi_v, new_theta_v);

        if (E_new == E_old) {
            phi_v[a] = new_phi_v;
            theta_v[b] = new_theta_v;
            break;  // one transfer per component per step
        }
    }
}
```

The "sign_from_curl" ensures the transferred theta velocity aligns
with the curl geometry (the rotation direction of the braid). This
preserves angular momentum as well as energy.

## Performance Considerations

### Integer arithmetic is FAST
- int16 × int16 = int32: 1 cycle on any CPU/GPU
- int32 × int32 = int64: 1 cycle
- No IEEE 754 rounding, no denormals, no NaN checks
- SIMD: 16 × int16 ops per 256-bit register (vs 4 × fp64)

### Memory is 4x smaller
- 24 bytes/voxel (int16) vs 96 bytes (fp64)
- N=384: 1.3 GB vs 5.2 GB
- Better cache utilization → higher throughput

### But: more evaluations per step
- Each "proposed move" requires 2 energy evaluations (before/after)
- 6 field components × 2 evals = 12 energy calls per voxel per step
- vs 1 force evaluation in the current PDE
- Net: ~6x more compute per step, but each eval is faster (int vs fp64)
- Rough estimate: similar wall-clock time, much better conservation

### GPU mapping
- Each voxel is one thread (same as current)
- Energy evaluation is local (7-point stencil, same as Laplacian)
- No global communication needed
- int16 arithmetic maps well to GPU tensor cores (via dp4a/dp2a)

## What This Gives Up

1. **No continuum limit.** The fields are discrete. The "dx → 0" limit
   doesn't recover a PDE. This is a lattice field theory, not a
   discretized continuum theory.

2. **Lattice artifacts.** Cubic lattice breaks rotational symmetry.
   Braids may align with grid axes. Mitigation: use a large enough
   grid that the lattice spacing is much smaller than the braid core.

3. **No analytic solutions.** Can't linearize, can't do perturbation
   theory, can't compute dispersion relations. Everything is numerical.

4. **Different physics at the margin.** The accept/reject dynamics
   are NOT equivalent to Hamiltonian flow. They're closer to
   microcanonical Monte Carlo. Long-time behavior may differ from
   the PDE at the edges of stability.

## What This Gains

1. **Exact conservation** — not to O(dt²), but to the last bit.
2. **No blowups** — moves that would increase energy beyond available
   kinetic are simply rejected. No NaN, no infinity, ever.
3. **No sign errors** — energy comparison is always well-defined.
4. **No dead zones** — discrete energy landscape has no flat regions
   (except exact degeneracies, which are measure-zero).
5. **Natural quantum structure** — the field comes in discrete units.
   The "quantum" of transfer is built into the representation.
6. **4x memory reduction** — enables larger grids or finer resolution.

## Implementation Path

Phase 1: Write a minimal int16 simulator for a SINGLE scalar field
(no theta, no curl). Verify that wave propagation works, that energy
is conserved exactly, and that localized structures (oscillons) form
and survive. This tests the core accept/reject dynamics.

Phase 2: Add the second field (theta) with curl coupling. Verify
that the curl interaction transfers energy between sectors correctly.
Test single-braid stability.

Phase 3: Add the transfer mechanism. Run two-proton binding test.
Compare to V48 fp64 results.

Phase 4: Optimize for GPU. Profile int16 throughput vs fp64 throughput.

## Open Questions

1. **Is the accept/reject dynamics ergodic?** In microcanonical MC,
   you need to visit all accessible states. The single-flip proposal
   (change one component by ±1) may be too local. May need multi-site
   updates or "cluster" moves.

2. **What's the effective "dt"?** Each step proposes ±1 changes.
   The time resolution is set by the velocity scale. If phi_v is
   typically ~100 (in int16 units), each step moves phi by ~1 unit.
   The effective dt ~ A_UNIT / (v_typical × A_UNIT) = 1/v_typical.

3. **How to set the integer scale factors?** The ratio of M2_FIXED
   to GRAD_FIXED to KIN_FIXED determines the physics. These need
   to be calibrated against the fp64 simulation to reproduce the
   same braid structure.

4. **Can the accept/reject model propagate waves at the correct speed?**
   The PDE has a well-defined dispersion relation (omega² = k² + m²).
   The integer model needs to reproduce this, at least approximately,
   for the physics to match.
