# v102 kernel step, K1 — pre-registration (the reference simulator)

**Written:** 2026-10-08, before any K1 result. Plan: `PLAN.md`.
Rulings: R8, R9 (working t = 0.1 J0), R10, D11.

Models: the leading-order effective lock (`../perm/cas/lo.py`, reading
V) at g/Δ = 0.1, energies in J0 = 2g²/Δ. The exact one-photon S-matrix
machinery is `../options/cas/e5_cool.py` and `../perm/cas/x4_form.py`.

---

## K1-1. Formation as a quantum jump

**Exact reference.**

- On the n = 4 lock (3+2 and W, leading order): the probability per
  resonant photon (ε_res ± 0.05, polarization averaged) that a hole is
  made. `../perm/` X4: 12.9% (3+2), 5.6% (W).
- The depleted-lock state left behind, ρ_exact. It is computed as the
  final P2 amplitude's flux into the space lead, resolved over the
  depleted lock's eigenstates.

**Jump model.**

- The event applies J_ν = Σ_μ P[μ, ν] b_0μ to the lock's ground state,
  normalized, on the cell next to the photon.
- Its rate is set by one number per photon: the exact probability
  above, the only fitted quantity.

**Gates.**

- (a) The jump state's overlap with ρ_exact (fidelity) ≥ 0.9.
- (b) The jump carries the photon's charges exactly: S^z, and I^z for
  W.
- (c) The energy balance closes: E(photon) + E(lock before) =
  E(lock after) + E(space quantum), to 1e-12. The space quantum takes
  whatever the lock and photon leave.

**Reported:** the jump state's weight in the dressed one-hole band (the
nucleon), and the remainder in nucleon + lock quanta.

## K1-2. Two holes in an isolated lock

- **Setup.**
  - An exact W ring, n = 6, T-floor at t = 0.1 J0. The leading-order
    mobile model has the floor exact.
  - Two holes are created by two jumps, at cells 0 and 3 (antipodal),
    with opposite photon charges, so the pair can be (0, 0).
  - Evolve to T = 200/J0. Also 3+2: ring n = 10, holes at 0 and 5.
- **Measured.**
  - The probability that the holes are in their bound pair state (the
    pair ground state of the (0, 0) sector) at t = 0, and its time
    average.
  - The energy above the pair ground state.
- **Prediction (registered):** the bound-pair weight does not grow. A
  closed lock cannot shed the binding energy: the gap is above the
  binding.
- **Gate:** the time-averaged bound weight ≤ its t = 0 value + 0.05
  (closed dynamics conserves the energy distribution). A violation means
  a coding error.

## K1-3. Capture by cold light (the sink)

- **Setup.**
  - The two-hole sector of the n = 4 lock (3+2 and W, exact), with one
    boundary link to a photon lead.
  - The `../options/` E5 multichannel S-matrix, now with the two-hole
    lock's eigenstates as channels.
  - Cold photons (the lowest δ = 0.01 of the band).
  - A master equation from the state in which the two holes are least
    bound (the highest-energy state of the two-hole manifold below the
    space-carrying levels) into the bound pair.
- **Reported:**
  - the photons needed to reach bound-pair population ≥ 0.9;
  - the per-photon capture probability;
  - whether any photon un-binds the pair (heats it).
- **Gate:** none (a report). It decides whether the kernel needs a
  boundary scattering rule for binding.

## K1-4. The classical field as quanta (a stated rule, not a gate)

The rule K2 will use, written here so K2 cannot choose it after seeing
results:

- A kernel cell's field holds `n = floor(Ee/ε)` whole atoms, with
  ε = A0·w1e/2π. An S-boundary event can fire only if n ≥ 1, and it
  consumes one atom.
- The photon's spin state is read from the two-plane amplitude as a
  Jones vector (fa1, fa2). Its S^z content is the circular projection,
  computed from (fa1, fa2) and the kernel's phase `th2`.
- Under R8 the isospin is the K2b partner pair (fa3, fa4). With that
  pair empty (the default), every photon has I^z = +½.
- The event's rate is the K1-1 probability per resonant atom.
- The event fires only if the neighbour cell's pitch puts the atom
  within the resonance window. The detuning goes into the space
  quantum's energy, so the ledger closes exactly.

## Endings

| ending | condition |
|---|---|
| **Go** | K1-1 passes and K1-2 behaves as registered. K2a may start with the jump rule, and the capture rule is set by K1-3 |
| **Revise** | K1-1 fidelity < 0.9: formation is not one jump. The event rule must carry more structure before K2 |
| **Bug** | K1-2 violates its gate |
