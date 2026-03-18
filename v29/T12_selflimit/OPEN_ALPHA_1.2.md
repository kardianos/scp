# Open Issue: Depletion Exponent α=1.2 (Need 1.0 for Newton)

## Status: PINNED — confirmed but incomplete

The M7 two-component model produces long-range power-law depletion
of the background field with exponent α≈1.21 (confirmed at L=30,60,100).
This is 20% steeper than the Newtonian 1/r required for gravity.

## Confirmed Facts
- Depletion is POWER-LAW (not Yukawa/exponential despite m=1.5)
- α≈1.2 is stable across domain sizes (not a finite-size effect)
- Depletion is broadband (all frequencies, good for universal gravity)
- M7 preserves braid structure (fc=0.91)

## Candidate Fixes (ordered by promise)

### Fix 1: Self-Sourcing Depletion ("Gravity Gravitates")
In GR, the gravitational field sources itself (nonlinear Einstein eqs).
Our depletion is linear — it doesn't self-deepen.

Modification: make coupling g depend on local depletion:
    g_eff(x) = g₀ × ρ₀ / ρ_B(x)

Where ρ_B is depleted (low), g_eff is STRONGER → well deepens itself.
This positive feedback pushes α toward 1.0.

Effort: Medium (modify M7 coupling). Likelihood: High.

### Fix 2: Asymmetric Coupling
M7's g×S²×B and g×B²×S cancel partially. The net depletion is
the difference of two ~1/r profiles → steeper.

Modification: different coupling strengths:
    S feels: -g_absorb × B² × S    (strong, braid absorbs)
    B feels: -g_radiate × S² × B   (weak, braid radiates poorly)

With g_absorb > g_radiate, net depletion is larger → α closer to 1.0.

Effort: Low (change one constant). Likelihood: Medium.

### Fix 3: Condensate Background
Current background is oscillating plane waves. A true condensate
(uniform |φ|=φ₀, no oscillation) has different depletion physics.
In superfluids, a vortex creates 1/r velocity field.

Modification: initialize B fields as uniform condensate, not plane waves.
    B_a(x,0) = φ₀ (constant everywhere)
    B_vel = 0

Effort: Low (change initial conditions). Likelihood: Medium.

### Fix 4: Measure Tensor Strain, Not Scalar Density
T5 found 63% quadrupolar strain. Maybe the TENSOR perturbation h_ij(r)
falls off as 1/r while the SCALAR ρ(r) falls off as 1/r^1.2. Other
braids respond to the tensor metric, not the energy density.

Modification: measure h_ij(r) profile of the B-field strain tensor
at various r, fit power law. If α_tensor < α_scalar: tensor is the
gravitational field.

Effort: Low (analysis only). Likelihood: Medium.

### Fix 5: Counter-Rotating Shell
The braid might need a structural shell (counter-rotating helical wrap)
that creates an impedance boundary between organized core and background.
The shell mediates absorption/radiation with different efficiency →
sharpens the depletion profile.

Effort: High (new initialization, unknown stability). Likelihood: Unknown.

### Fix 6: Binding Energy as Second Topology
A second field structure wrapping the braid (like electron cloud around
nucleus) could limit growth and set the mass. The binding energy between
inner braid and outer wrap IS the gravitational mass.

Effort: High. Likelihood: Speculative.

### Fix 7: Vacuum Lattice Structure
If the background has discrete structure (T9 substrate), the lattice
Green's function approaches 1/r at large r but may differ at intermediate
r. The α=1.2 could be a lattice correction.

Effort: High (requires T9 results). Likelihood: Unknown.

### Fix 8: Nonlinear Potential Mapping
Maybe the gravitational potential isn't ρ(r) but some nonlinear function
like Φ = -ln(1 - δρ/ρ₀). If braids respond to Φ rather than δρ, the
effective α for geodesics could differ from 1.2.

Effort: Low (analysis). Likelihood: Low.

## When to Revisit
- After T10G complex field extension (EMF path) is explored
- After T9 substrate gives constraints on background structure
- If two-braid attraction test shows force law F ~ 1/r^n with n≈2,
  that would confirm gravity works despite α=1.2 for the potential
