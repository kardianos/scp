# Quantum Solitons and the Gravitational Coupling

## The Problem

Six classical paths to gravity have been exhaustively explored:

| Path | Mechanism | Outcome |
|------|-----------|---------|
| 1. Finite-λ | ρ(r) in BLV metric | Dead (P/m=2 algebraic) |
| 2. L₆ sextic | Breaks P/m=2 | Nuclear only (0.55 fm) |
| 3. B⁰p coupling | Sourced Poisson | **1/r works**, g_top free |
| 4. Hopfion metric | Non-hedgehog BLV | Core-scale tensor only |
| 5. Constraint/WZW | e₀²=0 Gauss's law | Mechanism valid, g_top free |
| 6. Self-referential | Constant ε₀², constraint | Pion cloud (eigenvalue) |
| 3+4 combined | Scalar on aniso metric | 1/D³ tensor, not 1/D |

**Classical result**: The B⁰p coupling gives exact 1/r gravity (Path 3), but the
coupling constant g_top = 2.8×10⁻¹⁷ is not determined by the algebra. All classical
approaches to derive it have failed.

**The spec requires** (math/03_dynamics.md): "All particle masses, charges, and
coupling constants must emerge from (ρ₀, λ, e, μ, c) plus the topology of the
solutions." So g_top must come from somewhere.

**Key observation**: g_top = 0 at tree level. But the coupling B⁰p is allowed by
every symmetry of the theory (both are SU(2) singlets, both pseudoscalar under
parity). In quantum field theory, every coupling allowed by symmetry is generated
by quantum corrections, even if absent classically. The question is not WHETHER
g_top is generated, but WHAT VALUE the quantum theory produces.

---

## Why Quantum?

### Classical solitons don't exist

Physical nucleons are not classical Skyrmion profiles. They are quantum objects:
quantized collective coordinates give spin J = I = 1/2, vibrational modes have
zero-point energies, and the soliton fluctuates around the classical solution.

The classical profile f(r) is the saddle point of the functional integral:

    Z = ∫ DΨ exp(iS[Ψ]/ℏ)

The physical soliton is the full quantum state Ψ_quantum, which includes:
1. Translational zero-modes → momentum eigenstates
2. Rotational zero-modes → spin/isospin (J = I = 1/2 for nucleon)
3. Vibrational modes → zero-point energies (K=0 breathing, K=1 isorotation)
4. Degenerate sector fluctuations → p and j oscillations around zero

Channel (4) is where quantum gravity might hide.

### Discrete spectrum = specific knots at specific frequencies

The quantum soliton exists only in specific states. In each topological sector
(B, H), the allowed configurations form a discrete tower:

    |B=1, J=I=1/2, n=0⟩   (nucleon ground state)
    |B=1, J=I=3/2, n=0⟩   (Delta resonance)
    |B=1, J=I=1/2, n=1⟩   (Roper resonance N(1440))
    ...

Each state has definite energy, quantum numbers, and a specific spatial profile
that differs from the classical solution. The DEGENERATE SECTOR participates
in these quantum states — even if classically d = 0, quantum fluctuations
produce ⟨d²⟩ ≠ 0.

### The hierarchy problem is universal

Every theory of fundamental physics has the hierarchy problem: why is gravity
10³⁸× weaker than the strong force? Standard physics encodes this as
M_Planck/M_proton ≈ 10¹⁹. We encode it as g_top ≈ 10⁻¹⁷. The question is
whether the quantum theory PREDICTS this ratio.

---

## Four Avenues to Quantum Gravity

### Avenue A: One-Loop Effective Coupling

**Idea**: Pion loop corrections on the soliton background generate an effective
B⁰p vertex at one loop, even though the tree-level coupling is zero.

**Mechanism**: The constraint |q|² + ε₀²p² = ρ₀² creates a pion-pion-p-p
four-point vertex with coupling ~ ε₀²/ρ₀. A pion loop with a B⁰ insertion
(from the topological current) generates:

    g_top^{1-loop} ~ (ε₀²/ρ₀) × (1/16π²) × ∫ G_π(k) B̃⁰(k) d⁴k

where G_π is the pion propagator on the soliton background and B̃⁰ is the
Fourier-transformed baryon density.

The integral is finite (pion has mass m_π = 0.398 code⁻¹), giving:

    g_top^{1-loop} ~ ε₀² × m_π² / (16π² ρ₀) ~ ε₀² × 10⁻³

For G_Newton: g_top ≈ 2.8×10⁻¹⁷, so ε₀² ≈ 2.8×10⁻¹⁴. This is tiny but
finite — the gravitational weakness comes from the smallness of ε₀².

**Critical question**: Is ε₀² a free parameter, or is it fixed by the quantum
theory? If ε₀² is the sixth free parameter, we've merely shifted the problem.
But if ε₀² is determined by anomaly cancellation (Avenue C) or by the
requirement of UV finiteness, then g_top is PREDICTED.

**Advantage**: Most directly computable. Extends existing normal-mode
infrastructure (we have the pion propagator from Channel A, Phase 14).

**Limitation**: Requires knowing ε₀², which is currently free.


### Avenue B: WZW Anomaly Coefficient (Rotating Soliton)

**Idea**: The Wess-Zumino-Witten term has a topologically quantized coefficient
(N_c = 3 in QCD). For a ROTATING soliton (physical nucleon has J = 1/2), the
WZW term dynamically sources the pseudoscalar p. The coefficient is fixed by
topology, not by a free parameter.

**Mechanism**: Track A showed the WZW coupling B^μ w_μ vanishes on STATICS
(w₀ = ṗ = 0 for a static hedgehog). But the physical nucleon ROTATES:

    Ψ(x,t) = A(t) Ψ₀(x) Ã(t),  A(t) ∈ SU(2)

with angular velocity Ω ~ ℏJ/Λ. The rotating soliton has B^i ≠ 0
(spatial baryon current from rotation), and the WZW term becomes:

    L_WZW ⊃ (N/240π²) × B^i × ∂_i(degenerate sector)

The time-averaged source for p from a J = 1/2 nucleon:

    ⟨J=½|L_WZW|J=½⟩ ~ (N/240π²) × ℏ√(3/4)/Λ × B⁰(r)

This gives:

    g_top^{WZW} = N/(240π²) × ℏ/Λ

In code units (ℏ ≈ 0.039, Λ ≈ 142): g_top^{WZW} ≈ N × 1.2×10⁻⁷

For N = 1: g_top ≈ 10⁻⁷, giving gravity ~10²⁰× too strong.
For N = 3 (QCD-like): worse.

**Problem**: The WZW coefficient alone is too large by ~10¹⁰. Either:
- There are additional suppression factors (form factor, angular averaging)
- The WZW contribution is one of several that partially cancel
- The WZW route is not the dominant source of g_top

**Advantage**: No free parameters — coefficient is topologically quantized.
**Limitation**: Naive estimate gives wrong magnitude. Needs careful calculation
of rotating hedgehog matrix element (angular integrals may suppress it).


### Avenue C: Anomaly Cancellation / UV Consistency

**Idea**: The requirement that the quantum theory is consistent (no anomalies,
UV finite, or renormalizable) constrains the parameter space, potentially
fixing ε₀² and thus g_top.

**Mechanism**: In the Standard Model, anomaly cancellation (triangle diagrams
with gauge currents) fixes the hypercharge assignments and predicts the
number of quark colors. The condition is:

    Σ_f Q_f³ = 0  (gauge anomaly cancellation)

For Cl⁺(3,0,1), the analogous conditions would involve the topological
current B^μ and the degenerate sector:

1. **Mixed anomaly**: B^μ–p–p triangle. If B^μ is coupled to both the bulk
   sector (through topology) and the degenerate sector (through ε₀²), the
   triangle diagram must be finite or cancel between sectors.

2. **Gravitational anomaly**: In 4D, a chiral fermion coupled to gravity
   produces gravitational anomalies that must cancel. If the soliton's
   quantized modes include chiral excitations, anomaly cancellation
   constrains the mode spectrum.

3. **Index theorem**: The Atiyah-Singer index theorem relates the number
   of zero modes of the Dirac operator on a soliton background to the
   topological charge. If the degenerate sector has fermionic-like modes
   (half-integer spin from collective coordinates), the index theorem
   constrains ε₀².

**This is the most speculative avenue** but also the most powerful: if it
works, g_top (and thus G_Newton) is uniquely determined by the algebra.

**Process**: Classify all anomalies of the Cl⁺(3,0,1) quantum field theory.
Check which parameter combinations must satisfy cancellation conditions.
Determine if ε₀² is constrained.


### Avenue D: Spectral Quantization and Self-Consistency

**Idea**: The quantized soliton spectrum must be self-consistent: the
zero-point energies of all modes, including the degenerate sector, must
sum to a finite total. The CONDITION for this finiteness constrains the
parameters.

**Mechanism**: The quantum correction to the soliton mass is:

    M_quantum = M_classical + Σ_n (ℏω_n/2) + (ℏ²J(J+1))/(2Λ) + ...

where the sum is over ALL modes (bulk and degenerate). The sum diverges
unless regularized. The regularization prescription introduces counterterms
that depend on the parameters (ε₀² among them).

The requirement that:
1. The renormalized mass M_quantum is finite
2. The spectrum satisfies the Bohr-Sommerfeld quantization conditions
3. The partition function Z = Σ_states exp(-βE_n) converges

may single out specific values of ε₀². In particular, if the degenerate
sector's contribution to the zero-point energy DIVERGES unless ε₀² takes
a specific value, then ε₀² is fixed by quantum self-consistency.

**Connection to "specific knots at specific frequencies"**: The quantized
soliton exists only at specific energy levels. These energy levels depend
on ALL couplings, including g_top. The requirement that the spectrum
matches observations (nucleon mass, Delta-N splitting, pion mass) over-
constrains the five parameters — and the degenerate sector coupling is
the ADDITIONAL constraint needed to close the system.

Currently we have five parameters (ρ₀, λ, e, μ, c) and many observations:
- M_nucleon = 938 MeV
- M_Delta - M_nucleon = 293 MeV
- m_pion = 140 MeV
- r_proton = 0.84 fm
- g_A = 1.27 (axial coupling)
- ...

Five parameters vs. many observables → overconstrained → self-consistency
check. If the degenerate sector participates in the quantum spectrum, it
adds one parameter (ε₀²) but also new observables (gravitational effects).
The system might close.

---

## Research Process

### Phase 1: Rotating Soliton Matrix Elements (extends existing tools)

**Goal**: Compute the time-averaged degenerate sector source from a
quantized J = 1/2 soliton.

**Method**:
1. Use existing moment of inertia Λ = 142 (from normal_modes.c)
2. Quantize collective coordinates: A(t) ∈ SU(2), J = I = 1/2
3. Compute ⟨J=½| B^i(x) |J=½⟩ (spatial baryon current from rotation)
4. Evaluate the WZW matrix element ⟨J=½| B^i ∂_i p |J=½⟩
5. Extract effective g_top^{WZW}

**Numerical tool**: Extend normal_modes.c or write new solver for the
angular matrix elements. The Wigner D-functions are analytical; the
spatial integrals use the existing radial profiles.

**Expected outcome**: Either a definite prediction of g_top (if WZW
gives the right magnitude) or a definite upper bound (if it's too
large/small). Either way, this constrains the problem.

**Estimated effort**: Medium. Analytical derivation + modest numerical work.


### Phase 2: One-Loop B⁰p Vertex

**Goal**: Compute the one-loop correction to the B⁰p coupling from
pion fluctuations on the soliton background.

**Method**:
1. Expand the Lagrangian to quadratic order around the hedgehog:
   Ψ = Ψ₀ + δΨ, separate into pion (angular) and radial modes
2. Include the constraint coupling: |q|² + ε₀²p² = ρ₀²
3. Identify the pion-pion-p-p vertex from the constraint
4. Compute the one-loop diagram: pion loop with B⁰ insertion at one
   vertex and p external leg at the other
5. Extract g_top^{1-loop}(ε₀²)

**Numerical tool**: The pion propagator on the soliton background is
the Green's function of the K=1 angular mode equation (Phase 14,
Channel A). We already have the eigenvalues and eigenfunctions.

**Expected outcome**: g_top as a function of ε₀². Combined with the
Phase 1 result (WZW contribution), the total g_top is:

    g_top = g_top^{WZW} + g_top^{1-loop}(ε₀²)

Setting g_top = 2.8×10⁻¹⁷ gives ε₀². If ε₀² is small and positive,
the framework is self-consistent.

**Estimated effort**: Substantial. Requires careful treatment of the
soliton background propagator and regularization.


### Phase 3: Anomaly Classification

**Goal**: Determine if anomaly cancellation conditions constrain ε₀².

**Method**:
1. Classify the symmetries of the Cl⁺(3,0,1) Lagrangian:
   - SU(2)_L × SU(2)_R chiral symmetry (broken to SU(2)_V by soliton)
   - U(1)_B baryon number (topological)
   - Z₂ parity (p → -p)
   - Possible discrete symmetries from the Clifford algebra
2. Compute triangle diagrams: B^μ-B^ν-B^ρ, B^μ-p-p, etc.
3. Check for anomalies (divergences that cannot be removed by counterterms)
4. Determine if anomaly cancellation requires specific ε₀²

**This is primarily analytical work**, not numerical. It requires the
full structure of the Cl⁺(3,0,1) algebra including the degenerate sector.

**Expected outcome**: Either ε₀² is unconstrained (anomalies cancel for
all ε₀²) or ε₀² is fixed (specific value required for consistency).

**Estimated effort**: Large. Requires deep algebraic analysis.


### Phase 4: Full Quantum Spectrum

**Goal**: Compute the full quantum spectrum of the B=1 soliton including
degenerate sector, and check for self-consistency constraints.

**Method**:
1. Enumerate all modes: K=0 (breathing), K=1 (pion/rho), K=2, ...,
   plus degenerate sector modes (p and j fluctuations)
2. Compute zero-point energies: E_zp = Σ ℏω_n/2
3. Regularize using zeta-function or heat-kernel methods
4. Check: does the regularized sum depend on ε₀²? If so, does the
   requirement E_zp < ∞ (or renormalizability) fix ε₀²?
5. Compare predicted M_quantum with M_nucleon = 938 MeV

**Connection to existing work**: Phase 12 (normal modes) computed the
K=0 spectrum. Phase 14 computed K=1 channels A, B, C. The missing
piece is the degenerate sector's contribution.

**Expected outcome**: Either the quantum spectrum is consistent for
all ε₀² (no constraint) or specific ε₀² values are distinguished
(e.g., by finiteness, by matching M_nucleon, by stability).

**Estimated effort**: Very large. Full spectral computation.

---

## What Would Success Look Like?

### Strong success
G_Newton predicted from (ρ₀, λ, e, μ, c) with no free parameters.
This would mean: the gravitational constant is a DERIVED quantity,
determined by the same parameters that fix the nucleon mass, pion
mass, and nuclear forces. A genuine unification.

### Moderate success
g_top constrained to a discrete set of values (topological quantization).
For example: g_top = n × g₀ where n ∈ Z and g₀ is determined by the
algebra. The hierarchy g₀ ~ 10⁻¹⁷ would still need explanation, but
the discreteness would be a prediction.

### Partial success
ε₀² fixed by anomaly cancellation or UV consistency, giving g_top
as a function of the five parameters. The value might not match G_Newton
exactly (scheme dependence, higher-loop corrections), but the ORDER
OF MAGNITUDE would be a prediction.

### Null result (still informative)
g_top remains free at the quantum level. This would mean either:
- Gravity requires new physics beyond Cl⁺(3,0,1) (e.g., additional fields)
- The coupling is non-perturbative (instanton, not loop)
- The framework is incomplete and needs a sixth parameter

A null result at each phase still constrains the theory and narrows
the search space.

---

## Connection to Gravitational Wave Polarization (Open Problem A3)

The spec identifies gravitational wave polarization as a critical open
problem: GR predicts spin-2 (two polarizations, + and ×), while the
scalar B⁰p mechanism (Path 3) gives spin-0 (one polarization).

**Quantum effects could resolve this**: If the quantum soliton has
spin-dependent gravitational effects (from the WZW/rotation coupling),
the gravitational radiation from a binary system would carry angular
momentum information, potentially producing additional polarization
modes.

Specifically:
- Classical Path 3: scalar p → monopole radiation → one polarization
- Quantum WZW: rotating soliton → B^i source → dipole/quadrupole
  radiation with spin structure → additional polarizations

This is highly speculative but would be testable: the LIGO/Virgo
polarization measurements constrain the graviton spin.

---

## Recommended Starting Point

**Phase 1 (rotating soliton matrix elements)** is the recommended first
step because:
1. It extends existing infrastructure (moment of inertia, profiles)
2. It gives a definite numerical result (g_top^{WZW})
3. The WZW coefficient is topologically quantized → no free parameters
4. If the magnitude is wrong, it immediately constrains which avenues
   remain viable
5. It addresses the gravitational wave polarization question (spin
   structure from rotation)

The key computation is the matrix element:

    ⟨B=1, J=½, m| ∫ B^i(x) ∂_i p(x) d³x |B=1, J=½, m'⟩

in terms of the radial profile f(r) and the collective coordinate
wavefunctions D^{1/2}_{m,m'}(A).
