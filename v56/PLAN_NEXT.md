**Proposal: Density-Gradient + Frequency Self-Confinement Mechanisms for v56 Multivector Field**

**Goal**: Introduce interactions between local field density (or |M|²_bulk, rotor norm |q|²) and frequency content that can self-trap dispersing hedgehog or braid excitations without external boundaries. All mechanisms are first-principles additions to the existing 8-component Klein–Gordon + Higgs + minimal Skyrme Lagrangian on the Voronoi foam. They preserve the Hamiltonian structure and Lorentz invariance where possible. Implementation targets the `field_pga.h` plug-in and `foam_sim.c` diagnostics.

### 3. Resonance Locking
**Mechanism**: Local density ρ (or bulk norm deviation) sets a position-dependent natural frequency ω₀(ρ). When the carrier frequency content of the multivector field matches ω₀ at a density peak, phase-locking prevents radiation; energy is trapped in a resonant bound state analogous to a driven oscillator in a varying medium.

**Equations** (added to field_forces):
```
ω₀[i] = ω_base + α · (bulk_norm[i] - v²)          // α tunable
detuning = |∂t M[k] / M[k] - ω₀[i]|                // per-component instantaneous freq
force_res[i] = -β · detuning · M[k][i] * sign(detuning)   // restoring when locked
```
(β sets lock strength; only applied when |detuning| < ω_lock_width.)

**Implementation steps**:
1. Add `omega_base`, `alpha`, `beta`, `lock_width` to Config.
2. In `field_forces`, after Skyrme drive, compute instantaneous frequency via finite-difference velocity or Hilbert transform on recent frames (or use M_vel/M amplitude ratio).
3. Accumulate resonance force only on rotor components (k<4) or full M.
4. Update energy diagnostic with E_res = ½ β detuning² term.
5. Add config knobs to run_skyrme_L40.cfg and new resonance.cfg.

**Free parameters**: α (freq shift per density unit), β (lock gain), ω_base, lock_width.  
**Tuning plan**: Analytic derivation from linearized dispersion relation around vacuum (v56/EM_WAVE style Lean proof); then grid search {α,β} on small L=20 mesh with skyrme seed; CMA-ES on production L=40 if promising. 3–5 GPU runs to map lock window.

**Likelihood**: 45%. Strong analogy to polariton trapping and driven NLSE solitons; easy to add but requires precise freq estimation on unstructured mesh. Risk of over-damping or artificial locking.

### 4. Density-Dependent Effective Mass
**Mechanism**: Local density modulates the inertial term in the kinetic Lagrangian, making high-frequency (high-k) components heavier in low-density regions. High-k waves therefore cannot propagate out of high-ρ cores, self-confining the excitation (ponderomotive-like + variable-mass Klein–Gordon).

**Equations**:
```
m_eff[i] = m0 + γ · (bulk_norm[i] - v²)          // γ >0 → heavier in high density
L_kin = ½ m_eff · (∂t M)² - ½ |∇M|² + …
EOM contribution:  ∂t (m_eff ∂t M) term → extra force = - (∂t m_eff) (∂t M) in Verlet update
```
(Discretized carefully on foam faces to preserve symplecticity.)

**Implementation steps**:
1. Add `m0`, `gamma` to Config (m0 usually 0 for pure Higgs).
2. Pre-compute m_eff per cell in compute_forces or inside Verlet.
3. Modify velocity update: v ← v + (dt/2) a / m_eff; position update unchanged; second half-step uses new m_eff.
4. Add face-based correction for ∇m_eff terms to keep energy conserved.
5. Extend Diag with E_inertia = ½ (m_eff - m0) v².

**Free parameters**: γ (mass gradient strength), m0 baseline.  
**Tuning plan**: Derive γ from requirement that group velocity vg(ρ) drops below escape speed at core edge (analytic from polariton dispersion v56 already has). Grid γ on 2–3 short runs; full CMA-ES once combined with mechanism 3.

**Likelihood**: 55%. Very natural in variable-coefficient wave equations; proven in inhomogeneous media and optical trapping. Mesh discretization of variable-mass Verlet is the main technical risk but solvable.

### 5. Nonlinear Frequency Mixing (Four-Wave)
**Mechanism**: Density gradients provide the wave-vector kick needed for phase-matched four-wave mixing. High-frequency components scatter into lower-frequency bound modes that are trapped by the Higgs/Skyrme potential, transferring energy inward and building a self-consistent soliton.

**Equations** (phenomenological cubic interaction):
```
F_mix[k] = δ · Σ_{p,q,r} M[p]M[q]M[r] * exp(i Δk·x)   where Δk supplied by ∇ρ
```
(Practical discrete version: compute local triple-product of Fourier-filtered high-k bands and modulate by local |∇ bulk_norm|.)

**Implementation steps**:
1. Add `delta_mix` strength to Config.
2. In field_forces, after existing drives, add a slow cubic term on rotor components only (or full M).
3. To keep cost low: use two-band split (low-k carrier vs high-k tail) and modulate the low-k drive by local density gradient.
4. Verify phase-matching condition via diagnostic freq spectrum.
5. Add E_mix to energy (½ δ |M³|² term).

**Free parameters**: δ_mix (coupling), frequency cutoffs for bands.  
**Tuning plan**: Start with δ=0 (baseline), ramp in small steps while monitoring bound-mode growth in FFT of rotor_rms(t). Analytic estimate from χ⁽³⁾ susceptibility of the Higgs potential; numerical tuning via 5–10 short L=20 runs.

**Likelihood**: 35%. Powerful in optics (Kerr solitons) but requires careful band filtering on foam; risk of numerical instability or energy non-conservation if not derived from a potential. Promising if combined with 7.

### 6. Topological Gradient Force (Winding–Density Coupling)
**Mechanism**: The Skyrme winding number (or Hopf invariant on the rotor q) sources an emergent gauge field that couples to density gradients, producing a Lorentz-like force that pulls topological charge toward high-ρ regions and prevents unwinding.

**Equations**:
```
B_topo ~ ε_ijk ∂_j q · ∂_k q   (discrete curl of left-current L_μ)
F_topo = η · (∇ρ) × B_topo   on the rotor components
```
(η tunable; discretize via existing face gradients of bulk_norm and q.)

**Implementation steps**:
1. Add `eta_topo` to Config.
2. Extend field_forces to compute discrete left-invariant current L_μ = q⁻¹ ∇q on faces (quaternion multiplication on rotor 4-vector).
3. Form topological “magnetic” field from curl of L; multiply by ∇bulk_norm.
4. Add force only to k=0..3; update energy with ½ η |B_topo|² term (like EM energy).
5. Diagnostic: track integrated winding B = (1/24π²) ∫ Tr(L·L·L) per cluster.

**Free parameters**: η (topo–gradient coupling).  
**Tuning plan**: Derive η from requiring that the effective potential for a winding-1 hedgehog has a minimum at observed core density (match to v51 proton binding). Grid search after full Skyrme term is in; 3–4 runs.

**Likelihood**: 60%. Directly extends the existing Skyrme infrastructure; topological protection + gradient force is a classic confinement mechanism (QCD string tension analog). Highest near-term payoff.

### 7. Scale-Breaking Density-Dependent κ
**Mechanism**: The saturation parameter κ in the potential becomes κ(ρ) or κ(|∇ρ|). This explicitly breaks scale invariance in a density-dependent way, creating a preferred wavelength that is shorter (more confined) in high-density regions.

**Equations**:
```
κ_eff[i] = κ0 + ζ · bulk_norm[i] + ξ · |∇ bulk_norm[i]|
V = λ/4 (Mb - v²)² + (1/4) κ_eff (Mb - v²)²   // or on |q|² term
```
(ζ, ξ control linear and gradient dependence.)

**Implementation steps**:
1. Add `kappa0`, `zeta`, `xi` to Config (kappa0 re-uses existing saturation idea from v34–v43).
2. Compute local |∇ bulk_norm| from existing grad arrays in compute_grads_and_lap_all.
3. Modify the quartic drive in field_forces to use κ_eff.
4. Add E_kappa term to diagnostics.
5. Update all run_*.cfg files and the v56 PLAN.md.

**Free parameters**: ζ, ξ (and baseline κ0).  
**Tuning plan**: Analytic: require dV/dρ = 0 at equilibrium core density from v51 proton data. Numerical: CMA-ES on {zeta,xi} with skyrme seed; 10–15 GPU-hours total.

**Likelihood**: 50%. Very close to historical SCP parameter searches (V28 CMA-ES); easy to code and conservative. Risk of re-introducing the old |μ|/m² breathing instability if not careful.

### 8. Geometric (Berry) Phase Orbit
**Mechanism**: As the multivector state adiabatically follows a density gradient, it acquires a geometric phase (Berry curvature) that acts as an effective magnetic field, forcing the excitation into closed orbits or bound states around density peaks.

**Equations**:
```
A_Berry = i ⟨M| ∇_ρ M⟩   (connection in density space)
F_B = ∇_ρ × A_Berry
force_B = F_B × velocity_of_M   (Lorentz from Berry curvature)
```
(Computed via finite-difference overlap of neighboring cells’ M states.)

**Implementation steps**:
1. Add `berry_strength` knob.
2. In a post-force pass, compute local Berry connection using 6-neighbor overlaps (already have face topology).
3. Form discrete curl → Berry flux; apply Lorentz kick to velocity.
4. Add E_Berry = ½ |F_B|² diagnostic.
5. Requires small dt or higher-order integrator for adiabaticity.

**Free parameters**: berry_strength (overall scale of curvature).  
**Tuning plan**: Derive from the Hopf fibration geometry of the rotor S³ (or full 7-sphere vacuum manifold). Numerical: start with analytic value for hedgehog profile; refine with 4–6 runs monitoring orbit closure.

**Likelihood**: 25%. Elegant and fundamental (Berry phase is gauge-invariant), but computationally heavy (overlap integrals per cell per step) and adiabatic assumption may break in strongly nonlinear regime. Long-term high reward if implemented.

### 9. Chern–Simons-like Self-Interaction
**Mechanism**: The density gradient sources a Chern–Simons 3-form whose variation produces a self-force on the frequency current, effectively giving the multivector an anyonic or fractional statistics that favors compact, topologically protected configurations.

**Equations** (phenomenological):
```
L_CS = θ · (ρ · ∇ × A + (2/3) A³)   where A ~ rotor current or ∇q
```
Variation yields force term proportional to ρ · (∇ × v) + … (θ level controls statistics).

**Implementation steps**:
1. Add `theta_cs` (statistics parameter) to Config.
2. Define effective A from existing left-current L_μ or velocity field.
3. Add CS contribution to field_forces (3-form evaluated on faces).
4. Careful discrete exterior calculus to keep energy conservation.
5. Diagnostic: integrated CS invariant and anyonic phase accumulation.

**Free parameters**: θ_cs (anyonic angle).  
**Tuning plan**: Set θ_cs = 1/(2π) for semions or π for fermions; scan around these values. Analytic motivation from 2+1D anyon models projected to 3D foam. 3–5 exploratory runs.

**Likelihood**: 30%. Powerful for statistics-driven binding (fractional quantum Hall analog) but requires non-trivial discrete differential forms on Voronoi mesh; highest risk of breaking symplecticity.

### Compatibility & Recommended Combinations

- **High synergy (try together first)**: 6 (Topological gradient force) + 9 (Chern–Simons) — both act on the rotor current and reinforce topological protection; implement full Skyrme L_μ first, then layer both gauge-like terms.
- **Strong pairing**: 3 (Resonance locking) + 7 (Scale-breaking κ(ρ)) — frequency lock + preferred wavelength creates a narrow resonance band that is automatically confined.
- **Good stacking**: 4 (Density-dependent mass) + 6 (Topo force) — variable inertia slows escape while topology prevents unwinding.
- **Exploratory**: 5 (Four-wave mixing) + 7 — nonlinear scattering feeds the scale-selected modes.
- **Longer-term**: 8 (Berry) can be combined with any of the above once the mesh overlap machinery exists; it is largely orthogonal.

**Next actions**: Implement mechanisms 6 + 7 first (lowest risk, highest leverage on existing Skyrme code). Run skyrme seed with combined knobs on L=40; if bound state forms, add 3 and 4. Full backprop (CMA-ES) on the combined parameter set once any single mechanism shows >20% lifetime improvement.

Document version 1.0 — ready for v56/PLAN.md integration.
