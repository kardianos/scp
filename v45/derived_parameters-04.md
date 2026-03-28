# Derived Parameters for SCP Field Theory — Revision 04
## First Principles Recap
The theory rests on a single 3-component relativistic field φ with triple-product potential V(P=φ₀φ₁φ₂) plus Cosserat θ coupling. Stable structures must be true attractors satisfying virial identities derived from the Lagrangian. Free parameters break predictivity; they must emerge from consistency (virial balance, conservation laws, soliton equilibria). V45's null binding result plus virial mismatch (gradient energy >> binding term) shows current tuned values produce meta-stable "exploding" solitons rather than true bound states.
## Agreements with Prior Revisions
- **-01 & -02**: Strong agreement on reducing free parameters via virial theorem, self-tuning mechanisms, and η-sweep to find binding. Virial relation tying μ/κ to integrated energies is genuine first-principles constraint. Seed quality (gen_deuterium vs templates) is critical.
- **-03**: Excellent diagnosis of "exploding solitons" and phase-incoherence. Topological quench from compressed high-A_bg state + annealing is the right shift from kinematics to dynamics. Machian global A_bg constraint for automatic attraction is elegant.
## Disagreements
- -01 derivations are partly circular (P_core depends on μ/κ). 
- -03's specific jumps (m²=1.0, μ=-250, κ=20, Coriolis Ω term) lack derivation and would require kernel changes (forbidden without explicit request). Non-local or CA ideas are too radical and break Lorentz invariance.
- Treating phases as purely empirical is fine; exact 2π/3 may not be superior to CMA-ES values.
## What Might Work
- Virial-constrained μ/κ ratio computed from existing SFA data (E_grad, E_mass, E_pot) to check self-consistency.
- Compressed-spin quench without kernel mods: high A_bg=0.4 initial + imprinted velocity field via new seed generator, followed by annealing (damping ramp-down).
- η sweep on template-seeded D=40 runs; pre-converged baryons from V43.
- Dynamic but Lagrangian-preserving extensions (mode=3 density-κ globally).
## New Ideas
1. **Adiabatic effective potential**: Average fast carrier oscillations to derive a slow effective inter-baryon potential; fix μ from matching its minimum to nuclear scale.
2. **Topological charge conservation**: Treat carrier phases as topological invariants; derive δ exactly from winding numbers rather than optimization.
3. **Integrated-P constraint**: Add global Lagrange multiplier enforcing ∫P dV = constant; this creates Machian binding without local parameter changes.
4. **Coarse-graining**: Derive collective mediator-like force by integrating short-wavelength modes; small effective m_θ from nonlinear mixing could yield Yukawa attraction without explicit mass term.
5. **Breathing-phase locking**: Measure Adler-type synchronization energy between two baryons' breathing modes; set μ so locking energy exceeds thermal noise at observed periods (~150t).
## v46 Proposal
Implement virial checker on V45 data. Write quench seed (high-compression rotating background). Run η-varied template-seeded binding tests. Update CONCEPT.md to replace tuned constants with virial-derived relations and global constraints. Goal: true energy minimum at finite D without free parameters beyond overall scale.
This revision synthesizes the series into an actionable, kernel-safe path.