# CHPT — Unified Summary

Chiral Harmonic Particle Theory proposes that all of physics arises from a single continuous field in 3+1 dimensions. Particles are stable topological configurations ("knots") of this field, forces are density gradients, and radiation consists of propagating field oscillations ("null-rotors"). The theory aims to be a complete ontology: there is nothing else.

---

## The Field

The field $\Psi$ is valued in $Cl^+(3,0,1)$, the even subalgebra of Projective Geometric Algebra. This is isomorphic to the dual quaternions $\mathbb{H} + \varepsilon\mathbb{H}$ and decomposes as:

$$\Psi = \underbrace{(\rho + \mathbf{F})}_{q\;\text{(bulk)}} + e_0\underbrace{(\vec{J} + \tau I_3)}_{p\;\text{(weight)}}$$

The bulk quaternion $q$ carries 4 components: scalar density $\rho$ and spatial bivector $\mathbf{F}$ (field strength / internal twist). The weight part $p$ carries 4 more: flux $\vec{J}$ and pseudoscalar helicity $\tau$. Total: 8 real components.

The degenerate basis satisfies $e_0^2 = 0$. This is not a timelike dimension — time enters as an external evolution parameter, consistent with the theory's process ontology. Space is Euclidean; the geometric derivative $\nabla = \sum e_i \partial_i$ is purely spatial.

---

## Axioms

Six foundational postulates:

1. **One field**: A single continuous field over 3+1D spacetime is the sole fundamental entity.
2. **Density**: The field has a well-defined, non-negative local density $\rho$ at every point.
3. **Conservation**: Total integrated density is constant — density is redistributed, never created or destroyed.
4. **No dissipation**: The dynamics are fundamentally non-dissipative. All energy remains in the field.
5. **Speed limit**: There exists a finite maximum signal speed $c$. Topological non-locality (instantaneous global constraints on connected structures) is permitted, but superluminal signaling is not.
6. **Stable structures**: The field supports spatially localized, self-sustaining configurations (knots) stabilized by topology.

---

## Dynamics

The theory is governed by a four-term Lagrangian with five free parameters $(\rho_0, \lambda, e, \mu, c)$:

$$\mathcal{L} = \underbrace{\tfrac{1}{2}\eta^{\mu\nu}\langle\partial_\mu\Psi\,\widetilde{\partial_\nu\Psi}\rangle_0}_{\mathcal{L}_2} + \underbrace{\tfrac{1}{4e^2}\sum_{\mu<\nu}\langle[R_\mu, R_\nu]^2\rangle_0}_{\mathcal{L}_4} - \underbrace{\tfrac{\lambda}{4}(\langle\Psi\tilde\Psi\rangle_0 - \rho_0^2)^2}_{V} - \underbrace{\tfrac{\mu^2}{2}(|\vec{J}|^2 + \tau^2)}_{V_D}$$

- **$\mathcal{L}_2$** (kinetic): Quadratic gradient energy. Provides wave propagation.
- **$\mathcal{L}_4$** (Skyrme): Fourth-order derivative term from commutators of right-currents $R_\mu = \tilde\Psi\partial_\mu\Psi$. Prevents soliton collapse (evades Derrick's theorem).
- **$V$** (bulk potential): Mexican-hat potential constraining the bulk norm $|q|^2 \to \rho_0^2$ at spatial infinity.
- **$V_D$** (degenerate mass): Gives mass $\mu$ to the weight-sector modes $(\vec{J}, \tau)$, which are invisible to the standard Clifford norm since $e_0^2 = 0$.

The equation of motion is hyperbolic (wave-like) and Lorentz covariant:

$$\frac{1}{c^2}\partial_t^2\Psi - \nabla^2\Psi + \frac{\delta\mathcal{L}_4}{\delta\Psi} + \lambda\Psi(\langle\Psi\tilde\Psi\rangle_0 - \rho_0^2) = 0$$

---

## Particles as Knots

Particles are topologically stable field configurations. The vacuum manifold $\langle\Psi\tilde\Psi\rangle_0 = \rho_0^2$ is $S^3$, and maps from compactified space $S^3 \to S^3$ are classified by the homotopy group $\pi_3(S^3) = \mathbb{Z}$. The integer winding number $Q$ is identified with baryon number (or electric charge). It is exactly conserved because continuous deformations cannot change the homotopy class.

A knot is not a static object — it is a standing process of "knotting and unknotting" in the field. Motion is the propagation of this organizational pattern, not the movement of substance through a medium. This process ontology is why CHPT is not an aether theory: there is no object-medium duality, no preferred frame, and the Michelson-Morley null result is automatic.

**Mass** is the total integrated field energy of the configuration: $Mc^2 = E_2 + E_4 + E_V + E_D$. The virial theorem $E_2 - E_4 + 3E_\text{pot} = 0$ constrains the equilibrium, giving the mass formula $Mc^2 = 2E_4 - 2E_V - 2E_D$. Heavier particles correspond to more complex topological structures requiring more field energy to sustain.

---

## Linear Limit and Electromagnetism

Linearizing around the vacuum $\Psi = \rho_0 + \psi$ yields four propagating sectors:

| Mode | Mass | Identification |
|------|------|----------------|
| $\mathbf{F}$ (spatial bivector) | $0$ | Electromagnetic field |
| $S$ (scalar) | $\sqrt{2\lambda}\,\rho_0$ | Massive scalar (Higgs-like) |
| $P$ (pseudoscalar) | $\mu$ | Massive pseudoscalar |
| $\vec{J}$ (degenerate bivector) | $\mu$ | Massive flux (3 components) |

The massless bivector $\mathbf{F} = \vec{E} + I\vec{B}$ satisfies the free electromagnetic wave equation $\square\mathbf{F} = 0$, with the null condition $\mathbf{F}^2 = 0$ yielding $|\vec{E}| = |\vec{B}|$ and $\vec{E} \perp \vec{B}$. This is standard EM wave dynamics derived from the field structure, not postulated.

---

## Chirality, Charge, and Interactions

Two independent geometric properties classify knots:

- **Topological index $Q$** (winding number): Maps to electric charge. Quantized as integers from $\pi_3(S^3)$.
- **Chirality $\chi$** (handedness): Maps to weak interaction coupling. Independent of $Q$ — neutrinos prove this ($Q = 0$ but $\chi = L$).

Interaction rules follow from chirality:
- **Same chirality**: Repulsion (geometric incompatibility of density patterns). Maps to the Pauli exclusion principle.
- **Opposite chirality**: Attraction/nesting (constructive interference). Maps to matter-antimatter pairing and hadronic binding.
- **Achiral**: Weak interaction with other structures. Maps to photon-like transparency.

---

## Gravity

Gravity arises from density depletion zones. A knot concentrates density locally, creating a surrounding deficit (by conservation). Another knot in the vicinity experiences the pressure differential — higher-density surrounding field pushes it toward the depletion zone. This produces:

- **$1/r^2$ force law**: Geometric dilution of the deficit over a sphere.
- **Universality**: Every knot depletes the field, regardless of internal structure.
- **Attractive only**: Depletion always reduces local density below $\rho_0$.
- **Weakness**: An indirect, second-order geometric effect compared to direct knot-knot interactions (EM, strong).

The theory avoids the fatal problems of Le Sage push-gravity (heating, shielding, speed) because gravity comes from a static field gradient, not from absorption of discrete mediators.

For quantitative agreement with GR, CHPT must produce the correct post-Newtonian corrections (PPN parameters). The most promising path is showing that density gradients create an effective metric for null-rotor propagation that satisfies (or closely approximates) Einstein's field equations.

---

## Quantum Mechanics

CHPT is a deterministic field theory that addresses Bell's theorem via topological non-locality. Entangled particles are two ends of a single extended topological structure — changes to the global topology are instantaneous constraints that do not permit superluminal signaling. This makes CHPT a non-local realist theory, alongside de Broglie-Bohm mechanics.

- **Wave-particle duality**: A knot is localized (particle-like) but perturbs the surrounding field (wave-like). Detection is localized at whichever point the density patterns first interact.
- **Uncertainty**: Epistemic, not ontological — the knot has definite position and momentum, but the full field configuration is inaccessible.
- **Quantization**: Discrete energy levels correspond to the discrete set of stable resonant configurations.

The deepest challenge: deriving the Born rule ($P = |\psi|^2$) from field dynamics. This likely requires an additional postulate analogous to quantum equilibrium.

---

## Numerical Results

The bulk sector reduces exactly to the standard Skyrme model for static solitons (the "static decoupling theorem" — confirmed numerically). Key computed results:

- **$B=1$ Skyrmion**: $E/E_{FB} = 1.232$, virial $E_2/E_4 = 1.000$, topological charge $Q = 1.000$. Matches the established Skyrme literature exactly.
- **Higher-$B$ solitons** ($B = 2, 3, 4$): Computed via rational map ansatz. Multi-Skyrmions are bound: binding energy per baryon increases from 1.9% ($B=2$) to 7.7% ($B=4$).
- **Finite-$\lambda$ corrections**: Soliton mass decreases monotonically as $\lambda$ decreases (core softening). Below $\lambda \approx 8000$, the soliton collapses — the sigma-model soliton is a local minimum, with the global minimum at $\rho \to 0$.
- **Soliton scattering**: Head-on $B+B$ collision at $v = 0.5c$ shows perfect topology preservation ($Q = 1.9999$) through deep core interpenetration. Solitons decelerate $18\times$ under repulsive interaction. $B + \bar{B}$ annihilation produces inelastic pass-through with 48% mass radiated.
- **Degenerate sector**: Made dynamical via three coupling terms. Extensive bound-state searches (50+ parameter combinations, breathing modes, angular modes) yield no bound states. The degenerate sector's physical role remains unclear.
- **Parameter fitting**: At $e = 1$, $\rho_0 = 1$: 1 code energy = 9.1 MeV, 1 code length = 0.56 fm.

---

## Critical Open Problems

**What works:**
- Field equation established with proper time evolution and soliton stability
- Free-field Maxwell equations derived from the linear limit
- Lorentz invariance built in; special relativity automatic
- Topological charge conservation exact
- $B = 1$--$4$ soliton energies match Skyrme literature
- Equivalence principle holds by construction

**What does not yet work:**
- **SU(3) color charge**: No mechanism produces the three-fold color symmetry of QCD. This is the largest gap with the Standard Model.
- **Gauge symmetry**: The U(1) $\times$ SU(2) $\times$ SU(3) structure is absent. Whether geometric covariance can substitute is unproven.
- **Weak force**: All degenerate-sector bound-state searches are negative. No mechanism for parity violation or the W/Z bosons.
- **Particle spectrum**: Zero particle masses have been matched to experiment. The mapping from topological sectors to specific particles (electron, proton, etc.) is not established.
- **Sourced Maxwell equations**: Free EM waves are derived, but how knots source EM fields ($\nabla\mathbf{F} = J$) is not.
- **Gravitational wave polarization**: The depletion mechanism is scalar; it may not produce the two tensor polarizations observed by LIGO.
- **Born rule**: Not derived from field dynamics.
- **Three generations**: No explanation for why fermions come in three copies.

---

## Honest Assessment

CHPT's mathematical framework — a Skyrme-type Lagrangian on $Cl^+(3,0,1)$ — is well-defined, supports stable topological solitons, and naturally produces massless EM-like waves. The bulk sector is numerically equivalent to the standard Skyrme model, a respected (if incomplete) model of nuclear physics.

However, the theory remains largely promissory. The conceptual framework (knots as particles, topology as charge, density gradients as forces) is suggestive, but zero quantitative predictions match experiment beyond what the Skyrme model already achieves. The nuclear sector (color, confinement, weak force) has no working model. The quantum sector requires non-local hidden variables but has not derived the Born rule. The gravitational sector has qualitative plausibility but no post-Newtonian calculations.

The theory has five parameters and must produce the entire particle spectrum, all force strengths, and all coupling constants from those five inputs plus topology. Whether this is possible is an open mathematical question. The gap between the established equation and testable predictions is the primary challenge — it requires numerical computation of the full soliton spectrum and its mapping to observed particles.
