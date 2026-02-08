# Geometric-to-Physical Mapping: Verification Framework

**Question**: CHPT assigns geometric properties (winding number, chirality, orientation) to physical observables (charge, particle type, handedness of EM). How do we ensure these assignments are correct? What framework can verify the mapping?

---

## 1. The Current Assignments

The spec currently proposes:

| Geometric Property | Physical Observable | Source |
|---|---|---|
| Winding number $Q \in \pi_3(S^3) = \mathbb{Z}$ | Electric charge | `math/02_topology.md` |
| Chirality $\chi \in \{L, R\}$ | Particle type (matter/antimatter) | `05_chirality.md` |
| Pseudoscalar orientation (volume form $I = e_1 e_2 e_3$) | Right-hand rule / EM orientation | implicit in `math/04` |
| Topological current $B^\mu$ | Electric current $J^\mu$ | `08_electromagnetism.md` |

---

## 2. A Critical Internal Tension

**The Skyrme model identifies winding number with baryon number, not electric charge.**

This is the most important issue in the current spec. The standard Skyrme model — which CHPT's bulk sector has been *proven* to reduce to — maps:

| Skyrme Model | CHPT (current spec) |
|---|---|
| $B$ (winding number) = **baryon number** | $Q$ (winding number) = **electric charge** |
| Electric charge = $I_3 + B/2$ (derived) | Baryon number = ??? |
| Proton: $B=1, I_3 = +1/2, Q_{em} = +1$ | Proton: $Q=+1$? |
| Neutron: $B=1, I_3 = -1/2, Q_{em} = 0$ | Neutron: $Q=0$? |

In the Skyrme model:
- The **proton** and **neutron** are *the same soliton* ($B=1$), distinguished only by their orientation in isospin space (SU(2) collective coordinate quantization gives $I_3 = \pm 1/2$)
- Electric charge is a *derived* quantum number: $Q_{em} = I_3 + B/2$
- The winding number $B$ counts *baryons*, not electric charge

**CHPT identifies $Q$ with electric charge instead.** If this is correct, then:
- $Q=+1$ is a positron (or proton?) — which one?
- $Q=0$ is a neutral particle — but both the neutron and photon are neutral, and they're vastly different
- $Q=-1$ is an electron (or antiproton?) — which one?

The spec itself reveals the confusion: `math/02_topology.md` lists both "$Q=+1$: Proton/Positron?" and "$Q=-1$: Electron/Antiproton?" (with question marks).

**This is not a minor labeling issue.** It determines the entire physical interpretation of the theory.

### Three Resolutions

**Resolution A — Follow Skyrme (winding = baryon number)**

Accept the standard result: $B$ = baryon number. Electric charge emerges from collective coordinate quantization. This is well-established and has been verified experimentally (the Skyrme model reproduces nucleon properties to ~30%).

- Pro: Proven to work. Gives proton/neutron as states of a single $B=1$ soliton.
- Con: Requires quantizing collective coordinates (semiclassical). Leptons (electron, neutrino) don't emerge from the Skyrme sector.

**Resolution B — Multiple topological charges**

The full CHPT field $\Psi \in Cl^+(3,0,1)$ has 8 components. The bulk quaternion gives $\pi_3(S^3) = \mathbb{Z}$ (one integer invariant). But if the degenerate sector is made dynamical (Roadmap item 7), the vacuum manifold enlarges, potentially supporting *additional* independent topological charges:

- $B \in \mathbb{Z}$ from the bulk quaternion $\to$ baryon number
- $Q_{em}$ from a second topological invariant involving the degenerate sector $\to$ electric charge
- Possibly others $\to$ lepton number, etc.

This requires the full 8-component field to have a richer homotopy structure than $S^3$ alone.

- Pro: Could accommodate both baryon number and electric charge as independent topological quantities.
- Con: The degenerate sector is currently non-dynamical ($e_0^2 = 0$ kills its kinetic energy). A richer vacuum manifold requires new physics beyond the current Lagrangian.

**Resolution C — Reinterpret the target space**

Perhaps the $S^3$ winding number is *neither* baryon number nor electric charge, but a more fundamental quantity from which both are derived. For instance, if the fundamental unit is charge $e/3$ (the quark charge quantum), then:
- Electron: $Q = -3$ units $\to$ charge $-e$
- Up quark: $Q = +2$ units $\to$ charge $+2e/3$
- Down quark: $Q = -1$ unit $\to$ charge $-e/3$
- Proton (uud): $Q = 2+2-1 = +3$ units $\to$ charge $+e$

This requires the fundamental soliton to carry charge $e/3$, and baryons to be $B=3$ composite solitons.

- Pro: Explains charge quantization in units of $e/3$. Natural place for quarks.
- Con: The $B=1$ soliton would be a single quark — but quarks are confined. Need a confinement mechanism.

---

## 3. The Verification Framework

Regardless of which resolution is chosen, the mapping can be tested against five independent criteria. **All five must be satisfied simultaneously.** Failure of any one invalidates the mapping.

### Test 1: Topological Census

**Enumerate all independent topological invariants of the field.**

For a field $\Psi: S^3 \to \mathcal{M}$ (where $\mathcal{M}$ is the vacuum manifold):

- Compute $\pi_3(\mathcal{M})$ — this gives the classification of static solitons in 3D.
- Each independent $\mathbb{Z}$ factor in $\pi_3(\mathcal{M})$ gives one conserved integer charge.
- Each $\mathbb{Z}_n$ factor gives a charge conserved mod $n$.

**Required**: The number of independent topological charges must be $\geq$ the number of independently conserved charges in nature (electric charge, baryon number, lepton number, ...).

Currently: $\mathcal{M} = S^3$ (bulk quaternion), so $\pi_3(S^3) = \mathbb{Z}$ — only **one** independent charge. This is insufficient for the full Standard Model but may suffice if electric charge and baryon number are not independent (as in Resolution C).

**Action**: Determine the vacuum manifold of the full 8-component field (once the degenerate sector dynamics are specified) and compute its homotopy groups.

### Test 2: Discrete Symmetry Matching (C, P, T)

**Verify that the geometric operations corresponding to C, P, T have the correct effect on all mapped quantities.**

| Symmetry | Geometric operation | Effect on physical quantities |
|---|---|---|
| C (charge conjugation) | $q \to \bar{q}$ (quaternion conjugate) | Reverses electric charge, preserves mass, reverses baryon number |
| P (parity) | $\mathbf{x} \to -\mathbf{x}$ | Reverses chirality (handedness), preserves charge and mass |
| T (time reversal) | $t \to -t$ | Reverses velocities and magnetic field, preserves charge and mass |
| CPT (combined) | All three | Must be an exact symmetry of the Lagrangian |

**Specific checks**:

1. If winding number = electric charge, then $q \to \bar{q}$ must reverse the winding. Verify: for the hedgehog $q(r) = \cos f(r) + i\hat{n}\cdot\vec{\sigma}\sin f(r)$, conjugation gives $\bar{q}(r) = \cos f(r) - i\hat{n}\cdot\vec{\sigma}\sin f(r)$. The winding integral $B = -\frac{1}{2\pi^2}\int \epsilon_{ijk} \text{Tr}(A_i A_j A_k) d^3x$ where $A_i = \bar{q}\partial_i q$: under $q \to \bar{q}$, $A_i \to q\partial_i\bar{q} = -\bar{q}^{-1}A_i\bar{q}^{-1} \cdot (\bar{q}\bar{q})$... this needs careful computation.

2. The **chirality** operator $\chi = \text{sgn}(\langle \Psi^\dagger \nabla \Psi \rangle_4)$ must transform correctly under P: $\chi \to -\chi$ (parity flips handedness). Under C: the spec says C reverses both charge and chirality, making the positron ($Q=+1, \chi=R$) the CPT conjugate of the electron ($Q=-1, \chi=L$). This must be verified algebraically.

3. **CP violation**: The Standard Model violates CP (observed in kaon/B-meson decays). Can the CHPT Lagrangian accommodate CP violation? If the Lagrangian is CP-symmetric, CP violation must come from the vacuum or soliton dynamics.

### Test 3: Inter-Soliton Force Matching

**Compute the force between two solitons and verify the sign and scaling match the physical force between the corresponding particles.**

For well-separated solitons, the inter-soliton potential is determined by the tail of the soliton field. In the standard Skyrme model:

- The pion tail falls as $\sim e^{-m_\pi r}/r$ (Yukawa)
- The force between two $B=1$ solitons depends on their relative *isospin orientation*:
  - Attractive channel: deuteron forms (proton-neutron bound state)
  - Repulsive channel: proton-proton (after accounting for Coulomb)

**For the CHPT mapping to work**:

1. **Like-charge solitons must repel** (at long range, $\sim 1/r^2$). This means the inter-soliton force from the massless bivector ($\mathbf{F}$) tail must be repulsive for same-sign winding numbers.

2. **Unlike-charge solitons must attract** ($\sim 1/r^2$).

3. **The force must scale as $1/r^2$** at large separations (Coulomb law), not exponentially (which would indicate a massive mediator / short-range force).

4. **Moving solitons must produce magnetic-type forces** — a soliton moving past another should produce a velocity-dependent force perpendicular to both the velocity and the separation vector (the Lorentz force: $F = qv \times B$).

**This is the most computationally accessible test.** It can be checked numerically:
- Place two solitons at large separation
- Compute the energy as a function of separation and relative orientation
- Extract the force: $F(r) = -dE/dr$
- Verify sign, scaling, and angular dependence

### Test 4: Semiclassical Quantization (Spin-Statistics)

**Quantize the soliton's collective coordinates and verify the quantum numbers match physical particles.**

The $B=1$ Skyrmion has zero modes corresponding to:
- 3 spatial translations (give momentum)
- 3 spatial rotations (give spin)
- 3 isospin rotations (give isospin/flavor)

Quantizing these gives a tower of states with specific $(J, I)$ quantum numbers. In the standard Skyrme model:

- $B=1$: lowest state has $J = I = 1/2$ (proton and neutron, 4 states)
- $B=2$: lowest state has $J = I = 0$ or $J = I = 1$ (deuteron, etc.)

**The Finkelstein-Rubinstein constraint**: A $2\pi$ rotation of a $B=1$ Skyrmion in combined space-isospace is topologically equivalent to the identity but picks up a phase $(-1)^B$. For odd $B$: the soliton is a **fermion** (spin-1/2). For even $B$: **boson** (integer spin).

This is one of the deepest results in the Skyrme model — it connects topology to quantum statistics without any additional postulates.

**For CHPT**: This test is crucial because it determines whether the soliton has the correct spin. If the $B=1$ soliton must represent the electron (charge $Q=-1$, spin-1/2), then the F-R constraint must give half-integer spin for $B=1$. This is indeed what happens in the Skyrme model, regardless of which physical particle $B=1$ represents.

### Test 5: Asymptotic Field Identification (Long-Range Fields)

**Extract the long-range field of the soliton and verify it matches the known gauge field.**

A charged soliton at rest should produce, at large $r$:

$$\mathbf{F}(r) \sim \frac{Q}{4\pi r^2}\,\hat{r} \quad \text{(Coulomb field)}$$

where $\mathbf{F}$ is the bivector perturbation identified with the EM field.

A moving soliton (velocity $v$) should additionally produce:

$$\vec{B} \sim \frac{Q \vec{v} \times \hat{r}}{4\pi r^2} \quad \text{(Biot-Savart field)}$$

**This is where the right-hand rule enters.** The direction of $\vec{B}$ relative to $\vec{v}$ and $\hat{r}$ is determined by the cross product, which depends on the choice of orientation (volume form).

**The right-hand rule is a convention.** The physical content is:
- The *existence* of a velocity-dependent force perpendicular to both $v$ and $B$ (dynamical fact)
- The *relative* sign: two parallel currents attract (dynamical fact)
- The specific direction is fixed by the convention $\epsilon_{123} = +1$ (right-hand rule)

In CHPT, the orientation is set by the pseudoscalar $I = e_1 e_2 e_3$ in $Cl(3,0,1)$. Choosing $I^2 = -1$ (as standard) fixes the right-hand rule. Choosing $I^2 = +1$ would give the left-hand rule. The physics is the same either way — what matters is consistency.

**What IS physical (not conventional)**:
- Like charges repel: two solitons with same-sign $Q$ have positive interaction energy
- Parallel currents attract: two same-charge solitons moving in the same direction attract (weaker repulsion than at rest)
- The magnetic force is perpendicular to velocity: the force on a moving soliton in an external field has a component $\perp v$

**Action**: Compute the linearized field of a single static soliton. Verify it has $1/r^2$ falloff with the correct angular structure (monopole for electric, no magnetic monopole).

---

## 4. The Non-Local Behavior Question

The user specifically asks about "non-local behavior of particles and the behavior of EMF and a moving stream of electrons."

Several phenomena connect topology to non-local EM effects:

### 4.1 The Aharonov-Bohm Effect

A charged particle's wavefunction acquires a phase when encircling a region with magnetic flux, even if $\vec{B} = 0$ everywhere on the particle's path. The phase is:

$$\Delta\phi = \frac{q}{\hbar}\oint \vec{A}\cdot d\vec{l} = \frac{q\Phi_B}{\hbar}$$

In CHPT: the soliton is a topological object, and its interaction with the EM field $\mathbf{F}$ is through the *potential* (connection), not just the field strength. The topological winding of the soliton naturally couples to the gauge potential's topology (the holonomy of the connection around a loop). This is a deep connection between the topology of the soliton (charge) and the topology of the gauge field (flux), and it's *automatic* in any topological theory — it doesn't need to be imposed.

### 4.2 The Right-Hand Rule for Current-Carrying Wires

A stream of electrons (current $I$) produces a magnetic field that wraps around the wire:

$$\vec{B} = \frac{\mu_0 I}{2\pi r}\,\hat{\phi}$$

The direction (clockwise vs counterclockwise viewed from the direction of current) is the right-hand rule.

In CHPT: each moving soliton (electron) contributes a velocity-dependent field perturbation. The superposition of many moving solitons gives the macroscopic Biot-Savart field. The wrapping direction is determined by:
1. The sign of the topological charge ($Q = -1$ for electrons)
2. The orientation convention ($I = e_1 e_2 e_3$)
3. The direction of motion

**Key insight**: The magnetic field of a current is a *collective effect* of many moving topological objects. Its direction is determined by the SAME orientation convention that defines the cross product. Consistency is the requirement — not a specific handedness.

### 4.3 Electromagnetic Induction

A changing magnetic flux through a loop induces an EMF. The direction of the induced current opposes the change (Lenz's law / Faraday's law).

In CHPT: this should emerge from the coupled dynamics of the bivector field $\mathbf{F}$. Maxwell's equations in GA form: $\nabla \mathbf{F} = J$. The free-field version $\nabla \mathbf{F} = 0$ has been derived from linearization. The coupling of $\mathbf{F}$ to soliton currents (the sourced equation) would complete this picture.

---

## 5. A Specific Framework: The Consistency Ladder

I propose the following ordered verification protocol. Each step builds on the previous, and failure at any step requires revisiting earlier choices.

### Level 0 — Algebraic Consistency

Verify that the proposed mapping is **self-consistent at the algebraic level**:
- C reverses the mapped charge (check: $q \to \bar{q}$ reverses winding number)
- P reverses the mapped chirality (check: $\mathbf{x} \to -\mathbf{x}$ flips handedness)
- CPT is a symmetry of the Lagrangian
- Charge conjugation relates the proposed particle to its antiparticle

**Status**: Partially checked. The winding number does reverse under $q \to \bar{q}$. CPT of the Lagrangian has not been explicitly verified.

### Level 1 — Topological Counting

Verify that the **number of independent topological charges** matches the number of independent conserved charges in nature:
- At minimum: electric charge + baryon number = 2 independent charges
- Currently: $\pi_3(S^3) = \mathbb{Z}$ gives only 1

**Status**: INSUFFICIENT with current vacuum manifold. Either adopt Resolution A (winding = baryon number, electric charge derived) or Resolution B/C (richer manifold).

### Level 2 — Static Properties

Compute **static properties** of the soliton and compare:
- Mass (energy at rest) — match to one particle fixes the mass scale
- Topological charge — verify integer quantization
- Long-range field — verify $1/r^2$ Coulomb-like falloff of the bivector perturbation
- Angular structure — verify monopole (no preferred direction for a single charge at rest)

**Status**: Mass and charge computed for $B=1$. Long-range field not yet extracted.

### Level 3 — Semiclassical Quantum Numbers

Quantize collective coordinates:
- Spin: verify half-integer for the proposed fermion assignments
- Isospin: verify the correct multiplet structure (proton/neutron doublet, etc.)
- Magnetic moment: compare to measured values

**Status**: Not yet attempted. Requires understanding which collective coordinates the CHPT soliton has.

### Level 4 — Dynamic Properties

Compute velocity-dependent effects:
- Inter-soliton force at rest (Coulomb or Yukawa?)
- Moving soliton field (Biot-Savart structure?)
- Soliton scattering cross-sections

**Status**: Not yet attempted. Requires 3D time-dependent simulation (Roadmap item 5).

### Level 5 — Emergent Gauge Structure

Verify that the full Maxwell equations (including sources) emerge:
- Sourced equation $\nabla \mathbf{F} = J$ where $J$ is the topological current
- Lorentz force law for soliton motion in external field
- Gauge invariance as a consequence of the field structure

**Status**: Free-field equation derived (`math/04`). Sourced equation is open problem B4.

---

## 6. Recommendations

### Immediate (can be checked now)

1. **Resolve the winding-number assignment** — Is it baryon number (standard Skyrme) or electric charge (current spec)? This is the single most important theoretical decision. I recommend:

   **Adopt the standard Skyrme identification $B$ = baryon number** as the default, because:
   - The CHPT bulk sector is *proven* to reduce to the Skyrme model
   - The Skyrme identification has been verified against experiment
   - Electric charge can be derived from collective coordinate quantization
   - The alternative (winding = electric charge) creates the proton/neutron distinction problem

   However, note this makes the electron problematic: in the Skyrme model, there are no leptonic solitons. The electron must come from somewhere else (perhaps the degenerate sector, once it's made dynamical).

2. **Check C, P, T algebraically** — Explicitly compute the effect of $q \to \bar{q}$, $\mathbf{x} \to -\mathbf{x}$, and $t \to -t$ on the Lagrangian, the winding number, and the chirality $\chi$. This is a pure algebra exercise.

3. **Extract the long-range field** — From the existing $B=1$ soliton profile $f(r)$, compute the linearized bivector perturbation $\mathbf{F}$ at large $r$. Verify it's Coulomb-like ($1/r^2$) with the correct angular structure.

### Near-term (requires some development)

4. **Compute the inter-soliton potential** — Place two $B=1$ solitons at varying separations and orientations (using the product ansatz or the rational map results). Compute $E(r, \theta)$. Determine whether the force is attractive, repulsive, or orientation-dependent.

5. **Identify collective coordinates** — The $B=1$ CHPT soliton has spatial translations, spatial rotations, and isorotations. Count the zero modes. If the soliton is in the Skyrme model sector, this is already known: 6 translational/rotational + 3 isorotational = 9 zero modes total (for $B=1$).

### Long-term (foundational)

6. **Make the degenerate sector dynamical** — The current Lagrangian gives $P$ and $\vec{J}$ no kinetic energy ($e_0^2 = 0$ kills it). Extending the Lagrangian (Roadmap item 7) could enlarge the vacuum manifold, creating additional topological charges for electric charge and/or lepton number.

7. **Derive the sourced Maxwell equation** — Extract the source current $J$ from the nonlinear soliton and show that $\nabla \mathbf{F} = J$ holds, with $J$ proportional to the topological current.

---

## 7. Summary: What Determines the Correct Mapping?

The mapping is determined by **five independent consistency conditions** (Tests 1–5 above). Unlike a postulate, the mapping is *not a free choice* — it is forced by the mathematics once the Lagrangian is fixed. The procedure is:

1. Fix the Lagrangian and vacuum manifold.
2. Compute the homotopy groups $\to$ this determines how many independent topological charges exist.
3. Compute the static soliton $\to$ this determines its mass, charge, and long-range field.
4. Quantize collective coordinates $\to$ this determines spin, isospin, and the spin-statistics relation.
5. Compute the inter-soliton potential $\to$ this determines which charge maps to which force.

The mapping that satisfies all five simultaneously is the correct one. There is no freedom to "choose" — the geometry dictates the physics.

**The right-hand rule** specifically is resolved at step 5: the orientation convention $I = e_1 e_2 e_3$ propagates through the cross product to determine the direction of magnetic forces. The convention is arbitrary (right-hand vs left-hand), but once chosen, all electromagnetic phenomena must use the same convention consistently. The physical content (like charges repel, parallel currents attract) is independent of the convention.

---

## 8. Case Study: The Neutrino

The neutrino is the most instructive test case for the mapping framework, because it exposes the insufficiency of a single geometric property.

### 8.1 What We Know Experimentally

| Property | Value | Implication for CHPT |
|---|---|---|
| Electric charge | 0 | Winding number $Q = 0$ (if $Q$ = electric charge), or isospin cancels hypercharge |
| Spin | 1/2 | Must be a fermion — Finkelstein-Rubinstein constraint requires odd-$B$ (if in Skyrme sector) |
| Weak isospin | $I_W = 1/2$ (left-handed only) | Couples to W/Z bosons; has a non-trivial chirality |
| Chirality | Left-handed only | Only $\nu_L$ exists; $\nu_R$ has never been observed |
| Mass | ~0.05 eV (nonzero but tiny) | Neutrino oscillations prove $m \neq 0$; origin unknown even in SM |
| Interactions | Weak + gravity only | No EM, no strong; cross-section ~10^-38 cm² |
| Baryon number | 0 | It's a lepton, not a baryon |

### 8.2 The Achiral Neutrino (Rejected)

The original CHPT proposal mapped the neutrino to an achiral knot: $Q = 0$, $\chi = 0$. This picture explains why neutrinos pass through matter (no geometric interference). But it fails quantitatively:

- An achiral knot interacts only gravitationally (Rule 3).
- Gravitational cross-sections are ~10^-63 cm², but neutrino cross-sections are ~10^-38 cm².
- The discrepancy is **25 orders of magnitude**. This is not a small correction — the achiral picture is categorically wrong for neutrino interactions.

### 8.3 The Chiral Neutrino (Required)

The neutrino must be chiral — it must have some geometric property that enables weak-force coupling. But it must simultaneously lack whatever geometric property produces electromagnetic coupling.

This means CHPT needs **at least two independent geometric axes**:

```
                    Weak chirality (χ_W)
                         |
              ν_L        |        e_L
           (0, L)        |      (-1, L)
                         |
    ——————————————————————+———————————————————  Electric charge (Q)
                         |
              γ           |        e_R
           (0, 0)        |      (-1, 0)
                         |
```

- **Neutrino** ($\nu_L$): $Q = 0$, $\chi_W = L$ — no EM interaction, yes weak interaction
- **Left-handed electron** ($e_L$): $Q = -1$, $\chi_W = L$ — yes EM, yes weak
- **Right-handed electron** ($e_R$): $Q = -1$, $\chi_W = 0$ — yes EM, no weak (W-coupling)
- **Photon** ($\gamma$): $Q = 0$, $\chi_W = 0$ — no EM charge, no weak charge

The neutrino and the right-handed electron are the proof that these axes are independent.

### 8.4 Where Does Weak Chirality Live in the Algebra?

$Cl^+(3,0,1)$ has 8 components: 4 bulk ($s, f_1, f_2, f_3$) and 4 degenerate ($j_1, j_2, j_3, p$).

The bulk quaternion $q \in S^3$ provides the winding number (one topological charge). Weak chirality must come from somewhere else. The natural candidate is the **degenerate sector** $(J, P)$:

- Currently non-dynamical ($e_0^2 = 0$ kills its kinetic energy).
- If made dynamical (Roadmap item 7), could support its own topological structure.
- The degenerate sector has 4 components (a quaternion-like structure), potentially giving $\pi_3(S^3) = \mathbb{Z}$ again — a second independent topological charge.

**Speculative picture**: The full dynamical field $\Psi = q + (J + Ip)e_0$ could have a vacuum manifold $\mathcal{M} \cong S^3 \times S^3$ (bulk $\times$ degenerate), with:
- $\pi_3(S^3 \times S^3) = \mathbb{Z} \times \mathbb{Z}$ — **two independent integer charges**
- First charge from bulk: baryon number (or electric charge)
- Second charge from degenerate: weak isospin (or lepton number)

This is exactly the structure needed. But it requires the degenerate sector to have its own kinetic term and vacuum constraint, which the current Lagrangian does not provide.

### 8.5 Why Neutrinos Pass Through Matter

Even with chirality, neutrinos interact very weakly. The explanation is **range**, not **strength**:

- The electromagnetic force is mediated by the massless photon: infinite range, $1/r^2$ falloff.
- The weak force is mediated by the massive W/Z bosons (~80–91 GeV): range ~$\hbar/(m_W c) \approx 10^{-18}$ m.
- For a neutrino to interact weakly with a nucleus, it must pass within ~10^-18 m of the nucleus.
- A typical atom is ~10^-10 m across. The probability of passing within 10^-18 m is ~$(10^{-18}/10^{-10})^2 = 10^{-16}$.
- This geometric dilution explains the tiny cross-section without requiring the neutrino to be achiral.

In CHPT terms: the neutrino's weak-chiral knot structure is just as "real" and geometrically nontrivial as the electron's. But the weak interaction only operates when two knots' cores overlap (short range), while the electromagnetic interaction operates through the long-range bivector field $\mathbf{F}$ (infinite range). The neutrino is transparent not because it lacks structure, but because its structure couples only to a short-range field.

### 8.6 Applying the Five Tests to the Neutrino

| Test | Requirement | Status |
|---|---|---|
| 1. Topological Census | Need an independent charge for weak isospin | Requires dynamical degenerate sector |
| 2. C, P, T | C: $\nu_L \to \bar{\nu}_R$. P: $\nu_L \to \nu_R$ (doesn't exist!). P-violation built in | Must verify algebra accommodates maximal P-violation |
| 3. Inter-soliton forces | $\nu$–$e$ scattering via weak force only; no EM component | Requires short-range force from massive mediator |
| 4. Quantization | Must give spin-1/2, $I_W = 1/2$, $Q_{em} = 0$ | Requires collective coordinate analysis |
| 5. Asymptotic fields | No long-range EM field (Q=0); no long-range weak field (massive mediator) | Must verify soliton has no $1/r^2$ tail |

The neutrino is particularly constraining for Test 2: parity maps $\nu_L \to \nu_R$, but $\nu_R$ does not exist in nature. This means parity is **maximally violated** in the neutrino sector. The CHPT Lagrangian must either explicitly break P (Option A in `05_chirality.md`) or have a vacuum that breaks P (Option B).

---

## Appendix: Connections to Known Results

### Skyrme Model (established)
- $B$ = baryon number, verified to ~30% accuracy for nucleon properties
- Proton/neutron from $B=1$ collective coordinate quantization
- Deuteron from $B=2$ soliton
- Nuclear binding energies qualitatively correct
- No leptons in the theory

### 't Hooft–Polyakov Monopole
- Gauge theory soliton with magnetic charge
- Shows that topological charge can correspond to *magnetic* charge, not just electric
- Relevant if CHPT wants magnetic monopoles

### Faddeev-Niemi Hopfions
- Solitons in $S^2$-valued fields with Hopf invariant
- $\pi_3(S^2) = \mathbb{Z}$ (same as $\pi_3(S^3)$, by Hopf fibration)
- Energy scales as $E \propto Q^{3/4}$ (different from Skyrme's $E \propto Q$)
- Could be relevant if CHPT's effective target space is $S^2$ rather than $S^3$ in some limit

### Witten's Result (1983)
- In large-$N_c$ QCD, baryons are solitons of the pion field (= Skyrmions)
- The Skyrme model is the *correct* low-energy effective theory for baryons
- This provides the theoretical foundation for identifying winding number with baryon number
