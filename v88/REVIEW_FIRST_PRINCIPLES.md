# First-principles review of `DESIGN_CLEAN.md`

**Scope.** The design document alone. No project history, no earlier versions, no
appeal to other models or authority. Claims are sorted as **DERIVED**,
**POSTULATED**, or **GUESSED**. Numerical assertions that are not measured here
are not promoted to evidence.

**Implementation status.** No full 3D nonlinear simulation of this design was
run. Code created for this review: **none**. Reasons and a buildable
specification are in §D.

---

## Summary verdict

| Question | Verdict |
|----------|---------|
| Internally consistent as written? | **Partially.** Transfer-time and a speed *limit* of order \(a/\tau\) follow. Complete-cycle gating, isotropic \(c\), mass-as-cell-count with exterior distortion, and 3-from-7 do **not** follow from the stated premises. |
| Does it produce particles (§7)? | **Not shown.** Multi-cell trapping is *postulated* as the definition of a particle, not derived from dynamics. Discrete species and mass ratios are *hopes* pending a closed dynamical law. |
| Binding (the hard problem)? | **Not derived. Likely missing as stated.** Handedness gives two *distinct* channels, not an attractive one. Global \(\sum\theta_i=0\) does not supply local binding on a large fabric. Binding is a success criterion, not a consequence. |
| True experiment run? | **No.** Design underdetermines the equations of motion; any run would add postulates. Spec in §D. |

The single most valuable finding: **§7’s attractive, separation-vanishing binding channel does not follow from the design as written.** If that channel cannot be exhibited after the dynamics are completed, the design fails on its own terms.

---

## A. Coherence

### A.1 “Energy crosses only in complete cycles”

**What is written.**

- Cells do not exchange energy directly; every link carries its own harmonic.
- Energy is deposited into the link, carried through the link’s cycle, then
  delivered to the far cell.
- A partial cycle delivers nothing.
- Non-resonant drive never accepts energy; the gate is structural.

**What follows from link mediation alone (DERIVED).**

1. Transfer is not an instantaneous cell–cell term. There is an intermediate
   degree of freedom on the link.
2. Propagation of energy therefore has a characteristic time scale set by that
   intermediate dynamics (call it \(\tau\)).
3. Frequency-mismatched drive couples weakly if the coupling is near-linear and
   the link is a high-\(Q\) resonator. That is ordinary resonance filtering, not
   a structural all-or-nothing gate.

**What does not follow (gap).**

“Partial cycle delivers nothing” is **strictly stronger** than “transfer is
link-mediated and delayed.”

In ordinary coupled-oscillator or link-mediated lattice models, energy transfer
is continuous in time. Rabi exchange, nearest-neighbour hopping, and
wave-equation transport all move energy during a partial cycle. Completing a
cycle is not a prerequisite for nonzero net transfer; it is only a natural
period of the intermediate mode.

For “partial cycle \(\Rightarrow\) zero delivery” one needs an **extra structure**,
for example one of:

| Extra structure | What it does |
|-----------------|--------------|
| Stroboscopic / latch coupling | Cell–link interaction fires only at a phase mark (once per \(2\pi\)) |
| Cycle-integrated monodromy | Only the holonomy or \(\oint\) of the interaction deposits energy; the integrand is a total derivative mid-cycle |
| Strict action quantisation with selection rules | Only full quanta \(\hbar\omega_{\mathrm{link}}\) move, and only on resonance |
| Ideal rectifier on the link | Negative half-cycles are forbidden by construction |

None of these is implied by “the link has a harmonic.” A harmonic is a linear
(or weakly nonlinear) oscillator; that class *does* transfer energy mid-cycle.

Resonance gating (“cannot be driven resonantly \(\Rightarrow\) never accepts”) is
also incomplete as stated: it is approximately true for weak linear drive of a
sharp resonator, and false in general for strong or nonlinear coupling (virtual
exchange, multi-quantum processes, broadband kicks). Calling the gate
“structural” requires a definition of the allowed transitions of the link’s
configuration space, not only a frequency label.

**Conclusion (A.1).**

- **DERIVED:** link mediation \(\Rightarrow\) delayed, filtered transfer; a time
  scale \(\tau\) exists.
- **POSTULATED (not derived):** complete-cycle-only delivery; structural
  all-or-nothing acceptance.
- **Does the slogan follow from link-mediated transfer?** **No.** It needs a
  further gate, latch, monodromy rule, or quantisation rule. That further rule
  is currently a design hole, not a corollary.

---

### A.2 Does \(c\) emerge as claimed?

**What is written.**

\[
c \sim a \times (\text{rate of cycle coupling})
\]

is “a derived constant,” not postulated. The speed limit is “the fastest that a
link harmonic can complete a cycle and hand energy across.”

**What follows (DERIVED).**

If every energy hop across one cell spacing \(a\) requires a process of duration
at least \(\tau\), then the front speed of energy is bounded by

\[
v_{\mathrm{front}} \;\lesssim\; \frac{a}{\tau}.
\]

Identifying \(c\) with that bound gives \(c \sim a/\tau\). That much is elementary
and does not smuggle a continuum light cone.

**What is smuggled or unfinished.**

1. **The microscopic rate is still an input.** “Rate of cycle coupling” is not
   reduced to a more primitive, uniquely fixed quantity in the text. If the link
   natural frequency \(\omega_L\) (or coupling rate \(\gamma\)) is a free parameter of
   the fabric, then \(c \sim a\,\omega_L/(2\pi)\) merely *renames* that parameter.
   Derivation relative to continuum postulation of \(c\) is fair; derivation from
   nothing is not.

2. **Isotropy is not derived.** A discrete cell complex (especially a cubic
   index lattice with \(x(i)=a\cdot i+\xi(i)\)) generically has direction-dependent
   group velocity. A single scalar \(c\) requires either emergent continuum
   isotropy at long wavelength or an argument that cycle rates compensate the
   lattice anisotropy. Neither is given.

3. **Lorentz structure is not derived.** A speed limit is not special
   relativity. Frame independence of \(c\), light-cone causal structure for all
   signals (including nonlinear and multi-path), and a relativistic mass–energy
   relation are separate claims. The design does not produce them.

4. **Nonlinear front speed.** Even if linear waves have speed \(\sim a/\tau\),
   nonlinear or multi-cycle processes can in principle move energy differently.
   A proof that *all* energy transport respects the same \(c\) is missing.

**Conclusion (A.2).**

- **DERIVED:** a causal speed *limit of order* \(a/\tau\) from cycle-timed hops.
- **NOT DERIVED:** a unique isotropic Lorentzian \(c\); elimination of microscopic
  rate parameters; relativistic kinematics.
- Something *is* smuggled if one reads “\(c\) is derived” as “no free speed
  parameter remains.” The free parameter has been moved into the cycle rate.

---

### A.3 “Mass is a cell count” vs exterior distortion

**What is written.**

- \(\xi\) is the fabric’s degree of freedom; it defines the tessellation and
  “carries mass.” \(\theta=\nabla\cdot\xi\) is local density.
- A particle has (1) interior configuration, (2) **sharp exterior distortion**,
  (3) internal reflection from spectral mismatch.
- **Mass is a cell count.** A particle occupies an integer number of cells.
- On a closed fabric \(\sum_i\theta_i=0\) identically.

**The tension.**

Let \(R\) be the set of cells counted as “the particle” (the multi-cell support
with locked interior cycles). Let \(E\) be the exterior cells with nonzero
distortion (density/rate jump, decaying or not).

Two mass concepts are in play:

| Label | Definition | Status in text |
|-------|------------|----------------|
| \(M_{\mathrm{count}}\) | \(\lvert R\rvert\) (integer cell count) | Stated as *the* mass |
| \(M_\xi\) | Whatever functional of \(\xi\) (or of cyclic energy) “carries mass” | Implied by §1 |

**DERIVED conflict.**

1. If \(\xi\) carries mass everywhere it is nonzero, then cells in \(E\) contribute
   to mass. Then particle mass is not \(M_{\mathrm{count}}\) alone; it includes a
   contribution from the exterior cloud. That contribution is not an integer
   cell count of \(R\) and typically varies continuously with separation, box
   size, and mutual distortion — undermining “mass ratios fixed by configuration
   integers” unless the exterior contribution is proven identical for all
   isolated species or strictly zero.

2. If mass is *only* \(M_{\mathrm{count}}=\lvert R\rvert\), then exterior \(\xi\) does **not**
   carry particle mass. Then the slogan “the displacement … is the same object
   that carries mass” is false for the exterior, or “mass” has been redefined
   as a topological occupation number rather than a \(\xi\)-integral. The design
   does not state that redefinition.

3. \(\sum\theta_i=0\) makes densification a global zero-sum. The rarefaction that
   pays for a dense core is distributed over the fabric. If that rarefaction is
   part of the particle’s mass accounting, mass is nonlocal. If it is not, mass
   is only the core count while the compensating \(\theta>0\) lives “elsewhere”
   without mass — an asymmetry that needs a rule.

**Where does the mass of the distortion live?**

The design does **not** answer. The coherent options are:

| Option | Mass of distortion | Cost |
|--------|--------------------|------|
| A | Zero (pure kinematic mismatch, no \(\xi\)-mass) | Contradicts “\(\xi\) carries mass” uniformly |
| B | Included in particle mass via some \(M[\xi|_E]\) | Mass is not pure cell count; continuous pieces appear |
| C | Part of a separate radiation/field sector | Introduces a second mass-bearing account; risks §6.2 |
| D | Topological: only cells inside the mismatch surface count | Requires a sharp, invariant definition of that surface and a proof that exterior \(\xi\) is massless |

None of A–D is selected. **POSTULATE needed:** a mass functional \(M[\xi,\text{cycles}]\)
and a partition into particle vs environment.

**Conclusion (A.3).**

“Mass is a cell count” is **not consistent as stated** with both (i) \(\xi\) as the
mass-carrying DOF and (ii) a necessary exterior distortion. The mass of the
distortion is **unlocated**. This is an internal coherence failure of §1+§5,
not a wording quibble: §7’s “mass ratios fixed by configuration” depends on
which answer is chosen.

---

### A.4 The 3-from-7 claim (§3)

**Mathematical fact (external, used only as geometry, not as physics authority).**

A bilinear cross product \(\mathbb{R}^n\times\mathbb{R}^n\to\mathbb{R}^n\) with the usual
norm and orthogonality properties exists only for \(n=3\) and \(n=7\). That is a
theorem about real composition algebras / Hurwitz.

**What the design claims.**

Handedness requires a cross product \(\Rightarrow\) fabric has 7 degrees of freedom
from which 3 dimensions emerge — structural, not tunable.

**What the design’s own §2 actually uses.**

A link along \(\hat d\) has a **transverse plane** and two amplitudes \(u_1,u_2\):

\[
u_\pm=\frac{u_1\mp i u_2}{\sqrt{2}},\qquad
\chi=2\,\mathrm{Im}(u_1^*u_2).
\]

That construction needs:

- a distinguished longitudinal direction \(\hat d\) on the link, and
- an **oriented 2-plane** orthogonal to \(\hat d\).

In ordinary 3-dimensional space this is immediate: the orthogonal complement of a
line is 2-dimensional, and \(\chi\) is the signed area form on that plane. **No
7-dimensional cross product appears in the formulae.**

**Does “handedness requires a cross product” force 7 DOF and 3 emergent dimensions?**

**No.** Breakdown:

1. **Handedness of a link \(\neq\) ambient cross product algebra.** Circular
   polarisation on a 2-plane is enough. That is available as soon as space is
   3-dimensional (or, more minimally, as soon as each link has a 2D transverse
   fibre). It does not require \(\mathbb{R}^7\).

2. **Degrees of freedom \(\neq\) spatial dimensions.** “7 DOF per cell” and “3
   emergent spatial dimensions” are different kinds of count. A theory can have
   many internal components on a 3D complex without “dimension emerging from
   DOF.” The design never defines the map
   \(\{\text{7 DOF}\}\to\{\text{3 spatial directions}\}\).

3. **No mechanism for emergence.** For 3-from-7 to be load-bearing one would need
   something like: a 7-component object with a \(G_2\) (or related) structure;
   a dynamical preference that selects a 3-plane as “propagation space”; a proof
   that chiral links live only in that 3-plane. None is stated.

4. **7 is not used later.** Sections 4–7 never invoke seven components, \(G_2\),
   octonions, or a reduction \(7\to 3\). Species labels are mode indices,
   resonance vectors, handedness, cell count — not 7D representation labels.

**Classification.**

| Role | Assessment |
|------|------------|
| Load-bearing for transfer/handedness as written in §2 | **No.** §2 only needs a 2D transverse plane. |
| Decorative | **Yes** as currently written: a true geometric fact placed next to an unsupported architectural claim. |
| Unsupported implication | **Yes.** Cross product in 3 and 7 \(\not\Rightarrow\) this fabric has 7 DOF and 3 emergent dimensions. |

**Conclusion (A.4).** The 3-from-7 paragraph is **unsupported** as an implication
and **not load-bearing** for the rest of the design. Removing it does not damage
§2’s chiral link story. Keeping it requires a mechanism that is absent.

---

## B. Does it produce particles?

§7 lists five required products. For each: plausibility, derivability, and what
would demonstrate it. Status tags: **consequence**, **postulate**, **hope**.

### B.1 Multi-cell localised objects with sharp exterior

| | |
|--|--|
| **Status** | **Postulate** of what a particle *is* (§5), not a derived existence theorem. |
| **Plausible?** | Conditionally. Spectral mismatch + cycle gating *could* reflect energy at a boundary if the dispersion of interior and exterior cycle structures really do not share complete cycles. That is a standard idea (impedance mismatch, stop-band, trapped mode) once dynamics exist. |
| **Derivable from the text?** | **No.** There is no equation of motion, no spectrum of free-cell modes, no proof that a multi-cell locked arrangement is a stable invariant set. Tightening and saturation (§4) suggest a path to dense cores but do not construct a closed interior. |
| **Demonstration** | Full 3D nonlinear evolution (or a stationary variational problem on the full fabric) exhibiting a region \(R\) with: (i) \(\lvert R\rvert\ge 2\); (ii) interior cycle content that does not complete transfer cycles with exterior probes; (iii) a measurable jump in density/rate at \(\partial R\); (iv) lifetime \(\gg \tau\) against radiation into the bulk. Residual: leaked power through the boundary / interior cyclic energy \(\ll 1\). |

### B.2 Discrete spectrum of species (integer labels, recurring)

| | |
|--|--|
| **Status** | **Hope**, motivated by compact internal space (§1) and “existence is a threshold.” |
| **Plausible?** | Integer mode indices from a compact internal configuration space are **DERIVED** *if* that internal space is truly compact and the free motion is Hamiltonian with a discrete spectrum — that part of §1 is standard. That discrete labels *classify stable multi-cell particles* is not. |
| **Derivable?** | **No.** Threshold existence can still allow continuous moduli (positions, relative phases, soft shape modes). “No continuous family of particles” is **POSTULATED**. Recurrence across independent preparations is an empirical claim. |
| **Demonstration** | Many independent preparations (random or annealed initial data) on a large 3D fabric; cluster the resulting long-lived objects by integer invariants (cell count, mode content, net \(\chi\), resonance vector). Species are discrete only if clusters are isolated in invariant space with gaps that do not fill when the preparation measure is refined. Fit residual: within-cluster variance of continuous observables vs between-cluster gaps. |

### B.3 Binding (attractive channel vanishing with separation)

| | |
|--|--|
| **Status** | **Hope / requirement.** Not a consequence. See §C. |
| **Plausible?** | Unknown from the text. Handedness supplies channel *diversity*, not attraction. |
| **Derivable?** | **No** (detailed in §C). |
| **Demonstration** | Two multi-cell objects, controlled separation \(r\) and relative handedness; measure force or \(\Delta E(r)=E(r)-E(\infty)\). Binding requires some channel with \(\Delta E(r)<0\) for a range of \(r\) and \(\Delta E(r)\to 0\) as \(r\to\infty\). Report \(\Delta E(r)\) with numerical noise floor. |

### B.4 Mass ratios fixed by configuration

| | |
|--|--|
| **Status** | **Hope**, and **blocked** by the mass-accounting hole in §A.3. |
| **Plausible?** | Only after \(M\) is defined. If \(M=\lvert R\rvert\), ratios are rational by construction but may be trivial (all “masses” are small integers without hierarchy). If \(M\) includes exterior distortion energy, ratios need not be fixed by interior integers alone. |
| **Derivable?** | **No.** |
| **Demonstration** | Operational mass (inertia under a slow boost of the whole object, or total cyclic energy in \(R\) under a stated functional) for each species; ratios stable under box size and preparation; residuals of any claimed integer formula. |

### B.5 Derived \(c\) and genuinely cycle-gated transfer

| | |
|--|--|
| **Status** | **Partial consequence** for a speed limit of order \(a/\tau\); **postulate** for complete-cycle gating (§A.1–A.2). |
| **Plausible?** | Speed limit: yes. Hard gating: only with extra structure. |
| **Derivable?** | Limit \(v\lesssim a/\tau\): yes. Complete-cycle-only delivery: no. |
| **Demonstration** | (i) Inject a localised energy packet; measure signal arrival vs distance; fit \(v_{\mathrm{front}}\) and compare to \(a/\tau_{\mathrm{cycle}}\) (report both). (ii) For gating: drive a link off resonance and mid-cycle; measure delivered energy to the far cell as a function of drive duration and detuning. Complete-cycle claim requires delivered energy consistent with zero for drive windows that are not full cycles, within a stated noise floor — not merely a sinusoidal modulation. |

### B.6 Scoreboard

| §7 item | Consequence of design? | Best current status |
|---------|------------------------|---------------------|
| Multi-cell sharp objects | No — definitional postulate | Plausible target |
| Discrete species | Integer *labels* plausible; discrete *spectrum of particles* not derived | Hope |
| Binding | **No** | **Missing** |
| Mass ratios from configuration | No — mass undefined under distortion | Hope, blocked |
| Derived \(c\) + cycle gate | Speed *limit* yes; gate no | Partial |

---

## C. The hard problem: binding

### C.1 What is required

§7: an interaction between objects that

1. vanishes with separation, and
2. is **attractive in some channel**,

so composites can form.

§2: same-handed and opposite-handed configurations couple through **different
channels**. \(\chi=2\mathrm{Im}(u_1^*u_2)\) is a signed geometric quantity on links.

### C.2 What handedness actually gives

**DERIVED from §2.**

- Links carry a signed scalar \(\chi\) (or two circular amplitudes).
- Processes that depend on relative phase of transverse actuations can treat
  \(\chi>0\) and \(\chi<0\) differently.
- Therefore at least two coupling channels exist (schematically same-\(\chi\) vs
  opposite-\(\chi\), or more if relative orientation of two particles is included).

**NOT DERIVED.**

- The sign of the force (attractive vs repulsive) in either channel.
- The radial dependence \(F(r)\) or \(\Delta E(r)\).
- That either channel is attractive at intermediate range and \(\to 0\) at infinity.

Handedness is a **degeneracy-splitting / selection** structure. It is not an
attractive potential. Analogy only for clarity of logic (not an appeal to
another theory’s authority): polarisation-dependent scattering can repel in one
channel and attract in another; the existence of two channels never proves
attraction exists in either.

### C.3 Candidate mechanisms inside the design (examined)

#### (i) Overlap of exterior distortions

Two particles’ exterior density/rate fields overlap at finite separation.
Interaction energy is then some functional of the joint \(\xi\).

- **POSTULATE required:** the energy functional \(H[\xi,\text{cycles}]\).
- Without \(H\), the sign of \(\Delta E(r)\) is undefined.
- Even with a typical elastic energy \(\sim\lvert\nabla\xi\rvert^2\), same-sign density
  bumps often **repel** (positive definite quadratic forms give repulsion of
  like deformations). Attraction would need a specific nonlinear or
  sign-indefinite term the design does not supply.

**Verdict:** possible in principle; **not derived**; default quadratic guess
leans repulsive for like cores.

#### (ii) Global scarcity \(\sum_i\theta_i=0\)

On a closed fabric total dilation vanishes. Dense cores must be paid for by
rarefaction elsewhere.

**DERIVED:** densification is zero-sum on a closed manifold.

**NOT DERIVED:** local binding.

Argument: write \(\theta_i=\theta_i^{\mathrm{loc}}+\bar\theta\) with
\(\sum\theta^{\mathrm{loc}}=0\) enforced by a uniform shift \(\bar\theta\sim -Q/N\) for total
“charge” \(Q=\sum\theta^{\mathrm{loc}}\) of densified regions. On a fabric of \(N\) cells,
the constraint-mediated interaction is **\(O(1/N)\)** and **range-infinite**
(mean field). As \(N\to\infty\) at fixed particle content, this interaction
**vanishes**. It is not a local attractive channel that can form composites
whose binding is a bulk property independent of box size.

If the energy density \(f(\theta)\) is local and convex, isolated lumps cost bulk
self-energy; the constraint does not create a short-range well between two lumps.

**Verdict:** scarcity prevents a uniformly densified vacuum; it does **not**
provide §7 binding on a large fabric.

#### (iii) Cycle-mediated forces between exteriors

If complete cycles can run in the exterior and between two particles’ near
fields, energy exchange or virtual cycle processes might lower energy at finite
\(r\).

- Same-handed vs opposite-handed would modulate the rate of those processes
  (**channel diversity** — DERIVED as possibility).
- Lower energy at finite \(r\) (**attraction**) requires that the joint
  configuration with both objects present has less total cyclic + elastic energy
  than two isolated objects — **not shown**.

Spectral mismatch is what *defines* a particle (interior does not complete
cycles with exterior). That same idea makes the exterior a poor coupler into
the core; it does not by itself create an attractive *inter-particle* channel.
If anything, hard gating suppresses the very exchanges that might generate
forces — unless the force lives entirely in off-shell / incomplete virtual
processes, which the design’s “partial cycle delivers nothing” language
**forbids** if taken literally. That is a direct tension:

> Literal complete-cycle gating \(\Rightarrow\) no mid-cycle virtual transfer  
> \(\Rightarrow\) the usual field-theoretic sources of potential energy are
> restricted  
> \(\Rightarrow\) binding becomes **harder**, not easier, unless a separate
> mechanism (static elastic energy of \(\xi\), cycle-averaged stress) is
> introduced.

**Verdict:** under a literal gate, dynamical exchange forces are suppressed;
static distortion energy remains, with no derived attractive sign.

#### (iv) “Release into configuration” after tightening saturation (§4)

Saturation dumps energy into arrangement of internal harmonics and neighbours.
That could create multi-cell locked patterns (relevant to **existence** of
composites as single locked objects) but is not the same as a two-body potential
between well-separated particles that vanishes with separation. Fusion into one
object by violent overlap is not §7 binding.

### C.4 Can an attractive channel be derived?

**Short answer: no, not from the design as stated.**

What can be said carefully:

| Claim | Tag |
|-------|-----|
| Two handedness channels exist | **DERIVED** from §2 |
| Some channel is attractive | **Not derived**; no candidate with controlled sign |
| Interaction vanishes as \(r\to\infty\) | **Plausible** if energy is local (elastic + finite-range cycle coupling); **not proven** without \(H\) |
| Composites form | **Hope** |

**If the design cannot bind.**

As written, the design:

1. defines particles as spectrally isolated regions,
2. supplies chiral channel *labels*,
3. supplies a global densification constraint that does not bind locally,
4. does **not** supply an energy functional or a process that lowers energy when
   two isolable objects approach in some channel.

Therefore **binding is not a consequence of the design; it is an unpaid debt of
§7.** Stating that plainly is the correct first-principles outcome. Completing
the design requires either:

- an explicit local energy \(H\) whose two-body sector has a negative well in at
  least one handedness/orientation channel, with \(H_{\mathrm{int}}\to 0\) as
  \(r\to\infty\), or
- a derived effective interaction from the completed cycle dynamics that is
  measured to be attractive,

and a check that this mechanism is not a new field smuggled in violation of §6.2.

**GUESS (clearly labelled, not used as evidence):** the most natural place to
*look* for attraction after the dynamics are written is nonlinear elastic energy
of \(\xi\) with a sign-indefinite or saturation-induced term, or a
cycle-averaged stress from exterior link motion — **not** \(\sum\theta=0\) and
**not** \(\chi\) alone. That is a research direction, not a derivation.

### C.5 Consistency note: gate vs force

A design that both (a) forbids partial-cycle energy transfer and (b) needs soft
inter-particle potentials faces a structural choice:

- Soften the gate (allow cycle-averaged or virtual transfer) \(\Rightarrow\)
  “complete cycles only” is false as stated; or
- Keep a hard gate \(\Rightarrow\) inter-particle forces must come from **static**
  geometry (\(\xi\), locked configurations), not from dynamical energy exchange.

The text currently asserts (a) strongly and needs (b) for §7. That is an
internal design pressure, not yet a contradiction if static \(\xi\)-energy can
attract — but static attraction is **not written down**.

---

## D. Implementation

### D.1 Was a true experiment run?

**No.**

Constraints from the user and from the design:

- Full 3D, full nonlinear, no reduced dimension, no linear proxy, no substitute
  mechanism.
- Particles and forces must emerge from fabric DOF only (§6).

**Blocking reason (not budget alone):** the design document **underdetermines the
dynamics**. Missing pieces required before any admissible simulation:

1. State space per cell and per link (exact list of real degrees of freedom).
2. Equations of motion or a discrete variational update (energy \(H\) and
   symplectic/dissipation structure).
3. Precise meaning of “complete cycle delivery” as a term in those equations.
4. Internal compact configuration space (manifold, metric, free spectrum).
5. Coupling of cyclic energy to cell size (tightening law + saturation).
6. Operational definitions of particle mass, boundary, and force.

Any code written without fixing these would **add** physics. Reporting its
output as a test of *this* design would violate the “no substitute mechanisms”
and “true experiment” rules. I therefore did not build or run a simulator.

**Code created for this review: none.**

### D.2 Specification precise enough to implement

The following is a **completion sketch**: every item marked **OPEN** is a
postulate the design does not fix; an implementer must choose and document it.
Until OPEN items are fixed, results are tests of the *completion*, not of the
bare text.

#### State

- Lattice: cubic index set \(\Lambda\subset\mathbb{Z}^3\), spacing parameter \(a>0\).
  **OPEN:** topology (periodic box recommended so \(\sum\theta=0\) is meaningful).
- Per cell \(i\in\Lambda\):
  - displacement \(\xi_i\in\mathbb{R}^3\) (or scalar radial model only if proved
    equivalent — default **vector** \(\xi\));
  - internal configuration \(q_i\) on a compact manifold \(Q\). **OPEN:** \(Q\)
    (minimal chiral choice: \(Q=S^1\) per relevant mode is too small for
    multi-mode species; product of circles \(T^m\) or \(SU(2)\) etc. must be
    stated);
  - cyclic energy bookkeeping: either derived from \(\dot q_i\) or an explicit
    action variable \(J_i\). **OPEN.**
- Per oriented link \(\ell=\langle i,j\rangle\):
  - transverse amplitudes \((u_1,u_2)_\ell\in\mathbb{R}^2\) (or one complex \(u_\ell\));
  - phase / cycle coordinate \(\phi_\ell\in S^1\);
  - **OPEN:** conjugate momentum or second-order oscillator variables.
- Density: \(\theta_i=(\mathrm{Div}\,\xi)_i\) with a stated discrete divergence.
  Enforce \(\sum_i\theta_i=0\) identically if \(\xi\) is pure displacement on a
  closed complex (discrete Helmholtz projection). **OPEN:** whether \(\xi\)
  includes global translations (should not).

#### Geometry and handedness

- Link direction \(\hat d_\ell\).
- Transverse frame \((e_1,e_2)\) with orientation; \(\chi_\ell=2u_1 u_2\) after
  suitable phase convention, or \(\chi=2\mathrm{Im}(u_1^*u_2)\) in complex form.
- **Do not** implement 7D cross products unless a separate 7-component state is
  defined; §2 does not require them.

#### Dynamics (must be full nonlinear 3D)

A minimal admissible shape (illustrative, **all coefficients OPEN**):

1. **Internal / link oscillators** with natural frequencies from the free
   spectrum on \(Q\) and on the link fibre.
2. **Cell–link coupling** that realises the gate. To match the text literally,
   use a **latch**: net energy transfer across \(\ell\) updates only when
   \(\phi_\ell\) crosses a complete-cycle mark, by an amount determined by
   resonance of the two cells with the link mode; mid-cycle coupling to far-cell
   energy is zero. (Alternative monodromy form: transfer \(\propto\) increment of
   a cycle-averaged action, with instantaneous power a total derivative.)  
   **Document which.** Linear continuous coupling **without** a latch is a
   **different design** and must not be advertised as testing complete-cycle
   gating.
3. **Tightening:** effective cell volume or rest length decreases with cyclic
   energy \(J_i\) until saturation \(J_\star\); beyond \(J_\star\), excess goes into
   rearrangement of \(q_i\) and neighbours (configuration release), not further
   compression. **OPEN:** functions \(V_{\mathrm{cell}}(J)\), saturation rule.
4. **Elastic / geometric energy** of \(\xi\) so that exterior distortions cost
   energy and can mediate static forces. **OPEN:** \(H_\xi[\xi]\) — this is where
   binding must be shown or fail.
5. Integrator: structure-preserving if \(H\) is given; fixed \(\Delta t\) with
   documented conservation residuals (energy drift, \(\sum\theta\), discrete
   Gauss-like identities if any).

#### Measurements (admissible only in 3D nonlinear runs)

| ID | Measurement | Falsifies if |
|----|-------------|--------------|
| M1 | Front speed vs \(a/\tau\) | No finite speed limit of that order |
| M2 | Transfer vs drive duration (sub-cycle vs full cycle) | Sub-cycle delivery comparable to full-cycle (hard gate fails) |
| M3 | Search for multi-cell long-lived cores | None exist above noise for all tried completions |
| M4 | \(\Delta E(r)\) for pairs in same- and opposite-\(\chi\) channels | \(\Delta E(r)\ge 0\) all \(r\), all channels (binding fails) |
| M5 | Species clustering over preparations | Continuous smear of invariants (discrete spectrum fails) |
| M6 | Mass operational definition vs \(\lvert R\rvert\) and exterior energy | Exterior dominates or mass not integer-stable (mass claim fails) |

#### What must not be done

- 1D chains or 2D slabs (handedness transverse plane removed).
- Linearised acoustics as evidence of binding or species.
- Hand-placed “particle fields,” Mexican-hat stabilisers standing in for §5.1
  interior configuration, or imported species.
- Reporting a different gate than the one claimed.

#### Cost note

A honest first binding test (M3+M4) needs: enough cells for two multi-cell
objects plus separation scan; evolution time \(\gg\tau\); several handedness
orientations. That is a full nonlinear 3D code with OPEN choices fixed in a
written completion note. It was not executed here.

---

## E. What would falsify the design

### E.1 Cheapest *true* experiment that can kill it

**Target:** §7 binding (and secondarily multi-cell existence).

**Protocol (minimum admissible):**

1. Fix a completion of §D.2 (state, \(H\), gate rule) in a short companion note;
   freeze parameters; do not retune after seeing forces.
2. On a periodic 3D lattice large enough that two objects and a separation scan
   fit without forced overlap of exteriors at the maximum \(r\) (document \(N\) and
   minimum image distances).
3. Obtain at least one multi-cell long-lived object (by annealing, collision, or
   threshold construction **using only fabric DOF**). If none exists after a
   pre-registered search budget, the design completion **fails existence** — stop.
4. Place two copies at separation \(r\) along a lattice direction; for each
   same-handed and opposite-handed relative channel (and at least two relative
   orientations if objects are anisotropic); measure
   \(\Delta E(r)=E_{\mathrm{pair}}(r)-2E_1\) after relaxation of fast modes.
5. **Kill criterion (binding):** for all channels and orientations tested,
   \(\Delta E(r)\ge 0\) for all \(r\) down to contact, within the energy noise floor
   \(\varepsilon\), **or** \(\Delta E(r)\) does not tend to \(0\) as \(r\to\infty\)
   (non-local junk from constraints/box). Then §7 binding fails for that
   completion.
6. **Kill criterion (gate):** if M2 shows sub-cycle delivery of energy across
   links at a level comparable to full-cycle delivery, “energy crosses only in
   complete cycles” is false for that implementation; if the only way to restore
   the gate is a new field, §6.2 fails.

This is not a proxy: it is full 3D, nonlinear, and directly checks the design’s
own success list. It is the cheapest *decisive* true experiment because a single
negative binding survey (after existence is secured) ends the §7 claim without
needing mass-ratio spectroscopy.

### E.2 Cheaper *coherence* kills (no simulation)

These do not replace §E.1 but already constrain the text:

| Issue | Effect |
|-------|--------|
| Complete-cycle gate not following from link mediation | §2 slogan is incomplete until extra structure is added |
| Mass of exterior distortion unlocated | §5 mass claim incoherent until \(M\) is defined |
| 3-from-7 non sequitur | §3 drops out without damaging §2 kinematics |

A design can survive by **repair** (add gate law; define \(M\); drop or derive
3-from-7). It cannot survive a clean M4 failure after a fair completion.

### E.3 What would *not* count as falsification

- Failure of a 1D/2D model.
- Failure of a linearised model to bind.
- Failure of a substitute mechanism (e.g. ordinary Klein–Gordon Q-balls) that is
  not this design.
- Aesthetic dislike of discreteness or of complete-cycle language.

---

## F. Ledger: derive / postulate / guess

### Derived from the text’s premises

1. Energy transfer is delayed and mediated by link degrees of freedom.
2. A speed *limit* of order \(a/\tau\) exists if each hop needs time \(\ge\tau\).
3. A 2D transverse plane on a link supports a signed handedness \(\chi\).
4. Same- vs opposite-handed processes can couple differently (channel diversity).
5. Compact internal configuration space \(\Rightarrow\) discrete free mode indices
   (for the free internal problem).
6. Closed fabric + \(\theta=\mathrm{Div}\,\xi\) \(\Rightarrow\) \(\sum\theta_i=0\)
   (discrete Stokes), hence densification is zero-sum.
7. Global zero-sum constraint \(\Rightarrow\) only \(O(1/N)\) mean-field interaction
   from the constraint alone on large \(N\).

### Postulated (stated or required, not forced by prior lines)

1. Partial cycle delivers exactly nothing (hard gate).
2. Gate is structural rather than perturbative resonance.
3. Scalar isotropic \(c\) as *the* propagation constant.
4. Particle \(=\) spectrally mismatched multi-cell region with sharp exterior.
5. Interior stabilising arrangement exists and is never a field.
6. Mass \(=\lvert R\rvert\) cell count.
7. Existence is a hard threshold with no continuous particle families.
8. Attractive binding channel exists (§7).
9. Mass ratios fixed by configuration integers.
10. 7 DOF with 3 emergent dimensions from cross-product classification.

### Guesses (explicitly not evidence)

1. Static nonlinear \(\xi\)-energy is the least-implausible place to look for
   attraction once dynamics exist.
2. Literal hard gating will suppress exchange forces and push all binding burden
   onto static geometry.
3. Without a carefully chosen \(H_\xi\), like-signed dense cores will repel.

---

## G. Direct answers to the five brief questions

**A. Coherence.**

- Complete-cycle delivery: **does not** follow from link mediation alone; needs
  a latch/monodromy/quantisation rule.
- \(c\): a speed *limit* \(\sim a/\tau\) **does** follow; isotropic Lorentzian \(c\)
  does **not**; the microscopic rate remains an input.
- Mass as cell count vs exterior distortion: **inconsistent as written**;
  distortion mass is **unlocated**.
- 3-from-7: **unsupported** and **not load-bearing** for §2; decorative as
  currently stated.

**B. Particles.** Multi-cell trapping is a definition/postulate; discrete species
and mass ratios are hopes; binding is missing; speed limit is partial;
hard cycle gate is postulated. Nothing in §7 is demonstrated.

**C. Binding.** Handedness \(\Rightarrow\) two channels, **not** attraction. Global
\(\sum\theta=0\) \(\Rightarrow\) no local binding of useful strength on large
fabrics. **As stated, the design does not bind.** That is the central negative
result of this review.

**D. Implementation.** Not run. No code created. Specification in §D.2; OPEN
items must be fixed before any admissible experiment.

**E. Falsification.** Cheapest true kill: full 3D nonlinear pair survey of
\(\Delta E(r)\) in both handedness channels after multi-cell objects exist; all
channels non-attractive (or gate unrealisable without a new field) kills §7 / §6.

---

## H. What would make the design coherent next (not done here)

Repairs that stay inside first principles, each necessary:

1. **Write the gate** as a term in the update rule (latch or monodromy), not as
   a slogan.
2. **Write \(H\)** (or EOM) for \(\xi\), internal \(q\), and link \((u,\phi)\).
3. **Define mass** \(M\) and prove or drop “cell count only,” including the
   exterior.
4. **Either derive or delete** 3-from-7.
5. **Exhibit or abandon** an attractive channel: analytic effective
   \(\Delta E(r)\) from \(H\), then confirm with the §E.1 experiment.

Until (1)–(5) are done, §7 remains a wish list and binding remains the
unpaid hard problem.
