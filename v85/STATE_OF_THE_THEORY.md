# State of the theory — 2026-07-26
**Scope:** the standing synthesis after the v85 campaign (ANALYSIS → TOE
v85.3 → SLIP/TOPO/G2/SURV/QRK → X10 series → X11-PB). Four sections, per
request: what is understood, what is open-and-alive, ways forward, and
problems with the theory and numerics as currently expressed. Null/negative
results are cited only where they bound a claim; the registry of closures
lives in the individual results files.

---

## 1. The theory as currently understood

**Ontology.** One continuum, fundamentally in motion; particles, forces and
light are its modes. The vacuum reference is a motion pattern, not a stage.
[Program axiom; everything below is its measured content.]

**Matter.** A particle is a closed process (uptake = layment per cycle).
The unique stable sector is the gauged U(1) phase sector: Q-balls on a
continuous branch Q ∈ (86.6, 924) at g=0.05, with exact discrete Gauss law
(1e-13 floor) and lattice SR kinematics (E² = p² + E₀²). Persistence is
survivorship: what cannot close, dies — and survivorship discretizes **mode
indices at fixed charges** (winding, shell level), not charges or energies
(H5-weak).

**Energy and mass.** dE/dQ = ω exactly (Legendre, 8e-5); E = Qω(1+ε) with
ε = +0.9…+4.3% measured everywhere — Qω is the thin-wall skeleton, ∫T₀₀ the
inertial mass (D5′). For **bound linear modes E = ωQ is exact** (canonical
identity) — the load-bearing exception. Composite bookkeeping: ΔE = ∫ω dQ.
E = mc² is the effective SR identity, not an axiom.

**Substructure.** The nucleon-analog is a **multi-harmonic bound state, not
a composite of sub-particles** (QRK): three components locked by the
phase-blind product potential, confined by construction (no isolated
components, no mesons), carrying separately-conserved global flavor charges
(three U(1)s over one gauged total — the charge/lepton-number split,
code-verified). Spatially the harmonics are co-located (displaced patterns
collapse in ~50 t.u.); the substructure is real but lives in **mode space**,
measurable as split response lines (0.108/0.126 in the flavored ball) while
the clock partition itself is invisible to density probes (phase-blindness).

**Forces.** Gauge Coulomb (state-blind, long-range); phase-coherent contact
force cos(Δφ)e^{−κD} with the three-regime detuning law (coherent /
algebraic 2/(Δωτ) averaging / Adler locking); κ-knee saturation as the hard
core; capacity fold Q_max ≈ (2.1–2.6)/g² terminating the family; no
interior iron peak (proven). The **fabric is a weak superconductor**: London
screening measured ∝ g² to 2%, capped below 1% of the flux quantum in this
potential family. A second compact fiber (gauge2) exists and works: exact
integer-charge structure (color-Zeeman ∝ g₂², exact sign flips), conserved
color charge (C = −27.04, 0.3%/120 t.u.), imaged 1/r² color-Coulomb tail.

**Atomic physics.** Around a charged core, the same fabric supports discrete
gap-protected response modes (KG-Coulomb; α_f = g²Q_N/4π with a built-in
ceiling α_max ≈ 0.17–0.21 — no supercritical shells anywhere). Measured
across five GPU runs: shells form at the predicted (hard-core-corrected)
radii, are retained indefinitely when BC-clean and contact-clean, cascade
inward under radiation reaction, and breathe monopole-protected (a radial
charge oscillation cannot radiate). Smeared clouds cannot orbit — **orbits
require carriers** — and a genuine two-closure Coulomb orbit was achieved
(Kepler to 6%), bounded by carrier fragility: the only same-branch light
carrier is evaporation-metastable and runs away under sustained agitation.
Annihilation between conjugate closures is throttled ~50× by detuning +
flavor mismatch, surviving even full contact as a slow burn.

**Quantum boundary (mapped, not guessed).** Light quantization ⇔ charge
quantization (via exact linear-mode E = ωQ): universal ħ and quantized
charge are one fact in this fabric — but no classical mechanism found
supplies the quantum: slips quantize events not amounts; erosion selects
nothing; all classical integers are topological (winding, fiber ratios,
slip counts). Bell is closed by theorem (commutative algebra). The honest
ceiling: **structural analog up to the commutative wall**, with the wall's
location now measured from inside.

---

## 2. Open items that are alive (not null, not negative)

1. **The fiber (color) sector** — every measurement to date positive:
   integer charge structure exact, field imaged, benign to activate. Its
   dynamical force between objects remains unmeasured cleanly (probe-noise,
   not null). The CUDA kernel lacks gauge2 entirely (port = open item).
2. **Cloud-atoms (Stage 4, honest form)** — retention + shell localization
   are PASSED physics; a carbon-scale nucleus + self-limiting cloud at
   neutrality is designed (g=0.023, Q_N=650, D12-clean) and unrun.
3. **2B multi-center nuclei** — now load-bearing (three independent threads:
   carrier count, interference universality, constituent clocks). Untested
   at nucleus scale; the fiber can't hold spacing at g₂=0.05 but Coulomb +
   saturation (the v74 droplet mechanism) remains available.
4. **Second-sector carrier (TOE §D.2)** — the structural remedy for
   HIER/COUNT: a light-gap sector sharing the compact gauge U(1) (charge
   universality by construction if commensurate harmonics are used).
   Design-only; nothing yet contradicts it.
5. **Topological quantization in stiff regimes** — the g² screening law is
   confirmed physics; the option-C target (Σ|Φ|²a² ≥ 50, deficit only
   ×5–7) is concrete and untested.
6. **X10c/X11 late-time objects** — the breathing retained cloud and the
   throttled binary are new, characterized-but-unexplained states; their
   long-time fate (T ≫ 2000) is open.
7. **Mode spectroscopy as the fabric's "particle physics"** — QRK-1's split
   lines opened a working spectroscopy; the response-spectrum catalog
   (needed to define X1 alignment) is started, not finished.
8. **Cosmogonic selection (Stage 5)** — untested and *upgraded* by SURV:
   if charge selection exists it operates in agitated formation epochs;
   condensate-fragmentation runs would measure it directly.

## 3. Possible ways forward

- **(A) The carbon endgame on current physics** (recommended default):
  2B-style multi-nucleon assembly (droplet mechanism) + the D12-clean cloud
  at neutrality → a structural carbon atom: Z=6-mapped core, −Q_N cloud at
  shell radii. Two to three GPU runs. Accepts no-orbit, cloud-form
  electrons as the fabric's honest atom.
- **(B) Second-sector kernel design** — attack HIER directly: light gap
  m_l ≪ m sharing the gauge link (and optionally the fiber for commensurate
  charges). Kernel surgery with a real design doc first; the payoff is a
  stable light carrier and true orbital atoms.
- **(C) Fiber science** — port gauge2 to CUDA, measure the color force
  cleanly (scattering geometry), probe confinement and the g₂-scaling of
  pattern stabilization (QRK-2 at large g₂).
- **(D) Option-C stiff potential** — a vev-type modulus for real flux
  quantization; changes the whole branch structure; highest theory risk,
  highest topological payoff.
- **(E) Stage-5 fragmentation** — condensate collapse censuses: does
  formation select charges/objects? Directly tests cosmogonic Q-PIN and
  produces the abundance-inventory physics of the standing goal.
- **(F) Consolidation (v86)** — the "what the fabric is and is not" scoping
  document + CLAUDE.md stage-table update (2B mandatory; Stage 3 reframed);
  cheap, overdue, and orthogonal to all of the above.

## 4. Problems and issues, as currently expressed

**Theory-level.**
- **The Q-degeneracy stands.** Charge ≡ action ≡ (nearly) mass in one
  branch; HIER (carrier ≈ 0.8 nucleon masses vs 1/1836) and COUNT (no q₀
  divides Z-carbon into 6) are its two faces. Every atom result inherits it.
- **No quantum of exchange (T2/JUMP).** All transfer amounts are continuous;
  interrupted transitions legally trap fractional quanta; finite-carrier
  self-energy chirps lines ~29%. The quantum-light chain is conditional on
  carriers whose integrity mechanism does not exist in-kernel.
- **q₀ has no selector.** Erosion-floor failed dynamically; candidate scales
  (fold 90 / flux 126 / threshold 157) never coincide; "cosmogonic
  selection" is currently a hope with a test design, not a mechanism.
- **The atom's electron is unresolved in kind**: cloud (retained, breathing,
  no orbit) vs carrier (orbits, but dies). Nature's electron is both stable
  AND orbital; the fabric currently offers either/or.
- **E = Qω's ε-term** is booked but unexplained — no closed form for ε(Q);
  the "harmonic-native" axioms survive only as asymptotics.
- **Monopole-protected breathing** means relaxation to stationary states is
  generically obstructed for radial excitations — a structural annoyance
  for every future bound-state experiment (excitations linger ≫ runs).

**Numerics-level.**
- **Async SFA tail corruption**: final snapshot frames unreadable in two of
  five GPU runs (x11pb t≥1800). Data loss was masked by diag cadence; the
  writer needs a fix or a final-flush barrier (kernel I/O change, small).
- **Seed-shape mismatch**: ungauged shooter profiles relax visibly at high
  Q under g=0.05 (ball "swelling", X10b) — gauged-shooter seed profiles
  are needed for any high-Q precision claim.
- **Box/sponge systematics are load-bearing**: BC-limited decay (pathfinder),
  sponge-eaten evanescent tails, and the D7 4-axis certification (dt, dx,
  L, seed precision) has never actually been run (X2 remains undone —
  arguably the most overdue chore in the program).
- **Observable gaps**: per-component Q_a still absent from diag (SFA-frame
  analysis only); θ-probes useless at η=0; no in-run cluster charges
  (the c12/F6 lump census still BLOCKED on a 42 GB archive pull).
- **Renderer floor**: volview's transfer function hides low-amplitude
  objects (×0.15 satellite invisible) — fine for balls, misleading for
  perturbation-scale physics; trajectory tools (cpos/rprof/ccent, in
  scratchpad only) should be promoted into the repo's analysis toolset.
- **Resolution ceilings unaudited**: N=64 for all CPU campaigns and dx=0.287
  at N=384 for GPU runs; no refinement pass has confirmed that any of the
  new v85 phenomenology (slosh periods, burn rates, split lines) is
  dx-converged.
