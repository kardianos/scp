# v68 GAUGE DESIGN — Promoting the diagonal U(1) to a local gauge symmetry

**Date**: 2026-06-10. **Status**: design document (no code). This gates the next
kernel-modification authorization request (CLAUDE.md kernel policy).
**Problem being fixed** (v67 FINDINGS #1, FUTURE.md "Theta-Boundedness", route 3):
the complexified theory's Θ is massless + U(1)-charged + free — nature forbids this
combination; the measured η-drain bleeds every ball, and th1 (m_θ=1.6 > ω) closes the
drain 400× but kills the long-range mediator. Gauging restores BOTH: stable matter and
a long-range force.

**Verification**: `v68/theory/gauge_checks.mac` — **59/59 PASS**
(`maxima -b gauge_checks.mac`, output `gauge_checks.out`). Checks G1–G15 are cited
inline. Claim markers: [thm] exact, [verified] machine-checked, [estimate], [open].

Framework: v66/THEORY.md (12-field complexified Cosserat, Q-ball branch),
v67/theta_dynamics/THETA_DYNAMICS.md + DEBROGLIE.md. Standard parameters m²=2.25,
m_θ²=0 (to be revisited, §1.4), η=0.5, μ=−41.345, κ=50. Ball used for all numbers:
ω=1.39, Q=482, E₀=692, R_half=4.4, window ω∈(1.3087, 1.5).

---

## 0. Executive summary / decision

- **Field content (decision)**: add a NEW real abelian gauge field A_μ (option a);
  keep Θ as charged matter. Option (b) "Θ IS the connection" is structurally
  impossible as a reinterpretation — the η coupling is LINEAR in Φ, minimal coupling
  is QUADRATIC; no local field redefinition connects them [thm, G15]. Option (c)
  (composite phase-connection of Θ) is singular at Θ=0 and fixes nothing. The honest
  realization of the user's instinct is (a) + decoupling: A takes over Θ's mediator
  role; Θ becomes a massive Yukawa-dressed charged field (th1, measured).
- **Recommended package**: gauge the full diagonal U(1) (Φ AND Θ charged, the η
  coupling covariantized) + m_θ > ω. Gauging alone does NOT close the Θ-drain;
  the mass does (th1: 400×). Gauging supplies what the mass route destroyed:
  a massless NEUTRAL mediator → true 1/r Coulomb force, with zero charge drain
  [thm, G14].
- **Coupling strength**: g = 0.05 default (g=0.02 conservative). g=0.1 pushes the
  ω=1.39 ball's Coulomb self-shift δω ≈ +0.13 across most of the existence window
  and caps Q_max ≈ 780 [estimate, §5].
- **Lattice**: compact link angles θ_i = g·a·A_i + noncompact E_i (Kogut–Susskind
  Hamiltonian form), leapfrog-compatible; Gauss law conserved to machine rounding by
  construction (no projection step), seeded exactly at t=0 [standard + G12].
  Memory +6 N³ f32 arrays (+25%).
- **What reverses**: AB phase null (THETA_DYNAMICS §3.3) → 2πqn holonomy restored;
  "no chiral splitting" theorem void (fluctuation operator now contains A_bg);
  tb-style two-ball phenomenology gains a 1/r channel — same-charge pairs REPEL
  beyond r* ≈ 12 (g=0.05) [estimate].

---

## 1. Field content decision

### 1.1 Option (a) — new A_μ, Θ retained as charged matter (RECOMMENDED)

Add real A_μ (temporal gauge: 3 spatial components + 3 conjugate E_i; 2 physical
polarizations after the Gauss constraint). Minimal addition; the entire v66/v67
edifice (charge conservation proof, Q-ball branch, drain measurements, bath theory)
remains the g→0 limit, term-by-term [verified G13].

- Pro: smallest delta from the working kernel; exact validation ladder (g=0 must be
  bit-identical to v66); A is automatically U(1)-NEUTRAL [thm, G14] — the abelian
  connection carries no charge, so the new radiation channel cannot drain Q.
- Con: Θ is still massless+charged+free at m_θ=0 — gauging does not close the
  measured η-drain (the radiated Θ waves still carry ρ_Q^θ = ω|h|² > 0). The package
  therefore includes m_θ > ω (§1.4). With A present, the old objection to the mass
  route ("kills the long-range mediator", v64 H1) evaporates.

### 1.2 Option (b) — Θ reinterpreted as the connection (user's original instinct)

**DOF accounting**: Θ has 3 complex spatial components = 6 real fields (+6 momenta).
A_i in temporal gauge needs 3 real fields (+3 momenta), of which 2 are physical
after Gauss. A literal map (e.g. tu_i → A_i) strands the tv triplet (+ all momenta
mismatch): 4 of Θ's 6+6 phase-space pairs are leftover. Options for the remainder —
discard (then this is just option (a) with Θ deleted), or keep tv as a neutral
Proca triplet (needs a mass to be healthy; new unmotivated sector).

**Does η Θ·(∇×Φ) emerge as a limit of minimal coupling? NO** [thm, G15]:

- Minimal coupling cross-term: |D_iΦ|² = |∂Φ|² + 2g A·J^grad + g²A²|Φ|², with
  J^grad_i = Σ_a(u_a∂_iv_a − v_a∂_iu_a) = Im(Φ̄·∂_iΦ) [verified G3]. This vertex is
  QUADRATIC in matter (A·J_Q structure).
- The η coupling η Re[Θ̄·(∇×Φ)] is LINEAR in Φ — a two-point Proca-like mixing, not
  a connection vertex. Machine check G15a/b: ∂²/∂(matter)² of the η term vanishes
  identically; the same derivative of g A·J^grad is the nonzero constant gA. A local
  invertible field redefinition preserves the lowest matter power of a vertex, so no
  redefinition maps one into the other. This is the algebraic root of the
  THETA_DYNAMICS §3.3 [thm] "AB phase exactly zero" null.
- **Dual relation analysis**: integrating by parts, η Θ̄·(∇×Φ) = η (∇×Θ̄)·Φ — Θ
  enters only through ∇×Θ, i.e. Θ acts like a vector potential whose MAGNETIC field
  couples linearly (Zeeman/dipole-type) to Φ. The dual map A ~ curl⁻¹Θ is nonlocal
  AND still yields a matter-linear vertex — the obstruction is the matter power, not
  the curl. Θ is "B-field-like", never "connection-like". [thm-level]
- **What survives of the instinct** [estimate/open]: in a charged condensate
  background Φ_bg, minimal coupling linearizes to a two-point A–δΦ mixing
  g A_i·Im(Φ̄_bg ∂_i δΦ)-type — the same one-derivative two-point structure as the η
  coupling, with η_eff ~ g|Φ_bg|. The old η-phenomenology (polariton dressing,
  preferred frame) can re-emerge as IN-MEDIUM physics of the gauged theory. The η
  coupling is the "already-condensed" caricature of a gauge coupling — which is
  presumably why it felt connection-like.

Verdict: (b) as a reinterpretation is dead; (b) as an endpoint ("A is THE connection,
Θ decoupled") is reachable from (a) by config (η→0 or m_θ large).

### 1.3 Option (c) — composite gauge field from Θ's own phase

a^θ_i = Σ(tu∂_i tv − tv∂_i tu)/|Θ|² (the U(1) connection of Θ's phase; cf.
THETA_DYNAMICS §3.3: ∮a_θ·dl = 2πnB², quasi-quantized). Rejected: singular wherever
|Θ|→0 (generic — Θ radiates and passes through zero), no natural Maxwell term
(a Maxwell term for a^θ is a horrendous nonlocal Θ self-interaction), does not fix
Θ-boundedness (Θ remains the dynamical field), and the holonomy is amplitude-weighted
(2πnB²), not topological. Useful only as a diagnostic observable.

### 1.4 Recommended package (decision)

Gauge the FULL diagonal U(1) — both Φ and Θ minimally coupled with the same charge
g, η coupling covariantized (it must be: at η≠0 the term is gauge-invariant ONLY if
Θ co-rotates with Φ [verified G5]; a neutral-Θ variant requires deleting the η term,
i.e. option-(b)-endpoint). Θ-boundedness is closed by m_θ = 1.6 > ω (th1, measured
400×). Stage-0 config (η=0, Θ frozen) gives clean scalar-electrodynamics Q-ball
physics for validation.

---

## 2. The gauged Lagrangian, EOM, Gauss law

Conventions [all verified]: D_μ = ∂_μ + igA_μ on charge-(+1) fields; gauge
transformation Φ→e^{−igα}Φ, Θ→e^{−igα}Θ, A_μ→A_μ+∂_μα; covariant curl
(D×X)_k = ε_kij D_i X_j; E_i = ∂_i A_0 − ∂_t A_i. These signs keep the v66 charge
density ρ_Q = Σ_a(u_a v̇_a − v_a u̇_a + …) unchanged and give ∇·E = +gρ_Q.

    L =  Σ_a [ ½|D_tΦ_a|² − ½|D_iΦ_a|² − (m²/2)|Φ_a|² ]
       + Σ_a [ ½|D_tΘ_a|² − ½|D_iΘ_a|² − (m_θ²/2)|Θ_a|² ]
       − Vt(s)  +  η Re[ Θ̄ · (D×Φ) ]  +  ½|E|² − ½|∇×A|²

Gauge invariance [verified]: covariance of D_i (G1, 18 components), |DΦ|² invariant
(G2), covariant curl covariance (G4), gauged η term invariant (G5). Vt(s) depends
only on |Φ_a|² — trivially invariant; the WP-A "no bare mass from the potential"
theorem survives gauging, and the honest O(A²) seagull mass g²A²|Φ|² appears exactly
as predicted in THETA_DYNAMICS §6 [verified G3].

**EOM (temporal gauge A_0 = 0)** [thm; building blocks verified G7a–f]:

    Φ̈_a = D_iD_iΦ_a − m²Φ_a − 2Vt′(s)Φ_a Π_{b≠a}|Φ_b|² + η (D×Θ)_a
    Θ̈_a = D_iD_iΘ_a − m_θ²Θ_a + η (D×Φ)_a
    Ä_i  = −[∇×(∇×A)]_i + g J_i                       (E_i = −Ȧ_i)

with the single current (its three independent derivations agree):

    J = −Σ_a Im[ Φ̄_a D Φ_a + Θ̄_a D Θ_a ] + η Im[ Θ̄ × Φ ]
      = −J^grad(Φ) − J^grad(Θ) − gA(|Φ|²+|Θ|²) + η( u×tv + tu×v )

In components D_iD_iΦ reads: ∇²u − g(∇·A)v − 2gA·∇v − g²A²u (+ i↔ swap with signs)
[verified G7a,b]; the gauged η force on Φ is the covariant curl of Θ with the SAME
D [verified G7c–f].

**Gauss-law constraint** [verified G6, from the A_0 variation before gauge fixing]:

    ∇·E = g ρ_Q ,   ρ_Q = Σ_a Im[Φ̄_a D_tΦ_a + Θ̄_a D_tΘ_a]  (= v66 ρ_Q at A_0=0)

**Charge conservation is now LOCAL** [verified G8 + G8a]: substituting the full
gauged EOM, ∂_tρ_Q + ∇·J = 0 holds POINTWISE (potential pairing cancels pointwise by
U(1) invariance, G8a; gradient, gauge, seagull and η terms cancel exactly, G8). The
v66 result was global-up-to-flux; gauging upgrades it to a local continuity equation
whose current is precisely the v66 Noether current covariantized — the Ampère source
∂L/∂A_i = gJ_i is the SAME J [verified G9], and div(∇×∇×A) = 0 [G10] means the Gauss
constraint, once true, propagates under the EOM [thm].

**Static ball + Coulomb solves the constraint** [verified G12]: for Φ_a = f(r)e^{iωt},
ρ_Q = 3ωf², the radial field E_r(r) = (3gω/r²)∫₀^r f²s²ds satisfies ∇·E = gρ_Q
exactly (G12a); the ball carries no spatial current (J^grad = 0, G12b); Coulomb E is
curl-free so B stays zero (G12c). Asymptotically E_r → gQ/(4πr²): the ball wears its
charge as a measurable 1/r² field — **Q becomes readable from a boundary flux
integral alone** [thm, Gauss].

Temporal-gauge subtlety [thm]: a permanently static E forces A_L(t) = −E·t (pure
gauge, longitudinal); the exact stationary gauged ball instead carries the Coulomb
potential in its local phase rate, ω_local(r) = ω + g·a₀^equiv(r). On the lattice the
compact formulation absorbs this growth harmlessly (§4); seeding the ungauged profile
+ Gauss E gives an O(g²) settling transient [estimate].

---

## 3. What gauging fixes — item by item

**(i) Amplitude unboundedness of the connection = pure gauge** [thm, G11]. A → A+∇χ
shifts amplitudes arbitrarily at zero field strength (F[∇χ]=0, G11); the longitudinal
mode is eaten by the Gauss constraint, leaving 2 bounded transverse oscillators.
Contrast Θ: its amplitude is physical energy. CAVEAT [measured]: this fixes the
CONNECTION's boundedness; Θ as charged matter still needs m_θ > ω (th1). Under the
recommended package both statements hold simultaneously.

**(ii) Massless NEUTRAL mediator — long-range force WITHOUT charge drain** [thm,
G14]. Every A-dependent term in L (A·J^grad, A²|Φ|², the η seagull A·(u×tv+tu×v)) is
invariant under the global matter rotation with A held FIXED [verified G14a–c] ⟹ the
Noether charge density contains NO A or E contribution; a radiated A wave carries
energy and momentum but identically ZERO U(1) charge. The drain channel of v66 §5
(Θ radiation with ρ_Q^θ = ω|h|² > 0) has no analog in the A sector: **the drain
closes identically for the gauge mediator**. Moreover at η=0 a static ball does not
radiate A at all — its current is static (G12b) — so the gauged ball is exactly
non-radiating, while still sourcing a permanent 1/r Coulomb field. This is the
configuration th1 could not provide.

**(iii) AB phase / winding holonomy restored** [thm structure]. With minimal coupling
the matter phase acquires e^{−ig∮A·dl} around any loop: a flux tube of n quanta gives
holonomy 2πqn (q = g in our units). The THETA_DYNAMICS WINDING null ("encircling
phase = 0 exactly, at any order") was a theorem ABOUT the η-curl coupling (matter-
linear, G15); gauging changes the coupling class and the null REVERSES — winding
becomes genuinely topological for Φ. On compact links the flux through a plaquette is
quantized mod 2π/(ga), making the sectors exact lattice integers. The DEBROGLIE
fingerprint sharpens: interference of a moving ball around a flux line measures q per
charge quantum (k = p/Q chain, DEBROGLIE §2.4). Similarly the "no chiral splitting at
O(ηB)" theorem is void: the linear fluctuation operator now contains A_bg directly
(THETA_DYNAMICS §6 anticipated exactly this).

**(iv) True 1/r two-ball force — same-charge REPEL** [thm structure + estimate
numbers]. Two balls of charges g·Q₁, g·Q₂ interact by Coulomb F_C = g²Q₁Q₂/(4πD²)
(repulsive for like charges — sign fixed by the standard positive-definite ½|E|²
energy of the abelian field between like sources), ON TOP of the measured
phase-coherent contact force ~cos(Δφ)e^{−μ_t D}, μ_t = √(m²−ω²) = 0.564 (tb1 merger /
tb2 ejection taxonomy). Both forces now coexist:

- Coulomb (g=0.05, Q=482): F_C = 46.2/D². Contact-force band (DEBROGLIE §5 tail
  calibration, geometric mid): F_t ≈ 5.1, 1.3, 0.35, 0.10 at D = 8, 10, 12, 14.
- **Crossover radius r\*** (co-phase, same charge: attraction inside, repulsion
  outside — an UNSTABLE watershed: D < r\* → merger, D > r\* → escape):
  r\* ≈ 12 at g=0.05 (band 9.5–14.5 from the decade-wide tail calibration);
  r\* ≈ 9 at g=0.1 (7–11.5); r\* ≈ 16 at g=0.02. [estimate]
- Anti-phase same-charge: both channels repulsive — unconditional escape (tb2
  strengthened). Opposite-charge pair (Φ and Φ̄ balls): the relative phase beats at
  2ω, the contact force time-averages ≈ 0 [estimate], leaving a CLEAN attractive
  Coulomb orbit problem — see the positronium prediction, §5.

---

## 4. Lattice implementation

**Compact vs noncompact — recommend COMPACT links** (Kogut–Susskind Hamiltonian
U(1)): store link ANGLES θ_i(x) = g·a·A_i(x) wrapped to (−π,π] plus noncompact
conjugate E_i(x).

- Why compact: (1) the temporal-gauge secular growth A_L ∝ t around any static charge
  (§2) wraps harmlessly instead of losing f32 precision in a large-A representation
  (at g=0.05, Q=482, r=8: E≈0.03, A grows to ~9 by T=300 — gaA ≈ 0.06 rad and
  growing without bound; compact stays O(1) forever); (2) flux quantization gives exact
  integer winding sectors (the (iii) observable); (3) it is the standard scheme whose
  leapfrog conserves the discrete Gauss law to machine rounding.
- Cost of compact: forces involve sin(plaquette sums) — a few transcendentals per
  link per step (~10–15% step cost [estimate]); flux quantum 2π/(ga) is large at
  small g (IR artifact only for deliberately wound sectors — fine).
- Noncompact fallback: simpler stencils, identical continuum limit; acceptable if
  paired with periodic re-gauging (project ∇·A→0 every ~10³ steps via one FFT
  Poisson solve + matter phase rotation). More moving parts; not recommended.

**Update scheme** (leapfrog/Verlet, interleaved with the existing kernel):

    1. matter kick:  Φ̇,Θ̇ += dt·F[Φ,Θ,U]   — all spatial differences link-covariant:
                      (D_iΦ)(x) ≈ [U_i(x)Φ(x+î) − Φ(x)]/a,  U_i = e^{iθ_i}
                      (η covariant curl from the same transported differences)
    2. E kick:       E_i += dt·[ (1/(g a³))·Σ_staples sin(θ_plaq) + g·J_i^lat ]
    3. drift:        Φ += dt·Φ̇ ;  θ_i −= dt·g·a·E_i  (wrap to (−π,π])

**Gauss law maintenance**: NO projection/damping step. If (a) the lattice Hamiltonian
is exactly invariant under time-independent lattice gauge transformations (it is,
when ALL matter derivatives — including the η covariant curl — use link transport)
and (b) the seed satisfies the discrete constraint, then the symplectic update
conserves G(x) = (1/a)Σ_i[E_i(x)−E_i(x−î)] − gρ_Q(x) to rounding error [standard
real-time lattice-gauge result; the continuum analog is G9+G10]. Report
gauss_max(t), gauss_l2(t) in diag.tsv as integrator-health diagnostics (the analog of
the v66 Q-drift floor).

**Seeding (init=qball, complex_gauge=1)**: temporal gauge ⟹ seed E_i, not A_0.
Single ball (absorbing BC): exact radial Gauss integral E_r(r) = (3gω/r²)∫f²s²ds
[verified G12a], θ_i = 0. Multi-ball absorbing: superpose radial solutions (Gauss is
linear). Periodic BC: **a torus admits NO net charge** (∮E·dA = 0 ⟹ Σ gρ_Q = 0
exactly) — periodic runs must be net-neutral (±Q pairs) or use explicit jellium
subtraction (solve FFT Poisson for ρ_Q − ρ̄; document the fictitious background).
This is a real experimental-design constraint the v67 bath-ladder style (periodic,
one charged ball) inherits — those configs must move to absorbing BC or neutral
pairs when gauged.

**Memory**: +3 link angles +3 E_i = +6 N³ f32 arrays = +25% over the 24 matter arrays
(N=256: +0.40 GB; N=384: +1.36 GB → ~6.8 GB total, still fits V100-16GB).

**CFL/stability**: the A sector is a c=1 wave field — same dt ≤ a/√3 bound as Φ;
the seagull g²A²|Φ|² and current terms are small at g ≤ 0.1 (no new constraint)
[estimate]. f32 risk: per-step link increment g·a·E·dt ~ 10⁻⁵ rad accumulating over
10⁶ steps — marginal; mitigate with f64 (or Kahan-compensated) gauge-sector
accumulators [open].

---

## 5. Predictions table

Ball: ω=1.39, Q=482, E₀=692, R≈4.4 unless noted. "V#" = validation-ladder run (§6).

| # | observable | run | prediction | marker |
|---|---|---|---|---|
| P1 | g=0 mode | V1 | bit-identical to v66 kernel (forces reduce exactly) | [thm G13] |
| P2 | Gauss residual | V2 | flat at rounding floor for all T (no secular drift) | [standard] |
| P3 | Coulomb halo | V2 | E_r(r) = gQ_enc(r)/4πr²; at r=10: 0.019·(g/0.05) | [thm G12] |
| P4 | ball clock shift | V2 | δω ≈ +13.1 g² = +0.033 at g=0.05 (toward window top; factor-2 syst.) | [estimate] |
| P5 | charge drain at η=0 | V2 | ZERO beyond integrator floor; no A radiation (static current) | [thm] |
| P6 | A-wave dispersion | V3 | ω=k exactly, 2 transverse pols; no η mixing with A (A couples to J, not curl) | [thm] |
| P7 | full package η=0.5, m_θ=1.6, g=0.05 | V6 | drain stays closed (≤0.0025 Q/t.u., th1 level) AND Coulomb halo P3 present | [measured+thm] |
| P8 | same-charge co-phase pair | V4 | watershed r\*≈12: D=8 → merger (tb1-like), D=16 → Coulomb escape (tb1 REVERSED) | [estimate] |
| P9 | opposite-charge pair | V5b | clean Coulomb capture/orbit; "positronium" n=1 closure orbit D₁ = 8πM/(ω²g²Q²) ≈ 15.5 at g=0.05, v_rel≈0.093, T_orb≈1050 | [estimate, DEBROGLIE D13 + Coulomb] |
| P10 | AB holonomy | V5 | interference shift 2πn around an n-quantum flux line; WINDING null reversed | [thm structure] |
| P11 | maximum charge | scan | Coulomb self-energy caps the branch: Q_max ≈ (0.846/g²)^{3/2} ≈ 6200 (g=0.05), 780 (g=0.1); thin-wall balls Coulomb-forbidden | [estimate] |
| P12 | long-range drift | V4b | isolated pair at D=20, g=0.05: F=0.116, a=1.7×10⁻⁴ → Δx≈0.85 in T=100 (v51-scale, easily measurable) | [estimate] |

Falsifiers: P2 failing ⟹ the lattice η-curl term is not gauge-invariantly
discretized (implementation bug class #1). P4 with wrong SIGN ⟹ charge-convention
error. P8 with no watershed ⟹ tail-force calibration off by >decade — recalibrate
against tb data before trusting P9.

**Coupling-choice table** [estimate]:

| g | δω self | r\* | D₁ (positronium) | Q_max | verdict |
|---|---|---|---|---|---|
| 0.02 | +0.005 | ≈16 | 97 (outside box) | ~9.7×10⁴ | cleanest window, weakest signals |
| **0.05** | +0.033 | ≈12 | 15.5 (in box!) | ~6200 | **default** |
| 0.1 | +0.13 | ≈9 | 3.9 (< ball size) | ~780 | overbound; window nearly consumed |

---

## 6. Kernel-v2 architecture sketch (file-level)

Guarded mode `complex_gauge=1` (requires `complex_phi=1`; refuse with any of the
v66-refused options). New config keys: `g_gauge=0.05`, `gauge_compact=1`,
`gauss_check_every=100`. **All changes to `sfa/sim/scp_sim.cu`/`.c` require the
explicit user authorization this document is requesting.**

- `sfa/sim/scp_sim.cu` / `.c`: +6 arrays (θ_link[3], E[3]); link-covariant
  differences in the Φ/Θ force kernels (one fused pass; reuse halo logic); E-kick
  kernel (staples + matter current); link drift + wrap; Gauss-residual reduction.
  η covariant curl built from the SAME transported differences (gauge invariance of
  the discrete Hamiltonian is what makes Gauss exact — do not mix raw and
  transported stencils).
- `sfa/seed/gen_qball*` (and `gen_qball_pair`, `gen_qball_bath`): emit E_i seed
  (radial Gauss integral; superposition for pairs; FFT-Poisson + neutrality check
  for periodic). Refuse periodic + net charge ≠ 0.
- diag.tsv additions: `gauss_max`, `gauss_l2`, `E_em` (½E²+B² plaquette),
  `Q_flux_r8` (Gauss-readout of interior charge — independent cross-check of
  Q_core), optional Wilson-loop probe.
- SFA output: optional A/E columns behind a flag (viewer unaffected by default).
- **Pre-kernel theory task** (no authorization needed): gauged radial system —
  f″ + (2/r)f′ = (m² − ω̃(r)²)f + 2Vt′(f⁶)f⁵ with ω̃ = ω + g a₀(r),
  a₀″ + (2/r)a₀′ = −3gω̃f² — two-field shooter to recompute the existence window,
  Q(ω;g), Q_max(g), and the exact seed profiles (replaces the P4/P11 estimates).

**Validation ladder**: V1 g=0 bit-compare vs v66 → V2 single ball g=0.05 η=0
(P2–P5) → V3 vacuum A-wave dispersion (P6) → V4 same-charge pair watershed (P8,P12)
→ V5 flux-line AB + V5b opposite-charge orbit (P9,P10) → V6 full package η=0.5,
m_θ=1.6 (P7).

---

## 7. Open risks

1. **Existence-window erosion** [estimate→open]: Coulomb self-repulsion shifts ω up
   (P4) and caps Q (P11); at g=0.05 the usable window shrinks by ~20%. The gauged
   radial shooter (§6) must precede 3D runs; if the window closes faster than
   estimated, drop to g=0.02.
2. **Absorbing BC vs Gauss law** [open]: sponge-damping E/matter violates the
   constraint at the boundary. Mitigation: damp matter and the TRANSVERSE gauge
   sector only, monitor gauss residual confined to the sponge shell, keep all
   physics readouts interior (already standard practice). Periodic runs must be
   net-neutral (torus Gauss theorem, §4).
3. **η covariant curl discretization** [open]: ordering ambiguities O(a²); a
   non-gauge-invariant stencil silently breaks exact Gauss conservation (P2 is the
   tripwire). Budget a dedicated unit test (random gauge transformation on the
   lattice → all forces invariant).
4. **f32 precision in the gauge sector** [open]: tiny per-step link increments;
   consider f64 gauge arrays (+3% total memory at N=384) if gauss_l2 drifts.
5. **Θ at m_θ=0 stays sick** [measured]: gauging alone does NOT close the η-drain
   (th2-class behavior persists). The package is gauge + m_θ>ω; any m_θ=0 gauged
   run must be labeled as drain-afflicted.
6. **Parameter count grows** [open, theory-level]: g is a new genuine input
   (alongside a_lepton, α of the v59 audit). The v59 bridge numbers are untouched
   (the gauge sector adds no term to Vt or the mass bilinears), but a future
   derivation of g — e.g. from the holonomy quantum or a Θ-condensate matching
   η_eff ~ g|Φ_bg| (§1.2) — is the natural follow-up question.
7. **Phenomenology re-derivation debt** [open]: WP-A/B/C results (polariton speeds,
   Mathieu tongues, bath ladder predictions) are statements about the ungauged
   theory; each inherits O(g²) corrections and the chirality/AB no-gos flip. The
   v67 measured baselines remain valid as the g=0 control points.

## Files

- `v68/theory/gauge_checks.mac` — Maxima verification, **59/59 PASS**
  (G1–G15: covariance, invariance, EOM, Gauss law, local continuity, Ampère-current
  consistency, Coulomb-ball constraint, g→0 limit, A-neutrality, no-redefinition
  lemma). Output: `gauge_checks.out`.
- Inputs: `v66/THEORY.md`, `v67/FINDINGS.md`, `v67/theta_dynamics/THETA_DYNAMICS.md`,
  `v67/theta_dynamics/DEBROGLIE.md`, `FUTURE.md` (boundedness routes).
