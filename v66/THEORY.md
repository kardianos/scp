# v66 THEORY — Complexified 6-field Cosserat (Q-ball sector)

**Date**: 2026-06-09
**Status of each claim is marked**: [thm] = exact derivation, [verified] = checked
symbolically (Maxima, `theory/qball.mac`, 19/19 PASS) and/or numerically
(`theory/qball_checks.py`), [estimate] = numerical with known approximation,
[open] = unresolved.

**Motivation** (v65, `v65/FINDINGS_STABILITY.md`): the real 6-field theory provably has
no stable localized structure — static lumps are Derrick saddles, and the cubic source
P = φ₀φ₁φ₂ radiates at the third harmonic. The unique minimal fix is a Q-ball: complexify
the fields, make the potential U(1)-invariant; the conserved Noether charge Q removes the
Derrick instability and a stationary rotating ansatz Φ ∝ e^{iωt} makes |Φ| (and hence the
potential source) time-independent, killing the harmonic-radiation channel in the Φ sector.

---

## 1. Fields, Lagrangian, EOM

Complexify each of the six real fields:

    Φ_a = u_a + i v_a        (a = 0,1,2;  "phi" sector)
    Θ_a = tu_a + i tv_a      ("theta" sector)

12 real fields total. **Real limit**: v_a = tv_a = 0 reproduces the existing 6-field
system EXACTLY (term-by-term, verified below).

Lagrangian density [thm]:

    L =   Σ_a [ ½(u̇_a² + v̇_a²) − ½(|∇u_a|² + |∇v_a|²) − (m²/2)(u_a² + v_a²) ]
        + Σ_a [ ½(tu̇_a² + tv̇_a²) − ½(|∇tu_a|² + |∇tv_a|²) − (m_θ²/2)(tu_a² + tv_a²) ]
        − Vt(s)
        + η [ tu·(∇×u) + tv·(∇×v) ]                  (= η Re[ Θ̄ · (∇×Φ) ])

with

    s     = Π_a (u_a² + v_a²) = Π_a |Φ_a|²
    Vt(s) = (μ/2) s / (1 + κ s),     Vt′(s) = (μ/2) / (1 + κ s)²      [verified C2]

At v = 0:  s = (u₀u₁u₂)² = P², and Vt(P²) = (μ/2)P²/(1+κP²) = V(P) identically
[verified C1].

**Equations of motion (all 12 fields)** [thm; building blocks verified C3,C5,C6]:

    ü_a  = ∇²u_a  − m² u_a  − 2 Vt′(s) u_a Π_{b≠a}(u_b² + v_b²) + η (∇×tu)_a
    v̈_a  = ∇²v_a  − m² v_a  − 2 Vt′(s) v_a Π_{b≠a}(u_b² + v_b²) + η (∇×tv)_a
    tü_a = ∇²tu_a − m_θ² tu_a + η (∇×u)_a
    tv̈_a = ∇²tv_a − m_θ² tv_a + η (∇×v)_a

The curl coupling complexifies pairwise — (u,tu) and (v,tv) — exactly as the real system
couples (φ,θ); the v-sector is a second copy. Standard parameters: m² = 2.25 (m = 1.5),
m_θ² = 0, η = 0.5, μ = −41.345, κ = 50.

**Potential force** [verified C3a,b,c]:

    −∂Vt/∂u_a = −2 Vt′(s) u_a Π_{b≠a}(u_b² + v_b²)

and at v=0 this equals the existing real force −∂V/∂φ_a = −V′(P)·P/φ_a with
V′(P) = μP/(1+κP²)² (Maxima check C3b: exact symbolic equality). So a complex run with
all v, tv initialized to zero is bit-for-bit the same physics as the real kernel
(same forces, same kinetic normalization).

### 1.1 Optional terms — how they WOULD complexify (NOT in v1)

These are documentation only. **The v1 implementation REFUSES `complex_phi=1` combined
with any of: `alpha_cs≠0`, `beta_h≠0`, `kappa_h≠0`, `mode≠0`, `sigma_grad/sigma_cubic/
sigma_freq/sigma_cross≠0`, `lambda_self≠0`, `self_tune≠0`** (config error at startup).

- Cosserat strain `alpha_cs`: α|∇×φ/2 − θ|²  →  α( |∇×u/2 − tu|² + |∇×v/2 − tv|² )
  = α |∇×Φ/2 − Θ|² (complex modulus). Diagonal-U(1) invariant.
- Hardening `beta_h`: (β/2)|θ|²|∇×φ|²  →  (β/2)(Σ_a tu_a²+tv_a²)(|∇×u|² + |∇×v|²).
  Diagonal-U(1) invariant.
- Chiral helicity `kappa_h`: κ_h P² φ·(∇×φ)  →  κ_h s Re[Φ̄·(∇×Φ)]. Diagonal-U(1)
  invariant, but mixes u/v sectors; excluded from v1.
- `mode≠0` (inverse / density-κ mass couplings), the σ_* theta self-interactions,
  `lambda_self`, and `self_tune` (dynamical κ) each modify the potential sector in ways
  whose complexification has NOT been derived; refusing them keeps the v1 charge
  conservation proof airtight.

## 2. Symmetry and Noether charge

The potential Vt(s) is U(1)³-invariant (per-a phases), but the curl coupling mixes
a-components spatially [(∇×Φ)_x involves Φ_y, Φ_z]. The exact continuous symmetry is the
single **diagonal global U(1)** [thm; generator invariance verified C4]:

    Φ_a → e^{iα} Φ_a,   Θ_a → e^{iα} Θ_a   (same α for all a)
    i.e.  δu_a = −α v_a, δv_a = +α u_a, δtu_a = −α tv_a, δtv_a = +α tu_a

**Noether charge (sign convention: Q > 0 for Φ ∝ e^{+iωt}, ω > 0)** [thm]:

    Q = ∫ d³x  Σ_a ( u_a v̇_a − v_a u̇_a + tu_a tv̇_a − tv_a tu̇_a )

**Conservation** [verified C4, C5, C5b]: from the EOM,

    ∂_t ρ_Q = ∇·J_Q^(grad) − η ∇·( u×tv + tu×v ),
    J_Q^(grad) = Σ_a ( u_a∇v_a − v_a∇u_a + tu_a∇tv_a − tv_a∇tu_a )

- mass terms cancel pointwise (u_a v_a − v_a u_a = 0);
- potential terms cancel pointwise by U(1) invariance of s (Maxima C4);
- Laplacian terms are the total divergence ∇·(u∇v − v∇u) (Maxima C5b);
- curl cross-terms are EXACTLY the total divergence −∇·(u×tv + tu×v) (Maxima C5,
  via div(A×B) = B·∇×A − A·∇×B); no antisymmetry is lost at finite η.

So the conserved current is J_Q = −J_Q^(grad) + η(u×tv + tu×v) and dQ/dt = 0 in the
bulk for ANY field configuration, including η ≠ 0 and Θ ≠ 0.
**Boundary caveat** [thm]: with absorbing BC, dQ/dt = −(boundary flux); charge carried
to the sponge by radiation (Φ or Θ waves) is deliberately removed. The lattice
diagnostic should report Q(t) over the interior and treat the deficit as measured
radiated charge, not as a violation. (A discrete-exact conservation law on the lattice
stencil is not claimed here — [open], to be measured at η=0 as the integrator's Q-drift
floor.)

## 3. Symmetric Q-ball ansatz and radial ODE

Ansatz: Φ_a = f(r) e^{iωt} for all a, Θ = 0. Then s = f⁶ exactly [verified C6a] and the
potential force on every component factorizes with the SAME radial coefficient
[verified C6b,c], so all six (u,v) equations reduce to ONE radial ODE [thm]:

    f″ + (2/r) f′ = (m² − ω²) f + 2 Vt′(f⁶) f⁵,        f′(0)=0, f(∞)=0

Equivalently f″ + (2/r)f′ = U_ω′(f) with the effective potential [verified C7]

    U_ω(f) = ½(m² − ω²) f² + (1/3) Vt(f⁶)
           = ½(m² − ω²) f² + (μ/6) f⁶/(1 + κ f⁶)

(the 1/3 because the full potential is shared by three identical complex fields).
NOTE: Θ = 0 is NOT an exact solution when η ≠ 0 — see §5. The ansatz is exact for η = 0.

### 3.1 Existence window [thm + verified C8]

A localized solution requires (i) ω² < m² (exponential tail e^{−√(m²−ω²) r}) and
(ii) U_ω(f) < 0 for some f > 0, i.e.

    ω² > ω_min² = m² + min_f [ (2/3) Vt(f⁶)/f² ] = m² + min_f [ (μ/3) f⁴/(1+κf⁶) ]

Since μ < 0, the minimum is at the maximum of g(f) = f⁴/(1+κf⁶):
g′ ∝ f³(4 − 2κf⁶) = 0  ⟹  **κ f⁶ = 2**, f* = (2/κ)^{1/6}, g(f*) = (1/3)(2/κ)^{2/3},
giving the closed form

    ω_min² = m² + (μ/9) (2/κ)^{2/3}

Standard parameters (m²=2.25, μ=−41.345, κ=50):

    f* = 0.584804,  ω_min² = 1.712695,  ω_min = 1.308700
    **existence window: ω ∈ (1.3087, 1.5000)**   (width 0.191)

[verified: Maxima C8 closed form + Python brute-force minimization agree to 1e-6;
this confirms the prior estimate κf⁶=2, f*≈0.585, ω_min²≈1.712, ω_min≈1.3087.]

### 3.2 Charge, energy, stability criteria

On the ansatz [verified C10, C11]:

    ρ_Q = 3ω f²        ⟹  Q(ω) = 3ω ∫ f² dV
    e   = (3/2)(ω²f² + f′² + m²f²) + Vt(f⁶)   ⟹  E(ω) = ∫ e dV

Exact relation on the solution family: dE/dω = ω dQ/dω [thm, standard Q-ball identity;
verified numerically to ≤0.6% across the window, §6].

- **Classical (fission) stability**: dQ/dω < 0 [thm, standard criterion
  (Friedberg–Lee–Sirlin / Vakhitov–Kolokolov type); applied here, verified numerically].
- **Evaporation bound**: E < m·Q. Normalization check [verified C9]: a circularly
  polarized free plane wave u = A cos(ωt−kx), v = A sin(ωt−kx) with ω² = k²+m² has
  energy density ω²A² and charge density ωA², so E/Q = ω ≥ m for ANY free-quanta
  configuration. With our ½u̇² kinetic convention both E and Q scale as A², so the
  bound is E < mQ = 1.5·Q with no extra factors.

## 4. Numerical solutions (shooting, `theory/qball_checks.py`)

Profiles found by long-double RK4 shooting on f₀ ∈ (f₁, f_dip) of U_ω (under/over
bisection, 110 iterations); exponential tail attached at truncation; ODE residual
≤ 1e-4 (relative) except the thin-wall edge row. Table [estimate — rough shooting;
the dedicated C solver replaces this]:

| ω     | f(0)    | s₀=f⁶   | e(0)  | Q       | E       | E/Q    | E−mQ    | r_half | r_Q(rms) |
|-------|---------|---------|-------|---------|---------|--------|---------|--------|----------|
| 1.315*| 0.5893  | 0.0419  | 1.793 | 114174  | 151305  | 1.3252 | −19956  | 19.5   | 19.6     |
| 1.330 | 0.6004  | 0.0468  | 1.883 | 22207   | 29786   | 1.3413 | −3524   | 16.0   | 12.1     |
| 1.360 | 0.6238  | 0.0589  | 2.084 | 1745.0  | 2422.9  | 1.3884 | −194.7  | 6.86   | 5.31     |
| 1.390 | 0.6405  | 0.0690  | 2.253 | 482.2   | 691.9   | 1.4349 | −31.4   | 4.41   | 3.71     |
| 1.420 | 0.6366  | 0.0665  | 2.275 | 209.6   | 309.7   | 1.4778 | −4.66   | 3.29   | 3.14     |
| 1.440 | —       | —       | —     | —       | —       | —      | ~0      | —      | —        |
| 1.450 | 0.6121  | 0.0526  | 2.147 | 118.0   | 178.6   | 1.5128 | +1.51   | 2.67   | 3.09     |
| 1.470 | 0.5837  | 0.0395  | 1.979 | 91.9    | 140.5   | 1.5281 | +2.59   | 2.40   | 3.43     |
| 1.490 | 0.5326  | 0.0228  | 1.681 | 92.3    | 141.2   | 1.5287 | +2.65   | 2.26   | 5.03     |

(*thin-wall edge: bisection at the long-double precision floor, residual ~1e-3, 8% tail
correction — treat row as approximate. The ω=1.440 row is filled by the script run;
E−mQ crosses zero at ω ≈ 1.443 by interpolation of the 1.42/1.45 rows.)

Findings [estimate]:
- **dQ/dω < 0 over ω ≲ 1.47** (entire useful window) — classically stable branch.
  Near ω → m the branch turns (dQ/dω > 0 between 1.47 and 1.49): thick-wall unstable
  branch, **minimum charge Q_min ≈ 92**.
- **Evaporation bound E < 1.5·Q holds for ω ≲ 1.443**; violated on 1.45–1.49.
  ⟹ **absolutely stable Q-balls exist for ω ∈ (1.309, ≈1.443)** — both criteria met.
- Core amplitude is nearly flat across the window: f(0) ≈ 0.53–0.64 (s₀ = κf⁶/κ stays
  near the κf⁶ = 2…3.5 region); the ball grows by RADIUS, not amplitude, toward thin wall.
- Central energy density e(0) ≈ 1.7–2.3 for all ω.
- Suggested first 3D target: ω ≈ 1.39 (Q ≈ 480, E ≈ 690, r_half ≈ 4.4 — fits a
  64³–128³ box comfortably, well inside both stability criteria, far from both edges).
- dE/dω = ω·dQ/dω verified to ≤0.6% (finite differences) — consistency check on E, Q
  integrals.

## 5. Theta sector — honesty section

With Φ_a = f(r)e^{iωt} (same f for all a), the Φ-vector field is f(r)·(1,1,1)e^{iωt} and

    ∇×Φ = f′(r) ( r̂ × (1,1,1) ) e^{iωt}  ≠ 0

— a transverse, ℓ=1-type source oscillating at frequency ω. [thm]

(a) **Co-rotating Θ response is radiative.** Setting Θ_a = h_a(x) e^{iωt} with m_θ = 0:

    (∇² + ω²) h_a = −η (∇× f(r)(1,1,1))_a

is a Helmholtz problem with k = ω > 0. The Green's function is oscillatory
(e^{iωr}/r): the response has non-normalizable radiating tails, NOT a bound profile.
There is no localized stationary 12-field solution with m_θ = 0 and η ≠ 0. The U(1)
does NOT cure this channel: it forbids only the Φ-sector harmonic radiation (|Φ| static),
not the linear ω-frequency drive of a massless field. Θ radiation carries BOTH energy
and charge (ρ_Q^θ = ω|h|² > 0) out through the absorbing boundary. [thm]

(b) **Experimental design consequence.** The clean 3D stability experiment is staged
(config choice only — the kernel is always 12-field):
  1. **η = 0 run**: the radial solution is exact; expect a stationary Q-ball. Measures
     the integrator's Q-drift floor and lattice-discretization losses.
  2. **η = 0.5 run**: SAME seed; measure the Θ-radiation drain — dE/dt, dQ/dt through
     the sponge, charge retention half-life, and whether the ball settles to a finite-η
     quasi-stationary attractor (drain rate ∝ η² at leading order — perturbative
     prediction, [estimate]) or evaporates. Global Q remains exactly conserved
     (interior + measured outflow); this drain measurement IS the science.

(c) **Curl-source-killing ansatz options (future)** [estimate/open]:
  - **Gradient (curl-free) profiles**: Φ_a = ∂_a g(x) e^{iωt} gives ∇×Φ ≡ 0, removing
    the source entirely. Cost: for g = g(r), Φ = g′(r) r̂ — a hedgehog whose potential
    density acquires the angular factor x²y²z²/r⁶ (sphere average 1/105, zero on the
    coordinate planes), drastically weakening the attractive sextic term; the existence
    window must be recomputed with angular weighting and probably narrows severely.
    Also curl-freeness is not preserved exactly by the nonlinear evolution (the potential
    force of an anisotropic s is not curl-free), so the source is only suppressed, not
    forbidden. [open]
  - **Massive theta**: m_θ > ω makes the Θ response a bound Yukawa dressing instead of
    radiation (Helmholtz → modified Helmholtz). Departs from the locked m_θ² = 0 default;
    a parameter experiment, not a v1 change. [open]
  - **Phase-staggered ansatz** Φ_a = f(r) e^{i(ωt + δ_a)}: the relative phases δ_a do not
    affect s or the radial ODE, and they rotate but do not cancel ∇×Φ (the curl mixes
    different a with gradients of the same f). No cancellation found this way. [estimate]

## 6. Predicted observables for the 3D test (η = 0 first)

From §4, seeding ω = 1.39 (recommended): expect core amplitude |Φ_a| = f(0) ≈ 0.64
(u_a = f cos ωt, v_a = f sin ωt), s₀ ≈ 0.069, central energy density ≈ 2.25,
Q ≈ 482, E ≈ 692, half-amplitude radius ≈ 4.4, rms charge radius ≈ 3.7, internal
rotation period 2π/ω ≈ 4.52 t.u. Diagnostics to add: interior Q(t) (new), boundary
charge flux, E/Q ratio (should hold at 1.435 and stay < 1.5), |Φ| core stationarity
(the v65 P_max oscillation at the 3ω harmonic should be ABSENT — |Φ_a| is
time-independent; this is the smoking-gun signature distinguishing the Q-ball from the
failed real oscillons).

## Files

- `v66/theory/qball.mac` — Maxima verification, 19 checks, ALL PASS (`maxima -b qball.mac`)
- `v66/theory/qball_checks.py` — existence window cross-check + shooting solver + table
- `v66/theory/qball_checks.out` — script output
