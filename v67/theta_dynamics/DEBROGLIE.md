# DEBROGLIE — c as the rate of locality: de Broglie kinematics of the Q-ball

**Date**: 2026-06-10. Workpackage D of the v67 theta-dynamics series.
Framework: v66/THEORY.md §1–4 (complexified 6-field Cosserat, Q-ball branch),
THETA_DYNAMICS.md §1–2 (η≠0 medium, bath ladder). Standard parameters
m² = 2.25 (m = 1.5), μ = −41.345, κ = 50, η = 0.5, m_θ² = 0; c = 1 throughout
(the lattice signal speed of the φ sector).

**Verification**: `deBroglie.mac` — 28/28 PASS (`maxima -b deBroglie.mac`,
output `deBroglie.out`); `deBroglie.py` — 6/6 PASS (output `deBroglie_py.out`).
Claim markers: [thm] exact symbolic, [verified] machine-checked numeric,
[estimate], [open].

**Driving idea (user)**: c is the rate of locality of the medium; a particle is
a field pattern with an internal beat (the Q-ball's e^{iωt}). Then relativistic
AND quantum kinematics should both be *consequences* of how a beating pattern
must transform when it moves. This document derives that chain and finds one
structural surprise: **the ball's quantum of action is its charge** (§2).

---

## 1. The boosted Q-ball is an exact relativistic particle (η = 0)

The η = 0 φ sector is exactly Lorentz invariant (THETA_DYNAMICS §2.4, F11), so
boosts act on the v66 Q-ball Φ_a = f(r)e^{iωt} exactly.

**Exact moving solution** [thm, D1]: for any U(1) potential W(|Φ|²) (the v66
three-field potential is pointwise in the |Φ_a| and acts componentwise on the
symmetric ansatz), if Φ = F(x,y,z)e^{iωt} solves the EOM, then so does

    Φ_v(t,x) = F(γ(x−vt), y, z) · e^{iγω(t−vx)},    γ = 1/√(1−v²).

The envelope is the x-contracted profile; ALL the relativistic and de Broglie
kinematics live in the phase θ = γω(t − vx). Machine-checked by full chain rule:
the boosted residual reduces exactly to the rest profile equation (D1).

Phase kinematics [thm, D2–D5]:

| quantity | result | check |
|---|---|---|
| (a) de Broglie tilt k_dB = −∂θ/∂x | **k_dB = γωv** | D2 |
| (b) frequency at the moving center (x = vt) | **ω/γ** — time dilation: the internal beat slows | D3 |
| (c) frequency at a fixed lab point | **γω** — the quantum energy, not the clock | D4 |
| (d) phase velocity Ω/k = 1/v (>1, carries no energy); group velocity dΩ/dk = v | **v_ph·v_g = c² = 1** | D5c,d |
| dispersion of the packet | Ω(k) = √(ω² + k²); (γω, γωv) sits on it | D5a,b |

(b) vs (c) is exactly the twin/transverse-Doppler split: the comoving beat is
red by γ, the lab-frequency is blue by γ; both follow from the single phase θ.

**Energy–momentum** [thm, D6]: with E = γE₀ and p = γE₀v,

    E² = p² + E₀²   exactly  —  the ball IS a relativistic point particle of
    mass E₀ (= 692 at ω = 1.39), with internal (Compton-like) beat ω = 1.39.

**But ω ≠ E₀.** A textbook quantum particle has internal frequency = its rest
energy (Compton clock, ℏ=1). Here ω = 1.39 while E₀ = 692. The de Broglie
relation therefore reads [thm, D7]

    k_dB = γωv = (ω/E₀) · p  =  p / h_eff,      h_eff ≡ E₀/ω = 497.8 at ω=1.39.

The ball has an **effective ℏ that is a property of the ball, not of the
theory** — and h_eff turns out to be its charge (§2).

**In-sim consequence** [thm → testable]: to set a ball moving, stamp the seed
with the phase tilt e^{−iγωv x} and the contracted envelope. A moving ball's
snapshot must show phase stripes of wavelength 2π/(γωv) (= 45 code units at
v = 0.1, ω = 1.39) and a core beat slowed to ω/γ. This is a cheap, sharp test
of the whole §1 chain (one N=192 run with a tilted seed).

---

## 2. The action quantum is the charge: E = ωQ + surface term

### 2.1 Exact identity [thm, D8]

On the radial Q-ball (v66 §3), with e − ωρ_Q = 3(½f′² + U_ω(f)) pointwise
[D8a] and the Derrick/Pohozaev scaling identity K + 3P_U = 0 for the ODE
solution [D8b]:

    E − ωQ = ∫ f′² dV  ≡ G  >  0    (exact; D8c)

E exceeds ωQ by exactly the (positive) gradient/surface energy — a
surface-to-volume term that vanishes fractionally in the thin-wall (large-Q)
limit.

**Numerical verification** [verified]: computing E, ωQ and G = 4π∫f′²r²dr
independently from all 37 stored v66 profiles, the identity holds to
**3×10⁻⁵ relative** on 35/37 (the two ω ≤ 1.32 thin-wall edge rows sit at the
shooting bisection floor flagged in v66/THEORY.md §4 and are excluded;
`deBroglie_py.out` §2).

### 2.2 Free-quantum normalization [thm, D9]

A circularly polarized free wave has ρ_E/ρ_Q = ω exactly (D9): each unit of Q
is one quantum carrying E = ω·1, i.e. action 2π per beat period. **Q counts
quanta, and the field's bare action quantum is ℏ = 1** in code units.

### 2.3 The table — E/(ωQ) across the branch [verified]

From `v66/results/scan.tsv` (full table in `deBroglie_py.out`):

| ω | Q | E | h_eff = E/ω | **E/(ωQ)** |
|------|---------|---------|---------|--------|
| 1.3150 | 203356 | 269285 | 204779 | 1.0070 |
| 1.3300 | 22207 | 29786 | 22395 | 1.0085 |
| 1.3600 | 1745.0 | 2422.9 | 1781.5 | 1.0209 |
| 1.3900 | 482.2 | 691.9 | 497.8 | 1.0323 |
| 1.4200 | 209.6 | 309.7 | 218.1 | 1.0405 |
| 1.4350 | 153.0 | 229.0 | 159.6 | 1.0430 |
| 1.4800 | 86.7 | 132.8 | 89.7 | 1.0349 |

    E/(ωQ) ∈ [1.0046, 1.0430] on the stable branch (ω ∈ [1.315, 1.435]);
    max 1.0438 over the whole scan.  E = ωQ to ≤ 4.4% EVERYWHERE, and the
    deviation is exactly G/(ωQ) [thm], shrinking toward thin wall (large Q):
    0.6% at Q = 4.9×10⁴, 3.2% at Q = 482, 4.3% at Q = 153.

### 2.4 Statement of the result

**The Q-ball's quantum of action is its charge** [thm + verified]:

    h_eff = E₀/ω = Q · (1 + ε),     ε = G/(ωQ) ∈ [0.005, 0.044]  (→ 0 thin-wall)

    de Broglie relation:  k_dB = p / (Q(1+ε))  ≈  p/Q
    action per internal period:  S = E·T = 2π E/ω = 2π Q (1+ε)
    ⟹ each charge quantum carries action 2π(1+ε) ≈ h  (ℏ = 1).

The moving ball's phase tilt is the **per-quantum** de Broglie wavenumber
k_dB ≈ (p/Q)/ℏ: a packet of Q phase-locked quanta sharing ONE phase, each
obeying the standard ℏ = 1 relation, with the binding correction ε given in
closed form by the surface term.

**Honesty — what this is and is not**:
- This is classical field theory; nothing here quantizes Q. The statement is
  structural: IF Q comes in integers (as it would on quantization, since free
  quanta carry unit Q — D9), THEN the ball's action comes in units of h per
  quantum, with a computable O(ε) binding correction.
- The sign matters: a quantum-mechanical bound state of N constituents
  diffracts with k = p_total/ℏ (phases add across constituents). The classical
  ball instead has one COLLECTIVE phase, so it diffracts at the per-quantum
  scale k = p/Q — a factor-Q-longer wavelength. Measuring the interference
  scale of a moving ball distinguishes "coherent classical condensate"
  (k = p/Q, this theory) from "N-particle quantum composite" (k = p). This is
  the cleanest structural fingerprint the model makes about its own quantum
  limit. [thm given §1; interpretation [open]]
- ω < E₀/Q < m on the stable branch: the per-quantum energy E/Q is the
  internal frequency ω plus the surface share — consistent with E/(mQ) < 1
  (binding) and E/(ωQ) > 1 (gradient cost) simultaneously.

---

## 3. Medium corrections (η ≠ 0): de Broglie through the theta vacuum

With η ≠ 0 the φ sector is preferred-frame (THETA_DYNAMICS F11); the matter
branch is anisotropic: c_eff² = 1 + η²/m² for the two transverse
polarizations, 1 longitudinal (§2.2 there). Treating the moving ball as a
wavepacket on the matter branch Ω(k)² = ω₀² + c_eff²k² (linear-response proxy
for the nonlinear ball [estimate]):

**Lorentz-violating de Broglie relation** [thm on the dispersion, D10a,b]:

    k_dB(v) = ω₀ v / (c_eff √(c_eff² − v²))
            = γω₀v · [ 1 − (η²/m²) · (2−v²)/(2(1−v²)) + O(η⁴) ]

- At standard η = 0.5: η²/m² = 1/9, so the **transverse phase tilt is reduced
  11.1%** at low v (factor 0.8889, D10d); longitudinal tilt is uncorrected —
  the de Broglie wavelength of a moving ball is **birefringent/anisotropic**
  at O(η²), anchored to the theta frame (genuine Lorentz violation, not
  reciprocal relativity).
- In the medium v_ph·v_g = c_eff² ≠ c² [thm, D10c] — the §1 de Broglie product
  theorem picks up exactly the medium factor.
- The composite ball mixes polarizations; the observable tilt correction is an
  angular average between 0 and −11.1% depending on motion direction relative
  to the constituent-gradient geometry [estimate].

**Bath anchoring** [estimate]: in a theta bath of density e (bath frame), the
internal frequency runs: ω₀ → ω₀ + δω(e), δω(e) = +0.0119e − 0.0032e²
(THETA_DYNAMICS §2.5), and the O(e²) piece rescales with motion through the
bath as δm²(v) ≈ δm²(0)(1 + 8v²/3) (§2.4 there). So

    k_dB(v, x) = γ v [ ω₀ + δω(e(x)) ] / c_eff²(local) + …

— a position-dependent phase tilt. Its gradient is a force: §4.

**Caveat** [open]: at η ≠ 0 a moving ball is not an exact solution at all (the
massless-θ channel radiates, v66 §5); these are kinematic corrections assuming
the quasi-stationary finite-η attractor survives transport.

---

## 4. Clock-gradient force: weak-field-GR-like drift

Eikonal rays on Ω(x,k) = √(ω(x)² + k²) (internal clock rate ω(x) varying
slowly in space) give [thm, D11a,b]:

    dΩ/dt = 0 along rays,   and   a = −c² (1−v²) · d(ln ω)/dx
    v → 0:   a = −c² · d(δω/ω)/dx

— exactly the weak-field GR form a = −∇Φ with Φ/c² = δω/ω: **the ball falls
toward regions where its own clock runs slower.** Gravity-as-clock-gradient
emerges from pure phase kinematics; no new coupling is introduced.

**Number for the v67 bath ladder** [estimate]: δω(e) = +0.0119e − 0.0032e²
(blue with e), one rung across the box de/dx = 0.2/25 = 0.008, ω = 1.39:

    a = −(1/ω)(dδω/de)(de/dx) = −6.8×10⁻⁵ code units   (toward LOW e_bath)
    (−6.85/−6.48/−6.11 ×10⁻⁵ at e₀ = 0/0.1/0.2; D12a, deBroglie_py.out §3)

**V51 comparison** [estimate, order of magnitude only]: the measured V51
gradient drift (FUTURE.md F25, real-field proton in an A_bg gradient — a
DIFFERENT background and field content) was +6 code units in T = 400, i.e.
a_V51 ≈ 2x/T² = 7.5×10⁻⁵ (consistent with the quoted mean drift rate
0.015/t.u. = aT/2). Ratio a_pred/a_V51 = 0.91 — **the same order of
magnitude** (D12b). Not a fit; the agreement of scales across two different
backgrounds supports clock-gradient kinematics as the common mechanism behind
F25-type drifts. A direct test: run the bath ladder with a linear e_bath(x)
ramp and check both the magnitude and the 1/ω scaling.

---

## 5. Phase closure: orbit quantization [estimate]

Two balls in mutual orbit (separation D, relative velocity v_rel) feel a
phase-DEPENDENT force (v66 §5 taxonomy: co-phase attracts — tb1 merger;
anti-phase repels — tb2 ejection). On a circular relative orbit the
accumulated de Broglie phase around the orbit is ∮k_dB·dl with
k_rel = (ω/E₀)p_rel = ωv_rel/2 (reduced mass M/2, h_eff = E₀/ω from §2).
Single-valuedness of the relative phase — Bohr–Sommerfeld closure — requires
[thm arithmetic, D13a]:

    2πD · (ωv_rel/2) = 2πn   ⟹   D · v_rel = 2n/ω = 1.439 n   (ω = 1.39)

Physical meaning (classical, no quantization postulate): off closure, the
relative phase slips at beat rate Δω_rel, so the inter-ball force oscillates
between attraction and repulsion (tb1 vs tb2!) and time-averages away; ON
closure the force is stationary in the co-orbiting frame. Quantized orbits are
those where the phase-locked force can persist.

Dynamics estimate (n = 1, M = 692, tail decay μ_t = √(m²−ω²) = 0.564, tail
force band calibrated crudely from tb2 ejection energetics, U₀ ∈ [0.2, 2]):

| D | v_rel(n=1) | F_req = 2M/(ω²D³) | F_tail range | T_orb |
|----|-------|-------|---------------|------|
| 8 | 0.180 | 1.40 | 1.6 – 16 | 279 |
| 10 | 0.144 | 0.72 | 0.42 – 4.2 | 437 |
| 12 | 0.120 | 0.41 | 0.11 – 1.1 | 629 |
| 14 | 0.103 | 0.26 | 0.03 – 0.31 | 856 |

The n = 1 force balance crosses inside **D ≈ 9–14** — inside the
simulation-accessible window D ∈ [8, 20], L = 30 [verified arithmetic, P4].

**Experiment that would see it** (config/seed level, no kernel change):
`gen_qball_pair` co-phase seeds at D = 12 with TRANSVERSE boost tilts
e^{∓iγωv y}, v = v_rel/2 ≈ 0.06 (the §1 boost recipe — angular momentum
prevents the tb1 head-on merger). N = 256, L = 30, T ≈ 1300 (two orbits at
D = 12). Compare v_rel on/off closure (e.g. 0.120 vs 0.090 at D = 12):
prediction — near-closure initial conditions lock (D(t) oscillates about a
fixed radius, relative phase bounded); off-closure ones phase-slip, the
time-averaged binding dies, and the pair drifts apart or inspirals to merger.
Frame cadence ≤ 10 t.u. to resolve the relative phase (beat period
2π/Δω_rel ≈ 100 t.u. at 10% mismatch). Cost: ~2× a tb run. Caveats: the
tail-force calibration is order-of-magnitude; radiation drag (η = 0.5) adds a
secular inspiral on top — run η = 0 first [estimate/open].

---

## 6. Headline results

1. **Boosted ball** [thm, D1–D7]: Φ_v = f(r′)e^{iγω(t−vx)} is exact (η = 0);
   k_dB = γωv, comoving beat ω/γ, lab frequency γω, v_ph = 1/v, v_g = v,
   E² = p² + E₀² with E₀ = 692 at ω = 1.39. The ball is a relativistic
   particle whose quantum kinematics follow from its internal beat.
2. **Action quantum = charge** [thm + verified]: E − ωQ = ∫f′²dV exactly
   (verified 3×10⁻⁵ on 35 profiles); E/(ωQ) ∈ [1.005, 1.044] across the whole
   branch ⟹ h_eff = E/ω = Q(1+ε); k_dB = p/(Q(1+ε)); each charge quantum
   carries action 2π(1+ε) ≈ h. De Broglie per charge quantum — and a sharp
   classical-condensate vs quantum-composite fingerprint (k = p/Q vs p).
3. **Medium correction** [thm on dispersion/estimate for ball, D10]:
   k_dB = γωv[1 − (η²/m²)(2−v²)/(2(1−v²))]; 11.1% transverse reduction at
   η = 0.5, longitudinal unchanged — anisotropic, theta-frame-anchored
   (Lorentz-violating) de Broglie relation; v_ph·v_g = c_eff².
4. **Clock-gradient force** [thm kinematics, D11–D12]: a = −c²(1−v²)d(lnω)/dx;
   bath ladder number a = −6.8×10⁻⁵ (one rung across the box), V51 measured
   scale 7.5×10⁻⁵ — same order across different backgrounds.
5. **Orbit quantization** [estimate, D13]: closure D·v_rel = 1.44n; the n = 1
   orbit falls at D ≈ 9–14, inside the accessible range; phase-slip
   (attract/repel averaging) is the classical mechanism that selects quantized
   orbits; concrete two-ball experiment specified.

## Files

- `deBroglie.mac` / `deBroglie.out` — symbolic verification, 28/28 PASS
- `deBroglie.py` / `deBroglie_py.out` — action-quantum table, profile identity
  verification (35 profiles), force numbers, orbit table; 6/6 PASS
- Inputs: `v66/results/scan.tsv`, `v66/results/profile_omega*.txt`,
  `v66/THEORY.md`, `THETA_DYNAMICS.md`, FUTURE.md F25, tb1/tb2 cluster tables
