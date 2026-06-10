# v67 FINDINGS — Theta Sector, Bath, Fragmentation, and Emergent Kinematics

**Date**: 2026-06-10/11. **Status**: [measured] N=192–256 V100 campaign (11 runs, all rc=0)
+ Maxima-verified theory packages. Data: `v67/results/*_diag.tsv`, SFAs in `/space/scp/v67/`.
Theory: `v67/THEORY.md` (condensate/bath), `v67/theta_dynamics/THETA_DYNAMICS.md` (frame/
winding/forcing), `v67/theta_dynamics/DEBROGLIE.md` (kinematics). New instruments (v66
kernel diag upgrade, verified): `omega_core` (direct internal-frequency measurement),
`Q_core` (interior charge), θ probe columns (radiation spectra).

## 1. The θ-boundedness question (user-posed) — answered experimentally

The complexified theory's Θ is massless + U(1)-charged + free — a combination nature
forbids (photon neutral; gluons confined; W massive). Measured resolutions:

- **th1 (m_θ=1.6 > ω=1.39, η=0.5, T=300)**: after a one-time ~4% dressing transient the
  drain CLOSES — Q_core −0.0025/t.u. (~400× below massless baseline), no branch slide,
  core fully recovered. **A single parameter converts the theory from "all matter is
  soluble" to "matter absolutely stable as a Yukawa-dressed ball."**
- **th2 (m_θ=0.7 < ω)**: channel stays open (early drain identical to massless), late
  drain mildly slowed (Q=291 vs 253 at t=200) — kinematic phase-space reduction.
- Cost of the mass route: kills the long-range mediator (v64 H1). The structural fix
  that keeps BOTH stable matter and long-range fields is **gauging the diagonal U(1)**
  (θ → neutral connection; Gauss law pins charge; amplitude unboundedness becomes pure
  gauge) — design doc is the standing next step (task #10, FUTURE.md).

## 2. Edge-ball death mechanism — measured directly with omega_core

fd1 (ω=1.46 seed, η=0.5): omega_core reads 1.46000 at t=0, is kicked to **1.525 by t=5**
(past the VK turn 1.485 and above m=1.5, off the existence branch), core collapses
monotonically, gone by t=80. The v66 "b3 evaporation" is definitively reclassified:
**shock-driven frequency excursion off the solution branch.** The E=mQ line plays no
dynamical role (v66 c1 control: frozen at η=0 for T=300).

## 3. Condensate fragmentation — threshold theorem confirmed BOTH ways + new η-instability

Theory (`v67/THEORY.md`, Maxima 20/20): instability iff κA⁶ < 1/2, i.e. A < A* =
(2κ)^(−1/6) = 0.46416 [thm]; fastest mode k_max, σ_max from closed form; full 48×48
12-field eigenanalysis matches analytic to 6e−17.

| run | A | η | predicted | observed |
|---|---|---|---|---|
| cf2 | 0.40 | 0 | fragment, σ=0.290, t_frag~16–25 | **fragments**, growth ≈0.23–0.24/t.u., saturated t≈30 |
| cf3 | 0.585 | 0 | **must NOT fragment** (falsifier) | **flat for T=200** (s_max 0.044–0.049, zero growth) |
| cf4 | 0.585 | 0.5 | NEW first-order-in-η parametric instability, σ≈0.058–0.067 | **fragments from t≈100**, growth ≈0.078/t.u. (within ~16%) |
| cf1 | 0.40 | 0.5 | as cf2, (111)-direction unchanged | fragments (diag-level; spectra TBD) |

The η-only destabilization of a saturated condensate (cf3 stable ↔ cf4 fragments,
identical seed) is a **qualitatively new mechanism**: the torsion coupling CREATES
particles from dense stable matter — directly relevant to the star-interior program.
Fragment census (sfa_qball_track, cf2): interconnected web at t≈25 individualizes to
~63 distinct lumps by t=200 at thr=0.10 (spacing ≈12.6 vs λ_max=7.9 — right scale,
factor ~1.6, with visible merger/coarsening). Q conserved to 1e−4 through fragmentation.
Per-fragment charge spectrum: open (TSVs available).

## 4. Bath ladder — neutral bath ERODES; dressed states are real and strong

(η=0.5, periodic, T=300; Q_core = charge in r<8; theory ladder from `v67/THEORY.md` §3)

| e_bath | Q_core(299) | matched-Q clock shift δω (at Q=300) |
|---|---|---|
| 0 (control) | 267.9 | — (ω=1.4142) |
| 0.05 | 244.7 | **+0.012** (blue) |
| 0.2 | 223.5 | **−0.025** (red) |
| 0.8 | 122.9 | **−0.120** (deep red; ω=1.294 < ω_min=1.3087 — OFF the vacuum branch) |

- The self-bath alone (periodic control vs absorbing v66) re-feeds the ball ~25%.
- A charge-NEUTRAL bath cannot return charge (theory caveat confirmed): erosion
  accelerates with intensity — **detailed balance requires a charged (co-rotating)
  bath**: the sharp follow-up experiment (extend gen_qball_bath with coherent U(1)
  polarization).
- Clock shifts: blue at weak bath (predicted sign; magnitude ~20× the vacuum-perturbative
  estimate), red from e=0.2 (in-core mean-field well dominates — flagged as a systematic
  in WP-A, measured to be the leading effect). At e=0.8 the ball lives **below the vacuum
  existence window** — a strongly dressed state outside vacuum theory. [open: in-core
  dressed-state theory]

## 5. Neutral remnant (rem1, tb3 restart, T=1000) — identity solved

The v66 "persistent neutral two-lobed remnant" is a **standing-wave breather**: density
passes through ~0 every 2.2 t.u. ≈ π/ω (equal ±ω superposition — effectively a
real-sector oscillon, no charge protection). It radiates with strongly decelerating
decay: dE/dt = −0.76 → −0.06 across t=0–1000, E: 463 → 182, still re-forming cores
(s_max 0.066 at end). Long-lived (multi-1000 t.u.) but mortal — fully consistent with
the v65 real-sector theorems; only charged objects are absolutely protected.

## 6. Emergent relativistic-quantum kinematics (user's "c = rate of locality" idea)

`v67/theta_dynamics/DEBROGLIE.md` (Maxima 28/28 + numerics 6/6):

1. **Boosted ball = exact relativistic particle** [thm]: Φ = f(r′)e^{iγω(t−vx)} exact;
   k_dB = γωv, comoving clock ω/γ, v_ph·v_g = c², E² = p² + E₀² with E₀=692.
2. **The action quantum is the charge** [thm + verified]: Pohozaev identity
   **E − ωQ = ∫f′²dV > 0 exactly** (3e−5 on 35/37 profiles); E/(ωQ) ∈ [1.005, 1.044]
   on the branch → h_eff = E₀/ω = Q(1+ε). Each charge quantum carries action ≈ h;
   stability, mass, and quantization are the SAME conserved quantity. Fingerprint
   experiment: boosted-ball phase tilt k = p/Q (classical phase-locked) vs k = p
   (true quantum composite) — measurable, config-level.
3. **Clock-gradient force** [thm kinematics]: a = −c²(1−v²)∇(ln ω); one bath rung
   across the box gives a = 6.8e−5 vs V51 measured drift ≈ 7.5e−5 (**ratio 0.91**,
   order-of-magnitude, different backgrounds) — first mechanism candidate for the
   program's oldest unexplained result.
4. **Orbit quantization** [estimate]: phase closure D·v_rel = 2n/ω puts the n=1 orbit
   at D ≈ 9–14 — inside the simulatable window; experiment specced.

θ-dynamics package (`THETA_DYNAMICS.md`, 65 checks): unifying theorem — the η-coupling
is bilinear and V is Θ-free, so NO θ background enters the linear fluctuation operator;
all influence is via induced Φ condensate (relative mass O(η⁴B⁴), anisotropic,
bath-frame-anchored LV) or forcing (Mathieu tongues at W=2ω₀/n; in-band n=1 sum-beat
at W=ω). Winding: box-quantized Beltrami sectors; chiral splitting and AB phase are
exactly ZERO at this coupling (gauging would restore 2πqn holonomy). Honest nulls
included.

## 7. Open items / next experiments (priority order)

1. **Gauge the diagonal U(1)** (design doc) — resolves θ-boundedness while keeping
   long-range fields; restores AB/winding topology. [task #10]
2. **Charged (co-rotating) bath** — the actual detailed-balance test. [gen extension]
3. **Boost experiment** — k=γωv, ℏ_eff=Q fingerprint, time dilation. [config-level]
4. **Monochromatic drive at W=ω** — Mathieu n=1 tongue resonance curve. [gen config]
5. **In-core dressed-state theory** — explain the matched-Q clock shifts and the
   off-branch e=0.8 state. [theory]
6. **Quantized-orbit run** (D≈9–14, phase-tilted pair, T≈1300). [seed extension]
7. Fragment charge spectrum (Affleck–Dine distribution) from cf TSVs. [analysis]

## Infrastructure

Five Vast.ai V100 instances died during v66–v67 (4.5h/4.5h/35m/25m/—); v67c survived
its full 3h campaign. Standard practice now: self-serializing chained run scripts,
eager per-run diag pulls on a 6-min cron loop, SFA pulls at closeout, no reliance on
channel notifications or auto_download (script-mode outputs untracked).
