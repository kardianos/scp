# v71 QUARK — Is the Q-ball a quark or a proton? (Component-substructure verification)

**Date**: 2026-06-11. **Status**: [measured] — 5-run V100 campaign (q1–q5, all rc=0,
~25 min GPU total) + per-component analysis + flavored-baryon BVP solver.
**Question (user)**: verify whether the Q-ball is a quark rather than a proton,
and how a proton would be built from quarks.

New tools: `sfa/analysis/sfa_qcomp.c` (per-component Q/mass/centroid/rms),
`sfa/seed/gen_qball_quark.c` (component-displaced "quark" seeds),
`v71/analysis/flavored_qball.py` (3-field Newton BVP, asymmetric baryons),
`sfa_slice`/`render_slices.py` extended with per-component densities
rho2_0/1/2. Data: `/space/scp/v70/q*_*.sfa`, `v71/results/`.

## 0. Verdict

**The Q-ball is not a quark — it is the baryon.** The quark-level objects are
the three internal component fields Φ_a, and they pass every quark-style test
the theory admits:

| property | quark (QCD) | component Φ_a (this theory) | evidence |
|---|---|---|---|
| fractional charge | e/3 units | **exactly Q/3 each** (105.06×3 = 315.17) | sfa_qcomp on banked ball, stable in t |
| no isolated existence | confinement | **single-component lump cannot bind** (force ∝ Π_{b≠a}\|Φ_b\|² ≡ 0 [thm]; disperses as free KG) | q1: s ≡ 0 exactly, rms 3.4→13.2, 60% charge lost by t=100 |
| bound only in color-neutral combos | baryons + mesons | **only all-three combos bind** | q2: two components, s ≡ 0, same dispersal — **NO meson sector** (structural difference from QCD) |
| baryon = 3 bound quarks | proton = uud | symmetric ball = 3 co-located component lumps | q3/q4: assembly (below) |
| flavor structure | u vs d | **asymmetric (ω_0,ω_1,ω_2) stationary baryons exist**, Q_0/Q_1 down to 0.74 at Δω=0.04 | flavored_qball.py continuation (validated: symmetric case reproduces the v66 branch, Q=209.5 ✓) |

Caveats up front: the "color" here is a U(1)³ structure (only the diagonal is
gauged), NOT SU(3); there is no gluon sector and no linear-confinement string —
isolated components disperse rather than pull back with constant tension; and
the component index doubles as a spatial index in the η-curl term (η=0
throughout this campaign). The mapping is structural, not quantitative.

## 1. Free quark test (q1) and meson test (q2) — both die [measured + thm]

Force on Φ_a is −2Vt′(s)Φ_a Π_{b≠a}|Φ_b|²: with any component identically
zero, ALL components are free Klein–Gordon fields [thm, from the v66-verified
EOM]. Measured: a single-component lump of charge 105 (one quark's worth)
keeps s_max = 0 to machine zero for all t, spreads (rms 3.4 → 13.2 by t=100),
and loses 60% of its charge to the boundary; renders show pure dispersal
rings, the binding panel black, the Coulomb halo evaporating. Two co-located
components (q2) behave identically (phi_max 0.635 → 0.06). **Quarks and
"mesons" cannot exist; only three-component ("color-complete") objects bind.**

## 2. Proton assembly from three quarks (q3–q5) — the binding basin [measured]

Three single-component lumps at triangle radius d (each Q_a = 105, ω = 1.42
seeding, N=128, L=20, absorbing BC, full snapshot record):

| d | infall time | captured charge/quark | final object |
|---|---|---|---|
| 2 | ~15 t.u. | 75.3 of 105 (72%) | bound ball, Q ≈ 226, settled (drain → 0) |
| 4 | ~60 t.u. | 48 of 105 (46%) | bound ball, Q ≈ 144, still ringing at t=150 |
| 6 | ~120 t.u. | ~33 of 105 (31%), still falling | small remnant baryon amid dispersal fog |

Mechanism (visible frame by frame in the renders, `v71/results/q4.png`):
binding ignites ONLY in the central region where all three overlap; each
quark's core falls toward it while its outer envelope — not yet captured —
disperses as free waves. Assembly is a race between infall and dispersal:
the farther apart the quarks start, the smaller the baryon that survives.
Per-component charges remain exactly equal throughout every assembly
(e.g. 75.33/75.32/75.32 at q3 end) — capture is democratic.

**This is "how you make a proton from quarks" in this theory**: place three
component-lumps within roughly two core radii (d ≲ 4 ≈ r_half) and the baryon
self-assembles, radiating the excess as binding violence; beyond that, most
of each quark is lost to dispersal first.

## 3. Flavored baryons (proton vs neutron analog) [solved, stability open]

`flavored_qball.py`: Newton relaxation of the 3-component radial system with
unequal frequencies, Φ_a = f_a(r)e^{iω_a t} (s stays static, ansatz exact).
Continuation from the symmetric ω=1.42 solution converges smoothly to at
least Δω = 0.04 (ω = 1.38/1.42/1.42): localized stationary baryons with
**unequal charge partition** Q_a = (76.7, 103.7, 103.7) — the in-model analog
of uud vs udd flavor content. Energy rises along the branch (E 321 → 415).
Open: dynamical stability of the flavored branch (needs a 3D run with a
flavored seed — generator extension: per-component ω/profile).

## 4. Honest differences from QCD (do not oversell)

1. No mesons of any kind (the product potential binds only triples) — real
   QCD's richest sector is absent.
2. "Confinement" is non-existence (dispersal), not a linear string: pulling a
   quark out costs no growing energy; the quark simply stops being a particle.
3. Color is U(1)³/permutation, not SU(3); nothing rotates components into
   each other dynamically (at η=0), so per-component charge is separately
   conserved [verified: Q_a flat in every run].
4. The gauge charge is the diagonal U(1) — all three quarks carry the SAME
   sign of it; there is no analog of charge −1/3 vs +2/3 within one baryon.

## Files

- `v71/QUARK.md` (this), `v71/analysis/flavored_qball.py`,
  `v71/results/{q1,q4,q5}.png`, `v71/results/q*_qcomp.tsv`,
  `v71/results/flavored_profile_last.txt`
- `sfa/analysis/sfa_qcomp.c`, `sfa/seed/gen_qball_quark.c` (new tools)
- SFAs: `/space/scp/v70/q{1,2,3,4,5}_*.sfa` (per-frame visual record)
