# v71 NUCLEI — Composing He and Li from Q-ball nucleons

**Date**: 2026-06-11. **Status**: [measured] — 3-run V100 campaign (he2, li3, li3s;
all rc=0, ~40 min GPU), full snapshot record (121 frames each), renders in
`v71/results/{he2,li3,li3s}.png`. New tool: `sfa/seed/gen_qball_multi.c`
(N-ball seeds). Data: `/space/scp/v70/{he2,li3,li3s}*`. Branch reference:
`v69/theory/gscan.tsv` (g=0.05).

## What a nucleus IS in this theory

Same-charge balls have no non-merging bound state (co-phase contact attraction
is attractive to contact; quadrature/anti-phase configurations don't bind), so
nuclei here are **fused charge droplets** — the liquid-drop analogy is literal.
A nucleus of A nucleons = the fusion product, which forgets nucleon number and
relaxes onto the Q-ball branch at its surviving charge. Identity bookkeeping
("this is He") lives in the formation history, the charge, and the branch point.

## The three compositions (nucleon: ω=1.42, Q=311.5, E/Q=1.4649 unless noted)

| run | recipe | seed Q* | final (t≈272) | verdict |
|---|---|---|---|---|
| he2 | 2 nucleons, D=10, co-phase | 665.7 | **Q=604.6, E/Q=1.4421** | **He exists**: on-branch (pred 1.4406, Δ0.1%), binding ≈1.6% of rest mass, settled |
| li3 | 3 nucleons, triangle side 10 | 1050.9 | **Q=934, E/Q=1.4344, still shedding** | super-critical (Q_max=921): merges into ONE droplet that **evaporates toward the cap** — no fission, no fragmentation |
| li3s | 3 light nucleons (ω=1.46, Q=114.1, E/Q=1.5187) | 446.6 | **Q=356.7, E/Q=1.4652** | **Li exists** (light isotope): clean spherical droplet, binding ≈3.5% per charge, mid-branch |

(*seed Q exceeds A×Q_nucleon by 7–17%: constructive interference of overlapping
co-phase tails — the v66 "initial-charge interference" effect at small D.)

## What the pictures show (`he2.png`, `li3.png`, `li3s.png`)

- **he2**: two balls → contact by t≈30 → single elongated droplet that rings in
  its **quadrupole mode** (the bar visibly reorients prolate-x ↔ prolate-y
  across frames) while shedding radiation rings; one merged charge blob (rhoQ),
  one deformed Coulomb halo relaxing toward spherical.
- **li3**: trefoil → triangular merger intermediate (the binding field s shows
  the triangle filling in) → single central droplet by t≈50 → settled sphere
  with a strong Coulomb ring, surrounded by outgoing arc-shaped shells — the
  evaporation by which the super-critical droplet sheds toward Q_max.
- **li3s**: same sequence at one-third the charge; final frame is a clean
  spherical droplet — the stable Li-analog.

## Quantitative results

1. **Mass defect (fusion energy release), measured**: He lands at E/Q = 1.4421
   vs the free nucleon's 1.4649 — the captured charge is bound by **≈1.6% of
   rest mass**, released as φ radiation during the merger (E: 970.5 → 871.9).
   The light-Li releases ≈3.5% per charge (its loosely-bound ω=1.46 nucleons,
   E/Q=1.5187, fuse to a mid-branch droplet at 1.4652, still settling).
   Both final states sit on (he2: within 0.1%) or just above (still-hot) the
   independently-computed g=0.05 branch — fusion products are legitimate
   branch Q-balls, not new objects.
2. **The fission limit acts by evaporation, not explosion**: li3's merged
   droplet (Q up to 1050) sheds charge monotonically (1050.9 → 934 over
   300 t.u., decelerating) toward Q_max = 921, remaining a single coherent
   droplet throughout. Combined with v69's fiss1 (a COLD spherical
   super-critical ball survives unperturbed), the picture: super-critical
   charge is removed by radiation from hot/deformed states; binary fission has
   not yet been observed and may need the ℓ=2-seeded cold test (v69 next #1).
3. **Quark bookkeeping persists through fusion**: per-component charges remain
   exactly Q_total/3 at every stage (the composite has 3A quarks in 3 equal
   component pools — components, not nucleons, are the conserved structure).

## Caveats

- T=300 endpoint states are still draining slowly (he2 −0.05/t.u., li3
  −0.2/t.u.); "final" Q values are upper bounds on the parked values. li3's
  approach to exactly Q_max is extrapolated, not observed.
- The interference excess in seed Q means "A nucleons" is approximate at
  D=8–10; cleaner accounting would seed at larger D and pay longer infall.
- η=0 throughout (θ sector inert); the gauged drain channel is closed by
  m_θ=1.6 — these are the stable-matter conditions established in v68/v69.

## Files

- `v71/NUCLEI.md` (this), renders `v71/results/{he2,li3,li3s}.png`
- `sfa/seed/gen_qball_multi.c` — N-ball seed generator (new)
- SFAs + diags: `/space/scp/v70/{he2,li3,li3s}*` (121 frames each)
