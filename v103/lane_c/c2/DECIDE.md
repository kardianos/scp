# v103 lane c, step c2 — decision (two things: a sign, and something without one)

**Written:** 2026-10-11, after the registered runs. Gates: `PREREG.md`.
A measurement in models; it rules nothing.

- Tool: `cas/gap1d.c` (C). T0 passes; see `cas/out/t0.out`.
- Grid: `cas/run_c2.sh`, 32 runs (16 sets × L = 12, 16).
- Tables: `cas/tab.go` (Go), output in `cas/out/tab.out`.

**Frame (R14).** The medium is a ring of cells that are never empty.
Distance is counted in quanta. A gap is a missing quantum, read from a
cell's own count. Nothing is a container.

## Ending: **U** by the registered bars — but read the unlike leg

**The registered result.** In the random process with a flat level
(sets **d** and **bd**), gap pairs attract and excess pairs attract, at
both ring sizes. That meets the registered U.

**Read further, it is not "everything attracts everything".**
- The unlike leg (gap with excess) was reported, not scored. That was
  because in a *coherent* state, creation brings a gap and an excess
  together.
- In the random process there is no creation bias, and there the
  unlike pairs **repel**: ρ_ge = 1.19 and 1.13, V(contact) = +1.88.
- So d gives **like with like**: gaps cluster with gaps, excesses with
  excesses, and the two kinds avoid each other.
- That is surface tension. A flat level charges only for boundaries, so
  a cluster, which has fewer boundaries, is cheaper.
- **It acts at contact only.** V(D) is exactly constant beyond D = 0.

**The other attraction is the one you named: tweaked links with
inertia.** It is exactly Pretko's structure, in a saturated medium
(abc, abcd):
- an isolated gap cannot move;
- two gaps can only breathe between contact and three cells apart, and
  they cannot translate;
- gaps far apart are frozen and can never come together.

They bind by **B = −4t exactly** (stiff), and **−5.12t** (flat), at both
sizes. This too is contact range. In a saturated medium there are no
excesses, so its sign cannot be tested.

**Nothing reaches beyond contact in one dimension.** The only profiles
that extend across the ring are *repulsive*: hard-core gaps carried by
coherent ordinary re-forming.

## 1. The ledger (L = 12 | L = 16)

| set | B_gg | ρ_gg (W_min) | B_ee | ρ_ee (W_min) | verdict |
|---|---|---|---|---|---|
| a, ab (random, saturated, stiff) | 0 \| 0 | 1.000 (1) | — | — | nothing: counting only |
| ad, abd (random, saturated, flat) | −0.771 \| −0.616 | 0.463 \| 0.540 (1) | — | — | gaps cluster, at contact |
| none, b (random, headroom, stiff) | +0.283 \| +0.214 | 1.000 (0.62 \| 0.40) | same | same | nothing |
| **d, bd** (random, headroom, flat) | −0.115 \| −0.073 | 0.463 \| 0.540 (0.25 \| 0.08) | −0.115 \| −0.073 | 0.463 \| 0.540 | **like with like, at contact; unlike repel (ρ_ge 1.19 \| 1.13)** |
| ac, acd (coherent hop, saturated) | +0.273 \| +0.154; +0.203 \| +0.122 | 1.43 \| 1.42; 1.32 \| 1.34 | — | — | repel (hard-core exclusion) |
| **abc, abcd** (coherent tweaked links, saturated) | **−4.000 \| −4.000**; **−5.123 \| −5.123** | 0.440 \| 0.306; 0.333 \| 0.232 (1) | — | — | **bound by mobility (Mach), at contact** |
| c, cd (coherent hop, headroom) | +0.618 \| +0.471; +0.252 \| +0.199 | 1.27 \| 1.27; 1.12 \| 1.15 | +0.695 \| +0.515; +0.345 \| +0.252 | 1.30 \| 1.30; 1.22 \| 1.22 | both like pairs repel |
| bc, bcd (coherent tweaked links, headroom) | −0.871 \| −0.863; −0.769 \| −0.751 | 1.61 \| 1.70; 1.44 \| 1.66 (**0.013 \| 0.001**) | −0.428 \| −0.748; −0.276 \| −0.581 | 1.66 \| 1.70; 1.49 \| 1.60 | registered "repel"; really **the two gaps dissolve** (below) |

Reading the table:
- ρ = ⟨D⟩/⟨D⟩_ref, with D in quanta. The reference is counting alone.
- W_min is the share of the state in the bare composition: exactly two
  gaps, nothing else.
- "Random" means c off; "coherent" means c on. For the random process,
  B is measured against counting (PREREG §3.2, amended before any run).

**Candidate roles over the six sets where gaps attract.** No candidate
is present in all six. There are two separate routes:
1. **d in the random process:** clustering.
2. **b with c in a saturated medium:** mobility binding. d strengthens
   it, from −4t to −5.12t.

Three things never attract:
- **c alone:** coherence with ordinary re-forming repels, in all four
  a/d settings.
- **d in a coherent medium:** acd and cd repel.
- **b without c:** b has no static effect in the random process, which
  is checked exactly below.

## 2. The profiles (D in quanta, L = 16; reported, not scored)

| set | V_gg(D): D = 0, 2, 4, 6, … 14 |
|---|---|
| d, bd, ad, abd | −1.38, then +0.62 flat to the end: contact only |
| d, bd: gap with excess | +1.88, then −0.12 flat: unlike repel at contact |
| abc | −1.32 at D = 0 and at D = 4; no weight anywhere else. The breathing pair |
| abcd | −1.54 at D = 0, −1.04 at D = 4 |
| ac, c, cd (L = 12) | +2.10 → −0.61 (ac), +0.98 → −0.38 (c), +0.30 → −0.19 (cd): smooth repulsion across the ring |

## 3. Mach structure

This is read from the move graph and depends only on a and b.

| a, b | isolated gap moves | adjacent pair translates | far pair can join | fragment sizes (L = 16) |
|---|---|---|---|---|
| saturated, ordinary | yes | yes | yes | whole sector |
| **saturated, tweaked** | **no** (fragment 1) | **no** (fragment 2: it breathes) | **no** (fragment 1) | exact immobility |
| headroom, ordinary | yes | yes | yes | whole sector |
| **headroom, tweaked** | **no** | **no** | **no** | dressed fragments (309k, 271k) |

**With tweaked links, a gap's position is conserved.** With one gap
and every other cell full, P = const − x, so the moment conservation
pins the gap at x.
- With headroom, the gap can dress itself in gap–excess dipoles, but it
  still cannot move by one cell.
- That is the fracton rule, emerging here from the link tweak alone.

## 4. With headroom, coherent tweaked links dissolve the gaps (bc, bcd)

- The sector energies bind: B_gg = −0.87t, B_ee = −0.43 → −0.75t.
- **But the bare two-gap composition is 1.3 % of the state at
  L = 12, and 0.1 % at L = 16.**
- The medium becomes a fluctuating liquid of dipoles. "Two gaps" stop
  being two things, and the bare sliver that remains shows them apart
  (ρ ≈ 1.6–1.7).
- The registered verdict is "repel". The honest reading is that **no
  two-gap object survives** in that medium. This echoes lane b: a
  mobile medium dissolves whole voids.

## 5. What I expected (PREREG §5, written before running), against what happened

| expected | outcome |
|---|---|
| d makes like deviations cluster in the random process; unlike repel | **held**, exactly, at contact only |
| b with c binds by mobility: isolated gap frozen, the pair breathes between contact and three cells, B = −4t | **held**, exactly (−4.000t at both sizes) |
| b off, d off gives no attraction; a on with c on is free fermions, ρ > 1 | **held** (ρ = 1.43; gate T0.3 exact) |
| ending U, carried by d in the random process and by b with c in the coherent one | **half held.** U comes from d alone, and it is like-with-like, not universal. The b-with-c route is real but cannot be sign-tested where it holds, and where it could be tested (headroom) the gaps dissolve |

## 6. Checks

- **The c-off pairs agree to every digit,** across b on and off, in all
  four pairs and at both sizes, as required by construction.
- Worst residual 2.8e-10. No lowest manifold was truncated. Moment
  violations: 0. All 32 runs are present, and every .err file is empty.
- **Near-degenerate lowest levels in the coherent tweaked-headroom
  sectors:** the next state is 0.00065–0.018t above (smallest in bcd). The registered
  state is the exact lowest manifold. A finite-temperature reading would
  mix these in.

## 7. What lane c now says

1. **There are two kinds of sign-free attraction between gaps in one
   dimension, and both act at contact:**
   - clustering, from a level that charges only for boundaries;
   - mobility binding, from tweaked links with inertia.
2. **Your link tweak produces exactly the structure of fracton gravity
   (Pretko).**
   - A lone gap is immobile.
   - Gaps move only together.
   - Their mutual mobility binds them.
   - It does so whether or not the level is flat.

   In Pretko's analysis, the force becomes a Newtonian power law only
   under certain conditions and in more than one dimension. One
   dimension cannot show that.
3. **The like-with-like clustering is not gravity.**
   - Gravity would attract a gap to an excess too, and here the gap and
     the excess repel.
   - It is the inverse of the electric sign rule, not the absence of a
     sign.
4. **Coherence with ordinary re-forming only repels.** It is the one
   extended profile in the whole grid.

## 8. Limits

- **One dimension and small rings.** No range test is possible.
- **"The two gaps" is a conditioned composition.** It is reliable in
  saturated media (W_min = 1) and thin in coherent headroom media.
- **The unlike leg in coherent states** mixes force with creation.
- **The medium is still a ring of cells,** never empty and counted in
  quanta. That is frame-faithful in 1-D. In more dimensions, the medium
  must hold its own geometry (`../../FRAME.md`).

## Then

- **c3: tweaked links in two dimensions,** where the fracton mobility
  structure can carry a power law.
  - It needs a medium of never-empty cells, re-forming that conserves
    the first moment, and distance counted in quanta.
  - Measure whether the mobility binding reaches beyond contact, and
    how it falls off.
  - Its sign test: gaps against excesses with a little headroom, but
    not so much that the gaps dissolve (§4). The headroom is the
    registered cap, not a scan.
- **The 3-D geometry gate** (R14) before any 1/r test.
