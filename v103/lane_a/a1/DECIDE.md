# v103 lane a, step a1 — decision (the void as a quantum particle)

**Written:** 2026-10-11, after the registered runs. Gates: `PREREG.md`.
A measurement in models; it rules nothing.

- Tool: `../../tools/bx.c` (C), passing T0 at 1, 4, 8 and 16 threads.
- Runs: `cas/run_a1.sh`, 16 runs (8 sets × L = 6, 8). The last one, the
  2.2 × 10⁶-state sector, took 6 h 45 min.
- Tables: `cas/tab.go` (Go), output in `cas/out/tab.out`.

**Frame caveat (R14, ruled after this step was registered).**
- a1's medium is a fixed ring of cells, and a cell may be emptied: a
  double void is an empty cell. That is a container, "something in
  something".
- So a1 tests known quantum mechanics in the container frame (see
  `../../FRAME.md`, and §5 below for what the frame changes).

## Ending: **P(2)**

- **By the registered table: P(3), best c.**
- **But c's breathing-line pass is a defect in the bar.** Read as the
  bar states its intent, c fails it, and the ending is P(2). §2 shows
  this.
- **The decision is P(2).** Four sets tie at 2:
  - none, b and c (a wave and an antiparticle);
  - abc (a wave and spin ½).

| set | A1 wave | A2 exclusion | A4 spin ½ | A5 antiparticle | A6 breathing line | count |
|---|---|---|---|---|---|---|
| none | ✔ | ✘ | ✘ | ✔ | ✘ | 2 |
| a | ✔ | ✘ | ✘ | ✘ (no excess) | ✘ (cannot breathe) | 1 |
| b | ✔ | ✘ | ✘ (ferromagnetic medium) | ✔ | ✘ | 2 |
| c | ✔ | ✘ | ✘ | ✔ | ✘ (registered ✔; bar defect, §2) | 2 (registered 3) |
| ab | ✘ (degenerate medium) | ✘ | ✘ | ✘ | ✘ | 0 |
| ac | ✔ | ✘ | ✘ | ✘ | ✘ | 1 |
| bc | ✘ (degenerate, clumped medium) | ✘ | ✘ (medium gap 0 at L = 8) | ✘ | ✘ | 0 |
| **abc** | ✔ | ✘ (L = 8) | **✔** | ✘ | ✘ | 2 |

A property survives only if its bar holds at both L = 6 and L = 8.

**No set carries more than two of the five properties.** The tension
is the one predicted before running:
- the saturation that makes a void sharp (a) also forbids the
  antiparticle and the breathing;
- spin ½ needs all three ingredients at once: saturation, two planes
  and exchange.

## 1. The measurements (L = 6 | L = 8)

| set | medium | void W, Z | void S² | excess Z | breathing E / threshold, weight | P_double | two-void binding |
|---|---|---|---|---|---|---|---|
| none | gap 3.15 \| 2.43, singlet | 4.67, 0.982 \| 4.71, 0.985 | — | 0.984 \| 0.983 | 6.28 / 0.59, 0.26 \| 4.90 / 0.44, 0.09: above the continuum | 0.49 \| 0.51 | +0.47 \| +0.37 |
| a | frozen (one state) | 8.00, 1.000 | — | none | O\|G⟩ = 0 | 0.032 \| 0.014 | +0.86 \| +0.52 |
| b | ferromagnet, S² = 42 \| 72 (the full multiplet) | as none | 35.75 \| 63.75 | 0.979 \| 0.977 | as none | as none | as none |
| c | gap 0.15 \| 0.0018, number fluctuation 7.1 \| 12.0 | 0.34, 0.588 \| 0.057, 0.550 | — | 0.853 \| 0.849 | 1.88 / **−2.40**, 0.17 \| 0.28 / **−0.017**, 0.21 | **1.88 \| 2.78** | −0.83 \| −0.06 |
| ab | degenerate (gap 0) | 0.096, 0.001 \| 0.036, 0.000 | 35.75 \| 63.75 | none | O\|G⟩ = 0 | 0.032 \| 0.014 | +0.86 \| +0.52 |
| ac | frozen | 8.00, 1.000 | — | none | O\|G⟩ = 0 | 0.053 \| 0.051 | −0.91 \| −0.76 |
| bc | gap 0.0023 \| 0 (degenerate) | 0.0004, 0.33 \| 0.0000, 0.24 | 0.75 | 0.69 \| 0.001 | 2.21 / 2.85, 0.43 \| 4.77 / 2.00, 0.001 | 1.99 \| 2.00 | −3.09 \| −2.64 |
| **abc** | gap 7.21 \| 5.94, singlet | 1.46, 0.988 \| 1.23, 0.979 | **0.750000** | none | O\|G⟩ = 0 | 0.0001 \| 0.025 | **−7.50 \| −7.31** |

- |dAB| is ≥ 1e-8 in every A1-passing set. Values run from 5e-3 to
  0.54.
- Energies are in units of t.

## 2. The A6 bar defect, disclosed

- **The bar** was R = E_breath/threshold < 1, with the threshold being
  the bottom of the continuum of one void plus one excess.
- **It assumed a positive threshold.** In set c the threshold is
  negative (−2.40 and −0.017): E(N+1) + E(N−1) < 2E(N).
  - The energy is concave in N.
  - The pair-exchange medium wants to clump. Its number fluctuation of
    7–12 and its P_double of 1.9–2.8 empty cells show this.
- **With a negative threshold, R < 1 holds for any positive breathing
  energy, so the bar passes by sign alone.**
- **The bar's stated intent is a breathing line below the continuum.**
  For c the breathing pole sits *above* it: 1.88 > −2.40 and
  0.28 > −0.017. c fails.
- This reading is stricter, not softer, so it follows the ratchet.
  Future breathing bars compare E_breath with the threshold directly.

## 3. What I expected (PREREG §5, written before running), against what happened

| expected | outcome |
|---|---|
| saturation (a) gives sharp voids, but breaks the antiparticle and the breathing | **held**: W = 8, Z = 1; no excess; O\|G⟩ = 0 |
| headroom allows an excess, but softens the void | **half held.** The excess exists (Z ≈ 0.98), but the void stays sharp (Z ≈ 0.98) unless exchange is added (c: Z = 0.55, W → 0.06) |
| spin ½ needs b and c | **held, and it also needs a.** Only abc carries spin ½. Without saturation, bc's medium is degenerate and clumped |
| exclusion is not given by saturation at n̄ = 2, because a cell can hold zero | **held.** P_double is 0.014–0.05 under a, and no set passes A2 at both sizes |
| P(3) or less | **held**: P(2) |

**Unexpected:**
- **Two planes without exchange (b) make the medium a ferromagnet.**
  Its S² is the full multiplet, which is the known SU(2)-boson theorem.
  So the void cannot carry spin ½ against it.
- **Exchange without saturation (c, bc) makes the medium clump.** It is
  then no longer a uniform medium.
- **In abc, two voids bind hard:** −7.3t at L = 8, against a void
  bandwidth of 1.2t. The spin-½ voids pair.

## 4. Checks

- Every medium and void residual is ≤ 1.3e-7. The worst is bc, L = 8,
  whose medium is degenerate.
- Spectral convergence is ≤ 3.4e-13, except for bc at L = 8 (8e-6 and
  1e-4) and ab at L = 8 (2e-7). Both are degenerate media, and their
  bars fail by orders of magnitude, so the residuals do not decide
  anything.
- Every .err file is empty.

## 5. What lane a now says, and what the frame changes

1. **A void in a full medium is a good wave** (A1) in six of the eight
   sets.
   - The inversion survives that test: a missing quantum moves as a
     sharp quasiparticle, and it interferes around the ring.
2. **The other properties do not come together in one container
   medium.**
   - The antiparticle needs headroom.
   - Spin ½ needs saturation, two planes and exchange.
   - The headroom and the saturation exclude each other.
   - The breathing line never sits below the continuum.
3. **Exclusion (A2) fails everywhere here, because a cell can be
   emptied.**
   - Under R14 a cell cannot be emptied (a place never exists without
     medium). Then two gaps can never share a cell, and exclusion is
     structural: one gap per cell. That is your "exclusion principle in
     field over space".
   - c2 used exactly that rule, so in the gap frame A2 is not a property
     to be earned.
   - That is a frame statement, not a measurement here: a1 ran in the
     container frame.
4. **The spin-½ void (abc) is the a1 analogue of a hole in a singlet
   spin medium,** and its pairs bind. The nucleus-like reading, pairs of
   spin-½ voids bound well beyond their bandwidth, is left for a later
   lane-a step.

## Dimension note (2026-10-11, after the decision)

The user asked whether this was 3-D. **It was not: a1 ran on 1-D rings.**
Several of the five properties are not 1-D properties, so P(2) is the
1-D reading only.
- **A4 spin ½.** A ring has no rotations. What a1 measured is a doublet
  (S² = 3/4), not a sign flip under a 2π turn. Spin as rotation, and
  its link to statistics, need 3-D.
- **A2 exclusion (and the deferred A3).** In 1-D, hard-core bosons and
  fermions are equivalent, and quanta cannot pass one another, so
  exchange is not defined. A2 tested occupation only.
- **A1 wave.** On a long 1-D ring, an interacting medium's quasiparticle
  weight goes to zero. A void in a spin medium splits into a spin part
  and a charge part. So abc's spin-½ void at L = 6, 8 is likely a
  short-ring object in 1-D. In 3-D it can stay whole.
- **A6 and the void pairs.** In 1-D any net attraction binds. In 3-D
  binding needs a threshold. So abc's −7.3t pair is the generous case,
  and A6's failure, which happened even where binding is easiest, is
  the robust case.
- **The medium.** R6's n̄ = 2 is the ring value: two links per cell. A
  3-D medium has more links per cell, and its spin medium can order,
  where a 1-D one cannot.
- **Carried over to any dimension:** the per-cell facts. The excess
  needs headroom, and saturation forbids it.

A 3-D lane-a step needs the 3-D geometry gate first (R14: a fixed 3-D
grid would be a container). Exact enumeration cannot reach 3-D: the
smallest periodic cube without doubled bonds has 27 cells, about 10¹²
states at n̄ = 2.

## Limits

- **Rings of 6 and 8 cells.** Values near bars, such as c's Z = 0.55 and
  abc's P_double = 0.025, are finite-size readings.
- **Exchange J = J0/2 = 5t** came from R9's ratio, not from a scan. Under
  it, c and bc clump.
- **The container frame:** cells may empty, and distance is counted in
  cells.

## Then

- **a2, in the gap frame:** cells never empty, and distance counted in
  quanta.
  - Recheck A1, A4 and A5 with exclusion now structural.
  - Find whether any medium carries the antiparticle and spin ½ at
    once. That needs headroom without clumping.
- **The abc void pair:** its size and binding. It is the inversion's
  candidate for composite matter.
