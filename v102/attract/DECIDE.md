# v102 attract — decision (R8; decisions 2–3 by test; the search for an attraction)

**Written:** 2026-10-08, after the outputs in `cas/out/` (the analysis
is `cas/out/analyze_attract.out`). The rules and the endings are in
`PREREG.md`. Runs beyond the registration are marked *supplementary*.
**Status:** R8 is the user's ruling. Decisions 2 and 3 are settled here
by the registered rules. Item 4 is a measurement. No kernel or law
table was edited.

---

## 1. The answer

**Decision 2 → T-floor.**

- Transport inside a lock is allowed, with a bounded amplitude.
- Transport may not empty a lock cell.
- The measured window: t = 0.1 J0 passes and t = 1 J0 fails. So
  t_max lies between 0.1 and 1 J0, about 0.001–0.01 in rate units.

**Decision 3 → S-boundary.** Space may leave a lock only through a
link with one end outside it.

**Item 4: ending LN.** The fabric holds two objects together at a
distance, in two ways, but **not at the same transport strength**.

- **A nucleus (small t).** Two open slots bind into a pair of fixed
  size in the channel that is a singlet in everything:
  - W: (S, I) = (0, 0), B = −0.50 J0, mean separation 3.1 grains, the
    same at 16 and 24 cells;
  - 3+2: S = 0, B = −0.13 J0, separation 3.0, the same at 24 and 40.
  - Every other channel is two free holes.
- **A ladder (large t).** At t = 7 J0 the hole binds one lock quantum
  into an S = 3/2 (3+2) or (3/2, 3/2) (W) level below the threshold:
  −0.07 J0 (3+2, registered) and −0.06 J0 (W, supplementary),
  extrapolated. The ground state stays the hole's.
- At the t that keeps the hole firm and the pair bound (≲ 0.1 J0),
  there is no ladder. At the t that makes the ladder, the pairs come
  apart (supplementary):
  - 3+2 at t = 3, 5, 7 J0: B ≥ −0.005 J0, with a separation that grows
    with the chain;
  - W at 7 J0: B = +0.47 J0, repulsive.

**Where the attraction comes from.** No term was added; it is the
fiber's own correlations.

- Two holes on opposite sublattices couple through the gapped
  background. In the all-singlet channel they attract.
- Static, d = 3: −0.19 J0 (3+2), −0.75 J0 (W).
- The attraction falls off over about the correlation length, and it
  is close to zero for holes on the same sublattice.

## 2. Decision 2: transport inside the lock (B1, B2)

**B2: do two holes stay two?**

| option | two holes | source |
|---|---|---|
| T-none | pinned; the hole cannot walk | P4 |
| T-bound (a cell may empty) | **collapse** into one empty cell that cuts the fiber (−20.4 J in 3+2, −16.5 J in W) | `../perm/` X5, microscopic |
| **T-floor** | stay two, and bind in the all-singlet channel (§5) | the leading-order model, whose local space holds N = 2 or 1 only, so the floor is exact |

**B1: the moving W hole** (TeNPy, χ = 128):

| t | ground | profile, 16 → 24 → 32 cells | registered decay length (16 → 24) | 24 → 32 (supplementary) | 90% radius |
|---|---|---|---|---|---|
| 0.1 J0 | (½, ½) | identical to 3 decimals | 2.463 → 2.066 (16%) **✘** | 2.066 → 1.991 (3.6%) **✔** | 5, 5, 5 |
| 1 J0 | (½, ½) | identical; no longer staggered | 2.489 → 1.347 (46%) ✘ | — | 3, 3 |
| 3 J0 | (½, ½) | identical; a uniform polarization | 1.963 → 1.633 (17%) ✘ | — | 4, 5 |

**Verdict by the rule.**

- T-bound is out: holes collapse and cut the fiber.
- T-floor keeps the (½, ½) ground state and lets the hole walk (band
  ≈ 0.4 J0 at t = 0.1 J0, `../perm/` X3).
- Its N2 passes at t = 0.1 J0 on the two longest chains run (24 and 32;
  the registered pair 16/24 misses, as the static W hole's did, §8.1).
  It fails at 1 J0 and above.
- The 3+2 pair binds at 0.1 J0 and loosens at 1 J0 (§5).

**So: T-floor, t < t_max, with t_max between 0.1 and 1 J0.** It adds
one statement (the floor) and one bounded constant.

## 3. Decision 3: space leaving through the boundary

| condition | measured |
|---|---|
| the only formation route found | yes: `../fabric/` F2, `../perm/` X4. No other process makes a hole |
| a lock in cold light is untouched | yes: `../options/` E5 (no transitions, no destruction) |
| the holes it makes do not cut the fiber | under T-floor, yes (§2) |

**→ S-boundary.** Light makes a W hole with about 5% of resonant
photons, and the photon's (S^z, I^z) becomes the hole's.

## 4. A ladder (A1, A2, A4; W at t = 7 is supplementary)

The margin is E(next multiplet) − E(hole) − gap(L), extrapolated in
1/(L+1)².

| case | margins (shorter, longer) | extrapolated | verdict |
|---|---|---:|---|
| 3+2, t = 0.1, 1, 3, 5 J0 | +0.03 … +0.06 | +0.012 to +0.015 | not bound |
| **3+2, t = 7 J0** | −0.168, −0.107 | **−0.071** | **bound; ground S = ½** |
| W, t = 0.1 J0, (3/2, 3/2) | +0.055, +0.024 | −0.002 | not bound (bar −0.01) |
| W, t = 1 J0 | +0.042, +0.040 | +0.039 | not bound |
| W, t = 3 J0 | −0.014, −0.0005 | +0.011 | not bound |
| **W, t = 7 J0 (supplementary)** | −0.107, −0.084 | **−0.064** | **bound; ground (½, ½)** |
| 3+2 all orders, g/Δ = 0.1, 0.2, 0.3 (A2) | +0.10/+0.05, +0.10/+0.05, +0.09/+0.04 | +0.025, +0.022, +0.017 | not bound |
| W (3/2, ½) against (3/2, 3/2) (A4) | degenerate to 1e-6 at 24 and 32 | — | no lower channel |

- **The ladder comes from transport, not from the fourth-order terms.**
- At t = 7 J0 the moving hole binds one lock quantum: a spin(–isospin)
  polaron.
- This is the regime approaching the Nagaoka onset (9–10 J in 3+2). It
  is still inside the lock's transport margin, t/Δ ≈ 0.14 < 0.17.

## 5. A nucleus (A3; t ≥ 3 J0 supplementary)

**Static pairs, against d = 10** (J0):

| d | 3+2, S = 0 | 3+2, S = 1 | W (0,0) | W (1,0)/(1,1) |
|---:|---:|---:|---:|---:|
| 2 | +0.005 | +0.005 | −0.044 | −0.044 |
| **3** | **−0.188** | +0.168 | **−0.750** | +0.180 |
| 4 | −0.013 | −0.013 | −0.069 | −0.069 |
| 5 | **−0.124** | +0.113 | | |
| 6 | −0.012 | −0.012 | −0.043 | −0.043 |

**Mobile pairs (T-floor)**, B = E₂ + E₀ − 2E₁:

| fabric | t | channel | B (shorter / longer) | mean separation | verdict |
|---|---|---|---|---|---|
| 3+2 | 0.1 J0 | S = 0 | −0.130 / −0.129 | 2.99 / 2.99 (P(3) = 0.96) | **bound** |
| 3+2 | 0.1 J0 | S = 1 | +0.004 / +0.001 | 9.7 / 13.8 (grows) | free |
| 3+2 | 1 J0 | S = 0 | −0.073 / −0.064 | 6.9 / 8.0 (grows) | loosening, ✘ by the bar |
| **W** | 0.1 J0 | **(0, 0)** | **−0.504 / −0.503** | **3.06 / 3.06** | **bound** |
| W | 0.1 J0 | (1, 0), (1, 1) | +0.023 / +0.000 … +0.027 | 6.6 → 7.7 (grows) | free |
| 3+2 (supp.) | 3 J0 | S = 0 | +0.040 / −0.005 | 9.0 / 12.6 (grows) | free |
| 3+2 (supp.) | 5, 7 J0 | S = 0 | +0.14 / +0.06; +0.16 / +0.05 | 9.9 → 16.1; 10.4 → 16.9 | free (repulsive) |
| W (supp.) | 7 J0 | (0, 0) | +0.474 (16 cells) | 6.3 | repulsive |

- **The nucleus is a pair of open slots** about three grains apart,
  bound by the background's correlations. It is not two cells joined
  by a rule.
- In W it binds only as a double singlet. Nothing is mapped (the
  deuteron's channel is a different one; recorded, not read).

## 6. What is now true, and what is not

**True (measured in models):**

1. Under R8 + T-floor + S-boundary, one open slot is a firm (½, ½)
   object. That is p and n, both spin ½, one object at one size. It
   walks for t ≲ 0.1 J0, and light makes it.
2. **Two of them bind into a pair of fixed size** in the all-singlet
   channel, B ≈ −0.5 J0 (W). The force is the fiber's own, not an
   added term.
3. **A bound excited level (a ladder) exists only at large transport,**
   t ≈ 7 J0, where the hole binds a lock quantum. At that t, two holes
   no longer bind (supplementary, §5). **The ladder and the nucleus do
   not coexist at any t tried.**
4. The fourth-order terms at larger g do not make a ladder.
5. Transport without the floor collapses hole pairs into fiber cuts.

**Not true:**

- A ladder at the transport where the nucleon keeps its shape.
- Any mapping of these energies, sizes or channels onto measured nuclei.
- Any of it in a kernel.

## 7. What would move past this page (each is the user's call)

1. **Rule decisions 2 and 3 as settled here,** or overrule them:
   T-floor with t < t_max (0.1–1 J0), and S-boundary.
2. **The t tension.** The nucleon's shape and the bound pair want small
   t, and the ladder wants t ≈ 7 J0. Options:
   - accept no ladder at the nucleon's t;
   - look for a ladder source that works at small t: W at all orders,
     larger g with transport, or two-lock-quantum states;
   - or find whether the pair survives at the ladder's t (§5's
     supplementary rows).
3. **A kernel step** to put R8 + T-floor + S-boundary into a real
   simulation. That needs authorization.

## 8. Deviations from `PREREG.md`

1. **N2 for the moving W hole.** The registered lengths (16, 24) miss
   by 16%, while the profiles are identical. A supplementary 32-cell run
   gives 24 → 32 = 3.6% (pass). The same thing happened to the static W
   hole in `../perm/`: 16 → 24 = 7.5%, 24 → 32 = 2.8%. Both numbers are
   reported, and the decision uses the two longest chains run.
2. **Supplementary runs**, not in the registration:
   - W at t = 7 J0 (the ladder);
   - pairs at t = 3, 5, 7 J0 (3+2) and t = 7 J0 (W);
   - the moving W hole at L_bulk = 32.

   Each was added after seeing 3+2's ladder, or the moving hole's 16/24
   miss, and is labelled as such.
3. **Convergence.**
   - The mobile W excited states run at χ = 128 with truncation error
     up to 1e-3. The (3/2, ½) sector came out above (3/2, 3/2) at
     mobile W, which is impossible exactly, so it is unconverged by up to
     0.01–0.05 J0. The (3, 1) ladder rows are not used for verdicts.
   - The static W pair sectors disagree at the 1e-3 J0 level.
   - The 3+2 runs converge to 1e-5 or better.
4. **The all-order runs (A2)** use the in-tree compact-MPO DMRG. At
   g/Δ = 0.1 they reproduce `../fabric/` F4b exactly (−62.6553632054).

## Reproduce

From `v102/attract/cas/`. `PY` is `~/.venvs/scp-tenpy/bin/python`.

| step | command | output |
|---|---|---|
| registered batch | `./run_attract.sh` | `out/tp_*`, `out/a2_*` |
| supplementary | `./run_supp.sh`, `./run_supp2.sh`; the 32-cell moving W hole by `$PY -I ../../perm/cas/tp_chain.py 4 mob 32 1,1 0.1 128` | `out/tp_*t7.0*`, `out/tp_M4_mob_32_*` |
| analysis | `$PY -I analyze_attract.py` | `out/analyze_attract.out` |

Instruments:

- `../../perm/cas/tp_chain.py`: TeNPy chains. It now takes pairs
  (`pairD`) and two mobile holes (`mob2`), and reports the pair distance
  distribution.
- `a2_allorder.py`: 3+2 at all orders, at any g.
