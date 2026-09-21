# Phase 0 — decision

**Written:** 2026-09-21, after the run in `cas/out/profile.out`
and the algebra in `cas/out/virial.out`.
**Status:** this page closes phase 0. It does not open a mesh.
It does not edit a kernel.

Gates were locked in `PREREG.md` before the numbers. G2 was
written in `SPIN.md` before the profile existed. The calculation
is `PROFILE.md`.

---

## 1. What was asked

Two gates, both required (`PREREG.md` §5):

| gate | pass | fail |
|---|---|---|
| G1 | half-radius and double-radius both cost ≥ 5%, meters clean | a clean measurement under 5% |
| G2 | one odd spin level from the two-atom credit | spin \(1/2\) is an import |

Either failure stops the version before phase-1 code. A dirty
meter is inconclusive, which also does not open phase 1.

---

## 2. What came back

**G1 passes.** The hedgehog sector of the locked functional has
one minimizer. The meters are clean. The wall is the scaling
wall, not a 5% scrape.

| quantity | value |
|---|---:|
| \(E\) | 72.923 |
| \(E_2\) | 36.460 |
| \(E_4\) | 36.463 |
| \(E(1/2)/E\) | 1.25003 |
| \(E(2)/E\) | 1.25001 |
| virial \(\lvert E_2-E_4\rvert/E\) | \(3.4\times 10^{-5}\) |
| seed agreement | 0 |
| refinement, mid vs fine | \(2.7\times 10^{-5}\) |
| textbook \(\lvert E-72.9\rvert/72.9\) | \(3.2\times 10^{-4}\) |

Maxima: \(dE/d\lambda = E_2 - E_4/\lambda^2\), and both ratios are
\(5/4\) when \(E_2 = E_4\). The script re-integrated the scaled
profiles; it did not insert \(5/4\).

Half-fall of the profile, \(F=\pi/2\), sits at \(x = 0.876\) in
the unit \(R_0 = \sqrt{B/A}\).

**G2 fails.** `SPIN.md` §5. No odd level is fixed. A \(2\pi\)
rotation is a loop the double cover permits to act as \(-1\).
The two-atom credit books a clock cycle. It does not choose
that sign.

---

## 3. Instrument, stated once

`PREREG.md` named a softmax chart on monotone decrements.
L-BFGS on that chart, 200 nodes, box 20, stops at iteration
4000 with \(E = 72.942\) and \(\lVert\nabla_y E\rVert_\infty =
1.4\times 10^{-3}\) (`CHART` line). A constant shift of \(y\)
does not change the profile, and the tail decrements are
exponentially small. The chart does not reach the
stationarity line.

The reported minimizer is Newton on the same midpoint energy,
in the node values \(F\), with \(F(0)=\pi\) and \(F(X)=0\).
Same grids, same seeds, same functional, same numerical
thresholds. Stationarity is \(\lVert\partial E/\partial F\rVert_\infty
= 6\times 10^{-13}\).

A reader who counts the chart substitution as a dirty meter
should call G1 inconclusive. The measured ratios are still
1.250. That reading and the pass reading have the same
consequence under §5, because G2 fails either way.

The inertia integral on this box is 52.21. A 2024 evaluation
of the same integral is 53.38. The box wall bends the
\(1/x^2\) tail after \(x \approx 20\) (\(x^2 F\) falls from
2.15 toward 0). Energy does not feel that bend. Inertia does.
Inertia was not a gate.

---

## 4. Which ending this is

Phase 1 is not opened.

The advance rule is a conjunction. G1 is a pass on the
Newton minimizer. G2 is a fail. One fail is enough.

The classical fact is not quietly promoted into a mesh.
Inside the hedgehog ansatz the corner term does what the
plan asked: one stationary radius, and a 25% wall at half
and at double that radius. That fact is on record. It is
not an authorization.

---

## 5. What is now true, and what is not

**True:**

1. The functional in `PREREG.md` §1 has a hedgehog minimizer
   at \(E = 72.923\) in units \(f_\pi/e\).
2. At that point \(E_2 = E_4\), and a radial rescaling by
   2 or by \(1/2\) costs a quarter of the mass.
3. Two seeds and three grids agree past the meter floors.
4. The two-atom credit does not fix a spin level. Spin
   \(1/2\) of this lock is an import if it is asserted.
5. Handedness remains the sign of the degree.

**Not true:**

- Phase 1 is authorized.
- A mesh soliton has been shown. The ansatz is a sector,
  not the space of maps.
- \((E_\Delta - E_N)/E\) has been predicted. The rotor
  still carries the dimensionless Skyrme parameter.
- Equal quanta on a grain have been tested. That was the
  mesh's question.
- \(n = 6\), or any cell count, was read off this profile.
- A hedgehog can be the lock. User ruling, §6: the ansatz is
  incompatible with the theory, across the experiments. The
  wall does not lift that.

---

## 6. Hedgehog, refused

**User, 2026-09-21:** any hedgehog is incompatible with this
theory. True across the experiments.

The phase-0 minimizer is a hedgehog: \(F = F(r)\), a radial
profile on an origin. `v89/PRINCIPLE.md` §5 retires profiles.
A shape \(f(r)\) is a shape in a place. The 25% wall is real
inside that ansatz and lives in the construction the principle
retires.

The post-v89 record of the same object is the same refusal in
an experiment. `v89/KIMI_REVIEW_2026-08-02.md`: after relaxation
the index read \(0.0000\) on every shell, for every radius, and
the hedgehog was withdrawn with the rest of that wave. It had
also been built on a fixed lattice, which §4.1 forbids. v100
did not find a different object. It found the wall of the one
already set aside.

The split offered below this paragraph, before this ruling,
would have meshed this functional. That split is closed. A
hedgehog seed, a hedgehog relaxation, and a hedgehog score are
not a way past this page.

## 7. How to move past this page

To reverse the hedgehog refusal: the user withdraws the
sentence in §6. G1 and G2 do not withdraw it.

To reverse G1: show the Newton minimizer fails a meter that
was locked, or exhibit another critical point of this
functional whose clean wall is under 5%. That reversal would
still leave §6 in force.

To reverse G2: write the missing derivation. `SPIN.md` §4
lists the three readings that do not count. That reversal
would still leave §6 in force.
