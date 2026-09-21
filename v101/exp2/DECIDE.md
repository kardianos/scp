# Experiment 2 — decision

**Written:** 2026-09-21, after `cas/out/octave.out` and
`cas/out/monodromy.out`.
**Status:** this page closes the octave test. It does not edit
experiment 1. It does not edit a kernel.

Gates were locked in `PREREG.md` before either number existed.
Experiment 1's decision was already written.

---

## 1. What was asked

G2A: a \(1{:}1\) cycle returns at \(2\pi\), and the octave cycle
\((2{:}1,\ 1{:}2)\) returns at \(4\pi\) on \(\psi = 2\theta\).

G2B: the unison gate chain reproduces \(E = 7.432761\).

G2C: the \(2{:}1\) interval chain, gate variable
\(\Phi_i = \theta_i - 2\theta_{i-1}\), puts a localized core on the
closure boundary \(\Delta\theta = 2\pi\) at that same cost, within
5%. The odd boundary \(\Delta\theta = \pi\) must not be a
half-price copy of that core. A unison ring at twist \(2\pi/3\) is
a control, not a score.

## 2. Algebra

Maxima and the script agree.

| class | shift |
|---|---|
| \(1{:}1\) generator | \(2\pi\) |
| octave odd, \(n=(1,0)\) | \(\pi\) in \(\theta\), \(2\pi\) in \(\psi\) |
| octave even, \(n=(1,1)\) | \(2\pi\) in \(\theta\), \(4\pi\) in \(\psi\) |

**G2A passes.**

## 3. Instrument

The unison open chain, three seeds, lowest kept:

\[
E = 7.432761, \qquad \lVert\nabla E\rVert_\infty = 9.9\times 10^{-11}.
\]

Relative miss against 7.432761 is \(4.6\times 10^{-8}\). **G2B is
clean.**

The control ring at twist \(2\pi/3\) relaxes to \(E = 1.510164\)
with 2 non-vacuum sites. That is the withdrawn apparatus
(`CHARGE.md` §8.4), reproduced, and it is not a fractional core.

The analytic interval gradient matches a central difference to
relative \(8.3\times 10^{-8}\).

One repair, recorded so the floor is readable. The first descent
treated a \(10^{-18}\) wobble in \(E\) as a real rise and stopped
at \(\lVert\nabla E\rVert_\infty \approx 10^{-7}\). "Raises \(E\)"
is now a rise above \(10^{-12}\), which is still far under the
\(10^{-3}\) and 5% gates. Those gates were not moved. After the
repair the unison gradient is \(10^{-10}\).

## 4. The interval chain

| boundary | \(E\) | localized | \(\max\lvert\Delta\theta\rvert / \mathrm{mean}\) |
|---|---:|---|---:|
| vacuum \(\Delta\theta = 0\) | 0 | — | — |
| odd \(\Delta\theta = \pi\) | 1.601130 | yes | 191.5 |
| closure \(\Delta\theta = 2\pi\) | 6.225890 | yes | 185.2 |

The closure is a localized core. Its energy is 16.2% below
7.432761. The window was 5%. **That is a miss.**

The odd state also localizes. Half the closure would be 3.113;
half the unison core would be 3.716. The odd energy is 1.601,
outside 10% of either half. The chain does not restore the
withdrawn "half the integer" pattern. It also does not make the
odd class empty. Both boundaries hold a core, at two costs, and
neither cost is the unison number.

**G2C fails. Experiment 2 fails.**

## 5. What is now true, and what is not

**True:**

1. The corrected composition gives a \(4\pi\) return of the octave
   phase and a \(2\pi\) return of a \(1{:}1\) cycle.
2. The unison integer core and the unison fractional artifact both
   reproduce.
3. On \(\Phi = \theta - 2\theta_{\mathrm{prev}}\), a \(2\pi\) end
   shift is a localized core at \(E = 6.226\), not at \(7.433\).

**Not true:**

- The integer slip on this interval chain costs the gate's core
  energy.
- A \(4\pi\) return was measured as a dynamical period. G2A is the
  monodromy. The chain imposed the boundary; it did not discover it.
- The odd class is a fractional charge of one half.

## 6. How to disagree

To reverse G2A: show the corrected composition produces a different
group. The Maxima log is the witness.

To reverse G2C: on this functional, with these seeds, bring the
closure energy inside 5% of 7.432761 without changing \(A\), \(p\),
or the definition of \(\Phi\). Rescaling \(A\) to force the match
is the ratio this experiment does not add.
