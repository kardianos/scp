# Experiment 2 — pre-registration

**Written:** 2026-09-21, before the algebra and before the chain.
**Status:** commitment. Thresholds here are not revised after the
numbers exist. `DECIDE.md` applies them.

Experiment 1's counts are not an input. This file does not name
them.

---

## 1. Algebra (locked)

Corrected composition, \(\alpha = 0\):

\[
c \mapsto (p\, c + 2\pi n)/q, \qquad n\in\mathbb{Z}.
\]

Start at \(c = 0\).

- A pure \(1{:}1\) cycle (one or more edges \((p,q)=(1,1)\)) has
  monodromy in \(2\pi\mathbb{Z}\). The smallest positive element
  is \(2\pi\).
- The octave cycle \((2{:}1,\ 1{:}2)\) has monodromy
  \(\pi(n_1+n_2)\), so the group is \(\pi\mathbb{Z}\).
- The octave phase is \(\psi = 2\theta\). The odd class
  \(n_1+n_2\) odd shifts \(\psi\) by \(2\pi\) and is not a closure.
  The smallest positive even class shifts \(\theta\) by \(2\pi\) and
  \(\psi\) by \(4\pi\). That \(4\pi\) is the return this experiment
  means.

Maxima (`cas/monodromy.mac`) is the witness. The Python gate reads
the same composition.

**G2A pass** only if all three hold: \(1{:}1\) minimum \(2\pi\),
octave group \(\pi\mathbb{Z}\), even-class \(\psi\)-return \(4\pi\).
**G2A fail** if any one differs. No chain is required for G2A.

## 2. Unison instrument (locked)

The standing gate, copied from the repaired Frenkel–Kontorova
chain (`v89/charge_fk.c`, `v89/CHARGE.md` §8.1). Not a new functional.

\[
V(\theta) = 1 - \big((1+\cos\theta)/2\big)^{8},
\qquad
E = \sum_i V(\theta_i) + \tfrac{A}{2}\sum_i (\theta_{i+1}-\theta_i)^2
\]

with \(A = 1\), open chain, \(N = 400\), \(\theta_0 = 0\),
\(\theta_{N-1} = 2\pi\).

Overdamped descent. Initial step \(0.02\), halved when a step
raises \(E\), cap \(4\times 10^{5}\) steps, three seeds (uniform
ramp, site-centred step, site-centred step shifted by one), lowest
\(E\) kept.

**G2B clean** if \(\lvert E/7.432761 - 1\rvert \le 10^{-3}\) and
\(\lVert\nabla E\rVert_\infty \le 10^{-8}\) on the free sites.
**Dirty** otherwise. A dirty G2B makes G2C `INCONCLUSIVE`. It does
not flip G2A.

## 3. Interval chain (locked)

Same \(N\), \(p=8\), \(A=1\), same descent, same three seeds.
The gate variable is the comb phase of a \(2{:}1\) bond:

\[
\Phi_i = \theta_i - 2\,\theta_{i-1},
\qquad
E = \sum_{i=1}^{N-1} V(\Phi_i)
+ \tfrac{A}{2}\sum_{i=0}^{N-2}(\theta_{i+1}-\theta_i)^2.
\]

Ends pinned. The angle \(\theta\) is the monodromy angle of §1.

| boundary | \(\theta_{N-1}\) | what it is |
|---|---:|---|
| vacuum | \(0\) | \(\Phi = 0\) |
| odd | \(\pi\) | the odd octave class |
| closure | \(2\pi\) | even class, the \(4\pi\) return of \(\psi=2\theta\) |

Analytic \(\partial E/\partial\theta\) matches a central difference
to relative \(10^{-6}\) on one random interior step, or the chain
is dirty.

A state is **localized** when both hold: the fraction of bonds with
\(V(\Phi) > 0.01\) is at most \(0.05\), and
\(\max\lvert\Delta\theta\rvert \ge 3\,\mathrm{mean}\lvert\Delta\theta\rvert\).

**G2C pass** only if the closure state is localized and

\[
\big\lvert E(2\pi)\,/\,7.432761 - 1\big\rvert \le 0.05,
\]

and the odd state is not a half-price copy: it fails the
localization test, or it localizes and
\(\lvert E(\pi)/(E(2\pi)/2) - 1\rvert > 0.10\).

**G2C fail** if the meter is clean and any of those is false.

**G2C inconclusive** if G2B is dirty, the gradient check is dirty,
a reported state has \(\lVert\nabla E\rVert_\infty > 10^{-6}\), or
the vacuum energy exceeds \(10^{-6}\).

## 4. Negative control (not a score)

Periodic unison ring, \(N=400\), twist \(2\pi/3\), same potential
as §2. This is the withdrawn apparatus. It is printed so the run
shows it can see that artifact, and then it is not graded as G2C.

**Control clean** if non-vacuum sites \(\le 10\). Dirty control
makes G2C inconclusive. The control's energy is not a fractional
core.

## 5. Verdict

| | |
|---|---|
| **PASS** | G2A pass and G2C pass |
| **FAIL** | G2A fail, or G2C fail on a clean meter |
| **INCONCLUSIVE** | G2C inconclusive, and G2A did not fail |

G2A fail is a fail even if the chain was inconclusive.

## 6. What this experiment does not claim

A pass is an octave return and a core cost. It is not a radius,
not a hedgehog, and not a fractional charge of \(1/3\). A fail
closes this chain. It does not reopen experiment 1.
