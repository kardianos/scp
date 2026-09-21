# Phase 0 — the hedgehog

**Run:** 2026-09-21. Instrument `cas/profile.py`. Algebra
`cas/virial.mac`. Logs in `cas/out/`.

The functional is `PREREG.md` §1. No coefficient was scanned.

---

## 1. Algebra

A radial stretch \(F_\lambda(x) = F(x/\lambda)\) multiplies the
quadratic piece by \(\lambda\) and the corner piece by
\(1/\lambda\), up to the wall of a finite box. Maxima:

\[
\frac{dE}{d\lambda} = E_2 - \frac{E_4}{\lambda^2}.
\]

A critical point inside this family has \(E_2 = E_4\), and then

\[
\frac{E(\lambda)}{E(1)} = \frac{\lambda + \lambda^{-1}}{2}.
\]

Both \(\lambda = 1/2\) and \(\lambda = 2\) are \(5/4\). This is
an identity about scaling. It does not produce a profile.

---

## 2. Minimizer

Newton on the midpoint quadrature. Nodes \(F(0)=\pi\),
\(F(X)=0\), free in between. Seeds \(4\arctan e^{-x}\) and
\(4\arctan e^{-x/2}\), decrements rescaled to sum to \(\pi\).

The softmax chart named in the pre-registration stalls:
\(N=200\), \(X=20\), 4000 L-BFGS iterations, \(E=72.942\),
\(\lVert\nabla_y E\rVert_\infty = 1.4\times 10^{-3}\), not
successful. Production numbers below are the Newton solve.
`DECIDE.md` §3 records what that substitution does to the
meter.

Fine grid: 2000 intervals, box \(X=60\).

| | |
|---|---:|
| \(E\) | 72.923092 |
| \(E_2\) | 36.460319 |
| \(E_4\) | 36.462774 |
| \(\lVert\partial E/\partial F\rVert_\infty\) | \(6\times 10^{-13}\) |
| \(x\) where \(F=\pi/2\) | 0.876 |

Both seeds give the same \(E\) to the roundoff the script
printed (relative 0). Coarse \(X=30\) gives 72.924. Mid
\(X=60\), 1000 intervals, gives 72.921. Fine versus mid is
a relative \(2.7\times 10^{-5}\).

The profile falls. Its largest step upward is below 0.
Past \(x \approx 5\), \(x^2 F\) sits near 2.15 until the box
wall pulls it down (1.52 at \(x=40\), 0 at \(x=60\)).

---

## 3. The wall, re-integrated

\(F^\star(x/\lambda)\) was interpolated onto \([0,\lambda X]\)
and integrated again. The algebraic \(5/4\) was not written
into the ratio.

| \(\lambda\) | \(E\) | \(E/E(1)\) | \((\lambda+1/\lambda)/2\) |
|---:|---:|---:|---:|
| 0.50 | 91.156 | 1.25003 | 1.250 |
| 0.75 | 75.962 | 1.04168 | 1.04167 |
| 1.00 | 72.923 | 1.00000 | 1.000 |
| 1.25 | 74.747 | 1.02501 | 1.025 |
| 1.50 | 79.000 | 1.08334 | 1.08333 |
| 2.00 | 91.155 | 1.25001 | 1.250 |

The kill floor was 1.05. The measurement is 1.25.

---

## 4. Inertia, not a gate

\[
I = \frac{8\pi}{3}\int_0^X x^2\sin^2 F\Big(1 + F'^2 + \frac{\sin^2 F}{x^2}\Big)\,dx
= 52.211
\]

in units \(1/(e^3 f_\pi)\). A 2024 evaluation of the infinite-line
integral is 53.38. The difference tracks the wall in the tail.
Phase 0 does not turn \(I\) into a \(\Delta\)–nucleon ratio:
that ratio still contains the dimensionless parameter the
continuum unit scaled out of the profile only.

---

## 5. Reproduce

```
python3 v100/phase0/cas/profile.py
maxima --very-quiet -b v100/phase0/cas/virial.mac < /dev/null
```
