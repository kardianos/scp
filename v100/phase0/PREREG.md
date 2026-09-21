# Phase 0 — pre-registration

**Written:** 2026-09-21, before the profile run.
**Status:** commitment. Thresholds in this file are not revised
after the numbers exist. `DECIDE.md` applies them.

Phase 0 has two gates. The version proceeds to phase-1 code only if
both pass. A dirty meter is `INCONCLUSIVE`, which is not a pass.

---

## 1. The functional (locked)

Hedgehog profile \(F(0)=\pi\), \(F(\infty)=0\), energy in units
\(f_\pi/e = 1\), radius in units \(x = e f_\pi r\):

\[
\begin{aligned}
E_2[F] &= 2\pi\int_0^\infty\big(x^2 F'^2 + 2\sin^2 F\big)\,dx,\\
E_4[F] &= 2\pi\int_0^\infty\sin^2 F\Big(2 F'^2 + \frac{\sin^2 F}{x^2}\Big)\,dx,\\
E &= E_2 + E_4.
\end{aligned}
\]

No pion-mass term. No free coefficient. The split is the quadratic
piece and the corner piece of `../README.md`.

Off shell, a radial rescaling \(F_\lambda(x) = F(x/\lambda)\) obeys

\[
E_2[F_\lambda] = \lambda\, E_2[F],\qquad
E_4[F_\lambda] = \lambda^{-1}\, E_4[F],
\]

up to truncation of a finite box. At a stationary point \(E_2 = E_4\),
so \(E[F_\lambda]/E[F] = (\lambda + \lambda^{-1})/2\). In particular
both \(\lambda = 1/2\) and \(\lambda = 2\) sit at \(5/4\). Maxima
checks that algebra (`cas/virial.mac`). It does not prove a
minimizer exists.

---

## 2. G1 — stiffness

Let \(F^\star\) be the minimizer of \(E\) inside this ansatz.
**Pass** only if every meter in §4 is clean and

\[
\frac{E[F^\star_{1/2}]}{E[F^\star]} \ge 1.05
\quad\text{and}\quad
\frac{E[F^\star_{2}]}{E[F^\star]} \ge 1.05.
\]

**Fail** if every meter is clean and either ratio is below 1.05.
That means the corner piece is not load-bearing in this ansatz.

**Inconclusive** if any meter is dirty. Inconclusive does not open
phase 1, and it is not a physical failure of the corner term.

The 5% floor is the kill. The scaling identity predicts 25%. The
prediction is checked by a meter (§4), not by lowering the floor.

---

## 3. G2 — spin level

**Pass** only if a derivation, written in `SPIN.md`, produces one
odd integer \(k\) such that a \(2\pi\) rotation of a degree-1 frame
has phase \((-1)^k\), using only:

- a partial cycle converts nothing;
- the two-atom credit;
- \(\varepsilon = A_0\omega/2\pi\);
- the frame is the pair of transverse planes on a link.

**Fail** if that derivation inserts \(N_c\), inserts an odd integer,
appeals to QCD anomaly matching, or stops at "the double cover
permits a fermion." Topology permitting spin \(1/2\) is not a
derivation that nature uses it.

G2 does not use the profile numbers.

---

## 4. Meters (all required for a clean G1)

Instrument: `cas/profile.py`. Production value is the fine grid.

| meter | clean if |
|---|---|
| gradient | analytic \(\partial E/\partial y\) matches a central difference to relative \(10^{-6}\) on one random step |
| scaler | on seed A, before minimizing, re-integrated \(E[F_\lambda]\) matches \(\lambda E_2 + E_4/\lambda\) to relative \(10^{-3}\), at \(\lambda\in\{1/2, 2\}\) |
| seeds | seed A and seed B, on the fine grid, agree in \(E\) to relative \(5\times 10^{-3}\) |
| virial | \(\lvert E_2-E_4\rvert/(E_2+E_4) \le 10^{-3}\) on the reported minimizer |
| stationarity | optimizer reports success and \(\lVert\nabla E\rVert_\infty \le 10^{-6}\) |
| refinement | \(\lvert E_{\mathrm{fine}}-E_{\mathrm{mid}}\rvert/E_{\mathrm{fine}} \le 10^{-2}\) |
| textbook | \(\lvert E - 72.9\rvert/72.9 \le 0.02\). This is the published chiral-limit integral of this same functional. A miss is a normalization bug, not a new mass |
| wall | both measured ratios lie within \(2\%\) of \(5/4\). A clean well at 1.10 with a dirty tail is inconclusive |

Grids, fixed now:

| name | intervals \(N\) | box \(X\) |
|---|---:|---:|
| coarse | 500 | 30 |
| mid | 1000 | 60 |
| fine | 2000 | 60 |

Seeds, on every grid:

- A: \(F(x) = 4\arctan(e^{-x})\)
- B: \(F(x) = 4\arctan(e^{-x/2})\)

Decrements are rescaled to sum to \(\pi\) so the box endpoint is 0.
Parameterization is monotone: \(\Delta = \pi\,\mathrm{softmax}(y)\).

Reported \(E(\lambda)\) is a fresh midpoint integral of the linearly
interpolated \(F^\star(x/\lambda)\) on \([0,\lambda X]\), not the
algebraic \(\lambda E_2 + E_4/\lambda\). The algebraic value is only
the scaler meter.

The inertia integral is recorded and is not a gate. Phase 0 does
not claim \((E_\Delta - E_N)/E\). That ratio still contains the
dimensionless Skyrme parameter \(e\), which this continuum unit
has scaled out of the profile and not out of the rotor.

---

## 5. Advance rule

Phase-1 code is written only if G1 passes and G2 passes.

- G1 fail: the corner term does not stiffen the hedgehog. Stop.
- G2 fail: spin \(1/2\) is an import. Stop before any mesh.
- Inconclusive: stop. Do not widen the grid or loosen a meter
  inside this phase and then call it a pass. A later phase 0, if
  the user asks, is a new pre-registration.

`DECIDE.md` quotes the script's meter lines and the G2 verdict.
It does not retune them.
