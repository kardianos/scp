# Local Density and Proper Time in the 6-Field Cosserat Theory

## The density–time relation

The scalar density \(P = \phi_0\phi_1\phi_2\) sets the local rate of proper time. The
rate is read from the *dispersion relation of small field fluctuations* on a background
of density \(P\): a higher-density background slows the propagation of fluctuations,
which is the local lapse \(\alpha(x)\). The dynamically grounded form follows from the
saturating potential \(V(P) = \tfrac{\mu}{2}\,P^2/(1+\kappa P^2)\) that governs the
field (kernel `scp_sim.c`, with \(\kappa=50\) fixed), giving

\[
\alpha(x) = \frac{1}{\sqrt{1 + \kappa P(x)^2}} .
\]

Regions of higher \(P\) therefore experience slower proper time, exactly as required by
the equivalence principle. (The earlier candidate \(1/\sqrt{1+\kappa P}\) is the linear
reading of the same idea; the quadratic form is the one the field equations actually
carry, since the density enters the dynamics only through \(1+\kappa P^2\). Thread A of
the plan *measures* \(\alpha(P)\) from the linearized dispersion rather than positing
it — see [`PLAN.md`](PLAN.md).)

## Gradient ⇒ drift and spin

A spatial gradient \(\nabla\alpha\) produces differential aging across an extended body:
the higher-density side advances more slowly in coordinate time. The resulting mismatch
in spatial displacement per unit proper time generates a net drift together with a
torque. The drift obeys

\[
\ddot{x}^i \;\propto\; -\,\partial^i \ln\alpha ,
\]

which points **up** the density gradient, toward slower proper time — the *attractive*
sign of gravity. (The intuitive "differential aging ⇒ drift and spin" picture is correct
in mechanism; the direction is toward *slower*, not faster, time.)

Because the underlying 6-field dynamics couple the \(\phi\) triplet to the \(\theta\)
triplet through the \(\eta\,\mathrm{curl}\) terms, the drift is automatically
accompanied by a torsional response: the same density gradient that produces the drift
sources the spin (frame-dragging / spin-orbit) through the \(\theta\)-sector. The
torque is therefore not an add-on — it is the \(\theta\)-channel image of the lapse
gradient.

## The line element, and what sources the bending

The leading effective line element felt by test fluctuations is

\[
ds^2 = -\alpha^2(x)\,dt^2 + \delta_{ij}\,dx^i dx^j ,
\]

with the time-dilation arising from \(\partial_i\alpha\). This scalar (lapse-only) form
captures the drift but is **not** the full gravitational sector: a static, purely
\(g_{00}\) metric carries no transverse-traceless polarizations and is incompatible with
observed tensor gravitational waves (the v60 result). The physical bending must inherit
its two tensor degrees of freedom from the dynamical/torsional sector — the part of the
theory that is *nonlinear* (the \(\kappa\) saturation and the \(\eta\,\mathrm{curl}\,
\theta\) coupling). In this theory the gravitational response is sourced by the
nonlinear sector: a free/linear truncation is flat and inert, and the lapse gradient
itself only appears once \(\kappa\) is on.

## Regime — the lapse mechanism is nuclear, not gravitational

The lapse mechanism is real, but its *scale* settles its identity, and the scale is
nuclear. The effective potential is \(\Phi/c^2 = \ln\alpha\). At a representative core
density \(P=1\) this is \(\Phi/c^2 = -1.97\) (proper time runs \(\sim 7\times\) slower
at the core), against a *gravitational* \(\Phi/c^2 = Gm_p/(r_p c^2) \approx 1.5\times
10^{-39}\) at the proton scale — a factor \(\sim 1.3\times 10^{39}\) too strong. The
magnitude pre-check (`magnitude_precheck.py`) confirms this holds across every plausible
core density (\(P\in[0.1,320]\)) and for **both** lapse forms (\(1+\kappa P\) and
\(1+\kappa P^2\)): the smallest excess is \(\sim 10^{38}\). With a Yukawa range set by
the soliton core (\(\sim 0.5\)–\(0.8\) fm, not \(1/r\)), this exactly reproduces the
project's prior BLV / null-rotor result: a too-strong, short-range, nuclear-scale
effect.

The density–time relation above is therefore the correct *mechanism* — and it is the
mechanism of nuclear binding response, not of gravity. Gravity, by the v62/v63
number-type result and the "magic" reading, is not reachable from a static algebraic
function of the density; it would have to come from the dynamical/torsional sector and
at a vastly weaker coupling. The forward program (see [`PLAN.md`](PLAN.md)) is therefore
to pin the density–time *law* precisely (Thread A) and document *why* the lapse sector
is short-range, rather than to pursue it as a gravitational mechanism.

## Immediate next steps
1. ~~Magnitude pre-check (`magnitude_precheck.py`)~~ — **done**: nuclear wall confirmed
   (\(\sim 10^{38}\)–\(10^{39}\times\) too strong, both lapse forms; see Regime above).
2. Thread A: derive \(\alpha(P)\) from the linearized 6-field dispersion on a controlled
   density background; compare to \(1/\sqrt{1+\kappa P^2}\).
3. Thread B: evolve a wave packet in a `tanh` density ramp (real 6-field run, no imposed
   metric) and measure drift \(\ddot{x}(\partial_x P)\) and spin \(\omega(\partial_x P)\),
   checking the spin tracks \(\theta\)-rms growth.
4. Thread C: confirm the response is "magical" — zero drift in the linear limit
   (\(\kappa\to 0\), \(\eta\to 0\)); drift scaling with \(\kappa\) and \(\eta\).
