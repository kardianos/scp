# The Φ bar at g = 0 — corrected theorem (2026-08-13, PLAN step 3a′)

Replaces the fresh-eyes review's F2-as-written, per the grok consult
(`consult/REVIEW_grok_next.md` Q1). Go checks: `calc/phibar.go`,
tests green. **No bar moves; this is what the bars mean.**

## Statement A — H1's Φ clause is passable at g = 0

Under laws_V3a with `rad_clock=0` the per-cell radiance rate is
`O = k_rad·cap·x^p·comp` with `comp = w2/w2e = 1 + q_detune·x`, and
`Em = x·cap`. A uniformly converting motif with balanced books
(I = O) graded on snaps Δt apart therefore satisfies

    Φ_frame / Em = k_rad · x^{p−1} · (1 + q_detune·x) · Δt.

At the defaults (k=0.05, p=4, q=1.2, Δt=4) this exceeds the 0.015
frame bar for **x ≳ 0.373 — below ambient xb = 0.62**. At ambient
load the margin is ~3×. Two consequences:

- The H1 Φ clause is an **apparatus question** (field supply, workfn,
  atom grain, I/O pairing in time), not a stimulated-intake question.
  A small loaded motif with live balanced doors can pass at g = 0.
- The bar is **snap-dependent** (a stream snapped 4× finer fails 4×
  more easily). The meters now print the rate alongside; restating
  the bar with Δt is charter debt for the user, not this page.

The two arithmetic mistakes the original F2 made, kept as a record:
identifying the per-t.u. qgraph rate with the per-snap H1 integral
(smuggles Δt = 4), and dropping `comp` (at xb: 0.0119 bare vs 0.0208
with comp, either side of 0.015 — the "provably unpassable" claim
lived entirely in those two omissions).

## Statement B — what g > 0 actually buys is H2, not H1

At g = 0 the only balance point is x* = xb (qgraph Thm 3.4; `comp`
multiplies O and the coded dual I alike, so the box-fill root is
unchanged). The stimulated cross term (Thm 3.5) lifts the fixed
point: u = x*/xb solves u⁴ − g·u³ − 1 = 0, restoring, with u* > 1.
Holding 8% above ambient needs g ≈ 0.28 (`GForHold`).

So the g-term is the **H2 face** (a restoring root above ambient — a
law-table change with v92-closed-feedback priors to clear), and it
must never be introduced as "just the Φ bar in disguise."

## Statement C — census dilution (why the ~100-Em blobs still fail)

Φ/Em on a *census* is `k·Δt·Σ x_c^p(1+q·x_c) / Σ x_c` — the halo the
Em>0.05 floor drags in adds denominator and, at p = 4, almost no
numerator. The leads blobs fail the Φ bar mostly by ⟨x^{p−1}⟩ and
halo, not by missing g. Small, loaded, door-live motifs (the C6 at
x ≳ 0.4) are the right hunting ground — matching the A.5 arm choice.
