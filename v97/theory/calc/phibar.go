package calc

import "math"

// The corrected Φ-bar theorem (REVIEW_grok_next Q1; replaces the
// review's F2-as-written, whose algebra was right but whose H1
// identification smuggled a factor Δt and dropped the adopted law's
// clock compensation).
//
// Under laws_V3a with rad_clock=0 the radiance rate per cell is
//
//	O_rate = k_rad · cap · x^p · comp,   comp = w2/w2e = 1 + q_detune·x,
//
// and Em = x·cap, so a uniformly converting motif with balanced
// books (I = O) graded on snaps Δt apart satisfies
//
//	Φ_frame / Em = k_rad · x^{p−1} · (1 + q_detune·x) · Δt.
//
// STATEMENT A (H1): at the defaults (k=0.05, p=4, q=1.2, Δt=4) this
// exceeds the 0.015 frame bar for x ≳ 0.373 — BELOW ambient
// xb = 0.62. The H1 Φ clause is therefore passable at g = 0 by a
// small loaded motif with live, balanced doors. It is an apparatus
// question (field supply, workfn, atom grain, pairing), not a
// stimulated-intake question.
//
// STATEMENT B (H2): what g > 0 buys is not H1's Φ but a RESTORING
// hold above ambient. At g = 0 the only balance is x* = xb (qgraph
// Thm 3.4 — comp cancels between O and the coded dual I, so the
// box-fill root is unchanged). With the stimulated term (Thm 3.5)
// the balance u = x*/xb solves u⁴ − g·u³ − 1 = 0 with u* > 1 and
// restoring slope. Holding, e.g., 8% above ambient needs g ≈ 0.28.
// The g-term is an H2/law face, never an H1 patch.

// PhiFrameOverEm returns Φ_frame/Em for a uniformly loaded,
// balanced-book motif at load x.
func PhiFrameOverEm(k, p, q, x, dt float64) float64 {
	return k * math.Pow(x, p-1) * (1 + q*x) * dt
}

// PhiBarLoadThreshold solves PhiFrameOverEm(x) = bar for x by
// bisection on (0, 1.2].
func PhiBarLoadThreshold(k, p, q, dt, bar float64) float64 {
	lo, hi := 0.0, 1.2
	for i := 0; i < 200; i++ {
		mid := (lo + hi) / 2
		if PhiFrameOverEm(k, p, q, mid, dt) < bar {
			lo = mid
		} else {
			hi = mid
		}
	}
	return (lo + hi) / 2
}

// BalanceU solves u⁴ − g·u³ − 1 = 0 for the positive root u ≥ 1
// (qgraph Theorem 3.5: the load ratio x*/xb of the g-lifted fixed
// point).
func BalanceU(g float64) float64 {
	lo, hi := 1.0, 10.0+g
	f := func(u float64) float64 { return u*u*u*u - g*u*u*u - 1 }
	for i := 0; i < 200; i++ {
		mid := (lo + hi) / 2
		if f(mid) < 0 {
			lo = mid
		} else {
			hi = mid
		}
	}
	return (lo + hi) / 2
}

// GForHold inverts BalanceU: the stimulated fraction needed to hold
// u = x*/xb, g = (u⁴ − 1)/u³.
func GForHold(u float64) float64 {
	return (u*u*u*u - 1) / (u * u * u)
}
