# Phase 2 — size-spectrum bars

Pre-registered. Sharpen by measurement; do not soften to pass.
Implemented in `size/h0.go`. CLI: `go run ./cmd/sizespec -h0 CLASS FILE.fcs`.

A **line** (used by more than one class) is a compact species-shaped
peak:

```
top-bin fraction of motif-counts  ≥ 0.70
AND  mean(N_eff / N)              ≥ 0.50
```

`N` is the cell count of a connected live-Em motif (`Em > 0.05`,
`N ≥ 3`). `N_eff = (Σ Em)² / Σ Em²` on that motif.

---

## broad  (Face C must pass; a vacuum C6 must fail)

PASS if the stream is **not** a line, and at least one of:

- distinct `N` values ≥ 5
- top-bin < 0.40
- mean(`N_eff/N`) < 0.50   (graded clump: few components, internally spread)

FAIL if it is a line.

## line6  (vacuum C6 must pass; Face C must fail)

PASS if both:

- dominant motif has `N = 6` on ≥ 0.80 of frames
- among those `N = 6` frames with a chart, the unweighted
  Laplacian matches `{0,1,1,3,3,4}` to `1e-6` on ≥ 0.80

## ceiling  (tax-only must pass; a parked C6 must fail)

PASS if both:

- `Em_last / Em_first ≤ 0.85`  (inventory drops)
- the stream is **not** a line

A tax that drained *onto* a single compact `N` would fail this
class: that would be size selection, which MONO_FINAL said
radiance does not do.

---

Negative controls (same bars, opposite streams) are part of H0.
If a negative passes, the meter is wrong.
