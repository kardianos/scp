# v99 Phase 1 — execution (2026-08-19)

Paper only. Maxima and numpy. **1C is written:**
`phase1/DECIDE.md`. Calculation: `phase1/RELENT.md`.

---

## What was run

Gaussian relative entropy of a `W = 1` coherent wrap against the
vacuum of a free scalar on `C_n`, restricted to 1 site, 2 sites, or
a half-cycle. Identities for `E_w`, `ω`, `ε = ω/2π`, and the
discrete walks `ΔE_w + ε` and `ΔE_w + 1`.

## Result

Ending **C**. The selector does not force unique size. The
representation split is held. No kernel. Phase 3 not authorized.

Headline numbers (see `RELENT.md` for the rest):

- Half-cycle `S/K_{\mathrm{Casini}} ≈ 0.29–0.35` (`n = 4..24`).
  Not a 2π-flux theorem.
- `E_w + nε` hoards (`n* = 24`). `E_w/ε + n` sits at `n = 3`.
- `ΔE_w + ε < 0` on `n = 3..23`. Spreading always pays.
- Grain door at `m = 1`: `S_1 ≈ 1.55` independent of `n`.
- `ω_1 = 1` gives `n = 6`; `ω_1 = 1/(2π)` gives `n ≈ 39`.

## Reproduce

```
python3 v99/phase1/cas/srel.py
maxima --very-quiet -b v99/phase1/cas/srel.mac < /dev/null
```
