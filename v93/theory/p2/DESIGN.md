# P2 — dual intake kernel (authorized 2026-08-08)

**Status:** implemented behind knobs; not adopted into laws_V3a.
**Baseline snapshot:** `theory/p2/baseline/` (pre-change fab + freecell).

## Knobs

| key | default | role |
|---|---|---|
| `rad_dual` | **0** | scale of spontaneous dual intake. 0 = byte-inert. 1 = full dual of radiance. |
| `bath_xb` | 0 | apparatus: hold ambient field load on vacuum cells (`Em < bath_em`) at `Ee = bath_xb * cap`. |
| `bath_em` | 0.01 | Em threshold separating vacuum bath cells from matter. |

## Physics (qgraph UNIFIED_LAW §3.3)

```
O = k_rad * cap * x^p_rad          (existing radiance; x = Em/cap)
I = k_rad * rad_dual * cap * xb^p  (new; xb = Ee/cap)
```

At `rad_dual=1`, balance root is `x* = xb` (Theorem 3.4). Reverse
bookkeeping of radiance (field shrink + Es reverse-backs) so Em+Ee+Es
conserves. Gated on `Em>0` so vacuum does not mint matter.

## Battery

`rad_dual` defaults 0 → V3a surface unchanged. Full battery gates the change.

## Suite

`probe -suite p2`: xb ∈ {0.3, 0.6, 1.0, 1.5} × {rad_dual=0,1}, T≥500,
`bath_xb=xb`, `k_rad=0.05 p_rad=4`. Bar: dual on → x* tracks xb slope ~1;
dual off → x* bath-independent (slope ~0).
