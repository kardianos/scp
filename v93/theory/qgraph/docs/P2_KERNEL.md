# P2 kernel — dual intake (implemented 2026-08-08)

**Auth:** user explicit ("You have auth for P2 kernel... Go.")
**Baseline:** `v93/theory/p2/baseline/` (pre-change fab + freecell).
**Design:** `v93/theory/p2/DESIGN.md`.

## Knobs (defaults inert)

| key | default | meaning |
|---|---:|---|
| `rad_dual` | 0 | scale on dual intake `I = k_rad·rad_dual·cap·xb^p` |
| `bath_xb` | 0 | apparatus: hold vacuum-cell `Ee = bath_xb·cap` |
| `bath_em` | 0.01 | Em floor separating vacuum bath from matter |

## Where

- Go: `fab/step.go` pass-6 (after radiance), `fab/params.go`, ledger `radInTotal`
- C: `kernel/freecell.c` same placement (kernel of record)

## Suite

```
./theory/bin/probe -suite p2 -kernel fabrun -T 600 -j 2
```

9 arms: xb∈{0.3,0.6,1.0,1.5}×rad_dual∈{0,1} + k_rad=0 control.

## Battery

Full battery launched after the change; `rad_dual=0` must keep V3a surface green.
