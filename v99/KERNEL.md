# v99 kernel face — `pattern_only`

**Authorized:** user, 2026-08-19, “Do Tact 1.”
**Default 0.** Battery at defaults must stay green. Opening this
flag is not a law-table change.

Live kernels: `v97/kernel/freecell.c` and `v97/fab/` (A/B).

---

## What it is

The kernel had a **second traveling species**: scalar `Em` plus
`(dm1, dm2)` dense Givens, with beat-gated atoms copying `Ee ↔ Em`
and pulling `Es`. That is the Face C fluid.

Tact 1: when `pattern_only≥1`, pattern is `(fa1, fa2)` only.
Space is `Es`. Ledger `Σ(Es + |fa|² + |χ|²)`. `Em` is absorbed
into `fa` after seed and at the start of every step.

| value | meaning |
|---:|---|
| **0** | legacy. Byte-inert vs the pre-face kernel. |
| **1** | no dense hop; condensation is space↔pattern; clicks absorb `fa` into `Es` on `scond` |
| **2** | 1 plus a derived wrap-gate: a field hop is skipped if the link sits on a space-deficit triangle or 4-cycle whose `fa` phase winds (`|ΣΔθ| > 3`) |

`C_6` has girth 6. The wrap-gate does **not** lock a hexagon.
That is recorded, not papered over.

---

## Conversion (pattern_only ≥ 1)

- **Screen (`scond`):** field → space. `Ee -= d1`, `Es += d1`. Click
  recorded. The beam is absorbed.
- **Bulk:** space → pattern. `Es -= dsp`, `fieldInject(dsp)`.
  Triggered by the same beat/atom demand; traveling `fa` is not a
  source tank.
- **Over cap:** pattern → space. `Ee -= d2`, `Es += d2`.
- **Workfn presence:** space deficit `e_s0 - Es ≥ wf_floor`, not `Em`.
- **Door pair / radiance / dual intake:** skipped or inert (`Em = 0`).
  `door_g` settlement is not taken (no dense register to settle into).

Pitch load `x` is `(Ee + flload)/cap`.

---

## What this does not claim

Unique size is not forced. Slice 1 removes the false species.
Slice 2 is a local hop predicate, not `n*`. The v98 `ρ` gap stands.

Reproduce:

```
cd v97
make all
./battery                    # defaults, pattern_only=0
./fabrun exp=bath T=10 pattern_only=1 k_rad=0.05 wf_on=1
./fabrun exp=blob T=10 pattern_only=1 k_rad=0.05 wf_on=1
./fabrun exp=slit slit_clicks=1 T=20 pattern_only=1
```
