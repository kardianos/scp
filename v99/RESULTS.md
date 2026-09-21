# v99 — where the results live

Read `RESUME.md` first. This page is only an index.

---

## Phase 1 (paper + CAS)

| file | what it is |
|---|---|
| `PHASE1.md` | one-page execution summary |
| `phase1/RELENT.md` | the calculation: doors, identities, named readings |
| `phase1/DECIDE.md` | **1C verdict (ending C)** |
| `phase1/cas/srel.py` | Gaussian `S_rel` on `C_n` |
| `phase1/cas/srel.mac` | Maxima: `E_w`, `ε`, argmins, discrete `Δ` |
| `phase1/cas/out/` | logs and TSVs |

**1C in one sentence:** `S_rel` of a wrap against vacuum is
well-defined on a proper door and does not force unique size;
the representation split is held; no kernel.

## Kernel face

| file | what it is |
|---|---|
| `KERNEL.md` | `pattern_only` 0/1/2; conversion rules; what is not claimed |
| `v97/fab/pattern.go` | Go helpers (absorb, cond/evap, wrap-gate) |
| `v97/kernel/freecell.c` | C mirror |

**Battery at defaults (`pattern_only=0`):** ALL GREEN (87 bars),
including C==Go conserve and five A/B arms. 2026-08-19.

## Not results

- `n* = 6` is not derived.
- Gravity / 1/r was not obtained.
- Phase 3 is not authorized.
- The v98 size meter was not rerun (no new stream).
