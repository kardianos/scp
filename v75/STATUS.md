# v75 Status

**Updated**: 2026-07-15 — **F18 P2.1–P2.2–P2.5 smoke PASS** (Z2 class)

## Active goal

**`C12_ATOM_GOAL.md`** — multi-fabric Phases **1–3** → time-stable C₁₂.

| End state | Definition |
|-----------|------------|
| **Ideal** | Time-stable **C₁₂ analog** with firm Z/N |
| **Stretch** | Isotope **+N at fixed Z** + decay |

**Live phase:** P1 open; **P2.0–P2.2/P2.5 smoke PASS at Z=2**; scale Z=6 + P3 next.

## Closed baseline

| Item | Result |
|------|--------|
| F11–F16 | B1 isolation, B4 packaging, pair kinematics, Z6+L6 PASS_park (pre-P/N) |
| **F17** | P/N firm: B2 C-only n / C+Q p |
| **F18** | Z2N0/Z2N2 **PASS_nuc**; L=Z=2 massL/Ql hold; isotope Q_flux identical |

Data: `/space/scp/v75/pn/p2/` · instance `v75f16`.

## Phase checklist

| Phase | Focus | Status |
|-------|--------|--------|
| **P1** | Multi-rev H, shell-radius, binding | OPEN |
| **P2** | P/N + park + L-from-Z + isotope | **Z2 smoke PASS**; Z6 park + long-T open |
| **P3** | A≈12 C₁₂ package | NOT STARTED |

## Freeze (use this)

```
n_fabrics=3  mf_lock_CQ=0  q_C=0 q_Q=1 q_L=-1
gen_pn_core … nZ … nN … nL … [profL omegaL]
L count = Z (not A); same-sign ω; opposite EM via q_L
```

## Next

1. Scale park to **Z≈6 N≈6** (+ L=6) under F18 recipe  
2. Shell-radius diagnostic (P1.3) — COM D≈0 is not enough  
3. Optional P1 multi-rev on single-C + L  
4. Then P3 assemble + long-T visual C₁₂
