# v75 Status

**Updated**: 2026-07-15 — **F17 P/N PASS** (Z/N charge distinction firm)

## Active goal

**`C12_ATOM_GOAL.md`** — multi-fabric Phases **1–3** → time-stable C₁₂.

| End state | Definition |
|-----------|------------|
| **Ideal** | Time-stable **C₁₂ analog** with **firm Z/N** |
| **Stretch** | Isotope **+N at fixed Z** + decay calc/sim |

**Live phase:** P1 (bound hydrogenoid) open; **P2.0 DONE (F17)**; P2 packaging + P3 next.

## Closed baseline

| Item | Result |
|------|--------|
| B1 MF isolation + Coulomb (F11–F13) | PASS |
| Self-tune F15 (B4 packaging) | PASS |
| F16 pair kinematics + Z6+L6 PASS_park | PASS |
| **F17 P/N (B2 C-only n / C+Q p)** | **PASS** — \(Q_\mathrm{em}\propto Z\); flavor multiplet alone fails neutrality |

Data: `/space/scp/v75/pn/` · `PN_EXPERIMENT.md` · instance `v75f16`.

## Phase checklist

| Phase | Focus | Status |
|-------|--------|--------|
| **P1** | Retuned v_c, multi-rev B4, shell-radius, binding | OPEN |
| **P2** | **P/N firm** + parked multi-Z + L + isotope smoke | **P2.0 PASS**; rest open |
| **P3** | (Z,N) A≈12 + multi-L; visual C₁₂; +N decay | NOT STARTED |

## P/N freeze (use this)

- **p** = C bag + Q charge (`init_sfa_Q` / `gen_pn_core` Z balls)
- **n** = C only (Q empty at neutron sites)
- Config: `n_fabrics=3`, `mf_lock_CQ=0`, `q_C=0`, `q_Q=1`, `q_L=-1`
- Tool: `bin/gen_pn_core`

## Ops

1. Reuse or teardown `v75f16`
2. Next: P1 multi-rev **or** P2.1 parked templates with flavored p+n mixes + L from Z
3. Park-aware scorecard mandatory for multi-ball
