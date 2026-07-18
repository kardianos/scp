# NM campaign card — Phase 2

**Map:** `../../CAMPAIGN_MAP.md`  
**Mission:** **Lead R-compose** — RC1 co-field mandatory; RC2 stretch.

## Track

1. **CP0** — Keep R2 dual lite as regression only. **done**
2. **Wait CP-M1-NUM** — **ADOPTED** (O-011 / m1_claim=true).
3. **CP-RC1-SPEC** — **ADOPT** (NM-011); TM `rc1_joint_state_v0.md` + NE Maxwell2D.
4. **CP-RC1-NUM** — `sandbox_rc1_cofield.py` implemented; run `python3 run_rc1.py` → `rc1_claim`.
5. **CP-RC2-*** — Optional moving locks.

## Hard rules

- RC1 is **not** closed by Maxwell-lite \(\Phi\) alone.
- Must call NE dynamical Maxwell module (or embedded equivalent with TE/NE stamp).
- Sibling gates: neutral / same / opposite; shared \(c\); \(\psi\neq\Phi\).

## Co-agreement partners

| CP | Must co-agree |
|----|----------------|
| CP-RC1-SPEC/NUM | **NE, TE, TM** (+ TU) |
| CP-RC2 | **TD, ND, TM** |
