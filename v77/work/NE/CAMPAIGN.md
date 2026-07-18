# NE campaign card — Phase 2

**Map:** `../../CAMPAIGN_MAP.md`  
**Mission:** Comprehensive **dynamic** full Maxwell; Maxwell step API for RC1.

## Track

1. **CP0** — Keep M0 regression (`sandbox_full_maxwell_r2.py`).
2. **CP-M1-SPEC** — Co-agree TE gate list; design M1 suite.
3. **CP-M1-NUM** — Implement 2D beams, energy, dynamic Gauss, radiation, Ampère adversary; `outputs/m1_result.json`.
4. **CP-RC1-SPEC** — Freeze **API**: `step(rho_Q, J_Q) → E,B`; NM must call this (not reimplement Φ-only).
5. **CP-RC1-NUM** — Support/stamp NM co-field EM channel.
6. Optional **M2** 3D Yee; **M3** self-consistent charges.

## Hard rules

- M1 wave proof must **not** rely only on 1D TEM (1D may remain as regression).
- RC1 REJECT if NM never calls dynamical Maxwell.

## Co-agreement partners

| CP | Must co-agree |
|----|----------------|
| CP-M1-* | **TE** |
| CP-RC1-* | **NM, TE** (+ TM) |
| shared \(c\) | **ND** (note stamps) |
