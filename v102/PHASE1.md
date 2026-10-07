# v102 Phase 1 — the derived chain, as an instrument

**Opened:** 2026-10-07, under rulings R1, R2, R5, and R6
(`README.md`). Bars were locked in `phase1/PREREG.md` before any code.
**Decision:** `phase1/DECIDE.md`.

No kernel. Every calculator here is numpy/scipy and lives in
`phase1/cas/`. The DMRG is written in-tree (`dmrg.py`) and is
validated against the spin-1 chain (P2) before it is trusted.

| bar | what | script |
|---|---|---|
| P1 | microscopic ring with transport → spin-1 Heisenberg levels | `p1_micro.py` (plus diagnostics `p1_diag.py`, `p1_tc.py`, `p1_tc6.py`) |
| P2 | DMRG instrument: energy, gap, ξ, two-magnon check | `p2_validate.py` |
| P3/P5 | ξ from the microscopic chain; dependence on g/Δ and t/Δ | `p3_micro_dmrg.py` (did not converge; kept), `p3_cluster.py` (exact cluster effective Hamiltonian plus validated DMRG; used) |
| P4 | one line, no second line | `p1_micro.py` (P4a, P4b), `p2_validate.py` (open chain) |
| P6 | slit from the unlocked linearized field | `p6_slit.py` |

Reproduce:

```
cd v102/phase1/cas
python3 -I p2_validate.py
python3 -I p1_micro.py
python3 -I p1_diag.py
python3 -I p1_tc.py
python3 -I p1_tc6.py
python3 -I p6_slit.py
python3 -I p3_cluster.py 48 48
python3 -I p3_micro_dmrg.py eff 0 48 48
python3 -I p3_micro_dmrg.py 0.1 0 48 32        # did not converge; kept as record. also: 0.1 0 48 48 | 0.2 0 48 32 | 0.1 0.1 48 32 | 0.2 0.2 48 32 | 0.1 0 48 32 1
```
