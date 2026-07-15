# B4 orbit / binding kinematics (frozen θ)

**Frozen θ (F15 PASS):** n_C=1 ω_C=1.42, n_L=6 ω_L=1.46, R_shell=18, g=0.05, B1.

**Goal:** Show L cloud is not only rest-stable but follows **Coulomb kinematics**
(sub → inspiral, circ → arc, super → expand) without bag merge / L death.

## Runs

| ID | vt / v_circ | Expect |
|----|------------|--------|
| B4o_rest | 0 | Control (B4 rest redo short T optional) |
| B4o_sub | 0.7 | Shell COM radius shrink / inspiral |
| B4o_vc | 1.0 | Near-circular arc (~¼ period at T=400) |
| B4o_super | 1.3 | Shell expand |

v_circ ≈ g √(Qc Ql_one / (4π R mL)) ≈ **0.071** (code units).

## Pass

- massL loss ≤15%, Qc conserved (c_Q≈0), gauss floor  
- D_shell(t) (mean |r_L − r_C|) trends match sub/circ/super  
- Not full L absorb into C bag  

## Tool

`bin/gen_mf_shell_orbit` + `mf_pair_track` + cfg under `v75/cfg/b4o_*.cfg`
