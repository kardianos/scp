# Orbit2 — Gate O long-hold (N=128)

**Goal:** convert O SOFT → PASS/FAIL with kinematics-correct runs.  
**GPU:** V100-16GB (~2.6 GB VRAM).  
**Grid:** N=128 L=32 multi-fab B1.

| Job | \(v_t\) | \(T\) | snap | purpose |
|-----|---------|-------|------|---------|
| hard | 0.50 | 1200 | 10 | multi-rev (≥2) |
| soft | 0.03 | 4000 | 40 | long inspiral / capture |

Policy: track after each job; sparse snaps; f16.
