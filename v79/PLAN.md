# v79 GPU plan

## Instance
- GPU: Tesla V100 32GB (`min_ram=32` or 32GB V100)
- Disk: **80 GB** minimum (prefer 100) for multi SFA ~10–40GB each + scratch
- Name: `v79atom`

## Build
- Sources: `sfa/sim/scp_sim.cu` + `sfa/format/sfa.h` (multi-fab is in scp_sim.cu)

## Seeds (gen_pn_core)
- Profile: `v74/profiles/f_w146_g005.txt`
- Nuclear ω=1.46; light hierarchy same or 1.46
- Z=6 octa R=8; N=6 octa R=5.5; L=6 octa R=22

## Runs (priority order)
1. `z6_a_n6` T=800 — primary long-T atom (Z6N6+L6)
2. `z6_n6_nuc` T=800 — core-only control
3. Optional: `z6_a_n0` long-T if time

## Downloads
- diag.tsv, track if any, campaign logs first
- SFA f16 or selected snaps; full f32 to /space/scp/v79/
