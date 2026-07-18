# work/ND — Numeric Dynamics — J5, waves, dual-channel c

**Agent:** ND  
**Round 1:** J5 form-factor; free-wave \(c\); naïve kill.  
**Round 2:** Dual-channel shared \(c\); J5-β locked with numeric; joint TU export.

## Sandboxes

| Script | Demo ID | Round | Purpose |
|--------|---------|------:|---------|
| `sandbox_j5_formfactor.py` | `D-DYN-j5-formfactor` | 1 | J5 triad + form-factor scan |
| `sandbox_free_wave_c.py` | `D-DYN-free-wave-c` | 1 | Free capacity \(v=c\) |
| `sandbox_dual_channel_c.py` | `D-DYN-dual-channel-c` | 2 | ψ wave + EM wave same \(c\); dual-source independence |
| `j5_beta_numeric.md` | `D-DYN-J5` | 2 | J5-β documentation + shape-invariance |

## Run

```bash
cd /home/d/code/scp/v77/work/ND
python3 sandbox_j5_formfactor.py --skip-grid
python3 sandbox_free_wave_c.py --skip-leapfrog
python3 sandbox_dual_channel_c.py
```

## Round-2 key numbers

### Dual-channel shared \(c\) (DC1–DC5 **PASS**)

| Quantity | Value |
|----------|------:|
| \(C_{\mathrm{LOCAL}}\) | 1.0 |
| \(c_\psi = c_{\mathrm{em}}=1/\sqrt{\varepsilon\mu}\) | 1.0 |
| FD2 \(v/c\) (both channels) | **0.99973335** |
| rel_err | **2.67×10⁻⁴** |
| off-unit \(c\) (ε=4,μ=1) | **0.5** both track |
| dual-source cross-talk | **0** |

### J5-β locked

| Gate | Result |
|------|--------|
| ND-G1 naïve kill | **PASS** (rel_split=0.8586) |
| form factor baseline | **0.141421** |
| β1 per-shape \(s_*\) | **PASS** (\(s_*=7.071\)) |
| β1 shape-invariant | **FAIL** (σ=1 vs 1.5) |
| J5-α raw triad | **FAIL** |
| J5-β documentation | **LOCKED** |

## TU export

`outputs/joint_dual_channel_export.json` — V77-2 shared-\(c\) numeric support.
