# T5: Self-Consistent Metric Coupling — Results

## Setup
- Grid: N=96, L=20, dx=0.4211, dt=0.0842
- Mass: m=1.50, m^2=2.25
- T=300 evolution with absorbing damping
- Metric: g^{ij} = delta^{ij} - alpha_g * h_{ij}, where h_{ij} = d_i phi_j + d_j phi_i
- Far-field l=2 decomposition on R=10 shell (50 golden-ratio spiral points)
- Total wall time: 1007 s (~17 min)

## Summary Table

| alpha_g | fc     | trans_l2 | torsion | l2_h     | energy  | stable |
|---------|--------|----------|---------|----------|---------|--------|
| 0.0000  | 0.8337 | 0.0871   | 0.2229  | 0.6261   | 1462.3  | OK     |
| 0.0001  | 0.8337 | 0.0871   | 0.2229  | 0.6261   | 1462.3  | OK     |
| 0.0010  | 0.8337 | 0.0871   | 0.2228  | 0.6260   | 1462.3  | OK     |
| 0.0050  | 0.8337 | 0.0870   | 0.2222  | 0.6254   | 1462.3  | OK     |
| 0.0100  | 0.8337 | 0.0869   | 0.2216  | 0.6248   | 1462.3  | OK     |
| 0.0500  | 0.8338 | 0.0858   | 0.2162  | 0.6191   | 1462.2  | OK     |

## Key Findings

### 1. Stability across all alpha_g
All six coupling values (0 to 0.05) remained stable through T=300. No blowups.
The braid is robust against metric backreaction in this range.

### 2. l2_h (metric quadrupole fraction) is HIGH but nearly constant
The deviatoric (l>=2) fraction of h_{ij} on the R=10 shell is ~62% at all alpha_g values.
This is dominated by the **intrinsic** asphericity of the bimodal braid, not by
the metric coupling. The l2_h decreases very slightly with alpha_g (from 0.626 to 0.619),
meaning the metric backreaction marginally isotropizes the far field.

### 3. No amplification of spin-2 content
The critical test was whether self-consistent coupling amplifies the l=2 content.
It does not. The trend is:
- delta(l2_h) / delta(alpha_g) ~ -0.14 over the range [0, 0.05]
- The coupling **suppresses** rather than amplifies the quadrupolar metric component.

### 4. Core confinement unaffected
fc (core fraction) is essentially constant at 0.8337-0.8338 across all alpha_g.
The metric backreaction does not measurably affect the braid's self-confinement.

### 5. Torsion decreases with coupling
Torsion flux drops from 0.2229 (control) to 0.2162 at alpha_g=0.05, a 3% reduction.
The modified propagation slightly reduces the internal twist structure.

## Interpretation

The bimodal braid has a large intrinsic l=2 content (~62%) in h_{ij} at R=10 due to
its aspherical (elliptical) shape. However, self-consistent metric coupling does NOT
amplify this signal. The backreaction is perturbatively small even at alpha_g=0.05,
and its net effect is slightly suppressive. This is consistent with the negative results
from earlier gravity investigations: the strain tensor h_{ij} is a kinematic proxy
for asphericity, not a dynamical gravitational wave source.

The braid remains stable under metric backreaction, confirming it is a genuine
nonlinear attractor of the modified dynamics, not an artifact of the flat-metric EOM.
