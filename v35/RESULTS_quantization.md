# V35 Computation 3: Field Value Quantization

## Method

Modified the Cosserat simulation to round all field values to the nearest
multiple of epsilon after each Verlet drift step:

```c
phi[a][i] = round(phi[a][i] / eps) * eps;
theta[a][i] = round(theta[a][i] / eps) * eps;
```

The implied minimum action quantum is hbar_sim = eps * dx.
Rounding occurs AFTER drift, BEFORE force computation, so forces
are computed from the quantized field.

All runs: N=80, L=25, T=50, eta=0.5, m_theta^2=0, dx=0.6329, dt=0.06329.

Source: `src/v33_quant.c`

## Scan 1: Stability vs Epsilon

| eps    | hbar_sim   | E_drift % | theta_rms  | E_pot(0) | E_pot(end) | E_pot retain | survive | compress |
|--------|-----------|-----------|------------|----------|-----------|-------------|---------|----------|
| 0      | 0         | -0.271    | 2.371e-02  | -153.9   | -145.5    | 0.945       | YES     | 1.4x     |
| 0.0001 | 6.329e-05 | -0.226    | 2.371e-02  | -153.9   | -145.5    | 0.945       | YES     | --       |
| 0.001  | 6.329e-04 | +3.511    | 2.688e-02  | -153.9   | -145.7    | 0.946       | YES     | --       |
| 0.005  | 3.165e-03 | +18.816   | 2.539e-02  | -153.9   | -141.8    | 0.921       | YES     | 3.6x     |
| 0.01   | 6.329e-03 | +64.896   | 3.433e-02  | -153.9   | -132.2    | 0.859       | YES     | 4.7x     |
| 0.05   | 3.165e-02 | +1167     | 0          | -153.9   | -127.6    | 0.829       | NO      | --       |
| 0.1    | 6.329e-02 | +3413     | 0          | -153.9   | -427.0    | 2.774       | NO      | --       |

**Critical epsilon: 0.01** (largest eps where braid AND theta field survive to T=50)

At eps >= 0.05, the theta field is completely killed: all values round to zero.
theta_rms ~ 0.025, so eps=0.05 (twice the rms) wipes out most nonzero values.

## Energy Drift Scaling

Energy drift scales as |dE/E| ~ eps^1.4 (superlinear fit across all eps).

Drift rate comparison with expected rounding noise variance (eps^2/12 per value per step):
- eps=0.001: dE/dt = +4.2 (predicted: 4.0) -- quantitative match
- eps=0.005: dE/dt = +49 (predicted: 101) -- field dynamics absorb ~50%
- eps=0.01: dE/dt = +162 (predicted: 404) -- field dynamics absorb ~60%

At small eps, drift is purely stochastic (matches eps^2/12 prediction).
At larger eps, the dynamics partially cancel the rounding errors.

## Scan 2: Detailed Analysis at Critical eps=0.01

### Field Value Distribution

| eps   | Unique values (of 512000 per component) | Non-zero fraction | Value range         |
|-------|----------------------------------------|-------------------|---------------------|
| 0     | 512000 (all distinct)                  | 100%              | [-0.125, +0.126]    |
| 0.005 | 67-68                                  | 90.8%             | [-0.190, +0.170]    |
| 0.01  | 42-44                                  | 75.7%             | [-0.230, +0.220]    |

At eps=0.01, the theta field uses only ~42 distinct quantization levels
spanning [-0.23, +0.22]. 24.3% of grid points have theta rounded to exactly zero.

### Spatial Power Spectrum (3D FFT)

Radially-averaged power spectrum of theta[0] at t=50:

| eps   | S/S_max (spectral entropy) | Dominant k | Power fraction |
|-------|---------------------------|-----------|----------------|
| 0     | 0.1489                    | k=0.186   | 86.2%          |
| 0.01  | 0.1488                    | k=0.186   | 88.5%          |

**The spatial spectrum is essentially UNCHANGED by quantization.**

Both cases show the same dominant mode at k=0.186 (wavelength = 2*pi/0.186 = 33.8 code units,
close to 2L/pi ~ 15.9 -- a box standing wave). Quantization adds slight high-k noise:
modes at k > 1.5 are enhanced 2-25x, but these carry <1% of total power.

**The spectrum does NOT become discrete.** Field value quantization creates a
"staircase" in amplitude space but the spatial structure remains continuous.

### SFA Compression

Quantized fields compress significantly better:
- eps=0: 1.4x compression (nearly incompressible -- all unique values)
- eps=0.005: 3.6x
- eps=0.01: 4.7x

This confirms the quantization is genuinely reducing the field's information content.

## Quantum Implications

### hbar_sim at Critical Epsilon

    hbar_sim = eps * dx = 0.01 * 0.6329 = 6.33e-03 (code units)

### Hydrogen Comparison

From the V35 parameter sweep, hydrogen-like bound states require:
hbar^2 / m_eff ~ 1.1M to 5.5M, where M ~ 8300 is the soliton mass.

At the critical hbar_sim:

| C   | m_eff      | m_eff / M   |
|-----|-----------|-------------|
| 1.1 | 4.39e-09  | 5.29e-13    |
| 2.0 | 2.41e-09  | 2.91e-13    |
| 5.5 | 8.78e-10  | 1.06e-13    |

This m_eff is ~10^13 times smaller than the soliton mass.

### Physical hbar

Using code unit conversions (1 code E = 9.098 MeV, 1 code L = 0.5624 fm):

    hbar_physical = 197.3 MeV*fm / (9.098 * 0.5624) = 38.6 code units
    eps_physical = hbar_physical / dx = 38.6 / 0.633 = 61

**Physical hbar requires eps = 61, which is 76x the braid core amplitude (0.8).**
This would completely destroy all field structure.

Physical hbar is 6100x larger than hbar_sim at the critical eps.

## Conclusions

1. **Stability threshold**: eps_crit = 0.01 (10% of background amplitude,
   1.25% of core amplitude). Above this, quantization noise overwhelms the
   theta field. The phi braid persists to higher eps because its amplitude
   (0.8) is 8x larger than the critical eps.

2. **Spectrum unchanged**: Field value quantization does NOT discretize the
   spatial power spectrum. The theta field uses only ~42 amplitude levels at
   eps=0.01, but its spatial structure and dominant wavenumber are preserved.

3. **Energy non-conservation**: Rounding injects energy at dE/dt ~ eps^1.4.
   At small eps (0.001), the rate matches the stochastic prediction eps^2/12
   per value per step. At larger eps, dynamics partially absorb the noise.

4. **Physical hbar unreachable**: hbar_physical = 38.6 code units requires
   eps = 61, which is 6100x larger than the critical eps. Field value
   quantization at physical hbar would obliterate all structure.

5. **Implied m_eff is unphysical**: At the critical eps, m_eff ~ 10^-13 M.
   No known particle has a mass ratio this extreme relative to the proton.

6. **Bottom line**: Substrate granularity (field value quantization) is NOT
   a viable mechanism for introducing hbar. The stability limit is too
   restrictive (eps_crit = 0.01) and the resulting hbar_sim is too small
   by a factor of ~6000. This is analogous to real physics: the Planck
   scale sets the substrate granularity, but atomic quantization operates
   at 10^24 x larger scales. hbar must emerge from a different mechanism.

### Remaining paths to hbar

1. **Soliton de Broglie wavelength**: If theta forms solitons with mass m_eff,
   the effective hbar comes from the soliton's wave-particle duality, not
   the substrate. hbar_eff = m_eff * v * lambda_dB.

2. **Emergent from braid dynamics**: The braid's oscillation frequency
   omega=0.22 sets an action scale E/omega ~ 5000/0.22 ~ 23000. This is
   much larger than the substrate hbar and could be the relevant scale.

3. **Fine structure constant**: alpha = e^2/(hbar*c) ~ 1/137. If alpha is
   fixed by the dynamics and the charge quantum e (from winding number) is
   fixed, then hbar = e^2/(alpha*c) is determined. The free parameter
   shifts from hbar to alpha.

## Files

- `src/v33_quant.c` -- quantized Cosserat simulation
- `data/quant_eps_*_timeseries.tsv` -- energy/diagnostics timeseries at each eps
- `data/quant_eps_*_theta_final.bin` -- raw theta field for FFT analysis
- `analyze_quant.py` -- scan summary analysis
- `analyze_detailed.py` -- detailed quantum implications
- `fft_analysis.py` -- spatial FFT and field value distribution
- `drift_scaling.py` -- energy drift scaling with eps
