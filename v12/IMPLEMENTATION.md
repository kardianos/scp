# V12 Implementation Traceability

## Phase 1: BI-Modified Hartree (V9 with b_BI)

Phase 1 uses V9's existing `strong_geon.c` with the `-bBI` parameter enabled.
No code modifications — V9 already implements the scalar BI formula.

### Scalar BI Energy Density

**Equation (PLAN.md Section 3, Level 1):**
```
rho_M = N_l/(8pi r^2) [omega^2 u^2/f + f (u')^2]
rho_BI = b^2 (sqrt(1 + rho_M/b^2) - 1)      if rho_M > 0
```

**Code** (v9/src/strong_geon.c:214-239):
```c
double N_l = ell*(ell+1) / (2*ell+1);
double rho_M = N_l / (8.0 * M_PI * r * r) *
               (omega2 * u * u / f + f * up[i] * up[i]);

if (b_BI > 0.0 && rho_M > b_BI * b_BI) {
    rho[i] = b_BI * b_BI * (sqrt(1.0 + rho_M / (b_BI*b_BI)) - 1.0);
} else {
    rho[i] = rho_M;
}
```

**Mapping:**
- N_l = ℓ(ℓ+1)/(2ℓ+1) is the angular normalization factor (2/3 for ℓ=1)
- rho_M: Maxwell energy density from Regge-Wheeler mode (omega^2 u^2/f = electric, f u'^2 = magnetic)
- BI formula applied only when rho_M > b_BI^2 (condition is redundant but harmless)
- b_BI is the command-line parameter `-bBI` (default 0 = off)
- rho_BI >= 0 for all rho_M >= 0 (no sign issues)
- Weak field (rho_M << b^2): rho_BI = rho_M (Maxwell limit)
- Strong field (rho_M >> b^2): rho_BI ~ b*sqrt(rho_M) (square-root cap)

### Yukawa BVP

**Equation (PLAN.md Section 9):**
```
M'' - mu^2 M = -kappa r rho(r)
M(0) = 0,  M'(Rmax) = -mu M(Rmax)
```

**Code** (v9/src/strong_geon.c:242-296):
```c
/* Thomas algorithm for tridiagonal system */
/* Interior: M_{i-1} - (2 + mu^2*h^2)*M_i + M_{i+1} = -h^2 * kappa * r_i * rho_i */
```

**Mapping:**
- 2nd-order centered differences
- BC at origin: M(-h/2) = -M(h/2) (antisymmetry ghost point)
- BC at Rmax: M(Rmax+h) = M(Rmax)*exp(-mu*h) (Yukawa decay)
- Thomas algorithm: O(N) forward sweep + back substitution

### Pade Metric

**Equation (PLAN.md Section 3):**
```
f(r) = 1/(1 + 2M_eff(r)/r)
```

**Code** (v9/src/strong_geon.c:105-107, nl_metric=2):
```c
case 2: f[i] = 1.0 / (1.0 + 2.0 * M_eff[i] / r[i]); break;
```

**Mapping:**
- f > 0 for all M_eff > 0 (Pade denominator always positive for attractive well)
- Agrees with Schwarzschild f = 1 - 2M/r to O(M/r) in weak field
- No artificial floor needed (unlike linear metric which needs floor at 0.01)

---

## Phase 3: 3D BI-Hopfion in Yukawa Metric

Source: `v12/src/bi_geon.c`

### D-field BI Energy Density

**Equation (PLAN.md Section 9):**
```
R = 1 + (D^2 + B^2)/b^2 + |D x B|^2/b^4
rho_BI = b^2 (sqrt(R) - 1)
```

**Code** (bi_geon.c:223-238):
```c
double D2 = dx_f*dx_f + dy_f*dy_f + dz_f*dz_f;
double B2 = bx_f*bx_f + by_f*by_f + bz_f*bz_f;
double DdB = dx_f*bx_f + dy_f*by_f + dz_f*bz_f;

double DxBx = dy_f*bz_f - dz_f*by_f;
double DxBy = dz_f*bx_f - dx_f*bz_f;
double DxBz = dx_f*by_f - dy_f*bx_f;
double P2 = DxBx*DxBx + DxBy*DxBy + DxBz*DxBz;

double R = 1.0 + (D2 + B2) / b2 + P2 / b4;
double sqrtR = sqrt(R);
*rho_out = b2 * (sqrtR - 1.0);
```

**Mapping:**
- D^2 = |D|^2 → D2 = sum of component squares
- |D x B|^2 → P2 (explicit cross product components via BAC-CAB)
- R >= 1 always (sum of 1 + non-negative terms), so rho >= 0
- No NaN possible for any field strengths
- Poynting cross-term P2/b^4 vanishes for parallel D,B (Phase 1 standing waves)
  but is included for 3D knotted fields (Phase 3)

### Analytical BI Constitutive Relations

**Equation (PLAN.md Section 9):**
```
E = dH/dD = (D(1 + B^2/b^2) - B(B.D)/b^2) / sqrt(R)
H = dH/dB = (B(1 + D^2/b^2) - D(D.B)/b^2) / sqrt(R)
```

**Code** (bi_geon.c:241-254):
```c
double fac_D = (1.0 + B2 / b2) * inv_sqrtR;
double fac_B = (DdB / b2) * inv_sqrtR;
*ex_out = dx_f * fac_D - bx_f * fac_B;
*ey_out = dy_f * fac_D - by_f * fac_B;
*ez_out = dz_f * fac_D - bz_f * fac_B;

double fac_B2 = (1.0 + D2 / b2) * inv_sqrtR;
double fac_D2 = (DdB / b2) * inv_sqrtR;
*hx_out = bx_f * fac_B2 - dx_f * fac_D2;
*hy_out = by_f * fac_B2 - dy_f * fac_D2;
*hz_out = bz_f * fac_B2 - dz_f * fac_D2;
```

**Mapping:**
- BAC-CAB expansion already applied: B x (D x B) = D*B^2 - B*(B.D)
- fac_D = (1 + B^2/b^2)/sqrt(R) is the coefficient of D in E
- fac_B = (B.D)/(b^2 sqrt(R)) is the coefficient of B in E (subtracted)
- No Newton solver needed — purely algebraic (~30 flops per point)
- Weak-field limit: R → 1, E → D, H → B (Maxwell)
- Strong-D limit: E → D/|D| * b (capped at b)

### FFT Yukawa Solver

**Equation (PLAN.md Section 9):**
```
(nabla^2 - mu^2) Phi = -kappa rho_BI
=> Phi(k) = -kappa rho(k) / (k^2 + mu^2)
```

**Code** (bi_geon.c:277-310):
```c
double dk = 2.0 * M_PI / (N * h);
double mu2 = mu_param * mu_param;
double norm = 1.0 / (double)N3;

memcpy(fft_in, rho_BI, N3 * sizeof(double));
fftw_execute(fft_forward);

/* k-space Green's function multiplication */
double kx = ii * dk;  // ii = (i <= N/2) ? i : i - N
double ky = jj * dk;
double kz = kk * dk;
double k2 = kx*kx + ky*ky + kz*kz;

double G = kappa_param / (k2 + mu2) * norm;
fft_out[p][0] *= -G;
fft_out[p][1] *= -G;

fftw_execute(fft_backward);
memcpy(Phi, fft_in, N3 * sizeof(double));
```

**Mapping:**
- dk = 2π/(N*h) is the wavenumber spacing for the periodic box [-L, +L]
- Negative-frequency wrapping: ii = i - N for i > N/2 (FFTW convention)
- norm = 1/N^3 is FFTW's normalization (FFTW does unnormalized transforms)
- Sign: G = -kappa/(k^2 + mu^2) gives Phi < 0 for rho > 0 (attractive well)
- Zero mode (k=0): set to zero (no constant offset)
- The r2c/c2r transforms handle the Hermitian symmetry automatically

### Pade Lapse (Safeguard A)

**Equation (PLAN.md Section 5):**
```
alpha = 1 / sqrt(1 - 2*Phi)     (Pade lapse, always > 0)
```

**Code** (bi_geon.c:327-332):
```c
double arg = 1.0 - 2.0 * Phi[p];
if (arg < 0.01) arg = 0.01;  /* safety floor */
alpha_lapse[p] = 1.0 / sqrt(arg);
```

**Mapping:**
- alpha > 0 for all Phi (Pade form)
- Weak field: alpha ~ 1 + Phi (same as sqrt(1+2Phi) to O(Phi))
- Deep well: alpha → 0+ (extreme time dilation, never crosses zero)
- Safety floor at arg=0.01 → alpha_max = 10 (prevents numerical overflow)

### ADM Evolution Equations

**Equation (PLAN.md Section 9):**
```
dB/dt = -curl(alpha * E)
dD/dt = +curl(alpha * H)

curl(alpha E) = alpha curl(E) + grad(alpha) x E
                \_____________/   \_____________/
                time dilation     gravitational lensing
```

**Code** (bi_geon.c:349-393):
```c
/* Step 1: Form alpha*E and alpha*H products */
double al = alpha_lapse[p];
aEx[p] = al * ex;  aEy[p] = al * ey;  aEz[p] = al * ez;
aHx[p] = al * hx;  aHy[p] = al * hy;  aHz[p] = al * hz;

/* Step 2: Take curl of the PRODUCTS (alpha inside curl) */
curl_at(aHx, aHy, aHz, i, j, k, &ch_x, &ch_y, &ch_z);
curl_at(aEx, aEy, aEz, i, j, k, &ce_x, &ce_y, &ce_z);

rhs_dx[p] =  ch_x;   /* dD/dt = +curl(alpha*H) */
rhs_bx[p] = -ce_x;   /* dB/dt = -curl(alpha*E) */
```

**Mapping:**
- CRITICAL: alpha goes INSIDE the curl (multiply THEN differentiate)
- If alpha were OUTSIDE: curl(E) only gives time dilation, NO lensing
- The grad(alpha) x E term (from product rule) bends light toward well center
- This is the ENTIRE mechanism by which gravity confines the Hopfion
- The curl uses 4th-order central differences (same stencil as V11)

### 4th-Order Spatial Derivative

**Code** (bi_geon.c:261-268):
```c
double fm2 = f[idx(i - 2*di, j - 2*dj, k - 2*dk)];
double fm1 = f[idx(i - di,   j - dj,   k - dk)];
double fp1 = f[idx(i + di,   j + dj,   k + dk)];
double fp2 = f[idx(i + 2*di, j + 2*dj, k + 2*dk)];
return (-fp2 + 8.0*fp1 - 8.0*fm1 + fm2) / (12.0 * h);
```

**Mapping:**
- Standard 4th-order central difference: f'(x) = (-f(x+2h) + 8f(x+h) - 8f(x-h) + f(x-2h))/(12h)
- Periodic BCs via idx() modular arithmetic
- Same stencil as V11 (verified to give divB ~ O(h^4) truncation error)

### Hopfion Initialization

**Code** (bi_geon.c:154-202):
Identical to V11's `init_hopfion_B()`. Implements the Hopf map pullback:
```
B = C0 * n . (dn/dy x dn/dz, dn/dz x dn/dx, dn/dx x dn/dy)
```
where n: S^3 → S^2 is the CP^1 Hopf map, with C0 = b*a^2/16.

### RK4 Time Integration

Same 4-stage RK4 as V11, with two differences:
1. Uses analytical (D,B)→(E,H) instead of Newton inversion (~10x faster)
2. Forms alpha*E and alpha*H products before computing curl (ADM gravity)

CFL condition: dt = 0.4 * h / sqrt(3) (same as V11).
