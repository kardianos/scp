# v66 SPEC — Complexified 6-field Cosserat kernel (Q-ball sector), v1

**Date**: 2026-06-09
**Theory source**: `v66/THEORY.md` (all formulas below are copied from there; §-references
are to THEORY.md). Kernel line references are to the current working tree of
`sfa/sim/scp_sim.c`, `sfa/sim/scp_config.h`, `sfa/sim/scp_init.h`.
**Authorization**: the user has explicitly authorized modification of the simulation
kernel for this feature. Scope is limited to the files listed in §10. `sfa/format/sfa.h`
is NOT authorized for modification and this spec requires no change to it.

This spec is written for two independent implementation agents:

- **Agent K (kernel)**: sections 1–7, 10, and tests 9b–9e.
- **Agent R (radial solver)**: section 8 and test 9a (radial part). Agent R writes only
  `v66/radial_qball.c` and `v66/results/scan.tsv` + profile files; it does not touch the
  kernel. The ONLY contract between the two agents is the profile-file format in §8.3.

---

## 1. SCOPE

v1 is **CPU only**: `sfa/sim/scp_sim.c`, `sfa/sim/scp_config.h`, `sfa/sim/scp_init.h`.
The GPU port (`scp_sim.cu`) is deferred to phase 2 and MUST NOT be touched.

New opt-in config flag `complex_phi` (int, default 0):

- `complex_phi=0`: existing 6-field real kernel, **byte-for-byte unchanged behavior**
  (same memory plan, same code paths, same diag.tsv columns, same 12 SFA columns).
- `complex_phi=1`: 12-field complexified kernel per THEORY.md §1.

### 1.1 Hard incompatibilities (config validation)

Add `static void cfg_validate(const Config *c)` to `scp_config.h` (after `cfg_print`),
called from `main()` in `scp_sim.c` immediately after config load and before allocation.
When `c->complex_phi != 0`, it MUST `fprintf(stderr, "ERROR: complex_phi=1 is incompatible with %s\n", key); exit(1);`
for the FIRST violated condition among:

| condition | reason |
|---|---|
| `alpha_cs != 0` | Cosserat strain not complexified in v1 (THEORY §1.1) |
| `beta_h != 0` | hardening not complexified in v1 |
| `kappa_h != 0` | chiral helicity mixes u/v sectors, excluded (THEORY §1.1) |
| `mode != 0` | inverse / density-κ mass modes not derived for complex s |
| `self_tune != 0` | dynamical κ breaks the Q-conservation proof; Derrick saddle is cured by Q, not feedback |
| `sigma_grad != 0` | v54 theta self-interaction, not derived |
| `sigma_cubic != 0` | idem |
| `sigma_freq != 0` | idem |
| `sigma_cross != 0` | idem |
| `lambda_self != 0` | idem |
| `theta_vev != 0` | idem |
| `theta_sat != 0` | theta saturation/conversion not derived |
| `gamma_conf != 0` | mismatch confining potential not derived |
| `gamma_conv != 0` | idem |
| `sweep != 0` | sweep driver assumes real 6-field diagnostics |
| `tune_dt != 0` | v54 auto-tune (dynamical κ), same reason as self_tune |
| `vec_snap_dt != 0` | vector-frame output hardcodes 6 fields (scp_sim.c:966,1069-1085); deferred to phase 2 |

(`tune_dt` and `vec_snap_dt` are additions beyond the minimum list, included because
both are dynamical-κ / 6-field-hardcoded mechanisms; documented here as deliberate.)

Note: several of these parameters (`sigma_*`, `theta_sat`, `gamma_conf`, `gamma_conv`,
`lambda_self`, `theta_vev`, `tune_dt`) are parsed but UNUSED in the current CPU kernel
(GPU-only physics). The guard is still required: it future-proofs the phase-2 GPU port
and prevents a silently-meaningless config. `chi_chiral` is inert without `sigma_grad`
and is NOT guarded.

**Supported with complex_phi=1**: `m2` (`m`), `mtheta2` (`m_theta`), `eta`, `mu`,
`kappa`, all `bc_type` values (0 absorbing sphere, 1 gradient pinned, 2 periodic),
`bc_switch_time`, `damp_width`, `damp_rate`, snapshots (`snap_dt`, `burst_*`,
`precision`), diagnostics (`diag_dt`, `diag_file`), all init modes (behavior per §5).

## 2. FIELD LAYOUT

Keep `NFIELDS 3` (scp_config.h:15) — do NOT change it. The complex extension is a
**second copy** of the 6-field system, not a 12-component vector.

### 2.1 Naming convention (used throughout this spec and in code comments)

| THEORY.md symbol | Grid member | meaning |
|---|---|---|
| u_a   | `g->phi[a]`        | Re Φ_a (existing array, unchanged) |
| v_a   | `g->phi_im[a]`     | Im Φ_a (new) |
| tu_a  | `g->theta[a]`      | Re Θ_a (existing, unchanged) |
| tv_a  | `g->theta_im[a]`   | Im Θ_a (new) |

plus `phi_im_vel[a]`, `phi_im_acc[a]`, `theta_im_vel[a]`, `theta_im_acc[a]` mirroring
the existing `phi_vel/phi_acc/theta_vel/theta_acc`.

### 2.2 Grid struct changes (scp_sim.c:32-42)

Add to the Grid typedef:

```c
int complex_mode;                       /* copied from c->complex_phi at alloc */
double *phi_im[NFIELDS],   *phi_im_vel[NFIELDS],   *phi_im_acc[NFIELDS];
double *theta_im[NFIELDS], *theta_im_vel[NFIELDS], *theta_im_acc[NFIELDS];
double *pin_phi_im[NFIELDS], *pin_vel_im[NFIELDS];
double *pin_theta_im[NFIELDS], *pin_tvel_im[NFIELDS];
```

### 2.3 Allocation (grid_alloc, scp_sim.c:44-72)

Single flat `malloc`, block size N³ doubles:

- `complex_phi=0`: **exactly 18 blocks, identical offsets to today** (blocks 0-17 as in
  the current code). All `*_im` pointers set to NULL. This is the bit-compatibility
  guarantee: the real path allocates, indexes, and touches memory identically.
- `complex_phi=1`: 36 blocks. Blocks 0-17 keep their EXISTING assignment
  (phi 0-2, phi_vel 3-5, phi_acc 6-8, theta 9-11, theta_vel 12-14, theta_acc 15-17);
  the imaginary copy is appended:

| blocks | arrays |
|---|---|
| 18-20 | phi_im[0..2] (v) |
| 21-23 | phi_im_vel[0..2] |
| 24-26 | phi_im_acc[0..2] |
| 27-29 | theta_im[0..2] (tv) |
| 30-32 | theta_im_vel[0..2] |
| 33-35 | theta_im_acc[0..2] |

All blocks zero-initialized (use `calloc` or the existing memset pattern). Memory:
36·N³·8 B = 7.6 GB at N=288, 2.3 GB at N=192, 0.075 GB at N=64 — CPU runs in v1 are
N≤128, no problem.

**Real limit**: a `complex_phi=1` run whose v, tv, and their velocities are all zero
keeps them exactly zero under the EOM (every v/tv force term is homogeneous in the
v-sector; see §3), so it reproduces the real system up to floating-point reassociation
in the potential term (quantified tolerance in §9b).

## 3. FORCES (compute_forces_complex)

Implement a **separate function** `static void compute_forces_complex(Grid *g, const Config *c)`
placed directly after `compute_forces` (scp_sim.c:348-527). Do NOT thread `if (complex)`
branches through the existing hot loop — `compute_forces` stays untouched, which makes
the real-mode regression trivially exact. `verlet_step` dispatches:
`if (g->complex_mode) compute_forces_complex(g,c); else compute_forces(g,c);`

Since §1.1 forbids alpha_cs/beta_h/kappa_h/mode≠0, the complex function needs **no
two-pass scratch fields, no chiral block, no mode corrections** — it is the minimal
12-field loop.

### 3.1 Stencil conventions (identical to existing code, scp_sim.c:394-422)

- Voxel index `idx = (long)i*NN + j*N + k`; neighbors via periodic wrap
  `ip=(i+1)%N, im=(i-1+N)%N`, etc. (BC types 0/1 overwrite boundary behavior afterward,
  exactly as today).
- Laplacian: `lap = (f[n_ip]+f[n_im]+f[n_jp]+f[n_jm]+f[n_kp]+f[n_km] - 6*f[idx]) * idx2`,
  `idx2 = 1/(dx*dx)`.
- Curl, 2nd-order central, `idx1 = 1/(2*dx)`, reuse the existing `curl_component`
  helper on the field-pointer arrays:
  `curl(F)_0 = (F2[jp]-F2[jm]-F1[kp]+F1[km])*idx1`, cyclic.

### 3.2 Per-voxel computation (exact discrete expressions)

```c
/* amplitudes */
double u0=g->phi[0][idx], u1=g->phi[1][idx], u2=g->phi[2][idx];
double v0=g->phi_im[0][idx], v1=g->phi_im[1][idx], v2=g->phi_im[2][idx];
double s2_0 = u0*u0 + v0*v0;          /* |Phi_0|^2 */
double s2_1 = u1*u1 + v1*v1;
double s2_2 = u2*u2 + v2*v2;
double s    = s2_0 * s2_1 * s2_2;     /* THEORY eq. s = prod |Phi_a|^2 */
double den  = 1.0 + KAPPA*s;
double Vp   = 0.5*MU / (den*den);     /* Vt'(s) = (mu/2)/(1+kappa s)^2  [THEORY C2] */
/* product over b != a — computed as explicit pair products, NEVER s/s2_a */
double prod_rest[3] = { s2_1*s2_2, s2_0*s2_2, s2_0*s2_1 };
```

Then for each a = 0..2 (loops mirror scp_sim.c:468-525):

```c
/* u sector: pairs with tu via curl */
u_acc[a]  = lap(u[a])  - MASS2*u[a][idx]
          - 2.0*Vp*u[a][idx]*prod_rest[a]
          + ETA*curl(tu)[a];

/* v sector: pairs with tv via curl */
v_acc[a]  = lap(v[a])  - MASS2*v[a][idx]
          - 2.0*Vp*v[a][idx]*prod_rest[a]
          + ETA*curl(tv)[a];

/* tu sector */
tu_acc[a] = lap(tu[a]) - MTHETA2*tu[a][idx] + ETA*curl(u)[a];

/* tv sector */
tv_acc[a] = lap(tv[a]) - MTHETA2*tv[a][idx] + ETA*curl(v)[a];
```

This is THEORY.md §1 EOM verbatim (checks C3, C5, C6). The curl coupling pairs
**(u,tu)** and **(v,tv)** only — there is NO u↔tv or v↔tu coupling. OpenMP: one
`#pragma omp parallel for schedule(static)` over idx, same as the existing loop.

**Real-limit arithmetic note (binding for §9b)**: at v=0 the complex potential force is
`mu/(1+κ·u0²u1²u2²)² · u_a · (u_b² u_c²)` while the real kernel computes
`mu·(u0u1u2)/(1+κ·(u0u1u2)²)² · (u_b u_c)`. These are algebraically identical (Maxima
C3b) but NOT bit-identical (different multiplication order). Expected agreement:
~1 ULP per evaluation, i.e. ≤1e-15 relative per voxel per step.

## 4. CONFIG

### 4.1 New Config fields (scp_config.h, add after the `self_tune` block, ~line 69)

```c
/* v66 complexified Q-ball sector */
int    complex_phi;        /* 0=real 6-field (default), 1=complex 12-field */
char   qball_profile[512]; /* path to radial profile (r f) text file, init=qball */
double qball_omega;        /* internal rotation frequency omega */
double qball_x0, qball_y0, qball_z0;  /* qball center, physical coords; 1e30 = grid center */
```

### 4.2 Defaults (cfg_defaults, scp_config.h:102-134)

```c
c.complex_phi = 0;
c.qball_profile[0] = '\0';
c.qball_omega = 1.39;            /* THEORY §4: recommended first 3D target */
c.qball_x0 = 1e30;  c.qball_y0 = 1e30;  c.qball_z0 = 1e30;
```

Sentinel `1e30` means "grid center". The grid coordinate system is
`x = -L + i*dx`, `dx = 2L/(N-1)` (scp_init.h:23-26, scp_sim.c:49), so grid center is
(0,0,0); the init code replaces any coordinate ≥1e29 with 0.0.

### 4.3 Parsing (cfg_set, scp_config.h:136-213)

```c
else if (!strcmp(key,"complex_phi"))   c->complex_phi = atoi(val);
else if (!strcmp(key,"qball_profile")) strncpy(c->qball_profile, val, 511);
else if (!strcmp(key,"qball_omega"))   c->qball_omega = atof(val);
else if (!strcmp(key,"qball_x0"))      c->qball_x0 = atof(val);
else if (!strcmp(key,"qball_y0"))      c->qball_y0 = atof(val);
else if (!strcmp(key,"qball_z0"))      c->qball_z0 = atof(val);
```

### 4.4 cfg_print additions (scp_config.h:239-266)

After the `Physics:` line, add:

```c
if (c->complex_phi) {
    printf("Complex: 12-field U(1) kernel (complex_phi=1)\n");
    if (!strcmp(c->init, "qball"))
        printf("Q-ball:  omega=%.4f profile=%s center=(%.2f,%.2f,%.2f)\n",
               c->qball_omega, c->qball_profile,
               (c->qball_x0>=1e29?0.0:c->qball_x0),
               (c->qball_y0>=1e29?0.0:c->qball_y0),
               (c->qball_z0>=1e29?0.0:c->qball_z0));
}
```

The `Init:` line gains `if (!strcmp(c->init,"qball")) printf(" (%s, omega=%.4f)", c->qball_profile, c->qball_omega);`.

### 4.5 cfg_validate (new, scp_config.h after cfg_print)

Implements §1.1 plus:
- `init=qball` requires `complex_phi=1` AND non-empty `qball_profile` → error otherwise.
- `complex_phi=1` with `init=braid` → error ("braid init not defined for complex mode")
  — braid sets background phi velocities tied to the real dispersion; defining its
  complex extension is physics work outside v1 scope.

## 5. INIT

### 5.1 New mode `init=qball` (scp_init.h, new function + dispatch at lines 234-241)

Add `init_qball(Grid *g, const Config *c)` and a dispatch branch
`else if (!strcmp(c->init,"qball")) init_qball(g,c);`.

Profile loading:
1. Open `c->qball_profile` (text). Skip blank lines and lines starting with `#`.
   Each data line: two whitespace-separated doubles `r f` (extra trailing columns,
   if present, are ignored). r must be strictly increasing, first r ≥ 0.
   Store into dynamically grown arrays (realloc doubling); fail with
   `fprintf(stderr,...); exit(1)` on parse error or <2 points.
2. Report: `printf("Init: qball (%zu points, r=[%.3f,%.3f], f0=%.6f, omega=%.4f)\n", ...)`.
3. If the last tabulated f value satisfies `|f_last| > 1e-8`, print a WARNING (profile
   truncated while still sizable) but continue.

Grid fill (triple loop over i,j,k with `x=-L+i*dx` etc., same pattern as init_oscillon
scp_init.h:19-32):

```c
r  = sqrt((x-x0)^2 + (y-y0)^2 + (z-z0)^2);   /* x0,y0,z0 resolved per §4.2 */
f  = interp(r);   /* piecewise-linear in the table;
                     r < r_table[0]   -> f_table[0]   (flat core, f'(0)=0);
                     r > r_table[last]-> 0.0 */
for (a = 0..2) {
    u[a]      = f;                    /* g->phi[a][idx] */
    v[a]      = 0.0;                  /* g->phi_im[a][idx] */
    udot[a]   = 0.0;                  /* g->phi_vel[a][idx] */
    vdot[a]   = c->qball_omega * f;   /* g->phi_im_vel[a][idx] */
    tu[a] = tv[a] = tudot[a] = tvdot[a] = 0.0;
}
```

This is exactly Φ_a(t=0) = f(r), ∂_t Φ_a(0) = iωf(r) — the THEORY §3 ansatz at t=0,
giving Q = 3ω∫f²dV > 0 with the THEORY §2 sign convention. **`delta[]` is IGNORED by
init=qball** (the symmetric ansatz uses identical phase for all a; THEORY §5c shows
phase staggering buys nothing). Document this in the init printf.

### 5.2 Existing init modes under complex_phi=1 (v1 behavior — decided, documented)

- `init=oscillon`: real sector filled exactly as today (scp_init.h:19-32); the entire
  imaginary sector (v, tv, and their velocities) stays zero. This is the §9b
  regression-test configuration.
- `init=sfa` / `init=exec`: the column dispatch (scp_init.h:188-196) is extended with
  four new slot rules, active only when `complex_phi=1`:

  | semantic | component | target | slot |
  |---|---|---|---|
  | SFA_POSITION | 3-5 | `phi_im[comp-3]` | 12-14 |
  | SFA_ANGLE | 3-5 | `theta_im[comp-3]` | 15-17 |
  | SFA_VELOCITY | 6-8 | `phi_im_vel[comp-6]` | 18-20 |
  | SFA_VELOCITY | 9-11 | `theta_im_vel[comp-9]` | 21-23 |

  A legacy real 12-column file therefore loads real parts AND real velocities exactly
  as today; all imaginary slots remain zero. A v66 24-column file (per §7) round-trips
  fully. When `complex_phi=0` and a file contains the new components, skip them with
  a one-line WARNING. Missing-velocity warning (scp_init.h:115-116) unchanged.
- `init=template`: allowed; background + template stamping operate on the real sector
  exactly as today (scp_init.h:140-232); imaginary sector zero. (A real template in a
  complex run is a valid "real limit" seed.)
- `init=braid`: ERROR under complex_phi=1 (§4.5).

## 6. DIAGNOSTICS

### 6.1 New quantities (computed only when complex_phi=1)

Add `static void compute_charges(Grid *g, double *Q_phi, double *Q_theta, double *s_max, double *r_core)`
placed after `P_integrated` (scp_sim.c:726-731). One OpenMP reduction pass:

```
Q_phi   = sum_idx sum_a ( u_a*vdot_a - v_a*udot_a ) * dV          [THEORY §2]
Q_theta = sum_idx sum_a ( tu_a*tvdot_a - tv_a*tudot_a ) * dV
Q_total = Q_phi + Q_theta                                          (computed at call site)
s_max   = max_idx s(idx),  s = (u0²+v0²)(u1²+v1²)(u2²+v2²)
```

r_core (rms charge radius, matches THEORY table column r_Q(rms)): with amplitude
density `rho2(idx) = Σ_a (u_a²+v_a²)` (sign-definite; equals 3f² on the ansatz, so its
rms radius = r_Q since ρ_Q = ω·rho2),

```
M0 = Σ rho2*dV;   cx = Σ x*rho2*dV / M0  (likewise cy, cz)
r_core = sqrt( Σ ((x-cx)²+(y-cy)²+(z-cz)²) * rho2 * dV / M0 )
```

Guard: if `M0 < 1e-30`, set r_core = 0. (Centroid-based, so it tracks a drifting ball;
caveat: a large background inflates it — fine for Q-ball runs where background ≈ 0.)
This needs a second pass (centroid first) — acceptable, diagnostics run every
`diag_every` steps only. Velocities are full-step synchronized at diagnostic time
(diag runs after the second Verlet half-kick), so Q uses consistent (x, ẋ) pairs.

### 6.2 diag.tsv columns

`complex_phi=0`: header and rows **byte-identical to today** (scp_sim.c:1009-1010):

```
t  E_phi_kin  E_theta_kin  E_grad  E_mass  E_pot  E_tgrad  E_tmass  E_coupling  E_total  phi_max  P_max  P_int  theta_rms
```

`complex_phi=1`: **append 5 columns after theta_rms** (existing 14 keep position and
meaning, so column-index-based parsers keep working):

```
... theta_rms  Q_phi  Q_theta  Q_total  s_max  r_core
```

Edit points: the header fprintf at scp_sim.c:1009-1010 gains a conditional
`if (c.complex_phi) fprintf(fp, "\tQ_phi\tQ_theta\tQ_total\ts_max\tr_core");` before the
`\n` (restructure the single fprintf into base + conditional + newline). The three row
writers — initial (lines 1025-1034), in-loop (lines 1163-1178), final (lines 1202-1210)
— get the same conditional append, calling `compute_charges` once per diagnostic event.
Console output: append `Q=%.6e s_max=%.4e` to the in-loop printf when complex.

### 6.3 Existing column semantics under complex_phi=1

Computed in a new `compute_energy_complex` (mirror of compute_energy, scp_sim.c:637-716;
the real function stays untouched — same dispatch pattern as §3):

| column | complex definition | real limit |
|---|---|---|
| E_phi_kin | Σ_a ∫ ½(u̇_a²+v̇_a²) dV | identical |
| E_theta_kin | Σ_a ∫ ½(tu̇_a²+tv̇_a²) dV | identical |
| E_grad | Σ_a ∫ ½(\|∇u_a\|²+\|∇v_a\|²) dV | identical |
| E_mass | Σ_a ∫ ½m²(u_a²+v_a²) dV | identical |
| E_pot | ∫ Vt(s) dV = ∫ (μ/2)s/(1+κs) dV | = ∫ V(P) dV [THEORY C1] |
| E_tgrad, E_tmass | analogous with tu,tv | identical |
| E_coupling | −η Σ_a ∫ (u_a·curl(tu)_a + v_a·curl(tv)_a) dV | identical (matches existing sign, scp_sim.c:707-710) |
| E_total | sum of the 8 above | identical |
| phi_max | max_{a,idx} sqrt(u_a²+v_a²) (complex modulus — time-stationary on a Q-ball) | = max\|φ_a\| up to sqrt/square rounding |
| P_max | max_idx sqrt(s) (define P := √s = Π\|Φ_a\| ≥ 0) | = max\|P\| up to rounding |
| P_int | ∫ sqrt(s) dV | = ∫\|P\| dV up to rounding |
| theta_rms | sqrt( Σ_{a,idx}(tu²+tv²) / (6·N³) ) — **denominator 6**, not the hardcoded 3 (scp_sim.c:723) | = real value (tv=0 contributes 0, but /6 ≠ /3 — see note) |

**Note (deliberate)**: theta_rms in complex mode divides by 6 components even when
tv≡0, so its real-limit value is 1/√2 × the real-mode value. This is correct RMS over
the actual 6 theta components and is NOT an energy column; §9b compares energy columns
only. The smoking-gun observable (THEORY §6) is that P_max/phi_max are time-stationary
for a Q-ball (no 3ω harmonic), which these modulus definitions make directly visible.

## 7. OUTPUT (SFA)

No change to `sfa/format/sfa.h` (not authorized; none needed — `component` is uint8,
column count is uint32, names are ≤11 chars).

### 7.1 Column design (single frame group, 24 columns, complex_phi=1)

Columns 0-11: **identical names, semantics, components, and order as today**
(scp_sim.c:940-951) carrying the REAL parts — every existing analysis tool that looks
up `phi_x`/`SFA_POSITION comp<3` etc. keeps working on complex files. Columns 12-23
are appended:

| col | name | semantic | comp | array |
|---|---|---|---|---|
| 12-14 | `phiim_x` `phiim_y` `phiim_z` | SFA_POSITION | 3,4,5 | phi_im[0..2] |
| 15-17 | `thetaim_x` `thetaim_y` `thetaim_z` | SFA_ANGLE | 3,4,5 | theta_im[0..2] |
| 18-20 | `phiim_vx` `phiim_vy` `phiim_vz` | SFA_VELOCITY | 6,7,8 | phi_im_vel[0..2] |
| 21-23 | `thetaim_vx` `thetaim_vy` `thetaim_vz` | SFA_VELOCITY | 9,10,11 | theta_im_vel[0..2] |

(Longest name `thetaim_vx` = 10 chars < 12-char field. Re-uses only existing semantic
codes; component offset distinguishes imaginary parts, mirroring how theta velocities
already use SFA_VELOCITY comp 3-5.) When `complex_phi=0`, write exactly the existing
12 columns — output files are byte-compatible with today.

### 7.2 sfa_snap (scp_sim.c:821-847)

Make the pointer table and loop bounds count-driven:

```c
int nf = g->complex_mode ? 24 : 12;
double *arrays[24] = { phi[0..2], theta[0..2], phi_vel[0..2], theta_vel[0..2],
                       phi_im[0..2], theta_im[0..2], phi_im_vel[0..2], theta_im_vel[0..2] };
```

(first 12 entries in the existing order; the f32/f16 temp-buffer loops at lines 839-845
run to `nf`). Array order MUST match column registration order in §7.1.

### 7.3 Column registration (scp_sim.c:940-952)

After the existing 12 `sfa_add_column` calls, add the 12 new ones from §7.1 inside
`if (c.complex_phi) { ... }`, then `sfa_finalize_header` as today. Vector frames
(FRVD) are guarded off by §1.1 (`vec_snap_dt!=0` errors), so lines 966 and 1069-1085
need no change in v1.

## 8. RADIAL SOLVER — `v66/radial_qball.c`

Standalone C tool, no SFA dependency (1D text output is policy-compliant: SFA applies
to 3D field data). Build: `gcc -O3 -o v66/radial_qball v66/radial_qball.c -lm`.
Uses `long double` throughout the shooting integration (per THEORY §4 the bisection
hits the long-double precision floor near thin wall). Structure cloned from
`v65/radial_oscillon.c` (CLI pattern, diag cadence) but the integrator is RK4 shooting
on the STATIC profile ODE, not time evolution.

### 8.1 Equations (THEORY §3, verified C6-C8)

```
f'' + (2/r) f' = U'(f),
U'(f)  = (m² − ω²) f + 2·Vt'(f⁶)·f⁵,      Vt'(s) = (μ/2)/(1+κs)²
U(f)   = ½(m² − ω²) f² + (μ/6) f⁶/(1+κ f⁶)
```

Existence window (closed form, used for defaults and sanity checks):
`ω_min² = m² + (μ/9)(2/κ)^{2/3}`; standard parameters → ω ∈ (1.3087, 1.5000).

### 8.2 CLI

```
./radial_qball -omega 1.39 [options]                    # single profile
./radial_qball -omega_min 1.315 -omega_max 1.495 -omega_step 0.005 [options]   # scan
```

| flag | default | meaning |
|---|---|---|
| `-omega W` | — | single-ω mode; writes profile + one scan row |
| `-omega_min/-omega_max/-omega_step` | — | scan mode; profile written per ω |
| `-m2 X` | 2.25 | m² |
| `-mu X` | -41.345 | μ |
| `-kappa X` | 50.0 | κ |
| `-rmax X` | 60.0 | profile output extent |
| `-dr X` | 0.001 | RK4 step AND profile output spacing (output may be coarser, see `-out_dr`) |
| `-out_dr X` | 0.02 | profile file sampling interval |
| `-f0_lo X` | auto | bisection bracket low; auto = f₁·(1+1e-6), f₁ = smallest positive zero of U(f) (bisected numerically) |
| `-f0_hi X` | auto | bracket high; auto = f_dip = argmin U(f) (golden-section, or closed form κf⁶ near 2 as seed) |
| `-iters N` | 120 | max bisection iterations (stop early when bracket width < 1e-17·f0) |
| `-profile PATH` | `v66/results/profile_omega<W>.txt` (W with 4 decimals, e.g. `profile_omega1.3900.txt`) | profile output |
| `-scan PATH` | `v66/results/scan.tsv` | scan table output (appends header only if file absent) |

Errors: ω outside (ω_min, m) → stderr message with the computed window, exit 1.
The tool creates `v66/results/` (mkdir) if missing.

### 8.3 Method

1. **Series start** (r=0 regularity, f'(0)=0): take first step from r=0 to r=dr with
   `f(dr) = f0 + (dr²/6)·U'(f0)`, `f'(dr) = (dr/3)·U'(f0)`.
2. **RK4 integrate** the system (f, p=f') with the (2/r)p damping term to r_stop = 80
   (or until a bracket event). Bracket classification:
   - f crosses 0 → **overshoot** (f0 too large… or too small depending on branch; use
     the standard undamped-particle analogy: rolling from f0 in potential −U, crossing
     f=0 means too much initial "energy" → f0 too far past f₁; classify by which side
     produces which and bisect accordingly — match `theory/qball_checks.py` logic,
     under/over bisection).
   - f' becomes ≥ 0 while f > 0 (turns back up) → **undershoot**.
3. Bisect f0 within the bracket for `-iters` iterations (long double).
4. **Tail attach**: at r_t = the last radius where the solution is still monotonically
   decaying and `f < 1e-3·f0` (or where bracket noise exceeds f), replace the remainder
   with the analytic tail `f(r) = C·e^{−βr}/r`, `β = sqrt(m²−ω²)`, C matched to f(r_t).
   Record relative ODE residual over the kept numerical segment; report it.
5. **Integrals** (composite Simpson on the assembled profile, 4π r² dV weight):
   ```
   Q = 3ω ∫ f² 4πr² dr
   E = ∫ [ (3/2)(ω²f² + f'² + m²f²) + (μ/2) f⁶/(1+κf⁶) ] 4πr² dr     [THEORY C10,C11]
   r_half: largest r with f(r) ≥ f0/2
   r_Q   = sqrt( ∫ r²·f²·4πr²dr / ∫ f²·4πr²dr )
   ```

### 8.4 Outputs

**Profile file** (the init=qball contract, §5.1):

```
# radial_qball profile: omega=1.390000 m2=2.250000 mu=-41.345000 kappa=50.000000
# f0=0.640512  E=691.9  Q=482.2  E/(m*Q)=0.9566  residual=3.2e-05
# r f
0.000000 0.640512
0.020000 0.640509
...
60.000000 0.000000
```

Two whitespace-separated columns r f, `#` comments, r strictly increasing from 0.0 to
rmax inclusive, sampled every out_dr (linear interpolation of the RK4 trace + analytic
tail). Default rmax=60 covers any box up to L=60.

**Scan table** `v66/results/scan.tsv` (TSV, one header line):

```
omega	f0	E	Q	E_over_mQ	dQ_domega_sign	r_half	r_Q	residual
1.3900	0.640512	691.9	482.2	0.9566	-1	4.41	3.71	3.2e-05
```

- `E_over_mQ` = E/(m·Q) with m = sqrt(m2). Stability: classical needs
  `dQ_domega_sign = -1`; absolute needs `E_over_mQ < 1`.
- `dQ_domega_sign` = sign of centered finite difference of Q across adjacent scan rows
  (forward/backward at the ends; in single-ω mode compute Q at ω±0.002 internally).
- Sanity gates (stderr WARN, nonzero exit only on hard failure): ω=1.39 row should
  reproduce THEORY §4 within ~5%: f0≈0.6405, Q≈482, E≈692, r_half≈4.4.

## 9. TEST PLAN

All commands from repo root `/home/d/code/scp`. Test configs live in `v66/cfg/`,
outputs in `v66/results/` (small) — diag.tsv files are small; any .sfa goes to
`/space/scp/v66/` if >100 MB (these tests are N=64, so local is fine).

### (a) Build

```bash
gcc -O3 -march=native -fopenmp -o bin/scp_sim sfa/sim/scp_sim.c -lzstd -lm   # kernel
gcc -O3 -o v66/radial_qball v66/radial_qball.c -lm                            # radial
# regression that the unmodified-real path still builds via the canonical target:
make -C sfa install
```

Radial smoke: `./v66/radial_qball -omega 1.39` must emit the profile and a scan row
passing the §8.4 sanity gates. Also run the scan
`./v66/radial_qball -omega_min 1.315 -omega_max 1.495 -omega_step 0.005` and check
dQ_domega_sign flips from −1 to +1 somewhere in ω∈(1.46,1.50) and E_over_mQ crosses
1 near ω≈1.443 (THEORY §4).

### (b) Real-limit regression

Config `v66/cfg/reg_real.cfg`:

```
N=64
L=15
T=20
dt_factor=0.025
init=oscillon
A=0.8
sigma=3.0
eta=0.5
bc_type=0
diag_dt=1.0
snap_dt=0
diag_file=v66/results/reg_real_diag.tsv
output=/tmp/reg_real.sfa
precision=f32
complex_phi=0
```

`reg_cplx.cfg`: identical except `complex_phi=1`,
`diag_file=v66/results/reg_cplx_diag.tsv`, `output=/tmp/reg_cplx.sfa`. (oscillon init
zeroes the entire imaginary sector per §5.2.)

PASS criteria on the 9 energy columns (E_phi_kin … E_total):
- t=0 row: relative difference ≤ 1e-13 per column (same fields, reassociated potential
  formula only).
- all rows t ≤ 20: relative difference ≤ 1e-9 per column.

**Documented tolerance rationale**: bit-exactness (<1e-12 over the whole run via "same
arithmetic") is NOT claimed — the complex potential force evaluates
`(μ/2)/(1+κ·u0²u1²u2²)²·2u_a·u_b²u_c²` vs the real kernel's
`μP/(1+κP²)²·u_bu_c`, P=u0u1u2: algebraically equal, reassociated in FP (~1 ULP/step).
If 1e-9 at t=20 fails due to ULP Lyapunov growth in the strongly nonlinear oscillon,
the fallback gate is ≤1e-9 for t ≤ 5 plus monotonically-bounded divergence (no jump);
record the measured max divergence in v66/RESULTS. Also check: `theta_rms` complex
column = real column / √2 (×(1±1e-9)), and `P_max/P_int/phi_max` agree to ≤1e-9
(modulus-vs-abs rounding).

### (c) Q-conservation (integrator drift floor, THEORY §2 boundary caveat)

```
# v66/cfg/qcons.cfg
N=64
L=15
T=20
dt_factor=0.025
complex_phi=1
init=qball
qball_profile=v66/results/profile_omega1.3900.txt
qball_omega=1.39
eta=0
bc_type=2
diag_dt=0.5
snap_dt=0
diag_file=v66/results/qcons_diag.tsv
output=/tmp/qcons.sfa
precision=f32
```

(bc_type=2 periodic: no sponge, no boundary flux; box [-15,15] vs r_half≈4.4 — tails
~e^{−0.56·15} ≈ 2e-4, wrap negligible over T=20.) PASS: `max_t |Q_total(t) − Q_total(0)| / |Q_total(0)| < 1e-6`,
and Q_total(0) within 10% of the radial solver's Q≈482 (lattice-vs-continuum
discretization). Also record E_total drift (Verlet floor) for the same run.

### (d) Q-ball persistence smoke test (η=0, absorbing BC)

`v66/cfg/persist.cfg`: same as (c) but `T=50`, `bc_type=0`, `damp_width=3.0`,
`damp_rate=0.01`, `snap_dt=10`, `diag_file=v66/results/persist_diag.tsv`,
`output=/space/scp/v66/persist.sfa` (f32; ~50 MB at N=64, local /tmp also acceptable).
PASS (qualitative gates):
- s_max(t) stationary: stddev/mean < 5% over t∈[10,50] (no 3ω P_max oscillation — the
  THEORY §6 smoking gun vs v65 oscillons);
- Q_phi(t) retains > 99% over T=50 (only discretization radiation reaches the sponge);
- r_core stationary within ±10%;
- E_total/Q_phi ≈ ω·(E/Q)_theory consistency: E/Q ≈ 1.435 ± 5%.

### (e) η=0.5 drain measurement (the science run)

`v66/cfg/drain.cfg`: same seed as (d) but `eta=0.5`, `T=200`, `diag_dt=1.0`. Outputs to
`v66/results/drain_diag.tsv`. Measured quantities (analysis, not pass/fail):
dQ_phi/dt and dE/dt time series, Q_theta(t) (charge transferred to Θ before leaving),
charge half-life t_{1/2}, and a fit of drain rate vs η² (rerun at eta=0.25 for the
∝η² leading-order check, THEORY §5b). Compare s_max behavior: settles to finite-η
quasi-stationary attractor vs evaporates.

## 10. FILE-BY-FILE CHANGE LIST

### `sfa/sim/scp_config.h`

| location | edit |
|---|---|
| Config struct (after line 69, st_gain) | add `complex_phi`, `qball_profile[512]`, `qball_omega`, `qball_x0/y0/z0` (§4.1) |
| `cfg_defaults` (102-134) | defaults per §4.2 |
| `cfg_set` (136-213) | six new keys per §4.3 |
| `cfg_print` (239-266) | complex/qball lines per §4.4 |
| new `cfg_validate` (after cfg_print) | §1.1 incompatibility table + §4.5 init checks; signature `static void cfg_validate(const Config *c)` |

NFIELDS (line 15) unchanged. f16 helpers unchanged.

### `sfa/sim/scp_sim.c`

| function (lines) | edit |
|---|---|
| Grid struct (32-42) | add `complex_mode` + 12 `*_im` pointers + 4 `pin_*_im` pointer arrays (§2.2) |
| `grid_alloc` (44-72) | 18 vs 36 block allocation, real blocks at unchanged offsets, im pointers NULL when real (§2.3); set `g->complex_mode = c->complex_phi` (pass Config or flag in) |
| `grid_save_pinned` (84-95) | if complex: also allocate+memcpy pin_phi_im, pin_vel_im, pin_theta_im, pin_tvel_im |
| `compute_forces` (348-527) | **UNTOUCHED** |
| new `compute_forces_complex` (insert after 527) | §3.2 minimal 12-field loop, same stencils, no scratch/chiral/mode blocks |
| `apply_damping` (534-551) | inside the existing per-voxel damping, add `if (g->complex_mode)` block damping `phi_im_vel[a]` and `theta_im_vel[a]` by the same factor d |
| `apply_gradient_bc` (558-603) | x-slab pinning and accel-zeroing extended to the im arrays (using pin_*_im); y-extrapolation loop extended to phi_im/theta_im; guarded by complex_mode |
| `verlet_step` (609-631) | half-kick / drift / half-kick loops gain a `if (g->complex_mode)` twin loop over (phi_im, phi_im_vel, phi_im_acc, theta_im, theta_im_vel, theta_im_acc); force dispatch real/complex (§3) |
| `compute_energy` (637-716) | **UNTOUCHED**; new `compute_energy_complex` after it implementing §6.3 (8 energies + phi_max + P_max via √s) |
| `theta_rms` (718-724) | UNTOUCHED; new `theta_rms_complex` with /(6·N³) (§6.3) |
| `P_integrated` (726-731) | UNTOUCHED; new `P_integrated_complex` = ∫√s dV |
| new `compute_charges` (after 731) | §6.1: Q_phi, Q_theta, s_max, two-pass centroid r_core |
| `sfa_snap` (821-847) | arrays[24], `nf = complex ? 24 : 12`, loops to nf (§7.2) |
| `main` (853-…) | call `cfg_validate(&c)` after cfg_load/cfg_print; 12 extra `sfa_add_column` calls when complex (after 951, §7.3); diag header conditional append (1009-1010, §6.2); diag rows: dispatch energy/rms/Pint to *_complex variants and append charge columns at the initial (1025-1034), in-loop (1163-1178), and final (1202-1210) diagnostic sites |

Vector-frame code (lines 966, 1069-1085) untouched — unreachable under complex (guard).

### `sfa/sim/scp_init.h`

| location | edit |
|---|---|
| `init_oscillon` (19-32) | untouched (im sector calloc'd zero) |
| `init_braid` (34-56) | untouched; cfg_validate rejects it under complex |
| `init_from_sfa` (58-122) | extend slot dispatch (around 188-196 in scp_sim.c numbering — the semantic/component switch inside init_from_sfa) with the 4 imaginary rules of §5.2, active iff complex; warn+skip when real mode meets imaginary columns |
| `init_template` (140-232) | untouched (real sector only) |
| new `init_qball` (insert before do_init) | §5.1: profile load, linear interp, ansatz fill |
| `do_init` (234-241) | add `qball` dispatch branch |

### New files

- `v66/radial_qball.c` (§8) — Agent R.
- `v66/cfg/reg_real.cfg`, `reg_cplx.cfg`, `qcons.cfg`, `persist.cfg`, `drain.cfg` (§9).
- Generated, not committed: `v66/results/profile_omega*.txt`, `v66/results/scan.tsv`,
  `v66/results/*_diag.tsv`.

### Explicitly NOT modified

`sfa/format/sfa.h` (not authorized; not needed), `sfa/sim/scp_sim.cu` (phase 2),
seed generators in `sfa/seed/`, all analysis tools.

---

## Cross-agent contract recap

Agent R produces `v66/results/profile_omega1.3900.txt` in the §8.4 format; Agent K's
`init=qball` consumes exactly that format per §5.1. Neither agent needs anything else
from the other. Both must use Vt'(s) = (μ/2)/(1+κs)² with the SAME sign conventions
(μ negative, force −2Vt'(s)·f·∂s) — both formulas are copied verbatim from THEORY.md
§1/§3, which is the single source of truth on any discrepancy.
