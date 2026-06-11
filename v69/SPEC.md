# v69 SPEC — Gauged diagonal U(1) kernel (kernel-v3, `complex_gauge`)

**Date**: 2026-06-10
**Theory source**: `v68/GAUGE_DESIGN.md` (§2 Lagrangian/EOM/Gauss, §4 lattice scheme,
§5 predictions, §6 architecture sketch; Maxima `v68/theory/gauge_checks.mac` 59/59 PASS).
Extends the v66 complexified kernel (`v66/THEORY.md`, `v66/SPEC.md`) and the v67
diagnostics. Kernel line references are to the current working tree of
`sfa/sim/scp_sim.c` (1620 lines), `sfa/sim/scp_config.h` (374), `sfa/sim/scp_init.h` (402).
**Authorization**: the user has explicitly authorized kernel modification for this
feature (kernel-v3). Scope: the three files above plus new files in `v69/`.
`sfa/format/sfa.h` is NOT modified (none needed: semantics ANGLE/VELOCITY with
component offsets suffice; names ≤ 11 chars). `sfa/sim/scp_sim.cu` is NOT touched
(GPU port is a separate later phase, §7.6).

**Conventions** (fixed throughout; GAUGE_DESIGN §2 with one sign repair, see Open
Questions O1): D_μ = ∂_μ + igA_μ on charge-(+1) fields (both Φ_a and Θ_a);
temporal gauge A_0 = 0; **E_i = −∂_t A_i** (and in general E = −∇A_0 − Ȧ);
Gauss law ∇·E = +g ρ_Q with the v66 charge density
ρ_Q = Σ_a (u_a v̇_a − v_a u̇_a + tu_a tv̇_a − tv_a tu̇_a).
Compact Kogut–Susskind links: angles θ_i(x) = g·a·A_i(x+î/2) wrapped to (−π,π],
noncompact conjugate E_i(x) on the same link. a ≡ dx. Recommended defaults
g = 0.05, m_θ = 1.6 (full package); η = 0 stage-0 validation.

---

## 1. SCOPE

v1 is **CPU only**. New opt-in config flag `complex_gauge` (int, default 0) with
coupling `g_gauge` (double, default 0.05).

- `complex_gauge=0`: kernel behavior **byte-for-byte unchanged** for both
  `complex_phi=0` and `complex_phi=1` (no gauge allocation, no gauge code paths,
  no new diag/SFA columns — same dispatch structure as v66 SPEC §1).
- `complex_gauge=1`: gauged 12-field kernel per GAUGE_DESIGN §2/§4.

### 1.1 Hard incompatibilities (cfg_validate extensions)

When `c->complex_gauge != 0`, error (`fprintf(stderr,"ERROR: complex_gauge=1 ...");
exit(1)`) on the FIRST violated condition:

| condition | reason |
|---|---|
| `complex_phi == 0` | gauge sector is defined only for the complex 12-field theory |
| (entire v66 SPEC §1.1 table) | inherited automatically via the existing complex_phi guard (alpha_cs, beta_h, kappa_h, mode, self_tune, sigma_*, lambda_self, theta_vev, theta_sat, gamma_conf, gamma_conv, sweep, tune_dt, vec_snap_dt) |
| `bc_type == 1` | gradient-pinned BC + links not derived (pinning A is gauge-dependent) |
| `bc_switch_time != 0` | BC switching with a live constraint field not derived in v1 |
| `g_gauge < 0` | sign convention is fixed; flip charges via the seed, not g |
| `init == "braid"` | already refused under complex_phi |

Supported: `bc_type` 0 (absorbing) and 2 (periodic), all snapshot/diag options,
`init` ∈ {oscillon, sfa, exec, template, qball} (behavior per §5), `eta`, `m_theta`.

### 1.2 Net-charge check (runtime, after init — not in cfg_validate)

A periodic lattice (and the wrapped stencil of bc_type=0, see §4.1) admits no net
flux, so **bc_type=2 with net charge is refused**. In `main()` immediately after
`do_init()`, when `complex_gauge && g_gauge != 0 && bc_type == 2`: one reduction pass
computes `Q_net = Σ_x ρ_Q(x)·dV` and `Q_abs = Σ_x |ρ_Q(x)|·dV`; if
`|Q_net| > 1e-6 * fmax(Q_abs, 1.0)` →
`ERROR: bc_type=2 (periodic) requires net-neutral seed (Q_net=%.3e, Q_abs=%.3e); use +/- pairs or bc_type=0. jellium=1 is deferred (not implemented in v1).` and exit(1).
(v1: just error. A future `jellium=1` flag would instead document the fictitious
uniform background; the projection pass §5.4 already produces exactly that field,
so the flag is a one-line unlock later.)

### 1.3 g = 0 identity guarantee

When `complex_gauge=1 && g_gauge == 0.0` the kernel **dispatches to the existing
ungauged complex path** (`compute_forces_complex`, `compute_energy_complex`); links
and E are allocated, initialized to zero, and never updated (no E-kick, no θ-drift,
no projection pass); gauge diag columns are emitted as exact zeros (`Q_flux` = 0 by
the g==0 guard). **Expected tolerance: bit-identical** diag energy/charge columns and
matter SFA columns vs the same config with `complex_phi=1, complex_gauge=0` — it is
the same code path, not a numerically-reducing one.
The covariant-stencil reduction is validated separately at tiny g (test §7.1b):
with g = 1e-8 the gauged stencils run with θ ≡ O(g²)~1e-16, cos θ = 1.0 and
sin θ ≤ 1e-16 in f64, physics shifts are O(g²) ≈ 1e-16, and the remaining
difference is multiplication-order reassociation (~1 ULP/step): gate
**≤ 1e-9 relative per energy column at T = 20** (same rationale and fallback ladder
as v66 SPEC §9b).

---

## 2. FIELDS AND MEMORY

### 2.1 Precision decision: f64 gauge sector

GAUGE_DESIGN §4 sketched +6 N³ **f32** arrays (GPU framing) and flagged the
precision risk (§4 CFL note + open risk #4): the per-step link increment is
|Δθ| = g·a·E·dt ≈ 0.05·0.23·0.03·0.006 ≈ 2×10⁻⁷ rad, below f32 resolution at
θ ~ 0.1 (ε_f32·0.1 ≈ 1.2×10⁻⁸ is only 20× smaller — accumulation over 10⁵–10⁶
steps loses the increment). **Decision: the entire gauge sector is f64**, matching
the CPU matter sector (the kernel's 36 complex-mode blocks are already f64,
288 B/voxel — the "match the matter precision" choice is f64, not f32).

### 2.2 New Grid members and allocation layout

Grid struct additions (scp_sim.c:32–48):

```c
int    gauge_mode;            /* copied from c->complex_gauge at alloc */
double g_gauge;               /* cached coupling */
double G_offset;              /* frozen uniform Gauss offset, measured post-init (§4.1) */
double *th[3];                /* link angles theta_i(x) on link (x, x+i), wrapped (-pi,pi] */
double *Efield[3];            /* E_i(x) on the same link (noncompact) */
double *E_acc[3];             /* E-kick K_i (the gauge "acceleration", §3.4) */
double *link_c[3], *link_s[3];/* scratch: cos/sin theta_i, refreshed per force call */
double *plaq_s[3];            /* scratch: sin theta_P per plaquette plane (§3.4) */
```

`grid_alloc` (scp_sim.c:50–95) extends the single flat slab:

| complex_phi | complex_gauge | blocks | layout |
|---|---|---|---|
| 0 | 0 | 18 | unchanged (bit-compat) |
| 1 | 0 | 36 | unchanged v66 layout (blocks 0–35) |
| 1 | 1 | **54** | blocks 0–35 as v66; then th 36–38, Efield 39–41, E_acc 42–44, link_c 45–47, link_s 48–50, plaq_s 51–53 |

All zero-initialized. Gauge pointers NULL when `complex_gauge=0`. Memory:
54·8 = 432 B/voxel = 1.5× the v66 complex footprint (dynamic state is +9 blocks
= +25%; the other +9 are scratch). N=128: 0.91 GB; N=160: 1.77 GB; N=192: 3.06 GB
— all fine for CPU v1 (runs are N ≤ 160).

Link/site convention: `th[i][idx]` and `Efield[i][idx]` live on the link from voxel
`idx` to its +î neighbor (forward link). Periodic index wrap exactly as the matter
stencils; bc_type=0 keeps the wrapped stencil and damps in the sponge (§4.1/§3.6).

---

## 3. DYNAMICS — exact discrete expressions

Everything in this section is derived from ONE discrete Hamiltonian (per volume a³):

```
H = Σ_x a³ { ½Σ_a(u̇²+v̇²+tu̇²+tv̇²)
           + (1/2a²) Σ_i Σ_a [ |Φ_a(x+î)|² + |Φ_a(x)|² − 2 Re(Φ̄_a(x) U_i(x) Φ_a(x+î)) ]   (Φ and Θ alike, m²/m_θ² masses)
           + (m²/2)Σ|Φ_a|² + (m_θ²/2)Σ|Θ_a|² + Vt(s)
           − η Σ_a Re[ Θ̄_a(x) (D×Φ)^lat_a(x) ]                                            (covariant central curl, §3.3)
           + ½ Σ_i E_i(x)²  +  (1/(g²a⁴)) Σ_{p∈{xy,yz,zx}} (1 − cos θ_p(x)) }
```

with U_i(x) = e^{iθ_i(x)} and plaquette angle
`θ_ij(x) = θ_i(x) + θ_j(x+î) − θ_i(x+ĵ) − θ_j(x)`.
H is **exactly invariant** under arbitrary lattice gauge transformations
Φ→e^{−igα(x)}Φ, Θ→e^{−igα(x)}Θ, θ_i(x)→θ_i(x)+g[α(x+î)−α(x)] — every matter
difference (Laplacian AND η curl) uses link transport. That exact invariance is what
makes the discrete Gauss law a constant of the symplectic update (§3.5). **Do not
mix raw and transported stencils anywhere in the gauged force/current code**
(GAUGE_DESIGN risk #3; the tripwire is gauss_max, §4).

Canonical pairs: (Φ, Φ̇), (Θ, Θ̇), and (A_i, π_i = Ȧ_i = −E_i) with θ_i = g·a·A_i.
Hence: matter kicks from −∂H/∂(fields); **E-kick Ė_i = +(1/a³)·∂H/∂A_i(x)
= (g·a/a³)·∂H/∂θ_i(x)**; link drift θ̇_i = g·a·Ȧ_i = −g·a·E_i.

### 3.1 Transported neighbor values (the only new stencil primitive)

For any complex field F = (fu, fv) (one of Φ_a or Θ_a — same charge, same links),
direction i, voxel x, with `c⁺ = link_c[i][x]`, `s⁺ = link_s[i][x]`,
`c⁻ = link_c[i][x−î]`, `s⁻ = link_s[i][x−î]` and neighbors F⁺ = F(x+î), F⁻ = F(x−î):

```
/* forward:  U_i(x) · F(x+i)        (U·Phi rotates by +theta) */
TRp = c⁺*fu⁺ − s⁺*fv⁺ ;      TIp = c⁺*fv⁺ + s⁺*fu⁺
/* backward: U_i(x−i)† · F(x−i)     (conjugate link: rotate by −theta) */
TRm = c⁻*fu⁻ + s⁻*fv⁻ ;      TIm = c⁻*fv⁻ − s⁻*fu⁻
```

At θ = 0 these are exactly the raw neighbor values.

### 3.2 Covariant Laplacian (forces, both sectors)

```
LapU(F)(x) = (1/a²) Σ_{i=0..2} [ TRp_i + TRm_i − 2*fu(x) ]
LapV(F)(x) = (1/a²) Σ_{i=0..2} [ TIp_i + TIm_i − 2*fv(x) ]
```

(continuum D_iD_iΦ: ∇²u − g(∇·A)v − 2gA·∇v − g²A²u + i-swap, GAUGE_DESIGN G7a,b —
the seagull is inside the cos.)

### 3.3 Covariant curl (the η term)

Covariant **central** difference (transforms covariantly at x; reduces at θ=0 to the
existing `curl_component` stencil exactly):

```
DcU_i(F)(x) = (TRp_i − TRm_i) / (2a)        /* real part of D_i^c F */
DcV_i(F)(x) = (TIp_i − TIm_i) / (2a)        /* imag part */
```

Covariant curl, component a (same index pattern as curl_component, scp_sim.c:136–142):

```
(D×F)_0 = D_1^c F_2 − D_2^c F_1     i.e.  re: DcU_1(F_2) − DcU_2(F_1),  im: DcV_1(F_2) − DcV_2(F_1)
(D×F)_1 = D_2^c F_0 − D_0^c F_2
(D×F)_2 = D_0^c F_1 − D_1^c F_0
```

**Matter forces** (per voxel, a = 0..2; Vp and prod_rest exactly as
compute_forces_complex, scp_sim.c:592–602):

```
phi_acc[a]      (u̇̇):  LapU(Φ_a) − m²u_a − 2·Vp·u_a·prod_rest[a] + η·Re(D×Θ)_a
phi_im_acc[a]   (v̈):  LapV(Φ_a) − m²v_a − 2·Vp·v_a·prod_rest[a] + η·Im(D×Θ)_a
theta_acc[a]    (tü): LapU(Θ_a) − m_θ²tu_a + η·Re(D×Φ)_a
theta_im_acc[a] (tv̈): LapV(Θ_a) − m_θ²tv_a + η·Im(D×Φ)_a
```

Note the complex pairing: the gauged theory couples (D×Θ) into Φ̈ as a COMPLEX
equation (GAUGE_DESIGN §2 EOM) — at g=0 it reduces to the v66 (u,tu)/(v,tv) pairing
because Re(∇×Θ)_a = curl(tu)_a and Im(∇×Θ)_a = curl(tv)_a. These force expressions
are the exact −∂H/∂field of §3 H (the Φ-side η force was derived by lattice
summation-by-parts of the −η Re[Θ̄(D×Φ)] term; it comes out as +η(D×Θ) with the SAME
covariant central difference — this exact adjointness is required, do not
"simplify" one side to a forward difference).

### 3.4 E-kick: staples + lattice current (the core derivation)

Define plaquette-plane scratch (3 planes p: p0=(i,j)=(0,1), p1=(1,2), p2=(2,0)):

```
plaq_s[p][x] = sin( th[i][x] + th[j][x+î] − th[i][x+ĵ] − th[j][x] )    /* (i,j) per plane p */
```

**E-kick** stored in `E_acc[i][x]` (computed in the same force pass, after the
scratch passes):

```
E_acc[i][x] = (1/(g a³)) · Σ_{j≠i} [ S_ij(x) − S_ij(x−ĵ) ]   +   g · J_i^lat(x)
```

with S_ij = ±plaq_s by antisymmetry (S_ij = −S_ji); explicitly:

| i | staple sum |
|---|---|
| 0 (x) | +[plaq_s0(x) − plaq_s0(x−ĵ)] − [plaq_s2(x) − plaq_s2(x−k̂)] |
| 1 (y) | −[plaq_s0(x) − plaq_s0(x−î)] + [plaq_s1(x) − plaq_s1(x−k̂)] |
| 2 (z) | +[plaq_s2(x) − plaq_s2(x−î)] − [plaq_s1(x) − plaq_s1(x−ĵ)] |

(continuum limit +∂_j F_ij = +[∇×B]_i; the discrete divergence of the staple field
vanishes IDENTICALLY by S_ij = −S_ji telescoping — exact, no rounding caveat).

**The lattice current** — this is THE current, fixed uniquely as
J_i^lat ≡ (1/(g))·(g a/a³)·∂H_matter/∂θ_i(x), i.e. the lattice Hamiltonian's
A-derivative; any other choice breaks exact Gauss conservation:

```
J_i^lat(x) = (1/a) Σ_a Im[ Φ̄_a(x) U_i(x) Φ_a(x+î) ]                         (Φ gradient+gauge+seagull, all-in-one)
           + (1/a) Σ_a Im[ Θ̄_a(x) U_i(x) Θ_a(x+î) ]                         (Θ ditto)
           − (η/2) Σ_{j,k} ε_ijk { Im[ Θ̄_j(x)   U_i(x)  Φ_k(x+î) ]
                                  + Im[ Θ̄_j(x+î) U_i(x)† Φ_k(x)   ] }        (η seagull, symmetrized over the link ends)
```

Component forms (W ≡ U_i(x)·F(x+î) = (TRp,TIp) from §3.1; W′ ≡ U_i(x)†·F(x) with
components `(c⁺*fu + s⁺*fv, c⁺*fv − s⁺*fu)` — note: link angle AT x, field AT x):

```
Im[ Φ̄_a(x) U_i Φ_a(x+î) ]  =  u_a(x)*TIp_a − v_a(x)*TRp_a          /* at theta=0: u v⁺ − v u⁺, the standard Noether link current */
Im[ Θ̄_j(x) U_i Φ_k(x+î) ]  =  tu_j(x)*TIp_k − tv_j(x)*TRp_k
Im[ Θ̄_j(x+î) U_i† Φ_k(x) ]  =  tu_j(x+î)*W′I_k − tv_j(x+î)*W′R_k
```

ε_ijk written out — the η piece of J_i^lat for i=0 is
−(η/2)[ T(1,2) − T(2,1) ] with T(j,k) ≡ the {…} braces above; cyclic for i=1,2.

Continuum check (machine-checkable): J_i^lat → Im[Φ̄D_iΦ] + Im[Θ̄D_iΘ] − η(u×tv + tu×v)_i
= **−J_i** of GAUGE_DESIGN §2 (the Noether current). The sign is consistent:
Ė_i = ∂_jF_ij − gJ_i = staples + g·J_i^lat. Implementers: J^lat is minus the design's
J; the E-kick formula above is the ground truth.

Division-by-zero guard: the staple prefactor 1/(g a³) and the Q_flux/E_em magnetic
prefactor 1/(g²a⁴) are never evaluated when g==0 (g==0 short-circuits to the
ungauged path, §1.3).

### 3.5 Integration order (verlet_step, scp_sim.c:752–790)

The gauge sector is integrated EXACTLY like a matter pair (θ ↔ position,
E ↔ −velocity/(ga), E_acc ↔ force), interleaved at the same three stages:

```
1. half-kick:   vel += hdt*acc (all 12 matter, as today)
                AND  Efield[i] += hdt*E_acc[i]                      (gauge_mode && g!=0)
2. drift:       fields += dt*vel (all 12, as today)
                AND  th[i] += dt*(−g_gauge*a*Efield[i]); wrap th[i] to (−π,π]
                     /* wrap: th −= 2π*rint(th/(2π)) */
3. forces:      if (gauge_mode && g!=0) compute_forces_complex_gauge(g,c)
                else if (complex_mode)  compute_forces_complex(g,c)
                else                    compute_forces(g,c)
                /* _gauge fills all 12 matter acc arrays AND E_acc[3]:
                   pass A: link_c/link_s = cos/sin(th) (3 N³ sincos);
                   pass B: plaq_s[3];
                   pass C: fused per-voxel loop -> matter accs (§3.2/3.3) + E_acc (§3.4) */
4. half-kick:   same as 1 (matter + E)
5. BCs:         bc_type==0: apply_damping — damps the 12 matter velocity arrays
                (as today) AND Efield[0..2] by the same factor d (§3.6).
                bc_type==2: nothing. (bc_type==1 refused, §1.1)
```

`main()` already computes forces once after init (scp_sim.c:1283) — dispatch there
identically. **Why this conserves the discrete Gauss law to rounding** (no projection
or damping step in the loop): (kicks) E and matter velocities change simultaneously
using forces/currents from the SAME H at the SAME configuration; exact gauge
invariance of H gives the pointwise lattice Noether identity
g·Σ Im[Φ̄F_Φ + Θ̄F_Θ](x) = (g/a)·Σ_i[J_i^lat(x) − J_i^lat(x−î)] plus the identically-
divergence-free staples ⟹ ΔG(x) = 0 per half-kick. (drifts) ρ_Q = Im[Φ̄Φ̇+Θ̄Θ̇] is
invariant (ΔΦ ∝ Φ̇ ⟹ ΔIm[Φ̄Φ̇] = dt·Im[Φ̇̄Φ̇] = 0) and E, hence div E, untouched
(θ-drift doesn't enter G). So G(x) is constant at every sub-stage [standard
real-time KS result; continuum analog G9+G10].

### 3.6 Boundary conditions (decision)

bc_type=0 damps **matter velocities and E_i** with the same sponge factor (kills
outgoing A radiation, which would otherwise wrap — the stencil is periodic).
Consequences (documented, accepted for v1): (i) the static Coulomb tail inside the
sponge shell erodes toward 0; (ii) Gauss violations appear ONLY in the shell and are
frozen there (constraint violations do not propagate in temporal gauge); (iii) all
gauge diagnostics are therefore interior-only (§4) and `qdiag_radius`/flux radii
must satisfy R < L − damp_width. θ links are positions and are NOT damped.

---

## 4. GAUSS / GAUGE DIAGNOSTICS

### 4.1 Definitions

```
G(x)      = (1/a)·Σ_i [ Efield[i][x] − Efield[i][x−î] ]  −  g·ρ_Q(x)
ρ_Q(x)    = Σ_a ( u v̇ − v u̇ + tu tv̇ − tv tu̇ )            (velocities are full-step synchronized at diag time)
```

**Torus obstruction + frozen offset**: on the wrapped lattice Σ_x div E ≡ 0
(telescoping), so Σ_x G = −g·Q_total: a net-charged absorbing-BC run CANNOT have
G ≡ 0 — the best achievable is G(x) = uniform jellium offset
Ḡ = −g·Q_total(0)/(N³a³), which the integrator then freezes pointwise for all t.
The kernel measures `G_offset = (1/N³)Σ_x G(x)` ONCE, after init + projection
(§5.4), stores it in Grid, and all residual diagnostics are relative to it.
(bc_type=2 is net-neutral by §1.2 ⟹ G_offset ≈ 0 there.)

Interior domain Ω: voxels with r < L − damp_width when bc_type=0 (sponge excluded,
§3.6); the whole box when bc_type=2.

```
gauss_max = max_{x∈Ω} | G(x) − G_offset |
gauss_l2  = sqrt( (1/|Ω|) Σ_{x∈Ω} (G(x) − G_offset)² )
E_em      = Σ_x a³ [ ½ Σ_i Efield[i][x]²  +  (g!=0 ? (1/(g²a⁴)) Σ_{p=0..2} (1 − cos θ_p(x)) : 0) ]
Q_flux    = (g!=0) ? (1/g)·[ Σ_{x∈C} a³·(1/a)Σ_i(E_i(x)−E_i(x−î))  −  G_offset·|C|·a³ ] : 0
```

Q_flux: C is the centered **cube** of half-width `qdiag_radius` (default 8.0) —
a cube, not a sphere, because the volume-sum of the discrete divergence telescopes
EXACTLY to the surface link flux Σ_faces E_n·a², making the readout an E-only
boundary measurement with no interpolation error (GAUGE_DESIGN P3 "charge readable
from boundary flux"; the sphere wording is implemented as its exact lattice
analog). The −G_offset·|C|a³ term removes the known jellium bias exactly (without
it, a single ball in a 30³ box with a 20³ cube would read ~30% low). Expected:
Q_flux ≈ Q_enclosed = Q_core to few % (test §7.3).

### 4.2 diag.tsv columns

`complex_gauge=0`: byte-identical to today. `complex_gauge=1`: **append 4 columns
after the v67 column `thp2v`** (existing 25 columns keep position and meaning):

```
... thp1u thp1v thp2u thp2v  gauss_max  gauss_l2  E_em  Q_flux
```

Formats `%.6e` (gauss_*), `%.12e` (E_em, Q_flux). Implemented in a new
`compute_gauss(Grid*, const Config*, double *gauss_max, double *gauss_l2, double *e_em, double *q_flux)`
placed after `compute_charges` (scp_sim.c:~1082), called at the three diag sites
(initial 1413–1419, in-loop 1555–1561, final 1598–1611) under `if (c.complex_gauge)`.
Console major line gains `gauss=%.1e E_em=%.3e` when gauged.

### 4.3 E_total semantics (decision)

When `complex_gauge=1`, the `E_total` column = the 8 matter energies **+ E_em**
(it must be the conserved quantity for test §7.4; E_em is also printed separately
so the matter-only sum is recoverable). At g=0, E_em ≡ 0 and E_total is unchanged.
Energy bookkeeping in the gauged energy function (`compute_energy_complex_gauge`,
mirror of compute_energy_complex placed after it, scp_sim.c:~954):

- E_grad / E_tgrad: covariant **central** differences (§3.3 DcU/DcV, ½Σ|D_i^cF|²) —
  at g=0 these reduce bit-arithmetically to the v66 central-difference columns.
  (Known, accepted mismatch vs the forward-difference Hamiltonian §3 — inherited
  from the real/v66 kernel, diagnostic-only, shows up as the usual Verlet-floor
  drift, gated in §7.4.)
- E_coupling: −η Σ_a ∫ (u_a·Re(D×Θ)_a + v_a·Im(D×Θ)_a) dV — the v66 pairing with
  transported curls; g=0-identical to v66.
- E_mass, E_pot, E_kin columns: unchanged formulas (gauge-invariant point terms).
- E_em per §4.1, returned separately and added into E_total.

---

## 5. INIT

### 5.1 Shooter profile contract (coordination with the v69 shooter agent)

Files `v69/theory/profile_g<G>_omega<W>.txt` (e.g. `profile_g0.0500_omega1.3900.txt`),
text, `#` comments. Header comments MUST include `omega=`, `g=`, `m2=`, `mu=`,
`kappa=`. Data lines: **three** whitespace-separated columns

```
r   f   Er
```

r strictly increasing from 0; f the gauged radial profile (solves
f″+(2/r)f′ = (m²−ω̃²)f + 2Vt′(f⁶)f⁵, ω̃(r) = ω + g·a₀(r)); **Er(r) ≥ 0 for Q > 0**,
satisfying the radial Gauss law Er(r) = (g/r²)∫₀^r 3ω̃(s)f(s)²s²ds and
Er = −da₀/dr (a₀ > 0, Coulomb-like, a₀(∞)=0). Sign repair note: O1.
A legacy v66 2-column (r f) profile is also accepted (§5.3 fallback).

### 5.2 init=qball, gauged (extends scp_init.h:258–338)

Matter fill, single ball at (x0,y0,z0) (sentinel 1e30 → 0 as today):

```
ω̃(r)  = qball_omega + g·a0(r)
a0(r)  = ∫_r^{R_last} Er(s) ds  +  Er(R_last)·R_last        /* trapezoid on the table + exact 1/r² tail:
                                                               ∫_R^∞ gQ_∞/(4πs²)ds = Er(R)·R */
u_a = f(r);  v_a = 0;  udot_a = 0;  vdot_a = ω̃(r)·f(r)      /* a = 0,1,2; theta sector zero */
```

(Precompute the a0 table once from the (r,Er) table — cumulative backward
trapezoid — then interpolate a0 linearly alongside f. ρ_Q(0) = 3ω̃f² then matches
the shooter's Gauss source self-consistently. With a 2-column profile: ω̃ ≡ ω.)

**E seeding on links**: for each link (x,i), evaluate at the link **midpoint**
m = x + (a/2)î relative to the ball center, r_m = |m|, unit radial r̂ = m/r_m:

```
Efield[i][x] += Er_interp(r_m) · r̂_i
```

Er_interp: linear in the table; r < r₀: Er(r₀)·(r/r₀) (E→0 linearly at center);
r > R_last: **Er(R_last)·(R_last/r)²** (inverse-square continuation — do NOT zero;
the tail carries the flux). th[i] ≡ 0 (Coulomb is curl-free, G12c).

**Second ball** (new optional config keys, enables the §7.5 pair tests):
`qball2_x0/y0/z0` (sentinel 1e30 = disabled — all three must be set),
`qball2_sign` (int, +1 default; −1 = conjugate/opposite-charge ball),
`qball2_phase` (double, radians, default 0). With σ = qball2_sign, δ = qball2_phase,
f₂/ω̃₂ from the SAME profile evaluated at r₂ = |x − x₀₂|:

```
u_a    += f₂·cos(δ);            v_a    += f₂·sin(δ)
udot_a += −σ·ω̃₂·f₂·sin(δ);     vdot_a += +σ·ω̃₂·f₂·cos(δ)
Efield[i][x] += σ · Er_interp(r₂_mid) · r̂₂_i                 /* Gauss is linear: superpose */
```

(σ=−1 ⟹ Φ₂ ~ f e^{−iω̃t}: negative ρ_Q, inward E. Print both centers/signs in the
init line.) Tail overlap at D=16 is e^{−0.564·8} ≈ 1e-2 — superposition error is
absorbed by the projection pass and the O(g²) settling transient.

### 5.3 Other init modes under complex_gauge=1

- `init=oscillon` / `template`: real-sector seeding as today; ρ_Q(0)=0, E=0, links
  0 — Gauss already exact; projection pass is a no-op. Allowed.
- `init=sfa` (30 columns, §6): loads th/E per §6 slot rules; a 24-column v66 file
  (e.g. `gen_qball_pair` output) loads matter only, E=0 — the projection pass §5.4
  then **constructs the full periodic-Coulomb E from scratch** (this is the general
  fallback for any matter-only seed; for net-charge it yields the jellium field).
  Warn: `init=sfa: no gauge columns; E seeded by Gauss projection`.
- `init=exec`: via the sfa path, same rules.
- 2-column qball profile with g≠0: WARN
  `ungauged profile (no Er column): seeding vdot=omega*f, E from projection; O(g^2) settling transient expected` — then E=0 seed + projection.

### 5.4 Gauss projection pass (one-time, init only) — MANDATORY

Runs in `main()` after `do_init()` (and after the §1.2 net-charge check) whenever
`complex_gauge && g_gauge != 0`, for EVERY init mode (idempotent; also repairs f32
quantization on snapshot restarts — note: bit-exact restarts need `precision=f64`):

1. Compute G(x) (§4.1) from the seeded state; Ḡ = mean over the full box.
2. Solve the 7-point periodic lattice Poisson `∇²_lat χ = −(G − Ḡ)` by **conjugate
   gradient** (zero-mean RHS ⟹ consistent on the torus; keep χ zero-mean per
   iteration; no FFTW dependency — CG on N=128 converges in a few hundred stencil
   passes, seconds). Stop at `max|residual| < 1e-13·max(1, max|G−Ḡ|)` or 2·10⁴ iters
   (hard-error on non-convergence).
3. Correct the links' E by the forward gradient: `Efield[i][x] += (χ(x+î) − χ(x))/a`
   ⟹ div-correction = ∇²_lat χ exactly ⟹ G(x) → Ḡ uniformly, to the CG tolerance.
4. Measure and store `g->G_offset = mean_x G(x)` (≈ −g·Q_total/(N³a³); ≈0 net-neutral),
   print `Gauss projection: %d CG iters, gauss_max(0)=%.2e (offset %.3e)`.

Result: `gauss_max(0) ≲ 1e-12` (test §7.2 gate). The radial seed (§5.2) makes the
correction small (interpolation-level), keeping the field near the open-boundary
Coulomb solution rather than the periodic-image one.

---

## 6. OUTPUT (SFA)

`complex_gauge=1`: **30 columns** = the 24 v66 columns (names/order unchanged,
scp_sim.c:1299–1325) + 6 appended:

| col | name | semantic | comp | array |
|---|---|---|---|---|
| 24–26 | `th_x` `th_y` `th_z` | SFA_ANGLE | 6,7,8 | th[0..2] |
| 27–29 | `E_x` `E_y` `E_z` | SFA_VELOCITY | 12,13,14 | Efield[0..2] |

(Component offsets continue the v66 convention — ANGLE 0–2 real Θ, 3–5 imag Θ,
6–8 links; VELOCITY 0–11 matter, 12–14 E. Names ≤ 4 chars, fine. E_acc is NOT
written — recomputed at restart like matter accs, so 30 columns is an exact-restart
set at f64.) `sfa_snap` arrays[30] in registration order; `nf = gauge ? 30 :
complex ? 24 : 12`. **init=sfa** slot rules (scp_init.h:96–112, inside the
SCP_COMPLEX_FIELDS guard, active iff `complex_gauge`): SFA_ANGLE comp 6–8 → th[comp−6]
(slot 24–26); SFA_VELOCITY comp 12–14 → Efield[comp−12] (slot 27–29); when
`complex_gauge=0` meets these columns: skip with one-line WARNING (mirror of the
imaginary-sector warning). KVMD (scp_init.h:357–400): append keys
`complex_gauge`, `g_gauge` when complex_gauge=1 (29 keys) so `.sfa`-restart
reallocates the gauge sector.

---

## 7. TEST PLAN

Configs `v69/cfg/`, small outputs `v69/results/`, SFA > 100 MB → `/space/scp/v69/`.
Build: `gcc -O3 -march=native -fopenmp -o bin/scp_sim sfa/sim/scp_sim.c -lzstd -lm`
plus `make -C sfa install` regression. Ball: shooter profile g=0.05, ω=1.39
(Q≈482, E≈692, r_half≈4.4; exact numbers from the gauged shooter scan).

1. **(a) g=0 identity.** (a1) `complex_gauge=1, g_gauge=0` vs `complex_phi=1` on the
   v66 reg_cplx config (oscillon, N=64, L=15, T=20): diag energy/charge columns
   **bit-identical** (same code path; §1.3); gauge columns all exactly 0.
   (a2) same but `g_gauge=1e-8`, init=qball with the v66 2-col profile: every energy
   column within **1e-9 relative** of the complex_phi run at all t ≤ 20 (covariant-
   stencil reduction; fallback ladder as v66 §9b).
2. **(b) Gauss flatness.** Single gauged ball, g=0.05, η=0, bc_type=0, N=96, L=15,
   T=20, dt_factor=0.025, diag_dt=0.5. PASS: gauss_max(0) ≤ 1e-12;
   gauss_max(t) ≤ 1e-10 for all t with no secular trend (final ≤ 5× initial;
   growth is rounding accumulation only). FAIL ⟹ a non-gauge-invariant stencil
   (bug class #1, GAUGE_DESIGN P2 falsifier). Supporting unit check: config flag
   `test_gauge_xform=1` (new, default 0) — after init+forces, apply a deterministic
   random lattice gauge transformation (α(x) from a fixed-seed LCG), recompute
   forces and E_acc, verify all gauge-invariant bilinears (|acc| of each complex
   pair, E_acc, G) match to ≤ 1e-12 relative, print PASS/FAIL, exit.
3. **(c) Single gauged ball physics.** As (b) with T=100, qdiag_radius=10,
   damp_width=3 (R_damp=12 > 10 ✓). PASS: s_max stddev/mean < 5% over t∈[20,100]
   (stationary, non-radiating); Q_phi retention > 99%; **|Q_flux − Q_core|/|Q_core|
   ≤ 3%** (P3: charge read from E flux alone); E_em(t) settles to the static Coulomb
   value with < 1% drift over t∈[20,100]; ω_core ≈ ω + δω with δω > 0 (P4 sign check:
   wrong sign ⟹ charge-convention error).
4. **(d) Energy conservation incl. E_em.** Run (c)'s config at bc_type=2 with a
   ±pair (net-neutral; or the same single ball with bc_type=0 measuring interior
   only), T=50: |E_total(t) − E_total(0)|/E_total(0) < 1e-3, non-secular, and within
   2× the drift of the matched g=0 control run (Verlet floor inherited from v66).
5. **(e) Two-ball force sign.** N=160, L=20 (dx≈0.251), bc_type=0, damp_width=3,
   T=100, snap_dt=5, f32 snaps. Balls at x=±8 (D=16 > r*≈12), qball2 keys.
   (e1) same charge, co-phase (sign=+1, phase=0): centroid separation (offline:
   |Φ|²-weighted centroid in each half-box x≷0 from snapshots — no kernel change)
   INCREASES, ΔD ≥ +0.5 over T=100, monotonic trend after t≈20
   (Coulomb estimate: F = g²Q²/4πD² ≈ 0.18, a = F/E₀ ≈ 2.6e-4, ΔD ≈ 2·½at² ≈ 2.6 —
   well above the gate). (e2) opposite charge (qball2_sign=−1): ΔD ≤ −0.5
   (attraction; net-neutral ⟹ G_offset≈0, cleanest run). Both runs: gauss_max stays
   at the (b) floor.
6. **GPU port: separate later phase.** Only after (a)–(e) pass on CPU is a
   `scp_sim.cu` port authorized for request; this spec deliberately keeps all gauge
   code out of the shared headers' GPU paths (SCP_COMPLEX_FIELDS guard already
   isolates scp_init.h; the new init code goes inside that guard).

---

## 8. FILE-BY-FILE CHANGE LIST

### `sfa/sim/scp_config.h`

| location | edit |
|---|---|
| Config struct, after `qdiag_probe2` (line 79) | add `int complex_gauge; double g_gauge; double qball2_x0, qball2_y0, qball2_z0; int qball2_sign; double qball2_phase; int test_gauge_xform;` |
| `cfg_defaults` (after line 130) | `complex_gauge=0; g_gauge=0.05; qball2_x0/y0/z0=1e30; qball2_sign=1; qball2_phase=0; test_gauge_xform=0;` |
| `cfg_set` (after line 196) | keys `complex_gauge, g_gauge, qball2_x0, qball2_y0, qball2_z0, qball2_sign, qball2_phase, test_gauge_xform` |
| `cfg_print` complex block (272–282) | when complex_gauge: `printf("Gauge:   U(1) compact links, g=%.4f (complex_gauge=1)\n", c->g_gauge);` + second-ball line when qball2 enabled |
| `cfg_validate` (305–344) | §1.1 table: complex_gauge requires complex_phi; refuse bc_type==1, bc_switch_time!=0, g_gauge<0 under complex_gauge |

### `sfa/sim/scp_sim.c`

| function (lines) | edit |
|---|---|
| Grid struct (32–48) | §2.2 members (gauge_mode, g_gauge, G_offset, th/Efield/E_acc/link_c/link_s/plaq_s) |
| `grid_alloc` (50–95) | 54-block layout (§2.2); gauge pointers NULL unless complex_gauge; copy gauge_mode/g_gauge |
| `compute_forces` / `compute_forces_complex` (387–638) | **UNTOUCHED** |
| new `compute_forces_complex_gauge` (after 638) | §3: pass A cos/sin, pass B plaq_s, pass C fused matter accs + E_acc; OpenMP per pass |
| `apply_damping` (645–664) | inside `if (g->complex_mode)` block: also `if (g->gauge_mode) for(i) Efield[i][idx]*=d;` |
| `verlet_step` (752–790) | E half-kicks beside the matter im-sector kicks; θ drift+wrap beside the drifts; 3-way force dispatch (g==0 ⟹ ungauged complex path, §1.3) |
| `compute_energy_complex` (897–954) | UNTOUCHED; new `compute_energy_complex_gauge` after it (§4.3, extra out-param e_em) |
| new `compute_gauss` (after compute_charges, ~1082) | §4.1: G, gauss_max/l2 over Ω, E_em, Q_flux (cube, jellium-corrected) |
| new `init_gauss_project` (near compute_gauss) | §5.4 CG Poisson + E correction + G_offset; plus the §1.2 net-charge check helper |
| new `do_test_gauge_xform` | §7.2 unit check, runs+exits when test_gauge_xform=1 |
| `sfa_snap` (1172–1204) | arrays[30], `nf = gauge?30:complex?24:12` |
| `main` | force dispatch (1283–1284); after do_init: net-charge check + projection + optional xform test; SFA columns th_x..E_z after line 1325 (§6); diag header conditional after the v67 columns (1385–1386); the three diag row sites (1413–1419, 1555–1561, 1598–1611) call compute_gauss and append 4 columns; energy dispatch to *_gauge variant; console additions |

### `sfa/sim/scp_init.h`

| location | edit |
|---|---|
| `init_from_sfa` (96–112) | th/E slot rules under SCP_COMPLEX_FIELDS && complex_gauge (§6); skip+WARN otherwise |
| `init_qball` (258–338) | 3-column (r f Er) parsing (2-col fallback + WARN per §5.3); a0 cumulative integral + tail; ω̃ seed velocities; link-midpoint radial E seeding with 1/r² continuation; optional second ball (§5.2) |
| `sfa_embed_kvmd` (357–400) | append `complex_gauge`, `g_gauge` keys when complex_gauge=1 |

### New files

`v69/cfg/*.cfg` (tests §7), `v69/theory/profile_*` (shooter agent, contract §5.1).
Generated, not committed: `v69/results/*_diag.tsv`.

### Explicitly NOT modified

`sfa/format/sfa.h`, `sfa/sim/scp_sim.cu`, `sfa/seed/*` (gen_qball_pair works as-is
via the §5.3 init=sfa + projection fallback), all analysis tools.

---

## OPEN QUESTIONS / GAUGE_DESIGN ambiguities (resolved here, flagged for the record)

- **O1 (sign typo in GAUGE_DESIGN §2)**: it states "E_i = ∂_iA_0 − ∂_tA_i", but
  consistency with ∇·E = +gρ_Q (G12), with a₀″+(2/r)a₀′ = −3gω̃f² (§6, a₀ > 0
  Coulomb-like), and with E_r > 0 for Q > 0 requires **E = −∇A_0 − ∂_tA**. This spec
  adopts the standard sign; in temporal gauge E = −Ȧ is unaffected, so no equation
  in this spec changes — but the shooter must emit Er = −da₀/dr ≥ 0 (§5.1).
- **O2 (torus jellium for net-charge absorbing runs)**: GAUGE_DESIGN treats the
  periodic-BC neutrality theorem (§4) but not the fact that bc_type=0 still wraps
  the stencil ⟹ a single charged ball ALWAYS carries a uniform frozen Gauss offset
  Ḡ = −gQ/V_box. Resolved by G_offset bookkeeping + jellium-corrected Q_flux (§4.1);
  physically it is an exact uniform neutralizing background. Bigger boxes shrink the
  jellium field ∝ 1/V; revisit if precision Coulomb-halo fits (P3) need open-boundary
  solvers.
- **O3 (sponge vs gauge sector)**: design risk #2 suggests damping "the TRANSVERSE
  gauge sector only" — a Helmholtz split is overkill for v1; we damp all E in the
  sponge and accept shell-confined Gauss violation + shell Coulomb-tail erosion
  (§3.6). If two-ball runs show spurious boundary forces, implement the
  transverse-only damping as v2.
- **O4 (`gauge_compact` config key, design §6)**: v1 implements compact links only;
  the key is omitted (noncompact fallback not built). Flux quantum 2π/(ga) ≈ 550 at
  g=0.05, dx≈0.23 — winding sectors unreachable except by deliberate seeding; the
  AB/flux-line experiments (P10/V5) need a flux-tube seeder that is NOT specced here.
- **O5 (full package masses)**: m_θ = 1.6 is the recommended package but `mtheta2`
  default stays 0 (locked project default); P7 runs must set `m_theta=1.6`
  explicitly. Any gauged m_θ=0, η≠0 run is drain-afflicted (design risk #5) and must
  be labeled as such.
- **O6 (positronium orbit P9 / watershed P8 inner branch)**: D < r* merger runs and
  orbit-closure measurements need finer two-ball diagnostics (per-ball Q, momentum)
  than the half-box centroid analysis specced in §7.5 — deferred to the experiment
  plan, not a kernel gap.
