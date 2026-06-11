/*  scp_init.h — Shared initialization for CPU and CUDA simulation kernels
 *
 *  Requires: scp_config.h included first, SFA_IMPLEMENTATION defined,
 *  and a Grid struct with at minimum:
 *    double *phi[3], *phi_vel[3], *theta[3], *theta_vel[3]
 *    double *pin_phi[3], *pin_vel[3], *pin_theta[3], *pin_tvel[3]
 *    int N; long N3; double L, dx, dt;
 */

#ifndef SCP_INIT_H
#define SCP_INIT_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

static void init_oscillon(Grid *g, const Config *c) {
    const int N = g->N, NN = N * N;
    const double L = g->L, dx = g->dx;
    printf("Init: oscillon (A=%.3f sigma=%.3f)\n", c->A, c->sigma);
    for (int i = 0; i < N; i++) { double x = -L + i*dx;
    for (int j = 0; j < N; j++) { double y = -L + j*dx;
    for (int k = 0; k < N; k++) { double z = -L + k*dx;
        long idx = (long)i*NN + j*N + k;
        double r2 = x*x + y*y + z*z;
        double env = c->A * exp(-r2 / (2.0 * c->sigma * c->sigma));
        for (int a = 0; a < NFIELDS; a++)
            g->phi[a][idx] = env * cos(c->delta[a]);
    }}}
}

static void init_braid(Grid *g, const Config *c) {
    const int N = g->N, NN = N * N;
    const double L = g->L, dx = g->dx;
    const double kw = PI/L, omega = sqrt(kw*kw + c->m2);
    const double sx = 1+c->ellip, sy = 1-c->ellip;
    const double inv2R2 = 1.0/(2*c->R_tube*c->R_tube);
    const double k_bg = PI/L, omega_bg = sqrt(k_bg*k_bg + c->m2);
    printf("Init: braid (R=%.1f ellip=%.4f A=%.2f A_bg=%.2f)\n",
           c->R_tube, c->ellip, c->A, c->A_bg);
    for (int i = 0; i < N; i++) { double x = -L + i*dx;
    for (int j = 0; j < N; j++) { double y = -L + j*dx;
    for (int k = 0; k < N; k++) { double z = -L + k*dx;
        long idx = (long)i*NN + j*N + k;
        double r2e = x*x/(sx*sx) + y*y/(sy*sy);
        double env = exp(-r2e * inv2R2);
        for (int a = 0; a < NFIELDS; a++) {
            double ph = kw*z + c->delta[a];
            double ph_bg = k_bg*z + 2*PI*a/3.0;
            g->phi[a][idx] = c->A*env*cos(ph) + c->A_bg*cos(ph_bg);
            g->phi_vel[a][idx] = omega*c->A*env*sin(ph) + omega_bg*c->A_bg*sin(ph_bg);
        }
    }}}
}

static void init_from_sfa(Grid *g, const Config *c) {
    printf("Init: SFA file '%s' frame=%d\n", c->init_sfa, c->init_frame);
    SFA *sfa = sfa_open(c->init_sfa);
    if (!sfa) { fprintf(stderr, "FATAL: cannot open SFA '%s'\n", c->init_sfa); exit(1); }
    if ((int)sfa->Nx != g->N || (int)sfa->Ny != g->N || (int)sfa->Nz != g->N) {
        fprintf(stderr, "FATAL: SFA grid %ux%ux%u != sim grid %d^3\n",
                sfa->Nx, sfa->Ny, sfa->Nz, g->N);
        exit(1);
    }
    int frame = c->init_frame;
    if (frame < 0) frame = sfa->total_frames + frame;
    if (frame < 0 || frame >= (int)sfa->total_frames) {
        fprintf(stderr, "FATAL: frame %d out of range [0,%u)\n", frame, sfa->total_frames);
        exit(1);
    }
    printf("  Grid: %ux%ux%u, %d columns, %u frames, reading frame %d\n",
           sfa->Nx, sfa->Ny, sfa->Nz, sfa->n_columns, sfa->total_frames, frame);

    void *buf = malloc(sfa->frame_bytes);
    if (!buf) { fprintf(stderr, "FATAL: frame buffer alloc\n"); exit(1); }
    sfa_read_frame(sfa, frame, buf);

    int loaded[30] = {0};
    int warned_im = 0, warned_gauge = 0;
    uint64_t off = 0;
    for (uint32_t col = 0; col < sfa->n_columns; col++) {
        int dtype = sfa->columns[col].dtype;
        int sem   = sfa->columns[col].semantic;
        int comp  = sfa->columns[col].component;
        int es    = sfa_dtype_size[dtype];
        uint8_t *src = (uint8_t*)buf + off;

        double *target = NULL;
        int slot = -1;
        if (sem == SFA_POSITION && comp < 3) { target = g->phi[comp]; slot = comp; }
        else if (sem == SFA_ANGLE && comp < 3) { target = g->theta[comp]; slot = 3+comp; }
        else if (sem == SFA_VELOCITY && comp < 3) { target = g->phi_vel[comp]; slot = 6+comp; }
        else if (sem == SFA_VELOCITY && comp >= 3 && comp < 6) { target = g->theta_vel[comp-3]; slot = 9+comp-3; }
        else if ((sem == SFA_POSITION && comp >= 3 && comp < 6) ||
                 (sem == SFA_ANGLE    && comp >= 3 && comp < 6) ||
                 (sem == SFA_VELOCITY && comp >= 6 && comp < 12)) {
            /* v66 imaginary-sector columns (SPEC §5.2/§7.1) */
#ifdef SCP_COMPLEX_FIELDS
            if (c->complex_phi) {
                if      (sem == SFA_POSITION)             { target = g->phi_im[comp-3];       slot = 12+comp-3; }
                else if (sem == SFA_ANGLE)                { target = g->theta_im[comp-3];     slot = 15+comp-3; }
                else if (sem == SFA_VELOCITY && comp < 9) { target = g->phi_im_vel[comp-6];   slot = 18+comp-6; }
                else                                      { target = g->theta_im_vel[comp-9]; slot = 21+comp-9; }
            } else
#endif
            if (!warned_im) {
                printf("  WARNING: file has imaginary-sector columns; skipped (complex_phi=0)\n");
                warned_im = 1;
            }
        }
        else if ((sem == SFA_ANGLE    && comp >= 6  && comp < 9) ||
                 (sem == SFA_VELOCITY && comp >= 12 && comp < 15)) {
            /* v69 gauge-sector columns (SPEC §6): th links + E */
#ifdef SCP_COMPLEX_FIELDS
            if (c->complex_gauge) {
                if (sem == SFA_ANGLE) { target = g->th[comp-6];      slot = 24+comp-6; }
                else                  { target = g->Efield[comp-12]; slot = 27+comp-12; }
            } else
#endif
            if (!warned_gauge) {
                printf("  WARNING: file has gauge-sector columns; skipped (complex_gauge=0)\n");
                warned_gauge = 1;
            }
        }

        if (target) {
            long N3 = g->N3;
            if (dtype == SFA_F64)
                for (long i = 0; i < N3; i++) target[i] = ((double*)src)[i];
            else if (dtype == SFA_F32)
                for (long i = 0; i < N3; i++) target[i] = (double)((float*)src)[i];
            else if (dtype == SFA_F16)
                for (long i = 0; i < N3; i++) target[i] = f16_to_f64(((uint16_t*)src)[i]);
            loaded[slot] = 1;
            printf("  Loaded col %d '%s' (sem=%d comp=%d dtype=%d) -> slot %d\n",
                   col, sfa->columns[col].name, sem, comp, dtype, slot);
        }
        off += (uint64_t)g->N3 * es;
    }

    int n_fields = 0, n_vels = 0, n_gauge = 0;
    for (int i = 0; i < 6; i++) n_fields += loaded[i];
    for (int i = 6; i < 12; i++) n_vels += loaded[i];
    for (int i = 24; i < 30; i++) n_gauge += loaded[i];
    printf("  Loaded: %d/6 field arrays, %d/6 velocity arrays\n", n_fields, n_vels);
    if (n_vels == 0)
        printf("  WARNING: no velocity data — starting from rest (cold restart)\n");
#ifdef SCP_COMPLEX_FIELDS
    if (c->complex_gauge && c->g_gauge != 0 && n_gauge == 0)
        printf("  WARNING: init=sfa: no gauge columns; E seeded by Gauss projection\n");
#else
    (void)n_gauge;
#endif
    if (fabs(sfa->Lx - g->L) > 0.01)
        printf("  WARNING: SFA L=%.2f != config L=%.2f\n", sfa->Lx, g->L);

    free(buf);
    sfa_close(sfa);
}

static void init_from_exec(Grid *g, const Config *c) {
    printf("Init: exec '%s'\n", c->init_exec);
    char tmppath[512];
    snprintf(tmppath, sizeof(tmppath), "/tmp/scp_sim_init_%d.sfa", (int)getpid());
    char cmd[2048];
    snprintf(cmd, sizeof(cmd), "%s > %s", c->init_exec, tmppath);
    printf("  Running: %s\n", cmd);
    int ret = system(cmd);
    if (ret != 0) { fprintf(stderr, "FATAL: exec init failed (exit %d)\n", ret); exit(1); }
    Config tmp = *c;
    snprintf(tmp.init_sfa, sizeof(tmp.init_sfa), "%s", tmppath);
    tmp.init_frame = -1;
    init_from_sfa(g, &tmp);
    remove(tmppath);
}

static void init_template(Grid *g, const Config *c) {
    printf("Init: template '%s'\n", c->init_sfa);

    const int N = g->N, NN = N*N;
    const double L = g->L, dx = g->dx;
    const double k_bg = PI/L, omega_bg = sqrt(k_bg*k_bg + c->m2);

    for (long idx = 0; idx < g->N3; idx++) {
        int i = (int)(idx / NN), k = (int)(idx % N);
        double x = -L + i*dx, z = -L + k*dx;
        double A_bg_x = c->A_bg;
        if (c->bc_type == 1) {
            double frac = (x + L) / (2.0 * L);
            A_bg_x = c->gradient_A_high * (1.0 - frac) + c->gradient_A_low * frac;
        }
        for (int a = 0; a < NFIELDS; a++) {
            double ph_bg = k_bg*z + 2*PI*a/3.0;
            g->phi[a][idx] = A_bg_x * cos(ph_bg);
            g->phi_vel[a][idx] = omega_bg * A_bg_x * sin(ph_bg);
        }
    }

    SFA *tmpl = sfa_open(c->init_sfa);
    if (!tmpl) { fprintf(stderr, "FATAL: cannot open template '%s'\n", c->init_sfa); exit(1); }
    int TN = tmpl->Nx; double TL = tmpl->Lx, Tdx = 2.0*TL/(TN-1);
    long TN3 = (long)TN*TN*TN; int TNN = TN*TN;
    printf("  Template: %d^3, L=%.2f, dx=%.4f\n", TN, TL, Tdx);

    int frame = tmpl->total_frames - 1;
    void *buf = malloc(tmpl->frame_bytes);
    sfa_read_frame(tmpl, frame, buf);

    double *tphi[3], *tvel[3], *ttheta[3], *ttvel[3];
    for (int a = 0; a < 3; a++) {
        tphi[a]   = (double*)calloc(TN3, sizeof(double));
        tvel[a]   = (double*)calloc(TN3, sizeof(double));
        ttheta[a] = (double*)calloc(TN3, sizeof(double));
        ttvel[a]  = (double*)calloc(TN3, sizeof(double));
    }
    uint64_t off = 0;
    for (uint32_t col = 0; col < tmpl->n_columns; col++) {
        int dtype = tmpl->columns[col].dtype, sem = tmpl->columns[col].semantic, comp = tmpl->columns[col].component;
        int es = sfa_dtype_size[dtype]; uint8_t *src = (uint8_t*)buf + off;
        double *target = NULL;
        if (sem == SFA_POSITION && comp < 3) target = tphi[comp];
        else if (sem == SFA_ANGLE && comp < 3) target = ttheta[comp];
        else if (sem == SFA_VELOCITY && comp < 3) target = tvel[comp];
        else if (sem == SFA_VELOCITY && comp >= 3 && comp < 6) target = ttvel[comp-3];
        if (target) {
            if (dtype == SFA_F64) for (long i = 0; i < TN3; i++) target[i] = ((double*)src)[i];
            else if (dtype == SFA_F32) for (long i = 0; i < TN3; i++) target[i] = (double)((float*)src)[i];
            else if (dtype == SFA_F16) for (long i = 0; i < TN3; i++) target[i] = f16_to_f64(((uint16_t*)src)[i]);
        }
        off += (uint64_t)TN3 * es;
    }
    free(buf); sfa_close(tmpl);

    /* Map template voxels to destination by PHYSICAL coordinate (not index).
     * Template voxel (ti,tj,tk) is at physical position:
     *   tx = -TL + ti*Tdx, ty = -TL + tj*Tdx, tz = -TL + tk*Tdx
     * Find nearest destination voxel:
     *   gi = round((tx + L) / dx), etc.
     * This correctly handles different grid resolutions. */
    int placed = 0;
    for (int ti = 0; ti < TN; ti++) {
        double tx = -TL + ti*Tdx;
        int gi = (int)((tx + L) / dx + 0.5);
        if (gi < 0 || gi >= N) continue;
        for (int tj = 0; tj < TN; tj++) {
            double ty = -TL + tj*Tdx;
            int gj = (int)((ty + L) / dx + 0.5);
            if (gj < 0 || gj >= N) continue;
            for (int tk = 0; tk < TN; tk++) {
                double tz_phys = -TL + tk*Tdx;
                int gk = (int)((tz_phys + L) / dx + 0.5);
                if (gk < 0 || gk >= N) continue;
                long tidx = (long)ti*TNN + tj*TN + tk;
                long gidx = (long)gi*NN + gj*N + gk;
                for (int a = 0; a < NFIELDS; a++) {
                    double ph_bg_t = k_bg*tz_phys + 2*PI*a/3.0;
                    double bg_phi = c->A_bg*cos(ph_bg_t), bg_vel = omega_bg*c->A_bg*sin(ph_bg_t);
                    g->phi[a][gidx] += tphi[a][tidx] - bg_phi;
                    g->phi_vel[a][gidx] += tvel[a][tidx] - bg_vel;
                    g->theta[a][gidx] += ttheta[a][tidx];
                    g->theta_vel[a][gidx] += ttvel[a][tidx];
                }
                placed++;
            }
        }
    }
    for (int a = 0; a < 3; a++) { free(tphi[a]); free(tvel[a]); free(ttheta[a]); free(ttvel[a]); }
    printf("  Placed %d voxels from template\n", placed);
}

#ifdef SCP_COMPLEX_FIELDS
/* --- qball profile-table helpers (linear interp on strictly increasing r) --- */

/* generic table value: flat below rt[0], 'tail' above rt[np-1] */
static double qb_interp(const double *rt, const double *vt, size_t np,
                        double r, double tail) {
    if (r <= rt[0]) return vt[0];
    if (r >= rt[np-1]) return tail;
    size_t lo = 0, hi = np - 1;
    while (hi - lo > 1) {
        size_t mid = (lo + hi) / 2;
        if (rt[mid] <= r) lo = mid; else hi = mid;
    }
    double w = (r - rt[lo]) / (rt[hi] - rt[lo]);
    return vt[lo] + w * (vt[hi] - vt[lo]);
}

/* radial E: linear ramp to 0 at the center, inverse-square continuation
 * beyond the table (the tail carries the flux — SPEC §5.2) */
static double qb_er_interp(const double *rt, const double *et, size_t np, double r) {
    if (r <= rt[0])
        return (rt[0] > 0) ? et[0] * (r / rt[0]) : et[0];
    if (r >= rt[np-1]) {
        double R = rt[np-1];
        return et[np-1] * (R/r) * (R/r);
    }
    size_t lo = 0, hi = np - 1;
    while (hi - lo > 1) {
        size_t mid = (lo + hi) / 2;
        if (rt[mid] <= r) lo = mid; else hi = mid;
    }
    double w = (r - rt[lo]) / (rt[hi] - rt[lo]);
    return et[lo] + w * (et[hi] - et[lo]);
}

/* v66/v69 init=qball (v66 SPEC §5.1, v69 SPEC §5.2): symmetric Q-ball ansatz,
 *   u_a=f(r), v_a=0, udot_a=0, vdot_a=weff(r)*f(r); theta sector zero.
 * Profile columns: (r f) legacy, (r f Er) v69 contract, (r f Er weff) shooter
 * output (weff used directly). Gauged (complex_gauge && g!=0): E seeded on
 * radial link MIDPOINTS with 1/r^2 continuation; optional second ball via
 * qball2_* keys (matter + E superposed). delta[] is IGNORED. */
static void init_qball(Grid *g, const Config *c) {
    FILE *pf = fopen(c->qball_profile, "r");
    if (!pf) { fprintf(stderr, "FATAL: cannot open qball profile '%s'\n", c->qball_profile); exit(1); }
    size_t cap = 256, np = 0;
    int ncol = 0;
    double *rt = (double*)malloc(cap * sizeof(double));
    double *ft = (double*)malloc(cap * sizeof(double));
    double *et = (double*)malloc(cap * sizeof(double));
    double *wt = (double*)malloc(cap * sizeof(double));
    if (!rt || !ft || !et || !wt) { fprintf(stderr, "FATAL: qball profile alloc\n"); exit(1); }
    char line[1024];
    int lineno = 0;
    while (fgets(line, sizeof(line), pf)) {
        lineno++;
        char *p = line;
        while (*p == ' ' || *p == '\t') p++;
        if (*p == '#' || *p == '\n' || *p == '\r' || *p == '\0') continue;
        double r, f, er = 0, we = 0;
        int nf = sscanf(p, "%lf %lf %lf %lf", &r, &f, &er, &we);
        if (nf < 2) {
            fprintf(stderr, "FATAL: qball profile parse error at %s:%d\n",
                    c->qball_profile, lineno);
            exit(1);
        }
        if (ncol == 0) ncol = nf;                      /* first data line fixes layout */
        else if (nf < ncol) {
            fprintf(stderr, "FATAL: qball profile column count changed at %s:%d "
                    "(%d < %d)\n", c->qball_profile, lineno, nf, ncol);
            exit(1);
        }
        if (np == 0 && r < 0) {
            fprintf(stderr, "FATAL: qball profile first r must be >= 0 (got %f)\n", r);
            exit(1);
        }
        if (np > 0 && r <= rt[np-1]) {
            fprintf(stderr, "FATAL: qball profile r not strictly increasing at %s:%d\n",
                    c->qball_profile, lineno);
            exit(1);
        }
        if (np == cap) {
            cap *= 2;
            rt = (double*)realloc(rt, cap * sizeof(double));
            ft = (double*)realloc(ft, cap * sizeof(double));
            et = (double*)realloc(et, cap * sizeof(double));
            wt = (double*)realloc(wt, cap * sizeof(double));
            if (!rt || !ft || !et || !wt) { fprintf(stderr, "FATAL: qball profile realloc\n"); exit(1); }
        }
        rt[np] = r; ft[np] = f; et[np] = er; wt[np] = we; np++;
    }
    fclose(pf);
    if (np < 2) { fprintf(stderr, "FATAL: qball profile has <2 points\n"); exit(1); }

    const double om = c->qball_omega;
    const double GG = c->g_gauge;
    const int gauged = (c->complex_gauge && GG != 0.0);
    int have_er = (ncol >= 3);

    /* extra profile columns are honored ONLY when gauged: a 3/4-col profile
     * under complex_gauge=0 (or g==0) must reproduce the v66 seed exactly */
    if (!gauged || ncol == 2) {
        for (size_t q = 0; q < np; q++) wt[q] = om;   /* weff = omega everywhere */
        if (gauged && ncol == 2)
            printf("  WARNING: ungauged profile (no Er column): seeding "
                   "vdot=omega*f, E from projection; O(g^2) settling "
                   "transient expected\n");
        have_er = gauged && have_er;
    } else if (ncol == 3) {
        /* build weff = omega + g*a0(r) with a0(r) = -[int_r^R Er ds
         * + Er(R)*R] <= 0 (backward cumulative trapezoid + exact 1/r^2
         * tail). NOTE: SPEC §5.2 writes a0 with a + sign (a0>0), but that
         * is inconsistent with the kernel's own temporal-gauge dynamics:
         * stationarity of Phi=f(r)e^{i*Omega(r)t} under theta_dot=-g*a*E
         * requires Omega'(r)=+g*Er>=0 with Omega(inf)=omega, so the local
         * phase rate is BELOW omega inside the ball (Coulomb well a0<0).
         * GAUGE_DESIGN §2's original convention (E_i = d_i A_0 - d_t A_i,
         * A_0 a negative well) agrees, as does the shooter's 4-column
         * output (a0_GD(0)<0, weff0<omega). SPEC O1's E=-grad(A0)-dA/dt
         * "repair" is not invariant under SPEC §3's own gauge law and
         * §5.1/§5.2/§7.3(P4) inherited the flipped sign from it. */
        double a0 = -(et[np-1] * rt[np-1]);
        wt[np-1] = om + GG*a0;
        for (size_t q = np-1; q > 0; q--) {
            a0 -= 0.5*(et[q-1] + et[q])*(rt[q] - rt[q-1]);
            wt[q-1] = om + GG*a0;
        }
    }   /* ncol == 4 && gauged: weff straight from the shooter table */

    printf("Init: qball (%zu points, %d cols, r=[%.3f,%.3f], f0=%.6f, omega=%.4f, "
           "weff0=%.4f; delta[] ignored)\n",
           np, ncol, rt[0], rt[np-1], ft[0], om, wt[0]);
    if (fabs(ft[np-1]) > 1e-8)
        printf("  WARNING: profile truncated while still sizable (f_last=%.3e at r=%.3f)\n",
               ft[np-1], rt[np-1]);

    /* --- ball centers --- */
    const int N = g->N, NN = N * N;
    const double L = g->L, dx = g->dx;
    const double x0 = (c->qball_x0 >= 1e29) ? 0.0 : c->qball_x0;
    const double y0 = (c->qball_y0 >= 1e29) ? 0.0 : c->qball_y0;
    const double z0 = (c->qball_z0 >= 1e29) ? 0.0 : c->qball_z0;
    const int ball2 = (c->qball2_x0 < 1e29 && c->qball2_y0 < 1e29 && c->qball2_z0 < 1e29);
    const double x2 = ball2 ? c->qball2_x0 : 0, y2 = ball2 ? c->qball2_y0 : 0,
                 z2 = ball2 ? c->qball2_z0 : 0;
    const double sgn = (double)c->qball2_sign;
    const double cd = ball2 ? cos(c->qball2_phase) : 1.0;
    const double sd = ball2 ? sin(c->qball2_phase) : 0.0;
    if (ball2)
        printf("  ball1 center=(%.2f,%.2f,%.2f); ball2 center=(%.2f,%.2f,%.2f) "
               "sign=%+d phase=%.4f\n", x0, y0, z0, x2, y2, z2,
               c->qball2_sign, c->qball2_phase);

    /* --- matter fill --- */
    for (int i = 0; i < N; i++) { double x = -L + i*dx;
    for (int j = 0; j < N; j++) { double y = -L + j*dx;
    for (int k = 0; k < N; k++) { double z = -L + k*dx;
        long idx = (long)i*NN + j*N + k;
        double r = sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0) + (z-z0)*(z-z0));
        double f  = qb_interp(rt, ft, np, r, 0.0);
        double we = qb_interp(rt, wt, np, r, wt[np-1]);
        for (int a = 0; a < NFIELDS; a++) {
            g->phi[a][idx]        = f;        /* u_a = f(r) */
            g->phi_im[a][idx]     = 0.0;      /* v_a = 0 */
            g->phi_vel[a][idx]    = 0.0;      /* udot_a = 0 */
            g->phi_im_vel[a][idx] = we * f;   /* vdot_a = weff(r)*f(r) */
            /* theta sector (tu, tv, tudot, tvdot) stays zero (calloc'd) */
        }
        if (ball2) {
            double r2 = sqrt((x-x2)*(x-x2) + (y-y2)*(y-y2) + (z-z2)*(z-z2));
            double f2  = qb_interp(rt, ft, np, r2, 0.0);
            double we2 = qb_interp(rt, wt, np, r2, wt[np-1]);
            for (int a = 0; a < NFIELDS; a++) {
                g->phi[a][idx]        += f2*cd;
                g->phi_im[a][idx]     += f2*sd;
                g->phi_vel[a][idx]    += -sgn*we2*f2*sd;
                g->phi_im_vel[a][idx] +=  sgn*we2*f2*cd;
            }
        }
    }}}

    /* --- E seeding on radial link midpoints (gauged, Er available) --- */
    if (gauged && have_er) {
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < (long)N*NN; idx++) {
            int i=(int)(idx/NN), j=(int)((idx/N)%N), k=(int)(idx%N);
            double xc=-L+i*dx, yc=-L+j*dx, zc=-L+k*dx;
            for (int d = 0; d < 3; d++) {
                double mx = xc + (d==0 ? 0.5*dx : 0) - x0;
                double my = yc + (d==1 ? 0.5*dx : 0) - y0;
                double mz = zc + (d==2 ? 0.5*dx : 0) - z0;
                double rm = sqrt(mx*mx + my*my + mz*mz);
                double E = 0.0;
                if (rm > 1e-12) {
                    double Er = qb_er_interp(rt, et, np, rm);
                    double rhat = (d==0 ? mx : d==1 ? my : mz) / rm;
                    E += Er * rhat;
                }
                if (ball2) {
                    double m2x = xc + (d==0 ? 0.5*dx : 0) - x2;
                    double m2y = yc + (d==1 ? 0.5*dx : 0) - y2;
                    double m2z = zc + (d==2 ? 0.5*dx : 0) - z2;
                    double rm2 = sqrt(m2x*m2x + m2y*m2y + m2z*m2z);
                    if (rm2 > 1e-12) {
                        double Er2 = qb_er_interp(rt, et, np, rm2);
                        double rhat2 = (d==0 ? m2x : d==1 ? m2y : m2z) / rm2;
                        E += sgn * Er2 * rhat2;
                    }
                }
                g->Efield[d][idx] += E;     /* th links stay 0 (Coulomb curl-free) */
            }
        }
    }
    free(rt); free(ft); free(et); free(wt);
}
#endif /* SCP_COMPLEX_FIELDS */

static void do_init(Grid *g, const Config *c) {
    if      (!strcmp(c->init, "oscillon")) init_oscillon(g, c);
    else if (!strcmp(c->init, "braid"))    init_braid(g, c);
    else if (!strcmp(c->init, "sfa"))      init_from_sfa(g, c);
    else if (!strcmp(c->init, "exec"))     init_from_exec(g, c);
    else if (!strcmp(c->init, "template")) init_template(g, c);
#ifdef SCP_COMPLEX_FIELDS
    else if (!strcmp(c->init, "qball"))    init_qball(g, c);
#endif
    else { fprintf(stderr, "FATAL: unknown init mode '%s'\n", c->init); exit(1); }
}

/* ================================================================
   SFA KVMD metadata embedding — writes simulation parameters
   ================================================================ */

static void sfa_embed_kvmd(SFA *sfa, const Config *c) {
    char vN[32],vL[32],vT[32],vdt[32],vm[32],vmt[32],veta[32],vmu[32],vkappa[32];
    char vkh[32],vacs[32],vbh[32],vbcsw[32],vmode[32],via[32],vib[32],vkg[32];
    char vdw[32],vdr[32],vprec[32],vdelta[64],vbc[32],vgah[32],vgal[32],vgm[32];
    char vcplx[32],vqom[32];
    snprintf(vN,32,"%d",c->N); snprintf(vL,32,"%f",c->L);
    snprintf(vT,32,"%f",c->T); snprintf(vdt,32,"%f",c->dt_factor);
    snprintf(vm,32,"%f",sqrt(c->m2)); snprintf(vmt,32,"%f",sqrt(c->mtheta2));
    snprintf(veta,32,"%f",c->eta); snprintf(vmu,32,"%f",c->mu);
    snprintf(vkappa,32,"%f",c->kappa);
    snprintf(vkh,32,"%.6f",c->kappa_h);
    snprintf(vacs,32,"%.6f",c->alpha_cs); snprintf(vbh,32,"%.6f",c->beta_h);
    snprintf(vbcsw,32,"%.6f",c->bc_switch_time);
    snprintf(vmode,32,"%d",c->mode);
    snprintf(via,32,"%.6f",c->inv_alpha); snprintf(vib,32,"%.6f",c->inv_beta);
    snprintf(vkg,32,"%.6f",c->kappa_gamma);
    snprintf(vdw,32,"%.6f",c->damp_width); snprintf(vdr,32,"%.6f",c->damp_rate);
    snprintf(vprec,32,"%s", (const char*[]){"f16","f32","f64"}[c->precision]);
    snprintf(vdelta,64,"%.6f,%.6f,%.6f",c->delta[0],c->delta[1],c->delta[2]);
    snprintf(vbc,32,"%d",c->bc_type);
    snprintf(vgah,32,"%.6f",c->gradient_A_high);
    snprintf(vgal,32,"%.6f",c->gradient_A_low);
    snprintf(vgm,32,"%d",c->gradient_margin);
    /* v66: complex_phi MUST be embedded so the .sfa-restart path (scp_sim
     * <file>.sfa builds its config from KVMD) reallocates the imaginary
     * sector instead of silently dropping it; qball_omega for provenance.
     * The two keys are appended ONLY when complex_phi=1, keeping real-mode
     * output byte-identical to today (SPEC §7.1); files lacking the key
     * default to complex_phi=0, which is correct for all legacy files. */
    snprintf(vcplx,32,"%d",c->complex_phi);
    snprintf(vqom,32,"%.6f",c->qball_omega);
    /* v69: complex_gauge + g_gauge appended ONLY when complex_gauge=1 (29 keys)
     * so a .sfa-restart reallocates the gauge sector; legacy files default to 0. */
    char vcgauge[32],vgg[32];
    snprintf(vcgauge,32,"%d",c->complex_gauge);
    snprintf(vgg,32,"%.9f",c->g_gauge);
    const char *keys[] = {"N","L","T","dt_factor","m","m_theta","eta","mu","kappa",
                          "kappa_h","alpha_cs","beta_h","bc_switch_time",
                          "mode","inv_alpha","inv_beta","kappa_gamma",
                          "damp_width","damp_rate","precision","delta",
                          "bc_type","gradient_A_high","gradient_A_low","gradient_margin",
                          "complex_phi","qball_omega",
                          "complex_gauge","g_gauge"};
    const char *vals[] = {vN,vL,vT,vdt,vm,vmt,veta,vmu,vkappa,
                          vkh,vacs,vbh,vbcsw,
                          vmode,via,vib,vkg,vdw,vdr,vprec,vdelta,
                          vbc,vgah,vgal,vgm,
                          vcplx,vqom,
                          vcgauge,vgg};
    sfa_add_kvmd(sfa, 0, 0xFFFFFFFF, 0xFFFFFFFF, keys, vals,
                 c->complex_gauge ? 29 : (c->complex_phi ? 27 : 25));
}

#endif /* SCP_INIT_H */
