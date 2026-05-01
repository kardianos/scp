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

    int loaded[12] = {0};
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

    int n_fields = 0, n_vels = 0;
    for (int i = 0; i < 6; i++) n_fields += loaded[i];
    for (int i = 6; i < 12; i++) n_vels += loaded[i];
    printf("  Loaded: %d/6 field arrays, %d/6 velocity arrays\n", n_fields, n_vels);
    if (n_vels == 0)
        printf("  WARNING: no velocity data — starting from rest (cold restart)\n");
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

static void do_init(Grid *g, const Config *c) {
    if      (!strcmp(c->init, "oscillon")) init_oscillon(g, c);
    else if (!strcmp(c->init, "braid"))    init_braid(g, c);
    else if (!strcmp(c->init, "sfa"))      init_from_sfa(g, c);
    else if (!strcmp(c->init, "exec"))     init_from_exec(g, c);
    else if (!strcmp(c->init, "template")) init_template(g, c);
    else { fprintf(stderr, "FATAL: unknown init mode '%s'\n", c->init); exit(1); }
}

/* ================================================================
   SFA KVMD metadata embedding — writes simulation parameters
   ================================================================ */

static void sfa_embed_kvmd(SFA *sfa, const Config *c) {
    char vN[32],vL[32],vT[32],vdt[32],vm[32],vmt[32],veta[32],vmu[32],vkappa[32];
    char vkh[32],vacs[32],vbh[32],vbcsw[32],vmode[32],via[32],vib[32],vkg[32];
    char vdw[32],vdr[32],vprec[32],vdelta[64],vbc[32],vgah[32],vgal[32],vgm[32];
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
    const char *keys[] = {"N","L","T","dt_factor","m","m_theta","eta","mu","kappa",
                          "kappa_h","alpha_cs","beta_h","bc_switch_time",
                          "mode","inv_alpha","inv_beta","kappa_gamma",
                          "damp_width","damp_rate","precision","delta",
                          "bc_type","gradient_A_high","gradient_A_low","gradient_margin"};
    const char *vals[] = {vN,vL,vT,vdt,vm,vmt,veta,vmu,vkappa,
                          vkh,vacs,vbh,vbcsw,
                          vmode,via,vib,vkg,vdw,vdr,vprec,vdelta,
                          vbc,vgah,vgal,vgm};
    sfa_add_kvmd(sfa, 0, 0xFFFFFFFF, 0xFFFFFFFF, keys, vals, 25);
}

#endif /* SCP_INIT_H */
