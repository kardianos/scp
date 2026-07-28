/*  sfa_momentum.c — v86 momentum / stress-flux metrology for N6 and N7
 *
 *  Measures, per frame, the quantities the inertia lock needs WITHOUT ever
 *  assuming a relation between energy, charge and inertia:
 *
 *    Q, E                            (as scp_sim.c defines them)
 *    P_i = integral T^{0i} dV        (field momentum, matter + gauge Poynting)
 *    X_com                           (charge-weighted centre, per half-space)
 *    F_i = -integral T^{ji} n_j dA   (momentum flux through a plane: the
 *                                     MEASURED force on the region beyond it)
 *
 *  and reports them globally and for each half-space of a splitting plane, so
 *  a two-body run yields per-ball Q, P, X, and the interaction force, all from
 *  the fields.
 *
 *  WHY THIS EXISTS (v86/PART0_RESULTS.md, N7 design). The analytic Coulomb
 *  force g^2 Q1 Q2 / 4 pi D^2 is a continuum, static, point-monopole,
 *  infinite-volume formula. Its lattice form factor, finite-size and
 *  boundary-image errors are the SAME few percent as eps -- the very gap that
 *  separates the hypotheses M = E and M = Q omega. So the force must be
 *  measured, not assumed. Likewise the momentum: v70's "hbar_eff = p/k"
 *  used p = gamma E_0 v, which makes p/k = E_0/omega IDENTICALLY and is
 *  therefore not an independent leg of the hbar triple. Integrating T^{0i}
 *  makes it independent.
 *
 *  Conventions, taken from scp_sim.c and from T^{mu nu} = d^mu phi d^nu phi
 *  - g^{mu nu} L with signature (+,-,-,-):
 *      T^{00} = the kernel's energy density (verified against diag by
 *               sfa_radial.c to 1e-6)
 *      T^{0i} = -sum_a (phidot_a d_i phi_a)          [matter, both sectors]
 *             + (E x B)_i                            [gauge]
 *      T^{ji} =  sum_a (d_j phi_a d_i phi_a) + delta_ji * L_matter
 *             + [ -(E_j E_i + B_j B_i) + (1/2) delta_ji (E^2 + B^2) ]
 *      L_matter = sum_a (1/2)(phidot^2 - |D phi|^2) - m^2/2 |phi|^2 - Vt(s)
 *  Spatial derivatives are the kernel's covariant link differences when the
 *  gauge columns are present, plain central differences otherwise.
 *
 *  The SIGNS are not asserted -- they are VALIDATED at analysis time by the
 *  momentum-balance residual R_PF = |dP/dt - F| / |F| (grok's kill criterion).
 *  A sign error shows up as R_PF of order 2, not as a silent bias.
 *
 *  Handles f16 (viewing copies) and f32 (restartable output), gauged and
 *  ungauged files.
 *
 *  Build: gcc -O3 -march=native -fopenmp -o sfa_momentum sfa_momentum.c -lzstd -lm
 *  Usage: sfa_momentum file.sfa [--frames all|a,b,c] [--split X] [--out pfx]
 *                      [--g G] [--planes x1,x2,...]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#define SFA_IMPLEMENTATION
#include "../../sfa/format/sfa.h"

#define NF 3

static int find_col(const SFA *s, const char *name) {
    for (uint32_t c = 0; c < s->n_columns; c++)
        if (strncmp(s->columns[c].name, name, sizeof(s->columns[c].name)) == 0)
            return (int)c;
    return -1;
}

/* Decode any column (f16 or f32) into a float array. */
static float *getcol(const void *buf, const SFA *s, const char *name, int required) {
    int ci = find_col(s, name);
    if (ci < 0) {
        if (required) { fprintf(stderr, "error: column '%s' missing\n", name); exit(1); }
        return NULL;
    }
    uint64_t N3 = s->N_total, off = 0;
    for (int c = 0; c < ci; c++) off += N3 * sfa_dtype_size[s->columns[c].dtype];
    const uint8_t *src = (const uint8_t *)buf + off;
    float *a = (float *)malloc(N3 * sizeof(float));
    if (!a) { fprintf(stderr, "error: malloc\n"); exit(1); }
    if (s->columns[ci].dtype == SFA_F32) {
        memcpy(a, src, N3 * sizeof(float));
    } else if (s->columns[ci].dtype == SFA_F16) {
        const uint16_t *h = (const uint16_t *)src;
        #pragma omp parallel for schedule(static)
        for (uint64_t i = 0; i < N3; i++) {
            uint16_t v = h[i];
            uint16_t sign = v & 0x8000; int e = (v >> 10) & 0x1F; uint16_t mant = v & 0x3FF;
            float f;
            if (e == 0) { f = 0.0f; }
            else if (e == 31) { f = sign ? -1e30f : 1e30f; }
            else { uint32_t x = ((uint32_t)sign << 16) | ((uint32_t)(e - 15 + 127) << 23)
                               | ((uint32_t)mant << 13); memcpy(&f, &x, 4); }
            a[i] = f;
        }
    } else {
        fprintf(stderr, "error: column '%s' has unsupported dtype\n", name); exit(1);
    }
    return a;
}

typedef struct {
    double Q, E, rho2, m0;          /* m0 = charge magnitude, for the COM weight */
    double P[3], X[3];
    double n;
} Reg;

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "usage: %s file.sfa [--frames all|a,b,c] [--split X] "
                        "[--out pfx] [--g G] [--planes x1,x2,...]\n", argv[0]);
        return 1;
    }
    const char *path = argv[1], *outpref = "momentum";
    char frames_arg[4096] = "all", planes_arg[1024] = "";
    double split = 0.0, g_gauge = -1;
    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "--frames") && i + 1 < argc) snprintf(frames_arg, sizeof frames_arg, "%s", argv[++i]);
        else if (!strcmp(argv[i], "--split") && i + 1 < argc) split = atof(argv[++i]);
        else if (!strcmp(argv[i], "--out") && i + 1 < argc) outpref = argv[++i];
        else if (!strcmp(argv[i], "--g") && i + 1 < argc) g_gauge = atof(argv[++i]);
        else if (!strcmp(argv[i], "--planes") && i + 1 < argc) snprintf(planes_arg, sizeof planes_arg, "%s", argv[++i]);
        else { fprintf(stderr, "unknown arg %s\n", argv[i]); return 1; }
    }

    SFA *s = sfa_open(path);
    if (!s) { fprintf(stderr, "cannot open %s\n", path); return 1; }

    SFA_KVMDSet kv[16];
    int n_kv = sfa_read_kvmd(s, kv, 16);
    double m2 = 2.25, mt2 = 0.0, mu = -41.345, kappa = 50.0;
    for (int i = 0; i < n_kv; i++) for (int p = 0; p < kv[i].n_pairs; p++) {
        const char *k = kv[i].keys[p], *v = kv[i].values[p];
        if (!strcmp(k, "m")) m2 = atof(v) * atof(v);
        if (!strcmp(k, "m_theta")) mt2 = atof(v) * atof(v);
        if (!strcmp(k, "mu")) mu = atof(v);
        if (!strcmp(k, "kappa")) kappa = atof(v);
        if (!strcmp(k, "g_gauge") && g_gauge < 0) g_gauge = atof(v);
    }
    int N = (int)s->Nx, NN = N * N;
    long N3 = (long)N * N * N;
    double L = s->Lx, dx = 2.0 * L / (N - 1), dV = dx * dx * dx, idx1 = 1.0 / (2.0 * dx);

    /* planes for the stress-flux force */
    double planes[16]; int nplanes = 0;
    if (planes_arg[0]) {
        char *tok = strtok(planes_arg, ",");
        while (tok && nplanes < 16) { planes[nplanes++] = atof(tok); tok = strtok(NULL, ","); }
    } else { planes[nplanes++] = split; }

    int nfr = (int)s->total_frames;
    int *fl = malloc(sizeof(int) * nfr); int nsel = 0;
    if (!strcmp(frames_arg, "all")) { for (int i = 0; i < nfr; i++) fl[nsel++] = i; }
    else { char *t = strtok(frames_arg, ","); while (t) { int f = atoi(t); if (f >= 0 && f < nfr) fl[nsel++] = f; t = strtok(NULL, ","); } }

    uint64_t bufsz = 0;
    for (uint32_t c = 0; c < s->n_columns; c++) bufsz += s->N_total * sfa_dtype_size[s->columns[c].dtype];
    void *buf = malloc(bufsz);
    if (!buf) { fprintf(stderr, "cannot allocate %.1f GB\n", bufsz / 1e9); return 1; }

    char fn[1024]; snprintf(fn, sizeof fn, "%s_mom.tsv", outpref);
    FILE *fo = fopen(fn, "w");
    fprintf(fo, "frame\tt\tregion\tQ\tE\tPx\tPy\tPz\tXcom\tYcom\tZcom\trho2\tnvox\n");
    char fn2[1024]; snprintf(fn2, sizeof fn2, "%s_flux.tsv", outpref);
    FILE *ff = fopen(fn2, "w");
    fprintf(ff, "frame\tt\tplane_x\tFx\tFy\tFz\tFx_matter\tFx_gauge\n");

    int has_gauge = (find_col(s, "th_x") >= 0) && (find_col(s, "E_x") >= 0);
    fprintf(stderr, "sfa_momentum: N=%d L=%g dx=%.6f gauge=%s g=%.4f m2=%.4f mt2=%.4f\n",
            N, L, dx, has_gauge ? "yes" : "no", g_gauge, m2, mt2);

    for (int fi = 0; fi < nsel; fi++) {
        int f = fl[fi];
        double t = sfa_frame_time(s, f);
        if (sfa_read_frame(s, f, buf) != 0) { fprintf(stderr, "frame %d read failed\n", f); continue; }
        const char *ax = "xyz"; char nm[32];
        float *u[NF], *v[NF], *ud[NF], *vd[NF], *tu[NF], *tv[NF], *tud[NF], *tvd[NF];
        for (int a = 0; a < NF; a++) {
            snprintf(nm, sizeof nm, "phi_%c", ax[a]);      u[a]  = getcol(buf, s, nm, 1);
            snprintf(nm, sizeof nm, "phiim_%c", ax[a]);    v[a]  = getcol(buf, s, nm, 1);
            snprintf(nm, sizeof nm, "phi_v%c", ax[a]);     ud[a] = getcol(buf, s, nm, 1);
            snprintf(nm, sizeof nm, "phiim_v%c", ax[a]);   vd[a] = getcol(buf, s, nm, 1);
            snprintf(nm, sizeof nm, "theta_%c", ax[a]);    tu[a] = getcol(buf, s, nm, 0);
            snprintf(nm, sizeof nm, "thetaim_%c", ax[a]);  tv[a] = getcol(buf, s, nm, 0);
            snprintf(nm, sizeof nm, "theta_v%c", ax[a]);   tud[a]= getcol(buf, s, nm, 0);
            snprintf(nm, sizeof nm, "thetaim_v%c", ax[a]); tvd[a]= getcol(buf, s, nm, 0);
        }
        float *th[3] = {NULL,NULL,NULL}, *Ef[3] = {NULL,NULL,NULL};
        if (has_gauge) for (int d = 0; d < 3; d++) {
            snprintf(nm, sizeof nm, "th_%c", ax[d]); th[d] = getcol(buf, s, nm, 1);
            snprintf(nm, sizeof nm, "E_%c", ax[d]);  Ef[d] = getcol(buf, s, nm, 1);
        }

        int nthr = omp_get_max_threads();
        Reg *acc = calloc((size_t)3 * nthr, sizeof(Reg));   /* 0=all 1=left 2=right */
        double *pl_f = calloc((size_t)nplanes * 3 * nthr, sizeof(double));
        double *pl_mg = calloc((size_t)nplanes * 2 * nthr, sizeof(double));

        #pragma omp parallel
        {
            int tid = omp_get_thread_num();
            Reg *my = acc + 3 * tid;
            double *myf = pl_f + (size_t)nplanes * 3 * tid;
            double *mymg = pl_mg + (size_t)nplanes * 2 * tid;
            #pragma omp for schedule(static)
            for (long idx = 0; idx < N3; idx++) {
                int i = (int)(idx / NN), j = (int)((idx / N) % N), k = (int)(idx % N);
                double x = -L + i * dx, y = -L + j * dx, z = -L + k * dx;
                long np[3], nmm[3];
                np[0]=(long)((i+1)%N)*NN+(long)j*N+k;  nmm[0]=(long)((i-1+N)%N)*NN+(long)j*N+k;
                np[1]=(long)i*NN+(long)((j+1)%N)*N+k;  nmm[1]=(long)i*NN+(long)((j-1+N)%N)*N+k;
                np[2]=(long)i*NN+(long)j*N+(k+1)%N;    nmm[2]=(long)i*NN+(long)j*N+(k-1+N)%N;

                double DU[3][6], DV[3][6];
                for (int d = 0; d < 3; d++) {
                    double cP=1, sP=0, cM=1, sM=0;
                    if (has_gauge) { cP=cos(th[d][idx]); sP=sin(th[d][idx]);
                                     cM=cos(th[d][nmm[d]]); sM=sin(th[d][nmm[d]]); }
                    for (int a = 0; a < NF; a++) {
                        double upn=u[a][np[d]], vpn=v[a][np[d]];
                        double umn=u[a][nmm[d]], vmn=v[a][nmm[d]];
                        DU[d][a] = ((cP*upn - sP*vpn) - (cM*umn + sM*vmn))*idx1;
                        DV[d][a] = ((cP*vpn + sP*upn) - (cM*vmn - sM*umn))*idx1;
                        if (tu[a]) {
                            double tpn=tu[a][np[d]], wpn=tv[a][np[d]];
                            double tmn=tu[a][nmm[d]], wmn=tv[a][nmm[d]];
                            DU[d][3+a] = ((cP*tpn - sP*wpn) - (cM*tmn + sM*wmn))*idx1;
                            DV[d][3+a] = ((cP*wpn + sP*tpn) - (cM*wmn - sM*tmn))*idx1;
                        } else { DU[d][3+a]=0; DV[d][3+a]=0; }
                    }
                }
                double q=0, rho2=0, ekin=0, egrad=0, emass=0, sprod=1.0;
                double Pm[3]={0,0,0};
                double gradsq = 0.0, kinsq = 0.0;
                double Sij[3];                     /* T^{j x} matter, j=0,1,2 */
                double Six[3]={0,0,0}, Siy[3]={0,0,0}, Siz[3]={0,0,0};
                for (int a = 0; a < NF; a++) {
                    double U=u[a][idx], V=v[a][idx], UD=ud[a][idx], VD=vd[a][idx];
                    double TU = tu[a]?tu[a][idx]:0, TV = tv[a]?tv[a][idx]:0;
                    double TUD= tud[a]?tud[a][idx]:0, TVD= tvd[a]?tvd[a][idx]:0;
                    q += U*VD - V*UD;
                    double s2a = U*U + V*V; rho2 += s2a; sprod *= s2a;
                    ekin  += 0.5*(UD*UD + VD*VD) + 0.5*(TUD*TUD + TVD*TVD);
                    emass += 0.5*m2*s2a + 0.5*mt2*(TU*TU + TV*TV);
                    kinsq += UD*UD + VD*VD + TUD*TUD + TVD*TVD;
                    for (int d = 0; d < 3; d++) {
                        egrad += 0.5*(DU[d][a]*DU[d][a] + DV[d][a]*DV[d][a]
                                    + DU[d][3+a]*DU[d][3+a] + DV[d][3+a]*DV[d][3+a]);
                        gradsq += DU[d][a]*DU[d][a] + DV[d][a]*DV[d][a]
                                + DU[d][3+a]*DU[d][3+a] + DV[d][3+a]*DV[d][3+a];
                        Pm[d] -= UD*DU[d][a] + VD*DV[d][a]
                               + TUD*DU[d][3+a] + TVD*DV[d][3+a];
                    }
                    /* matter stress T^{j i} = sum d_j phi d_i phi + delta_ji L */
                    for (int d = 0; d < 3; d++) {
                        Six[d] += DU[d][a]*DU[0][a] + DV[d][a]*DV[0][a]
                                + DU[d][3+a]*DU[0][3+a] + DV[d][3+a]*DV[0][3+a];
                        Siy[d] += DU[d][a]*DU[1][a] + DV[d][a]*DV[1][a]
                                + DU[d][3+a]*DU[1][3+a] + DV[d][3+a]*DV[1][3+a];
                        Siz[d] += DU[d][a]*DU[2][a] + DV[d][a]*DV[2][a]
                                + DU[d][3+a]*DU[2][3+a] + DV[d][3+a]*DV[2][3+a];
                    }
                }
                double epot = (mu/2.0)*sprod/(1.0 + kappa*sprod);
                double Lm = 0.5*kinsq - 0.5*gradsq - emass - epot;
                Six[0] += Lm; Siy[1] += Lm; Siz[2] += Lm;

                double eelec=0, emag=0, B[3]={0,0,0};
                if (has_gauge) {
                    double ex=Ef[0][idx], ey=Ef[1][idx], ez=Ef[2][idx];
                    eelec = 0.5*(ex*ex+ey*ey+ez*ez);
                    const int pa[3]={0,1,2}, pb[3]={1,2,0};
                    double bb=0;
                    for (int p2=0;p2<3;p2++) {
                        int a2=pa[p2], b2=pb[p2];
                        double ang = th[a2][idx] + th[b2][np[a2]] - th[a2][np[b2]] - th[b2][idx];
                        bb += 1.0 - cos(ang);
                        B[(p2+2)%3] = ang/(g_gauge*dx*dx);
                    }
                    emag = bb/(g_gauge*g_gauge*dx*dx*dx*dx);
                    Pm[0] += ey*B[2]-ez*B[1];
                    Pm[1] += ez*B[0]-ex*B[2];
                    Pm[2] += ex*B[1]-ey*B[0];
                    double EE[3]={ex,ey,ez};
                    double half = 0.5*(ex*ex+ey*ey+ez*ez + B[0]*B[0]+B[1]*B[1]+B[2]*B[2]);
                    for (int d=0; d<3; d++) {
                        Six[d] += -(EE[d]*EE[0] + B[d]*B[0]) + (d==0?half:0.0);
                        Siy[d] += -(EE[d]*EE[1] + B[d]*B[1]) + (d==1?half:0.0);
                        Siz[d] += -(EE[d]*EE[2] + B[d]*B[2]) + (d==2?half:0.0);
                    }
                }
                double e = ekin + egrad + emass + epot + eelec + emag;

                int regs[2] = {0, (x < split ? 1 : 2)};
                for (int rr = 0; rr < 2; rr++) {
                    Reg *R = my + regs[rr];
                    R->Q += q; R->E += e; R->rho2 += rho2; R->n += 1;
                    R->m0 += fabs(q);
                    R->P[0] += Pm[0]; R->P[1] += Pm[1]; R->P[2] += Pm[2];
                    R->X[0] += fabs(q)*x; R->X[1] += fabs(q)*y; R->X[2] += fabs(q)*z;
                }
                /* stress flux through the x = plane surfaces: accumulate on the
                 * voxel slab immediately left of each plane (n_j = +xhat) */
                for (int p2 = 0; p2 < nplanes; p2++) {
                    if (x <= planes[p2] && x + dx > planes[p2]) {
                        myf[3*p2+0] += Six[0]; myf[3*p2+1] += Siy[0]; myf[3*p2+2] += Siz[0];
                        double mx = Six[0];
                        if (has_gauge) {
                            double ex=Ef[0][idx], ey=Ef[1][idx], ez=Ef[2][idx];
                            double half = 0.5*(ex*ex+ey*ey+ez*ez + B[0]*B[0]+B[1]*B[1]+B[2]*B[2]);
                            double gx = -(ex*ex + B[0]*B[0]) + half;
                            mymg[2*p2+1] += gx; mymg[2*p2+0] += mx - gx;
                        } else { mymg[2*p2+0] += mx; }
                    }
                }
            }
        }
        Reg tot[3]; memset(tot, 0, sizeof tot);
        for (int tid = 0; tid < nthr; tid++) for (int rr = 0; rr < 3; rr++) {
            Reg *a = acc + 3*tid + rr, *o = tot + rr;
            o->Q+=a->Q; o->E+=a->E; o->rho2+=a->rho2; o->n+=a->n; o->m0+=a->m0;
            for (int d=0; d<3; d++) { o->P[d]+=a->P[d]; o->X[d]+=a->X[d]; }
        }
        const char *rn[3] = {"all", "left", "right"};
        for (int rr = 0; rr < 3; rr++) {
            Reg *o = tot + rr;
            double w = (o->m0 > 1e-300) ? o->m0 : 1.0;
            fprintf(fo, "%d\t%.6f\t%s\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t"
                        "%.10g\t%.10g\t%.10g\t%.10g\t%.0f\n",
                    f, t, rn[rr], o->Q*dV, o->E*dV,
                    o->P[0]*dV, o->P[1]*dV, o->P[2]*dV,
                    o->X[0]/w, o->X[1]/w, o->X[2]/w, o->rho2*dV, o->n);
        }
        double dA = dx * dx;
        for (int p2 = 0; p2 < nplanes; p2++) {
            double Fx=0, Fy=0, Fz=0, Fm=0, Fg=0;
            for (int tid = 0; tid < nthr; tid++) {
                double *myf = pl_f + (size_t)nplanes*3*tid;
                double *mymg = pl_mg + (size_t)nplanes*2*tid;
                Fx += myf[3*p2+0]; Fy += myf[3*p2+1]; Fz += myf[3*p2+2];
                Fm += mymg[2*p2+0]; Fg += mymg[2*p2+1];
            }
            fprintf(ff, "%d\t%.6f\t%.6f\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\n",
                    f, t, planes[p2], -Fx*dA, -Fy*dA, -Fz*dA, -Fm*dA, -Fg*dA);
        }
        for (int a = 0; a < NF; a++) {
            free(u[a]); free(v[a]); free(ud[a]); free(vd[a]);
            free(tu[a]); free(tv[a]); free(tud[a]); free(tvd[a]);
        }
        if (has_gauge) for (int d = 0; d < 3; d++) { free(th[d]); free(Ef[d]); }
        free(acc); free(pl_f); free(pl_mg);
        fprintf(stderr, "frame %d/%d (t=%.2f)\n", fi + 1, nsel, t);
    }
    fclose(fo); fclose(ff);
    fprintf(stderr, "wrote %s_mom.tsv and %s_flux.tsv\n", outpref, outpref);
    free(buf); free(fl);
    return 0;
}
