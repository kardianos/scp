/*  sfa_radial.c — v86 radial decomposition of the gauged complex fabric
 *
 *  Emits, per frame and per radial shell about a chosen center, the FULL
 *  energy and charge decomposition, mirroring scp_sim.c's
 *  compute_energy_complex_gauge() / compute_charges() term by term so the
 *  shell sums reproduce the kernel's diag totals (verified by the --check
 *  flag, which prints the volume-integrated totals next to the diag columns).
 *
 *  Rungs served:
 *    N10  — cloud-only energy vs omega_cl*|Q_cl| needs energy AND charge
 *           resolved in radius with a core cut (protocol:
 *           v86/council/grok45/N_BATTERY_REVIEW.md §1 N10; "without core
 *           subtraction N10 is invalid").
 *    N5   — throughput vs charge needs E_kin resolved per object/region.
 *    EX-2 — sponge-flux spectroscopy needs the radial energy flux split into
 *           matter and gauge (Poynting) channels at a shell inside the sponge.
 *
 *  Conventions taken VERBATIM from sfa/sim/scp_sim.c (do not "improve" them):
 *    coordinates      x = -L + i*dx      (kernel compute_charges)
 *    charge density   rho_Q  = sum_a (u_a vdot_a - v_a udot_a)
 *    kinetic          0.5 sum_a (udot^2 + vdot^2)
 *    covariant grad   DcU[d][a] = ((cP*u(+d) - sP*v(+d)) - (cM*u(-d) + sM*v(-d)))/(2dx)
 *                     DcV[d][a] = ((cP*v(+d) + sP*u(+d)) - (cM*v(-d) - sM*u(-d)))/(2dx)
 *                     with cP=cos(th_d[idx]), sP=sin(th_d[idx]),
 *                          cM=cos(th_d[idx-d]), sM=sin(th_d[idx-d])
 *    gradient energy  0.5 sum_a sum_d (DcU^2 + DcV^2)
 *    mass             0.5 m^2 sum_a (u^2+v^2)
 *    potential        (mu/2) s/(1+kappa s),  s = prod_a (u_a^2+v_a^2)
 *    EM               0.5|E|^2 + (1/(g^2 dx^4)) sum_p (1 - cos(plaquette_p))
 *    theta sector     same forms with m_theta^2
 *    coupling         -eta*(u.Re(DxTheta) + v.Im(DxTheta))   (zero at eta=0)
 *  All neighbour indexing is periodic, exactly as the kernel does it.
 *
 *  Derived fields this tool adds (not in the kernel):
 *    rho_Q+ / rho_Q-  positive / negative parts of rho_Q (cloud vs ball)
 *    rho2             sum_a (u^2+v^2)               (density norm)
 *    omega_loc        rho_Q / rho2 per shell        (density-weighted local clock)
 *    S_r matter       -(sum_a udot_a DcU[r][a] + vdot_a DcV[r][a]) projected on rhat
 *    S_r gauge        (E x B)_r,  B_k = plaquette_angle_k/(g dx^2)
 *
 *  Build: gcc -O3 -march=native -fopenmp -o sfa_radial sfa_radial.c -lzstd -lm
 *  Usage: sfa_radial file.sfa [--frames a,b,c|all] [--nbins K] [--rmax R]
 *                    [--center x,y,z] [--out prefix] [--g G] [--check]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define SFA_IMPLEMENTATION
#include <omp.h>
#include "../../sfa/format/sfa.h"

#define NF 3            /* components per sector */

static int find_col(const SFA *s, const char *name) {
    for (uint32_t c = 0; c < s->n_columns; c++)
        if (strncmp(s->columns[c].name, name, sizeof(s->columns[c].name)) == 0)
            return (int)c;
    return -1;
}

/* column base pointer inside a frame buffer (all columns are f32 here) */
static const float *colptr(const void *buf, const SFA *s, int ci) {
    uint64_t off = 0;
    for (int c = 0; c < ci; c++)
        off += s->N_total * sfa_dtype_size[s->columns[c].dtype];
    return (const float *)((const uint8_t *)buf + off);
}

static const float *need(const void *buf, const SFA *s, const char *name) {
    int ci = find_col(s, name);
    if (ci < 0) { fprintf(stderr, "error: column '%s' missing\n", name); exit(1); }
    if (s->columns[ci].dtype != SFA_F32) {
        fprintf(stderr, "error: column '%s' is not f32 (this tool needs the "
                        "restartable f32 output, not an f16 viewing copy)\n", name);
        exit(1);
    }
    return colptr(buf, s, ci);
}

typedef struct {
    double vol;
    double q, qpos, qneg, qa[NF], rho2;
    double e_kin, e_grad, e_mass, e_pot;
    double e_tkin, e_tgrad, e_tmass, e_coup;
    double e_elec, e_mag;
    double flux_matter, flux_gauge;
    double n;
} Shell;

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "usage: %s file.sfa [--frames a,b,c|all] [--nbins K] "
                        "[--rmax R] [--center x,y,z] [--out prefix] [--g G] "
                        "[--check]\n", argv[0]);
        return 1;
    }
    const char *path = argv[1];
    const char *outpref = "radial";
    char frames_arg[4096] = "all";
    int nbins = 120, check = 0;
    double rmax = -1, cx = 0, cy = 0, cz = 0, g_gauge = -1;

    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "--frames") && i + 1 < argc) snprintf(frames_arg, sizeof frames_arg, "%s", argv[++i]);
        else if (!strcmp(argv[i], "--nbins") && i + 1 < argc) nbins = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--rmax") && i + 1 < argc) rmax = atof(argv[++i]);
        else if (!strcmp(argv[i], "--out") && i + 1 < argc) outpref = argv[++i];
        else if (!strcmp(argv[i], "--g") && i + 1 < argc) g_gauge = atof(argv[++i]);
        else if (!strcmp(argv[i], "--check")) check = 1;
        else if (!strcmp(argv[i], "--center") && i + 1 < argc)
            sscanf(argv[++i], "%lf,%lf,%lf", &cx, &cy, &cz);
        else { fprintf(stderr, "unknown arg %s\n", argv[i]); return 1; }
    }

    SFA *s = sfa_open(path);
    if (!s) { fprintf(stderr, "cannot open %s\n", path); return 1; }

    /* ---- physics parameters from KVMD (never guessed) ------------------- */
    SFA_KVMDSet kv[16];
    int n_kv = sfa_read_kvmd(s, kv, 16);
    double m2 = 2.25, mtheta2 = 1.6 * 1.6, mu = -41.345, kappa = 50.0, eta = 0.0;
    int have_m = 0, have_mt = 0;
    for (int i = 0; i < n_kv; i++) for (int p = 0; p < kv[i].n_pairs; p++) {
        const char *k = kv[i].keys[p], *v = kv[i].values[p];
        if (!strcmp(k, "m"))        { m2 = atof(v) * atof(v); have_m = 1; }
        if (!strcmp(k, "m_theta"))  { mtheta2 = atof(v) * atof(v); have_mt = 1; }
        if (!strcmp(k, "mu"))       mu = atof(v);
        if (!strcmp(k, "kappa"))    kappa = atof(v);
        if (!strcmp(k, "eta"))      eta = atof(v);
        if (!strcmp(k, "g_gauge") && g_gauge < 0) g_gauge = atof(v);
    }
    if (g_gauge < 0) {
        fprintf(stderr, "error: g_gauge not in KVMD and --g not given; the "
                        "magnetic energy 1/(g^2 dx^4)(1-cos) cannot be formed\n");
        return 1;
    }
    int N = (int)s->Nx;
    double L = s->Lx;
    double dx = 2.0 * L / (N - 1);          /* kernel convention: x = -L + i*dx */
    double dV = dx * dx * dx;
    double idx1 = 1.0 / (2.0 * dx);
    double BPRE = (g_gauge > 0) ? 1.0 / (g_gauge * g_gauge * dx * dx * dx * dx) : 0.0;
    if (rmax < 0) rmax = L;
    double dr = rmax / nbins;

    fprintf(stderr, "sfa_radial: N=%d L=%g dx=%.6f g=%.4f m2=%.4f%s mtheta2=%.4f%s "
                    "mu=%g kappa=%g eta=%g\n", N, L, dx, g_gauge, m2,
            have_m ? "" : "(default)", mtheta2, have_mt ? "" : "(default)",
            mu, kappa, eta);
    fprintf(stderr, "            center=(%g,%g,%g) nbins=%d rmax=%g dr=%.4f\n",
            cx, cy, cz, nbins, rmax, dr);

    /* ---- frame list ------------------------------------------------------ */
    int nfr = (int)s->total_frames;
    int *flist = malloc(sizeof(int) * nfr);
    int nsel = 0;
    if (!strcmp(frames_arg, "all")) { for (int i = 0; i < nfr; i++) flist[nsel++] = i; }
    else {
        char *tok = strtok(frames_arg, ",");
        while (tok) { int f = atoi(tok); if (f >= 0 && f < nfr) flist[nsel++] = f; tok = strtok(NULL, ","); }
    }

    uint64_t bufsz = 0;
    for (uint32_t c = 0; c < s->n_columns; c++)
        bufsz += s->N_total * sfa_dtype_size[s->columns[c].dtype];
    void *buf = malloc(bufsz);
    if (!buf) { fprintf(stderr, "cannot allocate %.1f GB frame buffer\n", bufsz / 1e9); return 1; }

    char fn[1024];
    snprintf(fn, sizeof fn, "%s_shells.tsv", outpref);
    FILE *fo = fopen(fn, "w");
    fprintf(fo, "frame\tt\tr\tvol\tnvox\tQ\tQpos\tQneg\tQ0\tQ1\tQ2\trho2\t"
                "E_kin\tE_grad\tE_mass\tE_pot\tE_tkin\tE_tgrad\tE_tmass\tE_coup\t"
                "E_elec\tE_mag\tS_matter\tS_gauge\n");

    const int NN = N * N;
    const long N3 = (long)N * N * N;

    for (int fi = 0; fi < nsel; fi++) {
        int f = flist[fi];
        double t = sfa_frame_time(s, f);
        if (sfa_read_frame(s, f, buf) != 0) {
            fprintf(stderr, "frame %d read failed\n", f); continue;
        }
        const float *u[NF], *v[NF], *ud[NF], *vd[NF];
        const float *tu[NF], *tv[NF], *tud[NF], *tvd[NF];
        const char *ax = "xyz";
        char nm[32];
        for (int a = 0; a < NF; a++) {
            snprintf(nm, sizeof nm, "phi_%c", ax[a]);        u[a]  = need(buf, s, nm);
            snprintf(nm, sizeof nm, "phiim_%c", ax[a]);      v[a]  = need(buf, s, nm);
            snprintf(nm, sizeof nm, "phi_v%c", ax[a]);       ud[a] = need(buf, s, nm);
            snprintf(nm, sizeof nm, "phiim_v%c", ax[a]);     vd[a] = need(buf, s, nm);
            snprintf(nm, sizeof nm, "theta_%c", ax[a]);      tu[a] = need(buf, s, nm);
            snprintf(nm, sizeof nm, "thetaim_%c", ax[a]);    tv[a] = need(buf, s, nm);
            snprintf(nm, sizeof nm, "theta_v%c", ax[a]);     tud[a]= need(buf, s, nm);
            snprintf(nm, sizeof nm, "thetaim_v%c", ax[a]);   tvd[a]= need(buf, s, nm);
        }
        const float *th[3], *Ef[3];
        for (int d = 0; d < 3; d++) {
            snprintf(nm, sizeof nm, "th_%c", ax[d]); th[d] = need(buf, s, nm);
            snprintf(nm, sizeof nm, "E_%c", ax[d]);  Ef[d] = need(buf, s, nm);
        }

        int nthr = 1;
        #ifdef _OPENMP
        #pragma omp parallel
        { }
        nthr = omp_get_max_threads();
        #endif
        Shell *acc = calloc((size_t)nbins * nthr, sizeof(Shell));

        #pragma omp parallel
        {
            int tid = 0;
            #ifdef _OPENMP
            tid = omp_get_thread_num();
            #endif
            Shell *my = acc + (size_t)nbins * tid;
            #pragma omp for schedule(static)
            for (long idx = 0; idx < N3; idx++) {
                int i = (int)(idx / NN), j = (int)((idx / N) % N), k = (int)(idx % N);
                double x = -L + i * dx - cx, y = -L + j * dx - cy, z = -L + k * dx - cz;
                double rr = sqrt(x * x + y * y + z * z);
                int b = (int)(rr / dr);
                if (b >= nbins) continue;

                long np[3], nmm[3];
                np[0] = (long)((i + 1) % N) * NN + (long)j * N + k;
                nmm[0]= (long)((i - 1 + N) % N) * NN + (long)j * N + k;
                np[1] = (long)i * NN + (long)((j + 1) % N) * N + k;
                nmm[1]= (long)i * NN + (long)((j - 1 + N) % N) * N + k;
                np[2] = (long)i * NN + (long)j * N + (k + 1) % N;
                nmm[2]= (long)i * NN + (long)j * N + (k - 1 + N) % N;

                double DcU[3][6], DcV[3][6];
                for (int d = 0; d < 3; d++) {
                    double cP = cos(th[d][idx]),      sP = sin(th[d][idx]);
                    double cM = cos(th[d][nmm[d]]),   sM = sin(th[d][nmm[d]]);
                    for (int a = 0; a < NF; a++) {
                        double upn = u[a][np[d]],  vpn = v[a][np[d]];
                        double umn = u[a][nmm[d]], vmn = v[a][nmm[d]];
                        double tpn = tu[a][np[d]], wpn = tv[a][np[d]];
                        double tmn = tu[a][nmm[d]],wmn = tv[a][nmm[d]];
                        DcU[d][a]     = ((cP*upn - sP*vpn) - (cM*umn + sM*vmn)) * idx1;
                        DcV[d][a]     = ((cP*vpn + sP*upn) - (cM*vmn - sM*umn)) * idx1;
                        DcU[d][3 + a] = ((cP*tpn - sP*wpn) - (cM*tmn + sM*wmn)) * idx1;
                        DcV[d][3 + a] = ((cP*wpn + sP*tpn) - (cM*wmn - sM*tmn)) * idx1;
                    }
                }

                double sprod = 1.0, rho2 = 0.0, q = 0.0, qa[NF];
                double ekin = 0, egrad = 0, emass = 0;
                double etkin = 0, etgrad = 0, etmass = 0, ecoup = 0;
                double smat = 0;                     /* matter radial flux (unprojected) */
                double sx = 0, sy = 0, sz = 0;
                const int ci1[3] = {1,2,0}, ci2[3] = {2,0,1};
                for (int a = 0; a < NF; a++) {
                    double U = u[a][idx], V = v[a][idx];
                    double UD = ud[a][idx], VD = vd[a][idx];
                    double TU = tu[a][idx], TV = tv[a][idx];
                    double TUD = tud[a][idx], TVD = tvd[a][idx];
                    qa[a] = U * VD - V * UD;
                    q += qa[a];
                    double s2a = U * U + V * V;
                    rho2 += s2a;
                    sprod *= s2a;
                    ekin  += 0.5 * (UD * UD + VD * VD);
                    etkin += 0.5 * (TUD * TUD + TVD * TVD);
                    emass += 0.5 * m2 * s2a;
                    etmass+= 0.5 * mtheta2 * (TU * TU + TV * TV);
                    for (int d = 0; d < 3; d++) {
                        egrad  += 0.5 * (DcU[d][a] * DcU[d][a] + DcV[d][a] * DcV[d][a]);
                        etgrad += 0.5 * (DcU[d][3+a]*DcU[d][3+a] + DcV[d][3+a]*DcV[d][3+a]);
                    }
                    double reDxT = DcU[ci1[a]][3+ci2[a]] - DcU[ci2[a]][3+ci1[a]];
                    double imDxT = DcV[ci1[a]][3+ci2[a]] - DcV[ci2[a]][3+ci1[a]];
                    ecoup -= eta * (U * reDxT + V * imDxT);
                    /* matter energy flux T^{0d} = -(udot D_d u + vdot D_d v), summed
                     * over the phi and theta sectors */
                    sx -= UD * DcU[0][a] + VD * DcV[0][a] + TUD * DcU[0][3+a] + TVD * DcV[0][3+a];
                    sy -= UD * DcU[1][a] + VD * DcV[1][a] + TUD * DcU[1][3+a] + TVD * DcV[1][3+a];
                    sz -= UD * DcU[2][a] + VD * DcV[2][a] + TUD * DcU[2][3+a] + TVD * DcV[2][3+a];
                }
                double epot = (mu / 2.0) * sprod / (1.0 + kappa * sprod);

                double ex = Ef[0][idx], ey = Ef[1][idx], ez = Ef[2][idx];
                double eelec = 0.5 * (ex * ex + ey * ey + ez * ez);
                const int pa[3] = {0,1,2}, pb[3] = {1,2,0};
                double bb = 0.0, B[3];
                for (int p = 0; p < 3; p++) {
                    int a = pa[p], bx = pb[p];
                    double ang = th[a][idx] + th[bx][np[a]] - th[a][np[bx]] - th[bx][idx];
                    bb += 1.0 - cos(ang);
                    /* plaquette in the (a,b) plane is the field component along the
                     * third axis; pa/pb = (0,1),(1,2),(2,0) -> axes 2,0,1 */
                    B[(p + 2) % 3] = ang / (g_gauge * dx * dx);
                }
                double emag = BPRE * bb;
                /* Poynting S = E x B */
                double px = ey * B[2] - ez * B[1];
                double py = ez * B[0] - ex * B[2];
                double pz = ex * B[1] - ey * B[0];

                double inv = (rr > 1e-12) ? 1.0 / rr : 0.0;
                smat = (sx * x + sy * y + sz * z) * inv;
                double sgau = (px * x + py * y + pz * z) * inv;

                Shell *sh = my + b;
                sh->n++;
                sh->q += q; sh->rho2 += rho2;
                if (q > 0) sh->qpos += q; else sh->qneg += q;
                for (int a = 0; a < NF; a++) sh->qa[a] += qa[a];
                sh->e_kin += ekin; sh->e_grad += egrad; sh->e_mass += emass;
                sh->e_pot += epot; sh->e_tkin += etkin; sh->e_tgrad += etgrad;
                sh->e_tmass += etmass; sh->e_coup += ecoup;
                sh->e_elec += eelec; sh->e_mag += emag;
                sh->flux_matter += smat; sh->flux_gauge += sgau;
            }
        }

        /* reduce threads */
        Shell *tot = calloc(nbins, sizeof(Shell));
        for (int tid = 0; tid < nthr; tid++)
            for (int b = 0; b < nbins; b++) {
                Shell *a = acc + (size_t)nbins * tid + b, *o = tot + b;
                o->n += a->n; o->q += a->q; o->qpos += a->qpos; o->qneg += a->qneg;
                o->rho2 += a->rho2;
                for (int c = 0; c < NF; c++) o->qa[c] += a->qa[c];
                o->e_kin += a->e_kin; o->e_grad += a->e_grad; o->e_mass += a->e_mass;
                o->e_pot += a->e_pot; o->e_tkin += a->e_tkin; o->e_tgrad += a->e_tgrad;
                o->e_tmass += a->e_tmass; o->e_coup += a->e_coup;
                o->e_elec += a->e_elec; o->e_mag += a->e_mag;
                o->flux_matter += a->flux_matter; o->flux_gauge += a->flux_gauge;
            }

        double T[16]; memset(T, 0, sizeof T);
        for (int b = 0; b < nbins; b++) {
            Shell *o = tot + b;
            double r = (b + 0.5) * dr;
            double vol = o->n * dV;
            fprintf(fo,
                "%d\t%.6f\t%.6f\t%.6g\t%.0f\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t"
                "%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t"
                "%.10g\t%.10g\t%.10g\t%.10g\n",
                f, t, r, vol, o->n,
                o->q*dV, o->qpos*dV, o->qneg*dV,
                o->qa[0]*dV, o->qa[1]*dV, o->qa[2]*dV, o->rho2*dV,
                o->e_kin*dV, o->e_grad*dV, o->e_mass*dV, o->e_pot*dV,
                o->e_tkin*dV, o->e_tgrad*dV, o->e_tmass*dV, o->e_coup*dV,
                o->e_elec*dV, o->e_mag*dV,
                /* fluxes are surface densities: report the shell-integrated
                 * value divided by the shell's mean area, i.e. the mean radial
                 * flux density on that sphere */
                (r > 0 && o->n > 0) ? o->flux_matter*dV / (4.0*M_PI*r*r*dr) : 0.0,
                (r > 0 && o->n > 0) ? o->flux_gauge *dV / (4.0*M_PI*r*r*dr) : 0.0);
            T[0]+=o->q*dV; T[1]+=o->e_kin*dV; T[2]+=o->e_grad*dV; T[3]+=o->e_mass*dV;
            T[4]+=o->e_pot*dV; T[5]+=o->e_tkin*dV; T[6]+=o->e_tgrad*dV;
            T[7]+=o->e_tmass*dV; T[8]+=o->e_coup*dV; T[9]+=o->e_elec*dV;
            T[10]+=o->e_mag*dV;
        }
        if (check) {
            double etot = T[1]+T[2]+T[3]+T[4]+T[5]+T[6]+T[7]+T[8]+T[9]+T[10];
            fprintf(stderr,
                "frame %2d t=%8.2f  (sums over r<%g only)\n"
                "   Q_phi=%.6e  E_phi_kin=%.6e E_grad=%.6e E_mass=%.6e E_pot=%.6e\n"
                "   E_theta_kin=%.6e E_tgrad=%.6e E_tmass=%.6e E_coupling=%.6e\n"
                "   E_em=%.6e (elec %.6e + mag %.6e)  E_total=%.6e\n",
                f, t, rmax, T[0], T[1], T[2], T[3], T[4], T[5], T[6], T[7], T[8],
                T[9]+T[10], T[9], T[10], etot);
        }
        free(acc); free(tot);
        fprintf(stderr, "frame %d/%d done (t=%.1f)\n", fi + 1, nsel, t);
    }
    fclose(fo);
    fprintf(stderr, "wrote %s_shells.tsv\n", outpref);
    free(buf); free(flist);
    return 0;
}
