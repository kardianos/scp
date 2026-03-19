/*  analyze.c — Braid detection and characterization from field snapshots
 *
 *  Reads binary field snapshots (N, L, t, phi[3][N³]) and produces:
 *  1. Braid detection (connected components of high energy density)
 *  2. Per-braid properties (position, mass, momentum, size, winding, binding)
 *  3. Field characterization (radial profiles, energy flux, background stats)
 *  4. Inter-braid properties (separations, relative velocity)
 *
 *  Build: gcc -O3 -march=native -fopenmp -o analyze src/analyze.c -lm
 *  Usage: ./analyze field_t0100.bin [field_t0200.bin ...]
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

#define NFIELDS 3
#define PI 3.14159265358979323846
#define MAX_BRAIDS 100

static double MU = -41.345, KAPPA = 50.0, MASS2 = 2.25;

/* ================================================================
   Field snapshot
   ================================================================ */

typedef struct {
    int N; double L, t, dx;
    long N3;
    double *phi[NFIELDS];
    /* Derived fields (computed once) */
    double *rho;       /* energy density */
    double *absP;      /* |triple product| */
    double *binding_w; /* 1/(1+|P|/P_thresh) */
    int *label;        /* braid label per cell (-1 = background) */
} Snapshot;

static Snapshot *load_snapshot(const char *path) {
    FILE *fp = fopen(path, "rb");
    if (!fp) { fprintf(stderr, "Cannot open %s\n", path); return NULL; }

    Snapshot *s = calloc(1, sizeof(Snapshot));
    fread(&s->N, sizeof(int), 1, fp);
    fread(&s->L, sizeof(double), 1, fp);
    fread(&s->t, sizeof(double), 1, fp);
    s->N3 = (long)s->N * s->N * s->N;
    s->dx = 2.0 * s->L / (s->N - 1);

    for (int a = 0; a < NFIELDS; a++) {
        s->phi[a] = malloc(s->N3 * sizeof(double));
        fread(s->phi[a], sizeof(double), s->N3, fp);
    }
    fclose(fp);

    /* Allocate derived fields */
    s->rho = calloc(s->N3, sizeof(double));
    s->absP = calloc(s->N3, sizeof(double));
    s->binding_w = calloc(s->N3, sizeof(double));
    s->label = malloc(s->N3 * sizeof(int));
    memset(s->label, -1, s->N3 * sizeof(int));

    return s;
}

static void free_snapshot(Snapshot *s) {
    for (int a = 0; a < NFIELDS; a++) free(s->phi[a]);
    free(s->rho); free(s->absP); free(s->binding_w); free(s->label);
    free(s);
}

/* ================================================================
   Compute derived fields
   ================================================================ */

static void compute_derived(Snapshot *s) {
    const int N = s->N, NN = N*N;
    const double dx = s->dx;
    double P_max = 0;

    #pragma omp parallel for reduction(max:P_max) schedule(static)
    for (long idx = 0; idx < s->N3; idx++) {
        int i=(int)(idx/NN), j=(int)((idx/N)%N), k=(int)(idx%N);

        /* Energy density (skip gradients for boundary cells) */
        double e = 0;
        for (int a = 0; a < NFIELDS; a++)
            e += 0.5 * MASS2 * s->phi[a][idx] * s->phi[a][idx];

        if (i > 0 && i < N-1 && j > 0 && j < N-1 && k > 0 && k < N-1) {
            for (int a = 0; a < NFIELDS; a++) {
                double gx = (s->phi[a][(long)(i+1)*NN+j*N+k] -
                             s->phi[a][(long)(i-1)*NN+j*N+k]) / (2*dx);
                double gy = (s->phi[a][(long)i*NN+(j+1)*N+k] -
                             s->phi[a][(long)i*NN+(j-1)*N+k]) / (2*dx);
                double gz = (s->phi[a][(long)i*NN+j*N+(k+1)] -
                             s->phi[a][(long)i*NN+j*N+(k-1)]) / (2*dx);
                e += 0.5 * (gx*gx + gy*gy + gz*gz);
            }
        }

        double P = s->phi[0][idx] * s->phi[1][idx] * s->phi[2][idx];
        e += (MU/2.0) * P*P / (1.0 + KAPPA*P*P);

        s->rho[idx] = e;
        s->absP[idx] = fabs(P);
        if (fabs(P) > P_max) P_max = fabs(P);
    }

    /* Binding weight */
    double P_thresh = P_max * 0.1;
    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < s->N3; idx++)
        s->binding_w[idx] = 1.0 / (1.0 + s->absP[idx] / (P_thresh + 1e-30));
}

/* ================================================================
   Braid detection: flood-fill connected components above threshold
   ================================================================ */

typedef struct {
    int id;
    /* Position */
    double cx, cy, cz;     /* centroid */
    /* Size */
    double Rxx, Ryy, Rzz;  /* second moments */
    double R_rms;           /* RMS radius */
    int n_cells;            /* number of cells */
    /* Energy/mass */
    double E_total;         /* integrated energy */
    double E_kin_proxy;     /* integrated |phi|² (proportional to KE) */
    double peak_rho;        /* max energy density */
    /* Triple product */
    double P_max;           /* peak |P| */
    double P_avg;           /* average |P| in braid */
    /* Binding structure */
    double w_avg;           /* average binding weakness */
    double w_min;           /* min w (most bound point) */
    /* Phase/winding */
    double winding_z;       /* phase winding along z through centroid */
    /* Momentum proxy */
    double px, py, pz;      /* Σ phi · ∇phi (momentum direction) */
} BraidInfo;

static int detect_braids(Snapshot *s, BraidInfo *braids) {
    const long N3 = s->N3;
    const int N = s->N, NN = N*N;
    const double dx = s->dx, dV = dx*dx*dx;

    /* Threshold: cells with rho > 5× average are "braid" */
    double avg_rho = 0;
    for (long idx = 0; idx < N3; idx++) avg_rho += s->rho[idx];
    avg_rho /= N3;
    double thresh = 5.0 * avg_rho;

    /* Flood-fill to find connected components */
    int n_braids = 0;
    long *stack = malloc(N3 * sizeof(long));

    for (long seed = 0; seed < N3; seed++) {
        if (s->rho[seed] < thresh || s->label[seed] >= 0) continue;
        if (n_braids >= MAX_BRAIDS) break;

        int bid = n_braids++;
        memset(&braids[bid], 0, sizeof(BraidInfo));
        braids[bid].id = bid;
        braids[bid].w_min = 1.0;

        /* BFS flood fill */
        int sp = 0;
        stack[sp++] = seed;
        s->label[seed] = bid;

        while (sp > 0) {
            long idx = stack[--sp];
            int i = (int)(idx/NN), j = (int)((idx/N)%N), k = (int)(idx%N);
            double L = s->L;
            double x = -L + i*dx, y = -L + j*dx, z = -L + k*dx;

            double rho = s->rho[idx];
            double p2 = 0;
            for (int a = 0; a < NFIELDS; a++)
                p2 += s->phi[a][idx] * s->phi[a][idx];

            /* Accumulate properties */
            braids[bid].cx += x * rho * dV;
            braids[bid].cy += y * rho * dV;
            braids[bid].cz += z * rho * dV;
            braids[bid].E_total += rho * dV;
            braids[bid].E_kin_proxy += p2 * dV;
            braids[bid].n_cells++;
            if (rho > braids[bid].peak_rho) braids[bid].peak_rho = rho;
            if (s->absP[idx] > braids[bid].P_max) braids[bid].P_max = s->absP[idx];
            braids[bid].P_avg += s->absP[idx];
            braids[bid].w_avg += s->binding_w[idx];
            if (s->binding_w[idx] < braids[bid].w_min)
                braids[bid].w_min = s->binding_w[idx];

            /* Expand to 6 neighbors */
            int di[] = {1,-1,0,0,0,0};
            int dj[] = {0,0,1,-1,0,0};
            int dk[] = {0,0,0,0,1,-1};
            for (int d = 0; d < 6; d++) {
                int ni = i+di[d], nj = j+dj[d], nk = k+dk[d];
                if (ni<0||ni>=N||nj<0||nj>=N||nk<0||nk>=N) continue;
                long nidx = (long)ni*NN + nj*N + nk;
                if (s->label[nidx] >= 0 || s->rho[nidx] < thresh) continue;
                s->label[nidx] = bid;
                stack[sp++] = nidx;
            }
        }

        /* Finalize centroid */
        if (braids[bid].E_total > 0) {
            braids[bid].cx /= braids[bid].E_total;
            braids[bid].cy /= braids[bid].E_total;
            braids[bid].cz /= braids[bid].E_total;
        }
        if (braids[bid].n_cells > 0) {
            braids[bid].P_avg /= braids[bid].n_cells;
            braids[bid].w_avg /= braids[bid].n_cells;
        }
    }

    /* Second pass: compute second moments and momentum */
    for (long idx = 0; idx < N3; idx++) {
        int bid = s->label[idx];
        if (bid < 0) continue;
        int i = (int)(idx/NN), j = (int)((idx/N)%N), k = (int)(idx%N);
        double L = s->L;
        double x = -L + i*dx - braids[bid].cx;
        double y = -L + j*dx - braids[bid].cy;
        double z = -L + k*dx - braids[bid].cz;
        double rho = s->rho[idx];

        braids[bid].Rxx += x*x * rho * dV;
        braids[bid].Ryy += y*y * rho * dV;
        braids[bid].Rzz += z*z * rho * dV;

        /* Momentum proxy: Σ phi_a × ∂phi_a/∂x_i */
        if (i>0 && i<N-1 && j>0 && j<N-1 && k>0 && k<N-1) {
            for (int a = 0; a < NFIELDS; a++) {
                double dpx = (s->phi[a][(long)(i+1)*NN+j*N+k]-s->phi[a][(long)(i-1)*NN+j*N+k])/(2*dx);
                double dpy = (s->phi[a][(long)i*NN+(j+1)*N+k]-s->phi[a][(long)i*NN+(j-1)*N+k])/(2*dx);
                double dpz = (s->phi[a][(long)i*NN+j*N+(k+1)]-s->phi[a][(long)i*NN+j*N+(k-1)])/(2*dx);
                braids[bid].px += s->phi[a][idx] * dpx * dV;
                braids[bid].py += s->phi[a][idx] * dpy * dV;
                braids[bid].pz += s->phi[a][idx] * dpz * dV;
            }
        }
    }

    /* Finalize sizes */
    for (int b = 0; b < n_braids; b++) {
        if (braids[b].E_total > 0) {
            braids[b].Rxx /= braids[b].E_total;
            braids[b].Ryy /= braids[b].E_total;
            braids[b].Rzz /= braids[b].E_total;
        }
        braids[b].R_rms = sqrt(braids[b].Rxx + braids[b].Ryy + braids[b].Rzz);
    }

    /* Compute winding number for each braid */
    for (int b = 0; b < n_braids; b++) {
        /* Find nearest grid point to centroid */
        int ic = (int)((braids[b].cx + s->L) / dx + 0.5);
        int jc = (int)((braids[b].cy + s->L) / dx + 0.5);
        if (ic < 0) ic = 0; if (ic >= N) ic = N-1;
        if (jc < 0) jc = 0; if (jc >= N) jc = N-1;

        double wind = 0;
        for (int kk = 0; kk < N-1; kk++) {
            long idx0 = (long)ic*NN + jc*N + kk;
            long idx1 = (long)ic*NN + jc*N + kk+1;
            double re0=s->phi[0][idx0], im0=s->phi[1][idx0];
            double re1=s->phi[0][idx1], im1=s->phi[1][idx1];
            wind += atan2(im1*re0 - re1*im0, re1*re0 + im1*im0);
        }
        braids[b].winding_z = wind / (2*PI);
    }

    free(stack);
    return n_braids;
}

/* ================================================================
   Radial profile around a braid
   ================================================================ */

static void radial_profile(Snapshot *s, BraidInfo *braid,
                           double *rho_bins, double *P_bins, double *w_bins,
                           int *counts, int nbins, double dr) {
    const int N = s->N, NN = N*N;
    const double dx = s->dx, L = s->L;

    for (int b = 0; b < nbins; b++) {
        rho_bins[b] = 0; P_bins[b] = 0; w_bins[b] = 0; counts[b] = 0;
    }

    for (int i = 1; i < N-1; i++) {
        double x = -L + i*dx - braid->cx;
        for (int j = 1; j < N-1; j++) {
            double y = -L + j*dx - braid->cy;
            double rp = sqrt(x*x + y*y);
            int b = (int)(rp / dr);
            if (b >= nbins) continue;
            for (int kk = 0; kk < N; kk++) {
                long idx = (long)i*NN + j*N + kk;
                rho_bins[b] += s->rho[idx];
                P_bins[b] += s->absP[idx];
                w_bins[b] += s->binding_w[idx];
                counts[b]++;
            }
        }
    }
    for (int b = 0; b < nbins; b++) {
        if (counts[b] > 0) {
            rho_bins[b] /= counts[b];
            P_bins[b] /= counts[b];
            w_bins[b] /= counts[b];
        }
    }
}

/* ================================================================
   Background field statistics
   ================================================================ */

typedef struct {
    double rho_mean, rho_std;
    double phi2_mean;
    double P_mean;
    long n_bg_cells;
} BackgroundStats;

static BackgroundStats compute_background(Snapshot *s) {
    BackgroundStats bg = {0};
    const long N3 = s->N3;

    /* Background = cells NOT labeled as any braid */
    double sum_rho = 0, sum_rho2 = 0, sum_phi2 = 0, sum_P = 0;
    long n = 0;
    for (long idx = 0; idx < N3; idx++) {
        if (s->label[idx] >= 0) continue;  /* skip braid cells */
        sum_rho += s->rho[idx];
        sum_rho2 += s->rho[idx] * s->rho[idx];
        for (int a = 0; a < NFIELDS; a++)
            sum_phi2 += s->phi[a][idx] * s->phi[a][idx];
        sum_P += s->absP[idx];
        n++;
    }
    if (n > 0) {
        bg.rho_mean = sum_rho / n;
        bg.rho_std = sqrt(sum_rho2/n - bg.rho_mean*bg.rho_mean);
        bg.phi2_mean = sum_phi2 / n;
        bg.P_mean = sum_P / n;
    }
    bg.n_bg_cells = n;
    return bg;
}

/* ================================================================
   Main: analyze one or more snapshots
   ================================================================ */

int main(int argc, char **argv) {
    if (argc < 2) {
        printf("Usage: %s field_t*.bin [field_t*.bin ...]\n", argv[0]);
        printf("Analyzes braid structure from field snapshots.\n");
        return 1;
    }

    omp_set_num_threads(16);
    BraidInfo braids[MAX_BRAIDS];

    /* Process each snapshot */
    for (int fi = 1; fi < argc; fi++) {
        printf("=== Analyzing: %s ===\n", argv[fi]);
        Snapshot *s = load_snapshot(argv[fi]);
        if (!s) continue;

        printf("  N=%d, L=%.1f, t=%.1f, dx=%.4f\n", s->N, s->L, s->t, s->dx);

        /* Compute derived fields */
        compute_derived(s);

        /* Detect braids */
        int nb = detect_braids(s, braids);
        printf("  Detected %d braid(s)\n\n", nb);

        /* Print braid properties */
        printf("  ID  cells   x       y       z       E_total  R_rms  peak_rho  |P|_max  |P|_avg  w_min  w_avg  winding  px      py      pz\n");
        printf("  --- ------  ------  ------  ------  -------  -----  --------  -------  -------  -----  -----  -------  ------  ------  ------\n");
        for (int b = 0; b < nb; b++) {
            BraidInfo *br = &braids[b];
            printf("  %2d  %6d  %+6.2f  %+6.2f  %+6.2f  %7.1f  %5.2f  %8.4f  %7.4f  %7.4f  %5.3f  %5.3f  %+6.3f  %+6.2f  %+6.2f  %+6.2f\n",
                   br->id, br->n_cells, br->cx, br->cy, br->cz,
                   br->E_total, br->R_rms, br->peak_rho,
                   br->P_max, br->P_avg, br->w_min, br->w_avg,
                   br->winding_z, br->px, br->py, br->pz);
        }

        /* Background stats */
        BackgroundStats bg = compute_background(s);
        printf("\n  Background: ρ_mean=%.4e ± %.4e, |φ|²_mean=%.4e, |P|_mean=%.4e, cells=%ld\n",
               bg.rho_mean, bg.rho_std, bg.phi2_mean, bg.P_mean, bg.n_bg_cells);

        /* Radial profiles for each braid */
        int nbins = 80;
        double dr = s->L / nbins;
        double *rho_bins = calloc(nbins, sizeof(double));
        double *P_bins = calloc(nbins, sizeof(double));
        double *w_bins = calloc(nbins, sizeof(double));
        int *counts = calloc(nbins, sizeof(int));

        /* Write per-braid profiles */
        for (int b = 0; b < nb && b < 5; b++) {
            radial_profile(s, &braids[b], rho_bins, P_bins, w_bins, counts, nbins, dr);

            char profname[512];
            snprintf(profname, sizeof(profname), "%s.braid%d_profile.tsv", argv[fi], b);
            FILE *fp = fopen(profname, "w");
            fprintf(fp, "r\trho\t|P|\tw\tcounts\n");
            for (int bn = 0; bn < nbins; bn++) {
                if (counts[bn] > 0)
                    fprintf(fp, "%.3f\t%.6e\t%.6e\t%.6f\t%d\n",
                            (bn+0.5)*dr, rho_bins[bn], P_bins[bn], w_bins[bn], counts[bn]);
            }
            fclose(fp);
            printf("  Wrote: %s\n", profname);
        }

        /* Write summary TSV */
        {
            char sumname[512];
            snprintf(sumname, sizeof(sumname), "%s.braids.tsv", argv[fi]);
            FILE *fp = fopen(sumname, "w");
            fprintf(fp, "id\tn_cells\tx\ty\tz\tE\tR_rms\tpeak_rho\tP_max\tP_avg\tw_min\tw_avg\twinding\tpx\tpy\tpz\n");
            for (int b = 0; b < nb; b++) {
                BraidInfo *br = &braids[b];
                fprintf(fp, "%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.6e\t%.6e\t%.6e\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",
                        br->id, br->n_cells, br->cx, br->cy, br->cz,
                        br->E_total, br->R_rms, br->peak_rho,
                        br->P_max, br->P_avg, br->w_min, br->w_avg,
                        br->winding_z, br->px, br->py, br->pz);
            }
            fclose(fp);
            printf("  Wrote: %s\n", sumname);
        }

        /* Inter-braid properties */
        if (nb >= 2) {
            printf("\n  Inter-braid distances:\n");
            for (int i = 0; i < nb; i++) {
                for (int j = i+1; j < nb; j++) {
                    double dx_ = braids[i].cx - braids[j].cx;
                    double dy_ = braids[i].cy - braids[j].cy;
                    double dz_ = braids[i].cz - braids[j].cz;
                    double dist = sqrt(dx_*dx_ + dy_*dy_ + dz_*dz_);
                    printf("    Braid %d — Braid %d: D=%.2f\n", i, j, dist);
                }
            }
        }

        printf("\n");
        free(rho_bins); free(P_bins); free(w_bins); free(counts);
        free_snapshot(s);
    }

    return 0;
}
