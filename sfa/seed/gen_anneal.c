/*  gen_anneal.c — Annealing seed: tiny high-amplitude chiral spots
 *
 *  Places small clusters of high-amplitude φ and θ spots on a background.
 *  Each spot is a few voxels wide with chirally-oriented θ.
 *  Spots are grouped in clusters of 3 (UUD or UDD chirality).
 *
 *  The idea: let the equations anneal these hot spots into particles
 *  via absorbing BC cooling, rather than constructing particles analytically.
 *
 *  Build: gcc -O3 -march=native -fopenmp -o gen_anneal gen_anneal.c -lzstd -lm
 *
 *  Usage:
 *    ./gen_anneal -N 128 -L 20 -o anneal_seed.sfa
 *    ./gen_anneal -N 384 -L 50 -n_clusters 6 -spot_A 2.0 -spot_R 1.0 -o big.sfa
 */

#define SFA_IMPLEMENTATION
#include "../format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define PI 3.14159265358979323846
#define NCOLS 12

/* A single chiral spot: high amplitude φ and θ at a point */
typedef struct {
    double x, y, z;     /* center position */
    double A;           /* amplitude */
    double R;           /* radius (Gaussian σ) */
    int chi;            /* chirality: +1 (Up) or -1 (Down) */
    int component;      /* which φ component (0,1,2) this spot dominates */
    double phase;       /* carrier phase offset: 0, 2π/3, or 4π/3 */
} Spot;

/* A cluster of 3 spots = one baryon seed */
typedef struct {
    double cx, cy, cz;  /* cluster center */
    int chirality[3];    /* per-component chirality: UUD=(+1,+1,-1), UDD=(+1,-1,-1) */
} Cluster;

int main(int argc, char **argv) {
    int N = 128;
    double L = 20.0;
    double A_bg = 0.1;
    double spot_A = 2.0;    /* spot amplitude (× A_bg) */
    double spot_R = 0.8;    /* spot radius in code units */
    double cluster_R = 2.5; /* separation between spots within a cluster */
    int n_clusters = 6;     /* total clusters (half UUD, half UDD) */
    double k_wave = 1.5;
    double delta[3] = {0.0, 3.0005, 4.4325};
    double Delta[3] = {0.0, 2*PI/3, 4*PI/3};
    double alpha_p = 0.1, eta_p = 0.5, beta_p = 0.5;
    char outpath[512] = "anneal_seed.sfa";
    int precision = 1; /* f32 */
    unsigned int seed = 42;

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-N") && i+1<argc) N = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-L") && i+1<argc) L = atof(argv[++i]);
        else if (!strcmp(argv[i], "-o") && i+1<argc) strncpy(outpath, argv[++i], 511);
        else if (!strcmp(argv[i], "-spot_A") && i+1<argc) spot_A = atof(argv[++i]);
        else if (!strcmp(argv[i], "-spot_R") && i+1<argc) spot_R = atof(argv[++i]);
        else if (!strcmp(argv[i], "-cluster_R") && i+1<argc) cluster_R = atof(argv[++i]);
        else if (!strcmp(argv[i], "-n_clusters") && i+1<argc) n_clusters = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-A_bg") && i+1<argc) A_bg = atof(argv[++i]);
        else if (!strcmp(argv[i], "-seed") && i+1<argc) seed = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-precision") && i+1<argc) {
            i++; precision = !strcmp(argv[i], "f16") ? 0 : 1;
        }
        else if (!strcmp(argv[i], "-h")) {
            fprintf(stderr, "Usage: %s [-N grid] [-L domain] [-spot_A amp] [-spot_R radius] "
                    "[-cluster_R sep] [-n_clusters N] [-seed N] [-o output.sfa]\n", argv[0]);
            return 0;
        }
    }

    srand(seed);
    double dx = 2.0 * L / N;
    double dt = 0.025 * dx;
    long N3 = (long)N * N * N;
    int NN = N * N;

    printf("gen_anneal: N=%d L=%.1f dx=%.4f\n", N, L, dx);
    printf("  spot_A=%.1f (×A_bg=%.1f → peak=%.2f) spot_R=%.2f cluster_R=%.2f\n",
           spot_A, A_bg, spot_A * A_bg, spot_R, cluster_R);
    printf("  n_clusters=%d (half UUD, half UDD), seed=%u\n", n_clusters, seed);

    /* Generate cluster positions — spread across the central region */
    double place_R = L * 0.5; /* place within 50% of domain */
    Cluster *clusters = calloc(n_clusters, sizeof(Cluster));

    for (int c = 0; c < n_clusters; c++) {
        /* Place in a rough grid pattern to avoid overlaps */
        double angle = 2 * PI * c / n_clusters;
        double radial = place_R * (0.3 + 0.5 * ((double)rand() / RAND_MAX));
        clusters[c].cx = radial * cos(angle) * (0.7 + 0.3 * ((double)rand() / RAND_MAX));
        clusters[c].cy = radial * sin(angle) * (0.7 + 0.3 * ((double)rand() / RAND_MAX));
        clusters[c].cz = place_R * (((double)rand() / RAND_MAX) - 0.5) * 0.6;

        /* Alternate UUD and UDD */
        if (c % 2 == 0) {
            clusters[c].chirality[0] = +1;
            clusters[c].chirality[1] = +1;
            clusters[c].chirality[2] = -1;
        } else {
            clusters[c].chirality[0] = +1;
            clusters[c].chirality[1] = -1;
            clusters[c].chirality[2] = -1;
        }
    }

    /* Generate spots: 3 per cluster, one per field component */
    int n_spots = n_clusters * 3;
    Spot *spots = calloc(n_spots, sizeof(Spot));

    for (int c = 0; c < n_clusters; c++) {
        for (int s = 0; s < 3; s++) {
            Spot *sp = &spots[c * 3 + s];
            /* Place each spot along a different axis relative to cluster center */
            double offset = cluster_R * (0.5 + 0.5 * ((double)rand() / RAND_MAX));
            double ox = 0, oy = 0, oz = 0;
            switch (s) {
                case 0: ox = offset; break;
                case 1: oy = offset; break;
                case 2: oz = offset; break;
            }
            /* Add some randomness */
            ox += cluster_R * 0.2 * (((double)rand() / RAND_MAX) - 0.5);
            oy += cluster_R * 0.2 * (((double)rand() / RAND_MAX) - 0.5);
            oz += cluster_R * 0.2 * (((double)rand() / RAND_MAX) - 0.5);

            sp->x = clusters[c].cx + ox;
            sp->y = clusters[c].cy + oy;
            sp->z = clusters[c].cz + oz;
            sp->A = spot_A * A_bg;
            sp->R = spot_R;
            sp->chi = clusters[c].chirality[s];
            sp->component = s;
            sp->phase = Delta[s];
        }
    }

    printf("  Generated %d spots in %d clusters\n", n_spots, n_clusters);
    for (int c = 0; c < n_clusters; c++) {
        printf("    Cluster %d: center=(%.1f,%.1f,%.1f) %s\n",
               c, clusters[c].cx, clusters[c].cy, clusters[c].cz,
               clusters[c].chirality[1] == 1 ? "UUD" : "UDD");
    }

    /* Allocate fields */
    double *fields[NCOLS];
    for (int c = 0; c < NCOLS; c++) {
        fields[c] = calloc(N3, sizeof(double));
    }

    /* Step 1: Initialize background carrier wave */
    double omega = sqrt(k_wave * k_wave + k_wave * k_wave); /* not physical, just for vel init */
    for (int iz = 0; iz < N; iz++)
    for (int iy = 0; iy < N; iy++)
    for (int ix = 0; ix < N; ix++) {
        long idx = (long)iz * NN + (long)iy * N + ix;
        double z = -L + (iz + 0.5) * dx;
        for (int a = 0; a < 3; a++) {
            double phase = k_wave * z + delta[a] + Delta[a];
            fields[a][idx] = A_bg * cos(phase);
            fields[6 + a][idx] = A_bg * omega * sin(phase); /* phi_vel */
        }
        /* θ and θ_vel start at zero (background has no curl) */
    }

    /* Step 2: Stamp spots — add localized bumps with chiral θ */
    double P_max = 0, theta_max = 0;

    for (int si = 0; si < n_spots; si++) {
        Spot *sp = &spots[si];
        double inv2sig2 = 1.0 / (2.0 * sp->R * sp->R);

        for (int iz = 0; iz < N; iz++)
        for (int iy = 0; iy < N; iy++)
        for (int ix = 0; ix < N; ix++) {
            long idx = (long)iz * NN + (long)iy * N + ix;
            double x = -L + (ix + 0.5) * dx;
            double y = -L + (iy + 0.5) * dx;
            double z = -L + (iz + 0.5) * dx;

            double rx = x - sp->x;
            double ry = y - sp->y;
            double rz = z - sp->z;
            double r2 = rx*rx + ry*ry + rz*rz;
            double env = exp(-r2 * inv2sig2);

            if (env < 1e-6) continue;

            /* Add amplitude to the dominant field component */
            int a = sp->component;
            double phase = k_wave * z + delta[a] + sp->phase;
            fields[a][idx] += sp->A * env * cos(phase);
            fields[6 + a][idx] += sp->A * env * omega * sin(phase);

            /* Add chirally-oriented θ:
             * θ_a gets a contribution proportional to chi × env
             * This creates alternating θ polarization between spots */
            double theta_amp = sp->A * 0.3 * sp->chi; /* 30% of φ amplitude, signed by chirality */
            /* θ contributes to the perpendicular components (curl-like) */
            double r = sqrt(r2 + 1e-30);
            /* Use a dipole-like pattern: θ ∝ cross(r_hat, axis) */
            double nx = rx/r, ny = ry/r, nz = rz/r;
            switch (sp->component) {
                case 0: /* x-axis spot → θ_y, θ_z */
                    fields[4][idx] += theta_amp * env * nz; /* θ_y ∝ nz */
                    fields[5][idx] -= theta_amp * env * ny; /* θ_z ∝ -ny */
                    break;
                case 1: /* y-axis spot → θ_x, θ_z */
                    fields[3][idx] -= theta_amp * env * nz; /* θ_x ∝ -nz */
                    fields[5][idx] += theta_amp * env * nx; /* θ_z ∝ nx */
                    break;
                case 2: /* z-axis spot → θ_x, θ_y */
                    fields[3][idx] += theta_amp * env * ny; /* θ_x ∝ ny */
                    fields[4][idx] -= theta_amp * env * nx; /* θ_y ∝ -nx */
                    break;
            }
        }
    }

    /* Measure diagnostics */
    for (long idx = 0; idx < N3; idx++) {
        double P = fabs(fields[0][idx] * fields[1][idx] * fields[2][idx]);
        if (P > P_max) P_max = P;
        for (int a = 0; a < 3; a++) {
            double t = fabs(fields[3+a][idx]);
            if (t > theta_max) theta_max = t;
        }
    }

    printf("\nDiagnostics:\n");
    printf("  P_max = %.6e\n", P_max);
    printf("  theta_max = %.6e\n", theta_max);

    /* Write SFA */
    uint8_t sfa_dtype = precision == 0 ? SFA_F16 : SFA_F32;
    SFA *sfa = sfa_create(outpath, N, N, N, L, L, L, dt);

    char buf[8][64];
    snprintf(buf[0], 64, "%d", n_clusters);
    snprintf(buf[1], 64, "%.2f", spot_A);
    snprintf(buf[2], 64, "%.2f", spot_R);
    snprintf(buf[3], 64, "%u", seed);
    const char *keys[] = {"generator", "n_clusters", "spot_A", "spot_R", "random_seed"};
    const char *vals[] = {"gen_anneal", buf[0], buf[1], buf[2], buf[3]};
    sfa_add_kvmd(sfa, 0, 0xFFFFFFFF, 0xFFFFFFFF, keys, vals, 5);

    sfa_add_column(sfa, "phi_x",    sfa_dtype, SFA_POSITION, 0);
    sfa_add_column(sfa, "phi_y",    sfa_dtype, SFA_POSITION, 1);
    sfa_add_column(sfa, "phi_z",    sfa_dtype, SFA_POSITION, 2);
    sfa_add_column(sfa, "theta_x",  sfa_dtype, SFA_ANGLE,    0);
    sfa_add_column(sfa, "theta_y",  sfa_dtype, SFA_ANGLE,    1);
    sfa_add_column(sfa, "theta_z",  sfa_dtype, SFA_ANGLE,    2);
    sfa_add_column(sfa, "phi_vx",   sfa_dtype, SFA_VELOCITY, 0);
    sfa_add_column(sfa, "phi_vy",   sfa_dtype, SFA_VELOCITY, 1);
    sfa_add_column(sfa, "phi_vz",   sfa_dtype, SFA_VELOCITY, 2);
    sfa_add_column(sfa, "theta_vx", sfa_dtype, SFA_VELOCITY, 3);
    sfa_add_column(sfa, "theta_vy", sfa_dtype, SFA_VELOCITY, 4);
    sfa_add_column(sfa, "theta_vz", sfa_dtype, SFA_VELOCITY, 5);
    sfa_finalize_header(sfa);

    if (precision == 1) {
        void *cols[NCOLS];
        for (int c = 0; c < NCOLS; c++) {
            float *p = malloc(N3 * sizeof(float));
            for (long i = 0; i < N3; i++) p[i] = (float)fields[c][i];
            cols[c] = p;
        }
        sfa_write_frame(sfa, 0.0, cols);
        for (int c = 0; c < NCOLS; c++) free(cols[c]);
    } else {
        void *cols[NCOLS];
        for (int c = 0; c < NCOLS; c++) {
            uint16_t *p = malloc(N3 * 2);
            for (long i = 0; i < N3; i++) {
                float f = (float)fields[c][i];
                uint32_t x; memcpy(&x, &f, 4);
                uint16_t sign = (x >> 16) & 0x8000;
                int exp = ((x >> 23) & 0xFF) - 127 + 15;
                uint16_t mant = (x >> 13) & 0x3FF;
                p[i] = (exp <= 0) ? sign : (exp >= 31) ? (sign|0x7C00) : (sign|(exp<<10)|mant);
            }
            cols[c] = p;
        }
        sfa_write_frame(sfa, 0.0, cols);
        for (int c = 0; c < NCOLS; c++) free(cols[c]);
    }

    sfa_close(sfa);
    printf("\nWrote %s (%d^3, %d clusters, %d spots)\n", outpath, N, n_clusters, n_spots);

    for (int c = 0; c < NCOLS; c++) free(fields[c]);
    free(clusters);
    free(spots);
    return 0;
}
