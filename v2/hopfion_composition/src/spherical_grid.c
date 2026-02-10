/*
 * spherical_grid.c — Masked-Cartesian spherical grid implementation
 */

#include "spherical_grid.h"
#include <stdio.h>

SphericalGrid *sg_alloc(int N, double L, double R, double sponge_width)
{
    SphericalGrid *g = calloc(1, sizeof(SphericalGrid));
    if (!g) return NULL;

    g->N = N;
    g->L = L;
    g->h = 2.0 * L / N;
    g->R = R;
    g->R_sponge = R - sponge_width;
    if (g->R_sponge < 0) g->R_sponge = 0;

    long total = (long)N * N * N;

    /* Allocate dense arrays */
    g->psi   = calloc(total, sizeof(Multivector));
    g->vel   = calloc(total, sizeof(Multivector));
    g->force = calloc(total, sizeof(Multivector));
    g->mask  = calloc(total, sizeof(unsigned char));
    g->sponge = calloc(total, sizeof(double));

    if (!g->psi || !g->vel || !g->force || !g->mask || !g->sponge) {
        fprintf(stderr, "sg_alloc: out of memory (N=%d, %.1f GB)\n",
                N, total * (3*sizeof(Multivector) + sizeof(unsigned char) + sizeof(double)) / 1e9);
        sg_free(g);
        return NULL;
    }

    /* Build mask and sponge */
    int n_active = 0;
    for (int i = 0; i < N; i++) {
        double x = -L + (i + 0.5) * g->h;
        for (int j = 0; j < N; j++) {
            double y = -L + (j + 0.5) * g->h;
            for (int k = 0; k < N; k++) {
                double z = -L + (k + 0.5) * g->h;
                double r = sqrt(x*x + y*y + z*z);
                int ix = sg_idx(N, i, j, k);

                if (r <= R) {
                    g->mask[ix] = 1;
                    n_active++;

                    /* Sponge damping: smooth ramp in [R_sponge, R] */
                    if (r > g->R_sponge) {
                        double s = (r - g->R_sponge) / sponge_width;
                        /* Smooth cubic: 3s² - 2s³ */
                        g->sponge[ix] = s * s * (3.0 - 2.0 * s);
                    }
                }
            }
        }
    }

    /* Build active cell list */
    g->active = malloc(n_active * sizeof(int));
    g->n_active = n_active;
    int n = 0;
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        int ix = sg_idx(N, i, j, k);
        if (g->mask[ix])
            g->active[n++] = ix;
    }

    /* Default vacuum: q = (rho0, 0, 0, 0), weight = 0 */
    g->vacuum = mv_zero();
    g->vacuum.s = 1.0;  /* will be set properly by sg_set_vacuum */

    return g;
}

void sg_free(SphericalGrid *g)
{
    if (!g) return;
    free(g->psi);
    free(g->vel);
    free(g->force);
    free(g->mask);
    free(g->sponge);
    free(g->active);
    free(g);
}

void sg_set_vacuum(SphericalGrid *g, double rho0)
{
    g->vacuum = mv_zero();
    g->vacuum.s = rho0;

    long total = (long)g->N * g->N * g->N;
    for (long ix = 0; ix < total; ix++) {
        g->psi[ix] = g->vacuum;
        g->vel[ix] = mv_zero();
        g->force[ix] = mv_zero();
    }
}

int sg_save(const SphericalGrid *g, const char *filename)
{
    FILE *fp = fopen(filename, "wb");
    if (!fp) return -1;

    /* Header */
    fwrite(&g->N, sizeof(int), 1, fp);
    fwrite(&g->L, sizeof(double), 1, fp);
    fwrite(&g->R, sizeof(double), 1, fp);
    fwrite(&g->R_sponge, sizeof(double), 1, fp);
    fwrite(&g->n_active, sizeof(int), 1, fp);

    /* Only write active cells (sparse) */
    for (int n = 0; n < g->n_active; n++) {
        int ix = g->active[n];
        fwrite(&ix, sizeof(int), 1, fp);
        fwrite(&g->psi[ix], sizeof(Multivector), 1, fp);
        fwrite(&g->vel[ix], sizeof(Multivector), 1, fp);
    }

    fclose(fp);
    return 0;
}

int sg_load(SphericalGrid *g, const char *filename)
{
    FILE *fp = fopen(filename, "rb");
    if (!fp) return -1;

    int N_file, n_active_file;
    double L_file, R_file, R_sponge_file;
    fread(&N_file, sizeof(int), 1, fp);
    fread(&L_file, sizeof(double), 1, fp);
    fread(&R_file, sizeof(double), 1, fp);
    fread(&R_sponge_file, sizeof(double), 1, fp);
    fread(&n_active_file, sizeof(int), 1, fp);

    if (N_file != g->N || fabs(L_file - g->L) > 1e-10) {
        fprintf(stderr, "sg_load: grid mismatch (file N=%d L=%.2f, grid N=%d L=%.2f)\n",
                N_file, L_file, g->N, g->L);
        fclose(fp);
        return -1;
    }

    /* Read active cells */
    for (int n = 0; n < n_active_file; n++) {
        int ix;
        Multivector psi_val, vel_val;
        fread(&ix, sizeof(int), 1, fp);
        fread(&psi_val, sizeof(Multivector), 1, fp);
        fread(&vel_val, sizeof(Multivector), 1, fp);
        if (ix >= 0 && ix < g->N * g->N * g->N) {
            g->psi[ix] = psi_val;
            g->vel[ix] = vel_val;
        }
    }

    fclose(fp);
    return 0;
}
