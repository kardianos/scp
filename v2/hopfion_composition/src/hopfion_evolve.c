/*
 * hopfion_evolve.c — Main driver for hopfion composition simulation
 *
 * Leapfrog (Störmer-Verlet) time evolution on a spherical grid.
 * Supports: initialization, gradient flow, time evolution, diagnostics.
 *
 * Usage:
 *   ./bin/hopfion_evolve [options]
 *
 * Options:
 *   -N <int>      Grid points per dimension (default: 128)
 *   -R <float>    Sphere radius (default: 8.0)
 *   -sponge <f>   Sponge layer width (default: 2.0)
 *   -e <float>    Skyrme coupling (default: 1.0)
 *   -lam <float>  Mexican hat coupling (default: 0)
 *   -lam6 <float> Sextic/BPS coupling (default: 0)
 *   -mpi <float>  Pion mass (default: 0)
 *   -rho0 <float> Vacuum density (default: 1.0)
 *   -dt <float>   Timestep (default: auto CFL)
 *   -T <float>    Max time (default: 10.0)
 *   -sigma        Enforce sigma-model |q|=ρ₀
 *   -init <mode>  Initialization: skyrmion, hopfion, hopfion_nm, composed, hybrid
 *   -profile <f>  Radial profile file for skyrmion init
 *   -a <float>    Hopfion size parameter (default: 2.0)
 *   -nm <n> <m>   Hopfion winding numbers (default: 1 1)
 *   -gflow        Run gradient flow (relax to minimum energy)
 *   -out <int>    Output interval in steps (default: 20)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "spherical_grid.h"
#include "hopfion_field.h"
#include "hopfion_init.h"
#include "topology.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static void print_usage(const char *prog)
{
    fprintf(stderr, "Usage: %s [options]\n", prog);
    fprintf(stderr, "  -N <int>       Grid points per dimension [128]\n");
    fprintf(stderr, "  -R <float>     Sphere radius [8.0]\n");
    fprintf(stderr, "  -sponge <f>    Sponge width [2.0]\n");
    fprintf(stderr, "  -e <float>     Skyrme coupling [1.0]\n");
    fprintf(stderr, "  -lam <float>   Mexican hat coupling [0=sigma]\n");
    fprintf(stderr, "  -lam6 <float>  Sextic/BPS coupling [0]\n");
    fprintf(stderr, "  -mpi <float>   Pion mass [0]\n");
    fprintf(stderr, "  -rho0 <float>  Vacuum density [1.0]\n");
    fprintf(stderr, "  -dt <float>    Timestep [auto]\n");
    fprintf(stderr, "  -T <float>     Max time [10.0]\n");
    fprintf(stderr, "  -sigma         Sigma model projection\n");
    fprintf(stderr, "  -init <mode>   skyrmion|hopfion|hopfion_nm|composed|hybrid\n");
    fprintf(stderr, "  -profile <f>   Radial profile for skyrmion\n");
    fprintf(stderr, "  -a <float>     Hopfion size [2.0]\n");
    fprintf(stderr, "  -nm <n> <m>    Hopfion winding [1 1]\n");
    fprintf(stderr, "  -gflow         Gradient flow (no dynamics)\n");
    fprintf(stderr, "  -out <int>     Output interval [20]\n");
}

int main(int argc, char *argv[])
{
    /* Default parameters */
    int N = 128;
    double R = 8.0;
    double sponge_width = 2.0;
    double e_skyrme = 1.0;
    double lambda = 0;
    double lambda6 = 0;
    double m_pi = 0;
    double rho0 = 1.0;
    double dt = 0;
    double T_max = 10.0;
    int sigma_model = 0;
    const char *init_mode = "hopfion";
    const char *profile_file = NULL;
    double hopfion_a = 2.0;
    int n_tor = 1, m_pol = 1;
    int gradient_flow = 0;
    int output_interval = 20;

    /* Parse arguments */
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-N") && i+1 < argc) N = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-R") && i+1 < argc) R = atof(argv[++i]);
        else if (!strcmp(argv[i], "-sponge") && i+1 < argc) sponge_width = atof(argv[++i]);
        else if (!strcmp(argv[i], "-e") && i+1 < argc) e_skyrme = atof(argv[++i]);
        else if (!strcmp(argv[i], "-lam") && i+1 < argc) lambda = atof(argv[++i]);
        else if (!strcmp(argv[i], "-lam6") && i+1 < argc) lambda6 = atof(argv[++i]);
        else if (!strcmp(argv[i], "-mpi") && i+1 < argc) m_pi = atof(argv[++i]);
        else if (!strcmp(argv[i], "-rho0") && i+1 < argc) rho0 = atof(argv[++i]);
        else if (!strcmp(argv[i], "-dt") && i+1 < argc) dt = atof(argv[++i]);
        else if (!strcmp(argv[i], "-T") && i+1 < argc) T_max = atof(argv[++i]);
        else if (!strcmp(argv[i], "-sigma")) sigma_model = 1;
        else if (!strcmp(argv[i], "-init") && i+1 < argc) init_mode = argv[++i];
        else if (!strcmp(argv[i], "-profile") && i+1 < argc) profile_file = argv[++i];
        else if (!strcmp(argv[i], "-a") && i+1 < argc) hopfion_a = atof(argv[++i]);
        else if (!strcmp(argv[i], "-nm") && i+2 < argc) {
            n_tor = atoi(argv[++i]); m_pol = atoi(argv[++i]);
        }
        else if (!strcmp(argv[i], "-gflow")) gradient_flow = 1;
        else if (!strcmp(argv[i], "-out") && i+1 < argc) output_interval = atoi(argv[++i]);
        else { print_usage(argv[0]); return 1; }
    }

    double L = R + 0.5;  /* small buffer beyond sphere */
    double h = 2.0 * L / N;

    /* Auto-set timestep from CFL */
    if (dt <= 0) {
        dt = 0.2 * h;  /* CFL factor 0.2 (conservative for L₂) */
        /* If Skyrme term active, need smaller dt */
        if (e_skyrme > 0 && e_skyrme < 100) {
            double dt_skyrme = 0.1 * h * h;  /* L₄ scales as h² */
            if (dt_skyrme < dt) dt = dt_skyrme;
        }
    }

    printf("==========================================================\n");
    printf(" Hopfion Composition Simulator\n");
    printf("==========================================================\n\n");
    printf("Grid: N=%d, R=%.1f, L=%.1f, h=%.4f\n", N, R, L, h);
    printf("Sponge: width=%.1f, R_inner=%.1f\n", sponge_width, R-sponge_width);
    printf("Parameters: rho0=%.1f, e=%.1f, lam=%.0f, lam6=%.2f, m_pi=%.3f\n",
           rho0, e_skyrme, lambda, lambda6, m_pi);
    printf("Init: %s", init_mode);
    if (!strcmp(init_mode, "hopfion") || !strcmp(init_mode, "hopfion_nm"))
        printf(" (a=%.2f, n=%d, m=%d)", hopfion_a, n_tor, m_pol);
    printf("\n");
    printf("Mode: %s, dt=%.6f, T=%.1f\n",
           gradient_flow ? "gradient flow" : (sigma_model ? "sigma-model dynamics" : "full dynamics"),
           dt, T_max);
    printf("Threads: %d\n\n", omp_get_max_threads());

    /* Allocate grid */
    SphericalGrid *g = sg_alloc(N, L, R, sponge_width);
    if (!g) { fprintf(stderr, "Grid allocation failed\n"); return 1; }

    printf("Grid: %d active cells (%.1f%% of cube, %.1f MB)\n",
           g->n_active, 100.0*g->n_active/((long)N*N*N),
           g->n_active * sizeof(Multivector) * 3.0 / 1e6);

    /* Set vacuum and initialize */
    sg_set_vacuum(g, rho0);

    HopfionParams params = {
        .rho0 = rho0,
        .e_skyrme = e_skyrme,
        .lambda = lambda,
        .lambda6 = lambda6,
        .m_pi_sq = m_pi * m_pi,
        .mu = 0,
        .c = 1.0
    };

    RadialProfile *prof = NULL;
    if (profile_file) {
        prof = profile_load(profile_file);
        if (!prof) { fprintf(stderr, "Failed to load profile\n"); sg_free(g); return 1; }
    }

    /* Initialize field */
    if (!strcmp(init_mode, "skyrmion")) {
        if (!prof) { fprintf(stderr, "Skyrmion init requires -profile\n"); sg_free(g); return 1; }
        init_skyrmion(g, prof, rho0, 0, 0, 0);
    } else if (!strcmp(init_mode, "hopfion")) {
        init_hopfion(g, rho0, hopfion_a, 0, 0, 0);
    } else if (!strcmp(init_mode, "hopfion_nm")) {
        init_hopfion_nm(g, rho0, hopfion_a, n_tor, m_pol, 0, 0, 0);
    } else if (!strcmp(init_mode, "composed")) {
        /* Two hopfions at offset positions */
        printf("Composing two hopfions...\n");
        init_hopfion(g, rho0, hopfion_a, 0, 0, +2.0);

        /* Save first hopfion, init second, compose */
        long total = (long)N * N * N;
        Multivector *h1 = malloc(total * sizeof(Multivector));
        memcpy(h1, g->psi, total * sizeof(Multivector));

        init_hopfion(g, rho0, hopfion_a, 0, 0, -2.0);

        /* Product ansatz */
        double inv_rho0 = 1.0 / rho0;
        #pragma omp parallel for schedule(static)
        for (int n = 0; n < g->n_active; n++) {
            int ix = g->active[n];
            Multivector a = h1[ix];
            Multivector b = g->psi[ix];
            g->psi[ix].s  = inv_rho0 * (a.s*b.s  - a.f1*b.f1 - a.f2*b.f2 - a.f3*b.f3);
            g->psi[ix].f1 = inv_rho0 * (a.s*b.f1 + a.f1*b.s  - a.f2*b.f3 + a.f3*b.f2);
            g->psi[ix].f2 = inv_rho0 * (a.s*b.f2 + a.f1*b.f3 + a.f2*b.s  - a.f3*b.f1);
            g->psi[ix].f3 = inv_rho0 * (a.s*b.f3 - a.f1*b.f2 + a.f2*b.f1 + a.f3*b.s);
        }
        free(h1);
    } else if (!strcmp(init_mode, "hybrid")) {
        if (!prof) { fprintf(stderr, "Hybrid init requires -profile\n"); sg_free(g); return 1; }
        init_skyrmion_hopfion(g, prof, rho0, hopfion_a, 3.0);
    } else {
        fprintf(stderr, "Unknown init mode: %s\n", init_mode);
        sg_free(g); return 1;
    }

    if (sigma_model) {
        hf_sigma_project(g, rho0);
    }

    /* Initial diagnostics */
    HopfionEnergy en = hf_energy(g, &params);
    double B = hf_baryon_charge(g);
    printf("\n--- Initial state ---\n");
    printf("E_pot    = %.6f  (E2=%.4f, E4=%.4f, E6=%.4f, EV=%.4f, Epi=%.4f)\n",
           en.Etotal, en.E2, en.E4, en.E6, en.EV, en.Epi);
    printf("B        = %.4f  (Skyrmion charge)\n", B);

    /* Compute Hopf charge (expensive — only at init and end) */
    printf("Computing Hopf charge (this may take a moment)...\n");
    double H = topo_hopf_charge(g);
    printf("H        = %.4f  (Hopf charge)\n", H);

    /* Time evolution */
    printf("\n--- %s ---\n",
           gradient_flow ? "Gradient flow" : "Time evolution");
    printf("%-6s %-8s %-12s %-12s %-12s %-8s\n",
           "step", "t", "E_pot", "E_kin", "E_tot", "B");

    double E_init = en.Etotal;
    int max_steps = (int)(T_max / dt + 0.5);

    for (int step = 1; step <= max_steps; step++) {
        double t = step * dt;

        /* Compute force */
        hf_force(g, &params);

        if (gradient_flow) {
            /* Gradient flow: dq/dtau = F (damped, no inertia) */
            double gf_dt = 0.01 * h * h;  /* much smaller for stability */

            #pragma omp parallel for schedule(static)
            for (int n = 0; n < g->n_active; n++) {
                int ix = g->active[n];
                g->psi[ix] = mv_add(g->psi[ix], mv_scale(gf_dt, g->force[ix]));
            }

            if (sigma_model) hf_sigma_project(g, rho0);
        } else {
            /* Leapfrog: v(t+dt/2) = v(t-dt/2) + dt*c²*F
             *           q(t+dt)   = q(t) + dt*v(t+dt/2) */
            double c2_dt = params.c * params.c * dt;

            #pragma omp parallel for schedule(static)
            for (int n = 0; n < g->n_active; n++) {
                int ix = g->active[n];

                /* Sponge damping */
                double damp = 1.0 - g->sponge[ix];

                g->vel[ix] = mv_add(mv_scale(damp, g->vel[ix]),
                                    mv_scale(c2_dt, g->force[ix]));
                g->psi[ix] = mv_add(g->psi[ix], mv_scale(dt, g->vel[ix]));
            }

            if (sigma_model) hf_sigma_project(g, rho0);
        }

        /* Diagnostics */
        if (step % output_interval == 0) {
            en = hf_energy(g, &params);
            double Ek = gradient_flow ? 0 : hf_kinetic_energy(g, &params);
            B = hf_baryon_charge(g);
            double Etot = en.Etotal + Ek;
            double dE = (Etot - E_init) / (fabs(E_init) + 1e-20);

            printf("%6d %8.3f %12.4f %12.4f %12.4f %8.4f",
                   step, t, en.Etotal, Ek, Etot, B);
            if (fabs(dE) > 0.01)
                printf("  dE=%.2e!", dE);
            printf("\n");

            /* Energy conservation check */
            if (!gradient_flow && fabs(dE) > 0.5) {
                printf("*** Energy conservation violated by >50%% — aborting ***\n");
                break;
            }
        }
    }

    /* Final diagnostics */
    en = hf_energy(g, &params);
    B = hf_baryon_charge(g);
    printf("\n--- Final state ---\n");
    printf("E_pot    = %.6f  (E2=%.4f, E4=%.4f, E6=%.4f, EV=%.4f)\n",
           en.Etotal, en.E2, en.E4, en.E6, en.EV);
    printf("B        = %.4f\n", B);

    printf("Computing final Hopf charge...\n");
    H = topo_hopf_charge(g);
    printf("H        = %.4f\n", H);

    /* Cleanup */
    if (prof) profile_free(prof);
    sg_free(g);

    return 0;
}
