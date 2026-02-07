/*
 * main.c — CHPT Soliton Search via Gradient Flow
 *
 * Finds minimum-energy static configurations of the CHPT field equation
 * by gradient descent (energy minimization) starting from a hedgehog ansatz.
 *
 * The full 8-component Cl+(3,0,1) field is evolved:
 *   dPsi/dtau = -dE/dPsi  (gradient flow in fictitious time tau)
 *
 * Output: binary field snapshots and ASCII energy log.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include "field.h"
#include "initial.h"

/* ---------- I/O ---------- */

static void write_energy_log(FILE *fp, int step, double tau, Energy en, double Q, double max_force)
{
    fprintf(fp, "%6d  %12.6e  %14.8e  %14.8e  %14.8e  %14.8e  %14.8e  %10.6f  %12.6e\n",
            step, tau, en.Etotal, en.E2, en.E4, en.EV, en.ED, Q, max_force);
    fflush(fp);
}

static void write_field_binary(const char *filename, const Field *f)
{
    FILE *fp = fopen(filename, "wb");
    if (!fp) { fprintf(stderr, "Cannot open %s\n", filename); return; }
    int N = f->N;
    double L = f->L;
    fwrite(&N, sizeof(int), 1, fp);
    fwrite(&L, sizeof(double), 1, fp);
    fwrite(f->psi, sizeof(Multivector), (size_t)N*N*N, fp);
    fclose(fp);
}

/* ---------- Gradient flow with adaptive step size ---------- */

static double max_force_norm(const Multivector *force, int N3)
{
    double maxf = 0;
    #pragma omp parallel for reduction(max:maxf) schedule(static)
    for (int i = 0; i < N3; i++) {
        double f2 = mv_dot(force[i], force[i]);
        if (f2 > maxf) maxf = f2;
    }
    return sqrt(maxf);
}

static void gradient_flow_step(Field *f, const Multivector *force, double dt)
{
    int N3 = f->N * f->N * f->N;
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < N3; i++) {
        f->psi[i] = mv_add(f->psi[i], mv_scale(dt, force[i]));
    }
}

/* Project quaternion part onto sigma model constraint |q| = rho0.
 * This enforces the strong-coupling limit (lambda -> infinity) and
 * protects topological charge during gradient flow. */
static void project_sigma_model(Field *f, double rho0)
{
    int N3 = f->N * f->N * f->N;
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < N3; i++) {
        double s = f->psi[i].s, b1 = f->psi[i].f1;
        double b2 = f->psi[i].f2, b3 = f->psi[i].f3;
        double r = sqrt(s*s + b1*b1 + b2*b2 + b3*b3);
        if (r > 1e-30) {
            double scale = rho0 / r;
            f->psi[i].s  *= scale;
            f->psi[i].f1 *= scale;
            f->psi[i].f2 *= scale;
            f->psi[i].f3 *= scale;
        }
    }
}

/* ---------- Main ---------- */

int main(int argc, char **argv)
{
    /* Default parameters */
    int N = 64;                /* grid points per dimension */
    double L = 8.0;           /* half-width of domain */
    double R_init = 1.5;      /* initial soliton size */
    int max_steps = 50000;    /* maximum gradient flow steps */
    double dt_init = 1e-4;    /* initial step size */
    double dt_max = 1e-2;     /* maximum step size */
    double tol = 1e-8;        /* convergence tolerance (relative energy change) */
    int log_interval = 100;   /* steps between log entries */
    int snap_interval = 5000; /* steps between field snapshots */
    double perturb_amp = 0.0; /* perturbation amplitude for degenerate sector */
    int sigma_model = 1;      /* project onto sigma model constraint (default: on) */

    /* Physical parameters (dimensionless units) */
    Params params;
    params.rho0 = 1.0;
    params.lambda = 100.0;    /* Strong potential keeps |q| ~ rho0, protects topology */
    params.e_skyrme = 4.0;    /* Standard Skyrme value ~4 */
    params.mu = 1.0;
    params.c = 1.0;

    /* Parse command line */
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-N") == 0 && i+1 < argc) N = atoi(argv[++i]);
        else if (strcmp(argv[i], "-L") == 0 && i+1 < argc) L = atof(argv[++i]);
        else if (strcmp(argv[i], "-R") == 0 && i+1 < argc) R_init = atof(argv[++i]);
        else if (strcmp(argv[i], "-steps") == 0 && i+1 < argc) max_steps = atoi(argv[++i]);
        else if (strcmp(argv[i], "-dt") == 0 && i+1 < argc) dt_init = atof(argv[++i]);
        else if (strcmp(argv[i], "-tol") == 0 && i+1 < argc) tol = atof(argv[++i]);
        else if (strcmp(argv[i], "-log") == 0 && i+1 < argc) log_interval = atoi(argv[++i]);
        else if (strcmp(argv[i], "-snap") == 0 && i+1 < argc) snap_interval = atoi(argv[++i]);
        else if (strcmp(argv[i], "-rho0") == 0 && i+1 < argc) params.rho0 = atof(argv[++i]);
        else if (strcmp(argv[i], "-lambda") == 0 && i+1 < argc) params.lambda = atof(argv[++i]);
        else if (strcmp(argv[i], "-e") == 0 && i+1 < argc) params.e_skyrme = atof(argv[++i]);
        else if (strcmp(argv[i], "-mu") == 0 && i+1 < argc) params.mu = atof(argv[++i]);
        else if (strcmp(argv[i], "-perturb") == 0 && i+1 < argc) perturb_amp = atof(argv[++i]);
        else if (strcmp(argv[i], "-sigma") == 0 && i+1 < argc) sigma_model = atoi(argv[++i]);
        else {
            fprintf(stderr, "Usage: %s [options]\n", argv[0]);
            fprintf(stderr, "  -N <int>       Grid size (default %d)\n", 64);
            fprintf(stderr, "  -L <float>     Half-width (default 8.0)\n");
            fprintf(stderr, "  -R <float>     Initial soliton size (default 1.5)\n");
            fprintf(stderr, "  -steps <int>   Max steps (default 50000)\n");
            fprintf(stderr, "  -dt <float>    Initial step size (default 1e-4)\n");
            fprintf(stderr, "  -tol <float>   Convergence tolerance (default 1e-8)\n");
            fprintf(stderr, "  -log <int>     Log interval (default 100)\n");
            fprintf(stderr, "  -snap <int>    Snapshot interval (default 5000)\n");
            fprintf(stderr, "  -rho0 <float>  Vacuum density (default 1.0)\n");
            fprintf(stderr, "  -lambda <float> Bulk coupling (default 100.0)\n");
            fprintf(stderr, "  -e <float>     Skyrme coupling (default 4.0)\n");
            fprintf(stderr, "  -mu <float>    Degenerate mass (default 1.0)\n");
            fprintf(stderr, "  -perturb <float> Perturbation amplitude (default 0.0)\n");
            fprintf(stderr, "  -sigma <0|1>    Sigma model projection (default 1)\n");
            return 1;
        }
    }

    printf("=== CHPT Soliton Search ===\n");
    printf("Grid: %d^3 = %d points\n", N, N*N*N);
    printf("Domain: [%.1f, %.1f]^3, h = %.4f\n", -L, L, 2.0*L/N);
    printf("Parameters: rho0=%.3f, lambda=%.3f, e=%.3f, mu=%.3f, c=%.3f\n",
           params.rho0, params.lambda, params.e_skyrme, params.mu, params.c);
    printf("Initial: hedgehog with R=%.3f\n", R_init);
    printf("Flow: max_steps=%d, dt_init=%.1e, tol=%.1e, sigma_model=%d\n",
           max_steps, dt_init, tol, sigma_model);
    printf("Threads: %d\n", omp_get_max_threads());
    printf("\n");

    /* Allocate field and force */
    Field *f = field_alloc(N, L);
    int N3 = N * N * N;
    Multivector *force = (Multivector *)malloc((size_t)N3 * sizeof(Multivector));
    if (!force) { fprintf(stderr, "malloc failed\n"); return 1; }

    /* Initialize */
    init_hedgehog(f, &params, R_init);

    /* Add perturbation to degenerate sector if requested */
    if (perturb_amp > 0) {
        init_perturb(f, perturb_amp, 42);
        printf("Added perturbation with amplitude %.3e\n", perturb_amp);
    }

    /* Open log file */
    FILE *logfp = fopen("energy.log", "w");
    if (!logfp) { fprintf(stderr, "Cannot open energy.log\n"); return 1; }
    fprintf(logfp, "# step  tau  E_total  E2  E4  EV  ED  Q  max_force\n");

    /* Initial diagnostics */
    Energy en = field_energy(f, &params);
    double Q = field_topological_charge(f, &params);
    printf("Initial: E=%.8e (E2=%.4e E4=%.4e EV=%.4e ED=%.4e) Q=%.4f\n",
           en.Etotal, en.E2, en.E4, en.EV, en.ED, Q);

    /* Save initial field */
    write_field_binary("field_0000.bin", f);

    /* Gradient flow loop */
    double dt = dt_init;
    double E_prev = en.Etotal;
    int snap_count = 1;
    double tau = 0;

    clock_t t_start = clock();

    for (int step = 1; step <= max_steps; step++) {
        /* Compute gradient (force = -dE/dPsi) */
        field_gradient(f, &params, force);
        double mf = max_force_norm(force, N3);

        /* Adaptive step size: Armijo-like backtracking */
        /* Save current field once, then try steps with decreasing dt */
        Multivector *psi_save = (Multivector *)malloc((size_t)N3 * sizeof(Multivector));
        memcpy(psi_save, f->psi, (size_t)N3 * sizeof(Multivector));

        int accepted = 0;
        for (int tries = 0; tries < 10; tries++) {
            /* Restore field to saved state */
            if (tries > 0)
                memcpy(f->psi, psi_save, (size_t)N3 * sizeof(Multivector));

            /* Take step */
            gradient_flow_step(f, force, dt);

            /* Check energy */
            Energy en_new = field_energy(f, &params);

            if (en_new.Etotal <= E_prev + 1e-15 * fabs(E_prev)) {
                /* Accept step */
                if (sigma_model)
                    project_sigma_model(f, params.rho0);
                en = field_energy(f, &params);  /* recompute after projection */
                E_prev = en.Etotal;
                tau += dt;
                /* Try increasing dt for next step */
                dt = fmin(dt * 1.1, dt_max);
                accepted = 1;
                break;
            } else {
                /* Reject: halve dt and retry */
                dt *= 0.5;
                if (tries == 9)
                    printf("WARNING: step %d failed after 10 attempts, dt=%.2e\n", step, dt);
            }
        }
        free(psi_save);

        if (!accepted) continue;

        /* Logging */
        if (step % log_interval == 0 || step == 1) {
            Q = field_topological_charge(f, &params);
            write_energy_log(logfp, step, tau, en, Q, mf);

            if (step % (log_interval * 10) == 0 || step == 1) {
                double elapsed = (double)(clock() - t_start) / CLOCKS_PER_SEC;
                printf("Step %6d: E=%.8e (E2=%.4e E4=%.4e EV=%.4e ED=%.4e) "
                       "Q=%.4f |F|=%.2e dt=%.2e [%.1fs]\n",
                       step, en.Etotal, en.E2, en.E4, en.EV, en.ED,
                       Q, mf, dt, elapsed);
            }
        }

        /* Snapshots */
        if (step % snap_interval == 0) {
            char fname[64];
            snprintf(fname, sizeof(fname), "field_%04d.bin", snap_count++);
            write_field_binary(fname, f);
        }

        /* Convergence check: relative energy change */
        if (step > 100 && mf < tol * params.rho0) {
            printf("\nConverged at step %d: max_force = %.2e < tol*rho0 = %.2e\n",
                   step, mf, tol * params.rho0);
            break;
        }
    }

    /* Final diagnostics */
    Q = field_topological_charge(f, &params);
    printf("\n=== Final Result ===\n");
    printf("E_total = %.10e\n", en.Etotal);
    printf("E2      = %.10e\n", en.E2);
    printf("E4      = %.10e\n", en.E4);
    printf("EV      = %.10e\n", en.EV);
    printf("ED      = %.10e\n", en.ED);
    printf("Q       = %.6f\n", Q);
    printf("Virial: E2 - E4 + 3*(EV+ED) = %.6e (should be ~0)\n",
           en.E2 - en.E4 + 3.0*(en.EV + en.ED));
    printf("Mass formula: Mc^2 = 2*E4 - 2*(EV+ED) = %.10e\n",
           2.0*en.E4 - 2.0*(en.EV + en.ED));

    /* Degenerate sector check */
    double max_J = 0, max_P = 0;
    for (int i = 0; i < N3; i++) {
        double jn = f->psi[i].j1*f->psi[i].j1 + f->psi[i].j2*f->psi[i].j2 + f->psi[i].j3*f->psi[i].j3;
        double pn = f->psi[i].p * f->psi[i].p;
        if (jn > max_J) max_J = jn;
        if (pn > max_P) max_P = pn;
    }
    printf("|J|_max = %.6e\n", sqrt(max_J));
    printf("|P|_max = %.6e\n", sqrt(max_P));
    if (sqrt(max_J) < 1e-10 && sqrt(max_P) < 1e-10) {
        printf("=> Degenerate sector confirmed to be zero (decoupled, as predicted).\n");
    }

    /* Save final field */
    write_field_binary("field_final.bin", f);
    printf("Field saved to field_final.bin\n");

    double elapsed = (double)(clock() - t_start) / CLOCKS_PER_SEC;
    printf("Total time: %.1f seconds\n", elapsed);

    /* Cleanup */
    fclose(logfp);
    free(force);
    field_free(f);

    return 0;
}
