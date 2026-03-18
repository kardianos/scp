/*  t10a.c — Wave-Wave Scattering via Triple Product
 *
 *  Test whether the triple-product potential V(P) = (μ/2)P²/(1+κP²)
 *  creates wave-wave scattering that depends on relative phase.
 *  This would be the EM signal at the field level.
 *
 *  Setup: Two Gaussian wave packets in the same field (φ₀) collide head-on.
 *  The triple-product coupling means they only interact when all 3 fields
 *  are nonzero simultaneously → need waves in all 3 fields overlapping.
 *
 *  Three test configurations:
 *  A) Same-component collision: both packets in φ₀. No triple-product interaction
 *     (since P = φ₀φ₁φ₂ and φ₁=φ₂=0). This is the NULL control.
 *  B) Cross-component collision: packet in φ₀ + packet in φ₁.
 *     Interaction requires φ₂ ≠ 0 → still no scattering in vacuum.
 *  C) Three-field collision: packets in all 3 fields overlapping.
 *     P ≠ 0 where all three overlap → nonlinear interaction.
 *     Vary relative phases to see if scattering is phase-dependent.
 *  D) On braid background: packets scatter off the structured field.
 *
 *  Build: gcc -O3 -fopenmp -o t10a src/t10a.c -lm
 */

#include "../../src/braid_core.h"
#include <float.h>

/* Measure total energy in a region */
static double measure_energy_region(Grid *g, const double *phys,
                                     double xmin, double xmax,
                                     double ymin, double ymax) {
    int N = g->N, NN = N*N;
    double dx = g->dx, L = g->L, dV = dx*dx*dx;
    double mu = phys[12], kappa = phys[13], mass2 = phys[14]*phys[14];
    double E = 0;

    for (int i = 1; i < N-1; i++) {
        double x = -L + i * dx;
        if (x < xmin || x > xmax) continue;
        for (int j = 1; j < N-1; j++) {
            double y = -L + j * dx;
            if (y < ymin || y > ymax) continue;
            for (int k = 0; k < N; k++) {
                int idx = i*NN + j*N + k;
                double ek = 0, eg = 0;
                for (int a = 0; a < NFIELDS; a++) {
                    ek += 0.5 * g->vel[a][idx] * g->vel[a][idx];
                    int kp = (k+1)%N, km = (k-1+N)%N;
                    double gx = (g->phi[a][idx+NN]-g->phi[a][idx-NN])/(2*dx);
                    double gy = (g->phi[a][idx+N]-g->phi[a][idx-N])/(2*dx);
                    double gz = (g->phi[a][i*NN+j*N+kp]-g->phi[a][i*NN+j*N+km])/(2*dx);
                    eg += 0.5*(gx*gx + gy*gy + gz*gz);
                }
                double p0 = g->phi[0][idx], p1 = g->phi[1][idx], p2 = g->phi[2][idx];
                double rho = p0*p0 + p1*p1 + p2*p2;
                double P = p0*p1*p2;
                double em = 0.5 * mass2 * rho;
                double ep = (mu/2.0)*P*P/(1.0+kappa*P*P);
                E += (ek + eg + em + ep) * dV;
            }
        }
    }
    return E;
}

/* Measure field amplitude in a region */
static double measure_amp_region(Grid *g, int comp,
                                  double xmin, double xmax,
                                  double ymin, double ymax) {
    int N = g->N, NN = N*N;
    double dx = g->dx, L = g->L;
    double max_amp = 0;
    for (int i = 1; i < N-1; i++) {
        double x = -L + i * dx;
        if (x < xmin || x > xmax) continue;
        for (int j = 1; j < N-1; j++) {
            double y = -L + j * dx;
            if (y < ymin || y > ymax) continue;
            for (int k = 0; k < N; k++) {
                int idx = i*NN + j*N + k;
                double a = fabs(g->phi[comp][idx]);
                if (a > max_amp) max_amp = a;
            }
        }
    }
    return max_amp;
}

/* Initialize a single wave packet */
static void add_wave_packet(Grid *g, int comp,
                            double x0, double y0, double z0,
                            double kx, double ky, double kz,
                            double sigma, double amplitude, double phase_offset) {
    int N = g->N, NN = N*N;
    double dx = g->dx, L = g->L;
    double inv_2s2 = 1.0 / (2.0 * sigma * sigma);
    double mass2 = 2.25;
    double k2 = kx*kx + ky*ky + kz*kz;
    double omega = sqrt(k2 + mass2);

    for (int i = 0; i < N; i++) {
        double x = -L + i * dx;
        for (int j = 0; j < N; j++) {
            double y = -L + j * dx;
            for (int kk = 0; kk < N; kk++) {
                double z = -L + kk * dx;
                int idx = i*NN + j*N + kk;
                double rx = x - x0, ry = y - y0, rz = z - z0;
                double r2 = rx*rx + ry*ry + rz*rz;
                double env = amplitude * exp(-r2 * inv_2s2);
                double phase = kx*rx + ky*ry + kz*rz + phase_offset;
                g->phi[comp][idx] += env * cos(phase);
                g->vel[comp][idx] += -omega * env * cos(phase);
            }
        }
    }
}

/* ================================================================
   Main experiment
   ================================================================ */

int main(int argc, char **argv) {
    int N = 96;
    double L = 25.0;  /* larger box for scattering */
    double T_run = 30.0;
    double amp = 0.3;   /* larger amplitude to see nonlinear effects */
    double sigma = 2.0;
    double k_wave = 3.0;

    double mass2 = 2.25;
    double omega = sqrt(k_wave*k_wave + mass2);
    double vg = k_wave / omega;

    printf("=== T10A: Wave-Wave Scattering via Triple Product ===\n\n");
    printf("Grid: N=%d, L=%.1f, dx=%.4f\n", N, L, 2.0*L/(N-1));
    printf("Wave: k=%.1f, omega=%.4f, vg=%.4f, amp=%.2f, sigma=%.1f\n",
           k_wave, omega, vg, amp, sigma);
    printf("Run time: T=%.0f\n\n", T_run);

    bimodal_init_params();
    double dt = 0.20 * (2.0*L/(N-1));
    int run_steps = (int)(T_run / dt);
    int meas_every = run_steps / 20;
    if (meas_every < 1) meas_every = 1;

    FILE *fout = fopen("data/t10a_scattering.tsv", "w");
    fprintf(fout, "test\tphase_offset\tt\tE_left\tE_right\tE_center\tmax_P\n");

    /* ============================================================
       Test 1: Single packet (control — no scattering expected)
       ============================================================ */
    printf("=== Test 1: Single packet (control) ===\n");
    {
        Grid *g = grid_alloc(N, L);
        grid_zero(g);
        add_wave_packet(g, 0, -10, 0, 0, k_wave, 0, 0, sigma, amp, 0);
        compute_forces(g, BIMODAL, 2.25);

        double E0 = measure_energy_region(g, BIMODAL, -L, L, -L, L);
        printf("  Initial energy: %.2f\n", E0);

        for (int s = 0; s < run_steps; s++) {
            verlet_full_step(g, BIMODAL, 2.25);
            apply_damping_xy(g);
            if (s % meas_every == 0) {
                double El = measure_energy_region(g, BIMODAL, -L, -2, -L, L);
                double Er = measure_energy_region(g, BIMODAL, 2, L, -L, L);
                double Ec = measure_energy_region(g, BIMODAL, -2, 2, -L, L);
                /* Max |P| */
                double maxP = 0;
                int NN = N*N;
                for (int idx = 0; idx < N*N*N; idx++) {
                    double P = fabs(g->phi[0][idx]*g->phi[1][idx]*g->phi[2][idx]);
                    if (P > maxP) maxP = P;
                }
                fprintf(fout, "single\t0\t%.3f\t%.2f\t%.2f\t%.2f\t%.4e\n",
                        s*dt, El, Er, Ec, maxP);
            }
        }
        double Ef = measure_energy_region(g, BIMODAL, -L, L, -L, L);
        printf("  Final energy: %.2f (loss from damping: %.1f%%)\n", Ef, 100*(1-Ef/E0));
        grid_free(g);
    }

    /* ============================================================
       Test 2: Head-on collision, same component (NULL control)
       Two packets in φ₀ moving toward each other. P=0 since φ₁=φ₂=0.
       ============================================================ */
    printf("\n=== Test 2: Same-component collision (P=0 control) ===\n");
    {
        Grid *g = grid_alloc(N, L);
        grid_zero(g);
        add_wave_packet(g, 0, -10, 0, 0, +k_wave, 0, 0, sigma, amp, 0);
        add_wave_packet(g, 0, +10, 0, 0, -k_wave, 0, 0, sigma, amp, 0);
        compute_forces(g, BIMODAL, 2.25);

        double E0 = measure_energy_region(g, BIMODAL, -L, L, -L, L);
        printf("  Initial energy: %.2f\n", E0);

        for (int s = 0; s < run_steps; s++) {
            verlet_full_step(g, BIMODAL, 2.25);
            apply_damping_xy(g);
            if (s % meas_every == 0) {
                double El = measure_energy_region(g, BIMODAL, -L, -2, -L, L);
                double Er = measure_energy_region(g, BIMODAL, 2, L, -L, L);
                double Ec = measure_energy_region(g, BIMODAL, -2, 2, -L, L);
                fprintf(fout, "same_comp\t0\t%.3f\t%.2f\t%.2f\t%.2f\t0\n",
                        s*dt, El, Er, Ec);
            }
        }

        /* After collision: measure energy in left vs right half */
        double El_f = measure_energy_region(g, BIMODAL, -L, 0, -L, L);
        double Er_f = measure_energy_region(g, BIMODAL, 0, L, -L, L);
        printf("  After: E_left=%.2f, E_right=%.2f (should be ~equal)\n", El_f, Er_f);
        printf("  → Pure superposition (no scattering), as expected\n");
        grid_free(g);
    }

    /* ============================================================
       Test 3: Three-field collision — the key test
       Packet in φ₀ from left, φ₁ from right, φ₂ from left.
       Where all three overlap: P ≠ 0 → nonlinear interaction.
       Vary relative phases.
       ============================================================ */
    printf("\n=== Test 3: Three-field collision (P≠0 → interaction!) ===\n");

    double phases[] = {0, PI/4, PI/2, PI, 3*PI/2};
    int nphases = 5;

    for (int ip = 0; ip < nphases; ip++) {
        double dphi = phases[ip];
        printf("  Phase offset = %.2f (%.0f deg): ", dphi, dphi*180/PI);
        fflush(stdout);

        Grid *g = grid_alloc(N, L);
        grid_zero(g);
        /* φ₀ from left */
        add_wave_packet(g, 0, -10, 0, 0, +k_wave, 0, 0, sigma, amp, 0);
        /* φ₁ from right */
        add_wave_packet(g, 1, +10, 0, 0, -k_wave, 0, 0, sigma, amp, dphi);
        /* φ₂ from left (co-propagating with φ₀) */
        add_wave_packet(g, 2, -10, 0, 0, +k_wave, 0, 0, sigma, amp, 0);
        compute_forces(g, BIMODAL, 2.25);

        double E0 = measure_energy_region(g, BIMODAL, -L, L, -L, L);

        /* Measure initial P overlap */
        double maxP_init = 0;
        for (int idx = 0; idx < N*N*N; idx++) {
            double P = fabs(g->phi[0][idx]*g->phi[1][idx]*g->phi[2][idx]);
            if (P > maxP_init) maxP_init = P;
        }

        for (int s = 0; s < run_steps; s++) {
            verlet_full_step(g, BIMODAL, 2.25);
            apply_damping_xy(g);
            if (s % meas_every == 0) {
                double El = measure_energy_region(g, BIMODAL, -L, -2, -L, L);
                double Er = measure_energy_region(g, BIMODAL, 2, L, -L, L);
                double Ec = measure_energy_region(g, BIMODAL, -2, 2, -L, L);
                double maxP = 0;
                for (int idx = 0; idx < N*N*N; idx++) {
                    double P = fabs(g->phi[0][idx]*g->phi[1][idx]*g->phi[2][idx]);
                    if (P > maxP) maxP = P;
                }
                fprintf(fout, "3field\t%.4f\t%.3f\t%.2f\t%.2f\t%.2f\t%.4e\n",
                        dphi, s*dt, El, Er, Ec, maxP);
            }
        }

        /* After collision: check for energy redistribution */
        double El_f = measure_energy_region(g, BIMODAL, -L, 0, -L, L);
        double Er_f = measure_energy_region(g, BIMODAL, 0, L, -L, L);
        double Ef = El_f + Er_f;
        double asym = (El_f - Er_f) / (Ef + 1e-30);

        /* Check for field mixing: did energy transfer between components? */
        double amp0_right = measure_amp_region(g, 0, 5, L, -5, 5);
        double amp1_left = measure_amp_region(g, 1, -L, -5, -5, 5);
        double amp2_right = measure_amp_region(g, 2, 5, L, -5, 5);

        printf("E0=%.1f, E_final=%.1f, L/R asym=%.4f, |P|_init=%.3e",
               E0, Ef, asym, maxP_init);
        printf(", φ₀ transmitted=%.4f, φ₁ transmitted=%.4f", amp0_right, amp1_left);
        if (fabs(asym) > 0.01)
            printf(" [PHASE-DEPENDENT!]");
        printf("\n");

        grid_free(g);
    }

    /* ============================================================
       Test 4: Amplitude dependence — does scattering scale as amp²?
       ============================================================ */
    printf("\n=== Test 4: Amplitude dependence ===\n");
    printf("  amp\tE_total\tL/R_asym\tmax|P|\n");

    double amps[] = {0.05, 0.1, 0.2, 0.5, 1.0, 2.0};
    int namps = 6;

    for (int ia = 0; ia < namps; ia++) {
        double a = amps[ia];
        Grid *g = grid_alloc(N, L);
        grid_zero(g);
        add_wave_packet(g, 0, -10, 0, 0, +k_wave, 0, 0, sigma, a, 0);
        add_wave_packet(g, 1, +10, 0, 0, -k_wave, 0, 0, sigma, a, 0);
        add_wave_packet(g, 2, -10, 0, 0, +k_wave, 0, 0, sigma, a, 0);
        compute_forces(g, BIMODAL, 2.25);

        double E0 = measure_energy_region(g, BIMODAL, -L, L, -L, L);
        double maxP = 0;
        for (int idx = 0; idx < N*N*N; idx++) {
            double P = fabs(g->phi[0][idx]*g->phi[1][idx]*g->phi[2][idx]);
            if (P > maxP) maxP = P;
        }

        for (int s = 0; s < run_steps; s++) {
            verlet_full_step(g, BIMODAL, 2.25);
            apply_damping_xy(g);
        }

        double El = measure_energy_region(g, BIMODAL, -L, 0, -L, L);
        double Er = measure_energy_region(g, BIMODAL, 0, L, -L, L);
        double asym = (El - Er) / (El + Er + 1e-30);

        printf("  %.2f\t%.1f\t%.6f\t%.4e\n", a, E0, asym, maxP);

        fprintf(fout, "ampscan\t%.4f\t%.3f\t%.2f\t%.2f\t%.2f\t%.4e\n",
                a, T_run, El, Er, 0.0, maxP);

        grid_free(g);
    }

    /* ============================================================
       Test 5: Transverse scattering (perpendicular collision)
       φ₀ from left (+x), φ₁ from bottom (+y), φ₂ from left (+x).
       Measure: does anything scatter at 90° or other angles?
       ============================================================ */
    printf("\n=== Test 5: Perpendicular collision ===\n");
    {
        Grid *g = grid_alloc(N, L);
        grid_zero(g);
        add_wave_packet(g, 0, -10, 0, 0, +k_wave, 0, 0, sigma, amp, 0);
        add_wave_packet(g, 1, 0, -10, 0, 0, +k_wave, 0, sigma, amp, 0);
        add_wave_packet(g, 2, -10, 0, 0, +k_wave, 0, 0, sigma, amp, 0);
        compute_forces(g, BIMODAL, 2.25);

        double E0 = measure_energy_region(g, BIMODAL, -L, L, -L, L);

        for (int s = 0; s < run_steps; s++) {
            verlet_full_step(g, BIMODAL, 2.25);
            apply_damping_xy(g);
        }

        /* Check four quadrants */
        double Epp = measure_energy_region(g, BIMODAL, 2, L, 2, L);    /* +x,+y */
        double Epm = measure_energy_region(g, BIMODAL, 2, L, -L, -2);  /* +x,-y */
        double Emp = measure_energy_region(g, BIMODAL, -L, -2, 2, L);  /* -x,+y */
        double Emm = measure_energy_region(g, BIMODAL, -L, -2, -L, -2);/* -x,-y */

        printf("  E0=%.1f\n", E0);
        printf("  Quadrants: (+x,+y)=%.2f, (+x,-y)=%.2f, (-x,+y)=%.2f, (-x,-y)=%.2f\n",
               Epp, Epm, Emp, Emm);
        double scattered = Epp + Emm; /* off-diagonal = scattered */
        double transmitted = Epm + Emp; /* along original directions */
        printf("  Scattered (diagonal): %.2f\n", scattered);
        printf("  Transmitted (straight): %.2f\n", transmitted);
        printf("  Scattering fraction: %.4f\n", scattered / (scattered + transmitted + 1e-30));

        grid_free(g);
    }

    fclose(fout);

    printf("\n=== T10A Complete ===\n");
    printf("Data in data/t10a_scattering.tsv\n");
    return 0;
}
