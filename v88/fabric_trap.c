/*  fabric_trap.c — v88: energy crosses cells only in COMPLETE CYCLES, so a
 *  region whose cycle rate differs from its surroundings cannot leak.
 *
 *  WHY THE PREVIOUS MODEL FAILED
 *    v88/fabric_mass.c made cell size dynamical with a double-well and the
 *    geometric zero-sum rule sum_i div xi = 0. It phase-separates (40% shrunk)
 *    but COARSENS: d=2 gives 78 lumps with sd/mean = 5.4, d=3 gives 2, d=4
 *    gives one blob of 26674 cells. No preferred size, no quantised mass.
 *    The reason is that it has NO CYCLES. Nothing sets a scale and nothing
 *    traps energy, so domains merge without limit.
 *
 *  THE MECHANISM HERE
 *    Each cell carries an internal cycle (complex amplitude psi = A e^{i phi}).
 *    Transfer between cells is RESONANT: the coupling moves energy only to the
 *    extent that neighbouring cycles keep step. The cycle rate is shifted by
 *    the cell's own amplitude,
 *        omega_i = omega_0 - gamma |psi_i|^2
 *    so a cell holding more cyclic energy runs at a different rate -- "more
 *    cyclic energy, tighter cell" in the size picture. Once a cluster's rate
 *    leaves the band its neighbours can support, no complete cycle connects it
 *    to the exterior and the energy is reflected back inside. That is the
 *    sharp exterior mismatch doing the trapping, and it is why the object has a
 *    finite size instead of coarsening.
 *
 *    Nothing is imported. There is one field, psi, and one nonlinearity. The
 *    trapping configuration is a pattern in psi, not an added scalar -- an
 *    intra-particle CONFIGURATION, emergent by construction.
 *
 *  THE BAND, WHICH IS WHAT MAKES IT DISCRETE
 *    Linear waves on a d-dim lattice occupy omega in [0, 4 d C]. A localised
 *    state must sit OUTSIDE that band or it resonates with the continuum and
 *    radiates. So the trapping condition is a threshold in amplitude, not a
 *    continuous dial -- and states either close their cycle or do not exist.
 *
 *  MEASURED
 *    From random initial conditions: does energy self-organise into localised
 *    objects, how many cells do they occupy (inverse participation ratio), and
 *    is that occupancy REPEATABLE (small spread) -- i.e. a quantised mass?
 *
 *  Build: gcc -O3 -march=native -fopenmp -o fabric_trap fabric_trap.c -lm
 *  Usage: fabric_trap [dmin] [dmax]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#define DMAX 7
static int D, L, N;
static int stride[DMAX];
static double *pr, *pi_, *kr, *ki;
static double C = 1.0, GAM = 4.0;

static unsigned long long rs = 0x9E3779B97F4A7C15ULL;
static double urand(void) {
    rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17;
    return (double)(rs >> 11) * (1.0 / 9007199254740992.0);
}
static inline int nbr(int idx, int a, int dir) {
    int c = (idx / stride[a]) % L;
    int nc = ((c + dir) % L + L) % L;
    return idx + (nc - c) * stride[a];
}

/* i dpsi/dt = -C sum_nbr (psi_j - psi_i) - gamma |psi_i|^2 psi_i
 * so dpsi/dt = i[ C sum(psi_j - psi_i) + gamma |psi|^2 psi ] */
static void deriv(const double *ar, const double *ai, double *br, double *bi) {
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < N; i++) {
        double sr = 0, si = 0;
        for (int a = 0; a < D; a++) {
            int p = nbr(i, a, +1), m = nbr(i, a, -1);
            sr += ar[p] + ar[m] - 2.0 * ar[i];
            si += ai[p] + ai[m] - 2.0 * ai[i];
        }
        double n2 = ar[i]*ar[i] + ai[i]*ai[i];
        double rr = C * sr + GAM * n2 * ar[i];
        double ri = C * si + GAM * n2 * ai[i];
        br[i] = -ri;          /* i*(rr + i ri) = -ri + i rr */
        bi[i] =  rr;
    }
}

int main(int argc, char **argv) {
    int dmin = (argc > 1) ? atoi(argv[1]) : 1;
    int dmax = (argc > 2) ? atoi(argv[2]) : 7;
    int Ls[8] = { 0, 4096, 200, 34, 14, 9, 6, 5 };

    printf("=====================================================================\n");
    printf("v88 cycle trapping: resonant transfer + amplitude-shifted cycle rate\n");
    printf("=====================================================================\n");
    printf("  C=%.2f gamma=%.2f   linear band omega in [0, %s]\n", C, GAM, "4dC");
    printf("  A cluster whose rate leaves the band cannot complete a cycle with\n");
    printf("  the exterior, so its energy is reflected inward. Trapping is a\n");
    printf("  THRESHOLD, not a dial -- states either close their cycle or do not\n");
    printf("  exist. That is where discreteness would come from.\n\n");
    printf("  %2s %6s %10s %11s %11s %11s %10s %10s\n",
           "d", "L", "sites", "band top", "peak |psi|^2", "cells/lump",
           "n_lumps", "sd/mean");

    for (D = dmin; D <= dmax && D <= DMAX; D++) {
        L = Ls[D];
        N = 1; for (int a = 0; a < D; a++) N *= L;
        stride[0] = 1;
        for (int a = 1; a < D; a++) stride[a] = stride[a-1] * L;
        pr = realloc(pr, sizeof(double)*N); pi_ = realloc(pi_, sizeof(double)*N);
        kr = realloc(kr, sizeof(double)*N); ki = realloc(ki, sizeof(double)*N);
        static double *tr, *ti; tr = realloc(tr, sizeof(double)*N);
        ti = realloc(ti, sizeof(double)*N);

        rs = 0x9E3779B97F4A7C15ULL + 1000u*D;
        double norm0 = 0;
        for (int i = 0; i < N; i++) {
            double a = urand(), b = urand();
            pr[i] = 0.5*(a*2-1); pi_[i] = 0.5*(b*2-1);
            norm0 += pr[i]*pr[i] + pi_[i]*pi_[i];
        }

        double dt = 0.004, T = (getenv("FT_T") ? atof(getenv("FT_T")) : 400.0);
        int NT = (int)(T/dt);
        for (int t = 0; t < NT; t++) {          /* RK2 midpoint, norm-preserving enough */
            deriv(pr, pi_, kr, ki);
            #pragma omp parallel for schedule(static)
            for (int i = 0; i < N; i++) { tr[i] = pr[i] + 0.5*dt*kr[i];
                                          ti[i] = pi_[i] + 0.5*dt*ki[i]; }
            deriv(tr, ti, kr, ki);
            #pragma omp parallel for schedule(static)
            for (int i = 0; i < N; i++) { pr[i] += dt*kr[i]; pi_[i] += dt*ki[i]; }
        }

        /* inverse participation ratio -> how many cells actually hold the energy */
        double s2 = 0, s4 = 0, pk = 0;
        for (int i = 0; i < N; i++) {
            double n2 = pr[i]*pr[i] + pi_[i]*pi_[i];
            s2 += n2; s4 += n2*n2; if (n2 > pk) pk = n2;
        }
        double ipr = (s2*s2) / (N * s4);          /* fraction of cells participating */
        double cells_tot = ipr * N;

        /* count peaks well above the background, and their spread */
        double thr = 0.25 * pk;
        int np = 0; double m = 0, v = 0;
        static double *pv; pv = realloc(pv, sizeof(double)*N);
        for (int i = 0; i < N; i++) {
            double n2 = pr[i]*pr[i] + pi_[i]*pi_[i];
            if (n2 < thr) continue;
            int loc = 1;
            for (int a = 0; a < D && loc; a++)
                for (int s = -1; s <= 1 && loc; s += 2) {
                    int mI = nbr(i,a,s);
                    double q = pr[mI]*pr[mI] + pi_[mI]*pi_[mI];
                    if (q > n2) loc = 0;
                }
            if (loc) { pv[np] = n2; m += n2; np++; }
        }
        if (np) m /= np;
        for (int i = 0; i < np; i++) v += (pv[i]-m)*(pv[i]-m);
        if (np) v = sqrt(v/np);
        printf("  %2d %6d %10d %11.2f %11.4f %11.1f %10d %10.3f\n",
               D, L, N, 4.0*D*C, pk, np ? cells_tot/np : 0.0, np,
               m > 0 ? v/m : 0.0);
        fflush(stdout);
    }
    printf("\n  READING\n");
    printf("  cells/lump small and REPEATABLE (sd/mean << 1) across many lumps\n");
    printf("  is the signature: localised objects of a definite occupancy, i.e.\n");
    printf("  a quantised mass set by a cycle-closure condition rather than by a\n");
    printf("  branch parameter. sd/mean ~ 1 means no quantisation.\n");
    return 0;
}
