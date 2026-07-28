/*  fabric_pn.c — v88: is the Kelvin (BCC 8+6) fabric's PINNING isotropic?
 *
 *  WHY THIS AND NOT lattice_iso.c
 *    lattice_iso measured the isotropy of the LINEAR DISPERSION -- the 4th
 *    moment of the neighbour set, which the tuned BCC 8+6 stencil zeroes
 *    exactly (A4 = 3, 1.49% spread at |k|a = 2 vs simple cubic's 22.69%).
 *    That is necessary but NOT sufficient. What killed simple cubic in
 *    v86/n_battery/fabric_aniso.c was the NONLINEAR Peierls-Nabarro pinning:
 *    the energy barrier a localised bump feels as it slides, whose directional
 *    spread sat at 40-79% at every width/spacing ratio. Dispersion isotropy and
 *    pinning isotropy are different tensors and one does not imply the other.
 *
 *  WHAT IS MEASURED
 *    A localised bump is slid along the shortest lattice translation in [100],
 *    [110] and [111], and the lattice energy
 *        E = sum_i V { (NF/2)(w^2 f^2 + m^2 f^2) + Vt(f^6)
 *                    + (NF/4) sum_j w_ij (f_j - f_i)^2 }
 *    is evaluated at each offset. PN barrier = (maxE - minE)/E over one period.
 *    SC and BCC 8+6 are compared at EQUAL NUMBER DENSITY (one site per unit
 *    volume), so "spacing" means the same for both.
 *
 *    Weight normalisation: sum_j w_j |d_j|^2 = 6 makes the graph Laplacian's
 *    symbol equal k^2 at long wavelength for BOTH fabrics, so they are compared
 *    at the same wave speed as well as the same density.
 *
 *  Build: gcc -O3 -march=native -fopenmp -o fabric_pn fabric_pn.c -lm
 *  Usage: fabric_pn [RW]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#define NF 3
static const double M2 = 2.25, MU = -41.345, KAP = 50.0;
static double AMP = 0.633, RW = 3.30, WFREQ = 1.42;

static inline double prof(double r) { double u = r / RW; return AMP / cosh(u * u); }
static inline double Vt(double s) { return 0.5 * MU * s / (1.0 + KAP * s); }

typedef struct {
    char name[32];
    int  nnb;
    double nb[16][3];      /* neighbour offsets */
    double w[16];          /* weights, normalised so sum w |d|^2 = 6 */
    double prim[3][3];     /* shortest lattice translations along 100,110,111 */
    int    two_site;       /* BCC has a body-centre sublattice */
    double a;
} Fab;

static void normalise(Fab *F) {
    double s = 0;
    for (int j = 0; j < F->nnb; j++) {
        double *d = F->nb[j];
        s += F->w[j] * (d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
    }
    for (int j = 0; j < F->nnb; j++) F->w[j] *= 6.0 / s;
}

/* energy of a bump centred at c, on fabric F, summed over a box of half-size R */
static double energy(const Fab *F, double R, const double c[3]) {
    double a = F->a;
    int n = (int)(R / a) + 2;
    double tot = 0.0;
    #pragma omp parallel for reduction(+:tot) schedule(dynamic)
    for (int i = -n; i <= n; i++) {
        for (int j = -n; j <= n; j++) {
            for (int k = -n; k <= n; k++) {
                for (int sub = 0; sub < (F->two_site ? 2 : 1); sub++) {
                    double o = sub ? 0.5 * a : 0.0;
                    double x = i * a + o, y = j * a + o, z = k * a + o;
                    double dx = x - c[0], dy = y - c[1], dz = z - c[2];
                    double r = sqrt(dx*dx + dy*dy + dz*dz);
                    if (r > R) continue;
                    double f0 = prof(r);
                    double s6 = f0*f0*f0*f0*f0*f0;
                    double e = 0.5 * NF * (WFREQ*WFREQ*f0*f0 + M2*f0*f0) + Vt(s6);
                    double g = 0.0;
                    for (int q = 0; q < F->nnb; q++) {
                        double nx = x + F->nb[q][0] - c[0];
                        double ny = y + F->nb[q][1] - c[1];
                        double nz = z + F->nb[q][2] - c[2];
                        double rn = sqrt(nx*nx + ny*ny + nz*nz);
                        double df = prof(rn) - f0;
                        g += F->w[q] * df * df;
                    }
                    e += 0.25 * NF * g;
                    tot += e;                     /* V = 1 per site by density */
                }
            }
        }
    }
    return tot;
}

static void scan(const Fab *F, double R, int dir, double *barrier, double *E0) {
    const int NS = 61;
    double mn = 1e300, mx = -1e300;
    for (int t = 0; t < NS; t++) {
        double u = (double)t / NS;
        double c[3] = { u * F->prim[dir][0], u * F->prim[dir][1], u * F->prim[dir][2] };
        double e = energy(F, R, c);
        if (t == 0) *E0 = e;
        if (e < mn) mn = e;
        if (e > mx) mx = e;
    }
    *barrier = (mx - mn) / fabs(*E0);
}

int main(int argc, char **argv) {
    if (argc > 1) RW = atof(argv[1]);

    /* ---- simple cubic, density 1 -> a = 1 ---- */
    Fab sc; memset(&sc, 0, sizeof(sc)); strcpy(sc.name, "SC 6"); sc.a = 1.0;
    { int s[6][3] = {{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
      for (int i = 0; i < 6; i++) {
          for (int c = 0; c < 3; c++) sc.nb[i][c] = s[i][c] * sc.a;
          sc.w[i] = 1.0 / (sc.a * sc.a);
      }
      sc.nnb = 6; sc.two_site = 0;
      /* shortest lattice translations */
      sc.prim[0][0] = sc.a;
      sc.prim[1][0] = sc.a; sc.prim[1][1] = sc.a;
      sc.prim[2][0] = sc.a; sc.prim[2][1] = sc.a; sc.prim[2][2] = sc.a;
      normalise(&sc); }

    /* ---- BCC 8+6 (Kelvin foam), density 1 -> 2 sites per a^3 -> a = 2^(1/3) */
    Fab bc; memset(&bc, 0, sizeof(bc)); strcpy(bc.name, "BCC 8+6 Kelvin");
    bc.a = cbrt(2.0);
    { int q = 0;
      for (int i = -1; i <= 1; i += 2)
        for (int j = -1; j <= 1; j += 2)
          for (int k = -1; k <= 1; k += 2) {
            bc.nb[q][0] = 0.5*i*bc.a; bc.nb[q][1] = 0.5*j*bc.a; bc.nb[q][2] = 0.5*k*bc.a;
            double d2 = 0; for (int c = 0; c < 3; c++) d2 += bc.nb[q][c]*bc.nb[q][c];
            bc.w[q] = 1.0 / d2;                        /* first shell, weight 1 */
            q++;
          }
      int s[6][3] = {{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
      for (int i = 0; i < 6; i++) {
          for (int c = 0; c < 3; c++) bc.nb[q][c] = s[i][c] * bc.a;
          double d2 = 0; for (int c = 0; c < 3; c++) d2 += bc.nb[q][c]*bc.nb[q][c];
          bc.w[q] = (2.0/3.0) / d2;                    /* TUNED: w2/w1 = 2/3 */
          q++;
      }
      bc.nnb = q; bc.two_site = 1;
      /* shortest BCC translations: [111] is a/2(1,1,1); [100] is a; [110] is a(1,1,0) */
      bc.prim[0][0] = bc.a;
      bc.prim[1][0] = bc.a; bc.prim[1][1] = bc.a;
      bc.prim[2][0] = 0.5*bc.a; bc.prim[2][1] = 0.5*bc.a; bc.prim[2][2] = 0.5*bc.a;
      normalise(&bc); }

    printf("=====================================================================\n");
    printf("v88 PN pinning isotropy: simple cubic vs Kelvin (BCC 8+6, w2/w1=2/3)\n");
    printf("=====================================================================\n");
    printf("  equal number density (1 site per unit volume) and equal wave speed\n");
    printf("  (sum_j w_j |d_j|^2 = 6 for both). Bump width RW = %.3f.\n", RW);
    printf("  PN barrier = (maxE-minE)/E over the SHORTEST lattice translation\n");
    printf("  in each direction.\n\n");
    printf("  %-16s %8s %11s %11s %11s %11s\n",
           "fabric", "width/a", "PN[100]", "PN[110]", "PN[111]", "anisotropy");

    double widths[] = { 2.00, 1.80, 1.60, 1.45, 1.30 };
    Fab *fabs_[2] = { &sc, &bc };
    for (unsigned q = 0; q < sizeof(widths)/sizeof(widths[0]); q++) {
        RW = widths[q];
        double R = 4.0 * RW;
        for (int fi = 0; fi < 2; fi++) {
            Fab *F = fabs_[fi];
            double pn[3], E0;
            for (int d = 0; d < 3; d++) scan(F, R, d, &pn[d], &E0);
            double mx = pn[0], mn = pn[0];
            for (int d = 1; d < 3; d++) { if (pn[d] > mx) mx = pn[d]; if (pn[d] < mn) mn = pn[d]; }
            printf("  %-16s %8.2f %11.3e %11.3e %11.3e %10.1f%%\n",
                   F->name, RW / F->a, pn[0], pn[1], pn[2],
                   mx > 0 ? 100.0 * (mx - mn) / mx : 0.0);
        }
        printf("\n");
    }
    printf("  READING: the anisotropy column is the decider. Simple cubic held\n");
    printf("  40-79%% at every ratio (v86/n_battery/fabric_aniso.c). If Kelvin's\n");
    printf("  is comparable, 4th-order dispersion isotropy did NOT buy pinning\n");
    printf("  isotropy and the design does not close.\n");
    return 0;
}
