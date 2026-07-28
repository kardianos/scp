/*  lattice_iso.c — v88: which CRYSTALLINE fabric is isotropic enough?
 *
 *  THE POINT
 *    fabric_aniso.c showed a simple-cubic fabric cannot be made discrete
 *    without becoming crystalline: 23% directional spread in wave speed at
 *    |k|a ~ 2. fabric_iso.c showed an amorphous fabric fixes that (<=5%) but
 *    pays with no translation symmetry at all, no momentum conservation,
 *    ensemble averaging, and irregular memory.
 *
 *    A periodic non-cubic fabric may get the isotropy WITHOUT those costs. The
 *    anisotropy of a lattice Laplacian sits in the 4th moment
 *        T4 = sum_j w_j d_j (x) d_j (x) d_j (x) d_j
 *    A neighbour set is 4th-order isotropic iff it is a spherical 4-design.
 *    The octahedron (SC, 6 nbrs) and the cube (BCC, 8 nbrs) are only
 *    3-designs -- hence SC's failure. But for a cubic-symmetric set T4 has just
 *    ONE independent anisotropy, measured here by
 *        A4 = sum_j w_j dx^4 / sum_j w_j dx^2 dy^2      (isotropic <=> A4 = 3)
 *    so ONE tunable shell-weight ratio can cancel it exactly, while keeping
 *    exact periodicity.
 *
 *    Note the Kelvin foam is not a separate option: its cell is the truncated
 *    octahedron, whose 14 faces are exactly BCC's 8 + 6 neighbour shells. So
 *    "foam lattice" and "two-shell BCC" are the same fabric.
 *
 *  All candidates are compared at EQUAL NUMBER DENSITY (one cell per unit
 *  volume), so "mean spacing" means the same thing for each.
 *
 *  Build: gcc -O3 -march=native -o lattice_iso lattice_iso.c -lm
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAXNB 64

typedef struct {
    char name[48];
    int  n;
    double d[MAXNB][3];
    double w[MAXNB];
    int shell[MAXNB];          /* 0 = first shell, 1 = second */
} Lat;

static void addv(Lat *L, double x, double y, double z, int sh) {
    L->d[L->n][0] = x; L->d[L->n][1] = y; L->d[L->n][2] = z;
    L->shell[L->n] = sh; L->w[L->n] = 1.0; L->n++;
}

/* A4 = <w dx^4> / <w dx^2 dy^2>; equals 3 for an isotropic 4th moment */
static double A4(const Lat *L) {
    double q = 0, p = 0;
    for (int j = 0; j < L->n; j++) {
        double x = L->d[j][0], y = L->d[j][1], z = L->d[j][2], w = L->w[j];
        q += w * (x*x*x*x + y*y*y*y + z*z*z*z) / 3.0;
        p += w * (x*x*y*y + y*y*z*z + z*z*x*x) / 3.0;
    }
    return (p != 0) ? q / p : INFINITY;
}

/* second-order (wave speed) coefficient: sum_j w_j |d_j|^2 / 3 */
static double C2(const Lat *L) {
    double s = 0;
    for (int j = 0; j < L->n; j++) {
        double *d = (double *)L->d[j];
        s += L->w[j] * (d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
    }
    return s / 3.0;
}

/* directional spread of Lambda(k)/|k|^2 */
static void spread(const Lat *L, double kmag, double *mean, double *spr) {
    const int ND = 20000;
    double mn = 1e300, mx = -1e300, sum = 0;
    for (int t = 0; t < ND; t++) {
        double u = 2.0 * ((t + 0.5) / ND) - 1.0;
        double phi = 2.399963229728653 * t;
        double sq = sqrt(fmax(0.0, 1.0 - u * u));
        double kx = kmag * sq * cos(phi), ky = kmag * sq * sin(phi), kz = kmag * u;
        double lam = 0;
        for (int j = 0; j < L->n; j++) {
            double kd = kx*L->d[j][0] + ky*L->d[j][1] + kz*L->d[j][2];
            lam += L->w[j] * (1.0 - cos(kd));
        }
        lam /= kmag * kmag;
        if (lam < mn) mn = lam;
        if (lam > mx) mx = lam;
        sum += lam;
    }
    *mean = sum / ND;
    *spr = (mx - mn) / *mean;
}

/* scale all vectors so number density is 1 cell per unit volume */
static void set_density(Lat *L, double cell_vol_in_a3, double a_for_unit) {
    (void)cell_vol_in_a3;
    for (int j = 0; j < L->n; j++)
        for (int c = 0; c < 3; c++) L->d[j][c] *= a_for_unit;
}

/* FV-like weight  w = 1/d^2  within a shell, with a tunable 2nd-shell factor */
static void set_weights(Lat *L, double shell2_factor) {
    for (int j = 0; j < L->n; j++) {
        double *d = (double *)L->d[j];
        double d2 = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
        L->w[j] = (1.0 / d2) * (L->shell[j] ? shell2_factor : 1.0);
    }
}

/* solve A4(factor) = 3 by bisection on the 2nd-shell weight */
static double tune(Lat *L) {
    double lo = -2.0, hi = 20.0;
    set_weights(L, lo); double flo = A4(L) - 3.0;
    set_weights(L, hi); double fhi = A4(L) - 3.0;
    if (flo * fhi > 0) return NAN;
    for (int it = 0; it < 200; it++) {
        double mid = 0.5 * (lo + hi);
        set_weights(L, mid);
        double f = A4(L) - 3.0;
        if (f * flo <= 0) { hi = mid; fhi = f; } else { lo = mid; flo = f; }
    }
    return 0.5 * (lo + hi);
}

int main(void) {
    Lat cand[8]; int nc = 0;

    /* ---- simple cubic: 1 site per a^3 -> a = 1 ; 6 nbrs at 1, 12 at sqrt2 */
    { Lat *L = &cand[nc++]; memset(L, 0, sizeof(*L)); strcpy(L->name, "SC 6");
      int s[6][3] = {{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
      for (int i = 0; i < 6; i++) addv(L, s[i][0], s[i][1], s[i][2], 0);
      set_density(L, 1, 1.0); }
    { Lat *L = &cand[nc++]; memset(L, 0, sizeof(*L)); strcpy(L->name, "SC 6+12");
      int s[6][3] = {{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
      for (int i = 0; i < 6; i++) addv(L, s[i][0], s[i][1], s[i][2], 0);
      for (int a = -1; a <= 1; a++) for (int b = -1; b <= 1; b++)
        for (int c = -1; c <= 1; c++) {
          int n2 = a*a+b*b+c*c; if (n2 != 2) continue; addv(L, a, b, c, 1); }
      set_density(L, 1, 1.0); }

    /* ---- BCC: 2 sites per a^3 -> a = 2^(1/3). 8 nbrs at a*sqrt3/2, 6 at a.
     *      This IS the Kelvin foam (truncated octahedron, 8 hex + 6 square). */
    { double a = cbrt(2.0);
      Lat *L = &cand[nc++]; memset(L, 0, sizeof(*L)); strcpy(L->name, "BCC 8");
      for (int i = -1; i <= 1; i += 2) for (int j = -1; j <= 1; j += 2)
        for (int k = -1; k <= 1; k += 2) addv(L, 0.5*i, 0.5*j, 0.5*k, 0);
      set_density(L, 2, a); }
    { double a = cbrt(2.0);
      Lat *L = &cand[nc++]; memset(L, 0, sizeof(*L));
      strcpy(L->name, "BCC 8+6 (Kelvin foam)");
      for (int i = -1; i <= 1; i += 2) for (int j = -1; j <= 1; j += 2)
        for (int k = -1; k <= 1; k += 2) addv(L, 0.5*i, 0.5*j, 0.5*k, 0);
      int s[6][3] = {{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
      for (int i = 0; i < 6; i++) addv(L, s[i][0], s[i][1], s[i][2], 1);
      set_density(L, 2, a); }

    /* ---- FCC: 4 sites per a^3 -> a = 4^(1/3). 12 nbrs at a/sqrt2, 6 at a.
     *      Voronoi cell = rhombic dodecahedron (12 faces). */
    { double a = cbrt(4.0);
      Lat *L = &cand[nc++]; memset(L, 0, sizeof(*L)); strcpy(L->name, "FCC 12");
      for (int i = -1; i <= 1; i++) for (int j = -1; j <= 1; j++)
        for (int k = -1; k <= 1; k++) {
          if (i*i + j*j + k*k != 2) continue; addv(L, 0.5*i, 0.5*j, 0.5*k, 0); }
      set_density(L, 4, a); }
    { double a = cbrt(4.0);
      Lat *L = &cand[nc++]; memset(L, 0, sizeof(*L)); strcpy(L->name, "FCC 12+6");
      for (int i = -1; i <= 1; i++) for (int j = -1; j <= 1; j++)
        for (int k = -1; k <= 1; k++) {
          if (i*i + j*j + k*k != 2) continue; addv(L, 0.5*i, 0.5*j, 0.5*k, 0); }
      int s[6][3] = {{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
      for (int i = 0; i < 6; i++) addv(L, s[i][0], s[i][1], s[i][2], 1);
      set_density(L, 4, a); }

    /* ---- icosahedral 12: a spherical 5-design, the ideal -- but NOT periodic */
    { Lat *L = &cand[nc++]; memset(L, 0, sizeof(*L));
      strcpy(L->name, "icosahedral 12 *");
      double g = 0.5 * (1.0 + sqrt(5.0));
      double v[12][3] = {{0,1,g},{0,1,-g},{0,-1,g},{0,-1,-g},
                         {1,g,0},{1,-g,0},{-1,g,0},{-1,-g,0},
                         {g,0,1},{-g,0,1},{g,0,-1},{-g,0,-1}};
      double nrm = sqrt(1.0 + g * g);
      for (int i = 0; i < 12; i++) addv(L, v[i][0]/nrm, v[i][1]/nrm, v[i][2]/nrm, 0);
      set_density(L, 1, 1.0); }

    printf("======================================================================\n");
    printf("v88 crystalline fabric isotropy -- equal number density throughout\n");
    printf("======================================================================\n");
    printf("  A4 = <w dx^4>/<w dx^2 dy^2>; ISOTROPIC 4th moment <=> A4 = 3\n");
    printf("  spread = directional variation of Lambda(k)/|k|^2 (wave speed)\n");
    printf("  * icosahedral is a spherical 5-design but is NOT compatible with\n");
    printf("    translational periodicity in 3D -- listed as the ideal bound.\n\n");

    printf("  %-24s %6s %8s %10s %10s %10s\n",
           "fabric", "nbrs", "A4", "spr@k=1.0", "spr@k=1.5", "spr@k=2.0");
    for (int i = 0; i < nc; i++) {
        set_weights(&cand[i], 1.0);
        double m, s1, s2, s3;
        spread(&cand[i], 1.0, &m, &s1);
        spread(&cand[i], 1.5, &m, &s2);
        spread(&cand[i], 2.0, &m, &s3);
        printf("  %-24s %6d %8.4f %9.3f%% %9.3f%% %9.3f%%\n",
               cand[i].name, cand[i].n, A4(&cand[i]), 100*s1, 100*s2, 100*s3);
    }

    printf("\n  TUNED TWO-SHELL STENCILS (2nd-shell weight solved for A4 = 3)\n");
    printf("  %-24s %10s %8s %10s %10s %10s\n",
           "fabric", "w2/w1", "A4", "spr@k=1.0", "spr@k=1.5", "spr@k=2.0");
    for (int i = 0; i < nc; i++) {
        int has2 = 0;
        for (int j = 0; j < cand[i].n; j++) if (cand[i].shell[j]) has2 = 1;
        if (!has2) continue;
        double f = tune(&cand[i]);
        if (isnan(f)) { printf("  %-24s   no root\n", cand[i].name); continue; }
        set_weights(&cand[i], f);
        double m, s1, s2, s3;
        spread(&cand[i], 1.0, &m, &s1);
        spread(&cand[i], 1.5, &m, &s2);
        spread(&cand[i], 2.0, &m, &s3);
        printf("  %-24s %10.5f %8.4f %9.3f%% %9.3f%% %9.3f%%\n",
               cand[i].name, f, A4(&cand[i]), 100*s1, 100*s2, 100*s3);
    }
    printf("\n  (a NEGATIVE w2/w1 means the 4th-order cancellation needs an\n");
    printf("   anti-bonded second shell, which is legitimate for a stencil but\n");
    printf("   makes the operator non-monotone -- flagged, not hidden.)\n");
    return 0;
}
