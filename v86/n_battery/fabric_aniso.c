/*  fabric_aniso.c — how discrete can a CUBIC fabric be before its own symmetry
 *  leaks into the physics?
 *
 *  WHY THIS MEASUREMENT DECIDES THE KERNEL REDESIGN
 *    The census's objects are fluid because the project always runs at
 *    width/spacing ~ 10-20, where the lattice is invisible (D8b measured the
 *    Peierls-Nabarro barrier at R/E ~ 1e-8 there). Making the fabric physically
 *    discrete means running at width/spacing ~ 1-3 instead. But a CUBIC lattice
 *    at that ratio does not just become discrete -- it becomes CUBIC, and the
 *    anisotropy shows up as an orientation-dependent energy, i.e. a particle
 *    whose mass depends on which way it points. That is fatal for emergent
 *    physics in a way that ordinary discretisation error is not.
 *
 *    So the design question "regular lattice vs irregular/amorphous cell
 *    complex" reduces to a number: at the ratio where PN pinning becomes
 *    physically significant, how large is the DIRECTIONAL SPREAD of that
 *    pinning? If the barrier is large but isotropic, a cubic fabric is usable.
 *    If barrier and anisotropy arrive together, it is not, and the fabric must
 *    be irregular (statistically isotropic) by construction.
 *
 *  WHAT IS MEASURED
 *    A Q-ball profile is slid through one lattice period along [100], [110] and
 *    [111]. The energy is the LATTICE energy whose variation gives the kernel's
 *    own 7-point Laplacian:
 *        E = sum_n { 1/2 sum_a [ w^2 f_a^2 + m^2 f_a^2 ]
 *                  + 1/2 sum_a sum_d ((f_a(n+d) - f_a(n))/dx)^2
 *                  + Vt(s) } dx^3
 *    The Peierls-Nabarro barrier is max E - min E over one period. Reported per
 *    direction, as a fraction of E, against the ratio r_half/dx.
 *
 *  Build: gcc -O3 -march=native -fopenmp -o fabric_aniso fabric_aniso.c -lm
 *  Usage: fabric_aniso profile.txt [omega] [rhalf]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#define NF 3
static const double M2 = 2.25, MU = -41.345, KAP = 50.0;

static double *pr, *pf;
static int pn_pts;

static void load_profile(const char *path) {
    FILE *f = fopen(path, "r");
    if (!f) { fprintf(stderr, "cannot open %s\n", path); exit(1); }
    int cap = 4096; pr = malloc(cap * sizeof(double)); pf = malloc(cap * sizeof(double));
    pn_pts = 0; char line[512];
    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '#') continue;
        double a, b;
        if (sscanf(line, "%lf %lf", &a, &b) != 2) continue;
        if (pn_pts >= cap) { cap *= 2; pr = realloc(pr, cap * sizeof(double));
                             pf = realloc(pf, cap * sizeof(double)); }
        pr[pn_pts] = a; pf[pn_pts] = b; pn_pts++;
    }
    fclose(f);
}

/* ANALYTIC bump. The table-interpolated profile is unusable for this
 * measurement: as the centre shifts, the interpolation nodes move relative to
 * the sample points, and that spurious variation scales like dx^2 -- swamping
 * the genuine Peierls-Nabarro barrier, which is EXPONENTIAL in width/dx. Using
 * a closed-form bump removes the table entirely, so the only offset dependence
 * left is real lattice pinning. */
static double AMP = 0.633, RW = 3.30, USE_ANALYTIC = 1;
static inline double prof_an(double r) {
    double u = r / RW;
    return AMP / cosh(u * u);          /* smooth, localised, exactly evaluable */
}
static inline double prof(double r) {
    if (USE_ANALYTIC) return prof_an(r);
    if (r >= pr[pn_pts - 1]) return 0.0;
    int lo = 0, hi = pn_pts - 1;
    while (hi - lo > 1) { int mid = (lo + hi) / 2; if (pr[mid] <= r) lo = mid; else hi = mid; }
    double t = (r - pr[lo]) / (pr[hi] - pr[lo]);
    return pf[lo] + t * (pf[hi] - pf[lo]);
}

static double Vt(double s) { return 0.5 * MU * s / (1.0 + KAP * s); }

/* lattice energy of a ball centred at (cx,cy,cz); the gradient term uses the
 * FORWARD link difference, which is the energy whose variation is the kernel's
 * 7-point Laplacian. */
static double energy(int N, double dx, double w, double cx, double cy, double cz) {
    double L = 0.5 * (N - 1) * dx;
    double tot = 0.0;
    #pragma omp parallel for reduction(+:tot) schedule(static)
    for (int i = 0; i < N; i++) {
        double x = -L + i * dx;
        for (int j = 0; j < N; j++) {
            double y = -L + j * dx;
            for (int k = 0; k < N; k++) {
                double z = -L + k * dx;
                double r0 = sqrt((x-cx)*(x-cx) + (y-cy)*(y-cy) + (z-cz)*(z-cz));
                double f0 = prof(r0);
                double s = f0 * f0 * f0 * f0 * f0 * f0;      /* three equal comps */
                double e = 0.5 * NF * (w * w * f0 * f0 + M2 * f0 * f0) + Vt(s);
                /* forward links in +x,+y,+z */
                double xp = sqrt((x+dx-cx)*(x+dx-cx) + (y-cy)*(y-cy) + (z-cz)*(z-cz));
                double yp = sqrt((x-cx)*(x-cx) + (y+dx-cy)*(y+dx-cy) + (z-cz)*(z-cz));
                double zp = sqrt((x-cx)*(x-cx) + (y-cy)*(y-cy) + (z+dx-cz)*(z+dx-cz));
                double gx = (prof(xp) - f0) / dx, gy = (prof(yp) - f0) / dx,
                       gz = (prof(zp) - f0) / dx;
                e += 0.5 * NF * (gx * gx + gy * gy + gz * gz);
                tot += e;
            }
        }
    }
    return tot * dx * dx * dx;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s profile.txt [omega] [rhalf]\n", argv[0]); return 1; }
    load_profile(argv[1]);
    double w = (argc > 2) ? atof(argv[2]) : 1.42;
    double rhalf = (argc > 3) ? atof(argv[3]) : 0.0;
    if (rhalf <= 0) {                      /* find r where f = f(0)/2 */
        double h = pf[0] * 0.5;
        for (int i = 1; i < pn_pts; i++) if (pf[i] <= h) { rhalf = pr[i]; break; }
    }

    printf("=========================================================================\n");
    printf("fabric anisotropy: does a CUBIC discrete fabric leak its symmetry?\n");
    printf("=========================================================================\n");
    printf("  profile %s   omega=%.4f   r_half=%.4f\n", argv[1], w, rhalf);
    printf("  PN barrier = (maxE-minE)/E over one lattice period, per direction\n\n");
    printf("  %7s %8s %12s %12s %12s %12s %11s\n",
           "dx", "width/dx", "E", "PN[100]/E", "PN[110]/E", "PN[111]/E", "anisotropy");

    /* directions and their primitive lattice period in units of dx */
    const double dirs[3][3] = { {1,0,0}, {1,1,0}, {1,1,1} };
    const char *dn[3] = { "100", "110", "111" };

    if (getenv("FA_TABLE")) USE_ANALYTIC = 0;
    if (getenv("FA_RW")) RW = atof(getenv("FA_RW"));
    printf("  mode: %s   bump width RW=%.3f\n\n",
           USE_ANALYTIC ? "ANALYTIC bump (no table interpolation)" : "table profile", RW);
    double dxs[] = { 2.20, 1.80, 1.50, 1.20, 1.00, 0.80, 0.60, 0.45 };
    for (unsigned q = 0; q < sizeof(dxs)/sizeof(dxs[0]); q++) {
        double dx = dxs[q];
        int N = (int)(2.0 * 4.0 * RW / dx);         /* box out to 4 widths */
        if (N % 2 == 0) N++;
        if (N > 221) N = 221;
        if (N < 21) N = 21;
        double pn[3], Ebase = 0;
        for (int d = 0; d < 3; d++) {
            double emin = 1e300, emax = -1e300;
            int NS = 17;
            for (int is = 0; is < NS; is++) {
                double u = (double)is / NS;         /* fraction of one period */
                double cx = u * dx * dirs[d][0];
                double cy = u * dx * dirs[d][1];
                double cz = u * dx * dirs[d][2];
                double e = energy(N, dx, w, cx, cy, cz);
                if (e < emin) emin = e;
                if (e > emax) emax = e;
                if (is == 0) Ebase = e;
            }
            pn[d] = (emax - emin) / fabs(Ebase);
        }
        double mx = pn[0], mn = pn[0];
        for (int d = 1; d < 3; d++) { if (pn[d] > mx) mx = pn[d]; if (pn[d] < mn) mn = pn[d]; }
        double aniso = (mx > 0) ? (mx - mn) / mx : 0.0;
        printf("  %7.3f %8.2f %12.4f %12.3e %12.3e %12.3e %10.1f%%\n",
               dx, RW / dx, Ebase, pn[0], pn[1], pn[2], 100.0 * aniso);
        (void)dn;
    }

    printf("\n  READING\n");
    printf("  * PN/E is the DISCRETENESS available: how strongly the fabric pins\n");
    printf("    and localises. At the census's operating point (width/dx ~ 10)\n");
    printf("    it is ~1e-8 -- the fabric is invisible and the physics is fluid.\n");
    printf("  * anisotropy is the PRICE: the spread of that pinning between\n");
    printf("    crystal directions. It is the fraction by which a particle's\n");
    printf("    energy depends on which way it is oriented / moving.\n");
    printf("  * If both rise together, a CUBIC fabric cannot be made discrete\n");
    printf("    without becoming crystalline, and the fabric must be irregular\n");
    printf("    (statistically isotropic) by construction.\n");
    return 0;
}
