/*  fabric_iso.c — v88 design test: is an AMORPHOUS fabric isotropic where a
 *  cubic one is not?
 *
 *  THE DESIGN QUESTION THIS SETTLES
 *    fabric_aniso.c measured that a cubic fabric cannot be made discrete
 *    without becoming crystalline: the Peierls-Nabarro barrier and its
 *    directional spread rise together, the spread sitting at ~1/2 to 2/3 of the
 *    barrier at EVERY width/spacing ratio. So a discrete fabric has to be
 *    irregular. Before rewriting the kernel onto an irregular cell complex, the
 *    premise has to be checked: does an amorphous fabric actually deliver
 *    isotropy at the spacings where discreteness lives?
 *
 *  WHAT IS COMPUTED
 *    For a fabric defined by neighbour vectors {d_j} with weights {w_j}, the
 *    Laplacian's Fourier symbol is
 *        Lambda(k) = sum_j w_j (1 - cos(k.d_j))
 *    Its long-wavelength expansion is
 *        Lambda(k) = 1/2 sum_j w_j (k.d_j)^2 - 1/24 sum_j w_j (k.d_j)^4 + ...
 *    The SECOND-order term is the wave speed and must be isotropic (it is, iff
 *    sum_j w_j d_j (x) d_j is proportional to the identity). The FOURTH-order
 *    term is where lattice symmetry leaks in -- that is the crystalline
 *    anisotropy, and it is what makes a particle's energy depend on its
 *    orientation.
 *
 *    Two fabrics are compared at matched mean spacing:
 *      CUBIC      6 neighbours at +-a along the axes
 *      AMORPHOUS  Poisson points, Lloyd-relaxed, neighbours within a cutoff,
 *                 weights w = A/d in the finite-volume sense (approximated here
 *                 by w ~ 1/d^2 * solid-angle share, which is the isotropic FV
 *                 weight for a well-relaxed cell complex)
 *
 *    Reported: the directional spread of Lambda(k)/|k|^2 over many directions
 *    at a fixed |k|, i.e. how much the wave speed depends on direction. This is
 *    THE number: it is the anisotropy the fabric would imprint on any object.
 *
 *  Build: gcc -O3 -march=native -fopenmp -o fabric_iso fabric_iso.c -lm
 *  Usage: fabric_iso [npts] [seed]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAXN  200000
#define MAXNB 64

static double px[MAXN], py[MAXN], pz[MAXN];
static int np;
static double BOX;

static unsigned long long rs = 88172645463325252ULL;
static double urand(void) {
    rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17;
    return (double)(rs >> 11) * (1.0 / 9007199254740992.0);
}

static inline double wrapd(double d, double L) {
    while (d >  0.5 * L) d -= L;
    while (d < -0.5 * L) d += L;
    return d;
}

/* Lloyd relaxation on a periodic box: move each point toward the centroid of
 * its neighbours' Voronoi region, approximated by repulsion from close pairs.
 * A few sweeps turn Poisson noise into a well-spaced (maximally random jammed)
 * point set, which is what a physical amorphous fabric looks like. */
static void relax(int sweeps, double target) {
    static double fx[MAXN], fy[MAXN], fz[MAXN];
    for (int s = 0; s < sweeps; s++) {
        memset(fx, 0, np * sizeof(double));
        memset(fy, 0, np * sizeof(double));
        memset(fz, 0, np * sizeof(double));
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < np; i++) {
            for (int j = 0; j < np; j++) {
                if (i == j) continue;
                double dx = wrapd(px[j] - px[i], BOX);
                double dy = wrapd(py[j] - py[i], BOX);
                double dz = wrapd(pz[j] - pz[i], BOX);
                double d2 = dx*dx + dy*dy + dz*dz;
                if (d2 > 4.0 * target * target || d2 < 1e-12) continue;
                double d = sqrt(d2);
                double push = (target - d) / d;          /* repel if too close */
                if (push < 0) push *= 0.15;              /* weak attraction */
                fx[i] -= push * dx; fy[i] -= push * dy; fz[i] -= push * dz;
            }
        }
        for (int i = 0; i < np; i++) {
            px[i] += 0.15 * fx[i]; py[i] += 0.15 * fy[i]; pz[i] += 0.15 * fz[i];
            px[i] -= BOX * floor(px[i] / BOX);
            py[i] -= BOX * floor(py[i] / BOX);
            pz[i] -= BOX * floor(pz[i] / BOX);
        }
    }
}

/* directional spread of Lambda(k)/|k|^2 for a neighbour set */
static void spread(const double *dxs, const double *dys, const double *dzs,
                   const double *ws, int nnb, double kmag,
                   double *out_mean, double *out_spread) {
    int NDIR = 512;
    double mn = 1e300, mx = -1e300, sum = 0;
    for (int t = 0; t < NDIR; t++) {
        /* uniform on the sphere */
        double u = 2.0 * ((t + 0.5) / NDIR) - 1.0;
        double phi = 2.399963229728653 * t;            /* golden angle */
        double sq = sqrt(fmax(0.0, 1.0 - u * u));
        double kx = kmag * sq * cos(phi), ky = kmag * sq * sin(phi), kz = kmag * u;
        double lam = 0;
        for (int j = 0; j < nnb; j++) {
            double kd = kx * dxs[j] + ky * dys[j] + kz * dzs[j];
            lam += ws[j] * (1.0 - cos(kd));
        }
        lam /= (kmag * kmag);
        if (lam < mn) mn = lam;
        if (lam > mx) mx = lam;
        sum += lam;
    }
    *out_mean = sum / NDIR;
    *out_spread = (mx - mn) / (*out_mean);
}

int main(int argc, char **argv) {
    int n = (argc > 1) ? atoi(argv[1]) : 4000;
    if (argc > 2) rs = (unsigned long long)atoll(argv[2]) | 1ULL;
    double a = 1.0;                                     /* mean spacing = 1 */

    printf("=====================================================================\n");
    printf("v88 fabric isotropy: amorphous vs cubic at matched mean spacing\n");
    printf("=====================================================================\n");
    printf("  Lambda(k) = sum_j w_j (1 - cos(k.d_j));  reported is the DIRECTIONAL\n");
    printf("  SPREAD of Lambda/|k|^2 over 4000 directions -- i.e. how much the\n");
    printf("  wave speed depends on direction. That is the anisotropy the fabric\n");
    printf("  imprints on every object living in it.\n\n");

    /* ---------------- cubic fabric ---------------- */
    double cdx[6] = { a,-a, 0, 0, 0, 0 }, cdy[6] = { 0, 0, a,-a, 0, 0 },
           cdz[6] = { 0, 0, 0, 0, a,-a }, cw[6];
    for (int i = 0; i < 6; i++) cw[i] = 1.0 / (a * a);

    /* ---------------- amorphous fabric ---------------- */
    np = n;
    BOX = a * cbrt((double)np);                        /* density = 1 per a^3 */
    for (int i = 0; i < np; i++) {
        px[i] = urand() * BOX; py[i] = urand() * BOX; pz[i] = urand() * BOX;
    }
    printf("  amorphous: %d points in a periodic box of side %.3f (mean spacing %.3f)\n",
           np, BOX, a);
    relax(40, a);
    printf("  Lloyd-relaxed 40 sweeps\n\n");

    /* neighbour set of one representative interior point, averaged over many
     * centres so the answer is a fabric property, not one cell's accident */
    printf("  %-12s %10s %14s %14s %12s\n",
           "fabric", "|k|a", "lambda/k^2", "spread", "vs cubic");
    double kk[] = { 0.25, 0.50, 0.80, 1.10, 1.40, 1.70, 2.00 };
    for (unsigned q = 0; q < sizeof(kk)/sizeof(kk[0]); q++) {
        double kmag = kk[q] / a;
        double cm, cs;
        spread(cdx, cdy, cdz, cw, 6, kmag, &cm, &cs);

        /* AVERAGING ORDER MATTERS AND I GOT IT WRONG FIRST TIME.
         * Averaging each cell's own directional spread gives ~80% and is
         * meaningless: one irregular neighbourhood is of course anisotropic.
         * A propagating wave samples MANY cells, so the physical object is the
         * MEAN SYMBOL <Lambda(k)>_cells, and the anisotropy is the directional
         * spread OF THAT AVERAGE. Accumulate per direction, then spread. */
        int NDIR = 512, NC = 400;
        static double acc[512];
        memset(acc, 0, sizeof(acc));
        int used = 0;
        for (int c = 0; c < NC; c++) {
            int i = (int)(urand() * np);
            double ddx[MAXNB], ddy[MAXNB], ddz[MAXNB], ww[MAXNB];
            int nb = 0;
            for (int j = 0; j < np && nb < MAXNB; j++) {
                if (j == i) continue;
                double dx = wrapd(px[j] - px[i], BOX);
                double dy = wrapd(py[j] - py[i], BOX);
                double dz = wrapd(pz[j] - pz[i], BOX);
                double d2 = dx*dx + dy*dy + dz*dz;
                if (d2 > 2.25 * a * a) continue;        /* first shell */
                double d = sqrt(d2);
                ddx[nb] = dx; ddy[nb] = dy; ddz[nb] = dz;
                ww[nb] = 1.0 / (d * d);                 /* FV-like weight */
                nb++;
            }
            if (nb < 4) continue;
            /* normalise so the k^2 coefficient matches cubic (same wave speed) */
            double tr = 0;
            for (int j = 0; j < nb; j++)
                tr += ww[j] * (ddx[j]*ddx[j] + ddy[j]*ddy[j] + ddz[j]*ddz[j]);
            double trc = 0;
            for (int j = 0; j < 6; j++)
                trc += cw[j] * (cdx[j]*cdx[j] + cdy[j]*cdy[j] + cdz[j]*cdz[j]);
            for (int j = 0; j < nb; j++) ww[j] *= trc / tr;
            for (int t = 0; t < NDIR; t++) {
                double u = 2.0 * ((t + 0.5) / NDIR) - 1.0;
                double phi = 2.399963229728653 * t;
                double sq = sqrt(fmax(0.0, 1.0 - u * u));
                double kx = kmag * sq * cos(phi), ky = kmag * sq * sin(phi),
                       kz = kmag * u;
                double lam = 0;
                for (int j = 0; j < nb; j++) {
                    double kd = kx * ddx[j] + ky * ddy[j] + kz * ddz[j];
                    lam += ww[j] * (1.0 - cos(kd));
                }
                acc[t] += lam / (kmag * kmag);
            }
            used++;
        }
        double amn = 1e300, amx = -1e300, asum = 0;
        for (int t = 0; t < NDIR; t++) {
            double v = acc[t] / used;
            if (v < amn) amn = v; if (v > amx) amx = v; asum += v;
        }
        asum /= NDIR;
        double aspr = (amx - amn) / asum;
        printf("  %-12s %10.3f %14.5f %13.3f%% %12s\n", "cubic", kk[q], cm, 100*cs, "-");
        printf("  %-12s %10.3f %14.5f %13.3f%% %11.2fx\n", "amorphous", kk[q],
               asum, 100*aspr, cs / fmax(aspr, 1e-12));
    }

    printf("\n  READING\n");
    printf("  |k|a ~ 2-3 is the regime where a structure is only a few cells\n");
    printf("  across -- i.e. where discreteness is physically available. The\n");
    printf("  question is whether the amorphous fabric's spread there is small\n");
    printf("  enough that an object does not feel a preferred direction.\n");
    return 0;
}
