/*  fabric_mass.c — v88: mass from the FABRIC'S OWN DENSITY, in d dimensions.
 *
 *  THE MECHANISM BEING TESTED
 *    Not a field living on a fixed grid. The dynamical variable IS the cell's
 *    size. Cells shrink or grow, shrinking is self-reinforcing and saturates
 *    against a hard core, and the shrunk regions must be paid for by expanded
 *    regions elsewhere -- so the ground state cannot be uniform. It must be
 *    LUMPY, and a lump is a definite NUMBER OF CELLS. That integer is the
 *    discreteness every soliton picture in this project has failed to supply:
 *    mass is a cell count, not a branch parameter.
 *
 *  HOW THE DENSITY ENTERS THE MAPPING
 *    This continues v88/fabric_proc.c, where the fabric was
 *        x(i) = a*i + A*delta(i)        with delta FROZEN and procedural.
 *    Here delta becomes the dynamical field xi(i), determined by minimising an
 *    energy, and the local fabric density is its divergence
 *        theta_i = div xi |_i           (theta < 0 = cell shrunk / densified)
 *    So "the density of the fabric changes in the mapping function" literally:
 *    the same displacement that defined the tessellation now carries the mass.
 *
 *  WHY THE LUMPINESS IS FORCED, NOT IMPOSED
 *    On a periodic lattice sum_i div xi = 0 IDENTICALLY -- a geometric fact,
 *    not a constraint anyone added. So the total dilation is pinned at zero:
 *    you cannot shrink everything. Any densification MUST be balanced by
 *    rarefaction elsewhere. That is exactly the frustration that turns a
 *    double-well into finite domains instead of one runaway phase. The
 *    conservation law does the work.
 *
 *  THE DIMENSION QUESTION
 *    Claim under test: this only works in roughly 3-7 dimensions. Two
 *    independent reasons to expect a LOWER bound near 3: the strain field of a
 *    dilational inclusion falls as r^{1-d}, so its self-energy integral
 *    int r^{d-1} * r^{2-2d} dr = int r^{1-d} dr converges at large r only for
 *    d > 2; and below that a lump's strain field is not normalisable so it
 *    cannot be a localised object. The upper bound is not derived here -- it is
 *    measured. d = 2..7 are all run.
 *
 *  W(theta): double well with minima at theta = 0 (normal fabric) and
 *  theta = -theta0 (shrunk, hard core). Degenerate, so neither phase wins
 *  globally and the split is set by the zero-sum rule.
 *
 *  Build: gcc -O3 -march=native -fopenmp -o fabric_mass fabric_mass.c -lm
 *  Usage: fabric_mass [dmin] [dmax]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#define DMAX 7

static int D, L, N;
static int stride[DMAX];
static double *xi;          /* displacement, D components, N sites */
static double *th;          /* dilation per site */
static double *fr;          /* force */

static double THETA0 = 1.0, LAM = 1.0, MU = 0.35;

static unsigned long long rs = 0x2545F4914F6CDD1DULL;
static double urand(void) {
    rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17;
    return (double)(rs >> 11) * (1.0 / 9007199254740992.0);
}

static inline int nbr(int idx, int a, int dir) {
    int c = (idx / stride[a]) % L;
    int nc = ((c + dir) % L + L) % L;
    return idx + (nc - c) * stride[a];
}

/* W(t) = (LAM/4)(t^2 - THETA0^2)^2 : minima at t = -THETA0 and +THETA0.
 *
 * THE FIRST VERSION OF THIS WAS WRONG AND GAVE A TRIVIAL ANSWER. It used
 * W = (LAM/4) t^2 (t+THETA0)^2, which has a minimum AT t = 0. The uniform
 * state theta == 0 then satisfies sum theta = 0 trivially AND sits at the
 * global minimum with E = 0, so the relaxation went straight to nothing:
 * 0.00% shrunk cells in every dimension.
 *
 * The zero-sum rule only forces lumpiness if the UNIFORM STATE IS DISFAVOURED.
 * Here t = 0 is a local MAXIMUM (W(0) = LAM*THETA0^4/4 > 0), so the fabric
 * wants every cell either shrunk or expanded and never neutral -- while the
 * geometry forbids shrinking them all. That conflict is the frustration, and
 * it is what makes phase separation compulsory rather than optional. */
static inline double Wp(double t) {
    return LAM * t * (t * t - THETA0 * THETA0);
}
static inline double W(double t) {
    double u = t * t - THETA0 * THETA0;
    return 0.25 * LAM * u * u;
}

static void compute_theta(void) {
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < N; i++) {
        double s = 0;
        for (int a = 0; a < D; a++)
            s += 0.5 * (xi[nbr(i, a, +1) * D + a] - xi[nbr(i, a, -1) * D + a]);
        th[i] = s;
    }
}

static double energy(void) {
    double e = 0;
    #pragma omp parallel for reduction(+:e) schedule(static)
    for (int i = 0; i < N; i++) {
        e += W(th[i]);
        for (int a = 0; a < D; a++)
            for (int b = 0; b < D; b++) {
                double g = 0.5 * (xi[nbr(i, a, +1) * D + b] - xi[nbr(i, a, -1) * D + b]);
                e += 0.5 * MU * g * g;
            }
    }
    return e;
}

static void forces(void) {
    memset(fr, 0, sizeof(double) * (size_t)N * D);
    /* dW/dxi : xi^a_j enters theta_{j-a} with +1/2 and theta_{j+a} with -1/2 */
    #pragma omp parallel for schedule(static)
    for (int j = 0; j < N; j++)
        for (int a = 0; a < D; a++)
            fr[j * D + a] -= 0.5 * (Wp(th[nbr(j, a, -1)]) - Wp(th[nbr(j, a, +1)]));
    /* elastic: mu * laplacian(xi) */
    #pragma omp parallel for schedule(static)
    for (int j = 0; j < N; j++)
        for (int b = 0; b < D; b++) {
            double lap = 0;
            for (int a = 0; a < D; a++)
                lap += xi[nbr(j, a, +1) * D + b] + xi[nbr(j, a, -1) * D + b]
                     - 2.0 * xi[j * D + b];
            fr[j * D + b] += MU * lap;
        }
}

/* connected clusters of shrunk cells (theta < -THETA0/2), 2D-connectivity */
static int clusters(int *sizes, int maxc, double *frac_shrunk) {
    static int *lab, *stack;
    lab = realloc(lab, sizeof(int) * N);
    stack = realloc(stack, sizeof(int) * N);
    for (int i = 0; i < N; i++) lab[i] = -1;
    int nc = 0, nsh = 0;
    for (int i = 0; i < N; i++) if (th[i] < -0.5 * THETA0) nsh++;
    for (int i = 0; i < N; i++) {
        if (th[i] >= -0.5 * THETA0 || lab[i] >= 0) continue;
        int sp = 0, cnt = 0;
        stack[sp++] = i; lab[i] = nc;
        while (sp) {
            int c = stack[--sp]; cnt++;
            for (int a = 0; a < D; a++)
                for (int s = -1; s <= 1; s += 2) {
                    int m = nbr(c, a, s);
                    if (th[m] < -0.5 * THETA0 && lab[m] < 0) { lab[m] = nc; stack[sp++] = m; }
                }
        }
        if (nc < maxc) sizes[nc] = cnt;
        nc++;
    }
    *frac_shrunk = (double)nsh / N;
    return nc;
}

int main(int argc, char **argv) {
    int dmin = (argc > 1) ? atoi(argv[1]) : 2;
    int dmax = (argc > 2) ? atoi(argv[2]) : 7;

    printf("==========================================================================\n");
    printf("v88 mass from fabric density: theta = div xi, W double-well, sum theta = 0\n");
    printf("==========================================================================\n");
    printf("  theta0=%.2f lambda=%.2f mu=%.2f   (minima at theta=0 and -theta0)\n",
           THETA0, LAM, MU);
    printf("  sum_i theta_i = 0 IDENTICALLY on a periodic lattice -- the zero-sum\n");
    printf("  rule is geometric, so lumpiness is forced, not imposed.\n");
    printf("  PREDICTION (stated before the runs): no stable finite lumps for d=2\n");
    printf("  (inclusion self-energy diverges); stable finite lumps for d>=3.\n\n");
    printf("  %2s %6s %10s %9s %10s %10s %10s %9s\n",
           "d", "L", "sites", "shrunk", "n_lumps", "mean size", "sd/mean", "E/site");

    int Ls[8] = { 0, 0, 256, 40, 16, 10, 7, 5 };
    for (D = dmin; D <= dmax && D <= DMAX; D++) {
        L = Ls[D];
        N = 1; for (int a = 0; a < D; a++) N *= L;
        stride[0] = 1;
        for (int a = 1; a < D; a++) stride[a] = stride[a-1] * L;

        xi = realloc(xi, sizeof(double) * (size_t)N * D);
        th = realloc(th, sizeof(double) * N);
        fr = realloc(fr, sizeof(double) * (size_t)N * D);
        rs = 0x2545F4914F6CDD1DULL + D;
        for (size_t q = 0; q < (size_t)N * D; q++) xi[q] = 0.35 * (urand() * 2 - 1);

        /* gradient descent with a modest step; the elastic term sets the scale */
        double dt = 0.10 / (1.0 + 2.0 * MU * D);
        int NIT = 6000;
        for (int it = 0; it < NIT; it++) {
            compute_theta();
            forces();
            #pragma omp parallel for schedule(static)
            for (size_t q = 0; q < (size_t)N * D; q++) xi[q] += dt * fr[q];
        }
        compute_theta();

        int *sz = malloc(sizeof(int) * 200000);
        double fs;
        int nc = clusters(sz, 200000, &fs);
        double m = 0, v = 0; int cap = nc < 200000 ? nc : 200000;
        for (int i = 0; i < cap; i++) m += sz[i];
        if (cap) m /= cap;
        for (int i = 0; i < cap; i++) v += (sz[i] - m) * (sz[i] - m);
        if (cap) v = sqrt(v / cap);
        printf("  %2d %6d %10d %8.2f%% %10d %10.1f %10.3f %9.5f\n",
               D, L, N, 100 * fs, nc, m, m > 0 ? v / m : 0.0, energy() / N);
        fflush(stdout);
        free(sz);
    }
    printf("\n  READING\n");
    printf("  n_lumps > 1 with a small sd/mean is the signature wanted: many\n");
    printf("  separate objects of a REPEATABLE cell count -- i.e. a quantised\n");
    printf("  mass. One giant lump (n_lumps=1) means coarsening won and there is\n");
    printf("  no preferred size. sd/mean ~ 1 means a broad distribution, i.e. no\n");
    printf("  quantisation.\n");
    return 0;
}
