/* shapelock.c — v89 — harmonics INSIDE a cell: the 3-axis shape spectrum.
 *
 * Standalone. Build: gcc -O2 -Wall -Wextra -o shapelock shapelock.c -lm
 *
 * WHY THIS IS A DIFFERENT QUESTION FROM harmfrust.c
 * -------------------------------------------------
 * Consonance was rejected as a bound BETWEEN cells: a degree-7 cell must
 * satisfy seven simultaneous commensurability conditions and generically
 * cannot, so only 5-13% of contacts ever locked, at every degree. That is
 * a network counting failure.
 *
 * A cell's OWN three axes are not a network. They are a closed triangle:
 *
 *      (a/b) * (b/c) * (c/a) = 1
 *
 * identically — the comb-closure condition of `charge_sym.mac`, with only
 * TWO of the three intervals free. One cell, one constraint set, no
 * neighbours to disagree with. So the mechanism that killed consonance
 * between cells is structurally absent within one.
 *
 * That is the literal reading of "harmonics within, dissonance between".
 *
 * WHAT IS BEING ASKED
 * -------------------
 * S1  Is there a discrete SPECTRUM of harmonically admissible shapes —
 *     integer-ratio ellipsoids whose three internal intervals all sit on
 *     low-Tenney-height rungs? How many, and does the count saturate or
 *     explode with the height ceiling?
 * S2  Under anisotropic (constant-volume) load, does a cell's shape LOCK
 *     to one of them and jump between them — a staircase in shape?
 * S3  Negative control: comb off must give a smooth response.
 * S4  Where are the LIMITS — at what load does the comb lose the shape,
 *     and what sets the ceiling?
 *
 * PARAMETRISATION
 * ---------------
 * x = ln a, y = ln b, z = ln c with x + y + z = 0, so the volume (and
 * hence the cell's energy) is fixed and only SHAPE moves. This is the
 * constant-volume soft mode a scalar radius does not have — the whole
 * point of going to three axes.
 *
 * Intervals u1 = log2(a/b), u2 = log2(b/c), u3 = log2(c/a); u1+u2+u3 = 0
 * identically, which is the closure, not an imposed constraint.
 *
 *   U = D(u1) + D(u2) + D(u3) - sum_i sigma_i x_i
 *   D(u) = min over rungs [ (u-u_k)^2/(2 s^2) + lam * H_k ]
 *
 * sigma is deviatoric (sum sigma = 0): pure shear, no volume change.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct { int p, q; double u, H; } Rung;
static Rung RUNGS[8192];
static int NRUNG = 0;

static int igcd(int a, int b) { while (b) { int t = a % b; a = b; b = t; } return a; }

static void build_rungs(double climb)
{
    NRUNG = 0;
    /* climb < 0 means "unison only" — the negative control. The height
     * test below rejects even 1/1 (0 > -1), leaving NRUNG = 0: energy
     * became 3e300, the backtracking comparison degenerated to equality
     * in floating point, relaxation exited after one iteration, and the
     * control read a flat response that looked like a pass. It had not
     * run. */
    if (climb < 0) {
        RUNGS[0].p = 1; RUNGS[0].q = 1; RUNGS[0].u = 0.0; RUNGS[0].H = 0.0;
        NRUNG = 1;
        return;
    }
    int lim = (int)(pow(2.0, climb > 0 ? climb : 0.0) + 0.5);
    if (lim < 1) lim = 1;
    for (int p = 1; p <= lim; p++)
        for (int q = 1; q <= lim; q++) {
            if (p * q > lim || igcd(p, q) != 1) continue;
            double H = log2((double)p * (double)q);
            if (H > climb + 1e-12 || NRUNG >= 8192) continue;
            RUNGS[NRUNG].p = p; RUNGS[NRUNG].q = q;
            RUNGS[NRUNG].u = log2((double)p / (double)q);
            RUNGS[NRUNG].H = H; NRUNG++;
        }
}

static double diss(double u, double s, double lam, int *kb)
{
    double best = 1e300; int b = 0;
    for (int k = 0; k < NRUNG; k++) {
        double d = u - RUNGS[k].u;
        double v = d * d / (2 * s * s) + lam * RUNGS[k].H;
        if (v < best) { best = v; b = k; }
    }
    if (kb) *kb = b;
    return best;
}
static double dgrad(double u, double s, double lam)
{ int k; diss(u, s, lam, &k); return (u - RUNGS[k].u) / (s * s); }

/* ------------------------------------------------- S1: the spectrum */

/* An admissible SHAPE is an integer triple (na:nb:nc), gcd 1, whose three
 * pairwise ratios ALL have Tenney height <= climb. Enumerated exactly in
 * integer arithmetic — no relaxation, no search, so this section cannot
 * be an artifact of a solver. */
static int height_ok(int p, int q, double climb)
{
    int g = igcd(p, q); p /= g; q /= g;
    return log2((double)p * (double)q) <= climb + 1e-12;
}

static int enumerate_shapes(double climb, int nmax, int *out, int maxout)
{
    int n = 0;
    for (int a = 1; a <= nmax; a++)
        for (int b = a; b <= nmax; b++)
            for (int c = b; c <= nmax; c++) {
                if (igcd(igcd(a, b), c) != 1) continue;
                if (!height_ok(b, a, climb)) continue;
                if (!height_ok(c, b, climb)) continue;
                if (!height_ok(c, a, climb)) continue;
                if (n < maxout) { out[3*n] = a; out[3*n+1] = b; out[3*n+2] = c; }
                n++;
            }
    return n;
}

/* --------------------------------------------------- the shape cell */

typedef struct { double x, y, z, s, lam, kE, tx, ty, tz; } Cell;

static void center(Cell *c)
{ double m = (c->x + c->y + c->z) / 3.0; c->x -= m; c->y -= m; c->z -= m; }

static double u1(Cell *c) { return (c->x - c->y) / M_LN2; }
static double u2(Cell *c) { return (c->y - c->z) / M_LN2; }
static double u3(Cell *c) { return (c->z - c->x) / M_LN2; }

/* QUADRATIC preference about a target shape, NOT a linear pull. A linear
 * drive is scale-free — unbounded benefit — so it always runs to the
 * outermost rung and the rung structure selects nothing. Diagnosed and
 * fixed once in harmgrow.c, then reintroduced here. */
static double energy(Cell *c)
{
    double dx = c->x - c->tx, dy = c->y - c->ty, dz = c->z - c->tz;
    return diss(u1(c), c->s, c->lam, NULL)
         + diss(u2(c), c->s, c->lam, NULL)
         + diss(u3(c), c->s, c->lam, NULL)
         + 0.5 * c->kE * (dx * dx + dy * dy + dz * dz);
}

static void grad(Cell *c, double g[3])
{
    double d1 = dgrad(u1(c), c->s, c->lam) / M_LN2;
    double d2 = dgrad(u2(c), c->s, c->lam) / M_LN2;
    double d3 = dgrad(u3(c), c->s, c->lam) / M_LN2;
    g[0] =  d1 - d3 + c->kE * (c->x - c->tx);
    g[1] = -d1 + d2 + c->kE * (c->y - c->ty);
    g[2] = -d2 + d3 + c->kE * (c->z - c->tz);
    double m = (g[0] + g[1] + g[2]) / 3.0;   /* project: volume fixed */
    g[0] -= m; g[1] -= m; g[2] -= m;
}

static double relax(Cell *c, int iters, double lr, double *gm)
{
    double U = energy(c), step = lr, g[3];
    for (int it = 0; it < iters; it++) {
        grad(c, g);
        double bx = c->x, by = c->y, bz = c->z, Un = 0; int ok = 0;
        for (int t = 0; t < 40; t++) {
            c->x = bx - step * g[0]; c->y = by - step * g[1];
            c->z = bz - step * g[2];
            center(c);
            Un = energy(c);
            if (isfinite(Un) && Un <= U) { ok = 1; break; }
            step *= 0.5;
        }
        if (!ok) { c->x = bx; c->y = by; c->z = bz; break; }
        if (U - Un < 1e-16 * (1 + fabs(U))) { U = Un; break; }
        U = Un; step *= 1.1; if (step > 1.0) step = 1.0;
    }
    grad(c, g);
    *gm = fmax(fabs(g[0]), fmax(fabs(g[1]), fabs(g[2])));
    return energy(c);
}

/* Rung-seeded exhaustive search: seed at every (rung, rung) pair for the
 * two free intervals. The charge work and harmgrow both showed random
 * restarts report the basin nearest the start as if it were the answer. */
static double solve(Cell *c, double *gm)
{
    double bx = 0, by = 0, bz = 0, Ub = 0, gb = 0; int have = 0;
    for (int k1 = 0; k1 < NRUNG; k1++)
        for (int k2 = 0; k2 < NRUNG; k2++) {
            double a1 = RUNGS[k1].u * M_LN2, a2 = RUNGS[k2].u * M_LN2;
            c->x = (2 * a1 + a2) / 3.0;
            c->y = c->x - a1;
            c->z = c->y - a2;
            center(c);
            double g2, U = relax(c, 500, 0.01, &g2);
            if (isfinite(U) && (!have || U < Ub)) {
                Ub = U; gb = g2; have = 1;
                bx = c->x; by = c->y; bz = c->z;
            }
        }
    c->x = bx; c->y = by; c->z = bz;
    *gm = gb;
    return Ub;
}

int main(void)
{
    setvbuf(stdout, NULL, _IOLBF, 0);
    printf("# shapelock — harmonics INSIDE a cell (3 axes, closed triangle)\n");
    printf("# (a/b)(b/c)(c/a) = 1 identically: closure, not a constraint.\n");
    printf("# volume fixed (x+y+z=0), so only SHAPE moves — the constant-\n");
    printf("# volume soft mode a scalar radius does not have.\n\n");

    /* ---------------- S1 ---------------- */
    printf("== S1. THE SHAPE SPECTRUM (exact integer enumeration) ==\n");
    printf("admissible = integer triple (a:b:c), gcd 1, with ALL THREE\n");
    printf("pairwise ratios at Tenney height <= ceiling. No solver here, so\n");
    printf("this section cannot be a search artifact.\n");
    printf(" ceiling   n_shapes(<=12)  n_shapes(<=32)   first few\n");
    {
        int out[3 * 4096];
        for (double cl = 0; cl <= 8.001; cl += 1.0) {
            int n12 = enumerate_shapes(cl, 12, out, 4096);
            int n32 = enumerate_shapes(cl, 32, out, 4096);
            printf("%8.0f  %14d  %14d   ", cl, n12, n32);
            int show = n32 < 5 ? n32 : 5;
            for (int i = 0; i < show; i++)
                printf("%d:%d:%d ", out[3*i], out[3*i+1], out[3*i+2]);
            printf("\n");
        }
    }
    printf("\nreading: a spectrum that GROWS with the ceiling but stays small\n");
    printf("and enumerable is a shape quantisation. If it explodes, the comb\n");
    printf("selects nothing; if it stays at 1 (the sphere), it forbids shape.\n\n");

    /* ---------------- S2 ---------------- */
    printf("== S2. DOES SHAPE LOCK UNDER SHEAR? ==\n");
    printf("Pure deviatoric load sigma = t*(2,-1,-1)/3 (volume preserved).\n");
    printf("comb_limit=4, s=0.05, lam=0.10. Rung-seeded exhaustive search.\n");
    printf("  t_pref  a:b ratio    b:c ratio    rung1   rung2    max|grad|\n");
    {
        build_rungs(4.0);
        for (double t = 0; t <= 3.001; t += 0.15) {
            Cell c = {0, 0, 0, 0.05, 0.10, 1.0, 2*t/3, -t/3, -t/3};
            double gm;
            solve(&c, &gm);
            int k1, k2;
            diss(u1(&c), c.s, c.lam, &k1);
            diss(u2(&c), c.s, c.lam, &k2);
            printf("%6.2f  %11.5f  %11.5f  %3d/%-3d %3d/%-3d  %10.1e\n",
                   t, exp(c.x - c.y), exp(c.y - c.z),
                   RUNGS[k1].p, RUNGS[k1].q, RUNGS[k2].p, RUNGS[k2].q, gm);
        }
    }
    printf("\nreading: flat runs = the shape is LOCKED to a rung pair and the\n");
    printf("shear is absorbed by the lock; jumps = the lock breaking. This is\n");
    printf("the staircase harmgrow found on a ring, but here there are NO\n");
    printf("neighbours to frustrate it — the closure is internal.\n\n");

    /* ---------------- S3 ---------------- */
    printf("== S3. NEGATIVE CONTROL — comb collapsed to the unison ==\n");
    printf("  t_pref  a:b ratio    b:c ratio\n");
    {
        build_rungs(-1.0);
        for (double t = 0; t <= 3.001; t += 0.375) {
            Cell c = {0, 0, 0, 0.05, 0.10, 1.0, 2*t/3, -t/3, -t/3};
            double gm;
            solve(&c, &gm);
            printf("%6.2f  %11.5f  %11.5f\n",
                   t, exp(c.x - c.y), exp(c.y - c.z));
        }
    }
    printf("\nreading: must be smooth. A staircase here would falsify S2.\n\n");

    /* ---------------- S4 ---------------- */
    printf("== S4. THE LIMITS — where does the comb lose the shape? ==\n");
    printf("Sweeping the ceiling at fixed strong shear: the outermost rung\n");
    printf("caps the axis ratio, so the comb span IS the aspect-ratio bound.\n");
    printf(" ceiling  nrungs   a:b at t=6   max rung ratio   locked?\n");
    {
        for (double cl = 1; cl <= 7.001; cl += 1.0) {
            build_rungs(cl);
            Cell c = {0, 0, 0, 0.05, 0.10, 1.0, 4.0, -2.0, -2.0};
            double gm;
            solve(&c, &gm);
            int k1; double off = 0;
            diss(u1(&c), c.s, c.lam, &k1);
            off = fabs(u1(&c) - RUNGS[k1].u);
            double hi = 0;
            for (int k = 0; k < NRUNG; k++) if (RUNGS[k].u > hi) hi = RUNGS[k].u;
            printf("%8.0f  %6d  %12.4f  %14.1f   %s (off %.1e)\n",
                   cl, NRUNG, exp(c.x - c.y), pow(2, hi),
                   off < 1e-3 ? "yes" : "NO ", off);
        }
    }
    printf("\nreading: a:b saturating at the outermost rung means the comb span\n");
    printf("sets the maximum aspect ratio — a bound on SHAPE made of\n");
    printf("intervals, with no length anywhere. That is the bound the size\n");
    printf("story could not deliver.\n");
    return 0;
}
