/* geopress.c — v89 — GEOMETRIC PRESSURE as the morphology bound.
 *
 * Standalone. Build: gcc -O2 -Wall -Wextra -o geopress geopress.c -lm
 *
 * THE PROPOSAL
 * ------------
 * Consonance was rejected as a morphology bound (`harmfrust.c`): a cell at
 * degree 7 must satisfy seven simultaneous commensurability conditions and
 * generically cannot, so only 5-13% of contacts ever locked, at every
 * degree, with no trend. That is a counting failure, not a tuning failure.
 *
 * A PRESSURE balance is one scalar equation per cell — internal drive
 * against the sum of contact forces — so there is nothing to frustrate.
 * This file asks whether that survives the same test consonance failed,
 * and whether it can hold a localised structure.
 *
 * THE MODEL
 * ---------
 * Cells in a periodic box. A cell's extent follows its own energy,
 *
 *      r_i = (E_i)^(1/3)
 *
 * which is "matter is converted space" as a size law: no absolute length
 * appears anywhere, only extents and separations, both of which are state.
 * Overlapping cells repel with contact energy
 *
 *      U_c = sum_contacts (k/2) delta_ij^2,  delta = r_i + r_j - d_ij > 0
 *
 * Energy is conserved (sum E fixed) and flows down the gradient of its own
 * chemical potential,
 *
 *      mu_i = dU_c/dE_i = (1/3) E_i^(-2/3) * P_i,   P_i = sum_j k delta_ij
 *
 * so a cell under high contact pressure SHEDS energy and one under low
 * pressure absorbs it. That is "pressure pushes, nothing reaches out" —
 * the standing space-transport law — arrived at from the geometry rather
 * than imposed.
 *
 * NOTHING IN U CONTAINS A MAXIMUM SIZE. If a bound appears it is the
 * balance between a cell's drive to grow and its neighbours' pressure.
 *
 * WHAT WOULD COUNT
 * ----------------
 * G1 a bounded, stable size distribution that does NOT degrade with
 *    coordination — the test consonance failed;
 * G2 a localised energy excess that grows to a size and STOPS, held by
 *    the surrounding pressure, rather than running away or dissolving;
 * G3 the negative control: with contact stiffness off, no bound and no
 *    localisation, so the effect is the pressure's.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    int N;
    double L, k, lam;
    double *x, *y, *z, *E, *r, *P, *mu;
} Box;

static unsigned long long RS = 20260802ULL;
static void rseed(unsigned long long s) { RS = s ? s : 1; }
static double rnd(void)
{ RS ^= RS << 13; RS ^= RS >> 7; RS ^= RS << 17;
  return (double)((RS >> 11) & 0xFFFFFFFFFFFFFULL) / 9007199254740992.0; }

static Box *box_new(int N, double L, double k, double lam)
{
    Box *b = calloc(1, sizeof(Box));
    b->N = N; b->L = L; b->k = k; b->lam = lam;
    b->x = calloc(N, sizeof(double)); b->y = calloc(N, sizeof(double));
    b->z = calloc(N, sizeof(double)); b->E = calloc(N, sizeof(double));
    b->r = calloc(N, sizeof(double)); b->P = calloc(N, sizeof(double));
    b->mu = calloc(N, sizeof(double));
    return b;
}
static void box_free(Box *b)
{ free(b->x); free(b->y); free(b->z); free(b->E);
  free(b->r); free(b->P); free(b->mu); free(b); }

static void box_init(Box *b, double E0, unsigned long long seed)
{
    rseed(seed);
    for (int i = 0; i < b->N; i++) {
        b->x[i] = b->L * rnd(); b->y[i] = b->L * rnd(); b->z[i] = b->L * rnd();
        b->E[i] = E0 * (0.8 + 0.4 * rnd());
    }
}

static double wrap(double d, double L)
{ while (d >  0.5 * L) d -= L; while (d < -0.5 * L) d += L; return d; }

static void radii(Box *b)
{ for (int i = 0; i < b->N; i++) b->r[i] = cbrt(b->E[i]); }

/* contact pressure per cell, mean contact number, total contact energy */
static double contacts(Box *b, double *zbar, double *Uc)
{
    int N = b->N; double L = b->L;
    memset(b->P, 0, N * sizeof(double));
    long nc = 0; double U = 0;
    for (int i = 0; i < N; i++)
        for (int j = i + 1; j < N; j++) {
            double dx = wrap(b->x[i] - b->x[j], L);
            double dy = wrap(b->y[i] - b->y[j], L);
            double dz = wrap(b->z[i] - b->z[j], L);
            double d = sqrt(dx * dx + dy * dy + dz * dz);
            double s = b->r[i] + b->r[j];
            if (d >= s || d < 1e-12) continue;
            double del = s - d;
            b->P[i] += b->k * del; b->P[j] += b->k * del;
            U += 0.5 * b->k * del * del;
            nc++;
        }
    if (zbar) *zbar = 2.0 * nc / N;
    if (Uc) *Uc = U;
    return U;
}

/* relax POSITIONS to mechanical equilibrium (energies fixed) */
static double relax_pos(Box *b, int iters, double step, double *fmax)
{
    int N = b->N; double L = b->L;
    double *fx = calloc(N, sizeof(double)), *fy = calloc(N, sizeof(double)),
           *fz = calloc(N, sizeof(double));
    double U = 0;
    for (int it = 0; it < iters; it++) {
        memset(fx, 0, N * sizeof(double)); memset(fy, 0, N * sizeof(double));
        memset(fz, 0, N * sizeof(double));
        U = 0;
        for (int i = 0; i < N; i++)
            for (int j = i + 1; j < N; j++) {
                double dx = wrap(b->x[i] - b->x[j], L);
                double dy = wrap(b->y[i] - b->y[j], L);
                double dz = wrap(b->z[i] - b->z[j], L);
                double d = sqrt(dx * dx + dy * dy + dz * dz);
                double s = b->r[i] + b->r[j];
                if (d >= s || d < 1e-12) continue;
                double del = s - d, f = b->k * del / d;
                fx[i] += f * dx; fy[i] += f * dy; fz[i] += f * dz;
                fx[j] -= f * dx; fy[j] -= f * dy; fz[j] -= f * dz;
                U += 0.5 * b->k * del * del;
            }
        for (int i = 0; i < N; i++) {
            b->x[i] += step * fx[i]; b->y[i] += step * fy[i];
            b->z[i] += step * fz[i];
            if (b->x[i] < 0) b->x[i] += L; if (b->x[i] >= L) b->x[i] -= L;
            if (b->y[i] < 0) b->y[i] += L; if (b->y[i] >= L) b->y[i] -= L;
            if (b->z[i] < 0) b->z[i] += L; if (b->z[i] >= L) b->z[i] -= L;
        }
    }
    *fmax = 0;
    for (int i = 0; i < N; i++) {
        double f = sqrt(fx[i] * fx[i] + fy[i] * fy[i] + fz[i] * fz[i]);
        if (f > *fmax) *fmax = f;
    }
    free(fx); free(fy); free(fz);
    return U;
}

/* one conserving energy-transport sweep: paired antisymmetric transfers
 * down the mu gradient. Conservation is exact by construction — the same
 * transfer-ownership discipline the local-clock suite established. */
static double flow_energy(Box *b, double dt)
{
    int N = b->N; double L = b->L;
    double zb, Uc;
    radii(b); contacts(b, &zb, &Uc);
    for (int i = 0; i < N; i++)
        b->mu[i] = (b->E[i] > 1e-12)
                 ? (1.0 / 3.0) * pow(b->E[i], -2.0 / 3.0) * b->P[i] : 0.0;
    double *dE = calloc(N, sizeof(double));
    double mx = 0;
    for (int i = 0; i < N; i++)
        for (int j = i + 1; j < N; j++) {
            double dx = wrap(b->x[i] - b->x[j], L);
            double dy = wrap(b->y[i] - b->y[j], L);
            double dz = wrap(b->z[i] - b->z[j], L);
            double d = sqrt(dx * dx + dy * dy + dz * dz);
            if (d >= b->r[i] + b->r[j]) continue;       /* only contacts carry */
            double q = dt * b->lam * (b->mu[i] - b->mu[j]);
            dE[i] -= q; dE[j] += q;                     /* paired: exact */
            if (fabs(q) > mx) mx = fabs(q);
        }
    for (int i = 0; i < N; i++) {
        b->E[i] += dE[i];
        if (b->E[i] < 1e-9) b->E[i] = 1e-9;             /* no negative energy */
    }
    free(dE);
    return mx;
}

static double totE(Box *b)
{ double s = 0; for (int i = 0; i < b->N; i++) s += b->E[i]; return s; }

static void rstats(Box *b, double *rmin, double *rmax, double *rmean, double *rsd)
{
    int N = b->N; double s = 0, s2 = 0;
    *rmin = 1e300; *rmax = 0;
    for (int i = 0; i < N; i++) {
        double r = b->r[i];
        if (r < *rmin) *rmin = r;
        if (r > *rmax) *rmax = r;
        s += r; s2 += r * r;
    }
    *rmean = s / N;
    *rsd = sqrt(fabs(s2 / N - (*rmean) * (*rmean)));
}

/* alternate mechanical and energetic relaxation to joint equilibrium */
static void equilibrate(Box *b, int cycles, double dt, double *mxq, double *fm)
{
    for (int c = 0; c < cycles; c++) {
        radii(b);
        relax_pos(b, 15, 0.02, fm);
        *mxq = flow_energy(b, dt);
    }
    radii(b);
    relax_pos(b, 15, 0.02, fm);
}

int main(void)
{
    setvbuf(stdout, NULL, _IOLBF, 0);
    printf("# geopress — geometric pressure as the morphology bound\n");
    printf("# r_i = E_i^(1/3); contact repulsion (k/2)delta^2; energy flows\n");
    printf("# down mu = (1/3)E^(-2/3) P, conserved by paired transfers.\n");
    printf("# NO maximum size anywhere in the model.\n\n");

    /* ---------------- G1 ---------------- */
    printf("== G1. BOUNDED SIZE vs COORDINATION (the test consonance failed) ==\n");
    printf("Packing fraction is swept by box size, which sweeps the mean\n");
    printf("contact number z. Consonance gave 5-13%% locked at EVERY degree\n");
    printf("with no trend; a working bound must instead stay sharp as z rises.\n");
    printf("    L      z    r_mean    r_sd/r_mean    r_max/r_mean   dE/E      max|f|\n");
    {
        double Ls[] = {22, 20, 18, 16, 14, 12};
        for (unsigned q = 0; q < sizeof(Ls) / sizeof(Ls[0]); q++) {
            Box *b = box_new(300, Ls[q], 1.0, 0.5);
            box_init(b, 1.0, 20260802);
            double e0 = totE(b), mq = 0, fm = 0;
            equilibrate(b, 90, 0.02, &mq, &fm);
            double zb, Uc; radii(b); contacts(b, &zb, &Uc);
            double rmn, rmx, rmean, rsd;
            rstats(b, &rmn, &rmx, &rmean, &rsd);
            printf("%6.1f %6.2f  %8.4f  %12.4f  %14.4f  %+8.1e  %8.1e\n",
                   Ls[q], zb, rmean, rsd / rmean, rmx / rmean,
                   (totE(b) - e0) / e0, fm);
            box_free(b);
        }
    }
    printf("\nreading: r_sd/r_mean is the width of the size distribution. If it\n");
    printf("stays small and does not blow up as z rises, pressure bounds size\n");
    printf("at every coordination — which is exactly what consonance could not\n");
    printf("do. dE/E must sit at the FP floor: the transport is paired.\n\n");

    /* ---------------- G2 ---------------- */
    printf("== G2. DOES A LOCALISED EXCESS HOLD? (the particle question) ==\n");
    printf("Seed a central sphere with excess energy, then relax. A structure\n");
    printf("that GROWS TO A SIZE AND STOPS — held by the surrounding pressure,\n");
    printf("with no bound in the potential — is the object the program wants.\n");
    printf("  boost   r_core/r_bulk t=0   after relax   retained   E_core frac\n");
    {
        double boosts[] = {2.0, 4.0, 8.0, 16.0, 32.0};
        for (unsigned q = 0; q < 5; q++) {
            Box *b = box_new(300, 20.0, 1.0, 0.5);
            box_init(b, 1.0, 20260802);
            double cx = 10, cy = 10, cz = 10, R = 4.0;
            int ncore = 0;
            for (int i = 0; i < b->N; i++) {
                double dx = wrap(b->x[i] - cx, b->L);
                double dy = wrap(b->y[i] - cy, b->L);
                double dz = wrap(b->z[i] - cz, b->L);
                if (dx * dx + dy * dy + dz * dz < R * R)
                    { b->E[i] *= boosts[q]; ncore++; }
            }
            radii(b);
            double rc0 = 0, rb0 = 0; int nb0 = 0;
            for (int i = 0; i < b->N; i++) {
                double dx = wrap(b->x[i] - cx, b->L);
                double dy = wrap(b->y[i] - cy, b->L);
                double dz = wrap(b->z[i] - cz, b->L);
                if (dx * dx + dy * dy + dz * dz < R * R) rc0 += b->r[i];
                else { rb0 += b->r[i]; nb0++; }
            }
            rc0 /= ncore; rb0 /= nb0;
            double e0 = totE(b), mq = 0, fm = 0;
            equilibrate(b, 120, 0.02, &mq, &fm);
            double rc = 0, rb = 0, ec = 0; int nb2 = 0;
            for (int i = 0; i < b->N; i++) {
                double dx = wrap(b->x[i] - cx, b->L);
                double dy = wrap(b->y[i] - cy, b->L);
                double dz = wrap(b->z[i] - cz, b->L);
                if (dx * dx + dy * dy + dz * dz < R * R)
                    { rc += b->r[i]; ec += b->E[i]; }
                else { rb += b->r[i]; nb2++; }
            }
            rc /= ncore; rb /= nb2;
            printf("%7.1f  %16.4f  %12.4f  %9.3f  %12.4f\n",
                   boosts[q], rc0 / rb0, rc / rb,
                   (rc / rb - 1) / (rc0 / rb0 - 1), ec / totE(b));
            box_free(b);
        }
    }
    printf("\nreading: 'retained' near 0 = the excess dissolved (no particle);\n");
    printf("near 1 = it held its contrast; >1 = runaway. A retained value that\n");
    printf("SATURATES as the boost grows is the bound doing its job: more\n");
    printf("energy stops buying proportionally more size.\n\n");

    /* ---------------- G3 ---------------- */
    printf("== G3. NEGATIVE CONTROL — contact stiffness off ==\n");
    printf("With k=0 there is no pressure, hence no mu, hence nothing to bound\n");
    printf("or localise. Any structure surviving here is not the pressure's.\n");
    printf("      k    r_sd/r_mean   r_max/r_mean   retained(boost 8)\n");
    {
        double ks[] = {0.0, 0.01, 0.1, 1.0};
        for (unsigned q = 0; q < 4; q++) {
            Box *b = box_new(300, 20.0, ks[q], 0.5);
            box_init(b, 1.0, 20260802);
            double cx = 10, cy = 10, cz = 10, R = 4.0;
            int ncore = 0;
            for (int i = 0; i < b->N; i++) {
                double dx = wrap(b->x[i] - cx, b->L);
                double dy = wrap(b->y[i] - cy, b->L);
                double dz = wrap(b->z[i] - cz, b->L);
                if (dx * dx + dy * dy + dz * dz < R * R)
                    { b->E[i] *= 8.0; ncore++; }
            }
            radii(b);
            double rc0 = 0, rb0 = 0; int nb0 = 0;
            for (int i = 0; i < b->N; i++) {
                double dx = wrap(b->x[i] - cx, b->L);
                double dy = wrap(b->y[i] - cy, b->L);
                double dz = wrap(b->z[i] - cz, b->L);
                if (dx * dx + dy * dy + dz * dz < R * R) rc0 += b->r[i];
                else { rb0 += b->r[i]; nb0++; }
            }
            rc0 /= ncore; rb0 /= nb0;
            double mq = 0, fm = 0;
            equilibrate(b, 120, 0.02, &mq, &fm);
            double rc = 0, rb = 0; int nb2 = 0;
            double rmn, rmx, rmean, rsd;
            for (int i = 0; i < b->N; i++) {
                double dx = wrap(b->x[i] - cx, b->L);
                double dy = wrap(b->y[i] - cy, b->L);
                double dz = wrap(b->z[i] - cz, b->L);
                if (dx * dx + dy * dy + dz * dz < R * R) rc += b->r[i];
                else { rb += b->r[i]; nb2++; }
            }
            rc /= ncore; rb /= nb2;
            rstats(b, &rmn, &rmx, &rmean, &rsd);
            printf("%7.2f  %12.4f  %13.4f  %18.3f\n",
                   ks[q], rsd / rmean, rmx / rmean,
                   (rc / rb - 1) / (rc0 / rb0 - 1));
            box_free(b);
        }
    }
    printf("\nreading: k=0 must show no bound and full retention (nothing moves\n");
    printf("at all, since mu is identically zero). The bound and any\n");
    printf("localisation must appear only as k rises.\n");
    return 0;
}
