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
    double L, k, lam, cap;
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
    b->N = N; b->L = L; b->k = k; b->lam = lam; b->cap = 1e12;
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

/* CAP SATURATION — the non-convexity.
 *
 * G2 showed a purely repulsive contact energy cannot localise: it is
 * CONVEX, a convex energy has one minimising phase, and dispersal is the
 * only behaviour available. Localisation needs mu to FALL with density
 * over some range so two phases can coexist.
 *
 * The v89-native way to get that is `cap`'s own mechanism, load flattens
 * pitch, read mechanically: a heavily loaded cell goes SOFT. A pair's
 * stiffness falls with the load the pair carries,
 *
 *     k_ij = k / (1 + (E_i + E_j) / (2 cap))
 *
 * which is symmetric, so contact forces stay paired and conservation is
 * untouched. It adds a NEGATIVE term to mu (softening lowers the cost of
 * holding more energy), and that term is the van der Waals loop. cap ->
 * infinity recovers the unsaturated model of G1/G2 exactly. */
static double kij(Box *b, int i, int j)
{ return b->k / (1.0 + (b->E[i] + b->E[j]) / (2.0 * b->cap)); }

/* d k_ij / d E_i  (negative: load softens) */
static double dkij_dE(Box *b, int i, int j)
{
    double u = 1.0 + (b->E[i] + b->E[j]) / (2.0 * b->cap);
    return -b->k / (2.0 * b->cap * u * u);
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
            double del = s - d, kk = kij(b, i, j);
            b->P[i] += kk * del; b->P[j] += kk * del;
            U += 0.5 * kk * del * del;
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
                double kk = kij(b, i, j);
                double del = s - d, f = kk * del / d;
                fx[i] += f * dx; fy[i] += f * dy; fz[i] += f * dz;
                fx[j] -= f * dx; fy[j] -= f * dy; fz[j] -= f * dz;
                U += 0.5 * kk * del * del;
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
    /* mu = dU/dE has TWO parts now:
     *   (a) growing costs overlap:      (1/3)E^(-2/3) * P    (positive)
     *   (b) load softens the contacts:  sum_j (delta^2/2) dk/dE  (negative)
     * (b) is the non-convexity. Where it wins, a denser cell has LOWER mu
     * and therefore draws energy instead of shedding it — coexistence. */
    double *soft = calloc(N, sizeof(double));
    for (int i = 0; i < N; i++)
        for (int j = i + 1; j < N; j++) {
            double dx = wrap(b->x[i] - b->x[j], L);
            double dy = wrap(b->y[i] - b->y[j], L);
            double dz = wrap(b->z[i] - b->z[j], L);
            double d = sqrt(dx * dx + dy * dy + dz * dz);
            double sm = b->r[i] + b->r[j];
            if (d >= sm) continue;
            double del = sm - d, t = 0.5 * del * del * dkij_dE(b, i, j);
            soft[i] += t; soft[j] += t;
        }
    for (int i = 0; i < N; i++)
        b->mu[i] = (b->E[i] > 1e-12)
                 ? (1.0 / 3.0) * pow(b->E[i], -2.0 / 3.0) * b->P[i] + soft[i]
                 : 0.0;
    free(soft);
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
    printf("localisation must appear only as k rises.\n\n");

    /* ---------------- G4: N16 ---------------- */
    printf("== G4 (N16). DOES cap SATURATION LET A CORE HOLD? ==\n");
    printf("k_ij = k/(1+(E_i+E_j)/2cap): a loaded pair goes soft, which puts a\n");
    printf("NEGATIVE term in mu. cap=1e12 is the unsaturated model and must\n");
    printf("reproduce G2 (retention about -0.25 at boost 8). If saturation\n");
    printf("opens a coexistence region, retention should RISE as cap falls.\n");
    printf("      cap   retained(boost 8)   r_sd/r_mean   r_max/r_mean   dE/E\n");
    {
        double caps[] = {1e12, 100, 30, 10, 3, 1, 0.3, 0.1};
        for (unsigned q = 0; q < sizeof(caps) / sizeof(caps[0]); q++) {
            Box *b = box_new(300, 20.0, 1.0, 0.5);
            box_init(b, 1.0, 20260802);
            b->cap = caps[q];
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
            double e0 = totE(b), mq = 0, fm = 0;
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
            printf("%9.3g  %17.3f  %13.4f  %14.4f  %+8.1e\n",
                   caps[q], (rc / rb - 1) / (rc0 / rb0 - 1),
                   rsd / rmean, rmx / rmean, (totE(b) - e0) / e0);
            box_free(b);
        }
    }
    /* ---------------- G5: the kinetic control ---------------- */
    printf("\n== G5. IS THE HELD CORE REAL, OR JUST SLOW? ==\n");
    printf("Low cap means SOFT contacts, hence slower relaxation, so retention\n");
    printf("could be incomplete dispersal rather than a held state. If the core\n");
    printf("is a true coexisting phase, retention is FLAT in relaxation time.\n");
    printf("If it decays toward the unsaturated -0.25, it was kinetics.\n");
    printf("     cap   cycles   retained   max|f|      max|dE per sweep|\n");
    {
        double caps[] = {0.1, 0.3, 1e12};
        int cyc[] = {120, 240, 480, 960};
        for (unsigned c = 0; c < 3; c++)
        for (unsigned q = 0; q < 4; q++) {
            Box *b = box_new(300, 20.0, 1.0, 0.5);
            box_init(b, 1.0, 20260802);
            b->cap = caps[c];
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
            equilibrate(b, cyc[q], 0.02, &mq, &fm);
            double rc = 0, rb = 0; int nb2 = 0;
            for (int i = 0; i < b->N; i++) {
                double dx = wrap(b->x[i] - cx, b->L);
                double dy = wrap(b->y[i] - cy, b->L);
                double dz = wrap(b->z[i] - cz, b->L);
                if (dx * dx + dy * dy + dz * dz < R * R) rc += b->r[i];
                else { rb += b->r[i]; nb2++; }
            }
            rc /= ncore; rb /= nb2;
            printf("%8.3g  %7d  %9.3f  %9.1e  %14.1e\n",
                   caps[c], cyc[q], (rc / rb - 1) / (rc0 / rb0 - 1), fm, mq);
            box_free(b);
        }
    }
    printf("\nreading: FLAT retention across a 8x range of relaxation time, with\n");
    printf("max|dE per sweep| falling toward zero, is a genuine coexisting\n");
    printf("phase. Retention that keeps sliding toward -0.25 is the soft-contact\n");
    printf("kinetics of low cap and nothing more.\n");

    printf("\nreading: retention rising above 0 as cap falls = the core HOLDS,\n");
    printf("and a value that saturates rather than growing without limit is a\n");
    printf("bounded object. Watch r_max/r_mean at the same time: if the bound\n");
    printf("blows up while retention rises, saturation bought localisation by\n");
    printf("destroying the bound, which is not a win. dE/E must stay at the\n");
    printf("FP floor throughout — the softening is symmetric, so transfers\n");
    printf("stay paired.\n");
    return 0;
}
