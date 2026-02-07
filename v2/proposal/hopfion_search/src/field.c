/*
 * field.c — Energy functional and gradient for CHPT soliton search
 *
 * KEY MATHEMATICAL INSIGHT:
 * The scalar extraction <M>_0 in Cl+(3,0,1) kills ALL terms containing e_0.
 * Therefore: <[R_i,R_j]^2>_0 depends ONLY on the quaternion parts of R_i,R_j.
 * Since R_d = Psi~ d_d Psi, and its quaternion part is A_d = q~ d_d q,
 * the Skyrme energy E4 depends ONLY on q = (s, f1, f2, f3).
 * Similarly, E2 = (1/2) sum_d <(d_d Psi)(d_d Psi)~>_0 only depends on q.
 * And EV depends only on |q|^2.
 * Only ED depends on the degenerate sector (j1, j2, j3, p).
 *
 * CONSEQUENCE: The static problem DECOUPLES.
 * - Bulk (s,f1,f2,f3): evolves under E2 + E4 + EV (= standard Skyrme model)
 * - Degenerate (j1,j2,j3,p): evolves under ED, trivially relaxes to zero.
 *
 * We implement all 8 components faithfully (no simplification).
 * The Skyrme force is computed ANALYTICALLY using the chain rule through
 * the quaternion right-currents.
 */

#define _USE_MATH_DEFINES
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <string.h>
#include <omp.h>
#include "field.h"

/* ---------- Quaternion helpers (4-component subset) ---------- */

typedef struct { double s, f1, f2, f3; } Quat;

static inline Quat q_from_mv(Multivector m) {
    return (Quat){m.s, m.f1, m.f2, m.f3};
}

/* Quaternion product in Cl+(3,0): basis {1, e23, e31, e12}
 * e23^2 = e31^2 = e12^2 = -1
 * e23*e31 = -e12,  e31*e12 = -e23,  e12*e23 = -e31
 * (opposite-handed compared to standard i,j,k) */
static inline Quat q_mul(Quat a, Quat b) {
    return (Quat){
        a.s*b.s  - a.f1*b.f1 - a.f2*b.f2 - a.f3*b.f3,
        a.s*b.f1 + a.f1*b.s  - a.f2*b.f3 + a.f3*b.f2,
        a.s*b.f2 + a.f1*b.f3 + a.f2*b.s  - a.f3*b.f1,
        a.s*b.f3 - a.f1*b.f2 + a.f2*b.f1 + a.f3*b.s
    };
}

static inline Quat q_rev(Quat a) {
    return (Quat){a.s, -a.f1, -a.f2, -a.f3};
}

static inline Quat q_add(Quat a, Quat b) {
    return (Quat){a.s+b.s, a.f1+b.f1, a.f2+b.f2, a.f3+b.f3};
}

static inline Quat q_sub(Quat a, Quat b) {
    return (Quat){a.s-b.s, a.f1-b.f1, a.f2-b.f2, a.f3-b.f3};
}

static inline Quat q_scale(double c, Quat a) {
    return (Quat){c*a.s, c*a.f1, c*a.f2, c*a.f3};
}

static inline Quat q_zero(void) {
    return (Quat){0, 0, 0, 0};
}

static inline double q_norm2(Quat a) {
    return a.s*a.s + a.f1*a.f1 + a.f2*a.f2 + a.f3*a.f3;
}

/* Clifford scalar product for quaternions: <AB>_0 = a.s*b.s - a.f1*b.f1 - ... */
static inline double q_scalar_prod(Quat a, Quat b) {
    return a.s*b.s - a.f1*b.f1 - a.f2*b.f2 - a.f3*b.f3;
}

/* Left multiplication by basis quaternions:
 * eps[0] = 1, eps[1] = e23, eps[2] = e31, eps[3] = e12 */
static inline Quat q_left_basis(int a, Quat q) {
    switch(a) {
    case 0: return q;
    case 1: return (Quat){-q.f1,  q.s,   q.f3, -q.f2}; /* e23 * q */
    case 2: return (Quat){-q.f2, -q.f3,  q.s,   q.f1}; /* e31 * q */
    case 3: return (Quat){-q.f3,  q.f2, -q.f1,  q.s};  /* e12 * q */
    default: return q_zero();
    }
}

/* ---------- Grid operations ---------- */

Field *field_alloc(int N, double L)
{
    Field *f = (Field *)malloc(sizeof(Field));
    if (!f) { fprintf(stderr, "field_alloc: malloc failed\n"); exit(1); }
    f->N = N;
    f->L = L;
    f->h = 2.0 * L / N;
    f->psi = (Multivector *)calloc((size_t)N * N * N, sizeof(Multivector));
    if (!f->psi) { fprintf(stderr, "field_alloc: calloc failed\n"); exit(1); }
    return f;
}

void field_free(Field *f)
{
    if (f) { free(f->psi); free(f); }
}

/* 4th-order central difference for quaternion part */
static inline Quat q_deriv(const Field *f, int i, int j, int k, int dir)
{
    int N = f->N;
    double inv12h = 1.0 / (12.0 * f->h);
    Multivector pm2, pm1, pp1, pp2;

    switch (dir) {
    case 0:
        pm2 = f->psi[idx(N,i-2,j,k)]; pm1 = f->psi[idx(N,i-1,j,k)];
        pp1 = f->psi[idx(N,i+1,j,k)]; pp2 = f->psi[idx(N,i+2,j,k)];
        break;
    case 1:
        pm2 = f->psi[idx(N,i,j-2,k)]; pm1 = f->psi[idx(N,i,j-1,k)];
        pp1 = f->psi[idx(N,i,j+1,k)]; pp2 = f->psi[idx(N,i,j+2,k)];
        break;
    default:
        pm2 = f->psi[idx(N,i,j,k-2)]; pm1 = f->psi[idx(N,i,j,k-1)];
        pp1 = f->psi[idx(N,i,j,k+1)]; pp2 = f->psi[idx(N,i,j,k+2)];
        break;
    }

    Quat r;
    r.s  = (-pp2.s  + 8*pp1.s  - 8*pm1.s  + pm2.s)  * inv12h;
    r.f1 = (-pp2.f1 + 8*pp1.f1 - 8*pm1.f1 + pm2.f1) * inv12h;
    r.f2 = (-pp2.f2 + 8*pp1.f2 - 8*pm1.f2 + pm2.f2) * inv12h;
    r.f3 = (-pp2.f3 + 8*pp1.f3 - 8*pm1.f3 + pm2.f3) * inv12h;
    return r;
}

/* Full 8-component derivative (for degenerate sector gradient energy, if needed) */
static inline Multivector mv_deriv(const Field *f, int i, int j, int k, int dir)
{
    int N = f->N;
    double inv12h = 1.0 / (12.0 * f->h);
    Multivector pm2, pm1, pp1, pp2;

    switch (dir) {
    case 0:
        pm2 = f->psi[idx(N,i-2,j,k)]; pm1 = f->psi[idx(N,i-1,j,k)];
        pp1 = f->psi[idx(N,i+1,j,k)]; pp2 = f->psi[idx(N,i+2,j,k)];
        break;
    case 1:
        pm2 = f->psi[idx(N,i,j-2,k)]; pm1 = f->psi[idx(N,i,j-1,k)];
        pp1 = f->psi[idx(N,i,j+1,k)]; pp2 = f->psi[idx(N,i,j+2,k)];
        break;
    default:
        pm2 = f->psi[idx(N,i,j,k-2)]; pm1 = f->psi[idx(N,i,j,k-1)];
        pp1 = f->psi[idx(N,i,j,k+1)]; pp2 = f->psi[idx(N,i,j,k+2)];
        break;
    }

    Multivector r;
    r.s  = (-pp2.s  + 8*pp1.s  - 8*pm1.s  + pm2.s)  * inv12h;
    r.f1 = (-pp2.f1 + 8*pp1.f1 - 8*pm1.f1 + pm2.f1) * inv12h;
    r.f2 = (-pp2.f2 + 8*pp1.f2 - 8*pm1.f2 + pm2.f2) * inv12h;
    r.f3 = (-pp2.f3 + 8*pp1.f3 - 8*pm1.f3 + pm2.f3) * inv12h;
    r.j1 = (-pp2.j1 + 8*pp1.j1 - 8*pm1.j1 + pm2.j1) * inv12h;
    r.j2 = (-pp2.j2 + 8*pp1.j2 - 8*pm1.j2 + pm2.j2) * inv12h;
    r.j3 = (-pp2.j3 + 8*pp1.j3 - 8*pm1.j3 + pm2.j3) * inv12h;
    r.p  = (-pp2.p  + 8*pp1.p  - 8*pm1.p  + pm2.p)  * inv12h;
    return r;
}

/* Consistent Laplacian for quaternion part: Σ_d D_d(D_d q)
 *
 * The E2 energy uses 4th-order central differences D_d for first derivatives.
 * The exact discrete gradient of E2 w.r.t. q(x) is Σ_d D_d(D_d q)(x), NOT
 * the standard 7-point Laplacian. These differ because D_d composed with
 * itself gives a 9-point stencil per direction:
 *
 * weights at offsets -4..+4: {1, -16, 64, 16, -130, 16, 64, -16, 1} / (144 h²)
 *
 * Derived from convolving the 4th-order stencil {1,-8,0,8,-1}/(12h) with itself.
 */
static inline Quat consistent_lap_q(const Field *f, int i, int j, int k)
{
    int N = f->N;
    double c = 1.0 / (144.0 * f->h * f->h);
    static const double w[9] = {1, -16, 64, 16, -130, 16, 64, -16, 1};

    double rs = 0, rf1 = 0, rf2 = 0, rf3 = 0;

    /* Direction 0 (x) */
    for (int m = -4; m <= 4; m++) {
        Multivector p = f->psi[idx(N, i+m, j, k)];
        double wm = w[m+4];
        rs += wm * p.s; rf1 += wm * p.f1; rf2 += wm * p.f2; rf3 += wm * p.f3;
    }
    /* Direction 1 (y) */
    for (int m = -4; m <= 4; m++) {
        Multivector p = f->psi[idx(N, i, j+m, k)];
        double wm = w[m+4];
        rs += wm * p.s; rf1 += wm * p.f1; rf2 += wm * p.f2; rf3 += wm * p.f3;
    }
    /* Direction 2 (z) */
    for (int m = -4; m <= 4; m++) {
        Multivector p = f->psi[idx(N, i, j, k+m)];
        double wm = w[m+4];
        rs += wm * p.s; rf1 += wm * p.f1; rf2 += wm * p.f2; rf3 += wm * p.f3;
    }

    return (Quat){c * rs, c * rf1, c * rf2, c * rf3};
}

/* ---------- Energy computation ---------- */

Energy field_energy(const Field *f, const Params *p)
{
    int N = f->N;
    double h3 = f->h * f->h * f->h;
    double e2_sum = 0, e4_sum = 0, ev_sum = 0, ed_sum = 0;
    double inv_4e2 = 1.0 / (4.0 * p->e_skyrme * p->e_skyrme);

    #pragma omp parallel for collapse(3) reduction(+:e2_sum,e4_sum,ev_sum,ed_sum) schedule(static)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        int ix = idx(N, i, j, k);
        Quat q = q_from_mv(f->psi[ix]);
        Quat qr = q_rev(q);

        /* Spatial derivatives */
        Quat dq[3];
        for (int d = 0; d < 3; d++)
            dq[d] = q_deriv(f, i, j, k, d);

        /* E2: (1/2) sum_d |d_d q|^2  (quaternion norm = Clifford bulk norm) */
        for (int d = 0; d < 3; d++)
            e2_sum += 0.5 * q_norm2(dq[d]);

        /* Right-currents A_d = q~ * d_d q (quaternion product) */
        Quat A[3];
        for (int d = 0; d < 3; d++)
            A[d] = q_mul(qr, dq[d]);

        /* E4: -(1/4e^2) sum_{d1<d2} <[A_d1, A_d2]^2>_0
         * [A,B]^2 is a quaternion whose scalar part is -(|commutator|^2)
         * So <[A,B]^2>_0 < 0, and E4 = -(1/4e^2)*negative > 0. */
        for (int d1 = 0; d1 < 3; d1++)
        for (int d2 = d1+1; d2 < 3; d2++) {
            Quat comm = q_sub(q_mul(A[d1], A[d2]), q_mul(A[d2], A[d1]));
            Quat comm2 = q_mul(comm, comm);
            e4_sum -= inv_4e2 * comm2.s;
        }

        /* EV: (lambda/4)(|q|^2 - rho0^2)^2 */
        double dev = q_norm2(q) - p->rho0 * p->rho0;
        ev_sum += 0.25 * p->lambda * dev * dev;

        /* ED: (mu^2/2)(j1^2 + j2^2 + j3^2 + P^2) */
        ed_sum += 0.5 * p->mu * p->mu * mv_weight_norm2(f->psi[ix]);
    }

    Energy en;
    en.E2 = e2_sum * h3;
    en.E4 = e4_sum * h3;
    en.EV = ev_sum * h3;
    en.ED = ed_sum * h3;
    en.Etotal = en.E2 + en.E4 + en.EV + en.ED;
    return en;
}

/* ---------- Skyrme force (analytical) ----------
 *
 * The Skyrme force on bulk component a of Psi at point x is:
 *
 *   F4_a(x) = (1/2e^2) sum_d { sigma_a * <eps_a * dq_d(x) , G_d(x)>_0
 *                               - D_d( <q~(x) eps_a, G_d(x)>_0 ) }
 *
 * where:
 *   A_d = q~ * d_d q           (quaternion right-current)
 *   C_{d1,d2} = [A_d1, A_d2]   (commutator)
 *   G_d = sum_{d'!=d} [A_{d'}, C_{d,d'}]  ("Skyrme tensor")
 *   eps_a = a-th basis quaternion (1, e23, e31, e12)
 *   sigma_a = sign under reversion (+1 for a=0, -1 for a=1,2,3)
 *   <X,Y>_0 = quaternion Clifford scalar product
 *   D_d = spatial derivative (finite diff)
 *
 * DERIVATION: From delta E4 = -(1/2e^2) sum_{x,d} <delta_R_d, G_d>_0
 * with delta_R_d = (delta q~)(d_d q) + q~ d_d(delta q).
 * Integration by parts on the second term gives the two contributions above.
 */

/* Precomputed data at each grid point */
typedef struct {
    Quat A[3];     /* right-currents */
    Quat G[3];     /* Skyrme G-tensor */
} SkyrmePre;

static void precompute_skyrme(const Field *f, const Params *p __attribute__((unused)), SkyrmePre *pre)
{
    int N = f->N;

    #pragma omp parallel for collapse(3) schedule(static)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        int ix = idx(N, i, j, k);
        Quat q = q_from_mv(f->psi[ix]);
        Quat qr = q_rev(q);

        /* Right-currents */
        for (int d = 0; d < 3; d++) {
            Quat dq = q_deriv(f, i, j, k, d);
            pre[ix].A[d] = q_mul(qr, dq);
        }

        /* Commutators C_{d1,d2} and G-tensor */
        Quat C[3]; /* C[0] = C_{01}, C[1] = C_{02}, C[2] = C_{12} */
        C[0] = q_sub(q_mul(pre[ix].A[0], pre[ix].A[1]),
                      q_mul(pre[ix].A[1], pre[ix].A[0]));
        C[1] = q_sub(q_mul(pre[ix].A[0], pre[ix].A[2]),
                      q_mul(pre[ix].A[2], pre[ix].A[0]));
        C[2] = q_sub(q_mul(pre[ix].A[1], pre[ix].A[2]),
                      q_mul(pre[ix].A[2], pre[ix].A[1]));

        /* G_d = sum_{d'!=d} [A_{d'}, C_{d,d'}]
         * C_{d,d'} with d < d': use C array
         * C_{d,d'} with d > d': = -C_{d',d}
         */
        for (int d = 0; d < 3; d++) {
            pre[ix].G[d] = q_zero();
            for (int dp = 0; dp < 3; dp++) {
                if (dp == d) continue;
                /* Get C_{d,dp} */
                Quat Cdd;
                if (d < dp) {
                    int ci = (d==0 && dp==1) ? 0 : (d==0 && dp==2) ? 1 : 2;
                    Cdd = C[ci];
                } else {
                    int ci = (dp==0 && d==1) ? 0 : (dp==0 && d==2) ? 1 : 2;
                    Cdd = q_scale(-1.0, C[ci]);
                }
                /* [A_{dp}, C_{d,dp}] */
                Quat term = q_sub(q_mul(pre[ix].A[dp], Cdd),
                                  q_mul(Cdd, pre[ix].A[dp]));
                pre[ix].G[d] = q_add(pre[ix].G[d], term);
            }
        }
    }
}

void field_gradient(const Field *f, const Params *p, Multivector *force)
{
    int N = f->N;
    int N3 = N * N * N;
    double inv_2e2 = 1.0 / (2.0 * p->e_skyrme * p->e_skyrme);

    memset(force, 0, (size_t)N3 * sizeof(Multivector));

    /* --- Step 1: Precompute Skyrme data --- */
    SkyrmePre *pre = (SkyrmePre *)malloc((size_t)N3 * sizeof(SkyrmePre));
    if (!pre) { fprintf(stderr, "field_gradient: malloc failed\n"); exit(1); }
    precompute_skyrme(f, p, pre);

    /* --- Step 2: Precompute pi_{a,d}(x) = <q~(x) * eps_a, G_d(x)>_0 --- */
    /* pi has 4 components (a=0..3) × 3 directions (d=0..2) per point */
    double *pi = (double *)calloc((size_t)N3 * 12, sizeof(double));
    if (!pi) { fprintf(stderr, "field_gradient: malloc failed\n"); exit(1); }

    #pragma omp parallel for collapse(3) schedule(static)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        int ix = idx(N, i, j, k);
        Quat qr = q_rev(q_from_mv(f->psi[ix]));
        for (int d = 0; d < 3; d++) {
            for (int a = 0; a < 4; a++) {
                /* q~ * eps_a */
                /* pi_{a,d} = scalar part of (q~ * eps_a * G_d)
                 * q~ * eps_a: right multiplication of q~ by basis element a */
                Quat qr_ea_r;
                switch (a) {
                case 0: qr_ea_r = qr; break;
                case 1: qr_ea_r = q_mul(qr, (Quat){0,1,0,0}); break;
                case 2: qr_ea_r = q_mul(qr, (Quat){0,0,1,0}); break;
                case 3: qr_ea_r = q_mul(qr, (Quat){0,0,0,1}); break;
                default: qr_ea_r = q_zero();
                }

                /* Scalar part of (qr_ea_r * G_d) */
                Quat prod = q_mul(qr_ea_r, pre[ix].G[d]);
                pi[12*ix + 3*a + d] = prod.s;
            }
        }
    }

    /* --- Step 3: Compute forces --- */
    double sigma[4] = {1.0, -1.0, -1.0, -1.0}; /* reversion signs */

    #pragma omp parallel for collapse(3) schedule(static)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        int ix = idx(N, i, j, k);
        Multivector psi_here = f->psi[ix];
        Quat q = q_from_mv(psi_here);

        /* --- E2 force: consistent Laplacian (only bulk) --- */
        Quat lap = consistent_lap_q(f, i, j, k);
        force[ix].s  = lap.s;
        force[ix].f1 = lap.f1;
        force[ix].f2 = lap.f2;
        force[ix].f3 = lap.f3;

        /* --- EV force: -lambda * (|q|^2 - rho0^2) * q --- */
        double dev = q_norm2(q) - p->rho0 * p->rho0;
        force[ix].s  -= p->lambda * dev * q.s;
        force[ix].f1 -= p->lambda * dev * q.f1;
        force[ix].f2 -= p->lambda * dev * q.f2;
        force[ix].f3 -= p->lambda * dev * q.f3;

        /* --- ED force: -mu^2 * (J, P) --- */
        force[ix].j1 = -p->mu * p->mu * psi_here.j1;
        force[ix].j2 = -p->mu * p->mu * psi_here.j2;
        force[ix].j3 = -p->mu * p->mu * psi_here.j3;
        force[ix].p  = -p->mu * p->mu * psi_here.p;

        /* --- E4 force (Skyrme, bulk only) --- */
        /* Term 1: sigma_a * <eps_a * d_d q, G_d>_0 summed over d */
        double f4[4] = {0, 0, 0, 0};
        for (int d = 0; d < 3; d++) {
            Quat dq = q_deriv(f, i, j, k, d);
            for (int a = 0; a < 4; a++) {
                Quat ea_dq = q_left_basis(a, dq);    /* eps_a * d_d q */
                Quat prod = q_mul(ea_dq, pre[ix].G[d]);
                f4[a] += sigma[a] * prod.s;           /* <eps_a*dq, G_d>_0 */
            }
        }

        /* Term 2: -D_d(pi_{a,d}) summed over d */
        /* Using 4th-order central difference on pi */
        double inv12h = 1.0 / (12.0 * f->h);
        for (int a = 0; a < 4; a++) {
            for (int d = 0; d < 3; d++) {
                /* pi at neighbors in direction d */
                int im2, im1, ip1, ip2;
                switch (d) {
                case 0:
                    im2 = idx(N,i-2,j,k); im1 = idx(N,i-1,j,k);
                    ip1 = idx(N,i+1,j,k); ip2 = idx(N,i+2,j,k);
                    break;
                case 1:
                    im2 = idx(N,i,j-2,k); im1 = idx(N,i,j-1,k);
                    ip1 = idx(N,i,j+1,k); ip2 = idx(N,i,j+2,k);
                    break;
                default:
                    im2 = idx(N,i,j,k-2); im1 = idx(N,i,j,k-1);
                    ip1 = idx(N,i,j,k+1); ip2 = idx(N,i,j,k+2);
                    break;
                }
                double pm2 = pi[12*im2 + 3*a + d];
                double pm1 = pi[12*im1 + 3*a + d];
                double pp1 = pi[12*ip1 + 3*a + d];
                double pp2 = pi[12*ip2 + 3*a + d];
                double dpi = (-pp2 + 8*pp1 - 8*pm1 + pm2) * inv12h;
                f4[a] -= dpi;
            }
        }

        /* Apply Skyrme force with prefactor */
        force[ix].s  += inv_2e2 * f4[0];
        force[ix].f1 += inv_2e2 * f4[1];
        force[ix].f2 += inv_2e2 * f4[2];
        force[ix].f3 += inv_2e2 * f4[3];
    }

    free(pi);
    free(pre);
}

/* ---------- Topological charge ---------- */

double field_topological_charge(const Field *f, const Params *p)
{
    int N = f->N;
    double h3 = f->h * f->h * f->h;
    double sum = 0.0;

    #pragma omp parallel for collapse(3) reduction(+:sum) schedule(static)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        int ix = idx(N, i, j, k);
        Quat qr = q_rev(q_from_mv(f->psi[ix]));

        /* Right-currents */
        Quat dq[3];
        for (int d = 0; d < 3; d++)
            dq[d] = q_deriv(f, i, j, k, d);

        Quat A[3];
        for (int d = 0; d < 3; d++)
            A[d] = q_mul(qr, dq[d]);

        /* Topological charge density: (1/(2 pi^2)) <A_0 A_1 A_2>_0 / rho0^4
         * The baryon number / degree of map. */
        Quat A01 = q_mul(A[0], A[1]);
        Quat A012 = q_mul(A01, A[2]);
        sum += A012.s;
    }

    /* B = -(1/(2 pi^2 rho0^4)) * sum * h^3 */
    double rho04 = p->rho0 * p->rho0 * p->rho0 * p->rho0;
    return -sum * h3 / (2.0 * M_PI * M_PI * rho04);
}
