/*
 * coupling.c — Degenerate sector coupling energy and forces
 *
 * Implements three coupling terms between bulk q=(s,f1,f2,f3) and
 * degenerate w=(j1,j2,j3,p) sectors of Cl+(3,0,1):
 *
 * 1. E_{2,D} = (1/2) sum_d |D_d w|^2
 * 2. E_{4,C} = (1/4e^2) sum_{i<j} |F^w_{ij}|^2
 *    where F^w_{ij} = [A_i, W_j] - [A_j, W_i] is the weight part
 *    of the commutator of full 8D right-currents.
 * 3. E_int = (g^2/2) |q|^2 |D_d w|^2
 *
 * For E_{4,C} gradient:
 * We use the direct variational approach. The energy density at point x is
 *   e4c(x) = (1/4e^2) sum_{i<j} |F^w_{ij}(x)|^2
 * where F^w depends on {q, D_d q, w, D_d w} at x. The force at x collects:
 *   (a) Local terms from variation of q(x) and w(x) through A, W at x
 *   (b) Non-local terms from D_d q(x) and D_d w(x) appearing in A, W at
 *       neighboring points — handled via pi-field divergence.
 */

#define _USE_MATH_DEFINES
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include "coupling.h"

/* ---------- Weight-sector quaternion type ---------- */

typedef struct { double j1, j2, j3, p; } Wquat;
typedef struct { double s, f1, f2, f3; } Quat;

static inline Quat q_from_mv(Multivector m) {
    return (Quat){m.s, m.f1, m.f2, m.f3};
}

static inline Wquat w_from_mv(Multivector m) {
    return (Wquat){m.j1, m.j2, m.j3, m.p};
}

static inline Quat q_rev(Quat a) {
    return (Quat){a.s, -a.f1, -a.f2, -a.f3};
}

static inline Quat q_mul(Quat a, Quat b) {
    return (Quat){
        a.s*b.s  - a.f1*b.f1 - a.f2*b.f2 - a.f3*b.f3,
        a.s*b.f1 + a.f1*b.s  - a.f2*b.f3 + a.f3*b.f2,
        a.s*b.f2 + a.f1*b.f3 + a.f2*b.s  - a.f3*b.f1,
        a.s*b.f3 - a.f1*b.f2 + a.f2*b.f1 + a.f3*b.s
    };
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

static inline Quat q_zero(void) { return (Quat){0,0,0,0}; }

static inline double q_norm2(Quat a) {
    return a.s*a.s + a.f1*a.f1 + a.f2*a.f2 + a.f3*a.f3;
}

/* Weight quaternion operations */
static inline Wquat w_zero(void) { return (Wquat){0,0,0,0}; }

static inline Wquat w_add(Wquat a, Wquat b) {
    return (Wquat){a.j1+b.j1, a.j2+b.j2, a.j3+b.j3, a.p+b.p};
}

static inline Wquat w_sub(Wquat a, Wquat b) {
    return (Wquat){a.j1-b.j1, a.j2-b.j2, a.j3-b.j3, a.p-b.p};
}

static inline Wquat w_scale(double c, Wquat a) {
    return (Wquat){c*a.j1, c*a.j2, c*a.j3, c*a.p};
}

static inline double w_norm2(Wquat a) {
    return a.j1*a.j1 + a.j2*a.j2 + a.j3*a.j3 + a.p*a.p;
}

static inline double w_dot(Wquat a, Wquat b) {
    return a.j1*b.j1 + a.j2*b.j2 + a.j3*b.j3 + a.p*b.p;
}

/* Weight reversion: (j1,j2,j3,p) -> (-j1,-j2,-j3,+p) */
static inline Wquat w_rev(Wquat a) {
    return (Wquat){-a.j1, -a.j2, -a.j3, a.p};
}

/* ---------- Mixed products ---------- */

/* qw_mul: q * w where q in Cl+(3,0), w in Cl-(3,0) -> Cl-(3,0) */
static inline Wquat qw_mul(Quat q, Wquat w) {
    return (Wquat){
        q.s*w.j1 - q.f1*w.p  - q.f2*w.j3 + q.f3*w.j2,
        q.s*w.j2 + q.f1*w.j3 - q.f2*w.p  - q.f3*w.j1,
        q.s*w.j3 - q.f1*w.j2 + q.f2*w.j1 - q.f3*w.p,
        q.s*w.p  + q.f1*w.j1 + q.f2*w.j2 + q.f3*w.j3
    };
}

/* wq_mul: w * q where w in Cl-(3,0), q in Cl+(3,0) -> Cl-(3,0) */
static inline Wquat wq_mul(Wquat w, Quat q) {
    return (Wquat){
        w.j1*q.s  - w.p*q.f1  + w.j3*q.f2 - w.j2*q.f3,
        w.j1*q.f3 + w.j2*q.s  - w.j3*q.f1 - w.p*q.f2,
       -w.j1*q.f2 + w.j2*q.f1 + w.j3*q.s  - w.p*q.f3,
        w.j1*q.f1 + w.j2*q.f2 + w.j3*q.f3 + w.p*q.s
    };
}

/* Commutator [A, W] = AW - WA where A bulk, W weight -> weight */
static inline Wquat qw_commutator(Quat a, Wquat w) {
    return w_sub(qw_mul(a, w), wq_mul(w, a));
}

/* Left multiplication by quaternion basis elements:
 * basis_mul(a, q) = eps_a * q */
static inline Quat q_left_basis(int a, Quat q) {
    switch(a) {
    case 0: return q;
    case 1: return (Quat){-q.f1,  q.s,   q.f3, -q.f2};
    case 2: return (Quat){-q.f2, -q.f3,  q.s,   q.f1};
    case 3: return (Quat){-q.f3,  q.f2, -q.f1,  q.s};
    default: return q_zero();
    }
}

/* ---------- Grid derivative helpers ---------- */

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

    return (Quat){
        (-pp2.s  + 8*pp1.s  - 8*pm1.s  + pm2.s)  * inv12h,
        (-pp2.f1 + 8*pp1.f1 - 8*pm1.f1 + pm2.f1) * inv12h,
        (-pp2.f2 + 8*pp1.f2 - 8*pm1.f2 + pm2.f2) * inv12h,
        (-pp2.f3 + 8*pp1.f3 - 8*pm1.f3 + pm2.f3) * inv12h
    };
}

static inline Wquat w_deriv(const Field *f, int i, int j, int k, int dir)
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

    return (Wquat){
        (-pp2.j1 + 8*pp1.j1 - 8*pm1.j1 + pm2.j1) * inv12h,
        (-pp2.j2 + 8*pp1.j2 - 8*pm1.j2 + pm2.j2) * inv12h,
        (-pp2.j3 + 8*pp1.j3 - 8*pm1.j3 + pm2.j3) * inv12h,
        (-pp2.p  + 8*pp1.p  - 8*pm1.p  + pm2.p)   * inv12h
    };
}

/* Consistent Laplacian for weight part */
static inline Wquat consistent_lap_w(const Field *f, int i, int j, int k)
{
    int N = f->N;
    double c = 1.0 / (144.0 * f->h * f->h);
    static const double wt[9] = {1, -16, 64, 16, -130, 16, 64, -16, 1};

    double rj1 = 0, rj2 = 0, rj3 = 0, rp = 0;

    for (int m = -4; m <= 4; m++) {
        double wm = wt[m+4];
        Multivector px = f->psi[idx(N, i+m, j, k)];
        rj1 += wm * px.j1; rj2 += wm * px.j2; rj3 += wm * px.j3; rp += wm * px.p;

        Multivector py = f->psi[idx(N, i, j+m, k)];
        rj1 += wm * py.j1; rj2 += wm * py.j2; rj3 += wm * py.j3; rp += wm * py.p;

        Multivector pz = f->psi[idx(N, i, j, k+m)];
        rj1 += wm * pz.j1; rj2 += wm * pz.j2; rj3 += wm * pz.j3; rp += wm * pz.p;
    }

    return (Wquat){c * rj1, c * rj2, c * rj3, c * rp};
}

/* ---------- Energy computation ---------- */

CouplingEnergy coupling_energy(const Field *f, const Params *p, double g_coupling)
{
    int N = f->N;
    double h3 = f->h * f->h * f->h;
    double inv_4e2 = 1.0 / (4.0 * p->e_skyrme * p->e_skyrme);
    double g2_half = 0.5 * g_coupling * g_coupling;

    double e2d_sum = 0, e4c_sum = 0, eint_sum = 0;

    #pragma omp parallel for collapse(3) reduction(+:e2d_sum,e4c_sum,eint_sum) schedule(static)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        int ix = idx(N, i, j, k);
        Quat q = q_from_mv(f->psi[ix]);
        Quat qr = q_rev(q);
        Wquat w = w_from_mv(f->psi[ix]);
        Wquat wr = w_rev(w);

        Quat dq[3];
        Wquat dw[3];
        for (int d = 0; d < 3; d++) {
            dq[d] = q_deriv(f, i, j, k, d);
            dw[d] = w_deriv(f, i, j, k, d);
        }

        /* E_{2,D} */
        double grad_w2 = 0;
        for (int d = 0; d < 3; d++)
            grad_w2 += w_norm2(dw[d]);
        e2d_sum += 0.5 * grad_w2;

        /* E_int */
        eint_sum += g2_half * q_norm2(q) * grad_w2;

        /* E_{4,C}: right-currents and coupling commutators */
        Quat A[3];
        Wquat W[3];
        for (int d = 0; d < 3; d++) {
            A[d] = q_mul(qr, dq[d]);
            W[d] = w_add(qw_mul(qr, dw[d]), wq_mul(wr, dq[d]));
        }

        for (int d1 = 0; d1 < 3; d1++)
        for (int d2 = d1+1; d2 < 3; d2++) {
            Wquat comm = w_sub(qw_commutator(A[d1], W[d2]),
                               qw_commutator(A[d2], W[d1]));
            e4c_sum += inv_4e2 * w_norm2(comm);
        }
    }

    CouplingEnergy en;
    en.E2D = e2d_sum * h3;
    en.E4C = e4c_sum * h3;
    en.Eint = eint_sum * h3;
    en.Etotal = en.E2D + en.E4C + en.Eint;
    return en;
}

/* ---------- Gradient (force) computation ----------
 *
 * Strategy: For E_{4,C}, the energy density at point x is
 *   e(x) = (1/4e^2) sum_{i<j} |F^w_{ij}(x)|^2
 * where F^w_{ij}(x) depends on {q(x), D_d q(x), w(x), D_d w(x)}.
 *
 * The variation delta e(x) w.r.t. field perturbation at point y:
 *   If y != x: only through D_d q(x) and D_d w(x) which feel delta at y
 *              via the stencil weights (handled by pi-field divergence)
 *   If y == x: through q(x), w(x), and D_d q(x), D_d w(x)
 *              (local terms + stencil self-weight terms)
 *
 * We separate: delta_F^w = (delta from local q,w) + (delta from D_d q, D_d w)
 *
 * For D_d terms: we define "response vectors"
 *   dF/d(D_d q_a) and dF/d(D_d w_alpha)
 * and build pi-fields that get divergenced.
 *
 * For local terms: we compute dF/d(q_a) and dF/d(w_alpha) directly.
 *
 * IMPLEMENTATION: We store pi-fields for q (4×3=12 scalars) and
 * w (4×3=12 scalars), then apply divergence in the force assembly.
 */

/* Precomputed data at each grid point */
typedef struct {
    Quat A[3];
    Wquat W[3];
    Wquat Fw[3]; /* F^w_{01}, F^w_{02}, F^w_{12} */
} CouplingPre;

static void precompute_coupling(const Field *f, CouplingPre *pre)
{
    int N = f->N;

    #pragma omp parallel for collapse(3) schedule(static)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        int ix = idx(N, i, j, k);
        Quat q = q_from_mv(f->psi[ix]);
        Quat qr = q_rev(q);
        Wquat w = w_from_mv(f->psi[ix]);
        Wquat wr = w_rev(w);

        for (int d = 0; d < 3; d++) {
            Quat dq = q_deriv(f, i, j, k, d);
            Wquat dw = w_deriv(f, i, j, k, d);
            pre[ix].A[d] = q_mul(qr, dq);
            pre[ix].W[d] = w_add(qw_mul(qr, dw), wq_mul(wr, dq));
        }

        pre[ix].Fw[0] = w_sub(qw_commutator(pre[ix].A[0], pre[ix].W[1]),
                               qw_commutator(pre[ix].A[1], pre[ix].W[0]));
        pre[ix].Fw[1] = w_sub(qw_commutator(pre[ix].A[0], pre[ix].W[2]),
                               qw_commutator(pre[ix].A[2], pre[ix].W[0]));
        pre[ix].Fw[2] = w_sub(qw_commutator(pre[ix].A[1], pre[ix].W[2]),
                               qw_commutator(pre[ix].A[2], pre[ix].W[1]));
    }
}

/* Get F^w_{d1,d2} from the precomputed array (handles d1>d2 via antisymmetry) */
static inline Wquat get_Fw(const CouplingPre *pre, int ix, int d1, int d2)
{
    if (d1 == d2) return w_zero();
    if (d1 < d2) {
        int ci = (d1==0 && d2==1) ? 0 : (d1==0 && d2==2) ? 1 : 2;
        return pre[ix].Fw[ci];
    } else {
        int ci = (d2==0 && d1==1) ? 0 : (d2==0 && d1==2) ? 1 : 2;
        return w_scale(-1.0, pre[ix].Fw[ci]);
    }
}

void coupling_gradient(const Field *f, const Params *p, double g_coupling,
                       Multivector *force)
{
    int N = f->N;
    int N3 = N * N * N;
    double inv_2e2 = 1.0 / (2.0 * p->e_skyrme * p->e_skyrme);
    double g2 = g_coupling * g_coupling;

    /* Precompute coupling data */
    CouplingPre *pre = (CouplingPre *)malloc((size_t)N3 * sizeof(CouplingPre));
    if (!pre) { fprintf(stderr, "coupling_gradient: malloc failed\n"); exit(1); }
    precompute_coupling(f, pre);

    /* === Build pi-fields for E_{4,C} ===
     *
     * The variation of E_{4,C} involves delta(D_d q) and delta(D_d w).
     * Through A_d = q~ * D_d q: delta A_d has a term q~ * delta(D_d q).
     * Through W_d = q~ * D_d w + w~ * D_d q:
     *   delta W_d has terms q~ * delta(D_d w) and w~ * delta(D_d q).
     *
     * The variation of F^w_{ij} w.r.t. delta(D_d q) propagates via:
     *   dA_d -> [dA_d, W_j] + [A_i, 0] - [dA_d,W_i] - [A_j,0] for (i,j) pairs involving d
     *   dW_d -> contribution when d is one of i,j
     *
     * For direction d, the "G-tensor" response:
     *   dE/d(D_d q) contributes through both A_d and W_d.
     *
     * Let's define the response vectors more carefully.
     *
     * E_{4,C} = (1/4e^2) sum_{i<j} |F^w_{ij}|^2
     * dE_{4,C} = (1/2e^2) sum_{i<j} F^w_{ij} . dF^w_{ij}
     *
     * F^w_{ij} = [A_i, W_j] - [A_j, W_i]
     *
     * dF^w_{ij} from dA_d (only when d=i or d=j):
     *   If d=i: + [dA_i, W_j] - 0 = [dA_i, W_j]    (and -[dA_j, W_i] = 0 since d!=j)
     *   If d=j: - [dA_j, W_i] + 0 = -[dA_j, W_i]
     *
     * dF^w_{ij} from dW_d (only when d=i or d=j):
     *   If d=j: + [A_i, dW_j]
     *   If d=i: - [A_j, dW_i]     (wait, check signs)
     *   More carefully: F^w_{ij} = [A_i, W_j] - [A_j, W_i]
     *   dF from dW_d=j: [A_i, dW_j]
     *   dF from dW_d=i: -[A_j, dW_i]
     *
     * So for pair (i,j) and perturbation in direction d:
     *   if d == i:  dF = [dA_i, W_j] - [A_j, dW_i]
     *   if d == j:  dF = [A_i, dW_j] - [dA_j, W_i]
     *
     * Now, dA_d = q~ * delta(D_d q), so through delta(D_d q_a):
     *   dA_d = q~ * eps_a * delta(D_d q_a)
     * Wait, that's not right. dA_d = d(q~) * D_d q + q~ * d(D_d q).
     * Through delta(D_d q): dA_d = q~ * delta(D_d q).
     * If delta(D_d q) = eps_a (unit perturbation), then dA_d = q~ * eps_a.
     *
     * Similarly, dW_d from delta(D_d q): dW_d = w~ * delta(D_d q) = w~ * eps_a
     * And dW_d from delta(D_d w): dW_d = q~ * delta(D_d w).
     *
     * So for a delta(D_d q_a) perturbation at direction d, the response on F^w is:
     * For each pair (i,j) with i<j:
     *   if d==i: F^w_{ij} += [q~*eps_a, W_j] - [A_j, w~*eps_a]
     *   if d==j: F^w_{ij} += [A_i, w~*eps_a] - [q~*eps_a, W_i]
     *
     * And for a delta(D_d w_alpha) perturbation:
     * For each pair (i,j) with i<j:
     *   if d==j: F^w_{ij} += [A_i, q~*e_alpha]
     *   if d==i: F^w_{ij} += -[A_j, q~*e_alpha]
     *
     * The pi-field for q: pi_q[4*d + a] = sum_{i<j with d=i or d=j} F^w_{ij} . (response to D_d q_a)
     * The pi-field for w: pi_w[4*d + alpha] = sum_{i<j with d=i or d=j} F^w_{ij} . (response to D_d w_alpha)
     */

    /* Pi-fields: pi_q[12*ix + 4*d + a] for q, pi_w[12*ix + 4*d + a] for w */
    double *pi_q = (double *)calloc((size_t)N3 * 12, sizeof(double));
    double *pi_w = (double *)calloc((size_t)N3 * 12, sizeof(double));
    if (!pi_q || !pi_w) { fprintf(stderr, "coupling_gradient: malloc failed\n"); exit(1); }

    #pragma omp parallel for collapse(3) schedule(static)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        int ix = idx(N, i, j, k);
        Quat qr = q_rev(q_from_mv(f->psi[ix]));
        Wquat wr = w_rev(w_from_mv(f->psi[ix]));

        for (int d = 0; d < 3; d++) {
            /* For each pair (d1,d2) with d1<d2 that includes d */
            /* Sum F^w_{d1,d2} . (response to delta D_d q_a) for each a */
            /* Sum F^w_{d1,d2} . (response to delta D_d w_alpha) for each alpha */

            double pq[4] = {0,0,0,0};
            double pw[4] = {0,0,0,0};

            for (int dp = 0; dp < 3; dp++) {
                if (dp == d) continue;
                int d1 = (d < dp) ? d : dp;
                int d2 = (d < dp) ? dp : d;
                Wquat Fwij = get_Fw(pre, ix, d1, d2);

                /* Response to delta(D_d q_a):
                 * if d==d1: dF = [q~*eps_a, W_{d2}] - [A_{d2}, w~*eps_a]
                 * if d==d2: dF = [A_{d1}, w~*eps_a] - [q~*eps_a, W_{d1}] */

                for (int a = 0; a < 4; a++) {
                    /* q~ * eps_a */
                    Quat qr_ea = q_mul(qr, (Quat){a==0,a==1,a==2,a==3});
                    /* w~ * eps_a: this is wq_mul(wr, eps_a) but eps_a is bulk
                     * Actually w~ * eps_a: w~ is weight, eps_a is bulk
                     * The result is wq_mul(wr, eps_a) which is in weight sector */
                    Wquat wr_ea = wq_mul(wr, (Quat){a==0,a==1,a==2,a==3});

                    Wquat response;
                    if (d == d1) {
                        /* dF = [qr_ea, W_{d2}] - [A_{d2}, wr_ea] */
                        response = w_sub(qw_commutator(qr_ea, pre[ix].W[d2]),
                                         qw_commutator(pre[ix].A[d2], wr_ea));
                    } else {
                        /* d == d2: dF = [A_{d1}, wr_ea] - [qr_ea, W_{d1}] */
                        response = w_sub(qw_commutator(pre[ix].A[d1], wr_ea),
                                         qw_commutator(qr_ea, pre[ix].W[d1]));
                    }
                    pq[a] += w_dot(Fwij, response);
                }

                /* Response to delta(D_d w_alpha):
                 * delta W_d from D_d(delta w) = q~ * delta(D_d w)
                 * So delta W_d = qw_mul(qr, e_alpha) for unit perturbation.
                 * if d==d1: dF = -[A_{d2}, qw_mul(qr, e_alpha)]   (d1=i in -[A_j,dW_i])
                 * if d==d2: dF = [A_{d1}, qw_mul(qr, e_alpha)]    (d2=j in [A_i,dW_j]) */
                for (int a = 0; a < 4; a++) {
                    Wquat qr_ew = qw_mul(qr, (Wquat){a==0,a==1,a==2,a==3});
                    Wquat response;
                    if (d == d1) {
                        response = w_scale(-1.0, qw_commutator(pre[ix].A[d2], qr_ew));
                    } else {
                        response = qw_commutator(pre[ix].A[d1], qr_ew);
                    }
                    pw[a] += w_dot(Fwij, response);
                }
            }

            for (int a = 0; a < 4; a++) {
                pi_q[12*ix + 4*d + a] = pq[a];
                pi_w[12*ix + 4*d + a] = pw[a];
            }
        }
    }

    /* E_int intermediate field */
    double *phi_int = NULL;
    if (g2 > 0) {
        phi_int = (double *)calloc((size_t)N3 * 12, sizeof(double));
        if (!phi_int) { fprintf(stderr, "coupling_gradient: malloc failed\n"); exit(1); }

        #pragma omp parallel for collapse(3) schedule(static)
        for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
        for (int k = 0; k < N; k++) {
            int ix = idx(N, i, j, k);
            double q2 = q_norm2(q_from_mv(f->psi[ix]));
            for (int d = 0; d < 3; d++) {
                Wquat dw = w_deriv(f, i, j, k, d);
                phi_int[12*ix + 4*d + 0] = q2 * dw.j1;
                phi_int[12*ix + 4*d + 1] = q2 * dw.j2;
                phi_int[12*ix + 4*d + 2] = q2 * dw.j3;
                phi_int[12*ix + 4*d + 3] = q2 * dw.p;
            }
        }
    }

    /* === Assemble all forces === */
    #pragma omp parallel for collapse(3) schedule(static)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        int ix = idx(N, i, j, k);
        Quat q = q_from_mv(f->psi[ix]);
        Quat qr = q_rev(q);
        Wquat w = w_from_mv(f->psi[ix]);
        Wquat wr = w_rev(w);

        /* --- E_{2,D} force on w: consistent Laplacian --- */
        Wquat lap = consistent_lap_w(f, i, j, k);
        force[ix].j1 += lap.j1;
        force[ix].j2 += lap.j2;
        force[ix].j3 += lap.j3;
        force[ix].p  += lap.p;

        /* --- E_{4,C} forces --- */

        /* Local terms from variation of q(x) and w(x) (not through derivatives).
         *
         * A_d = q~ * D_d q. Variation of q(x) -> d(q~) = -eps_a * q~ (for bulk comp a)
         * gives dA_d = -eps_a * q~ * D_d q = -eps_a * A_d.
         *
         * W_d = q~ * D_d w + w~ * D_d q.
         * Variation of q(x) -> d(q~) = sigma_a * eps_a * q~ gives
         *   dW_d(from dq) = sigma_a * eps_a * q~ * D_d w = sigma_a * qw_mul(eps_a * qr, dw_d)
         * Wait, that's not right either. q~ = rev(q), and d(q~)/dq_a...
         *
         * Let me be more careful. If q = sum_a q_a * eps_a, then
         * q~ = q_0 - q_1 e23 - q_2 e31 - q_3 e12 (reversion).
         * d(q~)/d(q_a) = sigma_a * eps_a where sigma = (+1,-1,-1,-1).
         *
         * So: dA_d/d(q_a)|_local = sigma_a * eps_a * D_d q
         *   (this is from the q~ factor; the D_d q factor is non-local)
         *
         * dW_d/d(q_a)|_local = sigma_a * qw_mul(eps_a, D_d w)
         *   (from q~ factor in first term)
         *
         * dW_d/d(w_alpha)|_local = wq_mul(d(w~)/d(w_alpha), D_d q)
         *   w~ = (-j1, -j2, -j3, +p), so d(w~)/d(w_alpha) = tau_alpha * e_alpha
         *   where tau = (-1,-1,-1,+1).
         *   = tau_alpha * wq_mul(e_alpha, D_d q)
         *
         * Now the full force from E_{4,C}:
         * F_q_a = -(1/2e^2) sum_{i<j} F^w_{ij} . dF^w_{ij}/d(q_a)
         * F_w_alpha = -(1/2e^2) sum_{i<j} F^w_{ij} . dF^w_{ij}/d(w_alpha)
         *
         * where dF^w_{ij}/d(q_a)|_local propagates through dA and dW:
         *   For pair (i,j): dF = [dA_i, W_j] + [A_i, dW_j] - [dA_j, W_i] - [A_j, dW_i]
         *   where dA_d = sigma_a * eps_a * D_d q (local q term)
         *         dW_d = sigma_a * qw_mul(eps_a, D_d w) (local q term through q~)
         *
         * And dF^w_{ij}/d(w_alpha)|_local propagates through dW:
         *   dW_d = tau_alpha * wq_mul(e_alpha, D_d q) (local w term through w~)
         */

        double sigma_q[4] = {1.0, -1.0, -1.0, -1.0};
        double tau_w[4] = {-1.0, -1.0, -1.0, 1.0};

        double f4_q[4] = {0, 0, 0, 0};
        double f4_w[4] = {0, 0, 0, 0};

        for (int d1 = 0; d1 < 3; d1++)
        for (int d2 = d1+1; d2 < 3; d2++) {
            Wquat Fwij = get_Fw(pre, ix, d1, d2);

            /* Derivatives at this point */
            Quat dq1 = q_deriv(f, i, j, k, d1);
            Quat dq2 = q_deriv(f, i, j, k, d2);
            Wquat dw1 = w_deriv(f, i, j, k, d1);
            Wquat dw2 = w_deriv(f, i, j, k, d2);

            /* Force on q_a (local) */
            for (int a = 0; a < 4; a++) {
                Quat ea = {a==0, a==1, a==2, a==3};

                /* dA_d from local q: sigma_a * eps_a * D_d q */
                Quat dA1 = q_scale(sigma_q[a], q_mul(ea, dq1));
                Quat dA2 = q_scale(sigma_q[a], q_mul(ea, dq2));

                /* dW_d from local q: sigma_a * qw_mul(eps_a, D_d w)
                 * Wait — need to think about this. q~ acts on D_d w.
                 * d(q~)/dq_a = sigma_a * eps_a.
                 * So d(qw_mul(q~, D_d w))/dq_a = sigma_a * qw_mul(eps_a, D_d w). */
                Wquat dW1_q = w_scale(sigma_q[a], qw_mul(ea, dw1));
                Wquat dW2_q = w_scale(sigma_q[a], qw_mul(ea, dw2));

                /* dF = [dA1, W2] + [A1, dW2_q] - [dA2, W1] - [A2, dW1_q] */
                Wquat resp = w_zero();
                resp = w_add(resp, qw_commutator(dA1, pre[ix].W[d2]));
                resp = w_add(resp, qw_commutator(pre[ix].A[d1], dW2_q));
                resp = w_sub(resp, qw_commutator(dA2, pre[ix].W[d1]));
                resp = w_sub(resp, qw_commutator(pre[ix].A[d2], dW1_q));

                f4_q[a] += w_dot(Fwij, resp);
            }

            /* Force on w_alpha (local) */
            for (int a = 0; a < 4; a++) {
                Wquat ea = {a==0, a==1, a==2, a==3};

                /* dW_d from local w: tau_alpha * wq_mul(e_alpha, D_d q) */
                Wquat dW1_w = w_scale(tau_w[a], wq_mul(ea, dq1));
                Wquat dW2_w = w_scale(tau_w[a], wq_mul(ea, dq2));

                /* dF = [A1, dW2_w] - [A2, dW1_w] */
                Wquat resp = w_sub(qw_commutator(pre[ix].A[d1], dW2_w),
                                   qw_commutator(pre[ix].A[d2], dW1_w));

                f4_w[a] += w_dot(Fwij, resp);
            }
        }

        /* Divergence terms: -D_d(pi_q) and -D_d(pi_w)
         * These come from the D_d(delta q) and D_d(delta w) contributions
         * in the right-currents, integrated by parts. */
        double inv12h = 1.0 / (12.0 * f->h);
        double div_q[4] = {0,0,0,0};
        double div_w[4] = {0,0,0,0};

        for (int a = 0; a < 4; a++) {
            for (int d = 0; d < 3; d++) {
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
                /* Divergence of pi_q */
                {
                    double pm2v = pi_q[12*im2 + 4*d + a];
                    double pm1v = pi_q[12*im1 + 4*d + a];
                    double pp1v = pi_q[12*ip1 + 4*d + a];
                    double pp2v = pi_q[12*ip2 + 4*d + a];
                    div_q[a] += (-pp2v + 8*pp1v - 8*pm1v + pm2v) * inv12h;
                }
                /* Divergence of pi_w */
                {
                    double pm2v = pi_w[12*im2 + 4*d + a];
                    double pm1v = pi_w[12*im1 + 4*d + a];
                    double pp1v = pi_w[12*ip1 + 4*d + a];
                    double pp2v = pi_w[12*ip2 + 4*d + a];
                    div_w[a] += (-pp2v + 8*pp1v - 8*pm1v + pm2v) * inv12h;
                }
            }
        }

        /* Apply E_{4,C} force: force = -dE/dpsi.
         * E_{4,C} = +(1/4e²) sum |F^w|² uses Euclidean dot (always positive),
         * so dE = +(1/2e²)(local - div), and force = -dE = -(1/2e²)(local - div).
         * This differs from field.c's E₄ which uses Clifford <C²>₀ = -|C|²,
         * giving dE₄ = -(1/2e²)<dR,G>₀ and force = +(1/2e²)(local - div). */
        force[ix].s  -= inv_2e2 * (f4_q[0] - div_q[0]);
        force[ix].f1 -= inv_2e2 * (f4_q[1] - div_q[1]);
        force[ix].f2 -= inv_2e2 * (f4_q[2] - div_q[2]);
        force[ix].f3 -= inv_2e2 * (f4_q[3] - div_q[3]);
        force[ix].j1 -= inv_2e2 * (f4_w[0] - div_w[0]);
        force[ix].j2 -= inv_2e2 * (f4_w[1] - div_w[1]);
        force[ix].j3 -= inv_2e2 * (f4_w[2] - div_w[2]);
        force[ix].p  -= inv_2e2 * (f4_w[3] - div_w[3]);

        /* --- E_int forces --- */
        if (g2 > 0) {
            /* Force on w: g^2 * sum_d D_d(|q|^2 * D_d w_alpha) */
            double f_int_w[4] = {0, 0, 0, 0};
            for (int a = 0; a < 4; a++) {
                for (int d = 0; d < 3; d++) {
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
                    double pm2v = phi_int[12*im2 + 4*d + a];
                    double pm1v = phi_int[12*im1 + 4*d + a];
                    double pp1v = phi_int[12*ip1 + 4*d + a];
                    double pp2v = phi_int[12*ip2 + 4*d + a];
                    f_int_w[a] += (-pp2v + 8*pp1v - 8*pm1v + pm2v) * inv12h;
                }
            }

            force[ix].j1 += g2 * f_int_w[0];
            force[ix].j2 += g2 * f_int_w[1];
            force[ix].j3 += g2 * f_int_w[2];
            force[ix].p  += g2 * f_int_w[3];

            /* Force on q from E_int: -g^2 * q_a * sum_d |D_d w|^2 */
            double grad_w2 = 0;
            for (int d = 0; d < 3; d++) {
                Wquat dw = w_deriv(f, i, j, k, d);
                grad_w2 += w_norm2(dw);
            }
            force[ix].s  -= g2 * q.s  * grad_w2;
            force[ix].f1 -= g2 * q.f1 * grad_w2;
            force[ix].f2 -= g2 * q.f2 * grad_w2;
            force[ix].f3 -= g2 * q.f3 * grad_w2;
        }

        (void)wr; (void)qr; (void)w;
    }

    free(pi_q);
    free(pi_w);
    free(pre);
    if (phi_int) free(phi_int);
}
