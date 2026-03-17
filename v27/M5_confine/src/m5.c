/*
 * m5.c — V27-M5: Topological Confinement via Chern-Simons Term
 *
 * Uses m=0, mu=-50, kappa=50 propagating braid (V27-M4 optimal config).
 *
 * Tests:
 *   M5a: Compute linking number (Chern-Simons density) at t=0,100,...,500
 *   M5b: Add CS potential L_CS_pot = lambda_CS * (Q_CS - Q0)^2
 *   M5c: Compare propagating braid vs static braid topology retention
 *
 * Compile: gcc -O3 -fopenmp -Wall -o m5 src/m5.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

/* ─── Parameters ─── */
static double mu_pot    = -50.0;
static double kappa     = 50.0;
static double mass      = 0.0;
static double A0        = 0.8;
static double R_tube    = 3.0;
static double lambda_pw = 0.0;
static double lambda_CS = 0.0;   /* CS potential strength */
static double Q_CS_0    = 0.0;   /* initial CS charge (set at t=0) */

static int    N         = 128;
static double L         = 20.0;
static double tfinal    = 500.0;
static double cfl_frac  = 0.20;
static char   outdir[512] = "data";

/* ─── Index helpers ─── */
#define IDX(i,j,k) ((long)(i)*N*N + (long)(j)*N + (long)(k))

/* ─── Globals ─── */
static double *phi[3], *vel[3], *acc[3];
static double *damp;
static double dx, dx2, m2, dt;
static long Ngrid;

/* ─── Periodic z wrap ─── */
static inline int wrap_z(int k)
{
    if (k < 0)   return k + N;
    if (k >= N)  return k - N;
    return k;
}

/* ─── Initialization: Propagating Helical Braid ─── */
static void init_braid_propagating(void)
{
    double Lz = 2.0 * L;
    double k_wave = 2.0 * M_PI / Lz;
    double omega = sqrt(k_wave * k_wave + m2);

    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < Ngrid; idx++) {
        int ii = idx / (N * N);
        int jj = (idx / N) % N;
        int kk = idx % N;
        double x = -L + ii * dx;
        double y = -L + jj * dx;

        double z = -L + kk * dx;

        double r_perp = sqrt(x * x + y * y);
        double envelope = A0 * exp(-r_perp * r_perp / (2.0 * R_tube * R_tube));

        for (int a = 0; a < 3; a++) {
            double phase = k_wave * z + 2.0 * M_PI * a / 3.0;
            phi[a][idx] = envelope * cos(phase);
            vel[a][idx] = omega * envelope * sin(phase);
        }
    }
}

/* ─── Initialization: Static Helical Braid (no propagation) ─── */
static void init_braid_static(void)
{
    double Lz = 2.0 * L;
    double k_wave = 2.0 * M_PI / Lz;

    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < Ngrid; idx++) {
        int ii = idx / (N * N);
        int jj = (idx / N) % N;
        int kk = idx % N;
        double x = -L + ii * dx;
        double y = -L + jj * dx;
        double z = -L + kk * dx;

        double r_perp = sqrt(x * x + y * y);
        double envelope = A0 * exp(-r_perp * r_perp / (2.0 * R_tube * R_tube));

        for (int a = 0; a < 3; a++) {
            double phase = k_wave * z + 2.0 * M_PI * a / 3.0;
            phi[a][idx] = envelope * cos(phase);
            vel[a][idx] = 0.0;  /* static: no initial velocity */
        }
    }
}

/* ─── Compute Chern-Simons density: rho_CS = eps_{abc} eps_{ijk} phi_a (d_i phi_b)(d_j phi_c) ─── */
/*
 * eps_{abc} has 6 terms (3! permutations, 3 even + 3 odd).
 * eps_{ijk} similarly. The full sum is:
 *   rho_CS = sum_{abc,ijk} eps_{abc} eps_{ijk} phi_a (d_i phi_b)(d_j phi_c)
 *
 * Expanding eps_{abc}: (0,1,2)=+1, (1,2,0)=+1, (2,0,1)=+1,
 *                      (0,2,1)=-1, (2,1,0)=-1, (1,0,2)=-1
 *
 * For a fixed (a,b,c) = (0,1,2) with eps=+1:
 *   contribution = phi_0 * [eps_{ijk} (d_i phi_1)(d_j phi_2)]
 *                = phi_0 * [(d_x phi_1)(d_y phi_2) - (d_y phi_1)(d_x phi_2)
 *                          +(d_y phi_1)(d_z phi_2) - (d_z phi_1)(d_y phi_2)
 *                          +(d_z phi_1)(d_x phi_2) - (d_x phi_1)(d_z phi_2)]
 *                = phi_0 * 2 * [curl-like cross of grad(phi_1) x grad(phi_2)]
 *
 * Actually eps_{ijk}(d_i phi_b)(d_j phi_c) = 2*(grad phi_b x grad phi_c)_k summed.
 * No, that's the k-component. We need the FULL contraction:
 *   eps_{ijk}(d_i phi_b)(d_j phi_c) sums over all i,j,k giving a SCALAR?
 * No - eps_{ijk} has 3 indices but we're contracting over i,j only with k free.
 * Wait, the expression is eps_{ijk} with all three summed:
 *   sum_{ijk} eps_{ijk} (d_i phi_b)(d_j phi_c) = 0 for ANY b,c
 * because eps_{ijk} is antisymmetric in i,j but (d_i phi_b)(d_j phi_c) has
 * no particular symmetry when b != c.
 *
 * Actually: sum_{ijk} eps_{ijk} (d_i phi_b)(d_j phi_c)
 *   = sum_k [sum_{ij} eps_{ijk} (d_i phi_b)(d_j phi_c)]
 *   = sum_k (grad phi_b x grad phi_c)_k   -- this is ZERO (dot product of cross with nothing)
 *
 * No! The k index in eps_{ijk} is also summed. So:
 *   sum_{ijk} eps_{ijk} X_i Y_j = sum_k (X x Y)_k -- but that's not right either.
 *
 * Let me be more careful. eps_{ijk} with all 3 indices summed against 2 vectors:
 *   sum_{ijk} eps_{ijk} A_i B_j = ?
 * This IS zero because for each k we get (A x B)_k, and summing over k gives nothing
 * (each k contributes independently, not summed to another vector).
 *
 * Wait, I think the actual CS density uses only the SPATIAL cross product components.
 * The correct form is:
 *
 * rho_CS = eps_{abc} phi_a (grad phi_b x grad phi_c) . e_z   [for a 1D linking along z]
 *        = eps_{abc} phi_a [(d_x phi_b)(d_y phi_c) - (d_y phi_b)(d_x phi_c)]
 *
 * But more generally, the topological charge density for three scalar fields in 3D:
 *   Q = eps_{abc} eps_{ijk} (d_i phi_a)(d_j phi_b)(d_k phi_c)
 *     = 6 det(d_i phi_a)  [the Jacobian determinant of the map (x,y,z) -> (phi1,phi2,phi3)]
 *
 * This is the DEGREE of the map. For our braid, this should be related to linking.
 *
 * So we compute TWO things:
 * 1) Q_jac = integral of det(J) where J_{ai} = d_i phi_a  (Jacobian density)
 * 2) Q_CS  = integral of eps_{abc} phi_a (d_x phi_b d_y phi_c - d_y phi_b d_x phi_c) (CS form)
 *
 * Both are topological invariants if the fields are smooth and well-behaved.
 */

typedef struct {
    double Q_jac;       /* Jacobian density integral: det(d_i phi_a) */
    double Q_CS_xy;     /* CS density with z-component of cross product */
    double Q_CS_yz;     /* CS density with x-component of cross product */
    double Q_CS_zx;     /* CS density with y-component of cross product */
    double H_12;        /* Helicity integral for pair (1,2) */
    double H_23;        /* Helicity integral for pair (2,3) */
    double H_31;        /* Helicity integral for pair (3,1) */
    double phase_wind;  /* Phase winding number along z */
} TopoCharges;

static TopoCharges compute_topo(void)
{
    TopoCharges tc;
    memset(&tc, 0, sizeof(tc));

    double Q_jac = 0, Q_CS_xy = 0, Q_CS_yz = 0, Q_CS_zx = 0;
    double H_12 = 0, H_23 = 0, H_31 = 0;

    #pragma omp parallel
    {
        double lQj = 0, lQxy = 0, lQyz = 0, lQzx = 0;
        double lH12 = 0, lH23 = 0, lH31 = 0;

        #pragma omp for schedule(static) nowait
        for (int i = 2; i < N-2; i++) {
            for (int j = 2; j < N-2; j++) {
                for (int k = 0; k < N; k++) {
                    long idx = IDX(i, j, k);
                    double dV = dx * dx * dx;

                    int km1 = wrap_z(k - 1);
                    int kp1 = wrap_z(k + 1);
                    double inv2dx = 1.0 / (2.0 * dx);

                    /* Gradients of all three fields */
                    double dphi[3][3]; /* dphi[a][dir]: d_{dir} phi_a */
                    for (int a = 0; a < 3; a++) {
                        dphi[a][0] = (phi[a][IDX(i+1,j,k)] - phi[a][IDX(i-1,j,k)]) * inv2dx;
                        dphi[a][1] = (phi[a][IDX(i,j+1,k)] - phi[a][IDX(i,j-1,k)]) * inv2dx;
                        dphi[a][2] = (phi[a][IDX(i,j,kp1)] - phi[a][IDX(i,j,km1)]) * inv2dx;
                    }

                    /* Jacobian determinant: det(J) where J_{ai} = d_i phi_a
                     * det = phi0_x(phi1_y*phi2_z - phi1_z*phi2_y)
                     *      -phi0_y(phi1_x*phi2_z - phi1_z*phi2_x)
                     *      +phi0_z(phi1_x*phi2_y - phi1_y*phi2_x)
                     */
                    double det = dphi[0][0] * (dphi[1][1]*dphi[2][2] - dphi[1][2]*dphi[2][1])
                               - dphi[0][1] * (dphi[1][0]*dphi[2][2] - dphi[1][2]*dphi[2][0])
                               + dphi[0][2] * (dphi[1][0]*dphi[2][1] - dphi[1][1]*dphi[2][0]);
                    lQj += det * dV;

                    /* CS densities: eps_{abc} phi_a (d_i phi_b)(d_j phi_c) for each ij-plane
                     * xy-plane: phi_a * (d_x phi_b * d_y phi_c - d_y phi_b * d_x phi_c)
                     */
                    double p0 = phi[0][idx], p1 = phi[1][idx], p2 = phi[2][idx];

                    /* xy cross: (d_x phi_b)(d_y phi_c) - (d_y phi_b)(d_x phi_c) */
                    double cross_xy_12 = dphi[1][0]*dphi[2][1] - dphi[1][1]*dphi[2][0];
                    double cross_xy_20 = dphi[2][0]*dphi[0][1] - dphi[2][1]*dphi[0][0];
                    double cross_xy_01 = dphi[0][0]*dphi[1][1] - dphi[0][1]*dphi[1][0];
                    /* eps_{012}=+1: phi0 * cross(1,2)_xy
                     * eps_{120}=+1: phi1 * cross(2,0)_xy
                     * eps_{201}=+1: phi2 * cross(0,1)_xy */
                    lQxy += (p0*cross_xy_12 + p1*cross_xy_20 + p2*cross_xy_01) * dV;

                    /* yz cross */
                    double cross_yz_12 = dphi[1][1]*dphi[2][2] - dphi[1][2]*dphi[2][1];
                    double cross_yz_20 = dphi[2][1]*dphi[0][2] - dphi[2][2]*dphi[0][1];
                    double cross_yz_01 = dphi[0][1]*dphi[1][2] - dphi[0][2]*dphi[1][1];
                    lQyz += (p0*cross_yz_12 + p1*cross_yz_20 + p2*cross_yz_01) * dV;

                    /* zx cross */
                    double cross_zx_12 = dphi[1][2]*dphi[2][0] - dphi[1][0]*dphi[2][2];
                    double cross_zx_20 = dphi[2][2]*dphi[0][0] - dphi[2][0]*dphi[0][2];
                    double cross_zx_01 = dphi[0][2]*dphi[1][0] - dphi[0][0]*dphi[1][2];
                    lQzx += (p0*cross_zx_12 + p1*cross_zx_20 + p2*cross_zx_01) * dV;

                    /* Helicity for pair (a,b): H_{ab} = int A_i B_i d3x
                     * where A_i = eps_{ijk} phi_a d_j phi_b
                     * and   B_i = eps_{ilm} d_l phi_a d_m phi_b
                     *
                     * A_x = phi_a*(d_y phi_b) - phi_a*(d_z phi_b)  ... no:
                     * A_x = eps_{x,j,k} phi_a d_j phi_b
                     *      = phi_a*(d_y phi_b*delta_{k=z} - d_z phi_b*delta_{k=y})
                     * Wait, A is a vector: A_i = eps_{ijk} phi_a (d_j phi_b)
                     * Actually this requires a third index. Let me use:
                     *   A = phi_a * grad(phi_b)  ... no, that's not eps.
                     *
                     * A = phi_a * curl(phi_b * e_k)?  The definition is unclear.
                     * Let me use the simpler helicity:
                     *   H_{ab} = int (phi_a * grad(phi_b) - phi_b * grad(phi_a)) . curl(phi_a grad(phi_b)) d3x
                     * This is too complex for a diagnostic. Use the simpler cross-helicity:
                     *   H_{ab} = int (grad phi_a x grad phi_b) . (phi_a grad phi_b - phi_b grad phi_a) d3x
                     */
                    /* Cross helicity: (grad phi_a x grad phi_b) . (phi_a grad phi_b) */
                    /* For pair (0,1): */
                    {
                        double cx = dphi[0][1]*dphi[1][2] - dphi[0][2]*dphi[1][1];
                        double cy = dphi[0][2]*dphi[1][0] - dphi[0][0]*dphi[1][2];
                        double cz = dphi[0][0]*dphi[1][1] - dphi[0][1]*dphi[1][0];
                        double vx = p0*dphi[1][0] - p1*dphi[0][0];
                        double vy = p0*dphi[1][1] - p1*dphi[0][1];
                        double vz = p0*dphi[1][2] - p1*dphi[0][2];
                        lH12 += (cx*vx + cy*vy + cz*vz) * dV;
                    }
                    /* For pair (1,2): */
                    {
                        double cx = dphi[1][1]*dphi[2][2] - dphi[1][2]*dphi[2][1];
                        double cy = dphi[1][2]*dphi[2][0] - dphi[1][0]*dphi[2][2];
                        double cz = dphi[1][0]*dphi[2][1] - dphi[1][1]*dphi[2][0];
                        double vx = p1*dphi[2][0] - p2*dphi[1][0];
                        double vy = p1*dphi[2][1] - p2*dphi[1][1];
                        double vz = p1*dphi[2][2] - p2*dphi[1][2];
                        lH23 += (cx*vx + cy*vy + cz*vz) * dV;
                    }
                    /* For pair (2,0): */
                    {
                        double cx = dphi[2][1]*dphi[0][2] - dphi[2][2]*dphi[0][1];
                        double cy = dphi[2][2]*dphi[0][0] - dphi[2][0]*dphi[0][2];
                        double cz = dphi[2][0]*dphi[0][1] - dphi[2][1]*dphi[0][0];
                        double vx = p2*dphi[0][0] - p0*dphi[2][0];
                        double vy = p2*dphi[0][1] - p0*dphi[2][1];
                        double vz = p2*dphi[0][2] - p0*dphi[2][2];
                        lH31 += (cx*vx + cy*vy + cz*vz) * dV;
                    }
                }
            }
        }

        #pragma omp critical
        {
            Q_jac += lQj; Q_CS_xy += lQxy; Q_CS_yz += lQyz; Q_CS_zx += lQzx;
            H_12 += lH12; H_23 += lH23; H_31 += lH31;
        }
    }

    tc.Q_jac = Q_jac;
    tc.Q_CS_xy = Q_CS_xy;
    tc.Q_CS_yz = Q_CS_yz;
    tc.Q_CS_zx = Q_CS_zx;
    tc.H_12 = H_12;
    tc.H_23 = H_23;
    tc.H_31 = H_31;

    /* Phase winding: measure how much the phase angle of (phi1+i*phi2) winds along z
     * at the tube center (x=0, y=0).
     * Total winding = (1/2pi) * integral of d(theta)/dz along z
     */
    {
        int ic = N / 2, jc = N / 2;
        double total_wind = 0.0;
        for (int k = 0; k < N; k++) {
            int kp = wrap_z(k + 1);
            long idx0 = IDX(ic, jc, k);
            long idx1 = IDX(ic, jc, kp);

            double re0 = phi[0][idx0], im0 = phi[1][idx0];
            double re1 = phi[0][idx1], im1 = phi[1][idx1];

            /* Phase difference via atan2 */
            double theta0 = atan2(im0, re0);
            double theta1 = atan2(im1, re1);
            double dtheta = theta1 - theta0;

            /* Wrap to [-pi, pi] */
            while (dtheta >  M_PI) dtheta -= 2.0 * M_PI;
            while (dtheta < -M_PI) dtheta += 2.0 * M_PI;

            total_wind += dtheta;
        }
        tc.phase_wind = total_wind / (2.0 * M_PI);
    }

    return tc;
}

/* ─── Compute CS density at a point (for CS force) ─── */
static double cs_density_at(int i, int j, int k)
{
    double inv2dx = 1.0 / (2.0 * dx);
    int km1 = wrap_z(k - 1);
    int kp1 = wrap_z(k + 1);

    double dphi[3][3];
    for (int a = 0; a < 3; a++) {
        dphi[a][0] = (phi[a][IDX(i+1,j,k)] - phi[a][IDX(i-1,j,k)]) * inv2dx;
        dphi[a][1] = (phi[a][IDX(i,j+1,k)] - phi[a][IDX(i,j-1,k)]) * inv2dx;
        dphi[a][2] = (phi[a][IDX(i,j,kp1)] - phi[a][IDX(i,j,km1)]) * inv2dx;
    }

    double p0 = phi[0][IDX(i,j,k)], p1 = phi[1][IDX(i,j,k)], p2 = phi[2][IDX(i,j,k)];

    /* Full CS density: eps_{abc} phi_a (grad phi_b x grad phi_c)
     * = phi_0*(grad1 x grad2) + phi_1*(grad2 x grad0) + phi_2*(grad0 x grad1)
     * summed over all 3 spatial components (dot with sum of unit vectors? No.)
     *
     * The scalar CS density is actually the Jacobian:
     * rho_CS = eps_{abc} eps_{ijk} (d_i phi_a)(d_j phi_b)(d_k phi_c) / 6
     * = det(J) where J_{ai} = d_i phi_a
     *
     * But the PROPOSAL says: rho_CS = eps_{abc} eps_{ijk} phi_a (d_i phi_b)(d_j phi_c)
     * which is different. Let me use the Jacobian determinant as the primary
     * topological charge, plus the proposal's version.
     */

    /* Jacobian determinant */
    double det = dphi[0][0] * (dphi[1][1]*dphi[2][2] - dphi[1][2]*dphi[2][1])
               - dphi[0][1] * (dphi[1][0]*dphi[2][2] - dphi[1][2]*dphi[2][0])
               + dphi[0][2] * (dphi[1][0]*dphi[2][1] - dphi[1][1]*dphi[2][0]);

    /* Proposal version: sum over k of eps_{abc} phi_a (grad_b x grad_c)_k
     * = sum_k [phi_0(grad1 x grad2)_k + phi_1(grad2 x grad0)_k + phi_2(grad0 x grad1)_k]
     * But sum_k of a vector = ... each component is independent.
     * Actually: sum_{ijk} eps_{ijk} phi_a d_i(phi_b) d_j(phi_c) with fixed a,b,c
     * and ijk all summed is NOT well-defined as written since eps has 3 free indices.
     *
     * I'll return the Jacobian det, which IS the proper topological density.
     */
    (void)p0; (void)p1; (void)p2;
    return det;
}

/* ─── Compute global CS charge (integral of Jacobian determinant) ─── */
static double compute_Q_CS_global(void)
{
    double Q = 0;
    #pragma omp parallel for schedule(static) reduction(+:Q)
    for (int i = 2; i < N-2; i++) {
        for (int j = 2; j < N-2; j++) {
            for (int k = 0; k < N; k++) {
                Q += cs_density_at(i, j, k) * dx * dx * dx;
            }
        }
    }
    return Q;
}

/* ─── Setup damping layer ─── */
static void setup_damping(void)
{
    double R_abs_inner = L * 0.70;
    double R_abs_outer = L * 0.95;
    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < Ngrid; idx++) {
        int ii = idx / (N * N);
        int jj = (idx / N) % N;
        double x = -L + ii * dx;
        double y = -L + jj * dx;
        double rperp = sqrt(x*x + y*y);
        if (rperp > R_abs_inner) {
            double f = (rperp - R_abs_inner) / (R_abs_outer - R_abs_inner);
            if (f > 1.0) f = 1.0;
            damp[idx] = 1.0 - 0.98 * f * f;
        } else {
            damp[idx] = 1.0;
        }
    }
}

/* ─── Compute acceleration ─── */
static void compute_acc(void)
{
    /* If CS potential is active, compute global Q_CS first */
    double Q_CS_now = 0.0;
    double cs_force_coeff = 0.0;
    if (lambda_CS > 0.0) {
        Q_CS_now = compute_Q_CS_global();
        /* Potential: V_CS = lambda_CS * (Q_CS - Q0)^2
         * Force on phi_a = -dV_CS/dphi_a = -2*lambda_CS*(Q_CS - Q0) * dQ_CS/dphi_a
         * dQ_CS/dphi_a at point x = d/dphi_a [det(J)] at x ... but det(J) doesn't
         * depend on phi_a directly, only on gradients of phi_a.
         *
         * Actually: the Jacobian density det(d_i phi_a) does NOT depend on phi_a
         * itself, only on its derivatives. So d(det)/d(phi_a) = 0 and there is
         * no direct local force.
         *
         * But: d(det)/d(d_i phi_a) != 0, so via integration by parts:
         * -dV_CS/dphi_a(x) = +2*lambda_CS*(Q_CS-Q0) * d_i [d(det)/d(d_i phi_a)]
         *
         * For det = eps_{abc} eps_{ijk}/6 * d_i phi_a d_j phi_b d_k phi_c
         * (using the standard Jacobian form):
         * d(det)/d(d_m phi_n) = cofactor_{nm} = (1/2) eps_{nbc} eps_{mjk} d_j phi_b d_k phi_c
         *
         * So: force on phi_n = +2*lambda_CS*(Q_CS-Q0) * d_m [cofactor_{nm}]
         *   = 2*lcs*(Q-Q0) * (1/2) eps_{nbc} eps_{mjk} d_m [d_j phi_b d_k phi_c]
         *   = lcs*(Q-Q0) * eps_{nbc} eps_{mjk} (d_m d_j phi_b)(d_k phi_c)
         *
         * This involves second derivatives. We'll compute it.
         */
        cs_force_coeff = 2.0 * lambda_CS * (Q_CS_now - Q_CS_0);
    }

    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (int i = 2; i < N-2; i++) {
            for (int j = 2; j < N-2; j++) {
                for (int k = 0; k < N; k++) {
                    long idx = IDX(i, j, k);

                    int km1 = wrap_z(k - 1);
                    int kp1 = wrap_z(k + 1);
                    double lapl = (phi[a][IDX(i+1,j,k)] + phi[a][IDX(i-1,j,k)]
                                 + phi[a][IDX(i,j+1,k)] + phi[a][IDX(i,j-1,k)]
                                 + phi[a][IDX(i,j,kp1)] + phi[a][IDX(i,j,km1)]
                                 - 6.0 * phi[a][idx]) / dx2;

                    /* Triple-product potential force */
                    double p0 = phi[0][idx], p1 = phi[1][idx], p2 = phi[2][idx];
                    double P = p0 * p1 * p2;
                    double P2 = P * P;
                    double denom2 = (1.0 + kappa * P2);
                    denom2 *= denom2;

                    double dP;
                    switch (a) {
                        case 0: dP = p1 * p2; break;
                        case 1: dP = p0 * p2; break;
                        default: dP = p0 * p1; break;
                    }
                    double dVdphi_triple = mu_pot * P * dP / denom2;

                    /* Pairwise coupling */
                    double dVdphi_pw = 0.0;
                    if (lambda_pw != 0.0) {
                        switch (a) {
                            case 0: dVdphi_pw = lambda_pw * (p1 + p2); break;
                            case 1: dVdphi_pw = lambda_pw * (p0 + p2); break;
                            default: dVdphi_pw = lambda_pw * (p0 + p1); break;
                        }
                    }

                    /* CS topological force */
                    double cs_acc = 0.0;
                    if (lambda_CS > 0.0 && fabs(cs_force_coeff) > 1e-30) {
                        /* force_a = cs_force_coeff * d_m [cofactor_{a,m}]
                         * cofactor_{a,m} = (1/2) eps_{abc} eps_{mjk} d_j phi_b d_k phi_c
                         *
                         * We need d_m of cofactor_{a,m} which requires second derivatives.
                         * Compute using central differences on the cofactor.
                         *
                         * For field a, the other two fields are b and c:
                         *   a=0: b=1, c=2, eps_{012}=+1
                         *   a=1: b=2, c=0, eps_{120}=+1
                         *   a=2: b=0, c=1, eps_{201}=+1
                         */
                        int b, c;
                        switch (a) {
                            case 0: b = 1; c = 2; break;
                            case 1: b = 2; c = 0; break;
                            default: b = 0; c = 1; break;
                        }

                        /* d_m cofactor_{a,m} = eps_{mjk} d_m(d_j phi_b * d_k phi_c)
                         * = eps_{mjk} [d_m d_j phi_b * d_k phi_c + d_j phi_b * d_m d_k phi_c]
                         *
                         * There are 6 nonzero terms from eps_{mjk}:
                         * (m,j,k) = (x,y,z): +[d_xx phi_b... no, d_m d_j = d_x d_y]
                         *
                         * Actually, since eps_{mjk} is antisymmetric in j,k:
                         * = eps_{mjk} d_m d_j phi_b * d_k phi_c + eps_{mjk} d_j phi_b * d_m d_k phi_c
                         * In second term, swap m<->j: = eps_{jmk} d_m phi_b * d_j d_k phi_c
                         *                            = -eps_{mjk} d_m phi_b * d_j d_k phi_c
                         * But d_j d_k phi_c is symmetric in j,k while eps_{mjk} is antisymmetric:
                         * => sum over j,k of eps_{mjk} d_j d_k phi_c = 0
                         *
                         * So: d_m cofactor = eps_{mjk} (d_m d_j phi_b)(d_k phi_c)
                         * (the second term vanishes!)
                         *
                         * This is: (curl of [d_j phi_b])_q dotted with d_q phi_c? No.
                         * Let me just enumerate the 6 terms:
                         * m=x,j=y,k=z: + (d_x d_y phi_b)(d_z phi_c)
                         * m=x,j=z,k=y: - (d_x d_z phi_b)(d_y phi_c)
                         * m=y,j=z,k=x: + (d_y d_z phi_b)(d_x phi_c)
                         * m=y,j=x,k=z: - (d_y d_x phi_b)(d_z phi_c)
                         * m=z,j=x,k=y: + (d_z d_x phi_b)(d_y phi_c)
                         * m=z,j=y,k=x: - (d_z d_y phi_b)(d_x phi_c)
                         *
                         * Since d_m d_j = d_j d_m, terms (1)+(4)=0, (2)+(5)=0, (3)+(6)=0.
                         * ALL CANCEL! d_m cofactor_{a,m} = 0 identically!
                         *
                         * This means the CS Jacobian density is a TOTAL DIVERGENCE and
                         * does not contribute any local force through the global penalty.
                         * The Q_CS is exactly conserved by any smooth dynamics.
                         *
                         * So the CS force is ZERO. We'll implement the ALTERNATIVE:
                         * use the PROPOSAL's rho_CS = eps_{abc} phi_a (d_i phi_b)(d_j phi_c)
                         * summed over a specific plane (xy), which is NOT a total divergence.
                         *
                         * V_CS = lambda_CS * (integral of rho_CS_xy - Q0)^2
                         * rho_CS_xy = eps_{abc} phi_a [(d_x phi_b)(d_y phi_c) - (d_y phi_b)(d_x phi_c)]
                         *
                         * d(rho_CS_xy)/d(phi_n) at point x:
                         *   = eps_{nbc}[(d_x phi_b)(d_y phi_c) - (d_y phi_b)(d_x phi_c)]   (from phi_a = phi_n)
                         *
                         * d(rho_CS_xy)/d(d_x phi_n) at point x:
                         *   = eps_{anc} phi_a (d_y phi_c) - eps_{acn} phi_a (d_y phi_c)   ...
                         * This gets complicated. Let's just do the direct force:
                         *
                         * force_n = -2*lcs*(Q-Q0) * [d(rho)/d(phi_n) - d_i(d(rho)/d(d_i phi_n))]
                         *
                         * Simpler: just use rho = eps_{abc} phi_a (d_x phi_b)(d_y phi_c) with eps summed.
                         * d(rho)/d(phi_n) = eps_{nbc}(d_x phi_b)(d_y phi_c) summed as:
                         *   n=0: +(dx phi1)(dy phi2) - (dx phi2)(dy phi1)  [from (b,c)=(1,2) and (2,1)]
                         *   n=1: +(dx phi2)(dy phi0) - (dx phi0)(dy phi2)
                         *   n=2: +(dx phi0)(dy phi1) - (dx phi1)(dy phi0)
                         *
                         * gradient part: d(rho)/d(d_x phi_n) = eps_{anc} phi_a d_y phi_c
                         *   n=0: eps_{a0c} phi_a dy phi_c
                         *     a=1,c=2: +phi1 dy phi2;  a=2,c=1: -phi2 dy phi1
                         *     = phi1 dy phi2 - phi2 dy phi1
                         *   Then -d_x of this contributes to force.
                         *
                         * Similarly d(rho)/d(d_y phi_n) = -eps_{abn} phi_a d_x phi_b
                         *   Then -d_y of this.
                         *
                         * FULL force = -2lcs(Q-Q0)*{direct + d_x[...] + d_y[...]}
                         *
                         * Implement numerical evaluation below.
                         */
                        /* Use the xy-plane CS density force */
                        double inv2dx_ = 1.0 / (2.0 * dx);

                        /* Direct term: eps_{nbc}(dx phi_b)(dy phi_c) */
                        double dphi_b[3], dphi_c[3];
                        /* First derivatives of phi_b and phi_c */
                        dphi_b[0] = (phi[b][IDX(i+1,j,k)] - phi[b][IDX(i-1,j,k)]) * inv2dx_;
                        dphi_b[1] = (phi[b][IDX(i,j+1,k)] - phi[b][IDX(i,j-1,k)]) * inv2dx_;
                        dphi_c[0] = (phi[c][IDX(i+1,j,k)] - phi[c][IDX(i-1,j,k)]) * inv2dx_;
                        dphi_c[1] = (phi[c][IDX(i,j+1,k)] - phi[c][IDX(i,j-1,k)]) * inv2dx_;

                        double direct = dphi_b[0]*dphi_c[1] - dphi_b[1]*dphi_c[0];

                        /* Gradient terms via numerical d_x and d_y of composite quantities */
                        /* F_x = phi_b * d_y phi_c - phi_c * d_y phi_b  at neighbors */
                        /* Force_grad_x = -d_x(F_x) */
                        double Fx_p = 0, Fx_m = 0;
                        {
                            /* At (i+1,j,k) */
                            double pb_p = phi[b][IDX(i+1,j,k)];
                            double pc_p = phi[c][IDX(i+1,j,k)];
                            double dyc_p = (phi[c][IDX(i+1,j+1,k)] - phi[c][IDX(i+1,j-1,k)]) * inv2dx_;
                            double dyb_p = (phi[b][IDX(i+1,j+1,k)] - phi[b][IDX(i+1,j-1,k)]) * inv2dx_;
                            Fx_p = pb_p * dyc_p - pc_p * dyb_p;
                        }
                        {
                            /* At (i-1,j,k) */
                            double pb_m = phi[b][IDX(i-1,j,k)];
                            double pc_m = phi[c][IDX(i-1,j,k)];
                            double dyc_m = (phi[c][IDX(i-1,j+1,k)] - phi[c][IDX(i-1,j-1,k)]) * inv2dx_;
                            double dyb_m = (phi[b][IDX(i-1,j+1,k)] - phi[b][IDX(i-1,j-1,k)]) * inv2dx_;
                            Fx_m = pb_m * dyc_m - pc_m * dyb_m;
                        }
                        double dxFx = (Fx_p - Fx_m) * inv2dx_;

                        /* F_y = -(phi_b * d_x phi_c - phi_c * d_x phi_b) at neighbors */
                        double Fy_p = 0, Fy_m = 0;
                        {
                            double pb_p = phi[b][IDX(i,j+1,k)];
                            double pc_p = phi[c][IDX(i,j+1,k)];
                            double dxc_p = (phi[c][IDX(i+1,j+1,k)] - phi[c][IDX(i-1,j+1,k)]) * inv2dx_;
                            double dxb_p = (phi[b][IDX(i+1,j+1,k)] - phi[b][IDX(i-1,j+1,k)]) * inv2dx_;
                            Fy_p = -(pb_p * dxc_p - pc_p * dxb_p);
                        }
                        {
                            double pb_m = phi[b][IDX(i,j-1,k)];
                            double pc_m = phi[c][IDX(i,j-1,k)];
                            double dxc_m = (phi[c][IDX(i+1,j-1,k)] - phi[c][IDX(i-1,j-1,k)]) * inv2dx_;
                            double dxb_m = (phi[b][IDX(i+1,j-1,k)] - phi[b][IDX(i-1,j-1,k)]) * inv2dx_;
                            Fy_m = -(pb_m * dxc_m - pc_m * dxb_m);
                        }
                        double dyFy = (Fy_p - Fy_m) * inv2dx_;

                        cs_acc = -cs_force_coeff * (direct - dxFx - dyFy);
                    }

                    acc[a][idx] = lapl - m2 * phi[a][idx] - dVdphi_triple - dVdphi_pw + cs_acc;
                }
            }
        }

        /* Boundary: zero acceleration for x,y edges */
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < Ngrid; idx++) {
            int ii = idx / (N * N);
            int jj = (idx / N) % N;
            if (ii < 2 || ii >= N-2 || jj < 2 || jj >= N-2) {
                acc[a][idx] = 0.0;
            }
        }
    }
}

/* ─── Velocity Verlet step ─── */
static void verlet_step(void)
{
    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < Ngrid; idx++)
            vel[a][idx] += 0.5 * dt * acc[a][idx];
    }

    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < Ngrid; idx++)
            phi[a][idx] += dt * vel[a][idx];
    }

    compute_acc();

    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < Ngrid; idx++)
            vel[a][idx] += 0.5 * dt * acc[a][idx];
    }

    /* Absorbing boundary damping (x,y only) */
    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < Ngrid; idx++) {
            vel[a][idx] *= damp[idx];
            phi[a][idx] *= damp[idx];
        }
    }
}

/* ─── Standard diagnostics ─── */
typedef struct {
    double Ek, Eg, Em, Ep, Epw, Et;
    double fc;
    double peak_P;
    double Pz;
} Diag;

static Diag compute_diag(double core_radius)
{
    Diag d;
    memset(&d, 0, sizeof(d));

    double Ek = 0, Eg = 0, Em = 0, Ep = 0, Epw = 0;
    double Ecore = 0, Eall = 0;
    double peak_P = 0;
    double Pz_total = 0;

    #pragma omp parallel
    {
        double lEk = 0, lEg = 0, lEm = 0, lEp = 0, lEpw = 0, lEc = 0, lEa = 0;
        double lpkP = 0;
        double lPz = 0;

        #pragma omp for schedule(static) nowait
        for (int i = 2; i < N-2; i++) {
            double x = -L + i * dx;
            for (int j = 2; j < N-2; j++) {
                double y = -L + j * dx;
                for (int k = 0; k < N; k++) {
                    long idx = IDX(i, j, k);
                    double dV = dx * dx * dx;
                    double e_loc = 0;

                    int km1 = wrap_z(k - 1);
                    int kp1 = wrap_z(k + 1);

                    for (int a = 0; a < 3; a++) {
                        double v2 = vel[a][idx] * vel[a][idx];
                        lEk += 0.5 * v2 * dV;
                        e_loc += 0.5 * v2;

                        double gx = (phi[a][IDX(i+1,j,k)] - phi[a][IDX(i-1,j,k)]) / (2*dx);
                        double gy = (phi[a][IDX(i,j+1,k)] - phi[a][IDX(i,j-1,k)]) / (2*dx);
                        double gz = (phi[a][IDX(i,j,kp1)] - phi[a][IDX(i,j,km1)]) / (2*dx);
                        double grad2 = gx*gx + gy*gy + gz*gz;
                        lEg += 0.5 * grad2 * dV;
                        e_loc += 0.5 * grad2;

                        double mass_e = 0.5 * m2 * phi[a][idx] * phi[a][idx];
                        lEm += mass_e * dV;
                        e_loc += mass_e;

                        lPz += -vel[a][idx] * gz * dV;
                    }

                    double p0 = phi[0][idx], p1 = phi[1][idx], p2 = phi[2][idx];
                    double P_val = p0 * p1 * p2;
                    double P2 = P_val * P_val;
                    double Vloc = 0.5 * mu_pot * P2 / (1.0 + kappa * P2);
                    lEp += Vloc * dV;
                    e_loc += Vloc;

                    if (lambda_pw != 0.0) {
                        double pw_e = lambda_pw * (p0*p1 + p1*p2 + p2*p0);
                        lEpw += pw_e * dV;
                        e_loc += pw_e;
                    }

                    double absP = fabs(P_val);
                    if (absP > lpkP) lpkP = absP;

                    double r = sqrt(x*x + y*y);
                    lEa += e_loc * dV;
                    if (r < core_radius) lEc += e_loc * dV;
                }
            }
        }

        #pragma omp critical
        {
            Ek += lEk; Eg += lEg; Em += lEm; Ep += lEp; Epw += lEpw;
            Ecore += lEc; Eall += lEa;
            Pz_total += lPz;
            if (lpkP > peak_P) peak_P = lpkP;
        }
    }

    d.Ek = Ek; d.Eg = Eg; d.Em = Em; d.Ep = Ep; d.Epw = Epw;
    d.Et = Ek + Eg + Em + Ep + Epw;
    d.fc = (fabs(Eall) > 1e-20) ? Ecore / Eall : 0.0;
    d.peak_P = peak_P;
    d.Pz = Pz_total;
    return d;
}

/* ═══════════════════════════════════════════════════════════════════
 * M5a: Evolve propagating braid, compute topological charges at intervals
 * ═══════════════════════════════════════════════════════════════════ */
static void run_m5a(void)
{
    printf("================================================================\n");
    printf("  M5a: Linking Number / CS Charge for Propagating Braid\n");
    printf("  m=%.1f, mu=%.1f, kappa=%.1f (M4 optimal)\n", mass, mu_pot, kappa);
    printf("================================================================\n\n");

    m2 = mass * mass;
    for (int a = 0; a < 3; a++) {
        memset(phi[a], 0, Ngrid * sizeof(double));
        memset(vel[a], 0, Ngrid * sizeof(double));
        memset(acc[a], 0, Ngrid * sizeof(double));
    }
    init_braid_propagating();
    setup_damping();

    double core_radius = 8.0;
    int Nt = (int)(tfinal / dt) + 1;
    int print_every = Nt / 20;
    if (print_every < 1) print_every = 1;

    /* Diagnostic time points */
    double diag_times[] = {0.0, 100.0, 200.0, 300.0, 400.0, 500.0};
    int n_diag_times = 6;
    int next_diag = 0;

    /* Output: detailed topology file */
    char path[600];
    snprintf(path, sizeof(path), "%s/m5a_topo_propagating.tsv", outdir);
    FILE *ftopo = fopen(path, "w");
    if (!ftopo) { fprintf(stderr, "Cannot open %s\n", path); return; }
    fprintf(ftopo, "time\tQ_jac\tQ_CS_xy\tQ_CS_yz\tQ_CS_zx\t"
                   "H_12\tH_23\tH_31\tphase_wind\t"
                   "E_total\tfc\tpeak_P\tPz\n");

    /* Fine-grained output every ~10 time units */
    int fine_every = (int)(10.0 / dt);
    if (fine_every < 1) fine_every = 1;
    snprintf(path, sizeof(path), "%s/m5a_topo_fine_prop.tsv", outdir);
    FILE *ffine = fopen(path, "w");
    if (ffine)
        fprintf(ffine, "time\tQ_jac\tQ_CS_xy\tphase_wind\tpeak_P\n");

    compute_acc();
    double wall_start = omp_get_wtime();

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        /* Check if we're at a diagnostic time */
        int do_full_diag = 0;
        if (next_diag < n_diag_times && t >= diag_times[next_diag] - 0.5*dt) {
            do_full_diag = 1;
            next_diag++;
        }
        int do_fine = (n % fine_every == 0);
        int do_print = (n % print_every == 0);

        if (do_full_diag || do_fine) {
            TopoCharges tc = compute_topo();
            Diag d = compute_diag(core_radius);

            if (do_full_diag) {
                fprintf(ftopo, "%.1f\t%.6e\t%.6e\t%.6e\t%.6e\t"
                               "%.6e\t%.6e\t%.6e\t%.4f\t"
                               "%.6e\t%.6f\t%.6e\t%.6e\n",
                        t, tc.Q_jac, tc.Q_CS_xy, tc.Q_CS_yz, tc.Q_CS_zx,
                        tc.H_12, tc.H_23, tc.H_31, tc.phase_wind,
                        d.Et, d.fc, d.peak_P, d.Pz);
                fflush(ftopo);

                printf("  t=%7.1f  Q_jac=%.4e  Q_CS_xy=%.4e  wind=%.3f  |P|=%.4f  fc=%.4f  E=%.1f\n",
                       t, tc.Q_jac, tc.Q_CS_xy, tc.phase_wind, d.peak_P, d.fc, d.Et);
            }

            if (do_fine && ffine) {
                fprintf(ffine, "%.2f\t%.6e\t%.6e\t%.4f\t%.6e\n",
                        t, tc.Q_jac, tc.Q_CS_xy, tc.phase_wind, d.peak_P);
            }
        } else if (do_print) {
            Diag d = compute_diag(core_radius);
            double elapsed = omp_get_wtime() - wall_start;
            double frac = (double)n / Nt;
            double eta_t = (frac > 0.001) ? elapsed * (1.0 - frac) / frac : 0;
            printf("  t=%7.1f  E=%.1f  fc=%.4f  |P|=%.4e  Pz=%.2f  [%.0fs, ETA %.0fs]\n",
                   t, d.Et, d.fc, d.peak_P, d.Pz, elapsed, eta_t);
            fflush(stdout);
        }

        if (n == Nt) break;
        verlet_step();
    }

    fclose(ftopo);
    if (ffine) fclose(ffine);
    printf("  M5a propagating done (%.1fs).\n\n", omp_get_wtime() - wall_start);
}

/* ═══════════════════════════════════════════════════════════════════
 * M5c: Same but for STATIC braid (no propagation) — compare topology loss
 * ═══════════════════════════════════════════════════════════════════ */
static void run_m5c_static(void)
{
    printf("================================================================\n");
    printf("  M5c-static: CS Charge for Static Braid (no propagation)\n");
    printf("  m=%.1f, mu=%.1f, kappa=%.1f\n", mass, mu_pot, kappa);
    printf("================================================================\n\n");

    m2 = mass * mass;
    for (int a = 0; a < 3; a++) {
        memset(phi[a], 0, Ngrid * sizeof(double));
        memset(vel[a], 0, Ngrid * sizeof(double));
        memset(acc[a], 0, Ngrid * sizeof(double));
    }
    init_braid_static();
    setup_damping();

    double core_radius = 8.0;
    int Nt = (int)(tfinal / dt) + 1;
    int print_every = Nt / 20;
    if (print_every < 1) print_every = 1;

    double diag_times[] = {0.0, 100.0, 200.0, 300.0, 400.0, 500.0};
    int n_diag_times = 6;
    int next_diag = 0;

    char path[600];
    snprintf(path, sizeof(path), "%s/m5c_topo_static.tsv", outdir);
    FILE *ftopo = fopen(path, "w");
    if (!ftopo) { fprintf(stderr, "Cannot open %s\n", path); return; }
    fprintf(ftopo, "time\tQ_jac\tQ_CS_xy\tQ_CS_yz\tQ_CS_zx\t"
                   "H_12\tH_23\tH_31\tphase_wind\t"
                   "E_total\tfc\tpeak_P\tPz\n");

    int fine_every = (int)(10.0 / dt);
    if (fine_every < 1) fine_every = 1;
    snprintf(path, sizeof(path), "%s/m5c_topo_fine_static.tsv", outdir);
    FILE *ffine = fopen(path, "w");
    if (ffine)
        fprintf(ffine, "time\tQ_jac\tQ_CS_xy\tphase_wind\tpeak_P\n");

    compute_acc();
    double wall_start = omp_get_wtime();

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        int do_full_diag = 0;
        if (next_diag < n_diag_times && t >= diag_times[next_diag] - 0.5*dt) {
            do_full_diag = 1;
            next_diag++;
        }
        int do_fine = (n % fine_every == 0);
        int do_print = (n % print_every == 0);

        if (do_full_diag || do_fine) {
            TopoCharges tc = compute_topo();
            Diag d = compute_diag(core_radius);

            if (do_full_diag) {
                fprintf(ftopo, "%.1f\t%.6e\t%.6e\t%.6e\t%.6e\t"
                               "%.6e\t%.6e\t%.6e\t%.4f\t"
                               "%.6e\t%.6f\t%.6e\t%.6e\n",
                        t, tc.Q_jac, tc.Q_CS_xy, tc.Q_CS_yz, tc.Q_CS_zx,
                        tc.H_12, tc.H_23, tc.H_31, tc.phase_wind,
                        d.Et, d.fc, d.peak_P, d.Pz);
                fflush(ftopo);

                printf("  t=%7.1f  Q_jac=%.4e  Q_CS_xy=%.4e  wind=%.3f  |P|=%.4f  fc=%.4f  E=%.1f\n",
                       t, tc.Q_jac, tc.Q_CS_xy, tc.phase_wind, d.peak_P, d.fc, d.Et);
            }

            if (do_fine && ffine) {
                fprintf(ffine, "%.2f\t%.6e\t%.6e\t%.4f\t%.6e\n",
                        t, tc.Q_jac, tc.Q_CS_xy, tc.phase_wind, d.peak_P);
            }
        } else if (do_print) {
            Diag d = compute_diag(core_radius);
            double elapsed = omp_get_wtime() - wall_start;
            double frac = (double)n / Nt;
            double eta_t = (frac > 0.001) ? elapsed * (1.0 - frac) / frac : 0;
            printf("  t=%7.1f  E=%.1f  fc=%.4f  |P|=%.4e  Pz=%.2f  [%.0fs, ETA %.0fs]\n",
                   t, d.Et, d.fc, d.peak_P, d.Pz, elapsed, eta_t);
            fflush(stdout);
        }

        if (n == Nt) break;
        verlet_step();
    }

    fclose(ftopo);
    if (ffine) fclose(ffine);
    printf("  M5c static done (%.1fs).\n\n", omp_get_wtime() - wall_start);
}

/* ═══════════════════════════════════════════════════════════════════
 * M5b: Propagating braid with CS potential — scan lambda_CS
 * ═══════════════════════════════════════════════════════════════════ */
static void run_m5b(double lcs_val)
{
    printf("================================================================\n");
    printf("  M5b: CS Potential, lambda_CS=%.4f\n", lcs_val);
    printf("  m=%.1f, mu=%.1f, kappa=%.1f\n", mass, mu_pot, kappa);
    printf("================================================================\n\n");

    lambda_CS = lcs_val;
    m2 = mass * mass;

    for (int a = 0; a < 3; a++) {
        memset(phi[a], 0, Ngrid * sizeof(double));
        memset(vel[a], 0, Ngrid * sizeof(double));
        memset(acc[a], 0, Ngrid * sizeof(double));
    }
    init_braid_propagating();
    setup_damping();

    /* Compute initial Q_CS for the penalty target */
    compute_acc();  /* need this before topo for clean state */
    {
        TopoCharges tc0 = compute_topo();
        /* Use Q_CS_xy as the conserved quantity for the penalty */
        Q_CS_0 = tc0.Q_CS_xy;
        printf("  Initial Q_CS_xy = %.6e (penalty target)\n", Q_CS_0);
        printf("  Initial Q_jac   = %.6e\n", tc0.Q_jac);
        printf("  Initial wind    = %.4f\n\n", tc0.phase_wind);
    }

    /* Recompute acc with CS force now that Q_CS_0 is set */
    compute_acc();

    double core_radius = 8.0;
    int Nt = (int)(tfinal / dt) + 1;
    int print_every = Nt / 20;
    if (print_every < 1) print_every = 1;

    double diag_times[] = {0.0, 100.0, 200.0, 300.0, 400.0, 500.0};
    int n_diag_times = 6;
    int next_diag = 0;

    char path[600], label[64];
    snprintf(label, sizeof(label), "lcs%.4f", lcs_val);
    snprintf(path, sizeof(path), "%s/m5b_topo_%s.tsv", outdir, label);
    FILE *ftopo = fopen(path, "w");
    if (!ftopo) { fprintf(stderr, "Cannot open %s\n", path); return; }
    fprintf(ftopo, "time\tQ_jac\tQ_CS_xy\tphase_wind\t"
                   "E_total\tfc\tpeak_P\tPz\tQ_CS_drift\n");

    int fine_every = (int)(10.0 / dt);
    if (fine_every < 1) fine_every = 1;
    snprintf(path, sizeof(path), "%s/m5b_fine_%s.tsv", outdir, label);
    FILE *ffine = fopen(path, "w");
    if (ffine)
        fprintf(ffine, "time\tQ_jac\tQ_CS_xy\tphase_wind\tpeak_P\n");

    double wall_start = omp_get_wtime();

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        int do_full_diag = 0;
        if (next_diag < n_diag_times && t >= diag_times[next_diag] - 0.5*dt) {
            do_full_diag = 1;
            next_diag++;
        }
        int do_fine = (n % fine_every == 0);
        int do_print = (n % print_every == 0);

        if (do_full_diag || do_fine) {
            TopoCharges tc = compute_topo();
            Diag d = compute_diag(core_radius);
            double drift = (fabs(Q_CS_0) > 1e-20) ?
                (tc.Q_CS_xy - Q_CS_0) / Q_CS_0 : tc.Q_CS_xy - Q_CS_0;

            if (do_full_diag) {
                fprintf(ftopo, "%.1f\t%.6e\t%.6e\t%.4f\t"
                               "%.6e\t%.6f\t%.6e\t%.6e\t%.6e\n",
                        t, tc.Q_jac, tc.Q_CS_xy, tc.phase_wind,
                        d.Et, d.fc, d.peak_P, d.Pz, drift);
                fflush(ftopo);

                printf("  t=%7.1f  Q_jac=%.4e  Q_xy=%.4e  wind=%.3f  |P|=%.4f  "
                       "drift=%.4e  E=%.1f\n",
                       t, tc.Q_jac, tc.Q_CS_xy, tc.phase_wind, d.peak_P, drift, d.Et);
            }

            if (do_fine && ffine) {
                fprintf(ffine, "%.2f\t%.6e\t%.6e\t%.4f\t%.6e\n",
                        t, tc.Q_jac, tc.Q_CS_xy, tc.phase_wind, d.peak_P);
            }
        } else if (do_print) {
            Diag d = compute_diag(core_radius);
            double elapsed = omp_get_wtime() - wall_start;
            double frac = (double)n / Nt;
            double eta_t = (frac > 0.001) ? elapsed * (1.0 - frac) / frac : 0;
            printf("  t=%7.1f  E=%.1f  fc=%.4f  |P|=%.4e  Pz=%.2f  [%.0fs, ETA %.0fs]\n",
                   t, d.Et, d.fc, d.peak_P, d.Pz, elapsed, eta_t);
            fflush(stdout);
        }

        if (n == Nt) break;
        verlet_step();
    }

    fclose(ftopo);
    if (ffine) fclose(ffine);
    printf("  M5b lambda_CS=%.4f done (%.1fs).\n\n", lcs_val, omp_get_wtime() - wall_start);

    lambda_CS = 0.0;  /* reset */
}

/* ─── Main ─── */
int main(int argc, char **argv)
{
    int test = -1;  /* -1=all, 0=M5a, 1=M5b, 2=M5c */
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-test") && i+1 < argc) { test = atoi(argv[++i]); }
        else if (!strcmp(argv[i], "-N")  && i+1 < argc) { N = atoi(argv[++i]); }
        else if (!strcmp(argv[i], "-L")  && i+1 < argc) { L = atof(argv[++i]); }
        else if (!strcmp(argv[i], "-o")  && i+1 < argc) { strncpy(outdir, argv[++i], sizeof(outdir)-1); }
        else if (!strcmp(argv[i], "-tfinal") && i+1 < argc) { tfinal = atof(argv[++i]); }
        else if (!strcmp(argv[i], "-cfl") && i+1 < argc) { cfl_frac = atof(argv[++i]); }
        else if (!strcmp(argv[i], "-mu")  && i+1 < argc) { mu_pot = atof(argv[++i]); }
        else if (!strcmp(argv[i], "-kappa") && i+1 < argc) { kappa = atof(argv[++i]); }
        else if (!strcmp(argv[i], "-mass") && i+1 < argc) { mass = atof(argv[++i]); }
    }

    dx  = 2.0 * L / (N - 1);
    dx2 = dx * dx;
    dt  = cfl_frac * dx;
    Ngrid = (long)N * N * N;

    printf("=== V27-M5: Topological Confinement via Chern-Simons ===\n");
    printf("  N=%d  L=%.1f  dx=%.4f  dt=%.5f  Ngrid=%.1fM\n",
           N, L, dx, dt, Ngrid/1e6);
    printf("  m=%.2f  mu=%.1f  kappa=%.1f\n", mass, mu_pot, kappa);
    printf("  Threads: %d\n", omp_get_max_threads());
    printf("  tfinal=%.0f  test=%d (-1=all)\n\n", tfinal, test);
    fflush(stdout);

    /* Allocate */
    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(Ngrid, sizeof(double));
        vel[a] = calloc(Ngrid, sizeof(double));
        acc[a] = calloc(Ngrid, sizeof(double));
        if (!phi[a] || !vel[a] || !acc[a]) {
            fprintf(stderr, "Allocation failed\n");
            return 1;
        }
    }
    damp = malloc(Ngrid * sizeof(double));
    if (!damp) { fprintf(stderr, "Damp allocation failed\n"); return 1; }

    /* ============================================================ */
    /* M5a: Topological charges of propagating braid (no CS force)  */
    /* ============================================================ */
    if (test == -1 || test == 0) {
        run_m5a();
    }

    /* ============================================================ */
    /* M5c: Static braid topology for comparison                     */
    /* ============================================================ */
    if (test == -1 || test == 2) {
        run_m5c_static();
    }

    /* ============================================================ */
    /* M5b: CS potential scan: lambda_CS = {0.01, 0.1, 1.0}        */
    /* ============================================================ */
    if (test == -1 || test == 1) {
        double lcs_vals[] = {0.01, 0.1, 1.0};
        int n_lcs = 3;

        char path[600];
        snprintf(path, sizeof(path), "%s/m5b_summary.tsv", outdir);
        FILE *fsum = fopen(path, "w");
        if (fsum)
            fprintf(fsum, "lambda_CS\tQ_jac_0\tQ_jac_f\tQ_CS_xy_0\tQ_CS_xy_f\t"
                          "wind_0\twind_f\tpeak_P_f\tfc_f\n");

        for (int il = 0; il < n_lcs; il++) {
            run_m5b(lcs_vals[il]);

            /* Read back final values from the fine file */
            if (fsum) {
                char fpath[600], label[64];
                snprintf(label, sizeof(label), "lcs%.4f", lcs_vals[il]);
                snprintf(fpath, sizeof(fpath), "%s/m5b_topo_%s.tsv", outdir, label);
                FILE *fr = fopen(fpath, "r");
                if (fr) {
                    char line[512];
                    double t_r, qj, qxy, pw, et, fc, pp, pz, dr;
                    double qj0 = 0, qxy0 = 0, pw0 = 0;
                    double qjf = 0, qxyf = 0, pwf = 0, ppf = 0, fcf = 0;
                    int first = 1;
                    while (fgets(line, sizeof(line), fr)) {
                        if (line[0] == 't') continue;  /* header */
                        if (sscanf(line, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
                                   &t_r, &qj, &qxy, &pw, &et, &fc, &pp, &pz, &dr) >= 7) {
                            if (first) { qj0 = qj; qxy0 = qxy; pw0 = pw; first = 0; }
                            qjf = qj; qxyf = qxy; pwf = pw; ppf = pp; fcf = fc;
                        }
                    }
                    fclose(fr);
                    fprintf(fsum, "%.4f\t%.6e\t%.6e\t%.6e\t%.6e\t%.4f\t%.4f\t%.6e\t%.6f\n",
                            lcs_vals[il], qj0, qjf, qxy0, qxyf, pw0, pwf, ppf, fcf);
                    fflush(fsum);
                }
            }
        }

        if (fsum) fclose(fsum);
    }

    /* Cleanup */
    for (int a = 0; a < 3; a++) {
        free(phi[a]); free(vel[a]); free(acc[a]);
    }
    free(damp);

    printf("\n=== V27-M5 Complete ===\n");
    return 0;
}
