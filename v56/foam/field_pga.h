/*  field_pga.h — V56 PGA field definition (relativistic quaternion, 8 comp)
 *
 *  Plug-in API contract (see field_iface.md):
 *    Constants block (always when included):
 *        NCOMP                     number of field components per cell
 *        g_metric[NCOMP]           diagonal metric for the bulk norm
 *        field_names[NCOMP]        SFA column names
 *        field_semantic[NCOMP]     SFA semantic codes
 *        field_component[NCOMP]    component index in semantic group
 *        FIELD_NAME                short label printed by the kernel
 *        N_SKYRME                  number of components with Skyrme L_4
 *                                  flux (always the first N_SKYRME of M)
 *
 *    Implementation block (when included with #define FIELD_IMPL,
 *    AFTER Sim and Config structs are defined):
 *        field_init(Sim*)          initialise s->M[*], s->M_vel[*]
 *                                  by dispatching on s->c.init
 *        field_forces(Sim*, lap[NCOMP], acc[NCOMP])
 *                                  combine Laplacian + self-coupling
 *                                  (potential + sigma constraint) into
 *                                  acc. Does NOT include Skyrme L_4
 *                                  divergence — that is added by the
 *                                  kernel after a face pass on J.
 *                                  Also writes s->bulk_norm[].
 *        field_skyrme_J(Sim*, grad_x, grad_y, grad_z, Jx, Jy, Jz)
 *                                  compute the Skyrme conjugate momentum
 *                                      J^{a,b} = ∂E_4/∂(∂_b q^a)
 *                                  per cell from grad_M, for a = 0..N_SKYRME-1
 *                                  and b = 0,1,2. Returns 1 if active,
 *                                  0 to skip the divergence pass.
 *
 *  Swap field_pga.h → field_cosserat.h / field_dirac.h to change physics.
 */

/* ===== Constants (always exposed) ===== */
#ifndef FIELD_PGA_H
#define FIELD_PGA_H

#define NCOMP 8
#define N_SKYRME 4              /* rotor sector M[0..3] carries the Skyrme L_4 */
#define FIELD_NAME "PGA-rel-quat (8-comp Higgs + Skyrme L_4)"

/* Bulk-norm metric: rotor part (timelike) → −1, spatial part → +1.
 * |M|²_bulk = Σ g_kk M[k]² has Minkowski (−,−,−,−,+,+,+,+) signature. */
static const double g_metric[NCOMP] = {-1, -1, -1, -1, +1, +1, +1, +1};

static const char *const field_names[NCOMP] = {
    "M0_e410", "M1_e420", "M2_e430", "M3_one",
    "M4_e1",   "M5_e2",   "M6_e3",   "M7_e321"
};
static const uint8_t field_semantic[NCOMP] = {
    SFA_CUSTOM, SFA_CUSTOM, SFA_CUSTOM, SFA_CUSTOM,
    SFA_POSITION, SFA_POSITION, SFA_POSITION, SFA_POSITION
};
static const uint8_t field_component[NCOMP] = { 0, 1, 2, 3, 0, 1, 2, 3 };

#endif  /* FIELD_PGA_H */


/* ===== Implementation (requires Sim to be defined first) ===== */
#ifdef FIELD_IMPL
#undef FIELD_IMPL

/* Initial conditions for the multivector field. */
static void field_init_vacuum(Sim *s) {
    Mesh *m = s->m;
    double v = sqrt(s->c.v2);
    /* Vacuum on the |M|²_bulk = v² hyperboloid: pick a SPACELIKE
     * component (M[4], basis e₁) so g_kk M[k]² = +v² as desired. */
    printf("[init] vacuum: M[4] = v = %.3f, others 0\n", v);
    #pragma omp parallel for schedule(static)
    for (uint32_t i = 0; i < m->N_cells; i++) {
        for (int k = 0; k < NCOMP; k++) {
            s->M[k][i] = (k == 4) ? v : 0.0;
            s->M_vel[k][i] = 0.0;
        }
    }
}

static void field_init_qball(Sim *s) {
    Mesh *m = s->m;
    Config *c = &s->c;
    double v = sqrt(c->v2);
    double inv2R2 = 1.0 / (2.0 * c->R_seed * c->R_seed);
    /* Q-ball: vacuum on M[4] + complex (M[5] + i M[6]) bump rotating
     * as exp(−i ω t). At t=0: M[5] = A·env, M[6] = 0,
     * d/dt M[5] = 0, d/dt M[6] = ω·A·env. */
    printf("[init] qball: vacuum M[4]=%.3f + Gaussian bump A=%.3f R=%.2f ω=%.3f\n",
           v, c->A, c->R_seed, c->omega_seed);
    #pragma omp parallel for schedule(static)
    for (uint32_t i = 0; i < m->N_cells; i++) {
        double x = m->cell_x[i], y = m->cell_y[i], z = m->cell_z[i];
        double r2 = x*x + y*y + z*z;
        double env = c->A * exp(-r2 * inv2R2);
        for (int k = 0; k < NCOMP; k++) {
            s->M[k][i] = (k == 4) ? v : 0.0;
            s->M_vel[k][i] = 0.0;
        }
        s->M[5][i] += env;
        s->M_vel[6][i] += c->omega_seed * env;
    }
}

/* Skyrme hedgehog: rotor projected onto S³ with winding 1.
 *
 *     q(r) = (cos f(r), n̂ · sin f(r))     |q| = 1 everywhere
 *     f(r) = π · exp(−r/R)               f(0)=π, f(∞)=0
 *
 * Bulk-Higgs vacuum: |M|²_bulk = v² with metric (−,−,−,−,+,+,+,+) and
 * |q|² = 1 → M[4]² = v² + 1. The seed is identically on the bulk-Higgs
 * vacuum hyperboloid AND on the rotor S³ at every cell, so both
 * constraints are simultaneously satisfied at t=0. Winding number = 1.
 */
static void field_init_skyrme(Sim *s) {
    Mesh *m = s->m;
    Config *c = &s->c;
    double v2 = c->v2;
    double M4 = sqrt(v2 + 1.0);     /* |q|=1 → M[4]² = v² + |q|² */
    double R = c->R_seed;
    printf("[init] skyrme: |q|=1 hedgehog winding-1, M[4]=%.4f R=%.2f\n", M4, R);
    #pragma omp parallel for schedule(static)
    for (uint32_t i = 0; i < m->N_cells; i++) {
        double x = m->cell_x[i], y = m->cell_y[i], z = m->cell_z[i];
        double r = sqrt(x*x + y*y + z*z) + 1e-12;
        double f = M_PI * exp(-r / R);
        double cf = cos(f), sf = sin(f);
        double nx = x / r, ny = y / r, nz = z / r;
        for (int k = 0; k < NCOMP; k++) {
            s->M[k][i] = 0.0;
            s->M_vel[k][i] = 0.0;
        }
        /* Rotor sphere: (q0, q⃗) = (cos f, n̂ sin f) */
        s->M[3][i] = cf;
        s->M[0][i] = nx * sf;
        s->M[1][i] = ny * sf;
        s->M[2][i] = nz * sf;
        /* Spatial vacuum: M[4]² = v² + 1 to satisfy bulk Higgs */
        s->M[4][i] = M4;
    }
}

/* ----- Braid seeds ------------------------------------------------------
 * Two mappings of the v55 6-field Cosserat braid into the 8-component
 * multivector. Both mirror v55's Gaussian-envelope tube along z with
 * three carrier waves at phase offsets δ_a — the difference is which
 * 4 components hold the braid amplitude.
 *
 * "braid_spatial":  carriers on M[5], M[6], M[7] (three SPACELIKE comp).
 *                   Mirrors v55's φ-as-structure intuition; the rotor
 *                   part stays at zero.
 *
 * "braid_rotor":    carriers on M[0], M[1], M[2] (three TIMELIKE comp).
 *                   Exercises the spinor part of the algebra. Different
 *                   sign of g_kk → different sign of the Higgs drive,
 *                   so dynamics will differ from braid_spatial.
 *
 * Both preserve |M|²_bulk = v² at t=0 by adjusting M[4] (the spacelike
 * Higgs vacuum component) to absorb the carrier amplitude.
 */

static void field_init_braid_common(Sim *s, int rotor_mode) {
    Mesh *m = s->m;
    Config *c = &s->c;
    double v2 = c->v2;
    double L = m->L;
    double kw = M_PI / L;
    double sx = 1.0 + c->ellip;
    double sy = 1.0 - c->ellip;
    double inv2R2 = 1.0 / (2.0 * c->R_tube * c->R_tube);
    double omega = sqrt(kw*kw + (c->mass2 > 0 ? c->mass2 : 1.0));   /* group velocity */

    int comp[3];
    if (rotor_mode) {
        comp[0] = 0; comp[1] = 1; comp[2] = 2;     /* M[0], M[1], M[2] — rotor */
    } else {
        comp[0] = 5; comp[1] = 6; comp[2] = 7;     /* M[5], M[6], M[7] — spatial */
    }
    printf("[init] braid_%s: v=%.3f A=%.3f A_bg=%.3f R_tube=%.2f ellip=%.3f δ=(%.2f,%.2f,%.2f)\n",
           rotor_mode ? "rotor" : "spatial",
           sqrt(v2), c->A, c->A_bg, c->R_tube, c->ellip,
           c->delta[0], c->delta[1], c->delta[2]);
    printf("[init]   carrier components: M[%d] M[%d] M[%d]; vacuum: M[4]\n",
           comp[0], comp[1], comp[2]);

    #pragma omp parallel for schedule(static)
    for (uint32_t i = 0; i < m->N_cells; i++) {
        double x = m->cell_x[i], y = m->cell_y[i], z = m->cell_z[i];
        double r2e = x*x/(sx*sx) + y*y/(sy*sy);
        double env = exp(-r2e * inv2R2);
        for (int k = 0; k < NCOMP; k++) {
            s->M[k][i] = 0.0;
            s->M_vel[k][i] = 0.0;
        }

        /* Three carriers — same form as v55 init_braid */
        double car[3];
        for (int j = 0; j < 3; j++) {
            double ph    = kw * z + c->delta[j];
            double ph_bg = kw * z + 2.0 * M_PI * j / 3.0;
            car[j] = c->A * env * cos(ph) + c->A_bg * cos(ph_bg);
            s->M[comp[j]][i] = car[j];
            s->M_vel[comp[j]][i] = omega * c->A * env * sin(ph)
                                 + omega * c->A_bg * sin(ph_bg);
        }

        /* Adjust M[4] (the Higgs vacuum component) so |M|²_bulk = v² at t=0.
         * Bulk norm so far: rotor_mode → -Σ car_j² (timelike contribution)
         *                  spatial_mode → +Σ car_j² (spacelike contribution)
         * To reach v² we set M[4]² so that M[4]² + (signed sum) = v².
         */
        double sum_sq = car[0]*car[0] + car[1]*car[1] + car[2]*car[2];
        double m4_sq;
        if (rotor_mode) m4_sq = v2 + sum_sq;     /* timelike carriers contribute −sum_sq */
        else            m4_sq = v2 - sum_sq;     /* spacelike carriers contribute +sum_sq */
        s->M[4][i] = (m4_sq > 0) ? sqrt(m4_sq) : 0.0;
    }
}

static void field_init_braid_spatial(Sim *s) { field_init_braid_common(s, 0); }
static void field_init_braid_rotor(Sim *s)   { field_init_braid_common(s, 1); }

static void field_init(Sim *s) {
    if      (!strcmp(s->c.init, "vacuum"))         field_init_vacuum(s);
    else if (!strcmp(s->c.init, "qball"))          field_init_qball(s);
    else if (!strcmp(s->c.init, "skyrme"))         field_init_skyrme(s);
    else if (!strcmp(s->c.init, "braid_spatial"))  field_init_braid_spatial(s);
    else if (!strcmp(s->c.init, "braid_rotor"))    field_init_braid_rotor(s);
    else { fprintf(stderr, "FATAL: unknown init '%s' (vacuum|qball|skyrme|braid_spatial|braid_rotor)\n", s->c.init); exit(1); }
}

/* Field-specific local force assembly. For PGA stage A + Stage C (Skyrme):
 *
 *     □ M[k] = lap[k] − m² M[k] − λ (|M|²_bulk − v²) g_kk M[k]
 *              − (k<4) · sigma_e2 (|q|² − 1) M[k]
 *              + (k<4) · ∂_b J^{k,b}     ← added by kernel via face pass
 *              + (k<4) · F_π^k           ← pion mass force, see below
 *
 * Pion mass: V_π = m_π² (1 − q₀), q₀ = M[3]. Tangent-projected force
 * on the rotor S³:
 *     F_π^k = m_π² (δ_{k3} − q₀ M[k])    for k = 0..3
 * (k=3: m_π² (1 − q₀²); k=0..2: −m_π² q₀ M[k]).
 * This penalizes deviations from the q₀=1 vacuum and gives the
 * Skyrmion a finite preferred size (Compton wavelength of the pion).
 *
 * The acceleration arrays are written; bulk_norm is also updated for
 * downstream diagnostics. The Skyrme L_4 divergence is added by the
 * kernel after this routine returns; here we only cover the local
 * (no-derivative) drives. */
static void field_forces(Sim *s,
                         double *const lap[NCOMP],
                         double *acc[NCOMP]) {
    uint32_t N = s->m->N_cells;
    double mass2 = s->c.mass2;
    double lambda = s->c.lambda;
    double v2 = s->c.v2;
    double sigma_e2 = s->c.sigma_e2;
    double m_pi2 = s->c.m_pi2;
    #pragma omp parallel for schedule(static)
    for (uint32_t i = 0; i < N; i++) {
        double Mb = 0;
        for (int k = 0; k < NCOMP; k++)
            Mb += g_metric[k] * s->M[k][i] * s->M[k][i];
        s->bulk_norm[i] = Mb;
        double drive = lambda * (Mb - v2);
        double q2 = 0.0;
        for (int j = 0; j < N_SKYRME; j++) q2 += s->M[j][i] * s->M[j][i];
        double q_drive = sigma_e2 * (q2 - 1.0);
        double q0 = s->M[3][i];
        for (int k = 0; k < NCOMP; k++) {
            double a = lap[k][i] - mass2 * s->M[k][i] - drive * g_metric[k] * s->M[k][i];
            if (k < N_SKYRME) a -= q_drive * s->M[k][i];
            if (m_pi2 > 0 && k < N_SKYRME) {
                /* Tangent-projected pion-mass force on S³. */
                if (k == 3) a += m_pi2 * (1.0 - q0 * q0);
                else        a -= m_pi2 * q0 * s->M[k][i];
            }
            acc[k][i] = a;
        }
    }
}

/* Skyrme conjugate momenta J^{a,b} = ∂E_4/∂(∂_b q^a) on cells.
 *
 * Energy density (rotor only, a = 0..N_SKYRME-1):
 *     e_4 = (c4/2) [ (Tr g)² − ‖g‖²_F ]
 *     g_ij = Σ_a (∂_i q^a)(∂_j q^a)            (3×3 symmetric)
 *
 * Derivative:
 *     J^{a,b} = c4 [ Tr(g)·(∂_b q^a) − Σ_j g_bj·(∂_j q^a) ]
 *
 * The kernel then computes ∂_b J^{a,b} via a face pass and adds it to
 * acc[a]. Energy-conservative: J on cells + finite-volume divergence
 * is the discrete adjoint of the cell-centered gradient operator that
 * defined g_ij in the first place (summation by parts).
 *
 * Returns 1 if Skyrme is active (c4 > 0), 0 to skip the divergence
 * pass entirely. */
static int field_skyrme_J(Sim *s,
                          double *const grad_x[NCOMP],
                          double *const grad_y[NCOMP],
                          double *const grad_z[NCOMP],
                          double *Jx[NCOMP],
                          double *Jy[NCOMP],
                          double *Jz[NCOMP]) {
    double c4 = s->c.skyrme_c4;
    if (c4 <= 0.0) return 0;
    uint32_t N = s->m->N_cells;
    #pragma omp parallel for schedule(static)
    for (uint32_t i = 0; i < N; i++) {
        /* Build 3×3 symmetric metric on the rotor gradients. */
        double gxx = 0, gyy = 0, gzz = 0, gxy = 0, gxz = 0, gyz = 0;
        for (int a = 0; a < N_SKYRME; a++) {
            double dx = grad_x[a][i], dy = grad_y[a][i], dz = grad_z[a][i];
            gxx += dx * dx;
            gyy += dy * dy;
            gzz += dz * dz;
            gxy += dx * dy;
            gxz += dx * dz;
            gyz += dy * dz;
        }
        double trg = gxx + gyy + gzz;
        double coef = 2.0 * c4;   /* J = 2·c4·[Tr(g)·∂q − g·∂q] from differentiating
                                   * e_4 = (c4/2)·[(Tr g)² − ‖g‖²_F]; the factor 2 comes
                                   * from chain rule on the squared trace. */
        for (int a = 0; a < N_SKYRME; a++) {
            double dx = grad_x[a][i], dy = grad_y[a][i], dz = grad_z[a][i];
            Jx[a][i] = coef * (trg * dx - (gxx * dx + gxy * dy + gxz * dz));
            Jy[a][i] = coef * (trg * dy - (gxy * dx + gyy * dy + gyz * dz));
            Jz[a][i] = coef * (trg * dz - (gxz * dx + gyz * dy + gzz * dz));
        }
    }
    return 1;
}

#endif  /* FIELD_IMPL */
