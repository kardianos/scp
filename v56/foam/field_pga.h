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
 *
 *    Implementation block (when included with #define FIELD_IMPL,
 *    AFTER Sim and Config structs are defined):
 *        field_init(Sim*)          initialise s->M[*], s->M_vel[*]
 *                                  by dispatching on s->c.init
 *        field_forces(Sim*, lap[NCOMP], acc[NCOMP])
 *                                  combine Laplacian + self-coupling
 *                                  into the per-component acceleration.
 *                                  May also write s->bulk_norm[].
 *
 *  Swap field_pga.h → field_cosserat.h / field_dirac.h to change physics.
 */

/* ===== Constants (always exposed) ===== */
#ifndef FIELD_PGA_H
#define FIELD_PGA_H

#define NCOMP 8
#define FIELD_NAME "PGA-rel-quat (8-comp Higgs)"

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

static void field_init_skyrme(Sim *s) {
    Mesh *m = s->m;
    Config *c = &s->c;
    double v = sqrt(c->v2);
    double v_q = c->A;
    double R = c->R_seed;
    printf("[init] skyrme: vacuum M[4]=%.3f + rotor hedgehog v_q=%.3f R=%.2f\n",
           v, v_q, R);
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
        s->M[3][i] = v_q * cf;
        s->M[0][i] = v_q * nx * sf;
        s->M[1][i] = v_q * ny * sf;
        s->M[2][i] = v_q * nz * sf;
        s->M[4][i] = v;
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

/* Field-specific force assembly: combines the kernel's pre-computed
 * Laplacian with this field's self-coupling. For PGA stage A:
 *
 *     □ M[k] = lap[k] − m² M[k] − λ (|M|²_bulk − v²) g_kk M[k]
 *
 * The acceleration arrays are written; bulk_norm is also updated for
 * downstream diagnostics. */
static void field_forces(Sim *s,
                         double *const lap[NCOMP],
                         double *acc[NCOMP]) {
    uint32_t N = s->m->N_cells;
    double mass2 = s->c.mass2;
    double lambda = s->c.lambda;
    double v2 = s->c.v2;
    #pragma omp parallel for schedule(static)
    for (uint32_t i = 0; i < N; i++) {
        double Mb = 0;
        for (int k = 0; k < NCOMP; k++)
            Mb += g_metric[k] * s->M[k][i] * s->M[k][i];
        s->bulk_norm[i] = Mb;
        double drive = lambda * (Mb - v2);
        for (int k = 0; k < NCOMP; k++) {
            acc[k][i] = lap[k][i]
                      - mass2 * s->M[k][i]
                      - drive * g_metric[k] * s->M[k][i];
        }
    }
}

#endif  /* FIELD_IMPL */
