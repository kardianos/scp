/*  field_qball.h — 6-field SO(3)×SO(3) Q-ball plug-in
 *
 *  Two SO(3)-valued sectors:
 *      φ = (M[0], M[1], M[2])   — "physical" triplet, rotates internally
 *      θ = (M[3], M[4], M[5])   — "torque"  triplet, static (or rotates
 *                                 in an orthogonal plane for a 2-charge
 *                                 Q-ball; toggle via init type)
 *
 *  Rotationally-invariant potential (no V(P), no determinant term):
 *      V(s, q²) = ½ m² s − ¼ a s² + ⅙ b s³ + ½ g q²
 *      s   = |φ|² + |θ|²
 *      q²  = |φ × θ|²   (only when both sectors nonzero & non-collinear)
 *
 *  Noether currents (exactly conserved by the symmetry):
 *      Q^i_φ = ε^{ijk} ∫ φ_j (∂_t φ_k) d³x
 *      Q^i_θ = ε^{ijk} ∫ θ_j (∂_t θ_k) d³x
 *
 *  Q-ball ansatz (φ rotates in xy-plane, θ static along z):
 *      φ_x(x,t) = R(r) cos(ωt),  φ_y = R(r) sin(ωt),  φ_z = 0
 *      θ_x = θ_y = 0,            θ_z = T(r)
 *      Q^z_φ = ω · ∫ R² d³x  (conserved)
 *
 *  Reference parameter set (from qball_profile.c, see qball_w0.7.dat):
 *      m²=1, a=2, b=1, g=0, ω=0.7  →  R₀=1.30, wall at r=4.24,
 *      M/Q = 0.80  (Q-ball 20% lighter than Q free quanta).
 *
 *  Build: gcc -O3 -march=native -fopenmp -o foam_sim_qball foam_sim.c \
 *             -DUSE_FIELD_QBALL -lzstd -lm
 *  (or edit the #include in foam_sim.c)
 */

/* ===== Constants (always exposed) ===== */
#ifndef FIELD_QBALL_H
#define FIELD_QBALL_H

#define NCOMP    6
#define N_SKYRME 0           /* no L_4 sector — disables Skyrme face pass */
#define FIELD_NAME "Q-ball SO(3)×SO(3) (6-comp, φ + θ)"

/* Bulk-norm metric: all spacelike. No mixed signature here. */
static const double g_metric[NCOMP] = {+1, +1, +1, +1, +1, +1};

static const char *const field_names[NCOMP] = {
    "phi_x", "phi_y", "phi_z",
    "theta_x", "theta_y", "theta_z"
};
static const uint8_t field_semantic[NCOMP] = {
    SFA_POSITION, SFA_POSITION, SFA_POSITION,
    SFA_ANGLE,    SFA_ANGLE,    SFA_ANGLE
};
static const uint8_t field_component[NCOMP] = { 0, 1, 2, 0, 1, 2 };

#endif  /* FIELD_QBALL_H */


/* ===== Implementation (requires Sim to be defined first) ===== */
#ifdef FIELD_IMPL
#undef FIELD_IMPL

/* Load qball_w0.7.dat (or whatever path) and store r, R(r), T(r) in
 * static arrays for radial interpolation. */
#define QBALL_PROFILE_PATH_DEFAULT "qball_w0.7.dat"

static int     qball_profile_loaded = 0;
static int     qball_profile_N = 0;
static double *qball_profile_r = NULL;
static double *qball_profile_R = NULL;
static double *qball_profile_T = NULL;
static double  qball_profile_r_max = 0;

static void qball_load_profile(const char *path) {
    if (qball_profile_loaded) return;
    FILE *fp = fopen(path, "r");
    if (!fp) {
        fprintf(stderr, "[qball] FATAL: cannot open profile '%s'\n", path);
        fprintf(stderr, "        run qball_profile and save to this path first.\n");
        exit(1);
    }
    int cap = 8192, N = 0;
    double *vr = malloc(sizeof(double)*cap);
    double *vR = malloc(sizeof(double)*cap);
    double *vT = malloc(sizeof(double)*cap);
    char line[1024];
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '#' || line[0] == '\n') continue;
        double r, R, Rp, T, Tp;
        if (sscanf(line, "%lf %lf %lf %lf %lf", &r, &R, &Rp, &T, &Tp) != 5) continue;
        if (N >= cap) { cap *= 2;
            vr = realloc(vr, sizeof(double)*cap);
            vR = realloc(vR, sizeof(double)*cap);
            vT = realloc(vT, sizeof(double)*cap); }
        vr[N] = r; vR[N] = R; vT[N] = T; N++;
    }
    fclose(fp);
    qball_profile_r = vr; qball_profile_R = vR; qball_profile_T = vT;
    qball_profile_N = N;
    qball_profile_r_max = vr[N-1];
    qball_profile_loaded = 1;
    printf("[qball] loaded profile from '%s': %d points, r ∈ [0, %.3f]\n",
           path, N, qball_profile_r_max);
}

/* Linear interpolation of R(r) and T(r) at arbitrary r ≥ 0. */
static inline void qball_profile_at(double r, double *R, double *T) {
    if (!qball_profile_loaded) { *R = 0; *T = 0; return; }
    if (r >= qball_profile_r_max) { *R = 0; *T = 0; return; }
    /* Profile is uniform-spaced; index directly. */
    double dr = qball_profile_r[1] - qball_profile_r[0];
    int i = (int)(r / dr);
    if (i < 0) i = 0;
    if (i >= qball_profile_N - 1) { *R = qball_profile_R[qball_profile_N-1];
                                    *T = qball_profile_T[qball_profile_N-1]; return; }
    double f = (r - qball_profile_r[i]) / dr;
    *R = (1 - f) * qball_profile_R[i] + f * qball_profile_R[i+1];
    *T = (1 - f) * qball_profile_T[i] + f * qball_profile_T[i+1];
}

/* Initial conditions ----------------------------------------------------- */

static void field_init_vacuum(Sim *s) {
    Mesh *m = s->m;
    printf("[init] vacuum: φ = θ = 0 everywhere\n");
    #pragma omp parallel for schedule(static)
    for (uint32_t i = 0; i < m->N_cells; i++) {
        for (int k = 0; k < NCOMP; k++) {
            s->M[k][i] = 0.0;
            s->M_vel[k][i] = 0.0;
        }
    }
}

/* Q-ball seed: φ rotates in xy-plane with frequency ω, θ static.
 *   φ_x(t=0) = R(r),         φ_y(t=0) = 0,           φ_z = 0
 *   ∂_t φ_x  = 0,            ∂_t φ_y  = ω·R(r)
 *   θ_z = T(r),  other θ = 0
 * Reads omega_seed for ω and the profile file from config (default path).
 */
static void field_init_qball(Sim *s) {
    Mesh *m = s->m;
    Config *c = &s->c;
    qball_load_profile(QBALL_PROFILE_PATH_DEFAULT);
    double w_phi   = c->omega_seed;
    double w_theta = c->omega_theta;
    /* Two-charge Q-ball: φ in (x,y), θ in (y,z). Different rotation
     * planes so that |φ × θ| is generically non-zero, exciting the
     * cross-coupling g·q² in V(s,q²). */
    printf("[init] qball: φ rotates at ω_φ=%.3f", w_phi);
    if (w_theta > 0.0) printf(", θ rotates at ω_θ=%.3f (in (y,z))", w_theta);
    else               printf(", θ static");
    printf(", profile from %s\n", QBALL_PROFILE_PATH_DEFAULT);
    /* When both sectors rotate, split the radial amplitude across them so
     * the total s = |φ|² + |θ|² ≈ R₀² (matching the single-sector profile
     * we solved for). Each sector then carries half the rotation kinetic
     * and half the gradient energy. With g = 0 this is an exact Q-ball
     * with two independent SO(3) charges; with g > 0 the equilibrium
     * shifts slightly and tangent damping is needed to find it. */
    double scale = (w_theta > 0.0) ? (1.0 / sqrt(2.0)) : 1.0;
    #pragma omp parallel for schedule(static)
    for (uint32_t i = 0; i < m->N_cells; i++) {
        double x = m->cell_x[i], y = m->cell_y[i], z = m->cell_z[i];
        double r = sqrt(x*x + y*y + z*z);
        double R, T;
        qball_profile_at(r, &R, &T);
        double Rs = R * scale;
        /* φ sector — rotates in (φ_x, φ_y) plane */
        s->M[0][i] = Rs;  s->M_vel[0][i] = 0.0;             /* φ_x */
        s->M[1][i] = 0;   s->M_vel[1][i] = w_phi * Rs;      /* φ_y */
        s->M[2][i] = 0;   s->M_vel[2][i] = 0.0;             /* φ_z */
        /* θ sector — rotates in (θ_y, θ_z) plane if w_theta > 0,
         * else static (T is zero unless the profile was generated
         * with T₀ ≠ 0). */
        if (w_theta > 0.0) {
            s->M[3][i] = 0;   s->M_vel[3][i] = 0.0;          /* θ_x */
            s->M[4][i] = Rs;  s->M_vel[4][i] = 0.0;          /* θ_y */
            s->M[5][i] = 0;   s->M_vel[5][i] = w_theta * Rs; /* θ_z */
        } else {
            s->M[3][i] = 0;  s->M_vel[3][i] = 0.0;
            s->M[4][i] = 0;  s->M_vel[4][i] = 0.0;
            s->M[5][i] = T;  s->M_vel[5][i] = 0.0;
        }
    }
}

static void field_init(Sim *s) {
    if      (!strcmp(s->c.init, "vacuum")) field_init_vacuum(s);
    else if (!strcmp(s->c.init, "qball"))  field_init_qball(s);
    else { fprintf(stderr, "FATAL: unknown init '%s' (vacuum|qball)\n", s->c.init); exit(1); }
}

/* Forces ----------------------------------------------------------------- */
/* EOM per component:
 *   □ M[k] = lap[k] − ∂V/∂M[k]
 *
 * V(s, q²) = ½m²s − ¼a s² + ⅙b s³ + ½g·q²
 *           with q² = |φ|²|θ|² − (φ·θ)²
 *
 * dV/d(M[k]) for k ∈ {0,1,2} (φ-component):
 *   ∂s/∂φ_k = 2 φ_k
 *   ∂(q²)/∂φ_k = 2 φ_k · |θ|² − 2 θ_k · (φ·θ)
 *   ∂V/∂φ_k = (m² − a·s + b·s² ) φ_k + g·(|θ|² φ_k − (φ·θ) θ_k)
 *
 * Symmetric for k ∈ {3,4,5} (θ-component) with φ↔θ swap on the g term.
 *
 * Effective mass terms a, b, g, mass2 are taken from config:
 *   a   ← s->c.lambda      (re-use existing knob)
 *   b   ← s->c.sigma_e2    (re-use existing knob — distinct from σ-penalty here)
 *   g   ← s->c.skyrme_c4   (re-use existing knob; cross-sector coupling)
 *   m²  ← s->c.mass2
 * We rename in the printf so the user can see what's mapped. The kernel's
 * Skyrme L_4 face pass is NOT triggered because we set N_SKYRME=0; skyrme_c4
 * is repurposed as `g` here.
 */
static void field_forces(Sim *s,
                         double *const lap[NCOMP],
                         double *acc[NCOMP]) {
    uint32_t N = s->m->N_cells;
    double m2 = s->c.mass2;
    double a  = s->c.lambda;      /* reused: quartic */
    double b  = s->c.sigma_e2;    /* reused: sextic  */
    double g  = s->c.skyrme_c4;   /* reused: cross   */

    #pragma omp parallel for schedule(static)
    for (uint32_t i = 0; i < N; i++) {
        double phx = s->M[0][i], phy = s->M[1][i], phz = s->M[2][i];
        double thx = s->M[3][i], thy = s->M[4][i], thz = s->M[5][i];
        double phi2  = phx*phx + phy*phy + phz*phz;
        double theta2 = thx*thx + thy*thy + thz*thz;
        double phith = phx*thx + phy*thy + phz*thz;
        double sums = phi2 + theta2;
        double poly = m2 - a*sums + b*sums*sums;     /* common factor */

        /* Bulk-norm diagnostic = s for this plug-in. */
        s->bulk_norm[i] = sums;

        /* ∂V/∂φ_k = poly·φ_k + g·(theta2·φ_k − phith·θ_k) */
        acc[0][i] = lap[0][i] - (poly * phx + g * (theta2 * phx - phith * thx));
        acc[1][i] = lap[1][i] - (poly * phy + g * (theta2 * phy - phith * thy));
        acc[2][i] = lap[2][i] - (poly * phz + g * (theta2 * phz - phith * thz));

        /* ∂V/∂θ_k = poly·θ_k + g·(phi2·θ_k − phith·φ_k) */
        acc[3][i] = lap[3][i] - (poly * thx + g * (phi2 * thx - phith * phx));
        acc[4][i] = lap[4][i] - (poly * thy + g * (phi2 * thy - phith * phy));
        acc[5][i] = lap[5][i] - (poly * thz + g * (phi2 * thz - phith * phz));
    }
}

/* No Skyrme L_4 sector — return 0 so the kernel skips the face pass. */
static int field_skyrme_J(Sim *s,
                          double *const grad_x[NCOMP],
                          double *const grad_y[NCOMP],
                          double *const grad_z[NCOMP],
                          double *Jx[NCOMP],
                          double *Jy[NCOMP],
                          double *Jz[NCOMP]) {
    (void)s; (void)grad_x; (void)grad_y; (void)grad_z;
    (void)Jx; (void)Jy; (void)Jz;
    return 0;
}

#endif  /* FIELD_IMPL */
