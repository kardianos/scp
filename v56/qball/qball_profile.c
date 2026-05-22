/*  qball_profile.c — 1D Q-ball ODE solver for the SO(3)×SO(3) 6-field
 *  potential V(s, q²) = ½m²·s − ¼a·s² + ⅙b·s³ + ½g·q²
 *  with s = R² + T² and q² = R²T² (the single-sector-rotation ansatz).
 *
 *  Coupled ODE system (radial, dimensionless code units):
 *      R″ + (2/r)R′ = (m² − ω² − a·s + b·s² + g·T²) R
 *      T″ + (2/r)T′ = (m²       − a·s + b·s² + g·R²) T
 *
 *  Method: 4th-order RK4 shoot from r=0 with central values (R₀, T₀),
 *  bisect on R₀ to satisfy R(R_max)=0 (Q-ball localization). T₀ is
 *  swept on a coarser grid since for g not too large the static θ
 *  sector mostly decouples from the rotating φ profile.
 *
 *  Output: stdout columns r R R' T T', and summary stats
 *      M_Q   = ∫ [½R′² + ½T′² + V(s) + ½ω²R²] · 4πr² dr     (energy)
 *      Q_φ  = ω · ∫ R² · 4πr² dr                          (Noether charge)
 *      R_eff = √(∫ R² r² dr / ∫ R² dr)                   (RMS radius)
 *
 *  Stability check: M_Q < m·Q (Q-ball lighter than Q free quanta).
 *
 *  Build:  gcc -O2 -o qball_profile qball_profile.c -lm
 *  Usage:  ./qball_profile <m²> <a> <b> <g> <ω> <T₀>
 *          (T₀ = 0 for pure-φ Q-ball; nonzero couples in the θ sector)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Potential and its R-, T-derivatives.
 * V = ½m²s − ¼a·s² + ⅙b·s³ + ½g·R²T²
 * dV/dR² = ½m² − ½a·s + ½b·s² + ½g·T²
 * Force on R (from RHS of ODE):
 *     F_R = 2·(dV/dR²)·R − ω²·R = (m² − ω² − a·s + b·s² + g·T²) R   ← matches docstring.
 */
static inline double force_R(double R, double T,
                             double m2, double a, double b, double g, double w2) {
    double s = R*R + T*T;
    return (m2 - w2 - a*s + b*s*s + g*T*T) * R;
}
static inline double force_T(double R, double T,
                             double m2, double a, double b, double g) {
    double s = R*R + T*T;
    return (m2 - a*s + b*s*s + g*R*R) * T;
}

/* One RK4 step on the radial ODE system (R, dR/dr, T, dT/dr) over dr.
 * For r=0, the centrifugal 2/r term diverges; handle the first step
 * with a Taylor expansion R(dr) ≈ R₀ + ½·(F_R/3)·dr², where the /3
 * accounts for the 3D regularity condition. */
static void rk4_step(double *R, double *Rp, double *T, double *Tp,
                     double r, double dr,
                     double m2, double a, double b, double g, double w2) {
    double k1R, k1Rp, k1T, k1Tp;
    double k2R, k2Rp, k2T, k2Tp;
    double k3R, k3Rp, k3T, k3Tp;
    double k4R, k4Rp, k4T, k4Tp;

    #define DERIVS(R_, Rp_, T_, Tp_, r_, dR, dRp, dT, dTp) do { \
        (dR)  = (Rp_); \
        (dT)  = (Tp_); \
        (dRp) = force_R((R_), (T_), m2, a, b, g, w2) - 2.0*(Rp_)/(r_); \
        (dTp) = force_T((R_), (T_), m2, a, b, g       ) - 2.0*(Tp_)/(r_); \
    } while (0)

    DERIVS(*R,                  *Rp,                  *T,                  *Tp,                  r,           k1R, k1Rp, k1T, k1Tp);
    DERIVS(*R + 0.5*dr*k1R,     *Rp + 0.5*dr*k1Rp,    *T + 0.5*dr*k1T,     *Tp + 0.5*dr*k1Tp,    r + 0.5*dr,  k2R, k2Rp, k2T, k2Tp);
    DERIVS(*R + 0.5*dr*k2R,     *Rp + 0.5*dr*k2Rp,    *T + 0.5*dr*k2T,     *Tp + 0.5*dr*k2Tp,    r + 0.5*dr,  k3R, k3Rp, k3T, k3Tp);
    DERIVS(*R +     dr*k3R,     *Rp +     dr*k3Rp,    *T +     dr*k3T,     *Tp +     dr*k3Tp,    r +     dr,  k4R, k4Rp, k4T, k4Tp);

    *R  += dr * (k1R  + 2.0*k2R  + 2.0*k3R  + k4R ) / 6.0;
    *Rp += dr * (k1Rp + 2.0*k2Rp + 2.0*k3Rp + k4Rp) / 6.0;
    *T  += dr * (k1T  + 2.0*k2T  + 2.0*k3T  + k4T ) / 6.0;
    *Tp += dr * (k1Tp + 2.0*k2Tp + 2.0*k3Tp + k4Tp) / 6.0;

    #undef DERIVS
}

/* Shoot from r=ε with (R, R'=0, T, T'=0), step out to r_max. Return the
 * value of R at r_max (the "miss"). Q-ball solution corresponds to
 * R(r_max) → 0 from above; we bisect on R₀. */
typedef struct {
    int    N;
    double dr;
    double *r;
    double *R, *Rp;
    double *T, *Tp;
    double R_at_max;     /* terminal value of R for the converged shoot */
    double T_at_max;
    int    diverged;      /* 1 if R blew up before r_max */
} Shoot;

static void shoot_alloc(Shoot *s, int N) {
    s->N = N;
    s->r  = malloc(sizeof(double) * N);
    s->R  = malloc(sizeof(double) * N);
    s->Rp = malloc(sizeof(double) * N);
    s->T  = malloc(sizeof(double) * N);
    s->Tp = malloc(sizeof(double) * N);
}
static void shoot_free(Shoot *s) {
    free(s->r); free(s->R); free(s->Rp); free(s->T); free(s->Tp);
}

static void shoot(Shoot *s, double R0, double T0,
                  double r_max,
                  double m2, double a, double b, double g, double w2) {
    int N = s->N;
    double dr = r_max / (N - 1);
    s->dr = dr;
    double eps = dr * 1e-6;

    /* Taylor seed at r=eps: R(eps) ≈ R₀ + ½(F_R₀/3)eps², R'(eps) ≈ (F_R₀/3)eps. */
    double FR0 = force_R(R0, T0, m2, a, b, g, w2);
    double FT0 = force_T(R0, T0, m2, a, b, g);
    double R  = R0 + 0.5 * FR0/3.0 * eps*eps;
    double Rp =          FR0/3.0 * eps;
    double T  = T0 + 0.5 * FT0/3.0 * eps*eps;
    double Tp =          FT0/3.0 * eps;

    s->r[0]  = eps;
    s->R[0]  = R;  s->Rp[0] = Rp;
    s->T[0]  = T;  s->Tp[0] = Tp;
    s->diverged = 0;

    for (int i = 1; i < N; i++) {
        double r = eps + (i-1)*dr;
        rk4_step(&R, &Rp, &T, &Tp, r, dr, m2, a, b, g, w2);
        s->r[i]  = r + dr;
        s->R[i]  = R;  s->Rp[i] = Rp;
        s->T[i]  = T;  s->Tp[i] = Tp;
        if (!isfinite(R) || fabs(R) > 1e6) { s->diverged = 1; for (int j=i+1;j<N;j++) { s->r[j]=s->r[i]+(j-i)*dr; s->R[j]=R; s->Rp[j]=Rp; s->T[j]=T; s->Tp[j]=Tp; } break; }
    }
    s->R_at_max = R;
    s->T_at_max = T;
}

/* Count zero crossings of R(r) over the integrated profile. */
static int zero_crossings_R(const Shoot *s) {
    int crossings = 0;
    for (int i = 1; i < s->N; i++) {
        if ((s->R[i-1] > 0 && s->R[i] < 0) || (s->R[i-1] < 0 && s->R[i] > 0))
            crossings++;
    }
    return crossings;
}

/* Find where R first deviates significantly from R₀ (the "wall location"). */
static double wall_location(const Shoot *s, double R0) {
    double thresh = 0.5 * R0;
    for (int i = 0; i < s->N; i++) {
        if (fabs(s->R[i]) < thresh) return s->r[i];
    }
    return s->r[s->N - 1];
}

/* Compute Q-ball aggregate diagnostics. Integrates only up to r_trunc,
 * the first r where |R| starts increasing again (bisection-precision
 * limit causes the asymptotic tail to drift upward eventually).
 * Returns r_trunc for transparency. */
static double aggregates(const Shoot *s,
                       double m2, double a, double b, double g, double w,
                       double *M_Q, double *Q_phi, double *R_eff_rms,
                       double *Norm_R, double *Norm_T) {
    double dr = s->dr;
    /* Find r_trunc: scan outward, stop at the first r where |R| has
     * dropped to <1% of R₀ AND |R(r+dr)| > |R(r)| (turn-around). */
    double R0 = s->R[0];
    int i_trunc = s->N - 1;
    double thresh = 0.01 * fabs(R0);
    for (int i = 1; i < s->N - 1; i++) {
        if (fabs(s->R[i]) < thresh && fabs(s->R[i+1]) > fabs(s->R[i])) {
            i_trunc = i;
            break;
        }
    }
    double r_trunc = s->r[i_trunc];

    double E  = 0, Q  = 0;
    double r2R2 = 0, r2T2 = 0, nR2 = 0;
    for (int i = 0; i <= i_trunc; i++) {
        double r = s->r[i];
        double R = s->R[i], Rp = s->Rp[i], T = s->T[i], Tp = s->Tp[i];
        double sums = R*R + T*T;
        double V = 0.5*m2*sums - 0.25*a*sums*sums + (1.0/6.0)*b*sums*sums*sums
                 + 0.5*g*R*R*T*T;
        double e_kin_t = 0.5 * w*w * R*R;          /* ω-rotation kinetic */
        double e_kin_grad = 0.5 * (Rp*Rp + Tp*Tp);  /* gradient kinetic   */
        double ed = e_kin_t + e_kin_grad + V;
        E  += ed * 4.0*M_PI*r*r * dr;
        Q  += w*R*R * 4.0*M_PI*r*r * dr;
        r2R2 += r*r * R*R * 4.0*M_PI*r*r * dr;
        nR2  +=        R*R * 4.0*M_PI*r*r * dr;
        r2T2 += r*r * T*T * 4.0*M_PI*r*r * dr;
    }
    *M_Q = E;
    *Q_phi = Q;
    *R_eff_rms = (nR2 > 0) ? sqrt(r2R2 / nR2) : 0;
    *Norm_R = nR2;
    *Norm_T = r2T2;
    return r_trunc;
}

/* Bisect R₀ to find the Q-ball. Standard shooting criterion: the
 * critical R₀ is the largest value such that R(r) has exactly 0 zero
 * crossings. Below it, R settles at a static eq (stays > 0 forever).
 * Above it, R overshoots and crosses zero (≥ 1 crossing). The Q-ball
 * lives at the boundary — R asymptotes to 0 without crossing. */
static int bisect_R0(Shoot *s,
                     double T0, double r_max,
                     double m2, double a, double b, double g, double w2,
                     double R_lo, double R_hi,
                     int max_iters,
                     double *R0_found) {
    for (int it = 0; it < max_iters; it++) {
        double R_mid = 0.5*(R_lo + R_hi);
        shoot(s, R_mid, T0, r_max, m2, a, b, g, w2);
        int nx = s->diverged ? 99 : zero_crossings_R(s);
        if (nx > 0) {
            /* R₀ too large — overshot through zero. */
            R_hi = R_mid;
        } else {
            /* R stayed positive — either at static eq (too small) or
             * still descending (just right, but we can't tell from one
             * shoot). Push R₀ up. */
            R_lo = R_mid;
        }
        if (R_hi - R_lo < 1e-12) { *R0_found = 0.5*(R_lo+R_hi); return it; }
    }
    *R0_found = 0.5*(R_lo+R_hi);
    return max_iters;
}

/* Coarse scan of R₀ to verify there IS a transition from 0 crossings
 * to >0 crossings somewhere in [R_lo, R_hi]. Returns the (rough) R₀
 * just above the transition (1 crossing), or -1 if none found. */
static double scan_for_qball(Shoot *s,
                             double T0, double r_max,
                             double m2, double a, double b, double g, double w2,
                             double R_lo, double R_hi, int n_steps) {
    double R_prev = R_lo;
    int nx_prev = -1;
    fprintf(stderr, "[scan] R₀ vs zero-crossings + R(r_max):\n");
    for (int i = 0; i <= n_steps; i++) {
        double R0 = R_lo + (R_hi - R_lo) * (double)i / n_steps;
        shoot(s, R0, T0, r_max, m2, a, b, g, w2);
        int nx = s->diverged ? 99 : zero_crossings_R(s);
        fprintf(stderr, "        R₀=%.5f  crossings=%d  R(rmax)=%+.3e  diverged=%d\n",
                R0, nx, s->R_at_max, s->diverged);
        if (nx_prev == 0 && nx > 0 && nx < 99) {
            /* Found the transition. */
            return R_prev;
        }
        R_prev = R0;
        nx_prev = nx;
    }
    return -1;
}

int main(int argc, char **argv) {
    if (argc < 7) {
        fprintf(stderr,
            "Usage: %s <m²> <a> <b> <g> <ω> <T₀>\n"
            "  m², a, b: V = ½m²s − ¼as² + ⅙bs³ + ½g·R²T²\n"
            "  ω: rotation frequency of the φ sector (ω<m for localized solution)\n"
            "  T₀: central value of the static θ sector (0 for pure-φ Q-ball)\n",
            argv[0]);
        return 1;
    }
    double m2 = strtod(argv[1], NULL);
    double a  = strtod(argv[2], NULL);
    double b  = strtod(argv[3], NULL);
    double g  = strtod(argv[4], NULL);
    double w  = strtod(argv[5], NULL);
    double T0 = strtod(argv[6], NULL);
    double w2 = w*w;

    /* Coleman bound on ω: ω² ≥ ω_min² = 2·min_{R>0}[V(R²)/R²]. */
    double R2_crit = 0.75 * a / b;
    double V_over_R2 = 0.5*m2 - 0.25*a*R2_crit + (1.0/6.0)*b*R2_crit*R2_crit;
    double w_min2 = 2.0 * V_over_R2;
    double w_min = (w_min2 > 0) ? sqrt(w_min2) : 0;
    double m = sqrt(m2);
    fprintf(stderr, "[qball] V parameters: m²=%.4f a=%.4f b=%.4f g=%.4f\n", m2,a,b,g);
    fprintf(stderr, "[qball] Coleman window: ω ∈ [%.4f, %.4f)  (m=%.4f)\n", w_min, m, m);
    fprintf(stderr, "[qball] ω requested: %.4f   (T₀=%.4f)\n", w, T0);
    if (w >= m) { fprintf(stderr, "[qball] FATAL: ω≥m, profile won't be localised\n"); return 1; }
    if (w_min > 0 && w < w_min) { fprintf(stderr, "[qball] WARNING: ω<ω_min — no Q-ball expected\n"); }

    /* Numerical setup */
    int N = 16000;
    double r_max = 50.0;
    Shoot s; shoot_alloc(&s, N);

    /* Bisection bracketing.
     *   V_eff(R) = V(R²) − ½ω²R² has extrema at R²_± = (a ± √Δ)/(2b)
     *   with Δ = a² − 4b(m²−ω²).
     * R_- (local max of V_eff) and R_+ (local min of V_eff, deep well).
     * Q-ball R₀ lies in (R_-, R_+) — for thick wall R₀ is well below R_+;
     * for thin wall R₀ approaches R_+ from below. We bracket the
     * full range and find the critical R₀ at the 0→1 zero-crossing
     * transition. */
    double disc = a*a - 4.0*b*(m2 - w2);
    if (disc < 0) { fprintf(stderr, "[qball] FATAL: no real static R amplitude — V_eff has no minimum\n"); return 1; }
    double R2_minus = (a - sqrt(disc)) / (2.0*b);
    double R2_plus  = (a + sqrt(disc)) / (2.0*b);
    double R_lo = sqrt(fabs(R2_minus)) * 1.001;
    double R_hi = sqrt(R2_plus);
    fprintf(stderr, "[qball] R₀ bracket: [%.4f, %.4f]  (R_-=%.4f, R_+=%.4f)\n", R_lo, R_hi, sqrt(R2_minus), sqrt(R2_plus));

    /* Adaptive coarse scan: 200 linear-spaced samples across [R_lo, R_hi]
     * find the first 0→1 zero-crossing flip. Q-ball can be anywhere in
     * the bracket — thick-wall (low ω) is mid-range, thin-wall (high ω)
     * is near R_+. */
    double R_seed = -1;
    int nx_prev = -1;
    double R_prev = R_lo;
    int n_scan = 200;
    fprintf(stderr, "[scan] linear-scan R₀ vs zero-crossings (showing transitions):\n");
    for (int i = 0; i <= n_scan; i++) {
        double R0 = R_lo + (R_hi - R_lo) * (double)i / n_scan;
        shoot(&s, R0, T0, r_max, m2, a, b, g, w2);
        int nx = s.diverged ? 99 : zero_crossings_R(&s);
        if (nx != nx_prev)
            fprintf(stderr, "        R₀=%.8f  crossings=%d  R(rmax)=%+.3e\n",
                    R0, nx, s.R_at_max);
        if (nx_prev == 0 && nx > 0 && nx < 99) {
            R_seed = R_prev;
            fprintf(stderr, "[scan] transition: 0→1 crossings between R₀=%.8f and %.8f\n",
                    R_prev, R0);
            break;
        }
        R_prev = R0;
        nx_prev = nx;
    }
    if (R_seed < 0) {
        fprintf(stderr, "[qball] FATAL: no zero-crossing transition found in [%.4f, %.4f].\n", R_lo, R_hi);
        return 1;
    }

    /* Fine bisection in a TIGHT bracket — the Q-ball lives within ~1e-6 of R_seed. */
    double R_bisect_lo = R_seed;
    double R_bisect_hi = R_seed + (R_hi - R_seed) * 0.01;
    double R0_found;
    int iters = bisect_R0(&s, T0, r_max, m2, a, b, g, w2,
                          R_bisect_lo, R_bisect_hi, 100, &R0_found);
    fprintf(stderr, "[qball] bisect converged after %d iters: R₀=%.6f, R(r_max)=%.3e\n",
            iters, R0_found, s.R_at_max);

    /* Final shoot to populate s with the converged solution. */
    shoot(&s, R0_found, T0, r_max, m2, a, b, g, w2);
    fprintf(stderr, "[qball] wall location (R drops below R₀/2): r ≈ %.3f\n",
            wall_location(&s, R0_found));

    /* Diagnostics */
    double M_Q, Q_phi, R_eff, nR, nT;
    double r_trunc = aggregates(&s, m2, a, b, g, w, &M_Q, &Q_phi, &R_eff, &nR, &nT);
    fprintf(stderr, "[qball] integration truncated at r=%.3f (after the exponential tail)\n", r_trunc);
    double mQ = m * Q_phi;
    fprintf(stderr, "\n[qball] Profile summary\n");
    fprintf(stderr, "        R₀ (centre)         = %.4f\n", R0_found);
    fprintf(stderr, "        R_eff (RMS radius)  = %.4f  code units\n", R_eff);
    fprintf(stderr, "        E (Q-ball energy)   = %.4e\n", M_Q);
    fprintf(stderr, "        Q (Noether charge)  = %.4e\n", Q_phi);
    fprintf(stderr, "        E/Q                 = %.4f   (need <m=%.4f for stability)\n", M_Q/Q_phi, m);
    fprintf(stderr, "        E/(m·Q)             = %.4f   (binding fraction; <1 ⇒ bound)\n", M_Q/mQ);
    fprintf(stderr, "        ∫R² dr              = %.4f\n", nR);
    fprintf(stderr, "        ∫T² dr              = %.4f\n", nT);

    /* Print profile to stdout */
    printf("# r\tR(r)\tR'(r)\tT(r)\tT'(r)\n");
    for (int i = 0; i < N; i++) {
        if (i % 4 != 0) continue;
        printf("%.5f\t%.6e\t%.6e\t%.6e\t%.6e\n",
               s.r[i], s.R[i], s.Rp[i], s.T[i], s.Tp[i]);
    }

    shoot_free(&s);
    return 0;
}
