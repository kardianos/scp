/* choir_pull.c — S2: the choir's correction, derived (v89)
 *
 * Two-voice micro-model of the dense exchange channel. The pair's
 * coincident partials (q*w0 ~ p*w1) exchange energy through the kernel's
 * own gated, comb-windowed, retarded channel — but as coupled AMPLITUDES
 * (phase-coherent exchange), not as two independent directed rates.
 * kappa_freq appears NOWHERE in twin B: any exchange bias there is
 * emergent.
 *
 *   Twin A — rate-level pair: kernel transfer law verbatim, WITH the
 *            posited choir's correction (kappa_freq from the law table).
 *   Twin B — amplitude-level pair: coherent gated exchange, NO posit.
 *
 * If twin B reproduces twin A's net-flow curve and healing, the bias is
 * derived: it is the odd (reactive) component of coherent exchange that
 * rate compression discards. The negative lemma motivating this: with
 * even gates on a rung separation, a locked rate-level pair has
 * g01 == g10 identically, so NO odd-in-det flow exists at rate level at
 * any entrainment strength — the correction cannot be a phase effect.
 *
 * All constants come from the law table (argv[1]). The one bridge
 * convention: kappa = 0.5 * kd * geo * gpl * res * head  (the kernel's
 * own channel conductance read as an amplitude coupling); sweep it with
 * --bridge to expose sensitivity.
 *
 * Modes:
 *   curve <p> <q> <m>   clamped spectroscopy: sweep det along the rung,
 *                       measure net flow into voice 0 in both twins
 *   heal  <p> <q> <m> <det0>   free run: det(t) healing, both twins
 *   dscan <p> <q> <m> <det>    off-rung d sweep at fixed det: tests the
 *                              mutual-gate-closure scaling of the bias
 *
 * build: gcc -O2 -o choir_pull choir_pull.c -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define TWO_PI 6.28318530717958647692

/* ---- law table (only the keys this channel uses) ---- */
static double L_C = 1, L_r0 = 0.85, L_w2 = 2.9, L_qdet = 1.2;
static double L_gres = 0.25, L_gresm = 0.10, L_pgate = 8;
static double L_lockf = 0.005, L_kdep = 1.2, L_kdepm = 2.0, L_cap = 2.5;
static double L_es0 = 1.0, L_esfloor = 0.05, L_spull = 0.15;
static double L_klock = 1.0, L_kfreq = 0.6, L_mobfloor = 0.004;

static void load_laws(const char *path)
{
    FILE *f = fopen(path, "r");
    if (!f) { fprintf(stderr, "cannot open laws %s\n", path); exit(1); }
    char line[256];
    while (fgets(line, sizeof line, f)) {
        char *h = strchr(line, '#'); if (h) *h = 0;
        char *eq = strchr(line, '='); if (!eq) continue;
        *eq = 0;
        char *k = line, *v = eq + 1;
        while (*k == ' ' || *k == '\t') k++;
        char *e = k + strlen(k);
        while (e > k && (e[-1] == ' ' || e[-1] == '\t')) *--e = 0;
        double x = atof(v);
        if (!strcmp(k, "C")) L_C = x;
        else if (!strcmp(k, "r0")) L_r0 = x;
        else if (!strcmp(k, "w2")) L_w2 = x;
        else if (!strcmp(k, "q_detune")) L_qdet = x;
        else if (!strcmp(k, "gamma_res")) L_gres = x;
        else if (!strcmp(k, "gamma_res_m")) L_gresm = x;
        else if (!strcmp(k, "p_gate")) L_pgate = x;
        else if (!strcmp(k, "lock_floor")) L_lockf = x;
        else if (!strcmp(k, "k_dep")) L_kdep = x;
        else if (!strcmp(k, "k_dep_m")) L_kdepm = x;
        else if (!strcmp(k, "cap")) L_cap = x;
        else if (!strcmp(k, "e_s0")) L_es0 = x;
        else if (!strcmp(k, "es_floor")) L_esfloor = x;
        else if (!strcmp(k, "s_pull")) L_spull = x;
        else if (!strcmp(k, "kappa_lock")) L_klock = x;
        else if (!strcmp(k, "kappa_freq")) L_kfreq = x;
        else if (!strcmp(k, "mob_floor")) L_mobfloor = x;
    }
    fclose(f);
}

/* ---- kernel functions, verbatim semantics ---- */
static double wrap_pi(double a)
{
    a = fmod(a + M_PI, TWO_PI);
    if (a < 0) a += TWO_PI;
    return a - M_PI;
}
static double gate_of(double dphi)
{
    double g = 0.5 * (1.0 + cos(dphi));
    int ip = (int)L_pgate;
    if ((double)ip == L_pgate && ip >= 1 && ip <= 8) {
        double r = 1.0;
        for (int q = 0; q < ip; q++) r *= g;
        return r;
    }
    return pow(g, L_pgate);
}

/* channel geometry as the kernel computes it, for a seeded pair:
 * both planes transverse to the link (gpl = 1), radii from the space
 * store after the seed's pull. */
static double geo_of(double d, double x0, double x1)
{
    double geo_r[2], x[2] = { x0, x1 };
    for (int a = 0; a < 2; a++) {
        double add = x[a] * L_cap / (1.0 + L_spull);
        double pull = L_spull * add;
        double avail = L_es0 - L_esfloor;
        if (pull > avail) pull = avail > 0 ? avail : 0;
        double Es = L_es0 - pull;
        geo_r[a] = L_r0 * cbrt(Es / L_es0 > 0 ? Es / L_es0 : 0);
    }
    double ri = geo_r[0], rj = geo_r[1], A = 0;
    if (d < ri + rj) {
        double t = d * d - rj * rj + ri * ri;
        double a2 = (4.0 * d * d * ri * ri - t * t) / (4.0 * d * d);
        if (a2 > 0) A = M_PI * a2;
        else if (d < fabs(ri - rj)) { double rm = ri < rj ? ri : rj; A = M_PI * rm * rm; }
    }
    double Aref = M_PI * L_r0 * L_r0, dref = 2.0 * L_r0;
    return (A / Aref) * (dref / d);
}

/* ---- shared channel state ---- */
typedef struct {
    double E0, E1;          /* dense stocks */
    double ph0, ph1;        /* partial-channel clocks (theta2) */
    double slot01, slot10;  /* in-flight accumulators (delivery timing) */
    double t01, t10;        /* flight cycle timers in [0,1) */
    double leak;            /* |d(E0+E1)| accumulated (twin B books) */
} Pair;

typedef struct {
    int p, q;               /* the lock */
    double d, tau, geo;     /* separation, flight time, geometry factor */
    double dt, bridge;
    int clamp;              /* 1: pitches+stocks frozen at targets */
    double xt0, xt1;        /* clamp targets */
} Chan;

/* per-step observables */
typedef struct {
    double net0;            /* energy flow into voice 0, per unit time */
    double g01, g10, gg, sn, det;
} Obs;

static double pitch_of(double x) { return L_w2 / (1.0 + L_qdet * x); }

static void voices(const Chan *ch, const Pair *P, double *w0, double *w1,
                   double *det)
{
    double fl = 0.5 * (P->slot01 + P->slot10);
    double x0 = ch->clamp ? ch->xt0 : (P->E0 + fl) / L_cap;
    double x1 = ch->clamp ? ch->xt1 : (P->E1 + fl) / L_cap;
    *w0 = pitch_of(x0);
    *w1 = pitch_of(x1);
    *det = ch->q * *w0 - ch->p * *w1;
}

/* one step of either twin. amp=0: rate twin (kernel law + posited
 * kappa_freq bias). amp=1: amplitude twin (coherent exchange, no posit). */
static Obs step_pair(const Chan *ch, Pair *P, int amp)
{
    Obs ob;
    double w0, w1, det;
    voices(ch, P, &w0, &w1, &det);
    double gw = L_gresm / (ch->p * ch->q), g2 = gw * gw;
    double res = (g2 / (g2 + det * det)) / (ch->p * ch->q);

    double u = ch->q * P->ph0, v = ch->p * P->ph1;
    double a = ch->q * w0 * ch->d / L_C, b = ch->p * w1 * ch->d / L_C;
    double ps01 = wrap_pi(u - a - v);
    double ps10 = wrap_pi(v - b - u);
    double g01 = gate_of(ps01), g10 = gate_of(ps10);

    double head0 = 1.0 - P->E0 / L_cap, head1 = 1.0 - P->E1 / L_cap;
    if (head0 < 0) head0 = 0; if (head0 > 1) head0 = 1;
    if (head1 < 0) head1 = 0; if (head1 > 1) head1 = 1;
    double fl0 = L_mobfloor * L_cap;
    double E0r = P->E0 > fl0 ? P->E0 : fl0;
    double E1r = P->E1 > fl0 ? P->E1 : fl0;
    double mi = sqrt(P->E0 * E1r), mj = sqrt(P->E1 * E0r);

    double kd = L_kdep * L_kdepm;
    double base = kd * ch->geo * res;
    /* kernel wants (per unit time) — also the entrainment's flux terms */
    double F01 = base * g01 * head1 * mi;
    double F10 = base * g10 * head0 * mj;

    double dE0 = 0, dE1 = 0;
    if (!amp) {
        /* twin A: the posited law */
        double Db = 2.0 * det * gw * g2 / ((det * det + g2) * (det * det + g2));
        double kb = L_kfreq * Db * g01 * g10;
        double w01 = F01 * (1.0 - kb); if (w01 < 0) w01 = 0;
        double w10 = F10 * (1.0 + kb); if (w10 < 0) w10 = 0;
        /* outflow limiter */
        double m01 = 0.98 * P->E0 / ch->dt, m10 = 0.98 * P->E1 / ch->dt;
        if (w01 > m01) w01 = m01;
        if (w10 > m10) w10 = m10;
        dE0 = w10 - w01;
        dE1 = w01 - w10;
        ob.net0 = dE0;
    } else {
        /* twin B: coherent exchange through the same gates and window.
         * bridge convention: kappa = 0.5*base (conductance as coupling).
         * Saturation throttles as the geometric mean of the two headrooms
         * so the coupling stays hermitian (books close on the rung). */
        double S = sqrt(P->E0 * P->E1);
        double kap = ch->bridge * 0.5 * base * sqrt(head0 * head1);
        double x01 = -2.0 * kap * g01 * S * sin(ps01);  /* into 1 */
        double x10 = -2.0 * kap * g10 * S * sin(ps10);  /* into 0 */
        dE1 = x01;
        dE0 = x10;
        ob.net0 = 0.5 * (dE0 - dE1);   /* antisymmetric part = net into 0 */
        P->leak += fabs((dE0 + dE1) * ch->dt);
        /* the same coupling pulls the partial clocks (injection pulling) */
        if (P->E0 > 1e-12 && P->E1 > 1e-12) {
            double pull0 = kap * g10 * sqrt(P->E1 / P->E0) * cos(ps10);
            double pull1 = kap * g01 * sqrt(P->E0 / P->E1) * cos(ps01);
            P->ph0 += pull0 / ch->q * ch->dt;
            P->ph1 += pull1 / ch->p * ch->dt;
        }
    }

    if (!ch->clamp) {
        /* stocks move; deposits ride the channel (flight load) */
        P->E0 += dE0 * ch->dt;
        P->E1 += dE1 * ch->dt;
        if (P->E0 < 0) P->E0 = 0;
        if (P->E1 < 0) P->E1 = 0;
    }

    /* the entrainment law (kernel pass 4/5), both twins: deposits and
     * completed deliveries entrain the receiver's clock toward the
     * retarded tail. Flux terms are the kernel wants. */
    double lockf = L_lockf * L_cap;
    double f01 = F01 * ch->dt, f10 = F10 * ch->dt;
    P->slot01 += f01;
    P->slot10 += f10;
    {   /* deposit-side entrainment */
        double mix1 = f01 / (f01 + P->E1 + lockf);
        double mix0 = f10 / (f10 + P->E0 + lockf);
        P->ph1 += L_klock * mix1 * ps01 / ch->p;
        P->ph0 += L_klock * mix0 * ps10 / ch->q;
    }
    double adv = ch->dt * L_C / ch->d;
    P->t01 += adv;
    if (P->t01 >= 1.0) {          /* delivery 0->1 completes a cycle */
        double take = P->slot01;
        double mix = take / (take + P->E1 + lockf);
        P->ph1 += L_klock * mix * ps01 / ch->p;
        P->slot01 = 0;
        P->t01 -= 1.0;
    }
    P->t10 += adv;
    if (P->t10 >= 1.0) {
        double take = P->slot10;
        double mix = take / (take + P->E0 + lockf);
        P->ph0 += L_klock * mix * ps10 / ch->q;
        P->slot10 = 0;
        P->t10 -= 1.0;
    }

    P->ph0 = fmod(P->ph0 + w0 * ch->dt, TWO_PI);
    P->ph1 = fmod(P->ph1 + w1 * ch->dt, TWO_PI);

    ob.g01 = g01; ob.g10 = g10; ob.gg = g01 * g10;
    ob.sn = sin(ps01); ob.det = det;
    return ob;
}

static void seed_pair(const Chan *ch, Pair *P, double x0, double x1)
{
    memset(P, 0, sizeof *P);
    P->E0 = x0 * L_cap;
    P->E1 = x1 * L_cap;
    P->ph0 = 1.234567;
    double w0 = pitch_of(x0);
    /* kernel pair_seedlock: open the 0->1 gate at seed */
    P->ph1 = fmod((ch->q * (P->ph0 - w0 * ch->d / L_C)) / ch->p
                  + 8.0 * TWO_PI, TWO_PI);
}

/* stay exactly on the rung m while sweeping det: q*w0+p*w1 = 2*pi*m/tau */
static void rung_targets(const Chan *ch, int m, double det,
                         double *x0, double *x1)
{
    double s = TWO_PI * m / ch->tau;
    double w0 = (s + det) / (2.0 * ch->q);
    double w1 = (s - det) / (2.0 * ch->p);
    *x0 = (L_w2 / w0 - 1.0) / L_qdet;
    *x1 = (L_w2 / w1 - 1.0) / L_qdet;
}

int main(int argc, char **argv)
{
    if (argc < 3) {
        fprintf(stderr, "usage: %s laws.cfg curve|heal|dscan p q m [args] "
                        "[--d D] [--bridge B] [--dt DT]\n", argv[0]);
        return 1;
    }
    load_laws(argv[1]);
    const char *mode = argv[2];
    Chan ch;
    memset(&ch, 0, sizeof ch);
    ch.p = argc > 3 ? atoi(argv[3]) : 1;
    ch.q = argc > 4 ? atoi(argv[4]) : 1;
    int m = argc > 5 ? atoi(argv[5]) : 1;
    double det0 = argc > 6 && argv[6][0] != '-' ? atof(argv[6]) : 0.15;
    ch.d = 1.25;
    ch.dt = 0.02;
    ch.bridge = 1.0;
    for (int i = 6; i < argc; i++) {
        if (!strcmp(argv[i], "--d") && i + 1 < argc) ch.d = atof(argv[++i]);
        else if (!strcmp(argv[i], "--bridge") && i + 1 < argc) ch.bridge = atof(argv[++i]);
        else if (!strcmp(argv[i], "--dt") && i + 1 < argc) ch.dt = atof(argv[++i]);
    }
    ch.tau = ch.d / L_C;

    printf("# choir_pull laws=%s mode=%s p=%d q=%d m=%d d=%g dt=%g bridge=%g\n",
           argv[1], mode, ch.p, ch.q, m, ch.d, ch.dt, ch.bridge);
    printf("# kfreq(lawtable, twin A only)=%g gamma_res_m=%g p_gate=%g "
           "kappa_lock=%g\n", L_kfreq, L_gresm, L_pgate, L_klock);

    double gw = L_gresm / (ch.p * ch.q);

    if (!strcmp(mode, "curve")) {
        /* clamped spectroscopy along the rung */
        ch.clamp = 1;
        printf("# det  net0_A  net0_B  kb_law  net0_law  ggA  ggB  snB  "
               "lockB  x0  x1  leakB\n");
        for (int k = -40; k <= 40; k++) {
            double det = 0.01 * k;      /* +-4*Gamma_b at 1:1 */
            double x0, x1;
            rung_targets(&ch, m, det, &x0, &x1);
            if (x0 < 0.01 || x1 < 0.01 || x0 > 0.9 || x1 > 0.9) continue;
            ch.xt0 = x0; ch.xt1 = x1;
            ch.geo = geo_of(ch.d, x0, x1);

            double warm = 200, T = 2000;
            int nw = (int)(warm / ch.dt), nT = (int)(T / ch.dt);
            double aA = 0, aB = 0, ggA = 0, ggB = 0, snB = 0;
            double slipB = 0, prevps = 0, unwrap = 0;
            Pair PA, PB;
            seed_pair(&ch, &PA, x0, x1);
            seed_pair(&ch, &PB, x0, x1);
            for (int s = 0; s < nw + nT; s++) {
                Obs oA = step_pair(&ch, &PA, 0);
                Obs oB = step_pair(&ch, &PB, 1);
                double ps = wrap_pi(ch.q * PB.ph0 - ch.p * PB.ph1);
                (void)oB;
                if (s >= nw) {
                    aA += oA.net0; ggA += oA.gg;
                    aB += oB.net0; ggB += oB.gg; snB += oB.sn;
                    if (s > nw) unwrap += wrap_pi(ps - prevps);
                }
                prevps = ps;
            }
            slipB = fabs(unwrap) / TWO_PI;
            aA /= nT; aB /= nT; ggA /= nT; ggB /= nT; snB /= nT;
            /* the posited law's prediction in absolute units, using the
             * unbiased symmetric flux at this det */
            double g2 = gw * gw;
            double res = (g2 / (g2 + det * det)) / (ch.p * ch.q);
            double Db = 2.0 * det * gw * g2
                        / ((det * det + g2) * (det * det + g2));
            double kb_law = L_kfreq * Db * ggA;
            double E0 = x0 * L_cap, E1 = x1 * L_cap;
            double S = sqrt(E0 * E1);
            double kd = L_kdep * L_kdepm;
            double head = 1.0 - 0.5 * (x0 + x1);
            double Fsym = kd * ch.geo * res * sqrt(ggA > 0 ? ggA : 0)
                          * head * S;
            double net_law = 2.0 * kb_law * Fsym;
            printf("%.4f %.6e %.6e %.5f %.6e %.4f %.4f %.4f %.1f %.4f %.4f %.2e\n",
                   det, aA, aB, kb_law, net_law, ggA, ggB, snB, slipB,
                   x0, x1, PB.leak);
        }
    } else if (!strcmp(mode, "heal")) {
        ch.clamp = 0;
        double x0, x1;
        rung_targets(&ch, m, det0, &x0, &x1);
        ch.geo = geo_of(ch.d, x0, x1);
        Pair PA, PB;
        seed_pair(&ch, &PA, x0, x1);
        seed_pair(&ch, &PB, x0, x1);
        printf("# t  detA  detB  E0A  E1A  E0B  E1B  leakB\n");
        double T = 4000;
        int nT = (int)(T / ch.dt), rep = (int)(5.0 / ch.dt);
        for (int s = 0; s <= nT; s++) {
            if (s % rep == 0) {
                double w0, w1, dA, dB;
                voices(&ch, &PA, &w0, &w1, &dA);
                voices(&ch, &PB, &w0, &w1, &dB);
                printf("%.1f %.6f %.6f %.5f %.5f %.5f %.5f %.2e\n",
                       s * ch.dt, dA, dB, PA.E0, PA.E1, PB.E0, PB.E1,
                       PB.leak);
            }
            /* refresh geometry as loads move (radii follow the stores) */
            double fl = 0.5 * (PB.slot01 + PB.slot10);
            ch.geo = geo_of(ch.d, (PB.E0 + fl) / L_cap,
                            (PB.E1 + fl) / L_cap);
            step_pair(&ch, &PA, 0);
            step_pair(&ch, &PB, 1);
        }
    } else if (!strcmp(mode, "dscan")) {
        /* fixed det, separations off the rung: bias vs gate closure */
        ch.clamp = 1;
        double dr = ch.d;
        printf("# d  K  net0_B  ggB  net0_A  kb_law*2Fsym\n");
        for (int k = -12; k <= 12; k++) {
            ch.d = dr * (1.0 + 0.02 * k);
            ch.tau = ch.d / L_C;
            double x0, x1;
            /* partition chosen at the REFERENCE rung: pitches fixed,
             * separation detunes the round trip (K != 0) */
            Chan cr = ch; cr.d = dr; cr.tau = dr / L_C;
            rung_targets(&cr, m, det0, &x0, &x1);
            ch.xt0 = x0; ch.xt1 = x1;
            ch.geo = geo_of(ch.d, x0, x1);
            double w0 = pitch_of(x0), w1 = pitch_of(x1);
            double K = wrap_pi((ch.q * w0 + ch.p * w1) * ch.tau - TWO_PI * m);
            double warm = 200, T = 2000;
            int nw = (int)(warm / ch.dt), nT = (int)(T / ch.dt);
            double aB = 0, ggB = 0, aA = 0, ggA = 0;
            Pair PA, PB;
            seed_pair(&ch, &PA, x0, x1);
            seed_pair(&ch, &PB, x0, x1);
            for (int s = 0; s < nw + nT; s++) {
                Obs oA = step_pair(&ch, &PA, 0);
                Obs oB = step_pair(&ch, &PB, 1);
                if (s >= nw) { aB += oB.net0; ggB += oB.gg;
                               aA += oA.net0; ggA += oA.gg; }
            }
            aB /= nT; ggB /= nT; aA /= nT; ggA /= nT;
            printf("%.4f %.4f %.6e %.4f %.6e %.4f\n",
                   ch.d, K, aB, ggB, aA, ggA);
        }
    } else {
        fprintf(stderr, "unknown mode %s\n", mode);
        return 1;
    }
    return 0;
}
