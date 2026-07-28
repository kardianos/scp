/*  fabric_link.c — v88: energy transfer MEDIATED BY LINK HARMONICS at rate c.
 *
 *  WHAT THIS FIXES
 *    The previous instrument (fabric_harmonic.c) moved energy with a bilinear
 *    cell-to-cell hop. Measured consequence (GROK_V88_SYMBOLIC.md M2):
 *        instantaneous transfer A_sec = 2*eps1*I = 0.075 EXACTLY,
 *        independent of detuning, CV = 0.000, flat fit.
 *    Detuning suppressed only the TIME AVERAGE; on-resonance/far ratio was 1.65,
 *    not a gate. So "energy crosses only in complete cycles" was not implemented
 *    at all -- with a direct cell-cell term the energy is already there the
 *    instant you write the Hamiltonian down, and detuning can only cancel it
 *    after the fact.
 *
 *  THE CORRECTED STRUCTURE
 *    Energy does not hop. It is MEDIATED. Every link (i,j) carries its own
 *    harmonic oscillator u_ij. Cells couple to links, never directly to each
 *    other:
 *
 *        cell i  <--g-->  link (i,j)  <--g-->  cell j
 *
 *    Energy must therefore be deposited into the link, the link must carry it
 *    through its own cycle, and only then can it be delivered. Transfer takes
 *    TIME, and it is gated STRUCTURALLY: a link whose frequency Omega cannot be
 *    driven resonantly by the cell modes never accepts the energy in the first
 *    place. There is no instantaneous channel left to average away.
 *
 *  WHERE c COMES FROM  [D, given the structure]
 *    A signal must traverse cell -> link -> cell, and each step is limited by
 *    how fast the coupling can complete a cycle. So the propagation speed is
 *        c ~ a * (rate at which link coupling completes)
 *    an EMERGENT constant set by g and Omega, not an input. c is a limit
 *    because that is the fastest cycles can couple, which is exactly the
 *    ontology's reading. This file measures c and checks it scales as the
 *    structure predicts.
 *
 *  MEASUREMENTS
 *    T1  propagation speed of a pulse, and its scaling with g and Omega.
 *        A well-defined finite c that scales as predicted is the first check.
 *    T2  transfer vs detuning: is the gate now STRUCTURAL (sharp) rather than
 *        the 1.65 ratio a bilinear hop gave? This is the decisive one.
 *    T3  two-lump static interaction: attractive or repulsive? The bilinear
 *        model gave repulsion (dE = +3.59 at sep 2, R^2 = 0.996) and therefore
 *        no bound object.
 *
 *  Build: gcc -O3 -march=native -fopenmp -o fabric_link fabric_link.c -lm
 *  Usage: fabric_link [test]        test = speed | gate | pair | all
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define LMAX 4096

/* 1D chain is enough for the structural questions; geometry comes later. */
static int NC = 512;                 /* cells */
static double OMC = 1.00;            /* cell mode frequency */
static double OMLK = 1.00;           /* link harmonic frequency  <- sets c */
static double G = 0.25;              /* cell<->link coupling     <- sets c */
static double DT = 0.002;

/* state: cell amplitudes z[i], link amplitudes u[i] on link (i,i+1) */
static double zr[LMAX], zi[LMAX], ur[LMAX], ui[LMAX];
static double dzr[LMAX], dzi[LMAX], dur[LMAX], dui[LMAX];
static double det[LMAX];             /* per-cell detuning of the mode */

/* H = sum_i (omc+det_i)|z_i|^2 + sum_i omlk|u_i|^2
 *   + g sum_i [ z_i* u_i + u_i* z_{i+1} + c.c. ]
 * i dz_i/dt = (omc+det_i) z_i + g( u_i + u_{i-1} )
 * i du_i/dt = omlk u_i     + g( z_i + z_{i+1} )
 * Cells are coupled ONLY through links. There is no z_i <-> z_j term. */
static void deriv(void) {
    for (int i = 0; i < NC; i++) {
        double sr = 0, si = 0;
        if (i < NC - 1) { sr += ur[i];   si += ui[i]; }
        if (i > 0)      { sr += ur[i-1]; si += ui[i-1]; }
        double hr = (OMC + det[i]) * zr[i] + G * sr;
        double hi = (OMC + det[i]) * zi[i] + G * si;
        dzr[i] =  hi;      /* dz/dt = -i H z  ->  re' = +im(H z) */
        dzi[i] = -hr;
    }
    for (int i = 0; i < NC - 1; i++) {
        double sr = zr[i] + zr[i+1], si = zi[i] + zi[i+1];
        double hr = OMLK * ur[i] + G * sr;
        double hi = OMLK * ui[i] + G * si;
        dur[i] =  hi;
        dui[i] = -hr;
    }
}

static void step(void) {
    static double z0r[LMAX], z0i[LMAX], u0r[LMAX], u0i[LMAX];
    memcpy(z0r, zr, sizeof(double)*NC); memcpy(z0i, zi, sizeof(double)*NC);
    memcpy(u0r, ur, sizeof(double)*NC); memcpy(u0i, ui, sizeof(double)*NC);
    deriv();
    for (int i = 0; i < NC; i++) { zr[i] = z0r[i] + 0.5*DT*dzr[i];
                                   zi[i] = z0i[i] + 0.5*DT*dzi[i];
                                   ur[i] = u0r[i] + 0.5*DT*dur[i];
                                   ui[i] = u0i[i] + 0.5*DT*dui[i]; }
    deriv();
    for (int i = 0; i < NC; i++) { zr[i] = z0r[i] + DT*dzr[i];
                                   zi[i] = z0i[i] + DT*dzi[i];
                                   ur[i] = u0r[i] + DT*dur[i];
                                   ui[i] = u0i[i] + DT*dui[i]; }
}

static void clear(void) {
    memset(zr,0,sizeof zr); memset(zi,0,sizeof zi);
    memset(ur,0,sizeof ur); memset(ui,0,sizeof ui);
    memset(det,0,sizeof det);
}

/* ---- T1: propagation speed. Excite one cell, track the leading edge. ---- */
static double measure_speed(double g, double omlk) {
    G = g; OMLK = omlk;
    clear();
    int src = NC/2;
    zr[src] = 1.0;
    double thr = 1e-4;
    double t = 0, tprev = 0; int rprev = 0;
    double vsum = 0; int vn = 0;
    for (int s = 0; s < 60000; s++) {
        step(); t += DT;
        if (s % 200) continue;
        int reach = 0;
        for (int i = src+1; i < NC; i++) {
            double a = zr[i]*zr[i] + zi[i]*zi[i];
            if (a > thr) reach = i - src; else break;
        }
        if (reach > rprev && rprev > 4 && reach < NC/2 - 8) {
            vsum += (double)(reach - rprev) / (t - tprev); vn++;
        }
        if (reach > rprev) { rprev = reach; tprev = t; }
        if (reach >= NC/2 - 8) break;
    }
    return vn ? vsum / vn : 0.0;
}

/* ---- T2: transfer vs detuning. Two cells, one link. THE decisive test. ---- */
static double transfer_at_detune(double delta) {
    clear();
    NC = 2;
    det[1] = delta;                       /* cell 1 detuned from cell 0 */
    zr[0] = 1.0;
    double maxI1 = 0;
    for (int s = 0; s < 200000; s++) {
        step();
        double a = zr[1]*zr[1] + zi[1]*zi[1];
        if (a > maxI1) maxI1 = a;
    }
    NC = 512;
    return maxI1;                          /* peak fraction delivered */
}

/* ---- T3: two localized excitations, static energy vs separation ---- */
static double pair_energy(int sep) {
    clear();
    int a = NC/2 - sep/2, b = a + sep;
    for (int k = -2; k <= 2; k++) {
        double w = exp(-0.5*k*k);
        zr[a+k] += w; zr[b+k] += w;
    }
    /* H expectation (real initial data) */
    deriv();
    double E = 0;
    for (int i = 0; i < NC; i++) {
        E += (OMC + det[i]) * (zr[i]*zr[i] + zi[i]*zi[i]);
        E += OMLK * (ur[i]*ur[i] + ui[i]*ui[i]);
    }
    for (int i = 0; i < NC-1; i++)
        E += 2.0 * G * (ur[i]*(zr[i]+zr[i+1]) + ui[i]*(zi[i]+zi[i+1]));
    return E;
}

int main(int argc, char **argv) {
    const char *test = (argc > 1) ? argv[1] : "all";
    printf("=====================================================================\n");
    printf("v88 link-mediated transfer: energy crosses via a link harmonic\n");
    printf("=====================================================================\n");
    printf("  cells couple ONLY to links, never to each other.\n");
    printf("  cell i <--g--> link(i,i+1) <--g--> cell i+1\n");
    printf("  transfer takes time; c is emergent from g and Omega_link.\n\n");

    if (!strcmp(test,"speed") || !strcmp(test,"all")) {
        printf("  T1 PROPAGATION SPEED (cells per unit time)\n");
        printf("  %8s %10s %12s %14s\n", "g", "Omega", "measured c", "c/g");
        for (double g = 0.10; g <= 0.85; g *= 1.6)
            for (double om = 1.0; om <= 2.1; om += 1.0) {
                double c = measure_speed(g, om);
                printf("  %8.3f %10.3f %12.4f %14.4f\n", g, om, c, g>0?c/g:0);
            }
        printf("  a finite, well-defined c that scales with g is the first check:\n");
        printf("  it means transfer is rate-limited by cycle coupling, not instant.\n\n");
    }

    if (!strcmp(test,"gate") || !strcmp(test,"all")) {
        G = 0.25; OMLK = 1.0;
        printf("  T2 TRANSFER vs DETUNING  -- the decisive test\n");
        printf("  bilinear hop gave on/off ratio 1.65 (no gate, CV=0.000 flat).\n");
        printf("  %10s %16s %12s\n", "detune", "peak I delivered", "vs zero");
        double base = transfer_at_detune(0.0);
        for (double d = 0.0; d <= 4.01; d += 0.25) {
            double v = transfer_at_detune(d);
            printf("  %10.3f %16.6e %12.4f\n", d, v, base>0 ? v/base : 0);
        }
        printf("  a SHARP fall means the gate is structural: a link that cannot\n");
        printf("  be driven resonantly never accepts the energy at all.\n\n");
    }

    if (!strcmp(test,"pair") || !strcmp(test,"all")) {
        G = 0.25; OMLK = 1.0;
        printf("  T3 TWO-LUMP STATIC ENERGY vs SEPARATION\n");
        printf("  bilinear model gave dE = +3.59 at sep 2 -> REPULSIVE, no bound object.\n");
        printf("  %8s %16s %14s\n", "sep", "E(sep)", "dE vs far");
        double efar = pair_energy(64);
        for (int s = 2; s <= 24; s += 2)
            printf("  %8d %16.6f %14.6f\n", s, pair_energy(s), pair_energy(s) - efar);
        printf("  NEGATIVE dE at small separation = attraction = a bound object\n");
        printf("  is possible. Positive = the link model inherits the same defect.\n");
    }
    return 0;
}
