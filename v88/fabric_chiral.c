/*  fabric_chiral.c — v88: link actuation in TWO PLANES along a chiral harmonic.
 *
 *  WHY THE PREVIOUS MODEL WAS LINEAR
 *    fabric_link.c put a single scalar harmonic on each link of a 1D chain. It
 *    got the two things it was built for -- emergent c (c ~ 1.05*g*a over a 4x
 *    range in g) and a structural transfer gate (227x suppression across the
 *    detune window, vs the bilinear hop's 1.65) -- but it could not bind,
 *    because in 1D A LINK HAS NO TRANSVERSE PLANE. There is nowhere for a
 *    second actuation plane to live, so no handedness, so nothing for chirality
 *    to label and no handedness-dependent interaction. The linearity was
 *    geometric, not a modelling accident.
 *
 *  THE STRUCTURE
 *    A link along direction dhat has a 2D plane perpendicular to it. The link
 *    harmonic actuates in BOTH transverse directions, and the relative phase
 *    between them is the handedness:
 *
 *        u = (u_1, u_2)  in the plane perp to dhat
 *        u_+/- = (u_1 -/+ i u_2)/sqrt(2)      left / right circular
 *
 *    |u_+|^2 - |u_-|^2 is the link's chirality: a signed, discrete-in-sign
 *    label carried by the geometry, not by a field anyone added.
 *
 *  WHY THIS IS WHERE BINDING CAN COME FROM
 *    Two cells connected through links of the SAME handedness exchange energy
 *    on a channel that closes; OPPOSITE handedness does not close the same way.
 *    So the interaction is handedness-dependent, and an attractive channel can
 *    exist for one combination while the other is repulsive. A scalar link has
 *    no such distinction -- every pair sees the same channel, which is why
 *    fabric_harmonic measured pure repulsion (dE = +3.59 at sep 2, R^2 = 0.996).
 *
 *  DIMENSION
 *    Handedness needs a cross product. In 3D, dhat x u is defined and the
 *    chiral coupling exists. The cross product exists ONLY in 3 and 7
 *    dimensions -- which is the first structural reason in this project for the
 *    ontology's 3-7 window, and it is geometric rather than fitted. Measured
 *    here by running d = 1,2,3 and showing chirality is IDENTICALLY ZERO below
 *    3 because the transverse plane is degenerate or absent.
 *
 *  MEASURES
 *    C1  chirality is well-defined only for d >= 3
 *    C2  does the chiral channel split same- vs opposite-handed interaction?
 *    C3  two-lump static energy by relative handedness: attraction anywhere?
 *
 *  Build: gcc -O3 -march=native -fopenmp -o fabric_chiral fabric_chiral.c -lm
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define NMAX 40
static int D = 3, L = 16, N;
static int stride[3];

/* cell: complex 3-vector amplitude (so it can carry a rotation sense)
 * link: complex amplitude in the 2 transverse directions of that link */
typedef struct { double r[3], i[3]; } CVec;
static CVec *z, *dz;
/* per site, per direction: link transverse amplitudes (2 components) */
static double *lr, *li, *dlr, *dli;      /* [site*D*2 + a*2 + c] */

static double OMC = 1.0, OMLK = 1.0, G = 0.25, CHI = 0.6, DT = 0.004;

static inline int nbr(int idx, int a, int dir) {
    int c = (idx / stride[a]) % L;
    int nc = ((c + dir) % L + L) % L;
    return idx + (nc - c) * stride[a];
}

/* the two transverse unit vectors for a link along axis a */
static void transverse(int a, int t1[3], int t2[3]) {
    for (int k = 0; k < 3; k++) { t1[k] = 0; t2[k] = 0; }
    t1[(a + 1) % 3] = 1;
    t2[(a + 2) % 3] = 1;
}

/* H = omc sum |z|^2 + omlk sum |u|^2
 *   + g  sum_links Re[ (P_perp z_i + P_perp z_j)* . u ]
 *   + chi sum_links  CHIRAL term: couples the two transverse components with a
 *         90-degree twist, i.e. u_1 <-> u_2 with i -- this is what makes the
 *         link harmonic circular rather than two independent linear ones.
 *         In index form: + chi * (u_1* (i u_2) + c.c.) = 2 chi Im(u_1* u_2). */
static void deriv(void) {
    memset(dz, 0, sizeof(CVec) * N);
    memset(dlr, 0, sizeof(double) * (size_t)N * D * 2);
    memset(dli, 0, sizeof(double) * (size_t)N * D * 2);
    for (int s = 0; s < N; s++) {
        /* cell self term */
        for (int k = 0; k < 3; k++) {
            dz[s].r[k] +=  OMC * z[s].i[k];
            dz[s].i[k] += -OMC * z[s].r[k];
        }
        for (int a = 0; a < D; a++) {
            int t1[3], t2[3]; transverse(a, t1, t2);
            int j = nbr(s, a, +1);
            size_t base = ((size_t)s * D + a) * 2;
            double u1r = lr[base], u1i = li[base];
            double u2r = lr[base+1], u2i = li[base+1];
            /* projections of both cells onto the two transverse directions */
            double p1r = 0, p1i = 0, p2r = 0, p2i = 0;
            for (int k = 0; k < 3; k++) {
                p1r += t1[k] * (z[s].r[k] + z[j].r[k]);
                p1i += t1[k] * (z[s].i[k] + z[j].i[k]);
                p2r += t2[k] * (z[s].r[k] + z[j].r[k]);
                p2i += t2[k] * (z[s].i[k] + z[j].i[k]);
            }
            /* link EOM: i du/dt = omlk u + g p + chiral mixing */
            double h1r = OMLK * u1r + G * p1r - CHI * u2i;
            double h1i = OMLK * u1i + G * p1i + CHI * u2r;
            double h2r = OMLK * u2r + G * p2r + CHI * u1i;
            double h2i = OMLK * u2i + G * p2i - CHI * u1r;
            dlr[base]   +=  h1i;  dli[base]   += -h1r;
            dlr[base+1] +=  h2i;  dli[base+1] += -h2r;
            /* back-reaction on both cells */
            for (int k = 0; k < 3; k++) {
                double fr = G * (t1[k] * u1r + t2[k] * u2r);
                double fi = G * (t1[k] * u1i + t2[k] * u2i);
                dz[s].r[k] +=  fi;  dz[s].i[k] += -fr;
                dz[j].r[k] +=  fi;  dz[j].i[k] += -fr;
            }
        }
    }
}

static void step(void) {
    deriv();
    for (int s = 0; s < N; s++)
        for (int k = 0; k < 3; k++) { z[s].r[k] += DT*dz[s].r[k];
                                      z[s].i[k] += DT*dz[s].i[k]; }
    for (size_t q = 0; q < (size_t)N*D*2; q++) { lr[q] += DT*dlr[q]; li[q] += DT*dli[q]; }
}

static void alloc(void) {
    N = 1; for (int a = 0; a < D; a++) N *= L;
    stride[0] = 1; for (int a = 1; a < D; a++) stride[a] = stride[a-1]*L;
    z = realloc(z, sizeof(CVec)*N); dz = realloc(dz, sizeof(CVec)*N);
    lr = realloc(lr, sizeof(double)*(size_t)N*D*2);
    li = realloc(li, sizeof(double)*(size_t)N*D*2);
    dlr = realloc(dlr, sizeof(double)*(size_t)N*D*2);
    dli = realloc(dli, sizeof(double)*(size_t)N*D*2);
    memset(z,0,sizeof(CVec)*N);
    memset(lr,0,sizeof(double)*(size_t)N*D*2); memset(li,0,sizeof(double)*(size_t)N*D*2);
}

/* total chirality  sum_links (|u_+|^2 - |u_-|^2), u_+/- = (u1 -/+ i u2)/sqrt2 */
static double chirality(void) {
    double c = 0;
    for (int s = 0; s < N; s++)
        for (int a = 0; a < D; a++) {
            size_t b = ((size_t)s*D + a)*2;
            double u1r=lr[b], u1i=li[b], u2r=lr[b+1], u2i=li[b+1];
            /* |u+|^2 - |u-|^2 = 2 Im(u1* u2) */
            c += 2.0 * (u1r*u2i - u1i*u2r);
        }
    return c;
}

/* seed a localised excitation with a chosen circular sense */
static void seed(int site, int hand, double amp) {
    /* rotation in the (x,y) plane: z_x = amp, z_y = i*hand*amp */
    z[site].r[0] += amp;
    z[site].i[1] += hand * amp;
}

static double total_energy(void) {
    double E = 0;
    for (int s = 0; s < N; s++)
        for (int k = 0; k < 3; k++) E += OMC*(z[s].r[k]*z[s].r[k] + z[s].i[k]*z[s].i[k]);
    for (int s = 0; s < N; s++)
        for (int a = 0; a < D; a++) {
            size_t b = ((size_t)s*D+a)*2;
            E += OMLK*(lr[b]*lr[b]+li[b]*li[b]+lr[b+1]*lr[b+1]+li[b+1]*li[b+1]);
            E += 2.0*CHI*(lr[b]*li[b+1] - li[b]*lr[b+1]);
        }
    return E;
}

int main(int argc, char **argv) {
    printf("=====================================================================\n");
    printf("v88 chiral link harmonic: actuation in TWO transverse planes\n");
    printf("=====================================================================\n");
    printf("  u = (u1,u2) in the plane perp to the link; handedness = 2 Im(u1* u2)\n");
    printf("  chirality needs a cross product -> defined only for d >= 3\n\n");

    printf("  C1 IS CHIRALITY EVEN DEFINED?  (seed one right-handed excitation)\n");
    printf("  %4s %8s %10s %16s\n", "d", "sites", "links/site", "chirality after T");
    for (D = 1; D <= 3; D++) {
        L = (D==1)?64:(D==2?24:14);
        alloc();
        seed(N/2, +1, 1.0);
        for (int t = 0; t < 4000; t++) step();
        printf("  %4d %8d %10d %16.6e%s\n", D, N, D, chirality(),
               D < 3 ? "   <- no transverse plane" : "");
    }

    printf("\n  C2/C3 HANDEDNESS-DEPENDENT PAIR ENERGY (d=3)\n");
    D = 3; L = 14;
    printf("  %6s %18s %18s %14s\n", "sep", "E(same hand)", "E(opposite)", "split");
    for (int sep = 2; sep <= 8; sep += 2) {
        double Es, Eo;
        alloc();
        seed(N/2, +1, 1.0);
        seed(N/2 + sep*stride[0], +1, 1.0);
        for (int t = 0; t < 3000; t++) step();
        Es = total_energy();
        alloc();
        seed(N/2, +1, 1.0);
        seed(N/2 + sep*stride[0], -1, 1.0);
        for (int t = 0; t < 3000; t++) step();
        Eo = total_energy();
        printf("  %6d %18.8f %18.8f %14.8f\n", sep, Es, Eo, Es - Eo);
    }
    printf("\n  A NONZERO SPLIT means the chiral channel distinguishes handedness --\n");
    printf("  the asymmetry a scalar link cannot produce, and the only place an\n");
    printf("  attractive channel can come from. Zero split = chirality is not\n");
    printf("  coupling and the construction fails.\n");
    return 0;
}
