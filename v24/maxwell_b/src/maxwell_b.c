/*
 * maxwell_b.c — Three complex scalars with U(1) gauge field (quark charges)
 *
 * Fields: Psi_a = phi_a + i*chi_a (a=1,2,3), gauge field A
 * Total: 7 real field arrays
 *
 * Charges: e_a = e_scale * {+2/3, +2/3, -1/3}  (UUD = proton)
 *   or     e_a = e_scale * {+2/3, -1/3, -1/3}  (UDD = neutron)
 *
 * Lagrangian:
 *   L = sum_a |D_mu Psi_a|^2 - m_a^2|Psi_a|^2 - V(|P|^2) - (1/4)F^2
 *   D_mu = d_mu - i*e_a*A_mu,  P = Psi_1*Psi_2*Psi_3
 *   V(|P|^2) = (mu/2)|P|^2/(1 + kappa*|P|^2)
 *
 * 1D temporal gauge (A_0=0). Velocity Verlet integrator.
 *
 * Compile: gcc -O3 -Wall -o maxwell_b src/maxwell_b.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Parameters */
static double mu       = -20.0;
static double kap      = 20.0;
static double mass_arr[3] = {1.0, 1.0, 1.0};  /* m1, m2, m3 */
static double e_scale  = 0.0;
static double charge_pattern[3] = {2.0/3, 2.0/3, -1.0/3}; /* UUD default */
static double A_init_amp = 0.8;
static double sigma_init = 3.0;
static double chi_seed = 1e-3;  /* seed amplitude for imaginary parts */
static int    Nx       = 4000;
static double xmax     = 100.0;
static double tfinal   = 10000.0;
static int    mode     = 0;  /* 0=e_scale scan, 1=UUD, 2=UDD, 3=single run */
static char   outdir[512] = "v24/maxwell_b/data";

static void parse_args(int argc, char **argv)
{
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mu"))      mu = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))   kap = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-m1"))      mass_arr[0] = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-m2"))      mass_arr[1] = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-m3"))      mass_arr[2] = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-e"))       e_scale = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))       A_init_amp = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sigma"))   sigma_init = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-seed"))    chi_seed = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-Nx"))      Nx = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-xmax"))    xmax = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-tfinal"))  tfinal = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mode"))    mode = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))       strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else { fprintf(stderr, "Unknown: %s\n", argv[i]); exit(1); }
    }
}

/*
 * V = (mu/2)|P|^2/(1+kappa|P|^2)
 * |P|^2 = |Psi_1|^2 |Psi_2|^2 |Psi_3|^2
 *
 * dV/dPsi_a* = (mu/2) * P * cofactor_a / (1+kappa|P|^2)^2
 * where cofactor_1 = Psi_2*Psi_3  (complex conjugate of product of others times P)
 *
 * More carefully: |P|^2 = P * P*,  P = Psi_1 Psi_2 Psi_3
 * d|P|^2/dPsi_a* = P * (product of Psi_b for b!=a)
 * So dV/dPsi_a* = (mu/2) * P * cof_a / (1+kappa|P|^2)^2
 *   where cof_1 = Psi_2*Psi_3*, etc.  (WRONG — need conjugate)
 *
 * Actually: |P|^2 = Psi_1 Psi_2 Psi_3 Psi_1* Psi_2* Psi_3*
 * d|P|^2/dPsi_1* = Psi_1 Psi_2 Psi_3 * Psi_2* Psi_3* = P * |Psi_2|^2|Psi_3|^2 / Psi_1*
 * Hmm, simpler: d|P|^2/dPsi_1* = P * Psi_2* Psi_3*  ... wait.
 * Let me be explicit. P = Psi_1 Psi_2 Psi_3, P* = Psi_1* Psi_2* Psi_3*.
 * |P|^2 = P Pbar.  d(P Pbar)/dPsi_1bar = P * dPbar/dPsi_1bar = P * Psi_2bar Psi_3bar.
 *
 * So: dV/dPsi_1* = (mu/2) / (1+kappa|P|^2)^2 * P * Psi_2* Psi_3*
 *
 * In real/imag decomposition, Psi_a = phi_a + i chi_a:
 *   force_phi_a = -Re(dV/dPsi_a*)
 *   force_chi_a = -Im(dV/dPsi_a*)
 */

typedef struct {
    double re, im;
} cplx;

static cplx cmul(cplx a, cplx b)
{
    return (cplx){a.re*b.re - a.im*b.im, a.re*b.im + a.im*b.re};
}

static cplx cconj(cplx a)
{
    return (cplx){a.re, -a.im};
}

/* Compute -dV/dPsi_a* (returns complex), given Psi_1, Psi_2, Psi_3 */
static cplx force_pot(cplx psi1, cplx psi2, cplx psi3, int a)
{
    cplx P = cmul(cmul(psi1, psi2), psi3);
    double P2 = P.re * P.re + P.im * P.im;
    double denom2 = (1.0 + kap * P2) * (1.0 + kap * P2);
    double coeff = -mu / (2.0 * denom2);  /* minus sign: force = -dV/d... , and mu<0 */

    /* cofactor = product of conjugates of the OTHER two fields */
    cplx cof;
    switch (a) {
        case 0: cof = cmul(cconj(psi2), cconj(psi3)); break;
        case 1: cof = cmul(cconj(psi1), cconj(psi3)); break;
        case 2: cof = cmul(cconj(psi1), cconj(psi2)); break;
        default: cof = (cplx){0,0};
    }

    /* dV/dPsi_a* = (mu/2) P * cof / (1+kap*P2)^2 */
    /* force = -dV/dPsi_a* */
    cplx Pcof = cmul(P, cof);
    return (cplx){coeff * Pcof.re, coeff * Pcof.im};
}

typedef struct {
    /* field values: phi[a], chi[a] for a=0,1,2, and gauge A */
    double *phi[3], *chi[3], *A;
    /* velocities */
    double *vphi[3], *vchi[3], *vA;
    /* accelerations */
    double *aphi[3], *achi[3], *aA;
} fields_t;

static fields_t alloc_fields(int N)
{
    fields_t f;
    for (int a = 0; a < 3; a++) {
        f.phi[a] = calloc(N, sizeof(double));
        f.chi[a] = calloc(N, sizeof(double));
        f.vphi[a] = calloc(N, sizeof(double));
        f.vchi[a] = calloc(N, sizeof(double));
        f.aphi[a] = calloc(N, sizeof(double));
        f.achi[a] = calloc(N, sizeof(double));
    }
    f.A = calloc(N, sizeof(double));
    f.vA = calloc(N, sizeof(double));
    f.aA = calloc(N, sizeof(double));
    return f;
}

static void free_fields(fields_t *f)
{
    for (int a = 0; a < 3; a++) {
        free(f->phi[a]); free(f->chi[a]);
        free(f->vphi[a]); free(f->vchi[a]);
        free(f->aphi[a]); free(f->achi[a]);
    }
    free(f->A); free(f->vA); free(f->aA);
}

static double *damp_arr;
static double dx_g, dx2_g;
static double charges[3];

/* Compute all accelerations */
static void compute_acc(fields_t *f)
{
    double dx = dx_g, dx2 = dx2_g;

    for (int a = 0; a < 3; a++) {
        f->aphi[a][0] = f->aphi[a][Nx-1] = 0;
        f->achi[a][0] = f->achi[a][Nx-1] = 0;
    }
    f->aA[0] = f->aA[Nx-1] = 0;

    for (int i = 1; i < Nx - 1; i++) {
        cplx psi[3];
        for (int a = 0; a < 3; a++)
            psi[a] = (cplx){f->phi[a][i], f->chi[a][i]};

        double Ai = f->A[i];
        double dAdx = (f->A[i+1] - f->A[i-1]) / (2.0 * dx);

        /* gauge current for Maxwell eq */
        double j_total = 0;

        for (int a = 0; a < 3; a++) {
            double ea = charges[a];
            double m2a = mass_arr[a] * mass_arr[a];

            double lapl_phi = (f->phi[a][i+1] - 2.0*f->phi[a][i] + f->phi[a][i-1]) / dx2;
            double lapl_chi = (f->chi[a][i+1] - 2.0*f->chi[a][i] + f->chi[a][i-1]) / dx2;

            double dchi_dx = (f->chi[a][i+1] - f->chi[a][i-1]) / (2.0 * dx);
            double dphi_dx = (f->phi[a][i+1] - f->phi[a][i-1]) / (2.0 * dx);

            /* Coupling potential force */
            cplx fp = force_pot(psi[0], psi[1], psi[2], a);

            /* EOM for phi_a:
             * d^2 phi/dt^2 = lapl_phi + 2*ea*A*dchi/dx + ea*(dA/dx)*chi
             *                - ea^2*A^2*phi - m^2*phi + Re(force_pot)
             * Wait — sign of covariant terms: |D Psi|^2 = |dPsi - i e A Psi|^2
             * Expanding spatial part:
             * |d_x Psi - i e A Psi|^2 = (d_x phi + e A chi)^2 + (d_x chi - e A phi)^2
             * EL for phi: d^2phi/dt^2 = d^2phi/dx^2 + e*(dA/dx)*chi + 2eA*dchi/dx
             *             + e^2 A^2 phi - m^2 phi - Re(dV/dPsi*)
             * Wait, careful with sign. Let me redo from L = |D_x Psi|^2 = (d_x phi + eA chi)^2 + (d_x chi - eA phi)^2
             * dL/dphi = 0 (no direct phi), d(dL/d(d_x phi))/dx = d/dx[2(d_x phi + eA chi)] = 2(d_xx phi + e dA/dx chi + eA d_x chi)
             * dL/dphi from second term: 2(d_x chi - eA phi)(-eA) = -2eA(d_x chi - eA phi)
             * Also from mass: -m^2 phi contribution.
             * EL: d^2phi/dt^2 = d^2phi/dx^2 + e(dA/dx)chi + eA(d_x chi) + eA(d_x chi) - e^2 A^2 phi - m^2 phi + force_pot_re
             * Hmm wait, let me just be careful. The spatial kinetic is (with D_x = d_x - i e A):
             * |D_x Psi|^2 = |(d_x phi + eA chi) + i(d_x chi - eA phi)|^2
             *             = (d_x phi + eA chi)^2 + (d_x chi - eA phi)^2
             * Vary w.r.t. phi (treating d_x phi as independent via integration by parts):
             * d/dx[2(d_x phi + eA chi)] - 2eA(-(d_x chi - eA phi)) = 0 from spatial kinetic
             * => 2(d_xx phi + e(dA/dx)chi + eA d_x chi) + 2eA(d_x chi - eA phi)
             * = 2 d_xx phi + 2e(dA/dx)chi + 4eA d_x chi - 2e^2A^2 phi
             * Hmm, that gives a factor 4. Let me recount.
             * S_spatial = integral |D_x Psi|^2 dx. L has +|D_x Psi|^2 kinetic (positive, since -(-|D_x|^2)).
             * Actually in the Lagrangian L = |D_t Psi|^2 - |D_x Psi|^2 - m^2|Psi|^2 - V
             * In temporal gauge D_t = d_t, so |D_t Psi|^2 = (d_t phi)^2 + (d_t chi)^2.
             * EL for phi: d^2phi/dt^2 = -dL_spatial/dphi + d/dx(dL_spatial/d(d_x phi))
             * where L_spatial = -|D_x Psi|^2 - m^2|Psi|^2 - V
             * -|D_x Psi|^2 term:
             * d[-|D_x|^2]/dphi = -d/dphi[(d_x chi - eA phi)^2] = 2eA(d_x chi - eA phi)
             * d[-|D_x|^2]/d(d_x phi) = -2(d_x phi + eA chi)
             * => d/dx[-2(d_x phi + eA chi)] + 2eA(d_x chi - eA phi)
             * = -2 d_xx phi - 2e(dA/dx)chi - 2eA(d_x chi) + 2eA d_x chi - 2e^2A^2 phi
             * = -2 d_xx phi - 2e(dA/dx)chi - 2e^2A^2 phi
             * From mass: d[-m^2|Psi|^2]/dphi = -2m^2 phi
             * But wait, the (1/2) convention. Usually L = (1/2)|D_t|^2 - (1/2)|D_x|^2 - (m^2/2)|Psi|^2
             * With the 1/2: everything divided by 2.
             * d^2phi/dt^2 = d^2phi/dx^2 + e(dA/dx)chi + e^2A^2 phi - m^2 phi + ...
             * Wait that gives +e^2A^2 phi? No:
             * From -(1/2)|D_x|^2, the EL gives:
             *   d/dx[(d_x phi + eAchi)] + eA(d_x chi - eAphi)      (flip sign for EL)
             * Hmm I keep getting confused. Let me just use the standard result.
             *
             * Standard: for L = (1/2)|D_mu Psi|^2,  D_mu = d_mu - i e A_mu
             * EOM: D_mu D^mu Psi + m^2 Psi + dV/dPsi* = 0
             * In 1+1 temporal gauge: d_t^2 Psi - (d_x - i e A)^2 Psi + m^2 Psi + dV/dPsi* = 0
             * (d_x - i e A)^2 Psi = d_xx Psi - 2ieA d_x Psi - ie(dA/dx)Psi - e^2A^2 Psi
             * So: d_t^2 Psi = d_xx Psi - 2ieA d_x Psi - ie(dA/dx)Psi - e^2A^2 Psi - m^2 Psi - dV/dPsi*
             *
             * Wait, -(d_x - ieA)^2 has wrong sign. Let me be more careful.
             * d_t^2 Psi = (d_x - ieA)^2 Psi - m^2 Psi - dV/dPsi*
             * (d_x - ieA)^2 = d_xx - 2ieA d_x - ie(dA/dx) - e^2A^2
             *
             * Re part: d_t^2 phi = d_xx phi + 2eA d_x chi + e(dA/dx)chi - e^2A^2 phi - m^2 phi - Re(dV/dPsi*)
             * Im part: d_t^2 chi = d_xx chi - 2eA d_x phi - e(dA/dx)phi - e^2A^2 chi - m^2 chi - Im(dV/dPsi*)
             */
            f->aphi[a][i] = lapl_phi
                + 2.0*ea*Ai*dchi_dx + ea*dAdx*f->chi[a][i]
                - ea*ea*Ai*Ai*f->phi[a][i]
                - m2a*f->phi[a][i]
                + fp.re;   /* force_pot already returns -dV/dPsi*, so + here */

            f->achi[a][i] = lapl_chi
                - 2.0*ea*Ai*dphi_dx - ea*dAdx*f->phi[a][i]
                - ea*ea*Ai*Ai*f->chi[a][i]
                - m2a*f->chi[a][i]
                + fp.im;

            /* Current: j_a = e_a * Im(Psi_a* D_x Psi_a)
             * = e_a * [phi_a (d_x chi_a - e_a A phi_a) - chi_a (d_x phi_a + e_a A chi_a)]
             * Hmm, more directly:
             * j_a = e_a*(phi_a d_x chi_a - chi_a d_x phi_a) - e_a^2 A |Psi_a|^2
             * Wait, from standard: j = e Im(Psi* D_x Psi) = e Im(Psi*(d_x Psi - ieA Psi))
             * = e Im(Psi* d_x Psi) + e^2 A |Psi|^2
             * = e(phi d_x chi - chi d_x phi) + e^2 A (phi^2 + chi^2)
             *
             * Maxwell: d_t^2 A - d_x^2 A = -j_total = -sum_a j_a
             */
            double mod2 = f->phi[a][i]*f->phi[a][i] + f->chi[a][i]*f->chi[a][i];
            j_total += ea*(f->phi[a][i]*dchi_dx - f->chi[a][i]*dphi_dx)
                     + ea*ea*Ai*mod2;
        }

        /* Maxwell equation: d_t^2 A = d_x^2 A - j_total */
        double lapl_A = (f->A[i+1] - 2.0*f->A[i] + f->A[i-1]) / dx2;
        f->aA[i] = lapl_A - j_total;
    }
}

typedef struct {
    double E_kin, E_grad, E_mass, E_pot, E_em, E_total;
    double Q_em;       /* total EM charge (Gauss law integral) */
    double peak[3];    /* peak |Psi_a| */
    double phi0[3];    /* phi_a at center */
    double chi0[3];    /* chi_a at center */
    double A0;         /* A at center */
    double f_core;     /* fraction of energy in core */
    double j_max;      /* peak current density */
} diag_t;

static diag_t diagnose(fields_t *f)
{
    double dx = dx_g;
    diag_t d = {0};
    int ic = Nx / 2;
    double core_r = 3.0 * sigma_init;
    double E_all = 0, E_core = 0;

    for (int i = 1; i < Nx - 1; i++) {
        double x = -xmax + i * dx;

        cplx psi[3];
        for (int a = 0; a < 3; a++)
            psi[a] = (cplx){f->phi[a][i], f->chi[a][i]};

        /* Kinetic energy */
        double ek = 0;
        for (int a = 0; a < 3; a++) {
            ek += 0.5*(f->vphi[a][i]*f->vphi[a][i] + f->vchi[a][i]*f->vchi[a][i]);
        }
        ek += 0.5*f->vA[i]*f->vA[i];  /* EM kinetic = (1/2)(d_t A)^2 = (1/2)E^2 */

        /* Gradient energy — use covariant derivative for scalars */
        double eg = 0;
        for (int a = 0; a < 3; a++) {
            double dphi = (f->phi[a][i+1] - f->phi[a][i-1]) / (2.0*dx);
            double dchi = (f->chi[a][i+1] - f->chi[a][i-1]) / (2.0*dx);
            double ea = charges[a];
            /* |D_x Psi|^2 = (dphi + eA chi)^2 + (dchi - eA phi)^2 */
            double Dx_re = dphi + ea*f->A[i]*f->chi[a][i];
            double Dx_im = dchi - ea*f->A[i]*f->phi[a][i];
            eg += 0.5*(Dx_re*Dx_re + Dx_im*Dx_im);
        }
        /* Magnetic energy: (1/2)(d_x A)^2 — but in 1D B=d_x A */
        double dA = (f->A[i+1] - f->A[i-1]) / (2.0*dx);
        eg += 0.5*dA*dA;

        /* Mass energy */
        double em = 0;
        for (int a = 0; a < 3; a++) {
            double mod2 = f->phi[a][i]*f->phi[a][i] + f->chi[a][i]*f->chi[a][i];
            em += 0.5*mass_arr[a]*mass_arr[a]*mod2;
        }

        /* Potential */
        cplx P = cmul(cmul(psi[0], psi[1]), psi[2]);
        double P2 = P.re*P.re + P.im*P.im;
        double V = 0.5*mu*P2/(1.0 + kap*P2);

        /* Peaks */
        for (int a = 0; a < 3; a++) {
            double mod = sqrt(f->phi[a][i]*f->phi[a][i] + f->chi[a][i]*f->chi[a][i]);
            if (mod > d.peak[a]) d.peak[a] = mod;
        }

        d.E_kin += ek * dx;
        d.E_grad += eg * dx;
        d.E_mass += em * dx;
        d.E_pot += V * dx;

        double etot = ek + eg + em + V;
        E_all += etot * dx;
        if (fabs(x) < core_r) E_core += etot * dx;

        /* Current for charge diagnostic */
        for (int a = 0; a < 3; a++) {
            double dphi = (f->phi[a][i+1] - f->phi[a][i-1]) / (2.0*dx);
            double dchi = (f->chi[a][i+1] - f->chi[a][i-1]) / (2.0*dx);
            double ea = charges[a];
            double j_a = ea*(f->phi[a][i]*dchi - f->chi[a][i]*dphi)
                       + ea*ea*f->A[i]*(f->phi[a][i]*f->phi[a][i]+f->chi[a][i]*f->chi[a][i]);
            if (fabs(j_a) > d.j_max) d.j_max = fabs(j_a);
        }
    }

    /* EM energy = kinetic(E^2/2) part of E_kin + magnetic(B^2/2) part of E_grad
       Already included above. Separate it: */
    d.E_em = 0;
    for (int i = 1; i < Nx - 1; i++) {
        double dA_l = (f->A[i+1] - f->A[i-1]) / (2.0*dx);
        d.E_em += (0.5*f->vA[i]*f->vA[i] + 0.5*dA_l*dA_l) * dx;
    }

    d.E_total = d.E_kin + d.E_grad + d.E_mass + d.E_pot;
    d.f_core = (E_all > 1e-20) ? E_core / E_all : 0.0;

    for (int a = 0; a < 3; a++) {
        d.phi0[a] = f->phi[a][ic];
        d.chi0[a] = f->chi[a][ic];
    }
    d.A0 = f->A[ic];

    /* Gauss law charge: Q = integral of charge density rho = -d_x E = -d_x(d_t A)
     * Equivalently Q = E(x=-inf) - E(x=+inf). In practice: E = d_t A at boundaries.
     * Or just integrate rho = sum_a j^0_a. In temporal gauge, j^0 isn't directly available,
     * but from Gauss: Q = -[d_x(d_t A)] integrated = d_t A at left - d_t A at right */
    d.Q_em = f->vA[1] - f->vA[Nx-2];  /* approximate */

    return d;
}

static void run_single(double es, const double cp[3], const double masses[3],
                       const char *label)
{
    for (int a = 0; a < 3; a++) {
        charges[a] = es * cp[a];
        mass_arr[a] = masses[a];
    }

    double dx = 2.0 * xmax / (Nx - 1);
    double dx2 = dx * dx;
    dx_g = dx; dx2_g = dx2;

    /* CFL: max speed is c=1, plus gauge contributions */
    double kmax = M_PI / dx;
    double m_max = fmax(mass_arr[0], fmax(mass_arr[1], mass_arr[2]));
    double dt = 0.4 * 2.0 / sqrt(kmax*kmax + m_max*m_max);
    int Nt = (int)(tfinal / dt) + 1;

    printf("\n=== %s ===\n", label);
    printf("  charges: (%.4f, %.4f, %.4f), total=%.4f\n",
           charges[0], charges[1], charges[2],
           charges[0]+charges[1]+charges[2]);
    printf("  masses: (%.4f, %.4f, %.4f)\n", mass_arr[0], mass_arr[1], mass_arr[2]);
    printf("  mu=%.1f kappa=%.1f A=%.3f sigma=%.2f seed=%.1e\n",
           mu, kap, A_init_amp, sigma_init, chi_seed);
    printf("  Nx=%d xmax=%.0f dx=%.5f dt=%.6f Nt=%d\n", Nx, xmax, dx, dt, Nt);

    fields_t f = alloc_fields(Nx);

    /* Damping */
    damp_arr = malloc(Nx * sizeof(double));
    double x_abs = xmax * 0.75;
    for (int i = 0; i < Nx; i++) {
        double x = -xmax + i * dx;
        double ax = fabs(x);
        if (ax > x_abs) {
            double frac = (ax - x_abs) / (xmax - x_abs);
            damp_arr[i] = 1.0 - 0.98 * frac * frac;
        } else {
            damp_arr[i] = 1.0;
        }
    }

    /* Initialize: real Gaussians + small chi seed */
    int ic = Nx / 2;
    for (int a = 0; a < 3; a++)
        for (int i = 0; i < Nx; i++) {
            double x = -xmax + i * dx;
            f.phi[a][i] = A_init_amp * exp(-x*x / (2.0*sigma_init*sigma_init));
            /* Seed chi with small odd perturbation (so it's localized and has structure) */
            if (es > 0)
                f.chi[a][i] = chi_seed * x * exp(-x*x / (2.0*sigma_init*sigma_init));
        }
    /* A starts at zero */

    compute_acc(&f);

    /* Output file */
    char tspath[600];
    snprintf(tspath, sizeof(tspath), "%s/%s_ts.tsv", outdir, label);
    FILE *fts = fopen(tspath, "w");
    if (!fts) { fprintf(stderr, "Cannot open %s\n", tspath); return; }
    fprintf(fts, "time\tphi1_0\tphi2_0\tphi3_0\tchi1_0\tchi2_0\tchi3_0\tA_0\t"
                 "peak1\tpeak2\tpeak3\t"
                 "E_kin\tE_grad\tE_mass\tE_pot\tE_em\tE_total\tQ_em\tf_core\tj_max\n");

    /* DFT storage */
    int max_dft = 50000;
    double *phi0_hist = malloc(max_dft * sizeof(double));
    double *t_hist = malloc(max_dft * sizeof(double));
    int n_dft = 0;

    int rec_every = Nt / 20000;
    if (rec_every < 1) rec_every = 1;
    int print_every = Nt / 20;
    if (print_every < 1) print_every = 1;
    int dft_every = Nt / max_dft;
    if (dft_every < 1) dft_every = 1;

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        if (n % dft_every == 0 && n_dft < max_dft) {
            phi0_hist[n_dft] = f.phi[0][ic];
            t_hist[n_dft] = t;
            n_dft++;
        }

        int do_rec = (n % rec_every == 0);
        int do_print = (n % print_every == 0);

        if (do_rec || do_print) {
            diag_t dg = diagnose(&f);

            if (do_rec)
                fprintf(fts, "%.4f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t"
                             "%.6e\t%.6e\t%.6e\t"
                             "%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.4f\t%.6e\n",
                        t, dg.phi0[0], dg.phi0[1], dg.phi0[2],
                        dg.chi0[0], dg.chi0[1], dg.chi0[2], dg.A0,
                        dg.peak[0], dg.peak[1], dg.peak[2],
                        dg.E_kin, dg.E_grad, dg.E_mass, dg.E_pot, dg.E_em,
                        dg.E_total, dg.Q_em, dg.f_core, dg.j_max);
            if (do_print)
                printf("  t=%7.0f  phi0=(%.3f,%.3f,%.3f)  chi0=(%.2e,%.2e,%.2e)  "
                       "A0=%.2e  E=%.4f  E_em=%.2e  Q=%.2e  fc=%.3f\n",
                       t, dg.phi0[0], dg.phi0[1], dg.phi0[2],
                       dg.chi0[0], dg.chi0[1], dg.chi0[2],
                       dg.A0, dg.E_total, dg.E_em, dg.Q_em, dg.f_core);
        }

        if (n == Nt) break;

        /* Velocity Verlet */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++) {
                f.vphi[a][i] += 0.5 * dt * f.aphi[a][i];
                f.vchi[a][i] += 0.5 * dt * f.achi[a][i];
            }
        for (int i = 1; i < Nx - 1; i++)
            f.vA[i] += 0.5 * dt * f.aA[i];

        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++) {
                f.phi[a][i] += dt * f.vphi[a][i];
                f.chi[a][i] += dt * f.vchi[a][i];
            }
        for (int i = 1; i < Nx - 1; i++)
            f.A[i] += dt * f.vA[i];

        compute_acc(&f);

        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++) {
                f.vphi[a][i] += 0.5 * dt * f.aphi[a][i];
                f.vchi[a][i] += 0.5 * dt * f.achi[a][i];
            }
        for (int i = 1; i < Nx - 1; i++)
            f.vA[i] += 0.5 * dt * f.aA[i];

        /* Absorbing boundary */
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++) {
                f.vphi[a][i] *= damp_arr[i];
                f.phi[a][i]  *= damp_arr[i];
                f.vchi[a][i] *= damp_arr[i];
                f.chi[a][i]  *= damp_arr[i];
            }
        for (int i = 0; i < Nx; i++) {
            f.vA[i] *= damp_arr[i];
            f.A[i]  *= damp_arr[i];
        }
    }

    fclose(fts);

    /* DFT of phi_1(0,t) — second half */
    int dft_start = n_dft / 2;
    double peak_om = 0, peak_pow = 0;
    if (n_dft - dft_start > 100) {
        char dftpath[600];
        snprintf(dftpath, sizeof(dftpath), "%s/%s_spectrum.tsv", outdir, label);
        FILE *fdft = fopen(dftpath, "w");
        fprintf(fdft, "omega\tpower\n");
        double T = t_hist[n_dft-1] - t_hist[dft_start];
        int nf = 500;
        for (int k = 0; k < nf; k++) {
            double omega = 3.0 * m_max * k / nf;
            double re = 0, im = 0;
            for (int j = dft_start; j < n_dft; j++) {
                double dtj = (j > dft_start) ?
                    (t_hist[j]-t_hist[j-1]) : (t_hist[dft_start+1]-t_hist[dft_start]);
                re += phi0_hist[j] * cos(omega * t_hist[j]) * dtj;
                im += phi0_hist[j] * sin(omega * t_hist[j]) * dtj;
            }
            double pw = (re*re + im*im) / (T*T);
            fprintf(fdft, "%.6f\t%.6e\n", omega, pw);
            if (pw > peak_pow) { peak_pow = pw; peak_om = omega; }
        }
        fclose(fdft);
    }
    printf("  Spectrum peak: omega=%.4f (mass=%.4f), oscillon=%s\n",
           peak_om, m_max, (peak_om > 0.01 && peak_om < m_max) ? "YES" : "NO");

    /* Final diagnostics */
    diag_t dg_final = diagnose(&f);
    printf("  FINAL: E_total=%.4f  E_em=%.6f  Q=%.6f  fc=%.3f\n",
           dg_final.E_total, dg_final.E_em, dg_final.Q_em, dg_final.f_core);
    printf("  Output: %s\n", tspath);

    free_fields(&f);
    free(damp_arr);
    free(phi0_hist);
    free(t_hist);
}

int main(int argc, char **argv)
{
    parse_args(argc, argv);

    if (mode == 0) {
        /* e_scale scan with UUD charges */
        double cp[3] = {2.0/3, 2.0/3, -1.0/3};
        double masses[3] = {1.0, 1.0, 1.0};
        double e_vals[] = {0.0, 0.1, 0.3, 1.0};
        int ne = 4;
        for (int ie = 0; ie < ne; ie++) {
            char label[64];
            snprintf(label, sizeof(label), "UUD_e%.2f", e_vals[ie]);
            run_single(e_vals[ie], cp, masses, label);
        }
    } else if (mode == 1) {
        /* UUD proton: m_U=1.0, m_D=0.95, Q=+1 */
        double cp[3] = {2.0/3, 2.0/3, -1.0/3};
        double masses[3] = {1.0, 1.0, 0.95};
        char label[64];
        snprintf(label, sizeof(label), "UUD_proton_e%.2f", e_scale);
        run_single(e_scale, cp, masses, label);
    } else if (mode == 2) {
        /* UDD neutron: m_U=1.0, m_D=0.95, Q=0 */
        double cp[3] = {2.0/3, -1.0/3, -1.0/3};
        double masses[3] = {1.0, 0.95, 0.95};
        char label[64];
        snprintf(label, sizeof(label), "UDD_neutron_e%.2f", e_scale);
        run_single(e_scale, cp, masses, label);
    } else if (mode == 3) {
        /* Single run with current e_scale */
        run_single(e_scale, charge_pattern, mass_arr, "single");
    } else {
        fprintf(stderr, "Unknown mode %d\n", mode);
        return 1;
    }

    printf("\nDone.\n");
    return 0;
}
