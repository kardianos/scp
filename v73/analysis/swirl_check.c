/*  swirl_check.c — Θ-sector swirl diagnostics of a 24-column complex SFA
 *  frame (kernel output or eta_qflow seed). The η curl coupling identifies
 *  the component index of Θ_a with a spatial index (Cosserat), so
 *  θ⃗_u = (theta_x, theta_y, theta_z) is a real vector field on space.
 *  This tool measures whether a ring carries a poloidal Θ circulation —
 *  the "smoke-ring swirl" — a pattern orthogonal to U(1) phase winding.
 *
 *  Global:  Q_phi, Q_theta, N_theta/N, E_int = -η∫(θu·∇×u + θv·∇×v),
 *           helicity H = ∫ θ⃗·(∇×θ⃗) for θu and θv,
 *           field momentum P⃗ = -∫ Σ_f ḟ ∇f  (all 12 pairs).
 *  Loops:   poloidal circulation Γ_pol(φ0,d) = ∮ θ⃗·t̂_χ dl around the tube
 *           cross-section at 8 azimuths × radii d, for θu, θv, and Φu
 *           (control: phase winding gives Φ no coherent vector swirl);
 *           toroidal circulation ∮ θ⃗·φ̂ dl at the peak circle.
 *
 *  Build: gcc -O3 -fopenmp -o swirl_check swirl_check.c -lzstd -lm
 *  Usage: swirl_check file.sfa eta [frame_idx=-1 (last)]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define SFA_IMPLEMENTATION
#include "../../sfa/format/sfa.h"

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s file.sfa eta [frame]\n", argv[0]); return 1; }
    SFA *s = sfa_open(argv[1]);
    if (!s) { fprintf(stderr, "cannot open\n"); return 1; }
    double eta = atof(argv[2]);
    int N = (int)s->Nx; long N3 = (long)N*N*N, NN = (long)N*N;
    double L = s->Lx, dx = 2.0*L/(N-1), dV = dx*dx*dx;
    int fr = (argc > 3) ? atoi(argv[3]) : (int)s->total_frames - 1;
    if (fr < 0) fr = (int)s->total_frames - 1;

    /* column map: POSITION 0-2 = u, 3-5 = v; ANGLE 0-2 = tu, 3-5 = tv;
     * VELOCITY 0-2 = u̇, 3-5 = tu̇, 6-8 = v̇, 9-11 = tv̇ */
    int cu[3],cv[3],ctu[3],ctv[3],cud[3],cvd[3],ctud[3],ctvd[3];
    for (int a=0;a<3;a++) cu[a]=cv[a]=ctu[a]=ctv[a]=cud[a]=cvd[a]=ctud[a]=ctvd[a]=-1;
    for (uint32_t c = 0; c < s->n_columns; c++) {
        int sem = s->columns[c].semantic, cmp = s->columns[c].component;
        if      (sem == SFA_POSITION && cmp < 3) cu[cmp] = c;
        else if (sem == SFA_POSITION && cmp < 6) cv[cmp-3] = c;
        else if (sem == SFA_ANGLE    && cmp < 3) ctu[cmp] = c;
        else if (sem == SFA_ANGLE    && cmp < 6) ctv[cmp-3] = c;
        else if (sem == SFA_VELOCITY && cmp < 3) cud[cmp] = c;
        else if (sem == SFA_VELOCITY && cmp < 6) ctud[cmp-3] = c;
        else if (sem == SFA_VELOCITY && cmp < 9) cvd[cmp-6] = c;
        else if (sem == SFA_VELOCITY && cmp < 12) ctvd[cmp-9] = c;
    }
    void *buf = malloc(s->frame_bytes);
    if (sfa_read_frame(s, fr, buf) != 0) { fprintf(stderr, "read fail\n"); return 1; }
    uint64_t *coff = malloc(s->n_columns * sizeof(uint64_t));
    { uint64_t off = 0; for (uint32_t c = 0; c < s->n_columns; c++) {
        coff[c] = off; off += (uint64_t)N3 * sfa_dtype_size[s->columns[c].dtype]; } }
    #define F(c, ii) ((float*)((uint8_t*)buf + coff[(c)]))[(ii)]
    #define IX(i,j,k) (((long)(i)*N+(j))*N+(k))
    #define DER(c, P, M) ((F((c),(P)) - F((c),(M)))*inv2dx)
    double inv2dx = 1.0/(2.0*dx);

    /* global sums */
    double Qphi=0, Qth=0, Nphi=0, Nth=0, Eint=0, Hu=0, Hv=0;
    double Px=0, Py=0, Pz=0, thmax=0;
    #pragma omp parallel for reduction(+:Qphi,Qth,Nphi,Nth,Eint,Hu,Hv,Px,Py,Pz) reduction(max:thmax)
    for (long ii = 0; ii < N3; ii++) {
        int i=(int)(ii/NN), j=(int)((ii/N)%N), k=(int)(ii%N);
        double th2 = 0;
        for (int a=0;a<3;a++) {
            Qphi += F(cu[a],ii)*F(cvd[a],ii) - F(cv[a],ii)*F(cud[a],ii);
            Qth  += F(ctu[a],ii)*F(ctvd[a],ii) - F(ctv[a],ii)*F(ctud[a],ii);
            Nphi += F(cu[a],ii)*F(cu[a],ii) + F(cv[a],ii)*F(cv[a],ii);
            Nth  += F(ctu[a],ii)*F(ctu[a],ii) + F(ctv[a],ii)*F(ctv[a],ii);
            th2  += F(ctu[a],ii)*F(ctu[a],ii) + F(ctv[a],ii)*F(ctv[a],ii);
        }
        if (th2 > thmax) thmax = th2;
        if (i<1||i>=N-1||j<1||j>=N-1||k<1||k>=N-1) continue;
        long Xp=IX(i+1,j,k),Xm=IX(i-1,j,k),Yp=IX(i,j+1,k),Ym=IX(i,j-1,k),Zp=IX(i,j,k+1),Zm=IX(i,j,k-1);
        long P3[3]={Xp,Yp,Zp}, M3[3]={Xm,Ym,Zm};
        for (int a=0;a<3;a++) {
            int b=(a+1)%3, c=(a+2)%3;
            /* (∇×w)_a = ∂_b w_c − ∂_c w_b */
            double curl_u  = DER(cu[c],  P3[b],M3[b]) - DER(cu[b],  P3[c],M3[c]);
            double curl_v  = DER(cv[c],  P3[b],M3[b]) - DER(cv[b],  P3[c],M3[c]);
            double curl_tu = DER(ctu[c], P3[b],M3[b]) - DER(ctu[b], P3[c],M3[c]);
            double curl_tv = DER(ctv[c], P3[b],M3[b]) - DER(ctv[b], P3[c],M3[c]);
            Eint += -eta*(F(ctu[a],ii)*curl_u + F(ctv[a],ii)*curl_v);
            Hu   += F(ctu[a],ii)*curl_tu;
            Hv   += F(ctv[a],ii)*curl_tv;
        }
        int call[12]; memcpy(call,   cu, 3*sizeof(int)); memcpy(call+3, cv, 3*sizeof(int));
        memcpy(call+6, ctu,3*sizeof(int)); memcpy(call+9, ctv,3*sizeof(int));
        int dall[12]; memcpy(dall,   cud,3*sizeof(int)); memcpy(dall+3, cvd,3*sizeof(int));
        memcpy(dall+6, ctud,3*sizeof(int)); memcpy(dall+9, ctvd,3*sizeof(int));
        for (int f=0; f<12; f++) {
            double fd = F(dall[f],ii);
            Px += -fd*DER(call[f],Xp,Xm);
            Py += -fd*DER(call[f],Yp,Ym);
            Pz += -fd*DER(call[f],Zp,Zm);
        }
    }
    Qphi*=dV; Qth*=dV; Nphi*=dV; Nth*=dV; Eint*=dV; Hu*=dV; Hv*=dV;
    Px*=dV; Py*=dV; Pz*=dV;

    /* tube geometry: peak circle radius in z=0, tube center in the
     * y=0 half-plane (amp^2-weighted), tube offset of the peak */
    int kz = N/2, jy = N/2;
    double rho_pk = 0, amp_pk = 0;
    for (int i = 1; i < N-1; i++) for (int j = 1; j < N-1; j++) {
        double x = -L+i*dx, y = -L+j*dx; long ii = IX(i,j,kz);
        double a2 = 0; for (int a=0;a<3;a++) a2 += F(cu[a],ii)*F(cu[a],ii)+F(cv[a],ii)*F(cv[a],ii);
        if (a2 > amp_pk) { amp_pk = a2; rho_pk = sqrt(x*x+y*y); }
    }
    double Rc_n = 0, Rc_d = 0;
    for (int i = N/2+1; i < N-1; i++) for (int k = 1; k < N-1; k++) {
        double x = -L+i*dx; long ii = IX(i,jy,k);
        double a2 = 0; for (int a=0;a<3;a++) a2 += F(cu[a],ii)*F(cu[a],ii)+F(cv[a],ii)*F(cv[a],ii);
        Rc_n += a2*x; Rc_d += a2;
    }
    double Rc = (Rc_d>0)? Rc_n/Rc_d : rho_pk;

    printf("frame %d (t=%.2f): Q_phi=%.4f  Q_theta=%.4f  N_theta/N=%.4f  max|Th|=%.4e\n",
           fr, sfa_frame_time(s, fr), Qphi, Qth, Nth/(Nphi+Nth), sqrt(thmax));
    printf("  E_int=%.6f  H_u=%.6f  H_v=%.6f  P=(%.3e, %.3e, %.3e)\n", Eint, Hu, Hv, Px, Py, Pz);
    printf("  tube: R_c=%.3f  rho_peak=%.3f\n", Rc, rho_pk);

    /* poloidal circulation Γ(φ0,d) = ∮ w⃗·t̂ dl, t̂ = −sinχ ρ̂(φ0) + cosχ ẑ,
     * around the circle of radius d centred on (R_c ρ̂(φ0), z=0). */
    int M = 256;
    double dlist[2] = {0.75, 1.5};
    for (int di = 0; di < 2; di++) {
        double d = dlist[di];
        double mu_tu=0, mu_tv=0, mu_u=0, mn_tu=1e30, mx_tu=-1e30;
        printf("  poloidal d=%.2f  Gamma[theta_u, theta_v, phi_u(ctrl)] per azimuth:\n", d);
        for (int p = 0; p < 8; p++) {
            double ph0 = 2.0*M_PI*p/8.0, cph = cos(ph0), sph = sin(ph0);
            double g_tu=0, g_tv=0, g_u=0;
            for (int t = 0; t < M; t++) {
                double ch = 2.0*M_PI*t/M, dl = 2.0*M_PI*d/M;
                double rho = Rc + d*cos(ch), z = d*sin(ch);
                double x = rho*cph, y = rho*sph;
                int i=(int)lround((x+L)/dx), j=(int)lround((y+L)/dx), k=(int)lround((z+L)/dx);
                if (i<1||i>=N-1||j<1||j>=N-1||k<1||k>=N-1) continue;
                long ii = IX(i,j,k);
                /* t̂ = −sinχ (cph, sph, 0) + cosχ (0,0,1) */
                double tx = -sin(ch)*cph, ty = -sin(ch)*sph, tz = cos(ch);
                g_tu += (F(ctu[0],ii)*tx + F(ctu[1],ii)*ty + F(ctu[2],ii)*tz)*dl;
                g_tv += (F(ctv[0],ii)*tx + F(ctv[1],ii)*ty + F(ctv[2],ii)*tz)*dl;
                g_u  += (F(cu[0],ii)*tx + F(cu[1],ii)*ty + F(cu[2],ii)*tz)*dl;
            }
            printf("    phi0=%3.0f  %+.5f  %+.5f  %+.5f\n", ph0*180/M_PI, g_tu, g_tv, g_u);
            mu_tu += g_tu/8; mu_tv += g_tv/8; mu_u += g_u/8;
            if (g_tu<mn_tu) mn_tu=g_tu; if (g_tu>mx_tu) mx_tu=g_tu;
        }
        printf("    mean: Gamma_tu=%+.5f [%+.5f..%+.5f]  Gamma_tv=%+.5f  Gamma_u_ctrl=%+.5f\n",
               mu_tu, mn_tu, mx_tu, mu_tv, mu_u);
    }

    /* toroidal Θ circulation at the peak circle */
    double gt_tu=0, gt_tv=0;
    for (int t = 0; t < M; t++) {
        double ph = 2.0*M_PI*t/M, dl = 2.0*M_PI*rho_pk/M;
        double x = rho_pk*cos(ph), y = rho_pk*sin(ph);
        int i=(int)lround((x+L)/dx), j=(int)lround((y+L)/dx);
        if (i<1||i>=N-1||j<1||j>=N-1) continue;
        long ii = IX(i,j,kz);
        double tx = -sin(ph), ty = cos(ph);
        gt_tu += (F(ctu[0],ii)*tx + F(ctu[1],ii)*ty)*dl;
        gt_tv += (F(ctv[0],ii)*tx + F(ctv[1],ii)*ty)*dl;
    }
    printf("  toroidal (rho_peak): Gamma_tu=%+.5f  Gamma_tv=%+.5f\n", gt_tu, gt_tv);

    sfa_close(s);
    return 0;
}
