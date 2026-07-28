/*  sfa_boost.c — v86 EX-1: boost an EXISTING SFA field configuration along x.
 *
 *  Every seed generator in the repo builds a boosted object from a radial
 *  PROFILE (gen_qball_boost, gen_qball_pair_boost, ...). EX-1 needs something
 *  none of them do: take a relaxed end-state that was produced by a simulation
 *  -- the X10c ball+cloud at t = 2000 -- and set it moving without disturbing
 *  the structure that took 2000 t.u. to form.
 *
 *  THE TRANSFORMATION, and its stated order of accuracy.
 *  A boost mixes space and time, so strictly it needs the solution at other
 *  times: at lab time t' = 0 the rest-frame time is t = gamma v x', which
 *  varies across the box. Expanding to first order in v with the frame's own
 *  time derivative (which the SFA stores):
 *
 *      Phi'(x')    = R(x') Phi(gamma x')
 *      Phidot'(x') = R(x') gamma [ Phidot(gamma x') - gamma v d_x Phi(gamma x') ]
 *
 *  where R(x') is rotation by the phase -gamma omega_loc v x'. The rotation is
 *  applied EXACTLY, not Taylor-expanded: the rest-time offset is gamma v x',
 *  so the phase it induces is gamma omega v x', which at v = 0.02 and |x| = 40
 *  is already ~1.1 radians. A first-order expansion of that rotation is not
 *  small, and in testing it lost exactly half the momentum.
 *
 *  omega_loc = (u vdot - v udot)/(u^2+v^2) is measured PER VOXEL AND PER
 *  COMPONENT, so a state carrying several clocks transforms correctly -- which
 *  matters here, because the X10c ball and its opposite-charge cloud rotate at
 *  different frequencies of opposite sign. A single global omega would boost
 *  one correctly and corrupt the other.
 *
 *  Validated against gen_qball_boost: boosting a static ball reproduces the
 *  purpose-built generator's momentum. Residual error is O(v^2) from the
 *  contraction (gamma - 1 = 2.0e-4 at v = 0.02) and from local
 *  non-monochromaticity.
 *
 *  THE GAUGE SECTOR is deliberately NOT transformed. Under a boost the links
 *  acquire an A_0 piece and the moving charge sources a magnetic field, both
 *  O(v); representing them correctly would break the kernel's temporal gauge.
 *  Instead the E columns are zeroed and the kernel's own init Gauss projection
 *  rebuilds E from the boosted charge density -- the same path every gauged
 *  seed in this program uses. The omitted initial B carries energy O(v^2)
 *  relative to the field energy (4e-4 at v = 0.02) and regenerates within a
 *  light-crossing time. This is a documented approximation, not an oversight:
 *  the alternative is a gauge-transformed seed that the kernel cannot ingest.
 *
 *  VALIDATION. The tool prints the momentum it expects the boosted state to
 *  carry (P ~ E v) so it can be checked against sfa_momentum before any GPU
 *  time is spent.
 *
 *  Build: gcc -O3 -march=native -fopenmp -o sfa_boost sfa_boost.c -lzstd -lm
 *  Usage: sfa_boost in.sfa --frame N --v VX --out out.sfa [--keep-gauge]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#define SFA_IMPLEMENTATION
#include "../../sfa/format/sfa.h"

static int find_col(const SFA *s, const char *name) {
    for (uint32_t c = 0; c < s->n_columns; c++)
        if (strncmp(s->columns[c].name, name, sizeof(s->columns[c].name)) == 0)
            return (int)c;
    return -1;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "usage: %s in.sfa --frame N --v VX --out out.sfa "
                        "[--keep-gauge]\n", argv[0]);
        return 1;
    }
    const char *in = argv[1], *out = NULL;
    int frame = 0, keep_gauge = 0;
    double vx = 0.0;
    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "--frame") && i + 1 < argc) frame = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--v") && i + 1 < argc) vx = atof(argv[++i]);
        else if (!strcmp(argv[i], "--out") && i + 1 < argc) out = argv[++i];
        else if (!strcmp(argv[i], "--keep-gauge")) keep_gauge = 1;
        else { fprintf(stderr, "unknown arg %s\n", argv[i]); return 1; }
    }
    if (!out) { fprintf(stderr, "--out required\n"); return 1; }
    if (fabs(vx) >= 0.5) { fprintf(stderr, "|v| must be < 0.5 (O(v^2) expansion)\n"); return 1; }

    SFA *s = sfa_open(in);
    if (!s) { fprintf(stderr, "cannot open %s\n", in); return 1; }
    int N = (int)s->Nx, NN = N * N;
    long N3 = (long)N * N * N;
    double L = s->Lx, dx = 2.0 * L / (N - 1);
    double gam = 1.0 / sqrt(1.0 - vx * vx);

    uint64_t bufsz = 0;
    for (uint32_t c = 0; c < s->n_columns; c++)
        bufsz += s->N_total * sfa_dtype_size[s->columns[c].dtype];
    void *buf = malloc(bufsz);
    if (!buf) { fprintf(stderr, "cannot allocate %.1f GB\n", bufsz / 1e9); return 1; }
    if (sfa_read_frame(s, frame, buf) != 0) { fprintf(stderr, "frame read failed\n"); return 1; }

    fprintf(stderr, "sfa_boost: N=%d L=%g dx=%.6f  frame %d (t=%.2f)  v=%.4f "
                    "gamma=%.8f  (contraction = %.2e of the box)\n",
            N, L, dx, frame, sfa_frame_time(s, frame), vx, gam, gam - 1.0);

    /* every column must be f32 -- a boosted seed is a restart, not a view */
    for (uint32_t c = 0; c < s->n_columns; c++)
        if (s->columns[c].dtype != SFA_F32) {
            fprintf(stderr, "error: column '%s' is %s, not f32; boost the "
                            "restartable output, not an f16 viewing copy\n",
                    s->columns[c].name, sfa_dtype_name(s->columns[c].dtype));
            return 1;
        }

    float *base[64];
    uint64_t off = 0;
    for (uint32_t c = 0; c < s->n_columns; c++) {
        base[c] = (float *)((uint8_t *)buf + off);
        off += s->N_total * sizeof(float);
    }

    /* The four complex sectors: each is (u, v, udot, vdot) and must be
     * transformed TOGETHER, because a boost acts on the complex phase. */
    struct { const char *u, *v, *ud, *vd; } sec[6] = {
        {"phi_x","phiim_x","phi_vx","phiim_vx"},
        {"phi_y","phiim_y","phi_vy","phiim_vy"},
        {"phi_z","phiim_z","phi_vz","phiim_vz"},
        {"theta_x","thetaim_x","theta_vx","thetaim_vx"},
        {"theta_y","thetaim_y","theta_vy","thetaim_vy"},
        {"theta_z","thetaim_z","theta_vz","thetaim_vz"},
    };

    float *tu = malloc(sizeof(float) * N3), *tv = malloc(sizeof(float) * N3);
    float *td = malloc(sizeof(float) * N3), *te = malloc(sizeof(float) * N3);
    if (!tu || !tv || !td || !te) { fprintf(stderr, "alloc\n"); return 1; }

    int done = 0;
    for (int p = 0; p < 6; p++) {
        int cu = find_col(s, sec[p].u), cv = find_col(s, sec[p].v);
        int cud = find_col(s, sec[p].ud), cvd = find_col(s, sec[p].vd);
        if (cu < 0 || cv < 0 || cud < 0 || cvd < 0) continue;
        float *U = base[cu], *V = base[cv], *UD = base[cud], *VD = base[cvd];
        memcpy(tu, U, sizeof(float) * N3); memcpy(tv, V, sizeof(float) * N3);
        memcpy(td, UD, sizeof(float) * N3); memcpy(te, VD, sizeof(float) * N3);
        /* Amplitude floor for the local clock. In the near-vacuum halo n2 is
         * ~machine noise and omega_loc = (u vdot - v udot)/n2 is meaningless.
         * Rotating by a NOISY omega_loc(x) is not harmless: the induced phase
         * is -gamma omega_loc v x, so a noisy omega_loc creates a spurious
         * phase GRADIENT proportional to x d_x(omega_loc), and integrated over
         * a large halo that injects real momentum. Measured: without the floor
         * the boost over-delivered momentum by 10%. Below the floor no rotation
         * is applied -- that region carries negligible energy, so leaving it
         * unboosted costs far less than boosting it wrongly. */
        double n2max = 0.0;
        #pragma omp parallel for reduction(max:n2max) schedule(static)
        for (long q = 0; q < N3; q++) {
            double nn = (double)tu[q]*tu[q] + (double)tv[q]*tv[q];
            if (nn > n2max) n2max = nn;
        }
        const double n2floor = 1e-6 * n2max;
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < N3; idx++) {
            int i = (int)(idx / NN), j = (int)((idx / N) % N), k = (int)(idx % N);
            double x = -L + i * dx;
            /* Lorentz-contracted sampling point (linear interpolation) */
            double fi = (gam * x + L) / dx;
            int i0 = (int)floor(fi);
            double w = fi - i0;
            int i0c = i0 < 0 ? 0 : (i0 > N - 1 ? N - 1 : i0);
            int i1c = i0 + 1 < 0 ? 0 : (i0 + 1 > N - 1 ? N - 1 : i0 + 1);
            if (i0c == i1c) w = 0.0;
            long a0 = (long)i0c * NN + (long)j * N + k;
            long a1 = (long)i1c * NN + (long)j * N + k;
            double u  = (1-w)*tu[a0] + w*tu[a1];
            double v  = (1-w)*tv[a0] + w*tv[a1];
            double ud = (1-w)*td[a0] + w*td[a1];
            double vd = (1-w)*te[a0] + w*te[a1];
            int im = i0c > 0 ? i0c - 1 : i0c, ip = i0c < N - 1 ? i0c + 1 : i0c;
            double denom = (ip - im) * dx;
            double dxu = (tu[(long)ip*NN + (long)j*N + k] - tu[(long)im*NN + (long)j*N + k]) / denom;
            double dxv = (tv[(long)ip*NN + (long)j*N + k] - tv[(long)im*NN + (long)j*N + k]) / denom;

            /* LOCAL clock. A boost evaluates the rest-frame field at the shifted
             * rest time t = -gamma v x, which for a locally monochromatic field
             * is an EXACT phase rotation by omega_loc * t. Expanding that
             * rotation to first order (as a naive Taylor boost does) fails at
             * large |x|, where gamma*v*x*omega is order 1 -- and it silently
             * loses half the momentum. Using the local clock also handles a
             * state with SEVERAL clocks (the X10c ball and its opposite-charge
             * cloud rotate at different, opposite-sign frequencies), which a
             * single global omega cannot. */
            double n2 = u*u + v*v;
            double wloc = (n2 > n2floor) ? (u*vd - v*ud) / n2 : 0.0;
            /* clamp to the physical band: bound states sit below m and free
             * radiation on this lattice runs to a few m; anything outside is
             * numerical noise, not a clock */
            if (wloc > 3.0) wloc = 3.0; else if (wloc < -3.0) wloc = -3.0;
            double ph = -gam * wloc * vx * x;
            double cph = cos(ph), sph = sin(ph);

            /* boosted velocity, before the phase rotation */
            double bu = gam * (ud - gam * vx * dxu);
            double bv = gam * (vd - gam * vx * dxv);

            U[idx]  = (float)(u * cph - v * sph);
            V[idx]  = (float)(u * sph + v * cph);
            UD[idx] = (float)(bu * cph - bv * sph);
            VD[idx] = (float)(bu * sph + bv * cph);
        }
        done++;
    }
    fprintf(stderr, "            boosted %d complex sectors (local-clock exact rotation)\n", done);

    if (!keep_gauge) {
        int nz = 0;
        const char *ec[3] = {"E_x", "E_y", "E_z"};
        for (int d = 0; d < 3; d++) {
            int c = find_col(s, ec[d]);
            if (c >= 0) { memset(base[c], 0, sizeof(float) * N3); nz++; }
        }
        if (nz) fprintf(stderr, "            zeroed %d E columns -- the kernel's "
                                "init Gauss projection rebuilds E from the\n"
                                "            boosted charge density (links left "
                                "untouched)\n", nz);
    } else {
        fprintf(stderr, "            --keep-gauge: E and links passed through "
                        "unchanged (only valid if the source was unmagnetised)\n");
    }

    /* expected momentum, for validation against sfa_momentum before spending GPU */
    fprintf(stderr, "            expected momentum of the boosted state: "
                    "P_x ~ E_rest * gamma * v\n"
                    "            -> run sfa_momentum on the output and check "
                    "P_x/E ~ %.5f before committing GPU time\n", vx);

    /* write */
    SFA *o = sfa_create(out, N, N, N, L, L, L, s->dt);
    for (uint32_t c = 0; c < s->n_columns; c++)
        sfa_add_column(o, s->columns[c].name, s->columns[c].dtype,
                       s->columns[c].semantic, s->columns[c].component);
    SFA_KVMDSet kv[16];
    int n_kv = sfa_read_kvmd(s, kv, 16);
    for (int i = 0; i < n_kv; i++) {
        const char *ks[128], *vs[128];
        for (int p = 0; p < kv[i].n_pairs; p++) { ks[p] = kv[i].keys[p]; vs[p] = kv[i].values[p]; }
        sfa_add_kvmd(o, kv[i].set_id, kv[i].first_frame, kv[i].frame_count,
                     ks, vs, kv[i].n_pairs);
    }
    sfa_finalize_header(o);
    /* sfa_write_frame takes an array of per-column pointers, not the flat
     * frame buffer that sfa_read_frame fills */
    void *cols[64];
    for (uint32_t c = 0; c < s->n_columns; c++) cols[c] = (void *)base[c];
    if (sfa_write_frame(o, 0.0, cols) != 0) {
        fprintf(stderr, "error: sfa_write_frame failed\n"); return 1;
    }
    sfa_close(o);
    fprintf(stderr, "wrote %s (1 frame, %d columns, f32)\n", out, s->n_columns);
    free(buf); free(tu); free(tv); free(td); free(te);
    return 0;
}
