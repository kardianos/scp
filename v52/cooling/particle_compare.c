/*  particle_compare.c — Compare soliton structure at two times in an SFA file.
 *
 *  Reads frames near t_a and t_b, computes cluster analysis on each,
 *  reports changes in mass, shape, binding energy, and centroid.
 *
 *  Build: gcc -O3 -fopenmp -o particle_compare particle_compare.c -lzstd -lm
 *  Usage: particle_compare file.sfa [t_a] [t_b]
 *         defaults: t_a = midpoint, t_b = last frame
 */
#define SFA_IMPLEMENTATION
#include "../../sfa/format/sfa.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

typedef struct {
    double mass;        /* sum |P| * dV */
    double E_pot;       /* sum V(P) * dV */
    double phi_max;     /* max |phi| */
    double P_peak;      /* max |P| */
    double cx, cy, cz;  /* centroid */
    double rms_r;       /* RMS radius from centroid */
    double theta_rms;   /* theta RMS within cluster */
    int nvox;           /* voxels above threshold */
    double phi_sq_sum;  /* sum phi² * dV */
    double theta_sq_sum;/* sum theta² * dV */
} ParticleStats;

static ParticleStats analyze_frame(SFA *s, uint32_t fi, double dx) {
    ParticleStats ps = {0};
    int N = s->Nx;
    long N3 = (long)N*N*N;
    int NN = N*N;
    double L = s->Lx;
    double dV = dx*dx*dx;
    double mu = -41.345, kappa = 50.0;

    void *buf = malloc(s->frame_bytes);
    if (sfa_read_frame(s, fi, buf) < 0) {
        free(buf);
        return ps;
    }

    /* Get column offsets */
    uint64_t col_off[12] = {0};
    for (int c = 1; c < (int)s->n_columns && c < 12; c++)
        col_off[c] = col_off[c-1] + N3 * sfa_dtype_size[s->columns[c-1].dtype];

    /* First pass: find P threshold and centroid */
    double P_max = 0;
    for (long i = 0; i < N3; i++) {
        float p[6];
        for (int a = 0; a < 6 && a < (int)s->n_columns; a++) {
            int es = sfa_dtype_size[s->columns[a].dtype];
            if (es == 2) {
                uint16_t h; memcpy(&h, (uint8_t*)buf + col_off[a] + i*2, 2);
                p[a] = sfa_f16_to_f32(h);
            } else {
                memcpy(&p[a], (uint8_t*)buf + col_off[a] + i*4, 4);
            }
        }
        double P = fabs((double)p[0] * p[1] * p[2]);
        if (P > P_max) P_max = P;
    }

    double thresh = P_max * 0.01;  /* 1% of peak */
    if (thresh < 1e-6) thresh = 1e-6;

    /* Second pass: compute stats for voxels above threshold */
    double sum_x=0, sum_y=0, sum_z=0, sum_w=0;
    for (long idx = 0; idx < N3; idx++) {
        int i = (int)(idx/NN), j = (int)((idx/N)%N), k = (int)(idx%N);
        float p[6];
        for (int a = 0; a < 6 && a < (int)s->n_columns; a++) {
            int es = sfa_dtype_size[s->columns[a].dtype];
            if (es == 2) {
                uint16_t h; memcpy(&h, (uint8_t*)buf + col_off[a] + idx*2, 2);
                p[a] = sfa_f16_to_f32(h);
            } else {
                memcpy(&p[a], (uint8_t*)buf + col_off[a] + idx*4, 4);
            }
        }
        double P = fabs((double)p[0] * p[1] * p[2]);
        if (P < thresh) continue;

        double x = -L + i*dx, y = -L + j*dx, z = -L + k*dx;
        double ps2 = p[0]*p[0]+p[1]*p[1]+p[2]*p[2];
        double ts2 = p[3]*p[3]+p[4]*p[4]+p[5]*p[5];
        double Vpot = (mu/2.0)*P*P/(1.0+kappa*P*P);

        ps.mass += P * dV;
        ps.E_pot += Vpot * dV;
        ps.phi_sq_sum += ps2 * dV;
        ps.theta_sq_sum += ts2 * dV;
        if (sqrt(ps2) > ps.phi_max) ps.phi_max = sqrt(ps2);
        if (P > ps.P_peak) ps.P_peak = P;
        sum_x += P * x; sum_y += P * y; sum_z += P * z;
        sum_w += P;
        ps.nvox++;
    }

    if (sum_w > 0) {
        ps.cx = sum_x / sum_w;
        ps.cy = sum_y / sum_w;
        ps.cz = sum_z / sum_w;
    }

    /* Third pass: RMS radius */
    double sum_r2 = 0;
    for (long idx = 0; idx < N3; idx++) {
        int i = (int)(idx/NN), j = (int)((idx/N)%N), k = (int)(idx%N);
        float p[3];
        for (int a = 0; a < 3; a++) {
            int es = sfa_dtype_size[s->columns[a].dtype];
            if (es == 2) {
                uint16_t h; memcpy(&h, (uint8_t*)buf + col_off[a] + idx*2, 2);
                p[a] = sfa_f16_to_f32(h);
            } else {
                memcpy(&p[a], (uint8_t*)buf + col_off[a] + idx*4, 4);
            }
        }
        double P = fabs((double)p[0] * p[1] * p[2]);
        if (P < thresh) continue;
        double x = -L + i*dx, y = -L + j*dx, z = -L + k*dx;
        double dr2 = (x-ps.cx)*(x-ps.cx)+(y-ps.cy)*(y-ps.cy)+(z-ps.cz)*(z-ps.cz);
        sum_r2 += P * dr2;
    }
    if (sum_w > 0) ps.rms_r = sqrt(sum_r2 / sum_w);
    ps.theta_rms = (ps.nvox > 0) ? sqrt(ps.theta_sq_sum / (ps.nvox * dV)) : 0;

    free(buf);
    return ps;
}

static uint32_t find_frame_near_time(SFA *s, double target_t) {
    uint32_t best = 0;
    double best_dt = 1e30;
    for (uint32_t fi = 0; fi < s->total_frames; fi++) {
        SFA_L2Entry entry;
        if (sfa_find_frame(s, fi, &entry) < 0) continue;
        double dt = fabs(entry.time - target_t);
        if (dt < best_dt) { best_dt = dt; best = fi; }
    }
    return best;
}

static void print_stats(const char *label, ParticleStats *ps, double t) {
    printf("  %s (t=%.1f):\n", label, t);
    printf("    mass=%.4f  E_pot=%.4f  phi_max=%.4f  P_peak=%.6f\n",
           ps->mass, ps->E_pot, ps->phi_max, ps->P_peak);
    printf("    centroid=(%.2f, %.2f, %.2f)  rms_r=%.3f\n",
           ps->cx, ps->cy, ps->cz, ps->rms_r);
    printf("    nvox=%d  phi²_sum=%.4f  θ²_sum=%.6f  θ_rms=%.6f\n",
           ps->nvox, ps->phi_sq_sum, ps->theta_sq_sum, ps->theta_rms);
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s file.sfa [t_a] [t_b]\n", argv[0]);
        return 1;
    }
    SFA *s = sfa_open(argv[1]);
    /* Fix total_frames by scanning JTOP chain (sfa_open only reads header) */
    {
        uint64_t jt = s->first_jtop_offset;
        uint32_t total = 0;
        while (jt) {
            fseek(s->fp, (long)jt + 12, SEEK_SET);
            uint32_t mx, cur; uint64_t nxt;
            fread(&mx, 4, 1, s->fp); fread(&cur, 4, 1, s->fp); fread(&nxt, 8, 1, s->fp);
            for (uint32_t i = 0; i < cur; i++) {
                uint64_t jo; uint32_t ff, fc;
                fread(&jo, 8, 1, s->fp); fread(&ff, 4, 1, s->fp); fread(&fc, 4, 1, s->fp);
                total += fc;
            }
            jt = nxt;
        }
        if (total > s->total_frames) s->total_frames = total;
    }
    if (!s) { fprintf(stderr, "Cannot open %s\n", argv[1]); return 1; }

    int N = s->Nx;
    double L = s->Lx;
    double dx = 2.0 * L / (N - 1);

    /* Find last valid frame time for default t_b */
    SFA_L2Entry last_entry;
    sfa_find_frame(s, s->total_frames - 1, &last_entry);
    double t_end = last_entry.time;

    double t_a = argc > 2 ? atof(argv[2]) : t_end / 2;
    double t_b = argc > 3 ? atof(argv[3]) : t_end;

    printf("Particle comparison: %s\n", argv[1]);
    printf("  Grid: %d³, L=%.1f, dx=%.4f, %d frames, t_end=%.1f\n",
           N, L, dx, s->total_frames, t_end);

    /* Find frames nearest to target times */
    uint32_t fi_a = find_frame_near_time(s, t_a);
    uint32_t fi_b = find_frame_near_time(s, t_b);

    SFA_L2Entry ea, eb;
    sfa_find_frame(s, fi_a, &ea);
    sfa_find_frame(s, fi_b, &eb);
    printf("  Comparing frame %d (t=%.2f) vs frame %d (t=%.2f)\n\n",
           fi_a, ea.time, fi_b, eb.time);

    ParticleStats pa = analyze_frame(s, fi_a, dx);
    ParticleStats pb = analyze_frame(s, fi_b, dx);

    print_stats("State A", &pa, ea.time);
    print_stats("State B", &pb, eb.time);

    /* Deltas */
    printf("\n  Changes (B - A):\n");
    if (pa.mass > 0) {
        printf("    mass:    %+.4f (%+.1f%%)\n", pb.mass - pa.mass, 100*(pb.mass-pa.mass)/pa.mass);
        printf("    E_pot:   %+.4f (%+.1f%%)\n", pb.E_pot - pa.E_pot,
               pa.E_pot != 0 ? 100*(pb.E_pot-pa.E_pot)/fabs(pa.E_pot) : 0);
        printf("    phi_max: %+.4f (%+.1f%%)\n", pb.phi_max - pa.phi_max, 100*(pb.phi_max-pa.phi_max)/pa.phi_max);
        printf("    P_peak:  %+.6f (%+.1f%%)\n", pb.P_peak - pa.P_peak, 100*(pb.P_peak-pa.P_peak)/pa.P_peak);
        printf("    rms_r:   %+.4f (%+.1f%%)\n", pb.rms_r - pa.rms_r, 100*(pb.rms_r-pa.rms_r)/pa.rms_r);
        double d_centroid = sqrt((pb.cx-pa.cx)*(pb.cx-pa.cx)+(pb.cy-pa.cy)*(pb.cy-pa.cy)+(pb.cz-pa.cz)*(pb.cz-pa.cz));
        printf("    centroid drift: %.4f\n", d_centroid);
        printf("    nvox:    %+d (%+.1f%%)\n", pb.nvox - pa.nvox, 100.0*(pb.nvox-pa.nvox)/pa.nvox);
    }

    sfa_close(s);
    return 0;
}
