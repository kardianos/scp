/*  v54/points.c — Point-based medium, v4 (strict locality)
 *
 *  Each point interacts ONLY with touching neighbors (r < 2*r_point).
 *  Forces propagate at finite speed through chains of grip connections.
 *  No action at a distance.
 *
 *  Build: gcc -O3 -march=native -fopenmp -o points points.c -lm
 *  Usage: ./points [N] [steps] [traj_file] [no_render] [conservative]
 *    traj_file    : path to write trajectory (default "traj.bin"; "-" disables)
 *    no_render    : "1" to skip PPM rendering (default 0)
 *    conservative : "1" to disable viscous drag + angular damping
 *                   (for frequency-resonance experiments) (default 0)
 *
 *  Trajectory file format (little-endian):
 *    magic[4]        "PNTS"
 *    version         u32 (=1)
 *    n_total         u32  total particle count
 *    n_tracked       u32  number of sampled particles
 *    n_samples       u32  number of sample rows
 *    steps_per_samp  u32
 *    dt              f64
 *    r_touch         f64
 *    L               f64
 *    tracked_ids     u32[n_tracked]
 *  For each sample (i = 0..n_samples-1):
 *    step            u32
 *    _pad            u32
 *    t               f64
 *    samples[n_tracked]:
 *       x,y,z           f32 x 3
 *       vx,vy,vz        f32 x 3
 *       n_neighbors     u8
 *       _pad[3]         u8 x 3
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <omp.h>

#define PI 3.14159265358979323846

typedef struct {
    double x, y, z;
    double vx, vy, vz;
    double dx, dy, dz;     /* unit direction */
    double wx, wy, wz;     /* angular velocity */
} Point;

/* ---- Spatial hash ---- */
#define HASH_SIZE (1 << 18)
#define MAX_PER_CELL 32
typedef struct { int count; int ids[MAX_PER_CELL]; } Cell;
typedef struct { Cell *cells; double cs; } SHash;

static inline int shash(int a, int b, int c) {
    return (int)((unsigned)(a*73856093u ^ b*19349663u ^ c*83492791u) & (HASH_SIZE-1));
}
static SHash *sh_new(double cs) {
    SHash *s = calloc(1, sizeof(SHash));
    s->cells = calloc(HASH_SIZE, sizeof(Cell));
    s->cs = cs;
    return s;
}
static void sh_build(SHash *s, Point *p, int N) {
    for (int i = 0; i < HASH_SIZE; i++) s->cells[i].count = 0;
    double inv = 1.0 / s->cs;
    for (int i = 0; i < N; i++) {
        int h = shash((int)floor(p[i].x*inv),(int)floor(p[i].y*inv),(int)floor(p[i].z*inv));
        if (s->cells[h].count < MAX_PER_CELL) s->cells[h].ids[s->cells[h].count++] = i;
    }
}
static void sh_free(SHash *s) { free(s->cells); free(s); }

static inline void normalize(double *x, double *y, double *z) {
    double m = sqrt(*x**x + *y**y + *z**z);
    if (m < 1e-12) { *x=1;*y=0;*z=0; return; }
    *x/=m; *y/=m; *z/=m;
}

/* Global flag: if nonzero, skip viscous drag terms (conservative run) */
static int g_conservative = 0;

/* ---- Force computation: STRICT LOCALITY ---- */
static void compute_forces(Point *pts, int N, SHash *sh,
                           double *fx, double *fy, double *fz,
                           double *tx, double *ty, double *tz,
                           double r_touch) {
    double inv_cs = 1.0 / sh->cs;
    double r_touch2 = r_touch * r_touch;

    memset(fx,0,N*8); memset(fy,0,N*8); memset(fz,0,N*8);
    memset(tx,0,N*8); memset(ty,0,N*8); memset(tz,0,N*8);

    #pragma omp parallel for schedule(dynamic, 128)
    for (int i = 0; i < N; i++) {
        Point *pi = &pts[i];
        int ci=(int)floor(pi->x*inv_cs), cj=(int)floor(pi->y*inv_cs), ck=(int)floor(pi->z*inv_cs);

        double fi_x=0,fi_y=0,fi_z=0;
        double ti_x=0,ti_y=0,ti_z=0;

        for (int a=-1;a<=1;a++) for (int b=-1;b<=1;b++) for (int c=-1;c<=1;c++) {
            Cell *cell = &sh->cells[shash(ci+a,cj+b,ck+c)];
            for (int n=0; n<cell->count; n++) {
                int j = cell->ids[n];
                if (j == i) continue;
                Point *pj = &pts[j];

                double rx=pj->x-pi->x, ry=pj->y-pi->y, rz=pj->z-pi->z;
                double r2 = rx*rx+ry*ry+rz*rz;
                if (r2 > r_touch2 || r2 < 1e-14) continue;

                double r = sqrt(r2), inv_r = 1.0/r;
                double ux=rx*inv_r, uy=ry*inv_r, uz=rz*inv_r;

                /* ======== 1. Contact repulsion/attraction ======== */
                /* Equilibrium at r_eq (half of r_touch).
                 * Below r_eq: repel. Above r_eq: attract.
                 * Both drop to zero at r=0 and r=r_touch. */
                double r_eq = r_touch * 0.5;
                double dr = r - r_eq;
                double f_contact = -2.0 * dr;  /* spring toward r_eq */
                /* Taper to zero at r_touch */
                double t_cut = r / r_touch;
                double taper = (1.0 - t_cut) * (1.0 - t_cut);
                f_contact *= taper;

                /* ======== 2. Grip: directional coupling ======== */
                double di_dot_r = pi->dx*ux + pi->dy*uy + pi->dz*uz;

                /* Chiral grip: (d_i × r̂) · d_j */
                double cx = pi->dy*uz - pi->dz*uy;
                double cy = pi->dz*ux - pi->dx*uz;
                double cz = pi->dx*uy - pi->dy*ux;
                double chiral = cx*pj->dx + cy*pj->dy + cz*pj->dz;

                double grip = 0.5*fabs(di_dot_r) + 0.5*fmax(chiral, 0.0);

                /* Grip modulates the spring: aligned neighbors hold tighter */
                f_contact *= (0.3 + 0.7 * grip);

                fi_x += f_contact * ux;
                fi_y += f_contact * uy;
                fi_z += f_contact * uz;

                /* ======== 3. Drag: velocity transmission through contact ======== */
                if (!g_conservative) {
                    double dvx=pj->vx-pi->vx, dvy=pj->vy-pi->vy, dvz=pj->vz-pi->vz;
                    double dv_along = dvx*ux + dvy*uy + dvz*uz;

                    /* Anisotropic: transmits better along direction */
                    double transmit = 0.3 + 0.7*fabs(di_dot_r);
                    double k_drag = 1.5 * grip * transmit * taper;

                    fi_x += k_drag * dv_along * ux;
                    fi_y += k_drag * dv_along * uy;
                    fi_z += k_drag * dv_along * uz;

                    /* Shear drag (transverse) */
                    double dvt_x=dvx-dv_along*ux, dvt_y=dvy-dv_along*uy, dvt_z=dvz-dv_along*uz;
                    double k_shear = 0.3 * grip * taper;
                    fi_x += k_shear * dvt_x;
                    fi_y += k_shear * dvt_y;
                    fi_z += k_shear * dvt_z;
                }

                /* ======== 4. Torque: align for better grip ======== */
                /* Maximize chiral: gradient w.r.t d_i is (r̂ × d_j) */
                double want_x = uy*pj->dz - uz*pj->dy;
                double want_y = uz*pj->dx - ux*pj->dz;
                double want_z = ux*pj->dy - uy*pj->dx;

                double tq_x = pi->dy*want_z - pi->dz*want_y;
                double tq_y = pi->dz*want_x - pi->dx*want_z;
                double tq_z = pi->dx*want_y - pi->dy*want_x;

                /* Parallel alignment torque */
                double par_x = pi->dy*pj->dz - pi->dz*pj->dy;
                double par_y = pi->dz*pj->dx - pi->dx*pj->dz;
                double par_z = pi->dx*pj->dy - pi->dy*pj->dx;

                double tq_str = 0.8 * taper;
                ti_x += tq_str * (0.7*tq_x + 0.3*par_x);
                ti_y += tq_str * (0.7*tq_y + 0.3*par_y);
                ti_z += tq_str * (0.7*tq_z + 0.3*par_z);
            }
        }
        fx[i]=fi_x; fy[i]=fi_y; fz[i]=fi_z;
        tx[i]=ti_x; ty[i]=ti_y; tz[i]=ti_z;
    }
}

/* ---- Angular integration ---- */
static void apply_angular(Point *pts, int N, double *tx, double *ty, double *tz, double dt) {
    for (int i = 0; i < N; i++) {
        pts[i].wx += dt * tx[i];
        pts[i].wy += dt * ty[i];
        pts[i].wz += dt * tz[i];
        if (!g_conservative) {
            /* Gentle angular damping */
            pts[i].wx *= 0.998; pts[i].wy *= 0.998; pts[i].wz *= 0.998;
        }
        /* Rotate direction */
        double w = sqrt(pts[i].wx*pts[i].wx+pts[i].wy*pts[i].wy+pts[i].wz*pts[i].wz);
        if (w > 1e-12) {
            double angle = w * dt;
            double kx=pts[i].wx/w, ky=pts[i].wy/w, kz=pts[i].wz/w;
            double c=cos(angle), s=sin(angle);
            double d=kx*pts[i].dx+ky*pts[i].dy+kz*pts[i].dz;
            double cx_=ky*pts[i].dz-kz*pts[i].dy;
            double cy_=kz*pts[i].dx-kx*pts[i].dz;
            double cz_=kx*pts[i].dy-ky*pts[i].dx;
            pts[i].dx = pts[i].dx*c + cx_*s + kx*d*(1-c);
            pts[i].dy = pts[i].dy*c + cy_*s + ky*d*(1-c);
            pts[i].dz = pts[i].dz*c + cz_*s + kz*d*(1-c);
            normalize(&pts[i].dx, &pts[i].dy, &pts[i].dz);
        }
    }
}

/* ---- Diagnostics ---- */
static void diag(Point *pts, int N, SHash *sh, int step, double t, double r_touch) {
    double ke=0, ke_r=0, r2=0, cx=0,cy=0,cz=0;
    for (int i=0;i<N;i++) {
        ke += pts[i].vx*pts[i].vx+pts[i].vy*pts[i].vy+pts[i].vz*pts[i].vz;
        ke_r += pts[i].wx*pts[i].wx+pts[i].wy*pts[i].wy+pts[i].wz*pts[i].wz;
        cx+=pts[i].x; cy+=pts[i].y; cz+=pts[i].z;
    }
    ke*=0.5/N; ke_r*=0.5/N; cx/=N; cy/=N; cz/=N;
    for (int i=0;i<N;i++) {
        double dx=pts[i].x-cx,dy=pts[i].y-cy,dz=pts[i].z-cz;
        r2+=dx*dx+dy*dy+dz*dz;
    }

    /* Sample chiral alignment + neighbor count */
    double inv_cs=1.0/sh->cs;
    double avg_ch=0; int n_pairs=0;
    double avg_nn=0; int n_samp=0;
    for (int i=0;i<N;i+=50) {
        int ci_=(int)floor(pts[i].x*inv_cs),cj_=(int)floor(pts[i].y*inv_cs),ck_=(int)floor(pts[i].z*inv_cs);
        int nn=0;
        for (int a=-1;a<=1;a++) for (int b=-1;b<=1;b++) for (int c=-1;c<=1;c++) {
            Cell *cell=&sh->cells[shash(ci_+a,cj_+b,ck_+c)];
            for (int n=0;n<cell->count;n++) {
                int j=cell->ids[n]; if(j==i) continue;
                double rx=pts[j].x-pts[i].x,ry=pts[j].y-pts[i].y,rz=pts[j].z-pts[i].z;
                double r=sqrt(rx*rx+ry*ry+rz*rz);
                if (r >= r_touch) continue;
                nn++;
                double inv_r=1.0/r;
                double ux=rx*inv_r,uy=ry*inv_r,uz=rz*inv_r;
                double crx=pts[i].dy*uz-pts[i].dz*uy;
                double cry=pts[i].dz*ux-pts[i].dx*uz;
                double crz=pts[i].dx*uy-pts[i].dy*ux;
                avg_ch += crx*pts[j].dx+cry*pts[j].dy+crz*pts[j].dz;
                n_pairs++;
            }
        }
        avg_nn += nn;
        n_samp++;
    }
    if (n_pairs>0) avg_ch /= n_pairs;
    avg_nn /= n_samp;

    printf("t=%7.2f KE=%.4f KEr=%.4f rms=%.2f nn=%.1f chiral=%.4f\n",
           t, ke, ke_r, sqrt(r2/N), avg_nn, avg_ch);
}

/* ---- Trajectory sampler ---- */
typedef struct {
    FILE *f;
    uint32_t n_tracked;
    uint32_t *ids;
    uint32_t n_written;
} TrajWriter;

/* Count neighbors of particle i within r_touch (uses existing spatial hash) */
static int count_neighbors(Point *pts, int i, SHash *sh, double r_touch) {
    double inv_cs = 1.0 / sh->cs;
    double r2max = r_touch * r_touch;
    int ci = (int)floor(pts[i].x*inv_cs);
    int cj = (int)floor(pts[i].y*inv_cs);
    int ck = (int)floor(pts[i].z*inv_cs);
    int n = 0;
    for (int a=-1;a<=1;a++) for (int b=-1;b<=1;b++) for (int c=-1;c<=1;c++) {
        Cell *cell = &sh->cells[shash(ci+a,cj+b,ck+c)];
        for (int m=0; m<cell->count; m++) {
            int j = cell->ids[m];
            if (j == i) continue;
            double rx = pts[j].x - pts[i].x;
            double ry = pts[j].y - pts[i].y;
            double rz = pts[j].z - pts[i].z;
            double r2 = rx*rx + ry*ry + rz*rz;
            if (r2 < r2max && r2 > 1e-14) n++;
        }
    }
    return n;
}

static TrajWriter *traj_open(const char *path, int n_total, int n_tracked_target,
                              int steps_per_sample, double dt,
                              double r_touch, double L) {
    TrajWriter *tw = calloc(1, sizeof(TrajWriter));
    tw->f = fopen(path, "wb");
    if (!tw->f) { fprintf(stderr, "traj_open: cannot open %s\n", path); free(tw); return NULL; }
    /* Evenly-spaced IDs */
    int stride = n_total / n_tracked_target;
    if (stride < 1) stride = 1;
    int nt = 0;
    tw->ids = malloc(n_tracked_target * sizeof(uint32_t));
    for (int i = 0; i < n_total && nt < n_tracked_target; i += stride) {
        tw->ids[nt++] = (uint32_t)i;
    }
    tw->n_tracked = (uint32_t)nt;

    /* Header (n_samples written at close) */
    char magic[4] = {'P','N','T','S'};
    uint32_t ver = 1;
    uint32_t n_total_u = (uint32_t)n_total;
    uint32_t n_samples_placeholder = 0;
    uint32_t sps = (uint32_t)steps_per_sample;
    fwrite(magic, 1, 4, tw->f);
    fwrite(&ver, 4, 1, tw->f);
    fwrite(&n_total_u, 4, 1, tw->f);
    fwrite(&tw->n_tracked, 4, 1, tw->f);
    fwrite(&n_samples_placeholder, 4, 1, tw->f);
    fwrite(&sps, 4, 1, tw->f);
    fwrite(&dt, 8, 1, tw->f);
    fwrite(&r_touch, 8, 1, tw->f);
    fwrite(&L, 8, 1, tw->f);
    fwrite(tw->ids, 4, tw->n_tracked, tw->f);
    return tw;
}

static void traj_sample(TrajWriter *tw, Point *pts, SHash *sh,
                        int step, double t, double r_touch) {
    if (!tw) return;
    uint32_t step_u = (uint32_t)step;
    uint32_t pad_u = 0;
    fwrite(&step_u, 4, 1, tw->f);
    fwrite(&pad_u, 4, 1, tw->f);
    fwrite(&t, 8, 1, tw->f);
    for (uint32_t k = 0; k < tw->n_tracked; k++) {
        int i = (int)tw->ids[k];
        float rec[6];
        rec[0]=(float)pts[i].x; rec[1]=(float)pts[i].y; rec[2]=(float)pts[i].z;
        rec[3]=(float)pts[i].vx;rec[4]=(float)pts[i].vy;rec[5]=(float)pts[i].vz;
        fwrite(rec, 4, 6, tw->f);
        int nn = count_neighbors(pts, i, sh, r_touch);
        if (nn > 255) nn = 255;
        uint8_t buf[4] = { (uint8_t)nn, 0, 0, 0 };
        fwrite(buf, 1, 4, tw->f);
    }
    tw->n_written++;
}

static void traj_close(TrajWriter *tw) {
    if (!tw) return;
    /* Seek back to n_samples field (offset 4+4+4+4 = 16 bytes) and write it */
    fseek(tw->f, 16, SEEK_SET);
    fwrite(&tw->n_written, 4, 1, tw->f);
    fclose(tw->f);
    free(tw->ids);
    free(tw);
}

/* ---- PPM renderer ---- */
#define IMG_W 800
#define IMG_H 800

static void render_ppm(const char *pfx, int frame, Point *pts, int N) {
    char fn[256]; snprintf(fn,256,"%s_%04d.ppm",pfx,frame);
    FILE *f=fopen(fn,"wb"); if(!f) return;
    fprintf(f,"P6\n%d %d\n255\n",IMG_W,IMG_H);

    unsigned char *img = calloc(IMG_W*IMG_H*3,1);
    float *depth = malloc(IMG_W*IMG_H*sizeof(float));
    for (int i=0;i<IMG_W*IMG_H;i++) depth[i]=1e30f;

    double ax=35.264*PI/180.0, ay=45.0*PI/180.0;
    double cax=cos(ax),sax=sin(ax),cay=cos(ay),say=sin(ay);

    double max_e=1.0;
    for (int i=0;i<N;i++) {
        double e=fabs(pts[i].x); if(e>max_e)max_e=e;
        e=fabs(pts[i].y); if(e>max_e)max_e=e;
        e=fabs(pts[i].z); if(e>max_e)max_e=e;
    }
    double scale=(IMG_W*0.4)/max_e;

    for (int i=0;i<N;i++) {
        double x0=pts[i].x,y0=pts[i].y,z0=pts[i].z;
        double x1=cay*x0+say*z0, z1=-say*x0+cay*z0, y1=y0;
        double y2=cax*y1-sax*z1, z2=sax*y1+cax*z1;

        double dx0=pts[i].dx,dy0=pts[i].dy,dz0=pts[i].dz;
        double dx1=cay*dx0+say*dz0, dz1=-say*dx0+cay*dz0, dy1=dy0;
        double dy2=cax*dy1-sax*dz1;

        int sx=(int)(IMG_W/2+x1*scale), sy=(int)(IMG_H/2-y2*scale);
        int rad=(int)(2.0+1.5*(1.0-z2/(max_e*2)));
        if(rad<1)rad=1; if(rad>5)rad=5;

        int cr=(int)(128+127*dx1), cg=(int)(128+127*dy2), cb=(int)(128+127*(-dz1));
        if(cr<0)cr=0;if(cr>255)cr=255;
        if(cg<0)cg=0;if(cg>255)cg=255;
        if(cb<0)cb=0;if(cb>255)cb=255;

        double bright=0.4+0.6*(1.0-(z2+max_e)/(2*max_e));
        cr=(int)(cr*bright); cg=(int)(cg*bright); cb=(int)(cb*bright);

        for (int dy=-rad;dy<=rad;dy++)
        for (int dx=-rad;dx<=rad;dx++) {
            if(dx*dx+dy*dy>rad*rad) continue;
            int px=sx+dx,py=sy+dy;
            if(px<0||px>=IMG_W||py<0||py>=IMG_H) continue;
            int idx=py*IMG_W+px;
            if(z2<depth[idx]) {
                depth[idx]=(float)z2;
                img[idx*3+0]=(unsigned char)cr;
                img[idx*3+1]=(unsigned char)cg;
                img[idx*3+2]=(unsigned char)cb;
            }
        }
    }
    fwrite(img,1,IMG_W*IMG_H*3,f);
    fclose(f); free(img); free(depth);
}

int main(int argc, char **argv) {
    int N_target = argc > 1 ? atoi(argv[1]) : 10000;
    int steps = argc > 2 ? atoi(argv[2]) : 200000;
    const char *traj_path = argc > 3 ? argv[3] : "traj.bin";
    int no_render = argc > 4 ? atoi(argv[4]) : 0;
    g_conservative = argc > 5 ? atoi(argv[5]) : 0;
    if (g_conservative) printf("CONSERVATIVE MODE: viscous drag + angular damping disabled\n");
    double L = 10.0;
    double r_touch = 1.2;   /* points interact only within this distance */
    double dt = 0.002;
    int render_every = 1000;
    int diag_every = 500;
    int traj_every = 10;        /* sample trajectory every 10 steps */
    int n_tracked_target = 200; /* track 200 evenly-spaced particles */

    int N_max = N_target * 2;
    Point *pts = calloc(N_max, sizeof(Point));
    double *Fx=malloc(N_max*8),*Fy=malloc(N_max*8),*Fz=malloc(N_max*8);
    double *Tx=malloc(N_max*8),*Ty=malloc(N_max*8),*Tz=malloc(N_max*8);

    /* Init: lattice at ~r_touch/2 spacing (touching neighbors), random dirs */
    double spacing = r_touch * 0.55;  /* slightly more than r_eq = r_touch/2 */
    int nside = (int)(2.0 * L / spacing);
    int N = 0;
    srand(42);
    for (int i=0;i<nside&&N<N_max;i++)
    for (int j=0;j<nside&&N<N_max;j++)
    for (int k=0;k<nside&&N<N_max;k++) {
        double px = -L+(i+0.5)*spacing;
        double py = -L+(j+0.5)*spacing;
        double pz = -L+(k+0.5)*spacing;
        /* Only place within sphere of radius L to avoid corners */
        if (px*px+py*py+pz*pz > L*L) continue;
        pts[N].x = px + 0.02*spacing*(2.0*rand()/RAND_MAX-1);
        pts[N].y = py + 0.02*spacing*(2.0*rand()/RAND_MAX-1);
        pts[N].z = pz + 0.02*spacing*(2.0*rand()/RAND_MAX-1);
        double th=acos(2.0*rand()/RAND_MAX-1.0), ph=2*PI*(double)rand()/RAND_MAX;
        pts[N].dx=sin(th)*cos(ph); pts[N].dy=sin(th)*sin(ph); pts[N].dz=cos(th);
        /* Small thermal velocity */
        double vth=0.15;
        th=acos(2.0*rand()/RAND_MAX-1.0); ph=2*PI*(double)rand()/RAND_MAX;
        double spd=vth*(0.5+(double)rand()/RAND_MAX);
        pts[N].vx=spd*sin(th)*cos(ph); pts[N].vy=spd*sin(th)*sin(ph); pts[N].vz=spd*cos(th);
        pts[N].wx=pts[N].wy=pts[N].wz=0;
        N++;
    }
    printf("v54 local v4: N=%d spacing=%.3f r_touch=%.2f dt=%.4f\n", N, spacing, r_touch, dt);

    /* Hash cell size = r_touch so 27-cell search covers all contacts */
    SHash *sh = sh_new(r_touch);
    sh_build(sh, pts, N);
    compute_forces(pts, N, sh, Fx, Fy, Fz, Tx, Ty, Tz, r_touch);

    /* Open trajectory writer (unless disabled with "-") */
    TrajWriter *tw = NULL;
    if (traj_path && strcmp(traj_path, "-") != 0) {
        tw = traj_open(traj_path, N, n_tracked_target, traj_every, dt, r_touch, L);
        if (tw) {
            printf("trajectory: %s  tracked=%u  stride=%d steps  (%d samples expected)\n",
                   traj_path, tw->n_tracked, traj_every, steps/traj_every);
        }
    }

    for (int step = 0; step < steps; step++) {
        if (step % diag_every == 0) diag(pts, N, sh, step, step*dt, r_touch);
        if (!no_render && step % render_every == 0) render_ppm("frame", step/render_every, pts, N);
        if (tw && step % traj_every == 0) traj_sample(tw, pts, sh, step, step*dt, r_touch);

        double hdt = 0.5*dt;
        for (int i=0;i<N;i++) {
            pts[i].vx+=hdt*Fx[i]; pts[i].vy+=hdt*Fy[i]; pts[i].vz+=hdt*Fz[i];
        }
        for (int i=0;i<N;i++) {
            pts[i].x+=dt*pts[i].vx; pts[i].y+=dt*pts[i].vy; pts[i].z+=dt*pts[i].vz;
        }
        apply_angular(pts, N, Tx, Ty, Tz, dt);

        sh_build(sh, pts, N);
        compute_forces(pts, N, sh, Fx, Fy, Fz, Tx, Ty, Tz, r_touch);

        for (int i=0;i<N;i++) {
            pts[i].vx+=hdt*Fx[i]; pts[i].vy+=hdt*Fy[i]; pts[i].vz+=hdt*Fz[i];
        }

        /* Soft confinement */
        for (int i=0;i<N;i++) {
            double r=sqrt(pts[i].x*pts[i].x+pts[i].y*pts[i].y+pts[i].z*pts[i].z);
            if (r > L) {
                double f=-0.1*(r-L)/r;
                pts[i].vx+=dt*f*pts[i].x;
                pts[i].vy+=dt*f*pts[i].y;
                pts[i].vz+=dt*f*pts[i].z;
            }
        }
    }
    if (!no_render) render_ppm("frame", steps/render_every, pts, N);
    diag(pts, N, sh, steps, steps*dt, r_touch);
    if (tw) {
        printf("trajectory: wrote %u samples to %s\n", tw->n_written, traj_path);
        traj_close(tw);
    }
    printf("Done.\n");

    sh_free(sh); free(pts);
    free(Fx);free(Fy);free(Fz);free(Tx);free(Ty);free(Tz);
    return 0;
}
