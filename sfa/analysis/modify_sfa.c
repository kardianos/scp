/*  modify_sfa.c — Read an SFA frame, apply modifications, write new seed SFA
 *
 *  Modifications:
 *    --add-braid cx cy cz kz_angle  — overlay a braid at position, rotated
 *    --add-oscillon cx cy cz A sigma — overlay a Gaussian pulse
 *    --perturb amplitude              — add random noise (fraction of phi_max)
 *    --rotate axis angle              — rotate all fields around axis (x/y/z)
 *    --scale factor                   — scale all field amplitudes
 *    --flip-chirality                 — reverse k_z → -k_z (flip z velocities)
 *
 *  Build: gcc -O3 -o modify_sfa modify_sfa.c -lzstd -lm
 *  Usage: ./modify_sfa input.sfa --frame -1 --add-braid 10 0 0 90 -o modified.sfa
 */

#define SFA_IMPLEMENTATION
#include "../format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#define PI 3.14159265358979323846
#define NFIELDS 3

static double MASS2 = 2.25;

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s input.sfa [--frame N] [--add-braid cx cy cz angle] "
                "[--perturb amp] [--flip-chirality] -o output.sfa\n", argv[0]);
        return 1;
    }

    const char *in_path = argv[1];
    const char *out_path = "modified.sfa";
    int frame_idx = -1;

    /* Parse modifications to apply */
    typedef struct { int type; double params[8]; } Mod;
    Mod mods[16]; int n_mods = 0;
    /* type: 0=add-braid, 1=add-oscillon, 2=perturb, 3=scale, 4=flip-chirality */

    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "--frame") && i+1<argc) frame_idx = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-o") && i+1<argc) out_path = argv[++i];
        else if (!strcmp(argv[i], "--add-braid") && i+4<argc) {
            Mod *m = &mods[n_mods++]; m->type = 0;
            m->params[0]=atof(argv[++i]); m->params[1]=atof(argv[++i]);
            m->params[2]=atof(argv[++i]); m->params[3]=atof(argv[++i]);
        }
        else if (!strcmp(argv[i], "--add-oscillon") && i+5<argc) {
            Mod *m = &mods[n_mods++]; m->type = 1;
            m->params[0]=atof(argv[++i]); m->params[1]=atof(argv[++i]);
            m->params[2]=atof(argv[++i]); m->params[3]=atof(argv[++i]);
            m->params[4]=atof(argv[++i]);
        }
        else if (!strcmp(argv[i], "--perturb") && i+1<argc) {
            Mod *m = &mods[n_mods++]; m->type = 2; m->params[0]=atof(argv[++i]);
        }
        else if (!strcmp(argv[i], "--scale") && i+1<argc) {
            Mod *m = &mods[n_mods++]; m->type = 3; m->params[0]=atof(argv[++i]);
        }
        else if (!strcmp(argv[i], "--flip-chirality")) {
            Mod *m = &mods[n_mods++]; m->type = 4;
        }
    }

    /* Open input SFA */
    SFA *in = sfa_open(in_path);
    if (!in) { fprintf(stderr, "Cannot open %s\n", in_path); return 1; }
    int N = in->Nx; double L = in->Lx;
    long N3 = (long)N*N*N;
    int NN = N*N;

    if (frame_idx < 0) frame_idx = in->total_frames + frame_idx;
    fprintf(stderr, "Loading %s frame %d (N=%d L=%.1f)\n", in_path, frame_idx, N, L);

    /* Read frame */
    void *buf = malloc(in->frame_bytes);
    sfa_read_frame(in, frame_idx, buf);

    /* Extract all 12 arrays as f64 */
    double *arrays[12]; memset(arrays, 0, sizeof(arrays));
    for (int i = 0; i < 12; i++) arrays[i] = (double*)calloc(N3, sizeof(double));

    uint64_t off = 0;
    for (int c = 0; c < in->n_columns; c++) {
        int dtype=in->columns[c].dtype, sem=in->columns[c].semantic, comp=in->columns[c].component;
        int es=sfa_dtype_size[dtype]; uint8_t *src=(uint8_t*)buf+off;
        int slot=-1;
        if(sem==SFA_POSITION&&comp<3) slot=comp;
        else if(sem==SFA_ANGLE&&comp<3) slot=3+comp;
        else if(sem==SFA_VELOCITY&&comp<3) slot=6+comp;
        else if(sem==SFA_VELOCITY&&comp>=3&&comp<6) slot=9+comp-3;
        if (slot >= 0) {
            if(dtype==SFA_F64) memcpy(arrays[slot],src,N3*8);
            else if(dtype==SFA_F32) for(long i=0;i<N3;i++) arrays[slot][i]=(double)((float*)src)[i];
            else if(dtype==SFA_F16) for(long i=0;i<N3;i++){
                uint16_t h=((uint16_t*)src)[i];int e=(h>>10)&0x1F;uint16_t m=h&0x3FF;
                if(e==0)arrays[slot][i]=0;else{float fv;
                uint32_t x=((uint32_t)(h&0x8000)<<16)|((uint32_t)(e-15+127)<<23)|((uint32_t)m<<13);
                memcpy(&fv,&x,4);arrays[slot][i]=(double)fv;}
            }
        }
        off += (uint64_t)N3*es;
    }
    free(buf);

    /* Apply modifications */
    double dx = 2.0*L/(N-1);
    srand(time(NULL));

    for (int mi = 0; mi < n_mods; mi++) {
        Mod *mod = &mods[mi];
        fprintf(stderr, "Applying mod type %d\n", mod->type);

        if (mod->type == 0) { /* add-braid */
            double bcx=mod->params[0], bcy=mod->params[1], bcz=mod->params[2];
            double angle_deg=mod->params[3];
            double angle = angle_deg * PI / 180.0;
            double ca=cos(angle), sa=sin(angle);
            double delta[3]={0, 3.0005, 4.4325};
            double A=0.4, R_tube=3.0, ellip=0.3325;
            double kw=PI/L, omega=sqrt(kw*kw+MASS2);
            double sx=1+ellip, sy=1-ellip;
            double inv2R2=1.0/(2*R_tube*R_tube);

            for (int i=0;i<N;i++) { double x=-L+i*dx;
            for (int j=0;j<N;j++) { double y=-L+j*dx;
            for (int k=0;k<N;k++) { double z=-L+k*dx;
                long idx=(long)i*NN+j*N+k;
                /* Rotate coordinates: rotate around y-axis by angle */
                double xr = (x-bcx)*ca + (z-bcz)*sa;
                double yr = y - bcy;
                double zr = -(x-bcx)*sa + (z-bcz)*ca;
                double r2e = xr*xr/(sx*sx) + yr*yr/(sy*sy);
                double env = exp(-r2e*inv2R2);
                for (int a=0;a<NFIELDS;a++) {
                    double ph = kw*zr + delta[a];
                    arrays[a][idx] += A*env*cos(ph);      /* phi */
                    arrays[6+a][idx] += omega*A*env*sin(ph); /* phi_vel */
                }
            }}}
            fprintf(stderr, "  Added braid at (%.1f,%.1f,%.1f) angle=%.0f deg A=%.2f\n",
                    bcx, bcy, bcz, angle_deg, A);
        }
        else if (mod->type == 1) { /* add-oscillon */
            double ocx=mod->params[0], ocy=mod->params[1], ocz=mod->params[2];
            double A=mod->params[3], sigma=mod->params[4];
            double delta[3]={0, 3.0005, 4.4325};
            for (int i=0;i<N;i++) { double x=-L+i*dx;
            for (int j=0;j<N;j++) { double y=-L+j*dx;
            for (int k=0;k<N;k++) { double z=-L+k*dx;
                long idx=(long)i*NN+j*N+k;
                double r2=(x-ocx)*(x-ocx)+(y-ocy)*(y-ocy)+(z-ocz)*(z-ocz);
                double env=A*exp(-r2/(2*sigma*sigma));
                for (int a=0;a<NFIELDS;a++) arrays[a][idx]+=env*cos(delta[a]);
            }}}
            fprintf(stderr, "  Added oscillon at (%.1f,%.1f,%.1f) A=%.2f sigma=%.1f\n",
                    ocx,ocy,ocz,A,sigma);
        }
        else if (mod->type == 2) { /* perturb */
            double amp = mod->params[0];
            double pm = 0;
            for(long i=0;i<N3;i++) for(int a=0;a<3;a++)
                if(fabs(arrays[a][i])>pm) pm=fabs(arrays[a][i]);
            for(long i=0;i<N3;i++) for(int a=0;a<3;a++)
                arrays[a][i] += amp*pm*((double)rand()/RAND_MAX*2-1);
            fprintf(stderr, "  Perturbed by %.1f%% of phi_max (%.4f)\n", amp*100, pm);
        }
        else if (mod->type == 3) { /* scale */
            double factor = mod->params[0];
            for(long i=0;i<N3;i++) for(int f=0;f<12;f++) arrays[f][i]*=factor;
            fprintf(stderr, "  Scaled all fields by %.3f\n", factor);
        }
        else if (mod->type == 4) { /* flip-chirality */
            /* Reverse z-velocities to flip helical winding direction */
            for(long i=0;i<N3;i++) {
                arrays[8][i] = -arrays[8][i];   /* phi_vz */
                arrays[11][i] = -arrays[11][i]; /* theta_vz */
            }
            fprintf(stderr, "  Flipped chirality (reversed z-velocities)\n");
        }
    }

    /* Write output SFA (f64 for maximum fidelity in seed) */
    SFA *out = sfa_create(out_path, N, N, N, L, L, L, in->dt);

    /* Copy KVMD from input */
    SFA_KVMDSet kv[16];
    int n_kv = sfa_read_kvmd(in, kv, 16);
    for (int s=0;s<n_kv;s++) {
        const char *keys[128], *vals[128];
        for(int p=0;p<kv[s].n_pairs;p++){keys[p]=kv[s].keys[p];vals[p]=kv[s].values[p];}
        sfa_add_kvmd(out, kv[s].set_id, kv[s].first_frame, kv[s].frame_count,
                     keys, vals, kv[s].n_pairs);
    }

    sfa_add_column(out,"phi_x",SFA_F64,SFA_POSITION,0);
    sfa_add_column(out,"phi_y",SFA_F64,SFA_POSITION,1);
    sfa_add_column(out,"phi_z",SFA_F64,SFA_POSITION,2);
    sfa_add_column(out,"theta_x",SFA_F64,SFA_ANGLE,0);
    sfa_add_column(out,"theta_y",SFA_F64,SFA_ANGLE,1);
    sfa_add_column(out,"theta_z",SFA_F64,SFA_ANGLE,2);
    sfa_add_column(out,"phi_vx",SFA_F64,SFA_VELOCITY,0);
    sfa_add_column(out,"phi_vy",SFA_F64,SFA_VELOCITY,1);
    sfa_add_column(out,"phi_vz",SFA_F64,SFA_VELOCITY,2);
    sfa_add_column(out,"theta_vx",SFA_F64,SFA_VELOCITY,3);
    sfa_add_column(out,"theta_vy",SFA_F64,SFA_VELOCITY,4);
    sfa_add_column(out,"theta_vz",SFA_F64,SFA_VELOCITY,5);
    sfa_finalize_header(out);

    sfa_write_frame(out, 0.0, (void**)arrays);
    sfa_close(out);
    sfa_close(in);

    fprintf(stderr, "Wrote %s (N=%d, %d modifications applied)\n", out_path, N, n_mods);
    for (int i=0;i<12;i++) free(arrays[i]);
    return 0;
}
