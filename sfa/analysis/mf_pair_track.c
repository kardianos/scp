/* mf_pair_track.c — track C vs L centroids in multi-fabric 54-col SFA
 *
 * Expects columns named phi_x/y/z (C) and lphi_x/y/z (L), optional im.
 * Writes TSV: frame t nC nL massC massL Qc Ql cxC cyC czC cxL cyL czL D
 *
 * Build: gcc -O3 -o bin/mf_pair_track sfa/analysis/mf_pair_track.c -lzstd -lm
 * Usage: mf_pair_track file.sfa [out.tsv]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define SFA_IMPLEMENTATION
#include "../format/sfa.h"

static int find_col(SFA *s, const char *name) {
    for (uint32_t c = 0; c < s->n_columns; c++)
        if (!strcmp(s->columns[c].name, name)) return (int)c;
    return -1;
}

static void *col_ptr(void *frame, SFA *s, int col) {
    if (col < 0) return NULL;
    uint64_t off = 0;
    for (int c = 0; c < col; c++)
        off += s->N_total * sfa_dtype_size[s->columns[c].dtype];
    return (uint8_t*)frame + off;
}

static double getf(void *base, int dtype, long i) {
    if (!base) return 0;
    if (dtype == SFA_F64) return ((double*)base)[i];
    if (dtype == SFA_F32) return (double)((float*)base)[i];
    if (dtype == SFA_F16) return (double)sfa_f16_to_f32(((const uint16_t*)base)[i]);
    return 0;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s file.sfa [out.tsv]\n", argv[0]);
        return 1;
    }
    const char *outpath = argc > 2 ? argv[2] : "mf_track.tsv";
    SFA *s = sfa_open(argv[1]);
    if (!s) { fprintf(stderr, "open failed\n"); return 1; }

    int c_px = find_col(s, "phi_x"), c_py = find_col(s, "phi_y"), c_pz = find_col(s, "phi_z");
    int c_ix = find_col(s, "phiim_x"), c_iy = find_col(s, "phiim_y"), c_iz = find_col(s, "phiim_z");
    int c_vx = find_col(s, "phi_vx"), c_vy = find_col(s, "phi_vy"), c_vz = find_col(s, "phi_vz");
    int c_ivx = find_col(s, "phiim_vx"), c_ivy = find_col(s, "phiim_vy"), c_ivz = find_col(s, "phiim_vz");
    int l_px = find_col(s, "lphi_x"), l_py = find_col(s, "lphi_y"), l_pz = find_col(s, "lphi_z");
    int l_ix = find_col(s, "lphiim_x"), l_iy = find_col(s, "lphiim_y"), l_iz = find_col(s, "lphiim_z");
    int l_vx = find_col(s, "lphi_vx"), l_vy = find_col(s, "lphi_vy"), l_vz = find_col(s, "lphi_vz");
    int l_ivx = find_col(s, "lphiim_vx"), l_ivy = find_col(s, "lphiim_vy"), l_ivz = find_col(s, "lphiim_vz");
    if (c_px < 0 || l_px < 0) {
        fprintf(stderr, "FATAL: need phi_* and lphi_* columns (got n_col=%u)\n", s->n_columns);
        return 1;
    }
    int N = (int)s->Nx;
    double L = s->Lx, dx = 2.0 * L / (N - 1), dV = dx*dx*dx;
    long N3 = s->N_total;
    void *frame = malloc(s->frame_bytes);
    FILE *out = fopen(outpath, "w");
    fprintf(out, "frame\tt\tmassC\tmassL\tQc\tQl\tcxC\tcyC\tczC\tcxL\tcyL\tczL\tD\n");

    for (uint32_t f = 0; f < s->total_frames; f++) {
        sfa_read_frame(s, f, frame);
        double t = 0; /* time from frame index * dt if available */
        if (s->dt > 0) t = f * s->dt; /* may not match snap_dt; ok for ordering */
        /* better: use embedded time if API provides - use frame index via snap */
        void *Cpx=col_ptr(frame,s,c_px), *Cpy=col_ptr(frame,s,c_py), *Cpz=col_ptr(frame,s,c_pz);
        void *Cix=col_ptr(frame,s,c_ix), *Ciy=col_ptr(frame,s,c_iy), *Ciz=col_ptr(frame,s,c_iz);
        void *Cvx=col_ptr(frame,s,c_vx), *Cvy=col_ptr(frame,s,c_vy), *Cvz=col_ptr(frame,s,c_vz);
        void *Civx=col_ptr(frame,s,c_ivx), *Civy=col_ptr(frame,s,c_ivy), *Civz=col_ptr(frame,s,c_ivz);
        void *Lpx=col_ptr(frame,s,l_px), *Lpy=col_ptr(frame,s,l_py), *Lpz=col_ptr(frame,s,l_pz);
        void *Lix=col_ptr(frame,s,l_ix), *Liy=col_ptr(frame,s,l_iy), *Liz=col_ptr(frame,s,l_iz);
        void *Lvx=col_ptr(frame,s,l_vx), *Lvy=col_ptr(frame,s,l_vy), *Lvz=col_ptr(frame,s,l_vz);
        void *Livx=col_ptr(frame,s,l_ivx), *Livy=col_ptr(frame,s,l_ivy), *Livz=col_ptr(frame,s,l_ivz);
        int dC = s->columns[c_px].dtype, dL = s->columns[l_px].dtype;

        double mC=0,mL=0,Qc=0,Ql=0,cxC=0,cyC=0,czC=0,cxL=0,cyL=0,czL=0;
        for (long idx=0; idx<N3; idx++) {
            int i=(int)(idx/(N*N)), j=(int)((idx/N)%N), k=(int)(idx%N);
            double x=-L+i*dx, y=-L+j*dx, z=-L+k*dx;
            double ux=getf(Cpx,dC,idx), uy=getf(Cpy,dC,idx), uz=getf(Cpz,dC,idx);
            double vx=getf(Cix,dC,idx), vy=getf(Ciy,dC,idx), vz=getf(Ciz,dC,idx);
            double rho2C = ux*ux+uy*uy+uz*uz+vx*vx+vy*vy+vz*vz;
            double udx=getf(Cvx,dC,idx), udy=getf(Cvy,dC,idx), udz=getf(Cvz,dC,idx);
            double vdx=getf(Civx,dC,idx), vdy=getf(Civy,dC,idx), vdz=getf(Civz,dC,idx);
            double qlocC = ux*vdx-vx*udx + uy*vdy-vy*udy + uz*vdz-vz*udz;
            mC += rho2C*dV; Qc += qlocC*dV;
            cxC += x*rho2C; cyC += y*rho2C; czC += z*rho2C;

            double Lx=getf(Lpx,dL,idx), Ly=getf(Lpy,dL,idx), Lz=getf(Lpz,dL,idx);
            double Lixv=getf(Lix,dL,idx), Liyv=getf(Liy,dL,idx), Lizv=getf(Liz,dL,idx);
            double rho2L = Lx*Lx+Ly*Ly+Lz*Lz+Lixv*Lixv+Liyv*Liyv+Lizv*Lizv;
            double Ludx=getf(Lvx,dL,idx), Ludy=getf(Lvy,dL,idx), Ludz=getf(Lvz,dL,idx);
            double Lvdx=getf(Livx,dL,idx), Lvdy=getf(Livy,dL,idx), Lvdz=getf(Livz,dL,idx);
            double qlocL = Lx*Lvdx-Lixv*Ludx + Ly*Lvdy-Liyv*Ludy + Lz*Lvdz-Lizv*Ludz;
            mL += rho2L*dV; Ql += qlocL*dV;
            cxL += x*rho2L; cyL += y*rho2L; czL += z*rho2L;
        }
        if (mC > 1e-30) { cxC*=dV/mC; cyC*=dV/mC; czC*=dV/mC; } else { cxC=cyC=czC=0; }
        if (mL > 1e-30) { cxL*=dV/mL; cyL*=dV/mL; czL*=dV/mL; } else { cxL=cyL=czL=0; }
        /* fix centroid: rho2 sum without dV in cx */
        /* actually cxC = sum x*rho2 / sum rho2 — I multiplied dV inconsistently */
        /* recompute: mass = sum rho2*dV, cx = sum x*rho2*dV / mass = sum x*rho2 / sum rho2 */
        double sC=0,sL=0, sxC=0,syC=0,szC=0,sxL=0,syL=0,szL=0;
        for (long idx=0; idx<N3; idx++) {
            int i=(int)(idx/(N*N)), j=(int)((idx/N)%N), k=(int)(idx%N);
            double x=-L+i*dx, y=-L+j*dx, z=-L+k*dx;
            double ux=getf(Cpx,dC,idx), uy=getf(Cpy,dC,idx), uz=getf(Cpz,dC,idx);
            double vx=getf(Cix,dC,idx), vy=getf(Ciy,dC,idx), vz=getf(Ciz,dC,idx);
            double r2=ux*ux+uy*uy+uz*uz+vx*vx+vy*vy+vz*vz;
            sC+=r2; sxC+=x*r2; syC+=y*r2; szC+=z*r2;
            double Lx=getf(Lpx,dL,idx), Ly=getf(Lpy,dL,idx), Lz=getf(Lpz,dL,idx);
            double Lixv=getf(Lix,dL,idx), Liyv=getf(Liy,dL,idx), Lizv=getf(Liz,dL,idx);
            double r2L=Lx*Lx+Ly*Ly+Lz*Lz+Lixv*Lixv+Liyv*Liyv+Lizv*Lizv;
            sL+=r2L; sxL+=x*r2L; syL+=y*r2L; szL+=z*r2L;
        }
        if (sC>1e-30) { cxC=sxC/sC; cyC=syC/sC; czC=szC/sC; }
        if (sL>1e-30) { cxL=sxL/sL; cyL=syL/sL; czL=szL/sL; }
        double D = sqrt((cxC-cxL)*(cxC-cxL)+(cyC-cyL)*(cyC-cyL)+(czC-czL)*(czC-czL));
        /* wall time from snap: use frame * unknown snap - print frame */
        fprintf(out, "%u\t%.6f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",
                f, (double)f, mC, mL, Qc, Ql, cxC, cyC, czC, cxL, cyL, czL, D);
        printf("frame %u D=%.4f massC=%.2f massL=%.2f Qc=%.2f Ql=%.2f\n", f, D, mC, mL, Qc, Ql);
    }
    fclose(out); free(frame); sfa_close(s);
    printf("Wrote %s\n", outpath);
    return 0;
}
