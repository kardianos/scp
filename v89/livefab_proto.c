/*
 * livefab_proto.c — livefab (S3) matter-insertion test (GLM, 2026-08-01).
 *
 * Standalone. Does NOT touch cellfab.c. Subordinate to LIVEFAB.md.
 *
 * Loads the S1 annealed substrate (sigma_d ~3%, the cleanest frozen
 * substrate the program has) and runs THE discriminator: insert matter
 * (convert space -> shrink radii in a central cluster) and compare a
 * FROZEN-link substrate against a LIVE-link (livefab) substrate.
 *
 *   frozen (S1):   link set fixed at the vacuum geometry. The shrunk cells'
 *                  old links go slack; surrounding cells cannot fill the
 *                  freed volume; sigma_d measured over frozen links grows.
 *   live  (S3):    links re-derive each sweep from current radii
 *                  (d < 1.15*(ri+rj)). Shrunk cells lose links; neighbours
 *                  re-jam into the freed volume; sigma_d over LIVE links
 *                  stays low. The contact graph reorganizes to account
 *                  for the matter — the capability a frozen substrate
 *                  cannot have, and the operational content of
 *                  PRINCIPLE §4.2-4.3 ("where there is matter, there is
 *                  less space; the remaining space must account for it").
 *
 * Build:  gcc -O2 -o livefab_proto livefab_proto.c -lm
 * Run:    ./livefab_proto [anneal.tsv] [matter_frac] [sweeps]
 *
 * Input: the cells_000000.tsv produced by cellfab's S1 anneal
 *        (cols: id x y z r Es Em Ee ...). Default path below.
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

typedef struct { double sigma_d, dbar; int nlinks; double deg; } Stats;

static int hg_h(int gx,int gy,int gz,int g){unsigned h=(unsigned)(gx*73856093)^(unsigned)(gy*19349663)^(unsigned)(gz*83492791);return h%(unsigned)(g*g*g);}

/* links where d < cfac*(ri+rj). */
static int build_links(int NC,const double*x,const double*y,const double*z,
                       const double*rad,double L,double cfac,
                       int**lio,int**ljo,int cap){
    double rmax=0; for(int i=0;i<NC;i++) if(rad[i]>rmax)rmax=rad[i];
    double rcut=cfac*2.0*rmax;
    int g=(int)(L/rcut); if(g<1)g=1; if(g>80)g=80;
    int nb=g*g*g; int*head=malloc(nb*sizeof(int)); int*cell=malloc(NC*sizeof(int));
    for(int b=0;b<nb;b++)head[b]=-1;
    for(int i=0;i<NC;i++){
        int gx=(int)(x[i]/L*g),gy=(int)(y[i]/L*g),gz=(int)(z[i]/L*g);
        if(gx>=g)gx=g-1; if(gy>=g)gy=g-1; if(gz>=g)gz=g-1; if(gx<0)gx=0; if(gy<0)gy=0; if(gz<0)gz=0;
        int b=hg_h(gx,gy,gz,g); cell[i]=head[b]; head[b]=i;
    }
    int*li=*lio,*lj=*ljo,ne=0;
    for(int i=0;i<NC;i++){
        int gx=(int)(x[i]/L*g),gy=(int)(y[i]/L*g),gz=(int)(z[i]/L*g);
        if(gx>=g)gx=g-1; if(gy>=g)gy=g-1; if(gz>=g)gz=g-1; if(gx<0)gx=0; if(gy<0)gy=0; if(gz<0)gz=0;
        for(int ax=gx-1;ax<=gx+1;ax++)
        for(int ay=gy-1;ay<=gy+1;ay++)
        for(int az=gz-1;az<=gz+1;az++){
            if(ax<0||ay<0||az<0||ax>=g||ay>=g||az>=g)continue;
            for(int q=head[hg_h(ax,ay,az,g)];q>=0;q=cell[q]){
                if(q<=i)continue;
                double dx=x[q]-x[i],dy=y[q]-y[i],dz=z[q]-z[i];
                double d2=dx*dx+dy*dy+dz*dz;
                double cut=cfac*(rad[i]+rad[q]);
                if(d2>=cut*cut||d2<=0)continue;
                if(ne>=cap){cap*=2;li=*lio=realloc(li,cap*sizeof(int));lj=*ljo=realloc(lj,cap*sizeof(int));}
                li[ne]=i;lj[ne]=q;ne++;
            }
        }
    }
    free(head);free(cell);
    return ne;
}

static Stats stats_region(int NC,const double*x,const double*y,const double*z,
                          const int*li,const int*lj,int ne,
                          double cx,double cy,double cz,double rshell,int region){
    Stats s={0,0,0,0}; double sm=0,sm2=0; int n=0;
    for(int e=0;e<ne;e++){
        int i=li[e],j=lj[e];
        if(region>=0){
            double mx=0.5*(x[i]+x[j])-cx,my=0.5*(y[i]+y[j])-cy,mz=0.5*(z[i]+z[j])-cz;
            double rr=sqrt(mx*mx+my*my+mz*mz);
            if(region==0 && rr>=rshell)continue;
            if(region==1 && rr<rshell*2.0)continue;
        }
        double dx=x[i]-x[j],dy=y[i]-y[j],dz=z[i]-z[j];
        double d=sqrt(dx*dx+dy*dy+dz*dz);
        sm+=d;sm2+=d*d;n++;
    }
    if(n<2)return s;
    double m=sm/n,var=sm2/n-m*m; if(var<0)var=0;
    s.dbar=m; s.sigma_d=sqrt(var)/m; s.nlinks=n; s.deg=2.0*ne/NC;
    return s;
}

/* pure repulsion below (ri+rj) over given edge set */
static void relax(int NC,double*x,double*y,double*z,const double*rad,
                  const int*li,const int*lj,int ne,double k){
    double*fx=calloc(NC,sizeof(double)),*fy=calloc(NC,sizeof(double)),*fz=calloc(NC,sizeof(double));
    for(int e=0;e<ne;e++){
        int i=li[e],j=lj[e];
        double dx=x[j]-x[i],dy=y[j]-y[i],dz=z[j]-z[i];
        double d=sqrt(dx*dx+dy*dy+dz*dz); if(d<=0)continue;
        double contact=rad[i]+rad[j];
        if(d>=contact)continue;
        double f=k*(contact-d)/d;
        fx[i]-=f*dx;fy[i]-=f*dy;fz[i]-=f*dz; fx[j]+=f*dx;fy[j]+=f*dy;fz[j]+=f*dz;
    }
    for(int i=0;i<NC;i++){x[i]+=fx[i];y[i]+=fy[i];z[i]+=fz[i];}
    free(fx);free(fy);free(fz);
}

static int load_tsv(const char*path,double**xo,double**yo,double**zo,double**ro,double*Lout){
    FILE*f=fopen(path,"r"); if(!f){perror(path);return -1;}
    char line[1024]; int NC=0, header=1;
    /* first pass: count, infer box */
    double xmin=1e30,xmax=-1e30,ymin=1e30,ymax=-1e30,zmin=1e30,zmax=-1e30;
    while(fgets(line,sizeof line,f)){
        if(line[0]=='#')continue;
        if(header){ /* skip header line (non-numeric starts) */
            char c; double v;
            if(sscanf(line,"%lf",&v)!=1){header=0;continue;}
            header=0;
        }
        double xx,yy,zz,rr; int id;
        if(sscanf(line,"%d %lf %lf %lf %lf",&id,&xx,&yy,&zz,&rr)>=4){
            NC++;
            if(xx<xmin)xmin=xx; if(xx>xmax)xmax=xx;
            if(yy<ymin)ymin=yy; if(yy>ymax)ymax=yy;
            if(zz<zmin)zmin=zz; if(zz>zmax)zmax=zz;
        }
    }
    double L=xmax; if(ymax>L)L=ymax; if(zmax>L)L=zmax;
    *Lout=L;
    *xo=malloc(NC*sizeof(double));*yo=malloc(NC*sizeof(double));
    *zo=malloc(NC*sizeof(double));*ro=malloc(NC*sizeof(double));
    rewind(f); NC=0; header=1;
    while(fgets(line,sizeof line,f)){
        if(line[0]=='#')continue;
        if(header){ double v; if(sscanf(line,"%lf",&v)!=1){header=0;continue;} header=0; }
        double xx,yy,zz,rr; int id;
        if(sscanf(line,"%d %lf %lf %lf %lf",&id,&xx,&yy,&zz,&rr)>=4){
            (*xo)[NC]=xx;(*yo)[NC]=yy;(*zo)[NC]=zz;(*ro)[NC]=rr; NC++;
        }
    }
    fclose(f);
    return NC;
}

int main(int argc,char**argv){
    const char*path=argc>1?argv[1]:"prestress/foam/anneal_snaps/cells_000000.tsv";
    double matter_frac=argc>2?atof(argv[2]):0.30;   /* r multiplier in matter region */
    int SW=argc>3?atoi(argv[3]):300;
    double cfac=1.15;                                 /* cellfab link rule */

    double L=0;
    double*x,*y,*z,*r; int NC=load_tsv(path,&x,&y,&z,&r,&L);
    if(NC<=0){fprintf(stderr,"no cells loaded from %s\n",path);return 1;}
    printf("# livefab_proto  substrate=%s NC=%d L=%.2f cfac=%.2f matter_r_mult=%.2f\n",
           path,NC,L,cfac,matter_frac);

    int cap=NC*30,*li=malloc(cap*sizeof(int)),*lj=malloc(cap*sizeof(int));

    /* baseline over the S1 substrate (vacuum) */
    int ne0=build_links(NC,x,y,z,r,L,cfac,&li,&lj,cap);
    Stats base_all=stats_region(NC,x,y,z,li,lj,ne0,0,0,0,0,-1);
    printf("# baseline (S1 vacuum): nlinks=%d deg=%.2f dbar=%.4f sigma_d=%.2f%%\n",
           ne0,base_all.deg,base_all.dbar,100*base_all.sigma_d);

    /* snapshot the VACUUM link set — the honest "frozen substrate" (S1)
     * model: candidate channels fixed at the vacuum geometry, THEN matter
     * is inserted and dynamics proceeds over the frozen candidate set. */
    int*li_vac=malloc(sizeof(int)*cap),*lj_vac=malloc(sizeof(int)*cap);
    int ne_vac=ne0;
    for(int e=0;e<ne_vac;e++){li_vac[e]=li[e];lj_vac[e]=lj[e];}

    double cx=L/2,cy=L/2,cz=L/2, Rm=0.16*L;   /* matter cluster radius */
    /* count cells in matter region */
    int nmat=0; for(int i=0;i<NC;i++){double dx=x[i]-cx,dy=y[i]-cy,dz=z[i]-cz; if(dx*dx+dy*dy+dz*dz<Rm*Rm)nmat++;}
    printf("# matter cluster: r=%.2f cells=%d (%.1f%%)  r_shrink %.2f->%.2f\n",
           Rm,nmat,100.0*nmat/NC,r[0],matter_frac*r[0]);

    /* helper: run one arm from a fresh copy of positions, with given link mode */
    #define ARM(title, LIVE) do {                                            \
        double*ax=malloc(NC*sizeof(double)),*ay=malloc(NC*sizeof(double)),   \
               *az=malloc(NC*sizeof(double)),*ar=malloc(NC*sizeof(double));  \
        for(int i=0;i<NC;i++){ax[i]=x[i];ay[i]=y[i];az[i]=z[i];ar[i]=r[i];}  \
        /* apply matter: shrink radii in core (space converted to pattern) */\
        for(int i=0;i<NC;i++){double dx=ax[i]-cx,dy=ay[i]-cy,dz=az[i]-cz;    \
            if(dx*dx+dy*dy+dz*dz<Rm*Rm) ar[i]=ar[i]*matter_frac;}            \
        /* FROZEN uses the vacuum graph (frozen candidate set); LIVE derives */\
        int ne_fixed = ne_vac;                                               \
        const int*li_use = (LIVE) ? NULL : li_vac;                           \
        const int*lj_use = (LIVE) ? NULL : lj_vac;                           \
        printf("\n== %s (links %s) ==\n",title,(LIVE)?"live":"frozen(vac)"); \
        for(int sw=0;sw<=SW;sw++){                                           \
            int ne; const int*li_r; const int*lj_r;                          \
            if(LIVE){ ne=build_links(NC,ax,ay,az,ar,L,cfac,&li,&lj,cap); li_r=li; lj_r=lj; } \
            else     { ne=ne_fixed; li_r=li_use; lj_r=lj_use; }              \
            Stats sa=stats_region(NC,ax,ay,az,li_r,lj_r,ne,cx,cy,cz,Rm,1);   \
            Stats sc=stats_region(NC,ax,ay,az,li_r,lj_r,ne,cx,cy,cz,Rm,0);   \
            if(sw%30==0||sw==SW)                                             \
                printf("  %-7s sw=%3d deg=%5.2f | core: n=%5d dbar=%.3f sd=%5.2f%% | " \
                       "far: n=%5d dbar=%.3f sd=%5.2f%%\n",                  \
                       title,sw,sc.deg, sc.nlinks,sc.dbar,100*sc.sigma_d,    \
                       sa.nlinks,sa.dbar,100*sa.sigma_d);                    \
            if(sw==SW)break;                                                 \
            relax(NC,ax,ay,az,ar,li,lj,ne,0.3);                              \
        }                                                                    \
        free(ax);free(ay);free(az);free(ar);                                 \
    } while(0)

    ARM("FROZEN", 0);
    ARM("LIVE",   1);

    free(x);free(y);free(z);free(r);free(li);free(lj);
    return 0;
}
