/*  volview_sfa.c — Interactive SFA volume viewer
 *
 *  Reads .sfa archives and steps through frames with animation.
 *  Colors: RED=bound(|P|), GREEN=fabric(Σφ²), BLUE=angle(Σθ²)
 *
 *  Controls:
 *    Mouse drag: rotate, Scroll/Up/Down: zoom
 *    Left/Right arrows or A/D: prev/next frame
 *    Space: play/pause animation
 *    1/2/3: toggle R/G/B, +/-: brightness
 *    Home: first frame, End: last frame
 *    S: screenshot, R: reset view, Q: quit
 *
 *  Build: gcc -O3 -fopenmp -o volview_sfa viewer/volview_sfa.c \
 *         $(sdl2-config --cflags --libs) -lzstd -lm
 */

#define SFA_IMPLEMENTATION
#include "sfa.h"

#include <SDL2/SDL.h>
#include <math.h>

#define WIN_W 900
#define WIN_H 700
#define MAX_STEPS 200

/* Volume data */
static int volN;
static float *vol_r, *vol_g, *vol_b;
static SFA *archive;
static int cur_frame = 0;
static int playing = 0;

/* Camera */
static float cam_theta = 0.7f, cam_phi = 0.3f, cam_dist = 2.5f;
static int show_r = 1, show_g = 1, show_b = 1;
static float brightness = 3.0f;
static float opacity = 2.0f;   /* lower = more translucent (default was 8) */
static int need_render = 1;

static void load_frame(int idx) {
    if (idx < 0) idx = 0;
    if (idx >= (int)archive->total_frames) idx = archive->total_frames - 1;
    cur_frame = idx;

    long N3 = (long)archive->N_total;
    void *buf = malloc(archive->frame_bytes);
    if (sfa_read_frame(archive, idx, buf) < 0) {
        fprintf(stderr, "Failed to read frame %d\n", idx);
        free(buf);
        return;
    }

    double *fields[6];
    uint64_t off = 0;
    for (int c = 0; c < (int)archive->n_columns && c < 6; c++) {
        fields[c] = (double *)((uint8_t *)buf + off);
        off += N3 * sfa_dtype_size[archive->columns[c].dtype];
    }
    for (int c = archive->n_columns; c < 6; c++)
        fields[c] = NULL;

    /* Compute channels */
    float maxP = 0, maxPhi2 = 0, maxTheta2 = 0;
    for (long i = 0; i < N3; i++) {
        double p0 = fields[0][i], p1 = fields[1][i], p2 = fields[2][i];
        float absP = (float)fabs(p0*p1*p2);
        float phi2 = (float)(p0*p0 + p1*p1 + p2*p2);
        float theta2 = 0;
        if (fields[3]) theta2 = (float)(fields[3][i]*fields[3][i] +
            fields[4][i]*fields[4][i] + fields[5][i]*fields[5][i]);
        if (absP > maxP) maxP = absP;
        if (phi2 > maxPhi2) maxPhi2 = phi2;
        if (theta2 > maxTheta2) maxTheta2 = theta2;
    }

    float normP = maxP * 0.3f; if (normP < 1e-20f) normP = 1;
    float normPhi2 = maxPhi2 * 0.15f; if (normPhi2 < 1e-20f) normPhi2 = 1;
    float normTheta2 = maxTheta2 * 0.3f; if (normTheta2 < 1e-20f) normTheta2 = 1;

    for (long i = 0; i < N3; i++) {
        double p0 = fields[0][i], p1 = fields[1][i], p2 = fields[2][i];
        float absP = (float)fabs(p0*p1*p2);
        float phi2 = (float)(p0*p0 + p1*p1 + p2*p2);
        float pn = absP/normP; if (pn > 1) pn = 1;
        vol_r[i] = pn;
        float unbound = phi2 * (1.0f - sqrtf(pn));
        vol_g[i] = unbound/normPhi2; if (vol_g[i] > 1) vol_g[i] = 1;
        float theta2 = 0;
        if (fields[3]) theta2 = (float)(fields[3][i]*fields[3][i] +
            fields[4][i]*fields[4][i] + fields[5][i]*fields[5][i]);
        vol_b[i] = theta2/normTheta2; if (vol_b[i] > 1) vol_b[i] = 1;
    }

    free(buf);
    need_render = 1;
}

static void sample_vol(float x, float y, float z, float *r, float *g, float *b) {
    float n = (float)(volN - 1);
    float fx = x*n, fy = y*n, fz = z*n;
    int ix = (int)fx, iy = (int)fy, iz = (int)fz;
    if (ix<0||ix>=volN-1||iy<0||iy>=volN-1||iz<0||iz>=volN-1) {*r=*g=*b=0; return;}
    float dx=fx-ix, dy=fy-iy, dz=fz-iz;
    int nn=volN, n2=nn*nn;
    #define S(a,i,j,k) (a)[(i)*n2+(j)*nn+(k)]
    float w[8]={(1-dx)*(1-dy)*(1-dz),(1-dx)*(1-dy)*dz,
                (1-dx)*dy*(1-dz),(1-dx)*dy*dz,
                dx*(1-dy)*(1-dz),dx*(1-dy)*dz,
                dx*dy*(1-dz),dx*dy*dz};
    *r=*g=*b=0;
    int di[8]={0,0,0,0,1,1,1,1}, dj[8]={0,0,1,1,0,0,1,1}, dk[8]={0,1,0,1,0,1,0,1};
    for(int q=0;q<8;q++){
        int ii=ix+di[q],jj=iy+dj[q],kk=iz+dk[q];
        *r+=w[q]*S(vol_r,ii,jj,kk);
        *g+=w[q]*S(vol_g,ii,jj,kk);
        *b+=w[q]*S(vol_b,ii,jj,kk);
    }
    #undef S
}

static void render(Uint32 *pixels) {
    float cx=cam_dist*sinf(cam_theta)*cosf(cam_phi);
    float cy=cam_dist*sinf(cam_theta)*sinf(cam_phi);
    float cz=cam_dist*cosf(cam_theta);
    float fwd[3]={-cx,-cy,-cz};
    float len=sqrtf(fwd[0]*fwd[0]+fwd[1]*fwd[1]+fwd[2]*fwd[2]);
    fwd[0]/=len;fwd[1]/=len;fwd[2]/=len;
    float up[3]={0,0,1};
    float right[3]={fwd[1]*up[2]-fwd[2]*up[1],fwd[2]*up[0]-fwd[0]*up[2],fwd[0]*up[1]-fwd[1]*up[0]};
    len=sqrtf(right[0]*right[0]+right[1]*right[1]+right[2]*right[2]);
    if(len<1e-6){right[0]=1;right[1]=0;right[2]=0;len=1;}
    right[0]/=len;right[1]/=len;right[2]/=len;
    float cup[3]={right[1]*fwd[2]-right[2]*fwd[1],right[2]*fwd[0]-right[0]*fwd[2],right[0]*fwd[1]-right[1]*fwd[0]};
    float step=1.4f/MAX_STEPS, fov=0.8f;

    #pragma omp parallel for schedule(dynamic,4)
    for(int py=0;py<WIN_H;py++){
        for(int px=0;px<WIN_W;px++){
            float u=(2.0f*px/WIN_W-1)*fov*WIN_W/WIN_H;
            float v=(1-2.0f*py/WIN_H)*fov;
            float dir[3]={fwd[0]+u*right[0]+v*cup[0],fwd[1]+u*right[1]+v*cup[1],fwd[2]+u*right[2]+v*cup[2]};
            len=sqrtf(dir[0]*dir[0]+dir[1]*dir[1]+dir[2]*dir[2]);
            dir[0]/=len;dir[1]/=len;dir[2]/=len;
            float ox=cx+.5f,oy=cy+.5f,oz=cz+.5f;
            float tmin=0,tmax=100;
            for(int d=0;d<3;d++){
                float o=(d==0)?ox:(d==1)?oy:oz;
                float di=(d==0)?dir[0]:(d==1)?dir[1]:dir[2];
                if(fabsf(di)<1e-8)continue;
                float t0=-o/di,t1=(1-o)/di;
                if(t0>t1){float tmp=t0;t0=t1;t1=tmp;}
                if(t0>tmin)tmin=t0;if(t1<tmax)tmax=t1;
            }
            float ar=0,ag=0,ab=0,aa=0;
            if(tmin<tmax&&tmax>0){
                if(tmin<0)tmin=0;
                for(float tt=tmin;tt<tmax&&aa<.98f;tt+=step){
                    float sr,sg,sb;
                    sample_vol(ox+tt*dir[0],oy+tt*dir[1],oz+tt*dir[2],&sr,&sg,&sb);
                    sr*=show_r;sg*=show_g;sb*=show_b;
                    float dens=(sr+sg+sb)*brightness;
                    float alpha=1-expf(-dens*step*opacity);
                    if(alpha<.001f)continue;
                    sr=powf(sr,.6f)*brightness;sg=powf(sg,.6f)*brightness;sb=powf(sb,.6f)*brightness;
                    ar+=(1-aa)*sr*alpha;ag+=(1-aa)*sg*alpha;ab+=(1-aa)*sb*alpha;aa+=(1-aa)*alpha;
                }
            }
            ar+=(1-aa)*.05f;ag+=(1-aa)*.05f;ab+=(1-aa)*.07f;
            int ir=(int)(fminf(ar,1)*255),ig=(int)(fminf(ag,1)*255),ib=(int)(fminf(ab,1)*255);
            pixels[py*WIN_W+px]=(ir<<16)|(ig<<8)|ib;
        }
    }
}

int main(int argc, char **argv) {
    if (argc < 2) {
        printf("Usage: %s <archive.sfa>\n", argv[0]);
        printf("Controls: drag=rotate, scroll=zoom, Left/Right=frame, Space=play\n");
        printf("  1/2/3=toggle R/G/B, +/-=brightness, S=screenshot, Q=quit\n");
        return 1;
    }

    archive = sfa_open(argv[1]);
    if (!archive) { fprintf(stderr, "Failed to open %s\n", argv[1]); return 1; }

    printf("SFA: %dx%dx%d, %d columns, %d frames\n",
           archive->Nx, archive->Ny, archive->Nz,
           archive->n_columns, archive->total_frames);
    for (uint32_t c = 0; c < archive->n_columns; c++)
        printf("  %s (%s)\n", archive->columns[c].name,
               sfa_dtype_name(archive->columns[c].dtype));

    volN = archive->Nx;  /* assumes cubic */
    long N3 = (long)archive->N_total;
    vol_r = malloc(N3 * sizeof(float));
    vol_g = malloc(N3 * sizeof(float));
    vol_b = malloc(N3 * sizeof(float));

    load_frame(0);

    SDL_Init(SDL_INIT_VIDEO);
    char title[256];
    snprintf(title, sizeof(title), "SFA Viewer — %s — frame 0/%d",
             argv[1], archive->total_frames);
    SDL_Window *win = SDL_CreateWindow(title, SDL_WINDOWPOS_CENTERED,
        SDL_WINDOWPOS_CENTERED, WIN_W, WIN_H, 0);
    SDL_Renderer *ren = SDL_CreateRenderer(win, -1, SDL_RENDERER_ACCELERATED);
    SDL_Texture *tex = SDL_CreateTexture(ren, SDL_PIXELFORMAT_RGB888,
        SDL_TEXTUREACCESS_STREAMING, WIN_W, WIN_H);
    Uint32 *pixels = malloc(WIN_W * WIN_H * sizeof(Uint32));

    int running = 1, dragging = 0;
    Uint32 last_anim = 0;
    int anim_interval = 200;  /* ms between frames when playing */

    while (running) {
        SDL_Event ev;
        while (SDL_PollEvent(&ev)) {
            switch (ev.type) {
            case SDL_QUIT: running=0; break;
            case SDL_KEYDOWN:
                switch (ev.key.keysym.sym) {
                case SDLK_q: case SDLK_ESCAPE: running=0; break;
                case SDLK_LEFT: case SDLK_a: load_frame(cur_frame-1); break;
                case SDLK_RIGHT: case SDLK_d: load_frame(cur_frame+1); break;
                case SDLK_HOME: load_frame(0); break;
                case SDLK_END: load_frame(archive->total_frames-1); break;
                case SDLK_SPACE: playing=!playing; break;
                case SDLK_1: show_r=!show_r; need_render=1; break;
                case SDLK_2: show_g=!show_g; need_render=1; break;
                case SDLK_3: show_b=!show_b; need_render=1; break;
                case SDLK_EQUALS: brightness*=1.3f; need_render=1; break;
                case SDLK_MINUS: brightness/=1.3f; need_render=1; break;
                case SDLK_UP: cam_dist*=0.9f; need_render=1; break;
                case SDLK_DOWN: cam_dist*=1.1f; need_render=1; break;
                case SDLK_r: cam_theta=.7f;cam_phi=.3f;cam_dist=2.5f;brightness=3;opacity=2;need_render=1; break;
                case SDLK_o: opacity*=0.7f; if(opacity<0.1f)opacity=0.1f; need_render=1; break;
                case SDLK_p: opacity*=1.4f; if(opacity>20)opacity=20; need_render=1; break;
                case SDLK_s: {
                    FILE *f=fopen("screenshot.ppm","wb");
                    fprintf(f,"P6\n%d %d\n255\n",WIN_W,WIN_H);
                    for(int i=0;i<WIN_W*WIN_H;i++){
                        unsigned char rgb[3]={(pixels[i]>>16)&0xff,(pixels[i]>>8)&0xff,pixels[i]&0xff};
                        fwrite(rgb,1,3,f);
                    }
                    fclose(f); printf("Saved screenshot.ppm\n");
                } break;
                }
                break;
            case SDL_MOUSEBUTTONDOWN:
                if(ev.button.button==SDL_BUTTON_LEFT)dragging=1; break;
            case SDL_MOUSEBUTTONUP:
                if(ev.button.button==SDL_BUTTON_LEFT)dragging=0; break;
            case SDL_MOUSEMOTION:
                if(dragging){
                    cam_phi+=ev.motion.xrel*.01f;
                    cam_theta-=ev.motion.yrel*.01f;
                    if(cam_theta<.1f)cam_theta=.1f;
                    if(cam_theta>3.04f)cam_theta=3.04f;
                    need_render=1;
                }
                break;
            case SDL_MOUSEWHEEL:
                cam_dist*=(ev.wheel.y>0)?.9f:1.1f;
                if(cam_dist<.5f)cam_dist=.5f;
                if(cam_dist>10)cam_dist=10;
                need_render=1;
                break;
            }
        }

        /* Animation */
        if (playing) {
            Uint32 now = SDL_GetTicks();
            if (now - last_anim > (Uint32)anim_interval) {
                int next = cur_frame + 1;
                if (next >= (int)archive->total_frames) next = 0;
                load_frame(next);
                last_anim = now;
            }
        }

        if (need_render) {
            render(pixels);
            SDL_UpdateTexture(tex, NULL, pixels, WIN_W*sizeof(Uint32));
            SDL_RenderCopy(ren, tex, NULL, NULL);
            /* HUD */
            double t = sfa_frame_time(archive, cur_frame);
            snprintf(title, sizeof(title), "SFA Viewer — frame %d/%d t=%.1f — R%d G%d B%d opac=%.1f %s",
                     cur_frame, archive->total_frames, t,
                     show_r, show_g, show_b, opacity, playing?"[PLAY]":"");
            SDL_SetWindowTitle(win, title);
            SDL_RenderPresent(ren);
            need_render = 0;
        }
        SDL_Delay(16);
    }

    free(pixels); free(vol_r); free(vol_g); free(vol_b);
    SDL_DestroyTexture(tex); SDL_DestroyRenderer(ren); SDL_DestroyWindow(win);
    SDL_Quit();
    sfa_close(archive);
    return 0;
}
