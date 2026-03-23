/*  volview.c — Interactive SFA volume viewer
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
#include "../format/sfa.h"

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
static int show_wireframe = 0;
static float *depth_buf = NULL;  /* per-pixel depth where volume becomes opaque */

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

    /* Read field values into f64 arrays, handling both f32 and f64 input */
    double *fld[6] = {NULL};
    for (int c = 0; c < 6; c++) fld[c] = (double*)calloc(N3, sizeof(double));

    uint64_t off = 0;
    for (int c = 0; c < (int)archive->n_columns && c < 6; c++) {
        int dtype = archive->columns[c].dtype;
        int es = sfa_dtype_size[dtype];
        uint8_t *src = (uint8_t*)buf + off;
        if (dtype == SFA_F64) {
            memcpy(fld[c], src, N3 * sizeof(double));
        } else if (dtype == SFA_F32) {
            float *fsrc = (float*)src;
            for (long i = 0; i < N3; i++) fld[c][i] = (double)fsrc[i];
        } else {
            /* Unsupported dtype — leave as zero */
        }
        off += N3 * es;
    }

    /* Compute channels */
    float maxP = 0, maxPhi2 = 0, maxTheta2 = 0;
    for (long i = 0; i < N3; i++) {
        double p0 = fld[0][i], p1 = fld[1][i], p2 = fld[2][i];
        float absP = (float)fabs(p0*p1*p2);
        float phi2 = (float)(p0*p0 + p1*p1 + p2*p2);
        float theta2 = 0;
        if (archive->n_columns >= 6)
            theta2 = (float)(fld[3][i]*fld[3][i] + fld[4][i]*fld[4][i] + fld[5][i]*fld[5][i]);
        if (absP > maxP) maxP = absP;
        if (phi2 > maxPhi2) maxPhi2 = phi2;
        if (theta2 > maxTheta2) maxTheta2 = theta2;
    }

    float normP = maxP * 0.3f; if (normP < 1e-20f) normP = 1;
    float normPhi2 = maxPhi2 * 0.15f; if (normPhi2 < 1e-20f) normPhi2 = 1;
    float normTheta2 = maxTheta2 * 0.3f; if (normTheta2 < 1e-20f) normTheta2 = 1;

    for (long i = 0; i < N3; i++) {
        double p0 = fld[0][i], p1 = fld[1][i], p2 = fld[2][i];
        float absP = (float)fabs(p0*p1*p2);
        float phi2 = (float)(p0*p0 + p1*p1 + p2*p2);
        float pn = absP/normP; if (pn > 1) pn = 1;
        vol_r[i] = pn;
        float unbound = phi2 * (1.0f - sqrtf(pn));
        vol_g[i] = unbound/normPhi2; if (vol_g[i] > 1) vol_g[i] = 1;
        float theta2 = 0;
        if (archive->n_columns >= 6)
            theta2 = (float)(fld[3][i]*fld[3][i] + fld[4][i]*fld[4][i] + fld[5][i]*fld[5][i]);
        vol_b[i] = theta2/normTheta2; if (vol_b[i] > 1) vol_b[i] = 1;
    }

    for (int c = 0; c < 6; c++) free(fld[c]);

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
            float vol_depth = tmax;  /* depth where volume becomes opaque */
            int depth_recorded = 0;
            if(tmin<tmax&&tmax>0){
                if(tmin<0)tmin=0;
                for(float tt=tmin;tt<tmax&&aa<.98f;tt+=step){
                    float sr,sg,sb;
                    sample_vol(ox+tt*dir[0],oy+tt*dir[1],oz+tt*dir[2],&sr,&sg,&sb);
                    sr*=show_r;sg*=show_g;sb*=show_b;
                    float dens=(sr+sg+sb)*brightness;
                    float alpha=1-expf(-dens*step*opacity);
                    if(alpha<.001f)continue;
                    if(!depth_recorded && aa > 0.15f) { vol_depth = tt; depth_recorded = 1; }
                    sr=powf(sr,.6f)*brightness;sg=powf(sg,.6f)*brightness;sb=powf(sb,.6f)*brightness;
                    ar+=(1-aa)*sr*alpha;ag+=(1-aa)*sg*alpha;ab+=(1-aa)*sb*alpha;aa+=(1-aa)*alpha;
                }
            }
            if(depth_buf) depth_buf[py*WIN_W+px] = vol_depth;
            ar+=(1-aa)*.05f;ag+=(1-aa)*.05f;ab+=(1-aa)*.07f;
            int ir=(int)(fminf(ar,1)*255),ig=(int)(fminf(ag,1)*255),ib=(int)(fminf(ab,1)*255);
            pixels[py*WIN_W+px]=(ir<<16)|(ig<<8)|ib;
        }
    }
}

/* Project 3D point (in [0,1]^3 volume space) to 2D pixel + depth */
static void project_point(float px, float py, float pz,
                          float cx, float cy, float cz,
                          float *fwd, float *right, float *cup,
                          float fov, int *sx, int *sy, float *depth) {
    float dx = px - cx, dy = py - cy, dz = pz - cz;
    float d = dx*fwd[0]+dy*fwd[1]+dz*fwd[2];
    *depth = d;
    if (d < 0.01f) { *sx = -1; *sy = -1; return; }
    float u = (dx*right[0]+dy*right[1]+dz*right[2]) / d;
    float v = (dx*cup[0]+dy*cup[1]+dz*cup[2]) / d;
    *sx = (int)((u / fov / (WIN_W/(float)WIN_H) + 1) * 0.5f * WIN_W);
    *sy = (int)((1 - v / fov) * 0.5f * WIN_H);
}

/* Draw a depth-tested line on the pixel buffer.
 * d0/d1 are the depths of the two endpoints (distance from camera along ray).
 * Only draws pixels where the edge is in front of the volume. */
static void draw_line(Uint32 *pixels, int x0, int y0, float d0,
                      int x1, int y1, float d1, Uint32 color) {
    int adx = abs(x1-x0), ady = abs(y1-y0);
    int sx = (x0<x1)?1:-1, sy = (y0<y1)?1:-1;
    int err = adx - ady;
    int steps = adx+ady+1;
    for (int i = 0; i < steps; i++) {
        if (x0>=0 && x0<WIN_W && y0>=0 && y0<WIN_H) {
            /* Interpolate depth along the line */
            float t = (steps > 1) ? (float)i / (steps-1) : 0;
            float edge_depth = d0 + t * (d1 - d0);
            float vol_d = depth_buf ? depth_buf[y0*WIN_W+x0] : 1e30f;
            if (edge_depth < vol_d) {
                /* Edge is in front of volume — draw bright */
                pixels[y0*WIN_W+x0] = color;
            } else {
                /* Edge is behind volume — draw dim, blended with existing */
                Uint32 existing = pixels[y0*WIN_W+x0];
                int er = (existing>>16)&0xff, eg = (existing>>8)&0xff, eb = existing&0xff;
                int cr = (color>>16)&0xff, cg = (color>>8)&0xff, cb = color&0xff;
                float blend = 0.25f;  /* how much of the edge shows through */
                int fr = (int)(er*(1-blend) + cr*blend);
                int fg = (int)(eg*(1-blend) + cg*blend);
                int fb = (int)(eb*(1-blend) + cb*blend);
                pixels[y0*WIN_W+x0] = (fr<<16)|(fg<<8)|fb;
            }
        }
        int e2 = 2*err;
        if (e2 > -ady) { err -= ady; x0 += sx; }
        if (e2 <  adx) { err += adx; y0 += sy; }
    }
}

/* Draw wireframe bounding box */
static void draw_wireframe(Uint32 *pixels,
                           float cx, float cy, float cz,
                           float *fwd, float *right, float *cup, float fov) {
    /* 8 corners of the unit cube [0,1]^3 */
    float corners[8][3] = {
        {0,0,0},{1,0,0},{1,1,0},{0,1,0},
        {0,0,1},{1,0,1},{1,1,1},{0,1,1}
    };
    int sx[8], sy[8]; float sd[8];
    for (int i = 0; i < 8; i++)
        project_point(corners[i][0], corners[i][1], corners[i][2],
                      cx, cy, cz, fwd, right, cup, fov, &sx[i], &sy[i], &sd[i]);

    /* 12 edges */
    int edges[12][2] = {
        {0,1},{1,2},{2,3},{3,0},  /* bottom face */
        {4,5},{5,6},{6,7},{7,4},  /* top face */
        {0,4},{1,5},{2,6},{3,7}   /* vertical edges */
    };
    Uint32 white = 0xCCCCCC;
    for (int e = 0; e < 12; e++) {
        int a = edges[e][0], b = edges[e][1];
        if (sx[a] >= 0 && sx[b] >= 0)
            draw_line(pixels, sx[a], sy[a], sd[a], sx[b], sy[b], sd[b], white);
    }
}

static void print_help(void) {
    printf("\n");
    printf("=== SFA Volume Viewer Controls ===\n");
    printf("  Mouse drag     Rotate camera\n");
    printf("  Scroll / Up/Dn Zoom in/out\n");
    printf("  Left/Right, A/D  Previous/next frame\n");
    printf("  Home / End     First/last frame\n");
    printf("  Space          Play/pause animation\n");
    printf("  1 / 2 / 3     Toggle Red/Green/Blue channel\n");
    printf("  + / -          Brightness up/down\n");
    printf("  O / P          Opacity down/up\n");
    printf("  W              Toggle wireframe bounding box\n");
    printf("  R              Reset view\n");
    printf("  S              Save screenshot (screenshot.ppm)\n");
    printf("  Q / Escape     Quit\n");
    printf("\n");
    printf("  Channels: RED=|P| (binding), GREEN=phi^2 (fabric), BLUE=theta^2 (angle)\n");
    printf("\n");
}

int main(int argc, char **argv) {
    if (argc < 2) {
        printf("Usage: %s <archive.sfa>\n", argv[0]);
        print_help();
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
    print_help();

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
    depth_buf = malloc(WIN_W * WIN_H * sizeof(float));
    for (int i = 0; i < WIN_W*WIN_H; i++) depth_buf[i] = 1e30f;

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
                case SDLK_w: show_wireframe=!show_wireframe; need_render=1; break;
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
            /* Wireframe overlay */
            if (show_wireframe) {
                float cx_w=cam_dist*sinf(cam_theta)*cosf(cam_phi)+.5f;
                float cy_w=cam_dist*sinf(cam_theta)*sinf(cam_phi)+.5f;
                float cz_w=cam_dist*cosf(cam_theta)+.5f;
                float fw[3]={-cx_w+.5f,-cy_w+.5f,-cz_w+.5f};
                float ln=sqrtf(fw[0]*fw[0]+fw[1]*fw[1]+fw[2]*fw[2]);
                fw[0]/=ln;fw[1]/=ln;fw[2]/=ln;
                float up_w[3]={0,0,1};
                float rt[3]={fw[1]*up_w[2]-fw[2]*up_w[1],fw[2]*up_w[0]-fw[0]*up_w[2],fw[0]*up_w[1]-fw[1]*up_w[0]};
                ln=sqrtf(rt[0]*rt[0]+rt[1]*rt[1]+rt[2]*rt[2]);
                if(ln<1e-6f){rt[0]=1;rt[1]=0;rt[2]=0;ln=1;}
                rt[0]/=ln;rt[1]/=ln;rt[2]/=ln;
                float cu[3]={rt[1]*fw[2]-rt[2]*fw[1],rt[2]*fw[0]-rt[0]*fw[2],rt[0]*fw[1]-rt[1]*fw[0]};
                draw_wireframe(pixels, cx_w, cy_w, cz_w, fw, rt, cu, 0.8f);
            }
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

    free(pixels); free(vol_r); free(vol_g); free(vol_b); free(depth_buf);
    SDL_DestroyTexture(tex); SDL_DestroyRenderer(ren); SDL_DestroyWindow(win);
    SDL_Quit();
    sfa_close(archive);
    return 0;
}
