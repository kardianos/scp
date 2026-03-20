/*  volview.c — Interactive volumetric viewer for 6-field Cosserat snapshots
 *
 *  Software raymarcher with SDL2 window.
 *  Colors:
 *    RED   = bound field  (|P| = |φ₀φ₁φ₂|)
 *    GREEN = unbound field (Σφ² where |P| small)
 *    BLUE  = angular field (Σθ²)
 *
 *  Controls:
 *    Mouse drag: rotate view
 *    Scroll / +/-: zoom
 *    1/2/3: toggle R/G/B channels
 *    T: cycle transfer function (density, MIP, emission)
 *    R: reset view
 *    S: save screenshot (screenshot.ppm)
 *    Q/ESC: quit
 *
 *  Build: gcc -O3 -o volview viewer/volview.c -lSDL2 -lm
 *  Usage: ./volview data/cosserat_mt0.0/field_t0300.bin
 */

#include <SDL2/SDL.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define WIN_W 800
#define WIN_H 600
#define MAX_STEPS 256
#define PI 3.14159265358979

/* ---- Volume data ---- */
static int VOL_N;
static float *vol_r, *vol_g, *vol_b;  /* normalized 0-1, N³ each */

/* ---- Camera ---- */
static float cam_theta = 0.6f, cam_phi = 0.3f, cam_dist = 2.5f;
static int show_r = 1, show_g = 1, show_b = 1;
static int mode = 0;  /* 0=emission, 1=MIP, 2=absorption */
static float brightness = 3.0f;

/* ---- Load snapshot ---- */
static int load_snapshot(const char *path) {
    FILE *fp = fopen(path, "rb");
    if (!fp) { perror("fopen"); return -1; }

    int N; double L, t;
    fread(&N, sizeof(int), 1, fp);
    fread(&L, sizeof(double), 1, fp);
    fread(&t, sizeof(double), 1, fp);

    /* Try to read nf marker */
    int nf = 3;
    int maybe_nf;
    if (fread(&maybe_nf, sizeof(int), 1, fp) == 1) {
        if (maybe_nf == 6 || maybe_nf == 3) {
            nf = maybe_nf;
        } else {
            /* Not a field count marker — rewind 4 bytes */
            fseek(fp, -4, SEEK_CUR);
        }
    }

    long N3 = (long)N * N * N;
    VOL_N = N;

    printf("Loading: N=%d L=%.1f t=%.1f nf=%d\n", N, L, t, nf);

    /* Read all fields */
    double *fields[6];
    for (int a = 0; a < nf; a++) {
        fields[a] = malloc(N3 * sizeof(double));
        fread(fields[a], sizeof(double), N3, fp);
    }
    for (int a = nf; a < 6; a++) {
        fields[a] = calloc(N3, sizeof(double));
    }
    fclose(fp);

    /* Compute channels */
    vol_r = malloc(N3 * sizeof(float));
    vol_g = malloc(N3 * sizeof(float));
    vol_b = malloc(N3 * sizeof(float));

    /* Find max values for normalization */
    double max_P = 0, max_phi2 = 0, max_theta2 = 0;
    for (long i = 0; i < N3; i++) {
        double p0 = fields[0][i], p1 = fields[1][i], p2 = fields[2][i];
        double absP = fabs(p0 * p1 * p2);
        double phi2 = p0*p0 + p1*p1 + p2*p2;
        double theta2 = fields[3][i]*fields[3][i] + fields[4][i]*fields[4][i] + fields[5][i]*fields[5][i];
        if (absP > max_P) max_P = absP;
        if (phi2 > max_phi2) max_phi2 = phi2;
        if (theta2 > max_theta2) max_theta2 = theta2;
    }

    printf("  max |P|=%.4f  max Σφ²=%.4f  max Σθ²=%.4f\n", max_P, max_phi2, max_theta2);

    /* Use percentile-based normalization */
    double norm_P = max_P * 0.3;
    double norm_phi2 = max_phi2 * 0.15;
    double norm_theta2 = (max_theta2 > 0) ? max_theta2 * 0.3 : 1.0;

    for (long i = 0; i < N3; i++) {
        double p0 = fields[0][i], p1 = fields[1][i], p2 = fields[2][i];
        double absP = fabs(p0 * p1 * p2);
        double phi2 = p0*p0 + p1*p1 + p2*p2;
        double P_norm = absP / norm_P;

        /* RED: bound field */
        vol_r[i] = (float)fmin(P_norm, 1.0);

        /* GREEN: unbound field (suppress where bound) */
        double unbound = phi2 * (1.0 - fmin(P_norm, 1.0));
        vol_g[i] = (float)fmin(unbound / norm_phi2, 1.0);

        /* BLUE: angular field */
        double theta2 = fields[3][i]*fields[3][i] + fields[4][i]*fields[4][i] + fields[5][i]*fields[5][i];
        vol_b[i] = (float)fmin(theta2 / norm_theta2, 1.0);
    }

    for (int a = 0; a < 6; a++) free(fields[a]);
    printf("  Volume loaded: %d³, channels normalized\n", N);
    return 0;
}

/* ---- Sample volume (trilinear) ---- */
static void sample_vol(float x, float y, float z, float *r, float *g, float *b) {
    /* x,y,z in [0,1] → grid coords */
    float fx = x * (VOL_N - 1), fy = y * (VOL_N - 1), fz = z * (VOL_N - 1);
    int ix = (int)fx, iy = (int)fy, iz = (int)fz;
    if (ix < 0 || ix >= VOL_N-1 || iy < 0 || iy >= VOL_N-1 || iz < 0 || iz >= VOL_N-1) {
        *r = *g = *b = 0; return;
    }
    float dx = fx - ix, dy = fy - iy, dz = fz - iz;
    int NN = VOL_N, N2 = NN*NN;

    /* Trilinear interpolation */
    #define S(arr, i, j, k) (arr)[(i)*N2 + (j)*NN + (k)]
    float w000 = (1-dx)*(1-dy)*(1-dz), w001 = (1-dx)*(1-dy)*dz;
    float w010 = (1-dx)*dy*(1-dz), w011 = (1-dx)*dy*dz;
    float w100 = dx*(1-dy)*(1-dz), w101 = dx*(1-dy)*dz;
    float w110 = dx*dy*(1-dz), w111 = dx*dy*dz;

    *r = w000*S(vol_r,ix,iy,iz) + w001*S(vol_r,ix,iy,iz+1)
       + w010*S(vol_r,ix,iy+1,iz) + w011*S(vol_r,ix,iy+1,iz+1)
       + w100*S(vol_r,ix+1,iy,iz) + w101*S(vol_r,ix+1,iy,iz+1)
       + w110*S(vol_r,ix+1,iy+1,iz) + w111*S(vol_r,ix+1,iy+1,iz+1);
    *g = w000*S(vol_g,ix,iy,iz) + w001*S(vol_g,ix,iy,iz+1)
       + w010*S(vol_g,ix,iy+1,iz) + w011*S(vol_g,ix,iy+1,iz+1)
       + w100*S(vol_g,ix+1,iy,iz) + w101*S(vol_g,ix+1,iy,iz+1)
       + w110*S(vol_g,ix+1,iy+1,iz) + w111*S(vol_g,ix+1,iy+1,iz+1);
    *b = w000*S(vol_b,ix,iy,iz) + w001*S(vol_b,ix,iy,iz+1)
       + w010*S(vol_b,ix,iy+1,iz) + w011*S(vol_b,ix,iy+1,iz+1)
       + w100*S(vol_b,ix+1,iy,iz) + w101*S(vol_b,ix+1,iy,iz+1)
       + w110*S(vol_b,ix+1,iy+1,iz) + w111*S(vol_b,ix+1,iy+1,iz+1);
    #undef S
}

/* ---- Render frame ---- */
static void render(Uint32 *pixels) {
    float cx = cam_dist * sinf(cam_theta) * cosf(cam_phi);
    float cy = cam_dist * sinf(cam_theta) * sinf(cam_phi);
    float cz = cam_dist * cosf(cam_theta);

    /* Camera basis vectors */
    float fwd[3] = {-cx, -cy, -cz};
    float len = sqrtf(fwd[0]*fwd[0]+fwd[1]*fwd[1]+fwd[2]*fwd[2]);
    fwd[0]/=len; fwd[1]/=len; fwd[2]/=len;

    float up[3] = {0, 0, 1};
    float right[3] = {
        fwd[1]*up[2] - fwd[2]*up[1],
        fwd[2]*up[0] - fwd[0]*up[2],
        fwd[0]*up[1] - fwd[1]*up[0]
    };
    len = sqrtf(right[0]*right[0]+right[1]*right[1]+right[2]*right[2]);
    if (len < 1e-6) { right[0]=1; right[1]=0; right[2]=0; len=1; }
    right[0]/=len; right[1]/=len; right[2]/=len;

    float cam_up[3] = {
        right[1]*fwd[2] - right[2]*fwd[1],
        right[2]*fwd[0] - right[0]*fwd[2],
        right[0]*fwd[1] - right[1]*fwd[0]
    };

    float step_size = 1.4f / MAX_STEPS;  /* volume is in [0,1]³ */
    float fov = 0.8f;

    #pragma omp parallel for schedule(dynamic, 4)
    for (int py = 0; py < WIN_H; py++) {
        for (int px = 0; px < WIN_W; px++) {
            float u = (2.0f * px / WIN_W - 1.0f) * fov * WIN_W / WIN_H;
            float v = (1.0f - 2.0f * py / WIN_H) * fov;

            float dir[3] = {
                fwd[0] + u*right[0] + v*cam_up[0],
                fwd[1] + u*right[1] + v*cam_up[1],
                fwd[2] + u*right[2] + v*cam_up[2]
            };
            len = sqrtf(dir[0]*dir[0]+dir[1]*dir[1]+dir[2]*dir[2]);
            dir[0]/=len; dir[1]/=len; dir[2]/=len;

            /* Ray origin: camera position, offset to volume center [0.5,0.5,0.5] */
            float ox = cx + 0.5f, oy = cy + 0.5f, oz = cz + 0.5f;

            /* Find intersection with [0,1]³ box */
            float tmin = 0, tmax = 100;
            for (int d = 0; d < 3; d++) {
                float o = (d==0)?ox:(d==1)?oy:oz;
                float di = (d==0)?dir[0]:(d==1)?dir[1]:dir[2];
                if (fabsf(di) < 1e-8) continue;
                float t0 = (0 - o) / di, t1 = (1 - o) / di;
                if (t0 > t1) { float tmp=t0; t0=t1; t1=tmp; }
                if (t0 > tmin) tmin = t0;
                if (t1 < tmax) tmax = t1;
            }

            float acc_r = 0, acc_g = 0, acc_b = 0, acc_a = 0;

            if (tmin < tmax && tmax > 0) {
                if (tmin < 0) tmin = 0;

                for (float tt = tmin; tt < tmax && acc_a < 0.98f; tt += step_size) {
                    float sx = ox + tt*dir[0];
                    float sy = oy + tt*dir[1];
                    float sz = oz + tt*dir[2];

                    float sr, sg, sb;
                    sample_vol(sx, sy, sz, &sr, &sg, &sb);

                    sr *= show_r; sg *= show_g; sb *= show_b;

                    /* Transfer function: emission-absorption */
                    float density = (sr + sg + sb) * brightness;
                    float alpha = 1.0f - expf(-density * step_size * 8.0f);
                    if (alpha < 0.001f) continue;

                    /* Apply gamma for visibility */
                    sr = powf(sr, 0.6f) * brightness;
                    sg = powf(sg, 0.6f) * brightness;
                    sb = powf(sb, 0.6f) * brightness;

                    acc_r += (1-acc_a) * sr * alpha;
                    acc_g += (1-acc_a) * sg * alpha;
                    acc_b += (1-acc_a) * sb * alpha;
                    acc_a += (1-acc_a) * alpha;
                }
            }

            /* Background: dark gray */
            acc_r += (1-acc_a) * 0.05f;
            acc_g += (1-acc_a) * 0.05f;
            acc_b += (1-acc_a) * 0.07f;

            int ir = (int)(fminf(acc_r, 1.0f) * 255);
            int ig = (int)(fminf(acc_g, 1.0f) * 255);
            int ib = (int)(fminf(acc_b, 1.0f) * 255);
            pixels[py * WIN_W + px] = (ir << 16) | (ig << 8) | ib;
        }
    }
}

static void save_screenshot(Uint32 *pixels) {
    FILE *fp = fopen("screenshot.ppm", "wb");
    fprintf(fp, "P6\n%d %d\n255\n", WIN_W, WIN_H);
    for (int i = 0; i < WIN_W*WIN_H; i++) {
        unsigned char rgb[3] = {
            (pixels[i]>>16)&0xff, (pixels[i]>>8)&0xff, pixels[i]&0xff
        };
        fwrite(rgb, 1, 3, fp);
    }
    fclose(fp);
    printf("Screenshot saved: screenshot.ppm\n");
}

int main(int argc, char **argv) {
    if (argc < 2) {
        printf("Usage: %s <snapshot.bin>\n", argv[0]);
        printf("Controls: drag=rotate, scroll=zoom, 1/2/3=toggle R/G/B, S=screenshot, Q=quit\n");
        return 1;
    }

    if (load_snapshot(argv[1]) < 0) return 1;

    SDL_Init(SDL_INIT_VIDEO);
    SDL_Window *win = SDL_CreateWindow("SCP Volume Viewer — RED:bound GREEN:fabric BLUE:angle",
        SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, WIN_W, WIN_H, 0);
    SDL_Renderer *ren = SDL_CreateRenderer(win, -1, SDL_RENDERER_ACCELERATED);
    SDL_Texture *tex = SDL_CreateTexture(ren, SDL_PIXELFORMAT_RGB888,
        SDL_TEXTUREACCESS_STREAMING, WIN_W, WIN_H);

    Uint32 *pixels = malloc(WIN_W * WIN_H * sizeof(Uint32));

    int running = 1, dragging = 0, need_render = 1;

    while (running) {
        SDL_Event ev;
        while (SDL_PollEvent(&ev)) {
            switch (ev.type) {
            case SDL_QUIT: running = 0; break;
            case SDL_KEYDOWN:
                switch (ev.key.keysym.sym) {
                case SDLK_q: case SDLK_ESCAPE: running = 0; break;
                case SDLK_1: show_r = !show_r; need_render = 1; break;
                case SDLK_2: show_g = !show_g; need_render = 1; break;
                case SDLK_3: show_b = !show_b; need_render = 1; break;
                case SDLK_s: save_screenshot(pixels); break;
                case SDLK_r: cam_theta=0.6; cam_phi=0.3; cam_dist=2.5; brightness=3.0; need_render=1; break;
                case SDLK_EQUALS: case SDLK_PLUS: brightness *= 1.3f; need_render=1; break;
                case SDLK_MINUS: brightness /= 1.3f; need_render=1; break;
                case SDLK_UP: cam_dist *= 0.9f; need_render=1; break;
                case SDLK_DOWN: cam_dist *= 1.1f; need_render=1; break;
                }
                break;
            case SDL_MOUSEBUTTONDOWN:
                if (ev.button.button == SDL_BUTTON_LEFT) dragging = 1;
                break;
            case SDL_MOUSEBUTTONUP:
                if (ev.button.button == SDL_BUTTON_LEFT) dragging = 0;
                break;
            case SDL_MOUSEMOTION:
                if (dragging) {
                    cam_phi += ev.motion.xrel * 0.01f;
                    cam_theta -= ev.motion.yrel * 0.01f;
                    if (cam_theta < 0.1f) cam_theta = 0.1f;
                    if (cam_theta > PI-0.1f) cam_theta = PI-0.1f;
                    need_render = 1;
                }
                break;
            case SDL_MOUSEWHEEL:
                cam_dist *= (ev.wheel.y > 0) ? 0.9f : 1.1f;
                if (cam_dist < 0.5f) cam_dist = 0.5f;
                if (cam_dist > 10.0f) cam_dist = 10.0f;
                need_render = 1;
                break;
            }
        }

        if (need_render) {
            render(pixels);
            SDL_UpdateTexture(tex, NULL, pixels, WIN_W * sizeof(Uint32));
            SDL_RenderCopy(ren, tex, NULL, NULL);
            SDL_RenderPresent(ren);
            need_render = 0;
            printf("  [render] θ=%.2f φ=%.2f d=%.2f bright=%.1f R%d G%d B%d\n",
                   cam_theta, cam_phi, cam_dist, brightness, show_r, show_g, show_b);
        }

        SDL_Delay(16);
    }

    free(pixels); free(vol_r); free(vol_g); free(vol_b);
    SDL_DestroyTexture(tex);
    SDL_DestroyRenderer(ren);
    SDL_DestroyWindow(win);
    SDL_Quit();
    return 0;
}
