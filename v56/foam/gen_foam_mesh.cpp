// gen_foam_mesh.cpp — generate a periodic 3D Voronoi mesh on [-L, L]^3
//
// Pipeline:
//   1. Bridson Poisson-disk sampling on [0, 2L)^3 with periodic distance
//   2. Hand the points to voro++ container_periodic
//   3. For each cell: extract volume, neighbor IDs, face areas, face deltas
//   4. Write binary foam_mesh.bin in the format described in DESIGN.md
//
// Build:
//   make -C voro_src libvoro++.a
//   g++ -O3 -Wall -Ivoro_src -o gen_foam_mesh gen_foam_mesh.cpp \
//       voro_src/libvoro++.a
//
// Usage:
//   ./gen_foam_mesh -L 20 -N 50000 -o foam_mesh.bin

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <ctime>
#include <vector>
#include <random>
#include <algorithm>

#include "voro++.hh"

using namespace voro;
using std::vector;

struct Vec3 { double x, y, z; };

// Wrap to [0, S) along each axis
static inline double wrap(double v, double S) {
    while (v < 0) v += S;
    while (v >= S) v -= S;
    return v;
}

// Periodic distance squared in cube [0, S)^3
static inline double dist2_periodic(const Vec3& a, const Vec3& b, double S) {
    double dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
    if (dx >  0.5*S) dx -= S; if (dx < -0.5*S) dx += S;
    if (dy >  0.5*S) dy -= S; if (dy < -0.5*S) dy += S;
    if (dz >  0.5*S) dz -= S; if (dz < -0.5*S) dz += S;
    return dx*dx + dy*dy + dz*dz;
}

// Bridson Poisson-disk in [0, S)^3 with periodic wrap.
// r_min: minimum distance between samples. k_tries: attempts per active.
static vector<Vec3> poisson_disk_periodic(
        double S, double r_min, int k_tries, std::mt19937_64& rng) {
    std::uniform_real_distribution<double> U01(0.0, 1.0);
    double gcell = r_min / std::sqrt(3.0);  // grid cell holds ≤ 1 sample
    int Ng = (int)std::ceil(S / gcell);
    gcell = S / Ng;  // exact fit so periodic wrap is clean
    long Ng3 = (long)Ng * Ng * Ng;
    vector<int> grid(Ng3, -1);

    auto gpos = [&](const Vec3& p) {
        int gx = (int)std::floor(p.x / gcell); if (gx < 0) gx = 0; if (gx >= Ng) gx = Ng-1;
        int gy = (int)std::floor(p.y / gcell); if (gy < 0) gy = 0; if (gy >= Ng) gy = Ng-1;
        int gz = (int)std::floor(p.z / gcell); if (gz < 0) gz = 0; if (gz >= Ng) gz = Ng-1;
        return (long)gx * Ng * Ng + gy * Ng + gz;
    };

    auto fits = [&](const Vec3& p, const vector<Vec3>& samples) {
        int gx = (int)std::floor(p.x / gcell);
        int gy = (int)std::floor(p.y / gcell);
        int gz = (int)std::floor(p.z / gcell);
        double r2 = r_min * r_min;
        for (int dx = -2; dx <= 2; dx++) {
            for (int dy = -2; dy <= 2; dy++) {
                for (int dz = -2; dz <= 2; dz++) {
                    int gi = ((gx + dx) % Ng + Ng) % Ng;
                    int gj = ((gy + dy) % Ng + Ng) % Ng;
                    int gk = ((gz + dz) % Ng + Ng) % Ng;
                    int s = grid[(long)gi * Ng * Ng + gj * Ng + gk];
                    if (s < 0) continue;
                    if (dist2_periodic(p, samples[s], S) < r2) return false;
                }
            }
        }
        return true;
    };

    vector<Vec3> samples;
    vector<int> active;

    Vec3 p0 { S * U01(rng), S * U01(rng), S * U01(rng) };
    samples.push_back(p0);
    grid[gpos(p0)] = 0;
    active.push_back(0);

    while (!active.empty()) {
        std::uniform_int_distribution<int> Ui(0, (int)active.size() - 1);
        int slot = Ui(rng);
        int i = active[slot];
        bool found = false;
        for (int t = 0; t < k_tries; t++) {
            double theta = U01(rng) * M_PI;
            double phi   = U01(rng) * 2.0 * M_PI;
            double r     = r_min * (1.0 + U01(rng));  // [r_min, 2 r_min]
            Vec3 p {
                samples[i].x + r * std::sin(theta) * std::cos(phi),
                samples[i].y + r * std::sin(theta) * std::sin(phi),
                samples[i].z + r * std::cos(theta)
            };
            p.x = wrap(p.x, S); p.y = wrap(p.y, S); p.z = wrap(p.z, S);
            if (fits(p, samples)) {
                samples.push_back(p);
                grid[gpos(p)] = (int)samples.size() - 1;
                active.push_back((int)samples.size() - 1);
                found = true;
                break;
            }
        }
        if (!found) {
            active[slot] = active.back();
            active.pop_back();
        }
    }
    return samples;
}

// Main mesh-build pipeline
struct FaceRec {
    uint32_t cell_a, cell_b;
    double area;
    double delta[3];
};

int main(int argc, char** argv) {
    double L = 20.0;
    int N_target = 20000;
    const char* output = "foam_mesh.bin";
    uint64_t seed = 42;

    for (int i = 1; i < argc; i++) {
        if (!std::strcmp(argv[i], "-L")) L = std::atof(argv[++i]);
        else if (!std::strcmp(argv[i], "-N")) N_target = std::atoi(argv[++i]);
        else if (!std::strcmp(argv[i], "-o")) output = argv[++i];
        else if (!std::strcmp(argv[i], "--seed")) seed = std::strtoull(argv[++i], nullptr, 10);
        else { std::fprintf(stderr, "unknown arg: %s\n", argv[i]); return 1; }
    }

    double S = 2.0 * L;       // box size
    double box_vol = S * S * S;
    double cell_vol = box_vol / N_target;
    // Bridson saturates near 0.65 V / r_min^3 samples in 3D, so to hit
    // N_target we set r_min = (0.65 × cell_vol)^(1/3).
    double r_min = std::pow(0.65 * cell_vol, 1.0/3.0);

    std::printf("Foam mesh generator\n");
    std::printf("  L = %.3f, S = box size = %.3f\n", L, S);
    std::printf("  N target = %d, cell_vol target = %.3f, r_min = %.3f\n",
                N_target, cell_vol, r_min);

    // ---- Poisson-disk sampling ----
    std::mt19937_64 rng(seed);
    std::time_t t0 = std::time(nullptr);
    vector<Vec3> samples = poisson_disk_periodic(S, r_min, 30, rng);
    int N = (int)samples.size();
    std::printf("  Poisson-disk: %d samples in %lds\n",
                N, (long)(std::time(nullptr) - t0));

    if (N < N_target * 0.5) {
        std::fprintf(stderr, "WARNING: only got %d/%d samples — reduce r_min or increase k_tries\n",
                     N, N_target);
    }

    // ---- Build periodic Voronoi via voro++ ----
    // Block grid: cube root of N for memory efficiency
    int n_blocks = std::max(3, (int)std::cbrt((double)N));
    container_periodic con(
        S, 0.0, S, 0.0, 0.0, S,             // bx, bxy, by, bxz, byz, bz
        n_blocks, n_blocks, n_blocks, 8);   // n_x, n_y, n_z, init_mem
    for (int i = 0; i < N; i++)
        con.put(i, samples[i].x, samples[i].y, samples[i].z);

    double vvol = con.sum_cell_volumes();
    std::printf("  Voronoi volumes computed: %.3f (expected %.3f)\n",
                vvol, box_vol);

    // ---- Extract per-cell data ----
    vector<double> volumes(N, 0.0);
    vector<FaceRec> faces;
    faces.reserve((size_t)N * 8);  // avg 14 neighbors → ~7 faces unique

    voronoicell_neighbor c;
    c_loop_all_periodic cl(con);
    if (cl.start()) {
        do {
            if (!con.compute_cell(c, cl)) continue;
            int id = cl.pid();
            double cx, cy, cz;
            cl.pos(cx, cy, cz);

            volumes[id] = c.volume();

            vector<int> neigh;
            vector<double> areas;
            c.neighbors(neigh);
            c.face_areas(areas);

            for (size_t f = 0; f < neigh.size(); f++) {
                int nb = neigh[f];
                if (nb < 0) continue;        // boundary (shouldn't happen periodic)
                if (nb == id) continue;      // self-neighbor (degenerate)
                if (nb < id) continue;       // dedupe: store each face once

                // delta: shortest-image vector from cell `id` to cell `nb`
                double dx = samples[nb].x - cx;
                double dy = samples[nb].y - cy;
                double dz = samples[nb].z - cz;
                if (dx >  0.5*S) dx -= S; if (dx < -0.5*S) dx += S;
                if (dy >  0.5*S) dy -= S; if (dy < -0.5*S) dy += S;
                if (dz >  0.5*S) dz -= S; if (dz < -0.5*S) dz += S;

                FaceRec fr;
                fr.cell_a = (uint32_t)id;
                fr.cell_b = (uint32_t)nb;
                fr.area = areas[f];
                fr.delta[0] = dx; fr.delta[1] = dy; fr.delta[2] = dz;
                faces.push_back(fr);
            }
        } while (cl.inc());
    }

    std::printf("  Faces (unique): %zu, avg degree = %.1f\n",
                faces.size(), 2.0 * (double)faces.size() / N);

    // ---- Build CSR cell→face index ----
    vector<uint32_t> cell_deg(N, 0);
    for (const auto& f : faces) {
        cell_deg[f.cell_a]++;
        cell_deg[f.cell_b]++;
    }
    vector<uint32_t> offsets(N + 1, 0);
    for (int i = 0; i < N; i++) offsets[i+1] = offsets[i] + cell_deg[i];
    vector<uint32_t> face_indices(offsets[N], 0);
    vector<uint32_t> cursor = offsets;  // copy first N entries as cursor
    for (size_t fid = 0; fid < faces.size(); fid++) {
        uint32_t a = faces[fid].cell_a, b = faces[fid].cell_b;
        face_indices[cursor[a]++] = (uint32_t)fid;
        face_indices[cursor[b]++] = (uint32_t)fid;
    }

    // ---- Write binary file ----
    // Convert sample positions back to [-L, L] for output
    FILE* fp = std::fopen(output, "wb");
    if (!fp) { std::fprintf(stderr, "cannot open %s\n", output); return 1; }

    // Header (32 bytes): magic(4) ver(4) L(8) N(4) Nf(4) reserved(8)
    std::fwrite("FOAM", 1, 4, fp);
    uint32_t version = 1; std::fwrite(&version, 4, 1, fp);
    std::fwrite(&L, 8, 1, fp);
    uint32_t N_u = (uint32_t)N; std::fwrite(&N_u, 4, 1, fp);
    uint32_t Nf  = (uint32_t)faces.size(); std::fwrite(&Nf, 4, 1, fp);
    uint64_t reserved = 0; std::fwrite(&reserved, 8, 1, fp);

    // Cell records (32 bytes each): pos[3] (24) + volume (8)
    for (int i = 0; i < N; i++) {
        double pos[4] = {
            samples[i].x - L,
            samples[i].y - L,
            samples[i].z - L,
            volumes[i]
        };
        std::fwrite(pos, 8, 4, fp);
    }

    // Face records (32 bytes each, padded for alignment):
    //   cell_a(4) cell_b(4) area(8) delta[3](24)  = 40 bytes; pad to 40
    // Note: DESIGN.md said 28 bytes — actually 40 is cleaner. Update.
    for (const auto& f : faces) {
        std::fwrite(&f.cell_a, 4, 1, fp);
        std::fwrite(&f.cell_b, 4, 1, fp);
        std::fwrite(&f.area, 8, 1, fp);
        std::fwrite(f.delta, 8, 3, fp);
    }

    // CSR offsets (4 × (N+1) bytes)
    std::fwrite(offsets.data(), 4, N + 1, fp);
    // Face index list (4 × 2*Nf bytes)
    std::fwrite(face_indices.data(), 4, face_indices.size(), fp);

    std::fclose(fp);

    long bytes = 32 + 32L*N + 40L*Nf + 4L*(N+1) + 4L*face_indices.size();
    std::printf("  Wrote %s (%.1f MB)\n", output, bytes / 1e6);

    // Volume distribution sanity
    double vmin = volumes[0], vmax = volumes[0], vsum = 0;
    for (double v : volumes) {
        if (v < vmin) vmin = v;
        if (v > vmax) vmax = v;
        vsum += v;
    }
    std::printf("  V_cell: mean=%.3f min=%.3f max=%.3f (target %.3f)\n",
                vsum / N, vmin, vmax, cell_vol);
    std::printf("  V_total = %.3f (box = %.3f)\n", vsum, box_vol);

    return 0;
}
