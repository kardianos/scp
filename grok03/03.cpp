// CMakeLists.txt
cmake_minimum_required(VERSION 3.15)
project(FieldSimulation)
find_package(FFTW3 REQUIRED)
find_package(HDF5 REQUIRED)
find_package(Open3D REQUIRED)  // Or VTK if preferred
add_executable(sim main.cpp)
target_link_libraries(sim FFTW3::fftw3 HDF5::HDF5 Open3D::Open3D)
target_compile_features(sim PRIVATE cxx_std_17)
target_compile_options(sim PRIVATE -O3 -fopenmp)  // Enable OpenMP

// main.cpp
#include <fftw3.h>
#include <hdf5.h>
#include <open3d/Open3D.h>  // For UI/rendering
#include <thread>
#include <mutex>
#include <condition_variable>
#include <vector>
#include <complex>
#include <omp.h>  // For multithreading

constexpr int N = 32;  // Grid size
constexpr double L = 10.0, dx = L / N, dt = 0.01, k_param = 1.0, lam = 0.1, A = 1.0;
constexpr double threshold = 0.2;

// 3D array type (use std::vector or Eigen for real code)
using Field = std::vector<std::vector<std::vector<double>>>;  // [N][N][N]

// Global state for threading
std::mutex mtx;
std::condition_variable cv;
bool recording_done = false;
std::vector<Field> frames;  // In-memory buffer (or load from HDF5 on demand)
int current_frame = 0;
bool auto_advance = true, reverse = false;

// FFTW plans (thread-safe)
fftw_plan fwd_plan, bwd_plan;

// Initialize grid and plans
Field init_phi() {
    Field phi(N, std::vector<std::vector<double>>(N, std::vector<double>(N, 0.0)));
    #pragma omp parallel for  // Multithread init
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            for (int k = 0; k < N; ++k) {
                double x = i * dx, y = j * dx, z = k * dx;
                phi[i][j][k] = A * sin(2 * M_PI * x / L) * sin(2 * M_PI * y / L) * sin(2 * M_PI * z / L);
            }
        }
    }
    return phi;
}

// Compute Laplacian via FFT (multithreaded FFTW)
Field compute_laplacian(const Field& phi) {
    // Allocate complex buffers
    std::complex<double> in[N*N*N], out[N*N*N];
    // Copy phi to in (parallel)
    #pragma omp parallel for
    for (int i = 0; i < N*N*N; ++i) in[i] = /* flatten phi */;
    
    fftw_execute(fwd_plan);  // Forward FFT (threaded)
    
    // Apply -k^2 in freq space (parallel)
    #pragma omp parallel for
    for (int idx = 0; idx < N*N*N; ++idx) {
        // Compute kx,ky,kz and multiply out[idx] *= -(kx^2 + ky^2 + kz^2)
    }
    
    fftw_execute(bwd_plan);  // Backward FFT
    Field lap(N /*...*/);  // Copy real part of out to lap (parallel)
    return lap;
}

// Time-step (leapfrog, multithreaded nonlinear term)
void step(Field& phi, Field& vel) {
    Field accel = compute_laplacian(phi);
    #pragma omp parallel for collapse(3)
    for (int i = 0; i < N; ++i) for (int j = 0; j < N; ++j) for (int k = 0; k < N; ++k) {
        accel[i][j][k] -= k_param * phi[i][j][k] - lam * pow(phi[i][j][k], 3);
    }
    // Update vel and phi (similar parallel loops)
    // ...
}

// Simulation thread: Run sim, record frames to HDF5/in-memory
void simulation_thread(hid_t file_id) {
    Field phi = init_phi(), vel(N /* zeros */);
    double t = 0.0;
    bool emitted = false;
    fftw_init_threads();  // Enable FFTW multithreading
    omp_set_num_threads(std::thread::hardware_concurrency());
    
    // Create FFT plans
    fwd_plan = fftw_plan_dft_3d(N, N, N, reinterpret_cast<fftw_complex*>(in), reinterpret_cast<fftw_complex*>(out), FFTW_FORWARD, FFTW_ESTIMATE);
    bwd_plan = fftw_plan_dft_3d(N, N, N, reinterpret_cast<fftw_complex*>(out), reinterpret_cast<fftw_complex*>(in), FFTW_BACKWARD, FFTW_ESTIMATE);
    
    while (!recording_done) {  // Run until stopped
        if (!emitted && t >= 1.0) {
            // Scale phi, add wave (parallel loops)
            emitted = true;
        }
        step(phi, vel);
        t += dt;
        
        {  // Record frame
            std::lock_guard<std::mutex> lock(mtx);
            frames.push_back(phi);  // Or write to HDF5 dataset
            cv.notify_one();
        }
    }
    fftw_destroy_plan(fwd_plan); fftw_destroy_plan(bwd_plan);
    fftw_cleanup_threads();
}

// Navigation thread: Handle reverse/seek (non-blocking)
void navigation_thread() {
    // Poll user input (e.g., via console or UI callbacks)
    while (true) {
        // If reverse toggled, set reverse = true;
        // If seek to frame X, set current_frame = X;
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
}

// UI thread: Render with Open3D, auto-advance from recording
void ui_thread() {
    auto viewer = std::make_shared<open3d::visualization::Visualizer>();
    viewer->CreateVisualizerWindow("Volumetric Sim", 800, 800);
    
    while (viewer->PollEvents()) {
        std::unique_lock<std::mutex> lock(mtx);
        cv.wait(lock, [&] { return current_frame < frames.size() || !auto_advance; });
        
        if (auto_advance) {
            current_frame = frames.size() - 1;  // Latest frame
        } else if (reverse) {
            current_frame = std::max(0, current_frame - 1);
        }
        
        // Create voxel grid from frames[current_frame]
        auto voxel_grid = std::make_shared<open3d::geometry::VoxelGrid>();
        // Loop over field, add voxels where abs(phi) > threshold, color by value
        #pragma omp parallel for
        for (/*...*/) {
            if (std::abs(frames[current_frame][i][j][k]) > threshold) {
                // Add voxel with color (red/green based on norm)
            }
        }
        viewer->ClearGeometries();
        viewer->AddGeometry(voxel_grid);
        viewer->UpdateGeometry();
        viewer->UpdateRender();
    }
}

int main() {
    // Open HDF5 file for recording (optional, if not in-memory)
    hid_t file_id = H5Fcreate("recording.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    
    std::thread sim_thread(simulation_thread, file_id);
    std::thread nav_thread(navigation_thread);
    std::thread ui_thread(ui_thread);  // Drives off recording
    
    // Run until user stops (e.g., Ctrl+C)
    ui_thread.join();  // Wait on UI
    recording_done = true;
    sim_thread.join();
    nav_thread.join();
    H5Fclose(file_id);
    return 0;
}