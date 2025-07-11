cmake_minimum_required(VERSION 3.15)
project(FieldSimulation)

find_package(FFTW3 REQUIRED)
find_package(HDF5 REQUIRED)
find_package(Open3D REQUIRED)  # Assumes Open3D is installed; use find_package(Open3D HINTS /usr/local/lib/cmake/Open3D) if custom path
find_package(OpenMP REQUIRED)  # For multithreading

add_executable(sim main.cpp)
target_link_libraries(sim PRIVATE FFTW3::fftw3 HDF5::HDF5 Open3D::Open3D OpenMP::OpenMP_CXX)
target_compile_features(sim PRIVATE cxx_std_17)
target_compile_options(sim PRIVATE -O3 -fopenmp)  # Enable OpenMP and optimization