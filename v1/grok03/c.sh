#!/bin/bash
g++ -std=c++17 -O3 -fopenmp 03.cpp -o 03 \
-I/usr/local/include/Open3D \
-I/usr/local/include/open3d/3rdparty \
-I/usr/include/hdf5/openmpi \
-I/usr/lib/x86_64-linux-gnu/openmpi/include \
-L/usr/local/lib \
-L/usr/lib/x86_64-linux-gnu/hdf5/openmpi \
-lOpen3D \
-lfftw3 \
-lfftw3_omp \
-lhdf5 \
-lz \
-I/usr/include/hdf5/openmpi -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi /usr/lib/x86_64-linux-gnu/hdf5/openmpi/libhdf5_hl.a /usr/lib/x86_64-linux-gnu/hdf5/openmpi/libhdf5.a -lcrypto -lcurl -lsz -lz -ldl -lm