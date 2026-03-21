#!/usr/bin/env python3
"""
Computation 1: θ Field Frequency Spectrum Around the Braid

Reads frames from V34 sfa_hires.sfa, extracts θ_φ(t) at various radii,
computes FFT power spectrum. Looks for discrete resonant frequencies.
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'sfa', 'format'))
import struct
import numpy as np
from numpy.fft import fft, fftfreq

SFA_PATH = "../v34/torsion_coupling/data/sfa_hires.sfa"

def read_sfa_header(path):
    """Read SFA header to get grid/frame info."""
    with open(path, 'rb') as f:
        type_tag = f.read(4)
        size = struct.unpack('<Q', f.read(8))[0]
        version = struct.unpack('<I', f.read(4))[0]
        flags = struct.unpack('<I', f.read(4))[0]
        Nx = struct.unpack('<I', f.read(4))[0]
        Ny = struct.unpack('<I', f.read(4))[0]
        Nz = struct.unpack('<I', f.read(4))[0]
        Lx = struct.unpack('<d', f.read(8))[0]
        Ly = struct.unpack('<d', f.read(8))[0]
        Lz = struct.unpack('<d', f.read(8))[0]
        dt_sim = struct.unpack('<d', f.read(8))[0]
        n_columns = struct.unpack('<I', f.read(4))[0]
        total_frames = struct.unpack('<I', f.read(4))[0]
        first_jtop = struct.unpack('<Q', f.read(8))[0]
        cdef_offset = struct.unpack('<Q', f.read(8))[0]
        jtop_max = struct.unpack('<I', f.read(4))[0]
        jmpf_max = struct.unpack('<I', f.read(4))[0]
    return {
        'Nx': Nx, 'Ny': Ny, 'Nz': Nz,
        'Lx': Lx, 'Ly': Ly, 'Lz': Lz,
        'dt': dt_sim, 'n_columns': n_columns,
        'total_frames': total_frames,
        'first_jtop': first_jtop, 'jmpf_max': jmpf_max,
        'flags': flags,
    }

def read_sfa_frame_times(path, header):
    """Read all frame times from the JMPF."""
    with open(path, 'rb') as f:
        # Read JTOP
        f.seek(header['first_jtop'])
        f.read(12)  # chunk header
        jtop_max = struct.unpack('<I', f.read(4))[0]
        jtop_cur = struct.unpack('<I', f.read(4))[0]
        next_jtop = struct.unpack('<Q', f.read(8))[0]

        # First JMPF offset
        jmpf_offset = struct.unpack('<Q', f.read(8))[0]
        f.read(8)  # first_frame + frame_count

        # Read JMPF
        f.seek(jmpf_offset)
        f.read(12)  # chunk header
        jmpf_max = struct.unpack('<I', f.read(4))[0]
        jmpf_cur = struct.unpack('<I', f.read(4))[0]

        times = []
        offsets = []
        sizes = []
        for i in range(jmpf_cur):
            t = struct.unpack('<d', f.read(8))[0]
            off = struct.unpack('<Q', f.read(8))[0]
            sz = struct.unpack('<Q', f.read(8))[0]
            cksum = struct.unpack('<I', f.read(4))[0]
            reserved = struct.unpack('<I', f.read(4))[0]
            times.append(t)
            offsets.append(off)
            sizes.append(sz)

    return times, offsets, sizes

def read_sfa_frame(path, offset, comp_size, N_total, n_columns, flags):
    """Read and decompress a single frame."""
    import zlib
    try:
        import pyzstd
        has_zstd = True
    except ImportError:
        try:
            import zstd as pyzstd
            has_zstd = True
        except ImportError:
            has_zstd = False

    with open(path, 'rb') as f:
        f.seek(offset + 12)  # skip FRMD chunk header
        compressed = f.read(comp_size)

    frame_bytes = n_columns * N_total * 8

    if has_zstd:
        raw = pyzstd.decompress(compressed)
    else:
        # Try raw zstd via subprocess
        import subprocess
        p = subprocess.run(['zstd', '-d', '--stdout'],
                          input=compressed, capture_output=True)
        raw = p.stdout

    # BSS decode if flag set
    codec = flags & 0xF
    if codec == 2:  # BSS + zstd
        raw_bytes = bytearray(raw)
        decoded = bytearray(len(raw_bytes))
        # BSS decode per column (each column is N_total float64 = N_total*8 bytes)
        col_offset = 0
        for c in range(n_columns):
            elem_size = 8  # float64
            n_values = N_total
            col_bytes = n_values * elem_size
            src = raw_bytes[col_offset:col_offset + col_bytes]
            dst_start = col_offset
            for b in range(elem_size):
                for i in range(n_values):
                    decoded[dst_start + i*elem_size + b] = src[b*n_values + i]
            col_offset += col_bytes
        raw = bytes(decoded)

    # Parse into columns
    data = np.frombuffer(raw, dtype=np.float64)
    columns = []
    for c in range(n_columns):
        columns.append(data[c*N_total:(c+1)*N_total].copy())

    return columns


def main():
    os.makedirs('data', exist_ok=True)

    print("=== θ Frequency Spectrum Analysis ===\n")

    # Check for the SFA file
    if not os.path.exists(SFA_PATH):
        print(f"SFA not found at {SFA_PATH}")
        print("Falling back to synthetic analysis from V34 timeseries data")
        analyze_from_timeseries()
        return

    header = read_sfa_header(SFA_PATH)
    N = header['Nx']
    L = header['Lx']
    n_col = header['n_columns']
    n_frames = header['total_frames']
    N_total = N * N * N
    dx = 2*L/(N-1)

    print(f"SFA: N={N}, L={L}, {n_col} columns, {n_frames} frames")
    print(f"Grid: dx={dx:.4f}")

    times, offsets, sizes = read_sfa_frame_times(SFA_PATH, header)
    print(f"Frame times: {times[0]:.2f} to {times[-1]:.2f}, dt_frame={times[1]-times[0]:.4f}\n")

    dt_frame = times[1] - times[0]

    # Read ALL frames, extract θ_φ(r, t) in cylindrical shells
    radii = [3, 5, 8, 12, 16, 20]  # radial shells to sample
    dr = 1.5  # shell width

    # Storage: theta_phi[r_idx][frame] = azimuthally-averaged θ_φ
    theta_phi_t = {r: [] for r in radii}
    theta_r_t = {r: [] for r in radii}
    theta_z_t = {r: [] for r in radii}
    phi2_t = {r: [] for r in radii}

    # First, find braid center from first frame
    cols0 = read_sfa_frame(SFA_PATH, offsets[0], sizes[0], N_total, n_col, header['flags'])
    phi2_all = cols0[0]**2 + cols0[1]**2 + cols0[2]**2
    avg = phi2_all.mean()
    thresh = 5 * avg
    mask = phi2_all > thresh

    coords = np.mgrid[0:N, 0:N, 0:N]
    w = phi2_all.reshape(N,N,N) * mask.reshape(N,N,N)
    ci = int(round((coords[0]*w).sum() / (w.sum()+1e-30)))
    cj = int(round((coords[1]*w).sum() / (w.sum()+1e-30)))
    ck = int(round((coords[2]*w).sum() / (w.sum()+1e-30)))
    cx = -L + ci*dx
    cy = -L + cj*dx
    cz = -L + ck*dx
    print(f"Braid center: ({cx:.1f}, {cy:.1f}, {cz:.1f})")

    # Pre-compute cylindrical radius for each grid point
    r_perp = np.zeros(N_total)
    angle = np.zeros(N_total)
    for ix in range(N):
        x = -L + ix*dx - cx
        for iy in range(N):
            y = -L + iy*dx - cy
            for iz in range(N):
                idx = ix*N*N + iy*N + iz
                r_perp[idx] = np.sqrt(x*x + y*y)
                angle[idx] = np.arctan2(y, x)

    print(f"\nReading {n_frames} frames...")

    # Read frames and extract radial time series
    for fi in range(n_frames):
        if fi % 50 == 0:
            print(f"  Frame {fi}/{n_frames} (t={times[fi]:.2f})")

        cols = read_sfa_frame(SFA_PATH, offsets[fi], sizes[fi],
                              N_total, n_col, header['flags'])

        # θ components: cols[3]=θ_x, cols[4]=θ_y, cols[5]=θ_z
        theta_x = cols[3] if n_col >= 6 else np.zeros(N_total)
        theta_y = cols[4] if n_col >= 6 else np.zeros(N_total)
        theta_z_arr = cols[5] if n_col >= 6 else np.zeros(N_total)

        # Cylindrical decomposition
        cos_a = np.cos(angle)
        sin_a = np.sin(angle)
        t_r = theta_x * cos_a + theta_y * sin_a       # radial
        t_phi = -theta_x * sin_a + theta_y * cos_a     # azimuthal
        t_z = theta_z_arr                                # axial

        phi2 = cols[0]**2 + cols[1]**2 + cols[2]**2

        for r_target in radii:
            shell = (r_perp > r_target - dr/2) & (r_perp < r_target + dr/2)
            n_pts = shell.sum()
            if n_pts > 0:
                theta_phi_t[r_target].append(t_phi[shell].mean())
                theta_r_t[r_target].append(t_r[shell].mean())
                theta_z_t[r_target].append(t_z[shell].mean())
                phi2_t[r_target].append(phi2[shell].mean())
            else:
                theta_phi_t[r_target].append(0)
                theta_r_t[r_target].append(0)
                theta_z_t[r_target].append(0)
                phi2_t[r_target].append(0)

    # FFT analysis
    print(f"\n=== Frequency Spectra ===\n")

    freq = fftfreq(n_frames, d=dt_frame)
    pos_mask = freq > 0
    freq_pos = freq[pos_mask]

    with open('data/theta_spectra.tsv', 'w') as f:
        f.write("radius\tfreq\tP_theta_phi\tP_theta_r\tP_theta_z\tP_phi2\n")

        for r_target in radii:
            sig = np.array(theta_phi_t[r_target])
            sig -= sig.mean()  # remove DC
            spectrum = np.abs(fft(sig))**2
            P_phi = spectrum[pos_mask]

            sig_r = np.array(theta_r_t[r_target])
            sig_r -= sig_r.mean()
            P_r = np.abs(fft(sig_r))**2
            P_r = P_r[pos_mask]

            sig_z = np.array(theta_z_t[r_target])
            sig_z -= sig_z.mean()
            P_z = np.abs(fft(sig_z))**2
            P_z = P_z[pos_mask]

            sig_p2 = np.array(phi2_t[r_target])
            sig_p2 -= sig_p2.mean()
            P_p2 = np.abs(fft(sig_p2))**2
            P_p2 = P_p2[pos_mask]

            # Find peaks
            top_idx = np.argsort(P_phi)[-5:][::-1]
            print(f"r={r_target}: Top θ_φ frequencies:")
            for idx in top_idx:
                print(f"  ω={freq_pos[idx]:.4f} (T={1/freq_pos[idx]:.2f})  "
                      f"P_φ={P_phi[idx]:.2e}  P_r={P_r[idx]:.2e}  P_z={P_z[idx]:.2e}")

            for i in range(len(freq_pos)):
                f.write(f"{r_target}\t{freq_pos[i]:.6f}\t{P_phi[i]:.6e}\t"
                        f"{P_r[i]:.6e}\t{P_z[i]:.6e}\t{P_p2[i]:.6e}\n")
            print()

    # Summary: are there discrete peaks shared across radii?
    print("=== Peak Frequency Summary ===\n")
    print(f"{'Radius':>6s}  {'ω₁':>8s}  {'ω₂':>8s}  {'ω₃':>8s}  {'ω₁/ω₂':>8s}")
    print("-" * 50)
    for r_target in radii:
        sig = np.array(theta_phi_t[r_target])
        sig -= sig.mean()
        spectrum = np.abs(fft(sig))**2
        P = spectrum[pos_mask]
        top3 = np.argsort(P)[-3:][::-1]
        omegas = freq_pos[top3]
        ratio = omegas[0]/omegas[1] if omegas[1] > 0 else 0
        print(f"{r_target:6d}  {omegas[0]:8.4f}  {omegas[1]:8.4f}  {omegas[2]:8.4f}  {ratio:8.3f}")

    print("\nIf ω₁/ω₂ ≈ constant across radii → discrete mode (like a drum)")
    print("If ω₁/ω₂ varies with radius → continuous (dispersive wave)")
    print("If ω₁ is the same at all radii → global oscillation (braid breathing)")


def analyze_from_timeseries():
    """Fallback: analyze θ_rms(t) from V34 timeseries."""
    import glob
    ts_files = glob.glob("../v34/torsion_coupling/data/cosserat_mt0.0/timeseries.tsv")
    if not ts_files:
        ts_files = glob.glob("../v34/torsion_coupling/data/cosserat_eta_0.5/timeseries.tsv")
    if not ts_files:
        print("No timeseries data found")
        return

    path = ts_files[0]
    print(f"Using timeseries: {path}")
    data = np.loadtxt(path, skiprows=1)
    t = data[:, 0]
    theta_rms = data[:, -1]

    # FFT of θ_rms(t)
    dt = t[1] - t[0]
    sig = theta_rms - theta_rms.mean()
    spectrum = np.abs(fft(sig))**2
    freq = fftfreq(len(sig), d=dt)
    pos = freq > 0

    print("\nTop frequencies in θ_rms(t):")
    P = spectrum[pos]
    top5 = np.argsort(P)[-5:][::-1]
    for idx in top5:
        print(f"  ω={freq[pos][idx]:.4f} (T={1/freq[pos][idx]:.2f})  P={P[idx]:.2e}")


if __name__ == '__main__':
    main()
