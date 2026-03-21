#!/usr/bin/env python3
"""SFA Compression Strategy Benchmark

Tests multiple compression pipelines on float64 volumetric field data
matching Cosserat simulation characteristics:
  - 64^3 grid, 6 float64 columns per frame
  - cos() background (A~0.1) + Gaussian braid perturbation (A~0.8)
  - Spatially smooth, temporally slowly evolving

Also loads real SFA data from a Cosserat simulation if available.
"""

import numpy as np
import struct
import time
import sys
import os

import pyzstd

# ============================================================
# 1. DATA GENERATION
# ============================================================

def generate_synthetic_frame(N=64, L=15.0, t=0.0, seed=42):
    """Generate one frame of 6 float64 fields mimicking Cosserat simulation.

    Fields: phi_x, phi_y, phi_z (position), theta_x, theta_y, theta_z (angle)
    Background: A_bg * cos(kz + phase)
    Braid core: Gaussian envelope * A_core * sin(...)
    """
    rng = np.random.RandomState(seed)
    dx = 2*L / (N-1)
    x = np.linspace(-L, L, N)
    y = np.linspace(-L, L, N)
    z = np.linspace(-L, L, N)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

    columns = []
    for col_idx in range(6):
        # Background: low-amplitude cosine oscillation
        A_bg = 0.1
        k = 2*np.pi / (2*L) * (1 + 0.3*col_idx)  # slightly different k per field
        phase = 0.5 * col_idx + 0.1*t
        bg = A_bg * np.cos(k*Z + phase)

        # Braid core: Gaussian envelope centered at origin, high amplitude
        A_core = 0.8
        sigma = 2.0  # core width
        r2 = X**2 + Y**2 + Z**2
        envelope = np.exp(-r2 / (2*sigma**2))
        core_phase = k*Z + 0.7*col_idx + 0.3*t
        core = A_core * envelope * np.sin(core_phase)

        # Small noise to prevent exact floating-point patterns
        noise = rng.normal(0, 1e-6, (N, N, N))

        field = bg + core + noise
        columns.append(field.astype(np.float64).ravel())

    return columns


def generate_frame_sequence(n_frames=10, N=64, L=15.0):
    """Generate a temporal sequence of frames with slow evolution."""
    frames = []
    for i in range(n_frames):
        t = i * 1.0  # dt=1.0
        cols = generate_synthetic_frame(N, L, t, seed=42+i)
        frames.append(cols)
    return frames


def load_sfa_frames(path, max_frames=11):
    """Load frames from a real SFA file (binary parsing)."""
    with open(path, 'rb') as f:
        # Read SFAH header
        chunk_type = f.read(4)
        chunk_size = struct.unpack('<Q', f.read(8))[0]

        version = struct.unpack('<I', f.read(4))[0]
        flags = struct.unpack('<I', f.read(4))[0]
        Nx, Ny, Nz = struct.unpack('<III', f.read(12))
        Lx, Ly, Lz = struct.unpack('<ddd', f.read(24))
        dt = struct.unpack('<d', f.read(8))[0]
        n_columns = struct.unpack('<I', f.read(4))[0]
        total_frames = struct.unpack('<I', f.read(4))[0]
        first_jtop_offset = struct.unpack('<Q', f.read(8))[0]
        cdef_offset = struct.unpack('<Q', f.read(8))[0]
        jtop_max = struct.unpack('<I', f.read(4))[0]
        jmpf_max = struct.unpack('<I', f.read(4))[0]

        codec = flags & 0xF
        N_total = Nx * Ny * Nz

        print(f"SFA: {Nx}x{Ny}x{Nz}, {n_columns} cols, {total_frames} frames, codec={codec}")

        # Read CDEF to get dtypes
        f.seek(cdef_offset)
        ct = f.read(4)
        cs = struct.unpack('<Q', f.read(8))[0]
        dtype_sizes = []
        for c in range(n_columns):
            name = f.read(12).rstrip(b'\x00').decode()
            dtype_code = struct.unpack('<B', f.read(1))[0]
            semantic = struct.unpack('<B', f.read(1))[0]
            component = struct.unpack('<B', f.read(1))[0]
            col_flags = struct.unpack('<B', f.read(1))[0]
            scale = struct.unpack('<d', f.read(8))[0]
            ds = [2,4,8,16,1,2,4,8,1,2,4,8][dtype_code]
            dtype_sizes.append(ds)

        frame_bytes = sum(ds * N_total for ds in dtype_sizes)

        # Read JTOP -> JMPF -> FRMD offsets
        f.seek(first_jtop_offset + 12)  # skip chunk header
        jtop_max_r = struct.unpack('<I', f.read(4))[0]
        jtop_cur = struct.unpack('<I', f.read(4))[0]
        next_jtop = struct.unpack('<Q', f.read(8))[0]

        # Read first L1 entry
        jmpf_offset = struct.unpack('<Q', f.read(8))[0]
        first_frame = struct.unpack('<I', f.read(4))[0]
        frame_count = struct.unpack('<I', f.read(4))[0]

        # Read JMPF entries
        f.seek(jmpf_offset + 12)  # skip chunk header
        jmpf_max_r = struct.unpack('<I', f.read(4))[0]
        jmpf_cur = struct.unpack('<I', f.read(4))[0]

        n_to_read = min(jmpf_cur, max_frames)
        frmd_entries = []
        for i in range(n_to_read):
            t_val = struct.unpack('<d', f.read(8))[0]
            offset = struct.unpack('<Q', f.read(8))[0]
            comp_size = struct.unpack('<Q', f.read(8))[0]
            checksum = struct.unpack('<I', f.read(4))[0]
            reserved = struct.unpack('<I', f.read(4))[0]
            frmd_entries.append((t_val, offset, comp_size, checksum))

        # Read and decompress frames
        frames = []
        for i, (t_val, offset, comp_size, checksum) in enumerate(frmd_entries):
            f.seek(offset + 12)  # skip FRMD chunk header (type + size)
            compressed = f.read(comp_size)

            if codec == 2:  # BSS + zstd
                bss_data = pyzstd.decompress(compressed)
                # Reverse BSS: grouped by byte position -> interleaved
                raw = bytearray(frame_bytes)
                col_offset = 0
                for c in range(n_columns):
                    es = dtype_sizes[c]
                    n_vals = N_total
                    for b in range(es):
                        for j in range(n_vals):
                            raw[col_offset + j*es + b] = bss_data[col_offset + b*n_vals + j]
                    col_offset += n_vals * es
            elif codec == 1:  # zstd only
                raw = pyzstd.decompress(compressed)
            else:
                f.seek(offset + 12)
                raw = f.read(frame_bytes)

            # Split into columns
            cols = []
            off = 0
            for c in range(n_columns):
                nbytes = N_total * dtype_sizes[c]
                arr = np.frombuffer(bytes(raw[off:off+nbytes]), dtype=np.float64)
                cols.append(arr.copy())
                off += nbytes
            frames.append(cols)

        return frames, (Nx, Ny, Nz, n_columns)


# ============================================================
# 2. COMPRESSION PRIMITIVES
# ============================================================

def bss_encode(data_bytes, n_values, elem_size):
    """Byte Stream Split: transpose bytes for float compression."""
    arr = np.frombuffer(data_bytes, dtype=np.uint8).reshape(n_values, elem_size)
    transposed = arr.T.copy()  # (elem_size, n_values) - groups by byte position
    return transposed.tobytes()


def bss_decode(data_bytes, n_values, elem_size):
    """Reverse Byte Stream Split."""
    arr = np.frombuffer(data_bytes, dtype=np.uint8).reshape(elem_size, n_values)
    return arr.T.copy().tobytes()


def delta_encode_f64(data):
    """Delta encoding on float64 (NOT bitwise lossless due to FP rounding).
    For lossless, use XOR encoding instead."""
    arr = np.frombuffer(data, dtype=np.float64)
    delta = np.empty_like(arr)
    delta[0] = arr[0]
    delta[1:] = arr[1:] - arr[:-1]
    return delta.tobytes()


def delta_decode_f64(data):
    """Reverse delta encoding. NOT bitwise lossless on float64."""
    arr = np.frombuffer(data, dtype=np.float64).copy()
    # cumsum is equivalent to sequential add, both suffer from FP rounding
    for i in range(1, len(arr)):
        arr[i] += arr[i-1]
    return arr.tobytes()


def idelta_encode(data):
    """Integer-domain delta: subtract successive values as int64 (LOSSLESS)."""
    arr = np.frombuffer(data, dtype=np.int64)
    delta = np.empty_like(arr)
    delta[0] = arr[0]
    delta[1:] = arr[1:] - arr[:-1]
    return delta.tobytes()


def idelta_decode(data):
    """Reverse integer delta (LOSSLESS)."""
    arr = np.frombuffer(data, dtype=np.int64).copy()
    for i in range(1, len(arr)):
        arr[i] += arr[i-1]
    return arr.tobytes()


def xor_encode_f64(data):
    """XOR successive values (as uint64). Similar float64s have many shared bits."""
    arr = np.frombuffer(data, dtype=np.uint64)
    xored = np.empty_like(arr)
    xored[0] = arr[0]
    xored[1:] = arr[1:] ^ arr[:-1]
    return xored.tobytes()


def xor_decode_f64(data):
    """Reverse XOR encoding."""
    arr = np.frombuffer(data, dtype=np.uint64).copy()
    for i in range(1, len(arr)):
        arr[i] ^= arr[i-1]
    return arr.tobytes()


def temporal_delta_encode(frames_bytes):
    """Temporal delta (float64): frame[n] - frame[n-1]. NOT bitwise lossless."""
    result = [frames_bytes[0]]
    for i in range(1, len(frames_bytes)):
        cur = np.frombuffer(frames_bytes[i], dtype=np.float64)
        prev = np.frombuffer(frames_bytes[i-1], dtype=np.float64)
        delta = cur - prev
        result.append(delta.tobytes())
    return result


def temporal_delta_decode(delta_frames_bytes):
    """Reverse temporal delta (float64). NOT bitwise lossless."""
    result = [delta_frames_bytes[0]]
    for i in range(1, len(delta_frames_bytes)):
        prev = np.frombuffer(result[i-1], dtype=np.float64)
        delta = np.frombuffer(delta_frames_bytes[i], dtype=np.float64)
        cur = prev + delta
        result.append(cur.tobytes())
    return result


def temporal_xor_encode(frames_bytes):
    """Temporal XOR: frame[n] ^ frame[n-1] as uint64. LOSSLESS."""
    result = [frames_bytes[0]]
    for i in range(1, len(frames_bytes)):
        cur = np.frombuffer(frames_bytes[i], dtype=np.uint64)
        prev = np.frombuffer(frames_bytes[i-1], dtype=np.uint64)
        xored = cur ^ prev
        result.append(xored.tobytes())
    return result


def temporal_xor_decode(xor_frames_bytes):
    """Reverse temporal XOR. LOSSLESS."""
    result = [xor_frames_bytes[0]]
    for i in range(1, len(xor_frames_bytes)):
        prev = np.frombuffer(result[i-1], dtype=np.uint64)
        xored = np.frombuffer(xor_frames_bytes[i], dtype=np.uint64)
        cur = prev ^ xored
        result.append(cur.tobytes())
    return result


def temporal_isub_encode(frames_bytes):
    """Temporal integer subtraction: frame[n] - frame[n-1] as int64. LOSSLESS."""
    result = [frames_bytes[0]]
    for i in range(1, len(frames_bytes)):
        cur = np.frombuffer(frames_bytes[i], dtype=np.int64)
        prev = np.frombuffer(frames_bytes[i-1], dtype=np.int64)
        delta = cur - prev  # wrapping subtraction in int64
        result.append(delta.tobytes())
    return result


def temporal_isub_decode(isub_frames_bytes):
    """Reverse temporal int64 subtraction. LOSSLESS."""
    result = [isub_frames_bytes[0]]
    for i in range(1, len(isub_frames_bytes)):
        prev = np.frombuffer(result[i-1], dtype=np.int64)
        delta = np.frombuffer(isub_frames_bytes[i], dtype=np.int64)
        cur = prev + delta
        result.append(cur.tobytes())
    return result


def gorilla_encode_f64(data):
    """Simplified Gorilla/XOR prediction encoding.

    For each value, XOR with previous, then store:
    - leading zero count (1 byte)
    - meaningful bits length (1 byte)
    - meaningful bits (variable bytes)

    This is a simplified byte-aligned version (not bit-packed like real Gorilla).
    """
    arr = np.frombuffer(data, dtype=np.uint64)
    n = len(arr)

    # Just do XOR + count leading zeros as metadata
    # Pack as: [xor_values_bytes] since full Gorilla needs bit-packing
    # Instead, use XOR + BSS which captures the same insight
    xored = np.empty_like(arr)
    xored[0] = arr[0]
    xored[1:] = arr[1:] ^ arr[:-1]

    # Count leading zero bytes per value (for analysis)
    lz_counts = np.zeros(n, dtype=np.uint8)
    for i in range(n):
        v = int(xored[i])
        lz = 0
        for b in range(7, -1, -1):
            if (v >> (b*8)) & 0xFF == 0:
                lz += 1
            else:
                break
        lz_counts[i] = lz

    return xored.tobytes(), lz_counts


def frame_to_bytes(columns):
    """Concatenate column arrays into single byte buffer."""
    parts = []
    for col in columns:
        parts.append(col.tobytes())
    return b''.join(parts)


def frame_to_columns_bytes(columns):
    """Return list of per-column byte buffers."""
    return [col.tobytes() for col in columns]


# ============================================================
# 3. COMPRESSION PIPELINES
# ============================================================

def compress_raw_zstd(data, level=3):
    """Plain zstd, no preprocessing."""
    return pyzstd.compress(data, level)

def decompress_raw_zstd(compressed, orig_size):
    return pyzstd.decompress(compressed)


def compress_bss_zstd(data, n_values, elem_size=8, level=3):
    """BSS + zstd (current SFA default)."""
    bss = bss_encode(data, n_values, elem_size)
    return pyzstd.compress(bss, level)

def decompress_bss_zstd(compressed, n_values, elem_size=8):
    bss = pyzstd.decompress(compressed)
    return bss_decode(bss, n_values, elem_size)


def compress_delta_bss_zstd(data, n_values, elem_size=8, level=3):
    """Delta + BSS + zstd."""
    delta = delta_encode_f64(data)
    bss = bss_encode(delta, n_values, elem_size)
    return pyzstd.compress(bss, level)

def decompress_delta_bss_zstd(compressed, n_values, elem_size=8):
    bss = pyzstd.decompress(compressed)
    delta = bss_decode(bss, n_values, elem_size)
    return delta_decode_f64(delta)


def compress_xor_bss_zstd(data, n_values, elem_size=8, level=3):
    """XOR + BSS + zstd."""
    xored = xor_encode_f64(data)
    bss = bss_encode(xored, n_values, elem_size)
    return pyzstd.compress(bss, level)

def decompress_xor_bss_zstd(compressed, n_values, elem_size=8):
    bss = pyzstd.decompress(compressed)
    xored = bss_decode(bss, n_values, elem_size)
    return xor_decode_f64(xored)


def compress_idelta_bss_zstd(data, n_values, elem_size=8, level=3):
    """Integer delta + BSS + zstd (LOSSLESS)."""
    delta = idelta_encode(data)
    bss = bss_encode(delta, n_values, elem_size)
    return pyzstd.compress(bss, level)

def decompress_idelta_bss_zstd(compressed, n_values, elem_size=8):
    bss = pyzstd.decompress(compressed)
    delta = bss_decode(bss, n_values, elem_size)
    return idelta_decode(delta)


def compress_f32_bss_zstd(data, level=3):
    """Float32 downcast + BSS + zstd (LOSSY)."""
    arr64 = np.frombuffer(data, dtype=np.float64)
    arr32 = arr64.astype(np.float32)
    n_values = len(arr32)
    bss = bss_encode(arr32.tobytes(), n_values, 4)
    return pyzstd.compress(bss, level), arr32

def decompress_f32_bss_zstd(compressed, n_values):
    bss = pyzstd.decompress(compressed)
    raw32 = bss_decode(bss, n_values, 4)
    arr32 = np.frombuffer(raw32, dtype=np.float32)
    return arr32.astype(np.float64).tobytes()


def compress_f16_bss_zstd(data, level=3):
    """Float16 downcast + BSS + zstd (VERY LOSSY)."""
    arr64 = np.frombuffer(data, dtype=np.float64)
    arr16 = arr64.astype(np.float16)
    n_values = len(arr16)
    bss = bss_encode(arr16.tobytes(), n_values, 2)
    return pyzstd.compress(bss, level), arr16

def decompress_f16_bss_zstd(compressed, n_values):
    bss = pyzstd.decompress(compressed)
    raw16 = bss_decode(bss, n_values, 2)
    arr16 = np.frombuffer(raw16, dtype=np.float16)
    return arr16.astype(np.float64).tobytes()


def compress_block_quant(data, block_size=256, bits=8):
    """GGUF-inspired block quantization.

    Divide data into blocks of block_size values.
    Per block: store (min, scale) as float32, then quantize to `bits`-bit integers.
    """
    arr = np.frombuffer(data, dtype=np.float64)
    n = len(arr)
    n_blocks = (n + block_size - 1) // block_size

    # Pad to full blocks
    padded = np.zeros(n_blocks * block_size, dtype=np.float64)
    padded[:n] = arr
    padded = padded.reshape(n_blocks, block_size)

    max_val = (1 << bits) - 1

    # Per-block min and scale
    block_min = padded.min(axis=1).astype(np.float32)
    block_max = padded.max(axis=1).astype(np.float32)
    block_scale = (block_max - block_min).astype(np.float32)
    block_scale[block_scale == 0] = 1.0  # avoid div by zero

    # Quantize
    normalized = (padded - block_min[:, None]) / block_scale[:, None]
    quantized = np.clip(np.round(normalized * max_val), 0, max_val).astype(np.uint8 if bits <= 8 else np.uint16)

    # Pack: [n_blocks × (4 + 4)] + [n_blocks × block_size × ceil(bits/8)]
    header = np.column_stack([block_min.view(np.uint32), block_scale.view(np.uint32)]).tobytes()
    qdata = quantized.tobytes()

    return header + qdata, n, n_blocks


def decompress_block_quant(packed, n_orig, n_blocks, block_size=256, bits=8):
    """Reverse block quantization."""
    max_val = (1 << bits) - 1
    header_size = n_blocks * 8  # 4 bytes min + 4 bytes scale per block

    header = np.frombuffer(packed[:header_size], dtype=np.uint32).reshape(n_blocks, 2)
    block_min = header[:, 0].view(np.float32).astype(np.float64)
    block_scale = header[:, 1].view(np.float32).astype(np.float64)

    elem_size = 1 if bits <= 8 else 2
    dtype = np.uint8 if bits <= 8 else np.uint16
    qdata = np.frombuffer(packed[header_size:], dtype=dtype).reshape(n_blocks, block_size)

    # Dequantize
    values = qdata.astype(np.float64) / max_val * block_scale[:, None] + block_min[:, None]
    return values.ravel()[:n_orig].tobytes()


# ============================================================
# 4. BENCHMARK HARNESS
# ============================================================

def measure_compression(name, compress_fn, decompress_fn, raw_data, n_values,
                        original_f64=None, n_iters=3):
    """Time compress/decompress, measure ratio, and (if lossy) error metrics."""
    raw_size = len(raw_data)

    # Compress
    t0 = time.perf_counter()
    for _ in range(n_iters):
        compressed = compress_fn(raw_data)
    t_comp = (time.perf_counter() - t0) / n_iters

    # Handle lossy methods that return tuples
    if isinstance(compressed, tuple):
        compressed = compressed[0]

    comp_size = len(compressed)

    # Decompress
    t0 = time.perf_counter()
    for _ in range(n_iters):
        decompressed = decompress_fn(compressed)
    t_decomp = (time.perf_counter() - t0) / n_iters

    ratio = raw_size / comp_size
    comp_speed = (raw_size / 1e6) / t_comp if t_comp > 0 else 0
    decomp_speed = (raw_size / 1e6) / t_decomp if t_decomp > 0 else 0

    # Error metrics (compare to original float64)
    ref = original_f64 if original_f64 is not None else raw_data
    ref_arr = np.frombuffer(ref, dtype=np.float64)
    dec_arr = np.frombuffer(decompressed, dtype=np.float64)

    if len(dec_arr) >= len(ref_arr):
        dec_arr = dec_arr[:len(ref_arr)]

    if original_f64 is None:
        # Lossless check
        is_lossless = np.array_equal(ref_arr, dec_arr)
        max_err = 0.0
        mean_err = 0.0
        psnr = float('inf')
    else:
        is_lossless = False
        diff = np.abs(ref_arr - dec_arr)
        max_err = np.max(diff)
        mean_err = np.mean(diff)
        signal_range = np.max(np.abs(ref_arr))
        if signal_range > 0 and max_err > 0:
            psnr = 20 * np.log10(signal_range / np.sqrt(np.mean(diff**2)))
        else:
            psnr = float('inf')

    return {
        'name': name,
        'raw_bytes': raw_size,
        'comp_bytes': comp_size,
        'ratio': ratio,
        'comp_MB_s': comp_speed,
        'decomp_MB_s': decomp_speed,
        'lossless': is_lossless if original_f64 is None else False,
        'max_err': max_err,
        'mean_err': mean_err,
        'psnr_dB': psnr,
    }


def run_single_frame_bench(columns, label="synthetic"):
    """Benchmark all compression methods on a single frame."""
    raw = frame_to_bytes(columns)
    n_values = sum(len(c) for c in columns)
    raw_size_mb = len(raw) / 1e6

    print(f"\n{'='*72}")
    print(f"SINGLE FRAME BENCHMARK ({label})")
    print(f"{'='*72}")
    print(f"Raw size: {len(raw)} bytes ({raw_size_mb:.2f} MB)")
    print(f"Grid: {int(round(len(columns[0])**(1/3)))}^3, {len(columns)} columns, float64")
    print()

    results = []

    # --- A. Raw zstd at various levels ---
    for level in [1, 3, 6, 9, 15, 22]:
        r = measure_compression(
            f"zstd-{level}",
            lambda d, lv=level: compress_raw_zstd(d, lv),
            lambda c: decompress_raw_zstd(c, len(raw)),
            raw, n_values
        )
        results.append(r)

    # --- B. BSS + zstd at various levels ---
    for level in [1, 3, 6, 9, 15, 22]:
        r = measure_compression(
            f"bss+zstd-{level}",
            lambda d, lv=level: compress_bss_zstd(d, n_values, 8, lv),
            lambda c: decompress_bss_zstd(c, n_values, 8),
            raw, n_values
        )
        results.append(r)

    # --- C. Delta + BSS + zstd ---
    for level in [3, 9, 22]:
        r = measure_compression(
            f"delta+bss+zstd-{level}",
            lambda d, lv=level: compress_delta_bss_zstd(d, n_values, 8, lv),
            lambda c: decompress_delta_bss_zstd(c, n_values, 8),
            raw, n_values
        )
        results.append(r)

    # --- D. XOR + BSS + zstd ---
    for level in [3, 9, 22]:
        r = measure_compression(
            f"xor+bss+zstd-{level}",
            lambda d, lv=level: compress_xor_bss_zstd(d, n_values, 8, lv),
            lambda c: decompress_xor_bss_zstd(c, n_values, 8),
            raw, n_values
        )
        results.append(r)

    # --- D2. Integer delta + BSS + zstd (LOSSLESS) ---
    for level in [3, 9, 22]:
        r = measure_compression(
            f"idelta+bss+zstd-{level}",
            lambda d, lv=level: compress_idelta_bss_zstd(d, n_values, 8, lv),
            lambda c: decompress_idelta_bss_zstd(c, n_values, 8),
            raw, n_values
        )
        results.append(r)

    # --- E. Float32 + BSS + zstd (LOSSY) ---
    for level in [3, 9]:
        r = measure_compression(
            f"f32+bss+zstd-{level}",
            lambda d, lv=level: compress_f32_bss_zstd(d, lv)[0],
            lambda c: decompress_f32_bss_zstd(c, n_values),
            raw, n_values, original_f64=raw
        )
        results.append(r)

    # --- F. Float16 + BSS + zstd (VERY LOSSY) ---
    r = measure_compression(
        "f16+bss+zstd-3",
        lambda d: compress_f16_bss_zstd(d, 3)[0],
        lambda c: decompress_f16_bss_zstd(c, n_values),
        raw, n_values, original_f64=raw
    )
    results.append(r)

    # --- G. Block quantization (GGUF-style) + zstd ---
    for bits in [8, 6, 4]:
        block_size = 256
        def comp_bq(d, b=bits, bs=block_size):
            packed, n_orig, n_blocks = compress_block_quant(d, bs, b)
            return pyzstd.compress(packed, 3)

        def decomp_bq(c, b=bits, bs=block_size):
            packed = pyzstd.decompress(c)
            n_orig = n_values
            n_blocks = (n_orig + bs - 1) // bs
            return decompress_block_quant(packed, n_orig, n_blocks, bs, b)

        r = measure_compression(
            f"bq{bits}+zstd-3",
            comp_bq, decomp_bq,
            raw, n_values, original_f64=raw
        )
        results.append(r)

    # --- H. Gorilla-style: XOR + BSS (analysis) ---
    xored_data, lz_counts = gorilla_encode_f64(raw)
    avg_lz = np.mean(lz_counts)
    print(f"Gorilla XOR analysis: avg leading zero bytes = {avg_lz:.2f}/8")
    # The actual compression is same as xor+bss+zstd above

    print_results_table(results)
    return results


def run_temporal_bench(frames, label="synthetic"):
    """Benchmark temporal compression on a sequence of frames."""
    n_frames = len(frames)
    n_values_per_frame = sum(len(c) for c in frames[0])

    # Convert to raw bytes per frame
    raw_frames = [frame_to_bytes(cols) for cols in frames]
    raw_total = sum(len(r) for r in raw_frames)

    print(f"\n{'='*72}")
    print(f"TEMPORAL SEQUENCE BENCHMARK ({label}, {n_frames} frames)")
    print(f"{'='*72}")
    print(f"Total raw: {raw_total} bytes ({raw_total/1e6:.2f} MB)")
    print()

    results = []

    # --- Baseline: compress each frame independently with BSS+zstd ---
    def compress_indep(frames_bytes, level=3):
        return [compress_bss_zstd(f, n_values_per_frame, 8, level) for f in frames_bytes]

    t0 = time.perf_counter()
    indep_compressed = compress_indep(raw_frames)
    t_comp = time.perf_counter() - t0
    indep_size = sum(len(c) for c in indep_compressed)

    t0 = time.perf_counter()
    indep_decompressed = [decompress_bss_zstd(c, n_values_per_frame, 8) for c in indep_compressed]
    t_decomp = time.perf_counter() - t0

    # Verify lossless
    lossless = all(d == r for d, r in zip(indep_decompressed, raw_frames))

    results.append({
        'name': f'indep bss+zstd-3 ({n_frames}f)',
        'raw_bytes': raw_total,
        'comp_bytes': indep_size,
        'ratio': raw_total / indep_size,
        'comp_MB_s': (raw_total/1e6) / t_comp if t_comp > 0 else 0,
        'decomp_MB_s': (raw_total/1e6) / t_decomp if t_decomp > 0 else 0,
        'lossless': lossless,
        'max_err': 0, 'mean_err': 0, 'psnr_dB': float('inf'),
    })

    # --- Temporal delta + BSS + zstd ---
    t0 = time.perf_counter()
    delta_frames = temporal_delta_encode(raw_frames)
    delta_compressed = [compress_bss_zstd(f, n_values_per_frame, 8, 3) for f in delta_frames]
    t_comp = time.perf_counter() - t0
    delta_size = sum(len(c) for c in delta_compressed)

    t0 = time.perf_counter()
    delta_decompressed = [decompress_bss_zstd(c, n_values_per_frame, 8) for c in delta_compressed]
    restored = temporal_delta_decode(delta_decompressed)
    t_decomp = time.perf_counter() - t0

    lossless = all(r == o for r, o in zip(restored, raw_frames))

    results.append({
        'name': f'temporal-delta bss+zstd-3 ({n_frames}f)',
        'raw_bytes': raw_total,
        'comp_bytes': delta_size,
        'ratio': raw_total / delta_size,
        'comp_MB_s': (raw_total/1e6) / t_comp if t_comp > 0 else 0,
        'decomp_MB_s': (raw_total/1e6) / t_decomp if t_decomp > 0 else 0,
        'lossless': lossless,
        'max_err': 0, 'mean_err': 0, 'psnr_dB': float('inf'),
    })

    # --- Temporal XOR + BSS + zstd (LOSSLESS) ---
    t0 = time.perf_counter()
    xor_frames = temporal_xor_encode(raw_frames)
    txor_compressed = [compress_bss_zstd(f, n_values_per_frame, 8, 3) for f in xor_frames]
    t_comp = time.perf_counter() - t0
    txor_size = sum(len(c) for c in txor_compressed)

    t0 = time.perf_counter()
    txor_decompressed = [decompress_bss_zstd(c, n_values_per_frame, 8) for c in txor_compressed]
    restored2 = temporal_xor_decode(txor_decompressed)
    t_decomp = time.perf_counter() - t0

    lossless = all(r == o for r, o in zip(restored2, raw_frames))

    results.append({
        'name': f'temporal-xor bss+zstd-3 ({n_frames}f)',
        'raw_bytes': raw_total,
        'comp_bytes': txor_size,
        'ratio': raw_total / txor_size,
        'comp_MB_s': (raw_total/1e6) / t_comp if t_comp > 0 else 0,
        'decomp_MB_s': (raw_total/1e6) / t_decomp if t_decomp > 0 else 0,
        'lossless': lossless,
        'max_err': 0, 'mean_err': 0, 'psnr_dB': float('inf'),
    })

    # --- Temporal int-sub + BSS + zstd (LOSSLESS) ---
    t0 = time.perf_counter()
    isub_frames = temporal_isub_encode(raw_frames)
    tisub_compressed = [compress_bss_zstd(f, n_values_per_frame, 8, 3) for f in isub_frames]
    t_comp = time.perf_counter() - t0
    tisub_size = sum(len(c) for c in tisub_compressed)

    t0 = time.perf_counter()
    tisub_decompressed = [decompress_bss_zstd(c, n_values_per_frame, 8) for c in tisub_compressed]
    restored2b = temporal_isub_decode(tisub_decompressed)
    t_decomp = time.perf_counter() - t0

    lossless = all(r == o for r, o in zip(restored2b, raw_frames))

    results.append({
        'name': f'temporal-isub bss+zstd-3 ({n_frames}f)',
        'raw_bytes': raw_total,
        'comp_bytes': tisub_size,
        'ratio': raw_total / tisub_size,
        'comp_MB_s': (raw_total/1e6) / t_comp if t_comp > 0 else 0,
        'decomp_MB_s': (raw_total/1e6) / t_decomp if t_decomp > 0 else 0,
        'lossless': lossless,
        'max_err': 0, 'mean_err': 0, 'psnr_dB': float('inf'),
    })

    # --- Temporal delta + f32 + BSS + zstd (LOSSY) ---
    t0 = time.perf_counter()
    delta_frames3 = temporal_delta_encode(raw_frames)
    # Compress deltas as f32 (except frame 0 which is full values)
    tf32_compressed = []
    tf32_compressed.append(compress_f32_bss_zstd(delta_frames3[0], 3)[0])
    for f in delta_frames3[1:]:
        tf32_compressed.append(compress_f32_bss_zstd(f, 3)[0])
    t_comp = time.perf_counter() - t0
    tf32_size = sum(len(c) for c in tf32_compressed)

    t0 = time.perf_counter()
    tf32_decompressed = [decompress_f32_bss_zstd(c, n_values_per_frame) for c in tf32_compressed]
    restored3 = temporal_delta_decode(tf32_decompressed)
    t_decomp = time.perf_counter() - t0

    # Error analysis
    ref_all = np.concatenate([np.frombuffer(r, dtype=np.float64) for r in raw_frames])
    dec_all = np.concatenate([np.frombuffer(r, dtype=np.float64) for r in restored3])
    diff = np.abs(ref_all - dec_all)
    max_err = np.max(diff)
    mean_err = np.mean(diff)
    sig_range = np.max(np.abs(ref_all))
    rms_err = np.sqrt(np.mean(diff**2))
    psnr = 20 * np.log10(sig_range / rms_err) if rms_err > 0 else float('inf')

    results.append({
        'name': f'temporal-delta f32+bss+zstd-3 ({n_frames}f)',
        'raw_bytes': raw_total,
        'comp_bytes': tf32_size,
        'ratio': raw_total / tf32_size,
        'comp_MB_s': (raw_total/1e6) / t_comp if t_comp > 0 else 0,
        'decomp_MB_s': (raw_total/1e6) / t_decomp if t_decomp > 0 else 0,
        'lossless': False,
        'max_err': max_err, 'mean_err': mean_err, 'psnr_dB': psnr,
    })

    # --- Concat all frames, single BSS+zstd pass ---
    all_raw = b''.join(raw_frames)
    for level in [3, 9]:
        t0 = time.perf_counter()
        all_bss = bss_encode(all_raw, n_values_per_frame * n_frames, 8)
        all_comp = pyzstd.compress(all_bss, level)
        t_comp = time.perf_counter() - t0

        t0 = time.perf_counter()
        all_decomp_bss = pyzstd.decompress(all_comp)
        all_decomp = bss_decode(all_decomp_bss, n_values_per_frame * n_frames, 8)
        t_decomp = time.perf_counter() - t0

        results.append({
            'name': f'concat bss+zstd-{level} ({n_frames}f)',
            'raw_bytes': raw_total,
            'comp_bytes': len(all_comp),
            'ratio': raw_total / len(all_comp),
            'comp_MB_s': (raw_total/1e6) / t_comp if t_comp > 0 else 0,
            'decomp_MB_s': (raw_total/1e6) / t_decomp if t_decomp > 0 else 0,
            'lossless': all_decomp == all_raw,
            'max_err': 0, 'mean_err': 0, 'psnr_dB': float('inf'),
        })

    print_results_table(results)
    return results


def print_results_table(results):
    """Print formatted results table."""
    print(f"{'Method':<42} {'Ratio':>6} {'Comp':>8} {'Decomp':>8} {'Loss':>5} {'MaxErr':>10} {'PSNR':>8}")
    print(f"{'':42} {'(x)':>6} {'(MB/s)':>8} {'(MB/s)':>8} {'':>5} {'':>10} {'(dB)':>8}")
    print('-' * 95)
    for r in results:
        loss_str = 'no' if r['lossless'] else 'YES'
        max_err_str = f"{r['max_err']:.2e}" if r['max_err'] > 0 else '-'
        psnr_str = f"{r['psnr_dB']:.1f}" if r['psnr_dB'] < 1e6 else 'inf'
        print(f"{r['name']:<42} {r['ratio']:>6.2f} {r['comp_MB_s']:>8.1f} {r['decomp_MB_s']:>8.1f} {loss_str:>5} {max_err_str:>10} {psnr_str:>8}")


def analyze_data_characteristics(columns, label=""):
    """Print statistics about the data to understand compression potential."""
    print(f"\n--- Data Characteristics ({label}) ---")
    for i, col in enumerate(columns):
        arr = col if isinstance(col, np.ndarray) else np.frombuffer(col, dtype=np.float64)
        raw = arr.tobytes()

        # Byte entropy
        byte_arr = np.frombuffer(raw, dtype=np.uint8)
        counts = np.bincount(byte_arr, minlength=256)
        probs = counts / len(byte_arr)
        probs = probs[probs > 0]
        entropy = -np.sum(probs * np.log2(probs))

        # BSS byte entropy (per stream)
        bss = bss_encode(raw, len(arr), 8)
        bss_arr = np.frombuffer(bss, dtype=np.uint8).reshape(8, len(arr))
        bss_entropies = []
        for b in range(8):
            bc = np.bincount(bss_arr[b], minlength=256)
            bp = bc / len(arr)
            bp = bp[bp > 0]
            bss_entropies.append(-np.sum(bp * np.log2(bp)))

        # Delta statistics
        delta = np.diff(arr)

        print(f"  col {i}: range=[{arr.min():.4f}, {arr.max():.4f}], "
              f"std={arr.std():.4f}, byte_H={entropy:.2f}, "
              f"BSS_H=[{min(bss_entropies):.2f}..{max(bss_entropies):.2f}], "
              f"delta_std={delta.std():.6f}")


# ============================================================
# 5. MAIN
# ============================================================

def main():
    print("SFA Compression Strategy Benchmark")
    print("=" * 72)

    # --- Generate synthetic data ---
    print("\nGenerating synthetic test data (64^3 x 6 fields)...")
    syn_frames = generate_frame_sequence(n_frames=10, N=64, L=15.0)
    analyze_data_characteristics(syn_frames[0], "synthetic frame 0")

    # --- Load real SFA data if available ---
    sfa_path = "/home/d/code/scp/v34/torsion_coupling/data/compress_test.sfa"
    real_frames = None
    if os.path.exists(sfa_path):
        print(f"\nLoading real SFA data from {sfa_path}...")
        try:
            real_frames, (Nx, Ny, Nz, ncols) = load_sfa_frames(sfa_path)
            print(f"Loaded {len(real_frames)} frames")
            analyze_data_characteristics(real_frames[0], "real frame 0")
        except Exception as e:
            print(f"Failed to load SFA: {e}")
            import traceback; traceback.print_exc()
            real_frames = None

    # --- Single frame benchmarks ---
    syn_single = run_single_frame_bench(syn_frames[0], "synthetic 64^3x6")

    if real_frames:
        real_single = run_single_frame_bench(real_frames[0], "real Cosserat 64^3x6")

    # --- Temporal benchmarks ---
    syn_temporal = run_temporal_bench(syn_frames, "synthetic")

    if real_frames:
        real_temporal = run_temporal_bench(real_frames, "real Cosserat")

    # --- Summary ---
    print(f"\n{'='*72}")
    print("SUMMARY")
    print(f"{'='*72}")

    # Find best lossless single frame
    all_lossless = [r for r in syn_single if r['lossless']]
    if all_lossless:
        best_ll = max(all_lossless, key=lambda r: r['ratio'])
        print(f"\nBest lossless (single frame, synthetic): {best_ll['name']} @ {best_ll['ratio']:.2f}x")

    # Find best lossy
    all_lossy = [r for r in syn_single if not r['lossless']]
    if all_lossy:
        best_lossy = max(all_lossy, key=lambda r: r['ratio'])
        print(f"Best lossy (single frame, synthetic): {best_lossy['name']} @ {best_lossy['ratio']:.2f}x, PSNR={best_lossy['psnr_dB']:.1f} dB")

    if real_frames:
        all_lossless_r = [r for r in real_single if r['lossless']]
        if all_lossless_r:
            best_ll_r = max(all_lossless_r, key=lambda r: r['ratio'])
            print(f"\nBest lossless (single frame, real): {best_ll_r['name']} @ {best_ll_r['ratio']:.2f}x")

        all_lossy_r = [r for r in real_single if not r['lossless']]
        if all_lossy_r:
            best_lossy_r = max(all_lossy_r, key=lambda r: r['ratio'])
            print(f"Best lossy (single frame, real): {best_lossy_r['name']} @ {best_lossy_r['ratio']:.2f}x, PSNR={best_lossy_r['psnr_dB']:.1f} dB")

    return syn_single, syn_temporal, real_single if real_frames else None, real_temporal if real_frames else None


if __name__ == '__main__':
    results = main()
