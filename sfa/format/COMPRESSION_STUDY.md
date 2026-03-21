# SFA Compression Study

Benchmark of compression pipelines for SFA (SCP Field Archive) volumetric
field data: float64 arrays from 3D physics simulations.

Test system: Linux x86-64, Python 3 + numpy 1.26 + pyzstd 0.19.
Benchmark script: `compression_bench.py` (in this directory).

---

## Data Characteristics

SFA frames store 3D volumetric fields: `N^3 x C` float64 values per frame
(N=grid size, C=column count). Typical: N=64, C=6, frame = 12.6 MB raw.

Tested on two datasets:

| Dataset | Description | Byte entropy | BSS byte entropy range |
|---------|-------------|-------------|----------------------|
| **Synthetic** | cos() bg (A=0.1) + Gaussian braid (A=0.8) + noise | 7.5 bits/byte | 1.0 -- 8.0 |
| **Real Cosserat** | Actual simulation, N=64, eta=0.5, 11 frames | 7.5 bits/byte | 1.0 -- 8.0 |

Key observation: the raw byte entropy is ~7.5 bits/byte (near maximum 8.0),
which is why plain zstd barely compresses. BSS (byte-stream split) separates
the high-entropy mantissa bytes from the low-entropy exponent/sign byte
(~1.0 bits/byte), enabling much better compression.

The real Cosserat test had 3 of 6 columns identically zero (theta fields
with mt=0.0, short T). The "non-zero only" column below isolates the actual
compression performance on live data.

---

## Single Frame Results

### Lossless Methods

Tested on a single 12.6 MB frame (64^3 x 6 columns, float64).

| Method | Ratio (synth) | Ratio (real, all) | Ratio (real, nz) | Comp MB/s | Decomp MB/s |
|--------|:---:|:---:|:---:|---:|---:|
| zstd-1 | 1.06 | 5.04 | 2.52 | 919 | 1059 |
| zstd-3 | 1.06 | 9.32 | 4.66 | 458 | 1030 |
| zstd-22 | 1.45 | 9.76 | 4.88 | 3.6 | 678 |
| **bss+zstd-1** | 1.53 | 10.15 | 5.08 | 440 | 519 |
| **bss+zstd-3** | 1.53 | 10.11 | 5.05 | 409 | 514 |
| bss+zstd-9 | 1.54 | 10.88 | 5.44 | 172 | 546 |
| bss+zstd-22 | 1.55 | 11.03 | 5.51 | 6.5 | 523 |
| xor+bss+zstd-3 | 1.53 | 10.19 | 5.10 | 246 | 31* |
| xor+bss+zstd-9 | 1.53 | 10.63 | 5.32 | 144 | 31* |
| **idelta+bss+zstd-3** | 1.52 | 10.87 | 5.43 | 261 | 40* |
| **idelta+bss+zstd-9** | 1.53 | 11.47 | 5.74 | 178 | 40* |
| idelta+bss+zstd-22 | 1.54 | 11.71 | 5.85 | 11 | 44* |

*XOR/idelta decompress speed limited by Python loop; C implementation
would match BSS speed (~500 MB/s).

**Key findings (lossless single frame):**

1. **BSS is essential.** On synthetic data, BSS+zstd-3 gets 1.53x vs zstd-3's
   1.06x -- a 44% improvement in compressed size. On real non-zero data: 5.05x
   vs 4.66x.

2. **Integer delta (idelta) adds ~5-7% ratio at higher zstd levels.** On real
   non-zero data: idelta+bss+zstd-9 = 5.74x vs bss+zstd-9 = 5.44x. The
   improvement is modest because spatial neighbors in the flattened array
   (ix varies slowest, iz fastest) are only neighbors along z -- not in x or y.

3. **XOR encoding is similar to idelta** but slightly worse ratio. Both need a
   sequential decode loop (not vectorizable as cumsum).

4. **zstd level barely matters above 3.** Level 3 vs 22 gains only 2-6% ratio
   but costs 60x in compression speed. Level 1 is only ~1% worse than 3.

5. **Float delta (not integer) is NOT bitwise lossless** due to FP rounding in
   cumsum. Must use integer-domain operations (XOR or int64 subtraction) for
   true lossless.

### Lossy Methods

| Method | Ratio (synth) | Ratio (real, nz) | Max Error | PSNR | Notes |
|--------|:---:|:---:|:---:|:---:|-------|
| f32+bss+zstd-3 | 4.81 | 17.61 | 3.0e-8 | 173 dB | ~7 decimal digits |
| f32+bss+zstd-9 | 4.92 | 19.68 | 3.0e-8 | 169 dB | |
| **f16+bss+zstd-3** | 89.3 | 56.9 | 2.4e-4 | 91-95 dB | ~3 decimal digits |
| bq8+zstd-3 | 163 | 119 | 1.9e-3 / 3.5e-3 | 65-70 dB | GGUF-style 8-bit |
| bq6+zstd-3 | 234 | 185 | 7.7e-3 / 1.4e-2 | 53-56 dB | GGUF-style 6-bit |
| bq4+zstd-3 | 261 | 216 | 3.3e-2 / 5.9e-2 | 41-45 dB | GGUF-style 4-bit |

**Key findings (lossy single frame):**

1. **f32 downcast is the best "safe lossy" option.** 3.0e-8 max error is
   negligible for any visualization or post-processing. Ratio improvement is
   dramatic: 4.8x (synth) to 17.6x (real). This is the clear choice for
   visualization codec.

2. **f16 is excellent for thumbnails/previews.** 56-89x compression with
   2.4e-4 max error. Fields in [-1,1] retain ~3 digits. Insufficient for
   computation but fine for visual inspection.

3. **Block quantization (GGUF-style) achieves extreme compression** (100-400x)
   but at significant precision cost. bq8 (8-bit) with max error ~3e-3 is
   acceptable for heatmap visualization. bq4 is too coarse for physics data.

---

## Temporal Sequence Results

Tested on 10-11 frame sequences (126-138 MB total raw).

| Method | Ratio (synth) | Ratio (real) | Lossless | Decomp MB/s |
|--------|:---:|:---:|:---:|---:|
| **indep bss+zstd-3** | 1.53 | 1.52 | YES | 509-530 |
| temporal-delta bss+zstd-3 | 1.45 | 1.52 | **NO** | 363-388 |
| temporal-xor bss+zstd-3 | 1.51 | 1.48 | YES | 392-396 |
| temporal-isub bss+zstd-3 | 1.52 | 1.50 | YES | 376-395 |
| temporal-delta f32+bss+zstd-3 | 3.84 | 3.93 | NO | 418-420 |
| concat bss+zstd-3 | 1.53 | 1.52 | YES | 335-343 |
| concat bss+zstd-9 | 1.54 | 1.53 | YES | 350-338 |

**Key findings (temporal):**

1. **Temporal differencing provides NO benefit** at this simulation cadence
   (dt=0.476, writing every 10 steps). The fields change too much between
   frames for temporal deltas to be small.

2. **Independent per-frame compression is optimal** when frames must be
   individually seekable (which SFA requires for random access). Temporal
   methods sacrifice seekability for negligible gain.

3. **Temporal-delta f32 (lossy)** does improve ratio to ~4x at the cost of
   accumulating FP32 rounding across frames (max error grows with frame count).

4. **Concatenation** (single zstd pass over all frames) gives no benefit
   vs independent -- the data doesn't have enough inter-frame redundancy
   at this temporal resolution.

**When temporal compression WOULD help:** If the simulation wrote every
timestep (dt=0.048), successive frames would differ by ~1/10th as much,
and temporal XOR/isub would see ~50% of bits unchanged, yielding ~2x extra.
This is relevant for checkpoint-level output (every step) vs diagnostic
output (every Nth step).

---

## GGUF Block Quantization Analysis

GGUF (used by llama.cpp) stores neural network weights using:
- **Block quantization**: divide tensor into blocks of 256 values
- **Per-block scale + min**: 2 x float32 = 8 bytes overhead per block
- **Quantized values**: 4-8 bits per value
- **K-quants**: hierarchical double-quantization (super-block scales quantized
  to 6-bit, further reducing metadata)

Techniques applicable to SFA:

1. **Block-adaptive precision**: Blocks in the low-amplitude background could
   use fewer bits than blocks in the high-gradient braid core. GGUF's
   per-block scale handles this automatically: small-range blocks get higher
   effective precision per bit.

2. **Super-block structure**: Grouping 16 blocks into a 4096-value super-block
   with a single FP16 scale factor reduces metadata from 3.1% to 0.2%.
   Worth adopting for bq8 and bq6.

3. **NOT directly applicable**: GGUF quantization is designed for weights that
   are read once and used many times. SFA data is written once and may be
   read for further computation, where bit-exact reproduction matters.
   Block quantization is appropriate only for visualization tiers.

---

## Recommendations

### Default Codec (lossless): `bss+zstd-3` (codec flag = 2)

**Keep the current default.** BSS+zstd-3 provides:
- 1.5x on dense float64 data, up to 10x when fields have limited dynamic range
- 400+ MB/s compress, 500+ MB/s decompress
- No decode loop (pure transpose + zstd)
- Well-tested, already implemented in sfa.h

The integer-delta variant (idelta+bss+zstd) gains only 5-7% ratio but adds a
sequential decode loop that halves decompression throughput in the current
Python prototype. In C the loop would be cheap, but the ratio gain is too small
to justify the added complexity as a default.

Increasing zstd level above 3 is not worthwhile: <2% ratio gain for 2-60x
slower compression.

### Visualization Codec (lossy OK): `f32+bss+zstd-3` (codec flag = 3, NEW)

Float32 downcast + BSS + zstd is the clear winner for visualization:
- 5-18x compression (data-dependent)
- Max error 3.0e-8 (7 decimal digits -- imperceptible)
- Faster than lossless (less data to process)
- Simple implementation: cast, BSS, zstd

This should be codec flag 3 in the SFA spec.

### Preview/Thumbnail Codec: `f16+bss+zstd-3` (codec flag = 4, NEW)

For quick previews or low-bandwidth streaming:
- 57-89x compression
- Max error 2.4e-4 (3 decimal digits)
- Sufficient for color-mapped volume rendering
- Extremely fast: 1.4+ GB/s compress

This should be codec flag 4.

### Block Quantization Codec: `bq8+zstd-3` (codec flag = 5, NEW)

For extreme compression (archival thumbnails, web previews):
- 119-163x compression
- Max error 1.9e-3 to 3.5e-3
- GGUF-inspired per-block (scale, min, 8-bit values)

### Temporal Compression: NOT RECOMMENDED for SFA

Temporal differencing conflicts with SFA's per-frame random-access design.
The compressed size gain is negligible at typical diagnostic output intervals.
If a temporal mode is ever needed, use **temporal XOR** (lossless) or
**temporal isub** (lossless), NOT float-domain delta (lossy).

---

## Updated Codec Flag Assignments

| Flag (bits 0-3) | Codec | Lossless | Description |
|:---:|---------|:---:|-------------|
| 0 | raw | yes | No compression |
| 1 | zstd | yes | zstd only (level 3) |
| 2 | bss+zstd | yes | Byte stream split + zstd (DEFAULT) |
| 3 | f32+bss+zstd | no | Float32 downcast + BSS + zstd |
| 4 | f16+bss+zstd | no | Float16 downcast + BSS + zstd |
| 5 | bq8+zstd | no | Block-quant 8-bit + zstd |
| 6 | idelta+bss+zstd | yes | Int64 spatial delta + BSS + zstd |
| 7-15 | reserved | -- | Future use |

Note: codecs 3-5 are lossy and store data at reduced precision. The
original dtype in CDEF remains f64; the codec flag tells the reader to
expect reduced-precision data and reconstruct to f64 on read.

---

## Raw Benchmark Data

### Single Frame, Synthetic (64^3 x 6 fields, 12.58 MB raw)

```
Method                                      Ratio     Comp   Decomp  Loss     MaxErr     PSNR
                                              (x)   (MB/s)   (MB/s)                      (dB)
-----------------------------------------------------------------------------------------------
zstd-1                                       1.06    919.0   1059.0    no          -      inf
zstd-3                                       1.06    458.0   1030.0    no          -      inf
zstd-22                                      1.45      3.6    678.0    no          -      inf
bss+zstd-1                                   1.53    440.0    519.0    no          -      inf
bss+zstd-3                                   1.53    409.0    514.0    no          -      inf
bss+zstd-9                                   1.54    172.0    546.0    no          -      inf
bss+zstd-22                                  1.55      6.5    523.0    no          -      inf
xor+bss+zstd-3                               1.53    246.0     31.0    no          -      inf
idelta+bss+zstd-3                            1.52    261.0     40.0    no          -      inf
idelta+bss+zstd-9                            1.53    141.0     39.0    no          -      inf
f32+bss+zstd-3                               4.81    686.0    644.0   YES   2.97e-08    173.1
f16+bss+zstd-3                              89.34   1435.0    678.0   YES   2.42e-04     95.1
bq8+zstd-3                                 163.18    294.0    893.0   YES   1.91e-03     70.0
bq4+zstd-3                                 260.92    301.0    880.0   YES   3.25e-02     45.7
```

### Single Frame, Real Cosserat Non-Zero Columns (64^3 x 3 fields, 6.29 MB raw)

```
Method                                      Ratio     Comp   Decomp  Loss     MaxErr     PSNR
                                              (x)   (MB/s)   (MB/s)                      (dB)
-----------------------------------------------------------------------------------------------
zstd-1                                       2.52    698.0   1451.0    no          -      inf
zstd-3                                       4.66    468.0   1537.0    no          -      inf
bss+zstd-1                                   5.08    894.0    889.0    no          -      inf
bss+zstd-3                                   5.05    650.0    955.0    no          -      inf
bss+zstd-9                                   5.44    189.0    936.0    no          -      inf
idelta+bss+zstd-3                            5.43    549.0     46.0    no          -      inf
idelta+bss+zstd-9                            5.74    215.0     46.0    no          -      inf
f32+bss+zstd-3                              17.61   1329.0   1076.0   YES   2.98e-08    168.6
f16+bss+zstd-3                              56.90   1591.0   1168.0   YES   2.44e-04     90.6
bq8+zstd-3                                 118.70    955.0   1733.0   YES   3.49e-03     65.4
```

### Temporal Sequence, Real Cosserat (11 frames, 138.4 MB raw)

```
Method                                      Ratio     Comp   Decomp  Loss     MaxErr     PSNR
                                              (x)   (MB/s)   (MB/s)                      (dB)
-----------------------------------------------------------------------------------------------
indep bss+zstd-3 (11f)                       1.52    371.0    530.0    no          -      inf
temporal-xor bss+zstd-3 (11f)                1.48    353.0    396.0    no          -      inf
temporal-isub bss+zstd-3 (11f)               1.50    341.0    395.0    no          -      inf
temporal-delta f32+bss+zstd-3 (11f)          3.93    495.0    420.0   YES   2.81e-07    164.2
concat bss+zstd-3 (11f)                      1.52    247.0    334.0    no          -      inf
```

---

## Methodology

- Each measurement averaged over 3 iterations
- Compression ratio = raw_bytes / compressed_bytes (higher is better)
- Speed measured on raw data size (not compressed)
- PSNR = 20*log10(max|signal|/RMS_error)
- "Real" data generated by: `cosserat_sfa -N 64 -L 15 -T 10 -eta 0.5 -mt 0.0`
- "Synthetic" data: cos() background + Gaussian braid + 1e-6 noise
- BSS implemented as numpy reshape+transpose (matches sfa.h logic)
- XOR/idelta decode uses Python loop (C would be ~10x faster)

---

## Extended Study: BSS-only and Alternative Codecs

Extended benchmark (`extended_bench.py`) testing BSS without compression (entropy
analysis) and six alternative compressor families applied with and without BSS
preprocessing. All tests on the same synthetic data (64^3 x 6 float64, 12.58 MB).

### BSS Entropy Analysis

BSS (byte stream split) transposes the float64 bytes so that all byte-0 values
are contiguous, all byte-1 values are contiguous, etc. On little-endian x86:
byte 7 is the MSB (sign + exponent high bits), byte 0 is the LSB (mantissa noise).

**BSS does NOT change the data size.** It is a pure transpose (ratio = 1.0).
Its value is in separating low-entropy streams from high-entropy streams so that
a downstream compressor can exploit the structure.

BSS transpose throughput: **1501 MB/s encode, 1120 MB/s decode** (numpy
reshape+transpose, dominated by memory bandwidth).

Per-stream entropy on the synthetic test data (all 6 columns concatenated):

| Stream | Byte position | Entropy (bits/byte) | Theo. ratio | Compressibility |
|:------:|---------------|--------------------:|:-----------:|-----------------|
| 7 | MSB: sign + exponent[10:4] | 1.00 | 8.05x | **excellent** |
| 6 | exponent[3:0] + mantissa[51:48] | 5.26 | 1.52x | moderate |
| 5 | mantissa[47:40] | 7.60 | 1.05x | marginal |
| 4 | mantissa[39:32] | 8.00 | 1.00x | none |
| 3 | mantissa[31:24] | 8.00 | 1.00x | none |
| 2 | mantissa[23:16] | 8.00 | 1.00x | none |
| 1 | mantissa[15:8] | 8.00 | 1.00x | none |
| 0 | LSB: mantissa[7:0] | 8.00 | 1.00x | none |

**Per-stream Shannon entropy theoretical limit: ratio 1.188x** (sum of
per-stream optimal encoding). Since BSS+zstd-22 achieves 1.55x, zstd is
exploiting sequential correlations (LZ77 matches) BEYOND what single-byte
entropy captures. The per-byte entropy underestimates the true compressibility
because adjacent values in each stream have correlated patterns that LZ77 can
match.

### Alternative Codec Comparison

Each codec tested both on raw (interleaved) data and after BSS preprocessing.
5 iterations per measurement.

#### Without BSS (raw data)

| Codec | Ratio | Compress MB/s | Decompress MB/s | Notes |
|-------|------:|:-------------:|:---------------:|-------|
| lz4 (default) | 1.00 | 1751 | 2611 | cannot compress random-looking bytes |
| snappy | 1.00 | 905 | 1128 | same problem |
| zstd-1 | 1.06 | 841 | 914 | |
| zstd-3 | 1.06 | 599 | 1008 | |
| brotli-1 | 1.06 | 322 | 198 | |
| lz4hc-9 | 1.08 | 55 | 2295 | slow compress, fast decompress |
| zstd-9 | 1.07 | 92 | 895 | |
| gzip-6 | 1.17 | 23 | 247 | |
| brotli-9 | 1.43 | 3 | 341 | needs high quality to find patterns |
| zstd-22 | 1.45 | 4 | 596 | |
| lzma-1 | 1.46 | 6 | 29 | |
| lzma-6 | 1.49 | 2 | 29 | highest ratio without BSS |

Without BSS, the float64 byte stream looks nearly random (7.5 bits/byte entropy).
Fast codecs (lz4, snappy) cannot compress it at all. Only codecs with large
dictionary windows and intensive matching (zstd-22, lzma, brotli-6+) achieve
meaningful compression, but at enormous speed cost.

#### With BSS Preprocessing

| Pipeline | Ratio | Compress MB/s | Decompress MB/s | Notes |
|----------|------:|:-------------:|:---------------:|-------|
| BSS + snappy | 1.46 | 295 | 465 | fastest compress |
| BSS + lz4 | 1.49 | 570 | 799 | **best speed** |
| BSS + brotli-1 | 1.50 | 587 | 681 | |
| BSS + gzip-1 | 1.51 | 62 | 394 | |
| BSS + brotli-4 | 1.52 | 222 | 656 | |
| BSS + idelta-z + zstd-3 | 1.52 | 261 | 40 | (from previous study) |
| BSS + zstd-1 | 1.53 | 494 | 670 | |
| BSS + zstd-3 | **1.53** | **552** | **676** | **current SFA default** |
| BSS + lz4hc-9 | 1.53 | 46 | 942 | **best decompress speed** |
| BSS + gzip-6 | 1.54 | 55 | 457 | |
| BSS + zstd-9 | 1.54 | 248 | 847 | |
| BSS + brotli-6 | 1.54 | 109 | 640 | |
| BSS + zstd-22 | 1.55 | 7 | 584 | |
| BSS + brotli-9 | 1.55 | 49 | 656 | |
| BSS + lzma-1 | 1.55 | 9 | 424 | |
| BSS + lzma-3 | 1.56 | 4 | 424 | |
| BSS + brotli-11 | 1.56 | 1 | 333 | |
| BSS + lzma-6 | 1.56 | 4 | 489 | highest ratio |

#### BSS Improvement Factor

BSS makes ALL codecs dramatically better on float64 data:

| Codec family | Ratio without BSS | Ratio with BSS | BSS factor |
|-------------|------------------:|---------------:|-----------:|
| lz4 (fast) | 1.00 | 1.49 | **1.49x** |
| snappy | 1.00 | 1.46 | **1.46x** |
| zstd-3 | 1.06 | 1.53 | 1.44x |
| gzip-6 | 1.17 | 1.54 | 1.32x |
| brotli-6 | 1.43 | 1.54 | 1.08x |
| lzma-6 | 1.49 | 1.56 | 1.05x |

BSS provides the largest improvement for fast codecs (lz4, snappy) that cannot
find patterns in random-looking interleaved bytes. For slow codecs that already
use large windows (lzma, brotli-9+, zstd-22), BSS still helps but by a smaller
margin since these codecs can partially discover the byte-position structure
themselves.

### Key Findings

1. **All codecs converge to ratio ~1.5-1.56 with BSS.** The compressible content
   (streams 6 and 7) is the same regardless of codec. The difference is how
   well each codec handles the incompressible streams (0-5). The spread between
   worst (BSS+snappy=1.46) and best (BSS+lzma-6=1.56) is only 7%.

2. **BSS+zstd-3 remains the optimal choice.** It achieves 1.53x ratio at
   552 MB/s compress / 676 MB/s decompress. No alternative codec beats this
   ratio-vs-speed tradeoff:
   - BSS+lz4 is faster (570/799 MB/s) but worse ratio (1.49x)
   - BSS+lzma-6 has better ratio (1.56x) but 140x slower compression
   - BSS+lz4hc-9 has the fastest decompression (942 MB/s) but only 46 MB/s compress

3. **If maximum decompression speed matters most:** BSS+lz4hc-9 achieves
   942 MB/s decompress (40% faster than zstd-3) with ratio 1.53x, but
   compress speed drops to 46 MB/s. BSS+lz4 (default) gives 799 MB/s
   decompress at 570 MB/s compress — a good "read-heavy" choice.

4. **gzip and brotli offer no advantage over zstd** for this data. gzip is
   universally worse (slower and lower ratio). brotli needs quality >= 6 to
   match zstd-3's ratio, at which point it's much slower.

5. **lzma achieves the highest ratio (1.56x)** but is impractical: 4 MB/s
   compress, and the ratio gain over zstd-3 is only 2%.

### Pareto-Optimal Pipelines

Pipelines not dominated in both ratio and decompress speed:

| Pipeline | Ratio | Compress MB/s | Decompress MB/s | Use case |
|----------|------:|:-------------:|:---------------:|----------|
| BSS + lz4hc-9 | 1.53 | 46 | **942** | read-heavy archival |
| BSS + zstd-9 | 1.54 | 248 | 847 | balanced |
| BSS + zstd-3 | 1.53 | **552** | 676 | **general purpose (default)** |
| BSS + lzma-6 | **1.56** | 4 | 489 | maximum ratio |

**Recommendation: no change.** BSS+zstd-3 remains the best default. For a
future "fast decompress" codec flag, BSS+lz4hc-9 or BSS+zstd-9 would be
the candidates to consider.

---

## Spatial Ordering Study

Does reordering the 3D grid data along space-filling curves improve compression?
The default SFA layout stores data in C row-major order (ix slowest, iz fastest),
which gives perfect locality along z but poor locality along x and y.

Tested spatial orderings and transforms, all with BSS+zstd-3 backend.
Benchmark script: `spatial_bench.py`.

### Spatial Correlation Analysis

For the synthetic test data (cos background + Gaussian braid), the mean absolute
difference between adjacent grid points varies by axis:

| Axis | Direction | Mean |delta| | Relative to z |
|------|-----------|---------------:|:-------------:|
| z | fastest-varying (row-major) | 0.00655 | 1.0x |
| y | middle | 0.00022 | **0.03x** |
| x | slowest-varying | 0.00022 | **0.03x** |

The data has 30x LESS variation along x and y than along z. This is because the
background field cos(kZ) oscillates primarily in z, making z-adjacent values
differ significantly while x/y-adjacent values (same z) are nearly identical.
This is an artifact of this specific test data — real simulations may have more
isotropic correlation structure.

### Results

| Pipeline | Ratio | Compress MB/s | Decompress MB/s | Lossless | Notes |
|----------|------:|:-------------:|:---------------:|:--------:|-------|
| **linear + BSS + zstd-3** | **1.533** | **497** | **614** | yes | **current SFA default** |
| plane-x + BSS + zstd-3 | 1.538 | 379 | 358 | yes | yz-planes stacked in x |
| idelta-z + BSS + zstd-3 | 1.522 | 484 | 424 | yes | int64 delta along z |
| idelta-z + BSS + zstd-22 | 1.542 | 5 | 447 | yes | z-delta, max zstd |
| idelta-y + BSS + zstd-3 | 1.423 | 374 | 397 | yes | int64 delta along y |
| idelta-x + BSS + zstd-3 | 1.425 | 364 | 396 | yes | int64 delta along x |
| morton + BSS + zstd-3 | 1.493 | 432 | 432 | yes | Z-order space-filling curve |
| hilbert + BSS + zstd-3 | 1.474 | 468 | 439 | yes | Hilbert space-filling curve |
| morton + idelta-z + BSS + zstd-3 | 1.369 | 258 | 349 | yes | Morton + z-delta combined |

### Key Findings

1. **Space-filling curves (Morton, Hilbert) HURT compression** on this data.
   Morton: 1.493x (2.6% worse), Hilbert: 1.474x (3.8% worse). The curves
   break the z-contiguity that the compressor exploits. Since this data has
   strong z-axis correlation (cos(kZ) background), any reordering that
   disrupts z-runs makes compression harder.

2. **Plane-major reordering has negligible effect** (1.538x vs 1.533x, +0.3%).
   For cubic grids, transposing axes barely matters because the data has the
   same number of contiguous runs in each direction.

3. **3D delta encoding along z (fastest axis) slightly hurts** (1.522x vs 1.533x,
   -0.7%). The delta values have higher byte entropy than the original values
   because the z-oscillation makes successive differences comparable in magnitude
   to the values themselves (delta_z mean = 0.00655 vs value range ~0.9). At
   high zstd levels (22), z-delta + BSS + zstd-22 reaches 1.542x, slightly
   beating plain BSS + zstd-22's 1.55x — so delta helps at extreme compression
   levels where zstd can exploit the reduced variance.

4. **Delta along x or y is much worse** (1.42x, -7%) because these take
   differences between values that are N^2 or N apart in memory, destroying
   the z-axis sequential correlation that zstd relies on.

5. **Combining Morton + z-delta is worst of all** (1.369x, -10.7%). Morton
   scrambles the z-ordering, then z-delta computes differences between
   Morton-adjacent (NOT z-adjacent) values, producing large noisy deltas.

### When Spatial Reordering WOULD Help

The test data has strongly anisotropic correlations (z-axis dominated). For data
with isotropic 3D correlations (e.g., a spherical soliton centered at origin),
Morton or Hilbert ordering could improve compression by ensuring ALL neighbor
pairs (not just z-neighbors) appear close in the byte stream. This would matter
most for:

- AMR patches with non-cubic aspect ratios (e.g., 4x4x256 slabs)
- Data with no preferred axis (spherical symmetry)
- Higher-order prediction filters that need multi-axis neighbors

For the current SFA use case (diagnostic output of 3D simulations), the
default linear row-major ordering is optimal. No spatial reordering is
recommended.

**Recommendation: keep linear (row-major) ordering.** The overhead of
computing and storing the permutation indices outweighs any marginal
compression benefit, and the benefit is data-dependent anyway.
