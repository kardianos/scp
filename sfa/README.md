# SFA — SCP Field Archive

Single-file compressed container for volumetric field simulation data.

## Directory Structure

```
sfa/
  format/
    SFA_SPEC.md     — Full specification (v3.1 with AMR support)
    sfa.h           — Single-header C library (write + read)
    sfa_test.c      — Roundtrip test
  viewer/
    volview.c       — Interactive SFA volume viewer (SDL2, raymarcher)
    volview_bin.c   — Legacy .bin file viewer (single snapshot)
```

## Build

```bash
# Library test
cd format && gcc -O2 -o sfa_test sfa_test.c -lzstd -lm && ./sfa_test

# SFA volume viewer (frame stepping, animation)
cd viewer && gcc -O3 -fopenmp -o volview volview.c \
    -I../format $(sdl2-config --cflags --libs) -lzstd -lm

# Legacy .bin viewer (single file)
cd viewer && gcc -O3 -fopenmp -o volview_bin volview_bin.c \
    $(sdl2-config --cflags --libs) -lm
```

## Usage

```bash
# View SFA archive with frame stepping
./volview simulation.sfa

# Controls:
#   Mouse drag: rotate
#   Scroll/Up/Down: zoom
#   Left/Right: prev/next frame
#   Space: play/pause animation
#   1/2/3: toggle RED(bound)/GREEN(fabric)/BLUE(angle)
#   O/P: decrease/increase opacity
#   +/-: brightness
#   S: screenshot
#   R: reset view
#   Q: quit
```

## Features

- Compressed (BSS + zstd, ~1.2-5× ratio on float64 field data)
- Seekable (two-level jump tables, O(2) reads for any frame)
- Columnar (per-column dtype, self-describing schema)
- Multi-resolution AMR (GRID patches at different scales, FGRP grouping)
- Mmap-friendly (fixed layout header + index)
- Appendable (frames written during simulation)
- Dependencies: libzstd only (viewer also needs SDL2)
