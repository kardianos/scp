#!/usr/bin/env python3
"""existence_spectrum.py — does the particle's EXISTENCE blink, or only its
components?

Reads a high-cadence diag.tsv and Fourier-analyzes:
  phi_max  — max over the box of the REAL phi component: for a rotating
             (charged) ball u = f cos(wt), this blinks at w even though the
             object is static.
  s_max    — max of the phase-invariant potential source s = prod|Phi_a|^2:
             this is the EXISTENCE density. Static for a charged ball,
             deeply modulated (through ~0) for a +/-w breather.
  Q_core, r_core — sanity.

Outputs modulation depth (max-min)/(max+min) over the analysis window and the
top spectral peaks of each series (plain periodogram of the detrended,
Hann-windowed series).

Usage: existence_spectrum.py diag.tsv [t_min] [t_max]
"""
import sys
import numpy as np


def peaks(t, y, nmax=5):
    dt = np.median(np.diff(t))
    y = y - np.mean(y)
    w = np.hanning(len(y))
    Y = np.abs(np.fft.rfft(y * w))
    f = np.fft.rfftfreq(len(y), dt) * 2 * np.pi  # angular frequency
    idx = np.argsort(Y[1:])[::-1] + 1
    out = []
    for i in idx:
        if all(abs(f[i] - fo) > 0.15 for fo, _ in out):
            out.append((f[i], Y[i]))
        if len(out) >= nmax:
            break
    return out, f, Y


def main(path, tmin=20.0, tmax=None):
    d = np.genfromtxt(path, names=True)
    t = d["t"]
    m = t >= tmin
    if tmax:
        m &= t <= tmax
    t = t[m]
    print(f"# {path}: {m.sum()} samples, t in [{t[0]:.1f},{t[-1]:.1f}], "
          f"dt={np.median(np.diff(t)):.3f}")
    print(f"# omega ref: ball internal w=1.42 -> expect phi_max line at w, "
          f"breather s line at 2w=2.84")
    for col in ("phi_max", "s_max", "Q_core", "r_core"):
        if col not in d.dtype.names:
            continue
        y = d[col][m]
        depth = (y.max() - y.min()) / (y.max() + y.min() + 1e-300)
        pk, _, _ = peaks(t, y)
        pks = "  ".join(f"w={f:.3f}(A={a:.3g})" for f, a in pk[:3])
        print(f"{col:8s} mean={np.mean(y):+.4e}  depth={depth:8.2%}  peaks: {pks}")


if __name__ == "__main__":
    main(sys.argv[1], *[float(a) for a in sys.argv[2:]])
