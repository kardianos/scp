#!/usr/bin/env python3
"""
analyze_candidates.py — analyze evolutionary search candidate .npz files.

Usage:
    python analyze_candidates.py [--dir candidates/] [--outdir analysis/]

Physics: V(P) = (mu/2)*P^2/(1+kappa*P^2)
    where P = phi_0 * phi_1 * phi_2
    mu = -41.345, kappa = 50.0, m^2 = 2.25
"""

import argparse
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import ndimage


# ── physics constants ────────────────────────────────────────────────────────

MU     = -41.345
KAPPA  = 50.0
M2     = 2.25          # mass^2 for E_mass term


# ── helpers ──────────────────────────────────────────────────────────────────

def load_candidate(path):
    """Return dict with phi (3,N,N,N float64), vel, N (int), L (float), dx."""
    d = np.load(path)
    phi = d["phi"].astype(np.float64)
    vel = d["vel"].astype(np.float64)
    N   = int(d["N"])
    L   = float(d["L"])
    dx  = 2.0 * L / (N - 1)
    return dict(phi=phi, vel=vel, N=N, L=L, dx=dx)


def compute_P(phi):
    """Product field P = phi_0 * phi_1 * phi_2."""
    return phi[0] * phi[1] * phi[2]


def compute_V(P):
    """Potential density V(P) = (mu/2)*P^2 / (1 + kappa*P^2)."""
    return (MU / 2.0) * P**2 / (1.0 + KAPPA * P**2)


def compute_gradient_energy_density(phi, dx):
    """
    E_grad density = 0.5 * sum_a |grad phi_a|^2
    Uses central differences with periodic roll (adequate for localized fields
    that vanish at the boundary).
    """
    density = np.zeros(phi.shape[1:], dtype=np.float64)
    for a in range(3):
        for axis in range(3):
            # shift along axis+1 because phi has leading component dimension
            gf = (np.roll(phi[a], -1, axis=axis) - np.roll(phi[a], 1, axis=axis)) / (2.0 * dx)
            density += 0.5 * gf**2
    return density


def make_grid(N, L):
    """Return (X, Y, Z) meshgrid and r array, all shape (N,N,N)."""
    x = np.linspace(-L, L, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    r = np.sqrt(X**2 + Y**2 + Z**2)
    return X, Y, Z, r


def parse_generation(filename):
    """Extract generation number from filename like 'ckpt_T10_pop128_gen100.npz'."""
    m = re.search(r"gen(\d+)", filename)
    return int(m.group(1)) if m else -1


def parse_series(filename):
    """
    Return a short series label for grouping, e.g. 'T10_pop128', 'T30_pop64'.
    Falls back to 'unknown'.
    """
    m = re.search(r"(T\d+_pop\d+)", filename)
    if m:
        return m.group(1)
    m = re.search(r"(T\d+)", filename)
    if m:
        return m.group(1)
    return "unknown"


# ── per-candidate analysis ───────────────────────────────────────────────────

def analyze_candidate(path):
    """
    Full analysis of one .npz candidate.
    Returns a flat dict of scalars + arrays for profiles.
    """
    fname = os.path.basename(path)
    c = load_candidate(path)
    phi, vel, N, L, dx = c["phi"], c["vel"], c["N"], c["L"], c["dx"]
    dV = dx**3

    # ── field products & densities ───────────────────────────────────────────
    P     = compute_P(phi)
    V_den = compute_V(P)           # potential energy density
    grad_den = compute_gradient_energy_density(phi, dx)
    mass_den = 0.5 * M2 * (phi**2).sum(axis=0)
    kin_den  = 0.5 * (vel**2).sum(axis=0)

    # ── integrated energies ──────────────────────────────────────────────────
    E_pot  = V_den.sum()  * dV
    E_grad = grad_den.sum() * dV
    E_mass = mass_den.sum() * dV
    E_kin  = kin_den.sum()  * dV
    E_total = E_pot + E_grad + E_mass + E_kin

    ratio_pot  = E_pot  / E_total  if E_total != 0 else 0.0
    ratio_gm   = E_grad / E_mass  if E_mass  != 0 else 0.0

    # ── grid ─────────────────────────────────────────────────────────────────
    X, Y, Z, r = make_grid(N, L)

    # ── centroid of |V| binding ──────────────────────────────────────────────
    weight  = np.abs(V_den)
    total_w = weight.sum()
    if total_w > 0:
        cx = (weight * X).sum() / total_w
        cy = (weight * Y).sum() / total_w
        cz = (weight * Z).sum() / total_w
    else:
        cx = cy = cz = 0.0

    # ── R_rms and R_half ─────────────────────────────────────────────────────
    r2_c = (X - cx)**2 + (Y - cy)**2 + (Z - cz)**2
    r_c  = np.sqrt(r2_c)

    if total_w > 0:
        R_rms = np.sqrt((weight * r2_c).sum() / total_w)
        # R_half: radius enclosing 50% of |E_pot|
        r_flat = r_c.ravel()
        w_flat = weight.ravel()
        sort_idx = np.argsort(r_flat)
        cumw     = np.cumsum(w_flat[sort_idx])
        half_idx = np.searchsorted(cumw, 0.5 * cumw[-1])
        half_idx = min(half_idx, len(r_flat) - 1)
        R_half   = r_flat[sort_idx[half_idx]]
    else:
        R_rms = R_half = 0.0

    # ── inertia tensor (|V|-weighted) ────────────────────────────────────────
    if total_w > 0:
        dx_ = X - cx; dy_ = Y - cy; dz_ = Z - cz
        Ixx = (weight * (dy_**2 + dz_**2)).sum() / total_w
        Iyy = (weight * (dx_**2 + dz_**2)).sum() / total_w
        Izz = (weight * (dx_**2 + dy_**2)).sum() / total_w
        Ixy = -(weight * dx_ * dy_).sum() / total_w
        Ixz = -(weight * dx_ * dz_).sum() / total_w
        Iyz = -(weight * dy_ * dz_).sum() / total_w
        I_tensor = np.array([[Ixx, Ixy, Ixz],
                              [Ixy, Iyy, Iyz],
                              [Ixz, Iyz, Izz]])
        evals = np.linalg.eigvalsh(I_tensor)
        evals = np.maximum(evals, 0.0)   # numerical noise guard
        aspect = (evals[-1] / evals[0]) if evals[0] > 1e-12 else np.inf
    else:
        evals  = np.zeros(3)
        aspect = 1.0

    # ── phase structure ───────────────────────────────────────────────────────
    absP      = np.abs(P)
    P_max     = absP.max()
    core_thr  = 0.1 * P_max

    core_mask = absP > core_thr
    n_core    = core_mask.sum()

    if n_core > 0:
        P_core    = P[core_mask]
        phi0_c    = phi[0][core_mask]
        phi1_c    = phi[1][core_mask]
        phi2_c    = phi[2][core_mask]
        frac_Ppos = (P_core > 0).sum() / n_core
        P_mean    = P_core.mean()
        P_std     = P_core.std()
        phase01_std = np.arctan2(phi1_c, phi0_c).std()
        phase02_std = np.arctan2(phi2_c, phi0_c).std()
        # Phase coherence: low std = coherent
        phase_coherence = 1.0 - min(phase01_std / np.pi, 1.0)
    else:
        frac_Ppos = 0.5
        P_mean = P_std = 0.0
        phase01_std = phase02_std = np.pi
        phase_coherence = 0.0

    # ── radial profiles ───────────────────────────────────────────────────────
    n_shells    = 20
    shell_edges = np.linspace(0.0, L, n_shells + 1)
    r_centers   = 0.5 * (shell_edges[:-1] + shell_edges[1:])

    phi2_sum  = np.zeros(n_shells)
    absP_sum  = np.zeros(n_shells)
    V_sum     = np.zeros(n_shells)
    counts    = np.zeros(n_shells)

    phi2_all  = (phi**2).sum(axis=0)
    shell_idx = (np.searchsorted(shell_edges, r.ravel(), side="right") - 1)
    shell_idx = np.clip(shell_idx, 0, n_shells - 1)

    np.add.at(phi2_sum, shell_idx, phi2_all.ravel())
    np.add.at(absP_sum, shell_idx, absP.ravel())
    np.add.at(V_sum,    shell_idx, V_den.ravel())
    np.add.at(counts,   shell_idx, 1.0)
    counts_safe = np.maximum(counts, 1)

    radial_phi2 = phi2_sum  / counts_safe
    radial_absP = absP_sum  / counts_safe
    radial_V    = V_sum     / counts_safe

    # ── cluster analysis ──────────────────────────────────────────────────────
    labeled, n_clusters = ndimage.label(core_mask)
    if n_clusters > 0:
        sizes = np.array(ndimage.sum(core_mask, labeled, range(1, n_clusters + 1)))
        largest_label = int(np.argmax(sizes)) + 1
        largest_vol   = float(sizes[np.argmax(sizes)]) * dV
        largest_mass  = float(np.abs(V_den[labeled == largest_label]).sum()) * dV
        lc_idx        = np.where(labeled == largest_label)
        largest_cx    = X[lc_idx].mean()
        largest_cy    = Y[lc_idx].mean()
        largest_cz    = Z[lc_idx].mean()
    else:
        n_clusters   = 0
        largest_vol  = 0.0
        largest_mass = 0.0
        largest_cx   = largest_cy = largest_cz = 0.0

    # ── generation / series metadata ─────────────────────────────────────────
    gen    = parse_generation(fname)
    series = parse_series(fname)

    return dict(
        # metadata
        fname        = fname,
        path         = path,
        gen          = gen,
        series       = series,
        N            = N,
        L            = L,
        dx           = dx,
        # energies
        E_pot        = E_pot,
        E_grad       = E_grad,
        E_mass       = E_mass,
        E_kin        = E_kin,
        E_total      = E_total,
        ratio_pot    = ratio_pot,
        ratio_gm     = ratio_gm,
        # spatial
        cx=cx, cy=cy, cz=cz,
        R_rms        = R_rms,
        R_half       = R_half,
        inertia_evals= evals,
        aspect_ratio = aspect,
        # phase
        P_max        = P_max,
        frac_Ppos    = frac_Ppos,
        P_mean       = P_mean,
        P_std        = P_std,
        phase_coh    = phase_coherence,
        phase01_std  = phase01_std,
        phase02_std  = phase02_std,
        # clusters
        n_clusters   = n_clusters,
        largest_vol  = largest_vol,
        largest_mass = largest_mass,
        largest_cx   = largest_cx,
        largest_cy   = largest_cy,
        largest_cz   = largest_cz,
        # radial profiles (arrays)
        r_centers    = r_centers,
        radial_phi2  = radial_phi2,
        radial_absP  = radial_absP,
        radial_V     = radial_V,
        # raw fields (for slices)
        _phi         = phi,
        _vel         = vel,
        _P           = P,
        _V_den       = V_den,
    )


# ── output helpers ────────────────────────────────────────────────────────────

def print_summary_table(results):
    """Print a formatted table to stdout."""
    # Sort by E_pot (most negative first)
    results = sorted(results, key=lambda r: r["E_pot"])

    header = (
        f"{'Filename':<35s}  {'gen':>5s}  "
        f"{'E_pot':>10s}  {'E_grad':>9s}  {'E_mass':>9s}  {'E_kin':>8s}  "
        f"{'E_pot/E_tot':>11s}  {'R_rms':>7s}  {'R_half':>7s}  "
        f"{'aspect':>7s}  {'N_clust':>7s}  {'PhCoh':>6s}"
    )
    sep = "-" * len(header)
    print(sep)
    print("CANDIDATE ANALYSIS SUMMARY")
    print(f"mu={MU}, kappa={KAPPA}, m^2={M2}")
    print(sep)
    print(header)
    print(sep)
    for r in results:
        print(
            f"{r['fname']:<35s}  {r['gen']:>5d}  "
            f"{r['E_pot']:>10.3f}  {r['E_grad']:>9.3f}  {r['E_mass']:>9.3f}  "
            f"{r['E_kin']:>8.3f}  {r['ratio_pot']:>11.4f}  "
            f"{r['R_rms']:>7.3f}  {r['R_half']:>7.3f}  "
            f"{r['aspect_ratio']:>7.3f}  {r['n_clusters']:>7d}  "
            f"{r['phase_coh']:>6.3f}"
        )
    print(sep)
    print()

    # Extra detail for top candidate
    top = results[0]
    print(f"TOP CANDIDATE: {top['fname']}")
    print(f"  N={top['N']}, L={top['L']}, dx={top['dx']:.4f}")
    print(f"  E_pot  = {top['E_pot']:.4f}   (binding energy — negative is good)")
    print(f"  E_grad = {top['E_grad']:.4f}")
    print(f"  E_mass = {top['E_mass']:.4f}")
    print(f"  E_kin  = {top['E_kin']:.4f}")
    print(f"  E_total= {top['E_total']:.4f}")
    print(f"  E_pot/E_total = {top['ratio_pot']:.4f}")
    print(f"  E_grad/E_mass = {top['ratio_gm']:.4f}")
    print(f"  Centroid = ({top['cx']:.3f}, {top['cy']:.3f}, {top['cz']:.3f})")
    print(f"  R_rms  = {top['R_rms']:.3f}   R_half = {top['R_half']:.3f}")
    print(f"  Inertia eigenvalues = {top['inertia_evals']}")
    print(f"  Aspect ratio (I_max/I_min) = {top['aspect_ratio']:.4f}")
    print(f"  P_max  = {top['P_max']:.4f}")
    print(f"  Fraction P>0 in core = {top['frac_Ppos']:.3f}  (0.5 = random, near 0 or 1 = coherent)")
    print(f"  P_mean = {top['P_mean']:.4f}  P_std = {top['P_std']:.4f}")
    print(f"  Phase coherence  = {top['phase_coh']:.3f}  (0=random, 1=locked)")
    print(f"  phase01_std = {top['phase01_std']:.3f} rad,  phase02_std = {top['phase02_std']:.3f} rad")
    print(f"  N_clusters = {top['n_clusters']}")
    print(f"  Largest cluster: vol={top['largest_vol']:.3f}, |E_pot|={top['largest_mass']:.3f}")
    print(f"  Largest cluster centroid = ({top['largest_cx']:.3f}, {top['largest_cy']:.3f}, {top['largest_cz']:.3f})")
    print()


def write_summary_tsv(results, outpath):
    """Write one TSV row per candidate."""
    results_sorted = sorted(results, key=lambda r: r["E_pot"])
    cols = [
        "fname", "series", "gen", "N", "L", "dx",
        "E_pot", "E_grad", "E_mass", "E_kin", "E_total",
        "ratio_pot", "ratio_gm",
        "cx", "cy", "cz", "R_rms", "R_half", "aspect_ratio",
        "P_max", "frac_Ppos", "P_mean", "P_std", "phase_coh",
        "phase01_std", "phase02_std",
        "n_clusters", "largest_vol", "largest_mass",
        "largest_cx", "largest_cy", "largest_cz",
    ]
    with open(outpath, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in results_sorted:
            row = "\t".join(
                str(r[c]) if not isinstance(r[c], float) else f"{r[c]:.6g}"
                for c in cols
            )
            fh.write(row + "\n")
    print(f"Wrote {outpath}")


# ── plotting ──────────────────────────────────────────────────────────────────

def plot_top_slices(top, outpath):
    """
    2D slice plots for the top candidate:
    xy-plane at z=0 and xz-plane at y=0.
    Columns: phi_0, phi_1, phi_2, P, V(P)
    """
    phi   = top["_phi"]
    P     = top["_P"]
    V_den = top["_V_den"]
    N     = top["N"]
    L     = top["L"]

    mid = N // 2
    extent = [-L, L, -L, L]

    slices_xy = [phi[0][:, :, mid], phi[1][:, :, mid], phi[2][:, :, mid],
                 P[:, :, mid],      V_den[:, :, mid]]
    slices_xz = [phi[0][:, mid, :], phi[1][:, mid, :], phi[2][:, mid, :],
                 P[:, mid, :],      V_den[:, mid, :]]
    titles    = [r"$\phi_0$", r"$\phi_1$", r"$\phi_2$", r"$P$", r"$V(P)$"]

    fig, axes = plt.subplots(2, 5, figsize=(18, 7))
    fig.suptitle(f"Top candidate: {top['fname']}\n"
                 f"E_pot={top['E_pot']:.3f}  R_rms={top['R_rms']:.3f}  "
                 f"N_clusters={top['n_clusters']}", fontsize=11)

    for col, (s_xy, s_xz, title) in enumerate(zip(slices_xy, slices_xz, titles)):
        vmax_xy = max(abs(s_xy.max()), abs(s_xy.min()), 1e-9)
        vmax_xz = max(abs(s_xz.max()), abs(s_xz.min()), 1e-9)

        cmap = "RdBu_r" if col < 4 else "viridis"
        cbar_label = title

        im0 = axes[0, col].imshow(
            s_xy.T, origin="lower", extent=extent,
            cmap=cmap, vmin=-vmax_xy, vmax=vmax_xy, aspect="equal"
        )
        axes[0, col].set_title(f"{title}  xy (z=0)")
        axes[0, col].set_xlabel("x"); axes[0, col].set_ylabel("y")
        fig.colorbar(im0, ax=axes[0, col], fraction=0.046)

        im1 = axes[1, col].imshow(
            s_xz.T, origin="lower", extent=extent,
            cmap=cmap, vmin=-vmax_xz, vmax=vmax_xz, aspect="equal"
        )
        axes[1, col].set_title(f"{title}  xz (y=0)")
        axes[1, col].set_xlabel("x"); axes[1, col].set_ylabel("z")
        fig.colorbar(im1, ax=axes[1, col], fraction=0.046)

    plt.tight_layout()
    plt.savefig(outpath, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {outpath}")


def plot_radial_profiles(results, outpath, n_top=3):
    """Radial profiles of |phi|^2, |P|, V(P) for top-n candidates."""
    results_sorted = sorted(results, key=lambda r: r["E_pot"])[:n_top]

    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    fig.suptitle(f"Radial profiles — top {n_top} candidates (ranked by E_pot)", fontsize=11)

    colors = plt.cm.tab10(np.linspace(0, 0.6, n_top))

    for ax, (ykey, ylabel) in zip(
        axes,
        [("radial_phi2", r"$\langle|\phi|^2\rangle$"),
         ("radial_absP", r"$\langle|P|\rangle$"),
         ("radial_V",    r"$\langle V(P)\rangle$")]
    ):
        for r, col in zip(results_sorted, colors):
            label = f"{r['fname']} (gen={r['gen']}, E_pot={r['E_pot']:.2f})"
            ax.plot(r["r_centers"], r[ykey], color=col, label=label, lw=1.5)
        ax.set_xlabel("r (code units)")
        ax.set_ylabel(ylabel)
        ax.axhline(0, color="k", lw=0.5, ls="--")
        ax.legend(fontsize=6, loc="upper right")
        ax.set_title(ylabel)

    plt.tight_layout()
    plt.savefig(outpath, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {outpath}")


def plot_evolution_trajectory(results, outpath):
    """E_pot vs generation for the pop128 series (checkpoint files) + finals."""

    # Group by series
    series_data = {}
    for r in results:
        s = r["series"]
        series_data.setdefault(s, []).append((r["gen"], r["E_pot"], r["fname"]))

    # Sort each series by generation
    for s in series_data:
        series_data[s].sort(key=lambda x: x[0])

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.set_title("Evolution trajectory: E_pot vs generation", fontsize=11)
    ax.set_xlabel("Generation")
    ax.set_ylabel("E_pot (code units, more negative = better binding)")
    ax.axhline(0, color="k", lw=0.5, ls="--")

    colors = plt.cm.tab10(np.linspace(0, 0.9, len(series_data)))
    for (s, data), col in zip(sorted(series_data.items()), colors):
        gens  = [d[0] for d in data]
        epots = [d[1] for d in data]
        fnames = [d[2] for d in data]
        # distinguish checkpoint vs best
        is_ckpt = ["ckpt" in fn for fn in fnames]
        ckpt_gens  = [g for g, c in zip(gens, is_ckpt) if c]
        ckpt_epots = [e for e, c in zip(epots, is_ckpt) if c]
        best_gens  = [g for g, c in zip(gens, is_ckpt) if not c]
        best_epots = [e for e, c in zip(epots, is_ckpt) if not c]

        if ckpt_gens:
            ax.plot(ckpt_gens, ckpt_epots, "o-", color=col, lw=1.5,
                    markersize=5, label=f"{s} (ckpts)")
        if best_gens:
            ax.plot(best_gens, best_epots, "D", color=col, markersize=9,
                    label=f"{s} (best)", zorder=5)

    ax.legend(fontsize=8)
    ax.invert_yaxis()   # more negative at top (better)
    ax.set_ylim(ax.get_ylim())   # let matplotlib decide, then freeze
    plt.tight_layout()
    plt.savefig(outpath, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {outpath}")


def plot_energy_decomposition(results, outpath):
    """Stacked bar chart of energy components for all candidates."""
    results_sorted = sorted(results, key=lambda r: r["E_pot"])

    fnames = [r["fname"].replace(".npz", "").replace("_", " ") for r in results_sorted]
    E_pot  = [r["E_pot"]  for r in results_sorted]
    E_grad = [r["E_grad"] for r in results_sorted]
    E_mass = [r["E_mass"] for r in results_sorted]
    E_kin  = [r["E_kin"]  for r in results_sorted]

    x = np.arange(len(fnames))
    width = 0.6

    fig, ax = plt.subplots(figsize=(max(10, 2 * len(fnames)), 5))
    ax.bar(x, E_grad, width, label="E_grad", color="steelblue")
    ax.bar(x, E_mass, width, bottom=E_grad, label="E_mass", color="orange")
    bottom2 = [g + m for g, m in zip(E_grad, E_mass)]
    ax.bar(x, E_kin,  width, bottom=bottom2, label="E_kin", color="green")
    ax.bar(x, E_pot,  width, label="E_pot",  color="red", alpha=0.85)

    ax.set_xticks(x)
    ax.set_xticklabels(fnames, rotation=30, ha="right", fontsize=7)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_ylabel("Energy (code units)")
    ax.set_title("Energy decomposition per candidate (ranked by E_pot)")
    ax.legend()
    plt.tight_layout()
    plt.savefig(outpath, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {outpath}")


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dir",    default="candidates/",
                        help="Directory containing .npz candidate files")
    parser.add_argument("--outdir", default="analysis/",
                        help="Directory for output files")
    args = parser.parse_args()

    cdir   = args.dir
    outdir = args.outdir

    if not os.path.isdir(cdir):
        sys.exit(f"ERROR: candidate directory not found: {cdir}")

    os.makedirs(outdir, exist_ok=True)

    # ── discover candidates ───────────────────────────────────────────────────
    npz_files = sorted(
        os.path.join(cdir, fn)
        for fn in os.listdir(cdir)
        if fn.endswith(".npz")
    )
    if not npz_files:
        sys.exit(f"ERROR: no .npz files found in {cdir}")

    print(f"Found {len(npz_files)} candidate files in '{cdir}'")

    # ── analyze each candidate ────────────────────────────────────────────────
    results = []
    for path in npz_files:
        fname = os.path.basename(path)
        print(f"  Analyzing {fname} ...", end=" ", flush=True)
        try:
            r = analyze_candidate(path)
            results.append(r)
            print(f"E_pot={r['E_pot']:.3f}  R_rms={r['R_rms']:.3f}  "
                  f"N_clust={r['n_clusters']}")
        except Exception as exc:
            print(f"FAILED: {exc}")

    if not results:
        sys.exit("ERROR: no candidates could be analyzed.")

    # ── sort by E_pot ─────────────────────────────────────────────────────────
    results_sorted = sorted(results, key=lambda r: r["E_pot"])
    top            = results_sorted[0]

    # ── print to stdout ───────────────────────────────────────────────────────
    print()
    print_summary_table(results_sorted)

    # ── write TSV summary ─────────────────────────────────────────────────────
    write_summary_tsv(results_sorted, os.path.join(outdir, "summary.tsv"))

    # ── plots ─────────────────────────────────────────────────────────────────
    print("Generating plots ...")

    plot_top_slices(
        top,
        os.path.join(outdir, "top_candidate_slices.png")
    )

    plot_radial_profiles(
        results_sorted,
        os.path.join(outdir, "radial_profiles.png"),
        n_top=min(3, len(results_sorted))
    )

    plot_evolution_trajectory(
        results_sorted,
        os.path.join(outdir, "evolution_trajectory.png")
    )

    plot_energy_decomposition(
        results_sorted,
        os.path.join(outdir, "energy_decomposition.png")
    )

    print()
    print(f"All output written to '{outdir}'")
    print(f"Top candidate by E_pot: {top['fname']}  (E_pot={top['E_pot']:.4f})")


if __name__ == "__main__":
    main()
