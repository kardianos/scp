"""
2D Token / Update-Budget Cellular Substrate (SCP v81 P3)

Ontology (shape 6 + 4 pure):
  - Tokens on directed lattice links are the sole currency (update budget).
  - c = hard per-tick hop cap on every directed bond (engine law, not a CFL note).
  - Charge = topological circulation / vortex handedness of the token flow.
  - Mass ≈ sequestered tokens cycling inside a coherent pattern.
  - No Maxwell field, no multiplet φ, no pairwise Coulomb, no scp_sim.

State:
  f[y, x, d]  non-negative int  — tokens at site ready to hop in direction d
              d ∈ {E=0, N=1, W=2, S=3}
  rest[y, x]  non-negative int  — optionally parked tokens (sequestered mass bucket)

Step:
  1. Collide (local, integer, mass + momentum conserving).
  2. Stream with hop cap c (excess stays; transfer ≤ c exact).
  3. Diagnostics: total tokens, max bond transfer, vorticity, vortex peaks.

All transfers are conservative: sum(f)+sum(rest) is bit-exact invariant.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Tuple, List, Dict, Any

import numpy as np

# Directions: E, N, W, S
E, N, W, S = 0, 1, 2, 3
DIRS = (E, N, W, S)
# unit vectors (dy, dx) in array index space (row=y, col=x)
EY = np.array([0, -1, 0, 1], dtype=np.int32)  # N decreases row
EX = np.array([1, 0, -1, 0], dtype=np.int32)


@dataclass
class TokenCA:
    """Integer directed lattice-gas with hard hop cap c."""

    ny: int
    nx: int
    c: int = 4  # hop cap per directed bond per tick
    omega: float = 1.0  # BGK relaxation toward integer equilibrium (0..1]
    seed: int = 0

    f: np.ndarray = field(init=False, repr=False)
    rest: np.ndarray = field(init=False, repr=False)
    rng: np.random.Generator = field(init=False, repr=False)

    # bookkeeping for hop-cap binding
    last_max_transfer: int = 0
    last_total_transfer: int = 0
    tick: int = 0

    def __post_init__(self) -> None:
        if self.ny < 8 or self.nx < 8:
            raise ValueError("grid too small")
        if self.c < 1:
            raise ValueError("hop cap c must be >= 1")
        self.f = np.zeros((self.ny, self.nx, 4), dtype=np.int64)
        self.rest = np.zeros((self.ny, self.nx), dtype=np.int64)
        self.rng = np.random.default_rng(self.seed)

    # ------------------------------------------------------------------ state
    def total_tokens(self) -> int:
        return int(self.f.sum() + self.rest.sum())

    def density(self) -> np.ndarray:
        return self.f.sum(axis=2) + self.rest

    def momentum(self) -> Tuple[np.ndarray, np.ndarray]:
        """Site momentum (mx, my) from directed occupations."""
        mx = self.f[:, :, E] - self.f[:, :, W]
        my = self.f[:, :, N] - self.f[:, :, S]  # N is -row in index, but +my chart
        # Chart: +mx = +x (east), +my = +y (north = decreasing row index in arrays)
        return mx.astype(np.float64), my.astype(np.float64)

    def velocity(self) -> Tuple[np.ndarray, np.ndarray]:
        rho = self.density().astype(np.float64)
        mx, my = self.momentum()
        inv = np.where(rho > 0, 1.0 / rho, 0.0)
        return mx * inv, my * inv

    def clear(self) -> None:
        self.f[:] = 0
        self.rest[:] = 0
        self.tick = 0
        self.last_max_transfer = 0
        self.last_total_transfer = 0

    def fill_uniform(self, rho: int) -> None:
        """Isotropic background density (rest + equal directed)."""
        if rho < 0:
            raise ValueError("rho >= 0")
        base, rem = divmod(int(rho), 5)
        self.f[:] = base
        self.rest[:] = base
        # dump remainder into rest
        self.rest += rem

    # --------------------------------------------------------------- collide
    def collide(self) -> None:
        """
        Integer BGK toward D2Q4 equilibrium; mass exact, momentum approx.

        feq_E = (ρ_dir + 2 mx)/4 etc. Supersonic sites (|m| large) are
        projected onto the non-negative orthant and renormalized so Σ feq = ρ.
        Rest tokens are sequestered (do not participate in collide/stream).
        """
        f = self.f
        rho = f.sum(axis=2)  # int64
        mx = f[:, :, E] - f[:, :, W]
        my = f[:, :, N] - f[:, :, S]

        r = rho.astype(np.float64)
        mxf = mx.astype(np.float64)
        myf = my.astype(np.float64)
        feq = np.empty_like(f, dtype=np.float64)
        feq[:, :, E] = 0.25 * r + 0.5 * mxf
        feq[:, :, W] = 0.25 * r - 0.5 * mxf
        feq[:, :, N] = 0.25 * r + 0.5 * myf
        feq[:, :, S] = 0.25 * r - 0.5 * myf

        # Non-negative mass-preserving projection (kills excess momentum if needed)
        np.maximum(feq, 0.0, out=feq)
        s = feq.sum(axis=2, keepdims=True)
        # where s==0 (rho==0): leave zeros; else renormalize to rho
        scale = np.divide(
            r[:, :, None],
            s,
            out=np.zeros_like(s),
            where=s > 0,
        )
        feq *= scale

        w = float(self.omega)
        target = (1.0 - w) * f.astype(np.float64) + w * feq
        np.maximum(target, 0.0, out=target)
        # re-normalize target to exact rho (float drift / clip)
        ts = target.sum(axis=2, keepdims=True)
        tscale = np.divide(
            r[:, :, None],
            ts,
            out=np.zeros_like(ts),
            where=ts > 0,
        )
        target *= tscale

        # Integerize: largest-remainder method (Hamilton) — exact mass
        floor = np.floor(target + 1e-9).astype(np.int64)
        residual = rho - floor.sum(axis=2)
        # Strip over-assigned (rare float edge)
        if np.any(residual < 0):
            for d in DIRS:
                take = np.minimum(floor[:, :, d], np.maximum(-residual, 0))
                floor[:, :, d] -= take
                residual += take
                if not np.any(residual < 0):
                    break
        frac = target - np.floor(target + 1e-9)
        # Award +1 to top-k fractional directions (k = residual ∈ {0,1,2,3,4})
        order = np.argsort(-frac, axis=2)
        for k in range(4):
            d_k = order[:, :, k]
            give = (residual > k).astype(np.int64)
            # scatter give into floor along direction d_k
            for d in DIRS:
                floor[:, :, d] += give * (d_k == d)
        # residual > 4 dump into E
        extra = residual - 4
        floor[:, :, E] += np.maximum(extra, 0)

        # HPP-style head-on scatter (mass + momentum conserving)
        fe = floor[:, :, E]
        fn = floor[:, :, N]
        fw = floor[:, :, W]
        fs = floor[:, :, S]
        pairs_ew = np.minimum(fe, fw)
        pairs_ns = np.minimum(fn, fs)
        scatter_ew = pairs_ew // 2
        scatter_ns = pairs_ns // 2
        do_ew = scatter_ew > scatter_ns
        do_ns = scatter_ns > scatter_ew
        se = np.where(do_ew, scatter_ew - scatter_ns, 0)
        sn = np.where(do_ns, scatter_ns - scatter_ew, 0)
        floor[:, :, E] = fe - se + sn
        floor[:, :, W] = fw - se + sn
        floor[:, :, N] = fn - sn + se
        floor[:, :, S] = fs - sn + se

        # Final hard assert path: repair if anything still drifted (paranoia)
        diff = int(floor.sum() - rho.sum())
        if diff != 0:
            # dump/remove at center cell
            cy, cx = self.ny // 2, self.nx // 2
            if diff > 0:
                # remove tokens
                for d in DIRS:
                    take = min(diff, int(floor[cy, cx, d]))
                    floor[cy, cx, d] -= take
                    diff -= take
                    if diff == 0:
                        break
            else:
                floor[cy, cx, E] += -diff

        self.f = floor

    # ---------------------------------------------------------------- stream
    def stream(self) -> None:
        """
        Stream with hard hop cap c.

        transfer = min(f[d], c); excess stays at origin in same direction.
        Periodic BC. Records max and total bond transfers this tick.
        """
        f = self.f
        c = int(self.c)
        ny, nx, _ = f.shape
        new_f = np.zeros_like(f)
        max_tr = 0
        tot_tr = 0

        # E: (y,x) -> (y, x+1)
        tr = np.minimum(f[:, :, E], c)
        stay = f[:, :, E] - tr
        new_f[:, :, E] += stay
        new_f[:, 1:, E] += tr[:, :-1]
        new_f[:, 0, E] += tr[:, -1]
        max_tr = max(max_tr, int(tr.max()) if tr.size else 0)
        tot_tr += int(tr.sum())

        # W: (y,x) -> (y, x-1)
        tr = np.minimum(f[:, :, W], c)
        stay = f[:, :, W] - tr
        new_f[:, :, W] += stay
        new_f[:, :-1, W] += tr[:, 1:]
        new_f[:, -1, W] += tr[:, 0]
        max_tr = max(max_tr, int(tr.max()) if tr.size else 0)
        tot_tr += int(tr.sum())

        # N: (y,x) -> (y-1, x)  (north = smaller row)
        tr = np.minimum(f[:, :, N], c)
        stay = f[:, :, N] - tr
        new_f[:, :, N] += stay
        new_f[:-1, :, N] += tr[1:, :]
        new_f[-1, :, N] += tr[0, :]
        max_tr = max(max_tr, int(tr.max()) if tr.size else 0)
        tot_tr += int(tr.sum())

        # S: (y,x) -> (y+1, x)
        tr = np.minimum(f[:, :, S], c)
        stay = f[:, :, S] - tr
        new_f[:, :, S] += stay
        new_f[1:, :, S] += tr[:-1, :]
        new_f[0, :, S] += tr[-1, :]
        max_tr = max(max_tr, int(tr.max()) if tr.size else 0)
        tot_tr += int(tr.sum())

        self.f = new_f
        self.last_max_transfer = max_tr
        self.last_total_transfer = tot_tr

    def step(self, n: int = 1) -> None:
        for _ in range(n):
            self.collide()
            self.stream()
            self.tick += 1

    # ----------------------------------------------------------- vorticity
    def vorticity(self) -> np.ndarray:
        """
        Discrete curl of velocity on cell-centered sites (periodic).

        ω ≈ ∂vx/∂y_chart - ∂vy/∂x but with array indices:
        chart y increases north (decreasing row), so
        ∂/∂y_chart = -∂/∂row.
        We report ω = d(vx)/d(row_neighbor_diff) style circulation density:
        ω[i,j] = 0.5*((vx[i,j+1]-vx[i,j-1]) - (vy[i+1,j]-vy[i-1,j])) with
        chart fix: vy_chart = +my component already.

        Simpler plaquette circulation (lower-left index i,j):
          C = vx[i,j] + vy[i, j+1] - vx[i+1, j] - vy[i, j]
        assigned to site (i,j). Sign: positive = CCW in chart (x east, y south-in-array...):
        with vy = N-S (positive when more N-bound tokens), positive C = CCW in
        (x right, y down-on-screen) which is CW in mathematical (x right, y up).
        We keep this sign convention consistently; charge = sign(C).
        """
        vx, vy = self.velocity()
        # plaquette circulation, periodic
        C = (
            vx
            + np.roll(vy, -1, axis=1)
            - np.roll(vx, -1, axis=0)
            - vy
        )
        return C

    def link_flux_curl(self) -> np.ndarray:
        """
        Integer-friendly circulation from directed occupations (not velocity).
        Uses mean east-bound and north-bound populations as edge proxies.
        """
        fe = self.f[:, :, E].astype(np.float64)
        fn = self.f[:, :, N].astype(np.float64)
        fw = self.f[:, :, W].astype(np.float64)
        fs = self.f[:, :, S].astype(np.float64)
        # east-edge flux proxy at site: fe - roll(fw from east neighbor)
        je = fe - np.roll(fw, -1, axis=1)
        jn = fn - np.roll(fs, 1, axis=0)  # N minus S from north neighbor
        C = je + np.roll(jn, -1, axis=1) - np.roll(je, -1, axis=0) - jn
        return C

    # -------------------------------------------------------------- seeding
    def _set_velocity_field(self, vx: np.ndarray, vy: np.ndarray) -> None:
        """
        Reassign directed occupations to match (rho_dir, vx, vy) without
        changing site directed-mass. Rest untouched. |v| clipped to < 1.
        """
        rho = self.f.sum(axis=2).astype(np.float64)
        # clip speed so D2Q4 equilibrium stays non-negative: |v| < 0.5 safe
        speed = np.sqrt(vx * vx + vy * vy)
        clip = np.maximum(speed / 0.45, 1.0)
        vx = vx / clip
        vy = vy / clip
        # feq with momentum m = rho * v
        mxf = rho * vx
        myf = rho * vy
        feq = np.empty_like(self.f, dtype=np.float64)
        feq[:, :, E] = 0.25 * rho + 0.5 * mxf
        feq[:, :, W] = 0.25 * rho - 0.5 * mxf
        feq[:, :, N] = 0.25 * rho + 0.5 * myf
        feq[:, :, S] = 0.25 * rho - 0.5 * myf
        np.maximum(feq, 0.0, out=feq)
        s = feq.sum(axis=2, keepdims=True)
        scale = np.divide(rho[:, :, None], s, out=np.zeros_like(s), where=s > 0)
        feq *= scale
        # integerize exact mass per site (vectorized largest-remainder)
        rho_i = rho.astype(np.int64)
        floor = np.floor(feq + 1e-9).astype(np.int64)
        residual = rho_i - floor.sum(axis=2)
        frac = feq - np.floor(feq + 1e-9)
        order = np.argsort(-frac, axis=2)
        for k in range(4):
            d_k = order[:, :, k]
            give = (residual > k).astype(np.int64)
            for d in DIRS:
                floor[:, :, d] += give * (d_k == d)
        floor[:, :, E] += np.maximum(residual - 4, 0)
        self.f = floor

    def add_vortex(
        self,
        cy: float,
        cx: float,
        gamma: float,
        core: float = 3.0,
        amp: float = 1.0,
        accumulate_v: Optional[Tuple[np.ndarray, np.ndarray]] = None,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Soft-core vortex velocity field (Lamb–Oseen / regularized).

        gamma > 0 → positive C in our plaquette convention after seeding.
        If accumulate_v is (vx,vy) arrays, add into them and do NOT assign f yet.
        Otherwise assign immediately onto current density.
        amp rescales gamma (kept for call-site compatibility).
        """
        yy, xx = np.mgrid[0 : self.ny, 0 : self.nx]
        dy = yy.astype(np.float64) - cy  # row offset (south positive in array)
        dx = xx.astype(np.float64) - cx
        dy -= self.ny * np.round(dy / self.ny)
        dx -= self.nx * np.round(dx / self.nx)
        r2 = dx * dx + dy * dy + core * core
        # Math frame: x east, y north (= -row). CCW Γ>0:
        #   v_x = -Γ (y_n - y0) / r2 = -Γ (-dy_row) / r2 = Γ dy_row / r2
        #   v_n =  Γ (x - x0) / r2 = Γ dx / r2
        g = float(gamma) * float(amp)
        vx = g * dy / r2
        vy = g * dx / r2
        if accumulate_v is not None:
            # in-place slice add (augmented assign on tuple subscript rebinds)
            accumulate_v[0][:] += vx
            accumulate_v[1][:] += vy
            return accumulate_v
        self._set_velocity_field(vx, vy)
        return vx, vy

    def add_vortex_pair(
        self,
        separation: float,
        gamma: float = 40.0,
        core: float = 3.0,
        amp: float = 1.0,
        opposite: bool = True,
        cy: Optional[float] = None,
        cx: Optional[float] = None,
    ) -> Tuple[Tuple[float, float], Tuple[float, float]]:
        """Place ± or ++ vortex pair along x; reassign f from superposed v."""
        if cy is None:
            cy = 0.5 * (self.ny - 1)
        if cx is None:
            cx = 0.5 * (self.nx - 1)
        half = 0.5 * separation
        c1 = (cy, cx - half)
        c2 = (cy, cx + half)
        vx = np.zeros((self.ny, self.nx), dtype=np.float64)
        vy = np.zeros((self.ny, self.nx), dtype=np.float64)
        self.add_vortex(c1[0], c1[1], +gamma, core=core, amp=amp, accumulate_v=(vx, vy))
        g2 = -gamma if opposite else +gamma
        self.add_vortex(c2[0], c2[1], g2, core=core, amp=amp, accumulate_v=(vx, vy))
        self._set_velocity_field(vx, vy)
        return c1, c2

    def centroid_vorticity(
        self,
        sign: int,
        guess_y: float,
        guess_x: float,
        window: int = 12,
    ) -> Optional[Tuple[float, float, float]]:
        """
        Signed-vorticity centroid in a window around a guess (minimal image).
        Returns (y, x, strength) or None if strength too small.
        """
        C = self.vorticity()
        yy, xx = np.mgrid[0 : self.ny, 0 : self.nx]
        dy = yy.astype(np.float64) - guess_y
        dx = xx.astype(np.float64) - guess_x
        dy -= self.ny * np.round(dy / self.ny)
        dx -= self.nx * np.round(dx / self.nx)
        mask = (np.abs(dy) <= window) & (np.abs(dx) <= window)
        field = sign * C
        w = np.where(mask & (field > 0), field, 0.0)
        mass = float(w.sum())
        if mass < 1e-8:
            return None
        # centroid in offset coords then map back
        cy_off = float((w * dy).sum() / mass)
        cx_off = float((w * dx).sum() / mass)
        y = (guess_y + cy_off) % self.ny
        x = (guess_x + cx_off) % self.nx
        return (y, x, sign * mass)

    def sequester_disk(self, cy: float, cx: float, radius: float, amount: int) -> None:
        """Park `amount` tokens into rest inside a disk (mass sequestration demo)."""
        yy, xx = np.mgrid[0 : self.ny, 0 : self.nx]
        dy = yy.astype(np.float64) - cy
        dx = xx.astype(np.float64) - cx
        dy -= self.ny * np.round(dy / self.ny)
        dx -= self.nx * np.round(dx / self.nx)
        mask = (dx * dx + dy * dy) <= radius * radius
        sites = np.argwhere(mask)
        if len(sites) == 0 or amount <= 0:
            return
        per = amount // len(sites)
        rem = amount - per * len(sites)
        self.rest[mask] += per
        for k in range(rem):
            y, x = sites[k]
            self.rest[y, x] += 1

    # ----------------------------------------------------------- diagnostics
    def find_vortices(
        self,
        n_peaks: int = 2,
        smooth: float = 1.5,
        min_abs: float = 0.05,
    ) -> Dict[str, Any]:
        """
        Find positive and negative vorticity peaks after Gaussian smooth.

        Returns dict with plus/minus lists of (y, x, strength) and separation
        between strongest + and strongest - (minimal image).
        """
        from math import exp

        C = self.vorticity()
        # separable Gaussian smooth
        sig = float(smooth)
        if sig > 0:
            rad = max(1, int(3 * sig))
            ax = np.arange(-rad, rad + 1, dtype=np.float64)
            ker = np.exp(-0.5 * (ax / sig) ** 2)
            ker /= ker.sum()
            # convolve periodic
            Cs = C
            for axis in (0, 1):
                pad = rad
                if axis == 0:
                    ext = np.pad(Cs, ((pad, pad), (0, 0)), mode="wrap")
                    out = np.zeros_like(Cs)
                    for i, k in enumerate(ker):
                        out += k * ext[i : i + self.ny, :]
                    Cs = out
                else:
                    ext = np.pad(Cs, ((0, 0), (pad, pad)), mode="wrap")
                    out = np.zeros_like(Cs)
                    for i, k in enumerate(ker):
                        out += k * ext[:, i : i + self.nx]
                    Cs = out
        else:
            Cs = C

        def peaks(field: np.ndarray, sign: int) -> List[Tuple[float, float, float]]:
            F = sign * field
            # non-maximum suppression in 3x3
            mx = F.copy()
            for dy in (-1, 0, 1):
                for dx in (-1, 0, 1):
                    if dy == 0 and dx == 0:
                        continue
                    mx = np.maximum(mx, np.roll(np.roll(F, dy, 0), dx, 1))
            is_peak = (F == mx) & (F >= min_abs)
            ys, xs = np.where(is_peak)
            vals = F[ys, xs]
            order = np.argsort(-vals)
            out: List[Tuple[float, float, float]] = []
            for idx in order[:n_peaks]:
                y, x = float(ys[idx]), float(xs[idx])
                out.append((y, x, float(sign * vals[idx])))
            return out

        plus = peaks(Cs, +1)
        minus = peaks(Cs, -1)

        sep = None
        mid = None
        if plus and minus:
            y1, x1, _ = plus[0]
            y2, x2, _ = minus[0]
            dy = y2 - y1
            dx = x2 - x1
            dy -= self.ny * round(dy / self.ny)
            dx -= self.nx * round(dx / self.nx)
            sep = float(np.hypot(dy, dx))
            mid = (0.5 * (y1 + y2), 0.5 * (x1 + x2))

        return {
            "plus": plus,
            "minus": minus,
            "separation": sep,
            "midpoint": mid,
            "C_max": float(Cs.max()),
            "C_min": float(Cs.min()),
            "C_abs_mean": float(np.mean(np.abs(Cs))),
            "C_smooth": Cs,
        }

    def kinetic_energy(self) -> float:
        mx, my = self.momentum()
        rho = self.density().astype(np.float64)
        # KE ~ |m|^2 / (2 rho)
        safe = np.maximum(rho, 1.0)
        return float(0.5 * np.sum((mx * mx + my * my) / safe))

    def snapshot(self) -> Dict[str, Any]:
        vinfo = self.find_vortices()
        return {
            "tick": self.tick,
            "total_tokens": self.total_tokens(),
            "max_transfer": self.last_max_transfer,
            "total_transfer": self.last_total_transfer,
            "KE": self.kinetic_energy(),
            "C_max": vinfo["C_max"],
            "C_min": vinfo["C_min"],
            "C_abs_mean": vinfo["C_abs_mean"],
            "separation": vinfo["separation"],
            "plus": vinfo["plus"],
            "minus": vinfo["minus"],
        }


def minimal_image_delta(
    y1: float, x1: float, y2: float, x2: float, ny: int, nx: int
) -> Tuple[float, float]:
    dy = y2 - y1
    dx = x2 - x1
    dy -= ny * round(dy / ny)
    dx -= nx * round(dx / nx)
    return dy, dx
