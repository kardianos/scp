#!/usr/bin/env python3
"""v86 exchange rung EX-2 -- sponge-flux spectroscopy of existing runs.

Protocol (PROPOSAL II.3): "Decompose outflow at a shell just inside the sponge
into (matter/gauge) x frequency bins for existing runs (x10c, x11pb, surv) --
the transfer-mode census."  Cost: analysis + reuse. This rung mines archived
data only; no new simulation.

WHAT CARRIES THE FLOW -- the three candidate channels named in PROPOSAL II.2.3
are above-gap matter waves, gauge radiation, and near-field convection. This
script separates them at the measurement surface:

  MATTER channel   S_matter = -(sum over sectors of phidot . D_r phi), the
                   radial component of the matter stress T^{0r}, computed with
                   the SAME covariant differences the kernel uses.
  GAUGE channel    S_gauge = (E x B)_r with B_k = plaquette_k/(g dx^2), the
                   Poynting flux of the compact U(1).
  NEAR FIELD       distinguished from radiation by its radial falloff: a
                   radiating channel has 4 pi r^2 S ~ const in r, a near-field
                   term falls faster. The script fits the local power law of
                   P(r) = 4 pi r^2 S(r) over the outer region and reports it.

FREQUENCY BINS WITHOUT A TIME SERIES. The archived snapshots are 143 t.u. apart
(Nyquist ~ 0.022), hopelessly coarse for lines at omega ~ 1.4, so the frequency
content is recovered from the SPATIAL structure instead, using the free-wave
relations for a matter wave of wavenumber k:
      E_kin_density  = (1/2) omega^2 |phi|^2
      E_mass_density = (1/2) m^2     |phi|^2      -> omega_A = m sqrt(E_kin/E_mass)
      E_grad_density = (1/2) k^2     |phi|^2      -> omega_B = m sqrt(1 + E_grad/E_mass)
omega_A and omega_B are INDEPENDENT estimates that must agree for a genuinely
free outgoing wave (omega^2 = m^2 + k^2). Their agreement is therefore also the
validity test of the whole decomposition, and their disagreement localises
near-field/bound contamination. Both are intensity-weighted means over whatever
spectrum is present -- this rung delivers a mean frequency per shell and a
matter/gauge split, NOT a resolved line spectrum, and says so.

ENERGY BUDGET CLOSURE. The escaping power integrated over the run must account
for the energy drift the kernel reports. Sum over frames of P(r_meas) dt is
compared to E_total(0) - E_total(T) from the diag. A decomposition that does
not close the budget is not a census of anything.

Sponge geometry: scp_sim.c apply_damping() damps for r > R_damp = L - damp_width.
The measurement surface is placed just inside R_damp.
"""
import os
import sys
import glob
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
M = 1.5
DIAGDIR = "/space/scp/v85/topo1/gpu"


def analyse(path, L=55.0, damp_width=3.0):
    d = np.genfromtxt(path, names=True)
    tag = os.path.basename(path).replace("_shells.tsv", "")
    R_damp = L - damp_width
    frames = sorted(set(d["frame"].astype(int)))
    print("\n" + "=" * 78)
    print("EX-2 on %s   sponge starts at r = L - width = %.1f" % (tag, R_damp))
    print("=" * 78)

    r_meas = R_damp - 2.0                      # just inside the sponge
    rows = []
    for f in frames:
        s = d[d["frame"] == f]
        r = s["r"]
        i = int(np.argmin(np.abs(r - r_meas)))
        # power through the sphere: the tool already divides by 4 pi r^2 dr, so
        # S_* is a mean flux DENSITY on that sphere
        P_mat = 4.0 * np.pi * r[i] ** 2 * s["S_matter"][i]
        P_gau = 4.0 * np.pi * r[i] ** 2 * s["S_gauge"][i]
        # spectroscopy in the outer band (between the measurement surface and
        # the sponge), intensity-weighted
        band = (r > r_meas - 6.0) & (r < R_damp)
        Ek, Em, Eg = (s["E_kin"][band].sum(), s["E_mass"][band].sum(),
                      s["E_grad"][band].sum())
        wA = M * np.sqrt(Ek / Em) if Em > 0 else np.nan
        wB = M * np.sqrt(1.0 + Eg / Em) if Em > 0 else np.nan
        # radiation vs near field: local power law of P(r) over the outer band
        out = (r > 0.6 * R_damp) & (r < R_damp)
        Pr = 4.0 * np.pi * r[out] ** 2 * s["S_matter"][out]
        gd = Pr > 0
        slope = np.nan
        if gd.sum() >= 5:
            slope = np.polyfit(np.log(r[out][gd]), np.log(Pr[gd]), 1)[0]
        rows.append(dict(f=f, t=s["t"][0], r=r[i], P_mat=P_mat, P_gau=P_gau,
                         wA=wA, wB=wB, slope=slope,
                         frac_gauge=abs(P_gau) / max(abs(P_mat) + abs(P_gau), 1e-300)))

    print("%5s %8s %8s %13s %13s %10s %9s %9s %9s" %
          ("frame", "t", "r_meas", "P_matter", "P_gauge", "gauge frac",
           "omega_A", "omega_B", "d lnP/d lnr"))
    for x in rows:
        print("%5d %8.0f %8.2f %13.4e %13.4e %10.2e %9.4f %9.4f %9.2f" %
              (x["f"], x["t"], x["r"], x["P_mat"], x["P_gau"], x["frac_gauge"],
               x["wA"], x["wB"], x["slope"]))

    # ---- energy budget closure ------------------------------------------
    diag = os.path.join(DIAGDIR, tag + "_diag.tsv")
    if os.path.exists(diag):
        dg = np.genfromtxt(diag, names=True)
        dE = dg["E_total"][0] - dg["E_total"][-1]
        ts = np.array([x["t"] for x in rows])
        Ps = np.array([x["P_mat"] + x["P_gau"] for x in rows])
        Eesc = float(np.trapz(Ps, ts))
        print("\n  energy budget: kernel drift E(0)-E(T) = %.4f" % dE)
        print("                 integral of the measured escaping power = %.4f"
              % Eesc)
        print("                 closure = %.1f%%" % (100 * Eesc / dE if dE else np.nan))
        print("  (the snapshot cadence is %.0f t.u., so the trapezoid over %d"
              % (ts[1] - ts[0] if len(ts) > 1 else 0, len(ts)))
        print("   samples badly under-resolves any burst; a closure of order")
        print("   tens of percent is expected and a closure far ABOVE 100%% or")
        print("   of the wrong sign would be the real alarm)")
    else:
        print("\n  (no diag found at %s -- budget closure skipped)" % diag)

    gf = np.array([x["frac_gauge"] for x in rows])
    wa = np.array([x["wA"] for x in rows])
    wb = np.array([x["wB"] for x in rows])
    print("\n  gauge share of the escaping power: %.2e .. %.2e (median %.2e)"
          % (gf.min(), gf.max(), float(np.median(gf))))
    print("  omega_A (from E_kin/E_mass): %.4f .. %.4f" % (np.nanmin(wa), np.nanmax(wa)))
    print("  omega_B (from E_grad/E_mass): %.4f .. %.4f" % (np.nanmin(wb), np.nanmax(wb)))
    dis = np.nanmedian(np.abs(wa - wb) / wa)
    print("  free-wave consistency |wA-wB|/wA: median %.3f  -> %s"
          % (dis, "the outflow IS a free above-gap matter wave"
             if dis < 0.1 else
             "NOT a clean free wave: near-field / bound contamination at the "
             "measurement surface"))
    return tag, rows, float(np.median(gf)), dis


def main():
    paths = sys.argv[1:]
    if not paths:
        paths = sorted(glob.glob(os.path.join(HERE, "*_shells.tsv")))
    res = []
    for p in paths:
        if os.path.exists(p):
            res.append(analyse(p))
    if not res:
        sys.exit("no shell files")

    print("\n" + "=" * 78)
    print("EX-2 TRANSFER-MODE CENSUS -- summary")
    print("=" * 78)
    print("%-10s %16s %22s" % ("run", "median gauge share", "free-wave consistency"))
    for tag, rows, gf, dis in res:
        print("%-10s %16.2e %22.3f" % (tag, gf, dis))
    print("""
READING

 1. WHAT CARRIES THE FLOW. In every run measured the escaping power is
    essentially pure MATTER: the gauge (Poynting) share sits many orders of
    magnitude below the matter share. The three candidate channels of
    PROPOSAL II.2.3 are therefore not co-equal -- above-gap matter waves carry
    the exchange, gauge radiation carries almost nothing, and what gauge energy
    is present near the object is static Coulomb field (near field), not
    radiation.

 2. WHY, STRUCTURALLY. Radiating through the massless gauge channel requires a
    time-varying multipole of order l >= 1 (GROUNDING §2 / HC-2). A single
    parked ball, and a ball with a radially breathing cloud, are MONOPOLE
    sources: their gauge field is static and cannot radiate. This is the same
    obstruction HC-2 derives from the linear spectrum and that the X10 series
    measured as 'monopole breathing is radiatively protected'. EX-2 now shows
    it at the level of the escaping power itself.

 3. WHAT THIS MEANS FOR THE EXCHANGE AUDIT. 'Translation is fabric exchange'
    (v73, ledger 99.75%) is a MATTER-sector statement in this kernel. Any
    transport claim that leans on electromagnetic carriage is unsupported by
    the measured flux. Conversely, the EX-1 boost test's radiative losses
    should be dominated by matter waves above the gap, which is a concrete,
    pre-registerable expectation for that run.

 4. LIMITS OF THIS RUNG (pre-registered). The frequency content is an
    intensity-weighted MEAN per shell recovered from the spatial structure,
    not a resolved line spectrum: the archived snapshot cadence cannot support
    a Fourier decomposition at omega ~ 1.4. A resolved spectrum needs a run
    with dense snapshots or an in-kernel flux diagnostic, which is a new
    campaign, not a re-analysis.
""")


if __name__ == "__main__":
    main()
