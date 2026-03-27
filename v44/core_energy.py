#!/usr/bin/env python3
"""Compute core-integrated energy from diag.tsv for mass defect analysis.

Instead of SFA (slow, 6GB per frame), use the diag.tsv time series with
smarter analysis:
1. Separate the breathing oscillation from the DC component
2. Use running averages over the breathing period
3. Compare matched time windows
"""
import csv, sys, numpy as np

def load_diag(path):
    with open(path) as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    data = {}
    for key in rows[0]:
        try:
            data[key] = np.array([float(r[key]) for r in rows])
        except:
            pass
    return data

def analyze(data, label, t_start=200, t_end=None):
    t = data['t']
    if t_end is None:
        t_end = t[-1]
    mask = (t >= t_start) & (t <= t_end)
    
    epot = data['E_pot'][mask]
    etot = data['E_total'][mask]
    pint = data['P_int'][mask]
    ekin = data['E_phi_kin'][mask]
    egrad = data['E_grad'][mask]
    emass = data['E_mass'][mask]
    etkin = data['E_theta_kin'][mask]
    etgrad = data['E_tgrad'][mask]
    ecoup = data['E_coupling'][mask]
    
    # The breathing period for a proton is ~150t (from CONCEPT.md)
    # Use a rolling average over 150t to smooth the oscillation
    window = 150  # samples (diag_dt=1)
    if len(epot) < window:
        window = len(epot) // 3
    
    def rolling_mean(x, w):
        if w <= 1: return x
        kernel = np.ones(w) / w
        return np.convolve(x, kernel, mode='valid')
    
    epot_smooth = rolling_mean(epot, window)
    etot_smooth = rolling_mean(etot, window)
    pint_smooth = rolling_mean(pint, window)
    
    print(f"\n{label} (t={t_start}–{t_end}, {len(epot)} samples):")
    print(f"  Raw:     <E_pot> = {epot.mean():.1f}  std = {epot.std():.1f}  (oscillation amplitude)")
    print(f"  Smooth:  <E_pot> = {epot_smooth.mean():.1f}  std = {epot_smooth.std():.1f}  (after {window}t rolling avg)")
    print(f"  Raw:     <E_total> = {etot.mean():.1f}  std = {etot.std():.1f}")
    print(f"  Smooth:  <E_total> = {etot_smooth.mean():.1f}  std = {etot_smooth.std():.1f}")
    print(f"  Raw:     <P_int> = {pint.mean():.1f}  std = {pint.std():.1f}")
    print(f"  Smooth:  <P_int> = {pint_smooth.mean():.1f}  std = {pint_smooth.std():.1f}")
    print(f"  Energy breakdown (means):")
    print(f"    E_phi_kin:  {ekin.mean():.1f}")
    print(f"    E_grad:     {egrad.mean():.1f}")
    print(f"    E_mass:     {emass.mean():.1f}")
    print(f"    E_pot:      {epot.mean():.1f}")
    print(f"    E_theta_kin:{etkin.mean():.1f}")
    print(f"    E_tgrad:    {etgrad.mean():.1f}")
    print(f"    E_coupling: {ecoup.mean():.1f}")
    print(f"    E_total:    {etot.mean():.1f}")
    
    # E_total drift rate (energy loss through absorbing BC)
    dt_range = t[mask][-1] - t[mask][0]
    etot_drift = (etot_smooth[-1] - etot_smooth[0]) / dt_range if dt_range > 0 else 0
    print(f"  E_total drift rate: {etot_drift:.1f} per time unit")
    
    return {
        'epot_mean': epot.mean(), 'epot_smooth': epot_smooth.mean(), 'epot_std_smooth': epot_smooth.std(),
        'etot_mean': etot.mean(), 'etot_smooth': etot_smooth.mean(), 'etot_std_smooth': etot_smooth.std(),
        'pint_mean': pint.mean(), 'pint_smooth': pint_smooth.mean(),
        'etot_drift': etot_drift,
    }

# Load data
uud = load_diag('/home/d/code/scp/v44/uud_isolated_diag.tsv')
udd = load_diag('/home/d/code/scp/v44/udd_isolated_diag.tsv')

# Analyze in matched windows
uud_r = analyze(uud, "UUD Proton", 200, 400)  # use 200-400 (both have SFA data here)
udd_r = analyze(udd, "UDD Neutron", 200, 400)

# Also analyze UDD in later window to check drift
udd_late = analyze(udd, "UDD Neutron (late)", 500, 900)

print("\n" + "="*60)
print("F24 MASS DEFECT COMPARISON (smoothed, t=200-400)")
print("="*60)
print(f"  UUD <E_pot> (smooth): {uud_r['epot_smooth']:.1f} ± {uud_r['epot_std_smooth']:.1f}")
print(f"  UDD <E_pot> (smooth): {udd_r['epot_smooth']:.1f} ± {udd_r['epot_std_smooth']:.1f}")
print(f"  Sum:                  {uud_r['epot_smooth'] + udd_r['epot_smooth']:.1f}")
print()
print(f"  UUD <E_total> (smooth): {uud_r['etot_smooth']:.1f} ± {uud_r['etot_std_smooth']:.1f}")
print(f"  UDD <E_total> (smooth): {udd_r['etot_smooth']:.1f} ± {udd_r['etot_std_smooth']:.1f}")
print(f"  Sum:                    {uud_r['etot_smooth'] + udd_r['etot_smooth']:.1f}")
print()
print(f"  UUD E_total drift: {uud_r['etot_drift']:.1f}/t")
print(f"  UDD E_total drift: {udd_r['etot_drift']:.1f}/t")
print()
print(f"  For V42 deuterium comparison, need <E_pot> from V42 at matched")
print(f"  time window (t=200-400 or later). V42 <E_pot>=-94.7 was likely")
print(f"  from a different time window.")

# Check if V42 diag exists
import os
v42_diag = "/home/d/code/scp/v42/results/deuterium_diag.tsv"
if os.path.exists(v42_diag):
    print(f"\n  Found V42 diag at {v42_diag}")
    v42 = load_diag(v42_diag)
    v42_r = analyze(v42, "V42 Deuterium", 200, 400)
    print(f"\n  V42 <E_pot> (smooth, t=200-400): {v42_r['epot_smooth']:.1f} ± {v42_r['epot_std_smooth']:.1f}")
    print(f"  V42 <E_total> (smooth): {v42_r['etot_smooth']:.1f}")
    print(f"\n  MASS DEFECT (E_pot):")
    print(f"    E_deut - (E_UUD + E_UDD) = {v42_r['epot_smooth']:.1f} - ({uud_r['epot_smooth']:.1f} + {udd_r['epot_smooth']:.1f})")
    print(f"                             = {v42_r['epot_smooth'] - uud_r['epot_smooth'] - udd_r['epot_smooth']:.1f}")
else:
    print(f"\n  V42 diag not found at {v42_diag}")
    # Check alternative locations
    for p in ["/home/d/code/scp/v42/deuterium_diag.tsv", 
              "/home/d/code/scp/v42/results/deuterium_diag.tsv"]:
        if os.path.exists(p):
            print(f"  Found at: {p}")
