#!/usr/bin/env python3
"""Build composite.tsv combining results, cluster analysis, and diag data
for the Lissajous chirality sweep.

Outputs one row per run, sorted by final P_int descending.

Columns:
  From results.tsv: d0, d1, d2, l0, l1, l2, E_init, E_final, phi_i, phi_f, Pint_i, Pint_f, trms_f, status
  From cluster_analysis.tsv: n_clusters, n_big_clusters, largest_mass, largest_Epot, largest_rms_r,
                              largest_P_peak, total_cluster_mass, second_mass
  From diag files: E_conservation, Pint_t0, Pint_t50, Pint_t100, trms_t0, trms_t50, trms_t100
"""

import os
import sys
import csv
from collections import defaultdict

BASE = "/space/scp/v52/lissajous"

def read_tsv(path):
    """Read TSV file, return (header_list, list_of_row_dicts)."""
    rows = []
    with open(path) as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        for row in reader:
            if not row or not row[0].strip():
                continue
            d = {}
            for i, h in enumerate(header):
                if i < len(row):
                    d[h] = row[i]
            rows.append(d)
    return header, rows


def parse_float(s, default=0.0):
    try:
        return float(s)
    except (ValueError, TypeError):
        return default


def get_diag_at_time(diag_rows, target_t):
    """Find the row closest to target_t and return it."""
    if not diag_rows:
        return None
    best = None
    best_dt = float('inf')
    for row in diag_rows:
        t = parse_float(row.get('t', ''), -1)
        dt = abs(t - target_t)
        if dt < best_dt:
            best_dt = dt
            best = row
    return best


def main():
    # 1. Read results.tsv
    _, results = read_tsv(os.path.join(BASE, "results.tsv"))
    print(f"Loaded {len(results)} runs from results.tsv", file=sys.stderr)

    # 2. Read cluster_analysis.tsv and aggregate per file
    _, clusters = read_tsv(os.path.join(BASE, "cluster_analysis.tsv"))

    # Group clusters by file name
    cluster_by_file = defaultdict(list)
    for c in clusters:
        fname = c.get('file', '')
        cluster_by_file[fname].append(c)

    print(f"Loaded {len(clusters)} cluster entries across {len(cluster_by_file)} files", file=sys.stderr)

    # 3. Process each run
    composite_rows = []

    for run in results:
        name = f"{run['l0']}_{run['l1']}_{run['l2']}"
        row = dict(run)  # copy all results fields
        row['name'] = name

        # -- Cluster analysis summary --
        file_clusters = cluster_by_file.get(name, [])
        n_total = len(file_clusters)

        # Filter to "significant" clusters (mass > 1.0, at least 100 voxels)
        big_clusters = []
        for c in file_clusters:
            mass = parse_float(c.get('mass', '0'))
            nvox = int(parse_float(c.get('n_voxels', '0')))
            if mass > 1.0 and nvox > 100:
                big_clusters.append(c)

        # Sort by mass descending
        big_clusters.sort(key=lambda c: parse_float(c.get('mass', '0')), reverse=True)

        row['n_clusters_raw'] = str(n_total)
        row['n_big_clusters'] = str(len(big_clusters))

        if big_clusters:
            largest = big_clusters[0]
            row['largest_mass'] = largest.get('mass', '0')
            row['largest_Epot'] = largest.get('E_pot', '0')
            row['largest_rms_r'] = largest.get('rms_r', '0')
            row['largest_P_peak'] = largest.get('P_peak', '0')
            row['largest_phi_max'] = largest.get('phi_max', '0')
            row['largest_nvox'] = largest.get('n_voxels', '0')
        else:
            row['largest_mass'] = '0'
            row['largest_Epot'] = '0'
            row['largest_rms_r'] = '0'
            row['largest_P_peak'] = '0'
            row['largest_phi_max'] = '0'
            row['largest_nvox'] = '0'

        total_mass = sum(parse_float(c.get('mass', '0')) for c in file_clusters)
        row['total_cluster_mass'] = f"{total_mass:.6e}"

        if len(big_clusters) >= 2:
            row['second_mass'] = big_clusters[1].get('mass', '0')
        else:
            row['second_mass'] = '0'

        # -- Diag file analysis --
        diag_path = os.path.join(BASE, f"{name}_diag.tsv")
        if os.path.exists(diag_path):
            _, diag_rows = read_tsv(diag_path)

            # Initial and final rows
            if diag_rows:
                first = diag_rows[0]
                last = diag_rows[-1]

                E_init_diag = parse_float(first.get('E_total', '0'))
                E_final_diag = parse_float(last.get('E_total', '0'))
                if abs(E_init_diag) > 1e-20:
                    row['E_conservation'] = f"{(E_final_diag - E_init_diag) / E_init_diag:.6e}"
                else:
                    row['E_conservation'] = '0'

                # P_int at t=0, t=50, t=100
                r0 = get_diag_at_time(diag_rows, 0.0)
                r50 = get_diag_at_time(diag_rows, 50.0)
                r100 = get_diag_at_time(diag_rows, 100.0)

                row['Pint_t0'] = r0.get('P_int', '0') if r0 else '0'
                row['Pint_t50'] = r50.get('P_int', '0') if r50 else '0'
                row['Pint_t100'] = r100.get('P_int', '0') if r100 else '0'

                row['trms_t0'] = r0.get('theta_rms', '0') if r0 else '0'
                row['trms_t50'] = r50.get('theta_rms', '0') if r50 else '0'
                row['trms_t100'] = r100.get('theta_rms', '0') if r100 else '0'

                # phi_max trajectory
                row['phi_max_t0'] = r0.get('phi_max', '0') if r0 else '0'
                row['phi_max_t50'] = r50.get('phi_max', '0') if r50 else '0'
                row['phi_max_t100'] = r100.get('phi_max', '0') if r100 else '0'
            else:
                for k in ['E_conservation', 'Pint_t0', 'Pint_t50', 'Pint_t100',
                           'trms_t0', 'trms_t50', 'trms_t100',
                           'phi_max_t0', 'phi_max_t50', 'phi_max_t100']:
                    row[k] = '0'
        else:
            print(f"Warning: missing diag file for {name}", file=sys.stderr)
            for k in ['E_conservation', 'Pint_t0', 'Pint_t50', 'Pint_t100',
                       'trms_t0', 'trms_t50', 'trms_t100',
                       'phi_max_t0', 'phi_max_t50', 'phi_max_t100']:
                row[k] = '0'

        composite_rows.append(row)

    # Sort by Pint_f descending
    composite_rows.sort(key=lambda r: parse_float(r.get('Pint_f', '0')), reverse=True)

    # Write composite TSV
    out_cols = [
        'name', 'd0', 'd1', 'd2', 'l0', 'l1', 'l2',
        'E_init', 'E_final', 'E_conservation',
        'phi_i', 'phi_f',
        'Pint_i', 'Pint_f',
        'Pint_t0', 'Pint_t50', 'Pint_t100',
        'trms_f', 'trms_t0', 'trms_t50', 'trms_t100',
        'phi_max_t0', 'phi_max_t50', 'phi_max_t100',
        'status',
        'n_clusters_raw', 'n_big_clusters',
        'largest_mass', 'largest_Epot', 'largest_rms_r', 'largest_P_peak',
        'largest_phi_max', 'largest_nvox',
        'total_cluster_mass', 'second_mass',
    ]

    out_path = os.path.join(BASE, "composite.tsv")
    with open(out_path, 'w') as f:
        f.write('\t'.join(out_cols) + '\n')
        for row in composite_rows:
            vals = [row.get(c, '') for c in out_cols]
            f.write('\t'.join(vals) + '\n')

    print(f"Wrote {len(composite_rows)} rows to {out_path}", file=sys.stderr)

    # -- Summary statistics --
    print("\n=== SUMMARY ===", file=sys.stderr)

    # Total distinct clusters
    total_raw = sum(int(parse_float(r.get('n_clusters_raw', '0'))) for r in composite_rows)
    total_big = sum(int(parse_float(r.get('n_big_clusters', '0'))) for r in composite_rows)
    print(f"Total raw clusters across all runs: {total_raw}", file=sys.stderr)
    print(f"Total significant clusters (mass>1, nvox>100): {total_big}", file=sys.stderr)

    # Multi-particle states
    multi = [r for r in composite_rows if int(parse_float(r.get('n_big_clusters', '0'))) >= 2]
    print(f"\nRuns with MULTIPLE bound clusters: {len(multi)}", file=sys.stderr)
    for r in multi[:10]:
        print(f"  {r['name']}: {r['n_big_clusters']} clusters, "
              f"largest={parse_float(r['largest_mass']):.1f}, "
              f"second={parse_float(r['second_mass']):.1f}, "
              f"Pint_f={r['Pint_f']}", file=sys.stderr)

    # Binding energy distribution (using E_pot of largest cluster)
    epots = [parse_float(r['largest_Epot']) for r in composite_rows if parse_float(r['largest_Epot']) != 0]
    if epots:
        print(f"\nBinding energy (E_pot of largest cluster):", file=sys.stderr)
        print(f"  min:    {min(epots):.2f}", file=sys.stderr)
        print(f"  max:    {max(epots):.2f}", file=sys.stderr)
        print(f"  median: {sorted(epots)[len(epots)//2]:.2f}", file=sys.stderr)
        print(f"  mean:   {sum(epots)/len(epots):.2f}", file=sys.stderr)

    # Correlation between phase offset sum and cluster count
    print(f"\nPhase offset vs cluster count:", file=sys.stderr)
    phase_groups = defaultdict(list)
    for r in composite_rows:
        d_sum = parse_float(r['d0']) + parse_float(r['d1']) + parse_float(r['d2'])
        phase_groups[f"{d_sum:.2f}"].append(int(parse_float(r.get('n_big_clusters', '0'))))

    # Group by delta sum ranges
    bins = [(0, 2), (2, 5), (5, 8), (8, 11), (11, 14)]
    for lo, hi in bins:
        counts = []
        for r in composite_rows:
            d_sum = parse_float(r['d0']) + parse_float(r['d1']) + parse_float(r['d2'])
            if lo <= d_sum < hi:
                counts.append(int(parse_float(r.get('n_big_clusters', '0'))))
        if counts:
            avg = sum(counts) / len(counts)
            multi_frac = sum(1 for c in counts if c >= 2) / len(counts) * 100
            print(f"  d_sum in [{lo},{hi}): n={len(counts)}, "
                  f"avg_clusters={avg:.2f}, multi_pct={multi_frac:.1f}%", file=sys.stderr)

    # P_int final distribution
    pint_fs = [parse_float(r.get('Pint_f', '0')) for r in composite_rows]
    if pint_fs:
        print(f"\nP_int final distribution:", file=sys.stderr)
        print(f"  min:    {min(pint_fs):.2f}", file=sys.stderr)
        print(f"  max:    {max(pint_fs):.2f}", file=sys.stderr)
        print(f"  median: {sorted(pint_fs)[len(pint_fs)//2]:.2f}", file=sys.stderr)

    # Top 10 by P_int final
    print(f"\nTop 10 runs by P_int_final:", file=sys.stderr)
    for i, r in enumerate(composite_rows[:10]):
        print(f"  {i+1}. {r['name']}: Pint_f={r['Pint_f']}, "
              f"n_big={r['n_big_clusters']}, mass={parse_float(r['largest_mass']):.1f}, "
              f"Epot={parse_float(r['largest_Epot']):.1f}", file=sys.stderr)

    # Bottom 10 by P_int final (most dissipated)
    print(f"\nBottom 10 runs by P_int_final (most dissipated):", file=sys.stderr)
    for i, r in enumerate(composite_rows[-10:]):
        print(f"  {121-9+i}. {r['name']}: Pint_f={r['Pint_f']}, "
              f"n_big={r['n_big_clusters']}, mass={parse_float(r['largest_mass']):.1f}", file=sys.stderr)

    # P_int retention ratio
    retentions = []
    for r in composite_rows:
        pi = parse_float(r.get('Pint_i', '0'))
        pf = parse_float(r.get('Pint_f', '0'))
        if pi > 1:
            retentions.append(pf / pi)
    if retentions:
        retentions.sort()
        print(f"\nP_int retention (Pint_f/Pint_i):", file=sys.stderr)
        print(f"  min:    {retentions[0]:.4f}", file=sys.stderr)
        print(f"  max:    {retentions[-1]:.4f}", file=sys.stderr)
        print(f"  median: {retentions[len(retentions)//2]:.4f}", file=sys.stderr)

    # Phase per-axis analysis
    print(f"\nPer-axis phase effect on P_int retention:", file=sys.stderr)
    for axis, key in [(0, 'd0'), (1, 'd1'), (2, 'd2')]:
        phase_vals = set(parse_float(r[key]) for r in composite_rows)
        for pv in sorted(phase_vals):
            subset = [r for r in composite_rows if abs(parse_float(r[key]) - pv) < 0.01]
            rets = []
            for r in subset:
                pi = parse_float(r.get('Pint_i', '0'))
                pf = parse_float(r.get('Pint_f', '0'))
                if pi > 1:
                    rets.append(pf / pi)
            if rets:
                print(f"  d{axis}={pv:.4f}: n={len(rets)}, "
                      f"mean_retention={sum(rets)/len(rets):.4f}", file=sys.stderr)


if __name__ == '__main__':
    main()
