#!/usr/bin/env python3
"""
predict.py — Experiment tracking and prediction for SCP braided soliton project.

Loads experiments.tsv, computes feature-outcome correlations,
trains a surrogate model, and suggests next experiments.

Usage:
    python3 predict.py                  # full analysis
    python3 predict.py --correlations   # just correlations
    python3 predict.py --suggest N      # suggest N next experiments
"""

import sys
import numpy as np
import pandas as pd
from pathlib import Path

# ============================================================
# Feature definitions
# ============================================================

# Continuous features
CONT_FEATURES = ['mass_value', 'mu', 'kappa', 'ellipticity']

# Categorical features (will be one-hot encoded)
CAT_FEATURES = ['field_type', 'potential', 'coupling', 'topology',
                'phase_struct', 'bc']

# Key outcomes to predict
OUTCOMES = ['stable', 'fc', 'trans_l2', 'torsion', 'winding_ok',
            'peak_P', 'thermal_eq']

# ============================================================
# Load and encode
# ============================================================

def load_data(path='experiments.tsv'):
    df = pd.read_csv(path, sep='\t')
    return df

def encode_features(df):
    """Convert categorical + continuous features into numeric matrix."""
    parts = []

    # Continuous
    for f in CONT_FEATURES:
        if f in df.columns:
            parts.append(df[f].astype(float).values.reshape(-1, 1))

    # Categorical (one-hot)
    for f in CAT_FEATURES:
        if f in df.columns:
            dummies = pd.get_dummies(df[f], prefix=f)
            parts.append(dummies.values)

    if not parts:
        return np.zeros((len(df), 0))
    return np.hstack(parts)

def get_feature_names(df):
    names = list(CONT_FEATURES)
    for f in CAT_FEATURES:
        if f in df.columns:
            for val in df[f].unique():
                names.append(f'{f}_{val}')
    return names

# ============================================================
# Correlation analysis
# ============================================================

def compute_correlations(df):
    """Compute Pearson correlation between each feature and each outcome."""
    print("\n=== Feature-Outcome Correlations ===\n")

    X = encode_features(df)
    names = get_feature_names(df)

    # Header
    header = f"{'Feature':<25s}"
    for o in OUTCOMES:
        header += f" {o:>10s}"
    print(header)
    print("-" * len(header))

    for i, name in enumerate(names):
        if i >= X.shape[1]:
            break
        row = f"{name:<25s}"
        for o in OUTCOMES:
            if o in df.columns:
                y = pd.to_numeric(df[o], errors='coerce').fillna(0).values
                if np.std(X[:, i]) > 1e-10 and np.std(y) > 1e-10:
                    r = np.corrcoef(X[:, i], y)[0, 1]
                    row += f" {r:+10.3f}"
                else:
                    row += f" {'---':>10s}"
            else:
                row += f" {'N/A':>10s}"
        print(row)

# ============================================================
# Suggestion engine
# ============================================================

# Define the UNEXPLORED hypothesis space
CANDIDATES = [
    # Phase B: Field studies
    {'id': 'T10A_scatter', 'phase': 'B', 'desc': 'Wave scattering in triple-product field',
     'mass_value': 1.5, 'mu': -41.3, 'kappa': 50, 'topology': 'wave_packet',
     'rationale': 'Test if V(P) creates phase-dependent scattering (→EM at field level)'},

    {'id': 'T10B_metric', 'phase': 'B', 'desc': 'Test wave propagation through braid background',
     'mass_value': 1.5, 'mu': -41.3, 'kappa': 50, 'topology': 'test_wave',
     'rationale': 'Does braid curve the field? (→gravity at field level)'},

    {'id': 'T10C_dispersion', 'phase': 'B', 'desc': 'Dispersion relation around condensate',
     'mass_value': 0.0, 'mu': -41.3, 'kappa': 50, 'topology': 'linear',
     'rationale': 'Can mass emerge from field configuration? Band structure?'},

    {'id': 'T10D_kibble', 'phase': 'B', 'desc': 'Cool random field, look for spontaneous defects',
     'mass_value': 1.5, 'mu': -41.3, 'kappa': 50, 'topology': 'random',
     'rationale': 'Do braids form naturally? Kibble-Zurek mechanism'},

    {'id': 'T10E_greens', 'phase': 'B', 'desc': 'Impulse response / Green\'s function',
     'mass_value': 1.5, 'mu': -41.3, 'kappa': 50, 'topology': 'impulse',
     'rationale': 'Spectral function reveals what propagates (spin-2? spin-1?)'},

    # Phase A: Substrate
    {'id': 'T9A_cubic', 'phase': 'A', 'desc': 'Cubic lattice + 3-body → coarse grain',
     'mass_value': 0, 'mu': 0, 'kappa': 0, 'topology': 'lattice',
     'rationale': 'Does cubic lattice produce V(det F) under Cauchy-Born?'},

    {'id': 'T9A_fcc', 'phase': 'A', 'desc': 'FCC lattice + 3-body → coarse grain',
     'mass_value': 0, 'mu': 0, 'kappa': 0, 'topology': 'lattice',
     'rationale': 'Close-packed lattice: different symmetry, different continuum limit'},

    {'id': 'T9C_expand', 'phase': 'A', 'desc': 'Dense field rapid expansion',
     'mass_value': 1.5, 'mu': -41.3, 'kappa': 50, 'topology': 'random',
     'rationale': 'Expansion + cooling → spontaneous braid formation?'},

    # Phase C: Particle refinements
    {'id': 'T7b_close', 'phase': 'C', 'desc': 'Two braids at D=10 (inside Yukawa range)',
     'mass_value': 1.5, 'mu': -41.3, 'kappa': 50, 'topology': 'helical_x2',
     'rationale': 'T7 repulsion was at D=30 (beyond range). D=10 might attract'},

    {'id': 'T_su2', 'phase': 'C', 'desc': 'Upgrade to SU(2) field (true topology)',
     'mass_value': 1.5, 'mu': -41.3, 'kappa': 50, 'topology': 'skyrmion',
     'rationale': 'π₃(SU(2))=Z gives REAL topological protection'},

    {'id': 'T_m_scan_bimodal', 'phase': 'C', 'desc': 'Fine mass scan m=0.5-2.0 with bimodal geometry',
     'mass_value': 'scan', 'mu': -41.3, 'kappa': 50, 'topology': 'helical',
     'rationale': 'Find optimal mass that maximizes BOTH trans + torsion'},

    {'id': 'T_4field', 'phase': 'C', 'desc': 'Four scalar fields with 4-product V(φ₀φ₁φ₂φ₃)',
     'mass_value': 1.5, 'mu': -41.3, 'kappa': 50, 'topology': 'helical',
     'rationale': '4 fields → richer topology, possible spin-2 from 4-index tensor'},
]

def suggest_experiments(n=5):
    """Score and rank candidate experiments."""
    print(f"\n=== Top {n} Suggested Experiments ===\n")

    # Simple scoring based on phase priority + information gain
    phase_priority = {'B': 3.0, 'A': 2.0, 'C': 1.0}  # Field > Substrate > Particle

    scored = []
    for c in CANDIDATES:
        score = phase_priority.get(c['phase'], 1.0)
        # Boost experiments that address known failures
        if 'gravity' in c['rationale'].lower() or 'metric' in c['rationale'].lower():
            score += 1.0  # gravity is our biggest gap
        if 'mass' in c['rationale'].lower() and 'emerge' in c['rationale'].lower():
            score += 0.5  # mass emergence is important
        if 'topology' in c['rationale'].lower():
            score += 0.3  # topology is a known weakness
        scored.append((score, c))

    scored.sort(key=lambda x: -x[0])

    for i, (score, c) in enumerate(scored[:n]):
        print(f"  {i+1}. [{c['phase']}] {c['id']} (score={score:.1f})")
        print(f"     {c['desc']}")
        print(f"     Why: {c['rationale']}")
        print()

# ============================================================
# Summary statistics
# ============================================================

def print_summary(df):
    print(f"\n=== Database Summary ===")
    print(f"Total experiments: {len(df)}")
    print(f"Versions: {', '.join(df['version'].unique())}")
    print(f"Phases: {dict(df['phase'].value_counts())}")
    print(f"Stable: {df['stable'].sum()}/{len(df)}")

    print(f"\nBest results:")
    for o in ['fc', 'trans_l2', 'torsion', 'peak_P']:
        if o in df.columns:
            idx = pd.to_numeric(df[o], errors='coerce').idxmax()
            val = df.loc[idx, o]
            name = df.loc[idx, 'id']
            print(f"  Best {o}: {val} ({name})")

# ============================================================
# Main
# ============================================================

def main():
    df = load_data(Path(__file__).parent / 'experiments.tsv')

    if '--correlations' in sys.argv:
        compute_correlations(df)
    elif '--suggest' in sys.argv:
        n = int(sys.argv[sys.argv.index('--suggest') + 1]) if len(sys.argv) > sys.argv.index('--suggest') + 1 else 5
        suggest_experiments(n)
    else:
        print_summary(df)
        compute_correlations(df)
        suggest_experiments(8)

if __name__ == '__main__':
    main()
