# Experiment Tracking and Prediction System

## Purpose
Track all experiments across V21-V29+ in a structured database.
Use feature-outcome correlations to predict which unexplored regions
of physics space are most promising.

## Architecture

```
tracking/
├── README.md          ← this file
├── experiments.tsv    ← master database (all experiments)
├── features.md        ← feature space definition
├── predict.py         ← surrogate model + suggestion engine
└── results/
    └── cycle_N.md     ← predictions per cycle
```

## The Feature Space (Hypothesis Axes)

Each experiment is described by features (inputs) and outcomes (outputs).

### Input Features (what we choose)
| Feature | Type | Values |
|---------|------|--------|
| n_fields | int | 1, 2, 3 |
| field_type | cat | real, complex, su2, o3 |
| mass_type | cat | zero, fixed, emergent |
| mass_value | float | 0.0 - 2.0 |
| potential | cat | none, phi4, triple, skyrme, gauge |
| mu | float | -100 to 0 |
| kappa | float | 0 to 100 |
| coupling | cat | none, pairwise, elastic, metric, gauge |
| topology_init | cat | gaussian, helical, borromean, ring, random |
| phase_structure | cat | symmetric_120, bimodal_085, 180, random |
| ellipticity | float | 0.0 - 1.0 |
| bc_type | cat | absorbing, periodic, open |
| dimension | int | 1, 2, 3 |
| N_grid | int | 48 - 256 |

### Output Observables (what we measure)
| Observable | Type | Target |
|------------|------|--------|
| stable | bool | true |
| fc | float | > 0.5 |
| trans_l2 | float | > 0.2 |
| torsion_flux | float | > 0.5 |
| winding_conserved | bool | true |
| peak_P | float | > 0.1 |
| long_range_force | cat | none, yukawa, 1/r |
| force_sign | cat | attract, repel, none |
| mass_emergent | bool | true |
| spin2_radiation | bool | true |
| gauge_invariance | bool | true |
| thermal_equilibrium | bool | true |
| radiation_rate | cat | zero, decaying, constant, growing |

## The Prediction Pipeline

### Step 1: Encode
Convert each experiment (V21-V29 tests) into a feature vector + outcome vector.
Store in experiments.tsv.

### Step 2: Correlate
Compute feature-outcome correlations (Pearson for continuous, chi² for categorical).
Identify which features drive which outcomes.

### Step 3: Predict
For unexplored feature combinations:
- k-NN: find nearest explored experiments, predict outcomes from their results
- Random Forest: train on explored experiments, predict fitness of unexplored
- Constraint propagation: which feature combos are IMPOSSIBLE given known results?

### Step 4: Suggest
Rank unexplored experiments by:
- Predicted fitness (exploit known-good regions)
- Uncertainty (explore unknown regions)
- Information gain (which experiment would most reduce uncertainty?)

Balance exploration vs exploitation using UCB (Upper Confidence Bound):
  score = predicted_fitness + β × uncertainty

### Step 5: Execute
Run top-ranked experiments, add results to database, repeat.

## Tools Needed

### Minimum viable (start here):
- `experiments.tsv`: hand-populated from V21-V29 results
- `predict.py`: Python script with pandas + scikit-learn
  - Load TSV, compute correlations, train random forest
  - Generate candidate experiments, score them, rank

### Growth path:
- SQLite database (replaces TSV when it gets large)
- Gaussian Process surrogate (better uncertainty quantification)
- Active learning loop (automated experiment selection)
- Symbolic regression (discover governing equations from data)
- Graph neural network (learn over hypothesis space structure)

### What we need to prep:
1. Python environment with: numpy, pandas, scikit-learn, matplotlib
2. Encode all V21-V29 experiments (~50 rows) into the feature schema
3. Write predict.py with basic correlation + random forest
4. Generate first round of predictions for Cycle 1
