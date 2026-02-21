# Field Knot Simulation

A 3D+time scalar field simulation featuring emergent topological structures through phase transitions.

## Concept

**Field knots** are topological defects that emerge spontaneously from a field undergoing a phase transition (Kibble-Zurek mechanism). Unlike placed parametric curves, these structures **self-organize** from random initial conditions.

## Features

- **Emergent knot formation** via phase transition (cooling)
- **Real-time visualization** of field density, phase, and defects
- **Density gradient dynamics** and field consumption
- **Boson-like transfers** between coupled structures
- **Stability analysis** through energy metrics

## Running

### Rust Visualizer (Real-time)
```bash
cd field_knot
cargo run --release
```

### Python Simulations (Analysis)
```bash
cd field_knot
python sim3_emergent_knots.py   # Emergent knot formation
python sim2_boson_transfer.py   # Energy transfer dynamics
python sim1_density_dynamics.py # Density gradients & consumption
python coupling_analysis.py     # Coupling stability analysis
```

## Controls (Rust Visualizer)

| Key | Action |
|-----|--------|
| Space | Pause/Resume |
| R | Reset simulation |
| D | Toggle field density |
| F | Toggle defects (knots) |
| W | Toggle phase winding |
| Mouse drag | Rotate camera |
| Scroll | Zoom |

## Visualization Modes

1. **Defects (Knots)**: Red/orange points showing phase singularities
2. **Field Density**: Blue/purple cloud showing |Φ|²
3. **Phase Winding**: Hot spots showing |∇φ| (where defects form)

## Mathematical Framework

### Phase Transition Dynamics
```
∂Φ/∂t = D·∇²Φ + (1 - |Φ|²)·Φ - T·Φ + noise
```

Where T is temperature that decreases during cooling.

### Defect Detection
Defects are detected where phase has non-zero winding:
```
W = |∇φ| > threshold  →  defect location
```

### Key Equations

| Quantity | Formula |
|----------|---------|
| Density | ρ = |Φ|² |
| Phase | φ = arg(Φ) |
| Winding | W = Σ |Δφ| |
| Energy | E = ∫ρ dV |

## Files

| File | Description |
|------|-------------|
| `src/main.rs` | Rust visualizer with emergent simulation |
| `sim3_emergent_knots.py` | Phase transition knot formation |
| `sim2_boson_transfer.py` | Energy transfer patterns |
| `sim1_density_dynamics.py` | Field consumption dynamics |
| `coupling_analysis.py` | Coupling stability analysis |
| `THEORY.md` | Mathematical foundations |

## Key Parameters

| Parameter | Effect |
|-----------|--------|
| Temperature (T) | Controls disorder → order transition |
| Diffusion | Smoothing and correlation length |
| Noise | Thermal fluctuations |
| Density threshold | Visualization cutoff |

## Observations

1. **Knot Formation**: Defects appear during cooling when T drops below critical value
2. **Defect Count**: Initially high (~350), stabilizes at ~80-100 persistent knots
3. **Self-Reinforcement**: Stable structures maintain themselves through phase coherence
4. **Energy Flow**: "Bosons" transfer energy between coupled structures
