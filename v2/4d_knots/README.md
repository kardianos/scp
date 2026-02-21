# 4D Knot Visualizer

Interactive visualizer for knots in 4D space where knots can slip through each other.

## Run

```bash
~/.cargo/bin/cargo run --release
```

## Controls

- **Drag**: Rotate camera
- **Scroll**: Zoom
- **Space**: Pause/Resume simulation
- **W-Slice slider**: View different "slices" of the 4th dimension

## Physics

Knots exist in 4D space but are visualized as 3D projections. The 4th dimension (W) allows knots to:
- Pass through each other when W-phases differ
- Become trapped when W-phases align (energy barrier)
- Slide freely along the field (binding energy well)

The energy landscape combines:
- **Energy barrier**: Sharp Gaussian barrier when knots overlap in W-space (activation energy to escape)
- **Binding strength**: Broader Gaussian well that traps knots at similar W-phases

This models topological solitons where the extra dimension provides the freedom for continuous deformation without intersection.

## Extending

- Add volumetric field visualization (energy density clouds)
- Add more knot types
- Implement proper 4D rotation matrices
- Add energy export for analysis
