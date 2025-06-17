// main.go
package main

import (
	"math/rand"
	"syscall/js"
)

// Field represents our 3D circular universe.
type Field struct {
	Width, Height, Depth int
	// We use two buffers for each property to calculate the next state
	// based on the current state without race conditions.
	Potential, PotentialNext []float64
	FlowX, FlowXNext         []float64
	FlowY, FlowYNext         []float64
	FlowZ, FlowZNext         []float64

	// Simulation parameters
	Dt            float64 // Time step
	DiffusionRate float64 // How fast potential spreads
	AdvectionRate float64 // How much flow carries potential
	PressureForce float64 // How much potential gradient creates flow
	Viscosity     float64 // How fast flow spreads/dampens
}

// NewField creates and initializes a new field.
func NewField(w, h, d int) *Field {
	size := w * h * d
	f := &Field{
		Width:  w,
		Height: h,
		Depth:  d,
		// State buffers
		Potential:     make([]float64, size),
		PotentialNext: make([]float64, size),
		FlowX:         make([]float64, size),
		FlowXNext:     make([]float64, size),
		FlowY:         make([]float64, size),
		FlowYNext:     make([]float64, size),
		FlowZ:         make([]float64, size),
		FlowZNext:     make([]float64, size),
		// Tunable constants
		Dt:            0.1,
		DiffusionRate: 0.02,
		AdvectionRate: 0.1,
		PressureForce: 0.5,
		Viscosity:     0.01,
	}
	f.randomize()
	return f
}

// getIndex calculates the 1D array index for a 3D coordinate.
func (f *Field) getIndex(x, y, z int) int {
	// Periodic (toroidal) boundary conditions
	x = (x + f.Width) % f.Width
	y = (y + f.Height) % f.Height
	z = (z + f.Depth) % f.Depth
	return (z*f.Height+y)*f.Width + x
}

// randomize sets up initial chaotic conditions.
func (f *Field) randomize() {
	for i := range f.Potential {
		if rand.Float64() < 0.01 { // Sparsely seed the field
			f.Potential[i] = rand.Float64()*2 - 1 // Random value in [-1, 1]
		}
	}
}

// Step advances the simulation by one time tick.
func (f *Field) Step() {
	size := f.Width * f.Height * f.Depth

	// --- Main Simulation Loop ---
	for z := 0; z < f.Depth; z++ {
		for y := 0; y < f.Height; y++ {
			for x := 0; x < f.Width; x++ {
				idx := f.getIndex(x, y, z)

				// Get potential and flow of neighbors
				p_c := f.Potential[idx]
				p_xp := f.Potential[f.getIndex(x+1, y, z)]
				p_xn := f.Potential[f.getIndex(x-1, y, z)]
				p_yp := f.Potential[f.getIndex(x, y+1, z)]
				p_yn := f.Potential[f.getIndex(x, y-1, z)]
				p_zp := f.Potential[f.getIndex(x, y, z+1)]
				p_zn := f.Potential[f.getIndex(x, y, z-1)]

				// 1. Update Flow vector (our "electroweak force")
				// a. Flow is pushed by the negative gradient of Potential ("pressure/gravity")
				gradP_x := (p_xp - p_xn) / 2.0
				gradP_y := (p_yp - p_yn) / 2.0
				gradP_z := (p_zp - p_zn) / 2.0

				// b. Flow diffuses (viscosity)
				laplacianFlowX := f.FlowX[f.getIndex(x+1, y, z)] + f.FlowX[f.getIndex(x-1, y, z)] + f.FlowX[f.getIndex(x, y+1, z)] + f.FlowX[f.getIndex(x, y-1, z)] + f.FlowX[f.getIndex(x, y, z+1)] + f.FlowX[f.getIndex(x, y, z-1)] - 6*f.FlowX[idx]
				laplacianFlowY := f.FlowY[f.getIndex(x+1, y, z)] + f.FlowY[f.getIndex(x-1, y, z)] + f.FlowY[f.getIndex(x, y+1, z)] + f.FlowY[f.getIndex(x, y-1, z)] + f.FlowY[f.getIndex(x, y, z+1)] + f.FlowY[f.getIndex(x, y, z-1)] - 6*f.FlowY[idx]
				laplacianFlowZ := f.FlowZ[f.getIndex(x+1, y, z)] + f.FlowZ[f.getIndex(x-1, y, z)] + f.FlowZ[f.getIndex(x, y+1, z)] + f.FlowZ[f.getIndex(x, y-1, z)] + f.FlowZ[f.getIndex(x, y, z+1)] + f.FlowZ[f.getIndex(x, y, z-1)] - 6*f.FlowZ[idx]

				f.FlowXNext[idx] = f.FlowX[idx] - f.PressureForce*gradP_x + f.Viscosity*laplacianFlowX
				f.FlowYNext[idx] = f.FlowY[idx] - f.PressureForce*gradP_y + f.Viscosity*laplacianFlowY
				f.FlowZNext[idx] = f.FlowZ[idx] - f.PressureForce*gradP_z + f.Viscosity*laplacianFlowZ

				// 2. Update Potential (our "density")
				// a. Potential diffuses
				laplacianPotential := (p_xp + p_xn + p_yp + p_yn + p_zp + p_zn) - 6*p_c

				// b. Potential is carried by flow (advection)
				// Divergence of (Potential * Flow)
				adv_x := (p_xp*f.FlowX[f.getIndex(x+1, y, z)] - p_xn*f.FlowX[f.getIndex(x-1, y, z)]) / 2.0
				adv_y := (p_yp*f.FlowY[f.getIndex(x, y+1, z)] - p_yn*f.FlowY[f.getIndex(x, y-1, z)]) / 2.0
				adv_z := (p_zp*f.FlowZ[f.getIndex(x, y, z+1)] - p_zn*f.FlowZ[f.getIndex(x, y, z-1)]) / 2.0
				advection := adv_x + adv_y + adv_z

				f.PotentialNext[idx] = p_c + f.Dt*(f.DiffusionRate*laplacianPotential-f.AdvectionRate*advection)

				// Clamp potential to [-1, 1] range to maintain stability
				if f.PotentialNext[idx] > 1.0 {
					f.PotentialNext[idx] = 1.0
				}
				if f.PotentialNext[idx] < -1.0 {
					f.PotentialNext[idx] = -1.0
				}
			}
		}
	}

	// --- Conservation and Buffer Swap ---
	var totalPotential float64 = 0
	for i := 0; i < size; i++ {
		totalPotential += f.PotentialNext[i]
	}
	// Renormalize to conserve total potential (corrects for float drift)
	var initialPotential float64 = 0 // should be ~0 for random start
	for i := 0; i < size; i++ {
		initialPotential += f.Potential[i]
	}
	if totalPotential != 0 {
		correctionFactor := initialPotential / totalPotential
		for i := 0; i < size; i++ {
			f.PotentialNext[i] *= correctionFactor
		}
	}

	// Swap buffers for next iteration
	f.Potential, f.PotentialNext = f.PotentialNext, f.Potential
	f.FlowX, f.FlowXNext = f.FlowXNext, f.FlowX
	f.FlowY, f.FlowYNext = f.FlowYNext, f.FlowY
	f.FlowZ, f.FlowZNext = f.FlowZNext, f.FlowZ
}

// --- WebAssembly Interface ---

var field *Field

func jsInitField(this js.Value, args []js.Value) interface{} {
	w, h, d := args[0].Int(), args[1].Int(), args[2].Int()
	field = NewField(w, h, d)
	return nil
}

func jsRunStep(this js.Value, args []js.Value) interface{} {
	if field != nil {
		field.Step()
	}
	return nil
}

func jsGetPotentialSlice(this js.Value, args []js.Value) interface{} {
	if field == nil {
		return nil
	}
	z := args[0].Int()
	sliceSize := field.Width * field.Height
	start := field.getIndex(0, 0, z)

	// Create a JS-managed buffer and copy the data into it
	jsBuffer := js.Global().Get("Uint8Array").New(sliceSize)
	js.CopyBytesToJS(jsBuffer, toUint8(field.Potential[start:start+sliceSize]))
	return jsBuffer
}

// Calculates the sum of all potential in the field.
func jsGetTotalPotential(this js.Value, args []js.Value) interface{} {
	if field == nil {
		return 0.0
	}
	var total float64
	for _, v := range field.Potential {
		total += v
	}
	return total
}

// Helper to convert float64 slice to byte slice for JS copy
func toUint8(slice []float64) []byte {
	out := make([]byte, len(slice))
	for i, v := range slice {
		// Map [-1, 1] to [0, 255]
		out[i] = byte((v + 1.0) * 0.5 * 255.0)
	}
	return out
}

func main() {
	c := make(chan struct{}, 0)
	js.Global().Set("initField", js.FuncOf(jsInitField))
	js.Global().Set("runStep", js.FuncOf(jsRunStep))
	js.Global().Set("getPotentialSlice", js.FuncOf(jsGetPotentialSlice))
	// Export the new function
	js.Global().Set("getTotalPotential", js.FuncOf(jsGetTotalPotential))
	<-c
}
