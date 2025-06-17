package main

import (
	"fmt"
	"log"
	"math/rand"
	"os"

	"github.com/hajimehoshi/ebiten/v2"
	"github.com/hajimehoshi/ebiten/v2/ebitenutil"
	"github.com/hajimehoshi/ebiten/v2/inpututil"
)

// --- SIMULATION CORE ---
// This section is a direct port of the previous simulation logic.

// Field represents our 3D circular universe.
type Field struct {
	Width, Height, Depth int
	// Double-buffered state arrays for stable calculation
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
		DiffusionRate: 0.01,
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
		if rand.Float64() < 0.005 { // Sparsely seed the field
			f.Potential[i] = rand.Float64()*2 - 1 // Random value in [-1, 1]
		}
	}
}

// Step advances the simulation by one time tick.
func (f *Field) Step() {
	size := f.Width * f.Height * f.Depth

	for z := 0; z < f.Depth; z++ {
		for y := 0; y < f.Height; y++ {
			for x := 0; x < f.Width; x++ {
				idx := f.getIndex(x, y, z)

				p_c := f.Potential[idx]
				i_xp := f.getIndex(x+1, y, z)
				i_xn := f.getIndex(x-1, y, z)
				i_yp := f.getIndex(x, y+1, z)
				i_yn := f.getIndex(x, y-1, z)
				i_zp := f.getIndex(x, y, z+1)
				i_zn := f.getIndex(x, y, z-1)
				p_xp := f.Potential[i_xp]
				p_xn := f.Potential[i_xn]
				p_yp := f.Potential[i_yp]
				p_yn := f.Potential[i_yn]
				p_zp := f.Potential[i_zp]
				p_zn := f.Potential[i_zn]

				// 1. Update Flow vector
				gradP_x := (p_xp - p_xn) / 2.0
				gradP_y := (p_yp - p_yn) / 2.0
				gradP_z := (p_zp - p_zn) / 2.0

				laplacianFlowX := f.FlowX[i_xp] + f.FlowX[i_xn] + f.FlowX[i_yp] + f.FlowX[i_yn] + f.FlowX[i_zp] + f.FlowX[i_zn] - 6*f.FlowX[idx]
				laplacianFlowY := f.FlowY[i_xp] + f.FlowY[i_xn] + f.FlowY[i_yp] + f.FlowY[i_yn] + f.FlowY[i_zp] + f.FlowY[i_zn] - 6*f.FlowY[idx]
				laplacianFlowZ := f.FlowZ[i_xp] + f.FlowZ[i_xn] + f.FlowZ[i_yp] + f.FlowZ[i_yn] + f.FlowZ[i_zp] + f.FlowZ[i_zn] - 6*f.FlowZ[idx]

				f.FlowXNext[idx] = f.FlowX[idx] - f.PressureForce*gradP_x + f.Viscosity*laplacianFlowX
				f.FlowYNext[idx] = f.FlowY[idx] - f.PressureForce*gradP_y + f.Viscosity*laplacianFlowY
				f.FlowZNext[idx] = f.FlowZ[idx] - f.PressureForce*gradP_z + f.Viscosity*laplacianFlowZ

				// 2. Update Potential
				laplacianPotential := (p_xp + p_xn + p_yp + p_yn + p_zp + p_zn) - 6*p_c

				adv_x := (p_xp*f.FlowX[i_xp] - p_xn*f.FlowX[i_xn]) / 2.0
				adv_y := (p_yp*f.FlowY[i_yp] - p_yn*f.FlowY[i_yn]) / 2.0
				adv_z := (p_zp*f.FlowZ[i_zp] - p_zn*f.FlowZ[i_zn]) / 2.0
				advection := adv_x + adv_y + adv_z

				next := p_c + f.Dt*(f.DiffusionRate*laplacianPotential-f.AdvectionRate*advection)

				if next > 1.0 {
					next = 1.0
				}
				if next < -1.0 {
					next = -1.0
				}
				f.PotentialNext[idx] = next
			}
		}
	}

	var totalPotential float64
	for i := 0; i < size; i++ {
		totalPotential += f.Potential[i]
	}
	var nextTotalPotential float64
	for i := 0; i < size; i++ {
		nextTotalPotential += f.PotentialNext[i]
	}

	if nextTotalPotential != 0 {
		correctionFactor := totalPotential / nextTotalPotential
		for i := 0; i < size; i++ {
			f.PotentialNext[i] *= correctionFactor
		}
	}

	f.Potential, f.PotentialNext = f.PotentialNext, f.Potential
	f.FlowX, f.FlowXNext = f.FlowXNext, f.FlowX
	f.FlowY, f.FlowYNext = f.FlowYNext, f.FlowY
	f.FlowZ, f.FlowZNext = f.FlowZNext, f.FlowZ
}

// TotalPotential calculates the sum of all potential in the field.
func (f *Field) TotalPotential() float64 {
	var total float64
	for _, v := range f.Potential {
		total += v
	}
	return total
}

// --- EBITENGINE APPLICATION ---

const (
	div          = 2
	screenWidth  = 640
	screenHeight = 640
	simWidth     = 640 / div
	simHeight    = 640 / div
	simDepth     = 640 / 10 / div
	axisX        = 0
	axisY        = 1
	axisZ        = 2
)

// Game implements ebiten.Game interface.
type Game struct {
	simField   *Field
	sliceImage *ebiten.Image // The image for the current slice
	sliceAxis  int
	sliceIndex int
}

// NewGame is the constructor for our application.
func NewGame() *Game {
	return &Game{
		simField:   NewField(simWidth, simHeight, simDepth),
		sliceAxis:  axisZ,
		sliceIndex: simDepth / 2,
	}
}

// Update proceeds the game state. It is called every tick.
func (g *Game) Update() error {
	// Handle Input to change slice view
	if inpututil.IsKeyJustPressed(ebiten.KeyZ) {
		g.sliceAxis = axisZ
		g.sliceIndex = g.simField.Depth / 2
	}
	if inpututil.IsKeyJustPressed(ebiten.KeyY) {
		g.sliceAxis = axisY
		g.sliceIndex = g.simField.Height / 2
	}
	if inpututil.IsKeyJustPressed(ebiten.KeyX) {
		g.sliceAxis = axisX
		g.sliceIndex = g.simField.Width / 2
	}
	if inpututil.IsKeyJustPressed(ebiten.KeyQ) {
		os.Exit(0)
	}

	shift := ebiten.IsKeyPressed(ebiten.KeyShift)
	if shift {
		switch {
		case ebiten.IsKeyPressed(ebiten.KeyArrowUp):
			g.sliceIndex++

		case ebiten.IsKeyPressed(ebiten.KeyArrowDown):
			g.sliceIndex--
		}
	} else {
		switch {
		case inpututil.IsKeyJustPressed(ebiten.KeyArrowUp):
			g.sliceIndex++
		case inpututil.IsKeyJustPressed(ebiten.KeyArrowDown):
			g.sliceIndex--
		}
	}

	// Clamp sliceIndex to be within bounds for the current axis
	maxIndex := 0
	switch g.sliceAxis {
	case axisZ:
		maxIndex = g.simField.Depth - 1
	case axisY:
		maxIndex = g.simField.Height - 1
	case axisX:
		maxIndex = g.simField.Width - 1
	}
	if g.sliceIndex > maxIndex {
		g.sliceIndex = 0
	}
	if g.sliceIndex < 0 {
		g.sliceIndex = maxIndex
	}

	// Run one step of the simulation
	g.simField.Step()
	return nil
}

// Draw draws the game screen. It is called every frame.
func (g *Game) Draw(screen *ebiten.Image) {
	// 1. Determine the dimensions of the slice we are viewing
	viewWidth, viewHeight := 0, 0
	switch g.sliceAxis {
	case axisZ:
		viewWidth, viewHeight = g.simField.Width, g.simField.Height
	case axisY:
		viewWidth, viewHeight = g.simField.Width, g.simField.Depth
	case axisX:
		viewWidth, viewHeight = g.simField.Height, g.simField.Depth
	}

	// 2. Create or resize the image buffer if necessary
	if g.sliceImage == nil || g.sliceImage.Bounds().Dx() != viewWidth || g.sliceImage.Bounds().Dy() != viewHeight {
		g.sliceImage = ebiten.NewImage(viewWidth, viewHeight)
	}
	pixels := make([]byte, viewWidth*viewHeight*4)

	// 3. Populate the pixel buffer based on the slice
	f := g.simField
	for j := 0; j < viewHeight; j++ {
		for i := 0; i < viewWidth; i++ {
			var potential float64
			switch g.sliceAxis {
			case axisZ:
				potential = f.Potential[f.getIndex(i, j, g.sliceIndex)]
			case axisY:
				potential = f.Potential[f.getIndex(i, g.sliceIndex, j)]
			case axisX:
				potential = f.Potential[f.getIndex(g.sliceIndex, i, j)]
			}

			// Map potential [-1, 1] to a color
			val := (potential + 1.0) * 0.5 // Map to [0, 1]
			idx := (j*viewWidth + i) * 4
			if val < 0.5 { // Negative potential -> Blue
				pixels[idx+0] = 0                          // R
				pixels[idx+1] = 0                          // G
				pixels[idx+2] = uint8(255 * (1.0 - val*2)) // B
			} else { // Positive potential -> Red
				pixels[idx+0] = uint8(255 * (val - 0.5) * 2) // R
				pixels[idx+1] = 0                            // G
				pixels[idx+2] = 0                            // B
			}
			pixels[idx+3] = 255 // Alpha
		}
	}
	g.sliceImage.WritePixels(pixels)

	// 4. Draw the slice image to the screen, scaled up
	op := &ebiten.DrawImageOptions{}
	scaleX := float64(screenWidth) / float64(viewWidth)
	scaleY := float64(screenHeight) / float64(viewHeight)
	scale := scaleX
	if scaleY < scaleX {
		scale = scaleY
	}
	op.GeoM.Scale(scale, scale)
	screen.DrawImage(g.sliceImage, op)

	// 5. Draw UI text
	axisStr := ""
	switch g.sliceAxis {
	case axisX:
		axisStr = "X"
	case axisY:
		axisStr = "Y"
	case axisZ:
		axisStr = "Z"
	}

	msg := fmt.Sprintf(
		"CONTROLS: Use X, Y, Z keys to change axis | Up/Down arrows to change slice\n"+
			"Slicing Axis: %s | Index: %d\n"+
			"Total Potential: %e",
		axisStr, g.sliceIndex, g.simField.TotalPotential(),
	)
	ebitenutil.DebugPrint(screen, msg)
}

// Layout takes the outside size (e.g. window size) and returns the logical screen size.
func (g *Game) Layout(outsideWidth, outsideHeight int) (int, int) {
	return screenWidth, screenHeight
}

func main() {
	game := NewGame()
	ebiten.SetWindowSize(screenWidth, screenHeight)
	ebiten.SetWindowTitle("Native Field F Simulator")
	if err := ebiten.RunGame(game); err != nil {
		log.Fatal(err)
	}
}
