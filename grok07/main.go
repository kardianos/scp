package main

import (
	"fmt"
	"log"
	"math"
	"math/rand"
	"os"

	sfft "github.com/davidkleiven/gosfft/sfft"
	"github.com/hajimehoshi/ebiten/v2"
	"github.com/hajimehoshi/ebiten/v2/ebitenutil"
	"github.com/hajimehoshi/ebiten/v2/inpututil"
	"gonum.org/v1/gonum/cmplxs"
)

// --- SIMULATION CORE ---

func sanityCheck(slice []float64, name string) {
	for i, v := range slice {
		if math.IsNaN(v) || math.IsInf(v, 0) {
			log.Fatalf("FATAL: Numerical instability detected. Field '%s' has an invalid value (%f) at index %d. Halting simulation.", name, v, i)
		}
	}
}

// FieldAnalysis holds the results of analyzing the field's state.
type FieldAnalysis struct {
	Total            float64
	Min              float64
	Max              float64
	MaxX, MaxY, MaxZ int
}

// Field represents our 3D circular universe.
type Field struct {
	Width, Height, Depth int
	// ... (rest of the fields are the same)
	Potential, PotentialNext                                                              []float64
	FlowX, FlowXNext                                                                      []float64
	FlowY, FlowYNext                                                                      []float64
	FlowZ, FlowZNext                                                                      []float64
	Ex, ExNext                                                                            []float64
	Ey, EyNext                                                                            []float64
	Ez, EzNext                                                                            []float64
	Bx, BxNext                                                                            []float64
	By, ByNext                                                                            []float64
	Bz, BzNext                                                                            []float64
	Dt, DiffusionRate, AdvectionRate, PressureForce, Viscosity, EmfCoupling, EmfWaveSpeed float64
	fft3                                                                                  *sfft.FFT3
}

// NewField creates and initializes a new field.
func NewField(w, h, d int) *Field {
	size := w * h * d
	f := &Field{
		Width: w, Height: h, Depth: d,
		Potential: make([]float64, size), PotentialNext: make([]float64, size),
		FlowX: make([]float64, size), FlowXNext: make([]float64, size),
		FlowY: make([]float64, size), FlowYNext: make([]float64, size),
		FlowZ: make([]float64, size), FlowZNext: make([]float64, size),
		Ex: make([]float64, size), ExNext: make([]float64, size),
		Ey: make([]float64, size), EyNext: make([]float64, size),
		Ez: make([]float64, size), EzNext: make([]float64, size),
		Bx: make([]float64, size), BxNext: make([]float64, size),
		By: make([]float64, size), ByNext: make([]float64, size),
		Bz: make([]float64, size), BzNext: make([]float64, size),
		Dt: 0.1, DiffusionRate: 0.005, AdvectionRate: 0.9, PressureForce: 0.6,
		Viscosity: 0.001, EmfCoupling: 0.4, EmfWaveSpeed: 1.0,
		fft3: sfft.NewFFT3(w, h, d),
	}
	f.randomize()
	return f
}

// getIndex calculates the 1D array index for a 3D coordinate.
func (f *Field) getIndex(x, y, z int) int {
	x = (x + f.Width) % f.Width
	y = (y + f.Height) % f.Height
	z = (z + f.Depth) % f.Depth
	return (z*f.Height+y)*f.Width + x
}

// randomize sets up initial chaotic conditions.
func (f *Field) randomize() {
	for i := range f.Potential {
		if rand.Float64() < 0.005 {
			f.Potential[i] = rand.Float64()*2 - 1
		}
	}
}

// fftSmooth applies 3D FFT low-pass filter to a 1D field slice to remove artifacts.
func (f *Field) fftSmooth(data []float64) []float64 {
	cdata := make([]complex128, len(data))
	for i, v := range data {
		cdata[i] = complex(v, 0)
	}
	coeff := f.fft3.FFT(cdata)

	// Low-pass: Zero high frequencies (e.g., keep lowest 1/2 in each dim)
	cutoff_kx := f.Width / 2
	cutoff_ky := f.Height / 2
	cutoff_kz := f.Depth / 2
	for z := 0; z < f.Depth; z++ {
		kz := z
		if kz > f.Depth/2 {
			kz -= f.Depth
		}
		if math.Abs(float64(kz)) > float64(cutoff_kz) {
			for y := 0; y < f.Height; y++ {
				for x := 0; x < f.Width; x++ {
					idx := (z*f.Height+y)*f.Width + x
					coeff[idx] = 0
				}
			}
			continue
		}
		for y := 0; y < f.Height; y++ {
			ky := y
			if ky > f.Height/2 {
				ky -= f.Height
			}
			if math.Abs(float64(ky)) > float64(cutoff_ky) {
				for x := 0; x < f.Width; x++ {
					idx := (z*f.Height+y)*f.Width + x
					coeff[idx] = 0
				}
				continue
			}
			for x := 0; x < f.Width; x++ {
				kx := x
				if kx > f.Width/2 {
					kx -= f.Width
				}
				if math.Abs(float64(kx)) > float64(cutoff_kx) {
					idx := (z*f.Height+y)*f.Width + x
					coeff[idx] = 0
				}
			}
		}
	}

	smoothed := f.fft3.IFFT(coeff)
	cmplxs.Scale(complex(1.0/float64(f.Width*f.Height*f.Depth), 0), smoothed)
	result := make([]float64, len(data))
	for i := range result {
		result[i] = real(smoothed[i])
	}
	return result
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

	// Apply FFT smoothing to remove voxel artifacts
	f.Potential = f.fftSmooth(f.Potential)
	f.FlowX = f.fftSmooth(f.FlowX)
	f.FlowY = f.fftSmooth(f.FlowY)
	f.FlowZ = f.fftSmooth(f.FlowZ)
}

// Step2 advances the simulation by one time tick with full EM.
func (f *Field) Step2() {
	c2 := f.EmfWaveSpeed * f.EmfWaveSpeed
	for z := 0; z < f.Depth; z++ {
		for y := 0; y < f.Height; y++ {
			for x := 0; x < f.Width; x++ {
				idx := f.getIndex(x, y, z)
				ixp, ixn := f.getIndex(x+1, y, z), f.getIndex(x-1, y, z)
				iyp, iyn := f.getIndex(x, y+1, z), f.getIndex(x, y-1, z)
				izp, izn := f.getIndex(x, y, z+1), f.getIndex(x, y, z-1)
				p_c := f.Potential[idx]
				p_xp, p_xn := f.Potential[ixp], f.Potential[ixn]
				p_yp, p_yn := f.Potential[iyp], f.Potential[iyn]
				p_zp, p_zn := f.Potential[izp], f.Potential[izn]
				curlE_x := (f.Ez[iyp] - f.Ez[iyn]) - (f.Ey[izp] - f.Ey[izn])
				curlE_y := (f.Ex[izp] - f.Ex[izn]) - (f.Ez[ixp] - f.Ez[ixn])
				curlE_z := (f.Ey[ixp] - f.Ey[ixn]) - (f.Ex[iyp] - f.Ex[iyn])
				f.BxNext[idx] = f.Bx[idx] - f.Dt*curlE_x
				f.ByNext[idx] = f.By[idx] - f.Dt*curlE_y
				f.BzNext[idx] = f.Bz[idx] - f.Dt*curlE_z
				curlB_x := (f.Bz[iyp] - f.Bz[iyn]) - (f.By[izp] - f.By[izn])
				curlB_y := (f.Bx[izp] - f.Bx[izn]) - (f.Bz[ixp] - f.Bz[ixn])
				curlB_z := (f.By[ixp] - f.By[ixn]) - (f.Bx[iyp] - f.Bx[iyn])
				f.ExNext[idx] = f.Ex[idx] + f.Dt*(c2*curlB_x-f.FlowX[idx])
				f.EyNext[idx] = f.Ey[idx] + f.Dt*(c2*curlB_y-f.FlowY[idx])
				f.EzNext[idx] = f.Ez[idx] + f.Dt*(c2*curlB_z-f.FlowZ[idx])
				gradP_x := p_xp - p_xn
				gradP_y := p_yp - p_yn
				gradP_z := p_zp - p_zn
				v_cross_B_x := f.FlowY[idx]*f.Bz[idx] - f.FlowZ[idx]*f.By[idx]
				v_cross_B_y := f.FlowZ[idx]*f.Bx[idx] - f.FlowX[idx]*f.Bz[idx]
				v_cross_B_z := f.FlowX[idx]*f.By[idx] - f.FlowY[idx]*f.Bx[idx]
				lorentz_x := p_c * (f.Ex[idx] + v_cross_B_x)
				lorentz_y := p_c * (f.Ey[idx] + v_cross_B_y)
				lorentz_z := p_c * (f.Ez[idx] + v_cross_B_z)
				laplacianFlowX := f.FlowX[ixp] + f.FlowX[ixn] + f.FlowX[iyp] + f.FlowX[iyn] + f.FlowX[izp] + f.FlowX[izn] - 6*f.FlowX[idx]
				laplacianFlowY := f.FlowY[ixp] + f.FlowY[ixn] + f.FlowY[iyp] + f.FlowY[iyn] + f.FlowY[izp] + f.FlowY[izn] - 6*f.FlowY[idx]
				laplacianFlowZ := f.FlowZ[ixp] + f.FlowZ[ixn] + f.FlowZ[iyp] + f.FlowZ[iyn] + f.FlowZ[izp] + f.FlowZ[izn] - 6*f.FlowZ[idx]
				f.FlowXNext[idx] = f.FlowX[idx] + f.Dt*(-f.PressureForce*gradP_x+f.EmfCoupling*lorentz_x+f.Viscosity*laplacianFlowX)
				f.FlowYNext[idx] = f.FlowY[idx] + f.Dt*(-f.PressureForce*gradP_y+f.EmfCoupling*lorentz_y+f.Viscosity*laplacianFlowY)
				f.FlowZNext[idx] = f.FlowZ[idx] + f.Dt*(-f.PressureForce*gradP_z+f.EmfCoupling*lorentz_z+f.Viscosity*laplacianFlowZ)
				laplacianPotential := (p_xp + p_xn + p_yp + p_yn + p_zp + p_zn) - 6*p_c
				advection := (p_xp*f.FlowX[ixp] - p_xn*f.FlowX[ixn] + p_yp*f.FlowY[iyp] - p_yn*f.FlowY[iyn] + p_zp*f.FlowZ[izp] - p_zn*f.FlowZ[izn])
				f.PotentialNext[idx] = p_c + f.Dt*(f.DiffusionRate*laplacianPotential-f.AdvectionRate*advection)
			}
		}
	}
	sanityCheck(f.PotentialNext, "PotentialNext")
	sanityCheck(f.FlowXNext, "FlowXNext")
	sanityCheck(f.FlowYNext, "FlowYNext")
	sanityCheck(f.FlowZNext, "FlowZNext")
	sanityCheck(f.ExNext, "ExNext")
	sanityCheck(f.EyNext, "EyNext")
	sanityCheck(f.EzNext, "EzNext")
	sanityCheck(f.BxNext, "BxNext")
	sanityCheck(f.ByNext, "ByNext")
	sanityCheck(f.BzNext, "BzNext")
	size := f.Width * f.Height * f.Depth
	minVal, maxVal := f.PotentialNext[0], f.PotentialNext[0]
	for i := 1; i < size; i++ {
		if f.PotentialNext[i] < minVal {
			minVal = f.PotentialNext[i]
		}
		if f.PotentialNext[i] > maxVal {
			maxVal = f.PotentialNext[i]
		}
	}
	if minVal < -1.0 || maxVal > 1.0 {
		valRange := maxVal - minVal
		if valRange < 1e-9 {
			valRange = 1e-9
		}
		for i := 0; i < size; i++ {
			f.PotentialNext[i] = -1.0 + 2.0*(f.PotentialNext[i]-minVal)/valRange
		}
	}
	var totalPotential, nextTotalPotential float64
	for i := 0; i < size; i++ {
		totalPotential += f.Potential[i]
		nextTotalPotential += f.PotentialNext[i]
	}
	if math.Abs(nextTotalPotential) > 1e-9 {
		correctionFactor := totalPotential / nextTotalPotential
		if correctionFactor < 0.5 || correctionFactor > 2 {
			fmt.Printf("CF: %0.16f\n", correctionFactor)
		} else {
			for i := 0; i < size; i++ {
				f.PotentialNext[i] *= correctionFactor
			}
		}
	}
	f.Potential, f.PotentialNext = f.PotentialNext, f.Potential
	f.FlowX, f.FlowXNext = f.FlowXNext, f.FlowX
	f.FlowY, f.FlowYNext = f.FlowYNext, f.FlowY
	f.FlowZ, f.FlowZNext = f.FlowZNext, f.FlowZ
	f.Ex, f.ExNext = f.ExNext, f.Ex
	f.Ey, f.EyNext = f.EyNext, f.Ey
	f.Ez, f.EzNext = f.EzNext, f.Ez
	f.Bx, f.BxNext = f.BxNext, f.Bx
	f.By, f.ByNext = f.ByNext, f.By
	f.Bz, f.BzNext = f.BzNext, f.Bz

	// Apply FFT smoothing to remove voxel artifacts
	f.Potential = f.fftSmooth(f.Potential)
	f.FlowX = f.fftSmooth(f.FlowX)
	f.FlowY = f.fftSmooth(f.FlowY)
	f.FlowZ = f.fftSmooth(f.FlowZ)
	f.Ex = f.fftSmooth(f.Ex)
	f.Ey = f.fftSmooth(f.Ey)
	f.Ez = f.fftSmooth(f.Ez)
	f.Bx = f.fftSmooth(f.Bx)
	f.By = f.fftSmooth(f.By)
	f.Bz = f.fftSmooth(f.Bz)
}

// AnalyzeField performs a full scan of the potential field for key metrics.
func (f *Field) AnalyzeField() FieldAnalysis {
	if len(f.Potential) == 0 {
		return FieldAnalysis{}
	}

	analysis := FieldAnalysis{
		Min: f.Potential[0],
		Max: f.Potential[0],
	}
	maxIndex1D := 0

	for i, v := range f.Potential {
		analysis.Total += v
		if v < analysis.Min {
			analysis.Min = v
		}
		if v > analysis.Max {
			analysis.Max = v
			maxIndex1D = i
		}
	}

	// Convert 1D index back to 3D coordinates
	w, h := f.Width, f.Height
	analysis.MaxZ = maxIndex1D / (w * h)
	analysis.MaxY = (maxIndex1D / w) % h
	analysis.MaxX = maxIndex1D % w

	return analysis
}

// --- EBITENGINE APPLICATION ---

const (
	div             = 4
	simWidth        = 640 / div
	simHeight       = 640 / div
	simDepth        = 640 / 10 / div
	infoPanelHeight = 80
	screenWidth     = simWidth * div
	screenHeight    = simHeight*div + infoPanelHeight
	axisX           = 0
	axisY           = 1
	axisZ           = 2
)

type Game struct {
	simField   *Field
	sliceImage *ebiten.Image
	sliceAxis  int
	sliceIndex int
	frameCount int64
	analysis   FieldAnalysis
}

func NewGame() *Game {
	return &Game{
		simField:   NewField(simWidth, simHeight, simDepth),
		sliceAxis:  axisZ,
		sliceIndex: simDepth / 2,
	}
}

func (g *Game) Update() error {
	// First, run analysis on the current state for use in input and drawing
	g.analysis = g.simField.AnalyzeField()

	// Handle Input
	if inpututil.IsKeyJustPressed(ebiten.KeyZ) {
		g.sliceAxis, g.sliceIndex = axisZ, g.simField.Depth/2
	}
	if inpututil.IsKeyJustPressed(ebiten.KeyY) {
		g.sliceAxis, g.sliceIndex = axisY, g.simField.Height/2
	}
	if inpututil.IsKeyJustPressed(ebiten.KeyX) {
		g.sliceAxis, g.sliceIndex = axisX, g.simField.Width/2
	}
	if inpututil.IsKeyJustPressed(ebiten.KeyQ) {
		os.Exit(0)
	}

	// 'H' key to go to the hottest slice
	if ebiten.IsKeyPressed(ebiten.KeyH) {
		switch g.sliceAxis {
		case axisX:
			g.sliceIndex = g.analysis.MaxX
		case axisY:
			g.sliceIndex = g.analysis.MaxY
		case axisZ:
			g.sliceIndex = g.analysis.MaxZ
		}
	}

	// Arrow keys to move slice
	if ebiten.IsKeyPressed(ebiten.KeyShift) {
		if ebiten.IsKeyPressed(ebiten.KeyArrowUp) {
			g.sliceIndex++
		}
		if ebiten.IsKeyPressed(ebiten.KeyArrowDown) {
			g.sliceIndex--
		}
	} else {
		if inpututil.IsKeyJustPressed(ebiten.KeyArrowUp) {
			g.sliceIndex++
		}
		if inpututil.IsKeyJustPressed(ebiten.KeyArrowDown) {
			g.sliceIndex--
		}
	}

	if inpututil.IsKeyJustPressed(ebiten.KeyR) {
		g.simField.randomize()
	}

	// Clamp sliceIndex to be within bounds
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

	// Run simulation step and increment frame counter
	g.simField.Step()
	g.frameCount++
	return nil
}

func (g *Game) Draw(screen *ebiten.Image) {
	// 1. Draw the info panel at the top
	axisStr := "Z"
	if g.sliceAxis == axisX {
		axisStr = "X"
	} else if g.sliceAxis == axisY {
		axisStr = "Y"
	}

	msg := fmt.Sprintf(
		"FPS: %.2f | Frame: %d\n"+
			"Slicing Axis: %s | Index: %d | Press H to jump to Max\n"+
			"Total Potential: %e | Min: %.4f | Max: %.4f\n"+
			"Hottest Voxel at (X: %d, Y: %d, Z: %d)\n"+
			"CONTROLS: Q to quit | X/Y/Z keys to change axis | Up/Down (Hold Shift for fast) to change slice",
		ebiten.ActualFPS(), g.frameCount,
		axisStr, g.sliceIndex,
		g.analysis.Total, g.analysis.Min, g.analysis.Max,
		g.analysis.MaxX, g.analysis.MaxY, g.analysis.MaxZ,
	)
	ebitenutil.DebugPrint(screen, msg)

	// 2. Prepare and draw the simulation slice below the info panel
	viewWidth, viewHeight := 0, 0
	switch g.sliceAxis {
	case axisZ:
		viewWidth, viewHeight = g.simField.Width, g.simField.Height
	case axisY:
		viewWidth, viewHeight = g.simField.Width, g.simField.Depth
	case axisX:
		viewWidth, viewHeight = g.simField.Height, g.simField.Depth
	}

	if g.sliceImage == nil || g.sliceImage.Bounds().Dx() != viewWidth || g.sliceImage.Bounds().Dy() != viewHeight {
		g.sliceImage = ebiten.NewImage(viewWidth, viewHeight)
	}
	pixels := make([]byte, viewWidth*viewHeight*4)

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

			val := (potential + 1.0) * 0.5
			idx := (j*viewWidth + i) * 4
			if val < 0.5 {
				pixels[idx+2] = uint8(255 * (1.0 - val*2))
			} else {
				pixels[idx+0] = uint8(255 * (val - 0.5) * 2)
			}
			pixels[idx+3] = 255
		}
	}
	g.sliceImage.WritePixels(pixels)

	op := &ebiten.DrawImageOptions{}
	// Calculate scaling to fit in the area below the info panel
	scaleW := float64(screenWidth) / float64(viewWidth)
	scaleH := float64(screenHeight-infoPanelHeight) / float64(viewHeight)
	scale := math.Min(scaleW, scaleH)
	op.GeoM.Scale(scale, scale)
	// Translate the image down to be below the info panel
	op.GeoM.Translate(0, infoPanelHeight)
	screen.DrawImage(g.sliceImage, op)
}

func (g *Game) Layout(outsideWidth, outsideHeight int) (int, int) {
	return screenWidth, screenHeight
}

func main() {
	game := NewGame()
	ebiten.SetWindowSize(screenWidth, screenHeight)
	ebiten.SetWindowTitle("Native Field F Simulator v2.2 (UI Enhanced)")
	if err := ebiten.RunGame(game); err != nil {
		log.Fatal(err)
	}
}
