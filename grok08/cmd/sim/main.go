package main

import (
	"bufio"
	"encoding/binary"
	"flag"
	"fmt"
	"log"
	"math"
	"math/rand"
	"os"

	sfft "github.com/davidkleiven/gosfft/sfft"
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
	Width, Height, Depth                                                                  int
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

	w, h := f.Width, f.Height
	analysis.MaxZ = maxIndex1D / (w * h)
	analysis.MaxY = (maxIndex1D / w) % h
	analysis.MaxX = maxIndex1D % w

	return analysis
}

const (
	simWidth  = 64
	simHeight = 64
	simDepth  = 64
)

func main() {
	steps := flag.Int("steps", 200, "number of simulation steps")
	flag.Parse()

	field := NewField(simWidth, simHeight, simDepth)

	file, err := os.Create("field.scpv")
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	buf := bufio.NewWriter(file)
	defer buf.Flush()

	// Write magic
	magic := []byte("SCPV\x00\x00\x00\x00")
	buf.Write(magic)

	// Write dimensions
	dim := func(tp int64, label string, value int32) {
		binary.Write(buf, binary.BigEndian, tp)
		binary.Write(buf, binary.BigEndian, int16(len(label)))
		buf.WriteString(label)
		binary.Write(buf, binary.BigEndian, int64(4))
		binary.Write(buf, binary.BigEndian, int32(value))
	}
	txt := func(tp int64, label string, value string) {
		binary.Write(buf, binary.BigEndian, tp)
		binary.Write(buf, binary.BigEndian, int16(len(label)))
		buf.WriteString(label)
		binary.Write(buf, binary.BigEndian, int64(len(value)))
		buf.WriteString(value)
	}
	dim(0, "height", int32(field.Height))
	dim(1, "width", int32(field.Width))
	dim(2, "length", int32(field.Depth))
	txt(3, "encoding", "float64")

	/*
		binary.Write(file, binary.BigEndian, int64(0)) // type height
		binary.Write(file, binary.BigEndian, int16(len("height")))
		file.Write([]byte("height"))
		binary.Write(file, binary.BigEndian, int64(4))
		binary.Write(file, binary.BigEndian, int32(field.Height))

		binary.Write(file, binary.BigEndian, int64(1)) // type width
		binary.Write(file, binary.BigEndian, int16(len("width")))
		file.Write([]byte("width"))
		binary.Write(file, binary.BigEndian, int64(4))
		binary.Write(file, binary.BigEndian, int32(field.Width))

		binary.Write(file, binary.BigEndian, int64(2)) // type length
		binary.Write(file, binary.BigEndian, int16(len("length")))
		file.Write([]byte("length"))
		binary.Write(file, binary.BigEndian, int64(4))
		binary.Write(file, binary.BigEndian, int32(field.Depth))

		binary.Write(file, binary.BigEndian, int64(3)) // type encoding
		binary.Write(file, binary.BigEndian, int16(len("encoding")))
		file.Write([]byte("encoding"))
		binary.Write(file, binary.BigEndian, int64(7))
		file.Write([]byte("float64"))
	*/

	data := make([]byte, len(field.Potential)*8)
	for step := 0; step < *steps; step++ {
		field.Step()

		// Write potential as example field data
		name := fmt.Sprintf("t%d", step)
		// data := make([]byte, len(field.Potential)*8)
		for i, v := range field.Potential {
			binary.LittleEndian.PutUint64(data[i*8:(i+1)*8], math.Float64bits(v))
		}

		binary.Write(buf, binary.BigEndian, int64(4)) // type
		binary.Write(buf, binary.BigEndian, int16(len(name)))
		buf.WriteString(name)
		binary.Write(buf, binary.BigEndian, int64(len(data)))
		buf.Write(data)
	}

	fmt.Println("Simulation complete, output to field.scpv")
}
