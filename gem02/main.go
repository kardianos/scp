package main

import (
	"fmt"
	"log"
	"math"
	"math/cmplx"
	"os"
	"runtime"
	"strconv"
	"sync"

	"github.com/jvlmdr/go-fftw/fftw"
)

// All parameters are unchanged
const (
	N               = 64
	L               = 20.0
	C               = 1.0
	G               = -10.0
	DT              = 0.005
	BeamWaist       = 1.5
	BeamSeparation  = L / 3.0
	BeamChargeL     = 2
	InitialVelocity = 15.0
	TotalTime       = 2.0
	OutputDir       = "output"
	VisEveryNSteps  = 5
)

// idx, createKVector are unchanged
func idx(i, j, k int) int { return i*N*N + j*N + k }
func createKVector(n int, l float64) []float64 { /* ... */
	k := make([]float64, n)
	factor := 2.0 * math.Pi / l
	for i := 0; i < n/2; i++ {
		k[i] = float64(i) * factor
	}
	for i := n / 2; i < n; i++ {
		k[i] = float64(i-n) * factor
	}
	return k
}

// CORRECTED createLaguerreGaussian function
func createLaguerreGaussian(x, y, z, l, p, w0, x0, y0, z0 float64) complex128 {
	xPrime, yPrime, zPrime := x-x0, y-y0, z-z0

	// THE FIX: Add a tiny epsilon to r to avoid Log(0)
	r := math.Sqrt(xPrime*xPrime+yPrime*yPrime) + 1e-15

	phi := math.Atan2(yPrime, xPrime)
	lgProfile := cmplx.Pow(complex(r/w0, 0), complex(math.Abs(l), 0)) *
		cmplx.Exp(complex(-r*r/(w0*w0), l*phi))
	zEnvelope := math.Exp(-zPrime * zPrime / (2 * w0 * w0))
	return lgProfile * complex(zEnvelope, 0)
}

// All other helper functions (applyKick, drift, scaleSlice) are unchanged
func applyKick(psi, v, laplacianPsi []complex128, wg *sync.WaitGroup, numWorkers int) []complex128 { /* ... */
	vNew := make([]complex128, len(v))
	chunkSize := len(psi) / numWorkers
	for i := 0; i < numWorkers; i++ {
		start := i * chunkSize
		end := start + chunkSize
		if i == numWorkers-1 {
			end = len(psi)
		}
		wg.Add(1)
		go func(start, end int) {
			defer wg.Done()
			dtHalf := complex(0.5*DT, 0)
			gComplex := complex(G, 0)
			c2Complex := complex(C*C, 0)
			for j := start; j < end; j++ {
				absPsi2 := cmplx.Abs(psi[j]) * cmplx.Abs(psi[j])
				nonlinearTerm := gComplex * complex(absPsi2, 0) * psi[j]
				acceleration := c2Complex*laplacianPsi[j] - nonlinearTerm
				vNew[j] = v[j] + acceleration*dtHalf
			}
		}(start, end)
	}
	wg.Wait()
	return vNew
}
func drift(psi, v []complex128, wg *sync.WaitGroup, numWorkers int) []complex128 { /* ... */
	psiNew := make([]complex128, len(psi))
	chunkSize := len(psi) / numWorkers
	for i := 0; i < numWorkers; i++ {
		start := i * chunkSize
		end := start + chunkSize
		if i == numWorkers-1 {
			end = len(psi)
		}
		wg.Add(1)
		go func(start, end int) {
			defer wg.Done()
			dtComplex := complex(DT, 0)
			for j := start; j < end; j++ {
				psiNew[j] = psi[j] + v[j]*dtComplex
			}
		}(start, end)
	}
	wg.Wait()
	return psiNew
}
func scaleSlice(slice []complex128, scaleFactor float64, wg *sync.WaitGroup, numWorkers int) { /* ... */
	chunkSize := len(slice) / numWorkers
	scaleComplex := complex(scaleFactor, 0)
	for i := 0; i < numWorkers; i++ {
		start := i * chunkSize
		end := start + chunkSize
		if i == numWorkers-1 {
			end = len(slice)
		}
		wg.Add(1)
		go func(start, end int) {
			defer wg.Done()
			for j := start; j < end; j++ {
				slice[j] *= scaleComplex
			}
		}(start, end)
	}
	wg.Wait()
}

func main() {
	// This is the version from your last message, which uses the correct
	// multi-dimensional plan functions.
	numWorkers := runtime.NumCPU()
	runtime.GOMAXPROCS(numWorkers)
	log.Printf("Starting simulation on %d CPU cores", numWorkers)
	dims := []int{N, N, N}
	totalPoints := N * N * N

	in := fftw.NewArrayN(dims)
	out := fftw.NewArrayN(dims)
	forwardPlan := fftw.NewPlanN(in, out, fftw.Forward, fftw.Estimate)
	backwardPlan := fftw.NewPlanN(out, in, fftw.Backward, fftw.Estimate)
	defer forwardPlan.Destroy()
	defer backwardPlan.Destroy()

	// The rest of the main function is correct and unchanged
	psi := make([]complex128, totalPoints)
	v := make([]complex128, totalPoints)
	kSquaredOp := make([]complex128, totalPoints)

	gridPts := make([]float64, N)
	for i := 0; i < N; i++ {
		gridPts[i] = -L/2 + float64(i)*(L/float64(N-1))
	}
	kxVec, kyVec, kzVec := createKVector(N, L), createKVector(N, L), createKVector(N, L)
	psi1 := make([]complex128, totalPoints)
	psi2 := make([]complex128, totalPoints)
	for i := 0; i < N; i++ {
		for j := 0; j < N; j++ {
			for k := 0; k < N; k++ {
				k2 := kxVec[i]*kxVec[i] + kyVec[j]*kyVec[j] + kzVec[k]*kzVec[k]
				kSquaredOp[idx(i, j, k)] = complex(-k2, 0)
				x, y, z := gridPts[i], gridPts[j], gridPts[k]
				psi1[idx(i, j, k)] = createLaguerreGaussian(x, y, z, BeamChargeL, 0, BeamWaist, -BeamSeparation/2, 0, 0)
				psi2[idx(i, j, k)] = createLaguerreGaussian(x, y, z, -BeamChargeL, 0, BeamWaist, BeamSeparation/2, 0, 0)
				psi[idx(i, j, k)] = psi1[idx(i, j, k)] + psi2[idx(i, j, k)]
			}
		}
	}

	copy(in.Elems, psi1)
	forwardPlan.Execute()
	psi1F := append(out.Elems[:0:0], out.Elems...)
	copy(in.Elems, psi2)
	forwardPlan.Execute()
	psi2F := append(out.Elems[:0:0], out.Elems...)

	vF := make([]complex128, totalPoints)
	for i := 0; i < N; i++ {
		for j := 0; j < N; j++ {
			for k := 0; k < N; k++ {
				kIndex := idx(i, j, k)
				kxComplex := complex(0, kxVec[i]*InitialVelocity)
				vF[kIndex] = (kxComplex * psi1F[kIndex]) + (-kxComplex * psi2F[kIndex])
			}
		}
	}
	copy(out.Elems, vF)
	backwardPlan.Execute()
	copy(v, in.Elems)

	scaleFactor := 1.0 / float64(totalPoints)
	wg := &sync.WaitGroup{}
	scaleSlice(v, scaleFactor, wg, numWorkers)

	if err := os.MkdirAll(OutputDir, 0755); err != nil {
		log.Fatalf("Failed to create output directory: %v", err)
	}
	log.Println("Initial conditions set. Writing initial state.")
	writeDataToFile(0, psi)

	numSteps := int(TotalTime / DT)
	log.Println("Starting simulation loop...")
	laplacianPsi := make([]complex128, totalPoints)
	for i := 1; i <= numSteps; i++ {
		copy(in.Elems, psi)
		forwardPlan.Execute()
		outData := out.Elems
		for j := range outData {
			outData[j] *= kSquaredOp[j]
		}
		backwardPlan.Execute()
		copy(laplacianPsi, in.Elems)
		scaleSlice(laplacianPsi, scaleFactor, wg, numWorkers)
		v = applyKick(psi, v, laplacianPsi, wg, numWorkers)
		psi = drift(psi, v, wg, numWorkers)
		copy(in.Elems, psi)
		forwardPlan.Execute()
		outData = out.Elems
		for j := range outData {
			outData[j] *= kSquaredOp[j]
		}
		backwardPlan.Execute()
		copy(laplacianPsi, in.Elems)
		scaleSlice(laplacianPsi, scaleFactor, wg, numWorkers)
		v = applyKick(psi, v, laplacianPsi, wg, numWorkers)

		if i%VisEveryNSteps == 0 {
			log.Printf("Step %d/%d (t=%.2f)", i, numSteps, float64(i)*DT)
			writeDataToFile(i, psi)
		}
	}
	log.Println("Simulation finished.")
}

// writeDataToFile and its helpers are unchanged and correct
func writeDataToFile(step int, psi []complex128) { /* ... */
	filePath := fmt.Sprintf("%s/output_%04d.vts", OutputDir, step)
	f, err := os.Create(filePath)
	if err != nil {
		log.Printf("Error creating file: %v", err)
		return
	}
	defer f.Close()
	density := make([]float64, len(psi))
	phase := make([]float64, len(psi))
	maxDensity := 0.0
	for i, val := range psi {
		d := cmplx.Abs(val) * cmplx.Abs(val)
		if d > maxDensity {
			maxDensity = d
		}
		density[i] = d
		phase[i] = cmplx.Phase(val)
	}
	log.Printf("Writing step %d, max density: %.4e", step, maxDensity)
	extentStr := "0 " + strconv.Itoa(N-1) + " 0 " + strconv.Itoa(N-1) + " 0 " + strconv.Itoa(N-1)
	f.WriteString("<VTKFile type=\"StructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n")
	f.WriteString("  <StructuredGrid WholeExtent=\"" + extentStr + "\">\n")
	f.WriteString("    <Piece Extent=\"" + extentStr + "\">\n")
	f.WriteString("      <Points>\n")
	coords := make([]float32, len(psi)*3)
	for i := 0; i < N; i++ {
		for j := 0; j < N; j++ {
			for k := 0; k < N; k++ {
				idx3 := idx(i, j, k) * 3
				coords[idx3+0] = float32(-L/2 + float64(i)*(L/float64(N-1)))
				coords[idx3+1] = float32(-L/2 + float64(j)*(L/float64(N-1)))
				coords[idx3+2] = float32(-L/2 + float64(k)*(L/float64(N-1)))
			}
		}
	}
	writeDataArray(f, "Coordinates", 3, coords)
	f.WriteString("      </Points>\n")
	f.WriteString("      <PointData Scalars=\"density\">\n")
	writeScalarData(f, "density", density)
	writeScalarData(f, "phase", phase)
	f.WriteString("      </PointData>\n")
	f.WriteString("    </Piece>\n")
	f.WriteString("  </StructuredGrid>\n")
	f.WriteString("</VTKFile>\n")
}
func writeScalarData(f *os.File, name string, data []float64) { /* ... */
	f.WriteString("        <DataArray type=\"Float64\" Name=\"" + name + "\" format=\"ascii\">\n")
	for _, v := range data {
		f.WriteString(fmt.Sprintf("          %.6e\n", v))
	}
	f.WriteString("        </DataArray>\n")
}
func writeDataArray(f *os.File, name string, numComponents int, data []float32) { /* ... */
	f.WriteString("        <DataArray type=\"Float32\" Name=\"" + name + "\" NumberOfComponents=\"" + strconv.Itoa(numComponents) + "\" format=\"ascii\">\n")
	for i := 0; i < len(data); i += 3 {
		f.WriteString(fmt.Sprintf("          %.6e %.6e %.6e\n", data[i], data[i+1], data[i+2]))
	}
	f.WriteString("        </DataArray>\n")
}
