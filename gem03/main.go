package main

import (
	"embed"
	"fmt"
	"image/color"
	"log"
	"math"
	"math/cmplx"
	"runtime"
	"sync"

	"github.com/hajimehoshi/ebiten/v2"
	"github.com/hajimehoshi/ebiten/v2/inpututil"
	"github.com/hajimehoshi/ebiten/v2/text"
	"github.com/jvlmdr/go-fftw/fftw"
	"golang.org/x/image/font"
	"golang.org/x/image/font/opentype"
)

//go:embed resources/GoRegular.ttf
var goFont embed.FS

// --- Parameters ---
const (
	N               = 128
	L               = 24.0
	C               = 1.0
	G               = -10.0
	DT              = 0.002
	BeamWaist       = 2.0
	BeamSeparation  = L / 3.0
	BeamChargeL     = 2
	InitialVelocity = 15.0
	FastSliceStep   = 10 // How many slices to jump when holding Shift
)

type Game struct {
	forwardPlan, backwardPlan *fftw.Plan
	in, out                   *fftw.ArrayN
	psi, v, laplacianPsi      []complex128
	kSquaredOp                []complex128
	wg                        *sync.WaitGroup
	pixels                    *ebiten.Image
	infoFont                  font.Face
	sliceZ                    int // The currently viewed Z-slice
}

func NewGame() (*Game, error) {
	log.Println("Setting up simulation...")
	numWorkers := runtime.NumCPU()
	runtime.GOMAXPROCS(numWorkers)
	totalPoints := N * N * N
	dims := []int{N, N, N}

	// Load font for UI
	tt, err := opentype.Parse(must(goFont.ReadFile("resources/GoRegular.ttf")))
	if err != nil {
		return nil, err
	}
	infoFont, err := opentype.NewFace(tt, &opentype.FaceOptions{
		Size:    16,
		DPI:     72,
		Hinting: font.HintingFull,
	})
	if err != nil {
		return nil, err
	}

	g := &Game{
		in:           fftw.NewArrayN(dims),
		out:          fftw.NewArrayN(dims),
		psi:          make([]complex128, totalPoints),
		v:            make([]complex128, totalPoints),
		laplacianPsi: make([]complex128, totalPoints),
		kSquaredOp:   make([]complex128, totalPoints),
		pixels:       ebiten.NewImage(N, N),
		wg:           &sync.WaitGroup{},
		infoFont:     infoFont,
		sliceZ:       N / 2, // Start in the middle
	}
	// Use 1D FFT plan as requested
	g.forwardPlan = fftw.NewPlanN(g.in, g.out, fftw.Forward, fftw.Estimate)
	g.backwardPlan = fftw.NewPlanN(g.out, g.in, fftw.Backward, fftw.Estimate)

	// --- Initial Conditions ---
	gridPts := make([]float64, N)
	for i := 0; i < N; i++ {
		gridPts[i] = -L/2 + float64(i)*(L/float64(N-1))
	}
	// We use a 1D k-vector for our 1D FFT
	kVec := createKVector(totalPoints, L*float64(N))
	for i := 0; i < totalPoints; i++ {
		g.kSquaredOp[i] = complex(-(kVec[i] * kVec[i]), 0)
	}

	psi1 := make([]complex128, totalPoints)
	psi2 := make([]complex128, totalPoints)
	for i := 0; i < N; i++ {
		for j := 0; j < N; j++ {
			for k := 0; k < N; k++ {
				x, y, z := gridPts[i], gridPts[j], gridPts[k]
				idx := idx(i, j, k)
				psi1[idx] = createLaguerreGaussian(x, y, z, BeamChargeL, 0, BeamWaist, -BeamSeparation/2, 0, 0)
				psi2[idx] = createLaguerreGaussian(x, y, z, -BeamChargeL, 0, BeamWaist, BeamSeparation/2, 0, 0)
				g.psi[idx] = psi1[idx] + psi2[idx]
			}
		}
	}

	copy(g.in.Elems, psi1)
	g.forwardPlan.Execute()
	psi1F := append(g.out.Elems[:0:0], g.out.Elems...)
	copy(g.in.Elems, psi2)
	g.forwardPlan.Execute()
	psi2F := append(g.out.Elems[:0:0], g.out.Elems...)

	vF := make([]complex128, totalPoints)
	// Apply initial velocity using a 3D k-vector for physical accuracy
	kxVec3D, _, _ := createKVector3D(N, L)
	for i := 0; i < totalPoints; i++ {
		kx := kxVec3D[i]
		kxComplex := complex(0, kx*InitialVelocity)
		vF[i] = (kxComplex * psi1F[i]) + (-kxComplex * psi2F[i])
	}
	copy(g.out.Elems, vF)
	g.backwardPlan.Execute()
	copy(g.v, g.in.Elems)

	scaleFactor := 1.0 / float64(totalPoints)
	scaleSlice(g.v, scaleFactor, g.wg, numWorkers)

	log.Println("Setup complete. Starting simulation loop.")
	return g, nil
}

func (g *Game) Update() error {
	// --- Handle Input for Slice Control ---
	step := 1
	if ebiten.IsKeyPressed(ebiten.KeyShift) {
		step = FastSliceStep
	}
	if inpututil.IsKeyJustPressed(ebiten.KeyUp) {
		g.sliceZ += step
	}
	if inpututil.IsKeyJustPressed(ebiten.KeyDown) {
		g.sliceZ -= step
	}
	// Clamp the slice value to be within bounds
	if g.sliceZ >= N {
		g.sliceZ = N - 1
	}
	if g.sliceZ < 0 {
		g.sliceZ = 0
	}

	// --- Physics Update ---
	numWorkers := runtime.NumCPU()
	scaleFactor := 1.0 / float64(N*N*N)
	for i := 0; i < 2; i++ {
		copy(g.in.Elems, g.psi)
		g.forwardPlan.Execute()
		outData := g.out.Elems
		for j := range outData {
			outData[j] *= g.kSquaredOp[j]
		}
		g.backwardPlan.Execute()
		copy(g.laplacianPsi, g.in.Elems)
		scaleSlice(g.laplacianPsi, scaleFactor, g.wg, numWorkers)
		g.v = applyKick(g.psi, g.v, g.laplacianPsi, g.wg, numWorkers)
		g.psi = drift(g.psi, g.v, g.wg, numWorkers)
		copy(g.in.Elems, g.psi)
		g.forwardPlan.Execute()
		outData = g.out.Elems
		for j := range outData {
			outData[j] *= g.kSquaredOp[j]
		}
		g.backwardPlan.Execute()
		copy(g.laplacianPsi, g.in.Elems)
		scaleSlice(g.laplacianPsi, scaleFactor, g.wg, numWorkers)
		g.v = applyKick(g.psi, g.v, g.laplacianPsi, g.wg, numWorkers)
	}
	return nil
}

func (g *Game) Draw(screen *ebiten.Image) {
	// --- Draw the 2D slice ---
	maxDensity := 1e-9
	for i := 0; i < N; i++ {
		for j := 0; j < N; j++ {
			density := cmplx.Abs(g.psi[idx(i, j, g.sliceZ)])
			if density > maxDensity {
				maxDensity = density
			}
		}
	}
	for i := 0; i < N; i++ {
		for j := 0; j < N; j++ {
			val := g.psi[idx(i, j, g.sliceZ)]
			density := cmplx.Abs(val)
			phase := cmplx.Phase(val)
			v := math.Sqrt(density / maxDensity)
			if v > 1.0 {
				v = 1.0
			}
			h := (phase + math.Pi) * (180 / math.Pi)
			g.pixels.Set(i, j, hsvToRgb(h, 1.0, v))
		}
	}
	op := &ebiten.DrawImageOptions{}
	sX := float64(screen.Bounds().Dx()) / float64(N)
	sY := float64(screen.Bounds().Dy()) / float64(N)
	op.GeoM.Scale(sX, sY)
	screen.DrawImage(g.pixels, op)

	// --- Draw the Informational Text ---
	infoStr := fmt.Sprintf("Slice Z = %d / %d", g.sliceZ, N-1)
	text.Draw(screen, infoStr, g.infoFont, 10, screen.Bounds().Dy()-10, color.White)
}

func (g *Game) Layout(outsideWidth, outsideHeight int) (screenWidth, screenHeight int) {
	return N, N // Our logical grid is N x N
}

func main() {
	game, err := NewGame()
	if err != nil {
		log.Fatal(err)
	}
	ebiten.SetWindowSize(N*5, N*5)
	ebiten.SetWindowTitle("3D Harmonics Simulation (Use Up/Down Arrows to Change Slice)")
	if err := ebiten.RunGame(game); err != nil {
		log.Fatal(err)
	}
}

// --- Helper Functions ---
func must[T any](v T, err error) T {
	if err != nil {
		log.Fatal(err)
	}
	return v
}
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
func createKVector3D(n int, l float64) ([]float64, []float64, []float64) {
	kx, ky, kz := make([]float64, n*n*n), make([]float64, n*n*n), make([]float64, n*n*n)
	k_base := createKVector(n, l)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			for k := 0; k < n; k++ {
				idx_ := idx(i, j, k)
				kx[idx_], ky[idx_], kz[idx_] = k_base[i], k_base[j], k_base[k]
			}
		}
	}
	return kx, ky, kz
}
func createLaguerreGaussian(x, y, z, l, p, w0, x0, y0, z0 float64) complex128 {
	xPrime, yPrime, zPrime := x-x0, y-y0, z-z0
	r := math.Sqrt(xPrime*xPrime+yPrime*yPrime) + 1e-15
	phi := math.Atan2(yPrime, xPrime)
	lgProfile := cmplx.Pow(complex(r/w0, 0), complex(math.Abs(l), 0)) *
		cmplx.Exp(complex(-r*r/(w0*w0), l*phi))
	zEnvelope := math.Exp(-zPrime * zPrime / (2 * w0 * w0))
	return lgProfile * complex(zEnvelope, 0)
}
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
			dtHalf, gComplex, c2Complex := complex(0.5*DT, 0), complex(G, 0), complex(C*C, 0)
			for j := start; j < end; j++ {
				absPsi2 := cmplx.Abs(psi[j]) * cmplx.Abs(psi[j])
				vNew[j] = v[j] + (c2Complex*laplacianPsi[j]-gComplex*complex(absPsi2, 0)*psi[j])*dtHalf
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
func hsvToRgb(h, s, v float64) color.RGBA {
	c, m := v*s, v-v*s
	x := c * (1 - math.Abs(math.Mod(h/60, 2)-1))
	var r, g, b float64
	switch {
	case h < 60:
		r, g, b = c, x, 0
	case h < 120:
		r, g, b = x, c, 0
	case h < 180:
		r, g, b = 0, c, x
	case h < 240:
		r, g, b = 0, x, c
	case h < 300:
		r, g, b = x, 0, c
	default:
		r, g, b = c, 0, x
	}
	return color.RGBA{uint8((r + m) * 255), uint8((g + m) * 255), uint8((b + m) * 255), 255}
}
