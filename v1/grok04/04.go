package main

import (
	"context"
	"fmt"
	"math"
	"sync"
	"time"

	sfft "github.com/davidkleiven/gosfft/sfft"
	rl "github.com/gen2brain/raylib-go/raylib"
	"gonum.org/v1/gonum/cmplxs"
)

const (
	N         = 32   // Grid size
	L         = 10.0 // Box length
	dx        = L / N
	dt        = 0.01 // Time step
	kParam    = 1.0  // Linear coefficient
	lam       = 0.1  // Nonlinear coefficient
	A         = 1.0  // Initial amplitude
	threshold = 0.2  // Voxel threshold
	voxelSize = float32(L / N)
)

type Field []complex128 // Flattened 3D field (N*N*N elements)

var (
	frames    []Field
	mu        sync.Mutex
	autoAdv   = true
	reverse   = false
	currFrame = 0
	emitted   = false
	t         = 0.0
	fft3      = sfft.NewFFT3(N, N, N)
)

// Init standing wave
func initPhi() Field {
	phi := make(Field, N*N*N)
	for i := 0; i < N; i++ {
		for j := 0; j < N; j++ {
			for k := 0; k < N; k++ {
				x, y, z := float64(i)*dx, float64(j)*dx, float64(k)*dx
				phi[idx(i, j, k)] = complex(A*math.Sin(2*math.Pi*x/L)*math.Sin(2*math.Pi*y/L)*math.Sin(2*math.Pi*z/L), 0)
			}
		}
	}
	return phi
}

func idx(i, j, k int) int { return i*N*N + j*N + k }

// Compute Laplacian via FFT
func laplacian(phi Field) Field {
	// Forward FFT
	coeff := fft3.FFT(phi)

	// Apply -k^2
	for i := 0; i < N; i++ {
		kx := 2 * math.Pi * (float64(i) - float64(N)/2) / L
		for j := 0; j < N; j++ {
			ky := 2 * math.Pi * (float64(j) - float64(N)/2) / L
			for k := 0; k < N; k++ {
				kz := 2 * math.Pi * (float64(k) - float64(N)/2) / L
				coeff[idx(i, j, k)] *= complex(-(kx*kx + ky*ky + kz*kz), 0)
			}
		}
	}

	// Inverse FFT
	result := fft3.IFFT(coeff)
	cmplxs.Scale(complex(1.0/float64(N*N*N), 0), result) // Normalize
	return result
}

// Time step (leapfrog)
func step(phi, vel *Field) {
	accel := laplacian(*phi)
	for i := range accel {
		realPhi := real((*phi)[i])
		accel[i] += complex(-kParam*realPhi-lam*math.Pow(realPhi, 3), 0)
	}

	// Update vel and phi (simple leapfrog)
	newVel := make(Field, len(*vel))
	copy(newVel, *vel)
	cmplxs.AddScaled(newVel, complex(0.5*dt, 0), accel)

	newPhi := make(Field, len(*phi))
	copy(newPhi, *phi)
	cmplxs.AddScaled(newPhi, complex(dt, 0), newVel)

	newAccel := laplacian(newPhi)
	for i := range newAccel {
		realPhi := real(newPhi[i])
		newAccel[i] += complex(-kParam*realPhi-lam*math.Pow(realPhi, 3), 0)
	}
	cmplxs.AddScaled(newVel, complex(0.5*dt, 0), newAccel)

	*phi = newPhi
	*vel = newVel
}

// Simulation goroutine: Evolve and record
func simulate(ctx context.Context) {
	phi := initPhi()
	vel := make(Field, N*N*N) // Zeros
	stepCount := 0

	for {
		select {
		case <-ctx.Done():
			return
		default:
			if !emitted && t >= 1.0 {
				cmplxs.Scale(complex(0.5, 0), phi) // Lower harmonic
				// Add wave
				for i := 0; i < N; i++ {
					for j := 0; j < N; j++ {
						for k := 0; k < N; k++ {
							x, y, z := float64(i)*dx, float64(j)*dx, float64(k)*dx
							phi[idx(i, j, k)] += complex(0.5*math.Sin(2*math.Pi*(x/L-(t-1.0)))*math.Exp(-((y-L/2)*(y-L/2)+(z-L/2)*(z-L/2))/(L/4*L/4)), 0)
						}
					}
				}
				emitted = true
			}

			step(&phi, &vel)
			t += dt

			mu.Lock()
			frames = append(frames, append(Field(nil), phi...)) // Copy frame
			mu.Unlock()

			time.Sleep(time.Millisecond * 10) // Throttle for real-time

			stepCount++
			if stepCount%10 != 0 {
				continue // Check ctx every 10 steps for performance
			}
		}
	}
}

// Navigation goroutine: Poll for changes (simplified; extend with channels if needed)
func navigate(ctx context.Context) {
	for {
		select {
		case <-ctx.Done():
			return
		default:
			time.Sleep(time.Millisecond * 100)
		}
	}
}

func main() {
	ctx, cancel := context.WithCancel(context.Background())
	defer cancel()

	go simulate(ctx)
	go navigate(ctx)

	rl.InitWindow(800, 800, "Volumetric Sim")
	rl.SetTargetFPS(60)

	camera := rl.NewCamera3D(rl.NewVector3(15, 15, 15), rl.NewVector3(5, 5, 5), rl.NewVector3(0, 1, 0), 45, rl.CameraPerspective)

	for !rl.WindowShouldClose() {
		rl.UpdateCamera(&camera, rl.CameraOrbital)

		if rl.IsKeyPressed(rl.KeyR) {
			reverse = !reverse
		}
		if rl.IsKeyPressed(rl.KeyA) {
			autoAdv = !autoAdv
		}

		mu.Lock()
		frameLen := len(frames)
		displayT := t
		mu.Unlock()

		if frameLen == 0 {
			rl.BeginDrawing()
			rl.ClearBackground(rl.Black)
			rl.DrawText("Simulating...", 10, 10, 20, rl.White)
			rl.EndDrawing()
			continue
		}

		mu.Lock()
		if autoAdv {
			currFrame = frameLen - 1
		} else if reverse && currFrame > 0 {
			currFrame--
		}
		if currFrame >= frameLen {
			currFrame = frameLen - 1
		}
		if currFrame < 0 {
			currFrame = 0
		}
		frame := frames[currFrame]
		mu.Unlock()

		rl.BeginDrawing()
		rl.ClearBackground(rl.Black)

		rl.BeginMode3D(camera)
		// Draw voxels
		for i := 0; i < N; i++ {
			for j := 0; j < N; j++ {
				for k := 0; k < N; k++ {
					val := real(frame[idx(i, j, k)])
					if math.Abs(val) > threshold {
						pos := rl.NewVector3(float32(i)*voxelSize, float32(j)*voxelSize, float32(k)*voxelSize)
						col := rl.NewColor(255, uint8(255*(1-math.Abs(val))), 128, 255) // Red-green by value
						rl.DrawCube(pos, voxelSize, voxelSize, voxelSize, col)
					}
				}
			}
		}
		rl.EndMode3D()

		text := fmt.Sprintf("Time: %.2f | Frame: %d | Auto: %t | Rev: %t", displayT, currFrame, autoAdv, reverse)
		rl.DrawText(text, 10, 10, 20, rl.White)
		rl.EndDrawing()
	}

	rl.CloseWindow()
}
