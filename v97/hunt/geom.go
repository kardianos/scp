package hunt

import "math"

func wrap(dx, L float64) float64 {
	if L <= 0 {
		return dx
	}
	return dx - L*math.Round(dx/L)
}

func wrapPi(a float64) float64 {
	const twoPi = 2 * math.Pi
	a = math.Mod(a+math.Pi, twoPi)
	if a < 0 {
		a += twoPi
	}
	return a - math.Pi
}

func hypot3(x, y, z float64) float64 {
	return math.Sqrt(x*x + y*y + z*z)
}

// wrapped centroid of weights at positions. origin is a reference
// point (typically the first live cell) so the average is unambiguous
// in a periodic box.
func centroid(xs, ys, zs, w []float64, L float64) (cx, cy, cz, wt float64) {
	if len(xs) == 0 {
		return
	}
	ox, oy, oz := xs[0], ys[0], zs[0]
	for i := range xs {
		if w[i] <= 0 {
			continue
		}
		cx += w[i] * wrap(xs[i]-ox, L)
		cy += w[i] * wrap(ys[i]-oy, L)
		cz += w[i] * wrap(zs[i]-oz, L)
		wt += w[i]
	}
	if wt <= 0 {
		return ox, oy, oz, 0
	}
	cx = ox + cx/wt
	cy = oy + cy/wt
	cz = oz + cz/wt
	if L > 0 {
		cx = cx - L*math.Floor(cx/L)
		cy = cy - L*math.Floor(cy/L)
		cz = cz - L*math.Floor(cz/L)
	}
	return cx, cy, cz, wt
}

func dist3(ax, ay, az, bx, by, bz, L float64) float64 {
	return hypot3(wrap(ax-bx, L), wrap(ay-by, L), wrap(az-bz, L))
}
