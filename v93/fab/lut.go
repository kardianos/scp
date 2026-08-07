package fab

// Trig kernel — port of v89/int_ledger/trig_lut.inc default path
// (lut_cos / lut_sin). Pure double arithmetic + fmod: bit-identical to
// the C instrument by construction (no libm involvement, no FMA — the
// C reference is built without FMA-capable codegen and gc does not
// contract on amd64).

import "math"

const (
	trigPio2 = 1.57079632679489661923
	trigPix2 = 6.28318530717958647692
)

var trigCosC = [6]float64{
	4.16666666666666019e-02,
	-1.38888888888741096e-03,
	2.48015872894767294e-05,
	-2.75573143513906633e-07,
	2.08757232129817483e-09,
	-1.13596475577881948e-11,
}

var trigSinS = [6]float64{
	-1.66666666666666324e-01,
	8.33333333332248946e-03,
	-1.98412698298579493e-04,
	2.75573137070700677e-06,
	-2.50507602534068634e-08,
	1.58969099521155010e-10,
}

func trigKernelCos(x float64) float64 {
	z := x * x
	r := trigCosC[0] + z*(trigCosC[1]+z*(trigCosC[2]+
		z*(trigCosC[3]+z*(trigCosC[4]+z*trigCosC[5]))))
	return 1.0 - 0.5*z + z*z*r
}

func trigKernelSin(x float64) float64 {
	z := x * x
	v := z * x
	r := trigSinS[0] + z*(trigSinS[1]+z*(trigSinS[2]+
		z*(trigSinS[3]+z*(trigSinS[4]+z*trigSinS[5]))))
	return x + v*r
}

func lutCos(a float64) float64 {
	x := math.Mod(a, trigPix2)
	if x < 0 {
		x += trigPix2
	}
	quad := int(x / trigPio2)
	t := x - float64(quad)*trigPio2
	var c, s float64
	if t <= trigPio2*0.5 {
		c = trigKernelCos(t)
		s = trigKernelSin(t)
	} else {
		u := trigPio2 - t
		c = trigKernelSin(u)
		s = trigKernelCos(u)
	}
	switch quad & 3 {
	case 0:
		return c
	case 1:
		return -s
	case 2:
		return -c
	default:
		return s
	}
}

func lutSin(a float64) float64 {
	return lutCos(a - trigPio2)
}
