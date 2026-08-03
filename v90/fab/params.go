package fab

// Config — LAW KEYS default to laws_V2g.cfg (2026-07-28) VERBATIM.
// Port of the Cfg struct + parser in v89/freecell.c. The law is the
// table; geometry/run keys are apparatus.

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"strconv"
	"strings"
)

const TwoPi = 2.0 * math.Pi

type Cfg struct {
	Seed uint64
	// --- laws_V2g (do not change defaults: they ARE the table) ---
	C, Dmin, R0, Rjit                       float64
	W1, W2, QDetune                         float64
	GammaRes, GammaResM, PGate, LockFloor   float64
	KDep, KDepM, Cap                        float64
	Es0, EsFloor, ECond                     float64
	FConv, FEvap, SPull                     float64
	KappaLock, KappaAlign, KappaFreq, KappaReac float64
	SK, SDisp, SigmaTumble                  float64
	CombLimit                               int
	RoughK, GammaRough                      float64
	MobSym                                  int
	MobFloor, FieldJ, QuantA0               float64
	QuantMode                               int
	// --- freecell geometry sector (apparatus, not law) ---
	Cfac      float64
	KRep      float64
	MobGeo    float64
	KappaBond float64
	JamSweeps int
	JamK      float64
	JamTol    float64
	FreezeGeo int
	BathFrac  float64
	// --- run + apparatus ---
	L, Dt, T             float64
	DiagEvery, SnapEvery int
	Exp                  string
	NoiseAmp             float64
	Amp, Sigma           float64
	Kx                   float64
	Px, Py, Pz           float64
	Bath                 int
	PairX0, PairX1, PairDoff float64
	PairM                int
	PairPP, PairQQ       int
	Seedlock             int
	RingN                int
	RingX, RingDoff      float64
	TriKind, TriBranch   int
	TriXU, TriXD, TriDoff float64
	Tri2Sep              float64
	Tri2K2               int
	OctX, OctDoff        float64
	ShearEps, ShearT     float64
	SnapDir              string
	SnapFile             string
}

func Defaults() Cfg {
	var p Cfg
	p.Seed = 20260802
	// laws_V2g.cfg verbatim
	p.C = 1.0
	p.Dmin = 1.0
	p.R0 = 0.85
	p.Rjit = 0.06
	p.W1 = 1.65
	p.W2 = 2.9
	p.QDetune = 1.2
	p.GammaRes = 0.25
	p.GammaResM = 0.10
	p.PGate = 8
	p.LockFloor = 0.005
	p.KDep = 1.2
	p.KDepM = 2.0
	p.Cap = 2.5
	p.Es0 = 1.0
	p.EsFloor = 0.05
	p.ECond = 0
	p.FConv = 0.25
	p.FEvap = 0.5
	p.SPull = 0.15
	p.KappaLock = 1.0
	p.KappaAlign = 0.5
	p.KappaFreq = 0
	p.KappaReac = 0
	p.SK = 0.06
	p.SDisp = 0.3
	p.SigmaTumble = 0.01
	p.CombLimit = 6
	p.RoughK = 0.35
	p.GammaRough = 0.5
	p.MobSym = 1
	p.MobFloor = 0.004
	p.FieldJ = 1.8
	p.QuantA0 = 1.15
	p.QuantMode = 2
	// freecell geometry
	p.Cfac = 1.15
	p.KRep = 1.0
	p.MobGeo = 1.0
	p.KappaBond = 1.0
	p.JamSweeps = 800
	p.JamK = 0.10
	p.JamTol = 1e-4
	p.FreezeGeo = 0
	p.BathFrac = 1.0
	// run
	p.L = 16.0
	p.Dt = 0.02
	p.T = 40.0
	p.DiagEvery = 100
	p.SnapEvery = 0
	p.Exp = "bath"
	p.NoiseAmp = 0
	p.Amp = 1.6
	p.Sigma = 2.2
	p.Kx = 0
	p.Px = -1
	p.Py = -1
	p.Pz = -1
	p.Bath = 1
	p.PairX0 = 0.325
	p.PairX1 = -1
	p.PairDoff = 0.05
	p.PairM = 1
	p.PairPP = 1
	p.PairQQ = 1
	p.Seedlock = 1
	p.RingN = 6
	p.RingX = 0.325
	p.RingDoff = 0.0
	p.TriKind = 0
	p.TriBranch = 0
	p.TriXU = 0.20
	p.TriXD = -1
	p.TriDoff = 0
	p.Tri2Sep = 2.6
	p.Tri2K2 = 0
	p.OctX = 0.325
	p.OctDoff = 0.0
	p.ShearEps = 0
	p.ShearT = 0
	p.SnapDir = ""
	p.SnapFile = ""
	return p
}

func atof(v string) float64 {
	f, _ := strconv.ParseFloat(strings.TrimSpace(v), 64)
	return f
}

func atoi(v string) int {
	// C atoi semantics: leading integer, 0 on garbage
	f, err := strconv.Atoi(strings.TrimSpace(v))
	if err != nil {
		g, _ := strconv.ParseFloat(strings.TrimSpace(v), 64)
		return int(g)
	}
	return f
}

func (p *Cfg) SetKV(k, v string) {
	switch k {
	case "seed":
		u, _ := strconv.ParseUint(strings.TrimSpace(v), 10, 64)
		p.Seed = u
	case "C":
		p.C = atof(v)
	case "dmin":
		p.Dmin = atof(v)
	case "r0":
		p.R0 = atof(v)
	case "rjit":
		p.Rjit = atof(v)
	case "w1":
		p.W1 = atof(v)
	case "w2":
		p.W2 = atof(v)
	case "q_detune":
		p.QDetune = atof(v)
	case "gamma_res":
		p.GammaRes = atof(v)
	case "gamma_res_m":
		p.GammaResM = atof(v)
	case "p_gate":
		p.PGate = atof(v)
	case "lock_floor":
		p.LockFloor = atof(v)
	case "k_dep":
		p.KDep = atof(v)
	case "k_dep_m":
		p.KDepM = atof(v)
	case "cap":
		p.Cap = atof(v)
	case "e_s0":
		p.Es0 = atof(v)
	case "es_floor":
		p.EsFloor = atof(v)
	case "e_cond":
		p.ECond = atof(v)
	case "f_conv":
		p.FConv = atof(v)
	case "f_evap":
		p.FEvap = atof(v)
	case "s_pull":
		p.SPull = atof(v)
	case "kappa_lock":
		p.KappaLock = atof(v)
	case "kappa_align":
		p.KappaAlign = atof(v)
	case "kappa_freq":
		p.KappaFreq = atof(v)
	case "kappa_reac":
		p.KappaReac = atof(v)
	case "s_k":
		p.SK = atof(v)
	case "s_disp":
		p.SDisp = atof(v)
	case "sigma_tumble":
		p.SigmaTumble = atof(v)
	case "comb_limit":
		p.CombLimit = atoi(v)
	case "rough_k":
		p.RoughK = atof(v)
	case "gamma_rough":
		p.GammaRough = atof(v)
	case "mob_sym":
		p.MobSym = atoi(v)
	case "mob_floor":
		p.MobFloor = atof(v)
	case "field_J":
		p.FieldJ = atof(v)
	case "quant_A0":
		p.QuantA0 = atof(v)
	case "quant_mode":
		p.QuantMode = atoi(v)
	case "cfac":
		p.Cfac = atof(v)
	case "k_rep":
		p.KRep = atof(v)
	case "mob_geo":
		p.MobGeo = atof(v)
	case "kappa_bond":
		p.KappaBond = atof(v)
	case "jam_sweeps":
		p.JamSweeps = atoi(v)
	case "jam_k":
		p.JamK = atof(v)
	case "jam_tol":
		p.JamTol = atof(v)
	case "freeze_geo":
		p.FreezeGeo = atoi(v)
	case "bath_frac":
		p.BathFrac = atof(v)
	case "L":
		p.L = atof(v)
	case "dt":
		p.Dt = atof(v)
	case "T":
		p.T = atof(v)
	case "diag_every":
		p.DiagEvery = atoi(v)
	case "snap_every":
		p.SnapEvery = atoi(v)
	case "exp":
		p.Exp = v
	case "noise_amp":
		p.NoiseAmp = atof(v)
	case "amp":
		p.Amp = atof(v)
	case "sigma":
		p.Sigma = atof(v)
	case "kx":
		p.Kx = atof(v)
	case "px":
		p.Px = atof(v)
	case "py":
		p.Py = atof(v)
	case "pz":
		p.Pz = atof(v)
	case "bath":
		p.Bath = atoi(v)
	case "pair_x0":
		p.PairX0 = atof(v)
	case "pair_x1":
		p.PairX1 = atof(v)
	case "pair_doff":
		p.PairDoff = atof(v)
	case "pair_m":
		p.PairM = atoi(v)
	case "pair_p":
		p.PairPP = atoi(v)
	case "pair_q":
		p.PairQQ = atoi(v)
	case "seedlock":
		p.Seedlock = atoi(v)
	case "ring_n":
		p.RingN = atoi(v)
	case "tri_kind":
		p.TriKind = atoi(v)
	case "tri_branch":
		p.TriBranch = atoi(v)
	case "tri_xU":
		p.TriXU = atof(v)
	case "tri_xD":
		p.TriXD = atof(v)
	case "tri_doff":
		p.TriDoff = atof(v)
	case "tri2_sep":
		p.Tri2Sep = atof(v)
	case "tri2_k2":
		p.Tri2K2 = atoi(v)
	case "ring_x":
		p.RingX = atof(v)
	case "ring_doff":
		p.RingDoff = atof(v)
	case "oct_x":
		p.OctX = atof(v)
	case "oct_doff":
		p.OctDoff = atof(v)
	case "shear_eps":
		p.ShearEps = atof(v)
	case "shear_t":
		p.ShearT = atof(v)
	case "snap_dir":
		p.SnapDir = v
	case "snap_file":
		p.SnapFile = v
	default:
		fmt.Fprintf(os.Stderr, "# WARN unknown key %s\n", k)
	}
}

func (p *Cfg) LoadCfg(path string) error {
	f, err := os.Open(path)
	if err != nil {
		return fmt.Errorf("cannot open %s", path)
	}
	defer f.Close()
	sc := bufio.NewScanner(f)
	for sc.Scan() {
		line := sc.Text()
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		eq := strings.IndexByte(line, '=')
		if eq < 0 {
			continue
		}
		k := strings.TrimSpace(line[:eq])
		v := line[eq+1:]
		if hs := strings.IndexByte(v, '#'); hs >= 0 {
			v = v[:hs]
		}
		v = strings.TrimRight(strings.TrimLeft(v, " "), " \r")
		if k != "" {
			p.SetKV(k, v)
		}
	}
	return sc.Err()
}
