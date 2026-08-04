package fab

// Sim state + allocation + small math — port of v89/freecell.c
// (state block, rng, gate, comb, atoms machinery, field_inject).
// The C globals live on the Sim struct; everything else is 1:1.

import (
	"fmt"
	"io"
	"math"
	"os"
)

const (
	sFree  = 0
	sAlive = 1
	sDying = 2
)

const hEmpty = int64(-1)
const hTomb = int64(-2)

type Sim struct {
	P Cfg

	Out  io.Writer // stdout of the instrument (log lines)
	Errw io.Writer // stderr (warnings, autopsies)

	rng uint64

	NC int
	// cells
	px, py, pz          []float64 // POSITIONS — dynamical here
	cr0, cr             []float64
	n1x, n1y, n1z       []float64
	n2x, n2y, n2z       []float64
	th2, cbeta          []float64
	w1e, w2e            []float64
	Es, Em, Ee          []float64
	fa1, fa2            []float64
	flload              []float64
	qcnvD, qcnvF        []float64
	roughq              []float64
	req1, scl1          []float64
	sprq, sscl          []float64
	fsum                []float64
	tag                 []uint8
	pin                 []uint8   // wall fixture: pass D skips (slit)
	scond               []uint8   // condensation-active override (clicks)
	fxb, fyb, fzb       []float64 // geometric force gather buffers
	rngbuf, nsnap, th2s []float64

	// channels — persistent slots with identity (birth/death ledger)
	NLMAX, NSLOT                int
	sli, slj                    []int
	sst                         []uint8
	sd                          []float64
	sux, suy, suz               []float64
	sA                          []float64
	slem, slph                  []float64 // [slot][dir]
	slp, slq                    []int8
	swant                       []float64 // [slot][dir]
	sflux                       []float64
	sldd                        []float64
	swl                         []float64
	sfluxd                      []float64 // cumulative DIRECTED dense flux [slot][dir]
	freelist                    []int
	nfree                       int
	births, deaths, betaReturns int64
	betaEnergy                  float64

	// pair -> slot hash (open addressing, tombstones)
	hkey  []int64
	hval  []int
	hsize int

	// cell -> incident slots (CSR, canonical)
	cls   []int
	clidx []int

	// spatial hash
	binHead []int
	binNext []int
	binG    int
	binSz   float64

	// conversion ledgers
	roughTotal, condTotal, evapTotal, backsTotal float64
	radTotal                                     float64 // v91 graded sub-cap radiance
	A0eff                                        float64
	qfireN                                       int64
	simT                                         float64
	E0                                           float64
	cenx, ceny, cenz                             float64
	ncandLast                                    int64
	fsT, fsR, fsY, fsZ                           [512]float64
	nfsamp                                       int

	// derived law constants
	gm, lockf, Aref, dref float64

	// the partial comb
	ncomb        int
	combp, combq [24]int

	// pair/truss registry
	repPairI, repPairJ int
	trussN             int
	trussE             [128][2]int
	trussDstar         [128]float64
	trussShape0        [3]float64
	// FCQ triangle registry
	triOn      bool
	triV       [2][3]int
	triP, triQ [2][3]int8
	ntri       int

	// COMPOSITE inter-ring boundary edges (directed, with charts)
	rintN      int
	rintE      [8][2]int
	rintP      [8]int8
	rintQ      [8]int8
	rintDstar  [8]float64
	ringsNring int
	ringsV0    [6]int
	ringsMode  [6]int
	ringsW     [6]float64
	ringsGrp   [6]int

	// DS screen (exp=slit): time-integrated Ee per y-bin + per cell
	dsI     [dsNBin]float64
	dsExpo  float64
	dsCellI []float64
	// DS tier-1 clicks
	clickT, clickY, clickE [clickMax]float64
	nclick                 int
	// blob2 per-blob COM samples
	b2ax, b2bx, b2tx [512]float64
	// XSEC sector meter + tag-split conversion ledger
	sectE, sectN2                    [sectMax]float64
	sectCx, sectCy                   float64
	ctRough, ctCond, ctEvap, ctBacks float64
	// p1_meter accumulators (cumulative first-moment flow per channel)
	p1sp, p1fl, p1fd, p1gm [3]float64

	// scratch (lazy)
	comp, stack []int
}

const (
	dsNBin   = 96
	clickMax = 8192
	sectMax  = 96
)

// ------------------------------------------------------------------
// rng + small math (cellfab-identical xorshift64)
// ------------------------------------------------------------------

func (s *Sim) xrand() uint64 {
	x := s.rng
	x ^= x << 13
	x ^= x >> 7
	x ^= x << 17
	s.rng = x
	return x
}

func (s *Sim) frand() float64 {
	return float64(s.xrand()>>11) * (1.0 / 9007199254740992.0)
}

func (s *Sim) grand() float64 {
	u1 := s.frand() + 1e-18
	u2 := s.frand()
	return math.Sqrt(-2.0*math.Log(u1)) * math.Cos(TwoPi*u2)
}

func wrapPi(a float64) float64 {
	a = math.Mod(a+math.Pi, TwoPi)
	if a < 0 {
		a += TwoPi
	}
	return a - math.Pi
}

// periodic box: minimum-image separation and position fold
func (s *Sim) wr(d float64) float64 {
	if d > 0.5*s.P.L {
		d -= s.P.L
	}
	if d < -0.5*s.P.L {
		d += s.P.L
	}
	return d
}

func (s *Sim) fold(x float64) float64 {
	x = math.Mod(x, s.P.L)
	if x < 0 {
		x += s.P.L
	}
	return x
}

func (s *Sim) gateOf(dphi float64) float64 { // cellfab.c:638 verbatim
	g := 0.5 * (1.0 + lutCos(dphi))
	ip := int(s.P.PGate)
	if float64(ip) == s.P.PGate && ip >= 1 && ip <= 8 {
		r := 1.0
		for q := 0; q < ip; q++ {
			r *= g
		}
		return r
	}
	return math.Pow(g, s.P.PGate)
}

func (s *Sim) randUnit() (x, y, z float64) {
	var a, b, c, n2 float64
	for {
		a, b, c = s.grand(), s.grand(), s.grand()
		n2 = a*a + b*b + c*c
		if n2 >= 1e-12 {
			break
		}
	}
	inv := 1.0 / math.Sqrt(n2)
	return a * inv, b * inv, c * inv
}

// the partial comb — cellfab.c:589 verbatim
func gcdI(a, b int) int {
	for b != 0 {
		a, b = b, a%b
	}
	return a
}

func (s *Sim) combBuild() {
	s.ncomb = 0
	for pp := 1; pp <= 6; pp++ {
		for qq := 1; qq <= 6; qq++ {
			if pp*qq > s.P.CombLimit {
				continue
			}
			if gcdI(pp, qq) != 1 {
				continue
			}
			s.combp[s.ncomb] = pp
			s.combq[s.ncomb] = qq
			s.ncomb++
		}
	}
	for a := 0; a < s.ncomb; a++ {
		for b := a + 1; b < s.ncomb; b++ {
			if s.combp[b]*s.combq[b] < s.combp[a]*s.combq[a] {
				s.combp[a], s.combp[b] = s.combp[b], s.combp[a]
				s.combq[a], s.combq[b] = s.combq[b], s.combq[a]
			}
		}
	}
}

// ------------------------------------------------------------------
// atoms at mode boundaries — cellfab.c:2739 verbatim
// ------------------------------------------------------------------

func (s *Sim) atomsEps(epsSrc, epsDst float64) float64 {
	if s.P.QuantMode >= 3 {
		return epsDst
	}
	return epsSrc
}

func (s *Sim) atomsW(wSrc, wDst float64) float64 {
	if s.P.QuantMode >= 3 {
		return wDst
	}
	return wSrc
}

func (s *Sim) atomsFire(demand, epsSrc, epsDst float64, cred *float64) float64 {
	if s.A0eff <= 0 || demand <= 0 {
		return demand
	}
	eps := s.atomsEps(epsSrc, epsDst)
	if eps <= 0 {
		return demand
	}
	if s.P.QuantMode == 2 || s.P.QuantMode == 4 {
		*cred += demand
		if *cred > 2.0*eps {
			*cred = 2.0 * eps
		}
		mv := math.Floor(*cred/eps) * eps
		*cred -= mv
		return mv
	}
	return math.Floor(demand/eps) * eps
}

func (s *Sim) atomsClamp(mv, ceilE, epsSrc, epsDst float64, cred *float64) float64 {
	if mv <= ceilE {
		return mv
	}
	if s.A0eff <= 0 {
		return ceilE
	}
	eps := s.atomsEps(epsSrc, epsDst)
	keep := ceilE
	if eps > 0 {
		keep = math.Floor(ceilE/eps) * eps
	}
	if (s.P.QuantMode == 2 || s.P.QuantMode == 4) && cred != nil {
		*cred += mv - keep
	}
	return keep
}

func (s *Sim) qatomDiag(fd int, w, e float64) {
	if s.A0eff <= 0 || e <= 0 {
		return
	}
	n := s.qfireN
	s.qfireN++
	if n%200 == 0 {
		dir := "DF"
		if fd != 0 {
			dir = "FD"
		}
		fmt.Fprintf(s.Out, "# QATOM t=%.2f dir=%s w=%.9g e=%.12g\n", s.simT, dir, w, e)
	}
}

func (s *Sim) fieldInject(i int, dE float64) { // cellfab.c:2711 verbatim
	if dE <= 0 {
		return
	}
	e := s.fa1[i]*s.fa1[i] + s.fa2[i]*s.fa2[i]
	if e > 1e-20 {
		fac := math.Sqrt((e + dE) / e)
		s.fa1[i] *= fac
		s.fa2[i] *= fac
	} else {
		s.fa1[i] = math.Sqrt(dE) * lutCos(s.th2[i])
		s.fa2[i] = math.Sqrt(dE) * lutSin(s.th2[i])
	}
	s.Ee[i] = e + dE
}

// ------------------------------------------------------------------
// allocation
// ------------------------------------------------------------------

func (s *Sim) allocAll(nc int) {
	s.NC = nc
	s.NLMAX = 80*nc + 1024
	s.px = make([]float64, nc)
	s.py = make([]float64, nc)
	s.pz = make([]float64, nc)
	s.cr0 = make([]float64, nc)
	s.cr = make([]float64, nc)
	s.n1x = make([]float64, nc)
	s.n1y = make([]float64, nc)
	s.n1z = make([]float64, nc)
	s.n2x = make([]float64, nc)
	s.n2y = make([]float64, nc)
	s.n2z = make([]float64, nc)
	s.th2 = make([]float64, nc)
	s.cbeta = make([]float64, nc)
	s.w1e = make([]float64, nc)
	s.w2e = make([]float64, nc)
	s.Es = make([]float64, nc)
	s.Em = make([]float64, nc)
	s.Ee = make([]float64, nc)
	s.fa1 = make([]float64, nc)
	s.fa2 = make([]float64, nc)
	s.flload = make([]float64, nc)
	s.qcnvD = make([]float64, nc)
	s.qcnvF = make([]float64, nc)
	s.roughq = make([]float64, nc)
	s.req1 = make([]float64, nc)
	s.scl1 = make([]float64, nc)
	s.sprq = make([]float64, nc)
	s.sscl = make([]float64, nc)
	s.fsum = make([]float64, nc)
	s.tag = make([]uint8, nc)
	s.pin = make([]uint8, nc)
	s.scond = make([]uint8, nc)
	s.fxb = make([]float64, nc)
	s.fyb = make([]float64, nc)
	s.fzb = make([]float64, nc)
	s.rngbuf = make([]float64, 6*nc)
	s.nsnap = make([]float64, 6*nc)
	s.th2s = make([]float64, nc)
	s.cls = make([]int, nc+1)
	s.clidx = make([]int, 2*s.NLMAX)

	s.sli = make([]int, s.NLMAX)
	s.slj = make([]int, s.NLMAX)
	s.sst = make([]uint8, s.NLMAX)
	s.sd = make([]float64, s.NLMAX)
	s.sux = make([]float64, s.NLMAX)
	s.suy = make([]float64, s.NLMAX)
	s.suz = make([]float64, s.NLMAX)
	s.sA = make([]float64, s.NLMAX)
	s.slem = make([]float64, 2*s.NLMAX)
	s.slph = make([]float64, 2*s.NLMAX)
	s.slp = make([]int8, s.NLMAX)
	s.slq = make([]int8, s.NLMAX)
	s.swant = make([]float64, 2*s.NLMAX)
	s.sflux = make([]float64, s.NLMAX)
	s.sldd = make([]float64, s.NLMAX)
	s.swl = make([]float64, s.NLMAX)
	s.sfluxd = make([]float64, 2*s.NLMAX)
	s.freelist = make([]int, s.NLMAX)
	s.nfree = 0
	for sl := s.NLMAX - 1; sl >= 0; sl-- {
		s.freelist[s.nfree] = sl
		s.nfree++
	}
	s.NSLOT = 0

	s.hsize = 1
	for s.hsize < 4*s.NLMAX {
		s.hsize <<= 1
	}
	s.hkey = make([]int64, s.hsize)
	s.hval = make([]int, s.hsize)
	for h := 0; h < s.hsize; h++ {
		s.hkey[h] = hEmpty
	}

	s.binHead = nil
	s.binNext = make([]int, nc)
	s.binG = 0
}

func (s *Sim) fatal(code int, format string, args ...any) {
	fmt.Fprintf(s.Errw, format, args...)
	os.Exit(code)
}
