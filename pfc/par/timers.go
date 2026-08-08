package par

import (
	"fmt"
	"io"
	"sync/atomic"
	"time"
)

// Coarse phases after the fuse rewrite (fewer barriers).
const (
	PhPullGeom = iota // pull ghosts + pitch + geom + spaceW + scale
	PhExch            // sscl exchange
	PhWork            // apply + area + field + hops + pub Fa
	PhBoundary        // serial boundary hops
	PhFinish          // merge Fa + force + motion + face publish
	PhCount
)

var (
	phaseOn atomic.Bool
	phaseNs [PhCount]int64
	phaseN  [PhCount]int64
)

// EnablePhaseTimers toggles accumulation of per-phase wall time in Step.
func EnablePhaseTimers(on bool) {
	phaseOn.Store(on)
	if on {
		for i := range phaseNs {
			phaseNs[i] = 0
			phaseN[i] = 0
		}
	}
}

func phaseAdd(id int, d time.Duration) {
	if !phaseOn.Load() {
		return
	}
	atomic.AddInt64(&phaseNs[id], d.Nanoseconds())
	atomic.AddInt64(&phaseN[id], 1)
}

// PrintPhaseTimers writes average phase costs to w.
func PrintPhaseTimers(w io.Writer, steps int) {
	names := [...]string{
		PhPullGeom: "pull+pitch+geom+space",
		PhExch:     "exch_sscl",
		PhWork:     "apply+field+hops",
		PhBoundary: "boundaryHops",
		PhFinish:   "force+motion+publish",
	}
	var total int64
	for i := 0; i < PhCount; i++ {
		total += phaseNs[i]
	}
	if total == 0 {
		fmt.Fprintln(w, "(no phase timer data)")
		return
	}
	fmt.Fprintf(w, "phase breakdown (sum over %d steps):\n", steps)
	for i := 0; i < PhCount; i++ {
		ns := phaseNs[i]
		pct := 100 * float64(ns) / float64(total)
		fmt.Fprintf(w, "  %-22s %8.1f ms total  %6.2f ms/step  %5.1f%%\n",
			names[i], float64(ns)/1e6, float64(ns)/1e6/float64(steps), pct)
	}
	fmt.Fprintf(w, "  %-22s %8.1f ms total\n", "SUM_phases", float64(total)/1e6)
}
