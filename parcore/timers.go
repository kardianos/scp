package parcore

import (
	"fmt"
	"io"
	"sync/atomic"
	"time"
)

// Coarse phase timers. Step accumulates wall time per phase index; the driver
// chooses how many phases to define, so the timer table is sized to a generous
// fixed cap (PhCap). The first len(phases) entries are meaningful; the rest
// stay zero. Phase PhBoundary is reserved as the last named slot.
const (
	PhBoundary = 63 // reserved slot for KindBoundary phases
	PhCap      = 64
)

var (
	phaseOn atomic.Bool
	phaseNs [PhCap]int64
	phaseN  [PhCap]int64
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
	if id < 0 || id >= PhCap {
		return
	}
	atomic.AddInt64(&phaseNs[id], d.Nanoseconds())
	atomic.AddInt64(&phaseN[id], 1)
}

// PhaseTimerReport writes the accumulated per-phase wall time to w. names[id]
// labels phase id; pass the same names your driver used when it built Step().
func PhaseTimerReport(w io.Writer, steps int, names []string) {
	var total int64
	for i := 0; i < PhCap; i++ {
		total += phaseNs[i]
	}
	if total == 0 {
		fmt.Fprintln(w, "(no phase timer data)")
		return
	}
	fmt.Fprintf(w, "phase breakdown (sum over %d steps):\n", steps)
	for i := 0; i < PhCap; i++ {
		ns := phaseNs[i]
		if ns == 0 {
			continue
		}
		label := ""
		if i < len(names) {
			label = names[i]
		}
		if i == PhBoundary {
			label = "boundary-hops"
		}
		if label == "" {
			label = fmt.Sprintf("phase%d", i)
		}
		pct := 100 * float64(ns) / float64(total)
		fmt.Fprintf(w, "  %-26s %9.1f ms total  %6.3f ms/step  %5.1f%%\n",
			label, float64(ns)/1e6, float64(ns)/1e6/float64(steps), pct)
	}
	fmt.Fprintf(w, "  %-26s %9.1f ms total\n", "SUM", float64(total)/1e6)
}
