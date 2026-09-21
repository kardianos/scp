package size

import (
	"fmt"
	"strings"
)

// H0Class is a pre-registered honesty class. Bars are in GradeH0.
// Sharpen by measurement; do not soften to pass.
type H0Class string

const (
	H0Broad   H0Class = "broad"   // Face C: hierarchy, not a line
	H0Line6   H0Class = "line6"   // vacuum C6: N=6 with C6 spectrum
	H0Ceiling H0Class = "ceiling" // tax-only: drain, still not a line
)

// H0Result is one graded stream.
type H0Result struct {
	Class   H0Class
	Path    string
	OK      bool
	Verdict string
	Detail  []string
	Report  *Report
}

func (r H0Result) String() string {
	var b strings.Builder
	mark := "FAIL"
	if r.OK {
		mark = "PASS"
	}
	fmt.Fprintf(&b, "H0 %s %s  %s\n", r.Class, mark, r.Path)
	fmt.Fprintf(&b, "  %s\n", r.Verdict)
	for _, d := range r.Detail {
		fmt.Fprintf(&b, "  %s\n", d)
	}
	if r.Report != nil {
		fmt.Fprintf(&b, "  %s\n", r.Report.dumpSummary())
	}
	return b.String()
}

// GradeH0 applies the named class to an already-built report.
//
// broad:
//
//	FAIL as a line if one N holds ≥ 0.70 of motif-counts AND
//	mean(N_eff/N) ≥ 0.50 (a compact species-shaped peak).
//	PASS as broad if not a line AND (distinct N ≥ 5
//	OR top-bin < 0.40 OR mean(N_eff/N) < 0.50).
//	The last clause is the Face C connected-clump case:
//	few components, internally graded.
//
// line6:
//
//	PASS if dominant N=6 on ≥ 0.80 of frames AND, among
//	those N=6 frames with a chart, Match6 ≥ 0.80.
//
// ceiling:
//
//	PASS if EmLast ≤ 0.85 EmFirst (inventory drops) AND
//	the stream is not a line under the broad/line test.
func GradeH0(class H0Class, rep *Report) H0Result {
	out := H0Result{Class: class, Path: rep.Path, Report: rep}
	switch class {
	case H0Broad:
		out.OK, out.Verdict, out.Detail = gradeBroad(rep)
	case H0Line6:
		out.OK, out.Verdict, out.Detail = gradeLine6(rep)
	case H0Ceiling:
		out.OK, out.Verdict, out.Detail = gradeCeiling(rep)
	default:
		out.Verdict = fmt.Sprintf("unknown class %q", class)
	}
	return out
}

func isLine(rep *Report) (bool, string) {
	top := rep.topBinFrac()
	ne := rep.meanNeffOverN()
	if top >= 0.70 && ne >= 0.50 {
		return true, fmt.Sprintf("line: top-bin=%.3f N_eff/N=%.3f (species-shaped)", top, ne)
	}
	return false, fmt.Sprintf("not-line: top-bin=%.3f N_eff/N=%.3f distinct=%d", top, ne, rep.distinctN())
}

func gradeBroad(rep *Report) (ok bool, verd string, det []string) {
	line, why := isLine(rep)
	det = append(det, why)
	det = append(det, fmt.Sprintf("hist N  %s", rep.histLine()))
	if line {
		return false, "read as a LINE (Face C must be broad)", det
	}
	ne := rep.meanNeffOverN()
	if rep.distinctN() >= 5 || rep.topBinFrac() < 0.40 || ne < 0.50 {
		return true, "broad (hierarchy / graded clump, not a line)", det
	}
	return false, "not a line, but too narrow to call broad", det
}

func gradeLine6(rep *Report) (ok bool, verd string, det []string) {
	f6 := rep.fracDomN(6)
	m6 := rep.fracMatch6()
	det = append(det, fmt.Sprintf("frac dominant N=6: %.3f (bar ≥ 0.80)", f6))
	det = append(det, fmt.Sprintf("C6 spectrum match among N=6: %.3f (bar ≥ 0.80)", m6))
	det = append(det, fmt.Sprintf("hist N  %s", rep.histLine()))
	if f6 < 0.80 {
		return false, "dominant motif is not N=6", det
	}
	if m6 < 0.80 {
		return false, "N=6 but Laplacian is not the C6 book", det
	}
	return true, "line at N=6 with C6 spectrum", det
}

func gradeCeiling(rep *Report) (ok bool, verd string, det []string) {
	var ratio float64
	if rep.EmFirst > 0 {
		ratio = rep.EmLast / rep.EmFirst
	}
	det = append(det, fmt.Sprintf("EmLast/EmFirst=%.3f (bar ≤ 0.85)", ratio))
	line, why := isLine(rep)
	det = append(det, why)
	det = append(det, fmt.Sprintf("hist N  %s", rep.histLine()))
	if ratio > 0.85 {
		return false, "inventory did not drop (not a tax ceiling)", det
	}
	if line {
		return false, "drained onto a line (tax selected a species — not the Face-C/MONO fact)", det
	}
	return true, "ceiling: drain without a line", det
}
