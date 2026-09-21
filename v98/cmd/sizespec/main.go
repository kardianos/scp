// sizespec — v98 Phase 2 size-spectrum meter (offline, no kernel).
//
//	sizespec FILE.fcs                 # per-frame TSV + histogram
//	sizespec -h0 broad|line6|ceiling FILE.fcs
//	sizespec -h0 broad A.fcs -h0 line6 B.fcs   (one class per file; use once)
package main

import (
	"flag"
	"fmt"
	"os"

	"scp/v97/hunt"
	"scp/v98/size"
)

func main() {
	h0 := flag.String("h0", "", "H0 class: broad|line6|ceiling")
	em := flag.Float64("em", 0.05, "live-Em floor")
	minN := flag.Int("minn", 3, "minimum motif size")
	tag := flag.Bool("tag", false, "restrict census to tagged cells")
	flag.Parse()
	if flag.NArg() < 1 {
		fmt.Fprintln(os.Stderr, "usage: sizespec [-h0 broad|line6|ceiling] [-em 0.05] [-minn 3] [-tag] FILE.fcs")
		os.Exit(2)
	}
	path := flag.Arg(0)
	rep, err := size.Analyze(path, size.Options{
		Census: hunt.CensusCfg{EmFloor: *em, MinN: *minN, TagOnly: *tag},
	})
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
	if *h0 != "" {
		g := size.GradeH0(size.H0Class(*h0), rep)
		fmt.Print(g.String())
		if !g.OK {
			os.Exit(1)
		}
		return
	}
	fmt.Println(rep.SummaryTSVHeader())
	for _, fs := range rep.Snaps {
		fmt.Println(fs.TSV())
	}
	fmt.Fprintf(os.Stderr, "# %s\n", rep.Dump())
}
