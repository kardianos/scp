// dipole — v97 conversion-dipole instruments (Go-only, offline).
//
// Reads an FCS v3 stream and writes a per-frame TSV of the five hunt
// meters on the dominant (or -id) motif: census, reconstructed J,
// Es surface flux, frame-diff I/O, I/O separation, winding.
//
//	dipole FILE.fcs                 # TSV on stdout
//	dipole -h0 sealed FILE.fcs      # H0 honesty grade
//	dipole -h0 leak|ring|bathring|current FILE.fcs
package main

import (
	"flag"
	"fmt"
	"os"

	"scp/v97/hunt"
)

func main() {
	h0 := flag.String("h0", "", "H0 class: sealed|leak|ring|bathring|current")
	em := flag.Float64("em", 0.05, "live-Em floor")
	minN := flag.Int("minn", 3, "minimum motif size")
	id := flag.Int("id", 0, "follow tracked motif id (0 = dominant)")
	tag := flag.Bool("tag", false, "restrict census to tagged cells")
	flag.Parse()
	if flag.NArg() < 1 {
		fmt.Fprintln(os.Stderr, "usage: dipole [-h0 class] [-em 0.05] [-minn 3] [-id N] file.fcs")
		os.Exit(2)
	}
	path := flag.Arg(0)
	s, err := hunt.Open(path)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
	defer s.Close()
	rows, err := hunt.Analyze(s, hunt.Options{
		Census:   hunt.CensusCfg{EmFloor: *em, MinN: *minN, TagOnly: *tag},
		FollowID: *id,
	})
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
	if *h0 != "" {
		rep := hunt.GradeH0(hunt.H0Class(*h0), path, rows)
		fmt.Print(rep.String())
		if !rep.OK {
			os.Exit(1)
		}
		return
	}
	fmt.Println(hunt.RowHeader())
	for _, r := range rows {
		fmt.Println(r.TSV())
	}
}
