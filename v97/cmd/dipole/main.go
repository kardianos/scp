// dipole — v97 conversion-dipole instruments (Go-only, offline).
//
// Reads an FCS v3 stream and writes a per-frame TSV of the five hunt
// meters on the dominant (or -id) motif: census, reconstructed J,
// Es surface flux, frame-diff I/O, I/O separation, winding.
//
//	dipole FILE.fcs                 # TSV on stdout
//	dipole -h0 sealed FILE.fcs      # H0 honesty grade
//	dipole -h0 leak|ring|bathring|current FILE.fcs
//	dipole -h1 FILE.fcs [FILE.fcs ...]   # inhabited / empty hunt
package main

import (
	"flag"
	"fmt"
	"os"

	"scp/v97/hunt"
)

func main() {
	h0 := flag.String("h0", "", "H0 class: sealed|leak|ring|bathring|current")
	h1 := flag.Bool("h1", false, "H1 inhabited/empty grade (accepts many files)")
	em := flag.Float64("em", 0.05, "live-Em floor")
	minN := flag.Int("minn", 3, "minimum motif size")
	id := flag.Int("id", 0, "follow tracked motif id (0 = dominant)")
	tag := flag.Bool("tag", false, "restrict census to tagged cells")
	flag.Parse()
	if flag.NArg() < 1 {
		fmt.Fprintln(os.Stderr, "usage: dipole [-h0 class|-h1] [-em 0.05] [-minn 3] [-id N] [-tag] file.fcs ...")
		os.Exit(2)
	}
	if *h1 {
		code := runH1(flag.Args(), hunt.Options{
			Census:   hunt.CensusCfg{EmFloor: *em, MinN: *minN, TagOnly: *tag},
			FollowID: *id,
		})
		os.Exit(code)
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

func runH1(paths []string, opt hunt.Options) int {
	var nIn, nEmpty, nShort, nErr int
	for _, path := range paths {
		s, err := hunt.Open(path)
		if err != nil {
			fmt.Fprintf(os.Stderr, "H1 ERROR  %s: %v\n", path, err)
			nErr++
			continue
		}
		rows, err := hunt.Analyze(s, opt)
		s.Close()
		if err != nil {
			fmt.Fprintf(os.Stderr, "H1 ERROR  %s: %v\n", path, err)
			nErr++
			continue
		}
		rep := hunt.GradeH1(path, rows)
		fmt.Print(rep.String())
		fmt.Println()
		switch {
		case rep.Inhabited:
			nIn++
		case rep.Span < 100:
			nShort++
		default:
			nEmpty++
		}
	}
	fmt.Printf("H1 SUMMARY  inhabited=%d  empty=%d  short=%d  err=%d  files=%d\n",
		nIn, nEmpty, nShort, nErr, len(paths))
	if nIn > 0 {
		return 0
	}
	return 1
}
