// fabrun — CLI for the v90 Go free-cell kernel (fab package).
// Same config surface as the C instrument:
//   fabrun [config.cfg] [key=value ...]
package main

import (
	"fmt"
	"os"
	"strings"

	"scp/v90/fab"
)

func main() {
	p := fab.Defaults()
	argi := 1
	args := os.Args
	if len(args) > 1 && !strings.Contains(args[1], "=") {
		if err := p.LoadCfg(args[1]); err != nil {
			fmt.Fprintln(os.Stderr, err)
			os.Exit(1)
		}
		argi = 2
	}
	for ; argi < len(args); argi++ {
		eq := strings.IndexByte(args[argi], '=')
		if eq < 0 {
			continue
		}
		p.SetKV(args[argi][:eq], args[argi][eq+1:])
	}
	s := &fab.Sim{P: p, Out: os.Stdout, Errw: os.Stderr}
	s.Run()
}
