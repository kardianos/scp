// run_lissajous.go — 5×5×5 Lissajous phase sweep runner
//
// go run run_lissajous.go
package main

import (
	"fmt"
	"math"
	"os"
	"os/exec"
	"path/filepath"
	"regexp"
	"strings"
)

const (
	baseSeed = "/home/d/code/scp/v52/chirality_test/base_seed.sfa"
	genTool  = "/home/d/code/scp/v52/chirality_test/gen_lissajous_seed"
	simTool  = "/home/d/code/scp/bin/scp_sim"
	outDir   = "/space/scp/v52/lissajous"
)

var phases = []float64{0, math.Pi / 3, 2 * math.Pi / 3, math.Pi, 4 * math.Pi / 3}
var labels = []string{"0", "p3", "2p3", "pi", "4p3"}

func extract(output, pattern string) string {
	re := regexp.MustCompile(pattern)
	m := re.FindStringSubmatch(output)
	if len(m) > 1 {
		return m[1]
	}
	return "?"
}

func main() {
	os.MkdirAll(outDir, 0755)

	resFile, _ := os.Create(filepath.Join(outDir, "results.tsv"))
	fmt.Fprintln(resFile, "d0\td1\td2\tl0\tl1\tl2\tE_init\tE_final\tphi_i\tphi_f\tPint_i\tPint_f\ttrms_f\tstatus")

	total := len(phases) * len(phases) * len(phases)
	run := 0

	for i, d0 := range phases {
		for j, d1 := range phases {
			for k, d2 := range phases {
				run++
				name := fmt.Sprintf("%s_%s_%s", labels[i], labels[j], labels[k])
				sfaOut := filepath.Join(outDir, name+".sfa")
				diagOut := filepath.Join(outDir, name+"_diag.tsv")
				seedFile := filepath.Join(outDir, name+"_seed.sfa")
				cfgFile := filepath.Join(outDir, name+".cfg")

				// Skip if diag already exists
				if _, err := os.Stat(diagOut); err == nil {
					fmt.Printf("[%d/%d] %s — exists, skipping\n", run, total, name)
					continue
				}

				fmt.Printf("[%d/%d] %s (%.2f,%.2f,%.2f) ... ", run, total, name, d0, d1, d2)

				// Generate seed
				genOut, err := exec.Command(genTool, baseSeed, seedFile,
					fmt.Sprintf("%.4f", d0), fmt.Sprintf("%.4f", d1), fmt.Sprintf("%.4f", d2)).CombinedOutput()
				if err != nil {
					fmt.Printf("SEED FAIL: %v\n", err)
					continue
				}
				_ = genOut

				// Write config
				cfg := fmt.Sprintf(`N=64
L=8.3
T=100
dt_factor=0.025
m=1.5
m_theta=0
eta=0.5
mu=-41.345
kappa=50
bc_type=2
damp_width=0
damp_rate=0
precision=f16
init=sfa
init_sfa=%s
output=%s
diag_file=%s
snap_dt=25
diag_dt=2
`, seedFile, sfaOut, diagOut)
				os.WriteFile(cfgFile, []byte(cfg), 0644)

				// Run simulation
				cmd := exec.Command(simTool, cfgFile)
				cmd.Env = append(os.Environ(), "OMP_NUM_THREADS=8")
				out, err := cmd.CombinedOutput()
				result := string(out)

				// Extract metrics from INIT line and COMPLETE block
				eInit := extract(result, `INIT: E_total=([0-9.e+-]+)`)
				phiInit := extract(result, `INIT:.*phi_max=([0-9.]+)`)
				pintInit := extract(result, `INIT:.*P_int=([0-9.e+-]+)`)

				// Find the final summary lines (after === COMPLETE ===)
				completeIdx := strings.LastIndex(result, "=== COMPLETE")
				tail := ""
				if completeIdx >= 0 {
					tail = result[completeIdx:]
				}
				eFinal := extract(tail, `E_total=([0-9.e+-]+)`)
				phiFinal := extract(tail, `phi_max=([0-9.]+)`)
				pintFinal := extract(tail, `P_int=([0-9.e+-]+)`)
				trmsFinal := extract(tail, `theta_rms=([0-9.e+-]+)`)

				status := "dissolved"
				if phiFinal != "?" {
					var pf float64
					fmt.Sscanf(phiFinal, "%f", &pf)
					if pf > 0.3 {
						status = "alive"
					}
				}

				fmt.Printf("%s (φ=%s P=%s)\n", status, phiFinal, pintFinal)

				fmt.Fprintf(resFile, "%.4f\t%.4f\t%.4f\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
					d0, d1, d2, labels[i], labels[j], labels[k],
					eInit, eFinal, phiInit, phiFinal, pintInit, pintFinal, trmsFinal, status)
				resFile.Sync()

				// Clean up seed and config
				os.Remove(seedFile)
				os.Remove(cfgFile)
			}
		}
	}

	resFile.Close()

	// Summary
	data, _ := os.ReadFile(filepath.Join(outDir, "results.tsv"))
	alive := strings.Count(string(data), "alive")
	dissolved := strings.Count(string(data), "dissolved")
	fmt.Printf("\n=== COMPLETE ===\n")
	fmt.Printf("Survivors: %d / %d\n", alive, total)
	fmt.Printf("Dissolved: %d / %d\n", dissolved, total)
	fmt.Printf("Results: %s/results.tsv\n", outDir)
}
