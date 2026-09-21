package main

import (
	"bufio"
	"os"
	"regexp"
	"strconv"
	"strings"
)

// LogResult holds freecell # RESULT lines of interest.
type LogResult struct {
	Drift, Cond, Evap, Rough, Rad, Backs float64
	QcnvD, QcnvF                         float64
	Fd, Fl, Sp, Gm                       float64 // p1x
	EmTag, BlobSpeed, BlobCos            float64
	FrontV                               float64
	HasP1, HasBlob                       bool
}

var (
	reKV   = regexp.MustCompile(`([A-Za-z_][A-Za-z0-9_]*)=([+\-0-9.eE]+)`)
	reNum  = regexp.MustCompile(`([+\-]?[0-9]*\.?[0-9]+(?:[eE][+\-]?[0-9]+)?)`)
)

func parseLog(path string) (LogResult, error) {
	var r LogResult
	f, err := os.Open(path)
	if err != nil {
		return r, err
	}
	defer f.Close()
	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 0, 64*1024), 4*1024*1024)
	for sc.Scan() {
		line := sc.Text()
		if !strings.Contains(line, "RESULT") {
			continue
		}
		kv := map[string]float64{}
		for _, m := range reKV.FindAllStringSubmatch(line, -1) {
			v, err := strconv.ParseFloat(m[2], 64)
			if err == nil {
				kv[m[1]] = v
			}
		}
		switch {
		case strings.Contains(line, "drift_rel"):
			// "# RESULT drift_rel -1.23e-15"
			nums := reNum.FindAllString(line, -1)
			if len(nums) > 0 {
				// last number is the value (first may be from RESULT alone — none)
				r.Drift, _ = strconv.ParseFloat(nums[len(nums)-1], 64)
			}
		case strings.Contains(line, "conv rough="):
			r.Rough = kv["rough"]
			r.Cond = kv["cond"]
			r.Evap = kv["evap"]
			r.Backs = kv["backs"]
			r.Rad = kv["rad"]
		case strings.Contains(line, "credit"):
			r.QcnvD = kv["qcnvD"]
			r.QcnvF = kv["qcnvF"]
		case strings.Contains(line, "p1x"):
			r.Sp = kv["sp"]
			r.Fl = kv["fl"]
			r.Fd = kv["fd"]
			r.Gm = kv["gm"]
			r.HasP1 = true
		case strings.Contains(line, "blob Em_tag"):
			r.EmTag = kv["Em_tag"]
			r.HasBlob = true
		case strings.Contains(line, "blob_drift"):
			r.BlobSpeed = kv["speed"]
			r.BlobCos = kv["cos_to_kdir"]
		case strings.Contains(line, "front_speed"):
			r.FrontV = kv["v"]
		}
	}
	return r, sc.Err()
}
