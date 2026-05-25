#!/usr/bin/env python3
"""
Small helper to extract a PythonReportSlice from a JSON report (or synthetic slice)
and print it in a simple key=value format that the Lean 7D_Algebra certificates
can parse easily.

Usage:
  python3 extract_report_slice.py <json_path>

Expected JSON shape (minimal):
{
  "fAmp": 0.55,
  "description": "...",
  "L_cross_observed": 13.333,
  "LF_cross_observed": 0.0,
  "L_degradation_min_pct": 0.12,
  "LF_degradation_max_pct": 0.04
}
"""
import json
import sys

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 extract_report_slice.py <json_path>", file=sys.stderr)
        sys.exit(1)

    path = sys.argv[1]
    with open(path) as f:
        data = json.load(f)

    # Support both the synthetic slice format and real report rows if needed
    fAmp = data.get("fAmp", data.get("lambda_nl", 0.0))
    desc = data.get("description", data.get("source", "Python report slice"))
    L_cross = float(data.get("L_cross_observed", data.get("cross_term_pct", 0.0)))
    LF_cross = float(data.get("LF_cross_observed", 0.0))
    L_deg_min = float(data.get("L_degradation_min_pct", 0.12))
    LF_deg_max = float(data.get("LF_degradation_max_pct", 0.04))

    print(f"fAmp={fAmp}")
    print(f"description={desc}")
    print(f"L_cross_observed={L_cross}")
    print(f"LF_cross_observed={LF_cross}")
    print(f"L_degradation_min_pct={L_deg_min}")
    print(f"LF_degradation_max_pct={LF_deg_max}")

if __name__ == "__main__":
    main()
