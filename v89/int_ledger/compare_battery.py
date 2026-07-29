#!/usr/bin/env python3
"""Run battery experiments under cellfab (FP64) and cellfabi ledger modes.

usage:
  python3 compare_battery.py --laws ../battery/laws_V2g.cfg \
      [--modes 0,1,2,3] [--only e1_conserve] [--jobs 4]

Writes logs under int_ledger/runs/<tag>/ and a SUMMARY.tsv.
Does NOT modify production battery/runs.
"""
import argparse
import concurrent.futures as cf
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
V89 = os.path.dirname(HERE)
BAT = os.path.join(V89, "battery")
EXPS = [
    "e1_conserve", "e2_pulse", "e3a_blob", "e3b_blob_tilt", "e4_curve",
    "e5_bell", "e6_pairs", "e7_tune", "e8_comma", "e9_fifth", "d1_slit",
    "t1_tonomura", "q2_eraser", "t4_hom", "qt_lo", "qt_hi",
    "p1_beam", "g1_footprint", "g3_shadow", "g4_throughput",
]


def build(bin_name, src_name):
    cmd = ["gcc", "-O2", "-march=native", "-fopenmp",
           "-o", os.path.join(V89, bin_name),
           os.path.join(V89, src_name), "-lm"]
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=V89)
    if r.returncode != 0:
        sys.exit(f"build {bin_name} failed:\n{r.stderr}")


def run_one(binary, laws_path, name, outdir, extra_lines, threads):
    os.makedirs(os.path.join(outdir, "cfg"), exist_ok=True)
    cfg = os.path.join(outdir, "cfg", name + ".cfg")
    log = os.path.join(outdir, name + ".log")
    with open(laws_path) as f:
        laws = f.read()
    with open(os.path.join(BAT, "apparatus", name + ".cfg")) as f:
        app = f.read()
    with open(cfg, "w") as f:
        f.write(laws)
        f.write("\n# --- apparatus ---\n")
        f.write(app)
        for line in extra_lines:
            f.write(line.rstrip() + "\n")
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    with open(log, "w") as f:
        r = subprocess.run([binary, cfg], stdout=f, stderr=subprocess.STDOUT,
                           env=env)
    return name, r.returncode, log


def parse_result(log_path):
    text = open(log_path, errors="replace").read()
    out = {"rc_ok": True, "rel_drift": None, "ledger_sum_err": None,
           "ledger_max_diff": None, "result_lines": []}
    m = re.search(r"# RESULT conservation .*rel_drift=([-\d.eE+]+)", text)
    if m:
        out["rel_drift"] = float(m.group(1))
    # prefer RESULT LEDGER; else last # LEDGER diag line
    for m in re.finditer(
            r"# RESULT LEDGER mode=(\d+) E0_int=([-\d]+) sum_int=([-\d]+) "
            r"max_sum_err=([-\d]+) max\|E-iEu\|=([-\d.eE+]+)", text):
        out["ledger_mode"] = int(m.group(1))
        out["ledger_sum_err"] = int(m.group(4))
        out["ledger_max_diff"] = float(m.group(5))
    if out["ledger_sum_err"] is None:
        for m in re.finditer(
                r"# LEDGER t=[\d.]+ mode=(\d+) sum_err=([-\d]+) "
                r"max\|E-iEu\|=([-\d.eE+]+) run_max_diff=([-\d.eE+]+) "
                r"run_max_sum_err=([-\d]+)", text):
            out["ledger_mode"] = int(m.group(1))
            out["ledger_sum_err"] = int(m.group(5))  # run max
            out["ledger_max_diff"] = float(m.group(4))
    out["result_lines"] = [ln for ln in text.splitlines()
                           if ln.startswith("# RESULT")]
    if "FATAL" in text:
        out["rc_ok"] = False
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--laws", default=os.path.join(BAT, "laws_V2g.cfg"))
    ap.add_argument("--modes", default="0,1,2,3",
                    help="cellfabi ledger_mode list; also always runs cellfab")
    ap.add_argument("--only", nargs="*", default=None)
    ap.add_argument("--jobs", type=int, default=2)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--skip-build", action="store_true")
    ap.add_argument("--skip-fp64", action="store_true",
                    help="do not re-run production cellfab")
    a = ap.parse_args()

    exps = a.only if a.only else EXPS
    modes = [int(x) for x in a.modes.split(",") if x.strip() != ""]

    if not a.skip_build:
        print("building cellfab + cellfabi...")
        build("cellfab", "cellfab.c")
        build("cellfabi", "cellfabi.c")

    bin_fab = os.path.join(V89, "cellfab")
    bin_i = os.path.join(V89, "cellfabi")
    root_out = os.path.join(HERE, "runs")
    os.makedirs(root_out, exist_ok=True)

    jobs = []
    # tag -> (binary, extra cfg lines)
    variants = {}
    if not a.skip_fp64:
        variants["fp64"] = (bin_fab, [])
    for m in modes:
        variants[f"i_m{m}"] = (bin_i, [f"ledger_mode={m}", "ledger_u=1e-12"])

    # run sequentially per variant to avoid OMP oversubscription; parallel exps inside
    summary_rows = []
    for tag, (binary, extra) in variants.items():
        outdir = os.path.join(root_out, tag)
        os.makedirs(outdir, exist_ok=True)
        print(f"=== variant {tag} ({len(exps)} exps) ===")
        thr = max(1, a.threads)
        # cap concurrent so jobs*threads ~ cores
        npar = max(1, min(a.jobs, max(1, 8 // thr)))
        results = {}
        with cf.ThreadPoolExecutor(max_workers=npar) as ex:
            futs = {
                ex.submit(run_one, binary, a.laws, name, outdir, extra, thr): name
                for name in exps
            }
            for fut in cf.as_completed(futs):
                name, rc, log = fut.result()
                pr = parse_result(log)
                pr["rc"] = rc
                results[name] = pr
                flag = "OK" if rc == 0 and pr["rc_ok"] else "FAIL"
                print(f"  {name}: {flag} drift={pr['rel_drift']} "
                      f"led_err={pr['ledger_sum_err']}")
        for name in exps:
            pr = results[name]
            summary_rows.append({
                "variant": tag,
                "exp": name,
                "rc": pr["rc"],
                "rel_drift": pr["rel_drift"],
                "ledger_sum_err": pr["ledger_sum_err"],
                "ledger_max_diff": pr["ledger_max_diff"],
                "n_result": len(pr["result_lines"]),
            })

    # write SUMMARY
    sum_path = os.path.join(root_out, "SUMMARY.tsv")
    with open(sum_path, "w") as f:
        f.write("variant\texp\trc\trel_drift\tledger_sum_err\tledger_max_diff\tn_result\n")
        for r in summary_rows:
            f.write(
                f"{r['variant']}\t{r['exp']}\t{r['rc']}\t{r['rel_drift']}\t"
                f"{r['ledger_sum_err']}\t{r['ledger_max_diff']}\t{r['n_result']}\n"
            )
    print(f"wrote {sum_path}")

    # drift comparison vs fp64 when present
    by = {}
    for r in summary_rows:
        by.setdefault(r["exp"], {})[r["variant"]] = r
    if "fp64" in variants:
        print("\n=== rel_drift vs fp64 ===")
        for name in exps:
            fp = by[name].get("fp64", {})
            base = fp.get("rel_drift")
            for tag in variants:
                if tag == "fp64":
                    continue
                d = by[name].get(tag, {}).get("rel_drift")
                print(f"  {name:16s} {tag:6s} drift={d}  fp64={base}")


if __name__ == "__main__":
    main()
