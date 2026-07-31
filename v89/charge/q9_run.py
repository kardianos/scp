#!/usr/bin/env python3
"""Q9 scan driver: run the dw scan 2-at-a-time, 4 threads each."""
import glob, os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, "..", "cellfab")
RUNS = os.path.join(HERE, "runs")


def main():
    cfgs = sorted(glob.glob(os.path.join(RUNS, "q9_dw*.cfg")))
    if len(sys.argv) > 1:
        cfgs = [c for c in cfgs if any(a in c for a in sys.argv[1:])]
    env = dict(os.environ, OMP_NUM_THREADS="4")
    procs = []
    for cfg in cfgs:
        log = cfg[:-4] + ".log"
        if os.path.exists(log) and "# RESULT slip" in open(log).read():
            print(f"skip (done): {os.path.basename(cfg)}")
            continue
        while len([p for p in procs if p.poll() is None]) >= 2:
            for p in procs:
                p.poll()
            os.wait()
        fh = open(log, "w")
        p = subprocess.Popen([BIN, cfg], stdout=fh, stderr=subprocess.STDOUT, env=env)
        procs.append(p)
        print(f"launch: {os.path.basename(cfg)} pid={p.pid}")
    for p in procs:
        p.wait()
    print("scan complete")


if __name__ == "__main__":
    main()
