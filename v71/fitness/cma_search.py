#!/usr/bin/env python3
"""
cma_search.py — separable CMA-ES over shape genomes (fitness-first soliton search).

Inverted methodology (FUTURE.md 2026-06-26): no Lagrangian/ansatz. Optimize a
free-form gauge-connected blob field to MAXIMIZE fitness
  J = viability * (1 + 0.5*min(nQFI,2)) - drift_penalty
where viability = S_mean * alive_frac (shape retention under kernel evolution).
The kernel is only the behavior oracle; this is seed-gen + scoring + search.

Genome (K blobs): per blob [cx,cy,cz, log-sigma, a0r,a0i,a1r,a1i,a2r,a2i] -> 10*K
real numbers in an UNBOUNDED search space, squashed to physical ranges on decode.
Pure stdlib (no numpy): separable (diagonal) CMA-ES — robust on noisy fitness.

Usage: cma_search.py --bin DIR --sim PATH --work DIR [--K 4 --gen 40 --pop 12 --seed 1]
"""
import os, sys, math, subprocess, argparse, random, time
from concurrent.futures import ThreadPoolExecutor

def decode(x, K, L):
    """unbounded vector -> genome text. tanh-squash to physical ranges."""
    lines = ["omega 1.45"]
    half = 0.45*L
    for b in range(K):
        o = b*10
        cx = half*math.tanh(x[o+0]); cy = half*math.tanh(x[o+1]); cz = half*math.tanh(x[o+2])
        sig = 1.0 + 2.5*0.5*(1+math.tanh(x[o+3]))      # sigma in [1.0, 3.5]
        amp = [math.tanh(x[o+4+i]) for i in range(6)]   # each in [-1,1]
        lines.append("blob %.4f %.4f %.4f %.4f %s" % (cx,cy,cz,sig,
                      " ".join("%.4f"%a for a in amp)))
    return "\n".join(lines)+"\n"

def evaluate(x, K, L, args, gid):
    gpath = os.path.join(args.work, "g_%s.gen"%gid)
    with open(gpath,"w") as f: f.write(decode(x,K,L))
    env = dict(os.environ, N=str(args.N), L=str(L), T=str(args.T),
               SNAP=str(args.snap), ETA=str(args.eta),
               OMP_NUM_THREADS=str(args.omp))
    try:
        out = subprocess.run(["bash", args.evalsh, gpath, gid, args.work, args.bin, args.sim],
                             capture_output=True, text=True, timeout=args.eval_timeout, env=env)
        line = out.stdout.strip().splitlines()[-1] if out.stdout.strip() else ""
    except Exception as e:
        return -1e9, "exc:%s"%e
    if not line.startswith("FITNESS"): return -1e9, line[:80]
    try: J = float(line.split()[1])
    except: J = -1e9
    return J, line

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True); ap.add_argument("--sim", required=True)
    ap.add_argument("--work", required=True); ap.add_argument("--evalsh", required=True)
    ap.add_argument("--K", type=int, default=4); ap.add_argument("--gen", type=int, default=40)
    ap.add_argument("--pop", type=int, default=12); ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--sigma0", type=float, default=1.0)
    ap.add_argument("--N", type=int, default=48); ap.add_argument("--L", type=float, default=12.0)
    ap.add_argument("--T", type=float, default=20.0); ap.add_argument("--snap", type=float, default=0.5)
    ap.add_argument("--eta", type=float, default=0.25); ap.add_argument("--eval_timeout", type=float, default=600)
    ap.add_argument("--workers", type=int, default=8); ap.add_argument("--omp", type=int, default=5)
    ap.add_argument("--warm", type=int, default=0)
    ap.add_argument("--warmgen", type=str, default="")
    args = ap.parse_args()
    os.makedirs(args.work, exist_ok=True)
    rng = random.Random(args.seed)
    D = 10*args.K
    L = args.L
    def atanh(y): y=max(-0.999,min(0.999,y)); return 0.5*math.log((1+y)/(1-y))
    def encode_genome(path, K, L):
        """inverse of decode(): a genome file -> CMA mean vector (pad/truncate to K)."""
        blobs=[]
        for ln in open(path):
            if ln.startswith("blob"):
                p=[float(v) for v in ln.split()[1:11]]
                if len(p)==10: blobs.append(p)
        half=0.45*L; m=[]
        for b in range(K):
            if b < len(blobs):
                cx,cy,cz,sig,*amp = blobs[b]
                m += [atanh(cx/half), atanh(cy/half), atanh(cz/half),
                      atanh((sig-1.0)/1.25 - 1.0)] + [atanh(a) for a in amp]
            else:
                m += [0.0]*10
        return m
    # sep-CMA-ES state
    if args.warmgen:
        mean = encode_genome(args.warmgen, args.K, L)
    elif args.warm:
        # warm start = the known-viable symmetric ball: blob0 carries equal real
        # amplitude in all 3 components (sigma~2.27, amp 0.61), blobs 1..K-1 start
        # at zero amplitude. CMA then climbs QFI by growing satellite structure
        # while the stability gate holds the core conserved.
        mean=[]
        for b in range(args.K):
            if b==0:
                mean += [0,0,0, atanh((2.27-1.0)/2.5*2-1),
                         atanh(0.61),0, atanh(0.61),0, atanh(0.61),0]
            else:
                mean += [0,0,0, 0, 0,0, 0,0, 0,0]
    else:
        mean = [rng.gauss(0,0.3) for _ in range(D)]
    C = [1.0]*D                     # diagonal variances
    sigma = args.sigma0
    lam = args.pop; mu = lam//2
    w = [math.log(mu+0.5)-math.log(i+1) for i in range(mu)]
    sw = sum(w); w = [wi/sw for wi in w]
    mueff = 1.0/sum(wi*wi for wi in w)
    cs = (mueff+2)/(D+mueff+5)
    ds = 1 + 2*max(0, math.sqrt((mueff-1)/(D+1))-1) + cs
    cc = (1+1.0/D+mueff/D)/(D+4+2*mueff/D)/ (1)  # approx
    ccov = (1.0/mueff)*2.0/((D+1.3)**2) + (1-1.0/mueff)*min(1, (2*mueff-1)/((D+2)**2+mueff))
    ccov = min(1.0, ccov*0.5)       # damp for diagonal/noisy
    ps = [0.0]*D; pc = [0.0]*D
    chiN = math.sqrt(D)*(1-1.0/(4*D)+1.0/(21*D*D))
    best = (-1e18, None, "")
    logp = os.path.join(args.work,"cma_log.tsv")
    with open(logp,"w") as lf: lf.write("gen\tbestJ\tmeanJ\tsigma\tline\n")
    evalc = 0
    for g in range(args.gen):
        # sample
        zs=[]; xs=[]
        for k in range(lam):
            z=[rng.gauss(0,1) for _ in range(D)]
            x=[mean[i]+sigma*math.sqrt(C[i])*z[i] for i in range(D)]
            zs.append(z); xs.append(x)
        # evaluate population in parallel (subprocesses release the GIL)
        def _ev(k): return evaluate(xs[k], args.K, L, args, "g%03d_%02d"%(g,k))
        with ThreadPoolExecutor(max_workers=args.workers) as ex:
            res=list(ex.map(_ev, range(lam)))
        fs=[r[0] for r in res]; evalc+=lam
        for k in range(lam):
            if fs[k]>best[0]: best=(fs[k], xs[k][:], res[k][1])
        # rank
        order=sorted(range(lam), key=lambda i:-fs[i])
        meanJ=sum(fs)/lam
        # recombination
        oldmean=mean[:]
        mean=[sum(w[m]*xs[order[m]][i] for m in range(mu)) for i in range(D)]
        zmean=[sum(w[m]*zs[order[m]][i] for m in range(mu)) for i in range(D)]
        # step-size path
        ps=[(1-cs)*ps[i]+math.sqrt(cs*(2-cs)*mueff)*zmean[i] for i in range(D)]
        psn=math.sqrt(sum(p*p for p in ps))
        sigma*=math.exp((cs/ds)*(psn/chiN-1))
        sigma=min(sigma, 3.0)
        # covariance path + diagonal update
        hsig=1.0 if psn/math.sqrt(1-(1-cs)**(2*(g+1)))/chiN < 1.4+2.0/(D+1) else 0.0
        pc=[(1-cc)*pc[i]+hsig*math.sqrt(cc*(2-cc)*mueff)*(mean[i]-oldmean[i])/(sigma*math.sqrt(C[i]) if C[i]>0 else 1) for i in range(D)]
        for i in range(D):
            rankone=pc[i]*pc[i]
            ranktmu=sum(w[m]*(zs[order[m]][i]**2) for m in range(mu))
            C[i]=(1-ccov)*C[i]+ccov*(rankone + ranktmu)
            C[i]=max(C[i],1e-6)
        with open(logp,"a") as lf:
            lf.write("%d\t%.5f\t%.5f\t%.4f\t%s\n"%(g,best[0],meanJ,sigma,best[2][:160]))
        # checkpoint best genome
        if best[1] is not None:
            with open(os.path.join(args.work,"best.gen"),"w") as bf: bf.write(decode(best[1],args.K,L))
        print("gen %d  bestJ=%.4f meanJ=%.4f sigma=%.3f evals=%d\n   %s"%(g,best[0],meanJ,sigma,evalc,best[2][:160]), flush=True)
    print("DONE bestJ=%.5f"%best[0], flush=True)
    print("BEST_GENOME:\n"+decode(best[1],args.K,L), flush=True)

if __name__=="__main__":
    main()
