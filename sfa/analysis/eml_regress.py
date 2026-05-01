#!/usr/bin/env python3
"""
eml_regress.py — Symbolic regression for SCP simulation data.

Brute-force search over RPN (Reverse Polish Notation) programs of growing
length K, matching the approach of the EML paper's RecognizeFunction.
Each program is a stack-machine sequence: constants/variables push to stack,
unary ops pop-apply-push, binary ops pop-pop-apply-push.

Data sources (SCP analysis tools):
  - cluster_profile JSON  ->  radial profile P(r), rho(r)
  - analyze_sfa JSON      ->  time series E(t), P_int(t)
  - diag TSV              ->  any two columns

Usage:
  python3 eml_regress.py profile CB15_profile.json
  python3 eml_regress.py profile CB15_profile.json --field rho
  python3 eml_regress.py timeseries analysis.json Ep
  python3 eml_regress.py tsv diag.tsv t E_pot
  python3 eml_regress.py tsv diag.tsv t E_pot --maxK 9

Requires: numpy.
"""

import argparse
import json
import math
import sys
import time
import itertools

import numpy as np

# ---------------------------------------------------------------------------
# Language definition
# ---------------------------------------------------------------------------
# Symbols grouped by arity (stack effect)
TERMS = [  # arity 0: push one value
    ("1",   0, lambda x: np.ones_like(x)),
    ("2",   0, lambda x: np.full_like(x, 2.0)),
    ("e",   0, lambda x: np.full_like(x, math.e)),
    ("x",   0, lambda x: x.copy()),
]

UNARY = [  # arity 1: pop one, push one
    ("exp",  1, lambda a: np.exp(np.clip(a, -500, 500))),
    ("ln",   1, lambda a: np.log(np.where(a > 0, a, np.nan))),
    ("sqrt", 1, lambda a: np.sqrt(np.where(a >= 0, a, np.nan))),
    ("sqr",  1, lambda a: a * a),
    ("neg",  1, lambda a: -a),
    ("inv",  1, lambda a: np.where(np.abs(a) > 1e-300, 1.0 / a, np.nan)),
]

BINARY = [  # arity 2: pop two, push one
    ("+",  2, lambda a, b: a + b),
    ("-",  2, lambda a, b: a - b),
    ("*",  2, lambda a, b: a * b),
    ("/",  2, lambda a, b: np.where(np.abs(b) > 1e-300, a / b, np.nan)),
    ("^",  2, lambda a, b: np.where(
        (a > 0) & np.isfinite(b) & (np.abs(b * np.log(np.maximum(a, 1e-300))) < 500),
        np.exp(b * np.log(np.maximum(a, 1e-300))), np.nan)),
]

ALL_SYMS = TERMS + UNARY + BINARY


def sym_name(s):  return s[0]
def sym_ar(s):    return s[1]


# ---------------------------------------------------------------------------
# RPN validation and evaluation
# ---------------------------------------------------------------------------
def valid_rpn(arities):
    """Check if a sequence of arities forms a valid RPN program (one result on stack)."""
    stack = 0
    for a in arities:
        if a == 0:
            stack += 1
        else:
            stack -= a
            if stack < 0:
                return False
            stack += 1
    return stack == 1


def eval_rpn(program, x):
    """Evaluate an RPN program on numpy array x. Returns result or None on failure."""
    stack = []
    with np.errstate(all="ignore"):
        for sym in program:
            name, arity, fn = sym
            if arity == 0:
                stack.append(fn(x))
            elif arity == 1:
                if len(stack) < 1:
                    return None
                a = stack.pop()
                stack.append(fn(a))
            elif arity == 2:
                if len(stack) < 2:
                    return None
                b = stack.pop()
                a = stack.pop()
                stack.append(fn(a, b))
    if len(stack) != 1:
        return None
    result = stack[0]
    if not np.all(np.isfinite(result)):
        return None
    if np.any(np.abs(result) > 1e15):
        return None
    return result


def rpn_to_str(program):
    """Convert RPN program to infix string."""
    stack = []
    for sym in program:
        name, arity, fn = sym
        if arity == 0:
            stack.append(name)
        elif arity == 1:
            a = stack.pop()
            if name == "neg":
                stack.append(f"(-{a})")
            elif name == "inv":
                stack.append(f"(1/{a})")
            elif name == "sqr":
                stack.append(f"({a})^2")
            else:
                stack.append(f"{name}({a})")
        elif arity == 2:
            b = stack.pop()
            a = stack.pop()
            if name == "^":
                stack.append(f"({a})^({b})")
            else:
                stack.append(f"({a} {name} {b})")
    return stack[0] if stack else "?"


def rpn_to_tree(program):
    """Convert RPN program to nested tuple tree (for EML compilation)."""
    stack = []
    for sym in program:
        name, arity, fn = sym
        if arity == 0:
            if name == "x":
                stack.append(("var",))
            else:
                val = {"1": 1.0, "2": 2.0, "e": math.e}[name]
                stack.append(("const", val, name))
        elif arity == 1:
            a = stack.pop()
            stack.append(("un", name, a))
        elif arity == 2:
            b = stack.pop()
            a = stack.pop()
            stack.append(("bin", name, a, b))
    return stack[0] if stack else None


# ---------------------------------------------------------------------------
# Enumeration engine
# ---------------------------------------------------------------------------
def enumerate_rpn(max_K, x_data, y_data, top_k=20):
    """Enumerate RPN programs of length 1..max_K, rank by fit to data.

    Uses numerical deduplication at test points to skip equivalent programs.
    Fits y = a * f(x) + b for each candidate (linear regression on top).
    """
    # Signature test points (distinct from data, irrational)
    sig_x = np.array([0.2713, 0.8517, 1.7321, 3.1416, 6.2832])
    y_var = np.var(y_data)
    if y_var < 1e-30:
        y_var = 1.0

    seen_sig = set()
    results = []  # (nmse, formula_str, tree)
    n_tested = 0
    n_valid = 0

    nsyms = len(ALL_SYMS)
    arity_list = [sym_ar(s) for s in ALL_SYMS]

    for K in range(1, max_K + 1):
        t0 = time.time()
        n_at_K = 0

        # Enumerate all length-K programs as index tuples
        for indices in itertools.product(range(nsyms), repeat=K):
            arities = [arity_list[i] for i in indices]
            if not valid_rpn(arities):
                continue
            n_valid += 1

            program = [ALL_SYMS[i] for i in indices]

            # Dedup via signature
            sig_result = eval_rpn(program, sig_x)
            if sig_result is None:
                continue
            sig_key = tuple(np.round(sig_result * 1e9).astype(np.int64))
            if sig_key in seen_sig:
                continue
            seen_sig.add(sig_key)
            n_tested += 1

            # Evaluate on actual data
            yp = eval_rpn(program, x_data)
            if yp is None:
                continue

            fstr = rpn_to_str(program)
            tree = rpn_to_tree(program)

            # Raw fit
            mse = np.mean((yp - y_data) ** 2) / y_var
            if mse < 0.5:
                results.append((mse, fstr, tree))
                n_at_K += 1

            # Scaled: y = a * f(x)
            ff = np.dot(yp, yp)
            if ff > 1e-30:
                a = np.dot(y_data, yp) / ff
                if np.isfinite(a) and 1e-10 < abs(a) < 1e10:
                    mse_s = np.mean((a * yp - y_data) ** 2) / y_var
                    if mse_s < 0.5:
                        results.append((mse_s, f"{a:.6g} * {fstr}", tree))
                        n_at_K += 1

            # Affine: y = a * f(x) + b
            n = len(x_data)
            sf = np.sum(yp); sy = np.sum(y_data); sfy = np.dot(yp, y_data)
            det = n * ff - sf * sf
            if abs(det) > 1e-30:
                a2 = (n * sfy - sf * sy) / det
                b2 = (sy - a2 * sf) / n
                if np.isfinite(a2) and np.isfinite(b2) and abs(a2) < 1e10:
                    mse_a = np.mean((a2 * yp + b2 - y_data) ** 2) / y_var
                    if mse_a < 0.5:
                        if b2 >= 0:
                            fs = f"{a2:.6g} * {fstr} + {b2:.6g}"
                        else:
                            fs = f"{a2:.6g} * {fstr} - {abs(b2):.6g}"
                        results.append((mse_a, fs, tree))
                        n_at_K += 1

        dt = time.time() - t0
        print(f"  K={K}: {n_tested} unique tested, {n_at_K} hits ({dt:.1f}s)")

        # Early stop if we already have very good fits
        results.sort(key=lambda r: r[0])
        if len(results) > 0 and results[0][0] < 1e-6:
            print(f"  [exact match found at K={K}, stopping]")
            break

    # Dedup by string, keep top_k
    results.sort(key=lambda r: r[0])
    seen_s = set()
    out = []
    for mse, s, t in results:
        if s not in seen_s:
            seen_s.add(s)
            out.append((mse, s, t))
            if len(out) >= top_k:
                break
    return out, n_tested


# ---------------------------------------------------------------------------
# EML compilation
# ---------------------------------------------------------------------------
def _eml(a, b):      return f"eml({a},{b})"
def _eml_exp(z):     return _eml(z, "1")
def _eml_log(z):     return _eml("1", _eml_exp(_eml("1", z)))
def _eml_sub(a, b):  return _eml(_eml_log(a), _eml_exp(b))
def _eml_neg(z):     return _eml_sub(_eml_log("1"), z)
def _eml_add(a, b):  return _eml_sub(a, _eml_neg(b))
def _eml_inv(z):     return _eml_exp(_eml_neg(_eml_log(z)))
def _eml_mul(a, b):  return _eml_exp(_eml_add(_eml_log(a), _eml_log(b)))
def _eml_div(a, b):  return _eml_mul(a, _eml_inv(b))
def _eml_pow(a, b):  return _eml_exp(_eml_mul(b, _eml_log(a)))
def _eml_sqrt(z):    return _eml_exp(_eml_mul(_eml_inv(_eml_add("1","1")), _eml_log(z)))


def to_eml(e):
    if e is None: return "?"
    if e[0] == "const":
        v = e[1]
        if v == 1.0:    return "1"
        if v == math.e:  return _eml_exp("1")
        if v == 2.0:    return _eml_add("1", "1")
        return f"C({v})"
    if e[0] == "var": return "x"
    if e[0] == "un":
        c = to_eml(e[2])
        m = {"exp": _eml_exp, "ln": _eml_log, "sqrt": _eml_sqrt,
             "sqr": lambda z: _eml_mul(z, z), "neg": _eml_neg, "inv": _eml_inv}
        return m[e[1]](c)
    if e[0] == "bin":
        l, r = to_eml(e[2]), to_eml(e[3])
        m = {"+": _eml_add, "-": _eml_sub, "*": _eml_mul, "/": _eml_div, "^": _eml_pow}
        return m[e[1]](l, r)
    return "?"


# ---------------------------------------------------------------------------
# Data loaders
# ---------------------------------------------------------------------------
def load_profile(path, field="P_abs", cluster_idx=0, frame_idx=0):
    with open(path) as f:
        data = json.load(f)
    frames = data.get("frames", [data])
    cluster = frames[frame_idx]["clusters"][cluster_idx]
    prof = cluster["radial_profile"]
    r = np.array([b["r"] for b in prof], dtype=np.float64)
    y = np.array([b[field] for b in prof], dtype=np.float64)
    mask = np.isfinite(y)
    return r[mask], y[mask], data.get("physics", {})


def load_timeseries(path, field):
    with open(path) as f:
        data = json.load(f)
    fmap = {"E":"E","Ep":"Ep","E_pot":"Ep","Pint":"Pint","P_int":"Pint",
            "phi_max":"phi_max","theta_rms":"theta_rms","aspect":"aspect","S":"S"}
    key = fmap.get(field, field)
    frames = data["frames"]
    t = np.array([f["t"] for f in frames], dtype=np.float64)
    y = np.array([f[key] for f in frames], dtype=np.float64)
    mask = np.concatenate([[True], np.diff(t) > 0])
    return t[mask], y[mask], data


def load_tsv(path, xcol, ycol):
    with open(path) as f:
        hdr = f.readline().strip().split("\t")
        rows = [l.strip().split("\t") for l in f if l.strip()]
    xi, yi = hdr.index(xcol), hdr.index(ycol)
    xv = np.array([float(r[xi]) for r in rows if len(r) > max(xi, yi)], dtype=np.float64)
    yv = np.array([float(r[yi]) for r in rows if len(r) > max(xi, yi)], dtype=np.float64)
    mask = np.isfinite(xv) & np.isfinite(yv)
    return xv[mask], yv[mask], {"columns": hdr}


# ---------------------------------------------------------------------------
# Display
# ---------------------------------------------------------------------------
def show(results, x, y, eml_out=True, n=10):
    ystd = np.std(y) if np.std(y) > 0 else 1.0
    print(f"\n{'='*72}")
    print(f"  TOP {min(len(results), n)} DISCOVERED FORMULAS")
    print(f"  Data: {len(x)} points, x in [{x.min():.4g}, {x.max():.4g}]")
    print(f"  y range: [{y.min():.6g}, {y.max():.6g}]")
    print(f"{'='*72}\n")

    for i, (nmse, fstr, tree) in enumerate(results[:n]):
        rmse = math.sqrt(max(nmse, 0)) * ystd
        pct = math.sqrt(max(nmse, 0)) * 100
        print(f"  #{i+1}  NMSE={nmse:.6g}  RMSE={rmse:.6g}  ({pct:.1f}% rel)")
        print(f"      {fstr}")
        if eml_out and tree is not None:
            try:
                em = to_eml(tree)
                if len(em) > 100:
                    em = em[:100] + "..."
                print(f"      EML: {em}")
            except Exception:
                pass
        print()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="Symbolic regression for SCP data (EML framework)")
    sp = ap.add_subparsers(dest="mode")

    p1 = sp.add_parser("profile", help="Fit radial profile")
    p1.add_argument("path")
    p1.add_argument("--field", default="P_abs")
    p1.add_argument("--cluster", type=int, default=0)
    p1.add_argument("--frame", type=int, default=0)

    p2 = sp.add_parser("timeseries", help="Fit time series")
    p2.add_argument("path"); p2.add_argument("field")

    p3 = sp.add_parser("tsv", help="Fit TSV columns")
    p3.add_argument("path"); p3.add_argument("x_col"); p3.add_argument("y_col")

    for p in [p1, p2, p3]:
        p.add_argument("--maxK", type=int, default=7,
                       help="Max RPN program length (default: 7)")
        p.add_argument("--top", type=int, default=10)
        p.add_argument("--no-eml", action="store_true")

    args = ap.parse_args()
    if not args.mode:
        ap.print_help(); sys.exit(1)

    if args.mode == "profile":
        x, y, meta = load_profile(args.path, args.field, args.cluster, args.frame)
        print(f"Profile: {len(x)} bins, field={args.field}")
        if meta: print(f"  Physics: {meta}")
    elif args.mode == "timeseries":
        x, y, meta = load_timeseries(args.path, args.field)
        print(f"Time series: {len(x)} frames, field={args.field}")
    else:
        x, y, meta = load_tsv(args.path, args.x_col, args.y_col)
        print(f"TSV: {len(x)} rows, {args.x_col} -> {args.y_col}")

    print(f"  x: [{x.min():.4g}, {x.max():.4g}]  y: [{y.min():.6g}, {y.max():.6g}]")
    print(f"\nSearching RPN programs K=1..{args.maxK}...")

    results, n_tested = enumerate_rpn(args.maxK, x, y, top_k=args.top)

    if not results:
        print(f"\nNo formulas found in {n_tested} candidates. Try --maxK {args.maxK + 2}.")
        sys.exit(1)

    show(results, x, y, eml_out=not args.no_eml, n=args.top)


if __name__ == "__main__":
    main()
