#!/usr/bin/env python3
"""
Deterministic Round-2 tables: evaluate L_fit at known-good / optimized seeds.
Writes round2_* without requiring interactive multi-hour search.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import congruence_score_r2 as R

OUT = Path(__file__).resolve().parent / "results"
OUT.mkdir(exist_ok=True)

OCC = R.Occam()


def fit_at(obj, seeds, bounds, step0):
    return R.multi_start(obj, seeds, bounds, step0)


def eval_map(rm: R.RayMap) -> dict:
    b, data_d, data_a = rm.b, rm.delay, rm.alpha
    Mled = rm.M_ledger
    scores = []

    def add(name, Lfit, pred_M, n_sec, link, softE, post, nogo, params, notes=""):
        scores.append(
            R.score_model(
                name,
                Lfit,
                pred_M,
                Mled,
                n_sectors=n_sec,
                has_link=link,
                softE=softE,
                postulated=post,
                long_range_local_claim=nogo,
                params=params,
                occ=OCC,
                notes=notes,
            )
        )

    # local GRIN
    def obj_lo(p):
        d, a, _ = R.local_grin_rays(p[0], p[1], b)
        return R.l_fit_rays(d, a, data_d, data_a)

    p_lo, L_lo = fit_at(
        obj_lo,
        [[0.35, 1.2], [0.15, 1.2], [0.5, 1.0], [0.2, 1.5], [0.4, 0.9], [0.25, 1.0]],
        [(1e-3, 0.94), (0.2, 4.0)],
        [0.04, 0.08],
    )
    _, _, M_lo = R.local_grin_rays(p_lo[0], p_lo[1], b)
    add(
        "monist_local_GRIN",
        L_lo,
        M_lo,
        1,
        True,
        False,
        False,
        False,
        {"A": p_lo[0], "sigma": p_lo[1]},
        "B-style n=ρ0/ρ_free",
    )

    # monist kernel
    def obj_mk(p):
        d, a, _ = R.monist_kernel_rays(p[0], p[1], p[2], b)
        return R.l_fit_rays(d, a, data_d, data_a)

    p_mk, L_mk = fit_at(
        obj_mk,
        [
            [Mled, 0.08, 0.15],
            [Mled, 0.08, 0.5],
            [Mled, 0.12, 0.85],
            [Mled, 0.05, 0.3],
            [Mled, 0.1, 1.0],
            [Mled * 1.1, 0.07, 0.2],
        ],
        [(0.1, 20.0), (1e-3, 2.0), (0.05, 4.0)],
        [0.15, 0.015, 0.08],
    )
    add(
        "monist_kernel_postulated",
        L_mk,
        p_mk[0],
        1,
        True,
        False,
        True,
        False,
        {"M": p_mk[0], "beta": p_mk[1], "eps": p_mk[2]},
        "postulated Class-C kernel",
    )

    if rm.meta.get("dynamical_free_response"):
        add(
            "monist_dynamical_free_response",
            L_mk,
            p_mk[0],
            1,
            True,
            False,
            False,
            False,
            {"M": p_mk[0], "beta": p_mk[1], "eps": p_mk[2]},
            "dynamical — no postulate penalty",
        )

    # dualist plummer
    def obj_pl(p):
        d, a, _ = R.dualist_plummer_rays(p[0], p[1], p[2], b)
        return R.l_fit_rays(d, a, data_d, data_a)

    p_pl, L_pl = fit_at(
        obj_pl,
        [[Mled, 0.08, 0.15], [Mled, 0.08, 0.5], [Mled, 0.12, 1.0], [Mled, 0.1, 0.3]],
        [(0.1, 20.0), (1e-3, 2.0), (0.05, 4.0)],
        [0.15, 0.015, 0.08],
    )
    add(
        "dualist_plummer_isomorphic",
        L_pl,
        p_pl[0],
        2,
        False,
        False,
        False,
        False,
        {"M": p_pl[0], "G": p_pl[1], "a": p_pl[2]},
        "isomorphic control",
    )

    # dualist log NON-iso
    def obj_log(p):
        d, a, _ = R.dualist_log_rays(p[0], p[1], p[2], b)
        return R.l_fit_rays(d, a, data_d, data_a)

    p_log, L_log = fit_at(
        obj_log,
        [
            [Mled, 0.05, 0.5],
            [Mled, 0.1, 1.0],
            [Mled, 0.2, 0.8],
            [Mled, 0.08, 0.3],
            [Mled, 0.3, 2.0],
            [2.0, 0.15, 1.2],
        ],
        [(0.1, 20.0), (1e-4, 2.0), (0.05, 5.0)],
        [0.15, 0.025, 0.12],
    )
    add(
        "dualist_log_2D",
        L_log,
        p_log[0],
        2,
        False,
        False,
        False,
        False,
        {"M": p_log[0], "G": p_log[1], "a": p_log[2]},
        "NON-isomorphic dualist",
    )

    # dualist yukawa
    def obj_yk(p):
        d, a, _ = R.dualist_yukawa_rays(p[0], p[1], p[2], p[3], b)
        return R.l_fit_rays(d, a, data_d, data_a)

    p_yk, L_yk = fit_at(
        obj_yk,
        [
            [Mled, 0.08, 0.0, 0.15],
            [Mled, 0.1, 0.3, 0.5],
            [Mled, 0.12, 0.8, 0.8],
            [Mled, 0.08, 0.1, 0.3],
            [Mled, 0.2, 1.0, 1.0],
            [Mled, 0.05, 0.05, 0.5],
        ],
        [(0.1, 20.0), (1e-4, 2.0), (0.0, 3.0), (0.05, 4.0)],
        [0.15, 0.02, 0.08, 0.08],
    )
    add(
        "dualist_yukawa",
        L_yk,
        p_yk[0],
        2,
        False,
        False,
        False,
        False,
        {"M": p_yk[0], "G": p_yk[1], "mu": p_yk[2], "a": p_yk[3]},
        "NON-iso (μ>0) or near-iso (μ→0)",
    )

    # free GM 1/b
    def obj_1b(p):
        d, a, _ = R.dualist_point_mass_1overb(p[0], b)
        return R.l_fit_rays(d, a, data_d, data_a)

    p_1b, L_1b = fit_at(
        obj_1b,
        [[0.25], [0.5], [0.1], [0.253], [0.08], [1.0], [0.4]],
        [(1e-4, 5.0)],
        [0.04],
    )
    add(
        "dualist_free_GM_1overb",
        L_1b,
        float("nan"),
        2,
        False,
        False,
        False,
        False,
        {"amp_GM": p_1b[0]},
        "free G*M; no ledger",
    )

    # softE
    sc = R.score_model(
        "softE_dualist_log_as_monist",
        L_log,
        p_log[0],
        Mled,
        n_sectors=1,
        has_link=True,
        softE=True,
        postulated=False,
        long_range_local_claim=False,
        params={"M": p_log[0], "G": p_log[1], "a": p_log[2]},
        occ=OCC,
        notes="log dualist labeled monist",
    )
    sc.L_occ = 0.0
    sc.L_M = 0.0
    sc.L_extra = 0.0
    sc.S = L_log + OCC.lambda_softE
    scores.append(sc)

    scores_sorted = sorted(scores, key=lambda s: s.S)
    winner = scores_sorted[0]
    fit_winner = min(scores, key=lambda s: s.L_fit)
    nc = R.nc_checklist(rm, winner, scores)
    return {
        "map_tag": rm.tag,
        "channel": rm.meta.get("channel"),
        "M_ledger": Mled,
        "free_deficit_core": rm.free_deficit_core,
        "score_winner": winner.name,
        "fit_winner": fit_winner.name,
        "scores": [R.asdict(s) for s in scores_sorted],
        "nc_checklist": nc,
    }


def main():
    maps = []
    maps.extend(R.load_b_local())
    maps.extend(R.load_b_kernel())
    maps.extend(R.load_b_dynamical())

    eval_maps = []
    for m in maps:
        if m.tag.startswith("lock_A0.35") or "kernel" in m.tag:
            eval_maps.append(m)
        elif "A0.15" in m.tag:
            eval_maps.append(m)

    report = {
        "round": 2,
        "occam": R.asdict(OCC),
        "maps_found": [
            {
                "tag": m.tag,
                "channel": m.meta.get("channel"),
                "M_ledger": m.M_ledger,
                "free_deficit_core": m.free_deficit_core,
                "dynamical": m.meta.get("dynamical_free_response"),
                "source": m.source,
            }
            for m in maps
        ],
        "evaluations": [],
        "dynamical_free_response_present": any(
            m.meta.get("dynamical_free_response") for m in maps
        ),
    }

    for rm in eval_maps:
        print(f"evaluating {rm.tag} ...")
        ev = eval_map(rm)
        report["evaluations"].append(ev)
        print(
            f"  score_winner={ev['score_winner']} fit_winner={ev['fit_winner']}"
        )
        for s in ev["scores"][:4]:
            print(f"    {s['name']}: L_fit={s['L_fit']:.4e} S={s['S']:.4e}")

    report["summary"] = R.make_summary(report)
    R.save(report)

    # model rank
    with open(OUT / "round2_model_rank.tsv", "w") as f:
        f.write("map\tchannel\tmodel\tL_fit\tS\trank_S\n")
        for ev in report["evaluations"]:
            for i, s in enumerate(ev["scores"]):
                f.write(
                    f"{ev['map_tag']}\t{ev['channel']}\t{s['name']}\t"
                    f"{s['L_fit']:.8e}\t{s['S']:.8e}\t{i+1}\n"
                )

    print("=== SUMMARY ===")
    print(json.dumps(report["summary"], indent=2))
    return report


if __name__ == "__main__":
    main()
