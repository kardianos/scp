#!/bin/bash
# v97 H0 honesty suite — reproducible regeneration (was manual for H0.md).
# Classes/streams per H0.md; `exact` class added by PLAN.md step 1'.
set -uo pipefail
cd "$(dirname "$0")/.."
D=./dipole
echo "v97 H0 honesty suite"
echo
$D -h0 sealed   ../v93/runs/faceC/c_deep2.fcs; echo
$D -h0 sealed   /space/scp/v93/fcs/hold_S00_uni_quiet.fcs; echo
$D -h0 leak     /space/scp/v93/fcs/rad_tax_S05.fcs; echo
$D -h0 leak     ../v91/runs/pad/m7_g4_k005.fcs; echo
$D -h0 ring     ../v93/runs/faceA/ringA_strang_live.fcs; echo
$D -h0 ring     ../v93/runs/quench/rv6_uni_m2.fcs; echo
$D -h0 bathring ../v93/runs/quench/rb6_uni_m2.fcs; echo
$D -h0 current  ../v93/runs/prog/e3b_lin.fcs; echo
$D -h0 exact    runs/leads/p2_space.fcs; echo
$D -h0 exact    runs/leads/p2_hot.fcs; echo
