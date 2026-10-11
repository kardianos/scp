#!/bin/bash
# T0 (../PREREG.md §4): exact enumeration against the sampler on small
# periodic graphs, every rule set whose enumeration fits; lattice and
# Green's-function gates on the production lattices.
cd "$(dirname "$0")"
mkdir -p out
T="seeds=32 sweeps=200000 warm=1000"
{
for r in c ac bc abc; do ./chain mode=t0 lat=square L1=2 L2=4 rules=$r $T; done
for r in bc abc; do ./chain mode=t0 lat=square L=4 rules=$r $T; done
for r in none a b ab; do ./chain mode=t0 lat=kagome L1=1 L2=2 rules=$r $T; done
for r in a b ab; do ./chain mode=t0 lat=kagome L=2 rules=$r $T; done
for r in ac bc abc; do ./chain mode=t0 lat=cubic L=2 rules=$r $T; done
for r in b ab; do ./chain mode=t0 lat=pyro L=1 rules=$r $T; done
echo "# not enumerable here (listed, not run): square 4x4 b-off; kagome 2x2 none; cubic 2x2x2 c; pyrochlore b-off"
for spec in "lat=ring L=128 rules=bc" "lat=square L=32 rules=bc" "lat=kagome L=18 rules=b" "lat=cubic L=12 rules=bc" "lat=pyro L=5 rules=b"; do
  ./chain mode=lat $spec | grep -E "^#|LATTICE|GREEN"
done
} > out/t0.out 2> out/t0.err
touch out/t0.done
