# Strict timelike strip through M/T=9: a geometric overlap plateau

This extends the [x≤8 thermal comparator](INDEPENDENT-HARD-VECTOR-NLO-X8-2026-11-26.md)
by a new exact native strip. The full measured bands are excluded
from the strict addition using their original windows.

At the 300 MeV lower temperature cut, the entire new x=8–9 mass
range for temperatures up to the old Born reference root is
[2.400, 3.036] GeV. This interval is strictly inside the first
BESIII window, [2.2324, 3.08] GeV. No new linear interpolant
contributes in the disjoint source over that rectangle. The conditional
root is unchanged to numerical roundoff in the linear comparison.
PCHIP reslopes adjacent intervals, but its root difference is only
1.7×10⁻⁹ GeV at αs=0.3.
This is a geometric masking plateau, not saturation or convergence
of full QCD beyond that strip.

With a 250 MeV lower cut, the new range is only partly masked:
the conditional root at αs=0.3 moves by −0.058 MeV from x≤8.
The absolute shifts from the same Born reference roots are
−2.005 and −0.879 MeV for the 250 and 300 MeV cuts; the removed
fractions of the strict optimistic count at the old roots are
19.13% and 33.34%. The benchmark still uses a fixed coupling,
bounded thermal support, and an optimistic all-convert count,
not a physical reheating floor.
For the 3 sampled native x=8.5 intermediate points, the absolute
strict coefficient residuals are at most 0.064% (PCHIP) and
0.116% (linear). These are not a whole-strip error estimate.

~~~sh
OPENBLAS_NUM_THREADS=1 python -W error partial_hard_vector_nlo_reheat_reference.py \
  EOS.dat KLOExSec.data BESIII-2112.11728.txt \
  results/hard-vector-nlo-through-x9-strips-2026-11-26.txt \
  --intermediate-check results/hard-vector-nlo-through-x9-checks-2026-11-26.txt
OPENBLAS_NUM_THREADS=1 python -W error hard_vector_nlo_overlap_audit_reference.py \
  EOS.dat results/hard-vector-nlo-through-x9-strips-2026-11-26.txt
OPENBLAS_NUM_THREADS=1 python -W error hard_vector_measured_window_saturation_reference.py \
  results/partial-hard-vector-nlo-through-x8-reheat-2026-11-26.json \
  results/partial-hard-vector-nlo-through-x9-reheat-2026-11-26.json
~~~
