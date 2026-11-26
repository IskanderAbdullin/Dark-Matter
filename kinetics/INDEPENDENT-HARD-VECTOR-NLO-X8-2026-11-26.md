# Strict timelike continuation through M/T=8 and two disjoint measured windows

This extends the [bounded x=7 result](INDEPENDENT-HARD-VECTOR-NLO-X7-2026-11-26.md)
by one complete native strip. It uses the same fixed-coupling vacuum
data patches, strict thermal comparator, and optimistic conversion
target; it is not a full finite-temperature QCD calculation. No
rates outside the provided mass and boost grids are included in
the strict correction.

The upper edge first intersects the lower BES window at 279.05 MeV
and the higher one at 425 MeV. The independent physical-mass
integration now separates each window. The second-window contribution
is exactly zero just below 425 MeV and positive just above it,
while the first window is already active. The summed partition
matches the production-code subtraction to about 1.1×10⁻⁷ relative
to the local unsplit rate. Only the first window intersects the strict
correction at the reported roots.

At illustrative αs=0.3, the conditional roots are 323.229 and
336.423 MeV for hard lower temperature cuts 250 and 300 MeV.
Shifts from the common measured + disjoint Born baselines are
−1.947 and −0.879 MeV. The cut excludes 8.56% and 22.25% of the
strict correction to optimistic counts evaluated at the Born
baseline roots. The second temperature selection changes its
conditional root by about −0.028 MeV versus x≤7 despite the
expanded strict support, because much of the newly sampled
mass range overlaps the measured first window. These shifts are
not a bound on omitted QCD terms or a physical reheating result.
The 3 native intermediate x=7.5 boost checks have maximum absolute
relative coefficient errors 0.136% (PCHIP) and 0.172% (linear);
they do not sample the whole strip.

~~~sh
OPENBLAS_NUM_THREADS=1 python -W error partial_hard_vector_nlo_reheat_reference.py   EOS.dat KLOExSec.data BESIII-2112.11728.txt   results/hard-vector-nlo-through-x8-strips-2026-11-26.txt   --intermediate-check results/hard-vector-nlo-through-x8-checks-2026-11-26.txt
OPENBLAS_NUM_THREADS=1 python -W error hard_vector_nlo_overlap_audit_reference.py   EOS.dat results/hard-vector-nlo-through-x8-strips-2026-11-26.txt
~~~
