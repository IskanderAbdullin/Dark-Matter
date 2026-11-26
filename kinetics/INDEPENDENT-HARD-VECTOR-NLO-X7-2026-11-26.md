# Strict timelike continuation through M/T=7, with a nonzero window exclusion

This extends the [bounded correction through x=6](INDEPENDENT-HARD-VECTOR-NLO-X6-2026-11-26.md)
with a complete native vector strip at x=M/T=7. It retains the same hard
boost range, visible radiation table and optimistic conversion
target. No region outside the native strips is extrapolated.

The top mass edge now crosses the first BESIII window at
T=318.9143 MeV, inside both reheating integrals. In an independent
integration in physical mass M, every panel is assigned either the
excluded or included contribution, before summation. The split result
reproduces the production-code comparator in x to within 1.61×10⁻⁷
of the unsplit local rate (counts and energy, both channels).
The independent 16-to-24 point comparison changes by at most
8.72×10⁻⁹ of the unsplit rate. At T just above the crossing, the
newly blocked heavy and light counts are positive; below it, exactly zero.
The blocked first window starts ending below the x=7 edge at T=440 MeV,
also tested by the local audit.

At αs=0.3, the shifted conditional roots are 323.520 and
336.452 MeV, shifts of −1.656 and −0.850 MeV relative to
the same measured+disjoint Born baseline. The mass-window cut
removes 0.208% and 3.282% of the strict correction to the
optimistic count at those **Born** baseline roots. These numbers
depend on the chosen temperature cuts of 250 and 300 MeV.
On the 2 independent x=6.5 boost probes included here, PCHIP
and linear relative interpolation errors are at most 0.210%
and 0.269%, respectively. These points check only parts
of the new strip.
The comparators do not determine a physical reheating threshold,
full NLO correction, or relic abundance.

~~~sh
OPENBLAS_NUM_THREADS=1 python -W error partial_hard_vector_nlo_reheat_reference.py \
  EOS.dat KLOExSec.data BESIII-2112.11728.txt \
  results/hard-vector-nlo-through-x7-strips-2026-11-26.txt \
  --intermediate-check results/hard-vector-nlo-through-x7-checks-2026-11-26.txt
OPENBLAS_NUM_THREADS=1 python -W error hard_vector_nlo_overlap_audit_reference.py \
  EOS.dat results/hard-vector-nlo-through-x7-strips-2026-11-26.txt
~~~
