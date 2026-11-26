# Strict timelike continuation through M/T=12 with the missing strip filled

The native main sweep jumped from x=10 to x=12. The x=11
strip was separately computed on all 16 matching native boost
nodes before either new interval was used. The same strict
coefficient, disjoint measured regions and fixed-coupling
optimistic count benchmark have been retained.

At αs=0.3 the resulting conditional roots through x=12 are
323.059 and 336.303 MeV for hard lower temperature cuts 250
and 300 MeV. These shift by −0.100 and −0.079 MeV relative
to x≤10, or roughly −2.117 and −0.999 MeV relative to the
measured-plus-Born reference. The x=11 intermediate controls
are 323.101 and 336.329 MeV. These are finite-strip changes;
they neither bound omitted corrections nor measure a physical
reheating condition.

The new support includes the gap between the two measured bands,
then reaches above the upper band once xT exceeds 3.671 GeV.
For the x=12 edge this first happens at 305.917 MeV. The
independent physical-mass integration maintains exact exclusions
for both windows. It does not fill any absent high strip or the
measured resonance intervals with an inferred source.
The four separate x=11.5 native probes have a maximum absolute
coefficient interpolation difference of 0.032% (PCHIP) and
0.076% (linear); this is not a bound across the whole strip.

~~~sh
OPENBLAS_NUM_THREADS=1 python -W error partial_hard_vector_nlo_reheat_reference.py \
  EOS.dat KLOExSec.data BESIII-2112.11728.txt \
  results/hard-vector-nlo-through-x12-strips-2026-11-26.txt \
  --intermediate-check results/hard-vector-nlo-through-x12-checks-2026-11-26.txt
OPENBLAS_NUM_THREADS=1 python -W error hard_vector_nlo_overlap_audit_reference.py \
  EOS.dat results/hard-vector-nlo-through-x12-strips-2026-11-26.txt
OPENBLAS_NUM_THREADS=1 python -W error hard_vector_measured_window_saturation_reference.py \
  results/partial-hard-vector-nlo-through-x10-reheat-2026-11-26.json \
  results/partial-hard-vector-nlo-through-x12-reheat-2026-11-26.json
~~~
