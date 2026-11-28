# Strict timelike strip through M/T=10: reopened support above the lower measured window

The [x=9 plateau](INDEPENDENT-HARD-VECTOR-NLO-X9-2026-11-26.md)
was due to exact overlap geometry at the 300 MeV lower temperature
cut, not convergence of the full thermal correction.

At x=9–10 the new strip enters the gap above the first BESIII band.
The upper new mass edge first passes 3.08 GeV at 308 MeV.
At the original Born roots, the x=9–10 portions between the two
measured bands span 3.08–3.252 and 3.08–3.373 GeV at the lower
temperature cuts 250 and 300 MeV, respectively. These are
**unmeasured in the adopted band selection**, not placeholders
for an independently measured continuum.

At illustrative αs=0.3 the conditional roots through x=10 are
323.159 and 336.382 MeV. Shifts from the original measured
plus disjoint Born roots are −2.017 and −0.920 MeV. Relative
to x≤9 the 300 MeV cut now moves by −0.0412 MeV instead
of an exact-overlap plateau; the 250 MeV cut moves by −0.0116 MeV.
This remains a bounded, fixed-coupling correction on the stated
boost grid and the same optimistic count. No full-QCD, thermal
threshold, or observational statement follows.
At the one independent x=9.5, u=0.25 native probe, the
coefficient residual is +0.027% (PCHIP) or
+0.050% (linear), not a full strip accuracy estimate.

~~~sh
OPENBLAS_NUM_THREADS=1 python -W error partial_hard_vector_nlo_reheat_reference.py \
  EOS.dat KLOExSec.data BESIII-2112.11728.txt \
  results/hard-vector-nlo-through-x10-strips-2026-11-26.txt \
  --intermediate-check results/hard-vector-nlo-through-x10-checks-2026-11-26.txt
OPENBLAS_NUM_THREADS=1 python -W error hard_vector_nlo_overlap_audit_reference.py \
  EOS.dat results/hard-vector-nlo-through-x10-strips-2026-11-26.txt
OPENBLAS_NUM_THREADS=1 python -W error hard_vector_measured_window_saturation_reference.py \
  results/partial-hard-vector-nlo-through-x9-reheat-2026-11-26.json \
  results/partial-hard-vector-nlo-through-x10-reheat-2026-11-26.json
~~~
