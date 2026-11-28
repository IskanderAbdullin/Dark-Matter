# Strict timelike continuation through M/T = 6

This extends the published
[bounded x=4–5 NLO correction](INDEPENDENT-PARTIAL-HARD-VECTOR-NLO-REHEAT-2026-11-26.md)
to the full native strip at x=M/T=6. Each complete strip contains
16 boost nodes over 0.0125≤(ω−M)/T≤32. The measured-vacuum regions remain
excluded from the hard complement. The temperature quadrature splits
at each crossing of a measured mass edge with an interpolation mass
edge; interpolation never extends outside the supplied grid. The independent x=5.5,
u=0.25 native value differs from the interpolated coefficient by
−0.149% (PCHIP) and −0.356% (linear), expressed as prediction/native−1.
One midpoint does not validate the full boost strip. An exact
regression on x=4,4.5,5 reproduced the previous numerical JSON fields
(controls, native midpoint residuals, and input hashes).

At illustrative αs=0.3, the PCHIP root shifts relative to the same
measured+disjoint Born baseline are −1.102 and −0.582 MeV for hard
lower-temperature cuts of 250 and 300 MeV, respectively, versus
−0.515 and −0.271 MeV through x=5. Refined linear/PCHIP shifts agree
within 0.001 MeV on both temperature windows. This is a sensitivity to
the range actually included, **not** a bound on the omitted spectral
regions and not a physical reheating threshold. The optimistic
all-particles-to-scalars count, fixed coupling, disjoint matching
patches and no-reaction thermal clock are unchanged.

Run with the same physical input files and native spectral interface
as the previous audit:

~~~sh
OPENBLAS_NUM_THREADS=1 python -W error partial_hard_vector_nlo_reheat_reference.py \
  EOS.dat KLOExSec.data BESIII-2112.11728.txt \
  results/hard-vector-nlo-through-x6-strips-2026-11-26.txt \
  --intermediate-check results/hard-vector-nlo-through-x6-checks-2026-11-26.txt
~~~
