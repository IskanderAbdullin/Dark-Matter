# Finite-band strict-NLO sensitivity of the disjoint hard source

Build on the [disjoint hard Born source](INDEPENDENT-DISJOINT-HARD-THERMAL-SOURCE-2026-11-26.md).
Add the external strict thermal vector coefficient **only** on
$4\le M/T\le5$, using complete fixed-mass strips at 4, 4.5 and 5
and 16 boost points on each. Interpolate the coefficient divided
by the one-flavor vacuum trace in $M/T$ and boost. Keep the same
Bose weight for the visible current, active invariant source
windows, radiation EOS, and optimistic counting convention.
The previously measured KLOE/BESIII vacuum contributions
are left unchanged.

The narrow correction remains between the measured mass regions
at the tested roots (here $M$ is at least 1 GeV and less than
1.69 GeV), so none of this increment overlaps the measured
BESIII windows or the KLOE pion band. Outside the narrow
thermal rectangle, use precisely the old comparator and no
assumption about the unknown NLO term.

| Lower T for hard source | Old disjoint Born root | Root with strict term only on sampled narrow band, fixed αs = 0.30 | Shift |
| ---: | ---: | ---: | ---: |
| 250 MeV | 325.1756 MeV | 324.6607 MeV | −0.5149 MeV |
| 300 MeV | 337.3019 MeV | 337.0306 MeV | −0.2714 MeV |

The NLO-only optimistic count increments at the old roots
are 0.00037897 and 0.00021086, respectively. The two
interpolation/quadrature controls shift the reported root
correction by less than 0.001 MeV. Six independent native
midpoint probes (M/T = 4.25, 4.75; u = 0.7, 2, 5) check the
interpolants. The largest strict-coefficient discrepancies
are 0.718% for linear and 0.318% for shape-preserving interpolation.
This check is not a bound elsewhere in the sampled rectangle.
At fixed illustrative
αs = 0.20–0.40 the 250 MeV cutoff shift is −0.344 to −0.685
MeV, **not** a QCD uncertainty band. In particular the shift
from the arbitrary 250 to 300 MeV hard-source boundary is
far larger in this comparator than the narrow correction
being tested here.

This is a numerical strict-NLO addition using the author's
[external v1.1-spf release](https://doi.org/10.5281/zenodo.10371215).
It does not determine the missing spectral coverage, matching near
the light cone, running coupling or genuine finite-temperature
QCD uncertainty. The 99%-of-0.120 target counts every new particle
as if it could become a scalar. It is not a physical reheating
constraint or relic prediction.

~~~sh
OPENBLAS_NUM_THREADS=1 python -W error partial_hard_vector_nlo_reheat_reference.py \
   EOS.dat KLOExSec.data BESIII-2112.11728.txt \
   results/hard-vector-nlo-three-strips-2026-11-26.txt \
   --intermediate-check results/hard-vector-nlo-interpolation-probes-2026-11-26.txt
~~~

Generate the strip inputs with the external adapter as in
[the two-strip test](INDEPENDENT-HARD-VECTOR-NLO-STRIPS-2026-11-26.md),
including M/T = 5 on the same boost grid for this extension.
