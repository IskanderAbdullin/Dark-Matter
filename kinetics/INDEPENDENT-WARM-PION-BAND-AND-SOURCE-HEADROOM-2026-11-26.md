# Test the measured two-pion window and budget omitted sources at 50 MeV

The [from-empty comparison](INDEPENDENT-WARM-FROM-EMPTY-RECYCLING-ENVELOPE-2026-11-26.md)
used positive vacuum charged-lepton production and elastic heating.
Fold the [published KLOE combination](https://www.lnf.infn.it/kloe/ppg/ppg_2017/KLOEppg2017/CrossSection/Combination/KLOExSec.data)
as a separate piecewise-constant **exclusive** `π⁺π⁻` vacuum
input on its 85 bins, s=.10–.95 GeV². Use the tabulated cross
section plus its tabulated per-bin error as a defined variant
of the source. This simultaneous shift is not a statistical
confidence band or a model of finite-temperature hadrons.

For a specified cross section in each bin, subdivide in s and
use positive upper products: powers, phase factors, and the
sub-pole propagator at the upper endpoints; the Bessel functions
at the lower endpoints. A separate bin quadrature at 50 MeV
checks that the finest positive cells lie about 1.3% above
the independently integrated input number and energy moments.
The normalization is obtained from the same vacuum electromagnetic
spectral response as the [earlier high-temperature band fold](INDEPENDENT-PARTIAL-HADRON-REHEAT-SOURCE-2026-11-26.md).
Credit two direct stable particles or four decay daughters per
produced pair, and allow all positive redistributed energy to
feed recycling. Repeat on three matched source and cooling grids.

## How much unknown energy would make this comparison inconclusive?

For the already included lepton and pion rates let `P(b)` be
the accumulated upper comoving source energy. Allow an extra
positive source with total budget `E_X` and no change to the
visible clock. Use the deliberately optimistic deposit
`P(b)+E_X` at every b in the recycling integral; the expansion
weight then gives a positive quadratic in `E_X`. For an
extra portal fermion source, grant up to `E_X/m₂` direct
particles including any prompt daughters. For extra elastic
heat omit that direct contribution. These controls do not
apply to injection of new light scalar particles at a much
lower energy cost per particle.

| 50 MeV start, central radiation table | Scalar Ωh² upper |
| --- | ---: |
| Vacuum charged leptons only | .02154169 |
| Charged leptons plus measured-band model at its own rate | .02717466 |
| Same prescribed pion energy credited immediately, direct count overestimated | .02785642 |

The prescribed pion energy input upper is 11.70% of the
lepton energy input upper. With no pion channel included,
the always-early extra energy needed for the comparison to
reach .060 is 57.26% of that lepton energy upper (59.88% if it
adds no direct particles). After adding the prescribed pion history and direct count,
its remaining budget for extra portal input is 1.153×10⁻⁹
GeV⁴ in start-normalized comoving energy (46.77% of the
lepton energy upper).
These are **sufficient budgets for the upper-bound
comparison**, not a measured requirement on all hadronic
sources and not a prediction of a successful relic abundance.
The selected band model stays below .060 if the total of
interest is .120. The uncovered `ππ` threshold strip, other
exclusive/inclusive hadronic channels, hadronic elastic heat,
medium corrections and new dark vertices are not bounded here.

~~~sh
export OPENBLAS_NUM_THREADS=1
D=results
python -W error warm_extra_injection_headroom_reference.py EOS.dat KLOExSec.data \
  "$D"/warm-energy-weighted-recycling-kernels-2026-11-26.json
~~~
