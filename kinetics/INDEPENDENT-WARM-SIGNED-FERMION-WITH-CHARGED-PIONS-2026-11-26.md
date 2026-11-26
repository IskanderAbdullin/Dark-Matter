# Include the selected charged-pion channel in the late signed fermion floor

Extend the [signed comparison without a scalar cap](INDEPENDENT-WARM-SIGNED-FERMION-TIME-COMPARISON-2026-11-26.md)
from the actual 10.442 MeV anchors. In addition to
the earlier vacuum charged-lepton terms, allow one
effective massive charged-pion pair source and sink.
Keep the premise that every light and heavy fermion
occupation per spin stays at most one half.
The scalar occupation is not capped.

For the charged pion spectral model use Bose-weighted
KLOE bins on .10–.95 GeV² and the declared vacuum
`|Fπ|²≤4` threshold cap from (4m_π²) to .10 GeV².
The [earlier massive-pion Bose and elastic control](INDEPENDENT-WARM-PION-STATISTICS-AND-THRESHOLD-2026-11-26.md)
is refolded at the actual late anchor, not reused at
50 MeV. Grant `|Fπ|²≤4` also for the independent
spacelike elastic heat comparison. The caps are
not inferred from the KLOE data. Apply finite positive
cooling cells and an analytic threshold-suppressed tail.

For the inverse visible loss into two massive pions,
bound each final Bose factor by (1/(1-e^{-m_π/T})).
The effective maximum spectral ratio is 8.873, but
the specified support ends below (sqrt{s}=.975) GeV.
The corresponding subpole rate ceiling is only
6.270×10⁻¹⁸ GeV⁻². It gains no new near-resonant loss
window. The earlier charged-lepton high-loss budget
can use the small extra pion injection and elastic heat
allowed here. No loss from unseen hadronic channels is
deduced from the measured window or from these limits.

| q32, 2,800-cell comparison | Value |
| --- | ---: |
| Added pion-pair comoving energy above anchor | ≤3.354×10⁻²² GeV⁴ |
| Extra elastic exponent | ≤1.391×10⁻⁷ |
| Fermion Ωh² floor, lepton-only | 1.1192275 |
| Floor with the selected charged-pion support | 1.1192254 |
| Floor allowing the same R at every subpole energy instead | 1.0985244 |

As a separate sensitivity, even assigning 10% more
gross available energy and the wider subpole R gives
a floor .3772. At 20% this conservative comparison
does not exclude .120. These are hypothetical changes
to input allowances, not upper estimates of missing
hadronic or photon rates. Neutral radiation, unknown
visible loss modes, changes to the visible clock, and
the dynamical fermion cap remain outside the result.

~~~sh
export OPENBLAS_NUM_THREADS=1
D=results
python -W error warm_signed_fermion_with_charged_pions_reference.py EOS.dat KLOExSec.data \
 "$D"/warm-isolated-fermion-survival-bound-2026-11-26.json \
 "$D"/warm-fermion-survival-with-leptons-2026-11-26.json \
 "$D"/warm-isolated-recycling-energy-envelope-2026-11-26.json \
 "$D"/warm-late-visible-recycling-envelope-2026-11-26.json
~~~
