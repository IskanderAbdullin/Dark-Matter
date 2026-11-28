# Direct prompt production couples the scalar and stable abundance budgets

A heavy particle or antiparticle decays into one stable
fermion or antifermion and one scalar. With unchanged
comoving counts and common entropy, the stable rest density
from exactly those parents is m_stable/m_scalar times the
scalar rest density. This relation does not depend on the
heating clock, portal normalization or daughter momenta.
For masses 20 MeV and 86 eV the ratio is 232,558.
Even if no other stable particles are produced, the scalar
mass share of the two daughter species is only about
4.30 parts per million under count conservation. Extra
independently produced stable particles lower that share.
A common entropy dilution cannot change this ratio.

| Start (MeV) | Stable omega from same parents | Separate central stable-pair omega | Fixed scalar: required stable loss | Complete conversion: scalar omega, both central |
|---:|---:|---:|---:|---:|
| 25 | 0.62809 | 5.6604 | 98.091800% | 2.9741e-05 |
| 40 | 12.586 | 26.75 | 99.695074% | 0.00022326 |
| 45 | 22.614 | 39.149 | 99.805866% | 0.00036282 |
| 50 | 36.971 | 54.724 | 99.869305% | 0.00055326 |
| 55.5 | 58.942 | 75.807 | 99.911133% | 0.00083287 |
| 60 | 82.799 | 96.337 | 99.933210% | 0.0011263 |
| 75 | 213.92 | 188.22 | 99.970388% | 0.002649 |
| 100 | 678.3 | 433.95 | 99.989473% | 0.0076994 |

There is a different upper count comparison if one also
allows stable-antistable to scalar pairs in the post-decay
population. In this restricted network these two-to-two
conversions conserve the combined scalar plus stable
number, so even perfect conversion grants at most one
additional scalar per stable produced. Use the central
direct vacuum charged-lepton source for stable pairs and
separate it from the stable daughters of heavy parents.
The stated upper comparison using an upper parent count
and **central** stable pair source is conditional on fixing
that second source; it is not an upper on its uncertainty
or on unspecified stable production. Even this favorable count
comparison at 100 MeV gives at most .01026, versus .060
for half the stated dark density. With both central sources
the complete-conversion value is .007699. This bound grants
all of the paired number even if the underlying momentum,
inverse rates and timing prevent it from occurring, and
in that complete-conversion limit there is no stable
rest density left. Reheating outside the table, parent
regeneration, non-two-to-two reactions, or other sources
can evade this restricted count budget.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error warm_direct_coproduction_budget_reference.py EOS.late.dat \
  results/warm-direct-vacuum-temperature-bracket-2026-11-26.json
~~~

[Full matched count ledger](results/warm-direct-coproduction-budget-2026-11-26.json).

The [weighted count upper](INDEPENDENT-WARM-DIRECT-RESTRICTED-COUNT-UPPER-2026-11-26.md)
uses positive envelopes for both light and heavy sources, including the stated
charged-pion input into each, instead of fixing the central direct light source.
