# Where cumulative vacuum production occurs after 50 MeV reheating

The early expanding controls stop at 45 MeV. Integrating the vacuum visible
source in log visible temperature, without a dark momentum grid, provides
an independent check on how much production that short run samples. The
assumptions are a fixed visible adiabatic background, heavy parents
promptly decaying, vacuum dark final states, and no inverse portal terms.
These are integrated number source terms, not final relic abundances.

The integrand for each pair channel is 2 R c_V/(3 H s²) in log(T), with
R the vacuum pair production rate, c_V = T ds/dT and s the visible entropy
density. The factor of 2 counts the two produced particles. For a heavy
pair the same yield sources both the scalar and the stable daughter.
Split the temperature integral at EOS knots. Cross-check the 45 MeV
yield against the separately implemented binned collisionless spectrum.

| Visible endpoint | Direct light count / cumulative to 12 MeV | Heavy-prompt count / cumulative to 12 MeV |
| ---: | ---: | ---: |
| 45 MeV | 0.2868 | 0.3709 |
| 35 MeV | 0.6899 | 0.8248 |
| 30 MeV | 0.8169 | 0.9316 |
| 25 MeV | 0.9035 | 0.9816 |
| 20 MeV | 0.9583 | 0.9974 |
| 15 MeV | 0.9896 | 0.9999 |

The denominator is the computed count at 12 MeV, not at zero
temperature; direct production after 12 MeV is not independently bounded
by this result. Dark-sector transport and occupation effects can change
the meaning of a cumulative **source** when projected into colder momentum
bins. Matching both number and energy while evolving through this
temperature range is a distinct follow-up.

```sh
cd kinetics
OPENBLAS_NUM_THREADS=1 python -W error visible_vacuum_source_cumulative_reference.py EOS.dat --spectra empty.json
```
