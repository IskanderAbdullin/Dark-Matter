# Conditional normalization of the retained species counts

This is a model-budget check, **not a relic prediction**. If the
retained particle numbers at the archived endpoint and the visible
comoving entropy both remain fixed thereafter, their rest masses
give a present-day density proxy. This does not assume every scalar
is cold during structure formation, nor that omitted reactions are
ineffective. It refers only to the stated benchmark masses and
already archived yields.

Compute s_today from the adopted standard-model late entropy degrees in
the [central EOS table](INDEPENDENT-BENCHMARK-SETUP.md) and the CMB
temperature in [Fixsen 2009](https://doi.org/10.1088/0004-637X/707/2/916).
The critical density uses H100=100 km/s/Mpc and the usual SI relation.
The resulting conversion is
Omega h² = 2.75928e8 (m/GeV) (n/s), using the same EOS definition of
visible entropy. The archive yields, originally divided by the
initial comoving visible entropy, are corrected by each independently
tracked final visible entropy ratio.

| 50 MeV reheating archive | scalar fixed-count Omega h² | stable fermion fixed-count Omega h² |
| --- | ---: | ---: |
| 12 MeV endpoint | 0.00014579 | 88.1776 |
| 9.3 MeV endpoint | 0.00014576 | 88.4445 |

The stable component alone is about 735 to 737 times the
[Planck base-Lambda-CDM cold-dark-matter reference](https://arxiv.org/abs/1807.06209)
of 0.120. If there were no extra dilution, and using the entire
observational budget for this stable species, more than 99.86% of the
retained stable number would have to be depleted relative to this
number-conserving benchmark. This quantifies the burden on omitted
annihilation/number change or an altered reheating history. It does
not show that such changes are impossible, or constrain arbitrary
portal parameters.

Separately integrate the full uncut vacuum pair source in the
prescribed tabulated radiation down to 1 MeV. Within that source-only,
number-conserving comparator, the stable proxy intersects 0.120 at
T_R about 9.02 MeV. The associated scalar proxy is only 2.36e-12.
Thus tuning this restricted source to the stable-count normalization
cannot simultaneously be read as a dominant scalar population.
This low-temperature root is a follow-up benchmark, not a
self-consistent solved kinetic reheating scenario. Two independent
moment/time quadrature settings agree at the root. At T_R=20 and
50 MeV the vacuum stable proxies are 2.77 and 88.65 respectively.

## Reproduce

Prepare [both EOS files](INDEPENDENT-BENCHMARK-SETUP.md). From
`kinetics`, run:

```sh
OPENBLAS_NUM_THREADS=1 python -W error conditional_relic_normalization_reference.py EOS.dat EOS.late.dat results/12-quantum-48-2026-11-26.json results/early-9p3-quantum-2026-11-26.json
```

Pass `--no-sweep` for the quick archived-number conversion only.
