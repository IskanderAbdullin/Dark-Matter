# A total-number ceiling on scalar dark matter after finite production

The archived stable fermion proxy can greatly exceed the observational
dark-matter budget, but its conversion into scalars does **not** carry
the same rest mass across. For the benchmark, mφ/mχ₂ = 4.3e-6.
χ₂χ̄₂ ↔ φφ interchanges two incoming and two outgoing particles.
Bose stimulation, Pauli blocking, momentum redistribution, and
Bose condensation do not increase total number for this reaction.
They can affect the time, particle distributions, and cold fraction.

Assume fixed total dark count per visible comoving entropy at the
specified checkpoint, or a subsequent decrease. Then even assigning
**every** available particle to the scalar leaves

`(Ωφ h²)max = (Ωφ h²)initial + (mφ/mχ₂)(Ωχ h²)initial`.

This is independent of the conversion matrix element or momentum
shape and is an upper bound on *rest mass*. Relativistic kinetic
energy is not part of the cold-dark-matter inventory measured by
[Planck](https://arxiv.org/abs/1807.06209). If some scalars remain
too warm or fail to cluster, the relevant bound is lower.

| Initial count basis | Largest scalar Ωh² | Fraction of 0.120 | Number factor needed for 99% of 0.120 |
| --- | ---: | ---: | ---: |
| 50 MeV history, archive at 12 MeV | 5.250e-4 | 0.437% | 226 |
| 50 MeV history, archive at 9.3 MeV | 5.261e-4 | 0.438% | 226 |
| 50 MeV prescribed full vacuum source through 1 MeV | 5.271e-4 | 0.439% | 225 |
| 9.019 MeV prescribed vacuum source root | 5.160e-7 | 0.000430% | 230,232 |

The full-vacuum-source rows are separate conditional input histories
from the [abundance comparison](INDEPENDENT-CONDITIONAL-RELIC-NORMALIZATION-2026-11-26.md),
not replacements for the archives or bounds on arbitrary interacting
production. For example, if arbitrary conversion of only the two
species in the 12 MeV archive were adjusted to leave total rest-mass
Ωh²=0.120, about 0.119476 would still be stable fermions and only
0.0005244 would be scalars. This is an algebraic allocation,
not a demonstration that the rate can reach it.

Thus further work on scalar dominance in this benchmark must quantify
a mechanism that increases particle number per visible entropy,
supplies an additional scalar population, or changes parameters or
initial conditions. Simply accelerating this number-conserving
conversion cannot meet that goal. The ceiling does not include new
visible injection after an archive or number-changing processes and
does not constrain different models under different assumptions.

```sh
OPENBLAS_NUM_THREADS=1 python -W error scalar_number_budget_reference.py results/conditional-relic-normalization-2026-11-26.json
```
