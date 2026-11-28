# Instantaneous elastic feed to the first scalar shell

Evaluate the nonlinear below-pole elastic flux on the populated
q32 archived distribution and its exactly cropped q16 interior.
Use one common H and the same surviving coefficients in both,
rather than mixing a background-normalization change into the
boundary difference. No time evolution is applied in this check.

The first-bin scalar number slope divided by its own number and H
is 0.27993 on q32 and 0.27990 on q16 in the finest shown quadrature.
The first-bin kinetic slopes are similarly positive. The change from
including the four exterior cells is smaller than the spread over
the three quadrature levels. This local cold gain can coexist with
a net loss of scalar kinetic energy over all cells: in the same
q32 check the bulk elastic scalar heating rate is -0.0430 times
the total dark kinetic energy per Hubble time. The opposite fermion
budget and both separate counts close to numerical precision.

The first bin is kinematically closed for direct on-shell coupling
on both supports. The positive elastic number slope here is a net
local redistribution on one initial shape, not an eigenvalue, a
constant extrapolatable source, or a condensate formation rate.

```sh
cd kinetics
OPENBLAS_NUM_THREADS=1 python -W error archived_cold_elastic_budget_reference.py EOS.dat EOS.late.dat results/wide32-embedded-vacuum-2026-11-26.json results/wide32-embedded-vacuum-shape-2026-11-26.json
```
Prepare the [external EOS inputs](INDEPENDENT-BENCHMARK-SETUP.md).
