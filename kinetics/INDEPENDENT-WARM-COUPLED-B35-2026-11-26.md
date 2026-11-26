# Warm benchmark: extend signed pair conversion to 3.43 MeV

Use the signed visible-source [50 MeV history](INDEPENDENT-50MEV-PAIR-FROM-EMPTY-2026-11-26.md)
as the common initial condition at 12.00 MeV. Continue the [coupled
warm benchmark](INDEPENDENT-WARM-COUPLED-LATE-PAIRS-2026-11-26.md)
to relative expansion (b=3.5). The radiation clock gives
(T_gamma=3.43387) MeV. Pair conversion is two signed midpoint
kicks around the conservative implicit prompt-plus-elastic collision
step. The newly escaping energy is tracked in the expansion ledger.

| Log-scale steps | Stable count at endpoint / starting retained count | Scalar count at endpoint / starting retained count | Conditional fixed-count stable Ωh² | Conditional fixed-count scalar Ωh² |
| ---: | ---: | ---: | ---: | ---: |
| 20 | 0.710208 | 1.723800 | 61.92487 | 0.000258724 |
| 32 | 0.708269 | 1.728643 | 61.75584 | 0.000259451 |

The count ratios include newly generated exterior particles. For the
fixed-count normalizations, also carry the exterior particle counts from
before 12 MeV unchanged, **not** as collision partners. The audit
reconstructs these contributions separately using archived total counts
and independently verifies each restart shape's retained count. The
inclusive fixed-count stable proxy declines about 29.2% relative to its
12 MeV value of 87.1873. It is still about 515 times .120, and not a
late-time prediction.

Refining the time grid from 20 to 32 intervals changes the fixed-count
stable and scalar proxies by about 0.27% and 0.28%, respectively.
This is a two-grid time test at fixed momentum support and collision
quadrature, not an error bound for removing those restrictions.
The relative energy ledger is below (6	imes10^{-13}) in both
runs. The reverse/forward gross midpoint ratio reaches about 0.85.
The maximum resolved scalar occupation is about (6.1	imes10^{-4})
and the scalar count-to-equivalent-saturation ratio is about
(4.85	imes10^{-4}) in the 32-step state. Its first cell contains
(3.16	imes10^{-8}) of the original retained scalar number, or
about (1.83	imes10^{-8}) of its own resolved population.

Neglect residual visible injection and feedback, heavy recycling,
number-changing scalar reactions and re-entry from exterior cohorts.
Earlier exterior cohorts are carried in the count audit but do not
affect the Hubble rate in this stage. Shape exports contain *retained*
cells only; passing them to another stage does not automatically carry
any particle tails. The conditional normalization assumes no change in
visible comoving entropy.

Reproduce with the official radiation table described in
[the common setup](INDEPENDENT-BENCHMARK-SETUP.md):
```sh
export OPENBLAS_NUM_THREADS=1
start=results/50MeV-from-empty-pair-quantum-shape-2026-11-26.json
for n in 20 32; do
  python -W error coupled_warm_pair_continuation_reference.py EOS.dat "$start" \
    --stop 3.5 --steps "$n" --save-snapshots "shape-${n}.json" > "rows-${n}.json"
done
python -W error warm_coupled_count_ledger_reference.py \
  results/50MeV-from-empty-pair-budget-2026-11-26.json "$start" \
  --results rows-20.json rows-32.json --bridges shape-20.json shape-32.json
```
