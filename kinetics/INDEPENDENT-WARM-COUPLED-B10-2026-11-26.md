# Warm signed pair continuation to 1.206 MeV

Restart both [3.43 MeV solutions](INDEPENDENT-WARM-COUPLED-B35-2026-11-26.md)
from their own retained distributions. Expand each by (10/3.5)
on the same radiation clock. Account for all exterior particle
production in the normalization without recycling those particles
as collision partners or entering their energy in the new stage's
Hubble rate. The stage starts at 3.433865 MeV and ends at
1.206135 MeV.

| Log-scale steps from 12 to 3.43 MeV, then to 1.206 MeV | Extra stable number loss after 3.43 MeV | Stable fixed-count Ωh² | Scalar fixed-count Ωh² |
| --- | ---: | ---: | ---: |
| 20 + 18 | 6.9690% | 57.60930 | 0.000277281 |
| 32 + 28 | 6.9768% | 57.44724 | 0.000277978 |

The final stable proxies differ by about 0.282%; the loss accrued
in this stage differs by only 0.0078 percentage points. As another
check, 48 steps in the *first* stage give stable Ωh²=61.6640
at 3.43 MeV, compared with 61.7558 at 32 steps. The 48-step
continuation below that temperature is not included here.

In the 32+28 case the total stable proxy has declined 34.1%
relative to its 12 MeV value and remains about 479 times .120.
The endpoint peak resolved scalar occupation is (7.06	imes10^{-4});
its scalar count-to-equivalent-saturation ratio is (6.54	imes10^{-4}).
The gross reverse/forward ratio peaks at 0.021 in the second stage.
The relative energy defect in the refined restarted run is
(7.34	imes10^{-13}), total particle count defect
(1.99	imes10^{-13}).

This comparison does not remove the fixed momentum and collision
quadrature limitations. The prescribed *vacuum forward* direct
visible source omitted below 12 MeV alone corresponds to a stable
proxy of about 0.42 if all its particles survive, a scale larger than
the present step refinement. It is neither a net-source prediction
nor an error bound on this nonlinear result. The second stage still
omits this residual source, visible feedback, heavy recycling,
number-changing scalar processes and exterior re-entry.
Do not identify the endpoint with a measured abundance or the
full asymptotic yield.

```sh
export OPENBLAS_NUM_THREADS=1
a=results/warm-coupled-pair-b35
for pair in '20 18' '32 28'; do
  read -r earlier later <<< "$pair"
  python -W error coupled_warm_pair_continuation_reference.py EOS.dat \
    "$a-$earlier-shape-2026-11-26.json" --stop 2.857142857142857 \
    --steps "$later" --save-snapshots "shape-$later.json" >"rows-$later.json"
done
python -W error warm_coupled_count_ledger_reference.py \
  results/50MeV-from-empty-pair-budget-2026-11-26.json \
  results/50MeV-from-empty-pair-quantum-shape-2026-11-26.json \
  --results "$a"-{20,32}-2026-11-26.json \
  --bridges "$a"-{20,32}-shape-2026-11-26.json \
  --follow-results rows-{18,28}.json --follow-bridges shape-{18,28}.json
```
