# Complete Gauss6 warm trajectory and comparison to Gauss5

The independent 255/448/488-node G6 trajectory now passes all five paired epochs from startup through x=.05. Every accepted fine restart and stored terminal state is bound through the uncompressed SHA chain. At x=.05 the fine conditional scalar yield is 0.0003716272903987, the paired relative yield drift is 3.18010e−5, and the largest sampled absolute continuity residual is 2.7741e−11. No chart or physical source is changed between the last two restarts.

At the corresponding native S4 measurement mass, G6/G5 pressure minus one is −0.080152%, −0.080490%, −0.080195%, and −0.081758% at x=.01,.02,.04,.05. At x=.05 the conditioned pressure contrast interval is [−0.00081766261, −0.00081750307]. The relative count change is −1.7310e−9. G6 terminal pressure remains between G4 and G5, and the targeted two-cell alternative also remains between G6 and G5. These are independently evolved warm-only chart choices, not an established extrapolation in quadrature order.

The x=.05 G5/G6 two-order pressure partition gives −0.082081 percentage points from changed cell counts and +0.000323 percentage points from changed pressure per particle. The same two leading warm cells contribute −0.084679 percentage points in total; smaller contributions offset part of the decrease. The auxiliary frozen stress sweep varies a formal measurement mass without evolving the distribution, and does not assign these changes to physical late times or a transfer function.

From `kinetics/`, check the chain and compare the unchanged mass chart:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results/early-pole-fine-charge2-scalar4; d=2026-11-28
python -W error early_pole_endpoint_gauss_warm_chain_receipt.py --order 6 \
  --reports "$p"-gauss6warm-{start,x01,x02,x04,x05}-$d.json.gz \
  --states "$p"-gauss6warm-{start,x01,x02,x04,x05}-state-$d.json.gz
python -W error early_pole_endpoint_warm_gauss_transition_attribution.py "$p-x05-state-$d.json.gz" \
  "$p"-gauss{5,6}warm-x05-state-$d.json.gz --levels 5 6
python -W error early_pole_endpoint_warm_gauss_order_contrast.py --levels 5 6 \
  --lower "$p-gauss5warm-comparison-x"{01,02,04,05}-$d.json \
  --higher "$p-gauss6warm-comparison-x"{01,02,04,05}-$d.json \
  --stress results/early-pole-frozen-stress-sweep-w124-gauss{5,6}-$d.json
```
