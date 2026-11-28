# Fourth-order warm Gauss history: the pressure change reverses

Hold the cold scalar chart, charge-d2 geometry, start condition and clock ladder fixed. On the same 22 complete warm cells, evolve four positive q²-volume Gauss nodes per cell independently. The G4 chart has 211 scalar, 448 daughter and 488 parent nodes. The archived five paired reports and five actual fine states form the complete restart chain from the same empty-daughter start to x=.05; SHA-256 binds each uncompressed restart. Replaying these receipts is not replaying the full evolution.

| x | G4/G3 pressure change | G4/W4 pressure change |
| --: | --: | --: |
| .01 | −0.12355% | −0.61403% |
| .02 | −0.12632% | −0.61516% |
| .04 | −0.12387% | −0.61412% |
| .05 | −0.12539% | −0.61509% |

At x=.05, G3/G2 was +0.08266%, but G4/G3 is negative. G4/G2 is −0.04284%. Each measurement is at the same W1 mass at that epoch. G4/G3 number changes by only −3.65e-9; two methods with the same number of warm nodes (G4 and geometric W4) give different pressures. On each chart separately retaining the node-conditioned count and mean energy still leaves the G4/G3 pressure change in [−0.1253977067%, −0.1253817588%]. These are reconstruction intervals, not bounds on dynamics. The G4 endpoint yield is 0.00037162728968835; its paired time yield drift on the last segment is 3.1801e-5. This receipt does not yet isolate the final-clock dependence of pressure.

The frozen-spectrum probe (hold terminal distributions fixed and vary the *measurement* mass) gives G4/G3 pressure −0.26558% at M=10 and −0.28253% at M=100, and pseudo-pressure −0.56658% and −0.61702%, respectively. These are numerical stress tests of fixed reconstructed spectra. The reversal with warm order is evidence against inferring convergence from the G2/G3 pair alone. It is not a late-time physical pressure history or a CLASS prediction.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
# Start, then use the previous fine state and the same clock steps for .01, .02, .04, .05.
python -W error early_pole_endpoint_higher_warm_gauss_history.py --warm-order 4 --steps 2048 --state /tmp/g4-start.json
p=results/early-pole-fine-charge2-scalar4; d=2026-11-28
python -W error early_pole_endpoint_gauss_warm_chain_receipt.py --order 4 \
  --reports "$p-gauss4warm-"{start,x01,x02,x04,x05}-$d.json.gz \
  --states "$p-gauss4warm-"{start,x01,x02,x04,x05}-state-$d.json.gz
python -W error early_pole_endpoint_warm_gauss_order_contrast.py --levels 3 4 \
  --lower "$p-gauss3warm-comparison-x"{01,02,04,05}-$d.json \
  --higher "$p-gauss4warm-comparison-x"{01,02,04,05}-$d.json \
  --stress results/early-pole-frozen-stress-sweep-w124-gauss3-$d.json \
           results/early-pole-frozen-stress-sweep-w124-gauss4-$d.json
```
