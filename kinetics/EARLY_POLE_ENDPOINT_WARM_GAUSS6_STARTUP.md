# Direct Gauss6 warm startup on the matched paired clock

The independently integrated Gauss6 startup is now checked and archived at x=.0075. It uses 255 scalar, 448 negative-charge and 488 positive-charge nodes, with 2048/4096 paired time steps. The higher-order chain validator checks this as one complete epoch with no predecessor and an uncompressed saved-state SHA. Through this epoch, the paired scalar-yield time drift is 5.12e−5 and the audited absolute continuity residual is below 2.55e−16. The later restart segments are separate work.

At the single native S4 startup mass, the whole-band pressure numerator changes as follows under the same moment-conditioned reconstruction:

| Compared startup charts | Pressure ratio minus one | Count ratio minus one |
| --- | ---: | ---: |
| G3 / G2 | +0.0923753% | +1.89e−7 |
| G4 / G3 | −0.1235425% | −4.15e−8 |
| G5 / G4 | +0.1120621% | +3.69e−8 |
| G6 / G5 | −0.0801770% | −2.48e−8 |
| targeted G4 / G6 | +0.0465629% | +1.22e−8 |

G6 lies between G4 and G5 at this epoch, and the targeted chart lies between G6 and G5. The signs survive the separate positive pressure envelopes conditional on each saved count and mean energy. This is an actual startup evolution comparison, not a frozen derivative and not a test of the later trajectories, continuum convergence, or physical pressure.

The direct cell-by-cell comparison also assigns most of the G6/G5 pressure difference to the same two warm cells used in the targeted rule. Together they contribute −0.082368 percentage points to the net −0.080177%. They slightly over-account for the decrease because other cells offset part of it. Dividing each cell pressure by its own count and symmetrizing the two replacement paths assigns −0.080497 percentage points to changed counts and +0.000320 percentage points to changed pressure per particle. This is a decomposition of two archived endpoint populations, not an intervention on the evolution.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28; p=results/early-pole-fine-charge2-scalar4
python -W error early_pole_endpoint_higher_warm_gauss_history.py \
 --warm-order 6 --x-end .0075 --steps 2048 --state /tmp/G6-start-state.json >/tmp/G6-start.json
python -W error early_pole_endpoint_gauss_warm_chain_receipt.py --prefix --order 6 \
 --reports /tmp/G6-start.json --states /tmp/G6-start-state.json
python -W error early_pole_endpoint_warm_startup_order_contrast.py --native "$p-start-state-$d.json.gz" \
 --native-chain "$p-chain-$d.json" \
 --reports "$p"-gauss{,3,4,5,6}warm-start-$d.json.gz "$p-targetedwarm-start-$d.json.gz" \
 --states "$p"-gauss{,3,4,5,6}warm-start-state-$d.json.gz "$p-targetedwarm-start-state-$d.json.gz"
```
