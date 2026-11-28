# The actual G7 startup gives a positive pressure change relative to G6

On the same cold and charged chart, G7 uses 277/448/488 nodes and the independently checked orthogonal-basis positive volume rule. Its x=.0075 startup completes the 2048/4096 paired clock, with yield 5.7823259971e−8, paired yield time drift 5.11556e−5 and absolute continuity residual below 1.78e−16. The trace and terminal bytes are bound by their SHA receipts.

On the common native S4 measurement mass, G7 relative to G6 has count +1.81e−8 and pressure +0.0590817%. The positive relative pressure interval at the respective stored cell counts and mean energies is [0.000590816397,0.000590816879]. The change is opposite to G6/G5 and follows the same sign as the earlier held-source test. This verifies only startup, not later trajectory stages, a continuum extrapolation or physical pressure.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results/early-pole-fine-charge2-scalar4; d=2026-11-28
python -W error early_pole_endpoint_warm_g7_startup_check.py \
 --native "$p-start-state-$d.json.gz" --native-chain "$p-chain-$d.json" \
 --reports "$p"-gauss{6,7}warm-start-$d.json.gz \
 --states "$p"-gauss{6,7}warm-start-state-$d.json.gz
```
