# Last-segment clock sensitivity of the third D8 fitted poles

Recover the 1024-step endpoint from the exact x=.04 restart used for the existing 2048-step terminal state. The recovered coarse trace exactly matches the archived five audited columns at every sampled time. Both actual terminal arrays then go through the same number/energy closure and transverse-pole refit at k=.0015. This compares a common initial state and differingly resolved **last** interval, not all accumulated clock error since startup.

| Quantity at x=.05 | Coarse/fine minus one |
| --- | ---: |
| Scalar count | −3.26824563e−5 |
| Fitted scalar gap | −1.08943804e−5 |
| Fitted transverse photon energy | −1.10873533e−11 |
| Fitted transverse residue | +1.10931264e−11 |

The photon fit hardly moves here, but the scalar gap changes at roughly a third of the count drift in relative magnitude. The [cold-pole conditioning test](INDEPENDENT-EARLY-POLE-ENDPOINT-COLD-CONDITION-NUMBERS-2026-11-28.md) was on a different chart, so cannot convert this optical drift into a rate error bound. Compare actual cold coefficients for the saved pair rather than carrying over a sensitivity from fine D2. The complete packed coarse state is retained alongside this receipt.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results/early-pole-cell-endpoint-third-d8; d=2026-11-28
python -c 'import json,sys; print(json.dumps(json.load(open(sys.argv[1]))["segments"][-1]))' \
  "$p-to-x05-$d.json" > /tmp/third-d8-final-report.json
python -W error early_pole_endpoint_charge_cold_clock.py \
  "$p-x04-state-$d.json.gz" "$p-x05-state-$d.json.gz" /tmp/third-d8-final-report.json \
  --coarse /tmp/third-d8-coarse-final.json
```
