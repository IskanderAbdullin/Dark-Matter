# Actual G8 continuation through x=.04

The G8 warm-order continuation from the accepted x=.02 state completed the 1024/2048 pair through x=.04. The scalar yield was 1.73209229123e-4, with relative paired time-step drift 1.27671e-4, maximum absolute continuity defect 7.21e-11, and maximum pair-number residual 1.12e-16. The receipt joins the four clocks to exact native state hashes and checks both G7 and G8 restarts.

At the common native measurement mass, the G8 pressure was 0.0482037% below G7, and total count differed by -1.164e-9 in relative units. The paired-step envelope was [-0.0482085%, -0.0481988%]. Redistribution among native momentum cells contributed -0.0458841 percentage points; the within-cell change contributed -0.00231955. The interval q=[3.96367894, 5.72884111] accounted for about 71.4% of the net pressure difference. The sign is stable over these four matched clocks; this does not establish the continuum limit or extend the actual G8 trajectory beyond x=.04.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results/early-pole-fine-charge2-scalar4; d=2026-11-28
python -W error early_pole_endpoint_warm_g7_prefix_check.py --levels 7 8 \
 --native "$p"-{start,x01,x02,x04}-state-$d.json.gz --native-chain "$p-chain-$d.json" \
 --reports6 "$p"-gauss7warm-{start,x01,x02,x04}-$d.json.gz --states6 "$p"-gauss7warm-{start,x01,x02,x04}-state-$d.json.gz \
 --reports7 "$p"-gauss8warm-{start,x01,x02,x04}-$d.json.gz --states7 "$p"-gauss8warm-{start,x01,x02,x04}-state-$d.json.gz
```
