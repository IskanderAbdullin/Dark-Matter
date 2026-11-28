# Actual G9 warm startup reverses the earlier G8 pressure trend

The positive G9 warm quadrature completes its own 2048/4096 startup to x=.0075 from the same empty-daughter problem, not by lifting a terminal G8 state. The paired scalar-yield drift is 5.116e-5, the maximum absolute continuity defect 2.31e-16, and maximum paired-number residual 6.52e-17.

At the same native S4 measurement mass G9 pressure is 0.0273012% higher than G8; its total count differs by only +7.96e-9 in relative units. About 84.9% of the pressure difference is in the original cell q=[3.96368,5.72884], which independently dominated the frozen source probe. This sign is opposite to the earlier G7→G8 pressure change at startup. A reversal at one clock is not a convergence result. The conditional pressure interval refers to original-cell count and energy, not time-stepping error.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28; p=results/early-pole-fine-charge2-scalar4
python -W error early_pole_endpoint_warm_g7_prefix_check.py --levels 8 9 \
 --native "$p-start-state-$d.json.gz" --native-chain "$p-chain-$d.json" \
 --reports6 "$p-gauss8warm-start-$d.json.gz" --states6 "$p-gauss8warm-start-state-$d.json.gz" \
 --reports7 "$p-gauss9warm-start-$d.json.gz" --states7 "$p-gauss9warm-start-state-$d.json.gz"
```

The clock family retains its optional exact substep checkpoint support. This is only the actual initial numerical prefix, not a fully evolved G9 history or a physical prediction.
