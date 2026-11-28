# Swap quadrature and stored histories separately at the warm terminal endpoint

Both G7 and mixed255 reach x=.05 with accepted paired traces. Freeze the
native S4 optical mass and clock. On each rule's low node set, evaluate
all combinations of the two positive scalar moment reconstructions
and the two actual charged occupancies. Separately redo all four
scalar/charge combinations on the common Volume160 rule. Both diagonal
Volume160 calculations agree with independently archived references.

The averaged log effects below hold at this one endpoint. The average
uses all six orders of low-cube substitution, retaining the range
of individual paths in the receipt.

| Swapped input | Warm forward log change | Warm inverse log change | Warm net log change |
|---|---:|---:|---:|
| Node rule | −1.21269e−4 | −1.98889e−6 | −1.21269e−4 |
| Represented scalar history | +4.59e−13 | −8.52132e−7 | +1.10933e−12 |
| Charged history | +6.94152e−8 | −6.99244e−7 | +6.94158e−8 |

The terminal warm net difference is therefore dominated by node rule
on this numerical control. For the inverse term, the history effects
matter appreciably. When the quadrature is held fixed to Volume160,
the net diagonal changes by only +6.96e−8. These counterfactuals neither
replay the coupled histories nor identify physical causal contributions.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results/early-pole-fine-charge2-scalar4; d=2026-11-28
python -W error early_pole_endpoint_warm_actual_source_factor_cube.py \
 --native "$p-x05-state-$d.json.gz" --states "$p"-{gauss7,mixed255}warm-x05-state-$d.json.gz \
 --control results/early-pole-warm-budgeted-gauss-probe-$d.json \
 --prefix results/early-pole-warm-mixed255-full-prefix-$d.json \
 --fine-reference results/early-pole-warm-actual-source-higher-volume-check-$d.json
```
