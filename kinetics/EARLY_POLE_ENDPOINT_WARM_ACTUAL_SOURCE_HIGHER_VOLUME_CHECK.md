# Higher volume refinement preserves the terminal source ranking

Re-evaluate the same frozen positive count/mean representatives of
the independently verified G7 and mixed255 terminal states at Volume96,
Volume128 and Volume160. Each uses its original charge occupations and
the same native measurement mass and clock. Independently paired
radial checks and the per-cell reconstructed moments remain accepted.

Relative to Volume160, G7 local net L1 residual is 0.0340991% and
mixed255 is 0.0253244%. In both charts the unsigned local difference
between Volume128 and Volume160 is 2.244e−6 of the net reference.
This higher-order comparison retains the ordering found using
Volume64 and Volume96. It does not prove monotone quadrature convergence
or physical time-integrated accuracy.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28; p=results/early-pole-fine-charge2-scalar4
python -W error early_pole_endpoint_warm_actual_source_reconstruction.py \
 --native "$p-x05-state-$d.json.gz" --g7 "$p-gauss7warm-x05-state-$d.json.gz" \
 --mixed "$p-mixed255warm-x05-state-$d.json.gz" --volumes 96 128 160 \
 --control results/early-pole-warm-budgeted-gauss-probe-$d.json \
 --history results/early-pole-warm-mixed255-full-prefix-$d.json
```
