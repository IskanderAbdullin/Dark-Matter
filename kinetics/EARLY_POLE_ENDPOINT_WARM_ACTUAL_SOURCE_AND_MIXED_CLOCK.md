# Warm source checks on the actual G7 and fixed mixed255 terminal states

Hold x=.05, the native S4 observation mass and clock, and each actual
archive's own retained charge occupations fixed. Compare the stored
warm occupations, a positive count/mean-energy reconstruction sampled
at the original nodes, and that same reconstruction evaluated on
32/64/96 volume rules.

| Actual chart | Local net L1 residual vs Volume96 | Signed residual vs Volume96 | Local L1 due to stored occupations vs fit at the same nodes |
|---|---:|---:|---:|
| G7 | 0.0338456% | +0.0210054% | 6.761e−10 as a fraction |
| mixed255 | 0.0250540% | +0.00887665% | 6.759e−10 as a fraction |

The mesh comparison uses paired radial orders for every volume and node
choice. The Volume64 to Volume96 unsigned residual is about 8.47e−6
of the reference in both charts. The two Volume96 net reference rates
agree relatively to 6.96e−8. The result separates quadrature on a
specified representative from changes in stored occupations; the
representative is not claimed to be the actual unresolved shape.

On the original finite node sets, hold each warm cell's count and mean
energy but vary its nonnegative nodal weights. The exact affine
per-cell net source extrema follow by two-node witnesses. These permit
a mixed255/G7 net-source relative contrast between −0.00012138673 and
−0.00012100118, compared with the actual −0.00012119212. Those are
fixed-node instantaneous source bounds, not bounds on possible
continuum shapes or time-integrated trajectories.

Independently recover the actual mixed255 coarse terminal branch from
the exact x=.04 snapshot. Its six archived trace columns match. The
mixed255/G7 fixed-mass pressure contrast is −0.000162255003 on 1024
steps and −0.000162258216 on 2048 steps. The last-segment clock changes
the relative contrast by only −3.21325e−9. This is a distinct observable
from the frozen instantaneous source and gives no full-history error bound.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28; p=results/early-pole-fine-charge2-scalar4
q=results/early-pole-warm-budgeted-gauss-probe-$d.json
r=results/early-pole-warm-mixed255-full-prefix-$d.json
python -W error early_pole_endpoint_warm_actual_source_reconstruction.py \
 --native "$p-x05-state-$d.json.gz" --g7 "$p-gauss7warm-x05-state-$d.json.gz" \
 --mixed "$p-mixed255warm-x05-state-$d.json.gz" --control "$q" --history "$r"
python -W error early_pole_endpoint_warm_actual_source_moment_hull.py \
 --native "$p-x05-state-$d.json.gz" --g7 "$p-gauss7warm-x05-state-$d.json.gz" \
 --mixed "$p-mixed255warm-x05-state-$d.json.gz" --control "$q" --mixed-chain "$r"
python -W error early_pole_endpoint_warm_fixed_budget_terminal_clock.py \
 --native "$p-x05-state-$d.json.gz" --prefix "$r" --control "$q" \
 --restart "$p-mixed255warm-x04-state-$d.json.gz" --report "$p-mixed255warm-x05-$d.json.gz" \
 --final "$p-mixed255warm-x05-state-$d.json.gz" \
 --g7pair results/early-pole-warm-terminal-clock-contrast-gauss6-gauss7-$d.json \
 --coarse /tmp/mixed-coarse-replay.json
```
