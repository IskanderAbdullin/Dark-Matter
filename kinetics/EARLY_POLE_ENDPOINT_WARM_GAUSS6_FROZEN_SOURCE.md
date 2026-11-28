# The sixth positive warm order changes the frozen source in the opposite direction again

Use each original S4 state at x=.01, .02, .04 and .05 and lift the **same** constant occupations to warm orders 2 through 6. The charge geometry and occupations, scalar cell counts, mass and collision clock are held. The sixth order has 255 scalar, 448 daughter and 488 parent nodes and remains positive and interior. Its known polynomial moment defect on the live cells is 1.24e-13 through degree eleven in the normalized coordinate.

| Frozen x | Change in total net source G6 minus G5 | Global L1 reference error G6 | Global L1 error of two-cell refinement (219 nodes) |
| ---: | ---: | ---: | ---: |
| .01 | −0.0803330% | 0.0422235% | 0.0161175% |
| .02 | −0.0810320% | 0.0417017% | 0.0170585% |
| .04 | −0.0733460% | 0.0480097% | 0.0249133% |
| .05 | −0.0859898% | 0.0529980% | 0.0370832% |

All comparisons use the corresponding Volume64 total net warm source. In these four frozen snapshots the G6/G5 change is almost all in the forward term, with largest contribution from q=[3.96368,5.72884]. The archived two-cell refinement remains better by this total-normalized L1 measure than changing to sixth order on every warm cell, with fewer nodes. These are fixed-state sensitivity checks, **not** a G6 trajectory, prediction of its pressure sign, or continuum guarantee. The G5-only default still reproduces its prior archive.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results/early-pole-fine-charge2-scalar4; w=results/early-pole-endpoint-warm-source-cells; d=2026-11-28
python -W error -m unittest test_early_pole_endpoint_warm_positive_gauss_ladder
for x in 01 02 04 05; do
 python -W error early_pole_endpoint_warm_frozen_gauss_transition.py \
  "$p"-x$x-state-$d.json.gz "$w"-x$x-$d.json --through 6
done
```
