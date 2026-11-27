# Why averaging the two orders has a dyadic limit

Fix the held elastic response and the 64 successive small elastic
maps. Formally multiply the entire slow field by a scalar
`epsilon`, with all other inputs fixed. At zero, the stationary cold
seed and the hard orbit at every small elastic boundary are identical
in all schedules. Each small elastic map leaves the cold count exactly
unchanged. In the first derivative at zero, cold injection is thus the
duration of each slow pulse times its raw cold RHS on this unperturbed
elastic orbit. Slow perturbations to hard and cold inputs first affect
a later slow injection at higher order in `epsilon`.

In any macro block, let `v0,v1/2,v1` be the cold RHS at its elastic
start, midpoint and end. Per unit macro exposure:

- Slow-elastic-slow gives `(v0+v1)/2`.
- Elastic-slow-elastic gives `v1/2`.
- The arithmetic average, `(v0+2*v1/2+v1)/4`, is exactly the
  composite trapezoid on two half-width macro blocks.

This identity holds separately per mode and per arm, on identical
elastic leaves. It concerns a derivative at zero slow coupling,
not the full slow-coupled trajectories, and does not imply a
second-order error bound or a valid mixed state. It also does not
apply unchanged if the elastic leaves or their graph differ.

## Full nonlinear held comparison

Both arms of each profile start from exactly the same archived hard
and cold data as their counterpart on the other time grids. Pair the
aperture-local daughter input against the native input. In all five
profiles, keep the same 64 small elastic backward solves and graph.

| Macro × elastic leaves | Slow-elastic-slow final mean | Elastic-slow-elastic final mean | Average minus next dyadic original mean |
|---|---:|---:|---:|
| 16 × 4 | +2.12757801e-4 | +4.80464418e-4 | -2.44596e-10 |
| 32 × 2 | +3.46611354e-4 | +4.22530365e-4 | -5.72459e-11 |
| 64 × 1 | +3.84570917e-4 | not tested | |

On the 17 matched recorded times, the peak paired count residual
between the average and the next original grid is 8.1532e-17
and 1.9082e-17, respectively. These are below 1e-6 of the peak
separation of the two coarse orderings, but are still 23.5 and 5.5
stored seed spacings. The individual 32 × 2 orders do not agree:
they use 956 and 508 slow RHS evaluations per arm and straddle
the 64 × 1 finite reference. Averaging two runs is not a cheaper
single admissible method; no runtime or cost advantage is claimed.

## Reproduction

From `kinetics/`, use the paired
`results/early-pole-joined-fine-216-memory-2026-11-27.json.gz`
archive. Pass each of `--case native` and `--case local-daughter`
to the held runner and retain the two arm files in that order.

- Original 16 and 32: run
  `early_pole_local_daughter_elastic_subcycle_reference.py`
  with `--macro-steps 16 --elastic-substeps 4` and
  `--macro-steps 32 --elastic-substeps 2`.
- Alternate 16 and 32: run
  `early_pole_daughter_alternate_outer_order.py` with
  `--macro-steps 16` and `--macro-steps 32`. Raw 32 × 2
  alternate arm ledgers are included as compressed JSON.
- Original 64: run
  `early_pole_prepared_hard_swap_trajectory_reference.py`
  with `--hold-screening --steps 64`.
- Verify the five pairs with
  `early_pole_daughter_dyadic_order_mixture.py` (native then local
  for each profile). The previous 16 comparison output replays
  byte-for-byte with the updated validation.

This is temporal method evidence on one held archived model.
It does not establish accuracy in the hard mesh, under refreshed
screening or for a physical relic prediction.

## Measured tangent on the two actual held elastic orbits

Compute both 64-leaf orbits with the slow stage absent, and
evaluate all eight raw cold-source modes at each of the 65
boundaries, keeping each initial cold mode fixed. The native and
local seed arrays and initial sources matched the accepted
nonlinear prepared snapshots exactly. The maximum absolute
energy defect of an individual elastic leaf was 3.331e-15
in both arms, below the 5e-15 gate. On these data the composite
trapezoid identity replays exactly for each stored block and
each mode of the paired source.

| Macro intervals | Alternate minus original mean: tangent | Full slow field | Full / tangent minus one |
|---|---:|---:|---:|
| 16 | +2.67703192e-4 | +2.67706617e-4 | +12.80 ppm |
| 32 | +7.59182687e-5 | +7.59190107e-5 | +9.77 ppm |

The comparison isolates the order contrast, not the total
contrast between the two physical inputs. For example, the
linear small-leaf trapezoid gives +3.84774506e-4 at 64
macro intervals while the full slow-coupled value is
+3.84570917e-4. Nonlinear dynamics remain necessary for the
latter value. The near agreement here does not provide a
bound for another exposure or regime.

Reproduce the fixed orbit via
`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python early_pole_daughter_fixed_orbit_quadrature.py results/early-pole-joined-fine-216-memory-2026-11-27.json.gz`.
Use `early_pole_daughter_order_linear_limit_compare.py` with
the fixed-orbit and the full dyadic-order JSONs in that order.
