# Warm-only refinement changes the evolved pressure while preserving the cold count

Repeat the fine S4 empty-start history with the same cold and charge charts, and split only complete native warm cells with lower edge ≥.01 into two subcells. The alternative has 167 scalar, 448 daughter and 488 parent nodes. Both use a 4096-step retained startup and identical later time-mesh counts. They are independent trajectories, not a static remap.

| x | Total count relative change | Reconstructed pressure-numerator relative change | Nodal pressure-numerator relative change | Conditional cold W₁ |
| --- | ---: | ---: | ---: | ---: |
| .01 | -2.7146e-06 | -6.2670% | -3.5460% | 1.667e-14 |
| .02 | -8.3157e-07 | -6.2743% | -3.5538% | 1.678e-14 |
| .04 | -3.6981e-07 | -6.3180% | -3.5978% | 1.730e-14 |
| .05 | -3.0948e-07 | -6.2348% | -3.5116% | 1.808e-14 |

At x=.05, the normalized warm W₁ within q≥.01 is 0.0329463. At the same fixed mass in this comparison, the geometric-node pressure bias is -3.7717% on the original warm mesh and -0.9770% on the refined one. The geometric-node comparison therefore combines the change in the actual computed distribution with a representation change. The earlier [fixed-shape fixture](INDEPENDENT-EARLY-POLE-ENDPOINT-WARM-NODAL-FIXTURE-2026-11-28.md) isolated only the latter. This is a conditional numerical pressure functional of a finite retained histogram, not a cosmological pressure or transfer prediction.

The [full warm chain](results/early-pole-fine-charge2-scalar4-warm2-chain-2026-11-28.json) passes the saved state and budget tests. The final fixed-start temporal difference is 3.18e-05 for total yield; this is not a measured pressure time error. The [full audit](results/early-pole-fine-charge2-scalar4-warm2-evolution-audit-2026-11-28.json) preserves exact and nodal moment comparisons. Replay with the scalar subshell driver at fixed `--charge-div 2 --scalar-div 4 --warm-floor .01 --warm-div 2`, using coarse step counts 2048 at startup, then 512, 1024, 1024, 1024 with exact saved-state restarts. The fixture does not establish warm-grid convergence or permit CLASS use.
