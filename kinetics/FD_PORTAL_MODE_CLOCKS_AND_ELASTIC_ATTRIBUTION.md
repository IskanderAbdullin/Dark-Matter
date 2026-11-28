# Separate local mode clocks and the conditional heat effect of elastic scattering

Numerically bracket the slow and fast eigenrates equal to the expansion rate in the tabulated two-bath equilibrium model. In the specified sampled brackets the rates increase through the crossing; both physical quadratures reproduce the four roots. Equality of a *local* linearized rate to the expansion rate is not a decoupling criterion.

| Collision model | Slow rate/H = 1 (GeV) | Fast rate/H = 1 (GeV) |
| --- | ---: | ---: |
| Timelike only | 0.556148 | 0.385182 |
| Timelike + charged-lepton elastic heat | 0.551876 | 0.382409 |

At the combined slow crossing, the fast mode is already about 198.4 times faster than expansion. At the fast crossing, slow/H is only 0.219. The frozen algebraic displacements at these crossings are not small enough to describe a nonlinear state.

The elastic effect is understood in force coordinates. In the symmetric two-by-two timelike rate matrix K, the heat response at zero number flux is `K11 - K01^2/K00`. A positive additional elastic heat conductance adds directly to that conditional value. The product of the two relaxation rates grows by exactly one plus the ratio of this added conductance to the conditional value; individual elastic flavor contributions add in this determinant ratio. The slow-rate ratio is this factor divided by the fast-rate ratio.

At 0.10 GeV the conditional heat is only 9.81% of the diagonal timelike heat rate. Adding elastic conductance equal to 16.23% of the diagonal raises the rate product by a factor 2.655; after dividing by the 1.1215 fast-mode factor, the slow-mode factor is 2.3675. At 0.68825 GeV the conditional heat is 0.750% of the diagonal, so even 0.01508% of the diagonal elastic increment increases the slow rate by about 1.986%. The numerical product identities are checked against both phase-space resolutions.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error fd_portal_local_mode_clock_crossings.py eos-range.dat
python -W error fd_portal_elastic_slow_mode_attribution.py eos-range.dat \
 results/fd-portal-combined-tabulated-modes-2026-11-28.json
```

Use the same locally normalized authors' EOS table as in the closed tabulated bath benchmark. These are local vacuum charged-lepton diagnostics only; no abundance or decoupling prediction is made.
