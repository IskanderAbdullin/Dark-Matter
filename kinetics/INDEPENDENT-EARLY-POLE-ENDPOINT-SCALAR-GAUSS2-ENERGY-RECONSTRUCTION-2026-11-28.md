# A positive reconstruction can preserve both moments of a scalar Gauss checkpoint

The [optional independent Gauss-node chart](INDEPENDENT-EARLY-POLE-ENDPOINT-SCALAR-GAUSS2-LIVE-FIXTURE-2026-11-28.md) stores separate nodal occupations, not one constant-in-cell occupation. A later diagnostic may need a positive continuous representative preserving both each cell's count and its mean energy at a specified mass. This does not define a unique physical shape; the [shape-envelope test](INDEPENDENT-EARLY-POLE-ENDPOINT-SCALAR-ENERGY-SHAPE-ENVELOPE-2026-11-28.md) explicitly gives alternatives.

Within one shell set s=(E−E_a)/(E_b−E_a), with s in [0,1], and take the numerical representative `f(q)=exp(α+β s(q))` relative to the phase-space measure q²dq/(2π²). Use the product form `E−E_a=(q−a)(q+a)/(E+E_a)` to resolve even very short shells. The mean s depends monotonically on β since its derivative is the variance under those positive normalized weights. This gives a unique fit for a strictly interior prescribed mean. Fix α from the count after computing β, and use log-sum-exp for stability. An empty stored cell remains empty. This positive exponential tilt is a reconstruction convention, not a physical Bose equilibrium prescription.

The [receipt](results/early-pole-endpoint-scalar-gauss2-energy-reconstruction-2026-11-28.json) uses the exact stored nodal smoke state at x=.005002, reached by a short startup and one short restart with only complete warm scalar cells above .01 on the two-node rule. All 86 nonempty scalar shells match their count and scaled energy to a maximum verification error 1.43e−14 across the independent 64- and 128-node quadratures. Each reconstructed pressure, velocity and inverse-energy moment stays in the separately derived shell envelope. Five synthetic specified tilts are additionally recovered. No physical spectrum or complete Gauss evolution is claimed.

From `kinetics/`, run:
```sh
OPENBLAS_NUM_THREADS=1 python -W error early_pole_endpoint_scalar_gauss2_energy_reconstruction.py \
  results/early-pole-endpoint-scalar-gauss2-energy-reconstruction-smoke-state-2026-11-28.json.gz
```
To replay the smoke state itself use `early_pole_endpoint_scalar_gauss2_history.py fine --gauss-floor .01 --steps 8 --x-end .005001 --state FIRST.json`, then the same chart to `--x-end .005002` with `--restart FIRST.json`. This postprocessing does not change the driver, default chart, or restart acceptance.
