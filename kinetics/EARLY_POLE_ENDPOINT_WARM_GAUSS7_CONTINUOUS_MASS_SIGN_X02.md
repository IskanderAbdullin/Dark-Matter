# Conditional G7/G6 pressure and pseudo-pressure signs on a continuous formal mass interval

Freeze the verified G6 and G7 states at x=.02. Retain each warm cell's
actual count and mean energy at the common native mass. Retain the common
piecewise-constant shape convention for cold cells, allowing their
different measured counts. Test pressure and pseudo-pressure at every
formal observation mass in [0,.1], without associating a mass with a time.

For each warm cell, brackets from the finite linear moment hull are
expanded by the analytic second derivative remainder in reference energy.
Combining terms in the curvature before bounding them prevents spurious
loss of accuracy at small observation mass. In cold cells, subtract the
counts before forming the common kernel bound rather than independently
varying two shared shapes. Both kernels decrease pointwise with observation
mass and increase with q. Brackets at the edges of each adaptive mass
interval therefore bound its entire interior as well.

The 103 certified intervals cover [0,.1] with no gaps. The worst lower
G7/G6 relative differences across all intervals are +3.65855e−4 in
pressure and +9.08720e−7 in pseudo-pressure, expressed as fractions.
These are conservative conditional bounds, not best estimates.
They do not include solver time error, uncertainty in the saved cell
moments, alternative independent cold shapes, or physical dynamics.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results/early-pole-fine-charge2-scalar4; d=2026-11-28
python -W error early_pole_endpoint_warm_g7_continuous_mass_sign.py \
 --native "$p-x02-state-$d.json.gz" --prefix "results/early-pole-warm-gauss7-prefix-x02-$d.json" \
 --states "$p"-gauss{6,7}warm-x02-state-$d.json.gz
```
