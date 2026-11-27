# Separate the initial Hubble-rate rescaling from the redshift path

While following the exact same log-scale history as the [joined expanding reference](INDEPENDENT-EARLY-POLE-JOINED-SCREENED-FRW-2026-11-26.md), multiply **all** collision rates by a single constant, chosen from the initial dark-to-visible energy ratio. The trace integral is not rescaled. Unlike rerunning a shorter expanding history, this control changes only the collision clock, not the redshift traversed by the hard and exterior particles.

On the fine five-step checkpoint fixture at exposure 2e−6, the initial energy ratio is .05135770288767686. Its constant clock factor is .9752697401968827. The archived dynamic full-H gain is 1.486745458481827e−7. Holding the initial clock factor gives 1.486745458438736e−7, a difference of −4.31e−18 in the signed gain. The two independently evaluated entropy and trace-subtracted energy ledgers retain the existing gates. This is a short reference comparison, not a statement that the conditional full-H clock can always be held constant.

## Recovering near-tolerance midpoint solves without loosening the gate

In the longer 96-step controls, a few full elastic midpoint solves reached a line-search floor just above the acceptance gate. The solver now detects that specific failure, retries the **elastic** midpoint operator as two consecutive half steps at the same scale and clock, and records any retries plus their entropy contributions. Other runtime errors propagate. The enclosing slow half steps and the atomic checkpoint granularity are unchanged; no energy, positivity or entropy gate is relaxed. The recovered fixed trajectory crossed macro step 70 and the expanding full-H trajectory crossed steps 80 and 83 after restarting from the last complete checkpoint. The full long comparisons and size of the local discretization difference are separate checks.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1
python -W error early_pole_joined_screened_frw_reference.py fine \
  --duration 2e-6 --steps 5 --expands --full-H
python -W error early_pole_joined_screened_frw_reference.py fine \
  --duration 2e-6 --steps 5 --expands --clock-scale .9752697401968827
python -W error -c 'from early_pole_joined_screened_frw_reference import retry_fixture; retry_fixture()'
```
Short outputs and accepted recovery diagnostics are in `results/early-pole-joined-frw-clock-and-recovery-2026-11-27.json.gz`.
