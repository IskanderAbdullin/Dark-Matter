# Independent warm scalar Gauss continuation

The fine S4 diagnostic holds the cold scalar and charged grids fixed. Change only the 22 complete warm scalar cells that start above q=.01 from one geometric node to two positive volume-Gauss nodes with independent occupations. The scalar/fermion/parent solver counts are 167/448/488, versus 145/448/488 in the original warm chart. The routine still rechecks all two-mesh pair-number/Pauli/time/continuity gates from the same restart. These are opt-in diagnostic histories, not an alternative official production solver.

The case starts from empty daughters at the same (x) as the matched W1/W2 histories. Use coarse step counts 2048 (through .0075), 512 (through .01), and 1024 for each subsequent segment. Each restart uses the prior *fine* nodal endpoint. The independent receipt checks stored reports and states separately, including all eight log-time checkpoints on both meshes, exact per-native-cell sums, terminal pole, and uncompressed hashes. It does not independently rerun every time step; it verifies the published traces. This is a count/collision ledger, not a proof of time convergence for pressure.

At each comparison epoch, evaluate all spectra at the same numerical mass from W1. For untouched cells retain the original constant occupation convention. In each doubled warm cell use a positive representative with its observed count and mean energy. The independent Jensen/chord interval bounds what any nonnegative within-cell distribution with the same two moments could do to pressure at that fixed mass. It excludes neither clock nor spatial discretization error. The per-cell count/shape decomposition is a finite spectral moment identity, not a reaction-history attribution.

Reproduction from `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_early_pole_endpoint_gauss_warm_chain_receipt
python -W error early_pole_endpoint_scalar_gauss2_history.py fine \
  --charge-div 2 --scalar-div 4 --gauss-floor .01 --steps 2048 \
  --state /tmp/warm-gauss-start-state.json
# For each continuation, add --restart /tmp/warm-gauss-<previous>-state.json,
# --x-end <target>, --steps <coarse> and a distinct --state path.
# The CLI writes each exact checkpoint state, and reports on stdout.
```

The archived reports and nodal states have names beginning `early-pole-fine-charge2-scalar4-gausswarm-` in `results/`. The independent receipt CLI accepts those gzipped archives directly, with `--reports` and `--states` each ordered as start, .01, .02, .04, .05. The comparison CLI takes the matching S4 state, warm Gauss state, and optional W2/W4 states at one explicit `--x`. Do not reinterpret the independent Gauss nodal occupations as piecewise-constant data.

## Resolved snapshot comparison

| x | Gauss pressure vs W1 | Gauss pressure vs W2 | Relative interval at fixed Gauss counts and energies |
| --- | ---: | ---: | ---: |
| .01 | −8.4891% | −2.3707% | 7.96e−10 |
| .02 | −8.4960% | −2.3705% | 7.57e−9 |
| .04 | −8.5415% | −2.3734% | 4.88e−8 |
| .05 | −8.4489% | −2.3614% | 7.98e−8 |

At the terminal snapshot, the total scalar number differs from W1 by only −4.25e−7 relative. The Gauss count/shape ledger accounts for −4.4015 percentage points through changed populations between the unchanged cells and −4.0474 points through changed mean pressure weights within them. The higher-order positive fit preserves both prescribed warm moments to 2.14e−14 or better over the four snapshots. The pressure interval is small even without selecting that particular reconstruction, because at fixed high momentum and positive mass the pressure weight is nearly linear in energy. It does not measure the uncertainty in the evolved nodal moments themselves.

The terminal independent trace receipt gives scalar yield 3.7162728866e−4. The clock difference for this yield using the same .04 restart is 3.1801e−5 relative; the pressure comparison above uses both finer saved snapshots and is not controlled by a count-only clock gate. These same-model results cannot be moved into CLASS, and W4 remains a separate diagnostic.
