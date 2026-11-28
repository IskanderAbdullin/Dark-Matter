# Short expanding FD benchmark with independent entropy and work ledgers

Use the tabulated visible EOS and the closed timelike FD collision operator. A symmetric substep advances the expanding geometry and both temperatures adiabatically, holding comoving dark number and separate bath entropies fixed. The central collision substep holds the scale factor and total energy fixed, changing both bath temperatures with the collision source. The second adiabatic substep closes the split. Hubble expansion is evaluated by an implicit midpoint solve, not imposed as an external clock.

For each adiabatic drift, two separate quadrature orders integrate pressure work along the recovered entropy chart. Across each full trajectory, verify both the change in comoving total energy plus pressure work and the change in comoving total entropy minus the collision ledger. An independent nonsplit DOP853 integration verifies both budgets and the state refinements. Across the three starts the error-halving ratios were 4.0000–4.0735; the largest of the sampled absolute normalized ledger defects was 4.65e-12. The longest simulated log-scale interval was only 0.01187.

Use the same local authors' EOS input as in [the tabulated-bath benchmark](FD_PORTAL_TABULATED_BATH_ENTROPY_AND_MODES.md).

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest -q test_fd_portal_tabulated_entropy_split_expansion
python -W error fd_portal_tabulated_entropy_split_expansion.py eos-range.dat
```

This is a short homogeneous numerical validation, not a full kinetic history, a complete portal interaction set, or a relic-density prediction. The only collision here is the restricted timelike vacuum charged-lepton FD term.
