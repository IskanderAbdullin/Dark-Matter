# Finite direct history with timelike and elastic FD rates

Start both baths at 0.90 GeV and zero chemical force, and evolve to log(a/a_start)=0.25 with no collision-rate interpolation. BDF and Radau use distinct direct phase-space resolutions. Verify the independent comoving number, total energy/pressure work, and total entropy budgets. The full reported coordinate discrepancy is below 3.38e-9 between these integrations, and the largest reported normalized budget defect is 7.32e-11.

At the end, visible temperature is about 0.71045 GeV. Adding elastic heat changes the terminal dark chemical parameter by +5.3313e-5 and the dark log-temperature by -7.113e-6. The final comoving number difference divided by the initial number is +2.8340e-5. Even though elastic scattering produces nonnegative entropy on a given instantaneous state, the combined solution's total integrated entropy production in this test is lower by 3.0495e-10 of initial entropy than in the separate timelike-only solution. Both cumulative entropy productions are positive. The two paths have different temperatures and chemical parameters; this is not a negative elastic entropy source.

```sh
cd kinetics
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error \
 fd_portal_direct_combined_expansion_tracking.py eos-range.dat
```

Use the same locally normalized tabulated EOS input as in the closed tabulated-bath benchmark. This is a finite, homogeneous, restricted vacuum charged-lepton history, not a complete reaction network or a prediction of relic abundance.
