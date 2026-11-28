# Sensitivity to moving the source-free switch from 12 to 9.3 MeV

The 50 MeV evolving visible/dark source calculation is continued to
9.3 MeV with both vacuum and signed state-dependent portals. Each
exported FD/BE shell is evolved without a thermal restart. The 12 MeV
controls switch off the source earlier. In both groups the subsequent
evolution uses tabulated radiation plus retained dark gravity and
the same shellwise below-pole and retained on-shell model.

Use the measured scale factors of the archived endpoints. Let
C=a_9.3/a_12. The 12 MeV run first advances source-free to b=C
in eight log intervals, then uses exactly the same absolute scale
factor mesh as the 9.3 MeV start: 36 intervals from a_9.3 to
410 a_12. The on-shell kinematic endpoint is also placed on the
mesh in each run. This controls for a strong time discretization
effect in the earlier unmatched-mesh comparison.

The switchover comparison includes residual portal production,
different prescription for visible backreaction between 12 and
9.3 MeV, and different pressure and collision scheduling over that
short precursor. It is a sensitivity test, not an isolated analytic
derivative of one of those changes. Below 9.3 MeV both runs omit
additional source and inverse portal reactions; the independent
small forward-source estimate is not a late inverse bound.

## Matched late result

For the quantum-source test the scalar N/N_sat changes from
0.00233819 (12 MeV switch) to 0.00234045 (9.3 MeV switch), a
0.096% relative shift. The scalar count fraction below q/T_R=0.005
changes from 4.13858e-8 to 4.13771e-8. The final equivalent
temperatures are 13.632/0.4013 keV versus 13.626/0.4012 keV.
Vacuum-source comparisons give nearly identical endpoint metrics.
This is much smaller than the earlier change from a thermal
reinitialization of the microscopic distributions; it does not
examine a substantially wider late momentum support.

Use the [benchmark EOS setup](INDEPENDENT-BENCHMARK-SETUP.md).
The early snapshots are reproducible with the preceding 12 MeV
source command replacing `--empty-stop .24` by `--empty-stop .186`.
Retrieve C from the ratio of `final_a` in the two early results.
For 12 MeV use `--steps 36 --preroll-stop C --preroll-steps 8
--stop 410`; for 9.3 MeV use `--steps 36 --stop (410/C)`
with the numerical value of 410/C. Both use `--prompt-until 8
--radiation-table EOS.late.dat --dark-gravity` and their matched
archived entropy-variable snapshots.
