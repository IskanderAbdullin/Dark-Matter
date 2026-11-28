# Fold one changed sub-shell profile through the frozen scalar pair kernel

A [moment feasibility test](INDEPENDENT-WARM-SHELL-OCCUPATION-FEASIBILITY-2026-11-26.md)
shows that the same retained light number and kinetic energy per
shell allow a dilute profile and profiles with tiny (f=1)
intervals. Test whether one such change can matter to the
instantaneous (chiarchileftrightarrowphiphi) fold at a
fixed saved q32 snapshot (visible temperature 10.391 MeV).

Replace only one light fermion shell in each control. Take the
same occupation for both charge sectors. Solve exactly in radial
volume for its number, and slide the full interval to match
kinetic energy. Run the same forward and reverse event integrals
on both profiles with the same refined radial and final-angle
cuts, including both interval edges. The scalar profile is
unchanged. No redistribution occurs outside the one chosen shell.

The shell used in the original cap witness holds only
(1.22	imes10^{-5}) of retained light number. Saturating a slice
of it makes no resolved difference to the pair fold at the
precision of the quadrature ladder. Instead saturating a slice
of shell 15, which holds 32.55% of retained light number, changes
the stimulated forward rate by about +0.247% and the signed net
scalar pair rate by about +1.07% on the fine quadrature. The
forward shift persists over radial/polar/azimuth/final-panel
orders (3,8,8,3), (4,10,10,4) and (6,14,12,5). The inverse
rate is much less affected, so equal number and kinetic energy
neither force equal rates nor make a sharp spike in an arbitrary
shell important to this integrated kernel.

These are separate mathematical perturbations of one frozen
reconstruction, not dynamical states, extremal rate bounds,
or evidence of a sustained percent-level yield shift. Other
operators and unresolved exterior distributions are untouched.

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error warm_shell_saturated_pair_fold_reference.py
~~~

[Diagnostics](results/warm-shell-saturated-pair-fold-2026-11-26.json).
