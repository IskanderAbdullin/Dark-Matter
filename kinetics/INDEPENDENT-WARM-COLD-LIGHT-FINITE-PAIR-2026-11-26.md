# The native scalar-pair quadrature misses reverse gain in the coldest light cell

Freeze the two [accepted 8.133100 MeV endpoints](INDEPENDENT-WARM-RECOIL-CROSSING-MATCHED-2026-11-26.md).
Compute the same scalar pair collision from the perspective of an outgoing
tagged light fermion, rather than counting hits to its first cell after
quadrature over two incident scalars. Use exact finite tagged momentum each
time and average the gain and loss with p²dp over the whole first shell
(upper edge 36.952 keV in momentum units). The partner is fully integrated.
The Pauli and Bose weights, invariant flux, identical scalar factor and
full azimuth-dependent matrix element are kept. The scalar lab-angle
integral splits at both outgoing scalar support/shell boundaries.
Contributions outside the retained baths use zero occupation, not a
thermal fill; spontaneous forward production into an empty scalar state
is allowed.

At q64 the converged independent first-cell gain is 1.40777e-5 and loss
2.53130e-5, giving net −1.12353e-5 per Hubble time. At exactly rest the
net is −1.10603e-5; finite cell averaging changes its magnitude by only
1.58%. The q32 shell yields −1.12362e-5. Refining the joint primary
quadrature changes the reported net by less than 1e-8 fractionally on
both supports. At three thermal tagged momenta, common temperature and
chemical potential give pointwise detailed balance to 2.1e-14 and
integrated balance to 5e-16. The p=0 evaluation also independently
reproduces the existing analytic rest-angle implementation to 1.4e-9
in both separate gross components.

At the same states the native loss agrees within 9e-7 relative. But its
reverse **gain** is about 2.0e-10 at the evolution quadrature and 7.40e-6
at the refined native quadrature, compared with 1.41e-5 independently.
The native signed pair loss is therefore 2.25× too large in magnitude
on the evolution quadrature and 1.59× too large on the refined one.
The discrepancy is not explained by finite momentum or the q32 prompt
support cutoff. Number/energy conservation checks of the native event
projection do not test whether it assigns enough reverse events to this
very small cell.

This is a frozen result. It does not correct the accepted FRW trajectories
or establish a corrected global pair rate. The natural next test is a
conservative partition of reverse phase space: use final-light variables
for the near-rest target and keep the native rule on its complement,
updating all four reactant/product moments for every event.

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error warm_cold_light_finite_pair_reference.py EOS.dat
~~~

[Detailed controls and gross components](results/warm-cold-light-finite-pair-2026-11-26.json).
