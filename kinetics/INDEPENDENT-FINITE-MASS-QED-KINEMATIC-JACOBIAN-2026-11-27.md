# Two independent reductions of the frozen decay geometry

Take the [continuous number-preserving dispersion inputs](INDEPENDENT-FINITE-MASS-QED-SHELL-EDGE-SLOPE-2026-11-27.md)
and include the kinematic Jacobian omitted from the unweighted aperture
test. Retain the frozen Q=4.1e−6 mode, fine and third archived shapes with
their own scalar poles, hard cuts 0.25 and 0.5, and g=0.005 and 0.01.
These are defined kinematic controls with a test amplitude of one.
They are not evaluations of a matched scalar decay amplitude.

Using parent momentum P, daughter momentum k, and fixed scalar
momentum Q, eliminating the angle gives a common factor 1/Q times
the two-momentum measure Pk dP dk. Use the energy delta function to
eliminate k: the remaining weight is
`Pk/(E_P^bare E_k^bare v_k) dP`, where v_k is the daughter group
velocity from the specified candidate energy. Independently eliminate P
and use `Pk/(E_P^bare E_k^bare v_P) dk`. In each chart invert the
energy condition separately and split at both sets of occupation edges
and smoothing joins. Replacing the bare denominators with the modeled
energies is an explicit additional convention check, not a residue
calculation. Both reductions also reproduce the known bare parent
energy integral, where k/(E_k v_k)=1.

At g=0.005 the real finite-mass logarithm changes the weighted
stimulated integral by −0.00459% to −0.01158% across the four
native grid/cut pairs. The change in the relative response compared
with the flat aperture is below 6.3e−10. At g=0.01 the largest
weighted reduction is 0.02271%, and the flat/weighted response
difference remains below 3.2e−9. Replacing just the denominator
convention changes the base stimulated integral by less than 9.0e−6
relative. The two chart reductions agree to 1.1e−16 absolute and
the vacuum control recovers the bare result to 1.4e−17. This verifies only that
the derivative and kinematic geometry do not amplify the small real
log effect in these frozen inputs. It does not constrain dielectric
resummation, widths, balance conditions for the full medium,
wavefunction residues, amplitudes, or evolved production.

Replay from kinetics:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_kinematic_jacobian_reference.py
```
[Full receipt](results/finite-mass-qed-kinematic-jacobian-2026-11-27.json).
