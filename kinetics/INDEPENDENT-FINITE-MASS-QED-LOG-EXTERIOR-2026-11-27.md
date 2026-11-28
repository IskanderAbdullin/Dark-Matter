# The frozen exterior contributes to the active log dispersion

The [archive pointwise check](INDEPENDENT-FINITE-MASS-QED-ARCHIVED-HARD-LOG-2026-11-27.md)
included all represented hard internal momenta. A finite active
collision domain may evolve only a subset of these momenta. Cutting
the same-species log interaction at its active boundary would
therefore discard not just exterior energy but a term in the
**active** energy gradient.

In the symmetric shell kernel K and its fixed-comoving derivative D,
split each archived occupation into disjoint x (active) and y
(exterior). The energy splits as
`xᵀKx + 2 xᵀKy + yᵀKy` with the same common prefactor.
The pure exterior term contributes to the scale derivative even
at fixed x; the cross term contributes to both scale work and
`2 Ky` in the unnormalized active first variation. An independently
differenced cross energy and the direct gradient agree to rounding.
Insert the hard floor and all active limits as exact subcell cuts
without altering the archive occupations.

Use both x=0.05 archive grids, active momenta [0.45,U] for
U=2.2, 5, 10, and exterior hard floors 0.1, 0.25, 0.45.
Probe the gradient on subcells in the decay band [1.40,1.46].
For parent fermions with hard floor 0.25 and U=2.2, the cross term
alone contributes about 31.2% of the represented parent log trace.
The exterior supplies 20.35–20.44% of the active log gradient probe
across both grids. With U=5 the exterior gradient fraction drops
to 5.83–6.08%, and with U=10 to 5.38–5.63%. The remaining floor
band [0.25,0.45] accounts for nearly all of that widest-domain
remainder: matching the hard floor to the active lower edge gives
only about 0.0017% for the widest domain. The daughter's exterior
gradient fraction is below 3e−6 in the U=2.2, hard-floor 0.25 test.

These fractions are relative to the small real log contribution,
not to the full energy or decay source. The frozen exterior has no
reaction flux or time-dependent occupation. This test does not
replace a live active/exterior variation, prove momentum convergence
for the physical trajectory, or extrapolate into the unmodeled
soft domain.

Replay from kinetics using the earlier fine and third archive snapshots:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_log_exterior_reference.py
```
[Raw receipt](results/finite-mass-qed-log-exterior-2026-11-27.json).
