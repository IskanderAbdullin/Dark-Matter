# Real finite-mass exchange in the frozen occupied aperture

Extend the [frozen equal-momentum check](INDEPENDENT-FINITE-MASS-QED-ARCHIVED-HARD-LOG-2026-11-27.md)
to the actual unequal parent and daughter momenta allowed by energy
and momentum conservation. Use the same archived fine and third endpoint
spectra at x=0.05, both frozen scalar poles crossed with both archive
grids, internal hard lower cuts 0.1, 0.25, 0.5, and the scalar momentum
Q=4.1e−6. Photon moment is a common thermal comparator at the selected
cut. Keep the distributions entering the energy correction piecewise
constant. For the occupancy factors in the aperture test separately
repeat using limited number-preserving reconstruction in p³/3.

For collinear faces, enforce the real energy difference of parent
momentum P and daughter momentum P±Q, evaluating the bare difference
in rationalized form. Add the finite-mass dispersion as an additive
order-g² energy term; solve the endpoints and separately invert the
daughter energy inside. Split integrals at parent edges and the
images of daughter edges under the equal-energy condition. On each
piece audit the constant occupancies at three distinct points. For the
limited reconstruction compare orders four and eight (absolute error
≤1.1e−18). Report only unweighted parent-momentum integrals of
f_parent(1−f_daughter) and f_parent−f_daughter. These are neither
invariant phase-space integrals nor physical rates: no dispersion
Jacobian, vertex, residues, widths, or screening are matched.

An independent vacuum closed-form check fixes the collinear roots,
and the exact effective-squared-mass comparator for the factorized
dispersion has the expected order-g⁴ difference from the additive
energy expansion. The endpoint change caused by the logarithm,
computed by subtraction of independently found roots, agrees with
the analytic log-only tangent within 1.76e−5 relative at g=0.001.
The interior energy defect is below 1.74e−16 absolute.

At g=0.005, the factorized dispersion moves the window center by
6.62–7.22% of the bare width. The added log correction moves each
face back by only about 0.0064–0.0072% of that width and changes the
width itself by less than 0.00005% relative. The stimulated integral
falls by at most 0.0142% across all grid/cut pairs, with similar
worst-case shifts for constant and limited occupations. At g=0.01
its largest fall is 0.0396%. None of the 36 tested frozen windows
closes or yields a negative integrated stimulus. The observed
response quantifies a specific kinematic control, not the uncertainty
in the actual finite-mass Boltzmann evolution.

Replay against the older archived snapshots in the results directory:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_archived_aperture_reference.py
```
[Full receipt](results/finite-mass-qed-archived-aperture-2026-11-27.json.gz).
