# The two charged cuts share an on-shell angular aperture

The [null photon cut](INDEPENDENT-FINITE-MASS-QED-THERMAL-PHOTON-VERTEX-CUT-2026-11-27.md)
does not settle the medium vertex. This test treats the
principal-value geometry of either one positive-energy charged line
on shell in the bare scalar vertex. It is a kinematic building block,
not a choice of inclusive real-time contour combination. For context,
the organization of one-loop thermal amplitudes into on-shell
forward cuts is described by [Brandt et al.](https://arxiv.org/abs/hep-th/9809177).

Write U for the external momentum of the cut species, V for the
other external momentum, W=U−V for the signed scalar transfer, R
for the cut on-shell loop momentum, and k=U−R. The remaining
charged denominator is exactly D_f=2 W·k. The projected scalar
numerator can thus be split as

    N / (k² D_f) = (C + 4 U·V/k²) / D_f + D / (2 k²),

where C and D are the mass-dependent coefficients in the
implementation. The remaining charged pole is linear in the cosine
to the scalar axis. The photon factor is separately averaged over
azimuth with a rapidity form that preserves the small difference of
near-collinear massive momenta. Integrate the angular pole both by
Cauchy-weight quadrature and by subtracting the root numerator
before adding its known logarithm.

Both charged cuts open in the same mapped parent-energy aperture.
At their simultaneous shell take R_d=R_P−q. They give the same
uncut photon average at paired energies to 7.0e−14 relative,
and the same unsigned energy Jacobian dE/(4 Q). The two principal-
value calculations agree within 7.4e−16 relative at the tested
noncoincident radii; the numerator projection is also checked
by independent spin traces. The orientation of the charged pole
reverses with the signed scalar transfer. Neither this orientation
nor the agreement of unsigned Jacobians determines the sum of
physical cuts with plasma occupations.

The exactly coincident radial point k=0 is not included in these
tests. Negative-energy cut branches, thermal distribution weights,
the external residues and the real-time/soft matching prescription
are not evaluated. No finite medium vertex, gauge-independent
thermal correction, or physical decay rate follows.

Replay with the archived spectra:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_thermal_charged_pv_geometry_reference.py
```
[Raw receipt](results/finite-mass-qed-thermal-charged-pv-geometry-2026-11-27.json).
