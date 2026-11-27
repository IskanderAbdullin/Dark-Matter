# Resolve the charged-cut principal value near the overlapping soft point

The [charged-cut geometry](INDEPENDENT-FINITE-MASS-QED-THERMAL-CHARGED-PV-GEOMETRY-2026-11-27.md)
compares the angular principal value at well separated radii.
Near equality of the internal and external charged energies,
the dressed interpretation is not known, but the **bare**
angular integral admits a local analytic subtraction.

Let p,r be the charged external and internal radial momenta,
c the external direction cosine relative to the scalar, and
u the internal cosine. Set d0=2m_cut² sinh²((eta_p−eta_r)/2),
where eta_p and eta_r are the charged rapidities. The
exact azimuthal discriminant of the uncut photon is

    X = (p r)² (u−c)² + 2 d0 p r (1−c u) + d0².

Subtract the analytic narrow profile
`Y=(p r)²((u−c)²+alpha²)`,
`alpha²=2 d0(1−c²)/(p r)`. The square-root inverse difference
is evaluated by rationalization with Y−X rather than direct
subtraction. Integrate the angular principal value of the
model exactly using u−c=alpha sinh(y), then integrate only the
remaining regularized kernel numerically.

For internal energy shift δ, both the separation of the angular
pole from c and alpha are linear in |δ|. In the odd part of the
angular principal value this gives a bare 1/δ term. The analytic
coefficient is reproduced at δ=1e−8 within 3.9e−10 relative.
It is approximately +13.27/−17.24 (parent/daughter)
on fine and +13.20/−17.14 on third in the declared units.
Thus the two unequal-mass bare coefficients are not exact
opposites. For δ≥1e−5 the subtraction reproduces direct
resolved angular integration. The symmetric even remainder
varies logarithmically with |δ|; the change over the last
decade agrees within 0.62% with the coefficient from the
azimuthal photon integral without the charged denominator.
This logarithm is integrable with respect to radial energy.

This isolates the singular terms of the **bare mathematical
angular integral**. It does not select a physical radial
principal-value prescription or supply the occupation factors,
negative-energy cuts, residues or real contributions needed
for a matched thermal calculation.

Replay with the archived spectra:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_charged_overlap_asymptotic_reference.py
```
[Raw receipt](results/finite-mass-qed-charged-overlap-asymptotic-2026-11-27.json).
