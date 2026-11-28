# The negative-energy charged cuts have no bare pole overlap

Complement the [positive charged-cut test](INDEPENDENT-FINITE-MASS-QED-THERMAL-CHARGED-PV-GEOMETRY-2026-11-27.md)
with the regular negative-frequency sector. For each species put
the internal charged four-vector R=(-E_r,r) on its past mass shell;
the external charged leg U is future on shell. The uncut photon
has k=U−R and satisfies k²≥4m_cut². For the two archived mass sets,
the remaining charged denominators also stay bounded away from
zero: the lower bound on their magnitudes is about 0.00258 on the
parent cut and 0.00250 on the daughter cut. These follow from
the signed timelike scalar transfer and hold for all loop angles
in this fixture, not just sampled radii.

Use the explicit projection from the positive cut, but its
azimuthal photon average now has positive sign. Resolve the sharp
antiforward region for the lighter fermion with an analytic
angular substitution. The angle integrals agree with independent
adaptive quadrature to 4.0e−14 relative. At fixed frozen archive
occupations, the radial quadratures agree within 5.0e−10 relative.
Under the declared common normalization, prior to any real-time
cut coefficient or statistical/coupling prefactor, the two archived
parent integrals are about +264 to +266 and the daughter integrals
about −34.0 to −33.6. For the antiparticle occupation we explicitly
use the same stored shape as for the particle in this control.
About 33% of the parent value versus 70–73% of the daughter value
comes from cells intersecting momentum [1,2]. This fraction
describes the represented cell grouping, not a cut on a continuous
distribution.

Regularity of these denominator factors is not a claim about
the complete thermal vertex. No thermal contour signs, external
residues, photon spectral resummation or inclusive rates are
specified or inferred. The positive-frequency overlap remains a
separate matching problem.

Replay with the two archived spectra:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_negative_charged_cut_reference.py
```
[Raw receipt](results/finite-mass-qed-negative-charged-cut-2026-11-27.json).
