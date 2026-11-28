# Boosted photon bath: clock convention and the first nonuniversal recoil

Use the [equal-charge photon-only scalar/chiral identity](INDEPENDENT-FINITE-MASS-QED-ANALYTIC-LOW-T-IDENTITY-2026-11-27.md),
now with the parent moving relative to an isotropic bath. Retain
bare photons, keep both fermion masses nonzero, and use the
low-temperature expansion at fixed boost. The hottest photon direction
in the parent rest frame must remain cold relative to the charged
masses and the available decay energy.

For bath temperature T and boost (beta,gamma), a parent-rest
photon direction has thermal occupation n_B(D k/T),
D=gamma(1+beta cos(theta)). The Planck moment of order n
is multiplied by the angular average of D^(-n). In particular,

    <D^(-2)> = 1,      <D^(-4)> = gamma² (1+beta²/3).

The complete angle-integrated real coefficient at order T² is
therefore unchanged in the parent proper clock. The external
phase-space pole shifts are unchanged scalar mass shifts. To verify
the external spin part independently, compute the parent and daughter
spin contractions in the parent rest frame, integrating over two polar
angles and a relative azimuth. The massive daughter logarithm and
parent contraction recover their rest-bath values to 1.2e−14.
Independently averaged Bose integrals and angular moments agree up to
about 9.1e−10 for boost controls through beta=.998.

The thermal mass shifts also change the clock. At fixed parent
spatial momentum in the bath, transform the rate by M_T/E_T and
normalize to the bare clock M/E. With delta(M²)=4h,
h=[g²/(4 pi²)] (pi²/6) T² in physical units, the clock ratio contributes
2h(1/M²−1/E²). The proper-clock relative leading term is
−g² T²/(12 M²); the bath-clock term relative to the bare rate at the
same **momentum** is −g² T²/(12 E²). This is not the derivative
at fixed velocity. At beta=.998 the latter normalized coefficient
in units [g²/(4 pi²)] (pi²/6) (T/M)² is −.007992.

Next remove the complete leading T² term and integrate the exact
signed real recoil against the angle-averaged Bose occupation.
Angular averaging can be performed with just a real logarithm.
Independent cap-aware numerical integrals and high-precision derivatives
of the signed-rate primitive verify the terms through T⁸.
The first higher recoil coefficient c4 multiplies the Planck moment
(pi⁴/15)(T/M)⁴ <D^(-4)> and g²/(4 pi²). On the archive-like scalar
mass fixtures c4 is approximately −4.4008. Its sign is not universal:
with finite mass ratios (m/M,s/M)=(.2,.65) the scalar coefficient
is +19.966, while the chiral coefficient at the same masses is −28.773.
The boost at beta=.998 multiplies the T⁴ Planck moment by about 333.33
at fixed T, and at the same time narrows the domain where the cold
expansion is valid. Neither this enhancement nor the T² result is
a hot-medium prediction for the archived production epoch.

Replay with the existing mass fixtures from kinetics:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_boosted_higher_bath_reference.py
```
[Raw numerical receipt](results/finite-mass-qed-boosted-higher-bath-2026-11-27.json).
