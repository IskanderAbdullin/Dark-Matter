# Fixed scalar momentum changes the sign of the leading-subtracted recoil test

The [plasma-bath control](INDEPENDENT-FINITE-MASS-QED-BOOSTED-SOFT-RECOIL-2026-11-27.md)
integrated all scalar angles and did not retain the cold selection.
Here hold the scalar momentum at Q=4.1e−6 and sample three
interior parent energies (20%, 50%, 80% of the appropriate
vacuum Born window) in each archive. Use the archived
parent and daughter Fermi occupations and an illustrative
isotropic plasma photon distribution with T=0.9.

For signed plasma photon energy k and angle z to the parent,
let ell=E_P−P z and p_R²=(P−k)²+2k(ell−ell_forward).
The on-shell daughter condition fixes the scalar angle to the
recoiling pair by

    cos(theta) = [A0 + k ell]/(Q p_R),
    A0 = (E_P−k) E_s − (M²+m_s²−m²)/2.

Solving |cos(theta)|≤1 in the signed coordinate k ell gives
a quadratic interval which can be clipped exactly to physical
photon angles. Integrate the exact spin-summed real amplitude
over the remaining azimuth, and use the corresponding 1/p_R
phase-space weight. An independent chart instead eliminates a
photon angle at fixed scalar direction and uses 1/|P−S|.
The two energy-delta reductions agree to 1.2e−14 absolute.
Resolving the daughter collinear region requires a log-distance
coordinate with the daughter mass scale, not just a log of
the parent photon denominator.

Multiply positive and negative real branches by their own
forward/reverse quantum factors, and subtract only the universal
leading-soft virtual term on the common Born event. For k≤1e−7
evaluate the even and odd differences through high-precision
local Taylor coefficients before the Bose factor. Above that
use extended-precision signed evaluations on exact
kinematic and occupation panels. Re-estimating with a higher
angle and energy quadrature changes the integrated residual
by at most 1.6e−9 absolute. Removing the final unresolved
infrared region below 1e−10 changes it by under 9e−11.

At fixed plasma photon cap 5e−6 the normalized partial
residual (relative to the local Born stimulus, **before**
the common g²/(4 pi²) factor) ranges from −0.00986 to −0.1193.
It is negative at every tested parent point, unlike the positive
all-scalar-angle recoil check. This shows that the scalar
selection changes this restricted signed contribution, not the
sign of physical production. In particular it excludes many
photon directions starting near 0.8e−6 at the 20% and 80%
points and near 2.0e−6 at the midpoint. Parent energies
outside the Born window can acquire real support and were
not integrated here. No finite virtual term, sum over the
full parent spectrum, hard/soft matching, in-medium width,
or physical inclusive rate is claimed.

Replay from kinetics:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_fixed_q_signed_soft_reference.py
```
[Raw receipt](results/finite-mass-qed-fixed-q-signed-soft-2026-11-27.json).
