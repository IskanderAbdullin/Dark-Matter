# Can the small finite-occupation rounding be seen at F583-4's sampling?

The [phase-space cap control](INDEPENDENT-HALO-BOSE-PRIMORDIAL-CAP-CORE-2026-11-26.md)
rounds a critical ideal Bose center. This forward check
uses F583-4 radii and velocity-error scales from
[SPARC Table 2](https://astroweb.case.edu/SPARC/MassModels_Lelli2016c.mrt)
and the distance/inclination from
[SPARC Table 1](https://astroweb.case.edu/SPARC/SPARC_Lelli2016c.mrt)
([catalog paper](https://arxiv.org/abs/1606.09251)).
It does **not** fit those velocities.

Keep the illustrative rotation scale $V_0=44.1$ km/s
and use the previously discussed $R=.6$ kpc at
49 Mpc as an angular benchmark, rescaled to
$R=.653$ kpc at SPARC's 53.3 Mpc. Treat this as
the reference critical model. With a finite
occupation cap at the same underlying scalar mass
and $V_0$, the central density changes and so must
the physical Poisson radius $R$; use
$R_F=R\sqrt{I_{\rm crit}(0)/I_F(0)}$.
Two equivalent integration conventions, either with
the corresponding $R_F$ or in one common physical
coordinate, agree to $5.4\times10^{-8}$ km/s at
the sampled radii for the $F=7000$ arm.
Holding both physical $R$ and $V_0$ fixed
instead is a distinct normalization comparison and
is not a fixed-particle-mass change.

At the first public radius, .22 kpc or about
.85 arcsec, the quoted velocity error is 4.9 km/s.
The ideal leading rounding scale for $F=7000$
in this comparison is 19.2 pc or about .074 arcsec.
The [original DensePak study](https://astroweb.case.edu/ssm/papers/ApJv676p920.pdf)
reports 3 arcsec fiber diameter and around 2 arcsec
effective resolution after dithering, not a measurement
of this 19 pc structure. The earlier
[long-slit data](https://astroweb.case.edu/ssm/papers/AJv122p2381.pdf)
also have seeing larger than this scale.

For a synthetic forward control, average the
projected line-of-sight velocity over a uniformly
weighted circular aperture before deprojection,
with an illustrative exponential emission profile.
Add the tabulated gas and .5 times the unit-disk
contribution in squared velocity, equally in both
models. Bracket several tracer radial scales.
This is neither the instrument fiber pipeline nor
a convolution of already reduced catalog velocities.
A solid-body velocity field with uniform light
provides a projection identity check.

For $F=7000$, the maximum fixed-mass model
difference at the sampled radii, even **before**
the aperture, is about $8.0\times 10^{-6}$
km/s in this narrow profile contrast. The sample's
diagonal error-weighted norm is about
$2.3\times10^{-6}$; the illustrative 3 arcsec
diameter aperture decreases it. Holding $R$
fixed without adjusting the central physical
normalization gives a considerably larger .071 km/s
maximum, still below the quoted error scale here,
but it is not the same physics comparison.
Small derivative projections against amplitude,
scale and disk mass-to-light changes illustrate
additional degeneracy rather than fitting data.
No error significance, measured beam response,
or population inference follows from this calculation.

~~~sh
python -W error halo_f583_4_core_projection_reference.py \
 results/halo-f583-4-sparc-table-extract-2026-11-26.json
~~~

[Data extract](results/halo-f583-4-sparc-table-extract-2026-11-26.json)
and [fixed-model projection ledger](results/halo-f583-4-core-projection-2026-11-26.json).
