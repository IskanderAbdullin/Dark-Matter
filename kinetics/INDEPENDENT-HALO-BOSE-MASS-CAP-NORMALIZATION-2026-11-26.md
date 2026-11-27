# Phase-space and distance effects on the F583-4 halo mass scale

The [upstream population setup](https://github.com/IskanderAbdullin/Dark-Matter/blob/f77e75a43e99f03c5ae4dc1d732731bbb9ba3fca/papers/04-hierarchical-galaxy-test.md)
gives a one-species ideal Bose normalization
$$
 m=\left[\sqrt{\pi/2}\,
 \frac{\hbar^3}{G\nu R^2 V_0}\right]^{1/4},
$$
where $\nu$ describes the chosen central density
in thermal units.
This is a normalization identity, not a fit or a
measurement of the particle mass.
Use the published-source benchmark $R=.6$ kpc
at 49 Mpc and $V_0=44.1$ km/s, as in the
[F583-4 projection](INDEPENDENT-HALO-F583-4-CORE-PROJECTION-2026-11-26.md).

At fixed angular core size and velocity scale,
the critical scale is 86.20 eV at 49 Mpc and
82.65 eV at the [SPARC catalog distance](https://astroweb.case.edu/SPARC/SPARC_Lelli2016c.mrt)
53.3 Mpc. Merely rescaling that angular core for
the quoted distance endpoints 42.6 and 64.0 Mpc
changes it to 92.45 and 75.43 eV respectively.
These endpoints are not a confidence interval:
$R,V_0,\nu$ are not estimated here. The
stored grid separately illustrates deprojection
at $55^\circ\pm10^\circ$, keeping the same
projected velocity amplitude.

At a fixed illustrative occupation cap, the
maximum central density is different for the
clipped-critical and ordinary smooth subcritical
Bose forms. The latter uses the **same** cap,
unlike the earlier deliberately matched static
comparison with an unconstrained peak.
At $F=7000$ the mass scales increase by only
.17% and .41%, respectively, relative to the
uncapped critical benchmark with the same $R,V_0$.
For the illustrative small cap .512 from the
held scalar production calculation, the increases
are 24.8% and 61.0%. A separate
[covariance-only requirement](INDEPENDENT-HALO-COLLISIONLESS-BOSE-COVARIANCE-CAP-2026-11-26.md)
using the critical model's density **and**
velocity moments, but otherwise arbitrary DF,
would give a minimum mass scale $1.469$ times
the critical benchmark at that low cap.
It is not the same constraint as matching just
the central density at fixed $V_0$.

The [finite lowered controls](INDEPENDENT-HALO-BOSE-FINITE-LOWERED-MIXING-2026-11-26.md)
at $F=7000$ and escape depths $W=2,5,8$
instead give about 88.83, 83.61 and 83.05 eV
at the same catalog distance, compared with
82.65 eV for the untruncated critical form.
Again these are **not** independent fits; freeing
$R,V_0$, changing the baryons, connecting
occupation to mass-dependent production, or
adding compact condensates would change the
inference.

~~~sh
python -W error halo_bose_mass_cap_normalization_reference.py \
 results/halo-f583-4-sparc-table-extract-2026-11-26.json
~~~

[Mass, distance, inclination and phase-space normalization ledger](results/halo-bose-mass-cap-normalization-2026-11-26.json).
