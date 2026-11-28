# An energy-selected tidal ejection need not reach the tidal sphere

The [tidal-versus-isolated energy ledger](INDEPENDENT-HALO-BOSE-ESCAPE-ENERGY-LEDGER-2026-11-26.md)
conditions on outgoing orbital energy above the selected
cutoff $E_t=\Phi(r_t)$ per unit mass.
At finite angular momentum that condition alone
does not make the radius $r_t$ accessible.
This check follows the subsequent collisionless
one-particle effective potential, not continued
scattering or escape in a time-dependent tidal field.

For any positive spherical density profile,
$$
 \frac{d}{dr}\left(\Phi+\frac{\ell^2}{2r^2}\right)
   =\frac{GM(r)r-\ell^2}{r^3}.
$$
The numerator is increasing, so the effective
potential has at most one minimum. Starting at an
allowed interior point $r$, the outer radius is
accessible precisely when
$E\ge E_t+\ell^2/(2r_t^2)$.
An initially inward particle may reach the boundary
after turning; this criterion is not a flight time.
At fixed outgoing speed $v$, rotational invariance
of the local bath makes its direction isotropic
relative to the radial direction. With
$s=2r_t^2(E-E_t)/(r^2v^2)$, the angular fraction
that cannot reach the boundary is
$\sqrt{1-s}$ for $0<s<1$, and zero for $s\ge1$.
The energy selection already excludes $s\le0$.
This evaluates the angular momentum correction
without approximating an effective barrier by
the potential value alone.

| Escape depth, peak occupation 7,000 | Fraction of energy-selected events blocked by angular momentum | Fraction of radius-accessible events still bound relative to infinity |
|---:|---:|---:|
| 2 | 8.42% | 85.6% |
| 5 | 4.53% | 73.6% |
| 8 | 1.17% | 51.7% |

The latter populations can cross the chosen
sphere but do not have positive isolated energy;
they can return in a vacuum exterior if no
external tides remove them and no further
collision intervenes. The test does not compute
how long they spend outside, actual tidal
stripping, or the evolved halo.
The largest grid/scramble effect on the
angular-momentum blocked rate is under .3%.

~~~sh
python -W error halo_bose_tidal_radius_crossing_reference.py \
 results/halo-bose-escape-energy-ledger-2026-11-26.json
~~~

[Geometric crossing and isolated-energy fractions](results/halo-bose-tidal-radius-crossing-2026-11-26.json).
