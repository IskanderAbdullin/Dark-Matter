# A gravitationally bound tidal crossing can remain outside for many dynamical times

The [angular-momentum check](INDEPENDENT-HALO-BOSE-TIDAL-RADIUS-CROSSING-2026-11-26.md)
finds many elastic products that can cross the tidal
radius despite having negative isolated orbital
energy. They cannot escape an isolated vacuum halo
to infinity, but they do not all return quickly.
Here compute only their **exterior** excursion, from
outward boundary crossing to its next inward crossing.
The time inside the halo before departure is not
included, nor are outside mass, tides or scattering.

In the vacuum exterior the Kepler parameter is
$\mu=GM$. For binding energy per mass $\epsilon>0$
and angular momentum per mass $\ell$, the semimajor
axis and eccentricity are $A=\mu/(2\epsilon)$ and
$e=(1-2\epsilon\ell^2/\mu^2)^{1/2}$. Let the
eccentric anomaly at crossing be
$u=\arccos[(1-r_t/A)/e]$. The exterior duration
is
$$
 t_{\rm ext}=2(\pi-u+e\sin u)\sqrt{A^3/\mu}.
$$
A direct radial quadrature provides an independent
check. The angular distribution of the high-energy
product is conditioned on reaching the boundary.
For a given time cut, use a restricted outgoing energy
band from the necessary condition that the **whole**
orbital period exceed the time cut. This improves
sampling of rare near-zero binding energies without
assigning all weakly bound particles the same period.

The output conditions the quoted duration fractions
on both negative isolated energy and geometric
boundary accessibility. Independent number-flux
checks agree with the prior angular integration.
A long excursion is not a lifetime or an actual
debris mass: the halo was held fixed and no
external perturbation or rescattering along the
flight is included.

| Escape depth | Fraction outside >1 Gyr | >10 Gyr | >100 Gyr |
|---:|---:|---:|---:|
| 2 | 51.1% | 8.39% | 1.59% |
| 5 | 85.2% | 22.3% | 4.25% |
| 8 | 99.0% | 76.8% | 19.9% |

The two finest independently scrambled arms differ
by less than .45% relative across these thresholds.
The coarse-to-fine difference is at most 1.72%.
A fixed potential over very long illustrative clocks
is not a physical extrapolation of the galaxy age.

~~~sh
python -W error halo_bose_exterior_return_clock_reference.py \
 results/halo-bose-tidal-radius-crossing-2026-11-26.json
~~~

[Ballistic exterior duration checks](results/halo-bose-exterior-return-clock-2026-11-26.json).
