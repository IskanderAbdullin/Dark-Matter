# Repeated elastic encounters depend on orbit, not just local density

The [upstream F583-4 comparison](https://github.com/IskanderAbdullin/Dark-Matter/blob/f77e75a43e99f03c5ae4dc1d732731bbb9ba3fca/continuation-02/10-collisionality-and-asymptotic-tail.md)
uses a local speed proxy and a single radial passage.
At the same illustrative scales $R=.6$ kpc,
$V_0=44.1\,{\rm km\,s}^{-1}$, and
$\sigma/m=1.25\,{\rm cm^2\,g^{-1}}$, this check
averages the vacuum elastic encounter kernel
along representative bound trajectories.
These are not independently fitted galaxy parameters;
the cross section is a benchmark, not a measurement.

Solve the spherical ideal Bose–Poisson problem
both at criticality and with a smooth central shift
corresponding to peak occupation 7,000. For a fixed
tagged speed $v$ and isotropic target speed $u$,
the average relative speed over the target direction
is $v+u^2/(3v)$ when $v\ge u$, and the
interchanged expression when $u>v$.
Average it over the full local Bose velocity distribution.
For each orbit with given turning radii use
$dt=dr/|v_r|$, where
$v_r^2=2[E-\Phi(r)]-L^2/r^2$.
A cosine change of variable resolves both turning points.
A circular orbit uses its local circular speed directly.

| Orbit (pericenter, apocenter) in $R$ | Vacuum encounters in 10 Gyr, shifted Bose profile |
|---|---:|
| Circular at $R$ | 9.66 |
| Circular at $10R$ | .657 |
| Circular at $20R$ | .138 |
| $(.1,\,10)$ | 3.48 |
| $(.1,\,20)$ | 1.80 |
| $(1,\,20)$ | 1.56 |
| $(5,\,20)$ | .659 |

The time is illustrative. The ratio between very
different orbits reflects residence time and the
actual relative speed; a local clock at apocenter
is not a surrogate for either. The same output
includes independently converged mean squared
velocity kicks and a tagged energy-exchange
moment, for equal masses and an isotropic CM
scattering angle. For the latter, the mean
single-event change of tagged kinetic energy per
mass is $(u^2-v^2)/4$. Central slow trajectories
typically gain energy while the more energetic
eccentric examples lose it in this held bath.

A large encounter count makes a fixed unperturbed
orbit less reliable as an actual trajectory. This is
a first interaction diagnostic under the fixed
potential, not a relaxation eigenvalue, proof of a
thermal halo, or an assembled halo simulation.
Bose enhancement of final states, assembly,
self-consistent redistribution and baryons have
not been included.

~~~sh
python -W error halo_bose_orbital_elastic_clock_reference.py
~~~

[Orbit periods, rate and kick controls](results/halo-bose-orbital-elastic-clock-2026-11-26.json).
