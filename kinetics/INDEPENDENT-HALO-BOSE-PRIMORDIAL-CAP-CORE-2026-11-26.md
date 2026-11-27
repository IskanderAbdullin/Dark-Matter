# Finite occupation rounds the ideal critical Bose halo center

The [upstream critical-core derivation](https://github.com/IskanderAbdullin/Dark-Matter/blob/f77e75a43e99f03c5ae4dc1d732731bbb9ba3fca/papers/01-critical-core-asymptotics.md)
finds a linear radial density correction for an
exactly critical ideal Bose halo. That endpoint has
unbounded occupation at the center of phase space.
Collisionless evolution transports the initial fine-grained
values, and coarse averaging cannot exceed their upper cap.
This is the weakest part of the more general
[phase-space mixing constraint](https://arxiv.org/abs/astro-ph/0504246).
It does not make a constructed static profile reachable.

Use the spherical isothermal local normalization
of the upstream derivation. For dimensionless potential
difference $y$, kinetic energy $t$, and a cap $F$,
construct the illustrative monotone energy DF
$$
 f_F(t+y)=\frac{1}{\exp[\max(t+y,\delta)]-1},
 \qquad \delta=\ln(1+F^{-1}).
$$
This clips the critical distribution rather than shifting
the chemical potential. Its density integral is
$I_F(y)=\int_0^\infty dt\,\sqrt t\, f_F(t+y)$.
At the center
$I_F'(0)=-\frac12\int_0^\infty dt\,t^{-1/2}f_F(t)$
is finite. The normalized Poisson solution has
$y=\xi^2/6+O(\xi^4)$, where $\xi=r/R$ with the
actual capped central density in $R$. Hence every
finite-$F$ version has quadratic central density,
not the exact critical linear slope. At $y\ge\delta$
the local density integrand is unchanged from the
uncapped gas; this does not make the entire potential
profile identical.

For $F\gg1$, use $w=y/\delta$. In the infrared,
$1/(\exp x-1)\sim1/x$. The inner-to-outer crossover
is explicitly
$$
 \frac{I_F(0)-I_F(\delta w)}{\sqrt{\delta}}\ \to\
 \begin{cases}
 2\sqrt w\,\arcsin \sqrt w+
 2\sqrt{1-w}-\frac23(1-w)^{3/2}-\frac43,&0\le w\le1,\\
 \pi\sqrt w-\frac43,&w\ge1.
 \end{cases}
$$
The first branch starts as $2w$, whereas the second
recovers the critical square-root dependence far
outside the rounding zone. The leading radius
scale $\sqrt{6\delta}\,R$ is a diagnostic, not an
observational radius.

| Occupation cap | Central density / uncapped critical density at fixed velocity dispersion | Numerical rounding radius $y=\delta$, in units of $R$ |
|---:|---:|---:|
| 100 | .94255 | .24743 |
| 5,000 | .99186 | .03470 |
| 7,000 | .99312 | .02932 |
| 200,000 | .99871 | .005479 |

The caps 5,000 and 7,000 are illustrative scales
from the separate [native input-shape ladder](INDEPENDENT-SCALAR-COLD-PACKING-CAP-LADDER-2026-11-26.md),
not measured halo caps. The integral, central
derivative, and Poisson shooting are checked
independently. This is a local single-component
core calculation, not a finite halo fit, a
prediction of collisionless mixing, or a claim
about observed rotation curves. A collisional
mechanism can change the initial fine-grained
occupation bound and needs a different analysis.

~~~sh
python -W error halo_bose_primordial_cap_core_reference.py
~~~

[Capped density, circular-velocity and matching controls](results/halo-bose-primordial-cap-core-2026-11-26.json).
