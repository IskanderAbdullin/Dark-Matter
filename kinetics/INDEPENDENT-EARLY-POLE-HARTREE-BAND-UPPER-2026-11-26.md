# Conditional production upper for any shape inside the early cold band

The [thin-shell feedback](INDEPENDENT-EARLY-POLE-SHELL-HARTREE-FEEDBACK-2026-11-26.md)
assumes uniform occupation in one shell. This
comparison keeps its positive Hartree mass and
narrow pole kinematics but removes the uniform
shell shape. All quantities below are in local
temperature units, as $T\propto a^{-1}$ with
constant visible entropy degrees of freedom.
There is no scattering or other source placing
new particles in the fixed $p/T<P_c$ band.
An arbitrarily shaped population elsewhere may
add positive mass feedback, but not reduce it.

Let $n$ be the band count in $T^3$ units.
For any distribution inside this band the positive
tadpole implies
$m^2\sqrt{m^2+P_c^2}\ge\lambda n/2$.
Choose $a>0$. After any first upward crossing
of $n_*=(2/\lambda)(aP_c)^2\sqrt{(aP_c)^2+P_c^2}$,
the mass must satisfy $m\ge aP_c$ whenever $n\ge n_*$.
The finite parent-energy window has length $2pD/m^2$.
Grant the entire interval the maximum FD parent
occupation $1/2$ and remove all daughter blocking.
The per-mode gain coefficient is bounded by
$g_P\Gamma M^3/(2 E m^2)$.
Using the tadpole inequality gives
$$
 \frac{d n}{d\ln a_{\rm FRW}}\le
 \frac{\kappa g_P B\sqrt{1+a^{-2}}}{\lambda}
 x^5[1+{\cal V}/n_*],\qquad x=M/T,
 \quad {\cal V}=P_c^3/(6\pi^2),
$$
where $\Gamma/H\le B\kappa x^2$.
The rate bounds spontaneous plus stimulated gain
for arbitrary scalar occupation in the band.
If the count never crosses $n_*$, that threshold
already bounds it. Otherwise integrate this upper
rate as if the crossing were at the initial $x_0$:
the excess is bounded by the constant prefactor
times $(x^5-x_0^5)/5$.
Later inverse loss cannot increase this upper.

For the illustrated scalar Yukawa rest width,
$B=1$ suffices at massless fermion daughter.
In general put $r=b/M$; at any nonnegative
scalar mass the width ratio is at most
$(1-r)(1+r)^3\le27/16$.
This follows from monotonically reducing the
rest momentum and spin sum with scalar mass,
and maximizing over $0\le r\le1$.
Both choices are tracked in the archive.

At $P_c=10^{-5}$, $x=.05$, and allowing
the generous $B=27/16$, the upper is
$Y_{\rm band}\le .00256$ for
$\lambda=7.65\times10^{-8}$, and
$Y_{\rm band}\le .000181$ for
$\lambda=1.085\times10^{-6}$.
These are 55.1% and 3.88% of the chosen
illustrative narrow target, respectively.
In the first case even the generous comparison
requires $x\gtrsim.0563$ before the bound
could reach that target, not that it will reach it.
For the smaller width envelope $B=1$ the
first upper drops to .00152.

This is a band-limited, time-limited bound **within
the specified Hartree pole model**. It is not an
upper on all scalar production or on arbitrary
dispersion, coherent dynamics, enhanced parent baths
or processes transferring already-produced particles
into the comoving band. The illustrative $P_c$
and count target do not turn it into a structure
formation or actual reheating constraint.

~~~sh
python -W error early_pole_hartree_band_upper_reference.py
~~~

[Finite band and width-envelope controls](results/early-pole-hartree-band-upper-2026-11-26.json).
