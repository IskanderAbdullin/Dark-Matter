# Coldest ideal thermal excited gas at fixed scalar abundance

The [chemical-equilibrium comparison](INDEPENDENT-SCALAR-CHEMICAL-EQUILIBRIUM-COLD-ORDER-2026-11-26.md)
fixes the scalar chemical potential $\mu=0$.
This extension asks whether kinetic equilibrium after
number-changing reactions shut off can evade its warmth.
Use ideal Bose statistics with the same constant
86 eV mass and arbitrary admissible noncondensed
chemical potential $\mu\le m$. Compare the distribution
in comoving momentum at the same comoving count,
not the free-streaming time of actual histories.

At fixed scale factor, increasing $\mu$ toward $m$
permits a lower thermal temperature while preserving
the count. The two occupation exponents
$(E-\mu)/T$ cross once: the saturated cooler gas has
more occupation at low $q$ and less at high $q$,
at equal count. Thus saturation maximizes the fraction
below every positive momentum cut among the ideal
noncondensed thermal states at that scale factor.
The limiting excited gas may have divergent occupation
at $q=0$; its radial count is finite.

## Ordering among saturated spectra

For $\mu=m$, set $r=m/T$, $\theta=aT$
and $x=\sqrt{q^2/\theta^2+r^2}-r$. At fixed
comoving count $\theta$ decreases with $r$:
$x$ decreases pointwise in $q/\theta$ as $r$
increases. Meanwhile $r\theta^2$ increases with $r$.
To see the latter, use $v=q/(\theta\sqrt{2r})$:
the exponent $x=\sqrt{r^2+2rv^2}-r$ increases
monotonically with $r$ at fixed $v$, toward $v^2$.
The radial number integral in $v$ therefore decreases,
so its scale $\theta\sqrt{2r}$ must increase at fixed
count.

Equivalently $q^2=\theta^2 x^2+2r\theta^2 x$. For a
larger $r$, the coefficient of $x^2$ is smaller
and that of $x$ larger; the two occupation slopes
therefore cross once. The larger-$r$ distribution
is enhanced on the low-$q$ side and depleted on the
high-$q$ side, with identical integrated count.
The nonrelativistic saturated limit gives the
coldest excited spectrum in this entire one-species
ideal thermal family:
$$
 f_{\rm lim}(q)=[\exp((q/q_T)^2)-1]^{-1},\quad
 n_c=\frac{q_T^3}{2\pi^2}\,\frac{\sqrt{\pi}}4\,\zeta(3/2).
$$
This is a bound on momentum fractions, not an
assertion that any physically available decoupling
reaches that limiting state.

## Numerical count targets

At $\Omega_\phi h^2=.060$, the limiting excited
spectrum has $\langle q/T_{\gamma,0}\rangle=.29876$
and minimum RMS .36904. Even this optimistic
spectrum puts at most **20.35%** of its count below
$q/T_{\gamma,0}=.1$, or **2.054%** below .01.
The explicit finite-mass quadrature approaches the
analytic limit from the warmer side.

A separate bookkeeping test allows a component at
exactly zero momentum and puts all other particles
in the most favorable limiting thermal gas, normalized
to their smaller count independently. Then at least
34.7% of the total count must be in the zero-momentum
component to have half the total below .1; to have
90% below .1 it must be at least 84.2%. At a cut of .01
the analogous fractions are 48.7% and 89.5%.
These are illustrative momentum targets, not
observational requirements or a physical model
of coherent condensate dynamics.

~~~sh
python -W error scalar_saturated_thermal_cold_order_reference.py \
 EOS.late.dat \
 results/warm-prompt-reheating-post-2mev-elastic-2026-11-26.json.gz
~~~

[Finite and limiting excited spectra and count brackets](results/scalar-saturated-thermal-cold-order-2026-11-26.json).
