# Early infrared production depends on the width of the parent-energy window

The [early-onset upstream control](https://github.com/IskanderAbdullin/Dark-Matter/blob/f77e75a43e99f03c5ae4dc1d732731bbb9ba3fca/kinetics/REHEATING-MATCHED-2026-11-25.md)
makes a very narrow infrared scalar population
at $M/T=.005$. The
[two-daughter geometry](https://github.com/IskanderAbdullin/Dark-Matter/blob/f77e75a43e99f03c5ae4dc1d732731bbb9ba3fca/kinetics/TWO-MASS-GEOMETRY-2026-11-25.md)
flags that its momentum could be smaller than
a thermal dispersion scale. Rather than treating
$m_\mathrm{pole}>p$ as a step-function veto,
compute the exact per-mode collision rate for
three **held** narrow Lorentz-invariant poles.

Write the pole masses as $M,m,b$ for parent,
scalar daughter and fermion daughter. For a scalar
with $E=\sqrt{p^2+m^2}$, define
$A=(M^2+m^2-b^2)/2$ and
$D=\sqrt{A^2-M^2m^2}=M p_*$.
If $M>m+b$ the parent-energy interval is
$$
 E_P^\pm=(AE\pm pD)/m^2.
$$
The lower root is rationalized for numerical stability.
The upper root is infinite for a massless scalar.
The collision coefficient per scalar occupation has
prefactor $g_P\Gamma M^2/(2 E p p_*)$ multiplying
the integral over that interval, with the actual
isotropic rest width $\Gamma$. The apparent
$1/p$ singularity at positive $m$ cancels
against the interval width $2pD/m^2$ at rest.
Integration over every scalar mode counts exactly
one daughter per time-dilated parent decay when
the two daughters are initially empty.

The quantum bracket is
$f_P(1-f_b)+(f_P-f_b)f_\phi$.
With parent and fermion daughter both FD at the
same temperature, parent chemical potential zero
and daughter chemical potential $\mu_b$, the
stimulated coefficient divided by the spontaneous
one is $1-\exp[(E+\mu_b)/T]$.
It changes sign at $\mu_b=-E$; a held thermal
fermion can suppress amplification without
removing spontaneous production. This identity
does not replace an evolving nonthermal daughter.

At $p/T=4.1\times10^{-6}$, with massless fermion
daughter, fixed rest width and initially empty
daughter distributions, the rate divided by the
chosen Hubble clock is:

| Held scalar pole mass over $T$ | Parent energy window $E_P/T$ | Seed/stimulated coefficient per Hubble time |
|---:|---:|---:|
| 0 | 1.524 to infinity | 127.16 |
| $4.1\times10^{-6}$ | 1.263 to 7.360 | 113.39 |
| $10^{-5}$ | .8385 to 1.8635 | 52.67 |
| $5.646\times10^{-5}$ | .2059 to .2381 | .6683 |
| $1.5\times10^{-4}$ | .0812 to .0857 | .03850 |

Thus a moderate held scalar pole can lower the
parent energy threshold while shrinking the
occupied window enough to reduce growth by
about 190-fold in the fourth row. That mass is
an illustrative scale, not a computed
self-energy of an initially empty scalar sector.
A fermion pole can also move the threshold or
close the channel. The retained data separate
fixed width from an illustrative scalar-Yukawa
Born width change. No real-time dispersion,
pole residue, damping, nonthermal parent
production or actual freezeout has been solved.

~~~sh
python -W error early_pole_decay_mode_clock_reference.py
~~~

[Mode scan and independent normalization controls](results/early-pole-decay-mode-clock-2026-11-26.json.gz).
