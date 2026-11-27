# Large Bose encounter counts do not imply fast kinetic equilibration

The [held orbit calculation](INDEPENDENT-HALO-BOSE-QUANTUM-ORBITAL-ELASTIC-2026-11-26.md)
gives strong stimulation of gross events in the center.
To see how much survives in a kinetic perturbation,
compute two local Rayleigh quotients of the linearized
elastic Uehling–Uhlenbeck operator for the same fixed
halo scales. This is an application of quantum collision
brackets, not a determination of the full
[transport spectrum](https://arxiv.org/abs/1207.6855).

Let $\delta f=f(1+f)\psi$ at fixed radius, and
$g=|\mathbf v_1-\mathbf v_2|$. The local rate for
the test shape is the collision bracket
$$
 \Gamma_\psi=
 \frac{\frac14\int d^3v_1\,d^3v_2\,d\Omega'\,
       (d\sigma/d\Omega')\,g\, f_1 f_2(1+f_3)(1+f_4)
       (\psi_1+\psi_2-\psi_3-\psi_4)^2}
      {\int d^3v\, f(1+f)\psi^2},
$$
with the same one-particle phase-space normalization
and isotropic vacuum total cross section as the orbit
comparison. Elastic particle number, energy and
momentum perturbations make the square vanish.
For other shapes the quotient gives initial linear
dissipation; it does not imply exponential evolution
of that shape. In a fixed local collision operator
it is a variational upper comparison on the
slowest rate in a subspace containing the trial,
not a guarantee that equilibration is at least this fast.

* **Shear:** $\psi=v_xv_y$. For fixed incoming/outgoing
  relative direction deflection $\gamma$, averaging the
  tensor over common rotations gives
  $(\Delta\psi)^2=g^4(1-\cos^2\gamma)/20$.
  Average the final Bose weights over the outgoing
  azimuth and integrate the remaining polar dependence.
* **Isotropic shape:** $\psi=v^4-c_2v^2-c_0$, projecting
  out particle number and kinetic energy in the
  $f(1+f)$ metric. With pair CM speed $V$,
  $\Delta\psi=2(Vg)^2(c_\mathrm{in}^2-c_\mathrm{out}^2)$
  up to the sign convention for incoming minus outgoing.
  This measures loss of a fourth-moment distortion, not a
  proof of convergence of arbitrary speed distributions.

For comparison use a Maxwell bath of the **same local
density and velocity variance**, with bare vacuum elastic
scattering. Its analytical coefficients are
$8\sigma_v/[5\Gamma(3/2)]$ and
$16\sigma_v/[15\Gamma(3/2)]$, respectively, in units
where the common prefactor is number density times
elastic cross section. Cartesian invariant and angular
checks and independent scrambled integrations check the
quoted values (worst combined relative resolution and
scramble spread about 1.1% for the shape test).

| At peak occupation 7,000 | Gross circular encounter rate / $\Omega_c$ | Shear / $\Omega_c$ | Isotropic shape / $\Omega_c$ |
|---:|---:|---:|---:|
| $.1R$ | 13.34 | .07284 | .02982 |
| $.3R$ | 1.684 | .06539 | .02812 |
| $R$ | .1843 | .04571 | .02277 |
| $10R$ | .00625 | .00517 | .00341 |

Here $\Omega_c=V_c/r$, not inverse orbital period.
At $.1R$, the gross encounter enhancement over its
vacuum comparator is about 478, while the local shear
and shape rates relative to their matched Maxwell
bath are only 2.08 and 1.27. The two test rates
are $3.12$ and $1.28\ {\rm Gyr}^{-1}$ there.
Thus gross stimulation cannot justify a fluid
closure on the orbital dynamical clock in this
held gas. Several Gyr can still give nonnegligible
local kinetic evolution; the calculation does
not settle a global halo steady state or formation,
since it includes neither spatial streaming in the
operator nor feedback of collisions on the halo.

~~~sh
python -W error halo_bose_linearized_shear_reference.py
python -W error halo_bose_linearized_shape_reference.py
python -W error halo_bose_kinetic_clock_comparison_reference.py \
 results/halo-bose-quantum-orbital-elastic-2026-11-26.json \
 results/halo-bose-linearized-shear-2026-11-26.json \
 results/halo-bose-linearized-shape-2026-11-26.json
~~~

[Shear](results/halo-bose-linearized-shear-2026-11-26.json),
[isotropic shape](results/halo-bose-linearized-shape-2026-11-26.json)
and [dynamical-clock comparison](results/halo-bose-kinetic-clock-comparison-2026-11-26.json).
