# Full finite lowered halos: where the mixing obstruction moves

The [local mixing check](INDEPENDENT-HALO-BOSE-HIGH-OCCUPATION-MIXING-2026-11-26.md)
covered only the central contribution and high occupation.
Here use the [upstream finite lowered Bose ansatz](https://github.com/IskanderAbdullin/Dark-Matter/blob/f77e75a43e99f03c5ae4dc1d732731bbb9ba3fca/papers/02-finite-halo-truncation.md)
with a smooth subcritical center and a genuine outer
escape boundary. The parameters are imposed, not fit.

For escape depth $W$, potential offset $y(\xi)$
and kinetic energy $t$, set the dimensionless scalar
occupation to
$$
 f=B(\delta+y+t)-B(\delta+W), \quad y+t<W,
 \qquad B(x)=(e^x-1)^{-1},
$$
and zero outside. Choose $\delta$ so that its actual
peak at $t=y=0$ equals $F$. Integrate the
spherical Poisson equation in logarithmic radius until
$y=W$, with density calculated from this lowered
distribution at every step and the actual central
density used for the radial scale $R$.
Examples at $F=7000$ and $W=2,5,8$ terminate
at $r_t/R=18.87,61.41,326.71$. Separate full-volume
integration closes the normalization to $1.2\times10^{-5}$
relative in the most extended example.
This does not test dynamical stability.

## Test both full distributions

Write $d_i(f_*)=\mathcal D_i(f_*)/N_i$ for a
finite initially homogeneous patch, and $d_h(f_*)$
for the whole final halo. The primordial momentum
example is the earlier smooth compact-tier logistic
shape with a nonzero tail (its peak is matched
exactly to the halo's). If the halo retains fraction
$\eta=N_h/N_i$, a necessary condition is
$\eta d_h(f_*)\le d_i(f_*)$ at each threshold.
This leaves possible debris out of the final
excess, in the conservative direction.

Across sampled $f_*/F=10^{-10}$ through .9 all
three $F=7000$ halos and an $F=.512,\ W=5$
control have no detected obstruction to unit
retention from those samples. This is not a
positive formation result. In fact exact full
retention is obstructed for this infinite-tail
progenitor at sufficiently **low** occupation,
unlike the clipped plateau's high-end obstruction.

To see it, write
$1-d(f_*)=N^{-1}\int\min(f,f_*)\,d^3x\,d^3p$.
For the finite halo this deficit is at most
$f_* V_h/N_h$, where $V_h$ is its finite
occupied phase volume. For a uniform primordial
patch the deficit is at least $f_* V_i(f>f_*)/N_i$.
The logistic tail makes the latter volume per
particle arbitrarily large as the threshold tends
to zero. For a given halo, take the threshold low
enough that the initial deficit exceeds the final
one; then $d_h>d_i$, contradicting exact full
retention. The sufficient $\log(f_*/F)$ thresholds
in the examples range from about $-106$ to
$-6485$, far below the directly sampled values.
This limiting statement need not imply any
appreciable missing mass or a meaningful
observational effect. Removing weakly occupied
debris or allowing a noncompact envelope changes
the full-retention problem.

Neither the local nor global test fits a galaxy,
and no gravitational evolution or collisions have
been simulated.

~~~sh
python -W error halo_bose_finite_lowered_mixing_reference.py
~~~

[Finite masses, phase volumes and normalized mixing checks](results/halo-bose-finite-lowered-mixing-2026-11-26.json).
