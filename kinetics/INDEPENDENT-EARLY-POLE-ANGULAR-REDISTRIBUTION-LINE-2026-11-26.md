# A conserving angular collision need not broaden the neutral transition

This is a fixed-radial-momentum response test, not a gauge rate or a new
cosmological trajectory. It isolates the consequence of keeping angular
gain and loss in a neutral, off-diagonal fermion bilinear driven by a weak
scalar probe. In the limit of identical velocities along identical paths,
the equal Abelian charges give opposite eikonal phases in the bilinear
and cancel. At nonzero momentum transfer, or for unequal masses, this is a
limit to check, not an assertion that the true thermal line is infinitely
sharp. A [single charged line can even decay nonexponentially](https://arxiv.org/abs/hep-ph/9607303).

Let the detuning on one radial shell be linear in direction, with half-span
$a$, so after centering it is $d-a\mu$. The toy collision operator
$-\nu(h-\langle h\rangle)$ redraws angles uniformly and exactly preserves
the isotropic component of the coherence. For $z=\epsilon+\nu-id$,
the loss-only result is $A=\langle(z+ia\mu)^{-1}\rangle$;
solving the same evolution with its gain term gives $A/(1-\nu A)$.
At $a=0$ this reduces to $1/(\epsilon-id)$, independent of the collision
rate. At large $\nu/a$ the conserving line has half-width $a^2/(3\nu)$,
from the isotropic component coupling to a damped dipole, while
loss-only gives a line with width of order $\nu$.
For uniform redraws the isolated decay pole at $d=-i\gamma$ exists on
the physical branch when $\nu>2a/\pi$, and
$\gamma=\nu-a\cot(a/\nu)$. The pole and measured half-width need not
coincide near this threshold.

| $\nu/a$ | Loss-only half-width / a | Redraw half-width / a | Smooth angular diffusion half-width / a |
| ---: | ---: | ---: | ---: |
| 1 | 1.4142 | .3344 | .4023 |
| 3 | 3.1623 | .1118 | .1142 |
| 10 | 10.0499 | .03336 | .03341 |
| 100 | 100.0050 | .00333336 | .00333341 |

Smooth diffusion uses Legendre damping eigenvalues
$\nu \ell(\ell+1)/2$, matching the dipole damping of the redraw model,
but not its higher multipoles. The independent full finite-angle solve
and Legendre recurrence reproduce the analytic redraw. The integrated
loss-only and redraw strengths both equal $\pi$; changes in peak and
width therefore do not imply changes in total oscillator strength.

This is not a replacement for the finite-temperature vertex. Real
changes in radial momenta, energy exchange, momentum-dependent dispersion,
nonlocal memory and unequal velocities are absent. Neither (
u) nor
the toy linewidth is inferred from the gauge coupling. A neutral
two-line coherence cannot be assigned a width by adding two independently
computed charged-particle lifetimes without checking correlated terms.

Reproduce from `kinetics/` with
`OPENBLAS_NUM_THREADS=1 python early_pole_angular_redistribution_line_reference.py`.
Raw data: `results/early-pole-angular-redistribution-line-2026-11-26.json`.
