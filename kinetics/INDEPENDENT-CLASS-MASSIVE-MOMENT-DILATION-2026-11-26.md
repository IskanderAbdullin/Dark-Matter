# Massive energy and pressure dilation on native quadrature nodes

The prior [number moment check](INDEPENDENT-CLASS-PSD-MANUAL-MOMENTA-2026-11-26.md)
does not measure how the same sampled derivative weights the
actual massive stress-energy across the relativistic transition.
For $D=qf'/f$, $\epsilon=\sqrt{q^2+(am/T_{\gamma,0})^2}$,
let $R=\int q^2\epsilon f$, $P=\int q^4 f/(3\epsilon)$,
and $K=\int q^6 f/(3\epsilon^3)$, all over the same domain.
Integrating by parts gives
$$
\begin{aligned}
 \int q^2\epsilon fD\,dq &= [q^3\epsilon f]-3(R+P),\\
 \int \frac{q^4 f D}{3\epsilon}\,dq
     &= \left[\frac{q^5 f}{3\epsilon}\right]-5P+K .
\end{aligned}
$$
The finite endpoint terms are included for finite-support
quadratures. The pressure identity tends to −4P for a
relativistic gas and −5P in the nonrelativistic limit with
vanishing endpoints. A synthetic exponential distribution
with a nonzero endpoint tests both identities. An independent
continuous quadrature over the archived spline panels also
reproduces the identities before inspecting the sparse nodes.

| Input on candidate CLASS grid | Max energy quadrature residual | Max pressure quadrature residual | Max energy dilation defect | Max pressure dilation defect |
| --- | ---: | ---: | ---: | ---: |
| Stable, raw 64 | 0.00690 | 0.00151 | 0.00760 | 0.00761 |
| Stable, raw 128 | 0.00209 | 0.000578 | 0.00541 | 0.00271 |
| Stable, raw 256 | 0.000742 | 0.0000743 | 0.01249 | 0.00201 |
| Stable, raw shifted 128 | 0.000798 | 0.000190 | 0.00134 | 0.00153 |
| Stable, pressure-balanced 256 | 0.0000702 | 0.0000137 | 0.0000400 | 0.0000308 |
| Stable, raw transformed 512 | 0.0000127 | 0.0000468 | 0.000398 | 0.000496 |
| Scalar, raw shifted 160 | 0.000145 | 0.000250 | 0.00270 | 0.00442 |
| Scalar, raw 384 | 0.00119 | 0.000597 | 0.000840 | 0.00176 |
| Scalar, pressure-balanced 128 | 0.000351 | 0.000270 | 0.000760 | 0.000571 |
| Scalar, pressure-balanced 256 | 0.000286 | 0.000118 | 0.000108 | 0.0000862 |

These are maxima over a set of scale factors including the
archived production surface, relativistic and nonrelativistic
epochs, and momentum-to-mass transition probes. The calculation
uses native-style sampling with an independent analytic
cubic-spline derivative, not the compiled solver's finite
difference. A grid can have small background residuals and
significantly larger errors in the perturbation's stress-energy
moments. The shifted and transformed schemes provide
independent numerical controls; none of the residuals by itself
is a bound on the transfer error or a physical prediction.

~~~sh
OPENBLAS_NUM_THREADS=1 python -W error archived_class_massive_moment_dilation_reference.py \
   fine stable-balanced.dat scalar-balanced.dat
~~~

Generate the fine phase-space table and the two positive smooth
comparators with the existing pressure-balanced and scalar-tracer
setups. The stable transformed control uses the native indefinite
trapezoid, q = 1/t − 1; the other controls use finite linear q.
