# Relax the hard packet stress with a direct-occupation prompt solve

The [two-packet tests](INDEPENDENT-WARM-HIGH-INVERSE-KNOWN-COHORT-SHAPE-BOUND-2026-11-26.md)
can activate a relatively enhanced high inverse source even
though each stored cohort's number and mean energy are
preserved. In a separate fixed-clock prompt-only solve, the
usual centered residual encountered an exponential overflow
on a nonlinear trial. This is not evidence of a physical
instability or exhausted computer resources.

Use the original prompt equations directly in finite FD/BE
occupation variables (A=f_L f_phi) and
(B=(1-f_L)(1+f_phi)), preserving the same phase-space
quadrature and parent stationarity. The alternate limited
trial path does not form the potentially enormous centered
(operatorname{expm1}) tilt. Its value and Jacobian-vector
product on the ordinary archive agree with the previous
formula to (5.2	imes10^{-13}) and (5.9	imes10^{-15})
relative, respectively; the independent finite-difference
Jacobian comparison is about (1.8	imes10^{-9}).

Both .52 and .58 GeV synthetic packet realizations complete
the fixed-clock evolution over (HDelta t=10^{-5}),
using two time subdivision levels. Conserved species
numbers and total energy agree with their starting values
to about (7	imes10^{-15}) at the sampled endpoints.
The initial normalized bare high inverse energy rates
are (1.49	imes10^{-20}) and (6.64	imes10^{-21}) in
units of retained dark energy times Hubble. By the last
sample both have fallen by at least five orders of magnitude
on **both** time grids. The fine late rate values are
however much smaller than the coarse ones; there is no
converged terminal high-tail rate to quote from this test.

This is a short isolated prompt relaxation. No coupled
radiation expansion, additional scattering/source, physical
unmeasured tail, or generic solver guarantee is inferred.
The original production integrator is unchanged.

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error warm_direct_prompt_stress_reference.py EOS.dat
~~~

[Diagnostics](results/warm-direct-prompt-stress-2026-11-26.json.gz).
