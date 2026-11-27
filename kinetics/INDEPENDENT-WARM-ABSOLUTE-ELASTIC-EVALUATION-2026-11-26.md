# Evaluate the reversible elastic panel without a large relative exponent

The [matched comoving expansion](INDEPENDENT-WARM-FINITE-PARENT-COMOVING-EXPANSION-2026-11-26.md)
normally uses elastic panels in coordinates relative to a thermal
reference. In a strongly displaced trial state the positive factor
exp(−z_scalar−z_light) can overflow while a separate thermal
factor would have made their product small. This is a failure of
the numerical factorization, not evidence of divergent physical
occupations or of a cloud capacity limit.

For the same below-pole panel and angular projection, form each
ordered rate from the bounded Fermi and physical Bose occupations
and the final-state blocking factors directly. Use the difference
of their log detailed-balance ratios with expm1 to evaluate the
antisymmetric net flux stably near balance. The analytic derivative
uses those same finite ordered rates. This optional evaluation
does not change the existing backward-Euler default.

At the exact .80 archive, on 9,901 elastic panels:

- Relative difference in flux versus the centered evaluation: 3.22e-14.
- Relative difference in the analytic directional derivative: 1.14e-15.
- Relative directional difference against an independent finite difference: 3.94e-10.
- Number residual relative to projected gross activity: 1.04e-15; the
  energy residual is zero at recorded precision.

A separate coordinate test uses the same physical occupations with
a deliberately cold synthetic comparison reference. The centered
presentation overflows, while the absolute rate and derivative stay
the same because they depend only on the absolute occupations.
This does not make implicit integrators universally admissible:
unphysical trial occupations and failed Newton solutions are still
rejected. No physical process or far-tail distribution is added.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MALLOC_ARENA_MAX=2
python -W error warm_absolute_elastic_audit_reference.py competing080/tapes/step-5.json
~~~

[Audit results](results/warm-absolute-elastic-evaluation-audit-2026-11-26.json).
