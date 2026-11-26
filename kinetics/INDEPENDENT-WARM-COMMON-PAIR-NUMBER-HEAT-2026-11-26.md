# Particle and heat Onsager response on the common pair chart

Use the [shared scalar-pair coordinates](INDEPENDENT-WARM-PAIR-COMMON-CHART-2026-11-26.md)
on explicitly thermal FD light and BE scalar baths, without dilute
or classical statistical factors. The equilibrium anchor is T = 17 MeV
and common chemical potential μ = −30 MeV. Perturb number affinity
x₁ and the dimensionless inverse-temperature difference x₂ with
the chemical offset chosen so the latter is conjugate to heat.
For an interior event of total energy E, where both thermal baths
are represented, exact energy conservation gives

    log(G/L) = 2 x₁ + x₂ (E − 2 μ)/T.

The transferred particle number and heat/T are the corresponding
event vector v = (2, (E − 2 μ)/T) per forward fermion gain event.
Thus the linear matrix is an integral of w v vᵀ, where w is the
nonnegative equilibrium forward event measure. This independently
predicts reciprocity, nonnegative dissipation, and the entire first
row from equilibrium event number and energy. At the artificial
finite support boundary the imposed vacuum is not exactly thermal;
the computed equilibrium drive from the boundary is below 8e-15
per equilibrium forward event for this test.

At the finer quadrature, central differences over two affinity step
sizes and an even-error extrapolation give, per equilibrium
forward event,

    [[4.00000000, 20.28866549],
     [20.28866549, 107.08186291]].

The eigenvalues of the symmetric part are approximately 0.15052 and
110.93134, both positive. The two independently measured off-diagonal
terms agree within 5.5e-12 on the finer rule and 7.1e-12 on the
coarser. The first row matches its independently accumulated analytic
moments. The heat diagonal changes by 2.3e-6 relatively between the
two rules. The collision projection and thermal finite-difference
tolerances are distinct from that angular/radial convergence.

This is a controlled two-bath number/heat test of scalar-pair
conversion. It is not a result for the full vacuum vector portal,
an evolving cosmology, or the other selected reactions.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error warm_common_pair_onsager_reference.py
python -W error warm_common_pair_onsager_reference.py --fine
~~~

[Full response and diagnostics](results/warm-common-pair-number-heat-onsager-2026-11-26.json).
