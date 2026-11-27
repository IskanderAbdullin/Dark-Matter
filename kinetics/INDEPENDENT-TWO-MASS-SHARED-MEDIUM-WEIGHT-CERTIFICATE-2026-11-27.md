# An admissible weight box for a fixed shared-medium step

This extends the [quadratic two-event test](INDEPENDENT-TWO-MASS-SHARED-MEDIUM-EVENT-WEIGHTS-2026-11-27.md).
It checks whether chosen bracketing shells can simultaneously balance
energy for an **already proposed** vector of signed event extents.
Such a check can be used before rebuilding shell support or retrying a
smaller step; it does not solve the extent/time-step problem.

Write v_j(w_j)=a_j+w_j d_j for the per-event deposit, H for the
constant energy Hessian, and ξ_j for the fixed signed extents. The
midpoint gradient is
`g(w)=grad E(n)+H sum_k ξ_k v_k(w_k)/2`.
At each probe vector z, the exact range of z·g(w) on the weight
box is its old-state value plus the sum of the individual range
endpoints of `ξ_k zᵀ H v_k(w_k)/2`. Each event coordinate is
independent, so its endpoint extrema can be summed with their
actual signs. This also bounds each coordinate of the updated state,
which must stay positive and below each Fermi capacity.

If g(w)·v_j at the left endpoint is strictly negative throughout
the box and its value at the right endpoint is strictly positive,
the frozen-gradient root map
`T_j(w)=−(a_j·g(w))/(d_j·g(w))` maps the box to its interior.
For positive lower slope s_j=min_w d_j·g(w), an infinity-norm
Lipschitz bound is the largest row sum of

    |ξ_k| max_{x at j's endpoints} |v_j(x)ᵀ H d_k| / (2 s_j).

If this bound is below one, the contraction theorem gives a unique
common energy root **in that box for the fixed extents**. Iterating
the weight map converges there. Conversely, positive lower slope
and a left endpoint strictly positive everywhere (or a right endpoint
strictly negative everywhere) for any one event exclude a common root.
Other inconclusive boxes are not automatically infeasible.

At the smallest tested step, cancellation from the opposing extents
makes the signed margin positive (7.80e−7) in the ±4e−5 box, while the
symmetric absolute bound is inconclusive (−3.21e−7). At the two longer
steps neither narrow box can contain a common root, while the original
[0,1] box works. Across the valid boxes the maximum contraction bound
is 5.20e−5; 4–5 iterations agree with the prior coupled weights to
5.28e−17 in absolute value. Direct corner enumeration independently
matches the computed energy and occupancy extrema within 4.90e−20.
The analytic bounds are evaluated in extended precision with a
1e−14 margin guard for this fixture, not general interval arithmetic.

Replay from the kinetics directory:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error two_mass_shared_medium_weight_certificate_reference.py
```
[Machine-readable receipt](results/two-mass-shared-medium-weight-certificate-2026-11-27.json).
No physical thermal dispersion, reaction rate, or kinetic convergence
is inferred from this manufactured medium.
