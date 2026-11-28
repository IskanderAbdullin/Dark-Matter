# Resolve the selective late-loss bound by parent shell

Use the same isolated 1.206 MeV state and assumptions as the
[parent-weighted Bose bootstrap](INDEPENDENT-WARM-PAIR-PARENT-WEIGHTED-BOSE-CEILING-2026-11-26.md).
The stable parents redshift and can only disappear. Both charge
populations have the same isotropic per-mode shape. There are no
inverse gains, redistributing collisions, imported older tails,
or other sources. The already omitted stable count is inert and
added once. The Bose cap covers all newly produced scalar
momenta in this specific comparison.

Two stable subpanels have anchor upper momenta (p_i,p_j).
The angular average of the squared CM momentum is at most
((p_i^2+p_j^2)/(4b^2)). Pointwise it is at most
((p_i^2+p_j^2)/(2b^2)), which permits the amplitude
coefficient (C_{ij}=C(\sqrt{(p_i^2+p_j^2)/2})) rather than
using the largest incident support for every pair. All
photons, spin factors and invariant-flux conventions are
the same as in the [previous bound](INDEPENDENT-WARM-PAIR-SELECTIVE-LOSS-BOUND-2026-11-26.md).

The opposite-charge occupation in subpanel j is never larger
than its initial positive radial number divided by (2b^3).
Thus the integrated upper exposure for each stable component i is

~~~text
tau_i = (1+F)^2 / (1536 pi m² H_gamma)
        * sum_j C_ij (p_i²+p_j²) n_j,upper .
~~~

This bounds loss even if hotter subpanels are preferentially
depleted. Within each subpanel the reconstructed Fermi
occupation is monotone in kinetic energy, so both initial
number bounds use just its exact endpoint extrema and radial
phase volume. To avoid subtracting the accuracy of the initial
number quadrature from the final floor, take the known full
stable initial number, then subtract the upper per-subpanel
losses (n_{i,\rm upper}(1-\exp[-\tau_i])).

At the weakened photon floor and the bootstrapped
(F\le 0.000719), refining to 256 subpanels per original shell
gives conditional inclusive stable Ωh² **at least 56.16**
(initial 57.45; reference target .120). The same global
moment-plus-support bound gives 56.00 at this cap. The finer
radial count brackets enclose the independent initial physical
number; the lower survival envelope strengthens monotonically
over the four tested subdivisions. No convergence of the
actual kinetic solution is inferred.

This is a tighter loss envelope, not a new physical prediction.
Parent scattering, inverse creation, an inherited
scalar exterior or a different expansion history invalidates
the stated comparison if it violates the assumptions above.

~~~sh
D=results
OPENBLAS_NUM_THREADS=1 python -W error warm_pair_shell_survival_bound_reference.py \
  "$D"/warm-coupled-pair-b10-28-shape-2026-11-26.json "$D"/warm-coupled-pair-b10-count-audit-2026-11-26.json \
  "$D"/warm-pair-parent-weighted-bose-ceiling-2026-11-26.json
~~~
