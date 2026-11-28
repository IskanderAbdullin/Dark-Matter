# Stable-count retargeting with a frozen conversion loss

The [finite cooling count benchmark](finite_cooling_reheat_portal_reference.py)
chose `TR` to obtain stable Ωh² = .120 **without subtracting** separate
stable-to-scalar conversion. Each net pair event removes and
adds the same number of particles of the two species. In rest-density
units a newly made scalar contribution ΔΩφ corresponds to
a positive *loss* `(m_stable/m_scalar) ΔΩφ` from the stable count,
not ΔΩφ itself. Here the masses are 20 MeV and 86 eV. Kinetic heat
and pressure are not included in either rest-density proxy.

Use the independently recorded early forward conversion through
1 MeV and the disjoint frozen spontaneous forward yield starting
at 1 MeV. At the old roots, the combined losses in stable Ωh²
are 7.083×10⁻⁵ (instantaneous) and 1.188×10⁻⁴
(finite cooling with Tmax/TR=5). Recompute the prompt production
function on two quadratures and solve for the stable-count target
with these losses held fixed:

| History | Old TR [MeV] | Fixed forward-loss increase [keV] | Alternate early finite-node loss increase [keV] |
| --- | ---: | ---: | ---: |
| Instantaneous | 9.019070 | 1.1111 | 1.1060 |
| Finite, Tmax/TR=5 | 6.236182 | 1.3071 | 1.3003 |

The alternate early calculation retains its own first-order
projected depletion iteration, not an additional loss on top of
the primary result. The post-1 MeV part is small and identical in
the two calculations for a given history. The exact root for
the frozen loss and the local production slope agree at the
0.05% level on these shifts. None of these results includes
a self-consistent change to pair conversion when TR changes, a
joint inverse solution, or a recalculated spectrum/transfer.
In particular, this is not a physical error bar or a preferred
observed reheating temperature.

~~~sh
export OPENBLAS_NUM_THREADS=1
python -W error finite_cooling_frozen_count_retarget_reference.py EOS.dat \
  results/finite-cooling-{reheat-portal,during-production-pair-shape,post-anchor-forward-pair,forward-depletion-refined}-2026-11-26.json
~~~
