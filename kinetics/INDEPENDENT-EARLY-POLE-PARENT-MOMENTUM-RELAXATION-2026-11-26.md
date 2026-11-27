# Momentum mixing is distinct from chemical replenishment

This independent comparison extends the
[resolved parent calculation](INDEPENDENT-EARLY-POLE-RESOLVED-PARENT-2026-11-26.md)
with zero chemical replenishment and two schematic redistribution
operators for the parent momentum cells. Both use rates per e-fold
rather than a physical scattering matrix element.

1. At fixed reference temperature, relax toward cell-average
   Fermi–Dirac occupations with chemical potential chosen to match
   instantaneous total parent count. This may exchange energy with
   the omitted thermal bath.
2. For isolated parent mixing, choose both multipliers in the
   discrete Fermi–Dirac target to match count and projected parent
   energy. The mixing source then exchanges neither moment; the
   separately discretized decay still has its own projection error.

The parent energy in the second control is the exact shell-averaged
one-particle energy at the given fixed $M/T$. The Jacobians use
implicit derivatives of the target moments. The static tests give
number and energy residuals at roundoff and directional derivatives
within (9	imes10^{-11}) relative. Both relaxations preserve Pauli
bounds in the differential equation by targeting an allowed FD state.

At the common endpoint $x=0.05$:

| Redistribution | Rate | Coarse $Y_\phi$ | Refined $Y_\phi$ |
|---|---:|---:|---:|
| None | 0 | 0.000366014 | 0.000356239 |
| Fixed reference temperature | 3 | 0.000369825 | 0.000359601 |
| Fixed reference temperature | 100 | 0.000422993 | 0.000409326 |
| Isolated, energy conserving | 3 | 0.000365331 | 0.000353249 |
| Isolated, energy conserving | 100 | 0.000402047 | 0.000398298 |
| Prescribed parent, different control | — | 0.000459463 | 0.000443716 |

The fast isolated operator recovers about 48% of the missing yield
between the resolved/no-mixing and held-parent refined runs. At the
same rate, the fixed-temperature version recovers about 61%. The slow
isolated operator slightly decreases the count in both comparisons;
no monotonic yield claim follows. The fast isolated run conserves
parent plus scalar count to within $7\times10^{-15}T^3$ at the
refined endpoint. Largest parent-cell deficits are still substantially
mesh dependent, so those extrema are not converged. This is a
sensitivity to redistribution and its energy budget, not a
derivation of particle scattering rates.

The first refined isolated rate-100 launch stopped at an overly
strict combined static check. The bare-collision static check now runs
in isolation and the elastic target has independent two-moment and
Jacobian checks. The full trajectory check was not weakened.

~~~sh
python -W error early_pole_parent_kinetic_relaxation_reference.py --static-only
python -W error early_pole_parent_elastic_relaxation_reference.py --static-only
# Both drivers, rate 3 or 100: default coarse resolution or the refined flags
python -W error early_pole_parent_elastic_relaxation_reference.py --scalar 65 --fermion 90 --parent 72 --rate 100
~~~

Raw parent occupations: [coarse](results/early-pole-parent-momentum-relaxation-coarse-2026-11-26.json.gz),
[refined](results/early-pole-parent-momentum-relaxation-refined-2026-11-26.json.gz).
