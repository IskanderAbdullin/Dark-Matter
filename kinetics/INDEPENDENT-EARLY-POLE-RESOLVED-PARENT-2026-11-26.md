# Resolving the finite parent momentum distribution changes early yield

This check complements the [count-only reservoir](INDEPENDENT-EARLY-POLE-PARENT-NUMBER-RESERVOIR-2026-11-26.md),
which preserved a prescribed parent shape. Unlike the upstream
[finite chemical bath](https://github.com/IskanderAbdullin/Dark-Matter/blob/f77e75a43e99f03c5ae4dc1d732731bbb9ba3fca/kinetics/FINITE-MOTHER-BATH-2026-11-25.md),
it does not enforce kinetic equilibrium. Here parent fermion occupation
evolves in each comoving momentum cell; the positive scalar Hartree
pole, massless fermion daughter, and fixed $T\propto1/a$ are unchanged.
There is initially a thermal parent and no daughter population.

For each scalar node, the overlap of the exact kinematic window with
one fermion-daughter cell and one parent cell gives a nonnegative
geometric multiplier $D_{ijk}$. The same signed event
$D_{ijk}[P_k(1+f_i-b_j)-f_i b_j]$ adds one scalar, adds one
fermion daughter and removes a parent. At both fermions' bounds
zero and one the flux points inward; the two daughter counts agree
and the total of parent plus either daughter is conserved when the
replenishment rate is zero. A fixed-gap directional derivative check
of the analytic collision Jacobian agrees within $3.3\times10^{-12}$
at the chosen check state. This construction is number conservative
but does not claim an exact discrete energy identity. The parent and
daughter occupancies are constant within each fermion cell.

| Cells (scalar, daughter, parent) | Held parent $Y$ | Live parent $Y$ | Relative decrease |
|---|---:|---:|---:|
| 37, 50, 74 | 0.00045946 | 0.00036601 | 20.34% |
| 67, 112, 122 | 0.00044372 | 0.00035624 | 19.72% |

At $x=0.05$, the refined live run consumes only 4.57% of the initial
parent count globally, but depletes the most affected shell
$1.308<p_P/T<1.361$ by 47.98%. The coarser shell with the largest
depletion loses 52.36% of its initial count. The refined live run
closes parent plus scalar to $3.2\times10^{-17}T^3$ at the endpoint;
both reported runs close paired daughter number within
$3\times10^{-17}T^3$. The roughly 20% reduction is relative to a
matched numerical held-parent run, not relative to a no-Hartree
comparison or to an observed relic abundance.

An additional pair of illustrative external chemical-relaxation
runs at fixed grid (47, 75, 86) gives yields 0.00036782 and 0.00042659
for rates 3 and 100 per e-fold respectively, with largest parent-cell
depletion 38.89% and 10.55%. This source changes total parent
number; the rates do not follow from a specified scattering model.
Momentum-only equilibration at fixed parent number is a separate
control. None of these calculations evolves the expansion, supplies
microscopic thermalization rates, or produces a freeze-out spectrum.

~~~sh
python -W error early_pole_resolved_parent_reference.py --static-only
python -W error early_pole_resolved_parent_reference.py --scalar 35 --fermion 40 --parent 40
python -W error early_pole_resolved_parent_reference.py --scalar 35 --fermion 40 --parent 40 --held
python -W error early_pole_resolved_parent_reference.py --scalar 65 --fermion 90 --parent 72
python -W error early_pole_resolved_parent_reference.py --scalar 65 --fermion 90 --parent 72 --held
python -W error early_pole_resolved_parent_reference.py --scalar 45 --fermion 60 --parent 48 --gamma 3
python -W error early_pole_resolved_parent_reference.py --scalar 45 --fermion 60 --parent 48 --gamma 100
~~~

Raw occupations: [coarse comparison](results/early-pole-resolved-parent-coarse-2026-11-26.json.gz),
[refined comparison](results/early-pole-resolved-parent-refined-2026-11-26.json.gz),
[chemical source test](results/early-pole-resolved-parent-chemical-2026-11-26.json.gz).
