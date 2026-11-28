# A reversible threshold control changes scalar number but preserves energy

The [vacuum 4-to-2 threshold check](INDEPENDENT-EARLY-POLE-SCALAR-NUMBER-CHANGE-2026-11-26.md)
motivates a reversible numerical test at fixed pole mass. Enumerate
unordered quadruplets of cold shells and use their permutation
multiplicities with the threshold phase space. For each event
deposit two equal-energy outgoing representatives at neighboring
hard nodes with positive energy weights. This is an equal-energy
quadrature closure, not the exact outgoing angular spectrum.
The four cold losses and the two hard gains preserve the first
energy moment. Reverse exactly the same graph with Bose factors,
so there is no hand-set irreversible energy leak.

A Bose distribution with zero chemical potential balances the
stoichiometric graph. Positive chemical potential gives net number
loss; negative chemical potential gives the opposite direction.
Empty receiving states have nonnegative flux. The entropy production
is nonnegative event by event through the stable gain-loss affinity.

On the three isolated frozen examples, at a normalized interval
equal to 0.03 times the inverse initial fractional loss clock,
about half of the represented scalar count remains, largely at
the new hard nodes. The endpoint does not follow from the initial
linear extrapolation: the incoming and outgoing occupation
factors evolve rapidly in this artificial threshold graph.
Energy drift is below 3 × 10⁻¹⁵; none of these runs assigns
this normalized interval a physical elapsed time. Expansion,
medium screening, scalar decay sourcing, angular broadening and
other number-changing channels are all omitted here.

~~~sh
python -W error early_pole_scalar_number_change_stencil_reference.py results/early-pole-spectrum-third-m06-2026-11-26.json.gz
~~~

[Three raw trajectories and static checks](results/early-pole-scalar-number-change-stencil-2026-11-26.json.gz).
