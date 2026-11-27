# Angular deposition and hot-shell refinement in the coupled threshold control

Separate the angular geometry from the still-unknown in-medium
4→2 amplitude. Use quasi-random isotropic directions for four
incoming cold particles and the exact two-body energies at each
sampled outgoing cosine. Retain the threshold matrix element,
but include the exact on-shell phase space and incoming 1/E
factors. Project each product separately between neighbors with
positive energy weights. Reverse the same graph, and re-solve
the Hartree mass as the live decay/chemical model evolves.

On the frozen original momentum grids, the angular quadrature
changes the formal total loss rate by less than 4 × 10⁻⁵
relative near x=0.05; most outgoing energies share the same
two receiving nodes. That agreement alone does not resolve the
product bandwidth.

At the nuisance multiplier 10⁻⁷, coupled runs to x=0.052 on
the two original grids both lose about 51% of the initial
scalar count to the chemical channel, before accounting for
new particles produced by decay. The equal-energy versions
yield 0.000235717 and 0.000239782; angular versions with four
samples yield 0.000235652 and 0.000239705. Eight samples on the
third grid gives 0.000235655. This is a reproducibility check
of the threshold comparator, not a measurement of the medium
rate. Local hot-shell refinement and Fermi gain remain
important distinct checks.

~~~sh
python -W error early_pole_scalar_number_change_angular_reference.py results/early-pole-spectrum-third-m06-2026-11-26.json.gz
python -W error early_pole_scalar_angular_chemical_coupled_reference.py results/early-pole-spectrum-third-m06-2026-11-26.json.gz --power 3 --strength 1e-7
~~~

[Frozen angular ledgers](results/early-pole-number-change-angular-2026-11-26.json.gz);
[coupled comparisons](results/early-pole-live-chemical-angular-comparison-2026-11-26.json).
