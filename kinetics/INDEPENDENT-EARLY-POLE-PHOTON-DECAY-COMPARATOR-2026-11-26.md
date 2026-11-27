# Shared pair chemistry can close the stimulated gate in a conditional evolution

Couple the finite photon heat reservoir to the conservative
scalar decay from the x=.05 ratio-.06 snapshot. The scalar
keeps its momentum distribution and Hartree pole. Both
fermions are replaced by joint instantaneous FD spectra
at the same starting number and energy moments; photons
keep a Bose shape at zero chemical potential in this first
test. The explicit visible reference radiation does not
exchange heat. All tracked charges share the same scale
factor and H. Pair phase space is fixed at the starting
fermion masses over the short interval; the unblocked central
normalization changes by below 6 × 10⁻⁵ if the parent mass is
updated from x=.05 to .052, before occupancy factors.
This is not a calculated thermal gauge spectral function.

| Shared gauge model | Final scalar n/Tref³ at x=.052 | Final μ_parent−μ_daughter over Tref |
|---|---:|---:|
| Pair source off, same fast-FD closure | .02279746 | 2.8110 |
| g=.001, hot photons | .02274323 | 2.3951 |
| g=.003, hot photons | .01902097 | .06136 |
| g=.003, colder initial photons | .02017513 | .36744 |
| g=.005, hot photons | .01730716 | 3.47 × 10⁻⁶ |

The strong illustrative case crosses the instantaneous
stimulated gate for q=4.1 × 10⁻⁶ around x=.050616. The
g=.004 hot case crosses around x=.05150; .0035 remains above
it at the endpoint. These thresholds apply only to the stated
hard kernel and instantaneous shape assumptions. The total
scalar count can still grow after the stimulated part changes
sign because spontaneous decays persist. The refined pair
quadrature at g=.003 changes the count by about 0.12%, while
the other archived momentum grid changes it by about 1.4%.
Strong branches required tightening ODE and independent ledger
quadrature tolerances; the earlier looser attempts did not
pass the prescribed count ledger.

**Separate photon-number budget.** Pure γγ↔fermion-pair
reactions conserve photon number plus both fermion counts.
A fixed-scale comparator keeps both this sum and total
energy, and lets photon fugacity vary rather than assuming
a process that fixes it to one. At exposure .8 in the hot
start, the photon temperature is 1.0285 with μγ/Tγ≈−.404
and daughter number .13001, versus daughter .15609 under
the zero-chemical-potential assumption. The independently
computed stationary μ/T≈−.399 and T≈.9751 are approached
in a longer fixed-scale run. No nontrivial zero-mode
population is needed in the examined states. This number
budget still assumes rapid thermal momentum redistribution,
not a study of true in-medium photon equilibration.

~~~sh
python -W error early_pole_joint_photon_decay_reference.py results/early-pole-spectrum-third-m06-2026-11-26.json.gz --g .004 --xend .052
python -W error early_pole_shared_photon_number_budget_reference.py results/early-pole-spectrum-third-m06-2026-11-26.json.gz --end .8
~~~

[Coupled counts, heat and gate](results/early-pole-joint-photon-decay-2026-11-26.json.gz);
[photon number/energy comparator](results/early-pole-shared-photon-number-budget-2026-11-26.json.gz).
