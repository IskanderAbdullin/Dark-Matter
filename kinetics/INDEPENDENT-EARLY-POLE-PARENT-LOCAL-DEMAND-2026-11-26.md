# Where the early held-parent source demands its particles

This independent audit examines the **local** supply implicit in the
[momentum-resolved paired daughter calculation](INDEPENDENT-EARLY-POLE-PAIRED-HARTREE-MATRIX-2026-11-26.md).
The parent is still imposed at its thermal FD distribution. The scalar
has the same self-generated positive Hartree pole, and the fermion
daughter evolves with Pauli blocking and inverse decay.

For each scalar momentum, split the parent energy aperture wherever
either a daughter-bin edge or a parent-shell edge is crossed. Parent
shells are defined in fixed comoving momentum $p_P/T$, not energy,
because $M/T$ changes. Integrate forward and inverse counts separately
inside each exact subinterval. Then integrate over each internal BDF
step and divide their signed difference by the amount held in that
parent shell initially. The sum is separately checked against the
growth of the scalar population. This is a **demand-to-inventory**
ratio along a held trajectory, not a physical depleted distribution.

For the coarse daughter grid at $x=0.05$, the net demand is 5.68%
of the initial parent count overall. The most affected parent shell
in the 64-edge scan, $1.086 < p_P/T < 1.252$, supplies 70.50% of its
own initial count. Changing just the number of diagnostic parent
edges from 36 to 64 changes this peak from 69.76% to 70.50%.
Neither of these scans exceeds the initial inventory in any
represented cell above the declared denominator floor. Refining
both daughter grids while retaining 64 diagnostic parent edges
gives 5.61% global demand and 69.55% in the same parent shell.
The sum of net withdrawals agrees with the final scalar count
within 3 × 10⁻⁷ relative for these cases.

This concentration means that small total demand does not imply
small departure from the assumed maintained shape. It does not show
that such a shape is impossible: parent momentum redistribution,
or replenishment, can change the local count. The rates of those
processes were not derived here.

~~~sh
python -W error early_pole_parent_momentum_demand_reference.py
python -W error early_pole_parent_momentum_demand_reference.py --bins 64
python -W error early_pole_parent_momentum_demand_reference.py --scalar 130 --fermion 190 --bins 64
~~~

[Raw shell inventories and signed demands](results/early-pole-parent-local-demand-2026-11-26.json.gz).
