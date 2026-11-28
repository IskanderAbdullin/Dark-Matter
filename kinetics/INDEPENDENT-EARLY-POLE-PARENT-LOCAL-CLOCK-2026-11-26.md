# Parent relaxation should be compared with the momentum-local drain

On each retained [spectral snapshot](INDEPENDENT-EARLY-POLE-SPECTRAL-SNAPSHOT-2026-11-26.md),
evaluate the bare decay/inverse-decay collision operator at the actual
instantaneous occupations. For each parent shell losing particles,
divide the signed loss by the number currently in that shell. Rates
are in units of the radiation-only reference H, not measured
microscopic scattering rates. The active filter excludes shells
below occupation 10⁻¹⁰, number 2 × 10⁻¹⁰ T³ or 10⁻⁴ of total
positive signed drain. The retained fraction is recorded.

At x=0.05 for mass ratio 0.06 without mixing, the highest retained
local drain is 22.9 reference H on the previous grid and 9.55 on
the third grid; the drain-weighted medians are 11.2 and 6.46.
This establishes that comparison to H alone is too weak for the
strongly affected band in these snapshots, not an accurate maximum
relaxation rate. Local derivatives are less converged than yields.
For illustration the second-grid continuation at x=0.10 reaches
a local rate near 100 reference H, explicitly not a converged rate.

Compare the collision direction to the energy- and number-matched
isolated target. In some drained cells the latter target is
already lower than the current occupation, so no positive constant
relaxation coefficient would replenish those specific cells at that
instant. Where the target instead restores population, the
ratio of drain to target deficit is retained only as an instantaneous
balance diagnostic. Its maximum can be dominated by a small
denominator. It is not a bound on a real elastic kernel.

~~~sh
python -W error early_pole_parent_local_collision_clock_reference.py results/early-pole-spectrum-third-m06-2026-11-26.json.gz
~~~

[Raw instantaneous rates](results/early-pole-parent-local-collision-clock-2026-11-26.json).
