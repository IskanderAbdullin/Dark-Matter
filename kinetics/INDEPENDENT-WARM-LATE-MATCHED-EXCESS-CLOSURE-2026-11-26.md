# Matching each daughter's excess is still insufficient to fix its total current

The [matched-total-excess control](INDEPENDENT-WARM-LATE-EQUAL-EXCESS-AND-PARENT-2026-11-26.md)
changes how the excess is divided between scalar and
light fermion. Here independently choose each species'
packet center so its own convex excess matches. Then
represent both shapes on the union of their two momentum
grids by restricting the original affine kinetic entropy
variable into every new cell, not reconstructing from
moments again. Both preserve the per-cohort birth number
and energy and have the same support and visible background.

At the unscaled earlier positive case's two target excesses,
all five tested combinations of widths have positive total
current, from +7.05e−10 to +1.10e−9 of initial dark energy
per Hubble time. This is a useful negative search result:
matching species individually removed the sign contrast
in this sample, not necessarily in other shapes.

Scale both specified target excesses by 1.29. The new
targets are 1.7389843e−16 for the scalar and 8.4101574e−17
for the light species, each divided by initial total dark
energy. The two widths give **opposite** signed currents:

| Half-width of both species (GeV) | Scalar / light centers (GeV) | Total current / (rho H) |
| --- | --- | ---: |
| .00250 | .7136187867 / .7136568960 | −2.64664e−10 |
| .00025 | .7149822250 / .7148030823 | +9.59003e−11 |

The two targets agree in each shape to better than
4.83e−13 relative and the measured initial species number
and energy stay within 2.23e−16 of the original common
phase. Both independent angular orders give the same
signs. These examples rule out predicting the current's
sign from just those moments to the tested tolerance.
This is a frozen prompt closure with closed 1.05 GeV
support, not a general sign boundary or an inference
about the true outside phase.

## Dynamic follow-up of the separate equal-total-W pair

For the original equal-*total*-W pair, advance both histories
using three fixed positive roots per each of the same six
intervals. The broader packet reaches a 1.54% sampled rise
above its initial W before the signed current turns between
Hdt=1e−8 and 3e−8. The narrower packet has negative current
at all recorded times. At Hdt=1e−7 the terminal W values
are 1.67911e−16 and 1.45746e−16 of initial dark energy, a
15.21% difference relative to the narrow value. Both are
below their common starting W. The history checks require
the exact input phase hashes and verify the resolved W and
invariants independently. No interpolated reversal time,
time-grid convergence, or finite-parent interpretation is
inferred for this separate pair.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MALLOC_ARENA_MAX=2
a=results/warm-corrected-anchor-signed-floor-2026-11-26.json
b=results/warm-multishell-corrected-q64-fine-stage-48-2026-11-26.json.gz
c=results/warm-late-equal-excess-shape-rate-2026-11-26.json
python -W error warm_late_equal_species_excess_shape_reference.py EOS.dat "$b" "$c"
python -W error warm_late_equal_species_excess_shape_reference.py EOS.dat "$b" "$c" --target-factor 1.29 --extremes-only
python -W error warm_late_near_crossing_shape_evolution_reference.py EOS.dat "$a" "$b" broad --hard .7134635061388191 --halfwidth .0025
python -W error warm_late_near_crossing_shape_evolution_reference.py EOS.dat "$a" "$b" narrow --hard .7148726662332907 --halfwidth .00025
python -W error warm_late_equal_excess_history_audit_reference.py "$c" broad narrow
~~~

[Species-matched initial comparison](results/warm-late-equal-species-excess-129-2026-11-26.json);
[unscaled controls](results/warm-late-equal-species-excess-base-2026-11-26.json.gz);
[matched-total-W history audit](results/warm-late-equal-excess-live-audit-2026-11-26.json);
[broad](results/warm-late-equal-excess-live-broad-2026-11-26.json.gz) and
[narrow](results/warm-late-equal-excess-live-narrow-2026-11-26.json.gz) saved histories.
