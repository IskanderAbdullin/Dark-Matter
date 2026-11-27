# Extend the physical scalar power check to k=2 h/Mpc

The [earlier exact-mode test](INDEPENDENT-WARM-DIRECT-NATIVE-EXACT-K-2026-11-26.md)
stopped at k=1. Keep the early common initialization and
the same physical scalar fraction and spectrum. Extend
RK to k support 2.12 h/Mpc and request the same points
exactly in the power tables for every model. The trace
includes the same scalar species at negligible density,
with unchanged PSD and numerical settings. Also collect
the completed independent NDF fine physical/trace pair:
128 momenta, lmax=32 and tolerance 8e−10. Its sampled grid
is different, so compare both cubic spline and PCHIP
in log k; the difference of their signed effects is
at most 1.13e−10 over the common comparison.

| At z=5 | Physical minus direct trace, Delta ln P at k=1 | Same at k=2 |
| --- | ---: | ---: |
| Original standard NDF | −5.05224e−6 | −7.85911e−6 |
| Wider early RK, q=96 and lmax=26 | −5.08819e−6 | −7.88597e−6 |
| Fine NDF physical/trace pair | −5.11556e−6 | −7.91875e−6 |

The three present estimates at k=2 span about .76% of the
effect, with the same sign at z=5 and 20. At z=20 the
resolved wide RK result there is −7.80726e−6. The resolved
held-versus-direct shape contrast at the same physical
fraction, with separate trace PSDs, is −9.10527e−7 at
k=2 and z=5. The 1e−4 amplified tangent predicts the
wide physical values within 7.31e−9 on the sampled grid.

## Vary the numerical axes individually

The compact wide RK run has q=70 and lmax=20. Raise q
to 96 at fixed lmax, raise lmax to 26 at fixed q, and
raise both. For k=2, z=5 the isolated q change is
−2.697e−8, while the isolated hierarchy-extension change
is only +1.329e−10; the summed shifts agree with both
changes together within 2e−13 there. Over the checked
grid the maximum nonadditive finite difference is
below 4e−13. Increasing angular hierarchy alone is
therefore not the efficient refinement for these
settings. This is a finite cross, not a proof of
continuum convergence or a universal hierarchy rule.

These are linear native tests of a prescribed, already
produced 86 eV distribution replacing an equal portion
of a CDM proxy. They do not include late interactions,
dynamical production perturbations, uncertain total
dark density or observational comparison. The full
higher-order amplified ladder is not needed for, or
claimed complete by, the separately completed fine
physical/trace pair.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
mkdir -p wide cross
python -W error warm_direct_vacuum_native_wide_rk_reference.py early wide > wide/spec.json
python -W error warm_direct_vacuum_native_hierarchy_cross_reference.py wide cross > cross/spec.json
# Run the generated CLASS inputs, storing a log per input name.
D=results
python -W error warm_direct_vacuum_native_wide_audit_reference.py wide "$D/warm-direct-vacuum-native-class-pilot-2026-11-26.json" > wide.json
python -W error warm_direct_vacuum_native_cross_audit_reference.py wide.json cross
python -W error warm_direct_vacuum_native_fine_pair_reference.py native wide.json "$D/warm-direct-vacuum-native-class-pilot-2026-11-26.json"
~~~

[Wide exact modes](results/warm-direct-native-wide-2026-11-26.json);
[individual numerical axes](results/warm-direct-native-cross-2026-11-26.json);
[fine NDF pair](results/warm-direct-native-fine-pair-2026-11-26.json).
Raw power tables and temporary-path-elided logs:
[compact](results/warm-direct-native-compact-raw-2026-11-26.json.gz),
[resolved](results/warm-direct-native-resolved-raw-2026-11-26.json.gz),
[held profiles](results/warm-direct-native-wide-held-raw-2026-11-26.json.gz),
[cross](results/warm-direct-native-cross-raw-2026-11-26.json.gz), and
[fine pair](results/warm-direct-native-fine-pair-raw-2026-11-26.json.gz).

## Resolve q more densely after isolating the numerical axes

Keep the compact lmax=20 and all exact-k RK/
common-start choices, but increase the radial
momentum nodes alone to 168. The z=5 physical power
difference is −5.11522e−6 at k=1 and −7.92412e−6 at
k=2. At z=20 it is −5.09125e−6 and −7.84514e−6.
Compared with the independent fine NDF result the
difference at k=2 is below 5.38e−9 absolute at both
redshifts, about .07% of the signal. This is agreement
of two numerically different tests, not a strict
continuum error bound.

Recreate the manual radial trapezoid on the exact
tabulated PSD. The present native normalization fixes
abundance even when the raw sampled count varies.
Over the 70, 96 and 128 node grids the normalized
q² errors are not monotone; on 168 the magnitude is
.071%, and on a sampling-only 336 grid .0056%.
An inverse-q moment is more sensitive to low momentum
and need not improve monotonically even when q² does.
Simple moment accuracy is a useful input diagnostic,
not a substitute for comparing the actual power.

~~~sh
python -W error warm_direct_vacuum_high_q_native_reference.py wide q168 --points 168
for arm in trace physical; do
 /path/to/class "q168/$arm.ini" >"q168/$arm.log"
done
python -W error warm_direct_vacuum_q_sampling_audit_reference.py native/direct-sig07.dat
~~~

[Exact q=168 power differences](results/warm-direct-vacuum-native-168-q-2026-11-26.json),
[raw q=168 tables](results/warm-direct-vacuum-native-168-q-raw-2026-11-26.json.gz),
and [manual radial quadrature audit](results/warm-direct-vacuum-native-q-sampling-2026-11-26.json).

At the same radial spacing, shorten only the high-q
range to about 18. This keeps over .9999998 of the
second-moment integral of the positive PSD. The k≥1
power comparisons change by less than 3.38e−9, much
less than the shift from refining radial spacing with
the original range. However, the initial all-mode
comparison target 1e−8 is **not met**: maxima are
1.10e−8 and 1.79e−8 on the two spacings, driven by
k=.1 and .3. Do not claim that tail chopping preserves
the small-mode response at that tolerance. These changes
may include numerical effects of the altered solver
setup; they are not evidence for a large physical high-q
population.

~~~sh
mkdir -p cap
python -W error warm_direct_vacuum_q_cap_spacing_control_reference.py wide cap > cap/spec.json
# Run each generated CLASS input.
python -W error warm_direct_vacuum_q_cap_native_audit_reference.py wide.json q168 cap
~~~

[Fixed-spacing outer cutoff control](results/warm-direct-vacuum-native-q-cap-2026-11-26.json)
and [raw power and elided-path logs](results/warm-direct-vacuum-native-q-cap-raw-2026-11-26.json.gz).
