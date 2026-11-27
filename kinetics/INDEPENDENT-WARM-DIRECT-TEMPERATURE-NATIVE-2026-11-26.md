# Reheating power response at matched and absolute scalar counts

The [positive temperature-dependent spectra](INDEPENDENT-WARM-DIRECT-TEMPERATURE-SHAPES-2026-11-26.md)
test 50 and 100 MeV with separate own-shape negligible
traces. Fix omega_dark=.12 with an independently
supplied cold component. Change just shape at the
50 MeV scalar yield, then switch to the 100 MeV
absolute direct scalar yield at that new shape.

| At z=5, compact exact mode | k=2 h/Mpc | k=10 h/Mpc |
|---|---:|---:|
| 50 MeV shape and 50 MeV direct abundance, fractional P change | −0.4961% | −0.8867% |
| 100 MeV shape at the 50 MeV abundance | −0.4733% | −0.8746% |
| 100 MeV shape and its own direct abundance | −8.3986% | −14.9776% |

At fixed density the 100 MeV spectrum gives *less*
suppression in this compact setup. Its own direct
yield, about 18.35 times the 50 MeV yield, dominates
that spectral difference in this exercise. At k=10
the shape-only difference is +.0001226 in ln P; changing
the abundance then contributes −.15347. The 100 MeV
trace alone differs from the 50 MeV trace by just about
8e−9 in ln P at k=10. The two 50 MeV density references (PSD fold and
independent abundance fold) differ by 1.90e−6 relative;
the resulting correction is negligible at the scale of
the grid differences below.

The denser 136-q/l40 native tier at the same q cap
gives −0.8781% and −0.8459% at the same reference
abundance and k=10, and −14.5233% at 100 MeV with
its own abundance. At k=2 the corresponding three
estimates are −0.4960%, −0.4543% and −8.0732%.
This supports the sign and hierarchy of count and shape
effects, but the shift in the small shape difference
itself (maximum 2.02e−4 in ln P) is not converged.
For example, the 100 MeV spectrum's native mean q
was 9.73% too high on the linear 78-point grid and
3.76% too high on the linear 136-point grid when both
used cap q=45. An independent numerical fold of the
same cubic spline reproduces the logged early ΔNeff.
The narrower-cap and transformed grids below
address the main momentum error in the broad grid.

## Momentum-resolved comparison

Two additional native configurations, a linear grid
with cap 18 and a transformed 100-point grid, agree
to within 1.49e−5 in the matched-density log-power
responses across all sampled modes and both redshifts.
The maximum discrepancy in the 100 MeV case's own
larger abundance response is 2.73e−4 in ln P, with the
worst agreement at the smallest k. These differences
are numerical checks, not error bars on the physical model.

| z=5, transformed 100-point q grid | k=2 h/Mpc | k=10 h/Mpc |
|---|---:|---:|
| 50 MeV shape, its own direct abundance | −0.4983% | −0.8808% |
| 100 MeV shape, at the 50 MeV abundance | −0.4466% | −0.8277% |
| 100 MeV shape, its own direct abundance | −7.9403% | −14.2333% |

The matched-density shape difference is now
+5.36e−4 in ln P at k=10, larger than the wide-grid
estimates but much smaller than the abundance effect.
The lower cap removes less than 1.25e−7 of the
prepared q² moment; its biggest mean-q or mean-q²
error is 0.251%. The transformed grid's equivalent
moment errors are below 3e−7. The background
relativistic growth of energy density, expressed as
ΔNeff, then agrees with the same momentum sums.
A still finer hierarchy cross can test residual error
in the transfer, which moment agreement alone cannot certify.

[Linear narrow-cap audit](results/warm-direct-vacuum-temperature-native-lowcap-2026-11-26.json),
[raw](results/warm-direct-vacuum-temperature-native-lowcap-raw-2026-11-26.json.gz);
[transformed audit](results/warm-direct-vacuum-temperature-native-transform-2026-11-26.json),
[raw](results/warm-direct-vacuum-temperature-native-transform-raw-2026-11-26.json.gz).

These are conditional linear transfer comparisons of
the chosen, already produced spectra. The 100 MeV
case does not predict the cold complement or the true
late-time density, nonlinear clustering, an observational
exclusion, or a full thermal hadronic source.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
# Prepare varied-spec.json and the native inputs as in the shapes note.
python -W error warm_direct_vacuum_temperature_native_tiers_reference.py varied dense linear136_cap45
python -W error warm_direct_vacuum_temperature_native_tiers_reference.py varied narrow linear136_cap18
python -W error warm_direct_vacuum_temperature_native_tiers_reference.py varied mapped transform100
# Run the five configurations on each tier. For example:
python -W error warm_direct_vacuum_temperature_native_audit_reference.py varied-spec.json varied 50 100
~~~

[Compact exact-mode audit](results/warm-direct-vacuum-temperature-native-compact-2026-11-26.json);
[compact raw](results/warm-direct-vacuum-temperature-native-compact-raw-2026-11-26.json.gz);
[dense audit](results/warm-direct-vacuum-temperature-native-dense-2026-11-26.json);
[dense raw](results/warm-direct-vacuum-temperature-native-dense-raw-2026-11-26.json.gz);
[separate grid moment diagnostic](results/warm-direct-vacuum-temperature-native-grid-moments-2026-11-26.json).
