# Elastic redistribution can reshape the peak without stopping production in these snapshots

Two independent counterfactuals separate redistribution of the
scalar from the simultaneous decay and expansion equations.

**Frozen screened shape evolution.** Apply the
[screened cold graph](INDEPENDENT-EARLY-POLE-SCALAR-VERTEX-RESPONSE-2026-11-26.md)
at fixed mass, fixed screening weights and a closed cold momentum
window, with no decays or cosmic expansion during the isolated
evolution. This window retains 99.9989% of the scalar count in
these inputs. Evolve for three units of the initial formal inverse
absolute flux; the chosen interval is only a normalization, not
a physical elapsed time. The 0.06 third-grid example moves about
50% of the number in total variation, while increasing discrete
entropy and preserving number and kinetic energy at about 10⁻¹⁵.
The previous momentum grid moves about 36–38%; the shape change
is not individually mesh-converged. Then evaluate the original
decay operator at the replacement distribution, with fermions
held at their snapshot values. On the third grid the signed
growth rate falls from 2.189 to about 1.948 per e-fold; on the
previous grid it changes only from 4.458 to about 4.460. The
self-consistent Hartree pole shifts by less than 4 × 10⁻¹¹
relatively for these replacements.

**Moment-constrained extremal check.** Hold the parent, fermion
daughter, scalar pole and epoch fixed. The total signed decay
source is then affine in the individual scalar occupations. Fix
the scalar number and kinetic energy in the represented IR band,
leaving occupations outside it untouched; optimize over every
nonnegative distribution in the band. An extremum can be
represented on at most two nodes. Enumerate those vertices and
cross-check with a separate linear program. On the third grid
the full represented q≤10⁻⁵ check for the 0.06 case at x=0.05
gives a net rate in [1.564, 2.262] per e-fold. At x=0.10, the
corresponding range is [3.117, 5.771]; restricting the lower
edge changes the upper endpoint but not the positive minimum.
All six snapshots and three reported lower cuts retain a
positive minimum. Re-solving the scalar gap at each selected
extreme gives effectively the same source in the tested cases.

This is a conditional non-freeze-out test under fixed fermion
states; it does not bound what happens after those states or
the total scalar particle count change. The dynamical screened
case is not a calibrated N=1 time evolution.

~~~sh
python -W error early_pole_scalar_screened_redistribution_reference.py results/early-pole-spectrum-third-m06-2026-11-26.json.gz --angle 12
python -W error early_pole_scalar_source_moment_envelope_reference.py --lower 0 results/early-pole-spectrum-third-m06-2026-11-26.json.gz
~~~

[Redistribution and source readout](results/early-pole-scalar-screened-redistribution-2026-11-26.json);
[moment envelopes](results/early-pole-scalar-source-moment-envelopes-2026-11-26.json.gz).
