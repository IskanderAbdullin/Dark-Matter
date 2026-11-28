# Two-mass geometry and conservative deposition below absolute-energy precision

The upstream [fast area reference](https://github.com/IskanderAbdullin/Dark-Matter/blob/f77e75a43e99f03c5ae4dc1d732731bbb9ba3fca/kinetics/TWO-MASS-AREA-2026-11-25.md)
uses absolute daughter energies. It explicitly reports zero-width cells
when their energy edges round to the same float64 value. This independent
control keeps momentum edges, daughter rest masses and the parent rest
excess δ=M−m−b as separate inputs. It tests whether the geometric and
conservative bookkeeping is still recoverable before any occupation or
collision rate is introduced. An excess already lost when rounding M
cannot be recovered and must be supplied separately.

## Geometry with two trace moments

For daughter rapidities u and t, the on-shell region satisfies
|u−t|≤Y≤u+t, where cosh(Y)−1=δ[2(m+b)+δ]/(2mb). At fixed t the
allowed u interval is [|Y−t|,Y+t]. Represent daughter kinetic offsets
as K=2m sinh²(u/2), split clipped rectangles at the rapidity
boundaries, and integrate energy differences via a factored
difference of hyperbolic cosines.

The reference records both geometrically weighted trace moments
m²/E and b²/E. It also records their deficits from the corresponding
rest energies. The scalar inner primitive uses log1p; a separate
series avoids cancellation in the trace deficit when K≪m.
All reported moments are normalized by the *whole* cell's energy
rectangle. The tests check the trace identities, daughter exchange,
and cell subdivision. An 85-digit adaptive reference uses the
unreduced absolute-energy formula E± on the same cells.

In the near-threshold test both float64 absolute energy widths vanish,
yet the underlying allowed area is .70875063 of the true small
rectangle. The two kinetic widths and both positive trace deficits
remain measurable in the offset formulation. The separate asymmetric
test also has a collapsed scalar width and a resolved fermion width.
These fractions do not mean the absolute phase-space measure is large.
At ordinary widths, 96 sample cells (23 partial) agree with the
pinned upstream area to 1.68e−14 in absolute area fraction. The
maximum relative discrepancy against the 85-digit reference in
the targeted suite is 2.70e−11.

## Positive event projection on momentum shells

Compute the boosted daughter three-momenta first, then recover each
kinetic offset as p²/(E+m). The daughter offsets must sum to the
parent kinetic energy plus its supplied rest excess. Deposit each
offset between its two neighboring shell-mean offsets using
nonnegative weights; their product preserves one count per daughter
and the represented event energy. If an event lies outside the convex
hull of neighboring shell representatives, reject it instead of
silently extrapolating or clipping.

Trace is not a per-event collision invariant. Its projection bias is
tracked separately, using K²/(m+K), since the constant and linear
parts of m²/(m+K) are already reproduced. Convexity implies a
nonnegative bias; it is bounded by the sum of the interpolation
variance and weighted within-shell variances, divided by m. The
three examples (ordinary, extremely small parent excess, and a
backward almost-at-rest daughter) check both the bound and
independent high-precision boosted energies.

This is a numerical reference for two positive on-shell masses,
not a replacement for the vectorized production kernel, a
thermal dispersion, an FRW collision update, or a derived relic
abundance. It does not impose conservation of the interaction trace.

Replay from the kinetics directory with NumPy, SciPy and mpmath:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error two_mass_kinetic_offset_trace_reference.py --check-upstream
python -W error two_mass_kinetic_offset_deposition_reference.py
```
Use the pinned upstream `two_mass_area.py` for the optional ordinary
comparison; omit `--check-upstream` to run without that file.
[Geometric receipt](results/two-mass-kinetic-offset-trace-2026-11-27.json)
and [deposition receipt](results/two-mass-kinetic-offset-deposition-2026-11-27.json).
