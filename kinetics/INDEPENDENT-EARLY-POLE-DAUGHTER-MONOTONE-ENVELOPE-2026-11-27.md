# Per-cell monotonicity is sufficient for a negative frozen contrast in this diagnostic

The [unrestricted fixed-moment envelope](INDEPENDENT-EARLY-POLE-DAUGHTER-FROZEN-MOMENT-ENVELOPE-2026-11-27.md)
allows both signs. As a separate shape prior, require the daughter
occupation to be nonincreasing **within each selected hard cell**;
allow jumps between cells. Do not assume the time-evolved physical
distribution meets this condition.

The native projected shape is constant within each selected cell.
Independently verify that its stored first two moments match that
constant to within 2.22e-18. With exact number and energy moments
of a constant, any monotone occupation must be constant almost
everywhere: the covariance of a strictly increasing energy and a
nonincreasing occupation has one sign, and it vanishes only for a
constant. This follows by writing the covariance as half the
average product of differences on two independent momenta.
The statement also holds with the opposite monotonic direction
if applied to a native cell separately.

A nonincreasing occupation is a positive mixture of cutoff
steps (1_{ple t}), with the zero and full-cell cases at the
two endpoints. The source depends linearly on that mixture.
Numerically optimize the retained moments over cutoffs and
bound the continuity between sampled cutoffs using the integrated
positive/negative parts of the dual residual density. The sign
check below also relaxes each native moment by (10^{-6}) of
the cell measure while keeping the local arm's targets fixed;
this is much larger than stored native moment discrepancies.
The reference fitted rate is reproduced independently within
3.23e-13, and the largest change in a dual bound between
32 and 64 quadrature nodes per interval is 5.61e-7.

| Hard depth × cold modes | Numerical upper endpoint, two moments | Also hold the reconstructed reciprocal-energy moment |
|---|---:|---:|
| 1 × 8 | -0.0010952 | -0.0031299 |
| 1 × 32 | -0.0049733 | -0.0057602 |
| 4 × 8 | -0.0113737 | -0.0115886 |
| 4 × 32 | -0.0119602 | -0.0121604 |

All tested upper estimates have the negative sign, even without the
extra reciprocal moment. This is a restriction on **shape** under fixed
initial geometry, not proof that the physical occupation is monotone or
that the later response is negative. The analytic cutoff-mixture
inequality bounds exact integrals if supplied; the reported envelopes
are numerical quadratures rather than interval-certified values.

Use `early_pole_daughter_monotone_envelope.py` on the four frozen
fixtures at `--subdivision 12` and 48, both without and with
`--fix-gauge`. Cross-check with
`early_pole_daughter_monotone_envelope_compare.py` against the
previous fitted-rate fixture.
