# How does the 50 MeV upper comparison respond to its couplings?

Vary the vacuum mixing parameter ε and the selected dark
Yukawa y while fixing masses, visible entropy/expansion, the
mediator width, and the list of reactions. Start with the
[from-empty lepton-only calculation](INDEPENDENT-WARM-FROM-EMPTY-RECYCLING-ENVELOPE-2026-11-26.md)
and the [charged-pion variant with the wider stated caps](INDEPENDENT-WARM-PION-STATISTICS-AND-THRESHOLD-2026-11-26.md).
No other hadronic source is added. Resolve direct
pair/decay counts separately from recycling. The source
and the charged elastic bounds scale as ε² in this
fixed-width setup. The energy-dependent recycling kernel
scales as y⁴. The prompt credit remains an upper credit at
every y; it does not presume that a decay is actually fast.

Write D, R and J for the baseline direct Ω proxy, recycling
Ω upper, and heat exponent. With `x=(ε/ε₀)²` and
`j=|y/y₀|`, the upper comparison is

`D x + R x² j⁴ exp[2J(x−1)]`.

The second x comes from the quadratic energy envelope.
The exponential adjusts the heating already included
in R. Recompute the measured-band direct and recycling
pieces from the existing code and recover its archived
baseline. Solve the half-of-.120 comparison independently
by a monotone numerical root and an algebraic y root.

| At 50 MeV | Baseline scalar Ωh² upper | ε/ε₀ where the upper first reaches .060, fixed y | y/y₀ where the upper first reaches .060, fixed ε |
| --- | ---: | ---: | ---: |
| Vacuum charged leptons | .02154169 | 1.29197 | 1.29751 |
| Same plus selected pions and stated caps | .02869973 | 1.20108 | 1.20595 |

References are ε₀=10⁻⁵ and y₀=10⁻³. Crossing this boundary
only removes this upper-bound exclusion. It is not evidence
of a particle population making half the target, a fit,
or a bound on the actual values of ε or y. Varying the
mediator width, temperature history, other vertices,
unmeasured sources, or the full cosmology requires
a separate run; the Yukawa rescaling only applies to the
specified tree-level vacuum heavy-forming channels.

~~~sh
export OPENBLAS_NUM_THREADS=1
D=results
python -W error warm_recycling_vacuum_coupling_sensitivity_reference.py EOS.dat KLOExSec.data \
  "$D"/warm-energy-weighted-recycling-kernels-2026-11-26.json \
  "$D"/warm-from-empty-recycling-envelope-2026-11-26.json \
  "$D"/warm-pion-bose-and-threshold-2026-11-26.json
~~~
