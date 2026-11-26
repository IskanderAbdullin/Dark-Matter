# Use the time dependence of the signed losses before debiting them

The [first signed comparison](INDEPENDENT-WARM-SIGNED-PAIR-FERMION-CAP-2026-11-26.md)
does not cap scalar occupation, but debits all linear-collision
losses pessimistically at the end. It therefore becomes
inconclusive for a moderately large fermion cap. Instead
retain the time-dependent number in every loss term. The
signed quadratic Bose term is still nonpositive if each
fermion occupation per spin in both masses and charges is
at most 1/2. No condition is imposed on the pointwise
scalar occupation; finite total energy and the other network
and expansion assumptions are unchanged.

With X the comoving number of all fermions and Z the
cumulative high-invariant charged-lepton loss, the comparison
for `Y=X+Z` has three monotone losses in Y: ordinary
pair losses proportional to `Y²`, soft single-Bose losses
to `Y^(2/3)`, and hard single-Bose losses to `Y^(4/3)`.
These powers follow from the Fermi-ball radial bound;
inverse gains and charged-lepton injection cannot worsen
the lower comparison. The bounded visible pair energy and
elastic heat are included. The corresponding integrals
are weighted by `b⁻³/H`, `b⁻⁵/H` and `b⁻⁴/H`.
The maximum possible high-invariant extraction is charged
once, after the solution for Y.

Each finite time cell has explicit upper integrals from a
lower radiation/then photon clock. The scalar comparison
trajectory is decreasing. Its loss over the cell is no
larger than using its entrance value in every rate.
Starting at the previous certified lower value, ordering
of the one-dimensional comparison gives a new lower value
by subtracting these endpoint rates. This is not an
error estimate for a physical Euler evolution. Intervals
are refined; a tiny numerical padding is applied, not
machine-directed interval rounding.

| q32, imposed per-spin fermion cap | Fermion Ωh² lower, 2,800 cells |
| --- | ---: |
| .001 | 22.8874 |
| .01 | 19.0954 |
| .10 | 8.9077 |
| .50 | 1.1192 |

At cap 1/2 the corresponding coarser bounds are .8182 and
1.0509, so the lower result already exceeds total .120 at
each grid. The refined bound therefore excludes .120
throughout the **entire signed-cap range tested**, while
allowing arbitrary scalar occupation. It does not
establish the fermion cap dynamically or extend the sign
argument above 1/2. Other hadronic or unknown reactions and
a different expansion history are outside this comparison.

~~~sh
export OPENBLAS_NUM_THREADS=1
D=results
python -W error warm_signed_fermion_time_comparison_reference.py \
 "$D"/warm-isolated-fermion-survival-bound-2026-11-26.json \
 "$D"/warm-fermion-survival-with-leptons-2026-11-26.json \
 "$D"/warm-isolated-recycling-energy-envelope-2026-11-26.json \
 "$D"/warm-signed-pair-fermion-cap-2026-11-26.json
~~~
