# A distinct lower bound on the warm fermion residue

The [scalar recycling upper](INDEPENDENT-WARM-ISOLATED-RECYCLING-ENVELOPE-2026-11-26.md)
does not test whether the **total** dark rest abundance can fall
from the warm 10.442 MeV snapshots to .120. This comparison
uses a different assumption: all outgoing scalar occupations
must be bounded by one uniform number `F_B` at all later times.
The final-state enhancement of a scalar pair is then no more
than `(1+F_B)^2`. No such cap is inferred for the full interacting
solution from the earlier ballistic forward-only check.

Let `X` be the combined comoving light and heavy **fermion**
count of both charges. A one-to-two heavy decay or reverse
fusion preserves X, as does two-to-two conversion between
the fermions. Elastic redistribution and reentry of all known
cohorts are allowed. In the stated isolated network, only
annihilation of a particle-antiparticle pair to scalars lowers
X, and any reverse gain only helps its lower bound. Visible
annihilation, a new loss process, and external heating or
injection are not included.

For light or heavy annihilation use the vacuum spin-trace
Cauchy ceiling `K_i=y⁴/(8π m_exchange²)` (the pair produces
two identical scalars). The largest coefficient is for heavy
annihilation. Since `2 n_+ n_-≤n_i²/2`, including arbitrary
charge ratios, the simple comparison obeys

`X_late ≥ X₀/[1 + K_max(1+F_B)² X₀ I/2]`,

where `I` is the conservative expansion integral of
`b⁻³/H` in d ln b from the earlier visible clock. Here X is
normalized as a number density at the anchor. The core
fermion inventory underestimates the fuller starting count.

An energy-assisted correction also imposes
`m₂ X + (m₁−m₂) X_heavy ≤ ρ_anchor,max`, for the comoving
energy bound under nonnegative pressure. At fixed X the
quadratic annihilation functional is convex in the heavy
fraction, so maximize it at either admissible endpoint.
Integrate its resulting worst loss; no assumed thermal shape
or known heavy abundance is used. The earlier archive gives
the upper starting energy including all known exterior tails,
which must not have been fed by a new external source.

| q32 anchor, uniform F_B | Fermion Ωh² lower, simple rate | Lower with energy constraint |
| --- | ---: | ---: |
| 0 | .99215 | 1.00620 |
| .1 | .82163 | .83125 |
| 1 | .25025 | .25113 |
| 2 | .11140 | .11158 |

Use the smaller stable mass for every surviving fermion, so
these are lower bounds on the fermion rest density even if
some heavy species remain late. To permit a total .120 in
this isolated comparison, the uniform cap would have to allow
at least about 1.89; a cap of 1 already rules out the target
total, regardless of the putative scalar fraction. A larger
allowed F_B is only loss of this exclusion. It does not show
that such occupancy is produced or that depletion succeeds.

~~~sh
export OPENBLAS_NUM_THREADS=1
D=results
python -W error warm_isolated_fermion_survival_bound_reference.py \
  "$D"/warm-isolated-recycling-energy-envelope-2026-11-26.json \
  "$D"/warm-known-exterior-recycling-envelope-2026-11-26.json \
  "$D"/warm-early-support-time-16-32-interaction-2026-11-26.json \
  "$D"/warm-q{16,32}-10p442-heavy-input-shape-2026-11-26.json
~~~
