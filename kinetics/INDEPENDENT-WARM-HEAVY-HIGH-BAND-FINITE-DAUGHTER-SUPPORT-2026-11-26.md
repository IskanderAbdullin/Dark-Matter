# Bound the high heavy band if the two dark daughter supports stay finite

The [finite-band heavy estimate](INDEPENDENT-WARM-HEAVY-COMOVING-BAND-ENERGY-2026-11-26.md)
bounds energy only below a comoving cut. Test a deliberately
restricted branch in which the two incoming dark species,
light fermion and scalar, stay below an equal, fixed
comoving momentum P throughout the post-10.442 MeV epoch.
Use the archived **retained** momentum edge as a candidate
for P, not as proof of the physical support. The actual
positive exterior, if allowed to feed reactions, cannot
silently use this assumption. Continue to assume the
heavy states initially empty and no unbounded elastic
transfer of heavy energy or momentum.

At expansion b, even ignoring other positive-energy
outgoing particles, a heavy particle created from
two incident masses (m_a,m_b) can have only
(b² p_H²≤2P² + (m_a²+m_b²−M²)b² +
2 sqrt[(P²+m_a²b²)(P²+m_b²b²)]).
This ceiling is concave in b², with negative derivative
already at the anchor for all three required initial
pairs. It therefore falls. For q32 the largest internal
cap is .60055 GeV, below a 1 GeV comoving cut.
Dark heavy production supplies no population outside
that cut in this restricted branch.

Thermal visible production is **not** given a hard
momentum cutoff. Reverse a tagged heavy pair via
Boltzmann detailed balance. Integrate the tagged
density for (p≥1/b) GeV against the exact
exponential moments of (p²(M+p)). Below invariant
mass 1 GeV use the subpole bound and supported
charged pions; above it keep half the Boltzmann
exponent to pay the invariant threshold and integrate
the rest at twice the bath temperature. The expansion
floor bounds the two time integrals. Drop heavy
decay in this high band, which only increases the
injection-based energy allowance.

For q32, the high band can gain at most
(2.596×10^{-27}) GeV⁴ in anchor-normalized
comoving energy through those visible paths. Combining
with the independently integrated lower band gives
a heavy energy upper of 2.8123% of the full budget.
This falls below the [4.8336% direct-light test criterion](INDEPENDENT-WARM-LIGHT-DIRECT-PAIR-DOSE-2026-11-26.md)
**under the stated support hypothesis**. It does not
prove a uniform light cap, bound gains from light
one↔two redistribution or elastic collisions, or
apply if the dark daughters populate the missing
higher momentum modes. The cosmological result
is not inferred from truncating a grid.

~~~sh
export OPENBLAS_NUM_THREADS=1
D=results
python -W error warm_heavy_high_band_finite_daughter_support_reference.py EOS.dat \
 "$D"/warm-signed-fermion-with-charged-pions-2026-11-26.json \
 "$D"/warm-isolated-recycling-energy-envelope-2026-11-26.json \
 "$D"/warm-heavy-band-aggregate-energy-bound-2026-11-26.json.gz \
 "$D"/warm-light-direct-pair-dose-2026-11-26.json \
 "$D"/warm-q{16,32}-10p442-heavy-input-shape-2026-11-26.json
~~~
