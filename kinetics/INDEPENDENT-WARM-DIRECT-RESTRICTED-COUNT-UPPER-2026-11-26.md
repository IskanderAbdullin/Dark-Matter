# A weighted particle-count ceiling for the selected direct vacuum network

A parent produced by a visible reaction is worth two
units of final scalar-or-stable count. Define the
comoving ledger as Q = N_scalar + N_light + 2 N_heavy,
counting both charges for each fermion. Decay of a heavy
fermion adds one scalar and one light fermion while
losing one heavy fermion, so Q is unchanged. This also
holds for inverse decay, quantum factors, and any
elastic momentum exchange. Stable-pair ↔ scalar-pair
conversion preserves Q. Heavy-pair loss to light or
scalar pairs lowers Q, but its dark **inverse** would
raise Q and is not included in this upper comparison.
Each heavy pair created from the visible bath increases
Q by four; each direct light pair increases it by two.
Inverse visible reactions only lower it. This does not
require instantaneous decay or assume any rate for the
two-to-two stable conversion.

Start empty and enclose both positive visible vacuum
sources through .1 MeV in the prescribed radiation
history. The [direct production bracket](INDEPENDENT-WARM-DIRECT-TEMPERATURE-BRACKET-2026-11-26.md)
is repeated at the existing heavy mass and independently
with the light mass in the pair threshold and kinematic
kernel. Use the same chosen vector coupling and two
charged-lepton/charged-pion sources for both, including
the KLOE bins plus the stated marginal allowance, separate
low positive dispersive moments, and R_pi≤1 above the
measured interval. Exclude other hadronic or thermal
channels. Because N_scalar≤Q, the final scalar rest
density cannot exceed m_scalar times the entire ledger.
The result grants all conversion even if the physical
rates make it impossible.

| Reheating (MeV) | Both source ceilings with all Q as scalar, Ωh² | Fraction of .060 |
|---:|---:|---:|
| 25 | 0.000036644 | 0.0611% |
| 50 | 0.00074844 | 1.247% |
| 75 | 0.00387454 | 6.458% |
| 100 | 0.01287463 | 21.46% |

At 100 MeV the parent counts contribute at most
0.00839110 and the direct-light count at most .00448353
in scalar units, including the charged-pion allowance
for the latter. The missing weighted count is at least
0.0471254 in scalar abundance units (about 1.99e−3 per
visible entropy) before one can reach .060. This is
stronger than just checking prompt decay but narrower
than a full kinetic-model result. A hidden number
source, heavy-pair regeneration, other hadrons, or
a different expansion/production history can violate
the stated ceiling. The printed figures are numerical
positive enclosures, not a confidence interval or an
observational exclusion.

The [separate co-production tally](INDEPENDENT-WARM-DIRECT-COPRODUCTION-BUDGET-2026-11-26.md)
shows how much stable mass survives if counts are held;
its central calculation cannot be used as the upper
for direct light production.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
D=results
python -W error warm_direct_restricted_count_upper_reference.py EOS.late.dat KLOExSec.data \
  "$D/warm-late-cold-pion-dispersive-recoil-2026-11-26.json" \
  "$D/warm-direct-vacuum-temperature-bracket-2026-11-26.json"
~~~

[Both invariant grids and all source terms](results/warm-direct-restricted-count-upper-2026-11-26.json).
