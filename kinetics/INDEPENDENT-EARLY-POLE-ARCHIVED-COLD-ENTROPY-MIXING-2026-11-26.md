# Distinguish intra-species reshaping from interspecies heat exchange

Earlier heat-exchange rows used fitted common or separate species
temperatures and imposed exchange rates. In the
[matched cold-window moment model](INDEPENDENT-EARLY-POLE-ARCHIVED-SHARED-COLD-WINDOW-2026-11-26.md),
construct a different, explicitly conservative comparison:
a projected entropy gradient on the *fermion* cell moments.
The positive block mobility in each cell is the local Fermi
occupation covariance of count and centered energy. Project out
the two species counts and either a single total energy gradient
(joint mixing) or one variational energy constraint for each
species separately (separate mixing). The same gauge energy used by
the collision solver defines these constraints. Photons and scalar
counts are untouched by the added reshaping operator.

For an entropy gradient g and constraint matrix C, the driver is
g−C(CᵀMC)⁻¹CᵀMg, using the cell mobility M. Its rate-scaled
quadratic form is nonnegative. This operator changes the hard
momentum shape even when it cannot exchange heat between species.
In the joint case individual energy changes are equal and opposite
at an instant. A manufactured Fermi fixture with fixed unequal
chemical potentials but common variational temperature gives
negligible response in both modes; choosing distinct temperatures
at the same chemical potentials activates only the joint mode.
The projection leaves both number constraints and the appropriate
energy constraint at numerical precision.

At the imposed rate 500 and reference exposure .0003, relative
scalar suppression by pair dynamics is:

| Archive | No extra mixing | Separate conserved-energy reshaping | Joint total-energy reshaping |
| --- | ---: | ---: | ---: |
| Fine | 12.23% | 8.90% | 9.33% |
| Third | 14.35% | 8.52% | 10.00% |

Each ratio uses a decay-only control with that same reshaping
operator. In this finite example intra-species reshaping alone
already changes the response materially. The additional difference
between joint and separate modes cannot be read as a measured
microscopic heat exchange. The rate, mobility and live hard
window are chosen for the comparison, not derived from the gauge
collision kernel. Absolute energy/count and entropy ledgers are
checked in the evolved test. No full production rate, physical
freeze-out or CLASS input follows.

Reproduce with early_pole_archived_cold_entropy_mixing_reference.py
--grid fine/third --mode separate/joint --rate 500.
Raw: results/early-pole-archived-cold-entropy-mixing-2026-11-26.json.gz.
