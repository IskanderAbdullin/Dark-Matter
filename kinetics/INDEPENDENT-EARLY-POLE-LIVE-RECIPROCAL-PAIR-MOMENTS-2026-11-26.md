# Live reciprocal hard photon-pair moments

The [evolving decay cells](INDEPENDENT-EARLY-POLE-LIVE-PAULI-GAUGE-MOMENT-DECAY-2026-11-26.md)
keep the photon moment fixed. This separate benchmark activates
the photon cells and both reversible gamma-gamma to Dirac-pair channels.
The hard photon distribution in each cell is reconstructed with a
Bose exponent a+b u, constrained by a>|b| on the full cell.
The reciprocal gauge energy derivative follows through the
Bose and Fermi fits together. The independent photon test checks
60 representable occupancies on four hard intervals; the gauge-moment
derivative and entropy derivative agree with finite differences.
Near-condensation states or arbitrary singular endpoint deposits
are not guaranteed to be represented by this two-moment closure.

A frozen pair geometry test uses the actual archived charged
populations and a held photon bath. The first-order dressed energies
define the allowed p,k interval at the fixed incoming energy and
total momentum. The angular Jacobian uses the outgoing group velocity.
The g=0 aperture agrees with the independent vacuum boost formula,
including the outgoing hard-support cut. The deposited outgoing
count, energy and vector momentum match the incident values to
roundoff. In a separate common-T test, choosing the same auxiliary
chemical potential for photon and charged populations gives
quantum forward/inverse equality at each sampled event to 1.5e−15.

The live calculation uses two sample points in incoming photon
momentum and angle, rather than integrating the full incident
spectrum. Outgoing momenta span the continuous on-shell kinematic
interval over the hard cells. At g=.005, with the same initialized
distributions, the dressed energy drift is zero to rounding,
versus +3.94e−9 when the gauge energy is stored but bare events are
used. The latter drift is explained by the integrated collisional
work. The combined photon plus fermion count holds to roundoff and
the independently integrated entropy increases. The photon exponent
stays above .6228, and the represented Fermi occupations remain
between .0102 and .2288. Refining 10 to 20 quadrature points changes
the combined outgoing count by far below 1e−11.

This is not a full collision rate. The calculation omits the
spin-dependent matrix element and the full incident integration,
and does not yet evolve decays, pair production and expansion
simultaneously. Its closest use is as a conservation and
reversibility validator for a future combined hard sector calculation.

Reproduce from kinetics/ with early_pole_hard_bose_cell_reconstruction_reference.py,
early_pole_gauge_pair_subcell_kinematics_reference.py and
early_pole_live_reciprocal_pair_moments_reference.py.
Raw: results/early-pole-live-reciprocal-pair-moments-2026-11-26.json.gz.
