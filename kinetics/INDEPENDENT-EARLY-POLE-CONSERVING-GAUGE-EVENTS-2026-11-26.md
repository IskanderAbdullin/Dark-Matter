# Reciprocal gauge energy in reversible collision and expansion tests

The [reciprocal energy functional](INDEPENDENT-EARLY-POLE-FROZEN-MASS-GATE-2026-11-26.md)
assigns consistent first-order hard quasiparticle shifts without
counting the interaction energy twice. Here we use its gradient in
three discrete decay/inverse and two pair channels. The outgoing
energy is reconstructed on two adjacent daughter nodes with
nonnegative weights. The event affinity is evaluated from the same
weighted occupancies and quantum statistical factors. These arbitrary
reference channels demonstrate energy and detailed balance; they do
not reproduce the exact momentum geometry or a gauge scattering rate.

At fixed scale, energy is quadratic in the nodal counts n, and an
instantaneous event of direction v exactly conserves it when its
gradient e(n) satisfies e(n)·v=0. For a finite step of size h, use the
gradient at n+h v/2. With the corresponding adjusted node weights,
this is the exact quadratic difference, not an approximate correction.
The reversible dressed equilibrium residual is at roundoff. At g=.005
the fixed-scale energy drift is −8.46e−9 using bare event energies,
versus +1.7e−16 with gradient-based channels. Both conserve the
auxiliary count N_gamma+N_parent+N_daughter; this is not electric
charge. The calculated entropy change agrees with the integrated
nonnegative entropy production.

In an expanding tagged hard support the single-particle momentum
redshifts and the field energy uses physical moments J(a,n).
The trace comes from an independent fixed-comoving-state scale
derivative. At g=.005 the continuity residual is −5.25e−9 for bare
events and +1.4e−16 for the gradient-based construction, with the
same initial distribution. The bare residual agrees with the
independently integrated spurious collisional energy source.
The expansion interval and event strengths are illustrative.

**Moving from energy to angle.** For fixed parent p and scalar q,
the outgoing magnitude k satisfies the exact gradient-energy equation
and the triangle interval |p−q|≤k≤p+q. For the monotone hard daughter
dispersion E_0+A/E_0, solving the quadratic for E_0 gives k directly.
The exact angular Jacobian k/(p q v_daughter) agrees with a finite
difference at better than 1.6e−7 relative. At the archived x=.05 and
g=.005, using the square-root shifted-mass proxy instead of the
variational dispersion moves the aperture boundaries by less than
2e−7 of a half-span under the same moments. This validates their
kinematic agreement through first order in this hard window, not the
matrix element, angular collision mesh or in-medium width.

In the archived early spectrum the estimated interaction energy is
4.18e−7 T_ref⁴ with no extra hard cut. It is negligible relative to
the reference visible radiation density (1.19e−8), but is about 3.7%
of the represented scalar energy. Restricting the reference moments
to p/T_ref≥.5 gives about 2.7% of initial scalar energy. The visible
Hubble comparison cannot substitute for auditing the dark energy
budget. The cutoff implementation still has to be matched to a
physical collision calculation; the tagged expansion witness does
not move its cutoff within the represented phase space.

Reproduce in kinetics/ with early_pole_conserving_gauge_events_reference.py,
early_pole_conserving_gauge_frw_events_reference.py,
early_pole_gauge_gradient_decay_geometry_reference.py, and
early_pole_snapshot_gauge_ledger_reference.py. Raw output:
results/early-pole-conserving-gauge-events-2026-11-26.json.gz.
