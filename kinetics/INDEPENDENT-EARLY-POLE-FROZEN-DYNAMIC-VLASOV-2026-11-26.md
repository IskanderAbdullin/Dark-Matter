# Dynamic collisionless response from frozen early spectra

The [static screening test](INDEPENDENT-EARLY-POLE-SNAPSHOT-HARD-LOOP-SCREENING-2026-11-26.md)
does not by itself determine a soft scattering correction. Here both
retarded and statistical terms are evaluated from the archived,
isotropic parent and daughter shapes at x=.05. The collisionless
longitudinal and transverse angular functions follow the relativistic
Vlasov construction (compare [Kneller](https://www.cambridge.org/core/journals/journal-of-plasma-physics/article/abs/relativistic-kinetic-theory-of-waves-in-isotropic-plasmas/3266B845EA48C0D431C04A5DAB233D6F)).
Each massive shell uses z=(omega/k)/v and the angular Landau cut
exists only when |z|<1. Its retarded radial weight is p E(−df/dp),
while fluctuations are reconstructed from independent occupation
variance f(1−f). No common-temperature fluctuation relation is imposed.
The massive response recovers the earlier static susceptibility
exactly for the declared cells. Smearing the edge derivatives inside
a fraction of the adjacent cell is a numerical comparison, not
kinetic broadening. The two raw momentum grids and three edge
smoothings are all shown.

**Spacelike soft kick comparison.** For a hard test line and the
spherical cutoff labelled 20 in the previous two-thermal convention:

| Match the upper cutoff by | Radial second moment / thermal | Transverse / thermal | Finite local phase / thermal |
| --- | ---: | ---: | ---: |
| Same k/T_ref | .5695–.5707 | .5664–.5676 | .5993–.6005 |
| Same k/m_D for each state | .5050–.5063 | .5050–.5063 | .5846–.5859 |

The phase is at one inverse angular half-span, for g=.005 in the
same reference geometry, and is the common-kick/local Doppler
quantity, not an additional effect. The two integration resolutions
change the dimensionless integrals by less than 3e−6 absolute.
For these particular integrated soft observables, a massless model
with both the actual static m_D and low-phase-velocity noise ratio
matched agrees with the massive result to within .06% relative.
This does not validate a single equilibrium temperature, the light-cone
dispersion, or a hard-transfer completion.

**Timelike real dispersion only.** The low-momentum plasma
frequency is about .333 g T_ref instead of .471 g T_ref, matching
the independent velocity-weighted static moment. In the two massive
raw snapshots the longitudinal timelike root reaches the light cone
near k/(g T_ref)=1.152–1.153 in this collisionless approximation;
the massless thermal control can approach the light cone arbitrarily
closely at higher momentum. The transverse pole is still present.
The absence of an identified timelike root beyond the endpoint does
not rule out a damped spacelike response. We do not assign timelike
collective modes a statistical occupation from the frozen Landau noise.

No scalar-production correction, equilibrated bath assumption, hard
scattering, or chemical-threshold revision is inferred. Reproduce in
kinetics/ with early_pole_frozen_dynamic_vlasov_reference.py,
early_pole_frozen_vlasov_soft_kick_reference.py (--nodes 8 or 12) and
early_pole_frozen_vlasov_dispersion_reference.py.
Raw output: results/early-pole-frozen-dynamic-vlasov-2026-11-26.json.gz.
