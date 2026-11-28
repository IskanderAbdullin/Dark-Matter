# A cold pair can redistribute within the infrared without making a hot pair

For two incoming scalar quasiparticles at the stored Hartree mass,
the two outgoing energies at fixed initial angle are
half the sum of the incoming energies plus or minus
the lab boost term times the center-of-mass scattering cosine.
This gives an analytic interval of outgoing angles for which both
outgoing momenta are inside a specified band. Average that interval
with the identical-scalar [dilute collision weight](INDEPENDENT-EARLY-POLE-SCALAR-DILUTE-CLOCK-2026-11-26.md).
We use the same momentum nodes as in the stored spectra and include
neither Bose enhancement nor net gain minus loss.

Across the three x=0.05 third-grid distributions, more than 99.994%
of the dilute event weight involves two incoming particles below
q=10⁻⁵. Conditional on this incoming subset, both outgoing particles
also lie below 10⁻⁵ in more than 99.9999% of the angular weight.
But only about 11–15% retain both outgoing particles in the single
peak occupation cell on that grid. This latter metric changes when
the bin width changes, so it is a geometry illustration, not a
converged redistribution fraction. Against any outgoing cutoff far
above the sum of the two incident kinetic energies, such as q=0.1
for this incoming subset, the absence of hot products follows
directly from two-body energy conservation.

The result does not bound cumulative heating across repeated
scatterings, energy exchange with the pre-existing tail, or the
effect of stimulated final-state factors. It motivates resolving
redistribution within the narrow infrared band rather than
treating all events as equivalent transfers into a hot component.

~~~sh
python -W error early_pole_scalar_elastic_angular_support_reference.py results/early-pole-spectrum-third-m06-2026-11-26.json.gz
~~~

[Raw angular support](results/early-pole-scalar-elastic-angular-support-2026-11-26.json).
