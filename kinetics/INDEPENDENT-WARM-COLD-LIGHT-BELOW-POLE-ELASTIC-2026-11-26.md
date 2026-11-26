# Dark elastic scattering changes the local light rate below the pole

The [fixed-light decay check](INDEPENDENT-WARM-NEAR-REST-LIGHT-INVERSE-DECAY-2026-11-26.md)
can leave out elastic contributions of comparable size. Evaluate
(Lphileftrightarrow Lphi) for a tagged light fermion at rest
on the same six warm snapshots. Use only invariant energies
(sqrt{s}/Mle0.70, 0.82, 0.90) in separate comparisons.
The resonant part is omitted deliberately rather than counted
again alongside decay and inverse decay.

For scalar partner lab energy (E), (s=m^2+m_phi^2+2mE)
at rest. The invariant elastic helicity trace is the same as
the prior scalar-rest calculation at fermion energy (mE/m_phi).
The lab boost has the opposite orientation relative to the
initial light in the center-of-mass frame: the outgoing scalar
energy is (A-Dcos	heta), so cut the angles at the actual
scalar and light energy edges with this orientation. Evaluate
direct and reverse together with exact FD and Bose factors.
Use the exact invariant relative velocity; at rest its
(p_phi/E_phi) cancels the (E_phi/p_phi) in converting
(dp_phi) to (dE_phi). Reverse occupation factors are allowed
to create a scalar final state outside the represented support
when the prescribed invariant cut permits it.

The vectorized scattering trace agrees with the independent
scalar-rest function and its Dirac-matrix trace to better than
(3	imes10^{-15}). Common-temperature equilibrium with separate
fermion and boson chemical potentials cancels event-by-event
within (5	imes10^{-15}), at rest and at two nonzero momenta.
The archived folds retain independent radial and scattering
angle refinement, and compare all three cuts.

At the first matched frame and cut 0.82, the q16 elastic
gain and loss are about (8.710	imes10^{-4}) and
(8.980	imes10^{-4}) in occupation per Hubble time.
Net is (-2.698	imes10^{-5}). The q32 net is
(-1.177	imes10^{-5}). The net stays negative at
rest through all six frames and for all three cuts.
The selected elastic contribution is large enough to matter
relative to the q32 cancellation of spontaneous gain and
inverse-decay loss; the individual elastic gain and loss
are much larger than their difference.

This is a stationary-snapshot comparison for the tagged mode.
It does not include collisions close to or above the pole,
other elastic channels, physical exterior populations,
or an evolved continuum history.

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error warm_cold_light_below_pole_elastic_reference.py
~~~

[Diagnostics](results/warm-cold-light-below-pole-elastic-2026-11-26.json.gz).
