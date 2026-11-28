# Equilibrium screening is too large on the frozen early snapshot

This is a frozen, isotropic, collisionless Abelian charge-response
check. Electromagnetic charges of each Dirac species are equal in
magnitude, with symmetric particles and antiparticles. Photons themselves
carry no Abelian charge. This is not a new cosmological trajectory.

For four spin and charge states, the static electric susceptibility
divided by g² is (2/π²) times the integral of p E (−df/dp).
Integrating by parts gives a positive, exact cell expression for the
archived piecewise-constant occupations: sum f times the difference of
p E across each cell. On a common massless FD bath it gives T²/3 per
Dirac species. For the spacelike Landau response per phase velocity
s=omega/k, the radial derivative weight is E² (−df/dp) for v>s.
The classical statistical factor for independent fermion occupations
is instead p E f(1−f) in the same domain. Their ratio recovers the
temperature for a single FD state at every s, but not necessarily for
a collection of nonthermal shapes. Here s=0 denotes the limiting cut
slope, not a division by a zero-frequency spectral density.

| Snapshot / fit | Debye mass squared vs two thermal Dirac species at T_ref | Ratio of noise to response at s→0, in T_ref |
| --- | ---: | ---: |
| Resolved x=.05, third grid | .49951 | 1.01347 |
| Resolved x=.05, fine grid | .49869 | 1.01268 |
| Separate within-species initial FD fits | .52702 | .95757 |
| Conditional g=.005 fast shared-heat endpoint | .67314 | .97551 |

The early resolved daughter contributes only about 9.1% of the
static susceptibility; the conditional endpoint with hard pair
production gives approximately equal contributions. The early raw
ratio remains close to 1.01 at s=.99, then rises to 1.065–1.068 at
s=.999. This last number is a comparison across represented cells,
not a calibrated prediction at the light cone; raw directional
threshold crossings are sensitive to the derivative reconstruction.

Replacing the equilibrium mass in the earlier path cutoff ratios
without also revisiting the noise and dispersion would be inconsistent.
For instance x=.05 would have m_D about .707 of the two-thermal
reference, although the low-phase-velocity effective ratio stays near
T_ref. No single temperature fits both inputs. We have not recomputed
the full dynamical spacelike and timelike propagators with these frozen
distributions, nor updated the kinetic trajectory or its scalar spectrum.

Reproduce from kinetics/ with early_pole_snapshot_hard_loop_screening_reference.py.
Raw diagnostic: results/early-pole-snapshot-hard-loop-screening-2026-11-26.json.gz.
