# Timelike soft poles and the finite angular-resolution time

Complement the [spacelike calculation](INDEPENDENT-EARLY-POLE-FINITE-WILSON-SPACELIKE-2026-11-26.md)
with the timelike poles of the leading soft equilibrium photon response.
This is still an external Abelian closed-current test, not a computed
neutral production vertex. Conventions for the soft thermal QED response
are compared with [Gorda et al.](https://doi.org/10.1103/PhysRevD.107.036012);
no next-to-leading terms are taken from that calculation.

The pole weights in the frequency-to-momentum ratio x are
2/(x |derivative of the real inverse propagator|) in the same convention
as the spacelike integral. Both branches tend to omega/mD = 1/sqrt(3)
as k vanishes. The longitudinal residue above k/mD=4 is exponentially
small and is omitted. The matched free thermal transverse comparator is
evaluated with the *same* closed current and cutoff. It is not a
prescription for subtracting part of an interacting production rate.

**Independent low-k calibration.** At fixed nonzero frequency and
vanishing k, the squared spatial current divided by (delta-v t)^2
tends to F(z)=1−2 sin(z)/z+2[1−cos(z)]/z², z=omega t. In the chosen
closed contour this is dominated at large z by the endpoint connector.
The finite-k test projects the Coulomb longitudinal term using the
Ward identity, so J0/k is retained even when J0 tends to zero. At the
physical comparison time, projected longitudinal and total transverse
geometric powers divided by t² approach 1.21437e−7 and 2.42873e−7
respectively. In contrast, taking omega/k fixed inside the cut makes
the squared current vanish like k². Treating these limits as
interchangeable would discard a finite timelike mode weight.

**Matched velocities and time.** The sample uses the parent pole-energy
midpoint at x=.05, q/T=4.1e−6. At g=.005, one inverse angular half-span
means t mD≈996 with spatial separation .6010/mD. Restoring the small
transverse daughter velocity moves the two spacelike cut terms by
less than 4e−5 *relative* in the tested coarse/fine azimuth scan. The
timelike results below use the collinear projection.

| Time in inverse angular half-spans | Closed-contour amplitude, cutoff 20 mD | Cutoff 40 mD |
| ---: | ---: | ---: |
| .25 | .99860 | .99829 |
| .5 | .99108 | .99040 |
| 1 | .94812 | .94673 |
| 2 | .75486 | .75257 |
| 4 | .29493 | .29312 |
| 8 | .01660 | .01640 |

These are exp(−χ) for the included classical equilibrium soft sectors;
the normalization matches the static leading HTL coefficient. At one
inverse half-span and 20 mD, spacelike longitudinal/transverse
contributions to χ are .03771/.01428; the two plasmon contributions sum to
.00128. The plasmon part is not negligible for every contour or time:
at one quarter of that time it supplies about 17% of the represented
exponent with the same cutoff. At g=.001, even the one-time total
exponent increases from .000154 to .000200 when the cutoff moves from
20 to 40 mD; a small absolute value is not evidence of cutoff
independence. The example hard cutoff, equilibrium screening,
straight trajectories, finite endpoint prescription and classical limit
remain assumptions. No production factor or CLASS handoff follows.

Reproduce from kinetics/ with OPENBLAS_NUM_THREADS=1:
early_pole_finite_wilson_poles_reference.py,
early_pole_matched_wilson_geometry_reference.py,
early_pole_noncollinear_wilson_cut_reference.py, and
early_pole_wilson_formation_time_reference.py. Use --level 0/1 or
--nodes 5/9 as appropriate. The compressed result bundle records both
resolutions, cutoffs and the separately executed small-k calibration.
