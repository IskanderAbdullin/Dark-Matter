# Parent support changes the sign of the cold signed-recoil control

The [fixed-parent diagnostic](INDEPENDENT-FINITE-MASS-QED-FIXED-Q-SIGNED-SOFT-2026-11-27.md)
showed negative residuals after holding three parent points
inside the Born window. For the signed real operator, the
allowed parent interval changes when a photon is emitted or
absorbed. It is not valid to integrate just those frozen
Born parents and subtract the virtual term there.

For each plasma photon energy and direction, solve both
collinear endpoints of the exact parent interval. Above
P=0.2 the endpoint equations are monotone for the tested
masses and |omega|≤1e−5; below that momentum the mass/energy
gap is much larger than the scalar-direction term. Split the
parent-energy integral at the native occupational edges of
both charged species, using the energy-translated daughter
edges for each sign. Also separate the parent region outside
the common Born window. Split photon angles where any such
edge enters or leaves the allowed interval, resolving the
daughter forward region in log angle. Both signed real
source integrals use the corresponding Pauli factors.
The virtual subtraction here is only its universal leading
soft part, normalized by the independently known integrated
Born aperture and the leading radiator.

At photon energy 1e−5, the inside-Born even difference and
the new outside-Born contribution nearly cancel. For example,
the latter is about +1.43e−6 on fine, while the signed residual
is about −1.86e−9. Dividing the remaining even part by
omega² gives −18.63 for fine and +5.14 for third at that
photon energy. Near 1e−7 the finite normalized even part is
positive on both grids (about 1.006 and 2.216). This directly
limits interpretation of the fixed-parent sign test.
Forward and backward modes and the stimulated factors are kept
separate until the subtraction. Substantial cancellation is
resolved with extended precision and checked under separate
photon-angle and parent-energy refinement.

An independent calculation integrates over parent energy first,
using the fixed-parent kernel and geometric panels resolving
the narrow daughter-collinear opening. The separate signed
stimulus integrals agree within 2.1e−17 absolute at signed
energies ±3e−6 and ±1e−5 (both archive grids). A single unsplit
parent panel did not resolve that opening.

Integrating the photon energy from zero to the plasma cap
5e−6 gives −0.02830 (fine) and +0.01389 (third) relative to each
integrated Born stimulus, **before** the common gauge factor.
The sign difference persists under caps 2.5e−6 and 1e−5.
A finite near-zero continuation and refined angular, parent and
photon panels change the integrated unnormalized components by
at most 2.9e−12 absolute. This is a sign difference between the
represented archive inputs, not between two numerical evaluations
of the same input or a sign of the full radiative correction.

This cutoff-dependent partial contribution does not include
finite virtual terms, an in-medium soft spectral function,
or matching to hard photons. It is not an inference about
the sign or convergence of the physical scalar rate.

Replay from kinetics with the archived spectra:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_parent_integrated_signed_soft_reference.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_parent_signed_soft_reverse_check.py
```
[Raw receipt](results/finite-mass-qed-parent-integrated-signed-soft-2026-11-27.json).
