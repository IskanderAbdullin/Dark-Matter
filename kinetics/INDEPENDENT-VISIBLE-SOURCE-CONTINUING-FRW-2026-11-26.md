# Continuing visible production in the two-moment FRW reference

This independent finite-support check adds ongoing input from the visible FD bath to the unscaled joint scalar/stable-fermion evolution. It is not a coupled visible/dark kinetic solution.

## Source construction and checks

The source separates direct light-pair production, light fermions from heavy-pair production and prompt decay, and scalars from the same decay. The direct lab source integrates the monopole and quadrupole angular weights analytically over each intersected shell. For each lab-frame heavy-particle energy node the scalar and stable-fermion daughters have uniform distributions over different lab energy intervals. Each overlaps the physical shell edges at p=q/a. Compute both positive particle counts and kinetic energies, maintaining the low and high tails separately. The sum of both daughter lab energies equals the parent lab energy.

At 45 MeV the direct number and energy and the prompt total energy agree with the reduced event moments within 6e-16 relative. The scalar number histogram agrees with the earlier angle-overlap source within 3.9e-15 in normalized L1. In a separate source-only integration from empty at 50 MeV to 45 MeV, summing the kept and omitted counts with a³ d(log a)/H reproduces the independent collisionless source yield. The maximum relative channel defect decreases from 1.46e-4 to 3.65e-5 under 6-to-12 midpoint refinement. Approximating each shell and tail by one energy when redshifting gives a maximum final channel energy-yield difference of 3.53e-5 at 12 steps. It is a separate cross-check, not an exact within-shell extrapolation.

## Warm-support unscaled evolution

The actual warm run starts at the same 45 MeV input snapshot as the earlier injection-disabled calculation, with fixed q/T_initial in approximately [0.1607, 15.8747]. At each midpoint apply the positive visible injection in two equal number/kinetic-energy increments bracketing the implicit collision solve. The bath EOS is imposed, and H is frozen at the pre-injection midpoint, including the current dark energy. Evolve to a=1.072508 relative to the 45 MeV input. Report each species' a³ number and a⁴ energy, the excluded low/high tails, the redshift trace contribution, and the collisional energy exchange separately. The source and collision results below use the actual source and collision amplitudes, not manufactured amplitude suppression.

| Substeps | Scalar energy from collisions, divided by initial total dark a⁴ energy |
| ---: | ---: |
| 3 | 0.1082123 |
| 6 | 0.1107116 |
| 12 | 0.1119231 |
| 24 | 0.1125195 |

At 24 steps the source adds 0.09032 (scalar) and 0.32521 (stable fermion) in the same energy units, and 43.99% and 45.58% of the respective initial counts. Number/energy ledgers close to 3.3e-13, including redshift work and the external source. Heat transfer has residual first-order time-step sensitivity; the 12-to-24 change is about 0.53% relative. Removing elastic collisions at 6 and 12 steps has no resolved effect on this warm-support bulk-energy statistic. It does not imply the colder states decouple or that the spectrum is insensitive.

The source keeps visible FD factors but uses vacuum dark final-state factors and omits inverse visible reactions. It assumes an immediate isotropic decay of each heavy particle with the specified benchmark kinematics; it does not evolve a heavy-state distribution. There is no feedback to visible entropy or temperature. The principal run omits both momentum tails, including the coldest scalar bins. Its seeded 45 MeV snapshot is a collisionless production calculation and is not claimed to be the correct initial state in the presence of collisions. The independent from-empty test verifies source integration only, not an interacting run from empty. This calculation establishes none of a relic yield, a condensate, or universal kinetic equilibrium. A next control is production from an actually empty initial state with positive first-source initialization, an expanding low-momentum support, and a dark final-state/portal detailed-balance check.

## Reproduce

Use the external author central EOS (see the EOS citation in the preceding benchmarks) and the verified empty-pair spectrum. Do not commit the external table.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1
python -W error empty_lepton_pair_spectrum_reference.py EOS.dat > empty.json
python -W error joint_visible_fd_source_reference.py EOS.dat --independent-spectra empty.json
python -W error joint_resonant_frw_split_reference.py EOS.dat empty.json --physical-amplitude --warm-source-init --visible-source --source-only
python -W error joint_resonant_frw_split_reference.py EOS.dat empty.json --physical-amplitude --warm-source-init --visible-source --steps 6 12
python -W error joint_resonant_frw_split_reference.py EOS.dat empty.json --physical-amplitude --warm-source-init --visible-source --with-elastic --steps 3 6 12 24
```
