# Heavy-parent recycling access on the represented archived spectra

The [upstream energy-moment upper bound](https://github.com/IskanderAbdullin/Dark-Matter/blob/90c141722376dab64e4cd10422d545474817787c/kinetics/PHI-BOSE-PORTAL-RECYCLING-ENERGY-BOUND-2026-11-24.md)
includes repeated production of χ₁/anti-χ₁ followed by decay,
even when one scalar count might eventually be converted to more.
Here evaluate the forward endothermic rates at fixed nonthermal
archives. Ignore final heavy Pauli blocking to obtain an optimistic
gross rate. If every produced heavy particle decays to a stable
fermion plus a scalar, each such parent adds one net dark particle.
Do not include later destruction of the parents or their decay
products when computing this instantaneous gain.

For φφ→χ₁ anti-χ₁ use the signed scalar-pair benchmark amplitude in
reverse with the actual heavier produced mass and lighter exchange
mass, including spin and identical-scalar normalization. For
opposite-charge χ₂ pairs, integrate the t-channel scalar transfer
rational function exactly over the two-body angle; independent
log-transfer quadrature verifies the integral to 1.6e-13. For the
two same-charge channels, Cauchy applied to the t-u terms and the
identical final-state factor bounds each cross section by twice
the corresponding opposite-charge result. Both same-charge initial
identical-pair factors are then included.

| Represented archive | φφ contribution to count growth / (N H) | Opposite-charge stable contribution | Upper incl. both same-charge channels |
| --- | ---: | ---: | ---: |
| 50 MeV reheat, 12 MeV archive | 3.645e-6 | 5.951e-5 | 1.822e-4 |
| 50 MeV reheat, 9.3 MeV archive | 8.159e-9 | 3.670e-7 | 1.109e-6 |
| 12 MeV populated q32 control | 4.137e-6 | 6.404e-5 | 1.963e-4 |
| Same q32 inner state cropped at q16 | 3.783e-6 | 6.045e-5 | 1.851e-4 |
| Low-reheat vacuum 3.157 MeV represented state | 0 | 0 | 0 |

Zero at low reheat is geometric: the maximum incoming invariant
on that retained support is below 2m₁. It does not constrain an
unrepresented hard tail. The scalar term changes more under
cropping than the stable term, even though both are small in
these snapshots. With pure momentum redshift and no redistribution
the 12 MeV q16 fermion channel closes at a/a_archive=1.783; the
populated q32 support allows it geometrically through 3.565.
These are thresholds of finite supports, not bounds for an
arbitrary interacting state. The rates are not integrated
cosmological multiplication ceilings. Earlier or harder states,
other number-changing channels, and changes to the source remain
outside this fixed-state test.

```sh
OPENBLAS_NUM_THREADS=1 python -W error archived_heavy_recycling_gate_reference.py EOS.late.dat results/{12-quantum-48-shape,early-9p3-quantum-shape,low-reheat-matched-vacuum-shape,wide32-embedded-vacuum-shape}-2026-11-26.json --crop-last 16
```
