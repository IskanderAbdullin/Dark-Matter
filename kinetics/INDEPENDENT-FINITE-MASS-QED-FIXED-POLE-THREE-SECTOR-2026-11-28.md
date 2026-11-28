# Charged energy has its own neutral occupation at the fixed transverse pole

Use the archived scalar number and energy matrices from the [low union](INDEPENDENT-FINITE-MASS-QED-FIXED-POLE-COLD-TO-SHOULDER-GAP-2026-11-28.md) and the [first future-transfer interval](INDEPENDENT-FINITE-MASS-QED-FIRST-FUTURE-SCALAR-BAND-2026-11-28.md). At one fixed photon frequency ω, if the scalar number and scalar energy currents are J_N and J_E, the photon energy current is −ω J_N and the total charged energy current is ω J_N−J_E. Thus the charged slope is ωD_N−D_E, and its neutral trial follows from the same affine matrix without any new quadrature.

| Archive | Charged-energy-neutral n | Scalar-energy-neutral n | Photon-energy / scalar-number-neutral n |
| --- | ---: | ---: | ---: |
| Fine | 4.3691983e−5 | 4.8430344e−5 | 6.7835248e−5 |
| Third | 4.2091704e−5 | 4.6748933e−5 | 6.5807624e−5 |

In the high band, below the first root the scalar sector loses energy to both other sectors. Between charged and scalar energy neutrality the photons gain while the scalar and charged sectors both lose. Between scalar and photon neutrality the scalar and photon sectors both gain while the charged sector loses. Above the highest root the scalar sector gains at the expense of the other two. At the charged root the photon gain equals the scalar loss (2.466e−7 and 2.417e−7 before the common normalization). Every sign test includes an independent closed energy sum.

The charged slope reverses between the two measured q regions, because their mean scalar energies differ relative to ω. Adding the measured low band to the high band shifts any of the latter roots by at most 1.40e−7 relative. This does not bound the missing interval [.0003,.004] or either tail. All roots are trial values in a frozen affine comparison, not actual photon occupations or an inclusive rate.

Run `OPENBLAS_NUM_THREADS=1 python -W error finite_mass_qed_fixed_pole_three_sector_neutrality.py` in `kinetics/` and compare the [receipt](results/finite-mass-qed-fixed-pole-three-sector-neutrality-2026-11-28.json).
