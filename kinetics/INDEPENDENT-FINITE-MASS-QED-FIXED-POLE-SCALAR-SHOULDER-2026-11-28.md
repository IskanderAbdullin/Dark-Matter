# Cold and shoulder scalar redistribution at one frozen transverse pole

Use the same k=.0015 and the same frozen occupations as in the [cold-band test](INDEPENDENT-FINITE-MASS-QED-FIXED-POLE-COLD-SCALAR-BAND-2026-11-28.md). Compare the cold [1e−6,1.5e−5] and shoulder [4e−5,3e−4] separately. They deliberately omit the gap [1.5e−5,4e−5] and all other momenta: the two-band union is not the scalar spectrum.

Each band has scalar number source nD−B, with n the same trial photon mode population, and scalar energy source nD_E−B_E. Each shoulder grid has two independently refined complete arm/cosine/q meshes. The number source is positive in the cold band and negative in the shoulder when n=1e−8, on both archives. The same n cannot neutralize both bands.

| Archive | Cold n (number) | Shoulder n (number) | Shoulder n (scalar energy) | Pair n (number) | Pair n (scalar energy) |
| --- | ---: | ---: | ---: | ---: | ---: |
| fine | 2.927279e-09 | 5.655159e-08 | 5.6648759e-08 | 1.5791237e-08 | 1.6233995e-08 |
| third | 2.3756479e-09 | 5.6792851e-08 | 5.6898145e-08 | 4.0805222e-08 | 4.1382232e-08 |

The two quadrature meshes agree within 4.77e−6 for the integrated shoulder coefficients on either fixed archive. This is not convergence between archives: the pair number-neutral occupation differs by a factor of about 2.58. The shoulder gross forward number coefficients differ by just .69%, but the corresponding cold coefficients differ by a factor of about 7.6. The sign contrast is reproduced; an archive-independent two-band magnitude is not established.

The two-band number-neutral value n_* balances the sum of the two number contributions X_i=n_*D_i−B_i. Write e_i^D=D_{E,i}/D_i and e_i^B=B_{E,i}/B_i. Its scalar-energy balance is
```text
energy(n_*) = X_c(e_c^D−e_s^D) + B_c(e_c^D−e_c^B) + B_s(e_s^D−e_s^B).
```
The first term measures redistribution between intervals with different mean energies; the last two are within-band asymmetries between forward and inverse terms. The decomposition is checked against direct full-energy differencing in extended precision.

In the cold interval at its own number-neutral n, both archives have a sequence of positive integrated panels followed by negative ones. The reducer compares each panel against its own last-mesh drift and floating-point uncertainty; these are panel signs, not pointwise sign certification. The absolute positive and negative partial integrals pair off even though the cold scalar energy is not neutral. All sources below omit the common factor g²y² k² dk/(64π⁵).

| Archive | First negative cold panel starts at q | Cold positive extent | Pair scalar energy at pair number-neutral n | Between-band part of pair scalar energy |
| --- | ---: | ---: | ---: | ---: |
| fine | 3.988178e-06 | 9.094593e-18 | -2.420625e-19 | -2.289037e-19 |
| third | 3.642568e-06 | 8.023284e-19 | -1.096422e-19 | -9.532814e-20 |

Run `OPENBLAS_NUM_THREADS=1 python -W error compare_fixed_pole_shoulder_pilot.py` from `kinetics/` to replay the [receipt](results/finite-mass-qed-fixed-pole-scalar-shoulder-2026-11-28.json) from the four archived shoulder meshes and existing cold grids. The raw shoulder panels can be regenerated with `python -W error fixed_pole_scalar_q_band_pilot.py GRID MESH shoulder OUTPUT`, for fine/third and meshes 0/1. No actual photon population, continuum collision term or inclusive rate is inferred.
