# Collisionless streaming diagnostics for the represented dark spectra

Take the actual represented shell occupations from the kinetic
snapshots, not moment-matched thermal shapes. Normalize the absolute
scale factor using the photon temperature and the visible entropy table.
After the selected snapshot, let these represented particles simply
redshift in a reference flat Planck-like background
(H₀=67.36 km/s/Mpc, Ωₘ=.3153; see [Planck 2018](https://arxiv.org/abs/1807.06209)).
Use the full EOS radiation contribution at MeV temperatures,
constant post-annihilation radiation factors below the last table
knot, and a massless-radiation approximation for the reference background.

For each lab momentum p at snapshot scale factor aᵢ, the comoving
momentum q=aᵢp is held fixed. Integrate v(q,a)/(a²H(a)) from aᵢ to
reference matter–radiation equality, with
v/c=q/√(q²+m²a²); then take number-weighted averages and quantiles
on the represented distribution. This is a **comoving path length**,
not a fitted power-suppression scale, halo prediction or Lyman-α
constraint. See, e.g., [Free-streaming integral for a nonthermal distribution](https://arxiv.org/abs/2112.14668)
and [Lesgourgues and Tram, distributions in CLASS](https://arxiv.org/abs/1104.2935)
for the separate collisionless perturbation computation.

| Input archive | Mean scalar distance | Mean stable fermion distance |
| --- | ---: | ---: |
| 9.019 MeV reheat, 3.157 MeV vacuum source only | 20.04 Mpc | 0.2403 kpc |
| Same, early scalar-pair conversion and vacuum visible source | 19.00 Mpc | 0.2403 kpc |
| Same, signed visible, prompt, elastic and pair evolution | 18.33 Mpc | 0.2403 kpc |
| 50 MeV reheat, 12 MeV archived source | 10.47 Mpc | 0.1992 kpc |
| Same high reheat, 9.3 MeV archived source | 10.20 Mpc | 0.2006 kpc |

For the low-reheat signed case, the fermion 10–90% range is about
0.106–0.398 kpc and its equality RMS speed is about 1.30e-7 c.
For scale only, a hypothetical 5.7 keV, two-spin thermal FD relic
adjusted to carry Ωh²=.120 in this same reference cosmology travels
about 55.6 kpc on average. That thermal comparator is not a
confidence limit derived in this work. The low-reheat stable component
is kinematically far colder than that comparator; it is not assumed
to be thermal. The scalar component is far more mobile and, in the
[early-plus-late paired abundance comparison](INDEPENDENT-SCALAR-PAIR-FROM-EMPTY-2026-11-26.md),
its total fixed-count mass fraction is tiny. The low input is
finite-support and does not resolve its coldest scalar modes.

As a separate homogeneous density diagnostic, if the same retained
low-reheat signed archive coasted collisionlessly to photon
temperature 1 MeV, its stable energy would be about 5.85e-7 of the
visible radiation, with the scalar contribution below 1e-10.
This is not a full BBN or neutrino-decoupling calculation.
The high-reheat snapshots are placed in the reference cosmology
only to compare their kinematics; their stable counts do not in
general provide the reference cosmological background.

Momentum and expansion quadrature are refined independently.
Do not infer a complete transfer function: the early creation
epoch, subsequent interaction history, perturbations and arbitrary
exterior modes are absent. Distances start at the specified archive.

```sh
OPENBLAS_NUM_THREADS=1 python -W error archived_collisionless_streaming_reference.py EOS.late.dat results/{low-reheat-matched-vacuum-shape,low-reheat-early-pair-vacuum-shape,low-reheat-early-pair-quantum-shape,12-quantum-48-shape,early-9p3-quantum-shape}-2026-11-26.json
```
