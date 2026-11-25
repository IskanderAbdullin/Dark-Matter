# Direct cubic branch: initial conditions and observational geometry

## Trigger and extra species

Keeping the background microphysical parameters fixed, varying χ's initial value from 0.1 to 10 times the reference gives first half-activation at z=4.161–4.548, versus 4.3453 at the reference. Near the reference, dz_half/dln(χ_initial)=0.085. The smallest seed gives extra homogeneous energy above 2% in the prior screen (2.045%). The kinetic exposure changes by −9.16% to +10.68% across the same range. These homogeneous trajectories assume an overdamped seed of fixed sign within each patch.

The χ seed in the reference is 9.05e16 eV. A free minimal light spectator would have fluctuation H_inf/(2π) per log scale. A 10% fluctuation relative to the required seed corresponds to H_inf≈5.69e7 GeV per log scale or 7.34e6 GeV in a 60-e-fold uncorrelated random-walk diagnostic. Neither figure is an inferred observational limit. An unbiased distribution can produce different signs in different patches. For a wall network with one wall per Hubble radius, the unmodified double-well tension would give ρ_wall/ρ_tot≈6.76e-4 near half-activation. Neither wall formation nor wall evolution was computed. The η quartic to χ is ≈1.96e-54 and the assumed η radiation population still needs a separate origin. See [stochastic spectator theory](https://arxiv.org/abs/astro-ph/9407016).

For the finite-field quartics, solve the homogeneous A and S spectator equations under radiation domination and match their late oscillatory action. A pure-quartic conformal solution gives an independent normalization, agreeing to ~0.15%. Refining the matching time changes the inferred result by 0.17%. The initial fixed root bracket missed the S solutions, so the final bracket is scaled using the independent conformal solution.

| Quantity at additional Ωh²=.0012 | Mediator S | Light A |
| --- | ---: | ---: |
| Initial field amplitude, eV | 9.58e23 | 4.27e21 |
| H_inf for 60-e-fold Gaussian RMS, GeV | 7.77e14 | 3.46e12 |
| H_inf for pure quartic stationary RMS, GeV | 7.25e13 | 4.80e7 |

The stationary light-scalar estimate is near the trigger's single-log-scale fluctuation diagnostic. It would require extremely long quartic stochastic relaxation, of order λ_A^−1/2 e-folds, and ignores the production portal. Moreover the homogeneous A and S solutions have a long potential self-resonance interval. For the 1% examples the pure-quartic peak Floquet estimate gives ~3.0e7 and 1.5e17 potential amplitude e-folds, respectively, before a crude quartic-to-quadratic crossover. This flags breakdown of the coherent relic extrapolation if modes remain in the unstable band. It does not compute fragmentation, number change or a final abundance. See [spectator dark matter](https://doi.org/10.1103/PhysRevD.98.123532), [the numerical diagnostics](direct_pole_initial_conditions.py), and [raw initial-condition results](results/direct-pole-initial-conditions-2026-11-25.json).

## Geometry with public DESI DR2 distances

Evaluate the original homogeneous direct-pole background with fixed standard densities and its prescribed η and χ evolution. Use all 13 compressed DESI DR2 distances and their covariance, profiling just one free common ruler. Compared with the fixed ΛCDM baseline, the direct-pole fit has χ²=16.30 versus 14.31, Δχ²≈+1.99. Predicted unscaled BAO distances shift by −0.45% to −0.92% and the fitted ruler shifts from 148.64 to 147.68 Mpc. An auxiliary model with η kept relativistic has χ²=14.31 with a free ruler. With an approximate ruler calculated at fixed drag redshift, even the sign of Δχ² changes. The free-ruler and fixed-redshift calculations both hold the original standard cosmology fixed, so this is a conditional consistency diagnostic.

At fixed z*=1090, D_M shifts by −0.472%, r_s by −0.338% and their angle by +0.134%. Recombination was not recomputed. Treating both extra components as smooth, the subhorizon growth amplitude D at z=0 is 1.65% smaller, at matched growth initialization z=100, with −0.23% at z=5. This deliberately smooth approximation leaves out η's force from χ, as well as CMB photon and matter perturbations. These differences motivate testing the coupled η–χ perturbations before an observational fit. Data provenance and calculation: [published DESI DR2](https://doi.org/10.1103/tr6y-kpc6), [pinned covariance and means](https://github.com/CobayaSampler/bao_data/tree/bb0c1c9009dc76d1391300e169e8df38fd1096db/desi_bao_dr2), [script](direct_pole_geometry.py) and [results](results/direct-pole-geometry-2026-11-25.json).
