# Degeneracy, daughter removal, and Yukawa-motivated rates

## Validated conditional candidate

The coarse width tuning at n=120, kappa2=.03 and x0=.1, with constant parent chemical relaxation/H=100:
- m2/M=1/sqrt(2): kappa=.0306074372, omega_phi h²=.119999984, p/T>.1 fraction=.062279.
- m2/M=.9: kappa=.00229534059, omega=.119999957, hot fraction=.00362262.

At the latter fixed width, n=240 and 320 (different IR cutoffs) give omega=.12110849/.12122412, hot fractions .00360959/.00360798, and mean p/T=.00240436/.00240276. The scalar width, not its spectrum normalization, was changed.

[Yukawa-motivated rates](yukawa_rate_bath.py) replace constant parent relaxation with time-dilated vacuum decay/inverse decay plus a relativistic scattering estimate gamma*y_i²*T, gamma=.01, tapered as 1/(1+(M_i/T)^2). The taper and equilibrium momentum/number relaxation are phenomenological. The same rates for chi2 include inverse creation and energy return. See Yin, JHEP05(2023)180, sec. IV.4, https://doi.org/10.1007/JHEP05(2023)180 for the RHN/ALP motivation. This is not a calculation of thermal self-energies or the full collision integrals of that theory.

For a parent of 86 TeV, m2/M=.9, kappa=.00224, kappa_SM1=4000, kappa2=.03, x0=.1 and xend=80:

| mesh | omega_phi h² | p/T>.1 fraction | mean p/T |
|---:|---:|---:|---:|
| 240 | .120112652 | .003516464 | .002336004 |
| 320 | .120223852 | .003515138 | .002334554 |

Initial effective parent rate ~99.72 H; chi2 at q=3 ~.000827 H. Number-source and continuity guards passed (final 320 continuity residual 2.07e-8). At xend the fermion energy is numerically negligible, and the bath remains above 1 TeV, consistent with the constant high-temperature SM g approximation. The parent is initialized in thermal equilibrium; whether it actually reaches it is a separate question. The model uses g1=g2=4 and a Dirac spin-averaged channel; it is not an identified Majorana-neutrino model. The illustrative coupling mapping has a nominal single-species seesaw mass contribution ~4.28 eV for kappa_SM1=4000; in a Majorana seesaw interpretation it would require cancellations or a different mechanism, and cannot be presented as a realistic fit.

## Spectral and numerical checks

The previous two-tracer export on the 240-cell constant-parent candidate failed its hot-tracer number quadrature by -6.56%; the 120-cell log-PCHIP reconstruction also failed the 1.5% moment guard (max bias 2.79%). The positive [three-tracer representation](../cosmology/three_scale_psd.py) partitions one PSD into cold, intermediate and hot coordinates. For the 320-cell Yukawa-motivated candidate, number/energy/q² reconstruction biases are -0.2171%, -0.3218%, -0.4957%; all transformed 120-point quadrature moment deviations <9.48e-5 and 240-point alternate-partition deviations <8.51e-6. No artificial rescaling.

[Free-streaming and CLASS interface](../cosmology/daughter_bath_class.py): exact shell-weighted homogeneous trajectories in the background including the scalar energy moment. To the equal-density CDM reference equality, mean path .007528 Mpc, median .001961 Mpc, RMS .102013 Mpc, fraction above .06 Mpc .003849. The tiny high-momentum component dominates RMS; a cold median is not itself a transfer constraint. The numerical integration is independently refined. The massive daughter has returned energy to the SM bath before neutrino decoupling, so no surviving chi2 dark radiation is added in this conditional model. The raw PSD already includes the entropy change.

The [public workflow](../.github/workflows/daughter-bath-class.yml) regenerates both kinetic meshes, runs synchronous CLASS at k=.003–20 Mpc^-1 and a full-hierarchy control up to 1.5 Mpc^-1. Its results are not claimed until logs are verified. Actual local spectra and summaries are in results/yukawa-degeneracy-2026-11-25.json. Next test empty-parent initial conditions and weaker parent coupling, then evaluate physical thermal dispersion and late self-interactions.
