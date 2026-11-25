# Velocity resonance and direct late activation: conditional consistency results

2026-11-25. Continuation from c32962f. Inputs: original independent parent-angle 240 and 320 bin scalar spectra at matched kappa=.00685. This is a diagnostic investigation, not a fit to data or proof of exclusion.

## Velocity-only interaction

For an isotropic frozen comoving PSD the pair relative-speed density follows from uniform relative cosine. We integrate the cosmological dilute exposure for the cold subset p<6 p_peak (about 80% of particles), restricting every selected particle to speeds below .03c. For any nonnegative cross section supported in the sampled interval, the integrated ratio of cosmological to halo collision counts is bounded below by the minimum of the ratio of the exposure kernels. The halo comparison assumes a constant 10 Gyr Maxwell distribution.

Independent grids 240/320 give minima 90.47/89.71 for rho=1e-24 g/cm3, v_T=100 km/s; 8.408/8.339 for illustrative critical F583-4 rho=6.70e-24, v_T=62.37 km/s; 2.69/2.67 for rho=1e-23,v_T=30; .087/.086 for rho=1e-22,v_T=10. The last case illustrates why this is not a universal no-go. The relative-speed normalization agrees with the squared cold number fraction at several ppm. These are positive dilute collision counts, NOT net kinetic relaxation.

For an infinitesimally narrow velocity resonance normalized to one dilute halo collision in 10 Gyr, direct NR angular integration at v_res=90 km/s gives cold all-four exposure 8.385/8.313 and large positive Bose entropy production on the frozen PSD. The isolated cold-subdomain entropy budget is exceeded by ~4.9e8 on n320; the entropy estimates differ by a factor of about two across the two quadratures, so do not quote them as converged. Even the dense slow-halo example (14 km/s) has nominal entropy production >2.6e6 times its cold entropy budget. These are frozen-trajectory inconsistency screens, NOT actual entropy changes or an exclusion. They assume instantaneous vacuum elastic scattering, which the next calculation shows is not appropriate here.

## Explicit intermediate state

For an identical-scalar single elastic s-wave Breit-Wigner resonance, matching the halo area implies vacuum lifetimes 120.3,52.2,11.9 Gyr for 40,90,140 km/s respectively, widths 1.73e-34,4.00e-34,1.75e-33 eV. For uniform production in a 10 Gyr halo, only about 4%,9%,32% decay during that interval. A local instantaneous cross section is not a valid description. Fixed-A inverse-decay evolution with the exact NR angular average predicts 1% conversion already at z of order several million. Its later unphysical fractions above one are a backreaction warning, not abundances.

A conservative expanding AA<->R operator uses positive geometric interpolation, preserves A+2R and represented kinetic plus excess energy, includes explicit expansion work, and passes detailed balance, analytic Jacobian and nonnegative entropy checks. Start z=2e6, qmax=20 peak, retains 97.6-97.7% of input number; omit hotter modes and medium self-energy. At 90 km/s BDF nA100 and independent nA140 Radau/320 input give maximum R-charge fractions .5118/.5062 and final .3909/.3943 at z~3.15-3.19e4. Charge residuals below 1e-15, energy residuals 1.7e-4/1.1e-9. Direct initial inverse-rate discrepancy on the finer run is -1.25%. The 40 and 140 km/s runs give final .3157/.3286. Some tighter BDF runs fail near strong IR pile-up, including n180; preserve their outputs. These are not predictions at z=0, and the tiny vacuum detuning and occupied medium require a dressed propagator and possibly coherent dynamics. No general statement about all velocity-dependent models is warranted.

## Direct quartic activation

V=m^2 A^2/2 + chi^2 A^4/(24 Lambda^2) + lambda_chi (chi^2-v^2)^2/4. At sigma/m=.1 cm2/g, z=3, lambda_late=1.082e-8. The illustrative MS-scale scalar double-bubble contribution lambda_late (m^2/(16 pi^2))^2/8 is 2.967e-6 eV^4, 3.87e5 times a 1% total-energy budget. A Gaussian density trigger is even more radiatively sensitive, and the induced lower-dimension operator cannot simply be discarded. This is a naturalness diagnostic, not an invariant finite prediction: counterterms and symmetries matter. A simple thermal bath with Delta N_eff=.2 also has a very small allowable loop quartic energy. The finite spherical screening test admits a narrow high-Lambda window (~1.83e17 to 6.4e17 eV for screening radius parameter <=1 and mu/H>=10 at z=3), so screening alone does not prove impossibility.

## Sources and priorities

- Chu, Garcia-Cely, Murayama, resonant SIDM: https://doi.org/10.1103/PhysRevLett.122.071103
- Symmetron: https://doi.org/10.1103/PhysRevD.84.103521
- Thermal dark energy: https://doi.org/10.1103/PhysRevD.101.023503 and https://doi.org/10.1103/PhysRevD.110.023533
- Two-loop effective potential: https://arxiv.org/abs/hep-ph/0111209
- ACT DR6 cosmological parameters: https://doi.org/10.1088/1475-7516/2025/11/063 (illustrative Delta N_eff choices here, not a likelihood)

Prefer a genuinely late interaction with a smaller radiative threshold, e.g. softly mixed light states, and a shift-symmetric production portal. The original constant-coupling nonlinear kinetic/IR problems remain open. Reproducible raw results in results/velocity-resonance-direct-2026-11-25.json; scripts in this directory; tests/test_resonance_activation.py.
