# Verified CLASS result: one real 86-eV thermal Bose relic

**Source of results:** successfully completed public [CLASS GitHub Actions run 44890003670](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/44890003670), using classy 3.4.3.0 on 21 November 2026. [The source program](thermal_bose_class.py), published [workflow](../.github/workflows/thermal-bose-class.yml), and raw zipped JSON artifact (GitHub retention 21 days) are reproducible. Crucially, the benchmark has not been compared to an observational likelihood.

## Precise assumptions

ONE real scalar thermal boson, mass m=86 eV, zero chemical potential, relativistically decoupled at gstarSdec=55, standard entropy dilution gstarS0=3.909, NO primordial condensate or later collision term. Photon T=2.7255 K; h=.674, omega_b=.02237, baseline massless neutrino contribution N_ur=3.044, n_s=.965, A_s=2.1e-9, tau_reio=.0544. Input PSD f_CLASS(q)=1/[(2pi)^3(exp(q)-1)] with deg_ncdm=1; m specified without manually forcing omega_ncdm. CLASS inferred omega_ncdm*h^2=0.11913012077745594, independently predicted 0.11912939486124709 (relative difference 6.09e-6).

## Actual CLASS background and transfer results

Early-time effective radiation: CLASS Neff=3.10881511794 versus CDM reference Neff=3.044, a difference 0.06481511794. z_eq=3381.43015 for thermal Bose model, z_eq=3382.07751 for cold reference. At z=0 and z=5 we calculated linear P(k), with the same primordial amplitudes and matched present non-baryonic omega (but not an identical early radiation density).

Define transfer-amplitude comparison T(k,z)=sqrt(P_Bose(k,z)/P_CDM(k,z)), where "Bose" here means the particular thermal relic. This is NOT the actual primordial production transfer of any nonthermal alpha scenario.

The first T=0.5 crossing:
- z=0: k=0.43268689336 inverse Mpc.
- z=5: k=0.43275181165 inverse Mpc.
- z=5, T normalized to its largest-scale value: k=0.43187632366 inverse Mpc.
- q-grid convergence: coarse PSD yields z5 T=0.5 crossing k=0.43291270113 inverse Mpc, difference only 0.037178%.

Sample z=5 CLASS results:

| k in 1/Mpc | Amplitude T | Power P_Bose/P_CDM | Coarse/fine relative P difference |
|---:|---:|---:|---:|
| 0.051026 | 0.997736 | 0.99548 | 0.13% |
| 0.097530 | 0.980891 | 0.96215 | 0.15% |
| 0.204496 | 0.881486 | 0.77702 | 0.074% |
| 0.296113 | 0.742952 | 0.55198 | 0.080% |
| 0.390874 | 0.574983 | 0.33060 | 0.107% |
| 0.428775 | 0.507183 | 0.25723 | 0.115% |
| 0.515960 | 0.361297 | 0.13054 | 0.135% |
| 0.620872 | 0.214849 | 0.04616 | 0.073% |
| 0.681075 | 0.168939 | 0.02854 | 0.498% |

At z=5, compare only modes with T>=0.05 (65 sampled points, k<=0.7471 inverse Mpc) for a robust relative P convergence statement: max coarse/fine error 1.70756%, and half-amplitude crossing stable to 0.0372%. At much higher k, the transfer-amplitude is tiny and P may contain weak oscillations/cancellations: the maximum relative coarse/fine disagreement over *all* 95 points reaches ~710.3 at k~4.34 inverse Mpc. Those deeply suppressed numerical tails are **not converged** and no exact high-k constraints should be extracted from them. The successful CI explicitly enforces convergence only on the non-negligible range and the half-amplitude location; it preserves the unresolved tail in the artifact for transparency.

## Physics interpretation — a conditional but severe structure-formation tension

The standard thermal 86-eV boson makes the linear matter power **only one quarter** of the matched-CDM power around k=0.433 inverse Mpc, and about 4.6% at k=0.621 inverse Mpc. These are very large comoving scales (2pi/k about 14.5 and 10.1 Mpc), far above an individual dwarf halo, indicating strong suppression of the seeds of structure in this *specific* thermally produced all-DM model.

For context, a 2024 hydrodynamical Lyman-alpha analysis reports a 95% lower bound of ~5.7 keV for its **standard thermal relic WDM fermion parameterization** (https://doi.org/10.1103/PhysRevD.109.043511). This bound is NOT a direct 86-eV boson mass limit and must not be transplanted blindly: species counting, PSD, phase-space temperature, effective primordial radiation and nonlinear IGM physics differ. The qualitative strong suppression found by CLASS motivates an explicit data confrontation, but no CMB/Lyman-alpha likelihood is evaluated in this project.

**Next falsifying target:** a physically produced *colder/nonthermal* Bose f(q) with fixed 86-eV mass and omega, with a concrete collision/production Lagrangian and validated CLASS transfer. A simple PSD-shape scan can illustrate dependence on momentum coldness, but cannot substitute for that physical production derivation.

## Original code references

[CLASS explanatory.ini](https://github.com/lesgourg/class_public/blob/master/explanatory.ini) for N_ncdm PSD file and T_ncdm; [source/background.c](https://github.com/lesgourg/class_public/blob/master/source/background.c) for default f0(q) containing the quantum-cell factor, and [python/classy.pyx](https://github.com/lesgourg/class_public/blob/master/python/classy.pyx) for the power-spectrum and background accessors.
