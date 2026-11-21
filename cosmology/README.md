# First full Einstein–Boltzmann transfer benchmark for an 86-eV Bose relic

**Model being tested:** ONE real thermal relativistic Bose scalar, decoupled at effective entropy degrees of freedom gSdec=55 with zero chemical potential and **no** primordial condensate, thereafter collisionlessly redshifted. This is not the only production history allowed by the original Abdullin–Popov JCAP 2021 galaxy phenomenology.

## Independent conservation and abundance

With T_CMB=2.7255 K, modern photon+neutrino gS0=3.909 and assumed gSdec=55:
- T_B/T_gamma=(3.909/55)^(1/3)=0.41422311185;
- relic scalar density today 14.59573862 particles per cubic centimetre;
- omega_B h² = 0.11912939486;
- early Delta Neff = 0.06481510, before the scalar becomes nonrelativistic.

These are conditional analytic calculations, not observed constraints. They differ slightly from rough earlier benchmarks using alternative entropy/critical-density conventions; the precise constants are documented in the code.

## Important normalization detail from the actual CLASS source

The official CLASS explanatory.ini permits a custom noncold distribution through use_ncdm_psd_files=1; source/background.c shows its built-in thermal Fermi distribution carries the (2pi)^(-3) quantum-cell factor. To describe ONE real bosonic state, not a neutrino+antineutrino pair or an overnormalised occupation, the input table contains

f_CLASS(q)= 1 / [ (2*pi)^3 * (exp(q)-1) ]

with q=p/T_ncdm, deg_ncdm=1, m_ncdm=86 eV, T_ncdm=0.41422311185 in units of photon temperature. The independent moments are 4pi int q²f dq = zeta(3)/pi² and 4pi int q³f dq = pi²/30. The locally tested finite-grid numerical moment errors are about 1.03e-5 and 7.8e-7, respectively. The code **does not** also prescribe omega_ncdm; passing both m and omega tells CLASS to renormalize the PSD, destroying the independent relic-abundance test.

Source: https://github.com/lesgourg/class_public/blob/master/explanatory.ini ; https://github.com/lesgourg/class_public/blob/master/source/background.c ; Python wrapper https://github.com/lesgourg/class_public/blob/master/python/classy.pyx .

## Actual test procedure in GitHub Actions

The scientific workflow thermal-bose-class.yml installs classy==3.4.3.0, constructs the Bose PSD numerically and runs the *linear* Einstein–Boltzmann solver at z=0 and z=5 over k=0.002 to 12 inverse Mpc. A cold reference is given the SAME physical non-baryonic omega today, but unlike the thermal model has no extra early radiation. The output therefore includes expansion-history as well as free-streaming differences.

Store both sqrt(P_Bose/P_CDM), a shape normalized at the smallest sampled k, a possible 0.5-amplitude crossing (if it exists **within** the simulated k range), and a second q-grid resolution check at z=5. Both the raw results and the independent analytic-only run are retained in an artifact for 21 days. It is critically important that this workflow produce **a real successful run** before making claims about CLASS-computed suppression.

## Limits and next follow-up

- This comparison does NOT include Lyman-alpha/CMB data likelihood or transfer-function half-mode observation limits.
- This baseline rules in neither a stable mixed Bose-gas galaxy nor an actual 86-eV particle detection.
- The primordial Bose-stimulated cold production, coherent condensate, self-interacting kinetic terms, and large-scale formation history need separate physical models; a standard thermal relic may have warm spectrum even at the same particle mass.
- A meaningful next paper must compare verified CLASS results with a physically derived colder f(q), keeping present abundance, temperature, chemical reactions and internal degrees-of-freedom consistent.

**Research status:** the [verified CLASS run](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/44890003670) completed successfully. See [actual transfer spectrum and numerical convergence](THERMAL-CLASS-VERIFIED.md), including the unresolved ultra-suppressed high-k numerical tail. An additional [nonthermal toy family](TOY-COLDER-PSD.md) explores dependence on primordial momentum shape without claiming any particle-production model.


### New verified nonthermal shape sensitivity scan

[Verified 86-eV toy colder PSD CLASS results](COLDER-CLASS-RESULTS.md) calculate four positive phase-space distributions with identical particle mass and number density but momentum scales alpha=1,0.5,0.3,0.1. They shift the linear z=5 half-amplitude k from 0.433 to 2.870 inverse Mpc, **without any claim** that such nonthermal distributions can be physically generated. The analytic distribution scaling was separately checked. Only the alpha=1 benchmark has a systematic fine/coarse PSD precision check at present.
