# Parent production and thermal-mass benchmark

At x0=0.01, initially empty parent, n=160, M1=86 TeV, M2/M1=0.9, kappa=0.00224, kappa2=0.03, the number-relaxation model returns:

| kappa_SM1 | initial rate/H | Omega h^2 | p/T>0.1 fraction | peak q (Tref) |
|---:|---:|---:|---:|---:|
|47|0.1181|0.0835428|0.00696694|9.346e-6|
|400|1.0054|0.1416198|0.00242402|3.432e-6|
|4000|10.0539|0.1572799|0.00193002|2.458e-6|

The first coupling maps to a nominal uncancelled 0.05024 eV seesaw contribution in the specified width convention. It underproduces at the previous width. Width bracketing is in progress. Initial empty occupation is not equivalent to an actual non-equilibrium momentum-resolved SM parent source: the parent is still constrained to a thermal momentum shape with chemical potential. Low initial rate/H makes this assumption particularly uncertain.

## Thermal and naturalness screen

Perez-Gonzalez, Sen, Tangarife, *Bose-enhanced Neutrino Decays in a Thermal Medium* (12 June 2026), https://arxiv.org/html/2606.14441v1, Eq. 18, use the simple-Yukawa equilibrium benchmarks delta m_phi^2=g^2 T^2/6 and delta m_f^2=g^2 T^2/16. We evaluate those coefficients as a diagnostic, not as the actual thermal self-energies of an ALP portal. The kinetic collision integral already includes statistical factors; multiplying by a published thermal width would double-count them. Equilibrium scalar thermal mass/T=1.09235e-6 in these examples. Its maximum ratio to the production-epoch peak p/T is 0.00354 for the thermal-parent x0=0.1 candidate, and 0.117, 0.316, 0.438 for the empty x0=0.01 cases. The earlier spectra are NOT robust against thermal dispersion. A frozen-PSD variance proxy is recorded separately, with no claim to be a time-dependent self-energy. The unprotected one-loop mass scale g M/(4 pi) is about 18 MeV versus m_phi=86 eV; a symmetry or other UV protection is required for this illustrative coupling map, not established by the calculation.

Data: results/empty-parent-thermal-screen-2026-11-25.json. Commands: see yukawa_rate_bath.py --initial-parent empty and thermal_consistency_screen.py.

## Numerical regression

Original test suite: 15 passed, one pre-existing subcritical halo integrator overflow on Python 3.9.7/SciPy 1.13.1 (unbounded automatic initial-step estimate). Explicit first_step=0.01 fixes it with the original physical parameters and tolerances; 16/16 tests passed.
