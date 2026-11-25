# Actual late weak-parent CLASS transfer and tail underflow (2026-11-25)

The late weak-coupling x0=0.1 initially empty parent candidate has M1=86 TeV, m2/M1=0.9, m_phi=86 eV, kappa_SM1=47, kappa2=.03, kappa=.005278180162600994. The old 480/960-node quadrature passed moment tests but CLASS failed with nonfinite state at x=2.036219e-5.

## Root cause and correction
The published classy-community 26.0.0 source (https://pypi.org/project/classy-community/26.0.0/), species/ncdm_base_species.cpp lines 237-266 and 296-326, computes q/f0*df0dq and extrapolates PSD from its last two rows. The old cold and middle tracer at the 480-point node q=480 had f0 ~2.73e-319: nonzero but q/f0 overflows; at 960 nodes there was also underflow to zero. A negligible floor proportional to 1e-100 exp(-u/50) instead of exp(-u) avoids this arithmetic pathology; added a reciprocal and extrapolation guard.

The old third scale selected a number-weighted shoulder at q~0.015. Select the peak of q^5 f (second momentum moment per log interval), option --energy-peak: hot scales 1.139,1.207; all four mesh/partition combinations pass at 120 nodes (max error 7.43e-5), and at 240 nodes (<4.01e-6). The 320-cell reconstructed number/first/second momentum biases remain -0.1905%/-0.3189%/-0.4961%. No abundance rescaling.

## Completed calculations
Windows portable Python 3.12.10 / classy-community 26.0.0 / NumPy 2.2.6 / SciPy 1.15.3. Command: python cosmology/daughter_bath_class.py --input outputs/yukawa_lateempty47_240.json outputs/yukawa_lateempty47_320.json --energy-peak --bins 120 --kmax 20 --output outputs/lateweak_transfer260_safe.json

Actual synchronous CLASS P/P_CDM(z=5,k=10/Mpc): n=240 0.985179558, n=320 0.984532708. Fine minimum 0.975671687; no half-power crossing up to 20/Mpc; k=1 ratio 0.999846335. Maximum independent partition/quadrature drift 0.235% (coarse), 0.0633% (fine). Cross-mesh drifts are in the JSON. Density, relative species fractions, radiation normalization and large-scale guards passed. These are conditional adiabatic collisionless linear calculations, not observational likelihoods or proof of physical viability. Full hierarchy control remains in progress.

[Raw output, failed run and source-based screen](results/late-weak-class-verified-2026-11-25.json). Targeted regressions 5/5 passed. Next: momentum-dependent parent with dark depletion and independent collision energy closure.
