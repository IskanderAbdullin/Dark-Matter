# Empty-sector charged-lepton injection, prompt decay and a fermion transport stress test

This benchmark is a controlled source-only calculation in the authors' tabulated visible adiabat. It starts without an occupied dark sector, uses the finite-width vector with direct FD electrons, muons and taus, and injects both neutral charge-symmetric Dirac pairs, m1=.1 GeV and m2=.02 GeV. The heavy state immediately decays into a stable .02 GeV partner and an 86 eV scalar, using the spin-averaged physical two-body kinematics. All products subsequently redshift without reinteraction in this injection control. No actual reheating model, reheating likelihood, hadron spectrum, physical freeze-out, or condensate is inferred.

## Production and distribution normalization

The empty-final-state source is evaluated directly, with visible Boltzmann weights factored from the exact FD denominators. The spin-correlated vacuum production angle has only the zeroth and second Legendre moments after averaging over the incoming relative azimuth; the individual stable-particle energy retains that angular weight. In the MB limit compare against independent Bessel moments with adaptive invariant-energy integration, for both the heavy and light masses. Check sum-energy normalization and independent squared energy under the production angles, then check fermion and scalar count and energy at the decay. The value labeled relative to FD equilibrium refers to density for the *same visible temperature*, not an actual FD dark bath.

At TR=.10 GeV and T_SM=.09 GeV the dimensionless stable Y is about 4.06e-5 (both stable parents and both charge sectors). The direct and decay-fed contributions are nearly equal and the total count is 8.21e-4 of neutral FD at the visible temperature. At TR=.05→.045 GeV it is 8.84e-5; at .02→.018 GeV it is 2.58e-6. These are source-side counts, not a two-body relaxation calculation. A matched FD fit to the .10→.09 stable spectrum, using two moments, gives T_dark/T_SM=1.043 and mu/T=-7.33, but still underpredicts its momentum-squared mean by 7.4% relative to the fitted value. The same fit at .02→.018 gives T_dark/T_SM=1.289 and mu/T=-13.78. Do not plug these fits into a thermal operator as though they establish the real shape.

## Very soft scalar tail and elastic vector proxy

The prompt scalar from .10→.09 GeV has rms p/T_SM about 2.52, using approximate shell compression that preserves parent count and energy before the decay. The numerical soft tail uses an exact interval probability in the decay angle, not a histogram zero. A separate analytic upper bound uses the smallest allowed parent energy at the largest elastic-free redshifted cutoff, a maximum vacuum sigma-v, the Maxwell-Boltzmann upper bound on initial occupancies, and the integrated visible expansion. For the **diagnostic** p/T_SM<.005, at TR=.05→.045 the upper yield divided by the numerical total scalar yield is about 10^-63. At .10→.09 the same scalar fraction is bounded by about 5.5e-5 within this source-only model and finite input cosmology. These broad bounds do not include scalar reprocessing. A zero output from a finite thermal quadrature is not treated as an exact absence of physical tail.

Separately calculate the l=2 CM transport clock from same-charge stable fermion vector scattering, as a dilute sensitivity test. It has exact t/u exchange interference and is checked against explicit spinors; the identical final-state phase-space factor is included. For the three above examples it is around 4.60e3, 77.4 and .185 times H respectively. This proxy is not a relaxation eigenvalue or evidence for instantaneous kinetic equilibrium. At sufficiently high reheating the source-only phase-space shell occupancy itself exceeds the dilute gate, so no dilute transport clock is reported for those rows. Pauli feedback and actual time-dependent kinetic redistribution are then mandatory.

## Reproduce

Use the authors' central [SM EOS table](https://member.ipmu.jp/satoshi.shirai/DM2020/eos2020.dat), not redistributed. Compute the spectra locally to regenerate the detailed normalized bins behind the compact published summary.
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error kinetics/empty_lepton_pair_spectrum_reference.py EOS_TABLE > raw-spectrum.json
python -W error kinetics/prompt_scalar_cold_tail_envelope_reference.py EOS_TABLE raw-spectrum.json
python -W error kinetics/dilute_dark_fermion_transport_clock_reference.py EOS_TABLE raw-spectrum.json
```
