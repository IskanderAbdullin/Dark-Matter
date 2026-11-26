# Finite-mediator 5D Pauli heat: first implementation gate (2026-11-22)

## Prior hypothesis and falsification criteria

For a spacelike t-channel massive vector with fixed tree-level vector-current numerator and contact coefficient G=gχgℓ/mV², the *eventwise* squared-amplitude correction is F(t)=[mV²/(mV²−t)]², t=−2p*²(1−cosθ*)≤0. In the mV→∞ limit at fixed G, the full five-dimensional Pauli heat must recover the previously tested contact integral; at equal χ/lepton temperature it must vanish by detailed balance; a light mediator must change non-equilibrium heat. Failure of any check refutes the implementation.

## Implementation and limits

The independently boosted 5D center-of-mass integrator in `kinetics/independent_massive_chie_5d.py` now accepts optional positive finite `mediator_mass` (same energy units as M,m,T). It multiplies the **full angular-dependent vector trace**, not an angularly constant surrogate, by F(t), retaining both initial/final Fermi factors and energy transfer E′−E. The old contact path is unchanged when the option is omitted. Tests: `python -m pytest -q tests/test_finite_mediator_5d.py`.

The tests exercise a massless lepton benchmark M/T=0.1, T=1, Tχ/T=.9 or 1, G=.001, incoming χ/e kinetic cut 18T, 12- and 16-point Gaussian radial/angular/azimuth quadrature, mV/T=1 and 10⁶. Tests assert contact recovery, finite-mediator nontriviality, equal-temperature cancellation, positive heat for colder χ, and invalid-mass rejection. These are **test specifications, not measured numerical outcomes** until CI actually runs. There is no independent second finite-mediator 5D implementation or demonstrated angular/tail convergence yet.

The mV/T=1 case can have important high-momentum and forward-scattering structure; the EFT is not valid uniformly in the thermal tail. This is a tree-level vector numerator, not a gauge-complete χ1/χ2/φ embedding, thermal propagator, plasma calculation, realistic ν decoupling or an eigenvalue of the nonlinear operator. No γel, ΩDM, BEC, frozen spectrum or CLASS claim follows.

## Next independently verifiable experiment

Run the above test in real public CI; independently derive a t-integrated angular kernel **with the nonconstant spin trace** and compare against 5D at n=12,16,24,32, cutoff 18/26/36T, mediator mass ratios 1/3/10/100 and both χe (me>0) and χν (m=0) with split bath temperatures. Quantify numerical uncertainty, energy-flow signs, finite-mass electron thresholds and the contact-limit thermal-tail bound before physical interpretation.
