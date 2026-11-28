# Density, pressure and pseudo-pressure on a frozen shape

For a decoupled distribution the background weights differ for energy density, pressure, and the higher moment called pseudo-pressure. See [Lesgourgues & Tram, CLASS IV, equations (14)–(15)](https://arxiv.org/html/1104.2935). This reference motivates the **choice of diagnostic moments**, not a use of CLASS or an inference about this project's physical relic.

Take the native W1, W2 and independent warm Gauss2 spectra at x=.05 with identical cold and charged charts. Retain each spectrum's shape while varying a formal positive mass parameter (M) in (E=sqrt{q^2+M^2}). The measures use weights (1, E, q^2/(3E), q^4/(3E^3)) against the same (q^2 f(q),dq). No mapping from (M) to this model's time or temperature is supplied; the large-(M) rows test the weights, not a physical propagation to late times. For Gauss, fit the positive representative once at the measured W1 mass, then hold it fixed for **all** test masses. It agrees at that mass with the separate terminal pressure result.

Using a centered independent dilation, (Mmapsto M e^t) and (a^{-4}) factors for the energy and stress, verify (dho/dln a=-3(ho+p)) and (dp/dln a=-5p+mathfrak p) for these formal frozen measures. Number uses (a^{-3}). These are background identities of the frozen test, not continuity or decoupling tests of the interacting histories.

At the original fixed reference mass (8.73183	imes10^{-4}), Gauss2 versus W1 changes energy by −0.1042%, pressure by −8.4489% and pseudo-pressure by −8.4609%. W1's q≥.01 tail is 9.288 particles per million, but carries 99.855% of pressure and 99.997% of pseudo-pressure at that mass. At a deliberately heavy test mass (M=100), Gauss2 versus W1 changes pressure by −12.559% and pseudo-pressure by −21.769% while the energy difference is still only −4.26e−7 relative. Similar abundance never certifies accuracy of higher stress moments. The 24/48-node integrals, independent dilation differences, and massless/heavy unit tests check the moment calculation; none bounds the numerical error in the input spectra.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results/early-pole-fine-charge2-scalar4; d=2026-11-28
python -W error -m unittest test_early_pole_endpoint_frozen_stress_sweep
python -W error early_pole_endpoint_frozen_stress_sweep.py \
  "$p-x05-state-$d.json.gz" "$p-warm2-x05-state-$d.json.gz" \
  --gauss "$p-gausswarm-x05-state-$d.json.gz"
```

The (M=0) limit is checked directly. In the high-mass limit, (3MP) and (3M^3mathfrak p) approach the (q^2) and (q^4) count-weighted moments; the geometric histograms also have independent analytic polynomial integrals. The Gauss2 limit contrasts against W1 are −12.566% and −21.818%. At the finite diagnostic (M=10), the adiabatic derivative drops from (1.24444	imes10^{-7}) to (1.09582	imes10^{-7}), about −11.94% relative but only (-1.4862	imes10^{-8}) absolute. This derivative is not an effective perturbation sound speed or a transfer prediction.
