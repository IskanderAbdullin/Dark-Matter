# Timelike pole number and photon heat on bounded momentum bands

Extend the [full-angle coefficients](INDEPENDENT-FINITE-MASS-QED-TIMELIKE-POLE-ANGULAR-2026-11-28.md) at fixed scalar q=4.1e−6 to two finite intervals in pole momentum. The transverse interval is [.0012,.003] and the longitudinal interval is [.00085,.003]. Both include the measured low-support openings, with verified zero panels below their respective past-null apertures. Neither interval covers all possible medium recoil or the full longitudinal branch. No distribution of pole quanta is assumed by the integrator.

## Measure and convergence

For each k, average both physical transverse polarizations or the conserved density projection in the Coulomb channel, including the positive-frequency retarded residue R in ω and the exact E-bin angular/energy support. Use dc/2 and the additional k² dk measure. The tabulated coefficients omit one common factor, 1/(32π³ E_s), as well as explicit g² y² and a separate charge-conjugate contribution. This factor follows from integrating the daughter momentum and the parent polar angle in the invariant Boltzmann phase space: the charged measure after azimuth integration is dE/(8πV), the photon pole measure is k² dk dc R/(4π²), and the fixed scalar projection supplies 1/(2 E_s). Unnormalizing a source divided by n_s also restores that factor of n_s. This is a single-charge orientation on the source side.

At each k, the real charged aperture opens when the separation of the on-shell daughter and parent spatial momenta equals V. The outer mesh splits the roots of |p_b−p_P| = k±q for both histograms and both ends of the cosine range, and the corresponding past-null conditions. Two independent bracketing grids must find the same outer roots. Each outer panel is complete, with an independently refined k/c/E/azimuth quadrature. The inverse and forward sides are accumulated separately even when one vanishes during opening. The raw receipts are retained as deterministic compressed JSONL; their indices, scope, contiguous coverage, nonnegative forward/inverse coefficients and signed ledgers are checked by the reducer.

## Conditional number and heat

Let A and B denote the integrated positive forward and inverse coefficients. For one constant trial occupation n inside the named band, the scalar number bracket normalized by its frozen n_s is n[(1+1/n_s)A−B]−B. Its sign changes at the number threshold shown below. The **photon bath** gains heat in the opposite direction to the analogous bracket with A and B separately weighted by ω(k). It has its own balancing occupation. This is not a constant trial temperature or a prescription for the actual population n(k).

No fitted heat model is needed for the bounds: the pole frequency is monotonically increasing in these frozen branches. For transverse modes the positive even series for Π_T and the factor 1/(1−β²) make ω² increasing in β=k/ω. For longitudinal modes ω²=Π_00(β)/β² is a positive power series in β². Multiplying each positive integrated panel by the two endpoint frequencies therefore brackets its exact photon-energy moment. A refinement check controls numerical panel error separately from these endpoint bounds. Since the heat threshold increases with the inverse coefficient and decreases with the forward coefficient, pairing opposite moment bounds gives a conservative interval in this model.

| Archive and mode | number-neutral n | photon-heat-neutral interval | maximum coefficient mesh difference |
| --- | ---: | ---: | ---: |
| fine T | 5.53982857e-08 | [5.49827302e-08, 5.64483041e-08] | 1.86e-06 |
| fine L | 7.18220553e-07 | [7.18589399e-07, 7.34571834e-07] | 5.43e-07 |
| third T | 5.36482656e-08 | [5.3507639e-08, 5.44553997e-08] | 1.86e-07 |
| third L | 6.87055473e-07 | [6.89291405e-07, 7.00811312e-07] | 1.91e-07 |

In the bounded longitudinal tests the lower heat-balance occupation exceeds the number-balance value on both archives by more than numerical mesh differences. At that **trial** number-neutral occupation, the scalar has no net number or scalar energy change at this fixed q, while the photon band still receives net energy. This does not establish a pole temperature or a physical heat flow for the unknown actual occupation. The corresponding transverse heat intervals overlap their number thresholds in the tested bounds.

## Reproduction

From `kinetics/` and with `OPENBLAS_NUM_THREADS=1`, run
```sh
mkdir -p pole-band-work
for g in fine third; do
  for p in T L; do
    python -W error finite_mass_qed_timelike_pole_momentum_band_reference.py \
      "$g" "$p" "pole-band-work/pole-band-$g-${p}0.jsonl"
    python -W error finite_mass_qed_timelike_pole_momentum_band_reference.py \
      "$g" "$p" "pole-band-work/pole-band-$g-${p}1.jsonl" --refine
  done
done
python -W error finite_mass_qed_timelike_pole_momentum_band_receipt.py pole-band-work
python -W error finite_mass_qed_timelike_pole_momentum_band_receipt.py results --gzip-input
```
The [summary](results/finite-mass-qed-timelike-pole-band-2026-11-28.json) reports coefficient units, forward-only support samples and the averaged photon-energy intervals. No virtual or spacelike matched term, all-q projection, cosmological solution, or populated collective spectrum is inferred.
