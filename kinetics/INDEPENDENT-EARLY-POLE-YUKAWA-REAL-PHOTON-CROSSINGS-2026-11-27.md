# Frozen real-photon crossings of the flavor-changing scalar Yukawa process

The existing cold decay and hard scattering references do not include the vacuum Yukawa process with an external real photon. For a neutral scalar and two equal-charge Dirac flavors, compute both diagrams for $\gamma+P\leftrightarrow b+\phi$ and $\gamma+b\leftrightarrow P+\phi$, and the crossed $\gamma+\phi\leftrightarrow b+\bar P$ (with charge conjugation). The same-mass topology agrees with the two diagrams in [Wang et al., Chin. Phys. C 45, 093001](https://doi.org/10.1088/1674-1137/ac0c0f); here both charged masses are distinct. The model uses a scalar Yukawa bilinear, not a pseudoscalar, and undressed Dirac and photon propagators.

The summed squared amplitude keeps interference. Independent on-shell invariant reduction reproduces the explicit gamma-matrix trace, with relative defects below 3e−10 over the fixture; the massless limit is $-4t^2/(su)$ before initial spin/polarization averages and coupling factors. Transverse sums, covariant polarization sums, explicit spinor Ward contractions, and crossed-spin traces all run independently. Lab and CM phase coordinates are checked by separately boosting four-vectors and comparing their invariant products. Forward and inverse quantum weights are evaluated on the same packets and satisfy common-temperature detailed balance in thermal fixtures, including Pauli blocking and photon stimulation.

For the two charged-photon absorption channels, $s+m_\phi^2-m_{\rm out}^2=2K\cdot q$ with $K$ the total incoming four-vector. On $.45<p_{P,b}<5$ and $.5<k_\gamma<4$, $q<10^{-5}$ requires $s<.0159$ or $s<.0184$, respectively. Thus the $s\ge1$ hard guard used elsewhere cannot measure **direct** injection into this cold range. The following finite-support calculation explicitly retains the small-$s$ space with a hard **plasma-frame** photon, including the near-collinear massive denominators. It does not generate a cutoff-independent soft-photon correction.

| Archive at q=4.1e−6 | P+γ stimulated | b+γ stimulated | Crossed pair stimulated | Sum relative to bare decay at g=.005 |
| --- | ---: | ---: | ---: | ---: |
| Fine | +25.9052 | −.59631 | −3.30868 | +6.78e−6 |
| Third | +28.4286 | −.60320 | −3.23947 | +8.83e−6 |

The three middle columns omit $g^2y^2$, in the same units as the matched decay coefficient with $y^2$ omitted before forming the last ratio. On both archives the sampled 10→16 node change is below 0.1% per stimulated channel with pole-aware panels. The two associated production channels independently have a positive combined source at $q=2\times10^{-6}$ and $8\times10^{-6}$ on fine; the crossed pair correction has only been computed at the central mode here. The finite external-photon contribution is small **in this frozen cut calculation**. Medium dispersion and residues, soft photon modes, virtual corrections, separate 1↔3 real emission, scattering-assisted coherent vertices, and a joint nonequilibrium bath were not calculated. Neither an inclusive correction, a freeze-out condition, nor a thermal bound follows from the quoted fraction.

From `kinetics/`:

```sh
export OPENBLAS_NUM_THREADS=1
python early_pole_hard_yukawa_photoproduction_reference.py
for g in fine third; do
  for n in 5 10 16; do
    python early_pole_frozen_yukawa_photoproduction_source_reference.py "$g" --order "$n"
    python early_pole_frozen_yukawa_photon_pair_crossing_reference.py "$g" --order "$n"
  done
  python early_pole_frozen_yukawa_photoproduction_source_reference.py "$g" --order 5 --trace
done
```

Raw results: `results/early-pole-frozen-yukawa-real-photon-crossings-2026-11-27.json.gz`.
