# Archived moments do not fix the sign of the initial bare decay source

This is a finite archived-spectrum counterexample to sign identification from three native hard moments, not a prediction for the later coupled trajectories. Take the published (x=.05), daughter mass ratio .06 fine and third snapshots, bare charged masses (.05 and .003), and each archived scalar Hartree mass. Hold the daughter and all scalar nodes fixed. Replace only one parent cell on fine, two on third, with explicit bounded subinterval occupations. Their number, bare energy and reciprocal gauge moment are equal cell by cell to archive values. A centered correction below (1.2\times10^{-12}) removes roundoff in the serialized profile; the unnormalized integral defects are below (4\times10^{-18}).

On each scalar node, integrate the exact parent-energy interval, split at both the parent edges and the daughter edges shifted by the scalar energy. If (N) and (V) are the archived number and shell phase volume, the contribution is ([N(I_P-I_b)+V I_{\rm seed}]/(E_s q)), omitting one common positive factor. Here (I_P) and (I_b) are the two individual occupation integrals and (I_{\rm seed}) integrates (f_P(1-f_b)). An independent radial integral, segmented on the same physical edges and using eight-point Gaussian quadrature for (p/E_P), agrees with the energy integrals within (1.8\times10^{-16}) absolute.

| Native archive | Step source | 85% mixed witness source | Zero along this mixture |
| --- | ---: | ---: | ---: |
| Fine | +6385.709 | −1345.206 | 70.210% |
| Third | +5217.779 | −1497.459 | 66.045% |

The spontaneous contribution is positive and stays near 0.0384 / 0.0376, respectively. The sign reversal is a change in the stimulated imbalance and persists when restricted to (q<10^{-5}). The witness is not the outcome of any evolution, and the zero is specific to this one-dimensional mixture, not the minimum perturbation needed for a reversal. Full screened dispersion, changing scalar and daughter distributions, finite width and elastic smoothing are not included. The previously observed five-channel positive finite integrated gain applied to the two tested initial reconstructions; agreement between them does not determine the sign over every profile sharing the same archived moments. Compare the parent and daughter aperture integrals separately when examining live or reprojected snapshots.

Reproduce from `kinetics/`:

```bash
OPENBLAS_NUM_THREADS=1 python early_pole_native_aperture_witness_reference.py \
  > /tmp/native-aperture-witness.json
cmp /tmp/native-aperture-witness.json results/early-pole-native-aperture-witness-2026-11-27.json
```
