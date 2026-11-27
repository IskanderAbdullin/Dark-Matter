# Crossing scalar refresh with two matched archived states

Complete the [screened scalar self cross-controls](INDEPENDENT-EARLY-POLE-SCALAR-SCREENED-TIME-REFRESH-2026-11-27.md) at 2e−8 exposure, eight macro steps and sixteen scalar momentum bins. In each archive the same prepared 97% parent fraction feeds the two hard screened terms and the finite slow channels. With the scalar self proxy disabled, the full held and refreshed output trajectories agree exactly, including the starting spectra and every ledger sample.

| Archive | No-self integrated gain | With held scalar term | With refreshed scalar term |
| --- | ---: | ---: | ---: |
| Fine | 9.691490653e−10 | 9.698804879e−10 (+0.07547%) | 9.698747006e−10 (+0.07487%) |
| Third | 1.608614961e−9 | 1.607208387e−9 (−0.08744%) | 1.607217374e−9 (−0.08688%) |

These two finite snapshots show that the sign difference persists with either scalar response choice. Updating the illustrative scalar mass and vertex changes the integrated correction by less than 0.00060 percentage points in these pairs. The final instantaneous rates, all positive, show the same small signs of the *corrections*. Neither these corrections nor the larger discrepancy induced by temporal refinement on fine are a precision determination of the physical elastic rate. The archive states, radial meshes, and initial decay apertures are different; this cross-control does not isolate which of those archive differences causes the opposite signs. The absolute energy and entropy ledger closures are below 1e−14 and 5e−13.

From `kinetics/`, with one BLAS thread:
```sh
export OPENBLAS_NUM_THREADS=1
for g in fine third; do
  for mode in held refresh; do
    opts=()
    if [ "$mode" = refresh ]; then opts=(--refresh); fi
    python early_pole_joined_scalar_self_proxy_reference.py "$g" \
      --duration 2e-8 --steps 8 --hard-slow --screened --mixes .97 "${opts[@]}"
  done
done
```

Raw: `results/early-pole-scalar-screened-cross-controls-2026-11-27.json.gz`.
