# Scalar rebinning, hard reprojection and dressed kinematics are distinct

A follow-up to the [bare native-profile witness](INDEPENDENT-EARLY-POLE-NATIVE-APERTURE-WITNESS-2026-11-27.md). Hold fixed the same physical profile, parent-only perturbation and archived daughter. Replace the original scalar nodes by the eight logarithmic cells used by the live split reference on [1e−6,1e−5]. Compare three initial source evaluations at that same scalar resolution:

1. Exact archived step/substep occupations integrated over the bare Lorentz-invariant aperture.
2. Projected/reconstructed hard occupations integrated over the same bare aperture.
3. The same reconstructed populations evaluated with the live variational gauge energy gradients and dressed aperture.

Use 1, 2, 4, 8 or 16 uniform subdivisions per targeted cell of the production aperture mesh. Projection retains integrated bare parent number and energy at better than 3e−16 absolute for this fixture. The original continuous subcell witness preserves the native reciprocal moment as well; the live two-moment reconstruction is not constrained to preserve it. The modeled hard-field energies of the baseline and mixture at 16 subdivisions differ by less than 3e−15. We stop the mixture at 97% so all live hard cells remain in the strict interior of their realizable Fermi moments.

| Initial mixture, 16 subdivisions | Bare raw, rebinned scalar | Bare projected | Dressed projected |
| --- | ---: | ---: | ---: |
| Fine, unchanged step | +.076759 | +.076759 | +.078635 |
| Fine, 85% witness | +.002290 | +.003079 | +.010028 |
| Fine, 97% witness | −.008224 | −.007259 | +.000328 |
| Third, unchanged step | +.063873 | +.063873 | +.072574 |
| Third, 85% witness | −.004958 | −.004156 | +.035483 |
| Third, 97% witness | −.014675 | −.014012 | +.030451 |

The original whole-grid bare sign reversal at 85% does **not** survive scalar rebinning alone on fine, but does survive on third. Under bare dynamics the native hard closure conceals the 97% reversal on both grids; increased hard resolution retains it for both at 8 and 16 subdivisions. In the dressed model all quoted cases are positive: the fine 97% example is particularly near cancellation, fully unresolved as to sign on further refinement. Third can even move in the wrong direction for two subdivisions, illustrating that the sampled agreement with the bare witness under refinement need not be monotone. These are instantaneous reconstructed rates only, not evidence about how the witness relaxes with elastic scattering, the sign throughout a joined evolution, or the sign in a complete in-medium collision operator.

From `kinetics/`:

```sh
for grid in fine third; do
  OPENBLAS_NUM_THREADS=1 python early_pole_live_aperture_witness_reprojection_reference.py "$grid" \
     --subdivisions 1 2 4 8 16 > "$grid-witness.json"
done
```

The raw results are concatenated in `results/early-pole-aperture-witness-reprojection-2026-11-27.json`.
