# First sampled sign change and a more finely resolved joint trajectory

Continue the [five-channel prepared fine witness](INDEPENDENT-EARLY-POLE-WITNESS-ELASTIC-SECTOR-ATTRIBUTION-2026-11-27.md) with held same-flavor screening and an active unlike response. The quoted numbers are reference exposures, not solver wall time.

On the three-subdivision, 16-scalar-bin mesh, run two increasingly short eight-step trajectories from the same moments. The longer one starts at −.00682279 and is already positive at 2.5e−10 (+.00160091). A fourfold time refinement measures −.00442342 at 6.25e−11, −.00221955 at 1.25e−10, and +.00169772 at 2.5e−10. Thus the earliest observed sign change is bracketed by the latter two sampled times on this one hard/soft projection. The spatial refinements below do not fix a hard-converged crossing time, and no interpolation between the two sign samples is needed for the bracket.

A separately refined joined trajectory uses eight hard subdivisions and 32 scalar bins (two enhancements at once), compared with six and 16 in the preceding finite control. At exposure 2e−8 the refined cumulative scalar gain is +8.7663e−10, about 2.14% above the six/16 result. Instantaneous values move from −.00659757 to +.0219049 by the first short sample at 2.5e−9 and end at +.0583905. The whole refined trajectory completes with energy and entropy residuals below 5e−16 and 5e−13, respectively. It was resumed from an exact measured macro boundary using the [checked restart protocol](INDEPENDENT-EARLY-POLE-SPLIT-RESTART-2026-11-27.md); all eight macro increments and both ledgers are present in the raw data.

These calculations support rapid sign relaxation in the tested prepared state. The impulse magnitude and first-crossing bracket are finite-model results; the two separate mesh sensitivities, role of scalar self-scattering, and physically populated initial spectra are not determined here.

From `kinetics/`:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MALLOC_ARENA_MAX=2
# Compare (duration, hard subdivisions, scalar bins):
# (2e-9,3,16), (5e-10,3,16), and (2e-8,8,32).
T=2e-8; H=8; B=32
python early_pole_cold_screened_split_reference.py fine --duration "$T" --count 8 \
  --same-soft-strength 1 --compton-strength 1 --conversion-strength 1 --pair \
  --screening-update full --snapshots --aperture-subdivision "$H" --scalar-bins "$B" \
  --witness-mix .97 --response-radial-order 32
```

Raw: `results/early-pole-witness-onset-higher-resolution-2026-11-27.json.gz`.
