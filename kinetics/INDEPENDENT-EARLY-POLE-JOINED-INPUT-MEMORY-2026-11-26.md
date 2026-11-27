# The longer joined fine-grid source loses much of its initial representation memory

The [joined finite model](INDEPENDENT-EARLY-POLE-HARD-FLAVOR-CONVERSION-2026-11-26.md)
contains photon pairs, hard Compton, hard pair conversion, unlike
screened elastic, and the same-flavor screened pole, together
with the cold scalar mode. At fixed scale and reference exposure
`8e-6`, compare on the *same fine hard mesh* the archived step
occupations and the count-conserving continuous reconstruction.
The unlike response is refreshed actively, with frozen exterior;
the same-flavor soft response remains held. This is a comparison
inside the same finite box, not a claim of a complete continuum rate.

| Fine archive input | Photon pairs, Compton, conversion only | With both screened elastic sectors | Relative change |
| --- | ---: | ---: | ---: |
| Native steps | 6.26338e-7 | 6.93480e-7 | +10.72% |
| Count-conserving continuous | 4.52824e-7 | 6.91367e-7 | +52.68% |

Both rows start with the same cold scalar counts. The absolute
accumulated gains with elastic included differ by about 0.30%,
even though the two corresponding no-elastic gains differ by
about 27.7% relative to the step result. Small variation in the
joined *absolute* response does not mean convergence of its
change relative to each representation's control.

The first joined mean scalar source per macro interval is
about .07234 and .06069, respectively; by the last interval it
is .096362 and .096383. The first two macro intervals already
produce 99.36% of the eventual accumulated difference. The
remaining rates have both signs of small representation
difference, and neither approaches zero in this interval.
This is a memory test of one observable on one spatially
homogeneous fine-grid model, not a proof that the full charged
spectra coincide or that the true collision dynamics has
erased early-universe memory.

After the first four macro steps, all subsequent interval
averages differ by less than 0.1% of the step-input rate.
The actual energy defect in the two joined runs stays below
`1.6e-13` and the entropy mismatch below `1.6e-12`.
A matched 24-step result changes the two absolute gains by
0.052% and 0.314%, respectively. The gap between the two initial
states is 0.570% at 24 steps and 0.306% at 48, while their two
last-interval rates differ by about 0.020% at both resolutions.
Thus the small absolute gap is qualitatively robust but its
precise sub-percent value is not temporally converged. No
independent hard-mesh convergence is claimed.

Reproduce from `kinetics/` on one BLAS thread, changing
`steps` to `cell-conservative` for the second initial state:

```sh
python early_pole_cold_screened_split_reference.py fine --duration 8e-6 --count 48 --same-soft-strength 1 --compton-strength 1 --conversion-strength 1 --pair --screening-update full --projection steps
python early_pole_cold_screened_split_reference.py fine --duration 8e-6 --count 4 --strength 0 --compton-strength 1 --conversion-strength 1 --pair --projection steps
```

The history comparison command takes two matched time series.
It compares finite-interval source averages; with optional
occupation snapshots, it integrates a relative (L^1) norm
of the actual reconstructed shapes over a common radial support.

Raw: `results/early-pole-joined-input-memory-2026-11-26.json.gz`.
24-step control: `results/early-pole-joined-input-memory-time-check-2026-11-26.json.gz`.
