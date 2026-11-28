# Third-archive time refinement and a matched Hubble-clock control

The [fixed-scale five-channel history](INDEPENDENT-EARLY-POLE-JOINED-LOCAL-SPECTRAL-MEMORY-2026-11-27.md) on the third archive is refined from 24 to 48 macro steps at exposure 8e−6. Unlike soft polarization is refreshed as before. The two accumulated gains change by +0.1453% (step) and −0.0144% (continuous) to 7.1256021e−7 and 7.1403226e−7. The gap relative to the step gain is 0.207%, down from 0.367%. The final interval difference stays about 0.021%. All interval means after 5e−7 differ by less than 0.1%. Against separately matched controls without either elastic channel, the finished gains are +23.46% and +4.10%, so those relative elastic effects still depend strongly on reconstruction. In the common aperture [1.39,1.48], the parent and daughter occupation differences fall from 5.61%/18.87% to 0.00992%/0.0301%; across their full common support the final differences are 0.201%/0.0385%. The two trajectories being close locally does not identify every allowed initial profile.

The independent [comoving energy closure](INDEPENDENT-EARLY-POLE-JOINED-SCREENED-FRW-2026-11-26.md) is also run on third. In this experiment **both** polarization responses are held at the archived values; do not compare the clock change to the active-response longer result as if all settings were fixed. The same causal rate supports and cold window are used along the short expansion. With exposure 2e−6, the cold gains are 1.6126491e−7 at fixed scale, 1.6126792e−7 under expansion using a radiation clock, and 1.5709055e−7 when the modeled dark energy is also included in H. Thus the last gain is 2.5903% lower, compared to 2.6603% on fine. Expansion alone at fixed radiation clock changes the third gain by +0.001864%. Trace subtraction leaves comoving energy residuals below 1.3e−13, with entropy residuals below 1.7e−12.

A fourth control integrates using the radiation clock only, choosing the shorter interval (2\times10^{-6}/\sqrt{1+\rho_{\rm dark,0}/\rho_{\rm vis,0}}=1.9504520007884388\times10^{-6}). Its accumulated gain differs from the full-H result by only 7.35e−14 absolute (4.68e−7 relative), and the maximum difference in their 13 stored history points at equal fractional exposure is below 7.35e−14. The shorter-reference comparison changes the redshift interval along with the collision clock; agreement shows the short test is well described by the initial rate rescaling, not that a constant H correction is justified over a long cosmological episode.

Reproduce from `kinetics/` using `OPENBLAS_NUM_THREADS=1`:

```sh
for p in steps cell-conservative; do
  python early_pole_cold_screened_split_reference.py third --duration 8e-6 --count 48 \
    --same-soft-strength 1 --compton-strength 1 --conversion-strength 1 \
    --screening-update full --pair --projection "$p" --snapshots > "$p.json"
done
python early_pole_joined_history_memory_reference.py steps.json cell-conservative.json
python early_pole_joined_snapshot_aperture_reference.py steps.json cell-conservative.json
python early_pole_joined_screened_frw_reference.py third
python early_pole_joined_screened_frw_reference.py third --expands
python early_pole_joined_screened_frw_reference.py third --expands --full-H
python early_pole_joined_screened_frw_reference.py third --expands --duration 1.9504520007884388e-6
```

Raw: `results/early-pole-joined-third-time-and-frw-2026-11-27.json.gz`.
