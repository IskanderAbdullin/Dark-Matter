# A simultaneous held-screening solve checks the split-width extrapolation

Evolve the same finite cold and hard state on the first original time interval, 3.3333333333333335e−7. Keep both elastic polarization tables fixed. Instead of alternating the slow and elastic flows, sum their instantaneous moment derivatives and integrate them together. Represent each charged cell by two FD logits; translate its moment derivative with the exact local susceptibility matrix. This keeps charged trial occupations in their admissible chart without requiring a physical claim about the held plasma model. The cold and photon variables remain in their original number/moment form.

The nominal DOP853 chart runs at rtol 3e−9 give mean cold rates .06649893898 (native) and .06552348262 (continuous), a 1.467% contrast relative to native. Both pass the unchanged total energy, entropy and hard-count closure checks; their energy defects are below 7e−16. They lie close to the previously reported [1.468% leading-linear extrapolation](INDEPENDENT-EARLY-POLE-FIRST-WINDOW-TIME-LOCAL-MEMORY-2026-11-27.md), which was estimated from time splits with active unlike response. The explicitly matched held-screening 64-step split paths give a direct finite-width control: the chart is -0.1690% and +0.0886% relative to the native and continuous matched split results, respectively. Tighter chart tolerances are a separate needed verification; good energy closure alone does not bound the error in cold gain or establish continuum convergence or physical validity.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
b=results/early-pole-joined-fine-216-memory-2026-11-27.json.gz
for p in steps cell-conservative; do
  python -W error early_pole_joined_simultaneous_chart_reference.py "$b" --projection "$p" > "$p-chart.json"
  python -W error early_pole_cold_screened_split_reference.py fine \
    --duration 3.3333333333333335e-7 --count 64 --pair --same-soft-strength 1 --compton-strength 1 \
    --conversion-strength 1 --screening-update frozen --projection "$p" --snapshots \
    --elastic-energy-budget 5e-15 > "$p-held64.json"
  python -W error early_pole_joined_frozen_split_chart_audit_reference.py "$p-chart.json" "$p-held64.json"
done
```

Raw: `results/early-pole-first-window-unsplit-chart-2026-11-27.json.gz`.
