# Most of the first count difference is generated in the first quarter of that window

On the same first original 3.3333333333333335e−7 interval and exact initial states, refine the [accepted 32-substep comparison](INDEPENDENT-EARLY-POLE-FIRST-WINDOW-DYADIC-32-REFINEMENT-2026-11-27.md) to 64 substeps, again with the 5e−15 elastic acceptance. The native and continuous mean cold rates are .06661149218 and .06546547131 at 64, a gap of 1.7205% relative to native, compared to 1.9720% at 32. The 16/32/64 observed powers are .985 and 1.007. Assuming the leading error is linear in width gives an indicative gap of 1.468%, close to the earlier 1.470% estimate. This is a numerical diagnostic, not an exact limit or an extrapolation in physical time.

In the 64-step history, the native branch's excess has accumulated 69.4% of its own terminal contrast within the first eighth of the interval, 81.1% within the first quarter and 91.4% within the first half. All subsequent relative rate gaps within this first window are below 1% after 9.375e−8 (the same threshold at 32 steps). The instantaneous end-rate gap is only .171%, much smaller than the mean across the whole interval. Refining the time grid therefore clarifies a short initial transient rather than establishing equality of the underlying spectra. This is not the cold fraction of a physical cosmological solution.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
for p in steps cell-conservative; do
  python -W error early_pole_cold_screened_split_reference.py fine \
    --duration 3.3333333333333335e-7 --count 64 --pair --same-soft-strength 1 --compton-strength 1 \
    --conversion-strength 1 --screening-update full --projection "$p" --snapshots \
    --elastic-energy-budget 5e-15 > "$p-first64.json"
done
python -W error early_pole_first_interval_memory_reference.py steps-first64.json cell-conservative-first64.json
```

The [dyadic comparison command](INDEPENDENT-EARLY-POLE-FIRST-WINDOW-DYADIC-32-REFINEMENT-2026-11-27.md) can use the 16-, 32- and 64-step histories in place of 8, 16 and 32. Raw: `results/early-pole-first-window-64-memory-2026-11-27.json.gz`.
