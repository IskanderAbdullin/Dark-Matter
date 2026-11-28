# Most of the measured fine past-support response lies near the upper cutoff

Hold the original fine endpoint and pole fixed. Use the exact 108 paired bridge panel edges, plus the previously verified low block, to cumulate the gross forward and inverse sources from q=1e−6. No integral is interpolated inside a panel. An arbitrary query cutoff is enclosed by the two neighboring exact cumulative totals.

For **number** on q=[1e−6,.0015], both the forward and inverse median lie in [.00119894,.00126568]. Both 90% crossing points lie in [.00145286,.00148614]; both 99% crossing points are only bounded by [.00148614,.0015]. The charged-energy gross medians are earlier, in [.00113148,.00117978], because they carry a smaller weight as the scalar energy approaches the fixed pole energy. At q=.001, at most 21.46% of the joined direct number response and 15.61% of the inverse number response have accumulated.

The two independent meshes give exactly the same native quantile brackets and a maximum absolute cumulative discrepancy of 2.42e−6. Three small-cutoff conditional roots fail the separate 6e−5 relative mesh check and are deliberately omitted from the root profile; the whole-band endpoint root passes. The brackets quantify the **verified truncated** band, and do not locate the global median, infer a peak outside the current support or supply a cosmological occupation. The end-of-band concentration makes the adjacent .0015–.004 transfer interval especially important to measure rather than extrapolate.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28; p=results
python -W error finite_mass_qed_native_past_cutoff_check.py --root "$p" \
 --low "$p/finite-mass-qed-fixed-pole-cold-to-shoulder-gap-$d.json" \
 --cold "$p/finite-mass-qed-fixed-pole-cold-scalar-band-$d.json" \
 --fixture "$p/finite-mass-qed-native-fine-past-bridge-$d.json"
```
