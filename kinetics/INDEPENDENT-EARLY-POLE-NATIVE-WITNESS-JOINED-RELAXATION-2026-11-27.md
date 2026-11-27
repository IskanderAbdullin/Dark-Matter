# The initial inverse-dominated witness becomes a source when the elastic sectors are included

Start from the fine 97% parent witness with 16 scalar bins and three targeted subdivisions of the hard aperture. At precisely these same initial moments and grids, compare the no-elastic operator (cold decay/inverse, photon pairs, hard Compton and hard flavor conversion) with the operator that also includes unlike screened elastic and the same-flavor pole. The latter condition refreshes the unlike response, with a held exterior and held same-flavor polarization. No scalar or hard profile is edited during the trajectory.

| Initial parent profile | No elastic gain, exposure 8e−8 | Both elastic sectors included |
| --- | ---: | ---: |
| Native steps | +6.3398834e−9 | +5.1492184e−9 |
| 97% witness | −5.4486888e−10 | +4.4366386e−9 |

At the initial witness moments, the scalar source is −.00682279 on both models. Without the elastic sectors all sampled sources remain negative, with final source −.00679894. With the elastic sectors the first stored sample, at 1e−8, is already positive (+.049067), and the final rate is +.062965, within 0.15% of the native-step rate +.062877. The witness still accumulates about 14% fewer scalars over the whole interval. On an independent short time grid with fourfold smaller macro steps, the witness source is positive already at the first sample 2.5e−9 (+.030263); we have not resolved the exact crossing time. Energy and independent entropy ledgers on these runs are below 2e−15 and 5e−12.

This comparison identifies the change upon adding the **two modeled elastic sectors together**, not their individual contributions or a complete plasma rate. A separate pure decay/inverse run with 32 hard subdivisions and 32 scalar bins still gives a negative net gain for this witness (−2.05976e−8 at 2e−6). The joint hard-mesh convergence of the faster sign reversal and its persistence on other prepared or physical inputs remain open.

The optional active-response radial quadrature is separate from the phase grid and ordinary moments. The original 12-node derivative integral can miss steep subcell logits. An independent integration-by-parts control on both 8-subdivision reconstructed witnesses compares derivative plus endpoint jumps to the static moment. Maximum absolute change discrepancies on fine fall from 1.77e−12 at 12 nodes to 3.62e−18 at 24. On third the same defect falls from 7.61e−8 at 12 nodes to 1.04e−18 at 64. The production static validation is unchanged; this run uses 32 nodes for the fine audited profiles.

Reproduce from `kinetics/`:

```sh
export OPENBLAS_NUM_THREADS=1
python early_pole_aperture_witness_static_response_reference.py
python early_pole_live_aperture_witness_trajectory_reference.py fine --subdivision 3 \
  --scalar-bins 16 --duration 8e-8 --mix 0 .97 --pair --hard-slow
for w in 0 .97; do
  python early_pole_cold_screened_split_reference.py fine --duration 8e-8 --count 8 \
    --same-soft-strength 1 --compton-strength 1 --conversion-strength 1 --pair \
    --screening-update full --snapshots --aperture-subdivision 3 --scalar-bins 16 \
    --witness-mix "$w" --response-radial-order 32
done
# Repeat the witness split with --duration 2e-8 for the short time control.
python early_pole_live_aperture_witness_trajectory_reference.py fine --subdivision 32 \
  --scalar-bins 32 --mix 0 .85 .97 --refine 2
```

Raw: `results/early-pole-native-witness-joined-relaxation-2026-11-27.json.gz`.
