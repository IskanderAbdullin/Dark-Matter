# An initial negative dressed source can persist in the smaller live operator

Initialize the [fine 97% dressed witness](INDEPENDENT-EARLY-POLE-APERTURE-WITNESS-TWO-MESH-2026-11-27.md), keeping 16 scalar bins and 16 hard subdivisions. The parent shape is projected into strictly admissible cell moments before starting the solve. The initial source matches the frozen fixture; the hard and scalar occupations, variational gauge energy and scalar Hartree gap now evolve. The reference exposure is 2e−6, with no elastic scattering, Compton or flavor pair conversion in this test.

- Cold decay/inverse with the archived step input adds 1.5845543e−7 scalar number; the 97% witness removes 2.0715291e−8. Its sampled instantaneous source starts at −.0103665 and ends at −.0103488, remaining negative at all 18 sampled times. Tightening tolerances and halving the maximum time step changes the timeline by less than 1e−15.
- Enabling the finite photon-pair reaction in the witness branch yields −2.0765187e−8, with final rate −.0103987. Photon-pair generation alone does not remove the opposite-sign behavior on this interval.
- The coupled energy drift is below 5e−16, with separately integrated entropy mismatch below 7e−12. Without photon pairs, the two opposite changes in fermion number mirror the scalar count; with photon pairs, the summed live hard count is conserved to 5e−13.

This is a finite early-time result for an explicitly prepared input, not a population estimate or a full collision hierarchy. It does not resolve what happens with screened elastic smoothing or later-time expansion.

Reproduce from `kinetics/`, with `OPENBLAS_NUM_THREADS=1`:

```sh
python early_pole_live_aperture_witness_trajectory_reference.py fine --mix 0 .97
python early_pole_live_aperture_witness_trajectory_reference.py fine --mix .97 --refine 2
python early_pole_live_aperture_witness_trajectory_reference.py fine --mix .97 --pair
```

Raw: `results/early-pole-native-witness-cold-trajectory-2026-11-27.json.gz`.
