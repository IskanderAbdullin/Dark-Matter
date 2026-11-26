# Fixed-state on-shell exchange across the old momentum boundary

The retained q/T_R=16 boundary suppresses all on-shell participation
after b=3.63867. Before then, not retaining an exterior daughter can
also conceal a large instantaneous bidirectional flux. We compute
that flux without evolving or thermally reconstructing either species.

One probe starts with the archived 12 MeV q/T_R=16 quantum-source
snapshot and vacuum outside it. Another starts with an independent
24-step q/T_R=32 vacuum-source archive, using the old q/T_R=16 edge
as an internal cut. At later b the snapshot is **collisionlessly
redshifted for the probe only**; it is not evolved with that outward or
inward rate. At every parent energy keep the full daughter angular
support, evaluate the archived BE/FD shape on each leg and use the
stationary parent occupation A/(A+B). Occupations outside the input
mesh remain vacuum. Angular integrations split at each shell and at
the tested boundary. Number and kinetic energy in the exterior are
reported separately for gains and losses.

The q16 archive with initially vacuum exterior gives an instantaneous
outgoing scalar count of about 5.67e3 times N_phi H at b=1. No inward
exterior count exists in that initial condition. The populated q32
archive has outward/inward rates of about 8002 and 8180 times N_phi H
through the same q16 edge; its instantaneous net is inward and much
smaller than either gross rate. These are rates averaged over the
retained scalar count. They must not be integrated for a Hubble time:
populating or draining the thin exterior can quickly change both rates.
The two archives are also not identical initial conditions. The
fixed-state check therefore diagnoses sensitivity and the need for
explicit exterior return, not a late-time abundance correction.

At b=5 in the redshift-only probe the q32 scalar gross rates are both
about 6.2e-5 times N_phi H. At b=8 the q32 input cannot form an on-shell
parent with two retained daughters, so the fixed-state source is zero.
An even higher unknown tail is not bounded by these statements.

## Reproduce

Prepare both [EOS inputs](INDEPENDENT-BENCHMARK-SETUP.md). The q32
source test uses a coarser time mesh and redistributes inner mesh edges;
it should not be read as the sole estimate of cutoff sensitivity.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1
python -W error joint_resonant_frw_split_reference.py EOS.dat --physical-amplitude --visible-source --matched-vacuum-order --visible-backreaction --trace-exterior --with-elastic --empty-pulse 1e-8 --empty-stop .24 --support-low .005 --empty-support-high 32 --empty-cells 22 --steps 24 --vacuum-scaled-residual --allow-empty-shells --include-zero-cell --save-snapshots wide.json
python -W error late_prompt_exterior_snapshot_reference.py EOS.dat EOS.late.dat results/wide32-early-vacuum-2026-11-26.json results/wide32-early-vacuum-shape-2026-11-26.json --boundary-q-over-TR 16
```
