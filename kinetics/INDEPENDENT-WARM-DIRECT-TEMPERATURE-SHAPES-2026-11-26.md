# Higher direct reheating changes both the count and the momentum spectrum

Repeat the central production calculation from the
[reheating ladder](INDEPENDENT-WARM-DIRECT-TEMPERATURE-BRACKET-2026-11-26.md)
retaining the parent energy. It includes only vacuum FD charged
leptons and the BE charged-pion measured window with central
bins. Propagate each prompt daughter through the exact uniform
decay-angle interval in lab energy, then map its physical
momentum to conserved q using the visible entropy clock.
Use the entire represented 6 MeV to reheating interval at
each temperature and keep the two production sources separate.

| Reheating (MeV) | Scalar omega from both included sources | Pion share of this count | Exact RMS q/T_gamma,0 | Shell-bounded mean path to equality (Mpc) |
|---:|---:|---:|---:|---:|
| 25 | 2.7007e−6 | 0.56% | 2.480 | [10.983, 11.912] |
| 50 | 1.58975e−4 | 8.24% | 1.931 | [8.453, 9.180] |
| 75 | 9.1986e−4 | 27.18% | 1.831 | [7.881, 8.561] |
| 100 | 2.91670e−3 | 41.89% | 1.708 | [7.404, 8.046] |

The independently sampled source-by-source CDFs change
by less than 3e−5 on shell boundaries when the angular,
invariant, lab-energy and time rules are refined. The
50-to-100 MeV change is about 0.105 in the combined CDF.
At 100 MeV the absolute dose is 18.35 times the 50 MeV
result, while RMS q/T_gamma,0 is about 11.5% smaller.
Both included source spectra move to smaller RMS q
as reheating rises; their different weights also change.
This refers to the redshifted spectrum, not to the local
energy of an individual daughter. The path comparison
starts all spectra at the same 6 MeV anchor. Monotonic
velocity bounds every shell's mean without assuming a
within-shell profile. The high-q tail is assigned the
causal matter-only upper path for the upper calculation.
The 50 and 100 MeV mixture path intervals do not overlap.

The data also prepare matched native inputs: one power
run keeps each absolute scalar count and another holds
the 50 MeV scalar count while varying only the shape.
Each has a negligible same-shape trace and fills the
rest of omega_dark=.12 with CDM. Neither this setup nor
the path result solves the stable production, nonlinear
recycling or a thermal hadronic spectrum.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
D=results
python -W error warm_direct_vacuum_temperature_psd_reference.py \
  EOS.late.dat KLOExSec.data "$D/warm-direct-vacuum-temperature-bracket-2026-11-26.json" \
  higherk/compact-trace.ini varied > varied-spec.json
python -W error warm_direct_vacuum_temperature_path_audit_reference.py \
  EOS.late.dat varied-spec.json
~~~

[Positive spectral counts, exact moments and PSD checksums](results/warm-direct-vacuum-temperature-psd-2026-11-26.json);
[ballistic path bounds by source](results/warm-direct-vacuum-temperature-path-audit-2026-11-26.json).

The [two-tier native power comparison](INDEPENDENT-WARM-DIRECT-TEMPERATURE-NATIVE-2026-11-26.md)
reports both density conventions and documents the remaining q-grid sensitivity.
