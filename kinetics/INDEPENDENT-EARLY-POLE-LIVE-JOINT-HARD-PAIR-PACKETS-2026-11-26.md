# Finite hard vacuum packets in the joint variational moment test

Replace the two incident samples in the
[joint channel witness](INDEPENDENT-EARLY-POLE-LIVE-JOINT-PAIR-DECAY-MOMENTS-2026-11-26.md)
with incident radial/relative-angle quadrature and the vacuum
Breit–Wheeler differential cross section. The selected hard support is
.5≤q_gamma/T≤2, s/T²≥1, |cos(theta*)|≤.8 and both fermion
momenta on the active hard support (.45–2.2). Neither the exterior
domains nor soft channels are represented. The integrated unblocked
angle weights before the outgoing cutoff agree with the independent
analytic central cross section. Identical photon events carry one
half the double-counted incident measure.

For each retained vacuum packet, keep the two incident momenta and
the first outgoing magnitude. Remap the second outgoing magnitude
so that the instantaneous variational *energy* closes, then verify
the momentum triangle on the remapped packet. Evaluate all quantum
occupation factors on the same on-shell packet. Weights remain the
vacuum coefficients; no in-medium amplitudes, gauge residues,
physical time, or hard/soft matching are inferred. The largest
checked remapping rejection in the runs is about 1.6e−6 of the
precomputed retained weight. The prior outgoing hard-support cut
is much larger: approximately 30% of the unblocked central weights
are outside that support.

A plain azimuth quadrature can sample nearly discontinuous cell
occupation or support changes. A second implementation solves the
azimuth angle at each crossed outgoing boundary and integrates each
side separately. This distinction matters for the weak coupled effect:

| Fixed-scale variant | Relative scalar gain change at pair multiplier 1 | At 4 |
| --- | ---: | ---: |
| Two old sampled angles | −13.7% | not compared |
| Full hard packets, unsplit coarsest | +.38% | −2.22% |
| Full hard packets, unsplit finest | −1.18% | −5.80% |
| Azimuth split, three refinements | −1.17% to −1.24% | −5.75% to −5.95% |

All rows compare to their own matched decay-only run and the
reported pair clocks are imposed dimensionless strengths.
The different support/rate normalization in the old angle example
makes its numerical size not directly comparable. The sign on the
coarse unsplit new grid is not a reliable conclusion. In the
refined hard test the joint occupation effect is negative, but
neither that sign nor its size is established for the cosmological
bath. The exact count, energy, and entropy ledgers remain valid.

Reproduce with early_pole_live_joint_hard_pair_packets_reference.py
(--level 0–2, optional --split-azimuth).
Raw output: results/early-pole-live-joint-hard-pair-packets-2026-11-26.json.gz.
