# Reevaluate the forward loss at the displaced low-reheat roots

The [fixed-loss retarget](INDEPENDENT-FINITE-COOLING-FROZEN-COUNT-RETARGET-2026-11-26.md)
increased TR by about 1.11 and 1.31 keV. This run regenerates
the prompt count, spectrum, and early forward rate at those
displaced temperatures, rather than assuming the loss stays
constant. It uses the same radiation/cooling clock and the
same two matched early reaction meshes, and treats the late
window separately from 1 MeV down to 30 keV.

| History | Increase in early forward loss | Stable Ωh² after the reevaluated loss | Additional local TR update |
| --- | ---: | ---: | ---: |
| Instantaneous | +0.13601% | 0.1199999044 | +1.50 eV |
| Finite, Tmax/TR=5 | +0.21118% | 0.1199997498 | +2.75 eV |

The absolute fractional changes between the two matched early
quadratures agree to 2.1×10⁻⁷ or better. The regenerated 1 MeV
prompt stable count agrees with the independent source-fold
normalization to about 2×10⁻¹⁰ in Ωh² units. The late forward
count is recomputed using the shifted frozen moment distribution
and both late quadrature choices. For the final number audit its
normalization alone is reduced by the early number loss (a
quadratic factor in the late forward count at fixed shape).
This is not a reconstruction of the momentum-dependent depletion.
The local update uses the prompt count slope checked in the
previous calculation; no third exact root or physical precision
at eV level is inferred.

The positive prompt calculation, spontaneous vacuum forward
conversion and prescribed radiation clock are still distinct
approximations. No joint inverse source, stimulation, coupled
depletion of the distribution, or measurement is included.
The frozen pair+prompt scalar RMS at 1 MeV changed from
7.37380 to 7.37345 (instantaneous), and 2.04729 to 2.04764
(finite), in units of p/T. That describes only the retained
mathematical proxy, not an updated transfer or viability bound.

## Reproduction

Use the original radiation table and the already archived roots.

~~~sh
export OPENBLAS_NUM_THREADS=1
P=results/finite-cooling-frozen-count-retarget-2026-11-26.json
E=results/finite-cooling-during-production-pair-shape-2026-11-26.json
python -W error finite_cooling_shifted_count_audit_reference.py prepare EOS.dat "$P" >input.json
python -W error finite_cooling_during_production_pair_reference.py EOS.dat input.json >early.json
python -W error finite_cooling_reheat_distribution_reference.py EOS.dat input.json >dist.json
python -W error finite_cooling_post_anchor_pair_diagnostic_reference.py EOS.dat dist.json >late.json
python -W error finite_cooling_shifted_count_audit_reference.py compare "$P" input.json "$E" early.json dist.json late.json
~~~
