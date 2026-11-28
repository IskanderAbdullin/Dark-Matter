# Bound the instantaneous gravity of the old omitted tails

The [matched q16/q32 runs](INDEPENDENT-WARM-EARLY-SUPPORT-TIME-INTERACTION-2026-11-26.md)
include gravity from new exterior cohorts but do not have the
original sampled momenta of the earlier production tails. Their
numbers enter the count audit once. Instead of inventing a
momentum distribution for them, use the archived lower and upper
energy for each old tail group, with its independently normalized
rest count. Assume those older particles are collisionless and
do not return.

For expansion (b), a particle's kinetic energy satisfies
(K/b^2 \le K(p/b)\le K/b).
This gives a late bracket from the initial energy limits, after
including exact rest-mass dilution. Compare adding that bracket
to the otherwise unchanged stage energy density: radiation,
retained dark particles, and the one-mean positive new exterior
cohorts actually used in the evolver. Direct positive quadrature
checks the retained contributions at both resolutions.

| Support | At the 12 MeV start: added ΔH/H | At 10.442 MeV, eight-block end |
| --- | ---: | ---: |
| q16 | [3.147, 3.149]×10⁻⁸ | [2.809, 3.187]×10⁻⁸ |
| q32 | [7.954, 7.955]×10⁻¹⁴ | [7.024, 8.012]×10⁻¹⁴ |

The four-block endpoints give virtually the same instantaneous
brackets. This substantially reduces the old *omitted energy*
when the physical momentum support is extended. It is **not**
a rerun with older gravity, not a bound on its induced abundance
or reaction-rate error, and not evidence that the old particles
could not stimulate reactions if re-entry were allowed.

~~~sh
export OPENBLAS_NUM_THREADS=1
D=results
python -W error warm_early_omitted_gravity_reference.py EOS.dat "$D"/50MeV-from-empty-pair-quantum-q32-audit-2026-11-26.json \
 --production "$D"/50MeV-from-empty-pair-quantum{,-q32}-2026-11-26.json \
 --seeds "$D"/50MeV-from-empty-pair-quantum{,-q32}-shape-2026-11-26.json \
 --q4 "$D"/warm-early-q{16,32}-serial-checkpoint-2026-11-26.json \
 --q8 "$D"/warm-early-q16-time8-checkpoint-2026-11-26.json.gz "$D"/warm-early-q32-time8-checkpoint-2026-11-26.json
~~~
