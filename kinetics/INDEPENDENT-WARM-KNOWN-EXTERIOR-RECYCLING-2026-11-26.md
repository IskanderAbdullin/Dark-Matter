# Admit the known exterior reservoirs to the late recycling network

The earlier [isolated energy envelope](INDEPENDENT-WARM-ISOLATED-RECYCLING-ENVELOPE-2026-11-26.md)
used the retained energy. The [passive old-tail bracket](INDEPENDENT-WARM-EARLY-OMITTED-GRAVITY-2026-11-26.md)
did not test reentry. Here allow **all known old and new exterior
particles to reenter and interact after 10.442 MeV**. A detailed
exterior momentum reconstruction is unnecessary: count every
old and new cohort exactly once, and add their upper energies
to the retained dark energy at the matched 16-block endpoint.
Because the entire enlarged post-anchor reservoir can redistribute,
do not impose a fixed support, charge abundance or negligible
stimulated emission on it.

For the old cohorts, reproduce the old-ledger energy interval
by ballistic dilution between the recorded 12 MeV start and
the anchor. A second, weaker pre-anchor condition assumes
only nonnegative pressure and no energy exchange between old
cohorts and the core in that interval. In that control the
old energy decreases at least as `a⁻³`; its momentum can change.
Use the positive restart archive and its distribution-independent
energy upper for **new** exterior cohorts. Check the endpoint
against both the fully hashed checkpoint and the earlier
instantaneous audit. Retain the larger inclusive number
normalization for the initial inventory, and use the
smaller core number in the denominator of the added-count
ratio to err on the side of overestimating recycling.

| Anchor input | Added energy over core, old ballistic plus new | Added energy over core, old positive pressure only plus new | Largest Ωφh² upper if total is .120 |
| --- | ---: | ---: | ---: |
| q16 | 6.4391×10⁻⁴ | 7.2745×10⁻⁴ | .002810111 |
| q32 | 1.5749×10⁻⁹ | 1.7936×10⁻⁹ | .002809227 |

The 2.35% upper fraction at .120 therefore also holds for this
enlarged known reservoir on the nominal clock. This does not
test earlier transfer between the core and the old omitted
particles, any **unrecorded** inflow or external heating, or
a different reaction network. It is not a new forward FRW
evolution; it is an energy/count bound on the enlarged
starting reservoir.

~~~sh
export OPENBLAS_NUM_THREADS=1
D=results
python -W error warm_known_exterior_recycling_envelope_reference.py EOS.dat \
 "$D"/50MeV-from-empty-pair-quantum-q32-audit-2026-11-26.json \
 "$D"/warm-isolated-recycling-energy-envelope-2026-11-26.json \
 "$D"/warm-early-omitted-gravity-2026-11-26.json \
 --cores "$D"/warm-q{16,32}-10p442-heavy-input-shape-2026-11-26.json \
 --checkpoints "$D"/warm-early-q{16,32}-time16-checkpoint-2026-11-26.json.gz \
 --seeds "$D"/50MeV-from-empty-pair-quantum{,-q32}-shape-2026-11-26.json \
 --producers "$D"/50MeV-from-empty-pair-quantum{,-q32}-2026-11-26.json
~~~
