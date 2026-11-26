# Fixed-state pair-rate separation at the broadened warm checkpoint

Evaluate the signed tree pair operator at the archived 12 MeV
surface of the [q32 production check](INDEPENDENT-WARM-SIGNED-Q32-PRODUCTION-2026-11-26.md).
Run three cases: the original q16 simulation state, the q32
simulation state, and the exact inner crop of the latter back
to q16. The last two use the same inner occupancies, so their
difference isolates the instantaneous effect of including
the resolved outer bins at the fixed state. It does not
separate the history of producing the two simulation shapes.

Normalize each rate to radiation H and the retained stable
number in that specific state. Fine angular controls:

| Archived state | Gross stable loss | Gross stable gain | Signed stable loss | Inverse / forward |
| --- | ---: | ---: | ---: | ---: |
| Original q16 | 1.382282 | 1.167399 | 0.214883 | 0.844545 |
| q32 | 1.382421 | 1.166891 | 0.215530 | 0.844092 |
| q32 cropped without changing inner occupancies | 1.382407 | 1.166833 | 0.215574 | 0.844059 |

Relative to its own crop, including the outer bins changes
gross loss by +1.07e-5 and gross gain by +4.99e-5.
The signed difference changes by −2.01e-4; its magnitude is
comparable to the coarse-to-fine angular change of 1.68e-4.
In contrast, the cropped q32 state and the separately evolved
original state differ in signed loss by 0.322%, despite sharing
the same grid. That includes changes to the archived inner
occupancies from the production history. The result preserves
the signed reverse term, which is about 84% of forward here.

The q32 grid represents nearly all the instantaneous scalar
final-state number in this test; the q16 crop represents about
99.9965% of it. This is a fixed-state pair test only. It does
not establish convergence of a resumed evolution, alter the
rest-frame tree kernel, or add omitted reactions.

At the same checkpoint the homogeneous expansion envelope
with the actual broader archived shape is about 5.04682e-5
in fractional H relative to dark-free radiation at the same
temperature. The original q16 estimate was about 5.04690e-5.
These are instantaneous calculations, not expansion likelihoods.

~~~sh
export OPENBLAS_NUM_THREADS=1
r=results/50MeV-from-empty-pair-quantum
python -W error warm_early_outer_pair_rates_reference.py EOS.dat \
    "$r-shape-2026-11-26.json" "$r-q32-shape-2026-11-26.json"
python -W error warm_archived_expansion_exterior_envelope_reference.py EOS.dat \
    "$r-q32-2026-11-26.json" "$r-q32-audit-2026-11-26.json" "$r-q32-shape-2026-11-26.json"
~~~
