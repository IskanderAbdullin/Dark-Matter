# Count outward and inward prompt transitions separately

The [short explicit-exterior transient](INDEPENDENT-WARM-PROMPT-REACTING-EXTERIOR-TRANSIENT-2026-11-26.md)
left very little material outside the original boundary despite
a large frozen first-crossing rate. Measure both directions
directly on the same coupled prompt trajectory.

At fixed virtual parent energy, write (A) for the inverse
weight of an occupied daughter pair and (B) for its decay
weight. Resolve inside and outside separately for each tagged
species. Elimination of the stationary parent gives effective
outward transitions proportional to (A_{m in} B_{m out}/
(A_{m all}+B_{m all})), and inward transitions with
inside and outside exchanged. A decay outside that follows an
inverse event outside for the same tagged species is not a
crossing. The difference agrees with the native projected
number rate without subtracting two separate trajectories.

Integrate the positive fluxes using each accepted backward
Euler step's actual endpoint. At the fine (HDelta t=5	imes10^{-6})
endpoint for initially empty q16 exterior, outward and inward
counts per initial scalar are .03859179 and .03857129. For
the light fermion they are .03040444 and .03036742.
The differences are (2.04921	imes10^{-5}) and
(3.70202	imes10^{-5}), matching the independently
tracked outside count to (3	imes10^{-17}) at that endpoint.
The two finest time subdivisions shift gross totals by less
than (8	imes10^{-6}) relative. Large positive crossing
counts do not imply persistent loss or a large net outward flow.

Starting with the actual q32 archive and measuring beyond q32
instead, the positive fluxes integrated over the same short
control are about (6.80	imes10^{-9}) and (1.10	imes10^{-8})
outward; the net resident fractions are only
(1.15	imes10^{-11}) and (2.22	imes10^{-11}).
This is still an isolated prompt model, not tagged particle
trajectories, full FRW evolution or permanent exit.

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error warm_prompt_gross_crossing_reference.py
~~~

[Diagnostics](results/warm-prompt-gross-crossing-2026-11-26.json.gz).
