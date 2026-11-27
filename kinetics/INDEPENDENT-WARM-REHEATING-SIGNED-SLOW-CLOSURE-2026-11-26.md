# A signed dilute slow closure for the hotter reheating anchor

The [fast return experiment](INDEPENDENT-WARM-REHEATING-FAST-PROMPT-RETURN-2026-11-26.md)
supports testing a closure that is distinct from
keeping both direct production spectra fixed. At
each initial anchor impose a common dark temperature
and enforce the fast decay affinity
eta_heavy=eta_light+eta_scalar. Match the **total**
dark energy and the two counts invariant under each
individual fast decay: L=N_light+N_heavy and
S=N_scalar+N_heavy. This allows a small, nonzero
heavy population rather than setting it to exactly
zero while also imposing inverse-decay equilibrium.

Consider three selected vacuum pair processes in
both directions: light-light ↔ heavy-heavy (opposite
charge and both same-charge channels), scalar-scalar ↔
heavy-antheavy, and light-antilight ↔ scalar-scalar.
The exact spin and identical-particle factors enter
the positive invariant mass integral; for each
common-temperature kernel the MB event rates just
multiply the respective squared fugacities. If
G_LH, G_SH and G_LS are signed net **event** rates,
then dL/dt=2(G_SH−G_LS), dS/dt=2(G_LH+G_LS).
A large positive gross heavy source therefore does
not determine even the sign of dS/dt.

Independently match the FD/BE thermal moments in
a convergent fugacity series at these dilute anchors.
For the resulting thermal state detailed balance
gives a constant reverse/forward ratio for each
collision, so use global occupation factors on
the positive forward rate and then apply its exact
signed affinity. The FD factors and the Bose final
factors have different bounds; neither was
substituted for the other. The brackets below
include these quantum factors on the sampled
fixed-T closure, not errors due to incomplete
physics or nonthermal distributions.

| Reheating → visible bath (MeV) | Net change of scalar-equivalent invariant S per Hubble time / visible entropy, thermal bracket |
|---:|---:|
| 50 → 35 | [−2.191e−5, −2.188e−5] |
| 100 → 85 | [+2.063e−4, +2.206e−4] |
| 100 → 75 | [−9.919e−5, −5.232e−5] |
| 100 → 60 | [−.002146, −.002039] |
| 100 → 50 | [−.004401, −.004268] |

At 100→60 MeV the heavy population is about .195% of
the three-species count. The reverse of heavy
regeneration is small on this closure, but
scalar-scalar → light-antilight more than cancels
the direct scalar-equivalent gain from light
regeneration. The forward light pair cross section
also matches an independent thermal lab-frame fold
within 3.8e−7 on the finer control. The heavy gross
forward total agrees within 0.7% with the independently
folded dilute quantum shape at the 100→60 MeV point.

A separate fixed-volume MB integration advances the two
invariants and re-solves temperature at fixed dark
energy. From this 100→60 MeV anchor, S first falls
from 5.8701e−5 to 5.6531e−5 by .003 Hubble times, then
rises to 7.5420e−5 at .02 Hubble times. Two quadratures
and time tolerances agree on this path. The quoted
time trajectory uses MB, while the initial sign was
also tested with the separate quantum bracket.
No continuing bath injection, expansion, elastic
feedback, visible reheating dynamics or other dark
channels are included in that short trajectory.
The calculation is not a cosmic yield or a final
mass-power prediction.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
D=results
python -W error warm_prompt_reheating_heavy_gain_probe_reference.py EOS.late.dat local-spectra > passive.json
python -W error warm_prompt_reheating_signed_slow_closure_reference.py EOS.late.dat local-spectra \
 "$D/warm-prompt-reheating-matched-entropy-gain-2026-11-26.json"
~~~

[All signed rates, closure checks and short history](results/warm-prompt-reheating-signed-slow-closure-2026-11-26.json).
