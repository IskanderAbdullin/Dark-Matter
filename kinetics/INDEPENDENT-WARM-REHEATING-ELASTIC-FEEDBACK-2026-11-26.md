# Feed the measured off-shell heating back into the two-temperature diagnostic

Retain the same 6 MeV handoff, signed MB pair rate,
quantum thermal moments and pressure work used in the
[isolated comparison](INDEPENDENT-WARM-REHEATING-TWO-TEMPERATURE-LATE-2026-11-26.md).
Add the [selected invariant elastic kernel](INDEPENDENT-WARM-REHEATING-ELASTIC-HEAT-AUDIT-2026-11-26.md)
strictly below s=.003 GeV², without reusing the on-shell
pole. Tabulate its positive thermal exchange coefficient
on the current path, update the histories of both temperatures
and their particle numbers, then recompute the coefficient
on the resulting path. The energy exchanged by the
elastic channel has opposite signs in the two species.
It changes no particle count directly.

| State at visible 2 MeV | Light dark T (MeV) | Scalar dark T (MeV) | Light yield | Scalar yield |
|---|---:|---:|---:|---:|
| Isolated, prompt heat only | .32661 | .88930 | 3.48875e−5 | 6.24403e−4 |
| With below-pole elastic exchange | .87562 | .88219 | 3.19786e−5 | 6.27312e−4 |

The extra exchange nearly removes the temperature split.
The light count is about 8.34% lower in this comparison,
not enough by itself to remove the stable overabundance
in the hypothetical frozen-count proxy. The 11- and
15-state calibrations change the terminal counts by
less than 6e−9 relative. Independently integrating
pressure work accounts for the total energy change
with a relative check of about 1.6e−11.

The direct held-state scattering also changes
particle **momenta** substantially. On the old
isolated 2 MeV endpoint, mean squared momentum kicks
per Hubble time, divided by the total species momentum
second moment, are about 340 for light particles and
42 for scalars. These positive moments are not
eigenmode damping rates. They do show why the
[conditional loss-only floor](INDEPENDENT-WARM-REHEATING-TWO-TEMPERATURE-LATE-2026-11-26.md#two-independent-conditional-endpoint-checks)
cannot be carried over when this redistribution remains
active. Nothing here proves the separate thermal shape
closure or a final relic abundance.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error warm_prompt_reheating_two_temperature_elastic_feedback_reference.py EOS.late.dat twoT.json
python -W error warm_prompt_reheating_late_kick_reference.py EOS.late.dat twoT.json
~~~

[Iterated thermal trajectories](results/warm-prompt-reheating-elastic-feedback-2026-11-26.json.gz);
[separate held-state momentum kick moments](results/warm-prompt-reheating-elastic-kicks-2026-11-26.json).

The [later shape-independent count inequality](INDEPENDENT-WARM-REHEATING-POST-2MEV-ELASTIC-2026-11-26.md)
does not require an assumption of ballistic
fermion momenta, but keeps an explicit bound on
every later scalar occupation and on the participating
number-changing interactions.
