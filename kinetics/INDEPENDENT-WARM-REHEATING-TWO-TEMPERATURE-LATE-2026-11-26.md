# Switch to two dark temperatures at the controlled heat-clock handoff

Use the [6 MeV handoff](INDEPENDENT-WARM-REHEATING-WIDE-QUANTUM-2026-11-26.md),
where the common-temperature description hits its
declared heat-clock ratio of 30. Compare independent
dark temperatures for the light fermion and the scalar.
Use FD/BE equilibrium moments separately, MB pair
reaction coefficients and energy moments in both
directions, and either retain or turn off the calculated
FD/BE prompt heat exchange. Each species also does
its own pressure work during expansion. This is a
selected two-temperature test, not a Boltzmann solution
for evolving nonthermal spectra.

| Endpoint at visible 2 MeV | Light dark T (MeV) | Scalar dark T (MeV) | Light yield | Scalar yield |
|---|---:|---:|---:|---:|
| Retain prompt heat | .32661 | .88930 | 3.48875e−5 | 6.24403e−4 |
| Omit prompt heat | .31892 | .88914 | 3.50606e−5 | 6.24230e−4 |

Both conserve the joined count for the isolated
light↔scalar conversion and the total energy including
pressure work. The two invariant-mass quadratures agree;
the independent energy-work audit has a relative defect
below 2e−14. The energy coefficient also matches a
temperature derivative of the number coefficient.
The approximation is explicitly limited:
heavy population at the handoff is negligible for these
tests; external visible inflow is omitted. A separate
unblocked lepton source folded on the retained path
would add about 2.2e−9 to the light yield and 1.7e−18
to the scalar-equivalent source count. It is not fed
into these runs. No other visible or dark channels are
included. The endpoint is a chosen 2 MeV probe,
not kinetic or chemical freeze-out.

## Two independent conditional endpoint checks

If those particular counts were frozen at 2 MeV,
the rest-density proxies would be 192.53 for the
light species and .014817 for the scalar in Ωh².
They cannot jointly be used as the fitted background
with a total .12. The scalar's ballistic shape
would travel on average 5.80 Mpc to equality in
a separate fixed reference background if it had
no further scattering. This is a normalized shape
diagnostic, not a background or transfer prediction
for the overabundant two-species mixture.

Independently, allow further forward χχ̄ loss on
ballistically cooling finite cells, with no inverse,
new particles, or redistribution. Use the same
positive shell-integrated [survival method](warm_pair_shell_survival_bound_reference.py),
a 35 MeV incoming momentum cutoff, and a photon-only
expansion floor. The assumption on scalar stimulation
is a *uniform supremum on all later occupations*.
An instantaneous thermal value is not such a bound.
With an assumed cap .5, at least .8931 of the retained
light population survives; its conditional Ωh² floor
is 171.95. A cap of 1 gives 157.67; cap 4 gives 58.84.
Cap 10 gives only a trivial floor. None of these
bounds applies when the stated loss-only and no-
redistribution assumptions fail.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error warm_prompt_reheating_two_temperature_late_reference.py EOS.late.dat cold.json > twoT.json
python -W error warm_prompt_reheating_late_ballistic_floor_reference.py EOS.late.dat twoT.json
python -W error warm_prompt_reheating_late_stream_diagnostic_reference.py EOS.late.dat twoT.json
~~~

[Two-temperature solutions](results/warm-prompt-reheating-two-temperature-late-2026-11-26.json.gz);
[conditional shellwise survival](results/warm-prompt-reheating-late-ballistic-floor-2026-11-26.json);
[conditional free-streaming and rest-density probes](results/warm-prompt-reheating-late-stream-diagnostic-2026-11-26.json).
