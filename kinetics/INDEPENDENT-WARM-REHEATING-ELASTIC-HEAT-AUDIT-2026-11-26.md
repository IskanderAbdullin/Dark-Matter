# Off-shell elastic heat is not negligible on the isolated two-temperature path

The [two-temperature 6→2 MeV comparison](INDEPENDENT-WARM-REHEATING-TWO-TEMPERATURE-LATE-2026-11-26.md)
omits elastic χφ→χφ scattering. Fold the selected
tree-level s/u amplitude on three held states of
that path. The scattering kinematics use the
initial and final directions in the center of mass
and boost the final energies back into the plasma
frame. Include the unpolarized stable occupation,
Pauli blocking and the scalar Bose factor in both
directions at the two different dark temperatures.
Use only s below .003 GeV², safely below the .01 GeV²
parent pole, so this window does not re-count the
prompt on-shell return already used in that path.

| Visible T (MeV) | Isotropic threshold control, Q_scalar/(Hρ_scalar) | Invariant s/u amplitude, Q_scalar/(Hρ_scalar) |
|---:|---:|---:|
| 4.558 | −4.730 | −2.696 |
| 3.463 | −4.998 | −3.128 |
| 2.000 | −3.000 | −2.190 |

Negative values transfer energy from the hotter
scalar population to the stable population.
The exact invariant result is smaller than the
threshold contact control but still comparable
to the scalar energy over a Hubble interval.
A same-temperature quantum fixture cancels
below 2e−18 of the gross event density times
temperature; the independent scalar-rest
trace agrees at sampled angles at 1e−15 relative.
Changing from 10 to 19 azimuth points has a
negligible effect on the held-state rate.
Lowering the cut from .003 to .0024 GeV²
changes the 4.558 MeV result by about .013%.
This does not evaluate the omitted higher-invariant
interval or any nonthermal medium correction.

The original endpoint is a valid **isolated**
diagnostic, but its temperature split and the
remaining stable-pair reactions cannot be
transferred to the coupled theory without
allowing this heat exchange to modify the
distributions and their densities.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error warm_prompt_reheating_two_temperature_contact_probe_reference.py EOS.late.dat twoT.json
python -W error warm_prompt_reheating_two_temperature_contact_probe_reference.py EOS.late.dat twoT.json --exact
python -W error warm_prompt_reheating_elastic_frame_audit_reference.py twoT.json
~~~

[Contact, invariant-amplitude and independent checks](results/warm-prompt-reheating-elastic-held-audit-2026-11-26.json).

An [iterated feedback test](INDEPENDENT-WARM-REHEATING-ELASTIC-FEEDBACK-2026-11-26.md)
then applies the below-cut heat to both temperature ledgers.
