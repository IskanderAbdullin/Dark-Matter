# Compare shifted spectra on the old fixed metric

This extends the [forward count reevaluation](INDEPENDENT-FINITE-COOLING-SHIFTED-COUNT-2026-11-26.md).
Here both prompt spectra and the actual early forward recoil bin shapes
are regenerated at the small displaced roots and exported with the
same positive number/second-moment reconstruction. There is no
recalculation of the gravitational driver. At z=20 the two comparisons
verify identical perturbation file hashes and conformal times, with the
same width 0.14 and exact modes k=4,10,22 h/Mpc.

| Shifted minus original scalar δ/δCDM | k=4 | k=10 | k=22 |
| --- | ---: | ---: | ---: |
| Instantaneous | +1.748×10⁻⁵ | +1.093×10⁻⁵ | +4.456×10⁻⁶ |
| Finite, Tmax/TR=5 | −2.982×10⁻⁵ | −4.812×10⁻⁵ | −5.533×10⁻⁵ |

This modest conditional sensitivity is much smaller than the preceding
[width stress test](INDEPENDENT-FINITE-COOLING-MATCHED-HISTORY-LOS-2026-11-26.md)
at these probes. It is a comparison of frozen scalar shapes,
not a claim about CDM, total matter, a continuous cutoff, or the
full perturbations of a coupled model.

## Independent local signed rate at 1 MeV

Reinsert the shifted prompt-plus-recoil scalar and stable prompt shape
in both reaction directions, at width 0.10 and 0.14, with 47 and
61 momentum panels. The inverse-to-quantum-forward ratios at
1 MeV are about 1.47–3.26×10⁻¹¹ in the instantaneous case and
3.96–9.81×10⁻¹² in the finite case. Refinement at fixed width gives
sub-percent changes in this diagnostic; the choice of width changes
the far tail enough to move the very small ratios by factors of order two.
The projected quantum and vacuum forward rates agree to about 9×10⁻¹⁵
in relative units at this anchor. The exact cubic maximum in the
reconstruction also gives a larger, shape-conditional stimulation
ceiling of about 3.4×10⁻⁶ in the finite case, without depending
on resolving the rate difference near roundoff.

The signed check is an **anchor-only** control. It does not bound the
earlier reaction history; early inverse can be significant locally.
Neither calculation evolves post-anchor recoil, a non-frozen stable
distribution, or a self-sourcing metric.

## Reproduction

Using the same released EOS and native CDM driver, with the
existing archive of shifted prompt counts and shapes:

~~~sh
export OPENBLAS_NUM_THREADS=1
D=shifted_spectrum
mkdir -p "$D"
REF=results/finite-cooling-shifted-count-input-2026-11-26.json
DIST=results/finite-cooling-shifted-dist-2026-11-26.json
python -W error finite_cooling_during_production_pair_reference.py EOS.dat "$REF" --shape >"$D/shape.json"
python -W error export_finite_cooling_class_psd_reference.py EOS.dat "$DIST" "$D/prompt" >"$D/prompt.json"
python -W error export_finite_cooling_pair_mixture_reference.py EOS.dat "$DIST" "$D/shape.json" "$D/mix" >"$D/mix.json"
python -W error finite_cooling_matched_history_los_reference.py recoil EOS.dat "$D/shape.json" "$D/recoil" >"$D/recoil.json"
for hist in instant ratio5; do
  python -W error finite_cooling_relativistic_los_reference.py EOS.dat \
    "$D/mix/$hist-scalar-pair-mix-sigma0.14.dat" "$D/prompt/$hist-scalar-sigma0.14.dat" \
    "$D/recoil/$hist-scalar-recoil-sigma0.14.dat" \
    cdm_native/cdm_los_00_perturbations_k{2,4,6}_s.dat --cdm-only --z 20 >"$D/$hist.json"
done
python -W error finite_cooling_shifted_spectrum_audit_reference.py "$D/"{prompt,mix,recoil}.json \
  --old results/{instant-reheat-cdm-los,finite-cooling-cdm-los}-z20-2026-11-26.json --new "$D/"{instant,ratio5}.json
python -W error finite_cooling_shifted_anchor_rate_reference.py EOS.dat "$D"
~~~
