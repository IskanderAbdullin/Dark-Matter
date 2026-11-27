# Late competing reactions for an already above-cut .80 GeV packet

This independently varies the reactivated hard-packet shape relative to
the [.68 GeV turning control](INDEPENDENT-WARM-LATE-COMPETING-FLUX-2026-11-26.md).
Use the same known exterior birth moments, the same cut 0.7148213914435242 GeV
and the exact no-floor .80 GeV prompt ladder. In this case the W observable,
the sum of positive energy excess per daughter above the cut, starts nonzero
and the signed prompt current is negative at all saved points.

At Hdt = 0, 1e-8, 1e-7 and 1e-6, the wide-lab-support below-pole
elastic and open-outgoing signed pair currents are negative under both
quadratures. Their summed gross W currents together are below 1.74e-6
of the absolute prompt current at the selected phases. The 32T lab
cap again does not intersect the energy range of the observable.
These ratios are comparisons on known phases, not a dynamical envelope.

From Hdt = 1e-7, evolve a matched interval Δlog a = 9e-7. With one
implicit collision block the physical W falls from roughly 6.622e-14
to 2.691e-14 of initial represented dark energy. Compared to matched
expanding prompt, the endpoint changes by −2.52e-6 after adding the
pair channel and −2.59e-6 after also adding the closed elastic channel.
Keep these as local numerical comparisons rather than a prediction of
the .80 GeV tail lifetime. Omitted channels and the unknown actual
exterior shape are not constrained by this test.

Halving the expansion block while keeping the exact same initial
phase and physical cut moves the endpoint by 17.40% in both the
prompt and fully coupled controls. Both still decline. The full
two-block endpoint is 2.223011614e-14 of initial dark energy,
2.63e-6 below its matched prompt control. The one- and two-block
differences thus have the same sign and order of magnitude, but
the overall endpoint is far from converged on this coarse clock.
The separate six-microstep non-expanding run is not the reference
solution for this matched-clock comparison.

The exact collisionless Liouville endpoint for this block would instead be
6.62183188e-14 (the initial value is 6.62190706e-14). The dominant decrease
comes from the reactions in this prescribed block. The initial pair
current is (−1.94 to −2.10)e-13 of initial dark energy per Hubble time,
and the initial elastic current is (−2.26 to −2.66)e-14. Independent
Liouville and reconstructed first-redshift W agree to 2.0e-15 relative.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MALLOC_ARENA_MAX=2
a=results/warm-corrected-anchor-signed-floor-2026-11-26.json
b=results/warm-multishell-corrected-q64-fine-stage-48-2026-11-26.json.gz
s=results/warm-late-aligned-hard080-prefix-2026-11-26.json.gz
mkdir -p competing080
python -W error warm_corrected_late_prompt_cut_aligned_compact_after_prefix_reference.py EOS.dat "$a" "$b" late080 --hard .80 --substeps 6 --predictor --prefix "$s"
python -W error warm_late_competing_legacy_tape_adapter_reference.py "$a" late080 competing080/tapes
d=competing080
t="$d/tapes"
for k in elastic pair; do
  python -W error warm_late_reactivated_competing_flux_reference.py EOS.dat "$d/$k-low.json" "$k" "$t"/step-{0,4,5,6}.json
  python -W error warm_late_reactivated_competing_flux_reference.py EOS.dat "$d/$k-high.json" "$k" "$t"/step-{0,4,5,6}.json --order-mode 1
done
run() { python -W error warm_late_reactivated_expanding_competing_reference.py EOS.dat "$@"; }
run "$t/step-5.json" "$t/step-6.json" "$d/expanding-prompt.json" prompt
run "$t/step-5.json" "$t/step-6.json" "$d/expanding-pair.json" prompt-pair
run "$t/step-5.json" "$t/step-6.json" "$d/expanding-all.json" all --extent 285.
python -W error warm_late_competing_flux_audit_reference.py "$t/step-5.json" "$d"
python -W error warm_late_competing_intermediate_guide_reference.py EOS.dat "$t/step-5.json" "$t/step-6.json" "$t/mid-5-to-6.json"
for mode in prompt all; do
  cp "$d/expanding-$mode.json" "$d/compare-$mode-coarse.json"
  run "$t/step-5.json" "$t/mid-5-to-6.json" "$d/compare-$mode-mid.json" "$mode" --extent 285.
  run "$t/mid-5-to-6.json" "$t/step-6.json" "$d/compare-$mode-fine.json" "$mode" --extent 285. --predecessor "$d/compare-$mode-mid.json"
done
python -W error warm_late_competing_clock_audit_reference.py "$t/step-5.json" "$d" --max-relative-delta 5e-6
~~~

[Two-order channel and one-block audit](results/warm-late-hard080-competing-flux-audit-2026-11-26.json),
[clock audit](results/warm-late-hard080-competing-clock-audit-2026-11-26.json),
[frozen and one-block data](results/warm-late-hard080-competing-flux-tapes-2026-11-26.json.gz),
[clock and intermediate guide data](results/warm-late-hard080-competing-clock-tapes-2026-11-26.json.gz).
