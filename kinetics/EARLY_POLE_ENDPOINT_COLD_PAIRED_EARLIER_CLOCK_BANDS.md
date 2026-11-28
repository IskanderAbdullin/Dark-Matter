# Both S2 and S4 whole bands respond to the earlier timing change

Two independently evolved S2/S4 endpoints at each earlier clock state each pass the complete same-band panel audit (57 and 71 panels respectively). Only x=.02–.04 changes from 2048 to 1024 steps, while all final .04–.05 intervals retain 2048. Both original and propagated earlier changes are tested on the same 17 display bins.

| Gross response | Whole-band change, S2 | Whole-band change, S4 | S2/S4 conditioned TV, original | S2/S4 conditioned TV, earlier-coarse |
| --- | ---: | ---: | ---: | ---: |
| forward number | -0.458954% | -0.4538384% | 2.022556% | 2.022196% |
| inverse number | -0.4986853% | -0.4931949% | 2.020395% | 2.020133% |
| forward excess | -0.4557899% | -0.4505841% | 1.87111% | 1.874391% |

The explicit scalar occupation shares of the integrated log shift, in table order, are S2: 1.11553%, 1.02575%, 1.09618%; S4: 1.10391%, 1.01306%, 1.0255%. Both replacement orders use exactly matching scalar-cell supports within each level; the remaining kernel is not further resolved.

The maximum paired profile and contrast discrepancy is 1.66e-13, and the two-order factor check discrepancy is 3.02e-13. The S2/S4 profile minima use independent signed-support bounds and positive witnesses. This is a discrete same-geometry earlier-time perturbation and not a global error estimate, full production model, or physical convergence result.

From `kinetics/`:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28; p=results/fine-charge2-scalar; t=results/early-pole-cold; s=results/early-pole-fine-charge2-scalar
a=results/finite-mass-qed-fixed-pole-cold-scalar-band-$d.json
args=(--s2-alt "$p"2-prior-clock-endpoint-continuum-cold-{1,2}.jsonl.gz --s4-alt "$p"4-prior-clock-endpoint-continuum-cold-{1,2}.jsonl.gz
      --s2-original "$p"2-endpoint-continuum-cold-{1,2}.jsonl.gz --s4-original "$p"4-endpoint-continuum-cold-{1,2}.jsonl.gz
      --states "$t"-s{2,4}-prior-clock-x05-state-$d.json.gz "$s"{2,4}-x05-state-$d.json.gz --archive "$a")
python -W error early_pole_endpoint_cold_prior_clock_band_comparison.py "${args[@]}" \
  --history "$t"-s{2,4}-prior-clock-report-$d.json --point-receipt "$t"-prior-clock-points-$d.json >/tmp/earlier-clock.json
python -W error early_pole_endpoint_cold_prior_clock_population_rectangle.py "${args[@]}" --fixture /tmp/earlier-clock.json
```
