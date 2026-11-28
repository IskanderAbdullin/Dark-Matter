# A third warm Gaussian node changes the measured stress direction

Independently evolve the positive three-node volume-Gauss chart from the same empty-daughter start, with the cold scalar chart and charge-d2 geometry unchanged. There are 189 scalar, 448 fermion and 488 parent nodes. Only the 22 complete warm cells change order. This is distinct from the frozen order ladder: all retained nodal occupations can evolve independently. The `gauss-volume-v2-general` state stamps its actual previous uncompressed restart with SHA-256. The chain receipt verifies all four links, chronological clocks, cell counts and paired/continuity/Pauli ledgers. It does not replace a full independent replay of the solver history.

| x | Gauss3/Gauss2 pressure change at held W1 mass | Gauss3/W4 pressure change |
| --- | ---: | ---: |
| .01 | +0.09232% | −0.49109% |
| .02 | +0.09469% | −0.48946% |
| .04 | +0.10536% | −0.49085% |
| .05 | +0.08266% | −0.49031% |

At x=.05 total number differs from Gauss2 by just +6.418e−9 relative. Given each warm cell's actual count and mean energy, the allowable end pressure contrast is [+0.0826486%, +0.0826646%] at the common measurement mass. That narrow interval concerns reconstruction only; it does not bound time stepping, momentum resolution or physical model error.

The exact same positive reconstructed distributions were held fixed while **formally** changing the probe mass M. The late terminal Gauss3-to-Gauss2 pressure contrast is positive at the archive mass, but **negative** at M=10 (−0.08008%) and M=100 (−0.12914%). The pseudo-pressure changes by −0.78698% and −1.34140%, respectively. These are opposite stress movements under different numerical observation weights at fixed distribution, not physical redshift evolution, perturbation sound speeds, or CLASS results. The old Gauss2 and geometric receipts reproduce unchanged after generalizing the postprocessor. A roundoff allowance of four parts in (10^{15}) in the massless physical inequality test only admits representational equality; the paired integration/number error gates are unchanged.

Reproduction from `kinetics/`, with the existing endpoint dependencies:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
# Startup and then the same exact time segments as the Gauss2 comparison.
python -W error early_pole_endpoint_higher_warm_gauss_history.py --warm-order 3 \
  --steps 2048 --state /tmp/g3-start.json
# Subsequent calls carry --restart previous-fine-state and --x-end .01/.02/.04/.05,
# with --steps 512/1024/1024/1024, and write each fine state with --state.
p=results/early-pole-fine-charge2-scalar4; d=2026-11-28
python -W error early_pole_endpoint_gauss_warm_chain_receipt.py --order 3 \
  --reports "$p-gauss3warm-"{start,x01,x02,x04,x05}-$d.json.gz \
  --states "$p-gauss3warm-"{start,x01,x02,x04,x05}-state-$d.json.gz
python -W error early_pole_endpoint_warm_gauss_order_contrast.py \
  --order2 "$p-warm4-gauss-comparison-x"{01,02,04,05}-$d.json \
  --order3 "$p-gauss3warm-comparison-x"{01,02,04,05}-$d.json \
  --stress results/early-pole-frozen-stress-sweep-w124-gauss-$d.json \
           results/early-pole-frozen-stress-sweep-w124-gauss3-$d.json
```

## Direct final-interval pressure check

The recovered 1024-step endpoint from the identical archived .04 Gauss3 restart reproduces every audited coarse trace entry. At the held W1 mass its pressure coarse/fine minus one is (+2.35630	imes10^{-7}), while count shifts by (-3.18010	imes10^{-5}). The Gauss3/Gauss2 pressure contrast changes by only (1.86219	imes10^{-8}) in relative fraction between coarse and fine time grids. This directly rules out the final time-step choice as the explanation of the +0.08266% spatial-order difference in this particular restart test. It is not a bound on time-error accumulation earlier in either history. Reproduce with `early_pole_endpoint_gauss_warm_pressure_clock.py` using the Gauss3 .04 state, .05 state and .05 paired trace, then pass its receipt as `--gauss3` along with the other four inputs to `early_pole_endpoint_warm_pressure_clock_contrast.py`.
