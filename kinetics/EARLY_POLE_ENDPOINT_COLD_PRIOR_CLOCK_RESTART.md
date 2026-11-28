# An earlier coarse interval changes the cold amplitude even at the same final clock

Recover the **actual** coarse branch on x=.02→.04 independently for each of S2 and S4, from the corresponding existing .02 state. The eight recovered ledger samples match the original coarse trace exactly. Then evolve to .05 with **2048** steps, the unchanged final clock of the original fine trajectory, starting from the recovered .04 state. The alternative terminal archive includes the SHA-256 of its actual prior state. The original path is (2048,2048) and the alternative here is (1024,2048); both hold charge partition D2 and the respective scalar geometry fixed.

| Change, alternative over original minus one | S2 | S4 |
| --- | ---: | ---: |
| Terminal scalar yield (fraction) | −5.15050e-5 | −5.09152e-5 |
| Fitted scalar gap (fraction) | −1.71689e-5 | −1.69723e-5 |
| Cold forward at q=3.5e-6 | −0.460785% | −0.454579% |
| Cold forward at q=4.7e-6 | −0.452604% | −0.446808% |

With both clocks preserving their cold cell geometry, the *count-normalized scalar occupation measure* across [1e-6,1.5e-5] differs in TV by only 3.74642e-6 for S2 and 5.17487e-6 for S4. These are occupation TVs, not collision-response TVs. The later point quadratures close to 7.41e-13 or better, and inverse responses decrease by roughly 0.48–0.50%. For S4 the ~0.45% earlier-window point effect exceeds that of changing only the final step at the same q, even though the final step is now held fixed. Prior timing therefore matters for absolute amplitudes; the two-point response cannot bound the integrated band or the full time history.

From `kinetics/`, regenerate either branch using `early_pole_endpoint_cold_prior_clock_restart.py --level 2` or `--level 4`, with `--inputs` pointing to the three native .02/.04/.05 states and the two native .04/.05 paired reports in order; specify `--earlier` and `--terminal` output paths. To reproduce the point receipt with stored artifacts:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
s=results/early-pole-fine-charge2-scalar; p=results/early-pole-cold-s; d=2026-11-28
python -W error early_pole_endpoint_cold_prior_clock_points.py \
 --states "$p"2-prior-clock-x{04,05}-state-$d.json.gz "$s"2-x05-state-$d.json.gz \
          "$p"4-prior-clock-x{04,05}-state-$d.json.gz "$s"4-x05-state-$d.json.gz \
 --reports "$p"{2,4}-prior-clock-report-$d.json
```
