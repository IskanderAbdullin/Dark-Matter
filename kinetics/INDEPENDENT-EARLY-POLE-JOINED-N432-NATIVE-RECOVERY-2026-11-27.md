# The native 432-step result passes with a tighter elastic acceptance budget

The [original native run](INDEPENDENT-EARLY-POLE-JOINED-N432-LEDGER-REJECTION-2026-11-27.md) completed all steps but failed the terminal energy gate by an order of magnitude. Reuse the same archived native hard reconstruction, channel supports, unlike active response, collision amplitudes, exposure 7.2e−5 and 432 macro steps. Require each backward elastic solve to meet an additional absolute energy defect limit of 5e−15. This does not relax the residual, terminal energy or terminal entropy gates, change a rate, or project the final state to enforce a preferred answer.

The reason to test the nonlinear acceptance separately is the [discrete-gradient chain rule](https://web.ma.utexas.edu/users/og/PUBLICATIONS/paper_DisHamSys.pdf). If the step increment is ΔY and the collision term is D, the energy defect is the sum of the discrete-energy gradient projected on the residual ΔY−hD and the energy work of hD. At matched frozen hard snapshots in both reconstructions, the unrefined late elastic steps have signed energy defects of 8.1e−14 to 2.3e−13. The independently reconstructed residual projection explains them within 1e−15, while collision work is below 1e−18 on those late snapshots. One extra Newton update at each of those snapshots reduces the energy defect below 7e−16.

| Native fine hard history | Terminal energy defect | Integrated scalar gain |
| --- | ---: | ---: |
| Accepted 216 macro steps | +2.3914e−13 | 7.845251615e−6 |
| Rejected default 432 macro steps | +2.2693e−11 | Not accepted |
| Accepted 432 with elastic budget | +2.0872e−14 | 7.847737198e−6 |

The tighter 432 arm passes the original 3e−12 energy and 3e−11 entropy gates (entropy defect −1.5224e−12). Its signed summed elastic defect is +2.1094e−14 and peak absolute elastic defect 4.8850e−15. At the **same** grid width, the tighter budget changes the rejected run's stored scalar count history by at most 1.05e−17. That is a numerical diagnostic only, not acceptance of the older run. Comparing the two accepted resolutions, the final gain rises by .03168% and the final instantaneous cold source by .0136%. The early interval rates differ more than the final integrated gain. Two accepted widths do not establish a convergence order, and no continuous-hard-input 432-step conclusion is included here.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
b=results/early-pole-joined-fine-216-memory-2026-11-27.json.gz
python -W error early_pole_cold_screened_split_reference.py fine \
  --duration 7.2e-5 --count 432 --same-soft-strength 1 --compton-strength 1 \
  --conversion-strength 1 --screening-update full --pair --projection steps --snapshots \
  --elastic-energy-budget 5e-15 --checkpoint native432-tight.json --checkpoint-every 12 > native432-tight-accepted.json
python -W error early_pole_joined_split_checkpoint_ledger_reference.py native432-tight.json
python -W error early_pole_joined_temporal_refinement_reference.py "$b" native432-tight-accepted.json --projection steps
for p in steps cell-conservative; do
  python -W error early_pole_elastic_energy_residual_projection_reference.py "$b" --projection "$p"
done
```

Raw: `results/early-pole-joined-n432-native-energy-refinement-2026-11-27.json.gz`.
