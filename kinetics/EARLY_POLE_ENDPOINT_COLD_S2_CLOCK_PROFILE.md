# S2's final time step moves the cold band's amplitude far more than its shape

From the exact same archived x=.04 S2 state and charge grid, recover the original 1024-step coarse branch and compare with the archived 2048-step fine branch at x=.05. The 1024-step recovery matches the stored time trace, and both held response archives use the recovered optical parameters. Pair full q-band integrals on [1e-6,1.5e-5]; the coarse raw panel archive has 57 panels per mesh. Use the exact shared 17 comparison bins to measure normalized shapes with no interpolation.

| Quantity | Coarse/fine change |
| --- | ---: |
| Total scalar yield | −3.16325e-5 |
| Forward cold number response | −0.282106% |
| Inverse cold number response | −0.306550% |
| Forward cold kinetic excess response | −0.280174% |

Normalizing each chart's profile by its own corresponding band total leaves total-variation distances of 8.248e-6 for forward number, 1.217e-5 for inverse number, and 1.139e-5 for forward excess. The paired quadratures reproduce these distances to better than 7.39e-14 absolute. The same numerical clock change has a visible amplitude effect and a much smaller normalized shape effect. This is within **one** S2 chart; it neither bounds earlier clock history nor substitutes for the matched S2-versus-S4 comparison.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results/fine-charge2-scalar2; s=results/early-pole-fine-charge2-scalar2; d=2026-11-28
python -W error early_pole_endpoint_cold_one_chart_clock_profile.py \
 --coarse "$p"-coarse-clock-endpoint-continuum-cold-{1,2}.jsonl.gz \
 --fine "$p"-endpoint-continuum-cold-{1,2}.jsonl.gz \
 --states "$s"-x04-state-$d.json.gz "$s"-coarse-x05-state-$d.json.gz "$s"-x05-state-$d.json.gz "$s"-x05-$d.json.gz \
 --recovery "$s"-coarse-x05-recovery-$d.json \
 --archive results/finite-mass-qed-fixed-pole-cold-scalar-band-$d.json
```
