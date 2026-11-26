# Match 12 and 24 corrected q32 steps at the same endpoint

Start from exactly the same two accepted q32 seeds at 8.133100 MeV
and use the same nine-cell reverse partition, supports, visible
table, inherited exterior gravity, and selected prompt and elastic
operators. Only the uniform expansion step changes. Both trajectories
end at 7.394049 MeV after relative expansion 1.10. Compare every
coarse frame against the corresponding pair of finer blocks.
The two raw input and step hash chains, old cohorts and total count
including new births are checked independently.

At the endpoint, the scalar/light number changes in normalized
shell L1 are 9.83e-4 and 7.31e-4. The corresponding energy
changes are 1.148e-3 and 7.24e-4, and pressure 1.148e-3 and
6.97e-4. Each metric is normalized by that initial species and
moment. This temporal sensitivity is substantially larger than
the previously observed selected-grid differences or the size
of the nine-cell correction at a fixed finer time grid.
Neither those comparisons nor this two-level ladder determine
an asymptotic convergence order.

The first cold cell declines at every accepted interval on both
time grids. Its final occupations are 4.1037515e-5 (12) and
4.0977127e-5 (24); the slopes remain between −5.08e-5 and
−4.71e-5 per ln(a) in the two runs. The three coarser fallbacks
were accepted with count conservation, and none was used on
the 24-step trajectory. This preserves the local cold-cell sign,
not the detailed corrected spectrum to the displayed digits.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
mkdir -p time
zcat results/warm-recoil-crossing-q32-terminal-2026-11-26.json.gz > time/seed.json
for n in 12 24; do
  s=$(python -c "print(1.10**(1/$n))")
  python -W error warm_multishell_corrected_continuation_reference.py EOS.dat \
   time/seed.json "time/q32-$n" --step "$s" --blocks "$n" --cells 9
done
python -W error warm_multishell_q32_time_audit_reference.py EOS.dat \
 time/seed.json time/q32-12 time/q32-24 --limit 12
~~~

[Matched time frames](results/warm-corrected-q32-time12-24-2026-11-26.json.gz);
[step 11](results/warm-corrected-q32-time12-step11-2026-11-26.json.gz),
[12-step endpoint](results/warm-corrected-q32-time12-step12-2026-11-26.json.gz)
and [24-step endpoint](results/warm-multishell-trajectory-lifted24-q32-2026-11-26.json.gz).
