# Separate time and outer support after the nine-cell reverse correction

The selected model is unchanged: the corrected scalar-pair sector,
selected prompt/elastic operators, tabulated visible radiation and the
inherited collisionless exterior used as gravity only. All runs start
at about 8.133100 MeV and end at 7.394049 MeV after expansion by 1.10.
The two comparisons below change one numerical control at a time.

## Time resolution on q32

Follow each trajectory's own accepted history, using the same input
seed and 12, 24 or 48 uniform ln(a) steps. Compare only identical
physical frames. Normalizing by the initial integral of each moment
and species, the largest endpoint difference in shell L1 falls from
0.1148% (12 vs 24) to 0.05786% (24 vs 48). The endpoint apparent
orders for number, full energy and pressure are 0.9856–0.9887 over
both species. Over all noninitial shared frames the range is
0.9380–0.9887; the two successive signed shell-error vectors have
cosines at least 0.9981. These observations support approximately
first-order behavior over the tested grids, not a continuum error
bound or a proof of the asymptotic order.

The first cold-cell occupation decreases in all intervals in all
three runs. The endpoints are 4.1037515e-5, 4.0977127e-5 and
4.0946258e-5, respectively. All 48 fine steps used their own history.
The fine run closes the inclusive count ledger to 3.28e-13.

## Outer support at 24 steps

Independently advance the larger q64-support state on the same
24-step time grid. The accepted q32 profiles enter only as nonlinear
initial guesses: q64 is not constrained to those values and accepts
three fallback steps without the narrow guess. Both cold cells
decrease in all 24 intervals. At the last checkpoint q32 has
geometrically zero first-cell prompt rate; the q64 representation
has 60 contributing samples and a first-cell prompt occupation rate
of −5.738e-9 per Hubble time. That difference does not restore the
native uncorrected positive step.

At the endpoint the scalar/light shared-grid number differences
are 3.49e-6 / 1.16e-6 of their respective initial q32 counts; the
L1 of the **change in the signed difference from its initial
offset** is 1.06e-6 / 2.94e-7. The full-energy and pressure metrics
in this support audit both use the initial total q32 energy, not
the per-species normalization of the time audit. q64-only cells
contain 7.35e-12 / 2.27e-11 of its retained scalar/light counts at
the endpoint. This is a test of the two represented supports,
not a bound on untracked incoming phase space or angular error.

## Reproduce from the exact published inputs

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
mkdir -p ladder
zcat results/warm-recoil-crossing-q32-terminal-2026-11-26.json.gz >ladder/seed32.json
for n in 12 24 48; do
  s=$(python -c "print(1.10**(1/$n))")
  python -W error warm_multishell_corrected_continuation_reference.py EOS.dat \
    ladder/seed32.json "ladder/q32-$n" --step "$s" --blocks "$n" --cells 9
done
python -W error warm_multishell_q32_three_level_time_reference.py EOS.dat \
  ladder/seed32.json ladder/q32-12 ladder/q32-24 ladder/q32-48
python -W error warm_q64_multishell_guided_continuation_reference.py EOS.dat \
  results/warm-recoil-crossing-q64-terminal-2026-11-26.json.gz ladder/q32-24 ladder/wide --blocks 24
python -W error warm_multishell_corrected_support_audit_reference.py EOS.dat \
  ladder/seed32.json ladder/q32-24 ladder/wide 24
~~~

[Three-level audit](results/warm-corrected-q32-time12-24-48-2026-11-26.json.gz);
[q32 step 47](results/warm-corrected-q32-time48-step47-2026-11-26.json.gz);
[q32 step 48](results/warm-corrected-q32-time48-step48-2026-11-26.json.gz);
[support audit](results/warm-multishell-corrected-support-audit-2026-11-26.json.gz);
[q64 step 23](results/warm-multishell-corrected-q64-stage-23-2026-11-26.json.gz);
[q64 step 24](results/warm-multishell-corrected-q64-stage-24-2026-11-26.json.gz).
