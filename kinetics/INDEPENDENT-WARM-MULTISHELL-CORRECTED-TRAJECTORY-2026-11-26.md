# Compare full matched q32 trajectories with and without the nine-cell lift

The [one-block replay](INDEPENDENT-WARM-COLD-PAIR-MULTISHELL-2026-11-26.md)
isolated the reverse projection error at the native cold spike. Here
independently restart **both** q32 branches at the 8.133100 MeV
accepted pre-spike state. Use the same two-states history initialization,
radiation table, selected operators, shared momentum support and
24 identical log-expansion intervals to relative expansion 1.10.
Only the reverse scalar-pair sector feeding the first nine light cells
changes. Initial profiles come from the same verified raw seed,
not from any later native snapshot. The corrected run uses its own
accepted predecessor for every block.

All stages check independent count closure inclusive of the new
collisionless exterior births. The inherited exterior is carried
as gravity only and its archived cohorts are verified without being
reactivated. The q32 prompt recoil window remains geometrically
closed in the first retained cell in both branches.

Both branches reach 7.394049 MeV. The native first-cell average rises
only at step 12; across all 24 steps the nine-cell version declines
smoothly between -5.080e-5 and -4.794e-5 per unit ln(a). All modified
nonlinear blocks use their own prior history without fallback.
The maximum normalized whole-grid number, energy or pressure shell
L1 change is 6.74e-07 during the run and 5.57e-07 at the endpoint,
using each initial species/moment as its denominator. The corrected
inclusive number and local energy ledgers close to
1e-13 and 2.5e-14,
respectively.

This is a genuine matched continuation, not a local replay near
the spike alone. It tests sensitivity to a known discretization
error, with the same native quadrature elsewhere. It is not a step-size
convergence claim for the entire track, a q64 trajectory or a final
cosmological prediction.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
mkdir -p trajectory
zcat results/warm-recoil-crossing-q32-terminal-2026-11-26.json.gz > trajectory/seed.json
STEP=$(python -c 'print(1.10**(1/24))')
python -W error warm_q64_history_serial_reference.py EOS.dat trajectory/seed.json trajectory/native --step-expansion "$STEP" --blocks 24
python -W error warm_multishell_corrected_continuation_reference.py EOS.dat trajectory/seed.json trajectory/lifted --step "$STEP" --blocks 24 --cells 9
python -W error warm_multishell_corrected_trajectory_audit_reference.py EOS.dat trajectory/seed.json trajectory/native trajectory/lifted --steps 24
~~~

[Full moment and rate comparison](results/warm-multishell-matched-trajectory-q32-2026-11-26.json.gz);
native [end at step 24](results/warm-multishell-trajectory-native24-q32-2026-11-26.json.gz);
lifted [spike-stage check](results/warm-multishell-trajectory-lifted12-q32-2026-11-26.json.gz),
[lifted step 23](results/warm-multishell-trajectory-lifted23-q32-2026-11-26.json.gz)
and [lifted endpoint](results/warm-multishell-trajectory-lifted24-q32-2026-11-26.json.gz).
