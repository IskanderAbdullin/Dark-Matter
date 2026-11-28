# Separate time resolution from the geometric shutoff of the native recoil gate

The [finer matched crossing](INDEPENDENT-WARM-RECOIL-CROSSING-MATCHED-2026-11-26.md)
started at two stored 8.568461 MeV states. The original coarser
native q32 and q64 continuations now extend through the same
endpoint. Starting from their exact shared states at that
temperature, compare eight original intervals with sixteen finer
ones down to 8.133100 MeV. On the shared cells the maximum
relative shell L1 time difference is 9.09e-4, versus 3.55e-6
for the tested support change. Forming the signed difference of
differences before taking the norm gives a mixed effect 6.87e-8.

At the matched 8.239831 MeV frame the q32 first-cell native
prompt rates are −1.491e-9 and −1.438e-9 per Hubble time.
At the 8.133100 MeV endpoint both q32 rates are exactly zero;
the matched q64 rates are −1.158e-7 and −1.111e-7.
The first light-cell average still declines at all the common
intervals in all four native trajectories. These are independent
time and support checks before the later corrected pair trajectory.

## What happens right at the finite cell boundary?

Freeze the archived fine q64 state at 8.239831 MeV and
vary only a common daughter momentum cut near the q32 gate.
Keep the original piecewise entropy reconstruction, truncating
its last cell without refitting. Hold the reference H fixed.
Write δ = E_cap − E_min(p₁), where p₁ is the outer edge
of the first light cell. The threshold is 0.2395511773 GeV
in scalar energy for this fixed snapshot. Below it there are
no exact first-cell events. At equality, the native floating-point
rate remains below 6e-30 per Hubble time.

Above the threshold the allowed scalar-energy interval has
length approximately δ−g(p₁−p), with g=−E_min'(p₁)>0.
Integrating this narrow wedge in light momentum and recoil
energy gives a leading factor δ²/(2g). The rest of the
integrand is evaluated at the **inner** one-sided value of the
first light cell, because adjacent reconstructed affine entropy
profiles need not join at a shell edge.

The resulting leading rate per H / δ² is −0.746705 GeV⁻².
At δ=2e-8 GeV the independent two-dimensional quadrature
recovers the leading prediction within 2.5e-5 relative; the
native sampled cell agrees with the independent rate within
2.3e-5. The larger finite-margin rows show nonquadratic
corrections. This explains a smooth shutoff in the exact
projected cold-cell prompt term, not a universal power law
for an evolving spectrum, an open exterior, or any other channel.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
mkdir -p coarse
s=$(python -c 'print(1.11**(1/16))')
for q in 32 64; do
  zcat "results/warm-history-lower-q$q-16-terminal-2026-11-26.json.gz" >"coarse/seed$q.json"
  python -W error warm_q64_history_serial_reference.py EOS.dat \
    "coarse/seed$q.json" "coarse/q$q" --step-expansion "$s" --blocks 16
done
# Use crossing/ from the linked finer-ladder instructions.
python -W error warm_native_recoil_time_support_cross_reference.py EOS.dat coarse crossing
python -W error warm_cold_recoil_gate_scaling_reference.py EOS.dat \
  results/warm-native-recoil-gate-q64-2026-11-26.json.gz
~~~

[Matched clocks and spectra](results/warm-native-recoil-time-support-cross-2026-11-26.json.gz);
[fixed-gate diagnostic](results/warm-native-recoil-gate-scaling-2026-11-26.json);
[exact gate snapshot](results/warm-native-recoil-gate-q64-2026-11-26.json.gz);
[coarse q32 terminal](results/warm-native-recoil-coarse-terminal-q32-2026-11-26.json.gz);
[coarse q64 terminal](results/warm-native-recoil-coarse-terminal-q64-2026-11-26.json.gz).
