# Advance the selected signed slow closure a short way into expansion

The [fixed-volume signed diagnostic](INDEPENDENT-WARM-REHEATING-SIGNED-SLOW-CLOSURE-2026-11-26.md)
starts from an imposed 100→60 MeV prompt-lepton seed.
Continue this same seed with the adiabatic visible
EOS, retaining the selected signed internal rates and
a common dark temperature with finite fast-decay
equilibrium population. Do not infer how the seed
itself would change during the earlier reheating
history. The Friedmann rate contains the computed
dark energy, but the visible entropy clock is prescribed:
dark inflow does not cool the visible bath in this test.

Use separate light and heavy FD lepton pair inputs
into the two fast-invariant yields. Dark blocking on
visible inflow, hadronic sources, elastic energy exchange
and other reactions are omitted. For dark energy per
visible entropy U the independent ledger is
dU/d ln a = J_E/(H s_visible) − 3 P_dark/s_visible.
The solver reconstructs pressure from the positive thermal
state. An independent numerical integral of input and pressure
work reproduces the energy change to within 2e−14 of
the initial energy on the refined grids.

| Selected continuation | Visible temperature at the comparison stop | S at the stop, per visible entropy | Dark temperature at the stop |
|---|---:|---:|---:|
| Internal signed reactions and expansion only | 58.1946 MeV | 8.6662e−5 | 24.1843 MeV |
| Same plus direct FD lepton inflow | 58.2185 MeV | 8.8400e−5 | 24.3503 MeV |

Both start with S=5.8701e−5 and an initial negative
scalar-equivalent slope, then turn upward on the
short followed interval. The stop is explicitly
triggered by the dilute-closure ceiling of .075 for
the maximum classical phase occupation, reached
by the scalar in these runs. It is **not** freeze-out,
a change in particle physics, or a statement that
any kinetic approximation is accurate past that
point. The small energy fraction of the dark sector
here (below .001) also does not imply a small internal
number or momentum response. There is no final
yield or dark-matter transfer result from this test.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error warm_prompt_reheating_heavy_gain_probe_reference.py EOS.late.dat local-spectra > passive.json
python -W error warm_prompt_reheating_expanding_slow_reference.py EOS.late.dat local-spectra/0.100-0.060-171.json
~~~

[The two arms, solver checks and signed histories](results/warm-prompt-reheating-expanding-slow-2026-11-26.json).

## Local quantum-occupation sensitivity

On the direct-inflow arm, independently reconstruct
FD/BE numbers, energy and pressure by the fugacity
series at each step. For each collision take the forward
occupation multiplier at the lower edge, midpoint, or
upper edge of the local positive bounds, and apply the
exact detailed-balance ratio to get its signed rate.
These choices are three sensitivity trajectories, **not**
global upper and lower solutions of the kinetic problem.

All three give a negative initial S slope and hit the
same dilute phase-factor threshold with visible T between
58.1775 and 58.2365 MeV. The terminal S is between 8.9032e−5
and 8.9338e−5 for these particular variants. The spread is
not a full physical error bar and these calculations do
not integrate the exact quantum collision operator.

~~~sh
python -W error warm_prompt_reheating_expanding_quantum_sensitivity_reference.py \
 EOS.late.dat local-spectra/0.100-0.060-171.json
~~~

[All three occupation-factor sensitivity runs](results/warm-prompt-reheating-expanding-quantum-sensitivity-2026-11-26.json).
