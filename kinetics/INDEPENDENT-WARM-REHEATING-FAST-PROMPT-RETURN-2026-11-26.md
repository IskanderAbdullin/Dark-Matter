# Short prompt-return relaxation at the hotter anchor

Start from the direct-lepton seed for 100 MeV reheating,
sampled when the visible bath reaches 60 MeV. Freeze
the expansion, direct sources and other reactions. Evolve
only the fast inverse-decay/prompt-return redistribution
on a shared momentum box, with separate counts and summed
energy conserved. This is a timescale-separation test,
not a model of the full reheating trajectory.

The 36-cell solution with a retained p/T cutoff of
16.918 transfers about 7.85% of the initial total
dark energy to scalars within 5e−9 Hubble times.
Halving the requested step changes that to 7.86%.
An independent 24-cell projection gives 7.75%.
The 36-cell fixed-box entropy gap is 95.1% closed,
but the spectra still differ from their own fixed-count
maximum-entropy comparison (total variation about
0.028 and 0.033 for scalars and stable fermions).
Separate counts and energy agree within numerical
roundoff across the tested steps.

Fold the same three gross regeneration kernels
on the reconstructed spectra before and after this
short relaxation without feeding them into the evolution.
At 36 cells the gross rate increases by about 22.1%
on these matched-support distributions. The agreement
under halved time steps tests the local conclusion,
not a single thermal closure or a full signed rate.
This reinforces the [entropy-only comparison](INDEPENDENT-WARM-REHEATING-MATCHED-ENTROPY-GAIN-2026-11-26.md),
but does not establish that a static equilibrium
follows over a Hubble time.

The fast redistribution and the much slower
advance of the number ledger should be handled with
distinct timesteps or an independently checked implicit
coupling. The box excludes about 9.7e−7 of initial
scalar counts and 1.85e−5 of stable counts.
No late source, Bose condensate, in-medium width,
coupled background or final dark-matter abundance is solved.

At 1.5e−8 Hubble times the corresponding 36-cell
control closes 99.1% of the fixed-box entropy gap
and transfers 8.64% of the original dark energy.
The gross regeneration rate folded on that terminal
spectrum is 25.9% above its own initial value.

[State, energy and rate controls](results/warm-reheating-fast-prompt-return-2026-11-26.json).

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error warm_prompt_reheating_heavy_gain_probe_reference.py EOS.late.dat local-spectra > passive.json
python -W error joint_prompt_two_moment_evolution_reference.py EOS.late.dat \
 local-spectra/0.100-0.060-171.json --TR .1 --TF .06 --step 4 --cap 17.7 \
 --height 5e-9 --dt 2.5e-10 --verify --record-terminal > terminal.json
python -W error warm_prompt_reheating_relaxed_regeneration_reference.py terminal.json EOS.late.dat
~~~

The [signed slow comparison](INDEPENDENT-WARM-REHEATING-SIGNED-SLOW-CLOSURE-2026-11-26.md)
adds the competing inverse channels on a separate
common-temperature closure, without extrapolating the
short shape relaxation over cosmological times.
