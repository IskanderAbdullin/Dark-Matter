# Conditional coherent-field lattice and a phase-mixing counterexample

The [coherent linear diagnostic](INDEPENDENT-EARLY-POLE-COHERENT-QUARTIC-FLOQUET-2026-11-26.md)
can be benchmarked beyond linear order without assuming it represents
the actual early production state. In the pure quartic limit evolve
a real field in a periodic line using a spectrally de-aliased cubic
force and velocity Verlet. Time is in homogeneous oscillator units
and field amplitude initially has a uniform part equal to one.
Expansion, fermions and Bose particle collision terms are absent.

A single resonant seed at k/Ω=1.272 grows at 0.03604–0.03605 Ω
before backreaction, consistent with the analytic Lamé band. On
512 and 1024 cells the off-uniform energy first exceeds 1% around
tΩ=308. The late sampled inhomogeneous fraction peaks near 29%
over the specified interval, and the quadratic resonant-band
contribution reaches about 27%. An out-of-band seed stays tiny
over the same interval. Changing the resonant seed by a factor
of ten shifts the 1% time by roughly 63 oscillator units, close
to ln(10)/0.03598. Relative energy drift in the refined aligned
run is below 3.2 × 10⁻⁴; no claim of a late steady state is made.

The large-box counterexamples add a deliberately chosen slow
amplitude modulation at k/Ω≈0.00497. A 1% modulation can make
nearly all energy *appear* to leave the uniform mode, while less
than 1.6% reaches k/Ω≥0.6 in the quadratic measure. At 10%
modulation the latter stays below 0.1% in the tested interval.
A no-spatial-gradient comparator made of independent exact
quartic oscillators already shows this cancellation of the
global mean and agrees within 0.1 percentage point during
the first 45 oscillator units. The decrease of mean-field
energy alone is therefore not evidence for short-scale power
or fragmentation. The added modulation is a hypothesis, not
a deduction of the phase from the kinetic spectra.

As a separate frozen check, the stored isotropic spectrum
implies a spatial two-point function proportional to
Σ(n/E) sinc(qr) under independent phases, and time
correlator proportional to Σ(n/E) exp(iEt) with a fixed pole.
Normalized correlation at one conditional resonant wavelength
exceeds 0.99992 in the tested snapshots. The narrow cold
frequency spread alone also allows a large temporal envelope
over 100 inverse conditional Floquet rates; a variance bound
verifies this without resolving high-frequency tail phases.
This does not assign a coherent mean field, a condensate
fraction, or phase correlation in the nonlinear medium.

These periodic field controls do not provide a frozen particle
spectrum, a chemical rate, or a structure prediction. Related
nonlinear scalar benchmarks are discussed by
[Micha and Tkachev](https://arxiv.org/abs/hep-ph/0403101) in a
different preheating setting.

~~~sh
python -W error early_pole_coherent_quartic_lattice_reference.py
python -W error early_pole_coherent_quartic_lattice_reference.py --cells 1024 --step .03
python -W error early_pole_coherent_quartic_lattice_reference.py --mode-count 256 --cells 4096 --step .05 --end 330 --modulation .01
python -W error early_pole_coherent_phase_mixing_reference.py results/early-pole-coherent-lattice-wide-1pct-2026-11-26.json.gz
python -W error early_pole_scalar_equal_time_coherence_reference.py results/early-pole-spectrum-third-m06-2026-11-26.json.gz
~~~

Raw runs: [aligned](results/early-pole-coherent-lattice-main-2026-11-26.json.gz),
[refined](results/early-pole-coherent-lattice-refined-2026-11-26.json.gz),
[stable](results/early-pole-coherent-lattice-stable-2026-11-26.json.gz),
[weak seed](results/early-pole-coherent-lattice-weakseed-2026-11-26.json.gz),
[strong seed](results/early-pole-coherent-lattice-strongseed-2026-11-26.json.gz),
[wide 1%](results/early-pole-coherent-lattice-wide-1pct-2026-11-26.json.gz),
[wide 1% finer time](results/early-pole-coherent-lattice-wide-1pct-timehalf-2026-11-26.json.gz),
[wide 10%](results/early-pole-coherent-lattice-wide-10pct-2026-11-26.json.gz),
[no-gradient check](results/early-pole-coherent-phase-mixing-2026-11-26.json.gz),
[stored-spectra correlations](results/early-pole-scalar-frozen-covariance-2026-11-26.json.gz).
