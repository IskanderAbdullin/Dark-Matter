# Two-dimensional field control with identical seeded Fourier modes

Extend the [one-dimensional coherent field control](INDEPENDENT-EARLY-POLE-COHERENT-QUARTIC-LATTICE-2026-11-26.md)
to a periodic plane, retaining the de-aliased quartic dynamics
and symplectic time integration. Initialize the same finite set
of plane-wave coefficients at each spatial resolution. The box,
initial spectrum and phases are fixed while numerical time and
space resolution are changed. Initial amplitudes are not inferred
from the early kinetic snapshot.

The angular multi-mode initial condition does not have to grow
at the peak single-mode exponent. Solve independent linear
oscillators for the actual seed weights against the exact
homogeneous quartic background, and apply the same growth fit.
For the matched runs the predicted amplitude exponent is about
0.03226–0.03227 in oscillator units; the nonlinear field evolves
initially at 0.03226–0.03229. It exceeds 1% inhomogeneous energy
at about 319 oscillator times and reaches a sampled maximum
near 30% by the specified 400-time endpoint. The stable seeded
band stays below 10⁻¹⁰ inhomogeneous energy in this test.
Halving the time step at 96 points per axis reduces the energy
drift approximately fourfold. A distinct 144-point mesh gives
a similar early slope and nonlinear energy readout.

An initial check comparing directly to the single-mode peak
rejected a perfectly plausible multi-mode run. Matching the
seed spectrum in the separate linear calculation is a more
appropriate control. This restricted coherent setting does not
show that the real dark scalar phase is initially coherent or
that its nonlinear or asymptotic spectrum follows the lattice.
No three-dimensional convergence is established here.

~~~sh
python -W error early_pole_coherent_quartic_lattice_2d_reference.py --end 400
python -W error early_pole_coherent_quartic_lattice_2d_reference.py --cells 144 --step .04 --end 400
python -W error early_pole_coherent_quartic_lattice_2d_reference.py --stable --end 400
~~~

Raw runs: [base](results/early-pole-coherent-lattice-2d-main-2026-11-26.json.gz),
[finer space](results/early-pole-coherent-lattice-2d-fine-2026-11-26.json.gz),
[halved time step](results/early-pole-coherent-lattice-2d-time-2026-11-26.json.gz),
[out-of-band control](results/early-pole-coherent-lattice-2d-stable-2026-11-26.json.gz).
