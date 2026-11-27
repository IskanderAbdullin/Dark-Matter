# An energy-matched Gaussian field is not a spatially uniform coherent field

Make a different, explicit initial-state hypothesis from the
[uniform coherent limit](INDEPENDENT-EARLY-POLE-COHERENT-QUARTIC-FLOQUET-2026-11-26.md).
Construct a zero-mean Gaussian random-phase ensemble with
⟨φ²⟩=Σn/E, ⟨π²⟩=ΣnE and ⟨|∇φ|²⟩=Σnq²/E.
Wick's fourth moment and the Hartree gap imply that its
canonical scalar energy exactly equals the saved Hartree
energy. The identity matches the energy, not higher
correlations or the actual quantum state. The reference
field uses only the retained excitation covariance.

Ignore gradients just for the local band calculation and
assign a separate quartic oscillator to each sampled pair
(φ,π). The mean local energy then omits the explicit gradient
part, about 0.5–0.9% of total for these snapshots. The associated
turning amplitudes and Floquet bands vary from patch to patch.
For the .06 case at x=.05, the homogeneous energy-matched
midband wavenumber lies in the local bands of patches carrying
about 8.2% of the initial energy. Maximizing over the tested
wavenumbers gives about 12.7% energy coverage, not simultaneous
resonance of the entire ensemble. The maximum volume coverage
is about 8.1%, a different weight. These are overlap diagnostics,
not the fraction ultimately transferred.

As an independent linear time-domain control, evolve each
local oscillator and two orthogonal infinitesimal perturbation
quadratures for a fixed probe momentum. After 240 uniform
oscillator times, about 7.15% of the energy in the 32,768
patch sample lies in patches amplifying the probe norm by
more than 100 at q/Ωuniform=1.272. At 1.573 it is about 11.72%.
This is seed-normalized **linear** gain without saturation.
Fields between patches are not coupled; gradients, transit
through patches, nonlinear backreaction, fermions and ongoing
production are omitted. No actual fragmented fraction, rate,
or observable spectrum is inferred.

~~~sh
python -W error early_pole_gaussian_field_energy_handoff_reference.py results/early-pole-spectrum-third-m06-2026-11-26.json.gz
python -W error early_pole_gaussian_local_field_growth_reference.py results/early-pole-spectrum-third-m06-2026-11-26.json.gz --refined
~~~

Raw results: [ratio-.06 handoff](results/early-pole-gaussian-handoff-third-m06-2026-11-26.json.gz),
[other momentum grid](results/early-pole-gaussian-handoff-fine-m06-2026-11-26.json.gz),
[ratio-.30 handoff](results/early-pole-gaussian-handoff-third-m30-2026-11-26.json.gz),
[x=.10 handoff](results/early-pole-gaussian-handoff-third-extend-2026-11-26.json.gz),
[linear probes and refinements](results/early-pole-gaussian-local-growth-2026-11-26.json.gz).
