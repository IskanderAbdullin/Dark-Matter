# Freezing pair chemistry restores drag at leading order in sound damping

Repeat the [static longitudinal control](INDEPENDENT-CHARGED-LEPTON-LOCAL-LONGITUDINAL-2026-11-27.md)
with **both** directions of the pair-changing reaction removed,
leaving only charged-lepton elastic heat and momentum exchange.
The background remains at common T and η_D=0, but pair density is
now conserved, including its perturbation. This is a limiting test
of the closure, not an inferred chemical decoupling history.

The slow manifold has three conserved coordinates: total energy,
total momentum and dark particle number. Its leading streamed
modes are an acoustic pair and a nonpropagating composition mode.
Even though the background η_D is zero, the conserved-number sound
perturbation need not follow δη_D=0. In particular the extra
n_D T δη_D contribution to dark pressure can split the two
accelerations δP_i/h_i. Relative velocity therefore contributes already
to the k² attenuation, unlike the [chemically equilibrated control](INDEPENDENT-CHARGED-LEPTON-LOCKED-SOUND-DAMPING-2026-11-27.md).

Project the full two-fluid entropy metric onto the three conserved
directions and compare its positive quadratic attenuation to actual
small-k eigenvectors using the independent entropy-loss quotient.
The sound speed also agrees with an independent constrained
thermodynamic slope at fixed total energy-to-number compression.
The direct damping coefficients match the projections to below 5e−5
relative across the tested controls.

| Common T (GeV), h_V/h_D=10 | Drag fraction of leading sound damping | Ratio of leading damping if setting K_p=K_Q |
|---:|---:|---:|
| .015 | .7640 | 2.2808 |
| .050 | .6435 | 2.0461 |
| .100 | .6251 | 1.9625 |
| .300 | .6177 | 1.9157 |
| .6882487 | .6232 | 1.8754 |

The nonpropagating composition mode has only vector damping at this
order in this model; that part changes inversely with the chosen
drag strength, as expected for diffusive damping after relative
velocities are eliminated. The displayed dimensionless damping
coefficients use Γ_Q=K_Q/h_D **for the active elastic** heat channel,
not the combined rate of the chemically equilibrated control.
The limiting calculations do not bound an evolving finite
annihilation rate, or include scalar free streaming and gravity.

Replay from the kinetics directory:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error \
  charged_lepton_frozen_chemical_acoustic_reference.py \
  results/charged-lepton-portal-combined-momentum-2026-11-27.json
```
[Full result](results/charged-lepton-frozen-chemical-acoustic-2026-11-27.json).
