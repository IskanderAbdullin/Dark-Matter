# Annihilation and elastic exchange separate differently in momentum and heat

At common visible and dark temperature and equilibrium dark fugacity,
linearize the charged-lepton vacuum portal with respect to a relative
bulk velocity and a fractional temperature offset. A timelike
annihilation or inverse event transports the entire lab-frame
four-momentum of the pair. With E the event's lab energy and u its
invariant mass, the average squared momentum along any given axis is
(E²−u²)/3. Its thermal weight is E². Hence the positive equilibrium
integral gives K_p/K_Q strictly below 1/3 for the timelike channel.
By contrast, in [elastic scattering](INDEPENDENT-CHARGED-LEPTON-FD-BULK-DRAG-2026-11-27.md)
the transported four-vector is the spacelike difference of particle
momenta, so the corresponding ratio is above 1/3. Angular isotropy
removes scalar-to-vector couplings in both channels.

The equilibrium forward/reverse weight, including Pauli blocking and
spin correlations, is shared within each channel. The timelike
implementation sums E², |P_pair|² and u² independently and checks
their exact relation. Compare with a separate adaptive integration
for both thermal and spatial moments. Both new source options are
restricted to common temperature and zero dark fugacity, leaving
the prior number/heat interfaces unmodified by default.

| Common T (GeV) | Timelike K_p/K_Q | Combined K_p/K_Q | Elastic increment relative to timelike momentum | Elastic increment relative to timelike heat |
|---:|---:|---:|---:|---:|
| .015 | .10593 | .41820 | +349.4% | +13.83% |
| .050 | .10502 | .41809 | +354.6% | +14.18% |
| .100 | .10769 | .44729 | +382.8% | +16.23% |
| .300 | .10270 | .39287 | +335.7% | +13.89% |
| .6882487 | .06195 | .06230 | +.5853% | +.01508% |

At .05 GeV, elastic scattering is a small addition relative to
annihilation in the energy response, yet supplies most of the momentum
response in this controlled vacuum comparison. That distinction
weakens at the warmer benchmark. These are local coefficients under
the same relative bulk and thermal forces, not an evolved relative
velocity, decoupling temperature, scalar drag from omitted channels,
or a cosmological transfer function.

Replay from the kinetics directory:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error charged_lepton_portal_combined_momentum_reference.py
```
[Full numerical receipt](results/charged-lepton-portal-combined-momentum-2026-11-27.json).
