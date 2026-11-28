# Direct light fermion source: dark FD blocking and inverse reaction

The existing production model is a vacuum-dark final-state rate. This
separate benchmark checks that approximation for the direct light fermion
pair channel at the 45 MeV source-only snapshot from 50 MeV reheating.
It does not change the FRW production operator.

For a visible equilibrium FD bath at zero lepton chemical potential, the
product of visible inverse Pauli factors equals the product of the
forward FD factors times exp(E_pair/T). At each outgoing angle the
light dark final-state factors therefore multiply the existing positive
forward weight by B = (1-f1)(1-f2), with reverse weight B exp(A), where
A = E_pair/T - theta1 - theta2 and theta = log((1-f)/f).
This shared geometry admits an expm1 evaluation of the net rate without
subtracting nearly equal integrated terms. The number and total-energy
moments use the same measure and respect partner exchange.

Two quadrature settings recover the vacuum source for an empty dark
sector, net production consistent with zero at equal visible-temperature
zero-chemical-potential FD equilibrium, and net removal at positive dark
chemical potential with the same temperature. In the actual collisionless
snapshot proxy the largest stable fermion bin occupation is 4.73e-4.
Both its binned shape and the independently moment-matched FD shape give
a forward correction of about 6.5e-6 of the vacuum count rate; the inverse
number correction is about 9e-9. The net energy correction is about 5.1e-6.
These are small effects in this particular direct-channel snapshot, not a
bound for the evolved distribution, heavy-state production, scalar
stimulation or later reheating histories.

From `kinetics`, use the same external central EOS and raw source spectrum:

```sh
OPENBLAS_NUM_THREADS=1 python -W error joint_visible_direct_detailed_balance_reference.py EOS.dat empty.json
```
